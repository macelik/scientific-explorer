// Merge rounds: the working segmentation recipe, the call-noise table and their exports.
// Spec: docs/superpowers/specs/2026-10-27-merge-rounds-design.md
import type { ClusterBreakpoint, ClusterSegment } from './types'
import { tsv, type KaryoMode, type RecallPair } from './cnRecall'

export interface WireRound { rules: unknown[]; allow_gap_crossing: boolean }
export interface ClusterNoise { segments: number; zigzags: number; zigzag_median_bins: number | null; gain_loss_adjacent: number }
export interface CellsNoise { n_cells: number; segments_mean: number; zigzags_mean: number; zigzag_median_bins: number | null; gain_loss_adjacent_mean: number }
export type Noise = Record<string, Record<string, { cluster: ClusterNoise | null; cells: CellsNoise | null }>>
export interface Round extends WireRound { label: string; merges: number; noise: Noise | null }
export interface BaseResponse {
  hash: string; n_rounds: number; breakpoints: ClusterBreakpoint[]; rounds: { pairs: RecallPair[]; merges: number }[]
  cluster_segments: ClusterSegment[]; noise: Noise; noise_original: Noise
  percell: { status: 'original' | 'ready' | 'running' | 'missing'; job: string | null }
}
export interface NoiseColumn { name: string; noise: Noise | null; cellsPending?: boolean }

export const wireRounds = (rs: WireRound[]): WireRound[] => rs.map((r) => ({ rules: r.rules, allow_gap_crossing: r.allow_gap_crossing }))
export const sameRounds = (a: WireRound[], b: WireRound[]) => JSON.stringify(wireRounds(a)) === JSON.stringify(wireRounds(b))

export function roundLabel(rules: any[]): string {
  const r0 = rules[0] ?? {}
  const users = rules.slice(1).filter((r) => r && r.enabled).length
  return `τ ${r0.tau} · ref ${r0.srdRef} · ${r0.combine}${users ? ` + ${users} user rule${users > 1 ? 's' : ''}` : ''}`
}
export const roundStripText = (rounds: Round[]) => ['Original', ...rounds.map((r, i) => `Round ${i + 1} (${r.label}, ${r.merges} merges)`)].join(' → ')
export const baseHeading = (n: number) => (n ? `Base: round ${n} (${n} round${n > 1 ? 's' : ''} of merges on the original)` : 'Base: original')

const MODE_LABELS: Record<KaryoMode, string> = { original: 'Original', before: 'Before', after: 'After', changed: 'Changed' }
export function karyoModes(nRounds: number, hasAfter: boolean): { mode: KaryoMode; label: string; enabled: boolean }[] {
  const modes: KaryoMode[] = nRounds ? ['original', 'before', 'after', 'changed'] : ['before', 'after', 'changed']
  return modes.map((mode) => ({ mode, label: MODE_LABELS[mode], enabled: mode === 'original' || mode === 'before' || hasAfter }))
}

const MEASURES: { level: 'cluster' | 'cells'; key: string; label: string; kind: 'count' | 'median' | 'mean' }[] = [
  { level: 'cluster', key: 'segments', label: 'segments', kind: 'count' },
  { level: 'cluster', key: 'zigzags', label: 'zig-zags', kind: 'count' },
  { level: 'cluster', key: 'zigzag_median_bins', label: 'zig-zag median length (bins)', kind: 'median' },
  { level: 'cluster', key: 'gain_loss_adjacent', label: 'gain next to loss', kind: 'count' },
  { level: 'cells', key: 'segments_mean', label: 'segments per cell', kind: 'mean' },
  { level: 'cells', key: 'zigzags_mean', label: 'zig-zags per cell', kind: 'mean' },
  { level: 'cells', key: 'zigzag_median_bins', label: 'zig-zag median length (bins)', kind: 'median' },
  { level: 'cells', key: 'gain_loss_adjacent_mean', label: 'gain next to loss per cell', kind: 'mean' },
]
const fmtValue = (v: number | null | undefined, kind: string) =>
  v === null || v === undefined ? '—' : kind === 'count' ? String(v) : kind === 'median' ? (Number.isInteger(v) ? String(v) : v.toFixed(1)) : v.toFixed(1)
const entry = (n: Noise | null, source: string, key: string) => n?.[source]?.[key]

export function noiseRows(cols: NoiseColumn[], source: string, key: string) {
  return MEASURES.map((m) => ({
    level: m.level, measure: m.label,
    values: cols.map((c) => {
      const e = entry(c.noise, source, key)
      const lv = m.level === 'cluster' ? e?.cluster : e?.cells
      if (!lv && m.level === 'cells' && c.cellsPending) return 're-calling…'
      return fmtValue(lv ? (lv as any)[m.key] : null, m.kind)
    }),
  }))
}
export function noiseCellsN(cols: NoiseColumn[], source: string, key: string): number | null {
  for (const c of cols) { const n = entry(c.noise, source, key)?.cells?.n_cells; if (n !== undefined) return n }
  return null
}
export const roundsHeader = (rounds: WireRound[]) => rounds.map((r, i) => `round ${i + 1} rules: ${JSON.stringify(r.rules)} · allow_gap_crossing: ${r.allow_gap_crossing}`)
export function noiseTsv(cols: NoiseColumn[], sources: string[], header: string[]): string {
  const rows: Record<string, unknown>[] = []
  for (const source of sources) {
    const keys: string[] = []
    for (const c of cols) for (const k of Object.keys(c.noise?.[source] ?? {})) if (!keys.includes(k)) keys.push(k)
    for (const key of keys) for (const m of MEASURES) {
      const r: Record<string, unknown> = { source, chromosome: key, level: m.level, measure: m.label }
      for (const c of cols) { const e = entry(c.noise, source, key); const lv = m.level === 'cluster' ? e?.cluster : e?.cells; r[c.name] = lv ? (lv as any)[m.key] : null }
      rows.push(r)
    }
  }
  return tsv(header, ['source', 'chromosome', 'level', 'measure', ...cols.map((c) => c.name)], rows)
}
export const roundsFromSession = (p: any): Round[] =>
  (Array.isArray(p?.rounds) ? p.rounds : [])
    .filter((r: any) => r && typeof r === 'object' && Array.isArray(r.rules) && typeof r.allow_gap_crossing === 'boolean')
    .map((r: any) => ({ ...r, label: typeof r.label === 'string' ? r.label : roundLabel(r.rules), merges: typeof r.merges === 'number' ? r.merges : 0, noise: r.noise ?? null }))

export const NOISE_DEFINITIONS = [
  'Calls are read per chromosome as the sequence of segment states. Segments with no call are skipped. Gain = state > 1 (weak gain included), loss = state < 1.',
  'Zig-zag: an internal segment whose state differs from both neighbours, while the two neighbours have the same state as each other.',
  'Gain next to loss: a boundary with a gain on one side and a loss on the other. This can be real (for example an unbalanced rearrangement), so it is counted but never merged automatically.',
  "Cell values: the mean per cell over all cells of the cluster. The zig-zag median length is taken over all of those cells' zig-zags. For all, the cells are every cell, each read under its own cluster's segmentation.",
]

async function postJson<T>(url: string, body: unknown): Promise<T> {
  const r = await fetch(url, { method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify(body) })
  if (!r.ok) { let msg = `${r.status}`; const status = r.status; try { const b = await r.json(); msg = typeof b.detail === 'string' ? b.detail : JSON.stringify(b.detail ?? b) } catch { /* keep status */ } throw Object.assign(new Error(msg), { status }) }
  return r.json() as Promise<T>
}
export const fetchBase = (rounds: WireRound[]) => postJson<BaseResponse>('/api/integration/base', { rounds })
export const pinRecallApi = (id: string, rounds: WireRound[]) => postJson<{ pinned: string }>(`/api/integration/recall/${encodeURIComponent(id)}/pin`, { rounds })
