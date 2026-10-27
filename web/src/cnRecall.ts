// Genome-wide merge + CN re-call: types, fetchers, karyogram pixel rule and TSV exports.
// Spec: docs/superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md
import type { ClusterSegment } from './types'
import { STORY_CN_LABELS } from './integration'

export const UNCHANGED = 252
export type KaryoMode = 'before' | 'after' | 'changed'
export interface RecallBody { rules: unknown[]; allow_gap_crossing: boolean }
export interface RecallPair { source: string; chrom: string; original: number; merged: number; merges: number; phase2: number }
export interface RecallSourceDiff { source: string; segments_before: number; segments_after: number; merges: number; phase2: number; best_s_before: number | null; best_s_after: number | null; bins_changed: number; bins_called: number; pct_changed: number }
export interface RecallTransition { before: number | null; after: number | null; bins: number }
export interface RecallCellRow { cell: string; wnn_cluster: number; best_s_before: number | null; best_s_after: number | null; bins_changed: number }
export interface RecallCells { n_cells: number; n_changed_cells: number; median_pct_changed: number; best_s_change: { median: number | null; min: number | null; max: number | null }; per_cell: RecallCellRow[] }
export interface RecallStatus {
  id: string; hash: string; status: 'running' | 'done' | 'error'; stage: 'merge' | 'cluster' | 'percell' | 'done'
  progress: { done: number; total: number } | null; error: string | null
  settings: { rules: unknown[]; allow_gap_crossing: boolean }; provenance: Record<string, unknown>
  merge: { pairs: RecallPair[]; merged_breakpoints: Record<string, unknown>[] } | null
  cluster_segments: ClusterSegment[] | null
  diff: { sources: RecallSourceDiff[]; cells: RecallCells | null; transitions: { clusters: RecallTransition[]; cells: RecallTransition[] | null } } | null
  percell_ready: boolean
}

async function asJson<T>(r: Response): Promise<T> {
  if (!r.ok) { let msg = `${r.status}`; try { const b = await r.json(); msg = typeof b.detail === 'string' ? b.detail : JSON.stringify(b.detail ?? b) } catch { /* keep status */ } throw new Error(msg) }
  return r.json() as Promise<T>
}
export const startRecallApi = async (body: RecallBody) =>
  asJson<{ id: string; hash: string; status: string }>(await fetch('/api/integration/recall', { method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify(body) }))
export const fetchRecall = async (id: string) => asJson<RecallStatus>(await fetch(`/api/integration/recall/${encodeURIComponent(id)}`))
export async function fetchRecallPercell(id: string): Promise<Uint8Array> {
  const r = await fetch(`/api/integration/recall/${encodeURIComponent(id)}/percell.bin`)
  if (!r.ok) throw new Error(`${r.status}`)
  return new Uint8Array(await r.arrayBuffer())
}

/** Pixel code: with a `before` value (Changed mode), a called bin whose state did not change is drawn as UNCHANGED. */
export function karyoCode(after: number, before?: number): number {
  return before !== undefined && before === after && after !== 255 ? UNCHANGED : after
}
export function clusterStateRow(segs: ClusterSegment[], group: string, nBins: number): Uint8Array {
  const a = new Uint8Array(nBins).fill(255)
  for (const g of segs) if (g.group === group) a.fill(Math.round(g.cn_call * 2), g.bin_start, g.bin_end + 1)
  return a
}
const stateText = (s: number | null) => (s === null ? 'not called' : `${s} (${STORY_CN_LABELS[s]})`)
export const stateLabel = (code: number) => stateText(code === 255 ? null : code / 2)
export const transitionText = (t: RecallTransition) => `${stateText(t.before)} → ${stateText(t.after)}: ${t.bins} bins`

const n = (v: number) => v.toLocaleString('en-US')
export function recallProgressText(st: RecallStatus | null): string {
  if (!st) return ''
  if (st.status === 'error') return `error: ${st.error}`
  const p = st.progress ? ` ${n(st.progress.done)} / ${n(st.progress.total)}` : ''
  if (st.status === 'done') return 'merging ✓ → cluster CN ✓ → per-cell ✓ · done'
  if (st.stage === 'merge') return `merging…${p}`
  if (st.stage === 'cluster') return `merging ✓ → cluster CN…${p}`
  return `merging ✓ → cluster CN ✓ → per-cell${p}`
}

export function tsv(header: string[], columns: string[], rows: Record<string, unknown>[]): string {
  const cell = (v: unknown) => (v === null || v === undefined ? '' : String(v))
  return [...header.map((h) => `# ${h}`), columns.join('\t'), ...rows.map((r) => columns.map((c) => cell(r[c])).join('\t'))].join('\n') + '\n'
}
export const recallHeader = (st: RecallStatus) => [`rules: ${JSON.stringify(st.settings.rules)}`, `allow_gap_crossing: ${st.settings.allow_gap_crossing}`, `provenance: ${JSON.stringify(st.provenance)}`]
const BREAKPOINT_COLUMNS = ['source', 'chromosome', 'absolute_bin', 'bp_genomic_start', 'bp_genomic_end', 'recursion_depth', 'ad_score', 'p_local', 'p_global', 'n_cells']
export const CLUSTER_SEGMENT_COLUMNS = ['group', 'n_cells', 'best_s', 'segment', 'chromosome', 'bin_start', 'bin_end', 'n_bins', 'continuous_cn', 'integer_cn', 'watson', 'cn_call']
export function mergedBreakpointsTsv(st: RecallStatus): string {
  const rows = st.merge?.merged_breakpoints ?? []
  return tsv(recallHeader(st), rows.length ? Object.keys(rows[0]) : BREAKPOINT_COLUMNS, rows)
}
export const recalledClusterSegmentsTsv = (st: RecallStatus) => tsv(recallHeader(st), CLUSTER_SEGMENT_COLUMNS, (st.cluster_segments ?? []) as unknown as Record<string, unknown>[])
export const recalledPercellTsv = (st: RecallStatus) => tsv(recallHeader(st), ['cell', 'wnn_cluster', 'best_s_before', 'best_s_after', 'bins_changed'], (st.diff?.cells?.per_cell ?? []) as unknown as Record<string, unknown>[])
export function changeSummaryTsv(st: RecallStatus): string {
  const cols = ['source', 'segments_before', 'segments_after', 'merges', 'phase2', 'best_s_before', 'best_s_after', 'bins_changed', 'bins_called', 'pct_changed']
  const tr = [...(st.diff?.transitions.clusters ?? []).map((t) => ({ scope: 'clusters', ...t })), ...(st.diff?.transitions.cells ?? []).map((t) => ({ scope: 'cells', ...t }))]
  return tsv(recallHeader(st), cols, (st.diff?.sources ?? []) as unknown as Record<string, unknown>[]) + '\n' + tsv([], ['scope', 'before', 'after', 'bins'], tr)
}
export function downloadText(text: string, name: string) {
  const url = URL.createObjectURL(new Blob([text], { type: 'text/tab-separated-values' })); const a = document.createElement('a')
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000)
}
