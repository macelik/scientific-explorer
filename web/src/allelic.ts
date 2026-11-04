// Allelic segmentation: types, API, settings, comparison tables.
// Spec: docs/superpowers/specs/2026-10-28-allelic-segmentation-design.md
import type { ClusterBreakpoint } from './types'

export type SignalId = 'dbaf' | 'abs_dbaf' | 'shrunk_abs_dbaf' | 'oriented_all'
export const SIGNALS: { id: SignalId; label: string; signed: boolean; formula: string }[] = [
  { id: 'dbaf', label: 'signed dBAF', signed: true, formula: '(HB − HA) / (2·AT)' },
  { id: 'abs_dbaf', label: '|dBAF|', signed: false, formula: '|HB − HA| / (2·AT)' },
  { id: 'shrunk_abs_dbaf', label: 'shrunk |dBAF|', signed: false, formula: '|HB − HA| / (2·(AT + 2))' },
  { id: 'oriented_all', label: 'oriented to all cells', signed: true, formula: 'sign(ΣHB_all − ΣHA_all) · (HB − HA) / (2·AT)' },
]
export const signalLabel = (id: string) => SIGNALS.find((s) => s.id === id)?.label ?? id

export type GlobalPool = 'genome' | 'chromosome'
export interface AdSettings { min_at: number; k: number; n_permutations: number; alpha: number; min_child_bins: number; min_terminal_bins: number; global_pool: GlobalPool }
export const DEFAULT_SETTINGS: AdSettings = { min_at: 10, k: 6, n_permutations: 1000, alpha: 0.001, min_child_bins: 20, min_terminal_bins: 50, global_pool: 'genome' }
export const POOL_NOTE = "The app's depth breakpoints now use the genome pool, as original pyEpi does. The chromosome pool is how the earlier depth run (raw_accepted_breakpoints.tsv in integration-story/data) sampled the global test."
export const SETTING_FIELDS: { key: Exclude<keyof AdSettings, 'global_pool'>; label: string; min: number; max: number; integer: boolean; minExclusive?: boolean }[] = [
  { key: 'min_at', label: 'min AT', min: 1, max: 1000, integer: true },
  { key: 'k', label: 'k (max depth)', min: 1, max: 8, integer: true },
  { key: 'n_permutations', label: 'permutations', min: 100, max: 5000, integer: true },
  { key: 'alpha', label: 'alpha', min: 0, max: 0.1, integer: false, minExclusive: true },
  { key: 'min_child_bins', label: 'min child (kept bins)', min: 1, max: 500, integer: true },
  { key: 'min_terminal_bins', label: 'min from chromosome end (kept bins)', min: 0, max: 500, integer: true },
]
export function settingsProblems(s: AdSettings): string[] {
  const out: string[] = []
  for (const f of SETTING_FIELDS) {
    const v = (s as any)[f.key]
    if (typeof v !== 'number' || !Number.isFinite(v)) { out.push(`${f.key}: not a number`); continue }
    if (f.integer && !Number.isInteger(v)) out.push(`${f.key}: must be a whole number`)
    if ((f.minExclusive ? v <= f.min : v < f.min) || v > f.max) out.push(`${f.key}: must be ${f.minExclusive ? '>' : '≥'} ${f.min} and ≤ ${f.max}`)
  }
  if (s.global_pool !== 'genome' && s.global_pool !== 'chromosome') out.push('global_pool: must be genome or chromosome')
  return out
}
export const sameRun = (run: { signal: string; settings: AdSettings } | null, signal: string, s: AdSettings) =>
  !!run && run.signal === signal && run.settings.global_pool === s.global_pool && SETTING_FIELDS.every((f) => run.settings[f.key] === s[f.key])

export interface AllelicStatus { available: boolean; state: 'idle' | 'building' | 'ready' | 'error'; problems: string[]; notes: string[]; allelic_dir: string; sources: string[] }
export interface AllelicTracks { chrom: string; signal: SignalId; min_at: number; n: number; tracks: { source: string; values: (number | null)[]; at: number[] }[] }
export interface AllelicBreakpoint { source: string; chromosome: string; absolute_bin: number; bp_genomic_start: number; bp_genomic_end: number; gap_start: number | null; gap_end: number | null; recursion_depth: number; ad_score: number; p_local: number; p_global: number }
export interface AllelicSegment { source: string; chromosome: string; bin_start: number; bin_end: number; n_bins: number; kept_bins: number; mean: number | null; median: number | null; at_sum: number; ha_sum: number; hb_sum: number }
export interface AllelicJob {
  id: string; status: 'queued' | 'running' | 'done' | 'error'; error: string | null; signal: SignalId; settings: AdSettings
  progress: { done: number; total: number; source: string | null; chromosome: string | null } | null
  result: { breakpoints: AllelicBreakpoint[]; audit: any[]; segments: AllelicSegment[] } | null
}

async function json<T>(r: Response): Promise<T> {
  if (!r.ok) { let msg = `${r.status}`; try { const b = await r.json(); msg = typeof b.detail === 'string' ? b.detail : JSON.stringify(b.detail ?? b) } catch { /* keep status */ } throw Object.assign(new Error(msg), { status: r.status }) }
  return r.json() as Promise<T>
}
const enc = encodeURIComponent
export const allelicApi = {
  status: () => fetch('/api/allelic/status').then((r) => json<AllelicStatus>(r)),
  tracks: (chrom: string, sources: string[], signal: SignalId, minAt: number) =>
    fetch(`/api/allelic/tracks?chrom=${enc(chrom)}&sources=${enc(sources.join(','))}&signal=${signal}&min_at=${minAt}`).then((r) => json<AllelicTracks>(r)),
  start: (signal: SignalId, settings: AdSettings) =>
    fetch('/api/allelic/segment', { method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify({ signal, settings }) }).then((r) => json<AllelicJob>(r)),
  job: (id: string) => fetch(`/api/allelic/segment/${enc(id)}`).then((r) => json<AllelicJob>(r)),
}

// ---------------- comparison (spec §5); all bins chromosome-local
export interface AllelicBp { bin: number; gap_start: number | null; gap_end: number | null }
export function bpDistance(depthBin: number, a: AllelicBp): number {
  if (a.gap_start !== null && a.gap_end !== null) {
    if (depthBin >= a.gap_start && depthBin <= a.gap_end) return 0
    return depthBin < a.gap_start ? a.gap_start - depthBin : depthBin - a.gap_end
  }
  return Math.abs(depthBin - a.bin)
}
export function matchBreakpoints(depth: number[], allelic: AllelicBp[], N: number) {
  const cand: { d: number; a: number; dist: number }[] = []
  depth.forEach((db, i) => allelic.forEach((ab, j) => { const dist = bpDistance(db, ab); if (dist <= N) cand.push({ d: i, a: j, dist }) }))
  cand.sort((x, y) => x.dist - y.dist || depth[x.d] - depth[y.d] || allelic[x.a].bin - allelic[y.a].bin)
  const usedD = new Set<number>(), usedA = new Set<number>(), pairs: typeof cand = []
  for (const c of cand) if (!usedD.has(c.d) && !usedA.has(c.a)) { usedD.add(c.d); usedA.add(c.a); pairs.push(c) }
  return { pairs, depthOnly: depth.map((_, i) => i).filter((i) => !usedD.has(i)), allelicOnly: allelic.map((_, j) => j).filter((j) => !usedA.has(j)) }
}
const median = (v: number[]) => { const s = [...v].sort((a, b) => a - b), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2 }

export interface AgreementRow { source: string; depth: number; allelic: number; matched: number; depthOnly: number; allelicOnly: number; medianOffset: number | null }
export function agreementRows(sources: string[], chroms: string[], depthBps: ClusterBreakpoint[], allelicBps: AllelicBreakpoint[], N: number): AgreementRow[] {
  return sources.map((source) => {
    const row: AgreementRow = { source, depth: 0, allelic: 0, matched: 0, depthOnly: 0, allelicOnly: 0, medianOffset: null }
    const offs: number[] = []
    for (const chrom of chroms) {
      const d = depthBps.filter((b) => b.source === source && b.chromosome === chrom).map((b) => b.absolute_bin)
      const a = allelicBps.filter((b) => b.source === source && b.chromosome === chrom).map((b) => ({ bin: b.absolute_bin, gap_start: b.gap_start, gap_end: b.gap_end }))
      const m = matchBreakpoints(d, a, N)
      row.depth += d.length; row.allelic += a.length; row.matched += m.pairs.length; row.depthOnly += m.depthOnly.length; row.allelicOnly += m.allelicOnly.length
      offs.push(...m.pairs.map((p) => p.dist))
    }
    row.medianOffset = offs.length ? median(offs) : null
    return row
  })
}

/** Depth cluster segment, chromosome-local half-open [s, e) as returned by clusterSegmentsLocal. */
export interface LocalDepthSeg { group: string; s: number; e: number; cn_call: number }
export interface DepthSegRow { source: string; bin_start: number; bin_end: number; cn_call: number; kept: number; mean: number | null; sd: number | null; allelicInside: number }
const round6 = (x: number) => Math.round(x * 1e6) / 1e6
export function depthSegmentRows(segs: LocalDepthSeg[], values: Record<string, (number | null)[]>, allelicBps: AllelicBreakpoint[], chrom: string, N: number): DepthSegRow[] {
  return segs.filter((g) => values[g.group]).map((g) => {
    const v = values[g.group].slice(g.s, g.e).filter((x): x is number => x !== null && Number.isFinite(x))
    const mean = v.length ? v.reduce((a, b) => a + b, 0) / v.length : null
    const sd = mean === null ? null : Math.sqrt(v.reduce((a, b) => a + (b - mean) ** 2, 0) / v.length)
    const inside = allelicBps.filter((b) => b.source === g.group && b.chromosome === chrom && b.absolute_bin - g.s > N && g.e - b.absolute_bin > N).length
    return { source: g.group, bin_start: g.s, bin_end: g.e - 1, cn_call: g.cn_call, kept: v.length, mean: mean === null ? null : round6(mean), sd: sd === null ? null : round6(sd), allelicInside: inside }
  })
}
export function callShares(seg: { bin_start: number; bin_end: number }, depthSegs: LocalDepthSeg[]): { call: number; frac: number }[] {
  const n = seg.bin_end - seg.bin_start + 1, acc = new Map<number, number>()
  for (const d of depthSegs) { const lo = Math.max(d.s, seg.bin_start), hi = Math.min(d.e - 1, seg.bin_end); if (hi >= lo) acc.set(d.cn_call, (acc.get(d.cn_call) ?? 0) + hi - lo + 1) }
  return [...acc].map(([call, k]) => ({ call, frac: k / n })).sort((a, b) => b.frac - a.frac || a.call - b.call)
}
export const sharesText = (sh: { call: number; frac: number }[]) => (sh.length ? sh.map((x) => `${x.call.toFixed(1)}: ${Math.round(x.frac * 100)}%`).join(' · ') : 'no depth call')
export interface AllelicSegRow { source: string; bin_start: number; bin_end: number; kept: number; mean: number | null; spanned: string }
export function allelicSegmentRows(segs: AllelicSegment[], depthSegs: LocalDepthSeg[], sources: string[], chrom: string): AllelicSegRow[] {
  return segs.filter((s) => s.chromosome === chrom && sources.includes(s.source)).map((s) => ({
    source: s.source, bin_start: s.bin_start, bin_end: s.bin_end, kept: s.kept_bins, mean: s.mean,
    spanned: sharesText(callShares(s, depthSegs.filter((d) => d.group === s.source))),
  }))
}

export const COMPARISON_DEFINITIONS = [
  'Matching is one-to-one, nearest first: all depth–allelic pairs within N bins are sorted by distance and accepted greedily, skipping any breakpoint already matched; ties go to the lower depth bin.',
  'An allelic breakpoint inside a gap of skipped bins is at distance 0 from any bin of that gap; outside it, the distance is to the nearer gap edge.',
  'Distances count genomic bins, not kept bins. Median offset is taken over matched pairs.',
  'Depth segments: kept bins have AT ≥ min AT; mean and sd (population) of the allelic signal over kept bins; "allelic bps inside" counts allelic breakpoints more than N bins from both segment edges.',
  'Allelic segments: each depth CN call with its share of the segment\'s bins, largest first.',
]
const num = (v: number | null, d = 4) => (v === null || !Number.isFinite(v) ? 'NA' : Number.isInteger(v) ? String(v) : v.toFixed(d))
export function runHeader(job: AllelicJob, N: number, depthSource: 'current' | 'original', nRounds: number, scope: string, identity: string | null | undefined): string[] {
  return [`signal: ${job.signal} (${signalLabel(job.signal)})`, `settings: ${JSON.stringify(job.settings)}`, `match N: ${N} bins`,
    `depth segments: ${depthSource === 'current' ? `current round (${nRounds} merge rounds)` : 'unmerged (no merge rounds)'}`, `depth identity: ${identity ?? 'unknown'}`, `table 1 scope: ${scope}`]
}
export function comparisonTsv(header: string[], t1: AgreementRow[], t2: DepthSegRow[], t3: AllelicSegRow[]): string {
  const lines = header.map((h) => `# ${h}`)
  lines.push('# table 1: breakpoint agreement', ['source', 'depth_bps', 'allelic_bps', 'matched', 'depth_only', 'allelic_only', 'median_offset'].join('\t'))
  for (const r of t1) lines.push([r.source, r.depth, r.allelic, r.matched, r.depthOnly, r.allelicOnly, num(r.medianOffset, 1)].join('\t'))
  lines.push('# table 2: depth segments seen by the allelic signal', ['source', 'bin_start', 'bin_end', 'cn_call', 'kept_bins', 'mean', 'sd', 'allelic_bps_inside'].join('\t'))
  for (const r of t2) lines.push([r.source, r.bin_start, r.bin_end, r.cn_call, r.kept, num(r.mean), num(r.sd), r.allelicInside].join('\t'))
  lines.push('# table 3: allelic segments seen by depth', ['source', 'bin_start', 'bin_end', 'kept_bins', 'mean', 'depth_calls_spanned'].join('\t'))
  for (const r of t3) lines.push([r.source, r.bin_start, r.bin_end, r.kept, num(r.mean), r.spanned].join('\t'))
  return lines.join('\n') + '\n'
}

// ---------------- session (spec §4)
export interface AllelicView { signal: SignalId; settings: AdSettings; chrom: string; sources: string[]; depthSource: 'current' | 'original'; matchN: number; scope: 'chrom' | 'genome' }
export const DEFAULT_VIEW: AllelicView = { signal: 'abs_dbaf', settings: DEFAULT_SETTINGS, chrom: 'chr8', sources: ['cluster0', 'cluster1', 'cluster2', 'cluster3', 'cluster4'], depthSource: 'current', matchN: 3, scope: 'chrom' }
export function allelicFromSession(p: any): { view: AllelicView; notice: string | null } {
  if (!p || typeof p !== 'object') return { view: DEFAULT_VIEW, notice: null }
  const bad: string[] = []
  const signal = SIGNALS.some((s) => s.id === p.signal) ? p.signal : (bad.push(`signal ${JSON.stringify(p.signal)}`), DEFAULT_VIEW.signal)
  const settings = { ...DEFAULT_SETTINGS, ...(p.settings && typeof p.settings === 'object' ? p.settings : {}) }
  const sp = settingsProblems(settings); if (sp.length) bad.push(...sp)
  const chrom = typeof p.chrom === 'string' ? p.chrom : (bad.push('chromosome'), DEFAULT_VIEW.chrom)
  const sources = Array.isArray(p.sources) && p.sources.every((s: unknown) => typeof s === 'string') ? p.sources : (bad.push('sources'), DEFAULT_VIEW.sources)
  const depthSource = p.depthSource === 'original' || p.depthSource === 'current' ? p.depthSource : (bad.push('depth segments'), DEFAULT_VIEW.depthSource)
  const matchN = Number.isInteger(p.matchN) && p.matchN >= 0 && p.matchN <= 50 ? p.matchN : (bad.push('match N'), DEFAULT_VIEW.matchN)
  const scope = p.scope === 'chrom' || p.scope === 'genome' ? p.scope : (bad.push('scope'), DEFAULT_VIEW.scope)
  if (bad.length) return { view: DEFAULT_VIEW, notice: `Saved allelic settings could not be applied: ${bad.join('; ')}` }
  return { view: { signal, settings: sp.length ? DEFAULT_SETTINGS : settings, chrom, sources, depthSource, matchN, scope }, notice: null }
}
