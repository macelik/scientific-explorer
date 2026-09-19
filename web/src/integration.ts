/** Pure helpers for the integration-story tabs (no store access). */
import type { ChromMeta, ClusterBreakpoint, ClusterSegment, FlankRow, IntegrationCell, ScreeningRow } from './types'

export const SOURCE_COLORS: Record<string, string> = { cluster0: '#1f77b4', cluster1: '#ff7f0e', cluster2: '#2ca02c', cluster3: '#d62728', cluster4: '#9467bd', all: '#4b5563' }
export const TAB10 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173', '#3182bd', '#e6550d']
export const STORY_CN: Record<number, string> = { 0: '#7B1FA2', 0.5: '#CE93D8', 1: '#4CAF50', 1.5: '#FF8A65', 2: '#AB0000' }
export const STORY_CN_LABELS: Record<number, string> = { 0: 'loss', 0.5: 'weak loss', 1: 'baseline', 1.5: 'weak gain', 2: 'gain' }
export const VIRIDIS = ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725']
export const DEPTH_COLORS = ['#440154', '#414487', '#2a788e', '#22a884', '#7ad151', '#fde725']
export const NOISE_CV_BINS = [0.6, 0.8, 1.0, 1.2]

export type Metric = 'delta_median' | 'delta_mean' | 'delta_trim' | 'log2fc' | 'srd' | 'srd_phi'
export type FlankRule = 'min_delta' | 'min_srd' | 'max_delta' | 'left' | 'right'
export const METRICS: { id: Metric; label: string; kind: 'effect' | 'stat'; bins: number[]; binLabels: string[]; source: string }[] = [
  { id: 'delta_median', label: 'log2 effect · median levels (delta)', kind: 'effect', bins: [0.2, 0.5, 0.75, 1.0], binLabels: ['< 0.20', '0.20–0.50', '0.50–0.75', '0.75–1.00', '≥ 1.00'], source: 'flank scores: delta_left/right (median per-bin levels)' },
  { id: 'delta_mean', label: 'log2 effect · mean levels', kind: 'effect', bins: [0.2, 0.5, 0.75, 1.0], binLabels: ['< 0.20', '0.20–0.50', '0.50–0.75', '0.75–1.00', '≥ 1.00'], source: 'flank scores: delta_mean_left/right' },
  { id: 'delta_trim', label: 'log2 effect · IQR-trimmed mean levels', kind: 'effect', bins: [0.2, 0.5, 0.75, 1.0], binLabels: ['< 0.20', '0.20–0.50', '0.50–0.75', '0.75–1.00', '≥ 1.00'], source: 'flank scores: delta_trim_left/right (levels trimmed to [Q1−1.5·IQR, Q3+1.5·IQR])' },
  { id: 'log2fc', label: 'log2 fold change · mean rate (SRD screen)', kind: 'effect', bins: [0.2, 0.5, 0.75, 1.0], binLabels: ['< 0.20', '0.20–0.50', '0.50–0.75', '0.75–1.00', '≥ 1.00'], source: 'screening: log2fc_left/right (summed-depth rate ratio)' },
  { id: 'srd', label: '|SRD| · signed root deviance (descriptive)', kind: 'stat', bins: [1.5, 3.3, 5.0, 8.0], binLabels: ['< 1.5', '1.5–3.3', '3.3–5.0', '5.0–8.0', '≥ 8.0'], source: 'screening: srd_left/right (Poisson-deviance contrast on summed depth)' },
  { id: 'srd_phi', label: '|SRD/√φ| · dispersion-scaled', kind: 'stat', bins: [1.5, 3.3, 5.0], binLabels: ['< 1.5', '1.5–3.3', '3.3–5.0', '≥ 5.0'], source: 'flank scores: srd_phi_left/right (pooled Pearson dispersion φ per cluster/chromosome)' },
]

export interface SmallSeg {
  key: string; source: string; chrom: string; segId: number; s: number; e: number; len: number; mbStart: number; mbEnd: number
  flank: FlankRow; screen: ScreeningRow | null
  left: Record<Metric, number | null>; right: Record<Metric, number | null>
  stats: { median: number; mean: number; trim: number; std: number; cv: number; flankStdL: number | null; flankStdR: number | null; flankMeanL: number | null; flankMeanR: number | null; flankCvL: number | null; flankCvR: number | null }
  reasonL: string; reasonR: string; proposal: string | null
}

const n = (v: any): number | null => (v === null || v === undefined || Number.isNaN(+v) ? null : +v)

/** Join the flank-score table with the SRD-screen table by (source, chromosome, bin_start). */
export function buildSmallSegments(flank: FlankRow[] | null, screen: ScreeningRow[] | null): SmallSeg[] {
  if (!flank) return []
  const sc = new Map<string, ScreeningRow>()
  for (const r of screen || []) sc.set(`${r.source}|${r.chromosome}|${r.bin_start}`, r)
  return flank.map((f) => {
    const key = `${f.source}|${f.chromosome}|${f.bin_start}`
    const s = sc.get(key) || null
    const side = (w: 'left' | 'right'): Record<Metric, number | null> => ({
      delta_median: n(f[`delta_${w}`]), delta_mean: n(f[`delta_mean_${w}`]), delta_trim: n(f[`delta_trim_${w}`]),
      log2fc: s ? n(s[`log2fc_${w}`]) : null, srd: s ? n(s[`srd_${w}`]) : null, srd_phi: n(f[`srd_phi_${w}`]),
    })
    return {
      key, source: f.source, chrom: f.chromosome, segId: f.orig_seg_id, s: f.bin_start, e: f.bin_end, len: f.seg_length, mbStart: f.Mb_start, mbEnd: f.Mb_end,
      flank: f, screen: s, left: side('left'), right: side('right'),
      stats: { median: f.median_seg, mean: f.mean_seg_perbin, trim: f.trimmean_seg_perbin, std: f.std_seg_perbin, cv: f.cv_seg_perbin,
        flankStdL: n(f.std_left_flank_perbin), flankStdR: n(f.std_right_flank_perbin), flankMeanL: n(f.mean_left_flank_perbin), flankMeanR: n(f.mean_right_flank_perbin), flankCvL: n(f.cv_left_flank_perbin), flankCvR: n(f.cv_right_flank_perbin) },
      reasonL: f.reason_left, reasonR: f.reason_right, proposal: s ? s.proposal : null,
    }
  })
}

/** Which flank the rule picks, given the effect metric used for "min/max delta". */
export function pickFlank(seg: SmallSeg, rule: FlankRule, effectMetric: Metric): 'left' | 'right' | null {
  const dl = seg.left[effectMetric], dr = seg.right[effectMetric]
  const sl = seg.left.srd_phi, sr = seg.right.srd_phi
  const abs = (v: number | null) => (v === null ? Infinity : Math.abs(v))
  if (rule === 'left') return dl === null && sl === null ? null : 'left'
  if (rule === 'right') return dr === null && sr === null ? null : 'right'
  if (rule === 'min_srd') { if (sl === null && sr === null) return null; return abs(sl) <= abs(sr) ? 'left' : 'right' }
  if (dl === null && dr === null) return null
  if (rule === 'max_delta') return dr === null || (dl !== null && Math.abs(dl) >= Math.abs(dr)) ? 'left' : 'right'
  return abs(dl) <= abs(dr) && dl !== null ? 'left' : 'right'
}
export function metricValue(seg: SmallSeg, metric: Metric, rule: FlankRule, effectMetric: Metric): { value: number | null; flank: 'left' | 'right' | null } {
  const f = pickFlank(seg, rule, effectMetric)
  if (!f) return { value: null, flank: null }
  const v = seg[f][metric]
  return { value: v === null ? null : (METRICS.find((m) => m.id === metric)!.kind === 'stat' ? Math.abs(v) : v), flank: f }
}
/** Consistency class used by the figures' hatching: SRD of the min-delta flank weak/strong, or the min-delta flank ≠ min-SRD flank. */
export function consistencyClass(seg: SmallSeg, effectMetric: Metric, srdRef = 3.3): { cls: 'weak' | 'strong' | 'inconsistent' | 'na'; label: string } {
  const fd = pickFlank(seg, 'min_delta', effectMetric), fs = pickFlank(seg, 'min_srd', effectMetric)
  if (!fd || !fs) return { cls: 'na', label: 'no flank metrics' }
  if (fd !== fs) return { cls: 'inconsistent', label: 'min-delta flank ≠ min-SRD flank (inconsistent, hatch /)' }
  const s = seg[fd].srd_phi
  if (s === null) return { cls: 'na', label: 'no SRD/√φ' }
  return Math.abs(s) < srdRef ? { cls: 'weak', label: `min-delta flank |SRD/√φ| ${Math.abs(s).toFixed(2)} < ${srdRef} (SRD weak, hatch x)` } : { cls: 'strong', label: `min-delta flank |SRD/√φ| ${Math.abs(s).toFixed(2)} ≥ ${srdRef} (SRD strong, hatch o)` }
}
/** Parallel diagonal segments of a "/" (or, mirrored, "\") hatch clipped to the unit square [0,1]x[0,1]. */
function unitDiag(step: number, back: boolean): { x0: number; y0: number; x1: number; y1: number }[] {
  const out: { x0: number; y0: number; x1: number; y1: number }[] = []
  const cMin = back ? 0 : -1, cMax = back ? 2 : 1
  for (let c = cMin + step / 2; c < cMax; c += step) {
    const x0 = back ? Math.max(0, c - 1) : Math.max(0, -c)
    const x1 = back ? Math.min(1, c) : Math.min(1, 1 - c)
    if (x1 - x0 > 1e-6) out.push({ x0, y0: back ? -x0 + c : x0 + c, x1, y1: back ? -x1 + c : x1 + c })
  }
  return out
}
/** Plotly `line` shapes hatching a rect (x0..x1 data units, y 0..1 of the given axis' domain): '/' single diagonal, 'x' crossed, 'none' nothing. */
export function hatchShapes(x0: number, x1: number, yref: string, style: 'none' | '/' | 'x', color = 'rgba(17,24,39,0.55)'): any[] {
  if (style === 'none') return []
  const step = 0.22, w = x1 - x0
  const lines = style === 'x' ? [...unitDiag(step, false), ...unitDiag(step, true)] : unitDiag(step, false)
  return lines.map((l) => ({ type: 'line', xref: 'x', yref: `${yref} domain`, x0: x0 + l.x0 * w, x1: x0 + l.x1 * w, y0: l.y0, y1: l.y1, line: { color, width: 1 }, layer: 'above' }))
}
export function binIndex(v: number, bins: number[]): number { let i = 0; while (i < bins.length && v >= bins[i]) i++; return i }
export function binColor(v: number | null, bins: number[]): string {
  if (v === null || !Number.isFinite(v)) return '#e5e7eb'
  const i = binIndex(Math.abs(v), bins); const pal = bins.length === 3 ? [VIRIDIS[0], VIRIDIS[1], VIRIDIS[2], VIRIDIS[4]] : VIRIDIS
  return pal[Math.min(i, pal.length - 1)]
}
export function viridisAt(t: number): string {
  const stops = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]]
  const x = Math.max(0, Math.min(1, t)) * (stops.length - 1); const i = Math.min(stops.length - 2, Math.floor(x)); const f = x - i
  const c = stops[i].map((a, k) => Math.round(a + (stops[i + 1][k] - a) * f)); return `rgb(${c[0]},${c[1]},${c[2]})`
}

/** pyEpi trimmed_mean_iqr(lb=False): keep [0, Q3+1.5·IQR], numpy linear percentiles. */
export function percentile(v: number[], q: number): number { const s = [...v].sort((a, b) => a - b); if (!s.length) return NaN; const p = (q / 100) * (s.length - 1); const lo = Math.floor(p), hi = Math.ceil(p); return s[lo] + (s[hi] - s[lo]) * (p - lo) }
export function trimmedMeanIqr(v: number[]): number { const q1 = percentile(v, 25), q3 = percentile(v, 75); const ub = q3 + 1.5 * (q3 - q1); const k = v.filter((x) => x >= 0 && x <= ub); return k.length ? k.reduce((a, b) => a + b, 0) / k.length : NaN }
export function median(v: number[]): number { const s = [...v].sort((a, b) => a - b); const m = s.length >> 1; return s.length ? (s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2) : NaN }
export function pearson(a: number[], b: number[]): number { const n = a.length; const ma = a.reduce((x, y) => x + y, 0) / n, mb = b.reduce((x, y) => x + y, 0) / n; let sab = 0, saa = 0, sbb = 0; for (let i = 0; i < n; i++) { const da = a[i] - ma, db = b[i] - mb; sab += da * db; saa += da * da; sbb += db * db } return sab / Math.sqrt(saa * sbb) }
export function spearman(a: number[], b: number[]): number { const rank = (v: number[]) => { const idx = v.map((x, i) => [x, i] as [number, number]).sort((p, q) => p[0] - q[0]); const r = new Array(v.length); let i = 0; while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++; const avg = (i + j) / 2 + 1; for (let k = i; k <= j; k++) r[idx[k][1]] = avg; i = j + 1 } return r }; return pearson(rank(a), rank(b)) }
export function skewness(v: number[]): number { const n = v.length; const m = v.reduce((a, b) => a + b, 0) / n; let s2 = 0, s3 = 0; for (const x of v) { const d = x - m; s2 += d * d; s3 += d * d * d } return (s3 / n) / Math.pow(s2 / n, 1.5) }
export function boxStats(v: number[]) { const q1 = percentile(v, 25), q3 = percentile(v, 75), med = percentile(v, 50); const iqr = q3 - q1; const inl = v.filter((x) => x >= q1 - 1.5 * iqr && x <= q3 + 1.5 * iqr); return { q1, q3, med, lo: Math.min(...inl), hi: Math.max(...inl), n: v.length } }

/** Level track of a pseudobulk profile on one chromosome. */
export function levelTrack(pb: Float32Array, chrom: ChromMeta, basis: 'chrmedian' | 'genome_trimmed' | 'raw', genomeRow: Float32Array, clip = 5): { y: (number | null)[]; ref: number; refLabel: string } {
  const seg = Array.from(pb.subarray(chrom.offset, chrom.offset + chrom.n))
  if (basis === 'raw') return { y: seg, ref: NaN, refLabel: 'pseudobulk X (mean over member cells)' }
  const ref = basis === 'chrmedian' ? median(seg) : trimmedMeanIqr(Array.from(genomeRow))
  const y = seg.map((v) => (ref > 0 && v > 0 ? Math.max(-clip, Math.min(clip, Math.log2(v / ref))) : v === 0 && ref > 0 ? -clip : null))
  return { y, ref, refLabel: basis === 'chrmedian' ? `chromosome median ${ref.toFixed(3)}` : `genome trimmed mean ${ref.toFixed(3)}` }
}
export function rollingMean(y: (number | null)[], w: number): (number | null)[] {
  if (w <= 1) return y
  const h = Math.floor(w / 2); const out: (number | null)[] = new Array(y.length)
  for (let i = 0; i < y.length; i++) { let s = 0, c = 0; for (let k = Math.max(0, i - h); k <= Math.min(y.length - 1, i + h); k++) { const v = y[k]; if (v !== null && Number.isFinite(v)) { s += v; c++ } } out[i] = c >= Math.max(1, Math.floor(w / 3)) ? s / c : null }
  return out
}
/** Break at omitted genomic intervals when ends are supplied; otherwise use step gaps. */
export function withGaps(mb: number[], y: (number | null)[], endsMb?: number[]): { x: (number | null)[]; y: (number | null)[] } {
  const d: number[] = []; for (let i = 1; i < mb.length; i++) d.push(mb[i] - mb[i - 1])
  const med = median(d); const X: (number | null)[] = [], Y: (number | null)[] = []
  for (let i = 0; i < mb.length; i++) { if (i > 0 && (endsMb ? mb[i] > endsMb[i - 1] + 1.01e-6 : mb[i] - mb[i - 1] > 1.5 * med)) { X.push(mb[i] - 1e-9); Y.push(null) } X.push(mb[i]); Y.push(y[i]) }
  return { x: X, y: Y }
}
/** Segments of a source on a chromosome from its accepted breakpoints (chromosome-local, half-open). */
export function segmentsFromBreakpoints(bps: ClusterBreakpoint[], source: string, chrom: string, n: number): { s: number; e: number }[] {
  const cuts = [0, ...bps.filter((b) => b.source === source && b.chromosome === chrom).map((b) => b.absolute_bin).sort((a, b) => a - b), n]
  const out: { s: number; e: number }[] = []; for (let i = 0; i + 1 < cuts.length; i++) if (cuts[i + 1] > cuts[i]) out.push({ s: cuts[i], e: cuts[i + 1] })
  return out
}
/** Cluster-level segments (genome-wide inclusive in the export) → chromosome-local half-open. */
export function clusterSegmentsLocal(segs: ClusterSegment[] | null, chroms: ChromMeta[]): (ClusterSegment & { s: number; e: number })[] {
  if (!segs) return []
  const off = new Map(chroms.map((c) => [c.name, c.offset]))
  return segs.map((g) => ({ ...g, s: g.bin_start - (off.get(g.chromosome) ?? 0), e: g.bin_end + 1 - (off.get(g.chromosome) ?? 0) }))
}
export function clusterOf(cell: IntegrationCell, key: string): string { const v = cell[key]; return v === undefined || v === null ? '' : String(v) }
export function clusterLevels(cells: IntegrationCell[], key: string): string[] { return [...new Set(cells.map((c) => clusterOf(c, key)))].filter((x) => x !== '').sort((a, b) => +a - +b) }
export function clusterColor(label: string, levels: string[]): string { const i = levels.indexOf(label); return i < 0 ? '#9ca3af' : TAB10[i % TAB10.length] }
/** Deterministic jitter (mulberry32, seed 42) for strip plots — display only. */
export function jitter(n: number, seed = 42, amp = 0.35): number[] { let a = seed >>> 0; const out: number[] = []; for (let i = 0; i < n; i++) { a += 0x6d2b79f5; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); out.push((((t ^ (t >>> 14)) >>> 0) / 4294967296 - 0.5) * 2 * amp) } return out }
