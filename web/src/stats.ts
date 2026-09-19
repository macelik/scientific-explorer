/**
 * Descriptive statistics for selected windows, computed in the browser from
 * fetched arrays. These are descriptive only (not significance tests) and are
 * verified for numerical parity against numpy/pandas in tests/.
 */
import type { AcfResult, NodeData, SummaryStats } from './types'
import { adIndices } from './coords'

export function mean(v: number[]): number | null { return v.length ? v.reduce((a, b) => a + b, 0) / v.length : null }
export function median(v: number[]): number | null {
  if (!v.length) return null
  const s = [...v].sort((a, b) => a - b)
  const m = s.length >> 1
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2
}
export function stdPop(v: number[]): number | null {
  if (!v.length) return null
  const mu = mean(v)!
  return Math.sqrt(v.reduce((a, b) => a + (b - mu) * (b - mu), 0) / v.length)
}
export function pearson(a: number[], b: number[]): { r: number | null; reason?: string } {
  const n = Math.min(a.length, b.length)
  if (n < 3) return { r: null, reason: `needs >= 3 paired observations (have ${n})` }
  const ma = mean(a.slice(0, n))!, mb = mean(b.slice(0, n))!
  let sab = 0, saa = 0, sbb = 0
  for (let i = 0; i < n; i++) { const da = a[i] - ma, db = b[i] - mb; sab += da * db; saa += da * da; sbb += db * db }
  if (saa === 0 || sbb === 0) return { r: null, reason: 'constant paired array' }
  return { r: sab / Math.sqrt(saa * sbb) }
}
/** Lag-k autocorrelation = Pearson(v[:-k], v[k:]). */
export function lagCorr(v: number[], k: number): { r: number | null; reason?: string } {
  if (v.length - k < 3) return { r: null, reason: `needs >= 3 pairs at lag ${k} (have ${Math.max(0, v.length - k)})` }
  return pearson(v.slice(0, v.length - k), v.slice(k))
}
export function summarize(values: number[]): SummaryStats {
  const finite = values.filter((x) => Number.isFinite(x))
  const nonfinite = values.length - finite.length
  const mu = mean(finite)
  const sd = stdPop(finite)
  let cv: number | null = null, cvReason: string | undefined
  if (mu === null) cvReason = 'empty'
  else if (mu === 0) cvReason = 'mean is zero'
  else cv = sd! / mu
  const l1 = lagCorr(finite, 1)
  return { n: finite.length, mean: mu, median: median(finite), std: sd, cv, cvReason, lag1: l1.r, lag1Reason: l1.reason, nonfinite }
}
export function acf(values: number[], maxLag = 10): AcfResult {
  const finite = values.filter((x) => Number.isFinite(x))
  const lags: number[] = [], vals: (number | null)[] = [], reasons: (string | null)[] = []
  for (let k = 1; k <= maxLag; k++) { const r = lagCorr(finite, k); lags.push(k); vals.push(r.r); reasons.push(r.reason ?? null) }
  return { lags, values: vals, reasons }
}
export function ecdf(values: number[]): { x: number[]; y: number[] } {
  const s = [...values].filter(Number.isFinite).sort((a, b) => a - b)
  return { x: s, y: s.map((_, i) => (i + 1) / s.length) }
}

/** AD samples of a node-local window [s,e): boundaries j with max(s,1)<=j<e, value d[j-1]. */
export function adSamples(node: NodeData, s: number, e: number): { j: number[]; values: number[] } {
  const j = adIndices(s, e, node.n)
  return { j, values: j.map((jj) => node.ad[jj - 1]) }
}
export function xSamples(node: NodeData, s: number, e: number): { i: number[]; values: number[] } {
  const i: number[] = []
  for (let k = Math.max(0, s); k < Math.min(e, node.n); k++) i.push(k)
  return { i, values: i.map((k) => node.x[k]) }
}

/** pandas Series.rolling(window, center=True, min_periods).mean() for odd windows. */
export function rollingMean(x: number[], window: number, minPeriods: number): (number | null)[] {
  const n = x.length
  const half = Math.floor((window - 1) / 2)
  const out: (number | null)[] = new Array(n)
  // prefix sums over finite values
  const ps = new Float64Array(n + 1), pc = new Int32Array(n + 1)
  for (let i = 0; i < n; i++) { const f = Number.isFinite(x[i]); ps[i + 1] = ps[i] + (f ? x[i] : 0); pc[i + 1] = pc[i] + (f ? 1 : 0) }
  for (let i = 0; i < n; i++) {
    const lo = Math.max(0, i - half), hi = Math.min(n - 1, i + (window - 1 - half))
    const c = pc[hi + 1] - pc[lo]
    out[i] = c >= minPeriods ? (ps[hi + 1] - ps[lo]) / c : null
  }
  return out
}

export interface Log2Track { values: (number | null)[]; status: ('ok' | 'zero' | 'nan' | 'undefined')[]; clip: number; undefinedReason: string | null }
/** log2(rolling_mean(x)/baseline); zero-signal windows are -inf, displayed clipped. */
export function log2Track(roll: (number | null)[], baseline: number | null, clip: number): Log2Track {
  const n = roll.length
  const values: (number | null)[] = new Array(n), status: Log2Track['status'] = new Array(n)
  if (baseline === null || !(baseline > 0)) {
    for (let i = 0; i < n; i++) { values[i] = null; status[i] = 'undefined' }
    return { values, status, clip, undefinedReason: baseline === 0 ? 'reference baseline is exactly zero: log2 ratio undefined (no substitution applied)' : 'reference baseline unavailable' }
  }
  for (let i = 0; i < n; i++) {
    const r = roll[i]
    if (r === null || !Number.isFinite(r)) { values[i] = null; status[i] = 'nan'; continue }
    if (r <= 0) { values[i] = -clip; status[i] = 'zero'; continue }
    const v = Math.log2(r / baseline)
    values[i] = Math.max(-clip, Math.min(clip, v)); status[i] = 'ok'
  }
  return { values, status, clip, undefinedReason: null }
}

export interface Peak { b: number; ad: number; isWinner: boolean; isEvent: boolean; relaxed: boolean; mb: number }
/**
 * Port of _peak_panel.find_spaced_peaks with the event-region rule from the plan:
 * winner first; then (if an event region is active) the strongest candidate inside
 * the region, even if spacing must be relaxed (flagged); then spaced local maxima.
 */
export function findCompetingPeaks(node: NodeData, winnerB: number, nPeaks: number, minSepMb: number,
  eventRange: { startMb: number; endMb: number } | null, mbOfBoundary: (b: number) => number): Peak[] {
  const d = node.ad, n = d.length
  if (n === 0) return []
  const isPeak = new Array(n).fill(false)
  for (let k = 1; k < n - 1; k++) isPeak[k] = d[k] > d[k - 1] && d[k] > d[k + 1]
  isPeak[0] = n > 1 ? d[0] > d[1] : true
  isPeak[n - 1] = n > 1 ? d[n - 1] > d[n - 2] : true
  const cand = new Set<number>()
  for (let k = 0; k < n; k++) if (isPeak[k]) cand.add(k)
  const w0 = winnerB - 1
  cand.add(w0)
  const sorted = [...cand].sort((a, b) => d[b] - d[a])
  const chosen: Peak[] = [{ b: winnerB, ad: d[w0], isWinner: true, isEvent: false, relaxed: false, mb: mbOfBoundary(winnerB) }]
  const chosenMb = [chosen[0].mb]
  if (eventRange) {
    let best = -1
    for (let k = 0; k < n; k++) {
      const mb = mbOfBoundary(k + 1)
      if (mb >= eventRange.startMb && mb <= eventRange.endMb && (best < 0 || d[k] > d[best])) best = k
    }
    if (best >= 0 && best !== w0) {
      const mb = mbOfBoundary(best + 1)
      const relaxed = chosenMb.some((m) => Math.abs(mb - m) < minSepMb)
      chosen.push({ b: best + 1, ad: d[best], isWinner: false, isEvent: true, relaxed, mb })
      chosenMb.push(mb)
    } else if (best === w0) {
      chosen[0].isEvent = true
    }
  }
  for (const k of sorted) {
    if (chosen.length >= nPeaks) break
    if (chosen.some((p) => p.b === k + 1)) continue
    const mb = mbOfBoundary(k + 1)
    if (chosenMb.every((m) => Math.abs(mb - m) >= minSepMb)) {
      chosen.push({ b: k + 1, ad: d[k], isWinner: false, isEvent: false, relaxed: false, mb })
      chosenMb.push(mb)
    }
  }
  return chosen
}

export function fmt(v: number | null | undefined, digits = 3): string {
  if (v === null || v === undefined || Number.isNaN(v)) return 'N/A'
  if (!Number.isFinite(v)) return v > 0 ? '+∞' : '−∞'
  return v.toFixed(digits)
}
export function fmtP(p: number, exact?: string): string {
  return `${p.toFixed(6)}${exact ? ` (${exact})` : ''}`
}
