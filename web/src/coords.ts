/**
 * Canonical coordinate conversions. One place for every mapping between
 * node-local slices, chromosome-local bins, genome-wide indices and genomic
 * coordinates. Genomic intervals are 1-based inclusive [start_bp, end_bp].
 */
import type { NodeData, SelectionCoords } from './types'

export const MB = 1e6

export function clampWindow(s: number, e: number, n: number): { s: number; e: number; clipped: { left: number; right: number } } {
  const cs = Math.max(0, Math.min(n, s))
  const ce = Math.max(0, Math.min(n, e))
  return { s: cs, e: Math.max(ce, cs), clipped: { left: cs - s, right: e - ce } }
}

/** Default peak window: [b-half, b+half) clipped to the node. */
export function peakWindow(b: number, n: number, half = 20) {
  return clampWindow(b - half, b + half, n)
}

/** AD membership: boundary j (value d[j-1]) is included when max(s,1) <= j < e. */
export function adIndices(s: number, e: number, n: number): number[] {
  const out: number[] = []
  for (let j = Math.max(s, 1); j < Math.min(e, n); j++) out.push(j)
  return out
}

export function selectionCoords(node: NodeData, s: number, e: number): SelectionCoords {
  const n = node.n
  const nBins = Math.max(0, e - s)
  const startBp = nBins ? node.start_bp[s] : NaN
  const endBp = nBins ? node.end_bp[e - 1] : NaN
  let retained = 0
  let gaps = 0
  for (let i = s; i < e; i++) {
    retained += node.end_bp[i] - node.start_bp[i] + 1
    if (i > s && node.start_bp[i] !== node.end_bp[i - 1] + 1) gaps++
  }
  const ad = adIndices(s, e, n)
  return {
    s, e, nBins,
    chromStart: node.start + s, chromEnd: node.start + e,
    genomeStart: node.offset_genome + node.start + s, genomeEnd: node.offset_genome + node.start + e,
    startBp, endBp,
    spanBp: nBins ? endBp - startBp + 1 : 0,
    retainedBp: retained,
    startMb: startBp / MB, endMb: endBp / MB,
    spanMb: nBins ? (endBp - startBp + 1) / MB : 0,
    retainedMb: retained / MB,
    gaps,
    adCount: ad.length, adFirst: ad.length ? ad[0] : null, adLast: ad.length ? ad[ad.length - 1] : null,
  }
}

/** Mb coordinate of a node-local boundary b: drawn at the start of its right-hand bin. */
export function boundaryMb(node: NodeData, b: number): number {
  if (b >= node.n) return node.end_bp[node.n - 1] / MB
  return node.start_bp[b] / MB
}

/** Bin whose [start,end] contains the Mb position, or nearest bin edge. */
export function mbToBin(node: NodeData, mb: number): number {
  const bp = mb * MB
  let lo = 0, hi = node.n - 1
  if (bp <= node.start_bp[0]) return 0
  if (bp >= node.end_bp[hi]) return hi
  while (lo < hi) {
    const mid = (lo + hi) >> 1
    if (node.end_bp[mid] < bp) lo = mid + 1
    else hi = mid
  }
  return lo
}

/** Snap an Mb coordinate to the nearest retained-bin *boundary* (0..n). */
export function mbToBoundary(node: NodeData, mb: number): number {
  const bp = mb * MB
  const i = mbToBin(node, mb)
  const mid = (node.start_bp[i] + node.end_bp[i]) / 2
  return bp < mid ? i : i + 1
}

/** Detect coordinate gaps (> 3x median step) to break plotted lines. */
export function gapMask(startBp: number[]): boolean[] {
  const n = startBp.length
  const diffs: number[] = []
  for (let i = 1; i < n; i++) diffs.push(startBp[i] - startBp[i - 1])
  const sorted = [...diffs].sort((a, b) => a - b)
  const med = sorted.length ? sorted[Math.floor(sorted.length / 2)] : 0
  const out = new Array(n).fill(false)
  for (let i = 1; i < n; i++) if (diffs[i - 1] > 3 * med) out[i] = true
  return out
}

/** Insert null after every gap so Plotly breaks the line. */
export function breakGaps(xs: number[], ys: (number | null)[], gapAt: boolean[]): { x: (number | null)[]; y: (number | null)[] } {
  const x: (number | null)[] = []
  const y: (number | null)[] = []
  for (let i = 0; i < xs.length; i++) {
    if (gapAt[i] && i > 0) { x.push(xs[i] - 1e-9); y.push(null) }
    x.push(xs[i]); y.push(ys[i])
  }
  return { x, y }
}
