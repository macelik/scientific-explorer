/** New exploratory hatching/merge: types, API client and classifiers.
 * Kept separate from integration.ts's existing archived-only consistency/
 * proposal classes (consistencyClass, pickFlank) and from the original merge
 * prototype. The classification rules here are a 1:1 port of
 * server/hatch_merge.py's classify_transition/classify_merge_proposal/
 * classify_ambiguous, so raw scores can be reclassified client-side as the
 * user adjusts a threshold, without a network round trip. */
import { hatchShapes } from './integration'

export type HatchMetric = 'srd' | 'srd_phi' | 'delta_log2fc'
export type Estimator = 'mean' | 'median' | 'iqr_mean'

export interface HatchScoreRow {
  start: number; end: number; phi: number | null; gap_left: boolean; gap_right: boolean
  srd_left: number | null; srd_right: number | null; srd_phi_left: number | null; srd_phi_right: number | null
  delta_mean_left: number | null; delta_mean_right: number | null
  delta_median_left: number | null; delta_median_right: number | null
  delta_iqr_mean_left: number | null; delta_iqr_mean_right: number | null
}
export interface HatchScoresResult {
  segments: { start: number; end: number }[]; phi: number | null; rows: HatchScoreRow[]
  source: string; chrom: string; start_bp: number[]; end_bp: number[]
}
export interface MergeStep {
  step: number; segments: { start: number; end: number }[]
  removed_boundary: number | null; merged_interval: [number, number] | null
  score: number | null; phi: number | null
}
export interface MergeResult {
  original: { start: number; end: number }[]; final: { start: number; end: number }[]
  steps: MergeStep[]; settings: Record<string, unknown>
  source: string; chrom: string; start_bp: number[]; end_bp: number[]; provenance: unknown
}
export interface MergeRequest {
  chrom: string; source: string; metric: HatchMetric; threshold: number; estimator?: Estimator
  veto_transition?: boolean; veto_ambiguous?: boolean; small_max_bins?: number; allow_gap_crossing?: boolean
}

async function postJson<T>(url: string, body: unknown): Promise<T> {
  const res = await fetch(url, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) })
  if (!res.ok) throw new Error((await res.json().catch(() => ({}))).detail ?? res.statusText)
  return res.json()
}
export const fetchHatchScores = (chrom: string, source: string) =>
  postJson<HatchScoresResult>('/api/integration/hatch-scores', { chrom, source })
export const runHatchMerge = (req: MergeRequest) => postJson<MergeResult>('/api/integration/hatch-merge', req)

const finite = (v: number | null | undefined): v is number => v !== null && v !== undefined && Number.isFinite(v)

export function weakerSide(left: number | null, right: number | null): 'left' | 'right' | 'tie' | null {
  if (left === null || right === null) return null
  const al = Math.abs(left), ar = Math.abs(right)
  if (al === ar) return 'tie'
  return al < ar ? 'left' : 'right'
}
export function classifyTransition(left: number | null, right: number | null): boolean {
  if (!finite(left) || !finite(right)) return false
  return left * right < 0
}
export function classifyMergeProposal(left: number | null, right: number | null, threshold: number) {
  const el = finite(left) && Math.abs(left) < threshold
  const er = finite(right) && Math.abs(right) < threshold
  let weakerSide_: 'left' | 'right' | 'both' | null = null
  if (el && er) { const w = weakerSide(left, right); weakerSide_ = w === 'tie' ? 'both' : (w as 'left' | 'right') }
  else if (el) weakerSide_ = 'left'
  else if (er) weakerSide_ = 'right'
  return { eligibleLeft: el, eligibleRight: er, weakerSide: weakerSide_ }
}
export function classifyAmbiguous(srdLeft: number | null, srdRight: number | null, fcLeft: number | null, fcRight: number | null): 'consistent' | 'ambiguous' | 'no_unique_preference' | 'undefined' {
  if (![srdLeft, srdRight, fcLeft, fcRight].every(finite)) return 'undefined'
  const srdSide = weakerSide(srdLeft, srdRight), fcSide = weakerSide(fcLeft, fcRight)
  if (srdSide === 'tie' || fcSide === 'tie') return 'no_unique_preference'
  return srdSide === fcSide ? 'consistent' : 'ambiguous'
}
export function scoreForMetric(row: HatchScoreRow, side: 'left' | 'right', metric: HatchMetric, estimator: Estimator): number | null {
  if (metric === 'delta_log2fc') return row[`delta_${estimator}_${side}` as const]
  return row[`${metric}_${side}` as const]
}
export { hatchShapes }
