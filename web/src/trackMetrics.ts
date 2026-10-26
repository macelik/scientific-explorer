/** Live per-flank colour metrics for the Chromosome tracks, plus mappings for
 * archived ids and saved sessions, and the candidate-audit reason derivation. */
import type { HatchScoreRow } from './hatchMerge'
import type { Metric } from './integration'

export type TrackMetric = 'delta_median' | 'delta_mean' | 'delta_iqr_mean' | 'srd' | 'srd_phi'
const EFFECT_BINS = { bins: [0.2, 0.5, 0.75, 1.0], binLabels: ['< 0.20', '0.20–0.50', '0.50–0.75', '0.75–1.00', '≥ 1.00'] }
export const TRACK_METRICS: { id: TrackMetric; label: string; kind: 'effect' | 'stat'; bins: number[]; binLabels: string[] }[] = [
  { id: 'delta_median', label: 'Δlog2FC · median levels', kind: 'effect', ...EFFECT_BINS },
  { id: 'delta_mean', label: 'Δlog2FC · mean levels', kind: 'effect', ...EFFECT_BINS },
  { id: 'delta_iqr_mean', label: 'Δlog2FC · two-sided IQR-mean levels', kind: 'effect', ...EFFECT_BINS },
  { id: 'srd', label: '|SRD| · signed root deviance', kind: 'stat', bins: [1.5, 3.3, 5.0, 8.0], binLabels: ['< 1.5', '1.5–3.3', '3.3–5.0', '5.0–8.0', '≥ 8.0'] },
  { id: 'srd_phi', label: '|SRD/√φ| · dispersion-scaled', kind: 'stat', bins: [1.5, 3.3, 5.0], binLabels: ['< 1.5', '1.5–3.3', '3.3–5.0', '≥ 5.0'] },
]

const finite = (v: unknown): v is number => typeof v === 'number' && Number.isFinite(v)

/** Value of metric m on its own smaller-|value| flank (missing side loses, ties left):
 * signed for Δlog2FC, absolute for SRD-type metrics. */
export function colourValue(row: HatchScoreRow, m: TrackMetric): { value: number | null; side: 'left' | 'right' | null } {
  const l = (row as any)[`${m}_left`], r = (row as any)[`${m}_right`]
  const fl = finite(l) ? l : null, fr = finite(r) ? r : null
  if (fl === null && fr === null) return { value: null, side: null }
  const side: 'left' | 'right' = fr === null ? 'left' : fl === null ? 'right' : Math.abs(fl) <= Math.abs(fr) ? 'left' : 'right'
  const v = (side === 'left' ? fl : fr) as number
  return { value: m === 'srd' || m === 'srd_phi' ? Math.abs(v) : v, side }
}

const ARCHIVED_TO_TRACK: Record<string, TrackMetric> = {
  delta_median: 'delta_median', delta_mean: 'delta_mean', log2fc: 'delta_mean', delta_trim: 'delta_iqr_mean', srd: 'srd', srd_phi: 'srd_phi',
}
export const archivedColourToTrack = (id: unknown): TrackMetric => (typeof id === 'string' && ARCHIVED_TO_TRACK[id]) || 'delta_mean'

/** The archived diagnostics panel still reads the archived tables; give it the matching effect column. */
export function diagnosticsEffectMetric(m: TrackMetric): Metric {
  return m === 'delta_median' ? 'delta_median' : m === 'delta_iqr_mean' ? 'delta_trim' : 'delta_mean'
}

/** Map a saved seg-settings object from before Rule 0 onto the current keys. */
export function migrateSegSettings(old: any): any {
  const o: any = { ...(old || {}) }
  if (o.colorMetric === undefined && 'metric' in o) o.colorMetric = archivedColourToTrack(o.metric)
  delete o.metric; delete o.effectMetric; delete o.flank
  if (o.labels && typeof o.labels === 'object') {
    const l: any = { ...o.labels }
    if (l.signs === undefined && typeof l.consistency === 'boolean') l.signs = l.consistency
    delete l.consistency; delete l.proposal
    o.labels = l
  }
  return o
}

const CUTOFF = 0.001
const NOTE = '(cutoff inferred: all accepted candidates have p < 0.001)'
const fmtP = (p: number) => String(Number(p.toPrecision(3)))

/** Candidate-audit reason: the recorded geometry reason, or a p-value reason
 * derived for tested-and-rejected candidates the archived table leaves blank. */
export function derivedRejectionReason(r: { selected_for_significance: boolean; accepted: boolean; rejection_reasons: string | null; p_local: number | null; p_global: number | null }): string | null {
  if (r.rejection_reasons) return r.rejection_reasons
  if (!r.selected_for_significance || r.accepted) return null
  if (finite(r.p_local) && r.p_local >= CUTOFF) return `local permutation p = ${fmtP(r.p_local)} ≥ ${CUTOFF} ${NOTE}`
  if (finite(r.p_global) && r.p_global >= CUTOFF) return `global permutation p = ${fmtP(r.p_global)} ≥ ${CUTOFF} ${NOTE}`
  return `not accepted (no recorded reason) ${NOTE}`
}
