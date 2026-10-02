/** User-defined merge rules: types, evaluator, presets and legacy migration.
 * The evaluator is a 1:1 port of server/merge_rules.py's condition_holds and
 * evaluate_rules; tests/fixtures/merge_rules_parity.json is asserted against
 * both so they cannot drift. Rules are exploratory diagnostics. */
import { classifyAmbiguous, classifyTransition, scoreForMetric, type Estimator, type HatchMetric, type HatchScoreRow, type HatchScoresResult } from './hatchMerge'

export type HatchPattern = '/' | '\\' | 'x' | '-' | '|' | '.'
export const HATCH_PATTERNS: HatchPattern[] = ['/', '\\', 'x', '-', '|', '.']
export type ConditionKind = 'weak_flank' | 'sign_mismatch' | 'ambiguous' | 'max_size'
export const CONDITION_KINDS: ConditionKind[] = ['weak_flank', 'sign_mismatch', 'ambiguous', 'max_size']
export type Condition =
  | { kind: 'weak_flank'; metric: HatchMetric; estimator: Estimator; threshold: number }
  | { kind: 'sign_mismatch'; metric: HatchMetric; estimator: Estimator }
  | { kind: 'ambiguous'; srdMetric: 'srd' | 'srd_phi'; estimator: Estimator }
  | { kind: 'max_size'; maxBins: number }
export interface MergeRule {
  id: string; name: string; enabled: boolean; conditions: Condition[]
  mergeInto: { metric: HatchMetric; estimator: Estimator }
  hatch: { pattern: HatchPattern; color: string }
}
export type PresetName = 'blank' | 'weak_flank' | 'sign_mismatch' | 'ambiguous' | 'small'
export interface SegmentEval {
  start: number; end: number; matched: string[]; claimedBy: string | null
  side: 'left' | 'right' | null; score: number | null; blocked: null | 'gap' | 'nonfinite'
}

const RULE_COLORS = ['#b91c1c', '#1d4ed8', '#047857', '#7c3aed', '#b45309', '#0f766e']
const METRIC_LABEL: Record<HatchMetric, string> = { srd: 'SRD', srd_phi: 'SRD/√φ', delta_log2fc: 'Δlog2FC' }
const ESTIMATOR_LABEL: Record<Estimator, string> = { mean: 'mean', median: 'median', iqr_mean: 'IQR mean' }

const finite = (v: number | null | undefined): v is number => v !== null && v !== undefined && Number.isFinite(v)
function rowScore(row: HatchScoreRow, side: 'left' | 'right', metric: HatchMetric, estimator: Estimator): number | null {
  const v = scoreForMetric(row, side, metric, estimator)
  return finite(v) ? v : null
}

export function metricLabel(metric: HatchMetric, estimator: Estimator): string {
  return metric === 'delta_log2fc' ? `Δlog2FC (${ESTIMATOR_LABEL[estimator]})` : METRIC_LABEL[metric]
}

export function describeCondition(c: Condition): string {
  switch (c.kind) {
    case 'weak_flank': return `weaker flank |${metricLabel(c.metric, c.estimator)}| < ${c.threshold}`
    case 'sign_mismatch': return `${metricLabel(c.metric, c.estimator)} flank signs disagree`
    case 'ambiguous': return `weaker ${METRIC_LABEL[c.srdMetric]} flank ≠ weaker Δlog2FC (${ESTIMATOR_LABEL[c.estimator]}) flank`
    case 'max_size': return `segment ≤ ${c.maxBins} bins`
  }
}

export function conditionHolds(c: Condition, row: HatchScoreRow): boolean {
  if (c.kind === 'max_size') return row.end - row.start <= c.maxBins
  if (c.kind === 'ambiguous') {
    return classifyAmbiguous(rowScore(row, 'left', c.srdMetric, c.estimator), rowScore(row, 'right', c.srdMetric, c.estimator),
      rowScore(row, 'left', 'delta_log2fc', c.estimator), rowScore(row, 'right', 'delta_log2fc', c.estimator)) === 'ambiguous'
  }
  const l = rowScore(row, 'left', c.metric, c.estimator), r = rowScore(row, 'right', c.metric, c.estimator)
  if (c.kind === 'sign_mismatch') return classifyTransition(l, r)
  return l !== null && r !== null && Math.min(Math.abs(l), Math.abs(r)) < c.threshold
}

export function conditionDetail(c: Condition, row: HatchScoreRow): string {
  const pass = conditionHolds(c, row) ? 'pass' : 'fail'
  const f = (v: number | null) => (v === null ? 'N/A' : v.toFixed(3))
  if (c.kind === 'max_size') return `${describeCondition(c)}: ${row.end - row.start} bins → ${pass}`
  if (c.kind === 'ambiguous') {
    return `${describeCondition(c)}: ${METRIC_LABEL[c.srdMetric]} L/R ${f(rowScore(row, 'left', c.srdMetric, c.estimator))} / ${f(rowScore(row, 'right', c.srdMetric, c.estimator))}, Δlog2FC L/R ${f(rowScore(row, 'left', 'delta_log2fc', c.estimator))} / ${f(rowScore(row, 'right', 'delta_log2fc', c.estimator))} → ${pass}`
  }
  return `${describeCondition(c)}: L/R ${f(rowScore(row, 'left', c.metric, c.estimator))} / ${f(rowScore(row, 'right', c.metric, c.estimator))} → ${pass}`
}

export function ruleProblems(rule: MergeRule): string[] {
  const problems: string[] = []
  if (!rule.conditions.length) problems.push('add at least one condition')
  for (const c of rule.conditions) {
    if (c.kind === 'weak_flank' && !(Number.isFinite(c.threshold) && c.threshold > 0)) problems.push('threshold must be a positive number')
    if (c.kind === 'max_size' && !(Number.isInteger(c.maxBins) && c.maxBins >= 1)) problems.push('size limit must be a whole number ≥ 1')
  }
  return problems
}

export const activeRules = (rules: MergeRule[]): MergeRule[] => rules.filter((r) => r.enabled === true && ruleProblems(r).length === 0)

export function canRunMerge(rules: MergeRule[]): boolean {
  const enabled = rules.filter((r) => r.enabled === true)
  return enabled.length > 0 && enabled.every((r) => ruleProblems(r).length === 0)
}

export function evaluateRules(rows: HatchScoreRow[], rules: MergeRule[], allowGapCrossing: boolean): SegmentEval[] {
  const act = activeRules(rules)
  return rows.map((row) => {
    const matched = act.filter((r) => r.conditions.every((c) => conditionHolds(c, row))).map((r) => r.id)
    const ev: SegmentEval = { start: row.start, end: row.end, matched, claimedBy: null, side: null, score: null, blocked: null }
    if (!matched.length) return ev
    const rule = act.find((r) => r.id === matched[0])!
    ev.claimedBy = rule.id
    const l = rowScore(row, 'left', rule.mergeInto.metric, rule.mergeInto.estimator)
    const r = rowScore(row, 'right', rule.mergeInto.metric, rule.mergeInto.estimator)
    if (l === null || r === null) { ev.blocked = 'nonfinite'; return ev }
    ev.side = Math.abs(l) <= Math.abs(r) ? 'left' : 'right'
    ev.score = ev.side === 'left' ? Math.abs(l) : Math.abs(r)
    if (!allowGapCrossing && (ev.side === 'left' ? row.gap_left : row.gap_right)) ev.blocked = 'gap'
    return ev
  })
}

export function evaluateForSources(hatchScores: Record<string, HatchScoresResult>, sources: string[], chrom: string, rules: MergeRule[], allowGapCrossing: boolean): Record<string, SegmentEval[]> {
  const out: Record<string, SegmentEval[]> = {}
  for (const s of sources) { const hs = hatchScores[`${s}|${chrom}`]; if (hs) out[s] = evaluateRules(hs.rows, rules, allowGapCrossing) }
  return out
}

export function ruleFlagCounts(evals: Record<string, SegmentEval[]>): Record<string, number> {
  const counts: Record<string, number> = {}
  for (const list of Object.values(evals)) for (const e of list) for (const id of e.matched) counts[id] = (counts[id] || 0) + 1
  return counts
}

export function newCondition(kind: ConditionKind): Condition {
  switch (kind) {
    case 'weak_flank': return { kind, metric: 'srd_phi', estimator: 'mean', threshold: 3.3 }
    case 'sign_mismatch': return { kind, metric: 'srd_phi', estimator: 'mean' }
    case 'ambiguous': return { kind, srdMetric: 'srd_phi', estimator: 'mean' }
    case 'max_size': return { kind, maxBins: 20 }
  }
}

export function defaultMergeInto(conditions: Condition[]): { metric: HatchMetric; estimator: Estimator } {
  for (const c of conditions) {
    if (c.kind === 'weak_flank' || c.kind === 'sign_mismatch') return { metric: c.metric, estimator: c.estimator }
    if (c.kind === 'ambiguous') return { metric: c.srdMetric, estimator: c.estimator }
  }
  return { metric: 'srd_phi', estimator: 'mean' }
}

let ruleSeq = 0
function newRuleId(existing: MergeRule[]): string {
  let id = ''
  do { id = `rule-${Date.now().toString(36)}-${++ruleSeq}` } while (existing.some((r) => r.id === id))
  return id
}

const PRESETS: Record<PresetName, { name: string; kinds: ConditionKind[]; pattern: HatchPattern }> = {
  blank: { name: 'new rule', kinds: [], pattern: '/' },
  weak_flank: { name: 'weak SRD/√φ flank', kinds: ['weak_flank'], pattern: 'x' },
  sign_mismatch: { name: 'flank signs disagree', kinds: ['sign_mismatch'], pattern: '/' },
  ambiguous: { name: 'ambiguous flank preference', kinds: ['ambiguous'], pattern: '\\' },
  small: { name: 'small segment', kinds: ['max_size'], pattern: '.' },
}

export function makeRule(preset: PresetName, existing: MergeRule[]): MergeRule {
  const p = PRESETS[preset]
  const conditions = p.kinds.map(newCondition)
  return { id: newRuleId(existing), name: p.name, enabled: true, conditions, mergeInto: defaultMergeInto(conditions),
    hatch: { pattern: p.pattern, color: RULE_COLORS[existing.length % RULE_COLORS.length] } }
}

/** Old sessions stored three fixed layers (`hatchLayers`); map each to the
 * equivalent preset rule, keeping its settings and on/off state. */
export function migrateHatchLayers(old: any): MergeRule[] {
  if (!old || typeof old !== 'object') return []
  const out: MergeRule[] = []
  const t = old.transition, p = old.proposal, a = old.ambiguous
  if (t) {
    const metric: HatchMetric = t.metric ?? 'srd_phi', estimator: Estimator = t.estimator ?? 'mean'
    out.push({ id: 'legacy-transition', name: 'transition zone (from saved session)', enabled: !!t.on,
      conditions: [{ kind: 'sign_mismatch', metric, estimator }], mergeInto: { metric, estimator }, hatch: { pattern: '/', color: RULE_COLORS[0] } })
  }
  if (p) {
    const metric: HatchMetric = p.metric ?? 'srd_phi', estimator: Estimator = p.estimator ?? 'mean'
    out.push({ id: 'legacy-proposal', name: 'merge proposal (from saved session)', enabled: !!p.on,
      conditions: [{ kind: 'weak_flank', metric, estimator, threshold: p.threshold ?? 3.3 }], mergeInto: { metric, estimator }, hatch: { pattern: 'x', color: RULE_COLORS[1] } })
  }
  if (a) {
    const estimator: Estimator = a.estimator ?? 'mean'
    out.push({ id: 'legacy-ambiguous', name: 'ambiguous flank preference (from saved session)', enabled: !!a.on,
      conditions: [{ kind: 'ambiguous', srdMetric: 'srd_phi', estimator }], mergeInto: { metric: 'srd_phi', estimator }, hatch: { pattern: '\\', color: RULE_COLORS[2] } })
  }
  return out
}
