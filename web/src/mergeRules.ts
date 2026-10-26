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
  id: string; type?: 'conditions'; name: string; enabled: boolean; conditions: Condition[]
  mergeInto: { metric: HatchMetric; estimator: Estimator }
  hatch: { pattern: HatchPattern; color: string }
}
export type PresetName = 'blank' | 'weak_flank' | 'sign_mismatch' | 'ambiguous' | 'small'
export interface SegmentEval {
  start: number; end: number; matched: string[]; claimedBy: string | null
  side: 'left' | 'right' | null; score: number | null; blocked: null | 'gap' | 'nonfinite'
  rule0: Rule0Eval | null
}

export type Rule0Criterion = 'delta_median' | 'delta_mean' | 'delta_iqr_mean' | 'srd_phi'
export const RULE0_CRITERIA: Rule0Criterion[] = ['delta_median', 'delta_mean', 'delta_iqr_mean', 'srd_phi']
export const CRITERION_LABEL: Record<Rule0Criterion, string> = {
  delta_median: 'min |Δlog2FC median| flank', delta_mean: 'min |Δlog2FC mean| flank',
  delta_iqr_mean: 'min |Δlog2FC IQR-mean| flank', srd_phi: 'min |SRD/√φ| flank',
}
export const RULE0_COLOR = '#b91c1c', SIGN_COLOR = '#d97706', FLANK_COLOR = '#7c3aed'
export interface Rule0 {
  id: 'rule0'; type: 'flank_threshold'; name: string; criteria: Rule0Criterion[]; combine: 'all' | 'any'
  tau: number; srdRef: number; maxBins: number; mergeSignMismatch: boolean; mergeFlankConflict: boolean
  hatch: { pattern: HatchPattern; color: string }
}
export type AnyRule = Rule0 | MergeRule
export const isRule0 = (r: AnyRule): r is Rule0 => (r as Rule0).type === 'flank_threshold'
export interface CriterionEval { criterion: Rule0Criterion; side: 'left' | 'right' | null; value: number | null; threshold: number; fraction: number | null; passed: boolean }
export interface Rule0Eval {
  start: number; end: number; evaluated: boolean; criteria: CriterionEval[]; passed: boolean
  signMismatch: boolean; flankConflict: boolean; heldBack: ('signs' | 'flanks')[]
  outcome: 'none' | 'would_merge' | 'held_back' | 'signs_only'; side: 'left' | 'right' | null; score: number | null; blocked: null | 'gap'
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

function criterionValue(row: HatchScoreRow, side: 'left' | 'right', c: Rule0Criterion): number | null {
  const v = (row as any)[`${c}_${side}`]
  return finite(v) ? v : null
}
function smallerSide(l: number | null, r: number | null): 'left' | 'right' | null {
  if (l === null && r === null) return null
  if (r === null) return 'left'
  if (l === null) return 'right'
  return Math.abs(l) <= Math.abs(r) ? 'left' : 'right'
}

export function evaluateRule0(row: HatchScoreRow, r0: Rule0, allowGapCrossing: boolean): Rule0Eval {
  const out: Rule0Eval = { start: row.start, end: row.end, evaluated: false, criteria: [], passed: false, signMismatch: false,
    flankConflict: false, heldBack: [], outcome: 'none', side: null, score: null, blocked: null }
  if (row.end - row.start > r0.maxBins) return out
  out.evaluated = true
  for (const c of r0.criteria) {
    const l = criterionValue(row, 'left', c), r = criterionValue(row, 'right', c)
    const threshold = c === 'srd_phi' ? r0.srdRef : r0.tau
    const side = smallerSide(l, r)
    const value = side === null ? null : Math.abs((side === 'left' ? l : r) as number)
    out.criteria.push({ criterion: c, side, value, threshold, fraction: value === null ? null : value / threshold, passed: value !== null && value < threshold })
  }
  const passing = out.criteria.filter((c) => c.passed)
  out.passed = r0.combine === 'all' ? passing.length === out.criteria.length : passing.length > 0
  const sl = criterionValue(row, 'left', 'srd_phi'), sr = criterionValue(row, 'right', 'srd_phi')
  out.signMismatch = sl !== null && sr !== null && sl * sr < 0
  if (!out.passed) { if (out.signMismatch) out.outcome = 'signs_only'; return out }
  out.flankConflict = new Set(passing.map((c) => c.side)).size > 1
  const fractions = passing.map((c) => c.fraction as number)
  out.score = r0.combine === 'all' ? Math.max(...fractions) : Math.min(...fractions)
  const leftFirst = (s: 'left' | 'right' | null) => (s === 'left' ? 0 : 1)
  const best = [...passing].sort((a, b) => (a.fraction as number) - (b.fraction as number) || leftFirst(a.side) - leftFirst(b.side))[0]
  out.side = out.flankConflict ? best.side : passing[0].side
  if (out.signMismatch && !r0.mergeSignMismatch) out.heldBack.push('signs')
  if (out.flankConflict && !r0.mergeFlankConflict) out.heldBack.push('flanks')
  if (out.heldBack.length) { out.outcome = 'held_back'; return out }
  out.outcome = 'would_merge'
  if (!allowGapCrossing && (out.side === 'left' ? row.gap_left : row.gap_right)) out.blocked = 'gap'
  return out
}

export function rule0Problems(r0: Rule0): string[] {
  const p: string[] = []
  if (!r0.criteria.length) p.push('check at least one criterion')
  if (!(Number.isFinite(r0.tau) && r0.tau > 0)) p.push('τ must be a positive number')
  if (!(Number.isFinite(r0.srdRef) && r0.srdRef > 0)) p.push('SRD/√φ ref must be a positive number')
  if (!(Number.isInteger(r0.maxBins) && r0.maxBins >= 1)) p.push('size limit must be a whole number ≥ 1')
  return p
}

export function describeRule0(r0: Rule0): string {
  const crit = r0.criteria.map((c) => `${CRITERION_LABEL[c]} < ${c === 'srd_phi' ? r0.srdRef : r0.tau}`).join(', ')
  return `${r0.combine} of: ${crit} · segments ≤ ${r0.maxBins} bins · signs disagree: ${r0.mergeSignMismatch ? 'merged' : 'held back'} · flanks disagree: ${r0.mergeFlankConflict ? 'merged' : 'held back'}`
}

export function rule0Lines(e: Rule0Eval): string[] {
  if (!e.evaluated) return ['Rule 0: segment longer than the size limit, not evaluated']
  const f = (v: number | null) => (v === null ? 'N/A' : v.toFixed(3))
  const lines = e.criteria.map((c) => `${CRITERION_LABEL[c.criterion]}: ${c.side ?? '–'} ${f(c.value)} vs ${c.threshold} (fraction ${f(c.fraction)}) → ${c.passed ? 'pass' : 'fail'}`)
  lines.push(`signs ${e.signMismatch ? 'disagree' : 'agree or N/A'} · criteria flanks ${e.flankConflict ? 'disagree' : 'agree'}`)
  const outcome = e.outcome === 'would_merge' ? `would merge ${e.side}${e.blocked === 'gap' ? ', blocked by gap' : ''} (score ${f(e.score)})`
    : e.outcome === 'held_back' ? `⚠ held back (${e.heldBack.map((h) => (h === 'signs' ? 'signs disagree' : 'flanks disagree')).join(', ')}); would merge ${e.side}`
    : e.outcome === 'signs_only' ? 'fails the test; signs disagree' : 'fails the test'
  lines.push(`Rule 0 → ${outcome}`)
  return lines
}

export function rule0FromSeg(s: { criteria: Rule0Criterion[]; combine: 'all' | 'any'; tau: number; srdRef: number; smallMax: number; mergeSignMismatch: boolean; mergeFlankConflict: boolean }): Rule0 {
  return { id: 'rule0', type: 'flank_threshold', name: 'Rule 0 — shown on the tracks', criteria: [...s.criteria], combine: s.combine,
    tau: s.tau, srdRef: s.srdRef, maxBins: s.smallMax, mergeSignMismatch: s.mergeSignMismatch, mergeFlankConflict: s.mergeFlankConflict,
    hatch: { pattern: 'x', color: RULE0_COLOR } }
}

export function rule0CaseCounts(evals: Record<string, SegmentEval[]>) {
  const c = { wouldMerge: 0, signs: 0, flanks: 0, heldBack: 0, signsOnly: 0, gap: 0 }
  for (const list of Object.values(evals)) for (const e of list) {
    const r = e.rule0; if (!r) continue
    if (r.outcome === 'would_merge') { c.wouldMerge++; if (r.signMismatch) c.signs++; if (r.flankConflict) c.flanks++; if (r.blocked === 'gap') c.gap++ }
    else if (r.outcome === 'held_back') c.heldBack++
    else if (r.outcome === 'signs_only') c.signsOnly++
  }
  return c
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

export function canRunMerge(rule0: Rule0, rules: MergeRule[]): boolean {
  return rule0Problems(rule0).length === 0 && rules.filter((r) => r.enabled === true).every((r) => ruleProblems(r).length === 0)
}

export function evaluateRules(rows: HatchScoreRow[], rules: AnyRule[], allowGapCrossing: boolean): SegmentEval[] {
  const first = rules[0]
  const r0 = first && isRule0(first) && rule0Problems(first).length === 0 ? first : null
  const act = activeRules(rules.filter((r): r is MergeRule => !isRule0(r)))
  return rows.map((row) => {
    const r0ev = r0 ? evaluateRule0(row, r0, allowGapCrossing) : null
    const matched = [...(r0ev && r0ev.outcome === 'would_merge' ? ['rule0'] : []), ...act.filter((r) => r.conditions.every((c) => conditionHolds(c, row))).map((r) => r.id)]
    const ev: SegmentEval = { start: row.start, end: row.end, matched, claimedBy: null, side: null, score: null, blocked: null, rule0: r0ev }
    if (!matched.length) return ev
    if (matched[0] === 'rule0' && r0ev) { ev.claimedBy = 'rule0'; ev.side = r0ev.side; ev.score = r0ev.score; ev.blocked = r0ev.blocked; return ev }
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

export function evaluateForSources(hatchScores: Record<string, HatchScoresResult>, sources: string[], chrom: string, rules: AnyRule[], allowGapCrossing: boolean): Record<string, SegmentEval[]> {
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
