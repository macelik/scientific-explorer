import { useIntegration } from '../integrationStore'
import type { Estimator, HatchMetric } from '../hatchMerge'
import { describeRule0, rule0FromSeg, rule0Problems, CONDITION_KINDS, HATCH_PATTERNS, defaultMergeInto, evaluateForSources, makeRule, newCondition, rule0CaseCounts, ruleFlagCounts, ruleProblems, type Condition, type ConditionKind, type HatchPattern, type MergeRule, type PresetName } from '../mergeRules'

const KIND_LABEL: Record<ConditionKind, string> = {
  weak_flank: 'weaker flank below threshold', sign_mismatch: 'flank signs disagree',
  ambiguous: 'SRD and Δlog2FC disagree on the weaker flank', max_size: 'segment size at most',
}
const PRESET_OPTIONS: [PresetName, string][] = [
  ['blank', 'blank rule'], ['weak_flank', 'preset: weak flank (old merge proposal)'],
  ['sign_mismatch', 'preset: flank signs disagree (old transition zone)'], ['ambiguous', 'preset: ambiguous flank preference'],
  ['small', 'preset: small segment'],
]
const PATTERN_LABEL: Record<HatchPattern, string> = { '/': '╱ diagonal', '\\': '╲ back-diagonal', x: '× cross', '-': '― horizontal', '|': '| vertical', '.': '· dots' }
const ESTIMATOR_OPTIONS = <><option value="mean">mean</option><option value="median">median</option><option value="iqr_mean">two-sided IQR mean</option></>

function MetricSelect({ label, metric, estimator, onChange }: { label: string; metric: HatchMetric; estimator: Estimator; onChange: (m: HatchMetric, e: Estimator) => void }) {
  return <>
    <select aria-label={`${label} metric`} value={metric} onChange={(e) => onChange(e.target.value as HatchMetric, estimator)}>
      <option value="srd">SRD</option><option value="srd_phi">SRD/√φ</option><option value="delta_log2fc">Δlog2FC</option>
    </select>
    {metric === 'delta_log2fc' && <select aria-label={`${label} estimator`} value={estimator} onChange={(e) => onChange(metric, e.target.value as Estimator)}>{ESTIMATOR_OPTIONS}</select>}
  </>
}

function ConditionRow({ label, c, onChange, onRemove }: { label: string; c: Condition; onChange: (c: Condition) => void; onRemove: () => void }) {
  return <div className="rule-cond">
    <select aria-label={`${label} kind`} value={c.kind} onChange={(e) => onChange(newCondition(e.target.value as ConditionKind))}>
      {CONDITION_KINDS.map((k) => <option key={k} value={k}>{KIND_LABEL[k]}</option>)}
    </select>
    {c.kind === 'weak_flank' && <>
      <MetricSelect label={label} metric={c.metric} estimator={c.estimator} onChange={(m, e) => onChange({ ...c, metric: m, estimator: e })} />
      <span>&lt;</span>
      <input aria-label={`${label} threshold`} type="number" step={0.1} value={Number.isFinite(c.threshold) ? c.threshold : ''}
        onChange={(e) => onChange({ ...c, threshold: e.target.value === '' ? NaN : +e.target.value })} style={{ width: 72 }} />
    </>}
    {c.kind === 'sign_mismatch' && <MetricSelect label={label} metric={c.metric} estimator={c.estimator} onChange={(m, e) => onChange({ ...c, metric: m, estimator: e })} />}
    {c.kind === 'ambiguous' && <>
      <select aria-label={`${label} SRD variant`} value={c.srdMetric} onChange={(e) => onChange({ ...c, srdMetric: e.target.value as 'srd' | 'srd_phi' })}><option value="srd">SRD</option><option value="srd_phi">SRD/√φ</option></select>
      <span>vs Δlog2FC</span>
      <select aria-label={`${label} estimator`} value={c.estimator} onChange={(e) => onChange({ ...c, estimator: e.target.value as Estimator })}>{ESTIMATOR_OPTIONS}</select>
    </>}
    {c.kind === 'max_size' && <>
      <input aria-label={`${label} max bins`} type="number" min={1} step={1} value={Number.isFinite(c.maxBins) ? c.maxBins : ''}
        onChange={(e) => onChange({ ...c, maxBins: e.target.value === '' ? NaN : +e.target.value })} style={{ width: 72 }} />
      <span>bins</span>
    </>}
    <button className="btn-xs" aria-label={`Remove ${label}`} onClick={onRemove}>×</button>
  </div>
}

export default function MergeRulesEditor() {
  const ig = useIntegration()
  const rules = ig.mergeRules
  const rule0 = rule0FromSeg(ig.seg)
  const r0Issues = rule0Problems(rule0)
  const evals = evaluateForSources(ig.hatchScores, ig.seg.sources, ig.seg.chrom, [rule0, ...rules], ig.allowGapCrossing)
  const counts = ruleFlagCounts(evals)
  const r0c = rule0CaseCounts(evals)
  const scoresLoaded = ig.seg.sources.some((s) => ig.hatchScores[`${s}|${ig.seg.chrom}`])
  const update = (id: string, patch: Partial<MergeRule>) => ig.setMergeRules(rules.map((r) => (r.id === id ? { ...r, ...patch } : r)))
  const setConditions = (r: MergeRule, conditions: Condition[]) =>
    update(r.id, r.conditions.length === 0 ? { conditions, mergeInto: defaultMergeInto(conditions) } : { conditions })
  const move = (i: number, d: -1 | 1) => {
    const j = i + d; if (j < 0 || j >= rules.length) return
    const next = [...rules]; [next[i], next[j]] = [next[j], next[i]]; ig.setMergeRules(next)
  }
  const duplicate = (r: MergeRule, i: number) => {
    const copy: MergeRule = { ...makeRule('blank', rules), name: `${r.name} (copy)`, enabled: r.enabled,
      conditions: r.conditions.map((c) => ({ ...c })), mergeInto: { ...r.mergeInto }, hatch: { ...r.hatch } }
    ig.setMergeRules([...rules.slice(0, i + 1), copy, ...rules.slice(i + 1)])
  }
  return <div className="merge-rules">
    <div className="controls-row">
      <select aria-label="Add rule" value="" onChange={(e) => { if (e.target.value) ig.setMergeRules([...rules, makeRule(e.target.value as PresetName, rules)]) }}>
        <option value="">+ Add rule…</option>
        {PRESET_OPTIONS.map(([k, l]) => <option key={k} value={k}>{l}</option>)}
      </select>
      <span className="muted small">A segment is flagged when ANY enabled rule matches; a rule matches when ALL its conditions hold. Rule order is merge priority.</span>
    </div>
    {!rules.length && <p className="muted small">No additional rules. Rule 0 below already defines what merges; add a rule to merge more segments.</p>}
    <div className="rule-card rule0-card">
      <div className="rule-head"><b>Rule 0 — shown on the tracks</b><span className="muted small">pinned first · edit it in section 1 above</span></div>
      <div className="small">{describeRule0(rule0)}</div>
      <div className="small" role="status" aria-label="Rule 0 status">
        {r0Issues.length ? <span className="notice small">{r0Issues.join('; ')}</span>
          : ig.hatchScoresError ? 'scores unavailable' : !scoresLoaded ? 'loading scores…'
          : `→ would merge now: ${r0c.wouldMerge} segment(s) on ${ig.seg.chrom} across the displayed rows${r0c.gap ? ` · blocked by gap: ${r0c.gap}` : ''}`}
      </div>
    </div>
    {rules.map((r, i) => {
      const n = i + 1, problems = ruleProblems(r)
      return <div key={r.id} className="rule-card">
        <div className="rule-head">
          <label className="chk-inline"><input type="checkbox" aria-label={`Rule ${n} enabled`} checked={r.enabled} onChange={(e) => update(r.id, { enabled: e.target.checked })} /><b>{n}</b></label>
          <input aria-label={`Rule ${n} name`} value={r.name} onChange={(e) => update(r.id, { name: e.target.value })} style={{ width: 240 }} />
          <label>hatch<select aria-label={`Rule ${n} pattern`} value={r.hatch.pattern} onChange={(e) => update(r.id, { hatch: { ...r.hatch, pattern: e.target.value as HatchPattern } })}>{HATCH_PATTERNS.map((p) => <option key={p} value={p}>{PATTERN_LABEL[p]}</option>)}</select></label>
          <input aria-label={`Rule ${n} colour`} type="color" value={r.hatch.color} onChange={(e) => update(r.id, { hatch: { ...r.hatch, color: e.target.value } })} />
          <span className="spacer" />
          <button className="btn-xs" aria-label={`Move rule ${n} up`} disabled={i === 0} onClick={() => move(i, -1)}>↑</button>
          <button className="btn-xs" aria-label={`Move rule ${n} down`} disabled={i === rules.length - 1} onClick={() => move(i, 1)}>↓</button>
          <button className="btn-xs" aria-label={`Duplicate rule ${n}`} onClick={() => duplicate(r, i)}>⧉</button>
          <button className="btn-xs" aria-label={`Delete rule ${n}`} onClick={() => ig.setMergeRules(rules.filter((x) => x.id !== r.id))}>🗑</button>
        </div>
        <div className="small muted">ALL of:</div>
        {r.conditions.map((c, k) => <ConditionRow key={k} label={`Rule ${n} condition ${k + 1}`} c={c}
          onChange={(nc) => setConditions(r, r.conditions.map((x, j) => (j === k ? nc : x)))}
          onRemove={() => setConditions(r, r.conditions.filter((_, j) => j !== k))} />)}
        <button className="btn-xs" aria-label={`Add condition to rule ${n}`} onClick={() => setConditions(r, [...r.conditions, newCondition('weak_flank')])}>+ condition</button>
        <div className="rule-cond">
          <span className="small">merge into the flank with the smaller |score| under</span>
          <MetricSelect label={`Rule ${n} merge-into`} metric={r.mergeInto.metric} estimator={r.mergeInto.estimator} onChange={(m, e) => update(r.id, { mergeInto: { metric: m, estimator: e } })} />
        </div>
        <div className="small" role="status" aria-label={`Rule ${n} status`}>
          {problems.length ? <span className="notice small">{problems.join('; ')}</span>
                        : r.enabled && ig.hatchScoresError ? 'scores unavailable'
            : r.enabled && !scoresLoaded ? 'loading scores…'
            : r.enabled ? `→ flags ${counts[r.id] || 0} segment(s) on ${ig.seg.chrom} across the displayed rows` : 'disabled'}
        </div>
      </div>
    })}
  </div>
}
