// tests/merge_rules_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs
//           node tests/_bundle/merge_rules_helpers.mjs
import fixture from './fixtures/merge_rules_parity.json'
import { canRunMerge, criterionSideLine, defaultMergeInto, evaluateRules, makeRule, makeRule0Copy, migrateHatchLayers, ruleProblems, rule0FromSeg, rule0Problems, rule0CaseCounts, describeRule0, type AnyRule, type MergeRule, type Rule0 } from '../web/src/mergeRules'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}
const close = (a: number | null, b: number | null) => (a === null || b === null ? a === b : Math.abs(a - b) < 1e-9)

const fx = fixture as any
for (const c of fx.cases) {
  const got = evaluateRules(fx.rows, c.rules as AnyRule[], c.allowGapCrossing)
  eq(got.length, c.expected.length, `parity ${c.name} row count`)
  c.expected.forEach((exp: any, i: number) => {
    const g = got[i]
    const same = g.start === exp.start && g.end === exp.end && JSON.stringify(g.matched) === JSON.stringify(exp.matched) &&
      g.claimedBy === exp.claimedBy && g.phase === exp.phase && g.side === exp.side && g.blocked === exp.blocked && close(g.score, exp.score)
    if (!same) { failures++; console.error(`FAIL parity ${c.name} row ${i}: got ${JSON.stringify(g)}, expected ${JSON.stringify(exp)}`) }
    else console.log(`ok parity ${c.name} row ${i}`)
    if (exp.rule0) {
      const r = g.rule0
      const r0same = !!r && r.evaluated === exp.rule0.evaluated && r.passed === exp.rule0.passed && r.signMismatch === exp.rule0.signMismatch &&
        r.flankConflict === exp.rule0.flankConflict && JSON.stringify(r.heldBack) === JSON.stringify(exp.rule0.heldBack) &&
        r.outcome === exp.rule0.outcome && r.side === exp.rule0.side && r.blocked === exp.rule0.blocked && close(r.score, exp.rule0.score)
      if (!r0same) { failures++; console.error(`FAIL parity ${c.name} row ${i} rule0: got ${JSON.stringify(r)}, expected ${JSON.stringify(exp.rule0)}`) }
      else console.log(`ok parity ${c.name} row ${i} rule0`)
    }
  })
}

const legacy = migrateHatchLayers({
  transition: { on: true, metric: 'delta_log2fc', estimator: 'median' },
  proposal: { on: false, metric: 'srd_phi', estimator: 'mean', threshold: 2.5 },
  ambiguous: { on: true, estimator: 'iqr_mean' },
})
eq(legacy.map((r) => [r.id, r.name, r.enabled, r.hatch.pattern]), [
  ['legacy-transition', 'transition zone (from saved session)', true, '/'],
  ['legacy-proposal', 'merge proposal (from saved session)', false, 'x'],
  ['legacy-ambiguous', 'ambiguous flank preference (from saved session)', true, '\\'],
], 'migration ids, names, on/off, patterns')
eq(legacy[0].conditions, [{ kind: 'sign_mismatch', metric: 'delta_log2fc', estimator: 'median' }], 'migration transition condition')
eq(legacy[0].mergeInto, { metric: 'delta_log2fc', estimator: 'median' }, 'migration transition mergeInto')
eq(legacy[1].conditions, [{ kind: 'weak_flank', metric: 'srd_phi', estimator: 'mean', threshold: 2.5 }], 'migration proposal keeps threshold')
eq(legacy[2].conditions, [{ kind: 'ambiguous', srdMetric: 'srd_phi', estimator: 'iqr_mean' }], 'migration ambiguous condition')
eq(legacy[2].mergeInto, { metric: 'srd_phi', estimator: 'iqr_mean' }, 'migration ambiguous mergeInto')
eq(migrateHatchLayers(undefined), [], 'no legacy layers -> no rules')

const base = makeRule('weak_flank', [])
eq(ruleProblems(base), [], 'preset rule is complete')
eq(ruleProblems({ ...base, conditions: [] }), ['add at least one condition'], 'empty rule')
eq(ruleProblems({ ...base, conditions: [{ kind: 'weak_flank', metric: 'srd_phi', estimator: 'mean', threshold: 0 }] }).length, 1, 'zero threshold')
eq(ruleProblems({ ...base, conditions: [{ kind: 'weak_flank', metric: 'srd_phi', estimator: 'mean', threshold: NaN }] }).length, 1, 'cleared threshold box (NaN)')
eq(ruleProblems({ ...base, conditions: [{ kind: 'max_size', maxBins: 2.5 }] }).length, 1, 'fractional size limit')
const r0 = rule0FromSeg({ criteria: ['delta_mean', 'srd_phi'], combine: 'all', tau: 0.2, srdRef: 3.3, smallMax: 100, mergeSignMismatch: false, mergeSignsOnly: false })
eq(canRunMerge(r0, []), true, 'Rule 0 alone can run')
eq(canRunMerge(r0, [{ ...base, enabled: false }]), true, 'disabled user rules do not block')
eq(canRunMerge(r0, [base, { ...base, id: 'b', conditions: [] }]), false, 'an incomplete enabled rule blocks the run')
eq(canRunMerge(r0, [base, { ...base, id: 'b', enabled: false, conditions: [] }]), true, 'an incomplete disabled rule does not block')
eq(canRunMerge({ ...r0, criteria: [] }, []), false, 'Rule 0 with no criterion blocks the run')
eq(evaluateRules(fx.rows, [{ ...base, conditions: [{ kind: 'weak_flank', metric: 'srd_phi', estimator: 'mean', threshold: NaN }] }], false).map((e) => e.matched),
  [[], [], [], [], []], 'an incomplete rule never hatches')
eq(defaultMergeInto([{ kind: 'max_size', maxBins: 3 }]), { metric: 'srd_phi', estimator: 'mean' }, 'size-only rule defaults to SRD/√φ')
eq(defaultMergeInto([{ kind: 'max_size', maxBins: 3 }, { kind: 'sign_mismatch', metric: 'delta_log2fc', estimator: 'median' }]),
  { metric: 'delta_log2fc', estimator: 'median' }, 'first metric condition sets mergeInto')
const r1 = makeRule('blank', []), r2 = makeRule('blank', [r1])
eq(r1.id !== r2.id, true, 'new rule ids are unique')
eq([makeRule('weak_flank', []).name, makeRule('sign_mismatch', []).name, makeRule('ambiguous', []).name, makeRule('small', []).name],
  ['weak SRD/√φ flank', 'flank signs disagree', 'ambiguous flank preference', 'small segment'], 'preset names')

eq([r0.id, r0.type, r0.maxBins, r0.hatch.pattern], ['rule0', 'flank_threshold', 100, 'x'], 'rule0FromSeg shape')
eq(rule0Problems(r0), [], 'default Rule 0 is valid')
eq(rule0Problems({ ...r0, criteria: [] }).length, 1, 'Rule 0 needs a criterion')
eq(rule0Problems({ ...r0, tau: NaN }).length, 1, 'cleared τ is a problem')
eq(rule0Problems({ ...r0, srdRef: 0 }).length, 1, 'zero SRD/√φ ref is a problem')
eq(rule0Problems({ ...r0, maxBins: 2.5 }).length, 1, 'fractional size limit is a problem')
const fx2 = fixture as any
eq(evaluateRules(fx2.rows, [{ ...r0, tau: NaN }], false).map((e) => e.rule0), [null, null, null, null, null], 'an invalid Rule 0 is ignored')
const anyCase = fx2.cases.find((c: any) => c.name === 'rule0_any_signs_allowed_gap_blocks')
const heldCase = fx2.cases.find((c: any) => c.name === 'rule0_all_mean_srd_held_back_by_signs')
eq(rule0CaseCounts({ a: evaluateRules(fx2.rows, anyCase.rules, false) }), { wouldMerge: 3, signs: 3, heldBack: 0, signsOnly: 0, gap: 1, signsMerge: 0, signsMergeBlocked: 0 }, 'case counts (any, signs allowed)')
eq(rule0CaseCounts({ a: evaluateRules(fx2.rows, heldCase.rules, false) }), { wouldMerge: 0, signs: 0, heldBack: 2, signsOnly: 1, gap: 0, signsMerge: 0, signsMergeBlocked: 0 }, 'case counts (held back)')
eq(describeRule0(r0), 'all of: min |Δlog2FC mean| flank < 0.2, min |SRD/√φ| flank < 3.3 · segments ≤ 100 bins · signs disagree: held back · fails Rule 0 with signs disagreeing: not merged', 'describeRule0')

const copy = makeRule0Copy({ criteria: ['srd_phi'], combine: 'all', tau: 0.2, srdRef: 1.0, smallMax: 100, mergeSignMismatch: true, mergeSignsOnly: true }, [])
eq([copy.type, copy.name, copy.enabled, copy.maxBins, copy.hatch.pattern, 'mergeSignsOnly' in copy, copy.id !== 'rule0'],
  ['flank_threshold', 'Rule 0 copy', true, 100, 'x', false, true], 'makeRule0Copy shape')
eq(canRunMerge(r0, [{ ...copy, tau: NaN }]), false, 'enabled copy with cleared τ blocks the run')
eq(canRunMerge(r0, [{ ...copy, tau: NaN, enabled: false }]), true, 'disabled invalid copy does not block')
eq(evaluateRules(fx.rows, [r0, { ...copy, srdRef: NaN }], false).every((e) => !e.matched.includes(copy.id)), true, 'invalid copy never hatches')
const row0 = fx.rows[0]
eq(criterionSideLine(row0, { criterion: 'delta_mean', side: 'left', value: 0.07, threshold: 0.2, fraction: 0.35, passed: true }),
  'min |Δlog2FC mean| flank < 0.2: L 0.070 ✓ · R −0.900 ✗ → left (fraction 0.35) pass', 'per-side tooltip line')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) }
console.log('all merge_rules_helpers checks passed')
