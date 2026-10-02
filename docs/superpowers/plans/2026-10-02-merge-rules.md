# Merge Rules Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the three fixed hatch layers with an ordered list of user-defined merge rules (OR of ANDs, per-rule hatching), and split the segmentation tab into a display-only Visualization section and a Merge rules section.

**Architecture:** Python stays the authority for the merge: a new `server/merge_rules.py` validates rules, evaluates them on per-segment flank-score rows, and runs the one-merge-per-step loop. A TypeScript port (`web/src/mergeRules.ts`) evaluates the same rules on the rows `hatch-scores` already returns, so hatches update live; one shared JSON fixture keeps the two evaluators in agreement. The frontend store holds `mergeRules` + `allowGapCrossing` in place of `hatchLayers`, migrating old sessions.

**Tech Stack:** Python 3.10, FastAPI, NumPy, pytest; React 18 + TypeScript + Vite + Zustand + Plotly.js; esbuild + Node for TS unit tests; Playwright (CommonJS script) for the browser walkthrough.

**Spec:** `docs/superpowers/specs/2026-10-02-merge-rules-design.md` — read it alongside this plan.

## Global Constraints

- Work from the repository root `scientific-explorer/`. Python commands assume `source .venv-scientific-explorer/bin/activate`.
- Do not change how scores are computed: `signed_srd`, `pooled_phi`, `flank_score`, `delta_log2fc`, `segment_summary`, `flank_view.summarize`, and the `POST /api/integration/hatch-scores` response stay byte-for-byte the same.
- No new dependencies (Python or npm). Rule reordering uses ↑/↓ buttons, not drag-and-drop.
- Exact semantics (spec §2): `weak_flank` is `min(|L|,|R|) < threshold` (strict) and needs both L and R finite; `sign_mismatch` is `L*R < 0` (zero is not a mismatch); `ambiguous` is false on a tie in either metric; `max_size` is `end - start <= maxBins`; the merge direction is the flank with the smaller |score| under the claiming rule's `mergeInto`, ties go left; non-finite `mergeInto` score → flagged but `blocked: 'nonfinite'`; gap on the chosen side with gap crossing off → `blocked: 'gap'`, never redirected; only internal segments are evaluated.
- Merge order (spec §4): claiming-rule index, then weaker-flank |score|, then segment start bin. φ is recomputed every step.
- JSON field names are camelCase as in the spec (`mergeInto`, `srdMetric`, `maxBins`, `claimedBy` on the TS side). Python evaluation results use snake_case keys (`claimed_by`); the parity test maps between them.
- Hatch patterns are exactly `'/'`, `'\\'`, `'x'`, `'-'`, `'|'`, `'.'`.
- UI copy calls rules exploratory diagnostics. Never write "validated", "proves", or "ground truth" about them.
- **Ruling (spec §5 deviation):** the API accepts `rules` as plain JSON objects and all validation lives in `merge_rules.validate_rules`, not in a Pydantic discriminated union. Reason: the spec requires disabled rules to be ignored, and a strict Pydantic model would reject a half-edited disabled rule. Cost if wrong: moving validation into Pydantic later is mechanical.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. **A threshold box cleared or set to 0** → no crash; the rule card shows a problem, the rule never hatches, and "Run merge" is disabled. Tests: Task 3 (`ruleProblems`/`canRunMerge`/incomplete-rule-never-hatches), Task 5 (browser).
2. **φ undefined (flat or all-zero profile)** → SRD/√φ conditions are false, an SRD/√φ `mergeInto` gives `blocked: 'nonfinite'`, and the merge ends with zero steps instead of looping. Test: Task 2 `test_run_rule_merge_terminates_when_every_flag_is_nonfinite`.
3. **A chromosome with no internal segment (≤ 2 segments)** → empty evaluation and a merge result with only step 0, no error. Test: Task 2 `test_run_rule_merge_no_internal_segments`.
4. **Every flagged segment is gap-blocked** → the merge ends without redirecting to the other neighbour. Test: Task 2 `test_run_rule_merge_gap_blocks_without_redirect`.
5. **An old saved session (only `hatchLayers`, or neither key)** → migrated preset rules or an empty list, no crash. Tests: Task 3 migration checks, Task 5 legacy-session browser check.

---

## File Structure

- Modify `server/hatch_merge.py` — extract `score_rows()`; remove the old layer-based loop.
- Create `server/merge_rules.py` — rule validation, condition evaluation, rule claiming/direction, `run_rule_merge()`.
- Modify `server/hatch_merge_api.py` — new `hatch-merge` request shape, calls `merge_rules`.
- Create `tests/test_merge_rules.py` — Python unit, loop, route and parity tests.
- Modify `tests/test_hatch_merge.py` — drop tests of the removed loop and old route shape.
- Create `tests/fixtures/merge_rules_parity.json` — shared Python/TS evaluator fixture.
- Create `web/src/mergeRules.ts` — TS types, evaluator, presets, migration, descriptions.
- Create `tests/merge_rules_helpers.ts` — TS unit + parity checks (esbuild + node).
- Modify `web/src/integration.ts` — `hatchShapes` gains `'-'`, `'|'`, `'.'`.
- Modify `web/src/hatchMerge.ts` — request/step/settings types for rules.
- Modify `web/src/integrationStore.ts` — `mergeRules`, `allowGapCrossing`, migration on restore.
- Create `web/src/components/MergeRulesEditor.tsx` — rule cards.
- Modify `web/src/components/HatchMergeView.tsx` — section 2 (editor + run + results).
- Modify `web/src/components/IntegrationSegmentationView.tsx` — Visualization heading, rule hatches, legend, tooltip.
- Modify `web/src/styles.css` — rule-card styles.
- Rewrite `tests/browser_hatch_merge_walkthrough.cjs`; append to `validation/hatch-merge/HANDOFF.md`.
- Modify docs: `README.md`, `docs/USAGE.md`, `docs/SCIENTIFIC_WORKFLOW.md`, `docs/ARCHITECTURE.md`, `tests/README.md`, `AGENTS.md`, `CONTRIBUTING.md`, `server/README.md`, `web/README.md`.

---

### Task 1: Python rule evaluator and validation

**Files:**
- Modify: `server/hatch_merge.py:181-204` (extract `score_rows`)
- Create: `server/merge_rules.py`
- Test: `tests/test_merge_rules.py`

**Interfaces:**
- Consumes: `hm.flat_partition`, `hm.segment_summary`, `hm.pooled_phi`, `hm.flank_score`, `hm.delta_log2fc`, `hm.is_gap`, `hm.classify_transition`, `hm.classify_ambiguous`, `hm.METRICS`, `hm.ESTIMATORS` (all existing in `server/hatch_merge.py`).
- Produces:
  - `hm.score_rows(profile, segs, var_start, var_end, chrom_offset) -> tuple[float | None, list[dict]]` — `segs` is a list of `(start, end)`; rows are exactly the dicts `hatch_scores` returns.
  - `mr.validate_rules(rules: list[dict]) -> list[dict]` — returns enabled rules in order; raises `ValueError`.
  - `mr.row_score(row: dict, side: str, metric: str, estimator: str = 'mean') -> float | None`
  - `mr.condition_holds(cond: dict, row: dict) -> bool`
  - `mr.evaluate_rules(rows: list[dict], rules: list[dict], allow_gap_crossing: bool = False) -> list[dict]` — `rules` must already be validated/enabled; each result is `{'start', 'end', 'matched': list[str], 'claimed_by': str | None, 'side': 'left' | 'right' | None, 'score': float | None, 'blocked': None | 'gap' | 'nonfinite'}`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_merge_rules.py`:

```python
from pathlib import Path
import json
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from server import hatch_merge as hm
from server import merge_rules as mr

FIXTURE = Path(__file__).resolve().parent / 'fixtures' / 'merge_rules_parity.json'


def _row(start=10, end=20, gap_left=False, gap_right=False, **scores):
    row = {'start': start, 'end': end, 'phi': 1.0, 'gap_left': gap_left, 'gap_right': gap_right}
    for key in ('srd', 'srd_phi', 'delta_mean', 'delta_median', 'delta_iqr_mean'):
        row[f'{key}_left'] = None
        row[f'{key}_right'] = None
    row.update(scores)
    return row


def _rule(rid, conditions, merge_metric='srd_phi', merge_estimator='mean', enabled=True, pattern='/'):
    return {'id': rid, 'name': rid, 'enabled': enabled, 'conditions': conditions,
            'mergeInto': {'metric': merge_metric, 'estimator': merge_estimator},
            'hatch': {'pattern': pattern, 'color': '#000000'}}


WEAK_SRD_PHI = {'kind': 'weak_flank', 'metric': 'srd_phi', 'estimator': 'mean', 'threshold': 1.0}
GOOD = _rule('ok', [WEAK_SRD_PHI])


def test_score_rows_matches_hatch_scores_rows():
    profile = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)])
    var_start = (np.arange(30) * 100).astype(float)
    var_end = var_start + 100
    phi, rows = hm.score_rows(profile, hm.flat_partition([10, 20], 30), var_start, var_end, 0)
    full = hm.hatch_scores(profile, [10, 20], 30, var_start, var_end, 0)
    assert rows == full['rows']
    assert phi == full['phi']


def test_weak_flank_is_strict_and_needs_both_flanks():
    row = _row(srd_phi_left=0.5, srd_phi_right=-3.0)
    assert mr.condition_holds(WEAK_SRD_PHI, row) is True
    assert mr.condition_holds({**WEAK_SRD_PHI, 'threshold': 0.5}, row) is False
    assert mr.condition_holds(WEAK_SRD_PHI, _row(srd_phi_left=0.5, srd_phi_right=None)) is False
    assert mr.condition_holds(WEAK_SRD_PHI, _row(srd_phi_left=float('nan'), srd_phi_right=0.1)) is False


def test_sign_mismatch_zero_same_sign_and_estimator_columns():
    c = {'kind': 'sign_mismatch', 'metric': 'delta_log2fc', 'estimator': 'median'}
    assert mr.condition_holds(c, _row(delta_median_left=0.2, delta_median_right=-0.1)) is True
    assert mr.condition_holds(c, _row(delta_median_left=0.0, delta_median_right=-0.1)) is False
    assert mr.condition_holds(c, _row(delta_median_left=0.2, delta_median_right=0.1)) is False
    assert mr.condition_holds(c, _row(delta_mean_left=0.2, delta_mean_right=-0.1)) is False


def test_ambiguous_condition_uses_selected_srd_variant_and_ties_are_not_ambiguous():
    c = {'kind': 'ambiguous', 'srdMetric': 'srd_phi', 'estimator': 'mean'}
    row = _row(srd_phi_left=1.0, srd_phi_right=3.0, srd_left=3.0, srd_right=1.0,
               delta_mean_left=3.0, delta_mean_right=1.0)
    assert mr.condition_holds(c, row) is True
    assert mr.condition_holds({**c, 'srdMetric': 'srd'}, row) is False
    assert mr.condition_holds(c, _row(srd_phi_left=2.0, srd_phi_right=2.0,
                                       delta_mean_left=3.0, delta_mean_right=1.0)) is False


def test_ambiguous_condition_on_real_scores():
    seg_c = np.array([3.0 + (0.1 if i % 2 == 0 else -0.1) for i in range(30)])
    profile = np.concatenate([np.full(2, 1.0), np.full(2, 2.0), seg_c])
    var_start = (np.arange(34) * 100).astype(float)
    _, rows = hm.score_rows(profile, hm.flat_partition([2, 4], 34), var_start, var_start + 100, 0)
    assert mr.condition_holds({'kind': 'ambiguous', 'srdMetric': 'srd_phi', 'estimator': 'mean'}, rows[0]) is True


def test_max_size_is_inclusive():
    c = {'kind': 'max_size', 'maxBins': 10}
    assert mr.condition_holds(c, _row(start=10, end=20)) is True
    assert mr.condition_holds(c, _row(start=10, end=21)) is False


def test_evaluate_rules_and_within_or_across():
    rows = [_row(start=10, end=20, srd_phi_left=0.5, srd_phi_right=-3.0),
            _row(start=20, end=25, srd_phi_left=2.0, srd_phi_right=-1.5)]
    both = _rule('and', [WEAK_SRD_PHI, {'kind': 'max_size', 'maxBins': 5}])
    assert [e['matched'] for e in mr.evaluate_rules(rows, [both])] == [[], []]
    size = _rule('size', [{'kind': 'max_size', 'maxBins': 5}])
    weak = _rule('weak', [WEAK_SRD_PHI])
    assert [e['matched'] for e in mr.evaluate_rules(rows, [size, weak])] == [['weak'], ['size']]


def test_claiming_rule_is_topmost_and_its_merge_into_sets_direction():
    row = _row(start=20, end=25, srd_phi_left=0.5, srd_phi_right=-3.0,
               delta_mean_left=-2.0, delta_mean_right=0.1)
    first = _rule('first', [{'kind': 'max_size', 'maxBins': 5}], merge_metric='delta_log2fc')
    second = _rule('second', [WEAK_SRD_PHI])
    ev = mr.evaluate_rules([row], [first, second])[0]
    assert ev['matched'] == ['first', 'second'] and ev['claimed_by'] == 'first'
    assert ev['side'] == 'right' and ev['score'] == pytest.approx(0.1) and ev['blocked'] is None
    ev = mr.evaluate_rules([row], [second, first])[0]
    assert ev['claimed_by'] == 'second' and ev['side'] == 'left' and ev['score'] == pytest.approx(0.5)


def test_merge_direction_tie_goes_left():
    ev = mr.evaluate_rules([_row(srd_phi_left=-2.0, srd_phi_right=2.0)],
                           [_rule('r', [{'kind': 'max_size', 'maxBins': 100}])])[0]
    assert ev['side'] == 'left' and ev['score'] == pytest.approx(2.0) and ev['blocked'] is None


def test_nonfinite_merge_into_flags_but_blocks():
    ev = mr.evaluate_rules([_row(srd_phi_left=None, srd_phi_right=1.0)],
                           [_rule('r', [{'kind': 'max_size', 'maxBins': 100}])])[0]
    assert ev['matched'] == ['r'] and ev['claimed_by'] == 'r'
    assert ev['blocked'] == 'nonfinite' and ev['side'] is None and ev['score'] is None


def test_gap_blocks_without_redirect_and_can_be_allowed():
    row = _row(srd_phi_left=0.2, srd_phi_right=5.0, gap_left=True)
    rules = [_rule('r', [WEAK_SRD_PHI])]
    ev = mr.evaluate_rules([row], rules)[0]
    assert ev['side'] == 'left' and ev['blocked'] == 'gap'
    ev = mr.evaluate_rules([row], rules, allow_gap_crossing=True)[0]
    assert ev['side'] == 'left' and ev['blocked'] is None


@pytest.mark.parametrize('rules', [
    [],
    [_rule('off', [WEAK_SRD_PHI], enabled=False)],
    [_rule('empty', [])],
    [_rule('k', [{'kind': 'bogus'}])],
    [_rule('m', [{'kind': 'weak_flank', 'metric': 'bogus', 'estimator': 'mean', 'threshold': 1.0}])],
    [_rule('e', [{'kind': 'sign_mismatch', 'metric': 'delta_log2fc', 'estimator': 'bogus'}])],
    [_rule('s', [{'kind': 'ambiguous', 'srdMetric': 'delta_log2fc', 'estimator': 'mean'}])],
    [_rule('t0', [{**WEAK_SRD_PHI, 'threshold': 0}])],
    [_rule('tn', [{**WEAK_SRD_PHI, 'threshold': float('nan')}])],
    [_rule('ti', [{**WEAK_SRD_PHI, 'threshold': float('inf')}])],
    [_rule('tb', [{**WEAK_SRD_PHI, 'threshold': True}])],
    [_rule('b0', [{'kind': 'max_size', 'maxBins': 0}])],
    [_rule('bf', [{'kind': 'max_size', 'maxBins': 2.5}])],
    [_rule('mi', [WEAK_SRD_PHI], merge_metric='bogus')],
    [_rule('p', [WEAK_SRD_PHI], pattern='#')],
    [GOOD, _rule('ok', [WEAK_SRD_PHI])],
    [{'name': 'no id', 'enabled': True, 'conditions': [WEAK_SRD_PHI]}],
    'not a list',
])
def test_validate_rules_rejects(rules):
    with pytest.raises(ValueError):
        mr.validate_rules(rules)


def test_validate_rules_ignores_disabled_and_returns_enabled_in_order():
    broken_but_off = _rule('off', [], enabled=False)
    second = _rule('two', [{'kind': 'max_size', 'maxBins': 3}])
    assert [r['id'] for r in mr.validate_rules([GOOD, broken_but_off, second])] == ['ok', 'two']


def _equivalence_rows():
    rows = []
    for profile, bounds in [
        (np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)]), [10, 20]),
        (np.concatenate([np.full(10, 8.0), np.full(10, 4.0), np.full(10, 2.0)]), [10, 20]),
        (np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(5, 8.0), np.full(15, 9.0)])
         + np.where(np.arange(40) % 2 == 0, 0.05, -0.05), [10, 20, 25]),
    ]:
        n = len(profile)
        vs = (np.arange(n) * 100).astype(float)
        rows += hm.score_rows(profile, hm.flat_partition(bounds, n), vs, vs + 100, 0)[1]
    return rows


@pytest.mark.parametrize('metric,estimator', [('srd', 'mean'), ('srd_phi', 'mean'), ('delta_log2fc', 'median')])
def test_presets_reproduce_old_layer_classifications(metric, estimator):
    for row in _equivalence_rows():
        l = mr.row_score(row, 'left', metric, estimator)
        r = mr.row_score(row, 'right', metric, estimator)
        assert mr.condition_holds({'kind': 'sign_mismatch', 'metric': metric, 'estimator': estimator}, row) \
            == hm.classify_transition(l, r)
        assert mr.condition_holds({'kind': 'weak_flank', 'metric': metric, 'estimator': estimator, 'threshold': 0.5}, row) \
            == (hm.classify_merge_proposal(l, r, 0.5)['weaker_side'] is not None)
        srd_l, srd_r = mr.row_score(row, 'left', 'srd_phi'), mr.row_score(row, 'right', 'srd_phi')
        fc_l = mr.row_score(row, 'left', 'delta_log2fc', estimator)
        fc_r = mr.row_score(row, 'right', 'delta_log2fc', estimator)
        assert mr.condition_holds({'kind': 'ambiguous', 'srdMetric': 'srd_phi', 'estimator': estimator}, row) \
            == (hm.classify_ambiguous(srd_l, srd_r, fc_l, fc_r) == 'ambiguous')
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest -q tests/test_merge_rules.py`
Expected: collection error `ImportError: cannot import name 'merge_rules' from 'server'`.

- [ ] **Step 3: Extract `score_rows` in `server/hatch_merge.py`**

Replace the whole `hatch_scores` function (lines 181-204) with:

```python
def score_rows(profile, segs, var_start, var_end, chrom_offset):
    """Pooled phi and one row per internal segment of the partition `segs`
    ((start, end) chromosome-local bins): SRD, SRD/sqrt(phi) and delta-log2FC
    (all three estimators) against both flanks, plus per-side genomic-gap flags."""
    summaries = [segment_summary(profile, s, e) for s, e in segs]
    phi = pooled_phi(profile, segs)
    rows = []
    for i in range(1, len(segs) - 1):
        seg, left, right = summaries[i], summaries[i - 1], summaries[i + 1]
        row = {
            'start': seg['start'], 'end': seg['end'], 'phi': phi,
            'gap_left': is_gap(var_start, var_end, chrom_offset + segs[i - 1][1] - 1, chrom_offset + segs[i][0]),
            'gap_right': is_gap(var_start, var_end, chrom_offset + segs[i][1] - 1, chrom_offset + segs[i + 1][0]),
        }
        for m in ('srd', 'srd_phi'):
            row[f'{m}_left'] = flank_score(m, seg, left, phi)
            row[f'{m}_right'] = flank_score(m, seg, right, phi)
        for est in ESTIMATORS:
            row[f'delta_{est}_left'] = delta_log2fc(seg, left, est)
            row[f'delta_{est}_right'] = delta_log2fc(seg, right, est)
        rows.append(row)
    return phi, rows


def hatch_scores(profile, boundaries, n_bins, var_start, var_end, chrom_offset):
    """Per internal segment: SRD, SRD/sqrt(phi) and delta-log2FC (all three
    estimators) against both flanks, plus the shared pooled phi and per-side
    genomic-gap flags. Computed genome-wide from live boundaries/profile, not
    restricted to chromosomes with archived diagnostic tables."""
    segs = flat_partition(boundaries, n_bins)
    phi, rows = score_rows(profile, segs, var_start, var_end, chrom_offset)
    return {'segments': [{'start': s, 'end': e} for s, e in segs], 'phi': phi, 'rows': rows}
```

- [ ] **Step 4: Create `server/merge_rules.py`**

```python
"""User-defined merge rules for the exploratory adjacent-segment merge.

A rule is a list of conditions on one internal segment's left/right flank
scores that must ALL hold; a segment is flagged when ANY enabled rule matches
(OR of ANDs). The topmost matching rule claims the segment, and that rule's
mergeInto metric picks the neighbour to merge into (smaller |score|, ties
left). Rows are the per-segment dicts produced by hatch_merge.score_rows, the
same rows /api/integration/hatch-scores returns, so web/src/mergeRules.ts can
evaluate identical rules client-side (see tests/fixtures/merge_rules_parity.json).

These are exploratory diagnostics, independent of the original
bootstrap/chain-guard merge prototype; nothing here recomputes CN.
"""
from __future__ import annotations

import math

from . import hatch_merge as hm

CONDITION_KINDS = ('weak_flank', 'sign_mismatch', 'ambiguous', 'max_size')
SRD_METRICS = ('srd', 'srd_phi')
HATCH_PATTERNS = ('/', '\\', 'x', '-', '|', '.')


def _number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _check_metric(metric, estimator, where):
    if metric not in hm.METRICS:
        raise ValueError(f'{where}: unknown metric {metric!r}')
    if estimator not in hm.ESTIMATORS:
        raise ValueError(f'{where}: unknown estimator {estimator!r}')


def _validate_condition(c, where):
    if not isinstance(c, dict):
        raise ValueError(f'{where}: each condition must be an object')
    kind = c.get('kind')
    if kind == 'weak_flank':
        _check_metric(c.get('metric'), c.get('estimator', 'mean'), where)
        t = c.get('threshold')
        if not (_number(t) and math.isfinite(t) and t > 0):
            raise ValueError(f'{where}: threshold must be a positive finite number')
    elif kind == 'sign_mismatch':
        _check_metric(c.get('metric'), c.get('estimator', 'mean'), where)
    elif kind == 'ambiguous':
        if c.get('srdMetric') not in SRD_METRICS:
            raise ValueError(f'{where}: srdMetric must be one of {SRD_METRICS}')
        if c.get('estimator', 'mean') not in hm.ESTIMATORS:
            raise ValueError(f'{where}: unknown estimator {c.get("estimator")!r}')
    elif kind == 'max_size':
        n = c.get('maxBins')
        if not (isinstance(n, int) and not isinstance(n, bool) and n >= 1):
            raise ValueError(f'{where}: maxBins must be an integer >= 1')
    else:
        raise ValueError(f'{where}: unknown condition kind {kind!r}')


def validate_rules(rules):
    """Validate and return the enabled rules, in list order.

    Disabled rules are ignored apart from the id checks, so a half-edited
    disabled rule never blocks a run. Raises ValueError on any problem."""
    if not isinstance(rules, list) or not all(isinstance(r, dict) for r in rules):
        raise ValueError('rules must be a list of objects')
    ids = [r.get('id') for r in rules]
    if any(not isinstance(i, str) or not i for i in ids):
        raise ValueError('every rule needs a non-empty string id')
    if len(set(ids)) != len(ids):
        raise ValueError('rule ids must be unique')
    enabled = [r for r in rules if r.get('enabled')]
    if not enabled:
        raise ValueError('At least one merge rule must be enabled')
    for r in enabled:
        where = f'rule {r["id"]!r}'
        conditions = r.get('conditions')
        if not isinstance(conditions, list) or not conditions:
            raise ValueError(f'{where}: needs at least one condition')
        for c in conditions:
            _validate_condition(c, where)
        merge_into = r.get('mergeInto') if isinstance(r.get('mergeInto'), dict) else {}
        _check_metric(merge_into.get('metric'), merge_into.get('estimator', 'mean'), f'{where} mergeInto')
        hatch = r.get('hatch') if isinstance(r.get('hatch'), dict) else {}
        if hatch.get('pattern') not in HATCH_PATTERNS:
            raise ValueError(f'{where}: unknown hatch pattern {hatch.get("pattern")!r}')
    return enabled


def row_score(row, side, metric, estimator='mean'):
    """The finite signed flank score for one side of a score row, else None."""
    key = f'delta_{estimator}_{side}' if metric == 'delta_log2fc' else f'{metric}_{side}'
    v = row.get(key)
    return float(v) if v is not None and math.isfinite(v) else None


def condition_holds(c, row):
    kind = c['kind']
    if kind == 'max_size':
        return (row['end'] - row['start']) <= c['maxBins']
    estimator = c.get('estimator', 'mean')
    if kind == 'ambiguous':
        return hm.classify_ambiguous(
            row_score(row, 'left', c['srdMetric']), row_score(row, 'right', c['srdMetric']),
            row_score(row, 'left', 'delta_log2fc', estimator), row_score(row, 'right', 'delta_log2fc', estimator),
        ) == 'ambiguous'
    left = row_score(row, 'left', c['metric'], estimator)
    right = row_score(row, 'right', c['metric'], estimator)
    if kind == 'sign_mismatch':
        return hm.classify_transition(left, right)
    if kind == 'weak_flank':
        return left is not None and right is not None and min(abs(left), abs(right)) < c['threshold']
    raise ValueError(f'unknown condition kind {kind!r}')


def evaluate_rules(rows, rules, allow_gap_crossing=False):
    """Per score row: which rules match, which claims it, and the merge
    direction/score or why it cannot merge. `rules` must be validated."""
    out = []
    for row in rows:
        matched = [r['id'] for r in rules if all(condition_holds(c, row) for c in r['conditions'])]
        ev = {'start': row['start'], 'end': row['end'], 'matched': matched,
              'claimed_by': None, 'side': None, 'score': None, 'blocked': None}
        if matched:
            rule = next(r for r in rules if r['id'] == matched[0])
            ev['claimed_by'] = rule['id']
            metric = rule['mergeInto']['metric']
            estimator = rule['mergeInto'].get('estimator', 'mean')
            left = row_score(row, 'left', metric, estimator)
            right = row_score(row, 'right', metric, estimator)
            if left is None or right is None:
                ev['blocked'] = 'nonfinite'
            else:
                side = 'left' if abs(left) <= abs(right) else 'right'
                ev['side'] = side
                ev['score'] = abs(left) if side == 'left' else abs(right)
                if not allow_gap_crossing and row[f'gap_{side}']:
                    ev['blocked'] = 'gap'
        out.append(ev)
    return out
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python3 -m pytest -q tests/test_merge_rules.py tests/test_hatch_merge.py`
Expected: all pass (the old `run_merge` tests in `test_hatch_merge.py` still pass too, since the old loop is untouched until Task 2).

- [ ] **Step 6: Commit**

```bash
git add server/hatch_merge.py server/merge_rules.py tests/test_merge_rules.py
git commit -m "feat: merge-rule evaluator and validation (OR of ANDs, per-segment)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Rule-driven merge loop and API; remove the layer-based loop

**Files:**
- Modify: `server/merge_rules.py` (append `run_rule_merge`)
- Modify: `server/hatch_merge.py:58-67` (docstring), `:207-353` (delete old loop)
- Modify: `server/hatch_merge_api.py`
- Modify: `tests/test_hatch_merge.py` (delete old loop/route tests)
- Test: `tests/test_merge_rules.py` (append)

**Interfaces:**
- Consumes: `mr.validate_rules`, `mr.evaluate_rules`, `hm.score_rows`, `hm.flat_partition` (Task 1).
- Produces:
  - `mr.run_rule_merge(profile, boundaries, n_bins, var_start, var_end, chrom_offset, rules, allow_gap_crossing=False, max_steps=10000) -> dict` with `{'original', 'final', 'steps', 'settings': {'rules', 'allow_gap_crossing'}}`. Each step: `{'step', 'segments', 'removed_boundary', 'merged_interval', 'segment', 'rule_id', 'rule_name', 'side', 'score', 'phi'}`; step 0 has `None` for everything after `segments`.
  - `POST /api/integration/hatch-merge` body `{chrom, source, rules: list[object], allow_gap_crossing: bool = false}`; 422 on any `ValueError`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_merge_rules.py`:

```python
def _staircase(gap_at=None):
    levels = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(5, 8.0), np.full(15, 9.0)])
    profile = levels + np.where(np.arange(40) % 2 == 0, 0.05, -0.05)
    var_start = (np.arange(40) * 100).astype(float)
    var_end = var_start + 100
    if gap_at is not None:
        var_end[gap_at - 1] += 500
    return profile, [10, 20, 25], 40, var_start, var_end


def _delta_rule(rid, conditions):
    return _rule(rid, conditions, merge_metric='delta_log2fc', merge_estimator='mean')


WEAK_DELTA = {'kind': 'weak_flank', 'metric': 'delta_log2fc', 'estimator': 'mean', 'threshold': 0.5}


def test_run_rule_merge_within_a_rule_smallest_score_first():
    profile, bounds, n, vs, ve = _staircase()
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [_delta_rule('r', [WEAK_DELTA])])
    assert [s['removed_boundary'] for s in out['steps'][1:]] == [10, 25]
    assert out['steps'][1]['segment'] == [10, 20] and out['steps'][1]['side'] == 'left'
    assert out['steps'][1]['score'] == pytest.approx(np.log2(4.2 / 4.0))
    assert out['steps'][1]['rule_id'] == 'r' and out['steps'][1]['rule_name'] == 'r'
    assert out['final'] == [{'start': 0, 'end': 20}, {'start': 20, 'end': 40}]


def test_run_rule_merge_rule_order_beats_score():
    profile, bounds, n, vs, ve = _staircase()
    first = _delta_rule('small-and-weak', [{'kind': 'max_size', 'maxBins': 5}, WEAK_DELTA])
    second = _delta_rule('weak', [WEAK_DELTA])
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [first, second])
    s1, s2 = out['steps'][1], out['steps'][2]
    assert (s1['removed_boundary'], s1['rule_id'], s1['side']) == (25, 'small-and-weak', 'right')
    assert (s2['removed_boundary'], s2['rule_id'], s2['side']) == (10, 'weak', 'left')


def test_run_rule_merge_equal_scores_merge_leftmost_first():
    profile = np.concatenate([np.full(15, 4.0), np.full(5, 20.0)])
    vs = (np.arange(20) * 100).astype(float)
    out = mr.run_rule_merge(profile, [5, 10, 15], 20, vs, vs + 100, 0,
                            [_delta_rule('r', [{**WEAK_DELTA, 'threshold': 1.0}])])
    assert out['steps'][1]['removed_boundary'] == 5


def test_run_rule_merge_recomputes_phi_every_step():
    profile, bounds, n, vs, ve = _staircase()
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [_delta_rule('r', [WEAK_DELTA])])
    for k in (1, 2):
        before = [(s['start'], s['end']) for s in out['steps'][k - 1]['segments']]
        assert out['steps'][k]['phi'] == pytest.approx(hm.pooled_phi(profile, before))
    assert out['steps'][1]['phi'] != pytest.approx(out['steps'][2]['phi'])


def test_run_rule_merge_gap_blocks_without_redirect():
    profile, bounds, n, vs, ve = _staircase(gap_at=10)
    rules = [_delta_rule('r', [{**WEAK_DELTA, 'threshold': 0.1}])]
    blocked = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, rules)
    assert len(blocked['steps']) == 1 and blocked['final'] == blocked['original']
    allowed = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, rules, allow_gap_crossing=True)
    assert [s['removed_boundary'] for s in allowed['steps'][1:]] == [10]


def test_run_rule_merge_no_internal_segments():
    profile, _, n, vs, ve = _staircase()
    out = mr.run_rule_merge(profile, [20], n, vs, ve, 0, [_delta_rule('r', [WEAK_DELTA])])
    assert len(out['steps']) == 1 and out['final'] == out['original']


def test_run_rule_merge_terminates_when_every_flag_is_nonfinite():
    vs = (np.arange(40) * 100).astype(float)
    out = mr.run_rule_merge(np.zeros(40), [10, 20, 30], 40, vs, vs + 100, 0,
                            [_rule('r', [{'kind': 'max_size', 'maxBins': 100}])])
    assert len(out['steps']) == 1


def test_run_rule_merge_echoes_settings_and_step_zero_shape():
    profile, bounds, n, vs, ve = _staircase()
    rules = [_delta_rule('r', [WEAK_DELTA]), _rule('off', [], enabled=False)]
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, rules)
    assert out['settings'] == {'rules': rules, 'allow_gap_crossing': False}
    assert out['steps'][0] == {'step': 0, 'segments': out['original'], 'removed_boundary': None,
                               'merged_interval': None, 'segment': None, 'rule_id': None,
                               'rule_name': None, 'side': None, 'score': None, 'phi': None}
    with pytest.raises(ValueError):
        mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [])


from fastapi import FastAPI
from fastapi.testclient import TestClient


class _FakeStore:
    def __init__(self):
        import pandas as pd
        self.identity = {'hash': 'fake', 'files': []}
        self.var_seq = np.array(['chr1'] * 30)
        self.var_start = (np.arange(30) * 100).astype(float)
        self.var_end = self.var_start + 100
        self.members = {'cluster0': np.arange(5), 'all': np.arange(5)}
        pb = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)])
        self.arrays = {'pb_x': np.tile(pb, (2, 1))}
        self.tables = {'breakpoints': pd.DataFrame({'source': ['cluster0', 'cluster0'],
                                                      'chromosome': ['chr1', 'chr1'],
                                                      'absolute_bin': [10, 20]})}
        self.available = True


def _client():
    from server.hatch_merge_api import register_hatch_merge
    app = FastAPI()
    store = _FakeStore()
    register_hatch_merge(app, lambda: store)
    return TestClient(app)


def test_hatch_merge_route_runs_rules_and_reports_rule_per_step():
    resp = _client().post('/api/integration/hatch-merge', json={
        'chrom': 'chr1', 'source': 'cluster0', 'rules': [_delta_rule('weak', [WEAK_DELTA])]})
    assert resp.status_code == 200
    body = resp.json()
    assert body['final'] == [{'start': 0, 'end': 20}, {'start': 20, 'end': 30}]
    assert body['steps'][1]['rule_id'] == 'weak' and body['steps'][1]['side'] == 'left'
    assert body['settings']['allow_gap_crossing'] is False
    assert body['provenance']['dataset'] == {'hash': 'fake', 'files': []}


@pytest.mark.parametrize('rules', [
    [],
    [_rule('off', [WEAK_DELTA], enabled=False)],
    [_rule('empty', [])],
    [_rule('t', [{**WEAK_DELTA, 'threshold': -1.0}])],
    [_rule('k', [{'kind': 'bogus'}])],
    [_rule('p', [WEAK_DELTA], pattern='#')],
    [_rule('dup', [WEAK_DELTA]), _rule('dup', [WEAK_DELTA])],
])
def test_hatch_merge_route_rejects_bad_rules(rules):
    resp = _client().post('/api/integration/hatch-merge', json={'chrom': 'chr1', 'source': 'cluster0', 'rules': rules})
    assert resp.status_code == 422


def test_hatch_merge_route_requires_rules_field():
    resp = _client().post('/api/integration/hatch-merge', json={'chrom': 'chr1', 'source': 'cluster0'})
    assert resp.status_code == 422
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest -q tests/test_merge_rules.py`
Expected: FAIL with `AttributeError: module 'server.merge_rules' has no attribute 'run_rule_merge'` and the route tests failing (old request shape).

- [ ] **Step 3: Append `run_rule_merge` to `server/merge_rules.py`**

```python
def run_rule_merge(profile, boundaries, n_bins, var_start, var_end, chrom_offset,
                   rules, allow_gap_crossing=False, max_steps=10000):
    """Greedy exploratory merge driven by user-defined rules, one merge per
    step. Each step recomputes every score (pooled phi included) for the
    current partition, evaluates the enabled rules, and merges the candidate
    with the lowest (claiming-rule index, weaker-flank |score|, start bin).
    Gap-blocked and non-finite candidates are skipped, never redirected.
    Stops when no candidate remains."""
    enabled = validate_rules(rules)
    order = {r['id']: k for k, r in enumerate(enabled)}
    names = {r['id']: (r.get('name') or r['id']) for r in enabled}
    segs = hm.flat_partition(boundaries, n_bins)
    original = [{'start': s, 'end': e} for s, e in segs]
    steps = [{'step': 0, 'segments': list(original), 'removed_boundary': None, 'merged_interval': None,
              'segment': None, 'rule_id': None, 'rule_name': None, 'side': None, 'score': None, 'phi': None}]
    for _ in range(max_steps + 1):
        phi, rows = hm.score_rows(profile, segs, var_start, var_end, chrom_offset)
        evals = evaluate_rules(rows, enabled, allow_gap_crossing)
        candidates = [(k + 1, ev) for k, ev in enumerate(evals)
                      if ev['claimed_by'] is not None and ev['blocked'] is None]
        if not candidates:
            break
        i, ev = min(candidates, key=lambda t: (order[t[1]['claimed_by']], t[1]['score'], t[1]['start']))
        a = i - 1 if ev['side'] == 'left' else i
        removed = segs[a][1]
        merged = (segs[a][0], segs[a + 1][1])
        segs = segs[:a] + [merged] + segs[a + 2:]
        steps.append({
            'step': len(steps), 'segments': [{'start': s, 'end': e} for s, e in segs],
            'removed_boundary': removed, 'merged_interval': [merged[0], merged[1]],
            'segment': [ev['start'], ev['end']], 'rule_id': ev['claimed_by'],
            'rule_name': names[ev['claimed_by']], 'side': ev['side'], 'score': ev['score'], 'phi': phi,
        })
    else:
        raise RuntimeError('Merge did not terminate within max_steps')
    return {'original': original, 'final': steps[-1]['segments'], 'steps': steps,
            'settings': {'rules': rules, 'allow_gap_crossing': allow_gap_crossing}}
```

- [ ] **Step 4: Remove the old layer-based loop from `server/hatch_merge.py`**

Delete these functions entirely: `_segment_status`, `_layer_eligibility`, `_eligible_boundaries`, `_validate_layer`, `run_merge` (everything after `hatch_scores`). In the `pooled_phi` docstring, replace `whenever that partition changes (see run_merge below): SRD/sqrt(phi) is` with `whenever that partition changes (see merge_rules.run_rule_merge): SRD/sqrt(phi) is`.

- [ ] **Step 5: Rewrite the merge route in `server/hatch_merge_api.py`**

Replace the imports and models block (lines 4-45) with:

```python
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from fastapi import HTTPException
from pydantic import BaseModel

from . import hatch_merge as hm
from . import merge_rules as mr
from .integration import IntegrationStore

_IMPLEMENTATION_FILES = (Path(__file__).with_name('merge_rules.py'), Path(__file__).with_name('hatch_merge.py'))


class HatchScoresRequest(BaseModel):
    chrom: str
    source: str


class HatchMergeRequest(BaseModel):
    """Rules are validated by merge_rules.validate_rules, not by this model, so
    disabled (possibly half-edited) rules can be sent and are ignored."""
    chrom: str
    source: str
    rules: List[Dict[str, Any]]
    allow_gap_crossing: bool = False
```

Replace the body of `hatch_merge` (the route) with:

```python
    @app.post('/api/integration/hatch-merge')
    def hatch_merge(req: HatchMergeRequest):
        ig = _ig()
        try:
            profile, boundaries, idx = _resolve(ig, req.chrom, req.source)
            result = mr.run_rule_merge(profile, boundaries, len(idx), ig.var_start, ig.var_end,
                                       chrom_offset=int(idx[0]), rules=req.rules,
                                       allow_gap_crossing=req.allow_gap_crossing)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e
        return {**result, 'source': req.source, 'chrom': req.chrom,
                'start_bp': ig.var_start[idx].tolist(), 'end_bp': ig.var_end[idx].tolist(),
                'provenance': {'dataset': ig.identity,
                               'implementation_sha256': hashlib.sha256(b''.join(p.read_bytes() for p in _IMPLEMENTATION_FILES)).hexdigest(),
                               'scope': 'User-defined exploratory merge rules; independent of the original bootstrap/chain-guard '
                                        'merge prototype and the archived SRD screen. No CN recomputation on merged boundaries.'}}
```

- [ ] **Step 6: Delete the obsolete tests in `tests/test_hatch_merge.py`**

Delete every test from `def test_run_merge_rejects_bad_settings` through the end of `test_run_merge_union_semantics_any_enabled_layer_can_trigger` (current lines 141-273), and the four route tests `test_hatch_merge_route_runs_and_reports_provenance`, `test_hatch_merge_route_rejects_bad_threshold`, `test_hatch_merge_route_rejects_no_layer_enabled`, `test_hatch_merge_route_accepts_transition_or_ambiguous_alone` (current lines 319-347). Keep `_toy_profile_and_bounds`, `_FakeStore`, `_client`, and the two `hatch-scores` route tests.

- [ ] **Step 7: Check nothing else uses the removed code**

Run: `grep -rn "run_merge\b\|_layer_eligibility\|small_max_bins\|veto_" server tests/*.py`
Expected: no output.

- [ ] **Step 8: Run the backend suite**

Run: `python3 -m pytest -q tests/test_merge_rules.py tests/test_hatch_merge.py tests/test_experiments.py tests/test_segment_estimators.py tests/test_integration_review.py`
Expected: all pass.

- [ ] **Step 9: Commit**

```bash
git add server/merge_rules.py server/hatch_merge.py server/hatch_merge_api.py tests/test_merge_rules.py tests/test_hatch_merge.py
git commit -m "feat: rule-driven merge loop and hatch-merge API; remove layer-based loop

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: TypeScript rule evaluator and Python/TS parity fixture

**Files:**
- Create: `web/src/mergeRules.ts`
- Create: `tests/fixtures/merge_rules_parity.json`
- Create: `tests/merge_rules_helpers.ts`
- Modify: `tests/test_merge_rules.py` (append parity test)

**Interfaces:**
- Consumes: from `web/src/hatchMerge.ts` (existing): `type HatchScoreRow`, `type HatchScoresResult`, `type HatchMetric`, `type Estimator`, `classifyAmbiguous`, `classifyTransition`, `scoreForMetric`. Python `mr.validate_rules`, `mr.evaluate_rules` (Task 1).
- Produces (exports of `web/src/mergeRules.ts`):
  - types `HatchPattern`, `ConditionKind`, `Condition`, `MergeRule`, `PresetName`, `SegmentEval { start; end; matched: string[]; claimedBy: string | null; side: 'left' | 'right' | null; score: number | null; blocked: null | 'gap' | 'nonfinite' }`
  - constants `HATCH_PATTERNS`, `CONDITION_KINDS`
  - `metricLabel(metric, estimator): string`, `describeCondition(c): string`, `conditionHolds(c, row): boolean`, `conditionDetail(c, row): string`
  - `ruleProblems(rule): string[]`, `activeRules(rules): MergeRule[]`, `canRunMerge(rules): boolean`
  - `evaluateRules(rows, rules, allowGapCrossing): SegmentEval[]`
  - `evaluateForSources(hatchScores: Record<string, HatchScoresResult>, sources: string[], chrom: string, rules, allowGapCrossing): Record<string, SegmentEval[]>`
  - `ruleFlagCounts(evals: Record<string, SegmentEval[]>): Record<string, number>`
  - `newCondition(kind): Condition`, `defaultMergeInto(conditions): { metric; estimator }`, `makeRule(preset: PresetName, existing: MergeRule[]): MergeRule`, `migrateHatchLayers(old: any): MergeRule[]`

- [ ] **Step 1: Create the shared fixture `tests/fixtures/merge_rules_parity.json`**

```json
{
  "rows": [
    {"start": 10, "end": 20, "phi": 1.0, "gap_left": false, "gap_right": false,
     "srd_left": 1.0, "srd_right": -3.0, "srd_phi_left": 0.5, "srd_phi_right": -1.5,
     "delta_mean_left": 0.07, "delta_mean_right": -0.9, "delta_median_left": 0.05, "delta_median_right": -0.8,
     "delta_iqr_mean_left": 0.06, "delta_iqr_mean_right": -0.85},
    {"start": 20, "end": 25, "phi": 1.0, "gap_left": false, "gap_right": true,
     "srd_left": 4.0, "srd_right": -0.4, "srd_phi_left": 2.0, "srd_phi_right": -0.2,
     "delta_mean_left": 0.9, "delta_mean_right": -0.17, "delta_median_left": 0.9, "delta_median_right": -0.2,
     "delta_iqr_mean_left": 0.9, "delta_iqr_mean_right": -0.18},
    {"start": 25, "end": 30, "phi": 1.0, "gap_left": false, "gap_right": false,
     "srd_left": 1.0, "srd_right": 1.2, "srd_phi_left": null, "srd_phi_right": null,
     "delta_mean_left": 0.3, "delta_mean_right": 0.2, "delta_median_left": 0.3, "delta_median_right": 0.2,
     "delta_iqr_mean_left": 0.3, "delta_iqr_mean_right": 0.2},
    {"start": 30, "end": 40, "phi": 1.0, "gap_left": false, "gap_right": false,
     "srd_left": 3.0, "srd_right": -3.0, "srd_phi_left": -2.0, "srd_phi_right": 2.0,
     "delta_mean_left": -0.5, "delta_mean_right": 0.5, "delta_median_left": -0.5, "delta_median_right": 0.5,
     "delta_iqr_mean_left": -0.5, "delta_iqr_mean_right": 0.5},
    {"start": 40, "end": 41, "phi": 1.0, "gap_left": false, "gap_right": false,
     "srd_left": 1.0, "srd_right": 3.0, "srd_phi_left": 1.0, "srd_phi_right": 3.0,
     "delta_mean_left": 3.0, "delta_mean_right": 1.0, "delta_median_left": 3.0, "delta_median_right": 1.0,
     "delta_iqr_mean_left": 3.0, "delta_iqr_mean_right": 1.0}
  ],
  "cases": [
    {
      "name": "weak_srd_phi_gap_blocked",
      "allowGapCrossing": false,
      "rules": [
        {"id": "a", "name": "a", "enabled": true,
         "conditions": [{"kind": "weak_flank", "metric": "srd_phi", "estimator": "mean", "threshold": 1.0}],
         "mergeInto": {"metric": "srd_phi", "estimator": "mean"}, "hatch": {"pattern": "x", "color": "#b91c1c"}}
      ],
      "expected": [
        {"start": 10, "end": 20, "matched": ["a"], "claimedBy": "a", "side": "left", "score": 0.5, "blocked": null},
        {"start": 20, "end": 25, "matched": ["a"], "claimedBy": "a", "side": "right", "score": 0.2, "blocked": "gap"},
        {"start": 25, "end": 30, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null},
        {"start": 30, "end": 40, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null},
        {"start": 40, "end": 41, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null}
      ]
    },
    {
      "name": "sign_mismatch_or_small_gap_allowed",
      "allowGapCrossing": true,
      "rules": [
        {"id": "b1", "name": "b1", "enabled": true,
         "conditions": [{"kind": "sign_mismatch", "metric": "delta_log2fc", "estimator": "median"}],
         "mergeInto": {"metric": "delta_log2fc", "estimator": "median"}, "hatch": {"pattern": "/", "color": "#1d4ed8"}},
        {"id": "b2", "name": "b2", "enabled": true,
         "conditions": [{"kind": "max_size", "maxBins": 5}],
         "mergeInto": {"metric": "srd_phi", "estimator": "mean"}, "hatch": {"pattern": ".", "color": "#047857"}}
      ],
      "expected": [
        {"start": 10, "end": 20, "matched": ["b1"], "claimedBy": "b1", "side": "left", "score": 0.05, "blocked": null},
        {"start": 20, "end": 25, "matched": ["b1", "b2"], "claimedBy": "b1", "side": "right", "score": 0.2, "blocked": null},
        {"start": 25, "end": 30, "matched": ["b2"], "claimedBy": "b2", "side": null, "score": null, "blocked": "nonfinite"},
        {"start": 30, "end": 40, "matched": ["b1"], "claimedBy": "b1", "side": "left", "score": 0.5, "blocked": null},
        {"start": 40, "end": 41, "matched": ["b2"], "claimedBy": "b2", "side": "left", "score": 1.0, "blocked": null}
      ]
    },
    {
      "name": "ambiguous_and_weak_srd_disabled_rule_ignored",
      "allowGapCrossing": false,
      "rules": [
        {"id": "c0", "name": "c0", "enabled": false,
         "conditions": [{"kind": "weak_flank", "metric": "srd", "estimator": "mean", "threshold": 100}],
         "mergeInto": {"metric": "srd", "estimator": "mean"}, "hatch": {"pattern": "|", "color": "#7c3aed"}},
        {"id": "c1", "name": "c1", "enabled": true,
         "conditions": [{"kind": "ambiguous", "srdMetric": "srd_phi", "estimator": "mean"},
                        {"kind": "weak_flank", "metric": "srd", "estimator": "mean", "threshold": 3.5}],
         "mergeInto": {"metric": "delta_log2fc", "estimator": "mean"}, "hatch": {"pattern": "\\", "color": "#b45309"}}
      ],
      "expected": [
        {"start": 10, "end": 20, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null},
        {"start": 20, "end": 25, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null},
        {"start": 25, "end": 30, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null},
        {"start": 30, "end": 40, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null},
        {"start": 40, "end": 41, "matched": ["c1"], "claimedBy": "c1", "side": "right", "score": 1.0, "blocked": null}
      ]
    }
  ]
}
```

- [ ] **Step 2: Append the Python parity test to `tests/test_merge_rules.py`**

```python
def test_python_evaluator_matches_parity_fixture():
    fixture = json.loads(FIXTURE.read_text())
    for case in fixture['cases']:
        got = mr.evaluate_rules(fixture['rows'], mr.validate_rules(case['rules']), case['allowGapCrossing'])
        assert len(got) == len(case['expected']), case['name']
        for g, exp in zip(got, case['expected']):
            assert (g['start'], g['end'], g['matched'], g['claimed_by'], g['side'], g['blocked']) == \
                (exp['start'], exp['end'], exp['matched'], exp['claimedBy'], exp['side'], exp['blocked']), case['name']
            assert (g['score'] is None) == (exp['score'] is None), case['name']
            if exp['score'] is not None:
                assert g['score'] == pytest.approx(exp['score']), case['name']
```

Run: `python3 -m pytest -q tests/test_merge_rules.py::test_python_evaluator_matches_parity_fixture`
Expected: PASS (the Python evaluator already exists; this pins the fixture's expected values).

- [ ] **Step 3: Write the failing TS test `tests/merge_rules_helpers.ts`**

```ts
// tests/merge_rules_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs
//           node tests/_bundle/merge_rules_helpers.mjs
import fixture from './fixtures/merge_rules_parity.json'
import { canRunMerge, defaultMergeInto, evaluateRules, makeRule, migrateHatchLayers, ruleProblems, type MergeRule } from '../web/src/mergeRules'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}
const close = (a: number | null, b: number | null) => (a === null || b === null ? a === b : Math.abs(a - b) < 1e-9)

const fx = fixture as any
for (const c of fx.cases) {
  const got = evaluateRules(fx.rows, c.rules as MergeRule[], c.allowGapCrossing)
  eq(got.length, c.expected.length, `parity ${c.name} row count`)
  c.expected.forEach((exp: any, i: number) => {
    const g = got[i]
    const same = g.start === exp.start && g.end === exp.end && JSON.stringify(g.matched) === JSON.stringify(exp.matched) &&
      g.claimedBy === exp.claimedBy && g.side === exp.side && g.blocked === exp.blocked && close(g.score, exp.score)
    if (!same) { failures++; console.error(`FAIL parity ${c.name} row ${i}: got ${JSON.stringify(g)}, expected ${JSON.stringify(exp)}`) }
    else console.log(`ok parity ${c.name} row ${i}`)
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
eq(canRunMerge([]), false, 'no rules -> cannot run')
eq(canRunMerge([{ ...base, enabled: false }]), false, 'only disabled rules -> cannot run')
eq(canRunMerge([base, { ...base, id: 'b', conditions: [] }]), false, 'an incomplete enabled rule blocks the run')
eq(canRunMerge([base, { ...base, id: 'b', enabled: false, conditions: [] }]), true, 'an incomplete disabled rule does not block')
eq(evaluateRules(fx.rows, [{ ...base, conditions: [{ kind: 'weak_flank', metric: 'srd_phi', estimator: 'mean', threshold: NaN }] }], false).map((e) => e.matched),
  [[], [], [], [], []], 'an incomplete rule never hatches')
eq(defaultMergeInto([{ kind: 'max_size', maxBins: 3 }]), { metric: 'srd_phi', estimator: 'mean' }, 'size-only rule defaults to SRD/√φ')
eq(defaultMergeInto([{ kind: 'max_size', maxBins: 3 }, { kind: 'sign_mismatch', metric: 'delta_log2fc', estimator: 'median' }]),
  { metric: 'delta_log2fc', estimator: 'median' }, 'first metric condition sets mergeInto')
const r1 = makeRule('blank', []), r2 = makeRule('blank', [r1])
eq(r1.id !== r2.id, true, 'new rule ids are unique')
eq([makeRule('weak_flank', []).name, makeRule('sign_mismatch', []).name, makeRule('ambiguous', []).name, makeRule('small', []).name],
  ['weak SRD/√φ flank', 'flank signs disagree', 'ambiguous flank preference', 'small segment'], 'preset names')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) }
console.log('all merge_rules_helpers checks passed')
```

Run: `web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs`
Expected: FAIL — `Could not resolve "../web/src/mergeRules"`.

- [ ] **Step 4: Create `web/src/mergeRules.ts`**

```ts
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

export const activeRules = (rules: MergeRule[]): MergeRule[] => rules.filter((r) => r.enabled && ruleProblems(r).length === 0)

export function canRunMerge(rules: MergeRule[]): boolean {
  const enabled = rules.filter((r) => r.enabled)
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
```

- [ ] **Step 5: Run the TS test to verify it passes**

Run:
```bash
web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs && node tests/_bundle/merge_rules_helpers.mjs
```
Expected: last line `all merge_rules_helpers checks passed`. Also re-run the existing TS helper test to confirm nothing it imports broke:
```bash
web/node_modules/.bin/esbuild tests/hatch_merge_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/hatch_merge_helpers.mjs && node tests/_bundle/hatch_merge_helpers.mjs
```
Expected: `all hatch_merge_helpers checks passed`.

- [ ] **Step 6: Commit** (the `.mjs` bundles are build output; do not commit them)

```bash
git add web/src/mergeRules.ts tests/fixtures/merge_rules_parity.json tests/merge_rules_helpers.ts tests/test_merge_rules.py
git commit -m "feat: TypeScript merge-rule evaluator with shared Python/TS parity fixture

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Frontend — rule state, hatches, editor, and merge panel

**Files:**
- Modify: `web/src/integration.ts:94-99` (`hatchShapes`)
- Modify: `web/src/hatchMerge.ts:24-48`
- Modify: `web/src/integrationStore.ts`
- Create: `web/src/components/MergeRulesEditor.tsx`
- Modify (rewrite): `web/src/components/HatchMergeView.tsx`
- Modify: `web/src/components/IntegrationSegmentationView.tsx`
- Modify: `web/src/styles.css` (append)

**Interfaces:**
- Consumes: everything exported by `web/src/mergeRules.ts` (Task 3); the API from Task 2.
- Produces: store fields `mergeRules: MergeRule[]`, `allowGapCrossing: boolean`, `setMergeRules(rules: MergeRule[])`, `setAllowGapCrossing(on: boolean)`; session keys `integration.mergeRules` and `integration.allowGapCrossing`. Accessible names used by Task 5: `Add rule`; `Rule N enabled|name|pattern|colour|status`; `Rule N condition K kind|threshold|max bins|metric|estimator|SRD variant`; `Move rule N up|down`; `Duplicate rule N`; `Delete rule N`; `Add condition to rule N`; `Rule N merge-into metric|estimator`; `Allow merging across a genomic gap`; button `Run merge`; heading text `2. Merge rules`; tracks heading text `1. Visualization`; history-table header `rule`.

- [ ] **Step 1: Add the new hatch patterns to `hatchShapes` in `web/src/integration.ts`**

Replace lines 94-99 (the whole `hatchShapes` function) with:

```ts
export type HatchStyle = 'none' | '/' | '\\' | 'x' | '-' | '|' | '.'
export function hatchShapes(x0: number, x1: number, yref: string, style: HatchStyle, color = 'rgba(17,24,39,0.55)'): any[] {
  if (style === 'none') return []
  const w = x1 - x0
  const line = (ax0: number, ax1: number, ay0: number, ay1: number, width = 1) =>
    ({ type: 'line', xref: 'x', yref: `${yref} domain`, x0: x0 + ax0 * w, x1: x0 + ax1 * w, y0: ay0, y1: ay1, line: { color, width }, layer: 'above' })
  const fractions = [0.1, 0.3, 0.5, 0.7, 0.9]
  if (style === '-') return fractions.map((y) => line(0, 1, y, y))
  if (style === '|') return fractions.map((x) => line(x, x, 0, 1))
  if (style === '.') return fractions.flatMap((x) => fractions.map((y) => line(x - 0.015, x + 0.015, y, y, 3)))
  const step = 0.22
  const lines = style === 'x' ? [...unitDiag(step, false), ...unitDiag(step, true)] : style === '\\' ? unitDiag(step, true) : unitDiag(step, false)
  return lines.map((l) => line(l.x0, l.x1, l.y0, l.y1))
}
```

- [ ] **Step 2: Update the merge types in `web/src/hatchMerge.ts`**

Add at the top, after the existing `import { hatchShapes } from './integration'` line:

```ts
import type { MergeRule } from './mergeRules'
```

Replace lines 24-48 (from `export interface MergeStep` through the end of `export interface MergeRequest { ... }`) with:

```ts
export interface MergeStep {
  step: number; segments: { start: number; end: number }[]
  removed_boundary: number | null; merged_interval: [number, number] | null
  segment: [number, number] | null; rule_id: string | null; rule_name: string | null
  side: 'left' | 'right' | null; score: number | null; phi: number | null
}
export interface MergeSettings { rules: MergeRule[]; allow_gap_crossing: boolean }
export interface MergeResult {
  original: { start: number; end: number }[]; final: { start: number; end: number }[]
  steps: MergeStep[]; settings: MergeSettings
  source: string; chrom: string; start_bp: number[]; end_bp: number[]; provenance: unknown
}
/** Rules are evaluated server-side by server/merge_rules.py; disabled rules are ignored. */
export interface MergeRequest { chrom: string; source: string; rules: MergeRule[]; allow_gap_crossing: boolean }
```

- [ ] **Step 3: Replace hatch-layer state with merge rules in `web/src/integrationStore.ts`**

1. Replace line 5 with:
```ts
import { fetchHatchScores, type HatchScoresResult } from './hatchMerge'
import { migrateHatchLayers, type MergeRule } from './mergeRules'
```
2. Delete lines 10-19 (`export interface HatchLayerSettings { ... }` and `const defaultHatchLayers ... }`).
3. In `IntegrationState`, replace `hatchLayers: HatchLayerSettings` with `mergeRules: MergeRule[]; allowGapCrossing: boolean`, and replace `setHatchLayers: (p: Partial<HatchLayerSettings>) => void` with:
```ts
  setMergeRules: (rules: MergeRule[]) => void
  setAllowGapCrossing: (on: boolean) => void
```
4. Replace `hatchLayers: defaultHatchLayers, hatchScores: {}, hatchScoresError: null,` with `mergeRules: [], allowGapCrossing: false, hatchScores: {}, hatchScoresError: null,`.
5. Replace `setHatchLayers: (p) => set({ hatchLayers: { ...get().hatchLayers, ...p } }),` with:
```ts
  setMergeRules: (rules) => set({ mergeRules: rules }),
  setAllowGapCrossing: (on) => set({ allowGapCrossing: on }),
```
6. In `serialize`, replace `hatchLayers: s.hatchLayers` with `mergeRules: s.mergeRules, allowGapCrossing: s.allowGapCrossing`.
7. In `restore`, replace `hatchLayers: { ...defaultHatchLayers, ...(p.hatchLayers || {}) }` with `mergeRules: Array.isArray(p.mergeRules) ? p.mergeRules : migrateHatchLayers(p.hatchLayers), allowGapCrossing: !!p.allowGapCrossing`.

- [ ] **Step 4: Create `web/src/components/MergeRulesEditor.tsx`**

```tsx
import { useIntegration } from '../integrationStore'
import type { Estimator, HatchMetric } from '../hatchMerge'
import { CONDITION_KINDS, HATCH_PATTERNS, defaultMergeInto, evaluateForSources, makeRule, newCondition, ruleFlagCounts, ruleProblems, type Condition, type ConditionKind, type HatchPattern, type MergeRule, type PresetName } from '../mergeRules'

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
  const counts = ruleFlagCounts(evaluateForSources(ig.hatchScores, ig.seg.sources, ig.seg.chrom, rules, ig.allowGapCrossing))
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
    {!rules.length && <p className="muted small">No merge rules yet. Add one to hatch segments on the Chromosome tracks and to enable a merge run.</p>}
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
            : r.enabled ? `→ flags ${counts[r.id] || 0} segment(s) on ${ig.seg.chrom} across the displayed rows` : 'disabled'}
        </div>
      </div>
    })}
  </div>
}
```

- [ ] **Step 5: Rewrite `web/src/components/HatchMergeView.tsx`**

Replace the whole file with:

```tsx
import { useRef, useState } from 'react'
import DeferredPlot from './DeferredPlot'
import MergeRulesEditor from './MergeRulesEditor'
import { useIntegration } from '../integrationStore'
import { useStore } from '../store'
import { runHatchMerge, hatchShapes, type MergeResult } from '../hatchMerge'
import { canRunMerge, describeCondition } from '../mergeRules'

function download(text: string, name: string, type: string) {
  const url = URL.createObjectURL(new Blob([text], { type })); const a = document.createElement('a')
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000)
}

export default function HatchMergeView() {
  const ig = useIntegration(), meta = ig.meta!
  const ex = useStore((s) => s.meta)!, presentation = useStore((s) => s.presentation)
  const [source, setSource] = useState('cluster4')
  const [result, setResult] = useState<MergeResult | null>(null)
  const [stepIdx, setStepIdx] = useState(0)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const token = useRef(0)
  const chrom = ig.seg.chrom
  const runnable = canRunMerge(ig.mergeRules)
  const request = { chrom, source, rules: ig.mergeRules, allow_gap_crossing: ig.allowGapCrossing }
  const signature = JSON.stringify(request)
  const [frozen, setFrozen] = useState('')
  const stale = !!result && frozen !== signature

  const run = async () => {
    if (!runnable) return
    const id = ++token.current; setBusy(true); setError('')
    try {
      const data = await runHatchMerge(request)
      if (token.current !== id) return
      setResult(data); setFrozen(signature); setStepIdx(data.steps.length - 1)
    } catch (e: any) { if (token.current === id) setError(String(e.message)) }
    finally { if (token.current === id) setBusy(false) }
  }
  const step = result?.steps[stepIdx]
  // Describe the displayed run from result.settings, never from the live editor state.
  const ruleSummary = result ? result.settings.rules.filter((r) => r.enabled).map((r, i) => `${i + 1}. ${r.name}: ${r.conditions.map(describeCondition).join(' AND ')}`).join(' | ') : ''

  const layout = { height: 260, margin: { l: 60, r: 20, t: 30, b: 45 }, paper_bgcolor: '#fff', plot_bgcolor: '#fff', font: { family: 'system-ui', size: 12 }, hovermode: presentation ? false : 'closest' as const }
  const segmentTrace = (segs: { start: number; end: number }[], startBp: number[], endBp: number[], color: string, name: string) => ({
    type: 'bar', orientation: 'h', name, y: segs.map(() => name), base: segs.map((s) => startBp[s.start] / 1e6),
    x: segs.map((s) => (endBp[s.end - 1] - startBp[s.start]) / 1e6), marker: { color }, hovertemplate: '%{x:.3f} Mb<extra></extra>',
    customdata: segs.map((s) => `${s.start},${s.end}`),
  })
  const shapesForSegments = (segs: { start: number; end: number }[], startBp: number[], endBp: number[]) =>
    segs.map((s, i) => ({ type: 'rect', xref: 'x', yref: 'paper', x0: startBp[s.start] / 1e6, x1: endBp[s.end - 1] / 1e6, y0: 0, y1: 1,
      fillcolor: i % 2 === 0 ? 'rgba(150,73,29,0.18)' : 'rgba(40,120,160,0.18)', line: { width: 1, color: '#334155' } }))

  return <section className="panel hatch-merge">
    <h3>2. Merge rules — exploratory adjacent-segment merge</h3>
    <p className="lead">Build an ordered list of rules. A rule is a set of conditions on a segment's left/right flank scores that must ALL hold; a segment is flagged when ANY enabled rule matches, and each rule's hatch on the Chromosome tracks above shows exactly the segments it flags. A merge run repeatedly joins one flagged segment to its more similar neighbour (the flank with the smaller |score| under that rule's "merge into" metric): topmost rule first, then smallest score, then leftmost, recomputing every score, including the pooled dispersion φ, after each merge. These are user-defined exploratory diagnostics, not validated decision rules. This is not the original bootstrap/chain-guard merge prototype, and CN calls are not recomputed on merged boundaries.</p>
    <MergeRulesEditor />
    <div className="controls-row">
      <label>Source<select aria-label="Merge source" value={source} onChange={(e) => setSource(e.target.value)}>{meta.sources.map((s) => <option key={s}>{s}</option>)}</select></label>
      <label>Chromosome<select aria-label="Merge chromosome" value={chrom} onChange={(e) => ig.setSeg({ chrom: e.target.value })}>{ex.chromosomes.map((c) => <option key={c.name}>{c.name}</option>)}</select></label>
      <label className="chk-inline"><input type="checkbox" aria-label="Allow merging across a genomic gap" checked={ig.allowGapCrossing} onChange={(e) => ig.setAllowGapCrossing(e.target.checked)} />allow merging across a genomic gap (off by default)</label>
      <button className="btn primary" disabled={busy || !runnable} onClick={run}>{busy ? 'Merging…' : 'Run merge'}</button>
    </div>
    {!runnable && <p className="muted small">Enable at least one rule, and complete every enabled rule (each needs at least one condition with valid values), to run a merge.</p>}
    {error && <p role="alert">{error}</p>}
    {result && step && <>
      <p className={stale ? 'notice' : 'muted'} role="status">{stale ? 'Settings changed; showing the previous run. ' : ''}{result.source} {result.chrom}: {result.original.length} original segments → {result.final.length} after {result.steps.length - 1} merge(s). Rules: {ruleSummary}{result.settings.allow_gap_crossing ? ' · gap crossing allowed' : ''}.</p>
      <div className="controls-row">
        <label>Step<input aria-label="Merge step" type="range" min={0} max={result.steps.length - 1} value={stepIdx} onChange={(e) => setStepIdx(+e.target.value)} /></label>
        <span className="muted small">step {stepIdx} of {result.steps.length - 1}{step.removed_boundary !== null ? ` — removed boundary at bin ${step.removed_boundary}: rule "${step.rule_name}" merged segment [${step.segment?.[0]}, ${step.segment?.[1]}) ${step.side}${step.score !== null ? `, |score| ${step.score.toFixed(3)}` : ''}` : ' — original segmentation'}</span>
        <button className="btn-sm" onClick={() => setStepIdx(0)} disabled={stepIdx === 0}>Reset to original</button>
        <button className="btn-sm" onClick={() => setStepIdx((i) => Math.max(0, i - 1))} disabled={stepIdx === 0}>Undo one step</button>
        <button className="btn-sm" onClick={() => setStepIdx(result.steps.length - 1)} disabled={stepIdx === result.steps.length - 1}>Jump to final</button>
        <button className="btn-sm" onClick={() => download(JSON.stringify(result, null, 2), `hatch-merge-${result.chrom}-${result.source}.json`, 'application/json')}>Download JSON</button>
        <button className="btn-sm" onClick={() => {
          const rows = [['step', 'source', 'chromosome', 'merged_start', 'merged_end', 'removed_boundary', 'segment_start', 'segment_end', 'rule_id', 'rule_name', 'side', 'score', 'phi'],
            ...result.steps.filter((s) => s.removed_boundary !== null).map((s) => [s.step, result.source, result.chrom, s.merged_interval?.[0], s.merged_interval?.[1], s.removed_boundary, s.segment?.[0], s.segment?.[1], s.rule_id, JSON.stringify(s.rule_name ?? ''), s.side, s.score ?? '', s.phi ?? ''])]
          download(rows.map((r) => r.join(',')).join('\n') + '\n', `hatch-merge-history-${result.chrom}-${result.source}.csv`, 'text/csv')
        }}>Download merge-history CSV</button>
      </div>
      <DeferredPlot data={[segmentTrace(result.original, result.start_bp, result.end_bp, '#94a3b8', 'Original')]}
        layout={{ ...layout, title: { text: 'Original segmentation', font: { size: 12 } }, xaxis: { title: `${result.chrom} position (Mb)` }, yaxis: { visible: false }, shapes: shapesForSegments(result.original, result.start_bp, result.end_bp) }} />
      <DeferredPlot data={[segmentTrace(step.segments, result.start_bp, result.end_bp, '#96491d', `Step ${stepIdx}`)]}
        layout={{ ...layout, title: { text: `Segmentation at step ${stepIdx} (${step.segments.length} segments)`, font: { size: 12 } }, xaxis: { title: `${result.chrom} position (Mb)` }, yaxis: { visible: false }, shapes: [...shapesForSegments(step.segments, result.start_bp, result.end_bp), ...(step.merged_interval ? hatchShapes(result.start_bp[step.merged_interval[0]] / 1e6, result.end_bp[step.merged_interval[1] - 1] / 1e6, 'y', '/', '#111827') : [])] }} />
      <div className="table-scroll"><table className="tbl"><thead><tr><th>step</th><th>source</th><th>chromosome</th><th>merged interval (bins)</th><th>removed boundary</th><th>segment</th><th>rule</th><th>side</th><th>score</th><th>φ</th></tr></thead><tbody>
        {result.steps.filter((s) => s.removed_boundary !== null).map((s) => <tr key={s.step} className={s.step === stepIdx ? 'sel' : 'clickable'} onClick={() => setStepIdx(s.step)}>
          <td>{s.step}</td><td>{result.source}</td><td>{result.chrom}</td><td>[{s.merged_interval?.[0]}, {s.merged_interval?.[1]})</td><td>{s.removed_boundary}</td>
          <td>[{s.segment?.[0]}, {s.segment?.[1]})</td><td>{s.rule_name}</td><td>{s.side}</td><td>{s.score !== null ? s.score.toFixed(3) : '–'}</td><td>{s.phi != null ? s.phi.toFixed(3) : '–'}</td>
        </tr>)}
      </tbody></table></div>
      <p className="small muted">The rule hatches on the Chromosome tracks are computed from the original accepted-breakpoint segmentation, not from this panel's selected step. This panel shows the merged segmentation's boundaries (the two plots and the history table); it does not display recomputed per-segment scores for the merged partition. CN calls are not recomputed on merged boundaries.</p>
      <details><summary>Merge run provenance and settings</summary><pre style={{ whiteSpace: 'pre-wrap', overflowWrap: 'anywhere' }}>{JSON.stringify({ settings: result.settings, provenance: result.provenance }, null, 2)}</pre></details>
    </>}
  </section>
}
```

- [ ] **Step 6: Update `web/src/components/IntegrationSegmentationView.tsx`**

1. Replace line 13 (`import { classifyAmbiguous, classifyMergeProposal, classifyTransition, scoreForMetric, type Estimator, type HatchMetric } from '../hatchMerge'`) with:
```tsx
import { activeRules, conditionDetail, evaluateForSources, ruleFlagCounts } from '../mergeRules'
```
2. Replace lines 35-39 (from `const hatchLayers = ig.hatchLayers` through the closing `}, [...])` of that `useEffect`) with:
```tsx
  const rulesOn = activeRules(ig.mergeRules).length > 0
  useEffect(() => {
    if (!meta || !rulesOn) return
    for (const s of seg.sources) ig.ensureHatchScores(s, seg.chrom)
  }, [rulesOn, seg.chrom, seg.sources, meta])
  const ruleEvals = useMemo(() => evaluateForSources(ig.hatchScores, seg.sources, seg.chrom, ig.mergeRules, ig.allowGapCrossing), [ig.hatchScores, seg.sources, seg.chrom, ig.mergeRules, ig.allowGapCrossing])
  const ruleCounts = useMemo(() => ruleFlagCounts(ruleEvals), [ruleEvals])
```
3. Replace the block from the comment `// new exploratory hatch layers (genome-wide; independent of the archived-only consistency/proposal layers above)` through the closing `}` of its `for (const row of hs.rows) { ... }` loop (current lines 144-183) with:
```tsx
      // merge-rule hatches (computed client-side from hatch-scores; independent of the archived layers above)
      const hs = ig.hatchScores[`${src}|${seg.chrom}`]
      const evals = ruleEvals[src]
      if (hs && evals && rulesOn) hs.rows.forEach((row, ri) => {
        const ev = evals[ri]
        const x0 = mb[row.start], x1 = row.end < chromMeta.n ? mb[row.end] : bins[2 * (chromMeta.offset + chromMeta.n - 1) + 1] / MB
        for (const id of ev.matched) { const rule = ig.mergeRules.find((r) => r.id === id)!; shapes.push(...hatchShapes(x0, x1, yref, rule.hatch.pattern, rule.hatch.color)) }
        if (row.gap_left) shapes.push({ type: 'line', xref: 'x', yref: `${yref} domain`, x0, x1: x0, y0: 0, y1: 1, line: { color: '#b91c1c', width: 2, dash: 'dot' }, opacity: 0.7 })
        if (row.gap_right) shapes.push({ type: 'line', xref: 'x', yref: `${yref} domain`, x0: x1, x1, y0: 0, y1: 1, line: { color: '#b91c1c', width: 2, dash: 'dot' }, opacity: 0.7 })
        if (!ev.matched.length) return
        const active = new Set(activeRules(ig.mergeRules).map((r) => r.id))
        const ruleLines = ig.mergeRules.map((r, k) => ({ r, n: k + 1 })).filter(({ r }) => active.has(r.id))
          .map(({ r, n }) => `${n}. ${r.name}${ev.matched.includes(r.id) ? ' — MATCH' : ''}<br>${r.conditions.map((c) => `&nbsp;&nbsp;${conditionDetail(c, row)}`).join('<br>')}`)
        const claimed = ig.mergeRules.find((r) => r.id === ev.claimedBy)
        const outcome = ev.blocked === 'gap' ? `would merge ${ev.side}, blocked by gap` : ev.blocked === 'nonfinite' ? 'not mergeable: the merge-into score is missing' : `merges ${ev.side} (|score| ${fmt(ev.score, 3)})`
        data.push({ type: 'scatter', mode: 'markers', x: [(x0 + x1) / 2], y: [seg.basis === 'raw' ? 0 : yr[0] + 0.2], marker: { size: 7, color: claimed?.hatch.color ?? '#155e75', line: { color: '#111827', width: 1 }, symbol: 'diamond' }, xaxis: 'x', yaxis: yref, showlegend: false, name: 'merge rule',
          hovertemplate: `<b>${src} ${seg.chrom} · merge rules</b><br>bins [${row.start}, ${row.end}) · ${row.end - row.start} bins<br>claimed by: ${claimed?.name ?? '–'} → ${outcome}<br>${ruleLines.join('<br>')}<extra></extra>` })
      })
```
4. In the `trackFig` dependency array, replace `hatchLayers, ig.hatchScores]` with `ig.mergeRules, ig.hatchScores, ruleEvals, rulesOn]`.
5. Immediately before the first `<div className="controls-row">` inside the `#cluster-tracks` panel (the one starting with the `chromosome` select), insert:
```tsx
        <h4 className="vis-heading">1. Visualization <span className="muted small">— display only; nothing here affects merging (merge rules are in section 2 below the tracks)</span></h4>
```
6. Replace the whole third controls row (the `<div className="controls-row">` that starts with `new exploratory hatch layers (genome-wide, computed here — separate from the archived layers above):` and ends after the `{ig.hatchScoresError && ...}` span) with:
```tsx
        {ig.hatchScoresError && <div className="notice small">{ig.hatchScoresError}</div>}
```
7. In the `seg-legend` div, replace the four lines starting `{hatchLayers.transition.on && <span>new: hatch /`, `{hatchLayers.proposal.on && ...`, `{hatchLayers.ambiguous.on && ...`, and `{(hatchLayers.transition.on || hatchLayers.proposal.on || hatchLayers.ambiguous.on) && ...` with:
```tsx
          {activeRules(ig.mergeRules).map((r) => <span key={r.id}><b className="rule-swatch" style={{ color: r.hatch.color }}>{r.hatch.pattern}</b>rule {ig.mergeRules.indexOf(r) + 1} “{r.name}” · flags {ruleCounts[r.id] || 0}</span>)}
          {rulesOn && <span><i className="ln" style={{ borderColor: '#b91c1c' }} />dotted red = a genomic gap sits at that flank boundary</span>}
```

- [ ] **Step 7: Append rule styles to `web/src/styles.css`**

```css
.rule-card { border: 1px solid #e2e8f0; border-radius: 8px; padding: 8px 10px; margin: 8px 0; background: #fafafa; }
.rule-head, .rule-cond { display: flex; flex-wrap: wrap; gap: 6px; align-items: center; margin: 4px 0; }
.rule-cond { margin-left: 16px; }
.rule-swatch { display: inline-block; min-width: 16px; margin-right: 4px; text-align: center; font-weight: 700; }
.vis-heading { margin: 6px 0 4px; font-size: 13px; }
```

- [ ] **Step 8: Check for leftovers and build**

Run: `grep -rn "hatchLayers\|setHatchLayers\|HatchLayerSettings\|TransitionLayerConfig\|ProposalLayerConfig\|AmbiguousLayerConfig\|small_max_bins\|\.fired" web/src`
Expected: only the `migrateHatchLayers` function and its call in `integrationStore.ts`.

Run: `cd web && npm run build`
Expected: `tsc --noEmit` passes and `vite build` prints `✓ built in …`.

- [ ] **Step 9: Commit**

```bash
git add web/src/integration.ts web/src/hatchMerge.ts web/src/integrationStore.ts web/src/components/MergeRulesEditor.tsx web/src/components/HatchMergeView.tsx web/src/components/IntegrationSegmentationView.tsx web/src/styles.css
git commit -m "feat: merge-rules editor, per-rule hatches, and visualization/merge split

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Browser walkthrough and validation handoff

**Files:**
- Rewrite: `tests/browser_hatch_merge_walkthrough.cjs`
- Modify: `validation/hatch-merge/HANDOFF.md` (append a section)
- Generated (commit): `validation/hatch-merge/merge_rules_check_results.json`, `validation/hatch-merge/merge-rules.png`

**Interfaces:**
- Consumes: the accessible names listed in Task 4's Produces block; session JSON key `integration.mergeRules` (found by recursive search, so nesting changes don't break the test); the existing SessionBar button `download JSON` and file input `input[type="file"][accept="application/json"]`.
- Produces: a reproducible end-to-end check of the UI; closes the UI-verification gap flagged on commit `3c86329`.

- [ ] **Step 1: Replace `tests/browser_hatch_merge_walkthrough.cjs` with:**

```js
// Merge-rules browser walkthrough (2026-10-02): visualization/merge split, rule editor,
// live per-rule hatches, reorder, run, stale notice, session restore, legacy migration.
// Requires a running server (python3 run.py --port 8766) and Playwright + Chromium, e.g.:
//   PLAYWRIGHT_MODULE=/path/to/node_modules/playwright CHROMIUM_PATH=/path/to/chrome \
//   node tests/browser_hatch_merge_walkthrough.cjs
// Writes validation/hatch-merge/merge_rules_check_results.json and validation/hatch-merge/merge-rules.png.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const fs=require('fs'),path=require('path'),os=require('os'),assert=require('assert/strict');
function findOwner(obj,key){if(!obj||typeof obj!=='object')return null;if(Object.prototype.hasOwnProperty.call(obj,key))return obj;for(const v of Object.values(obj)){const f=findOwner(v,key);if(f)return f}return null}
(async()=>{
 const browser=await chromium.launch({headless:process.env.HEADED!=='1',executablePath:process.env.CHROMIUM_PATH});
 try {
 const page=await browser.newPage({viewport:{width:1400,height:1100}}),errors=[];
 page.on('pageerror',e=>errors.push(String(e)));
 page.on('console',msg=>{if(msg.type()==='error')errors.push('console: '+msg.text())});
 await page.goto('http://127.0.0.1:8766');await page.getByText('300 matching cells',{exact:true}).waitFor();
 await page.getByRole('button',{name:'Integration · segmentation & CN',exact:true}).click();
 await page.getByText(/Showing 760 of 3460 rows/).waitFor({timeout:60000});

 const tracks=page.locator('#cluster-tracks');await tracks.scrollIntoViewIfNeeded();
 const plot=tracks.locator('.js-plotly-plot');await plot.waitFor();
 await page.waitForFunction(()=>document.querySelector('#cluster-tracks .js-plotly-plot')?._fullLayout);
 const shapeCount=()=>plot.evaluate(el=>el.layout.shapes.length);
 const visHeading=await tracks.getByText('1. Visualization').count();
 const shapesBaseline=await shapeCount();

 const panel=page.locator('.hatch-merge');await panel.scrollIntoViewIfNeeded();
 await page.getByRole('heading',{name:/2\. Merge rules/}).waitFor();
 const runBtn=panel.getByRole('button',{name:'Run merge'});
 const runDisabledWithNoRules=await runBtn.isDisabled();

 // rule 1: weak-flank preset; drive the threshold high, then near zero
 await panel.getByLabel('Add rule').selectOption('weak_flank');
 const status1=panel.getByLabel('Rule 1 status');
 const flags=async(label)=>{const m=(await panel.getByLabel(label).innerText()).match(/flags (\d+)/);return m?+m[1]:null};
 const thr1=panel.getByLabel('Rule 1 condition 1 threshold');
 await thr1.fill('1000');
 await page.waitForFunction(()=>/flags [1-9]\d*/.test(document.querySelector('[aria-label="Rule 1 status"]')?.textContent||''),{timeout:30000});
 const flagsHigh=await flags('Rule 1 status');
 await page.waitForTimeout(600);
 const shapesWithRule=await shapeCount();
 await thr1.fill('0.0001');await page.waitForTimeout(400);
 const flagsLow=await flags('Rule 1 status');
 await thr1.fill('');await page.waitForTimeout(300);
 const emptyThresholdProblem=(await status1.innerText()).includes('threshold');
 const runDisabledWhenIncomplete=await runBtn.isDisabled();
 await thr1.fill('3.3');await page.waitForTimeout(400);

 // rule 2: sign-mismatch preset; legend lists both; move it to the top
 await panel.getByLabel('Add rule').selectOption('sign_mismatch');await page.waitForTimeout(400);
 const legendText=await tracks.locator('.seg-legend').last().innerText();
 const legendHasBothRules=legendText.includes('weak SRD/√φ flank')&&legendText.includes('flank signs disagree');
 await panel.getByLabel('Move rule 2 up').click();
 const rule1NameAfterMove=await panel.getByLabel('Rule 1 name').inputValue();

 // run
 await panel.getByLabel('Allow merging across a genomic gap').check();
 await runBtn.click();
 await page.waitForFunction(()=>document.querySelectorAll('.hatch-merge .js-plotly-plot').length>=2,{timeout:60000});
 await page.waitForTimeout(500);
 const summaryText=await panel.locator('p[role="status"]').first().innerText();
 const headers=await panel.locator('table.tbl thead th').allInnerTexts();
 const historyRows=await panel.locator('table.tbl tbody tr').count();
 const firstRuleCell=historyRows?await panel.locator('table.tbl tbody tr').first().locator('td').nth(6).innerText():null;
 const out=path.resolve(__dirname,'../validation/hatch-merge');fs.mkdirSync(out,{recursive:true});
 await panel.screenshot({path:path.join(out,'merge-rules.png')});

 // stale notice after editing a rule (the weak-flank rule is now rule 2)
 await panel.getByLabel('Rule 2 condition 1 threshold').fill('2.0');await page.waitForTimeout(300);
 const staleText=await panel.locator('p[role="status"]').first().innerText();

 // session export -> change -> restore
 const enabled1=panel.getByLabel('Rule 1 enabled');
 const [download]=await Promise.all([page.waitForEvent('download'),page.getByRole('button',{name:'download JSON',exact:true}).click()]);
 const dir=fs.mkdtempSync(path.join(os.tmpdir(),'merge-rules-session-'));const sessionPath=path.join(dir,'session.json');
 await download.saveAs(sessionPath);
 await enabled1.uncheck();const enabledBeforeRestore=await enabled1.isChecked();
 const fileInput=page.locator('input[type="file"][accept="application/json"]');
 await fileInput.setInputFiles(sessionPath);await page.waitForTimeout(600);
 const enabledAfterRestore=await panel.getByLabel('Rule 1 enabled').isChecked();
 const rulesAfterRestore=await panel.locator('.rule-card').count();

 // legacy session: only hatchLayers present -> three migrated preset rules
 const session=JSON.parse(fs.readFileSync(sessionPath,'utf8'));
 const owner=findOwner(session,'mergeRules');delete owner.mergeRules;delete owner.allowGapCrossing;
 owner.hatchLayers={transition:{on:true,metric:'srd_phi',estimator:'mean'},proposal:{on:false,metric:'srd_phi',estimator:'mean',threshold:2.5},ambiguous:{on:true,estimator:'median'}};
 const legacyPath=path.join(dir,'legacy.json');fs.writeFileSync(legacyPath,JSON.stringify(session));
 await fileInput.setInputFiles(legacyPath);await page.waitForTimeout(600);
 const nCards=await panel.locator('.rule-card').count();const legacyNames=[],legacyEnabled=[];
 for(let i=1;i<=nCards;i++){legacyNames.push(await panel.getByLabel(`Rule ${i} name`).inputValue());legacyEnabled.push(await panel.getByLabel(`Rule ${i} enabled`).isChecked())}
 fs.rmSync(dir,{recursive:true,force:true});

 const result={generated_at:new Date().toISOString(),visHeading,shapesBaseline,shapesWithRule,runDisabledWithNoRules,flagsHigh,flagsLow,emptyThresholdProblem,runDisabledWhenIncomplete,legendHasBothRules,rule1NameAfterMove,summaryText,headers,historyRows,firstRuleCell,staleText,enabledBeforeRestore,enabledAfterRestore,rulesAfterRestore,legacyNames,legacyEnabled,errors};
 fs.writeFileSync(path.join(out,'merge_rules_check_results.json'),JSON.stringify(result,null,2));
 console.log(JSON.stringify(result,null,2));
 assert.ok(visHeading>0,'visualization section heading missing');
 assert.equal(runDisabledWithNoRules,true,'Run merge should be disabled with no rules');
 assert.ok(flagsHigh>flagsLow,'a high threshold should flag more segments than a near-zero one');
 assert.ok(shapesWithRule>shapesBaseline,'an enabled rule should add hatch shapes to the tracks');
 assert.equal(emptyThresholdProblem,true,'a cleared threshold should show a problem on the rule card');
 assert.equal(runDisabledWhenIncomplete,true,'an incomplete enabled rule should block the run');
 assert.equal(legendHasBothRules,true,'legend should list both rules');
 assert.equal(rule1NameAfterMove,'flank signs disagree','moving rule 2 up should put it first');
 assert.ok(headers.includes('rule'),'history table needs a rule column');
 if(historyRows>0)assert.ok(firstRuleCell&&firstRuleCell.trim().length>0,'history rows should name the rule');
 assert.ok(staleText.includes('Settings changed'),'editing a rule after a run should mark the result stale');
 assert.equal(enabledBeforeRestore,false);
 assert.equal(enabledAfterRestore,true,'session restore should bring back the saved rule state');
 assert.equal(rulesAfterRestore,2);
 assert.deepEqual(legacyNames,['transition zone (from saved session)','merge proposal (from saved session)','ambiguous flank preference (from saved session)'],'legacy hatchLayers should migrate to three preset rules');
 assert.deepEqual(legacyEnabled,[true,false,true]);
 assert.deepEqual(errors,[],'console/page errors present');
 console.log('PASS merge rules: split sections, editor, live hatches, reorder, run, stale, session restore, legacy migration');
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1)});
```

- [ ] **Step 2: Start the server with the freshly built frontend**

Run (background): `python3 run.py --port 8766`
Then: `curl -s http://127.0.0.1:8766/api/health`
Expected: a JSON health response. If a server is already running on 8766 from before Task 4, stop it first so the new `web/dist` and backend are served.

- [ ] **Step 3: Run the walkthrough**

These local paths were verified to exist on 2026-10-02:
```bash
PLAYWRIGHT_MODULE=/home/mami/Desktop/EGA/job_search/career-ops/node_modules/playwright \
CHROMIUM_PATH=/home/mami/.cache/ms-playwright/chromium-1234/chrome-linux64/chrome \
node tests/browser_hatch_merge_walkthrough.cjs
```
Expected: last line `PASS merge rules: ...`. If Chromium fails to launch because of a version mismatch with that Playwright module, retry with `CHROMIUM_PATH=/home/mami/.cache/ms-playwright/chromium-1217/chrome-linux64/chrome`. If an assertion fails, fix the code, not the assertion (use superpowers:systematic-debugging).

- [ ] **Step 4: Append the outcome to `validation/hatch-merge/HANDOFF.md`**

Append a section titled `## Merge rules redesign (2026-10-02)` stating, using the actual values from `validation/hatch-merge/merge_rules_check_results.json`: the command that was run; `flagsHigh` and `flagsLow` for the weak-flank rule; `historyRows` and the summary text of the merge run; that session restore and legacy migration passed; that `errors` was empty; and the backend test count from the Task 6 suite run. Close with one sentence noting this supersedes the union-of-layers model from commit `3c86329`, whose UI-level check this run also covers.

- [ ] **Step 5: Commit**

```bash
git add tests/browser_hatch_merge_walkthrough.cjs validation/hatch-merge/HANDOFF.md validation/hatch-merge/merge_rules_check_results.json validation/hatch-merge/merge-rules.png
git commit -m "test: browser walkthrough for merge rules; record validation results

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Documentation and final verification

**Files:**
- Modify: `README.md`, `docs/USAGE.md`, `docs/SCIENTIFIC_WORKFLOW.md`, `docs/ARCHITECTURE.md`, `tests/README.md`, `AGENTS.md`, `CONTRIBUTING.md`, `server/README.md`, `web/README.md`

**Interfaces:**
- Consumes: the final behavior from Tasks 1-5.
- Produces: docs that describe merge rules, not hatch layers.

- [ ] **Step 1: README API paragraph**

In `README.md`, replace the paragraph that begins `` `POST /api/integration/hatch-merge` `{chrom, source, transition?, proposal?, ambiguous?, small_max_bins?, allow_gap_crossing?}` `` with:

```markdown
`POST /api/integration/hatch-merge` `{chrom, source, rules, allow_gap_crossing?}` runs the exploratory greedy merge driven by user-defined rules. Each rule is `{id, name, enabled, conditions, mergeInto, hatch}`; its conditions (`weak_flank`, `sign_mismatch`, `ambiguous`, `max_size`) must all hold for an internal segment, and a segment is flagged when any enabled rule matches. The topmost matching rule claims the segment and its `mergeInto` metric picks the neighbour to merge into (smaller |score|, ties left). One merge happens per step, ordered by rule, then score, then position, with every score and the pooled φ recomputed after each merge; gap-blocked segments are skipped, not redirected. Disabled rules are ignored; invalid enabled rules return 422. The response returns the full step history (original → final), naming the rule behind each merge. See `server/merge_rules.py` and `tests/test_merge_rules.py`.
```

Also, in the README's "What the explorer can do" Integration · segmentation & CN bullet, replace `the exploratory hatch-layer/adjacent-segment merge tool` with `user-defined merge rules with per-rule hatching and an exploratory adjacent-segment merge`.

- [ ] **Step 2: USAGE section**

In `docs/USAGE.md`, replace the whole `### Exploratory hatching and adjacent-segment merging` section (heading plus its two paragraphs) with:

```markdown
### Visualization and merge rules

The Chromosome tracks panel starts with **1. Visualization**: colouring, flank choice, scale, track level, smoothing, size and reference settings, and the archived label layers. These settings only change what is drawn; none of them affect merging.

**2. Merge rules**, below the tracks, holds an ordered list of rules you build yourself. Add a blank rule or a preset (weak flank, flank signs disagree, ambiguous flank preference, small segment). Each rule has:

- **Conditions, all of which must hold** for an internal segment: *weaker flank below threshold* (min(|L|, |R|) < threshold, for SRD, SRD/√φ or Δlog2FC with an estimator), *flank signs disagree*, *SRD and Δlog2FC disagree on the weaker flank*, and *segment size at most N bins*.
- **Merge into**: the metric that decides which neighbour a flagged segment joins (the flank with the smaller |score|; ties go left).
- **A hatch pattern and colour.** Each enabled rule hatches exactly the segments it flags, live, and its card shows how many segments it flags on the displayed rows. Click a hatched segment's diamond marker for each condition's values and the merge outcome.

A segment is flagged when any enabled rule matches. **Run merge** then merges one flagged segment per step: rules higher in the list go first, then the smallest score, then the leftmost segment, recomputing every score, including the pooled dispersion φ, after each merge. Genomic gaps are never crossed unless you allow it; a segment whose chosen neighbour sits across a gap is skipped, not redirected. Use the step slider (or Reset/Undo) to inspect intermediate segmentations and download the run (JSON) or the merge history (CSV), which names the rule behind each merge. Rules are saved in sessions; sessions saved before this change are converted into the equivalent preset rules. These rules are exploratory diagnostics, not the original bootstrap/chain-guard merge prototype, and CN calls are not recomputed on merged boundaries.
```

- [ ] **Step 3: SCIENTIFIC_WORKFLOW worked example 2**

In `docs/SCIENTIFIC_WORKFLOW.md`, worked example 2, replace the **Diagnostic** bullet with:

```markdown
- **Diagnostic:** user-defined merge rules on the Chromosome tracks. Each rule
  combines conditions on a segment's flank scores (e.g. "weaker SRD/√φ flank
  below 3 AND segment ≤ 20 bins", or "Δlog2FC flank signs disagree") and
  hatches the segments it flags, followed by the exploratory adjacent-segment
  merge, which acts on exactly those flagged segments.
```

and replace the **Recompute** bullet with:

```markdown
- **Recompute:** one segment merges per step, into its more similar
  neighbour under the claiming rule's metric; pooled dispersion φ is
  recomputed, unconditionally, after every merge, since it is shared across
  the whole chromosome's segmentation.
```

- [ ] **Step 4: ARCHITECTURE**

In `docs/ARCHITECTURE.md`:
- Replace the backend table row for `server/hatch_merge.py` / `server/hatch_merge_api.py` with two rows:
```markdown
| `server/hatch_merge.py` | Flank scoring for the exploratory tools: SRD, SRD/√φ, Δlog2FC, pooled φ, gap flags, and `score_rows()` / `hatch_scores()`. Independent of `flank_view.py` and the archived merge prototype. |
| `server/merge_rules.py` / `server/hatch_merge_api.py` | User-defined merge rules: validation, per-segment evaluation (OR of ANDs), and the one-merge-per-step loop with φ recomputed every step; the `hatch-scores` and `hatch-merge` routes. See [SCIENTIFIC_WORKFLOW.md](SCIENTIFIC_WORKFLOW.md)'s worked example 2. |
```
- Replace the frontend row for `hatchMerge.ts` with two rows:
```markdown
| `hatchMerge.ts` | Score/result types and the API client for `hatch-scores` / `hatch-merge`. |
| `mergeRules.ts` | Rule types, presets, legacy-session migration, and a 1:1 TypeScript port of `server/merge_rules.py`'s evaluator, so hatches update live while rules are edited. `tests/fixtures/merge_rules_parity.json` is checked against both implementations. |
```
- In the `components/` row, add `MergeRulesEditor` to the list of components.
- In "Data flow: one representative request", replace `a Zustand store action reads the relevant settings` with `HatchMergeView reads the merge rules from the integration store`.

- [ ] **Step 5: Test commands and directory READMEs**

- In `tests/README.md`: add `test_merge_rules.py` to the pytest list; add `merge_rules_helpers.ts` to the TypeScript-helpers list with its run command (`web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs && node tests/_bundle/merge_rules_helpers.mjs`); add a line that `fixtures/merge_rules_parity.json` is shared by the Python and TypeScript evaluator tests; and add `tests/test_merge_rules.py` to the "Minimal sanity check" pytest command.
- In `AGENTS.md`, `CONTRIBUTING.md`, and `README.md` (Quick start → Tests), add `tests/test_merge_rules.py` to the pytest command.
- In `server/README.md`: replace the `hatch_merge.py` / `hatch_merge_api.py` bullet with `hatch_merge.py` (exploratory flank scoring) and `merge_rules.py` / `hatch_merge_api.py` (user-defined merge rules and their routes).
- In `web/README.md`: add `MergeRulesEditor` to the components list, and add one line under Stores: `mergeRules.ts` holds the rule types and the client-side evaluator (a port of `server/merge_rules.py`).

- [ ] **Step 6: Check for stale wording and links**

Run: `grep -rn "hatch layer\|hatchLayers\|transition zone\|veto\|small_max_bins" README.md docs/USAGE.md docs/ARCHITECTURE.md docs/SCIENTIFIC_WORKFLOW.md docs/EXTENDING.md server/README.md web/README.md tests/README.md`
Expected: no matches except (a) docs/EXTENDING.md's worked example and history sentence about the feature's past designs and (b) USAGE's mention of "preset" names. Reword anything else that still describes the old layers as current.

Then re-run the repository's Markdown link check (from `docs/superpowers/plans/2026-09-30-portfolio-repository-polish.md` Task 14 Step 1, with `docs/superpowers/plans` excluded).
Expected: `All links in N Markdown files resolve`.

- [ ] **Step 7: Final verification**

Run:
```bash
python3 -m pytest -q tests/test_merge_rules.py tests/test_hatch_merge.py tests/test_experiments.py tests/test_segment_estimators.py tests/test_integration_review.py
web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs && node tests/_bundle/merge_rules_helpers.mjs
web/node_modules/.bin/esbuild tests/hatch_merge_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/hatch_merge_helpers.mjs && node tests/_bundle/hatch_merge_helpers.mjs
cd web && npm run build && cd ..
git diff --stat 3c86329..HEAD -- server/science.py server/flank_view.py server/pyepi_adapter.py reference/
```
Expected: pytest all pass; both TS helpers print `all … checks passed`; build succeeds; the final `git diff --stat` prints nothing (scoring and reference code untouched). Record the pytest count in the HANDOFF.md section from Task 5 if it was left pending.

- [ ] **Step 8: Commit**

```bash
git add README.md docs/USAGE.md docs/SCIENTIFIC_WORKFLOW.md docs/ARCHITECTURE.md tests/README.md AGENTS.md CONTRIBUTING.md server/README.md web/README.md validation/hatch-merge/HANDOFF.md
git commit -m "docs: describe user-defined merge rules and the visualization/merge split

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
