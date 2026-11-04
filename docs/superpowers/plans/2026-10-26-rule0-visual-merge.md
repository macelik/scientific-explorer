# Rule 0 (merge what the tracks hatch) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make section 1 of the segmentation tab define "Rule 0": a checklist of flank criteria with thresholds that hatches segments on the tracks and is always the first rule of the merge. Also switch the tracks to live scores, redefine the consistency hatch as SRD/√φ sign disagreement, add p-value reasons to the candidate audit, and drop the SRD-screen proposal layer.

**Architecture:** Rule 0 is a new rule type (`type: "flank_threshold"`, id `rule0`) in the existing rule engine, implemented in Python (`server/merge_rules.py`, authority for the merge) and ported 1:1 to TypeScript (`web/src/mergeRules.ts`, live hatching), with the shared parity fixture extended. The tracks draw small segments from the live `hatch-scores` rows instead of the archived tables, using new pure helpers in `web/src/trackMetrics.ts`. Every merge request sends Rule 0 first, built from the section-1 settings.

**Tech Stack:** Python 3.10, FastAPI, NumPy, pytest; React 18 + TypeScript + Vite + Zustand + Plotly.js; esbuild + Node for TS unit tests; Playwright (CommonJS) for the browser walkthrough.

**Spec:** `docs/superpowers/specs/2026-10-26-rule0-visual-merge-design.md` — read it alongside this plan. It builds on `docs/superpowers/specs/2026-10-02-merge-rules-design.md`; where they disagree, the 2026-10-26 spec wins.

## Global Constraints

- Work from the repository root `scientific-explorer/`, branch `merge-rules`. Python: `source .venv-scientific-explorer/bin/activate`.
- Do not change score computation (`server/hatch_merge.py`), the `hatch-scores` response, or the archived diagnostics panel's own logic (figs 4–8, merge-prototype summary).
- No new dependencies.
- Never stage with `git add -A`, `git add .` or a directory path. Stage files by name. Two untracked files in `docs/superpowers/plans/` (`2026-09-21-typesafe-jev-judgment-layers.md`, `2026-09-30-portfolio-repository-polish.md`) must stay untracked and untouched. Do not commit `web/dist` or `tests/_bundle/*.mjs`.
- Never stop or restart the process on port 8766 (it belongs to the user). Tasks that need a server start their own on port 8767 and stop only that one.
- Rule 0 object, exactly: `{"id": "rule0", "type": "flank_threshold", "name": <string>, "criteria": <non-empty distinct subset of ["delta_median","delta_mean","delta_iqr_mean","srd_phi"]>, "combine": "all"|"any", "tau": >0 finite, "srdRef": >0 finite, "maxBins": int ≥ 1, "mergeSignMismatch": bool, "mergeFlankConflict": bool, "hatch": {"pattern": <one of '/','\\','x','-','|','.'>, "color": <string>}}`. It must be `rules[0]`; no other rule may have id `rule0` or a `type` other than `"conditions"` (default).
- Rule 0 semantics (spec §3):
  - Only internal segments with `end - start ≤ maxBins` are evaluated.
  - **Criterion flank:** each criterion picks the side with the smaller |value|. A missing side loses, and ties go left. A criterion passes when `|value| < threshold`, using τ for Δlog2FC criteria and srdRef for `srd_phi`. Its fraction is `|value| / threshold`.
  - **Combining:** `all` requires every checked criterion to pass; `any` requires at least one.
  - **Signs and flanks:** a sign mismatch means `srd_phi_left * srd_phi_right < 0`. A flank conflict means the passing criteria pick different sides.
  - **Direction:** when there is no conflict, use the shared side. When there is a conflict, use the passing criterion with the smallest fraction (ties go left).
  - **Score:** `all` uses the max fraction over the checked criteria; `any` uses the min fraction over the passing ones.
  - **Held back:** a segment is held back on a sign mismatch with `mergeSignMismatch` false, and/or on a flank conflict with `mergeFlankConflict` false.
  - **Gap blocking:** a would-merge segment whose direction side has a gap, while gap crossing is off, is blocked (never redirected).
- Merge order is unchanged: (rule index, score, segment start), Rule 0 is index 0, and φ and every criterion are recomputed after each merge.
- Colours (spec §4): Rule 0 would-merge `×` in `#b91c1c`; sign mark amber `|` in `#d97706`; flank-conflict mark violet `―` (pattern `'-'`) in `#7c3aed`.
- UI copy calls Rule 0 and rules exploratory diagnostics. Never "validated", "proves", or "ground truth". The candidate-audit cutoff is labelled "inferred".
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. **τ or SRD/√φ ref cleared or 0, or no criterion checked** → the tracks still colour segments; Rule 0 draws no hatches, its card and section 1 show a problem, and "Run merge" is disabled. Tests: Task 2 (`rule0Problems`, `evaluateRules` ignoring an invalid Rule 0, `canRunMerge`), Task 5 (browser).
2. **φ undefined (constant levels)** → `srd_phi` is missing, so the SRD/√φ criterion fails, no sign mismatch is reported, and Δlog2FC criteria still work. Test: Task 1 `test_run_rule0_merges_by_fraction_with_constant_levels` (φ = 0 there) and `test_rule0_missing_values`.
3. **Every checked criterion missing on both sides** → not passed, no crash, outcome `none`. Test: Task 1 `test_rule0_missing_values`.
4. **A non-pilot chromosome** → tracks show live colouring and Rule 0 hatches, where they used to be empty. Test: Task 5 (browser, non-pilot chromosome).
5. **An old saved session** (`seg.metric` of `log2fc` or `delta_trim`, `labels.proposal` / `labels.consistency`, `seg.flank` / `seg.effectMetric`) → it migrates without a crash. Tests: Task 2 `migrateSegSettings`, Task 5 legacy-session restore.

---

## File Structure

- Modify `server/merge_rules.py` — Rule 0 validation and `evaluate_rule0`; wire into `validate_rules`, `evaluate_rules`, `run_rule_merge`.
- Modify `tests/test_merge_rules.py` — Rule 0 unit/run/route tests; existing validate/run/route tests pass Rule 0 first.
- Modify `tests/fixtures/merge_rules_parity.json` — inert Rule 0 in existing cases; four Rule 0 cases.
- Modify `web/src/mergeRules.ts` — Rule 0 types, evaluator, problems, description, tooltip lines, case counts; `canRunMerge(rule0, rules)`.
- Create `web/src/trackMetrics.ts` — live colour metrics, colour value, archived-id mapping, session migration, derived candidate-audit reason.
- Modify `tests/merge_rules_helpers.ts` — parity incl. Rule 0; new `canRunMerge` semantics.
- Create `tests/track_helpers.ts` — tests for `trackMetrics.ts`.
- Create `tests/compare_archived_live.py`, `validation/rule0/archived_vs_live.json`, `validation/rule0/README.md`.
- Modify `web/src/hatchMerge.ts`, `web/src/integrationStore.ts`, `web/src/components/IntegrationSegmentationView.tsx`, `web/src/components/HatchMergeView.tsx`, `web/src/components/MergeRulesEditor.tsx`, `web/src/styles.css`.
- Modify `tests/browser_hatch_merge_walkthrough.cjs`; `validation/hatch-merge/HANDOFF.md`; docs (`README.md`, `docs/USAGE.md`, `docs/SCIENTIFIC_WORKFLOW.md`, `docs/ARCHITECTURE.md`, `tests/README.md`, `web/README.md`).

---

### Task 1: Python Rule 0

**Files:**
- Modify: `server/merge_rules.py`
- Test: `tests/test_merge_rules.py`

**Interfaces:**
- Consumes (existing): `hm.flat_partition`, `hm.score_rows`, `HATCH_PATTERNS`, `_number`, `condition_holds`, `validate_rules`, `evaluate_rules`, `run_rule_merge`.
- Produces:
  - constants `RULE0_ID = 'rule0'`, `RULE0_TYPE = 'flank_threshold'`, `RULE0_CRITERIA = ('delta_median', 'delta_mean', 'delta_iqr_mean', 'srd_phi')`, `RULE0_COMBINE = ('all', 'any')`
  - `evaluate_rule0(row: dict, r0: dict, allow_gap_crossing: bool = False) -> dict` with keys `start, end, evaluated, criteria (list of {criterion, side, value, threshold, fraction, passed}), passed, sign_mismatch, flank_conflict, held_back (list of 'signs'/'flanks'), outcome ('none'|'would_merge'|'held_back'|'signs_only'), side, score, blocked (None|'gap')`
  - `validate_rules(rules)` now requires Rule 0 at `rules[0]` and returns `[rule0, *enabled_user_rules]`. A list with only Rule 0 is valid.
  - `evaluate_rules(rows, rules, allow_gap_crossing)`: if `rules[0]` is a Rule 0, each result gains `'rule0': <evaluate_rule0 dict>` and `matched` starts with `'rule0'` when its outcome is `would_merge`. When Rule 0 claims, `side`/`score`/`blocked` come from it. Without a leading Rule 0, `'rule0'` is `None` and behaviour is unchanged.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_merge_rules.py`, directly after the `GOOD = _rule('ok', [WEAK_SRD_PHI])` line:

```python
def _rule0(**kw):
    r = {'id': 'rule0', 'type': 'flank_threshold', 'name': 'Rule 0', 'criteria': ['delta_mean'],
         'combine': 'all', 'tau': 0.5, 'srdRef': 3.3, 'maxBins': 100,
         'mergeSignMismatch': False, 'mergeFlankConflict': False,
         'hatch': {'pattern': 'x', 'color': '#b91c1c'}}
    r.update(kw)
    return r


# Never matches: only 1-bin segments are evaluated, and |SRD/√φ| < 1e-9 essentially never holds.
INERT_RULE0 = _rule0(criteria=['srd_phi'], tau=1e-9, srdRef=1e-9, maxBins=1)


def _with_r0(rules):
    return [INERT_RULE0, *rules]
```

Append to the end of `tests/test_merge_rules.py`:

```python
# ---------------- Rule 0

def _r0_row(start=10, end=15, gap_left=False, gap_right=False, **scores):
    return _row(start=start, end=end, gap_left=gap_left, gap_right=gap_right, **scores)


def test_rule0_criterion_picks_smaller_flank_and_strict_threshold():
    row = _r0_row(delta_mean_left=-0.4, delta_mean_right=0.1)
    ev = mr.evaluate_rule0(row, _rule0(tau=0.5))
    c = ev['criteria'][0]
    assert (c['criterion'], c['side'], c['value'], c['threshold']) == ('delta_mean', 'right', pytest.approx(0.1), 0.5)
    assert c['fraction'] == pytest.approx(0.2) and c['passed'] is True
    assert ev['passed'] is True and ev['outcome'] == 'would_merge' and ev['side'] == 'right'
    assert ev['score'] == pytest.approx(0.2)
    assert mr.evaluate_rule0(row, _rule0(tau=0.1))['criteria'][0]['passed'] is False


def test_rule0_tie_goes_left_and_missing_side_loses():
    assert mr.evaluate_rule0(_r0_row(delta_mean_left=-0.2, delta_mean_right=0.2), _rule0())['side'] == 'left'
    ev = mr.evaluate_rule0(_r0_row(delta_mean_left=None, delta_mean_right=0.3), _rule0())
    assert ev['criteria'][0]['side'] == 'right' and ev['outcome'] == 'would_merge'


def test_rule0_missing_values():
    ev = mr.evaluate_rule0(_r0_row(), _rule0(criteria=['delta_mean', 'srd_phi'], combine='any'))
    assert ev['evaluated'] is True and ev['passed'] is False and ev['outcome'] == 'none'
    assert [c['side'] for c in ev['criteria']] == [None, None]
    assert ev['sign_mismatch'] is False and ev['score'] is None


def test_rule0_all_vs_any_and_scores():
    row = _r0_row(delta_mean_left=0.1, delta_mean_right=0.9, srd_phi_left=2.0, srd_phi_right=5.0)
    both = ['delta_mean', 'srd_phi']
    ev = mr.evaluate_rule0(row, _rule0(criteria=both, combine='all', tau=0.5, srdRef=3.3))
    assert ev['passed'] is True and ev['score'] == pytest.approx(max(0.1 / 0.5, 2.0 / 3.3))
    ev = mr.evaluate_rule0(row, _rule0(criteria=both, combine='any', tau=0.5, srdRef=3.3))
    assert ev['score'] == pytest.approx(min(0.1 / 0.5, 2.0 / 3.3))
    ev = mr.evaluate_rule0(row, _rule0(criteria=both, combine='all', tau=0.5, srdRef=1.0))
    assert ev['passed'] is False and ev['outcome'] == 'none'
    ev = mr.evaluate_rule0(row, _rule0(criteria=both, combine='any', tau=0.5, srdRef=1.0))
    assert ev['passed'] is True and ev['score'] == pytest.approx(0.2)


def test_rule0_sign_mismatch_held_back_unless_switch_on():
    row = _r0_row(delta_mean_left=0.1, delta_mean_right=0.9, srd_phi_left=1.0, srd_phi_right=-4.0)
    ev = mr.evaluate_rule0(row, _rule0())
    assert ev['sign_mismatch'] is True and ev['held_back'] == ['signs'] and ev['outcome'] == 'held_back'
    assert ev['side'] == 'left' and ev['score'] == pytest.approx(0.2)
    ev = mr.evaluate_rule0(row, _rule0(mergeSignMismatch=True))
    assert ev['held_back'] == [] and ev['outcome'] == 'would_merge'
    assert mr.evaluate_rule0(_r0_row(srd_phi_left=0.0, srd_phi_right=-4.0, delta_mean_left=0.1, delta_mean_right=0.9),
                             _rule0())['sign_mismatch'] is False


def test_rule0_signs_only_when_test_fails():
    ev = mr.evaluate_rule0(_r0_row(delta_mean_left=0.9, delta_mean_right=0.9, srd_phi_left=1.0, srd_phi_right=-1.0),
                           _rule0(tau=0.5))
    assert ev['passed'] is False and ev['outcome'] == 'signs_only' and ev['score'] is None


def test_rule0_flank_conflict_direction_and_switch():
    # delta_mean picks right (0.1/0.5 = 0.2), srd_phi picks left (0.33/3.3 = 0.1): conflict, smaller fraction is srd_phi -> left
    row = _r0_row(delta_mean_left=0.4, delta_mean_right=0.1, srd_phi_left=0.33, srd_phi_right=3.0)
    both = ['delta_mean', 'srd_phi']
    ev = mr.evaluate_rule0(row, _rule0(criteria=both))
    assert ev['flank_conflict'] is True and ev['held_back'] == ['flanks'] and ev['outcome'] == 'held_back'
    assert ev['side'] == 'left' and ev['score'] == pytest.approx(0.2)
    ev = mr.evaluate_rule0(row, _rule0(criteria=both, mergeFlankConflict=True))
    assert ev['outcome'] == 'would_merge' and ev['side'] == 'left'


def test_rule0_both_switches_needed_when_both_disagree():
    row = _r0_row(delta_mean_left=0.4, delta_mean_right=0.1, srd_phi_left=0.33, srd_phi_right=-3.0)
    both = ['delta_mean', 'srd_phi']
    assert mr.evaluate_rule0(row, _rule0(criteria=both))['held_back'] == ['signs', 'flanks']
    assert mr.evaluate_rule0(row, _rule0(criteria=both, mergeSignMismatch=True))['held_back'] == ['flanks']
    ev = mr.evaluate_rule0(row, _rule0(criteria=both, mergeSignMismatch=True, mergeFlankConflict=True))
    assert ev['outcome'] == 'would_merge'


def test_rule0_size_limit_and_gap_block():
    big = _r0_row(start=10, end=30, delta_mean_left=0.1, delta_mean_right=0.9)
    ev = mr.evaluate_rule0(big, _rule0(maxBins=10))
    assert ev['evaluated'] is False and ev['outcome'] == 'none' and ev['criteria'] == []
    gapped = _r0_row(delta_mean_left=0.1, delta_mean_right=0.9, gap_left=True)
    ev = mr.evaluate_rule0(gapped, _rule0())
    assert ev['outcome'] == 'would_merge' and ev['side'] == 'left' and ev['blocked'] == 'gap'
    assert mr.evaluate_rule0(gapped, _rule0(), allow_gap_crossing=True)['blocked'] is None


def test_evaluate_rules_rule0_claims_before_user_rules():
    row = _r0_row(start=10, end=15, delta_mean_left=0.1, delta_mean_right=0.9, srd_phi_left=None, srd_phi_right=1.0)
    user = _rule('small', [{'kind': 'max_size', 'maxBins': 10}], merge_metric='delta_log2fc')
    ev = mr.evaluate_rules([row], [_rule0(), user])[0]
    assert ev['matched'] == ['rule0', 'small'] and ev['claimed_by'] == 'rule0'
    assert ev['side'] == 'left' and ev['score'] == pytest.approx(0.2) and ev['rule0']['outcome'] == 'would_merge'
    held = _r0_row(start=10, end=15, delta_mean_left=0.1, delta_mean_right=0.9, srd_phi_left=1.0, srd_phi_right=-1.0)
    ev = mr.evaluate_rules([held], [_rule0(), user])[0]
    assert ev['matched'] == ['small'] and ev['claimed_by'] == 'small' and ev['rule0']['outcome'] == 'held_back'
    assert mr.evaluate_rules([row], [user])[0]['rule0'] is None


@pytest.mark.parametrize('rules', [
    [GOOD],
    [GOOD, INERT_RULE0],
    [INERT_RULE0, {**GOOD, 'id': 'x', 'type': 'flank_threshold'}],
    [INERT_RULE0, {**GOOD, 'id': 'rule0'}],
    [_rule0(criteria=[])],
    [_rule0(criteria=['bogus'])],
    [_rule0(criteria=['delta_mean', 'delta_mean'])],
    [_rule0(criteria='delta_mean')],
    [_rule0(combine='some')],
    [_rule0(tau=0)],
    [_rule0(tau=float('nan'))],
    [_rule0(tau=True)],
    [_rule0(srdRef=-1.0)],
    [_rule0(maxBins=0)],
    [_rule0(maxBins=2.5)],
    [_rule0(maxBins=True)],
    [_rule0(mergeSignMismatch='false')],
    [_rule0(mergeFlankConflict=None)],
    [_rule0(hatch={'pattern': '#', 'color': '#000'})],
    [{**_rule0(), 'type': 'conditions'}],
])
def test_validate_rules_rejects_bad_rule0(rules):
    with pytest.raises(ValueError):
        mr.validate_rules(rules)


def test_validate_rules_rule0_alone_is_valid_and_first():
    assert mr.validate_rules([INERT_RULE0]) == [INERT_RULE0]
    off = _rule('off', [], enabled=False)
    assert [r['id'] for r in mr.validate_rules([INERT_RULE0, off, GOOD])] == ['rule0', 'ok']


def _two_small_bumps(gap_at=None):
    # constant levels: phi = 0, so srd_phi is missing everywhere (Review Focus 2)
    profile = np.concatenate([np.full(10, 4.0), np.full(5, 4.5), np.full(10, 6.0), np.full(5, 6.3), np.full(10, 9.0)])
    var_start = (np.arange(40) * 100).astype(float)
    var_end = var_start + 100
    if gap_at is not None:
        var_end[gap_at - 1] += 500
    return profile, [10, 15, 25, 30], 40, var_start, var_end


def test_run_rule0_merges_by_fraction_with_constant_levels():
    profile, bounds, n, vs, ve = _two_small_bumps()
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [_rule0(tau=0.5, maxBins=5)])
    assert [s['removed_boundary'] for s in out['steps'][1:]] == [25, 10]
    s1 = out['steps'][1]
    assert (s1['rule_id'], s1['rule_name'], s1['side'], s1['segment']) == ('rule0', 'Rule 0', 'left', [25, 30])
    assert s1['score'] == pytest.approx(np.log2(6.3 / 6.0) / 0.5)
    assert out['final'] == [{'start': 0, 'end': 15}, {'start': 15, 'end': 30}, {'start': 30, 'end': 40}]


def test_run_rule0_gap_blocks_without_redirect():
    profile, bounds, n, vs, ve = _two_small_bumps(gap_at=25)
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [_rule0(tau=0.5, maxBins=5)])
    assert [s['removed_boundary'] for s in out['steps'][1:]] == [10]


def test_hatch_merge_route_runs_rule0_alone():
    resp = _client().post('/api/integration/hatch-merge', json={
        'chrom': 'chr1', 'source': 'cluster0', 'rules': [_rule0(tau=0.5)]})
    assert resp.status_code == 200
    body = resp.json()
    assert body['steps'][1]['rule_id'] == 'rule0'
    assert body['final'] == [{'start': 0, 'end': 20}, {'start': 20, 'end': 30}]
```

Then update the existing tests so each still tests what it was written for:

1. `test_validate_rules_rejects`: in its parametrize list, delete the case `[_rule('off', [WEAK_SRD_PHI], enabled=False)]` (a list with only a disabled user rule is now valid once Rule 0 is present). Keep `[]` and `'not a list'` as they are. Wrap every other case with `_with_r0(...)`, e.g. `_with_r0([_rule('empty', [])])` and `_with_r0([GOOD, _rule('ok', [WEAK_SRD_PHI])])`.
2. `test_validate_rules_ignores_disabled_and_returns_enabled_in_order`: call `mr.validate_rules(_with_r0([GOOD, broken_but_off, second]))` and expect `['rule0', 'ok', 'two']`.
3. `test_non_boolean_enabled_is_ignored`: wrap its rule list with `_with_r0(...)` and prepend `'rule0'` to the expected ids.
4. Every `mr.run_rule_merge(...)` call in the existing `test_run_rule_merge_*` tests: wrap its rules argument with `_with_r0(...)`. In `test_run_rule_merge_echoes_settings_and_step_zero_shape`, build `rules = _with_r0([_delta_rule('r', [WEAK_DELTA]), _rule('off', [], enabled=False)])`, keep `assert out['settings'] == {'rules': rules, 'allow_gap_crossing': False}`, and keep the `mr.run_rule_merge(..., [])` → `ValueError` check.
5. `test_hatch_merge_route_runs_rules_and_reports_rule_per_step`: wrap the rules with `_with_r0(...)`.
6. `test_hatch_merge_route_rejects_bad_rules`: delete the case `[_rule('off', [WEAK_DELTA], enabled=False)]`, keep `[]`, wrap every other case with `_with_r0(...)`.
7. Do not change `test_python_evaluator_matches_parity_fixture` here; Task 2 updates the fixture and that test.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest -q tests/test_merge_rules.py`
Expected: FAIL. The new tests fail with `AttributeError: module 'server.merge_rules' has no attribute 'evaluate_rule0'`, and the wrapped existing tests fail because Rule 0 is unknown to `validate_rules`. Also, `test_python_evaluator_matches_parity_fixture` passes for now: its fixture has no Rule 0, and Task 2 changes it.

- [ ] **Step 3: Implement Rule 0 in `server/merge_rules.py`**

After the `HATCH_PATTERNS` line, add:

```python
RULE0_ID = 'rule0'
RULE0_TYPE = 'flank_threshold'
RULE0_CRITERIA = ('delta_median', 'delta_mean', 'delta_iqr_mean', 'srd_phi')
RULE0_COMBINE = ('all', 'any')
```

After `_validate_condition`, add:

```python
def _positive(v):
    return _number(v) and math.isfinite(v) and v > 0


def _check_hatch(r, where):
    hatch = r.get('hatch') if isinstance(r.get('hatch'), dict) else {}
    if hatch.get('pattern') not in HATCH_PATTERNS:
        raise ValueError(f'{where}: unknown hatch pattern {hatch.get("pattern")!r}')


def _validate_rule0(r):
    where = 'Rule 0'
    criteria = r.get('criteria')
    if not isinstance(criteria, list) or not criteria:
        raise ValueError(f'{where}: check at least one criterion')
    if len(set(criteria)) != len(criteria) or any(c not in RULE0_CRITERIA for c in criteria):
        raise ValueError(f'{where}: criteria must be distinct values from {RULE0_CRITERIA}')
    if r.get('combine') not in RULE0_COMBINE:
        raise ValueError(f'{where}: combine must be one of {RULE0_COMBINE}')
    for key in ('tau', 'srdRef'):
        if not _positive(r.get(key)):
            raise ValueError(f'{where}: {key} must be a positive finite number')
    n = r.get('maxBins')
    if not (isinstance(n, int) and not isinstance(n, bool) and n >= 1):
        raise ValueError(f'{where}: maxBins must be an integer >= 1')
    for key in ('mergeSignMismatch', 'mergeFlankConflict'):
        if not isinstance(r.get(key), bool):
            raise ValueError(f'{where}: {key} must be true or false')
    _check_hatch(r, where)
```

Replace the whole `validate_rules` function with:

```python
def validate_rules(rules):
    """Validate and return [Rule 0, *enabled user rules], in list order.

    rules[0] must be Rule 0 (always active). Disabled user rules are ignored
    apart from the id checks, so a half-edited disabled rule never blocks a
    run. Raises ValueError on any problem."""
    if not isinstance(rules, list) or not all(isinstance(r, dict) for r in rules):
        raise ValueError('rules must be a list of objects')
    ids = [r.get('id') for r in rules]
    if any(not isinstance(i, str) or not i for i in ids):
        raise ValueError('every rule needs a non-empty string id')
    if len(set(ids)) != len(ids):
        raise ValueError('rule ids must be unique')
    if not rules or rules[0].get('id') != RULE0_ID or rules[0].get('type') != RULE0_TYPE:
        raise ValueError('rules[0] must be Rule 0 (id "rule0", type "flank_threshold")')
    for r in rules[1:]:
        if r.get('type', 'conditions') != 'conditions':
            raise ValueError(f'rule {r["id"]!r}: only rules[0] may be Rule 0')
    _validate_rule0(rules[0])
    enabled = [r for r in rules[1:] if r.get('enabled') is True]
    for r in enabled:
        where = f'rule {r["id"]!r}'
        conditions = r.get('conditions')
        if not isinstance(conditions, list) or not conditions:
            raise ValueError(f'{where}: needs at least one condition')
        for c in conditions:
            _validate_condition(c, where)
        merge_into = r.get('mergeInto') if isinstance(r.get('mergeInto'), dict) else {}
        _check_metric(merge_into.get('metric'), merge_into.get('estimator', 'mean'), f'{where} mergeInto')
        _check_hatch(r, where)
    return [rules[0], *enabled]
```

After `condition_holds`, add:

```python
def _criterion_value(row, side, criterion):
    v = row.get(f'{criterion}_{side}')
    return float(v) if v is not None and math.isfinite(v) else None


def _smaller_side(left, right):
    if left is None and right is None:
        return None
    if right is None:
        return 'left'
    if left is None:
        return 'right'
    return 'left' if abs(left) <= abs(right) else 'right'


def evaluate_rule0(row, r0, allow_gap_crossing=False):
    """Rule 0 on one score row: per-criterion flank/threshold checks, the
    all/any test, sign and flank-conflict checks, direction, ordering score
    and outcome. See docs/superpowers/specs/2026-10-26-rule0-visual-merge-design.md §3."""
    out = {'start': row['start'], 'end': row['end'], 'evaluated': False, 'criteria': [], 'passed': False,
           'sign_mismatch': False, 'flank_conflict': False, 'held_back': [], 'outcome': 'none',
           'side': None, 'score': None, 'blocked': None}
    if row['end'] - row['start'] > r0['maxBins']:
        return out
    out['evaluated'] = True
    for c in r0['criteria']:
        left, right = _criterion_value(row, 'left', c), _criterion_value(row, 'right', c)
        threshold = r0['srdRef'] if c == 'srd_phi' else r0['tau']
        side = _smaller_side(left, right)
        value = None if side is None else abs(left if side == 'left' else right)
        out['criteria'].append({'criterion': c, 'side': side, 'value': value, 'threshold': threshold,
                                'fraction': None if value is None else value / threshold,
                                'passed': value is not None and value < threshold})
    passing = [c for c in out['criteria'] if c['passed']]
    out['passed'] = len(passing) == len(out['criteria']) if r0['combine'] == 'all' else bool(passing)
    sl, sr = _criterion_value(row, 'left', 'srd_phi'), _criterion_value(row, 'right', 'srd_phi')
    out['sign_mismatch'] = sl is not None and sr is not None and sl * sr < 0
    if not out['passed']:
        if out['sign_mismatch']:
            out['outcome'] = 'signs_only'
        return out
    out['flank_conflict'] = len({c['side'] for c in passing}) > 1
    fractions = [c['fraction'] for c in passing]
    out['score'] = max(fractions) if r0['combine'] == 'all' else min(fractions)
    best = min(passing, key=lambda c: (c['fraction'], 0 if c['side'] == 'left' else 1))
    out['side'] = best['side'] if out['flank_conflict'] else passing[0]['side']
    if out['sign_mismatch'] and not r0['mergeSignMismatch']:
        out['held_back'].append('signs')
    if out['flank_conflict'] and not r0['mergeFlankConflict']:
        out['held_back'].append('flanks')
    if out['held_back']:
        out['outcome'] = 'held_back'
        return out
    out['outcome'] = 'would_merge'
    if not allow_gap_crossing and row[f'gap_{out["side"]}']:
        out['blocked'] = 'gap'
    return out
```

(For `all`, `passing` is every checked criterion, so the max over `passing` is the max over the checked criteria.)

Replace the whole `evaluate_rules` function with:

```python
def evaluate_rules(rows, rules, allow_gap_crossing=False):
    """Per score row: which rules match, which claims it, and the merge
    direction/score or why it cannot merge. `rules` must be validated; a
    leading Rule 0 is evaluated first and claims before user rules."""
    r0 = rules[0] if rules and rules[0].get('type') == RULE0_TYPE else None
    user = rules[1:] if r0 is not None else rules
    out = []
    for row in rows:
        r0ev = evaluate_rule0(row, r0, allow_gap_crossing) if r0 is not None else None
        matched = ([RULE0_ID] if r0ev is not None and r0ev['outcome'] == 'would_merge' else []) + \
            [r['id'] for r in user if all(condition_holds(c, row) for c in r['conditions'])]
        ev = {'start': row['start'], 'end': row['end'], 'matched': matched,
              'claimed_by': None, 'side': None, 'score': None, 'blocked': None, 'rule0': r0ev}
        if matched and matched[0] == RULE0_ID:
            ev.update(claimed_by=RULE0_ID, side=r0ev['side'], score=r0ev['score'], blocked=r0ev['blocked'])
        elif matched:
            rule = next(r for r in user if r['id'] == matched[0])
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

`run_rule_merge` needs no change (it orders by `order[claimed_by]`, and Rule 0 is index 0 of `validate_rules`' result). Update the module docstring's first paragraph to mention that `rules[0]` is Rule 0, the always-active rule built from the tracks' section-1 settings.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python3 -m pytest -q tests/test_merge_rules.py tests/test_hatch_merge.py`
Expected: all pass except `test_python_evaluator_matches_parity_fixture`, which now fails because its fixture rules lack Rule 0 (Task 2 fixes it). Confirm that this is the only failure.

- [ ] **Step 5: Commit**

```bash
git add server/merge_rules.py tests/test_merge_rules.py
git commit -m "feat: Rule 0 (flank-threshold rule) in the Python merge-rule engine

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: TypeScript Rule 0, track helpers, parity fixture

**Files:**
- Modify: `web/src/mergeRules.ts`, `tests/fixtures/merge_rules_parity.json`, `tests/merge_rules_helpers.ts`, `tests/test_merge_rules.py` (parity test only)
- Create: `web/src/trackMetrics.ts`, `tests/track_helpers.ts`

**Interfaces:**
- Consumes: Python `evaluate_rules` with `'rule0'` detail (Task 1); `HatchScoreRow` and `scoreForMetric` from `web/src/hatchMerge.ts`; `type Metric` from `web/src/integration.ts`.
- Produces (`web/src/mergeRules.ts`):
  - `type Rule0Criterion = 'delta_median'|'delta_mean'|'delta_iqr_mean'|'srd_phi'`, `RULE0_CRITERIA`, `CRITERION_LABEL: Record<Rule0Criterion, string>`
  - `interface Rule0 { id: 'rule0'; type: 'flank_threshold'; name: string; criteria: Rule0Criterion[]; combine: 'all'|'any'; tau: number; srdRef: number; maxBins: number; mergeSignMismatch: boolean; mergeFlankConflict: boolean; hatch: { pattern: HatchPattern; color: string } }`
  - `type AnyRule = Rule0 | MergeRule`, `isRule0(r): r is Rule0`
  - `interface CriterionEval { criterion; side: 'left'|'right'|null; value: number|null; threshold: number; fraction: number|null; passed: boolean }`
  - `interface Rule0Eval { start; end; evaluated; criteria: CriterionEval[]; passed; signMismatch; flankConflict; heldBack: ('signs'|'flanks')[]; outcome: 'none'|'would_merge'|'held_back'|'signs_only'; side; score; blocked: null|'gap' }`
  - `SegmentEval` gains `rule0: Rule0Eval | null`
  - `RULE0_COLOR = '#b91c1c'`, `SIGN_COLOR = '#d97706'`, `FLANK_COLOR = '#7c3aed'`
  - `evaluateRule0(row, r0, allowGap): Rule0Eval`, `rule0Problems(r0): string[]`, `describeRule0(r0): string`, `rule0Lines(e: Rule0Eval): string[]`
  - `rule0FromSeg(s: { criteria; combine; tau; srdRef; smallMax; mergeSignMismatch; mergeFlankConflict }): Rule0`
  - `rule0CaseCounts(evals: Record<string, SegmentEval[]>): { wouldMerge; signs; flanks; heldBack; signsOnly; gap }`
  - `evaluateRules(rows, rules: AnyRule[], allowGap)` and `evaluateForSources(..., rules: AnyRule[], ...)` (a leading Rule 0 is evaluated; an invalid Rule 0 is ignored, `rule0` null)
  - **changed:** `canRunMerge(rule0: Rule0, rules: MergeRule[]): boolean`
- Produces (`web/src/trackMetrics.ts`):
  - `type TrackMetric = 'delta_median'|'delta_mean'|'delta_iqr_mean'|'srd'|'srd_phi'`, `TRACK_METRICS` (id, label, kind, bins, binLabels)
  - `colourValue(row: HatchScoreRow, m: TrackMetric): { value: number|null; side: 'left'|'right'|null }`
  - `archivedColourToTrack(id: unknown): TrackMetric`, `diagnosticsEffectMetric(m: TrackMetric): Metric`
  - `migrateSegSettings(old: any): any`
  - `derivedRejectionReason(r: { selected_for_significance; accepted; rejection_reasons; p_local; p_global }): string | null`

- [ ] **Step 1: Update the parity fixture**

In `tests/fixtures/merge_rules_parity.json`:

(a) In each of the three existing cases, insert this object as the first element of `"rules"`:

```json
{"id": "rule0", "type": "flank_threshold", "name": "Rule 0", "criteria": ["srd_phi"], "combine": "all", "tau": 1e-9, "srdRef": 1e-9, "maxBins": 1, "mergeSignMismatch": false, "mergeFlankConflict": false, "hatch": {"pattern": "x", "color": "#b91c1c"}}
```

Their `"expected"` entries stay unchanged.

(b) Append these four cases to `"cases"`. Expectations were hand-derived from the five fixture rows. Each expected entry carries `"rule0"`; the parity tests compare the listed Rule 0 fields.

```json
{
  "name": "rule0_all_mean_srd_held_back_by_signs",
  "allowGapCrossing": false,
  "rules": [
    {"id": "rule0", "type": "flank_threshold", "name": "Rule 0", "criteria": ["delta_mean", "srd_phi"], "combine": "all", "tau": 0.2, "srdRef": 3.3, "maxBins": 100, "mergeSignMismatch": false, "mergeFlankConflict": false, "hatch": {"pattern": "x", "color": "#b91c1c"}}
  ],
  "expected": [
    {"start": 10, "end": 20, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": true, "passed": true, "signMismatch": true, "flankConflict": false, "heldBack": ["signs"], "outcome": "held_back", "side": "left", "score": 0.35, "blocked": null}},
    {"start": 20, "end": 25, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": true, "passed": true, "signMismatch": true, "flankConflict": false, "heldBack": ["signs"], "outcome": "held_back", "side": "right", "score": 0.85, "blocked": null}},
    {"start": 25, "end": 30, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": true, "passed": false, "signMismatch": false, "flankConflict": false, "heldBack": [], "outcome": "none", "side": null, "score": null, "blocked": null}},
    {"start": 30, "end": 40, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": true, "passed": false, "signMismatch": true, "flankConflict": false, "heldBack": [], "outcome": "signs_only", "side": null, "score": null, "blocked": null}},
    {"start": 40, "end": 41, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": true, "passed": false, "signMismatch": false, "flankConflict": false, "heldBack": [], "outcome": "none", "side": null, "score": null, "blocked": null}}
  ]
},
{
  "name": "rule0_any_signs_allowed_gap_blocks",
  "allowGapCrossing": false,
  "rules": [
    {"id": "rule0", "type": "flank_threshold", "name": "Rule 0", "criteria": ["delta_mean", "srd_phi"], "combine": "any", "tau": 0.2, "srdRef": 3.3, "maxBins": 100, "mergeSignMismatch": true, "mergeFlankConflict": false, "hatch": {"pattern": "x", "color": "#b91c1c"}}
  ],
  "expected": [
    {"start": 10, "end": 20, "matched": ["rule0"], "claimedBy": "rule0", "side": "left", "score": 0.15151515151515152, "blocked": null,
     "rule0": {"evaluated": true, "passed": true, "signMismatch": true, "flankConflict": false, "heldBack": [], "outcome": "would_merge", "side": "left", "score": 0.15151515151515152, "blocked": null}},
    {"start": 20, "end": 25, "matched": ["rule0"], "claimedBy": "rule0", "side": "right", "score": 0.06060606060606061, "blocked": "gap",
     "rule0": {"evaluated": true, "passed": true, "signMismatch": true, "flankConflict": false, "heldBack": [], "outcome": "would_merge", "side": "right", "score": 0.06060606060606061, "blocked": "gap"}},
    {"start": 25, "end": 30, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": true, "passed": false, "signMismatch": false, "flankConflict": false, "heldBack": [], "outcome": "none", "side": null, "score": null, "blocked": null}},
    {"start": 30, "end": 40, "matched": ["rule0"], "claimedBy": "rule0", "side": "left", "score": 0.6060606060606061, "blocked": null,
     "rule0": {"evaluated": true, "passed": true, "signMismatch": true, "flankConflict": false, "heldBack": [], "outcome": "would_merge", "side": "left", "score": 0.6060606060606061, "blocked": null}},
    {"start": 40, "end": 41, "matched": ["rule0"], "claimedBy": "rule0", "side": "left", "score": 0.30303030303030304, "blocked": null,
     "rule0": {"evaluated": true, "passed": true, "signMismatch": false, "flankConflict": false, "heldBack": [], "outcome": "would_merge", "side": "left", "score": 0.30303030303030304, "blocked": null}}
  ]
},
{
  "name": "rule0_all_median_srd_flank_conflict_held",
  "allowGapCrossing": false,
  "rules": [
    {"id": "rule0", "type": "flank_threshold", "name": "Rule 0", "criteria": ["delta_median", "srd_phi"], "combine": "all", "tau": 2.0, "srdRef": 3.3, "maxBins": 100, "mergeSignMismatch": true, "mergeFlankConflict": false, "hatch": {"pattern": "x", "color": "#b91c1c"}}
  ],
  "expected": [
    {"start": 10, "end": 20, "matched": ["rule0"], "claimedBy": "rule0", "side": "left", "score": 0.15151515151515152, "blocked": null,
     "rule0": {"evaluated": true, "passed": true, "signMismatch": true, "flankConflict": false, "heldBack": [], "outcome": "would_merge", "side": "left", "score": 0.15151515151515152, "blocked": null}},
    {"start": 20, "end": 25, "matched": ["rule0"], "claimedBy": "rule0", "side": "right", "score": 0.1, "blocked": "gap",
     "rule0": {"evaluated": true, "passed": true, "signMismatch": true, "flankConflict": false, "heldBack": [], "outcome": "would_merge", "side": "right", "score": 0.1, "blocked": "gap"}},
    {"start": 25, "end": 30, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": true, "passed": false, "signMismatch": false, "flankConflict": false, "heldBack": [], "outcome": "none", "side": null, "score": null, "blocked": null}},
    {"start": 30, "end": 40, "matched": ["rule0"], "claimedBy": "rule0", "side": "left", "score": 0.6060606060606061, "blocked": null,
     "rule0": {"evaluated": true, "passed": true, "signMismatch": true, "flankConflict": false, "heldBack": [], "outcome": "would_merge", "side": "left", "score": 0.6060606060606061, "blocked": null}},
    {"start": 40, "end": 41, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": true, "passed": true, "signMismatch": false, "flankConflict": true, "heldBack": ["flanks"], "outcome": "held_back", "side": "left", "score": 0.5, "blocked": null}}
  ]
},
{
  "name": "rule0_conflict_allowed_size_limit_user_rule_after",
  "allowGapCrossing": false,
  "rules": [
    {"id": "rule0", "type": "flank_threshold", "name": "Rule 0", "criteria": ["delta_median", "srd_phi"], "combine": "all", "tau": 2.0, "srdRef": 3.3, "maxBins": 5, "mergeSignMismatch": true, "mergeFlankConflict": true, "hatch": {"pattern": "x", "color": "#b91c1c"}},
    {"id": "u", "name": "u", "enabled": true, "conditions": [{"kind": "max_size", "maxBins": 5}], "mergeInto": {"metric": "srd_phi", "estimator": "mean"}, "hatch": {"pattern": ".", "color": "#047857"}}
  ],
  "expected": [
    {"start": 10, "end": 20, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": false, "passed": false, "signMismatch": false, "flankConflict": false, "heldBack": [], "outcome": "none", "side": null, "score": null, "blocked": null}},
    {"start": 20, "end": 25, "matched": ["rule0", "u"], "claimedBy": "rule0", "side": "right", "score": 0.1, "blocked": "gap",
     "rule0": {"evaluated": true, "passed": true, "signMismatch": true, "flankConflict": false, "heldBack": [], "outcome": "would_merge", "side": "right", "score": 0.1, "blocked": "gap"}},
    {"start": 25, "end": 30, "matched": ["u"], "claimedBy": "u", "side": null, "score": null, "blocked": "nonfinite",
     "rule0": {"evaluated": true, "passed": false, "signMismatch": false, "flankConflict": false, "heldBack": [], "outcome": "none", "side": null, "score": null, "blocked": null}},
    {"start": 30, "end": 40, "matched": [], "claimedBy": null, "side": null, "score": null, "blocked": null,
     "rule0": {"evaluated": false, "passed": false, "signMismatch": false, "flankConflict": false, "heldBack": [], "outcome": "none", "side": null, "score": null, "blocked": null}},
    {"start": 40, "end": 41, "matched": ["rule0", "u"], "claimedBy": "rule0", "side": "left", "score": 0.5, "blocked": null,
     "rule0": {"evaluated": true, "passed": true, "signMismatch": false, "flankConflict": true, "heldBack": [], "outcome": "would_merge", "side": "left", "score": 0.5, "blocked": null}}
  ]
}
```

- [ ] **Step 2: Extend the Python parity test**

In `tests/test_merge_rules.py`, replace the body of `test_python_evaluator_matches_parity_fixture` with:

```python
    fixture = json.loads(FIXTURE.read_text())
    keys = (('evaluated', 'evaluated'), ('passed', 'passed'), ('sign_mismatch', 'signMismatch'),
            ('flank_conflict', 'flankConflict'), ('held_back', 'heldBack'), ('outcome', 'outcome'),
            ('side', 'side'), ('blocked', 'blocked'))
    for case in fixture['cases']:
        got = mr.evaluate_rules(fixture['rows'], mr.validate_rules(case['rules']), case['allowGapCrossing'])
        assert len(got) == len(case['expected']), case['name']
        for g, exp in zip(got, case['expected']):
            assert (g['start'], g['end'], g['matched'], g['claimed_by'], g['side'], g['blocked']) == \
                (exp['start'], exp['end'], exp['matched'], exp['claimedBy'], exp['side'], exp['blocked']), case['name']
            assert (g['score'] is None) == (exp['score'] is None), case['name']
            if exp['score'] is not None:
                assert g['score'] == pytest.approx(exp['score']), case['name']
            if 'rule0' in exp:
                for py_key, js_key in keys:
                    assert g['rule0'][py_key] == exp['rule0'][js_key], (case['name'], g['start'], py_key)
                assert (g['rule0']['score'] is None) == (exp['rule0']['score'] is None), case['name']
                if exp['rule0']['score'] is not None:
                    assert g['rule0']['score'] == pytest.approx(exp['rule0']['score']), case['name']
```

Run: `python3 -m pytest -q tests/test_merge_rules.py::test_python_evaluator_matches_parity_fixture`
Expected: PASS. Task 1 already implemented the Python side, so this pins the fixture. If it fails, recheck the fixture values against the hand derivation in `evaluate_rule0`. Fix the fixture only if the code matches the spec and the fixture value was mis-derived, and state which in the report.

- [ ] **Step 3: Write the failing TS tests**

In `tests/merge_rules_helpers.ts`:

(a) Change the import line to:

```ts
import { canRunMerge, defaultMergeInto, evaluateRules, makeRule, migrateHatchLayers, ruleProblems, rule0FromSeg, rule0Problems, rule0CaseCounts, describeRule0, type AnyRule, type MergeRule, type Rule0 } from '../web/src/mergeRules'
```

(b) In the parity loop, cast `c.rules as AnyRule[]`, and after the existing `same` comparison add a Rule 0 comparison when `exp.rule0` is present:

```ts
    if (exp.rule0) {
      const r = g.rule0
      const r0same = !!r && r.evaluated === exp.rule0.evaluated && r.passed === exp.rule0.passed && r.signMismatch === exp.rule0.signMismatch &&
        r.flankConflict === exp.rule0.flankConflict && JSON.stringify(r.heldBack) === JSON.stringify(exp.rule0.heldBack) &&
        r.outcome === exp.rule0.outcome && r.side === exp.rule0.side && r.blocked === exp.rule0.blocked && close(r.score, exp.rule0.score)
      if (!r0same) { failures++; console.error(`FAIL parity ${c.name} row ${i} rule0: got ${JSON.stringify(r)}, expected ${JSON.stringify(exp.rule0)}`) }
      else console.log(`ok parity ${c.name} row ${i} rule0`)
    }
```

(c) Replace the four existing `canRunMerge` lines with:

```ts
const r0 = rule0FromSeg({ criteria: ['delta_mean', 'srd_phi'], combine: 'all', tau: 0.2, srdRef: 3.3, smallMax: 100, mergeSignMismatch: false, mergeFlankConflict: false })
eq(canRunMerge(r0, []), true, 'Rule 0 alone can run')
eq(canRunMerge(r0, [{ ...base, enabled: false }]), true, 'disabled user rules do not block')
eq(canRunMerge(r0, [base, { ...base, id: 'b', conditions: [] }]), false, 'an incomplete enabled rule blocks the run')
eq(canRunMerge(r0, [base, { ...base, id: 'b', enabled: false, conditions: [] }]), true, 'an incomplete disabled rule does not block')
eq(canRunMerge({ ...r0, criteria: [] }, []), false, 'Rule 0 with no criterion blocks the run')
```

(d) Before the final `if (failures)` line add:

```ts
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
eq(rule0CaseCounts({ a: evaluateRules(fx2.rows, anyCase.rules, false) }), { wouldMerge: 4, signs: 3, flanks: 0, heldBack: 0, signsOnly: 0, gap: 1 }, 'case counts (any, signs allowed)')
eq(rule0CaseCounts({ a: evaluateRules(fx2.rows, heldCase.rules, false) }), { wouldMerge: 0, signs: 0, flanks: 0, heldBack: 2, signsOnly: 1, gap: 0 }, 'case counts (held back)')
eq(describeRule0(r0), 'all of: min |Δlog2FC mean| flank < 0.2, min |SRD/√φ| flank < 3.3 · segments ≤ 100 bins · signs disagree: held back · flanks disagree: held back', 'describeRule0')
```

Create `tests/track_helpers.ts`:

```ts
// tests/track_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/track_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/track_helpers.mjs
//           node tests/_bundle/track_helpers.mjs
import { archivedColourToTrack, colourValue, derivedRejectionReason, diagnosticsEffectMetric, migrateSegSettings } from '../web/src/trackMetrics'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}
const row: any = { start: 10, end: 15, phi: 1, gap_left: false, gap_right: false,
  srd_left: -4, srd_right: 1.5, srd_phi_left: null, srd_phi_right: 2.0,
  delta_mean_left: -0.3, delta_mean_right: 0.1, delta_median_left: 0.2, delta_median_right: -0.2,
  delta_iqr_mean_left: null, delta_iqr_mean_right: null }

eq(colourValue(row, 'delta_mean'), { value: 0.1, side: 'right' }, 'effect colour: signed value on smaller flank')
eq(colourValue(row, 'delta_median'), { value: 0.2, side: 'left' }, 'effect colour tie goes left')
eq(colourValue(row, 'srd'), { value: 1.5, side: 'right' }, 'stat colour: absolute value on smaller flank')
eq(colourValue(row, 'srd_phi'), { value: 2.0, side: 'right' }, 'missing side loses')
eq(colourValue(row, 'delta_iqr_mean'), { value: null, side: null }, 'both missing -> N/A')

eq(['delta_median', 'delta_mean', 'log2fc', 'delta_trim', 'srd', 'srd_phi', 'bogus', undefined].map(archivedColourToTrack),
  ['delta_median', 'delta_mean', 'delta_mean', 'delta_iqr_mean', 'srd', 'srd_phi', 'delta_mean', 'delta_mean'], 'archived colour ids map to live')
eq(['delta_median', 'delta_mean', 'delta_iqr_mean', 'srd', 'srd_phi'].map((m) => diagnosticsEffectMetric(m as any)),
  ['delta_median', 'delta_mean', 'delta_trim', 'delta_mean', 'delta_mean'], 'diagnostics effect metric mapping')

const old = { chrom: 'chr9', metric: 'log2fc', effectMetric: 'delta_trim', flank: 'min_srd', tau: 0.3, labels: { consistency: false, proposal: true, cn: true } }
eq(migrateSegSettings(old), { chrom: 'chr9', tau: 0.3, labels: { cn: true, signs: false }, colorMetric: 'delta_mean' }, 'old session migrates')
eq(migrateSegSettings({ colorMetric: 'srd', labels: { signs: true } }), { colorMetric: 'srd', labels: { signs: true } }, 'new session passes through')
eq(migrateSegSettings(undefined), {}, 'missing seg -> empty')

const base = { selected_for_significance: true, accepted: false, rejection_reasons: null, p_local: 0.002997, p_global: null }
eq(derivedRejectionReason(base), 'local permutation p = 0.003 ≥ 0.001 (cutoff inferred: all accepted candidates have p < 0.001)', 'local p reason')
eq(derivedRejectionReason({ ...base, p_local: 0.000999, p_global: 0.0123 }), 'global permutation p = 0.0123 ≥ 0.001 (cutoff inferred: all accepted candidates have p < 0.001)', 'global p reason')
eq(derivedRejectionReason({ ...base, p_local: 0.0005, p_global: 0.0005 }), 'not accepted (no recorded reason) (cutoff inferred: all accepted candidates have p < 0.001)', 'no reason found')
eq(derivedRejectionReason({ ...base, rejection_reasons: 'left_child_too_short' }), 'left_child_too_short', 'recorded reason kept')
eq(derivedRejectionReason({ ...base, accepted: true }), null, 'accepted -> no reason')
eq(derivedRejectionReason({ ...base, selected_for_significance: false }), null, 'not tested and no recorded reason -> null')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) }
console.log('all track_helpers checks passed')
```

Run: `web/node_modules/.bin/esbuild tests/track_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/track_helpers.mjs`
Expected: FAIL — `Could not resolve "../web/src/trackMetrics"`.

- [ ] **Step 4: Implement Rule 0 in `web/src/mergeRules.ts`**

1. Change `MergeRule` to carry an optional type: `export interface MergeRule { id: string; type?: 'conditions'; name: string; ... }` (other fields unchanged).
2. Add `rule0: Rule0Eval | null` to `SegmentEval`.
3. After the `SegmentEval` interface add:

```ts
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
```

4. After `conditionDetail`, add:

```ts
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
```

5. Replace `canRunMerge`, `evaluateRules` and `evaluateForSources` with:

```ts
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
```

- [ ] **Step 5: Create `web/src/trackMetrics.ts`**

```ts
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
```

Check: in the `migrateSegSettings` test, the expected object key order is `chrom, tau, labels, colorMetric`. JSON comparison is order-sensitive. The code adds `colorMetric` after the spread and replaces `labels` in place, so `labels` keeps its original position. With `old = {chrom, metric, effectMetric, flank, tau, labels}`, the spread yields `chrom, metric, effectMetric, flank, tau, labels`. Then `colorMetric` is appended, the three keys are deleted, and `labels` is reassigned in place. The result order is `chrom, tau, labels, colorMetric`, which matches. Inside labels, `cn` is kept and `signs` is appended, giving `{cn, signs}`, which also matches.

- [ ] **Step 6: Run all TS helper tests and the parity test**

Run:
```bash
web/node_modules/.bin/esbuild tests/track_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/track_helpers.mjs && node tests/_bundle/track_helpers.mjs
web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs && node tests/_bundle/merge_rules_helpers.mjs
web/node_modules/.bin/esbuild tests/hatch_merge_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/hatch_merge_helpers.mjs && node tests/_bundle/hatch_merge_helpers.mjs
python3 -m pytest -q tests/test_merge_rules.py tests/test_hatch_merge.py
```
Expected: the three helper scripts print `all … checks passed`, and pytest passes in full.

`web` does not compile yet. `HatchMergeView.tsx` and `MergeRulesEditor.tsx` still call the old `canRunMerge(rules)`, and Task 4 fixes that. Do not run `npm run build` in this task. Do run `cd web && npx tsc --noEmit`, and confirm that every error it prints is in `HatchMergeView.tsx`, `MergeRulesEditor.tsx`, or `IntegrationSegmentationView.tsx` and is about `canRunMerge`/`SegmentEval.rule0`. Paste the error list into the report.

- [ ] **Step 7: Add the generated bundle to `.gitignore` and commit**

Append `tests/_bundle/track_helpers.mjs` to `.gitignore`, next to the two helper bundles already listed there.

```bash
git add .gitignore web/src/mergeRules.ts web/src/trackMetrics.ts tests/fixtures/merge_rules_parity.json tests/merge_rules_helpers.ts tests/track_helpers.ts tests/test_merge_rules.py
git commit -m "feat: TypeScript Rule 0 evaluator, track metric helpers, parity fixture cases

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Record archived-vs-live agreement on the pilot chromosomes

**Files:**
- Create: `tests/compare_archived_live.py`
- Create (generated, commit): `validation/rule0/archived_vs_live.json`
- Create: `validation/rule0/README.md`

**Interfaces:**
- Consumes: `GET /api/integration/meta` (`pilot_chromosomes`, `sources`, `flank_scores`, `screening`), `POST /api/integration/hatch-scores`. Archived `bin_start`/`bin_end` are chromosome-local and half-open, matching live `start`/`end` (verified: `bin_end - bin_start == seg_length` for all 162 archived rows).
- Produces: the agreement record the spec §2 requires before switching the tracks to live scores.

- [ ] **Step 1: Write `tests/compare_archived_live.py`**

```python
"""Compare the archived small-segment flank tables with live hatch-scores on
the pilot chromosomes (spec 2026-10-26 §2). Needs a running server:
    BASE_URL=http://127.0.0.1:8767 python3 tests/compare_archived_live.py
Writes validation/rule0/archived_vs_live.json and prints a summary."""
import json
import math
import os
import sys
import urllib.request
from pathlib import Path

BASE = os.environ.get('BASE_URL', 'http://127.0.0.1:8766')
OUT = Path(__file__).resolve().parents[1] / 'validation' / 'rule0' / 'archived_vs_live.json'

# (archived table, archived column stem, live column stem)
PAIRS = [
    ('flank_scores', 'delta', 'delta_median'),
    ('flank_scores', 'delta_mean', 'delta_mean'),
    ('flank_scores', 'delta_trim', 'delta_iqr_mean'),
    ('flank_scores', 'srd_phi', 'srd_phi'),
    ('screening', 'log2fc', 'delta_mean'),
    ('screening', 'srd', 'srd'),
]


def get(path):
    with urllib.request.urlopen(BASE + path, timeout=300) as r:
        return json.loads(r.read())


def post(path, body):
    req = urllib.request.Request(BASE + path, data=json.dumps(body).encode(), headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.loads(r.read())


def num(v):
    return None if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)


def main():
    meta = get('/api/integration/meta')
    pilots, sources = meta['pilot_chromosomes'], meta['sources']
    live = {}
    for chrom in pilots:
        for src in sources:
            res = post('/api/integration/hatch-scores', {'chrom': chrom, 'source': src})
            for row in res['rows']:
                live[(src, chrom, row['start'], row['end'])] = row
    report = {'base_url': BASE, 'pilot_chromosomes': pilots, 'pairs': []}
    for table, a_stem, l_stem in PAIRS:
        rows = meta.get(table) or []
        matched = unmatched = compared = 0
        max_diff, worst = 0.0, None
        for r in rows:
            key = (r['source'], r['chromosome'], int(r['bin_start']), int(r['bin_end']))
            lr = live.get(key)
            if lr is None:
                unmatched += 1
                continue
            matched += 1
            for side in ('left', 'right'):
                a, b = num(r.get(f'{a_stem}_{side}')), num(lr.get(f'{l_stem}_{side}'))
                if a is None or b is None:
                    continue
                compared += 1
                d = abs(a - b)
                if d > max_diff:
                    max_diff, worst = d, {'key': list(key), 'side': side, 'archived': a, 'live': b}
        report['pairs'].append({'archived': f'{table}.{a_stem}', 'live': l_stem, 'archived_rows': len(rows),
                                'matched_segments': matched, 'unmatched_segments': unmatched,
                                'compared_values': compared, 'max_abs_diff': max_diff, 'worst': worst})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + '\n')
    for p in report['pairs']:
        print(f"{p['archived']:>26} vs {p['live']:<15} matched {p['matched_segments']:>4}/{p['archived_rows']:<4} "
              f"values {p['compared_values']:>4} max|diff| {p['max_abs_diff']:.3g}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
```

- [ ] **Step 2: Run it against a fresh server on port 8767**

```bash
source .venv-scientific-explorer/bin/activate
nohup python3 run.py --port 8767 > /tmp/claude-1000/rule0-compare-8767.log 2>&1 &
echo $! > /tmp/claude-1000/rule0-compare-8767.pid
until curl -sf http://127.0.0.1:8767/api/integration/status | grep -q '"done"'; do sleep 3; done
BASE_URL=http://127.0.0.1:8767 python3 tests/compare_archived_live.py
kill "$(cat /tmp/claude-1000/rule0-compare-8767.pid)"
```

Expected: six summary lines, and `validation/rule0/archived_vs_live.json` written. If the `status` poll never returns `"done"`, read the log and report rather than proceeding. Stop only the PID you started, never the process on port 8766.

- [ ] **Step 3: Write `validation/rule0/README.md`**

Record the following, using the real numbers from the JSON:
- the command
- one line per pair: archived column, live column, matched/archived rows, compared values, and max |diff|
- **Conclusion:** state which pairs agree to floating-point noise (max |diff| < 1e-6) and which differ.
  - If any pair differs, state that track values on pilot chromosomes changed when the tracks moved to live scores, and give that pair's `worst` example.
  - If segments are unmatched, state that the archived tables cover segments the current live breakpoints do not reproduce, and give the counts.

Keep the wording descriptive ("agree", "differ by up to"), not "validated".

- [ ] **Step 4: Commit**

```bash
git add tests/compare_archived_live.py validation/rule0/archived_vs_live.json validation/rule0/README.md
git commit -m "test: record archived vs live flank-score agreement on pilot chromosomes

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Frontend — live tracks, Rule 0 controls and hatching, pinned Rule 0 card

**Files:**
- Modify: `web/src/hatchMerge.ts`, `web/src/integrationStore.ts`, `web/src/components/IntegrationSegmentationView.tsx`, `web/src/components/HatchMergeView.tsx`, `web/src/components/MergeRulesEditor.tsx`, `web/src/styles.css`

**Interfaces:**
- Consumes (Task 2): `Rule0`, `Rule0Criterion`, `RULE0_CRITERIA`, `CRITERION_LABEL`, `RULE0_COLOR`, `SIGN_COLOR`, `FLANK_COLOR`, `AnyRule`, `isRule0`, `evaluateForSources`, `rule0FromSeg`, `rule0Problems`, `describeRule0`, `rule0Lines`, `rule0CaseCounts`, `canRunMerge(rule0, rules)`, `ruleFlagCounts`, `conditionDetail`, `activeRules` from `mergeRules.ts`; `TRACK_METRICS`, `TrackMetric`, `colourValue`, `diagnosticsEffectMetric`, `migrateSegSettings`, `derivedRejectionReason` from `trackMetrics.ts`.
- Produces: accessible names Task 5 uses:
  - heading text `1. Visualization and Rule 0`
  - checkboxes `Rule 0 criterion min |Δlog2FC median| flank`, `Rule 0 criterion min |Δlog2FC mean| flank`, `Rule 0 criterion min |Δlog2FC IQR-mean| flank`, `Rule 0 criterion min |SRD/√φ| flank` (aria-label is `Rule 0 criterion ` + `CRITERION_LABEL[c]`)
  - select `Rule 0 combine` (options `all`, `any`)
  - inputs `τ` (aria-label `Rule 0 tau`) and `SRD/√φ ref` (aria-label `Rule 0 SRD/√φ ref`)
  - checkboxes `Merge sign-disagreeing segments` and `Merge flank-conflict segments`
  - status `Rule 0 status` (text `→ flags N segment(s) on <chrom> across the displayed rows`, or the problems), a legend element with class `rule0-legend`, and the archived layers list without the SRD-screen proposal entry.
  - The existing user-rule names (`Add rule`, `Rule N …`, `Run merge`, `Allow merging across a genomic gap`) are unchanged.

- [ ] **Step 1: Types — `web/src/hatchMerge.ts`**

Change `import type { MergeRule } from './mergeRules'` to `import type { AnyRule } from './mergeRules'`. Set `MergeSettings.rules: AnyRule[]` and `MergeRequest.rules: AnyRule[]`.

- [ ] **Step 2: Store — `web/src/integrationStore.ts`**

- Imports:
  - Replace `type FlankRule, type Metric` in the `./integration` import with nothing (keep `buildSmallSegments`, `SmallSeg`).
  - Add `import { migrateSegSettings, type TrackMetric } from './trackMetrics'`.
  - Add `type Rule0Criterion` to the `./mergeRules` import.
- `LabelLayers`: remove `proposal`, rename `consistency` to `signs`.
- `SegSettings`: replace `metric: Metric; effectMetric: Metric; flank: FlankRule` with `colorMetric: TrackMetric; criteria: Rule0Criterion[]; combine: 'all' | 'any'; mergeSignMismatch: boolean; mergeFlankConflict: boolean`.
- `defaultSeg`: replace `metric: 'delta_mean', effectMetric: 'delta_mean', flank: 'min_delta'` with `colorMetric: 'delta_mean', criteria: ['delta_mean', 'srd_phi'], combine: 'all', mergeSignMismatch: false, mergeFlankConflict: false`. In `labels`, remove `proposal: false` and replace `consistency: true` with `signs: true`.
- `restore`: replace the `seg:` expression with `seg: (() => { const ms = migrateSegSettings(p.seg); return { ...defaultSeg, ...ms, labels: { ...defaultSeg.labels, ...(ms.labels || {}) } } })()`.

- [ ] **Step 3: Segmentation view — `web/src/components/IntegrationSegmentationView.tsx`**

Make these changes in order. Locate each one by its quoted content.

1. **Imports:**
   - From `../integration`, remove `consistencyClass` and `metricValue`.
   - Delete the `PROPOSAL_LABEL` and `PROPOSAL_SHORT` constants. Keep `REASON_LABEL` and `REASON_SHORT`.
   - Replace the `../mergeRules` import with `import { CRITERION_LABEL, FLANK_COLOR, RULE0_CRITERIA, SIGN_COLOR, activeRules, conditionDetail, evaluateForSources, rule0CaseCounts, rule0FromSeg, rule0Lines, rule0Problems } from '../mergeRules'`.
   - Add `import { TRACK_METRICS, colourValue, derivedRejectionReason, diagnosticsEffectMetric } from '../trackMetrics'`.
2. **Hatch-score fetching and evaluation.** Replace the block from `const rulesOn = activeRules(ig.mergeRules).length > 0` through `const ruleCounts = useMemo(...)` with:
```tsx
  useEffect(() => { if (meta) for (const s of seg.sources) ig.ensureHatchScores(s, seg.chrom) }, [seg.chrom, seg.sources, meta])
  const rule0 = useMemo(() => rule0FromSeg(seg), [seg.criteria, seg.combine, seg.tau, seg.srdRef, seg.smallMax, seg.mergeSignMismatch, seg.mergeFlankConflict])
  const rule0Issues = rule0Problems(rule0)
  const allRules = useMemo(() => [rule0, ...ig.mergeRules], [rule0, ig.mergeRules])
  const ruleEvals = useMemo(() => evaluateForSources(ig.hatchScores, seg.sources, seg.chrom, allRules, ig.allowGapCrossing), [ig.hatchScores, seg.sources, seg.chrom, allRules, ig.allowGapCrossing])
  const r0Counts = useMemo(() => rule0CaseCounts(ruleEvals), [ruleEvals])
  const userRulesOn = activeRules(ig.mergeRules).length > 0
```
3. **Metric definitions.** Replace `const metricDef = METRICS.find((m) => m.id === seg.metric)!` and the `effectMetric` line with:
```tsx
  const trackDef = TRACK_METRICS.find((m) => m.id === seg.colorMetric) ?? TRACK_METRICS[1]
  const effectMetric: Metric = diagnosticsEffectMetric(trackDef.id)
```
   `effectMetric` keeps feeding the archived diagnostics panel unchanged. Replace every remaining `metricDef` use inside the tracks figure and the legend with `trackDef`. The diagnostics panel never used `metricDef`.
4. **Track drawing.** Delete the `pilot` constant and the `// small segments` block (the whole `if (seg.labels.small) for (const s of smallHere...) { ... }`). Also delete the `// merge-rule hatches` block (the whole `if (hs && evals && rulesOn) hs.rows.forEach(...)`). Replace both with this single block, placed where the small-segments block was:
```tsx
      // small internal segments from live hatch-scores: colour, Rule 0 cases, user-rule hatches
      const hs = ig.hatchScores[`${src}|${seg.chrom}`]
      const evals = ruleEvals[src]
      if (hs && evals) hs.rows.forEach((row, ri) => {
        if (row.end - row.start > seg.smallMax) return
        const ev = evals[ri], r0 = ev.rule0
        const key = `${src}|${seg.chrom}|${row.start}`
        const x0 = mb[row.start], x1 = row.end < chromMeta.n ? mb[row.end] : bins[2 * (chromMeta.offset + chromMeta.n - 1) + 1] / MB
        const cv = colourValue(row, trackDef.id)
        const color = seg.scale === 'binned' ? binColor(cv.value, trackDef.bins) : (cv.value === null ? '#e5e7eb' : viridisAt(Math.min(1, Math.abs(cv.value) / (trackDef.kind === 'stat' ? 8 : 2))))
        const selected = seg.selectedSegment === key
        if (seg.labels.small) shapes.push({ type: 'rect', xref: 'x', yref: `${yref} domain`, x0, x1, y0: 0, y1: 1, fillcolor: color, opacity: selected ? 0.75 : 0.5, line: { width: selected ? 2.5 : 0.6, color: '#111827' }, layer: 'below' })
        const merging = !!r0 && r0.outcome === 'would_merge'
        const held = !!r0 && r0.outcome === 'held_back'
        if (merging) shapes.push(...hatchShapes(x0, x1, yref, 'x', rule0.hatch.color))
        if (r0 && r0.signMismatch && (merging || held || (r0.outcome === 'signs_only' && seg.labels.signs))) shapes.push(...hatchShapes(x0, x1, yref, '|', SIGN_COLOR))
        if (r0 && r0.flankConflict && (merging || r0.heldBack.includes('flanks'))) shapes.push(...hatchShapes(x0, x1, yref, '-', FLANK_COLOR))
        if (merging && r0!.blocked === 'gap') shapes.push({ type: 'rect', xref: 'x', yref: `${yref} domain`, x0, x1, y0: 0.02, y1: 0.98, fillcolor: 'rgba(0,0,0,0)', line: { width: 1.5, color: rule0.hatch.color, dash: 'dot' }, layer: 'above' })
        for (const id of ev.matched) if (id !== 'rule0') { const rule = ig.mergeRules.find((r) => r.id === id); if (rule) shapes.push(...hatchShapes(x0, x1, yref, rule.hatch.pattern, rule.hatch.color)) }
        if (row.gap_left) shapes.push({ type: 'line', xref: 'x', yref: `${yref} domain`, x0, x1: x0, y0: 0, y1: 1, line: { color: '#b91c1c', width: 2, dash: 'dot' }, opacity: 0.7 })
        if (row.gap_right) shapes.push({ type: 'line', xref: 'x', yref: `${yref} domain`, x0: x1, x1, y0: 0, y1: 1, line: { color: '#b91c1c', width: 2, dash: 'dot' }, opacity: 0.7 })
        const archived = ig.smallSegs.find((x) => x.key === key)
        const labels: string[] = []
        if (seg.labels.reason && archived) labels.push(`${REASON_SHORT[archived.reasonL] || archived.reasonL}|${REASON_SHORT[archived.reasonR] || archived.reasonR}`)
        if (seg.labels.segLabels) labels.push(`${row.end - row.start}b ${fmt(cv.value, 2)}`)
        if (labels.length) annotations.push({ xref: 'x', yref, x: (x0 + x1) / 2, y: seg.basis === 'raw' ? 0 : yr[1] - 0.15, yanchor: 'top', xanchor: 'center', textangle: -90, showarrow: false, font: { size: 9, color: '#111827' }, bgcolor: 'rgba(255,255,255,0.7)', text: labels.join(' · ') })
        const active = new Set(activeRules(ig.mergeRules).map((r) => r.id))
        const userLines = ig.mergeRules.map((r, k) => ({ r, n: k + 1 })).filter(({ r }) => active.has(r.id) && ev.matched.includes(r.id))
          .map(({ r, n }) => `rule ${n} “${r.name}” matches<br>${r.conditions.map((c) => `&nbsp;&nbsp;${conditionDetail(c, row)}`).join('<br>')}`)
        const claimNote = ev.claimedBy && ev.claimedBy !== 'rule0' ? `<br>merge claimed by rule “${ig.mergeRules.find((r) => r.id === ev.claimedBy)?.name ?? ev.claimedBy}”${ev.blocked ? ` (${ev.blocked === 'gap' ? 'blocked by gap' : 'merge-into score missing'})` : ` → ${ev.side}`}` : ''
        data.push({ type: 'scatter', mode: held ? 'markers+text' : 'markers', text: held ? ['⚠'] : undefined, textposition: 'top center', textfont: { size: 12, color: SIGN_COLOR },
          x: [(x0 + x1) / 2], y: [seg.basis === 'raw' ? 0 : yr[1] - 0.35], marker: { size: 9, color, line: { color: merging ? rule0.hatch.color : '#111827', width: merging ? 2 : 1 }, symbol: merging || held ? 'diamond' : 'square' },
          customdata: [key], xaxis: 'x', yaxis: yref, showlegend: false, name: 'small segment',
          hovertemplate: `<b>${src} ${seg.chrom} · bins [${row.start}, ${row.end}) · ${row.end - row.start} bins</b><br>${trackDef.label}: ${fmt(cv.value, 3)} (${cv.side ?? '–'} flank)<br>${(r0 ? rule0Lines(r0) : [rule0Issues.length ? `Rule 0 inactive: ${rule0Issues.join('; ')}` : 'Rule 0: no score row']).join('<br>')}${userLines.length ? `<br>${userLines.join('<br>')}` : ''}${claimNote}${row.gap_left || row.gap_right ? `<br>genomic gap: ${row.gap_left ? 'left ' : ''}${row.gap_right ? 'right' : ''}` : ''}<extra></extra>` })
      })
```
   (`rule0` is in scope from item 2. The marker hover replaces both old hover carriers.)
5. **Figure dependencies and title.**
   - Remove `smallHere`, `effectMetric`, `metricDef` and `rulesOn` from the `trackFig` `useMemo` dependency array. Add `trackDef`, `rule0`, `rule0Issues` and `ig.smallSegs`.
   - Change the title tail from `` coloured by ${metricDef.label} (${seg.flank.replace('_', ' ')} flank) `` to `` coloured by ${trackDef.label} (smaller-|value| flank) · Rule 0 hatched ``.
6. **Candidate audit.** In the candidate-audit `hovertemplate`, replace `${r.rejection_reasons || 'no rejection reason'}` with `${derivedRejectionReason(r) ?? 'no rejection reason'}`.
7. **Track-panel heading and lead.**
   - In the `#cluster-tracks` panel, change the `<h3>` text to `Chromosome tracks — pseudobulk depth, accepted breakpoints and small segments`.
   - Change the `vis-heading` `<h4>` text to `1. Visualization and Rule 0`. Keep its span, but change its text to `— colouring is display only; Rule 0 below hatches the segments that the merge in section 2 starts from`.
   - Replace the `<p className="lead">` text with this sentence, and drop the pilot-chromosome sentence and the archived hatch descriptions: "Grey line: per-bin pseudobulk depth, broken at genomic gaps. Vertical lines: accepted breakpoints by recursion depth. Shaded spans: internal segments ≤ {seg.smallMax} bins, coloured by live flank scores (any chromosome) on the colour metric's smaller-|value| flank. Rule 0 hatches × the segments whose checked flank criteria are below τ / the SRD/√φ ref; an amber | marks flank SRD/√φ signs that disagree and a violet ― marks criteria that pick different flanks."
8. **Controls row 1.**
   - Replace the `colour by` select's options with `TRACK_METRICS`, bound to `seg.colorMetric`: `onChange={(e) => ig.setSeg({ colorMetric: e.target.value as any })}`.
   - Delete the `effect metric for flank choice` select and the `flank` select.
   - Rename the `SRD ref` label text to `SRD/√φ ref` and give its input `aria-label="Rule 0 SRD/√φ ref"`. Give the `τ` input `aria-label="Rule 0 tau"`.
   - Allow clearing: for both inputs use `value={Number.isFinite(seg.X) ? seg.X : ''}` and `onChange={(e) => ig.setSeg({ X: e.target.value === '' ? NaN : +e.target.value })}`.
9. **New Rule 0 controls row.** Add this row directly after row 1:
```tsx
        <div className="controls-row rule0-controls">
          <span className="small"><b>Rule 0</b> — hatch and merge segments where</span>
          <select aria-label="Rule 0 combine" value={seg.combine} onChange={(e) => ig.setSeg({ combine: e.target.value as 'all' | 'any' })}><option value="all">all</option><option value="any">any</option></select>
          <span className="small">of these hold:</span>
          {RULE0_CRITERIA.map((c) => <label key={c} className="chk-inline"><input type="checkbox" aria-label={`Rule 0 criterion ${CRITERION_LABEL[c]}`} checked={seg.criteria.includes(c)}
            onChange={(e) => ig.setSeg({ criteria: RULE0_CRITERIA.filter((x) => (x === c ? e.target.checked : seg.criteria.includes(x))) })} />{CRITERION_LABEL[c]} &lt; {c === 'srd_phi' ? 'SRD/√φ ref' : 'τ'}</label>)}
          <label className="chk-inline"><input type="checkbox" aria-label="Merge sign-disagreeing segments" checked={seg.mergeSignMismatch} onChange={(e) => ig.setSeg({ mergeSignMismatch: e.target.checked })} />merge would-merge segments whose flank signs disagree</label>
          <label className="chk-inline"><input type="checkbox" aria-label="Merge flank-conflict segments" checked={seg.mergeFlankConflict} onChange={(e) => ig.setSeg({ mergeFlankConflict: e.target.checked })} />merge segments whose criteria pick different flanks, toward the one furthest below its threshold</label>
          {rule0Issues.length > 0 && <span className="notice small">Rule 0 inactive: {rule0Issues.join('; ')}</span>}
        </div>
```
10. **Layers list.**
    - Remove the `['proposal', ...]` entry.
    - Replace `['consistency', 'flank-consistency hatch (×/‒/) + x/o/ label']` with `['signs', 'signs-disagree hatch (amber |) on segments that do not pass Rule 0']`.
    - Change `['reason', 'retained reason (effect ≥ τ / terminal)']` to `['reason', 'retained reason (archived prototype, where available)']`.
11. **Legend.**
    - In the `seg-legend` div, replace `metricDef` with `trackDef`, and delete the two lines `{seg.labels.proposal && ...}` and `{seg.labels.consistency && !seg.labels.proposal && ...}`.
    - Replace the existing per-user-rule legend line with `activeRules(ig.mergeRules).map(...)`, keeping its content. Keep the gap legend line, shown when `userRulesOn || rule0Issues.length === 0`.
    - Add the Rule 0 legend:
```tsx
          {rule0Issues.length === 0 && <span className="rule0-legend"><b className="rule-swatch" style={{ color: rule0.hatch.color }}>×</b>Rule 0 would merge: {r0Counts.wouldMerge}{r0Counts.signs ? ` · signs disagree: ${r0Counts.signs}` : ''}{r0Counts.flanks ? ` · flanks disagree: ${r0Counts.flanks}` : ''}{r0Counts.gap ? ` · blocked by gap: ${r0Counts.gap}` : ''} · <b style={{ color: SIGN_COLOR }}>|</b>/<b style={{ color: FLANK_COLOR }}>―</b> ⚠ held back: {r0Counts.heldBack} · <b style={{ color: SIGN_COLOR }}>|</b> signs disagree only: {r0Counts.signsOnly}</span>}
```
12. **Selected-segment card.**
    - Replace the `seg.selectedSegment && (() => { const s = ig.smallSegs.find(...) ... })()` block with a version that parses `const [ssrc, schrom, sstart] = seg.selectedSegment.split('|')`.
    - It finds `const row = ig.hatchScores[`${ssrc}|${schrom}`]?.rows.find((r) => r.start === +sstart)`, the matching eval `ruleEvals[ssrc]?.[index of row]`, and `const archived = ig.smallSegs.find((x) => x.key === seg.selectedSegment)`. It renders nothing if `row` is missing.
    - The card shows:
      - a title `{ssrc} {schrom} · bins [{row.start}, {row.end}) · {row.end - row.start} bins`, with the existing ✕ button
      - a table with one row per live metric (SRD, SRD/√φ, Δlog2FC median, mean, IQR-mean), columns left / right using `fmt(row[`${id}_left`])`, highlighting the row whose id equals `seg.colorMetric`
      - a line `φ = {fmt(row.phi)} · genomic gap: {left/right/none}`
      - the `rule0Lines(eval.rule0)` lines as a list, when available
      - archived extras when `archived` exists: the existing "Original median bootstrap fraction" row, the "eligible / merged / reason" row, and the SRD-screen bootstrap SE row (moved from the old card unchanged, under a sub-heading "Archived prototype (this segment)")
13. **Segment table.**
    - In the `Segments on …` table, replace the `small-segment metric (...)`, `consistency` and `proposal` columns with `colour value ({trackDef.label})` and `Rule 0`.
    - For each cluster segment `g`, find the live row with `start === g.s && end === g.e` for `g.group`, plus its eval. Show `fmt(colourValue(row, trackDef.id).value)` (or blank).
    - Show the Rule 0 outcome label: `would merge` / `held back` / `signs disagree` / blank for `none` or no row.
    - Clicking a row with a live row selects `${g.group}|${seg.chrom}|${g.s}`, setting `selectedSegment` without changing sources.
    - Change the `<summary>` count text to `{segTable.length} cluster segments`.
14. **Clean up.**
    - Delete `smallHere`. Its only remaining use was the summary count, which item 13 replaces.
    - Make sure `selectSeg` (used by the diagnostics click) still sets `selectedSegment: s.key`; archived keys use the same `source|chrom|start` format.
    - Delete any now-unused imports or constants so `tsc` reports no unused-symbol errors, if the project's tsconfig enforces them.

- [ ] **Step 4: Merge panel — `web/src/components/HatchMergeView.tsx`**

- Imports: replace `import { canRunMerge, describeCondition } from '../mergeRules'` with `import { canRunMerge, describeCondition, describeRule0, isRule0, rule0FromSeg } from '../mergeRules'`.
- After `const chrom = ig.seg.chrom`, add `const rule0 = rule0FromSeg(ig.seg)`. Replace `const runnable = canRunMerge(ig.mergeRules)` with `const runnable = canRunMerge(rule0, ig.mergeRules)`, and the request's `rules: ig.mergeRules` with `rules: [rule0, ...ig.mergeRules]`.
- `ruleSummary`: map `result.settings.rules.filter((r) => isRule0(r) || r.enabled)` to `isRule0(r) ? `Rule 0: ${describeRule0(r)}` : `${r.name}: ${r.conditions.map(describeCondition).join(' AND ')}`` and join with ` | `. Drop the old numbering.
- Not-runnable message: `Fix Rule 0's settings in section 1 (at least one criterion; τ and SRD/√φ ref positive; size limit a whole number) and complete every enabled rule (each needs at least one condition with valid values) to run a merge.`
- Lead paragraph:
  - Replace its first two sentences with: "Rule 0 (pinned first, built from section 1) merges the segments it hatches on the Chromosome tracks; your own rules below add to it with lower priority. A user rule is a set of conditions on a segment's left/right flank scores that must ALL hold."
  - Keep the rest, including the sentence about hatches showing what is flagged at the start.

- [ ] **Step 5: Rule cards — `web/src/components/MergeRulesEditor.tsx`**

- Import `describeRule0`, `rule0FromSeg` and `rule0Problems`.
- Compute `const rule0 = rule0FromSeg(ig.seg)` and `const r0Issues = rule0Problems(rule0)`. Change the counts call to `evaluateForSources(..., [rule0, ...rules], ig.allowGapCrossing)`.
- Before `{rules.map(...)}`, render the pinned card:
```tsx
    <div className="rule-card rule0-card">
      <div className="rule-head"><b>Rule 0 — shown on the tracks</b><span className="muted small">pinned first · edit it in section 1 above</span></div>
      <div className="small">{describeRule0(rule0)}</div>
      <div className="small" role="status" aria-label="Rule 0 status">
        {r0Issues.length ? <span className="notice small">{r0Issues.join('; ')}</span>
          : ig.hatchScoresError ? 'scores unavailable' : !scoresLoaded ? 'loading scores…'
          : `→ flags ${counts['rule0'] || 0} segment(s) on ${ig.seg.chrom} across the displayed rows`}
      </div>
    </div>
```
  (`scoresLoaded` already exists in this component; move its declaration above this block if needed.)
- Change the empty-list message to `No additional rules. Rule 0 above already defines what merges; add a rule to merge more segments.`

- [ ] **Step 6: Styles — `web/src/styles.css`**

Append:
```css
.rule0-card { border-color: #fca5a5; background: #fff7f7; }
.rule0-controls { border-top: 1px dashed #e2e8f0; padding-top: 4px; }
.rule0-legend { white-space: normal; }
```

- [ ] **Step 7: Leftover check and build**

Run: `grep -rn "seg\.metric\b\|seg\.effectMetric\|seg\.flank\b\|labels\.proposal\|labels\.consistency\|consistencyClass\|PROPOSAL_" web/src`
Expected: no output.

Run: `cd web && npm run build`
Expected: `tsc --noEmit` passes, and `vite build` prints `✓ built in …`.

Run all three TS helper tests and `python3 -m pytest -q tests/test_merge_rules.py tests/test_hatch_merge.py` again.
Expected: all pass.

- [ ] **Step 8: Commit**

```bash
git add web/src/hatchMerge.ts web/src/integrationStore.ts web/src/components/IntegrationSegmentationView.tsx web/src/components/HatchMergeView.tsx web/src/components/MergeRulesEditor.tsx web/src/styles.css
git commit -m "feat: live-score tracks with Rule 0 hatching, criteria checklist, pinned Rule 0 card

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Browser walkthrough

**Files:**
- Modify: `tests/browser_hatch_merge_walkthrough.cjs`
- Modify: `validation/hatch-merge/HANDOFF.md` (append a section)
- Generated (commit): `validation/hatch-merge/merge_rules_check_results.json`, `validation/hatch-merge/merge-rules.png`

**Interfaces:**
- Consumes: the accessible names from Task 4's Produces block, plus the existing user-rule names.

- [ ] **Step 1: Update the walkthrough**

Keep the script's structure: `BASE_URL`, `findOwner`, the existing user-rule steps, session export/restore, and the legacy `hatchLayers` migration check. Make these changes:

1. Add a helper `const r0Flags = async () => { const m = (await panel.getByLabel('Rule 0 status').innerText()).match(/flags (\d+)/); return m ? +m[1] : null }`.
2. **Rule 0 alone can run.** Replace the `runDisabledWithNoRules` check with `runEnabledWithRule0Only = !(await runBtn.isDisabled())`, recorded before any user rule is added, and assert it is `true`.
3. **Rule 0 responds to settings.** Before adding user rules:
   - set the `Rule 0 combine` select to `any` and the `Rule 0 tau` input to `5`
   - wait until `Rule 0 status` matches `/flags \d+/`
   - record `r0High`
   - set the tau input to `0.0001` and uncheck `Rule 0 criterion min |SRD/√φ| flank`
   - record `r0Low`
   - assert `r0High > r0Low` and `r0Low === 0`
   - restore: tau `0.2`, check the SRD/√φ criterion again, combine `all`
4. **Switches.** Set tau to `5`, combine `any`, and turn both merge switches off. Then:
   - Read the legend text from `.rule0-legend` and parse `would merge: (\d+)` and `held back: (\d+)` into `wmOff` and `hbOff`.
   - Turn on both `Merge sign-disagreeing segments` and `Merge flank-conflict segments`, then parse again into `wmOn` and `hbOn`.
   - Assert `hbOn === 0` and `wmOn === wmOff + hbOff`.
   - Turn both switches off again and restore tau `0.2` and combine `all`.
5. **Invalid Rule 0.** Clear the tau input.
   - Assert `Rule 0 status` contains `τ` and that `Run merge` is disabled.
   - Set tau back to `0.2`.
6. **Non-pilot chromosome.** Switch the tracks' chromosome select to `chr3` and wait for the plot title to include `chr3`. Then:
   - Assert `#cluster-tracks .js-plotly-plot` has at least one data trace named `small segment` (`plot.evaluate(el => el.data.filter(d => d.name === 'small segment').length)` > 0).
   - Switch back to `chr8`.
7. **Candidate audit.** Check the `candidate audit (rejected / blocked)` layer checkbox and wait 800 ms. Then:
   - Assert `plot.evaluate(el => el.data.some(d => [].concat(d.hovertemplate || []).some(h => String(h).includes('permutation p ='))))` is `true`.
   - Uncheck the layer.
8. **Merge history names Rule 0.**
   - Before adding user rules, set tau to `5` and combine to `any`, click `Run merge`, and wait for the history table.
   - Assert that the first history row's rule cell (td index 6) contains `Rule 0` when `historyRows > 0`. Also record `historyRows`.
   - Then restore tau `0.2` and combine `all`.
   - Keep the existing user-rule run later in the script. Its stale-notice check is unchanged.
9. **Assertions that change because of Rule 0:**
   - `shapesWithRule > shapesBaseline` still holds (adding a user rule with threshold 1000 adds hatches), so keep it.
   - The legacy-session check expects the three migrated user rules: `legacyNames` unchanged.
   - Add an old `seg` payload to the legacy session: before writing `legacy.json`, set `owner.seg = { ...owner.seg, metric: 'log2fc', flank: 'min_srd', labels: { ...owner.seg.labels, proposal: true, consistency: false } }`, with `colorMetric`, `criteria` and `signs` deleted.
   - After restore, assert the `colour by` select's value is `delta_mean` and the `signs-disagree hatch` layer checkbox is unchecked.
10. Record every new value in the results JSON (`runEnabledWithRule0Only, r0High, r0Low, wmOff, hbOff, wmOn, hbOn, chr3SmallSegments, auditDerivedReason, rule0HistoryRows`), alongside the existing ones.

- [ ] **Step 2: Build and start a fresh server on port 8767**

```bash
cd web && npm run build && cd ..
source .venv-scientific-explorer/bin/activate
nohup python3 run.py --port 8767 > /tmp/claude-1000/rule0-browser-8767.log 2>&1 &
echo $! > /tmp/claude-1000/rule0-browser-8767.pid
until curl -sf http://127.0.0.1:8767/api/integration/status | grep -q '"done"'; do sleep 3; done
```

- [ ] **Step 3: Run the walkthrough**

```bash
BASE_URL=http://127.0.0.1:8767 \
PLAYWRIGHT_MODULE=/home/mami/Desktop/EGA/job_search/career-ops/node_modules/playwright \
CHROMIUM_PATH=/home/mami/.cache/ms-playwright/chromium-1234/chrome-linux64/chrome \
node tests/browser_hatch_merge_walkthrough.cjs
```
Expected: the last line is the PASS line. If an assertion fails because the UI is wrong, fix the UI, rebuild, restart your 8767 server, and re-run. Fix the script only when its assumption cannot match a correct UI, and explain why in the report. Lazy plots (DeferredPlot) need scrolling into view. When finished, stop only your server: `kill "$(cat /tmp/claude-1000/rule0-browser-8767.pid)"`.

- [ ] **Step 4: Append to `validation/hatch-merge/HANDOFF.md`**

Add a section `## Rule 0 (2026-10-26)` with:
- the command run
- the key values from the results JSON (`r0High`/`r0Low`, `wmOff`/`hbOff`/`wmOn`/`hbOn`, `chr3SmallSegments`, `auditDerivedReason`, `rule0HistoryRows`)
- the PASS line
- the backend test count from `python3 -m pytest -q tests/test_*.py`
- a pointer to `validation/rule0/README.md` for the archived-vs-live comparison

- [ ] **Step 5: Commit**

```bash
git add tests/browser_hatch_merge_walkthrough.cjs validation/hatch-merge/HANDOFF.md validation/hatch-merge/merge_rules_check_results.json validation/hatch-merge/merge-rules.png
git commit -m "test: browser walkthrough covers Rule 0, live tracks and audit reasons

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Documentation and final verification

**Files:**
- Modify: `README.md`, `docs/USAGE.md`, `docs/SCIENTIFIC_WORKFLOW.md`, `docs/ARCHITECTURE.md`, `tests/README.md`, `web/README.md`

**Interfaces:**
- Consumes: final behaviour from Tasks 1–5 and `validation/rule0/README.md`.

- [ ] **Step 1: Update docs**

- **`docs/USAGE.md`**: rewrite the `### Visualization and merge rules` section to describe:
  - Section 1, **"1. Visualization and Rule 0"**:
    - colour by the live metrics (any chromosome), on each metric's smaller-|value| flank
    - the Rule 0 criteria checklist with all/any, τ and SRD/√φ ref
    - the two switches
    - the hatching cases and the legend
    - the ⚠ held-back marker
    - hovering a segment's marker for the per-criterion detail
    - the signs-disagree layer toggle
  - Section 2: the pinned Rule 0 card, user rules adding to it, Rule 0 merging first, and "Run merge" working with Rule 0 alone.
  - The candidate audit's derived p-value reasons, with the inferred cutoff.
  - That the SRD-screen proposal layer was removed.
  - That old sessions map their colour choice (table from spec §2).
  - Keep the existing sentences about recompute-after-each-merge, gaps, history, CSV/JSON, φ in history and the scientific framing.
- **`README.md`**:
  - In the `POST /api/integration/hatch-merge` paragraph, state that `rules[0]` must be Rule 0, `{id: "rule0", type: "flank_threshold", criteria, combine, tau, srdRef, maxBins, mergeSignMismatch, mergeFlankConflict, hatch}`, built from section 1, and that user rules follow it.
  - Update the "What the explorer can do" Integration · segmentation & CN bullet to mention live-score tracks with Rule 0 hatching.
- **`docs/SCIENTIFIC_WORKFLOW.md`** worked example 2: update the Diagnostic bullet to say the tracks' Rule 0 criteria (chosen flank below τ / SRD/√φ ref) hatch candidate segments, with amber/violet marks for sign and flank disagreement, and user rules adding further conditions. Keep the "starts from those flagged segments and re-evaluates every rule after each merge" wording.
- **`docs/ARCHITECTURE.md`**: add a frontend row for `trackMetrics.ts` (live colour metrics, session migration, candidate-audit reason derivation). Extend the `merge_rules.py` and `mergeRules.ts` rows with Rule 0.
- **`tests/README.md`**: add `track_helpers.ts` with its esbuild/node command, and `compare_archived_live.py` (needs a running server; `BASE_URL`).
- **`web/README.md`**: add `trackMetrics.ts` under Stores (`src/`).

- [ ] **Step 2: Check for stale wording and links**

Run: `grep -rn "flank-consistency hatch\|SRD-screen proposal: hatch\|effect metric for flank choice\|min |effect| flank\|pilot chromosomes)" README.md docs/USAGE.md docs/SCIENTIFIC_WORKFLOW.md docs/ARCHITECTURE.md docs/EXTENDING.md web/README.md tests/README.md`
Expected: no hits that describe current behaviour. Historical mentions must be framed as history. Then run the link check:

```bash
python3 - <<'EOF'
import re, os, sys
files=[]
for root, dirs, fs in os.walk('.'):
    if any(p in root for p in ('/.git','/node_modules','/.venv','/cache','docs/superpowers/plans','/.superpowers')): continue
    files += [os.path.join(root,f) for f in fs if f.endswith('.md')]
broken=[]
for f in files:
    for m in re.finditer(r'\[[^\]]*\]\(([^)]+)\)', open(f,encoding='utf-8',errors='ignore').read()):
        t=m.group(1).split('#')[0]
        if not t or t.startswith(('http://','https://','mailto:')): continue
        if not os.path.exists(os.path.normpath(os.path.join(os.path.dirname(f),t))): broken.append((f,t))
print('\n'.join(f'{a} -> {b}' for a,b in broken) or f'All links in {len(files)} Markdown files resolve.')
sys.exit(1 if broken else 0)
EOF
```
Expected: `All links in N Markdown files resolve.`

- [ ] **Step 3: Final verification**

```bash
python3 -m pytest -q tests/test_*.py
for t in merge_rules_helpers hatch_merge_helpers track_helpers; do web/node_modules/.bin/esbuild tests/$t.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/$t.mjs --log-level=warning && node tests/_bundle/$t.mjs | tail -1; done
cd web && npm run build && cd ..
git diff --stat 59c8303..HEAD -- server/hatch_merge.py server/science.py server/flank_view.py server/pyepi_adapter.py reference/
git status --short
```
Expected:
- pytest passes in full
- each helper prints `all … checks passed`
- the build succeeds
- the `git diff --stat` prints nothing (scoring and reference code untouched)
- `git status --short` shows only the two untracked plan files

- [ ] **Step 4: Commit**

```bash
git add README.md docs/USAGE.md docs/SCIENTIFIC_WORKFLOW.md docs/ARCHITECTURE.md tests/README.md web/README.md
git commit -m "docs: describe Rule 0, live-score tracks and derived audit reasons

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
