# Rule 0 revision 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the always-true genomic-gap flag. Always merge toward the criterion furthest below its threshold. Add a last-phase merge for signs-disagree segments that fail Rule 0. Allow editable Rule 0 copies as user rules. Always draw breakpoints by depth. Show per-side pass/fail in hovers.

**Architecture:**
- Python `server/merge_rules.py` stays the authority.
- `web/src/mergeRules.ts` ports it for live hatching, and the shared parity fixture pins the two together.
- The UI changes live in `IntegrationSegmentationView.tsx` (tracks and section 1) and `MergeRulesEditor.tsx` (section 2).

**Tech Stack:** FastAPI/Python 3 + numpy, pytest; React 18 + TypeScript + Vite + Zustand + Plotly; esbuild+node for TS helper tests; Playwright browser walkthrough.

**Spec:** `docs/superpowers/specs/2026-10-26-rule0-revision-2-design.md` (builds on `docs/superpowers/specs/2026-10-26-rule0-visual-merge-design.md`).

## Global Constraints

- Gap test: `is_gap` ⇔ `var_start[right] > var_end[left] + 1`.
- Rule 0 object: `{"id":"rule0","type":"flank_threshold","name","criteria","combine","tau","srdRef","maxBins","mergeSignMismatch":bool,"mergeSignsOnly":bool,"hatch":{"pattern","color"}}`.
  - `mergeFlankConflict` is ignored if present and never required.
- Rule 0 copy (user rule): `{"id": <not "rule0">, "type":"flank_threshold", "name", "enabled":bool, "criteria","combine","tau","srdRef","maxBins","mergeSignMismatch":bool,"hatch"}`.
  - `mergeSignsOnly` is ignored on copies.
- Direction for a passing segment: the flank of the passing criterion with the smallest fraction (ties → left). `held_back` ⊆ `["signs"]`.
- Outcomes: `none | would_merge | held_back | signs_only | signs_merge`.
- Phases: Rule 0 `would_merge` and every user-rule claim are phase 1. Rule 0 `signs_merge` claimed by no user rule is phase 2.
  - Candidate order: `(phase, rule index, score, start)`.
  - Steps carry `phase` (`null` on step 0).
  - Phase-2 `rule_name` = `f"{name} · signs disagree, fails Rule 0 (last phase)"`.
- Colours:
  - Rule 0 `#b91c1c`; signs-merge black `#111827`; amber `#d97706`.
  - The violet flank-conflict colour is removed.
- UI strings, verbatim:
  - `"also merge signs-disagree segments that fail Rule 0 (last, black ×)"`
  - menu entry `"Rule 0 copy (current section-1 settings)"`
  - copy default name `"Rule 0 copy"`
- All subagents run on sonnet or haiku.
- Never kill the user's server on port 8766; use a separate port (e.g. 8799) for browser runs.
- Stage files by name only.
- Leave `docs/superpowers/plans/2026-09-21-typesafe-jev-judgment-layers.md` and `docs/superpowers/plans/2026-09-30-portfolio-repository-polish.md` untracked and untouched.
- Commit trailer: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

- **Saved sessions** with `mergeFlankConflict`, `labels.breakpoints` or `labels.depthColor` must restore without error and without stray keys. Covered in Task 3 `track_helpers.ts`.
- **A Rule 0 copy with a cleared τ** (NaN) must not hatch, must block "Run merge" while enabled, and must not block it while disabled. Covered in Task 3.
- **A phase-2 merge that creates a phase-1 candidate** must merge that candidate next. Covered in Task 2 `test_phase1_candidate_created_by_phase2_goes_next`.
- **A segment that is `signs_merge` for Rule 0 and also matches a user rule** is claimed by the user rule in phase 1. Covered in Task 2 and in the parity fixture.
- **The real-data gap count** must be small: about 549 of 25,336 boundaries genome-wide, not every one. Covered in Task 1 (recorded) and Task 5 (browser).

---

### Task 1: Genomic-gap fix

**Files:**
- Modify: `server/hatch_merge.py:118-121`
- Test: `tests/test_hatch_merge.py` (near `test_is_gap_detects_bp_discontinuity`, line ~70)
- Create: `validation/rule0/gap_fix.md`

**Interfaces:** Produces the corrected `is_gap(var_start, var_end, left_end_bin, right_start_bin) -> bool`. No signature change.

- [ ] **Step 1: Write the failing test** (append after the existing `is_gap` test)

```python
def test_is_gap_one_based_inclusive_coordinates():
    # pyEpiAneufinder bins: [800001, 900000], [900001, 1000000], then a 100 kb hole
    var_start = np.array([800001, 900001, 1100001], dtype=float)
    var_end = np.array([900000, 1000000, 1200000], dtype=float)
    assert hm.is_gap(var_start, var_end, 0, 1) is False
    assert hm.is_gap(var_start, var_end, 1, 2) is True
```

Check that the file already imports `numpy as np`, and add the import if not.

- [ ] **Step 2: Run it and check that it fails**

Run: `python3 -m pytest tests/test_hatch_merge.py -k is_gap -v`
Expected: `test_is_gap_one_based_inclusive_coordinates` FAILS at the first assert. The old test passes.

- [ ] **Step 3: Implement**

```python
def is_gap(var_start, var_end, left_end_bin, right_start_bin):
    """True if the two retained bins straddling a boundary are not genomically
    contiguous (a filtered/omitted interval sits between them). Works for
    1-based inclusive bins (next start = end + 1) and half-open bins (next
    start = end)."""
    return bool(float(var_start[right_start_bin]) > float(var_end[left_end_bin]) + 1)
```

- [ ] **Step 4: Run the full suite**

Run: `python3 -m pytest -q`
Expected: all pass (148 + 1). If a test asserted a gap flag on contiguous half-open data, it still passes, because `start == end` gives no gap.

- [ ] **Step 5: Record real-data counts**

Run from the workspace root (`..` of the repo):

```bash
python3 - <<'EOF'
import h5py, numpy as np
h = h5py.File('pyepi_results/count_matrix.h5ad', 'r')
s = h['var/start'][...].astype(np.int64); e = h['var/end'][...].astype(np.int64)
seq = h['var/seqnames'][...] if 'seqnames' in h['var'] else None
d = s[1:] - e[:-1]
print('boundaries', len(d), 'old rule (start != end)', int((s[1:] != e[:-1]).sum()), 'new rule (start > end+1)', int((d > 1).sum()))
EOF
```

Write `validation/rule0/gap_fix.md`. It records the command and the three numbers, which should be about 25336 / 25336 / 549. Note that the count includes chromosome transitions, which are never internal segment boundaries. Add one sentence: "Before this fix every internal boundary was flagged as a genomic gap; revision 1's 'blocked by gap: 13' on chr8 came from this bug."

- [ ] **Step 6: Commit**

```bash
git add server/hatch_merge.py tests/test_hatch_merge.py validation/rule0/gap_fix.md
git commit -m "fix: genomic-gap test for 1-based inclusive bin coordinates

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Python engine — fraction direction, signs phase, Rule 0 copies

**Files:**
- Modify: `server/merge_rules.py` (`_validate_rule0`, `validate_rules`, `evaluate_rule0`, `evaluate_rules`, `run_rule_merge`, module docstring)
- Modify: `tests/test_merge_rules.py`
- Modify: `tests/fixtures/merge_rules_parity.json`

**Interfaces:**
- Produces:
  - `evaluate_rule0(row, r0, allow_gap_crossing=False) -> dict` with keys `start, end, evaluated, criteria, passed, sign_mismatch, flank_conflict, held_back, outcome, side, score, blocked`.
    - `outcome` adds `'signs_merge'`.
    - `blocked` may be `'gap'` or `'nonfinite'`.
  - Each `evaluate_rules` row adds a `phase` key (`1`, `2` or `None`).
  - Each `run_rule_merge` step adds `phase`.
  - User rules may have `type: 'flank_threshold'`.
- Task 3 ports these exactly.

- [ ] **Step 1: Update test helpers and the parity fixture's Rule 0 objects**

In `tests/test_merge_rules.py`, change `_rule0` so the dict has `'mergeSignMismatch': False, 'mergeSignsOnly': False` (drop `'mergeFlankConflict'`). Add:

```python
def _copy(rid='c1', enabled=True, **kw):
    r = _rule0(**kw)
    r.pop('mergeSignsOnly')
    r.update(id=rid, name=rid, enabled=enabled)
    return r
```

In `tests/fixtures/merge_rules_parity.json`, replace every `"mergeFlankConflict": false` / `true` with `"mergeSignsOnly": false`. Case `rule0_conflict_allowed_size_limit_user_rule_after` had `mergeFlankConflict: true`. That is now the default behaviour, so its expectations stay the same apart from `heldBack` (see Step 5).

- [ ] **Step 2: Rewrite the flank-conflict tests and add new failing tests**

Replace the existing tests that use `mergeFlankConflict` (around lines 440–452) with the tests below. Keep the existing sign-mismatch tests.

```python
def test_flank_conflict_always_merges_toward_smallest_fraction():
    both = ['delta_mean', 'srd_phi']
    # delta_mean picks left (0.1/0.5 = 0.2); srd_phi picks right (0.33/3.3 = 0.1)
    row = _row(delta_mean_left=0.1, delta_mean_right=0.4, srd_phi_left=3.0, srd_phi_right=0.33)
    ev = mr.evaluate_rule0(row, _rule0(criteria=both))
    assert ev['flank_conflict'] is True
    assert ev['held_back'] == []
    assert ev['outcome'] == 'would_merge'
    assert ev['side'] == 'right'


def test_merge_flank_conflict_key_is_ignored():
    row = _row(delta_mean_left=0.1, delta_mean_right=0.4)
    rules = mr.validate_rules([_rule0(mergeFlankConflict=False)])
    assert mr.evaluate_rule0(row, rules[0])['outcome'] == 'would_merge'


def test_signs_only_display_when_option_off():
    row = _row(delta_mean_left=0.9, delta_mean_right=-0.8, srd_phi_left=-2.0, srd_phi_right=2.0)
    ev = mr.evaluate_rule0(row, _rule0())
    assert (ev['outcome'], ev['side'], ev['score']) == ('signs_only', None, None)


def test_signs_merge_direction_and_score_from_failing_criteria():
    both = ['delta_mean', 'srd_phi']
    # fails: delta_mean min 0.8 ≥ 0.5; srd_phi min 2.0 < 3.3 passes, but combine=all needs both
    row = _row(delta_mean_left=0.9, delta_mean_right=-0.8, srd_phi_left=-2.0, srd_phi_right=2.5)
    ev = mr.evaluate_rule0(row, _rule0(criteria=both, mergeSignsOnly=True))
    assert ev['passed'] is False
    assert ev['outcome'] == 'signs_merge'
    # fractions: delta_mean right 0.8/0.5 = 1.6; srd_phi left 2.0/3.3 ≈ 0.606 → left
    assert ev['side'] == 'left'
    assert ev['score'] == pytest.approx(2.0 / 3.3)
    assert ev['blocked'] is None


def test_signs_merge_gap_and_nonfinite():
    row = _row(gap_left=True, delta_mean_left=0.9, delta_mean_right=-1.0, srd_phi_left=-2.0, srd_phi_right=2.0)
    assert mr.evaluate_rule0(row, _rule0(mergeSignsOnly=True))['blocked'] == 'gap'
    assert mr.evaluate_rule0(row, _rule0(mergeSignsOnly=True), allow_gap_crossing=True)['blocked'] is None
    row2 = _row(srd_phi_left=-2.0, srd_phi_right=2.0)  # delta_mean missing on both sides
    ev = mr.evaluate_rule0(row2, _rule0(mergeSignsOnly=True))
    assert (ev['outcome'], ev['side'], ev['blocked']) == ('signs_merge', None, 'nonfinite')


def test_held_back_is_not_signs_merge():
    row = _row(delta_mean_left=0.1, delta_mean_right=0.4, srd_phi_left=-1.0, srd_phi_right=1.0)
    ev = mr.evaluate_rule0(row, _rule0(mergeSignsOnly=True))
    assert (ev['outcome'], ev['held_back']) == ('held_back', ['signs'])


def test_evaluate_rules_phases():
    r0 = _rule0(mergeSignsOnly=True)
    merge_row = _row(start=10, end=20, delta_mean_left=0.1, delta_mean_right=0.4)
    signs_row = _row(start=20, end=30, delta_mean_left=0.9, delta_mean_right=-0.8, srd_phi_left=-2.0, srd_phi_right=2.0)
    evs = mr.evaluate_rules([merge_row, signs_row], mr.validate_rules([r0]))
    assert [(e['claimed_by'], e['phase'], e['matched']) for e in evs] == [('rule0', 1, ['rule0']), ('rule0', 2, ['rule0'])]


def test_user_rule_claims_signs_merge_segment_in_phase_1():
    r0 = _rule0(mergeSignsOnly=True)
    user = _rule('u', [{'kind': 'sign_mismatch', 'metric': 'srd_phi', 'estimator': 'mean'}])
    signs_row = _row(start=20, end=30, delta_mean_left=0.9, delta_mean_right=-0.8, srd_phi_left=-2.0, srd_phi_right=2.5)
    ev = mr.evaluate_rules([signs_row], mr.validate_rules([r0, user]))[0]
    assert (ev['claimed_by'], ev['phase'], ev['matched']) == ('u', 1, ['u'])
    assert ev['rule0']['outcome'] == 'signs_merge'


def test_rule0_copy_matches_on_would_merge_and_claims_with_its_side():
    copy = _copy(criteria=['srd_phi'], srdRef=1.0)
    row = _row(srd_phi_left=0.9, srd_phi_right=0.2, delta_mean_left=5.0, delta_mean_right=5.0)
    ev = mr.evaluate_rules([row], mr.validate_rules([_rule0(), copy]))[0]
    assert (ev['matched'], ev['claimed_by'], ev['side'], ev['phase']) == (['c1'], 'c1', 'right', 1)
    assert ev['score'] == pytest.approx(0.2)


def test_disabled_copy_is_ignored_even_if_invalid():
    rules = mr.validate_rules([_rule0(), _copy(enabled=False, criteria=[])])
    assert [r['id'] for r in rules] == ['rule0']


def test_copy_ignores_merge_signs_only():
    copy = _copy(criteria=['delta_mean'])
    copy['mergeSignsOnly'] = True
    row = _row(delta_mean_left=0.9, delta_mean_right=-0.8, srd_phi_left=-2.0, srd_phi_right=2.0)
    ev = mr.evaluate_rules([row], mr.validate_rules([_rule0(), copy]))[0]
    assert ev['matched'] == [] and ev['claimed_by'] is None


def _step_profile(levels, width=5):
    return np.concatenate([np.full(width, float(v)) for v in levels])


def test_run_merges_phase1_before_lower_score_phase2():
    # segments of 5 bins each; boundaries every 5. Levels chosen so the middle segment (A)
    # passes Rule 0 and a later segment (B) only has disagreeing signs.
    profile = _step_profile([10, 10.2, 10, 30, 10, 30])
    n = len(profile)
    vs = np.arange(n, dtype=float) * 100
    ve = vs + 100
    out = mr.run_rule_merge(profile, [5, 10, 15, 20, 25], n, vs, ve, 0,
                            [_rule0(tau=0.2, maxBins=5, mergeSignMismatch=False, mergeSignsOnly=True)])
    phases = [s['phase'] for s in out['steps'][1:]]
    assert phases, 'expected at least one merge'
    assert phases == sorted(phases), f'phase 2 ran before phase 1 finished: {phases}'
    assert out['steps'][0]['phase'] is None
    if 2 in phases:
        k = phases.index(2) + 1
        assert out['steps'][k]['rule_name'] == 'Rule 0 · signs disagree, fails Rule 0 (last phase)'
```

The exact levels in `test_run_merges_phase1_before_lower_score_phase2` may not produce both phases. The implementer must print the run's `(rule_id, phase, segment)` steps and adjust the levels until the run contains at least one phase-1 and at least one phase-2 step. Then replace the soft `if 2 in phases` with hard asserts: `assert 1 in phases and 2 in phases`.

Use the same approach (adjust levels until it holds, then assert hard) for:

```python
def test_phase1_candidate_created_by_phase2_goes_next():
    # Construct (by adjusting levels) a run in which some phase-2 step is immediately
    # followed by a phase-1 step; assert that such a pair exists and that every step
    # after a phase-1-eligible state is phase 1 (i.e. the phase list is not simply sorted).
    ...
```

If no profile yields a phase-2 → phase-1 transition after reasonable search (≤ 30 minutes), replace this test with a unit test on the candidate key. Add a small helper `_candidate_key(order, ev)` to `merge_rules.py`, used by `run_rule_merge`, returning `(ev['phase'], order[ev['claimed_by']], ev['score'], ev['start'])`. Assert that a phase-1 ev sorts before a phase-2 ev with a lower score. Report which variant you used.

Add validation cases to the existing parametrized invalid list (near line 495). Replace `[_rule0(mergeFlankConflict=None)]` with:

```python
    [_rule0(mergeSignsOnly='yes')],
    [_rule0(), _copy(rid='rule0')],
    [_rule0(), _copy(tau=0)],
    [_rule0(), {**_copy(), 'type': 'bogus'}],
```

- [ ] **Step 3: Run the tests and check that they fail**

Run: `python3 -m pytest tests/test_merge_rules.py -q`
Expected: the new tests FAIL (KeyError `phase`, outcome mismatches, copy rejected by validation).

- [ ] **Step 4: Implement in `server/merge_rules.py`**

Replace `_validate_rule0` with a shared validator:

```python
def _validate_flank_rule(r, where, signs_only):
    criteria = r.get('criteria')
    if not isinstance(criteria, list) or not criteria:
        raise ValueError(f'{where}: check at least one criterion')
    if not all(isinstance(c, str) for c in criteria) or len(set(criteria)) != len(criteria) or any(c not in RULE0_CRITERIA for c in criteria):
        raise ValueError(f'{where}: criteria must be distinct values from {RULE0_CRITERIA}')
    if r.get('combine') not in RULE0_COMBINE:
        raise ValueError(f'{where}: combine must be one of {RULE0_COMBINE}')
    for key in ('tau', 'srdRef'):
        if not _positive(r.get(key)):
            raise ValueError(f'{where}: {key} must be a positive finite number')
    n = r.get('maxBins')
    if not (isinstance(n, int) and not isinstance(n, bool) and n >= 1):
        raise ValueError(f'{where}: maxBins must be an integer >= 1')
    keys = ('mergeSignMismatch', 'mergeSignsOnly') if signs_only else ('mergeSignMismatch',)
    for key in keys:
        if not isinstance(r.get(key), bool):
            raise ValueError(f'{where}: {key} must be true or false')
    _check_hatch(r, where)
```

In `validate_rules`:
- Allow `type` in `('conditions', RULE0_TYPE)` for `rules[1:]`. Keep rejecting any other id equal to `RULE0_ID`. The existing unique-id check already covers that, but add an explicit message: `f'rule {r["id"]!r}: only rules[0] may use id "rule0"'`. Raise `f'rule {r["id"]!r}: unknown rule type {t!r}'` for other types.
- Call `_validate_flank_rule(rules[0], 'Rule 0', signs_only=True)`.
- For each enabled user rule of type `flank_threshold`, call `_validate_flank_rule(r, where, signs_only=False)`. Otherwise run the existing conditions validation.

Replace `evaluate_rule0` from `out['flank_conflict'] = ...` onward, and its not-passed branch:

```python
    if not out['passed']:
        if out['sign_mismatch']:
            out['outcome'] = 'signs_only'
            if r0.get('mergeSignsOnly') is True:
                out['outcome'] = 'signs_merge'
                valued = [c for c in out['criteria'] if c['fraction'] is not None]
                if not valued:
                    out['blocked'] = 'nonfinite'
                    return out
                best = min(valued, key=lambda c: (c['fraction'], 0 if c['side'] == 'left' else 1))
                out['side'], out['score'] = best['side'], best['fraction']
                if not allow_gap_crossing and row[f'gap_{out["side"]}']:
                    out['blocked'] = 'gap'
        return out
    out['flank_conflict'] = len({c['side'] for c in passing}) > 1
    fractions = [c['fraction'] for c in passing]
    out['score'] = max(fractions) if r0['combine'] == 'all' else min(fractions)
    best = min(passing, key=lambda c: (c['fraction'], 0 if c['side'] == 'left' else 1))
    out['side'] = best['side']
    if out['sign_mismatch'] and not r0['mergeSignMismatch']:
        out['held_back'].append('signs')
        out['outcome'] = 'held_back'
        return out
    out['outcome'] = 'would_merge'
    if not allow_gap_crossing and row[f'gap_{out["side"]}']:
        out['blocked'] = 'gap'
    return out
```

When `evaluate_rules` evaluates a copy, it calls `evaluate_rule0(row, {**rule, 'mergeSignsOnly': False}, allow_gap_crossing)`, so copies never produce `signs_merge`. Update the docstring to cite both specs.

Replace `evaluate_rules`:

```python
def evaluate_rules(rows, rules, allow_gap_crossing=False):
    """Per score row: which rules match, which claims it (and in which phase),
    and the merge direction/score or why it cannot merge. `rules` must be
    validated; a leading Rule 0 is evaluated first and claims before user rules.
    Phase 1: Rule 0 would-merge and every user-rule claim. Phase 2: Rule 0
    signs_merge (fails the test, signs disagree) when no user rule matches."""
    r0 = rules[0] if rules and rules[0].get('id') == RULE0_ID else None
    user = rules[1:] if r0 is not None else rules
    out = []
    for row in rows:
        r0ev = evaluate_rule0(row, r0, allow_gap_crossing) if r0 is not None else None
        user_evals = {}
        matched = [RULE0_ID] if r0ev is not None and r0ev['outcome'] == 'would_merge' else []
        for r in user:
            if r.get('type') == RULE0_TYPE:
                e = evaluate_rule0(row, {**r, 'mergeSignsOnly': False}, allow_gap_crossing)
                if e['outcome'] == 'would_merge':
                    user_evals[r['id']] = e
                    matched.append(r['id'])
            elif all(condition_holds(c, row) for c in r['conditions']):
                matched.append(r['id'])
        ev = {'start': row['start'], 'end': row['end'], 'matched': matched, 'claimed_by': None, 'phase': None,
              'side': None, 'score': None, 'blocked': None, 'rule0': r0ev}
        if matched and matched[0] == RULE0_ID:
            ev.update(claimed_by=RULE0_ID, phase=1, side=r0ev['side'], score=r0ev['score'], blocked=r0ev['blocked'])
        elif matched:
            rule = next(r for r in user if r['id'] == matched[0])
            ev.update(claimed_by=rule['id'], phase=1)
            if rule['id'] in user_evals:
                e = user_evals[rule['id']]
                ev.update(side=e['side'], score=e['score'], blocked=e['blocked'])
            else:
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
        elif r0ev is not None and r0ev['outcome'] == 'signs_merge':
            ev.update(matched=[RULE0_ID], claimed_by=RULE0_ID, phase=2, side=r0ev['side'],
                      score=r0ev['score'], blocked=r0ev['blocked'])
        out.append(ev)
    return out
```

In `run_rule_merge`:
- Change the candidate key to `(ev['phase'], order[ev['claimed_by']], ev['score'], ev['start'])`, via `_candidate_key` if you added it.
- Add `'phase': None` to step 0.
- Add `'phase': ev['phase']` to each step.
- Set `rule_name` to `names[...]` when phase is 1, else `f"{names[RULE0_ID]} · signs disagree, fails Rule 0 (last phase)"`.
- Update the docstring.
- Update the module docstring: mention Rule 0 copies and the last signs phase.

- [ ] **Step 5: Extend the parity fixture**

In every case's `expected` rows, add `"phase"` (`1` where `claimedBy` is non-null, else `null`). Where a case's rule0 expectation had `"heldBack": ["flanks"]` or `["signs","flanks"]`, recompute: flanks no longer hold back. Fix `outcome`, `side`, `score`, `claimedBy`, `matched` and `blocked` to match. Python is the authority: compute each expectation by reasoning from the rules, then confirm Python agrees. If they disagree, find which is wrong. Do not just copy Python's output.

Add two cases over the existing `rows`:
1. `rule0_signs_merge_last_phase`: Rule 0 with `criteria: ["delta_mean"], combine: "all", tau: 0.2, srdRef: 3.3, maxBins: 100, mergeSignMismatch: false, mergeSignsOnly: true`, and `allowGapCrossing: false`. Row `[30,40)` (srd_phi −2 / 2, delta_mean ∓0.5) fails the test and has signs disagreeing. Expect `signs_merge`, phase 2, side `left` (tie → left), score `2.5`.
2. `rule0_copy_and_user_rule`: inert Rule 0, then a copy `{id:"c", type:"flank_threshold", enabled:true, criteria:["srd_phi"], combine:"all", tau:0.2, srdRef:1.0, maxBins:100, mergeSignMismatch:true, hatch:{pattern:"x",color:"#1d4ed8"}}`, then the existing user rule `a` (weak srd_phi < 1.0). Row `[10,20)` has srd_phi min 0.5 on the left, so both match, and `c` claims it with side left, score 0.5, phase 1, and matched `["c","a"]`. Work out the other rows the same way.

Extend `test_python_evaluator_matches_parity_fixture` to compare `g['phase'] == exp['phase']`. Drop `flank_conflict`/`flankConflict` from the compared rule0 keys only if no case asserts it; otherwise keep it.

- [ ] **Step 6: Run all Python tests**

Run: `python3 -m pytest -q`
Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add server/merge_rules.py tests/test_merge_rules.py tests/fixtures/merge_rules_parity.json
git commit -m "feat: Rule 0 merges toward smallest fraction, last-phase signs merge, Rule 0 copies

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: TypeScript engine port, store and migration

**Files:**
- Modify: `web/src/mergeRules.ts`
- Modify: `web/src/hatchMerge.ts` (`MergeStep` gains `phase: 1 | 2 | null`; rule types)
- Modify: `web/src/integrationStore.ts`
- Modify: `web/src/trackMetrics.ts` (`migrateSegSettings`)
- Modify: `tests/merge_rules_helpers.ts`, `tests/track_helpers.ts`

**Interfaces:**
- Consumes: the Task 2 semantics and fixture (`phase`, `signs_merge`, copies).
- Produces, for Task 4:
  - `interface FlankRule { id: string; type: 'flank_threshold'; name: string; enabled?: boolean; criteria: Rule0Criterion[]; combine: 'all' | 'any'; tau: number; srdRef: number; maxBins: number; mergeSignMismatch: boolean; mergeSignsOnly?: boolean; hatch: { pattern: HatchPattern; color: string } }`
  - `type Rule0 = FlankRule & { id: 'rule0'; mergeSignsOnly: boolean }`
  - `type UserRule = MergeRule | FlankRule`; `type AnyRule = Rule0 | UserRule`
  - `isRule0(r) = r.id === 'rule0'`; `isFlankRule(r): r is FlankRule = r.type === 'flank_threshold'`
  - `Rule0Eval.outcome` adds `'signs_merge'`; `Rule0Eval.blocked: null | 'gap' | 'nonfinite'`; `heldBack: 'signs'[]`
  - `SegmentEval.phase: 1 | 2 | null`
  - `flankRuleProblems(r: FlankRule): string[]` (renamed from `rule0Problems`; keep `rule0Problems` as an alias)
  - `userRuleProblems(r: UserRule): string[]` (dispatches to `ruleProblems` or `flankRuleProblems`)
  - `activeRules(rules: UserRule[]): UserRule[]`
  - `canRunMerge(rule0: Rule0, rules: UserRule[]): boolean`
  - `rule0FromSeg(s: { criteria; combine; tau; srdRef; smallMax; mergeSignMismatch; mergeSignsOnly }): Rule0`
  - `makeRule0Copy(s: same shape as rule0FromSeg input, existing: UserRule[]): FlankRule` (name `'Rule 0 copy'`, pattern `'x'`, colour from `RULE_COLORS[existing.length % RULE_COLORS.length]`, `enabled: true`, `maxBins: s.smallMax`, no `mergeSignsOnly`)
  - `criterionSideLine(row: HatchScoreRow, c: CriterionEval): string`. It returns, e.g., `min |Δlog2FC mean| flank < 0.2: L −0.123 ✓ · R 0.401 ✗ → left (fraction 0.62) pass`, with signed values to 3 decimals, `N/A` for a missing side, `✓` if |v| < threshold, and `→ – fail` when no side has a value. The minus sign is U+2212.
  - `rule0Lines(e: Rule0Eval, row: HatchScoreRow): string[]` (now takes the row; per-criterion lines use `criterionSideLine`)
  - `rule0CaseCounts(evals)` returns `{ wouldMerge, gap, signs, heldBack, signsOnly, signsMerge, signsMergeBlocked }` (no `flanks`)
  - `SIGNS_MERGE_COLOR = '#111827'`; remove `FLANK_COLOR`.
  - The `SegSettings` field `mergeFlankConflict` is replaced by `mergeSignsOnly: boolean` (default false). `LabelLayers` drops `breakpoints` and `depthColor`.
  - `IntegrationState.mergeRules: UserRule[]`.

- [ ] **Step 1: Update the TS tests first**

In `tests/merge_rules_helpers.ts`:
- In the parity loop, also compare `g.phase === exp.phase`.
- Replace `mergeFlankConflict: false` in the `rule0FromSeg` call with `mergeSignsOnly: false`.

Add:

```ts
const copy = makeRule0Copy({ criteria: ['srd_phi'], combine: 'all', tau: 0.2, srdRef: 1.0, smallMax: 100, mergeSignMismatch: true, mergeSignsOnly: true }, [])
eq([copy.type, copy.name, copy.enabled, copy.maxBins, copy.hatch.pattern, 'mergeSignsOnly' in copy, copy.id !== 'rule0'],
  ['flank_threshold', 'Rule 0 copy', true, 100, 'x', false, true], 'makeRule0Copy shape')
eq(canRunMerge(r0, [{ ...copy, tau: NaN }]), false, 'enabled copy with cleared τ blocks the run')
eq(canRunMerge(r0, [{ ...copy, tau: NaN, enabled: false }]), true, 'disabled invalid copy does not block')
eq(evaluateRules(fx.rows, [r0, { ...copy, srdRef: NaN }], false).every((e) => !e.matched.includes(copy.id)), true, 'invalid copy never hatches')
const row0 = fx.rows[0]
eq(criterionSideLine(row0, { criterion: 'delta_mean', side: 'left', value: 0.07, threshold: 0.2, fraction: 0.35, passed: true }),
  'min |Δlog2FC mean| flank < 0.2: L 0.070 ✓ · R −0.900 ✗ → left (fraction 0.35) pass', 'per-side tooltip line')
```

The existing `describeRule0` assertion, if any, must change to expect the new summary text: `… · signs disagree: held back · fails Rule 0 with signs disagreeing: not merged`. With `mergeSignsOnly` true, the last part reads `merged last`.

In `tests/track_helpers.ts`, add a migration test. `migrateSegSettings({ mergeFlankConflict: true, labels: { breakpoints: false, depthColor: false, signs: true } })` must return an object without `mergeFlankConflict`, with `mergeSignsOnly` absent or false, and with `labels` lacking `breakpoints`/`depthColor` but keeping `signs: true`.

- [ ] **Step 2: Bundle and run them; check that they fail**

```bash
web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs && node tests/_bundle/merge_rules_helpers.mjs
web/node_modules/.bin/esbuild tests/track_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/track_helpers.mjs && node tests/_bundle/track_helpers.mjs
```

Expected: bundle errors (missing exports) or FAIL lines.

- [ ] **Step 3: Implement `mergeRules.ts`**

Port Task 2's Python 1:1.

`evaluateRule0(row, r0: FlankRule, allowGapCrossing)`:
- After computing criteria, `passed` and `signMismatch`:
  - If the segment does not pass: on a sign mismatch, set `outcome = 'signs_only'`. If `r0.mergeSignsOnly === true`, set `outcome = 'signs_merge'` and pick the side and score from the criteria with a non-null fraction (smallest fraction, left on ties). With none, set `blocked = 'nonfinite'`. Otherwise apply the gap check.
  - If it passes: `flankConflict` = more than one distinct side among passing criteria; `side` = side of the passing criterion with the smallest fraction (left on ties); `score` as before; on a sign mismatch with `!r0.mergeSignMismatch`, set `heldBack = ['signs']` and `outcome = 'held_back'`; otherwise `would_merge` plus the gap check.

`evaluateRules(rows, rules: AnyRule[], allowGapCrossing)` mirrors Python's `evaluate_rules`:
- Rule 0 is used only if `rules[0]` is Rule 0 and `flankRuleProblems(rules[0])` is empty.
- User rules = `activeRules(rest)`. A copy is evaluated with `{ ...r, mergeSignsOnly: false }`.
- Phase 2 fallback for `signs_merge`.
- `ev.phase` as in Python.

Other changes:
- `activeRules` keeps rules with `enabled === true && userRuleProblems(r).length === 0`.
- `canRunMerge` uses `flankRuleProblems(rule0)` and `userRuleProblems` over enabled rules.
- `describeRule0(r0)` → `${combine} of: ${crit} · segments ≤ ${maxBins} bins · signs disagree: ${mergeSignMismatch ? 'merged' : 'held back'} · fails Rule 0 with signs disagreeing: ${mergeSignsOnly ? 'merged last' : 'not merged'}`. For a copy (`!isRule0`), omit the last part.
- `rule0Lines(e, row)`:
  - one `criterionSideLine` per criterion;
  - `SRD/√φ flank signs ${e.signMismatch ? 'disagree' : 'agree or N/A'}`;
  - if `e.flankConflict`: `criteria pick different flanks → toward ${e.side}`;
  - the outcome line:

| Outcome | Line |
|---|---|
| `would_merge` | `Rule 0 → merges into ${side}` (+ `, blocked by genomic gap` when gap-blocked) |
| `held_back` | `Rule 0 → ⚠ held back: signs disagree (would merge into ${side})` |
| `signs_merge` | `Rule 0 → fails the test; signs disagree, merges last into ${side}` (or `… blocked by genomic gap` / `… no flank value`) |
| `signs_only` | `Rule 0 → fails the test; signs disagree (not merged)` |
| `none` | `Rule 0 → fails the test` |

  - Segments over the limit: `Rule 0: segment longer than the size limit, not evaluated`.
- `rule0CaseCounts`:
  - `would_merge`: blocked → `gap++`, else `wouldMerge++`; plus `signs++` on a sign mismatch.
  - `held_back` → `heldBack++`.
  - `signs_only` → `signsOnly++`.
  - `signs_merge`: blocked → `signsMergeBlocked++`, else `signsMerge++`.
- `makeRule0Copy` per Interfaces.
- `ruleFlagCounts` is unchanged.

Keep `PresetName` and `makeRule` as they are. The menu entry for the copy is handled in Task 4 through `makeRule0Copy`.

- [ ] **Step 4: Store and migration**

- `integrationStore.ts`:
  - `SegSettings`: `mergeSignsOnly: boolean` replaces `mergeFlankConflict`.
  - `defaultSeg.mergeSignsOnly = false`.
  - `LabelLayers` and its defaults drop `breakpoints` and `depthColor`.
  - `mergeRules: UserRule[]`; `setMergeRules(rules: UserRule[])`.
- `migrateSegSettings`: delete `mergeFlankConflict`, and delete `labels.breakpoints` and `labels.depthColor` if present. Keep the existing migrations.
- `hatchMerge.ts`: `MergeStep` gains `phase: 1 | 2 | null`; request and settings types use `AnyRule`.

- [ ] **Step 5: Run the TS tests and the type-check build**

```bash
web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs && node tests/_bundle/merge_rules_helpers.mjs
web/node_modules/.bin/esbuild tests/track_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/track_helpers.mjs && node tests/_bundle/track_helpers.mjs
web/node_modules/.bin/esbuild tests/hatch_merge_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/hatch_merge_helpers.mjs && node tests/_bundle/hatch_merge_helpers.mjs
cd web && npx tsc --noEmit -p . ; cd ..
```

Expected: the three helper suites print `all … checks passed`. `tsc` errors are allowed only in `IntegrationSegmentationView.tsx`, `MergeRulesEditor.tsx` and `HatchMergeView.tsx`, which are fixed in Task 4. To keep this commit building, make the minimal edits those files need to compile:
- rename `mergeFlankConflict` → `mergeSignsOnly` in the `rule0FromSeg` memo deps and checkbox;
- remove `FLANK_COLOR` uses;
- pass `row` to `rule0Lines`;
- type `ig.mergeRules` as `UserRule`.

Then confirm `cd web && npm run build` succeeds.

- [ ] **Step 6: Commit**

```bash
git add web/src/mergeRules.ts web/src/hatchMerge.ts web/src/integrationStore.ts web/src/trackMetrics.ts tests/merge_rules_helpers.ts tests/track_helpers.ts web/src/components/IntegrationSegmentationView.tsx web/src/components/MergeRulesEditor.tsx web/src/components/HatchMergeView.tsx
git commit -m "feat: port Rule 0 revision 2 to TypeScript; store and session migration

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: UI — tracks, section 1, section 2

**Files:**
- Modify: `web/src/components/IntegrationSegmentationView.tsx`
- Modify: `web/src/components/MergeRulesEditor.tsx`
- Modify: `web/src/components/HatchMergeView.tsx`
- Modify: `web/src/styles.css` (only if needed)

**Interfaces:** Consumes everything Task 3 produces.

- [ ] **Step 1: Tracks — breakpoints** (around line 109)

Remove the `seg.labels.breakpoints` guard. Always colour lines with `DEPTH_COLORS[Math.min(bp.recursion_depth, DEPTH_COLORS.length - 1)]`. Remove `['breakpoints', …]` and `['depthColor', …]` from the layer checkbox list. Show the depth legend `<span>depth: …</span>` unconditionally, with the label `accepted breakpoints by recursion depth:`.

- [ ] **Step 2: Tracks — drawing per segment** (around lines 126–135)

Replace the Rule 0 drawing block with the spec §6 table:

```tsx
const merging = !!r0 && r0.outcome === 'would_merge'
const held = !!r0 && r0.outcome === 'held_back'
const signsMerge = !!r0 && r0.outcome === 'signs_merge'
const outline = (c: string) => shapes.push({ type: 'rect', xref: 'x', yref: `${yref} domain`, x0, x1, y0: 0.02, y1: 0.98, fillcolor: 'rgba(0,0,0,0)', line: { width: 1.5, color: c, dash: 'dot' }, layer: 'above' })
if (merging) { if (r0!.blocked === null) shapes.push(...hatchShapes(x0, x1, yref, 'x', rule0.hatch.color)); else outline(rule0.hatch.color) }
if (signsMerge) { if (r0!.blocked === null) shapes.push(...hatchShapes(x0, x1, yref, 'x', SIGNS_MERGE_COLOR)); else outline(SIGNS_MERGE_COLOR) }
if (r0 && r0.signMismatch && (merging || held || (r0.outcome === 'signs_only' && seg.labels.signs))) shapes.push(...hatchShapes(x0, x1, yref, '|', SIGN_COLOR))
```

- Remove the violet flank-conflict line.
- For user rules, the hatch loop over `ev.matched` must skip `'rule0'` and look rules up in `ig.mergeRules` (now `UserRule[]`). Both kinds have `hatch`.
- The marker `line.color` is `rule0.hatch.color` when merging, `SIGNS_MERGE_COLOR` when signsMerge, else `#111827`. Its width is 2 when merging or signsMerge. Its symbol is a diamond when merging, held or signsMerge.

- [ ] **Step 3: Tracks — hover text**

Build the hover from `rule0Lines(r0, row)`.
- For user rules that match: a conditions rule keeps `conditionDetail` lines. A copy (`isFlankRule`) shows `rule ${n} “${name}” matches` followed by `rule0Lines(evaluateRule0(row, { ...r, mergeSignsOnly: false }, ig.allowGapCrossing), row)`, with the `Rule 0 →` prefix replaced by `→`. Import `evaluateRule0`.
- The `claimNote` stays.
- Change the gap text to `genomic gap on the ${[row.gap_left && 'left', row.gap_right && 'right'].filter(Boolean).join(' and ')} flank`, shown only when a side has a gap.
- The click card (around line 308) uses `rule0Lines(ev.rule0, row)`.

- [ ] **Step 4: Section 1 controls** (around lines 259–270)

- Remove the `Merge flank-conflict segments` checkbox.
- After the sign checkbox, add:
  `<label className="chk-inline"><input type="checkbox" aria-label="Merge signs-only segments" checked={seg.mergeSignsOnly} onChange={(e) => ig.setSeg({ mergeSignsOnly: e.target.checked })} />also merge signs-disagree segments that fail Rule 0 (last, black ×)</label>`
- Update the `rule0` memo deps.
- Keep the M=0 & G>0 notice as is.
- Update the lead paragraph (line ~248):
  - remove the violet sentence;
  - say breakpoints are coloured by recursion depth;
  - say merges go toward the flank of the criterion furthest below its threshold;
  - say a black × marks signs-disagree segments that fail Rule 0 and merge last, when that option is on.

- [ ] **Step 5: Legend** (around lines 281–285)

Replace the Rule 0 legend span:

```tsx
{rule0Issues.length === 0 && <span className="rule0-legend"><b className="rule-swatch" style={{ color: rule0.hatch.color }}>×</b>Rule 0 would merge now: {r0Counts.wouldMerge}{r0Counts.gap ? ` · blocked by gap: ${r0Counts.gap} (merge only with gap crossing)` : ''}{r0Counts.signs ? ` · signs disagree: ${r0Counts.signs}` : ''} · <b style={{ color: SIGN_COLOR }}>|</b> ⚠ held back: {r0Counts.heldBack}{seg.mergeSignsOnly
  ? <> · <b className="rule-swatch" style={{ color: SIGNS_MERGE_COLOR }}>×</b> signs disagree, fails Rule 0, merges last: {r0Counts.signsMerge}{r0Counts.signsMergeBlocked ? ` (+${r0Counts.signsMergeBlocked} blocked by gap)` : ''}</>
  : <> · <b style={{ color: SIGN_COLOR }}>|</b> signs disagree only: {r0Counts.signsOnly}</>}</span>}
```

Show the dotted-red gap legend entry only when some displayed row has `gap_left || gap_right`. Compute that with `useMemo` over `ig.hatchScores` for the displayed sources and chromosome.

- [ ] **Step 6: Section 2 editor**

In `MergeRulesEditor.tsx`:
- **Menu.** The "+ Add rule…" select gets a first option `<option value="rule0copy">Rule 0 copy (current section-1 settings)</option>` before the presets. On change, `'rule0copy'` appends `makeRule0Copy(ig.seg, rules)`; any other value uses `makeRule` as now.
- **Pinned Rule 0 card.** It uses `describeRule0(rule0)`. Its status line is `→ would merge now: ${r0c.wouldMerge}` + gap part + (when `rule0.mergeSignsOnly`) ` · merges last: ${r0c.signsMerge}`.
- **Card per rule.** Branch on `isFlankRule(r)`.
  - **Copy card:**
    - head: enabled checkbox, number, name input, move/duplicate/delete, as on conditions cards;
    - a combine `<select aria-label={`Rule ${n} combine`}>` all/any;
    - one checkbox per `RULE0_CRITERIA`, with aria-label `Rule ${n} criterion ${CRITERION_LABEL[c]}`;
    - number inputs `Rule ${n} tau`, `Rule ${n} SRD/√φ ref` and `Rule ${n} size limit`, where an empty value becomes NaN, as Rule 0 does;
    - a checkbox `Rule ${n} merge sign-disagreeing` (label "merge segments whose flank signs disagree");
    - hatch pattern and colour controls as on conditions cards;
    - problems from `flankRuleProblems(r)`, shown like `ruleProblems`;
    - count `flags ${counts[r.id] || 0}`.
  - **Duplicate a copy:** `{ ...r, id: <new unique id>, name: `${r.name} (copy)`, criteria: [...r.criteria], hatch: { ...r.hatch } }`. Generate the id the way `makeRule0Copy` does; export a small `newUserRuleId(existing)` from `mergeRules.ts` if needed.
  - **Conditions cards** are unchanged.
- **Merge request.** `canRunMerge(rule0, rules)` gates the run.

In `HatchMergeView.tsx`:
- The step description shows `phase 2 (last): ` before `rule "${step.rule_name}"` when `step.phase === 2`.
- The rules summary (`ruleSummary`) names copies as `“${name}” (Rule 0 copy)`.

- [ ] **Step 7: Build and manual check**

```bash
cd web && npm run build && cd ..
```

Start a separate server, without touching 8766: `PORT=8799 python3 run.py` (check `run.py` for the real port flag/env name), in the background. Open `http://127.0.0.1:8799` with the Playwright MCP or Chrome tools and check on chr8:
- no breakpoint checkboxes;
- depth-coloured lines;
- gap lines only at real gaps (expect few or none on chr8 internal segments);
- the hover shows L/R values with ✓/✗;
- toggling the new option draws black ×;
- "+ Add rule…" → "Rule 0 copy" adds an editable card;
- no console errors.

Stop the 8799 server afterwards.

- [ ] **Step 8: Commit**

```bash
git add web/src/components/IntegrationSegmentationView.tsx web/src/components/MergeRulesEditor.tsx web/src/components/HatchMergeView.tsx web/src/styles.css web/src/mergeRules.ts
git commit -m "feat: breakpoints by depth, signs-merge hatch, per-side hovers, Rule 0 copy cards

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

(Stage only the files you actually changed.)

---

### Task 5: Docs and browser walkthrough

**Files:**
- Modify: `tests/browser_hatch_merge_walkthrough.cjs`
- Modify: `validation/hatch-merge/merge_rules_check_results.json`, `validation/hatch-merge/merge-rules.png` (regenerated by the walkthrough)
- Modify: `README.md`, `docs/USAGE.md`, `docs/SCIENTIFIC_WORKFLOW.md`, `docs/ARCHITECTURE.md`, `docs/EXTENDING.md`, `tests/README.md`, `validation/rule0/README.md`

**Interfaces:** Consumes the Task 4 UI aria-labels: `Merge signs-only segments`, `Add rule`, `Rule N combine`, `Rule N tau`, `Rule N criterion …`.

- [ ] **Step 1: Update the walkthrough** (it takes `BASE_URL`)

- Remove the steps that toggle the flank-conflict switch, and any `labels.breakpoints`/`depthColor` handling.
- Remove the gap-crossing workaround the old script used because every segment was gap-adjacent. Instead, assert that on chr8 the number of segments whose hover or tooltip text contains `genomic gap on the` is less than the number of evaluated small segments. Read this from the plot's data with `page.evaluate` over the Plotly `hovertemplate`s.
- Assert that no element has aria-label `Merge flank-conflict segments`, and that the layer list has no "accepted breakpoints" checkbox.
- Tick `Merge signs-only segments` at a setting where `signs disagree only` was > 0 before ticking. Assert that the legend's `merges last:` count is > 0, or `+N blocked` is shown. Run the merge, and assert that the history (step slider walk or result JSON) contains a step with `phase === 2`. If chr8 has no such segment at any tried setting, try chr1–chr22 in turn and record which chromosome was used.
- Pick "Rule 0 copy (current section-1 settings)" from `Add rule`. Assert that a card with `Rule 1 combine` appears. Set `Rule 1 tau` to empty and assert that "Run merge" is disabled. Restore it and assert that it is enabled.
- Keep the existing checks (switch isolation for `mergeSignMismatch`, all/any, audit reason, no console errors). Drop the flank-switch rows from the isolation table.
- Write results to `validation/hatch-merge/merge_rules_check_results.json` and the screenshot to `merge-rules.png`, as now.

- [ ] **Step 2: Run it against a fresh server on port 8799** (never 8766)

```bash
(PORT=8799 python3 run.py > /tmp/claude-sx-8799.log 2>&1 &) ; sleep 8
BASE_URL=http://127.0.0.1:8799 node tests/browser_hatch_merge_walkthrough.cjs
```

Use the real port mechanism from `run.py`. Expected: all checks pass, with zero console errors. Then stop only the 8799 process (find it by its port, e.g. `fuser -k 8799/tcp` or `kill` of the PID you started).

- [ ] **Step 3: Docs**

Update every place that mentions `mergeFlankConflict`, the violet ― marker, the "flanks disagree" legend entry, the breakpoint layer checkboxes, or gap-blocked counts caused by the bug: `grep -rn "mergeFlankConflict\|violet\|flanks disagree\|depthColor\|colour by recursion depth" README.md docs tests/README.md validation web/README.md`. Leave untouched: `docs/superpowers/` (specs and plans are history) and the two untracked plan files.

State:
- direction toward the criterion furthest below its threshold;
- the `mergeSignsOnly` option and its last phase (black ×). Include the scientific-framing sentence: "This phase merges segments that fail every displayed threshold, solely because their flank SRD/√φ signs disagree; it is exploratory.";
- Rule 0 copies in section 2;
- `phase` on steps;
- the gap fix, linking `validation/rule0/gap_fix.md` from `validation/rule0/README.md`.

In README's API paragraph, update the Rule 0 object to the Global Constraints shape. Mention copies (`type: "flank_threshold"` user rules) and the `(phase, rule, score, start)` order.

- [ ] **Step 4: Full verification**

```bash
python3 -m pytest -q
for t in merge_rules_helpers hatch_merge_helpers track_helpers; do web/node_modules/.bin/esbuild tests/$t.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/$t.mjs && node tests/_bundle/$t.mjs | tail -1; done
cd web && npm run build && cd ..
git status --short
```

Expected:
- all pytest tests pass;
- three `all … checks passed`;
- the build succeeds;
- `git status` shows only your intended changes plus the two untracked plan files. `tests/_bundle/` is git-ignored; check that it is, and don't stage it.

- [ ] **Step 5: Commit**

```bash
git add tests/browser_hatch_merge_walkthrough.cjs validation/hatch-merge/merge_rules_check_results.json validation/hatch-merge/merge-rules.png README.md docs/USAGE.md docs/SCIENTIFIC_WORKFLOW.md docs/ARCHITECTURE.md docs/EXTENDING.md tests/README.md validation/rule0/README.md
git commit -m "docs, test: Rule 0 revision 2 walkthrough and documentation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

(Stage only files that changed.)
