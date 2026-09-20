# Integration hatching and exploratory merging — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add three new, genome-wide, user-configurable hatch layers (transition zone, merge proposal, ambiguous flank preference) and a separate exploratory adjacent-segment merge tool to Integration · segmentation & CN → Chromosome tracks, reusing existing scientific formulas and clearly distinguishing new results from the archived merge prototype / SRD screen.

**Architecture:** A new, small, pure-Python numerical module (`server/hatch_merge.py`) computes per-segment flank scores (SRD, SRD/√φ, delta log2FC with a selectable estimator) and runs the greedy exploratory merge; it reuses `flank_view.summarize` for level estimators and ports (with agreement tests) the existing signed-root-deviance and pooled-dispersion formulas from `reference/integration-prototype/src/srd_effect`. Two new FastAPI routes expose raw per-segment scores (for live, threshold-adjustable hatch classification) and the full merge step history (for undo/reset by array indexing, no server-side session state). The frontend gets a TypeScript port of the same three classification rules (transition / merge-proposal / ambiguous) for cheap client-side reclassification as thresholds change, a small `integrationStore.ts` extension for settings persistence and result caching, three new toggle+control rows wired into the existing Chromosome tracks layer controls (drawing hatches on the existing `trackFig` via the existing `hatchShapes` helper), and a new self-contained `HatchMergeView.tsx` panel (modeled on `AlternativeFlankView.tsx`) for the merge tool's own original-vs-merged plot, history table, reset/undo and exports.

**Tech Stack:** Python 3.10 / NumPy / SciPy (`xlogy`) / FastAPI / Pydantic on the backend; React 18 / TypeScript / Zustand / Plotly.js (via the existing `Plot`/`DeferredPlot` wrapper) on the frontend; pytest for backend tests.

**Spec:** `../../../INTEGRATION_HATCHING_AND_MERGING_TASK.md` (this plan implements it in full; read both together)

## Global Constraints

- Modify `scientific-explorer` only. Do not modify `interactive-explorer`, original scientific inputs, or archived reference outputs.
- Reuse existing scientific formulas (`flank_view.summarize` for mean/median/two-sided-IQR-mean levels; the exposure-aware signed-root-deviance and pooled-Pearson-dispersion formulas from `reference/integration-prototype/src/srd_effect/{srd_pruning,flank_dispersion}.py`). Do not introduce a replacement statistic or an unlabelled pseudocount.
- Keep thresholds separate per metric (SRD, SRD/√φ, delta log2FC never share a numeric threshold).
- Absolute-magnitude comparisons for "weaker"/"lower" flank selection; signed values only for direction.
- Zero is not a valid "opposite sign"; missing/undefined scores are never treated as zero or as eligible.
- Threshold comparisons use strict `<`.
- Compute genome-wide (any chromosome), not restricted to the seven pilot chromosomes with archived tables.
- Use original normalized pseudobulk X (`pb_x`), never smoothed/clipped/log-transformed display curves.
- Chromosome-local half-open retained-bin intervals internally; actual bin/bp coordinates for plotting.
- Never silently merge across a genomic gap; default-disallow, expose an explicit opt-in.
- The merge rule operates within one selected source/chromosome at a time, never across clusters or chromosomes.
- Transition/ambiguous classifications must not silently veto merge eligibility; any veto is an explicit, separately labelled setting.
- Do not silently inherit the original prototype's small-segment size cutoff; any size restriction here is a new, explicitly labelled, off-by-default setting.
- Never present the new merge tool's output as the original bootstrap/chain-guard merge prototype's result; preserve the original prototype's tables/UI untouched.
- Do not recompute or display CN calls on merged boundaries (out of scope here — flag as a stated limitation).
- Keep `DeferredPlot` and the existing rendering fixes; do not eagerly initialize every off-screen plot.
- Avoid iterative scientific computation inside React render; expensive numeric work (SRD/φ) happens once per backend fetch, cheap threshold comparisons happen client-side.
- Do not rewrite unrelated tabs, `flank_view.py`, or the original merge prototype code.

---

## File Structure

- Create `server/hatch_merge.py` — pure numerical module: segment building, level summaries (via `flank_view.summarize`), SRD/pooled-φ (ported), the three classification rules, and the greedy merge algorithm. No FastAPI/h5py/pandas imports.
- Create `tests/test_hatch_merge.py` — all of task section 6's numerical test list, plus agreement tests against the reference implementation.
- Create `server/hatch_merge_api.py` — two FastAPI routes (`/api/integration/hatch-scores`, `/api/integration/hatch-merge`) translating `IntegrationStore` data into `hatch_merge.py` calls.
- Modify `server/app.py` — register the new routes (one import line + one call, mirroring `register_integration`).
- Create `web/src/hatchMerge.ts` — types matching the backend responses; a thin `fetch` client; TypeScript ports of the three classification rules (`classifyTransition`, `classifyMergeProposal`, `classifyAmbiguous`) for cheap client-side reclassification; hatch-style mapping reusing `hatchShapes` from `integration.ts`.
- Create `tests/hatch_merge_helpers.ts` — unit tests for the TS classifiers (mirrors the Python test cases), run the same way `tests/experiment_overlays.ts` is (see that file for the harness pattern).
- Modify `web/src/integrationStore.ts` — new `hatchLayers` settings slice (persisted via `serialize`/`restore`), a `hatchScores` cache keyed by `source|chrom`, and `ensureHatchScores(source, chrom)`.
- Modify `web/src/components/IntegrationSegmentationView.tsx` — three new layer checkboxes + per-layer metric/threshold controls in the existing Chromosome tracks `controls-row`; hatch shapes for the new layers folded into the existing `trackFig` shapes array; legend text; clear "(new, computed here)" labelling distinct from the existing archived-only "flank-consistency hatch"/"SRD-screen proposal" layers.
- Create `web/src/components/HatchMergeView.tsx` — the exploratory merge tool: source/chromosome/metric/estimator/threshold/veto/size-restriction controls, "Run merge" action, original-vs-merged segment plot, merge-history table, step slider (reset/undo/inspect), JSON/CSV export. Mounted inside `IntegrationSegmentationView.tsx` near `AlternativeFlankView`.
- Modify `scientific-explorer/START_HERE.md` and `scientific-explorer/README.md` — document the new controls (mirrors how the recursion-depth feature was documented).
- Create `scientific-explorer/validation/hatch-merge/HANDOFF.md` — the required handoff/validation record (task section "Verify before reporting completion").

---

### Task 1: Backend numerical module — scoring

**Files:**
- Create: `server/hatch_merge.py`
- Test: `tests/test_hatch_merge.py`

**Interfaces:**
- Consumes: `server.flank_view.summarize(values) -> dict` (existing; keys `n, mean, median, iqr_mean, retained, std, cv, zero_fraction`).
- Produces (used by Task 2 and by Task 1's own merge code in Task 1 continued below):
  - `flat_partition(boundaries: list[int], n_bins: int) -> list[tuple[int,int]]`
  - `segment_summary(profile: np.ndarray, s: int, e: int) -> dict` (adds `start`, `end`, `sum` to `summarize`'s keys)
  - `signed_srd(y_t: float, e_t: float, y_r: float, e_r: float) -> float`
  - `pooled_phi(profile: np.ndarray, seg_bounds: list[tuple[int,int]]) -> float | None`
  - `estimator_level(summary: dict, estimator: str) -> float`
  - `delta_log2fc(seg_summary: dict, flank_summary: dict, estimator: str) -> float | None`
  - `flank_score(metric: str, seg_summary: dict, flank_summary: dict, phi: float | None, estimator: str = 'mean') -> float | None`
  - `is_gap(var_start: np.ndarray, var_end: np.ndarray, left_end_bin: int, right_start_bin: int) -> bool`

- [ ] **Step 1: Write the failing tests for segment building and level summaries**

```python
# tests/test_hatch_merge.py
from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from server import hatch_merge as hm


def test_flat_partition_half_open_and_drops_empty():
    assert hm.flat_partition([3, 3, 7], 10) == [(0, 3), (3, 7), (7, 10)]
    assert hm.flat_partition([], 5) == [(0, 5)]


def test_segment_summary_adds_start_end_sum_to_summarize():
    profile = np.array([2., 2., 4., 4.])
    s = hm.segment_summary(profile, 1, 3)
    assert s['start'] == 1 and s['end'] == 3
    assert s['sum'] == pytest.approx(6.0)
    assert s['mean'] == pytest.approx(3.0)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'server.hatch_merge'`

- [ ] **Step 3: Implement segment building and summaries**

```python
# server/hatch_merge.py
"""New exploratory hatching and adjacent-segment merging for the segmentation
tab's Chromosome tracks. Deliberately separate from the original bootstrap/
chain-guard merge prototype (analysis/depth_segment_prototype/merge_core.py)
and from the archived SRD screen: this is a new, independent, exploratory
rule, not a replacement or recomputation of either.

Level estimators (mean/median/two-sided-IQR-mean) are reused unmodified from
flank_view.summarize. The signed-root-deviance and pooled-Pearson-dispersion
formulas are ported from reference/integration-prototype/src/srd_effect
(srd_pruning.exposure_aware_srd, flank_dispersion.pooled_phi); tests in
test_hatch_merge.py assert numerical agreement with those reference functions
on shared examples so the two never silently drift apart.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.special import xlogy

from .flank_view import summarize

METRICS = ('srd', 'srd_phi', 'delta_log2fc')
ESTIMATORS = ('mean', 'median', 'iqr_mean')


def flat_partition(boundaries, n_bins):
    """Half-open chromosome-local segments from internal boundary bin offsets."""
    cuts = [0, *sorted(int(b) for b in boundaries), int(n_bins)]
    return [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1) if cuts[i + 1] > cuts[i]]


def segment_summary(profile, s, e):
    """summarize() over profile[s:e], plus start/end/sum (sum needed for SRD)."""
    x = np.asarray(profile, dtype=float)[s:e]
    return {**summarize(x), 'start': int(s), 'end': int(e), 'sum': float(x.sum())}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add server/hatch_merge.py tests/test_hatch_merge.py
git commit -m "feat(hatch-merge): segment partitioning and level summaries"
```

---

### Task 2: Backend numerical module — SRD, pooled φ, agreement tests

**Files:**
- Modify: `server/hatch_merge.py`
- Test: `tests/test_hatch_merge.py`

**Interfaces:**
- Consumes: Task 1's `segment_summary`.
- Produces: `signed_srd`, `pooled_phi`, `estimator_level`, `delta_log2fc`, `flank_score`, `is_gap` (used by Task 3's classifiers and Task 4's merge loop).

- [ ] **Step 1: Write failing tests, including agreement with the reference implementation**

```python
def test_signed_srd_matches_reference_examples():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'reference' / 'integration-prototype' / 'src' / 'srd_effect'))
    from srd_pruning import exposure_aware_srd  # reference-only import, test comparison only
    cases = [(10., 5., 10., 5.), (20., 5., 5., 5.), (0., 4., 0., 4.), (7., 3., 2., 9.)]
    for y_t, e_t, y_r, e_r in cases:
        assert hm.signed_srd(y_t, e_t, y_r, e_r) == pytest.approx(exposure_aware_srd(y_t, e_t, y_r, e_r))


def test_signed_srd_zero_pair_and_validation():
    assert hm.signed_srd(0., 4., 0., 4.) == 0.0
    with pytest.raises(ValueError):
        hm.signed_srd(1., 0., 1., 4.)
    with pytest.raises(ValueError):
        hm.signed_srd(-1., 4., 1., 4.)


def test_pooled_phi_matches_reference_and_handles_short_segments():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'reference' / 'integration-prototype' / 'src' / 'srd_effect'))
    from flank_dispersion import pooled_phi as ref_pooled_phi
    rng = np.random.default_rng(0)
    profile = rng.poisson(5, size=40).astype(float)
    bounds = [(0, 5), (5, 6), (6, 20), (20, 40)]  # one segment with only 1 bin
    seg_ids = np.zeros(40, dtype=int)
    for i, (s, e) in enumerate(bounds):
        seg_ids[s:e] = i
    assert hm.pooled_phi(profile, bounds) == pytest.approx(ref_pooled_phi(profile, seg_ids))


def test_pooled_phi_none_when_no_segment_qualifies():
    assert hm.pooled_phi(np.array([0., 0., 0.]), [(0, 1), (1, 3)]) is None


def test_delta_log2fc_undefined_on_zero_level_no_pseudocode():
    seg = hm.segment_summary(np.array([0., 0.]), 0, 2)
    flank = hm.segment_summary(np.array([2., 2.]), 0, 2)
    assert hm.delta_log2fc(seg, flank, 'mean') is None
    flank2 = hm.segment_summary(np.array([4., 4.]), 0, 2)
    assert hm.delta_log2fc(hm.segment_summary(np.array([2., 2.]), 0, 2), flank2, 'mean') == pytest.approx(-1.0)


def test_flank_score_srd_phi_none_when_phi_undefined():
    seg = hm.segment_summary(np.array([5., 5.]), 0, 2)
    flank = hm.segment_summary(np.array([2., 2.]), 0, 2)
    assert hm.flank_score('srd_phi', seg, flank, None) is None
    assert hm.flank_score('srd_phi', seg, flank, 0.0) is None
    assert hm.flank_score('srd', seg, flank, None) is not None


def test_is_gap_detects_bp_discontinuity():
    var_start = np.array([0, 100, 250])
    var_end = np.array([100, 200, 350])
    assert hm.is_gap(var_start, var_end, 1, 2) is True   # 200 != 250
    assert hm.is_gap(var_start, var_end, 0, 1) is False  # 100 == 100
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: FAIL — `AttributeError: module 'server.hatch_merge' has no attribute 'signed_srd'` (and similar for the others)

- [ ] **Step 3: Implement the ported formulas and score dispatch**

```python
# append to server/hatch_merge.py

def signed_srd(y_t, e_t, y_r, e_r):
    """Signed root deviance between two Poisson-like exposures.

    Ported from reference/integration-prototype/src/srd_effect/srd_pruning.py
    (exposure_aware_srd): D = 2*[xlogy(y_t,y_t/mu_t) + xlogy(y_r,y_r/mu_r)],
    lambda_hat = (y_t+y_r)/(e_t+e_r), mu_* = e_* * lambda_hat. Both-zero pair
    returns 0.0. Raises on nonpositive exposure or negative depth.
    """
    if e_t <= 0 or e_r <= 0:
        raise ValueError('exposure must be strictly positive')
    if y_t < 0 or y_r < 0:
        raise ValueError('depth values must be nonnegative')
    y_sum = y_t + y_r
    if y_sum == 0:
        return 0.0
    lam = y_sum / (e_t + e_r)
    mu_t, mu_r = e_t * lam, e_r * lam
    d = 2.0 * (xlogy(y_t, y_t / mu_t) + xlogy(y_r, y_r / mu_r))
    return float(np.sign(y_t / e_t - y_r / e_r) * np.sqrt(max(float(d), 0.0)))


def pooled_phi(profile, seg_bounds):
    """Pooled Pearson dispersion over segments with >=2 bins and mean>0.

    Ported from reference/integration-prototype/src/srd_effect/flank_dispersion.py
    (pooled_phi): sum_S sum_b (p-m_S)^2/m_S / sum_S (E_S-1), taking (start,end)
    bounds instead of a seg_ids array. This is pooled across EVERY segment of
    the current partition for one source/chromosome, so it must be recomputed
    whenever that partition changes (see run_merge below): SRD/sqrt(phi) is
    not a per-segment-pair quantity.
    """
    profile = np.asarray(profile, dtype=float)
    num, df = 0.0, 0
    for s, e in seg_bounds:
        p = profile[s:e]
        if len(p) < 2:
            continue
        m = p.mean()
        if m <= 0:
            continue
        num += float(((p - m) ** 2 / m).sum())
        df += len(p) - 1
    return float(num / df) if df > 0 else None


def estimator_level(summary, estimator):
    if estimator not in ESTIMATORS:
        raise ValueError(f'Unknown estimator {estimator!r}')
    return summary[estimator if estimator != 'iqr_mean' else 'iqr_mean']


def delta_log2fc(seg_summary, flank_summary, estimator):
    """Signed log2(segment level / flank level) for the selected estimator.

    Reuses flank_view.summarize's estimators via estimator_level. No
    pseudocount: undefined (None) whenever either level is not > 0.
    """
    a, b = estimator_level(seg_summary, estimator), estimator_level(flank_summary, estimator)
    return float(np.log2(a / b)) if a > 0 and b > 0 else None


def flank_score(metric, seg_summary, flank_summary, phi, estimator='mean'):
    """Signed score of seg relative to flank for the selected metric.

    metric='delta_log2fc' ignores phi. metric='srd' returns the signed root
    deviance directly. metric='srd_phi' divides by sqrt(phi) (never by phi
    itself) and is None when phi is missing, nonfinite or non-positive.
    """
    if metric not in METRICS:
        raise ValueError(f'Unknown metric {metric!r}')
    if metric == 'delta_log2fc':
        return delta_log2fc(seg_summary, flank_summary, estimator)
    srd = signed_srd(seg_summary['sum'], seg_summary['end'] - seg_summary['start'],
                      flank_summary['sum'], flank_summary['end'] - flank_summary['start'])
    if metric == 'srd':
        return srd
    if phi is None or not np.isfinite(phi) or phi <= 0:
        return None
    return float(srd / np.sqrt(phi))


def is_gap(var_start, var_end, left_end_bin, right_start_bin):
    """True if the two retained bins straddling a boundary are not genomically
    contiguous (a filtered/omitted interval sits between them)."""
    return float(var_start[right_start_bin]) != float(var_end[left_end_bin - 1])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: PASS (all tests so far)

- [ ] **Step 5: Commit**

```bash
git add server/hatch_merge.py tests/test_hatch_merge.py
git commit -m "feat(hatch-merge): port signed SRD and pooled phi with reference agreement tests"
```

---

### Task 3: Backend numerical module — classification rules

**Files:**
- Modify: `server/hatch_merge.py`
- Test: `tests/test_hatch_merge.py`

**Interfaces:**
- Consumes: Task 2's `flank_score`, `delta_log2fc`.
- Produces: `classify_transition(left, right) -> bool`, `weaker_side(left, right) -> str | None` (`'left' | 'right' | 'tie' | None`), `classify_merge_proposal(left, right, threshold) -> dict`, `classify_ambiguous(srd_left, srd_right, fc_left, fc_right) -> str` (`'consistent' | 'ambiguous' | 'no_unique_preference' | 'undefined'`). Used by Task 4 (merge vetoes), Task 5 (API), and mirrored in TS by Task 6.

- [ ] **Step 1: Write the failing tests (covers task section 6's classification bullets)**

```python
def test_classify_transition_opposite_same_zero_missing():
    assert hm.classify_transition(1.0, -1.0) is True
    assert hm.classify_transition(1.0, 1.0) is False
    assert hm.classify_transition(0.0, -1.0) is False    # zero is not an opposite sign
    assert hm.classify_transition(None, -1.0) is False
    assert hm.classify_transition(float('nan'), 1.0) is False


def test_weaker_side_uses_absolute_magnitude_when_both_negative():
    assert hm.weaker_side(-1.0, -3.0) == 'left'
    assert hm.weaker_side(-3.0, -1.0) == 'right'
    assert hm.weaker_side(-2.0, 2.0) == 'tie'
    assert hm.weaker_side(None, 1.0) is None


def test_classify_merge_proposal_strict_threshold_and_sides():
    assert hm.classify_merge_proposal(0.5, 5.0, 1.0) == {'eligible_left': True, 'eligible_right': False, 'weaker_side': 'left'}
    assert hm.classify_merge_proposal(1.0, 5.0, 1.0)['eligible_left'] is False  # strict <, not <=
    assert hm.classify_merge_proposal(0.5, 0.5, 1.0)['weaker_side'] == 'both'
    assert hm.classify_merge_proposal(None, 5.0, 1.0) == {'eligible_left': False, 'eligible_right': False, 'weaker_side': None}


def test_classify_merge_proposal_terminal_segment_one_neighbor():
    # a terminal segment only ever supplies one side; the other stays ineligible/None, never invented
    assert hm.classify_merge_proposal(0.2, None, 1.0) == {'eligible_left': True, 'eligible_right': False, 'weaker_side': 'left'}


def test_classify_ambiguous_differs_ties_and_undefined():
    assert hm.classify_ambiguous(1.0, 3.0, 3.0, 1.0) == 'ambiguous'      # weaker SRD=left, weaker log2fc=right
    assert hm.classify_ambiguous(1.0, 3.0, 1.0, 3.0) == 'consistent'
    assert hm.classify_ambiguous(2.0, 2.0, 1.0, 3.0) == 'no_unique_preference'  # exact SRD tie
    assert hm.classify_ambiguous(1.0, 3.0, None, 3.0) == 'undefined'
    assert hm.classify_ambiguous(1.0, float('nan'), 1.0, 3.0) == 'undefined'
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: FAIL — `AttributeError: module 'server.hatch_merge' has no attribute 'classify_transition'`

- [ ] **Step 3: Implement the classifiers**

```python
# append to server/hatch_merge.py

def _finite(v):
    return v is not None and np.isfinite(v)


def classify_transition(left_score, right_score):
    """True iff both scores are finite and have strictly opposite signs.

    Zero never counts as an opposite sign; a missing/nonfinite score never
    counts as a transition. Descriptive only, not a merge decision.
    """
    if not (_finite(left_score) and _finite(right_score)):
        return False
    return (left_score * right_score) < 0


def weaker_side(left_score, right_score):
    """Side with the smaller absolute score. 'tie' on exact equality (not
    treated as evidence either way); None if either side is missing."""
    if left_score is None or right_score is None:
        return None
    al, ar = abs(left_score), abs(right_score)
    if al == ar:
        return 'tie'
    return 'left' if al < ar else 'right'


def classify_merge_proposal(left_score, right_score, threshold):
    """Eligible sides under abs(score) < threshold (strict), and the weaker
    eligible side when both qualify. A terminal segment's missing side is
    never invented as eligible."""
    el = _finite(left_score) and abs(left_score) < threshold
    er = _finite(right_score) and abs(right_score) < threshold
    side = None
    if el and er:
        w = weaker_side(left_score, right_score)
        side = 'both' if w == 'tie' else w
    elif el:
        side = 'left'
    elif er:
        side = 'right'
    return {'eligible_left': bool(el), 'eligible_right': bool(er), 'weaker_side': side}


def classify_ambiguous(srd_left, srd_right, fc_left, fc_right):
    """'ambiguous' iff the weaker SRD-variant flank differs from the weaker
    delta-log2FC flank; 'consistent' if they agree; 'no_unique_preference' on
    an exact tie in either metric; 'undefined' if any input is missing or
    nonfinite. Requires an internal segment (both sides present by construction)."""
    vals = (srd_left, srd_right, fc_left, fc_right)
    if any(v is None or not np.isfinite(v) for v in vals):
        return 'undefined'
    srd_side, fc_side = weaker_side(srd_left, srd_right), weaker_side(fc_left, fc_right)
    if srd_side == 'tie' or fc_side == 'tie':
        return 'no_unique_preference'
    return 'consistent' if srd_side == fc_side else 'ambiguous'
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: PASS (all tests so far)

- [ ] **Step 5: Commit**

```bash
git add server/hatch_merge.py tests/test_hatch_merge.py
git commit -m "feat(hatch-merge): transition, merge-proposal and ambiguous classifiers"
```

---

### Task 4: Backend numerical module — genome-wide scoring and the greedy merge algorithm

**Files:**
- Modify: `server/hatch_merge.py`
- Test: `tests/test_hatch_merge.py`

**Interfaces:**
- Consumes: everything from Tasks 1-3.
- Produces:
  - `hatch_scores(profile: np.ndarray, boundaries: list[int], n_bins: int, var_start: np.ndarray, var_end: np.ndarray, chrom_offset: int) -> dict` — used by Task 5's API.
  - `run_merge(profile: np.ndarray, boundaries: list[int], n_bins: int, var_start: np.ndarray, var_end: np.ndarray, chrom_offset: int, metric: str, threshold: float, estimator: str = 'mean', veto_transition: bool = False, veto_ambiguous: bool = False, small_max_bins: int | None = None, allow_gap_crossing: bool = False) -> dict` — used by Task 5's API. Raises `ValueError` on bad `metric`/`estimator`/`threshold`.

- [ ] **Step 1: Write the failing tests**

```python
def _toy_profile_and_bounds():
    # chr-local bins 0..30; boundaries at 10 and 20 give three segments
    profile = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)])
    var_start = (np.arange(30) * 100).astype(float)
    var_end = var_start + 100
    return profile, [10, 20], 30, var_start, var_end


def test_hatch_scores_reports_every_internal_segment_with_phi_and_gap_flags():
    profile, boundaries, n_bins, var_start, var_end = _toy_profile_and_bounds()
    out = hm.hatch_scores(profile, boundaries, n_bins, var_start, var_end, chrom_offset=0)
    assert [tuple(s.values()) for s in out['segments']] == [(0, 10), (10, 20), (20, 30)]
    assert len(out['rows']) == 1  # only the middle segment has both flanks
    row = out['rows'][0]
    assert row['start'] == 10 and row['end'] == 20
    assert row['gap_left'] is False and row['gap_right'] is False
    assert row['phi'] == out['phi']
    assert row['delta_mean_left'] == pytest.approx(np.log2(4.2 / 4.0))
    assert row['srd_left'] is not None and row['srd_phi_right'] is not None


def test_hatch_scores_flags_a_genomic_gap():
    profile, boundaries, n_bins, var_start, var_end = _toy_profile_and_bounds()
    var_end = var_end.copy()
    var_end[9] += 500  # open a gap between bin 9 and bin 10 (the first boundary)
    out = hm.hatch_scores(profile, boundaries, n_bins, var_start, var_end, chrom_offset=0)
    assert out['rows'][0]['gap_left'] is True


def test_run_merge_rejects_bad_settings():
    profile, boundaries, n_bins, var_start, var_end = _toy_profile_and_bounds()
    with pytest.raises(ValueError):
        hm.run_merge(profile, boundaries, n_bins, var_start, var_end, 0, 'bogus', 1.0)
    with pytest.raises(ValueError):
        hm.run_merge(profile, boundaries, n_bins, var_start, var_end, 0, 'delta_log2fc', 1.0, estimator='bogus')
    with pytest.raises(ValueError):
        hm.run_merge(profile, boundaries, n_bins, var_start, var_end, 0, 'srd', -1.0)


def test_run_merge_merges_the_closest_pair_and_recomputes_before_next_step():
    # three near-identical segments (0-10, 10-20) close, (20-30) far: expect exactly one merge
    profile, boundaries, n_bins, var_start, var_end = _toy_profile_and_bounds()
    out = hm.run_merge(profile, boundaries, n_bins, var_start, var_end, 0, 'delta_log2fc', 0.5, estimator='mean')
    assert out['original'] == [{'start': 0, 'end': 10}, {'start': 10, 'end': 20}, {'start': 20, 'end': 30}]
    assert out['final'] == [{'start': 0, 'end': 20}, {'start': 20, 'end': 30}]
    assert len(out['steps']) == 2  # step 0 = original, step 1 = the one merge
    assert out['steps'][1]['removed_boundary'] == 10
    assert out['steps'][1]['merged_interval'] == [0, 20]
    assert out['steps'][0]['segments'] == out['original']


def test_run_merge_stops_when_nothing_eligible():
    profile, boundaries, n_bins, var_start, var_end = _toy_profile_and_bounds()
    out = hm.run_merge(profile, boundaries, n_bins, var_start, var_end, 0, 'delta_log2fc', 0.01, estimator='mean')
    assert out['final'] == out['original']
    assert len(out['steps']) == 1


def test_run_merge_deterministic_tie_break_leftmost():
    profile = np.concatenate([np.full(5, 4.0), np.full(5, 4.0), np.full(5, 4.0), np.full(5, 20.0)])
    var_start = (np.arange(20) * 100).astype(float); var_end = var_start + 100
    out = hm.run_merge(profile, [5, 10, 15], 20, var_start, var_end, 0, 'delta_log2fc', 1.0, estimator='mean')
    assert out['steps'][1]['removed_boundary'] == 5  # both (5,10)-pair and (10,15)-pair tie at score 0; leftmost wins


def test_run_merge_respects_gap_by_default_and_allows_opt_in():
    profile, boundaries, n_bins, var_start, var_end = _toy_profile_and_bounds()
    var_end = var_end.copy(); var_end[9] += 500  # gap right at the only cheap boundary
    blocked = hm.run_merge(profile, boundaries, n_bins, var_start, var_end, 0, 'delta_log2fc', 0.5, estimator='mean')
    assert blocked['final'] == blocked['original']
    allowed = hm.run_merge(profile, boundaries, n_bins, var_start, var_end, 0, 'delta_log2fc', 0.5, estimator='mean', allow_gap_crossing=True)
    assert allowed['final'] != allowed['original']


def test_run_merge_small_max_bins_is_off_by_default_and_explicit_when_set():
    profile = np.concatenate([np.full(50, 4.0), np.full(50, 4.2), np.full(50, 4.1)])
    var_start = (np.arange(150) * 100).astype(float); var_end = var_start + 100
    unrestricted = hm.run_merge(profile, [50, 100], 150, var_start, var_end, 0, 'delta_log2fc', 1.0, estimator='mean')
    assert unrestricted['final'] != unrestricted['original']  # all three segments are 50 bins; still merges
    restricted = hm.run_merge(profile, [50, 100], 150, var_start, var_end, 0, 'delta_log2fc', 1.0, estimator='mean', small_max_bins=10)
    assert restricted['final'] == restricted['original']  # neither side of the only boundary is <=10 bins


def test_run_merge_veto_transition_is_opt_in_and_does_not_veto_by_default():
    # middle segment is a transition zone (opposite-sign flanks) under delta_log2fc
    profile = np.concatenate([np.full(10, 8.0), np.full(10, 4.0), np.full(10, 8.4)])
    var_start = (np.arange(30) * 100).astype(float); var_end = var_start + 100
    default = hm.run_merge(profile, [10, 20], 30, var_start, var_end, 0, 'delta_log2fc', 3.0, estimator='mean')
    assert default['final'] != default['original']
    vetoed = hm.run_merge(profile, [10, 20], 30, var_start, var_end, 0, 'delta_log2fc', 3.0, estimator='mean', veto_transition=True)
    assert vetoed['final'] == vetoed['original']
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: FAIL — `AttributeError: module 'server.hatch_merge' has no attribute 'hatch_scores'`

- [ ] **Step 3: Implement genome-wide scoring and the merge loop**

```python
# append to server/hatch_merge.py

def hatch_scores(profile, boundaries, n_bins, var_start, var_end, chrom_offset):
    """Per internal segment: SRD, SRD/sqrt(phi) and delta-log2FC (all three
    estimators) against both flanks, plus the shared pooled phi and per-side
    genomic-gap flags. Computed genome-wide from live boundaries/profile, not
    restricted to chromosomes with archived diagnostic tables."""
    segs = flat_partition(boundaries, n_bins)
    summaries = [segment_summary(profile, s, e) for s, e in segs]
    phi = pooled_phi(profile, segs)
    rows = []
    for i in range(1, len(segs) - 1):
        seg, left, right = summaries[i], summaries[i - 1], summaries[i + 1]
        row = {
            'start': seg['start'], 'end': seg['end'], 'phi': phi,
            'gap_left': is_gap(var_start, var_end, chrom_offset + segs[i - 1][1], chrom_offset + segs[i][0]),
            'gap_right': is_gap(var_start, var_end, chrom_offset + segs[i][1], chrom_offset + segs[i + 1][0]),
        }
        for m in ('srd', 'srd_phi'):
            row[f'{m}_left'] = flank_score(m, seg, left, phi)
            row[f'{m}_right'] = flank_score(m, seg, right, phi)
        for est in ESTIMATORS:
            row[f'delta_{est}_left'] = delta_log2fc(seg, left, est)
            row[f'delta_{est}_right'] = delta_log2fc(seg, right, est)
        rows.append(row)
    return {'segments': [{'start': s, 'end': e} for s, e in segs], 'phi': phi, 'rows': rows}


def _segment_status(idx, segs, summaries, phi, metric, estimator):
    """Transition/ambiguous status of the internal segment at idx, or None for
    a terminal segment (fewer than two neighbors)."""
    if idx <= 0 or idx >= len(segs) - 1:
        return None
    seg, left, right = summaries[idx], summaries[idx - 1], summaries[idx + 1]
    l_score = flank_score(metric, seg, left, phi, estimator)
    r_score = flank_score(metric, seg, right, phi, estimator)
    srd_l = flank_score('srd_phi', seg, left, phi)
    srd_r = flank_score('srd_phi', seg, right, phi)
    fc_l = delta_log2fc(seg, left, estimator)
    fc_r = delta_log2fc(seg, right, estimator)
    return {
        'transition': classify_transition(l_score, r_score),
        'ambiguous': classify_ambiguous(srd_l, srd_r, fc_l, fc_r) == 'ambiguous',
    }


def _vetoed(i, segs, summaries, phi, metric, estimator, veto_transition, veto_ambiguous):
    if not (veto_transition or veto_ambiguous):
        return False
    for idx in (i, i + 1):
        status = _segment_status(idx, segs, summaries, phi, metric, estimator)
        if status is None:
            continue
        if veto_transition and status['transition']:
            return True
        if veto_ambiguous and status['ambiguous']:
            return True
    return False


def _eligible_boundaries(segs, summaries, phi, metric, estimator, threshold,
                          veto_transition, veto_ambiguous, small_max_bins,
                          allow_gap_crossing, var_start, var_end, chrom_offset):
    out = []
    for i in range(len(segs) - 1):
        left_s, left_e = segs[i]
        right_s, right_e = segs[i + 1]
        if not allow_gap_crossing and is_gap(var_start, var_end, chrom_offset + left_e, chrom_offset + right_s):
            continue
        if small_max_bins is not None and (left_e - left_s) > small_max_bins and (right_e - right_s) > small_max_bins:
            continue
        score = flank_score(metric, summaries[i], summaries[i + 1], phi, estimator)
        if score is None or not np.isfinite(score) or abs(score) >= threshold:
            continue
        if _vetoed(i, segs, summaries, phi, metric, estimator, veto_transition, veto_ambiguous):
            continue
        out.append((i, score))
    return out


def run_merge(profile, boundaries, n_bins, var_start, var_end, chrom_offset, metric, threshold,
              estimator='mean', veto_transition=False, veto_ambiguous=False, small_max_bins=None,
              allow_gap_crossing=False, max_steps=10000):
    """Greedy exploratory adjacent-segment merge, operating on one source's
    one-chromosome profile/boundaries. Repeatedly merges the eligible boundary
    (abs(score) < threshold, strict) with the smallest absolute score,
    tie-breaking on the leftmost boundary bin, recomputing every score
    (including the shared pooled phi) after each merge, until none remain
    eligible. Transition/ambiguous vetoes are off by default and explicit when
    set; a genomic-gap boundary is never crossed unless allow_gap_crossing is
    set. This is a new exploratory rule, independent of the original
    bootstrap/chain-guard merge prototype.
    """
    if metric not in METRICS:
        raise ValueError(f'Unknown merge metric {metric!r}')
    if estimator not in ESTIMATORS:
        raise ValueError(f'Unknown log2FC estimator {estimator!r}')
    if not (isinstance(threshold, (int, float)) and np.isfinite(threshold) and threshold > 0):
        raise ValueError('Threshold must be a positive finite number')
    segs = flat_partition(boundaries, n_bins)
    original = [{'start': s, 'end': e} for s, e in segs]
    steps = [{'step': 0, 'segments': list(original), 'removed_boundary': None,
              'merged_interval': None, 'score': None, 'phi': None}]
    guard = 0
    while True:
        guard += 1
        if guard > max_steps:
            raise RuntimeError('Merge did not terminate within max_steps')
        summaries = [segment_summary(profile, s, e) for s, e in segs]
        phi = pooled_phi(profile, segs) if metric != 'delta_log2fc' else None
        candidates = _eligible_boundaries(segs, summaries, phi, metric, estimator, threshold,
                                           veto_transition, veto_ambiguous, small_max_bins,
                                           allow_gap_crossing, var_start, var_end, chrom_offset)
        if not candidates:
            break
        i, score = min(candidates, key=lambda t: (abs(t[1]), segs[t[0]][1]))
        removed_boundary = segs[i][1]
        merged = (segs[i][0], segs[i + 1][1])
        segs = segs[:i] + [merged] + segs[i + 2:]
        steps.append({
            'step': len(steps), 'segments': [{'start': s, 'end': e} for s, e in segs],
            'removed_boundary': removed_boundary, 'merged_interval': [merged[0], merged[1]],
            'score': score, 'phi': phi,
        })
    return {
        'original': original, 'final': steps[-1]['segments'], 'steps': steps,
        'settings': {'metric': metric, 'estimator': estimator, 'threshold': threshold,
                     'veto_transition': veto_transition, 'veto_ambiguous': veto_ambiguous,
                     'small_max_bins': small_max_bins, 'allow_gap_crossing': allow_gap_crossing},
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: PASS (all tests)

- [ ] **Step 5: Commit**

```bash
git add server/hatch_merge.py tests/test_hatch_merge.py
git commit -m "feat(hatch-merge): genome-wide hatch scoring and greedy exploratory merge"
```

---

### Task 5: Backend API routes

**Files:**
- Create: `server/hatch_merge_api.py`
- Modify: `server/app.py:431-433` (add one import + one registration call)
- Test: `tests/test_hatch_merge.py` (route-level tests using FastAPI's TestClient against a minimal fake store)

**Interfaces:**
- Consumes: `hatch_merge.hatch_scores`, `hatch_merge.run_merge` (Task 4); `IntegrationStore` (`server/integration.py`, existing: `.arrays['pb_x']`, `.var_seq`, `.var_start`, `.var_end`, `.members`, `.tables['breakpoints']`, `.identity`) exactly as `flank_view.generate_view` already consumes it.
- Produces: `register_hatch_merge(app, get_integration)` — same signature convention as `register_integration` in `integration_api.py`. Routes:
  - `POST /api/integration/hatch-scores` — body `{chrom: str, source: str}` → `hatch_merge.hatch_scores(...)` result, plus `source`, `chrom`, `start_bp`/`end_bp` arrays (mirrors `flank_view.generate_view`'s bp arrays) so the frontend can map bins to Mb without a second round trip.
  - `POST /api/integration/hatch-merge` — body `{chrom, source, metric, threshold, estimator?, veto_transition?, veto_ambiguous?, small_max_bins?, allow_gap_crossing?}` → `hatch_merge.run_merge(...)` result, plus the same `start_bp`/`end_bp` arrays and a `provenance` object (mirrors `flank_view.generate_view`'s provenance dict).

- [ ] **Step 1: Write the failing route tests**

```python
# append to tests/test_hatch_merge.py
from fastapi import FastAPI
from fastapi.testclient import TestClient


class _FakeStore:
    """Minimal stand-in for IntegrationStore covering what the routes touch."""
    def __init__(self):
        self.identity = {'hash': 'fake', 'files': []}
        self.var_seq = np.array(['chr1'] * 30)
        self.var_start = (np.arange(30) * 100).astype(float)
        self.var_end = self.var_start + 100
        self.members = {'cluster0': np.arange(5), 'all': np.arange(5)}
        pb = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)])
        self.arrays = {'pb_x': np.tile(pb, (2, 1))}
        import pandas as pd
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


def test_hatch_scores_route_returns_rows_and_bp_arrays():
    resp = _client().post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'cluster0'})
    assert resp.status_code == 200
    body = resp.json()
    assert body['source'] == 'cluster0' and body['chrom'] == 'chr1'
    assert len(body['start_bp']) == 30
    assert len(body['rows']) == 1


def test_hatch_scores_route_rejects_unknown_source():
    resp = _client().post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'nope'})
    assert resp.status_code == 422


def test_hatch_merge_route_runs_and_reports_provenance():
    resp = _client().post('/api/integration/hatch-merge', json={
        'chrom': 'chr1', 'source': 'cluster0', 'metric': 'delta_log2fc', 'threshold': 0.5, 'estimator': 'mean'})
    assert resp.status_code == 200
    body = resp.json()
    assert body['final'] != body['original']
    assert 'provenance' in body and body['provenance']['dataset'] == {'hash': 'fake', 'files': []}


def test_hatch_merge_route_rejects_bad_threshold():
    resp = _client().post('/api/integration/hatch-merge', json={
        'chrom': 'chr1', 'source': 'cluster0', 'metric': 'srd', 'threshold': -1.0})
    assert resp.status_code == 422
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'server.hatch_merge_api'`

- [ ] **Step 3: Implement the routes**

```python
# server/hatch_merge_api.py
"""New exploratory hatching/merge routes. Kept separate from
integration_api.py's flank-view/bootstrap routes, which back the existing
(unmodified) Alternative flank view."""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Optional

import numpy as np
from fastapi import HTTPException
from pydantic import BaseModel, Field

from . import hatch_merge as hm
from .integration import IntegrationStore


class HatchScoresRequest(BaseModel):
    chrom: str
    source: str


class HatchMergeRequest(BaseModel):
    chrom: str
    source: str
    metric: str
    threshold: float
    estimator: str = 'mean'
    veto_transition: bool = False
    veto_ambiguous: bool = False
    small_max_bins: Optional[int] = Field(default=None, ge=1)
    allow_gap_crossing: bool = False


def _resolve(ig: IntegrationStore, chrom: str, source: str):
    if ig.arrays is None:
        raise HTTPException(503, 'Pseudobulk index is not ready')
    if source not in ig.members:
        raise HTTPException(422, f'Unknown source {source!r}')
    idx = np.flatnonzero(ig.var_seq == chrom)
    if not len(idx):
        raise HTTPException(422, f'Unknown chromosome {chrom!r}')
    if not np.array_equal(idx, np.arange(idx[0], idx[-1] + 1)):
        raise HTTPException(422, 'Chromosome bins are not contiguous')
    si = list(ig.members).index(source)  # ig.members is a dict {source: cell-index-array}; list() gives its keys in insertion order
    profile = ig.arrays['pb_x'][si, idx]
    bp = ig.tables['breakpoints']
    boundaries = bp.loc[(bp.source == source) & (bp.chromosome == chrom), 'absolute_bin'].tolist()
    return profile, boundaries, idx


def register_hatch_merge(app, get_integration):
    def _ig() -> IntegrationStore:
        ig = get_integration()
        if ig is None or not getattr(ig, 'available', False):
            raise HTTPException(503, 'integration data not loaded')
        return ig

    @app.post('/api/integration/hatch-scores')
    def hatch_scores(req: HatchScoresRequest):
        ig = _ig()
        try:
            profile, boundaries, idx = _resolve(ig, req.chrom, req.source)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e
        result = hm.hatch_scores(profile, boundaries, len(idx), ig.var_start, ig.var_end, chrom_offset=int(idx[0]))
        return {**result, 'source': req.source, 'chrom': req.chrom,
                'start_bp': ig.var_start[idx].tolist(), 'end_bp': ig.var_end[idx].tolist()}

    @app.post('/api/integration/hatch-merge')
    def hatch_merge(req: HatchMergeRequest):
        ig = _ig()
        try:
            profile, boundaries, idx = _resolve(ig, req.chrom, req.source)
            result = hm.run_merge(profile, boundaries, len(idx), ig.var_start, ig.var_end, chrom_offset=int(idx[0]),
                                   metric=req.metric, threshold=req.threshold, estimator=req.estimator,
                                   veto_transition=req.veto_transition, veto_ambiguous=req.veto_ambiguous,
                                   small_max_bins=req.small_max_bins, allow_gap_crossing=req.allow_gap_crossing)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e
        return {**result, 'source': req.source, 'chrom': req.chrom,
                'start_bp': ig.var_start[idx].tolist(), 'end_bp': ig.var_end[idx].tolist(),
                'provenance': {'dataset': ig.identity,
                               'implementation_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                               'scope': 'New exploratory merge rule; independent of the original bootstrap/chain-guard '
                                        'merge prototype and the archived SRD screen. No CN recomputation on merged boundaries.'}}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd scientific-explorer && python3 -m pytest -q tests/test_hatch_merge.py -v`
Expected: PASS (all tests)

- [ ] **Step 5: Register the routes in the running app**

Edit `server/app.py` at the block shown below (currently lines 431-433):

```python
from .experiment_api import register_experiments
register_experiments(app, _store, lambda: jobs, cfg)
register_integration(app, lambda: integration)
```

Change to:

```python
from .experiment_api import register_experiments
register_experiments(app, _store, lambda: jobs, cfg)
register_integration(app, lambda: integration)

from .hatch_merge_api import register_hatch_merge
register_hatch_merge(app, lambda: integration)
```

This stays before the `StaticFiles`/catch-all SPA mount later in the file, matching the existing routes (the handoff doc's lesson 1: API routes must be registered before the catch-all).

- [ ] **Step 6: Smoke-test against the real running server**

Run (from `scientific-explorer/`, with the app already serving — see `README.md` for `python3 run.py --port 8766`):
```bash
curl -s -X POST http://127.0.0.1:8766/api/integration/hatch-scores -H 'Content-Type: application/json' -d '{"chrom":"chr8","source":"cluster4"}' | python3 -m json.tool | head -30
```
Expected: 200 with `rows`, `segments`, `phi`, `start_bp`/`end_bp` populated (or a clear 503/422 if the integration store isn't loaded/available in that environment — in which case retry after `/api/integration/status` reports `"available": true`).

- [ ] **Step 7: Commit**

```bash
git add server/hatch_merge_api.py server/app.py tests/test_hatch_merge.py
git commit -m "feat(hatch-merge): expose hatch-scores and hatch-merge API routes"
```

---

### Task 6: Frontend types, API client, and TypeScript classifiers

**Files:**
- Create: `web/src/hatchMerge.ts`
- Create: `tests/hatch_merge_helpers.ts`

**Interfaces:**
- Consumes: `hatchShapes` from `web/src/integration.ts` (existing, signature `(x0: number, x1: number, yref: string, style: 'none' | '/' | 'x', color?: string) => any[]`).
- Produces (used by Task 7's store and Task 8/9's components):
  - `type HatchMetric = 'srd' | 'srd_phi' | 'delta_log2fc'`
  - `type Estimator = 'mean' | 'median' | 'iqr_mean'`
  - `interface HatchScoreRow { start: number; end: number; phi: number | null; gap_left: boolean; gap_right: boolean; srd_left: number | null; srd_right: number | null; srd_phi_left: number | null; srd_phi_right: number | null; delta_mean_left: number | null; delta_mean_right: number | null; delta_median_left: number | null; delta_median_right: number | null; delta_iqr_mean_left: number | null; delta_iqr_mean_right: number | null }`
  - `interface HatchScoresResult { segments: { start: number; end: number }[]; phi: number | null; rows: HatchScoreRow[]; source: string; chrom: string; start_bp: number[]; end_bp: number[] }`
  - `interface MergeStep { step: number; segments: { start: number; end: number }[]; removed_boundary: number | null; merged_interval: [number, number] | null; score: number | null; phi: number | null }`
  - `interface MergeResult { original: { start: number; end: number }[]; final: { start: number; end: number }[]; steps: MergeStep[]; settings: Record<string, unknown>; source: string; chrom: string; start_bp: number[]; end_bp: number[]; provenance: unknown }`
  - `fetchHatchScores(chrom: string, source: string): Promise<HatchScoresResult>`
  - `runHatchMerge(req: { chrom: string; source: string; metric: HatchMetric; threshold: number; estimator?: Estimator; veto_transition?: boolean; veto_ambiguous?: boolean; small_max_bins?: number; allow_gap_crossing?: boolean }): Promise<MergeResult>`
  - `weakerSide(left: number | null, right: number | null): 'left' | 'right' | 'tie' | null`
  - `classifyTransition(left: number | null, right: number | null): boolean`
  - `classifyMergeProposal(left: number | null, right: number | null, threshold: number): { eligibleLeft: boolean; eligibleRight: boolean; weakerSide: 'left' | 'right' | 'both' | null }`
  - `classifyAmbiguous(srdLeft: number | null, srdRight: number | null, fcLeft: number | null, fcRight: number | null): 'consistent' | 'ambiguous' | 'no_unique_preference' | 'undefined'`
  - `scoreForMetric(row: HatchScoreRow, side: 'left' | 'right', metric: HatchMetric, estimator: Estimator): number | null`

- [ ] **Step 1: Write the failing helper tests, mirroring the Python cases**

```typescript
// tests/hatch_merge_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/hatch_merge_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/hatch_merge_helpers.mjs
//           node tests/_bundle/hatch_merge_helpers.mjs
import { classifyTransition, classifyMergeProposal, classifyAmbiguous, weakerSide } from '../web/src/hatchMerge'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}

eq(classifyTransition(1, -1), true, 'transition opposite sign')
eq(classifyTransition(1, 1), false, 'transition same sign')
eq(classifyTransition(0, -1), false, 'transition zero is not opposite')
eq(classifyTransition(null, -1), false, 'transition missing left')
eq(classifyTransition(NaN, 1), false, 'transition NaN')

eq(weakerSide(-1, -3), 'left', 'weaker side both negative picks smaller magnitude')
eq(weakerSide(-2, 2), 'tie', 'weaker side exact tie')
eq(weakerSide(null, 1), null, 'weaker side missing')

eq(classifyMergeProposal(0.5, 5, 1), { eligibleLeft: true, eligibleRight: false, weakerSide: 'left' }, 'merge proposal one side')
eq(classifyMergeProposal(1, 5, 1), { eligibleLeft: false, eligibleRight: false, weakerSide: null }, 'merge proposal strict threshold')
eq(classifyMergeProposal(0.5, 0.5, 1), { eligibleLeft: true, eligibleRight: true, weakerSide: 'both' }, 'merge proposal tie -> both')
eq(classifyMergeProposal(0.2, null, 1), { eligibleLeft: true, eligibleRight: false, weakerSide: 'left' }, 'merge proposal terminal segment')

eq(classifyAmbiguous(1, 3, 3, 1), 'ambiguous', 'ambiguous differing weaker flanks')
eq(classifyAmbiguous(1, 3, 1, 3), 'consistent', 'ambiguous consistent weaker flanks')
eq(classifyAmbiguous(2, 2, 1, 3), 'no_unique_preference', 'ambiguous exact SRD tie')
eq(classifyAmbiguous(1, 3, null, 3), 'undefined', 'ambiguous missing input')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) }
console.log('all hatch_merge_helpers checks passed')
```

- [ ] **Step 2: Run to verify it fails**

Run:
```bash
cd scientific-explorer
web/node_modules/.bin/esbuild tests/hatch_merge_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/hatch_merge_helpers.mjs
```
Expected: esbuild error — `Could not resolve "../web/src/hatchMerge"`

- [ ] **Step 3: Implement `hatchMerge.ts`**

```typescript
// web/src/hatchMerge.ts
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
```

- [ ] **Step 4: Run to verify it passes**

Run:
```bash
cd scientific-explorer
web/node_modules/.bin/esbuild tests/hatch_merge_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/hatch_merge_helpers.mjs
node tests/_bundle/hatch_merge_helpers.mjs
```
Expected: `all hatch_merge_helpers checks passed`, exit code 0

- [ ] **Step 5: Commit**

```bash
git add web/src/hatchMerge.ts tests/hatch_merge_helpers.ts
git commit -m "feat(hatch-merge): frontend types, API client and ported classifiers"
```

---

### Task 7: `integrationStore.ts` — settings, caching, session persistence

**Files:**
- Modify: `web/src/integrationStore.ts`

**Interfaces:**
- Consumes: `fetchHatchScores`, `runHatchMerge`, `HatchScoresResult`, `MergeResult`, `HatchMetric`, `Estimator` from `web/src/hatchMerge.ts` (Task 6).
- Produces (used by Task 8 and Task 9):
  - New state: `hatchLayers: HatchLayerSettings`, `hatchScores: Record<string, HatchScoresResult>` (key `` `${source}|${chrom}` ``), `hatchScoresError: string | null`.
  - New actions: `setHatchLayers(p: Partial<HatchLayerSettings>): void`, `ensureHatchScores(source: string, chrom: string): Promise<void>`.
  - `HatchLayerSettings` shape (all persisted): `{ transition: { on: boolean; metric: HatchMetric; estimator: Estimator }; proposal: { on: boolean; metric: HatchMetric; estimator: Estimator; threshold: number }; ambiguous: { on: boolean; estimator: Estimator } }`.
  - `serialize()`/`restore()` include `hatchLayers` (settings only — fetched `hatchScores` results are not session-persisted, matching the existing `AlternativeFlankView` convention of local/regenerable results).

- [ ] **Step 1: Extend the store**

Read the current file first (`web/src/integrationStore.ts`, 82 lines) to get exact surrounding context before editing — the snippets below show the target end state of each touched region.

Add to the imports:
```typescript
import { fetchHatchScores, type Estimator, type HatchMetric, type HatchScoresResult } from './hatchMerge'
```

Add a new exported type near `LabelLayers`:
```typescript
export interface HatchLayerSettings {
  transition: { on: boolean; metric: HatchMetric; estimator: Estimator }
  proposal: { on: boolean; metric: HatchMetric; estimator: Estimator; threshold: number }
  ambiguous: { on: boolean; estimator: Estimator }
}
const defaultHatchLayers: HatchLayerSettings = {
  transition: { on: false, metric: 'srd_phi', estimator: 'mean' },
  proposal: { on: false, metric: 'srd_phi', estimator: 'mean', threshold: 3.3 },
  ambiguous: { on: false, estimator: 'mean' },
}
```

Add to `IntegrationState` interface:
```typescript
  hatchLayers: HatchLayerSettings
  hatchScores: Record<string, HatchScoresResult>
  hatchScoresError: string | null
  setHatchLayers: (p: Partial<HatchLayerSettings>) => void
  ensureHatchScores: (source: string, chrom: string) => Promise<void>
```

Add to the store body (initial state):
```typescript
  hatchLayers: defaultHatchLayers, hatchScores: {}, hatchScoresError: null,
```

Add to the store body (actions, alongside `ensureCandidates`):
```typescript
  setHatchLayers: (p) => set({ hatchLayers: { ...get().hatchLayers, ...p } }),
  ensureHatchScores: async (source, chrom) => {
    const k = `${source}|${chrom}`
    if (get().hatchScores[k]) return
    try { const r = await fetchHatchScores(chrom, source); set({ hatchScores: { ...get().hatchScores, [k]: r } }) }
    catch (e: any) { set({ hatchScoresError: `hatch-scores: ${e.message || e}` }) }
  },
```

Modify `serialize()` to include `hatchLayers`:
```typescript
  serialize: () => { const s = get(); return { clusterKey: s.clusterKey, colorBy: s.colorBy, selectedCluster: s.selectedCluster, highlight: s.highlight, focusCell: s.focusCell, dimOthers: s.dimOthers, pointSize: s.pointSize, compareKey: s.compareKey, seg: s.seg, hatchLayers: s.hatchLayers } },
```

Modify `restore()` to include `hatchLayers`:
```typescript
  restore: (p) => { if (!p) return; set({ clusterKey: p.clusterKey || 'wnn_leiden_0.3', colorBy: p.colorBy || 'cluster', selectedCluster: p.selectedCluster ?? null, highlight: p.highlight || [], focusCell: p.focusCell ?? null, dimOthers: !!p.dimOthers, pointSize: p.pointSize || 3, compareKey: p.compareKey || 'wnn_leiden_1.0', seg: { ...defaultSeg, ...(p.seg || {}), labels: { ...defaultSeg.labels, ...((p.seg || {}).labels || {}) } }, hatchLayers: { ...defaultHatchLayers, ...(p.hatchLayers || {}) } }) },
```

- [ ] **Step 2: Typecheck**

Run: `cd scientific-explorer/web && npx tsc --noEmit`
Expected: no new errors (existing unrelated errors, if any, are pre-existing — confirm with `git stash` + rerun if unsure)

- [ ] **Step 3: Commit**

```bash
git add web/src/integrationStore.ts
git commit -m "feat(hatch-merge): store settings, fetch/cache for hatch scores"
```

---

### Task 8: Wire the three new hatch layers into the Chromosome tracks panel

**Files:**
- Modify: `web/src/components/IntegrationSegmentationView.tsx`

**Interfaces:**
- Consumes: `ig.hatchLayers`, `ig.setHatchLayers`, `ig.hatchScores`, `ig.ensureHatchScores` (Task 7); `classifyTransition`, `classifyMergeProposal`, `classifyAmbiguous`, `weakerSide`, `scoreForMetric`, `hatchShapes`, `type HatchMetric`, `type Estimator` from `hatchMerge.ts` (Task 6).
- Produces: no new exports; extends the existing `trackFig` shapes/annotations and the existing layer-toggle `controls-row`.

- [ ] **Step 1: Fetch scores when a new layer is enabled**

Add near the existing `useEffect` that calls `ig.ensureCandidates` (around the current line `useEffect(() => { if (seg.labels.candidates && meta) ... }, [...])`):

```typescript
const hatchLayers = ig.hatchLayers
useEffect(() => {
  if (!meta) return
  if (hatchLayers.transition.on || hatchLayers.proposal.on || hatchLayers.ambiguous.on) for (const s of seg.sources) ig.ensureHatchScores(s, seg.chrom)
}, [hatchLayers.transition.on, hatchLayers.proposal.on, hatchLayers.ambiguous.on, seg.chrom, seg.sources, meta])
```

- [ ] **Step 2: Add the hatch shapes into `trackFig`'s per-source loop**

Inside the `srcs.forEach((src, k) => { ... })` block of `trackFig`'s `useMemo`, immediately after the existing "small segments" `if (seg.labels.small) for (const s of smallHere...) { ... }` block (i.e. after its closing `}`), add:

```typescript
      // new exploratory hatch layers (genome-wide; independent of the archived-only consistency/proposal layers above)
      const hs = ig.hatchScores[`${src}|${seg.chrom}`]
      if (hs && (hatchLayers.transition.on || hatchLayers.proposal.on || hatchLayers.ambiguous.on)) for (const row of hs.rows) {
        const x0 = mb[row.start], x1 = row.end < chromMeta.n ? mb[row.end] : bins[2 * (chromMeta.offset + chromMeta.n - 1) + 1] / MB
        const styles: ('none' | '/' | 'x')[] = []
        if (hatchLayers.transition.on) {
          const l = scoreForMetric(row, 'left', hatchLayers.transition.metric, hatchLayers.transition.estimator)
          const r = scoreForMetric(row, 'right', hatchLayers.transition.metric, hatchLayers.transition.estimator)
          if (classifyTransition(l, r)) styles.push('/')
        }
        if (hatchLayers.proposal.on) {
          const l = scoreForMetric(row, 'left', hatchLayers.proposal.metric, hatchLayers.proposal.estimator)
          const r = scoreForMetric(row, 'right', hatchLayers.proposal.metric, hatchLayers.proposal.estimator)
          const cls = classifyMergeProposal(l, r, hatchLayers.proposal.threshold)
          if (cls.weakerSide) styles.push('x')
        }
        if (hatchLayers.ambiguous.on) {
          const srdL = scoreForMetric(row, 'left', 'srd_phi', hatchLayers.ambiguous.estimator)
          const srdR = scoreForMetric(row, 'right', 'srd_phi', hatchLayers.ambiguous.estimator)
          const fcL = scoreForMetric(row, 'left', 'delta_log2fc', hatchLayers.ambiguous.estimator)
          const fcR = scoreForMetric(row, 'right', 'delta_log2fc', hatchLayers.ambiguous.estimator)
          if (classifyAmbiguous(srdL, srdR, fcL, fcR) === 'ambiguous') styles.push('x')
        }
        for (const style of styles) shapes.push(...hatchShapes(x0, x1, yref, style, 'rgba(21,94,117,0.6)'))
        if (row.gap_left || row.gap_right) shapes.push({ type: 'line', xref: 'x', yref: `${yref} domain`, x0: row.gap_left ? x0 : x1, x1: row.gap_left ? x0 : x1, y0: 0, y1: 1, line: { color: '#b91c1c', width: 2, dash: 'dot' }, opacity: 0.7 })
      }
```

- [ ] **Step 3: Add `hs`/`hatchLayers`/`ig.hatchScores` to the `trackFig` dependency array**

Change the existing dependency array (last line of the `trackFig` `useMemo`):
```typescript
  }, [meta, ig.pb, seg, smallHere, chromMeta, bins, effectMetric, metricDef, clSegs, ig.candidates, presentation, nBins])
```
to:
```typescript
  }, [meta, ig.pb, seg, smallHere, chromMeta, bins, effectMetric, metricDef, clSegs, ig.candidates, presentation, nBins, hatchLayers, ig.hatchScores])
```

- [ ] **Step 4: Add the three new imports**

At the top of the file, extend the existing `integration.ts` import line with nothing new (unchanged), and add a new import line right after it:
```typescript
import { classifyAmbiguous, classifyMergeProposal, classifyTransition, scoreForMetric, type Estimator, type HatchMetric } from '../hatchMerge'
```

- [ ] **Step 5: Add the three new toggle rows with per-layer controls**

In the existing layer-toggle `controls-row` (the one built from the array literal `[['breakpoints', ...], ..., ['allelic', ...]]`), add three new inline controls in a new `controls-row` directly below it (keeping the existing archived-layer row untouched, per the requirement to visibly distinguish new from archived layers):

```tsx
        <div className="controls-row">
          <span className="small muted">new exploratory hatch layers (genome-wide, computed here — separate from the archived layers above):</span>
          <label className="chk-inline"><input type="checkbox" checked={hatchLayers.transition.on} onChange={(e) => ig.setHatchLayers({ transition: { ...hatchLayers.transition, on: e.target.checked } })} />transition zone (/)</label>
          <select aria-label="Transition metric" value={hatchLayers.transition.metric} onChange={(e) => ig.setHatchLayers({ transition: { ...hatchLayers.transition, metric: e.target.value as HatchMetric } })}><option value="srd">SRD</option><option value="srd_phi">SRD/√φ</option><option value="delta_log2fc">delta log2FC</option></select>
          {hatchLayers.transition.metric === 'delta_log2fc' && <select aria-label="Transition estimator" value={hatchLayers.transition.estimator} onChange={(e) => ig.setHatchLayers({ transition: { ...hatchLayers.transition, estimator: e.target.value as Estimator } })}><option value="mean">mean</option><option value="median">median</option><option value="iqr_mean">two-sided IQR mean</option></select>}
          <label className="chk-inline"><input type="checkbox" checked={hatchLayers.proposal.on} onChange={(e) => ig.setHatchLayers({ proposal: { ...hatchLayers.proposal, on: e.target.checked } })} />merge proposal (×)</label>
          <select aria-label="Merge-proposal metric" value={hatchLayers.proposal.metric} onChange={(e) => ig.setHatchLayers({ proposal: { ...hatchLayers.proposal, metric: e.target.value as HatchMetric } })}><option value="srd">SRD</option><option value="srd_phi">SRD/√φ</option><option value="delta_log2fc">delta log2FC</option></select>
          {hatchLayers.proposal.metric === 'delta_log2fc' && <select aria-label="Merge-proposal estimator" value={hatchLayers.proposal.estimator} onChange={(e) => ig.setHatchLayers({ proposal: { ...hatchLayers.proposal, estimator: e.target.value as Estimator } })}><option value="mean">mean</option><option value="median">median</option><option value="iqr_mean">two-sided IQR mean</option></select>}
          <label>threshold<input aria-label="Merge-proposal threshold" type="number" step={0.1} value={hatchLayers.proposal.threshold} onChange={(e) => ig.setHatchLayers({ proposal: { ...hatchLayers.proposal, threshold: +e.target.value } })} style={{ width: 64 }} /></label>
          <label className="chk-inline"><input type="checkbox" checked={hatchLayers.ambiguous.on} onChange={(e) => ig.setHatchLayers({ ambiguous: { ...hatchLayers.ambiguous, on: e.target.checked } })} />ambiguous flank preference (×)</label>
          <select aria-label="Ambiguous estimator" value={hatchLayers.ambiguous.estimator} onChange={(e) => ig.setHatchLayers({ ambiguous: { ...hatchLayers.ambiguous, estimator: e.target.value as Estimator } })}><option value="mean">mean</option><option value="median">median</option><option value="iqr_mean">two-sided IQR mean</option></select>
          {ig.hatchScoresError && <span className="notice small">{ig.hatchScoresError}</span>}
        </div>
```

Place this immediately after the existing `</div>` that closes the layer-toggle `controls-row` (the one containing `rows:` and `layers:`), i.e. right before the `<div className="seg-legend">` block.

- [ ] **Step 6: Extend the legend with the new layers' meaning**

Immediately after the existing legend lines for `seg.labels.proposal`/`seg.labels.consistency`, add:
```tsx
          {hatchLayers.transition.on && <span>new: hatch / where left/right {hatchLayers.transition.metric === 'delta_log2fc' ? `delta log2FC (${hatchLayers.transition.estimator})` : hatchLayers.transition.metric} scores have opposite signs (descriptive, not a merge decision)</span>}
          {hatchLayers.proposal.on && <span>new: hatch × where at least one flank's |{hatchLayers.proposal.metric === 'delta_log2fc' ? `delta log2FC (${hatchLayers.proposal.estimator})` : hatchLayers.proposal.metric}| &lt; {hatchLayers.proposal.threshold}</span>}
          {hatchLayers.ambiguous.on && <span>new: hatch × where the weaker SRD/√φ flank differs from the weaker delta log2FC ({hatchLayers.ambiguous.estimator}) flank</span>}
          <span><i className="ln" style={{ borderColor: '#b91c1c' }} />dotted red = a genomic gap sits at that flank boundary</span>
```

- [ ] **Step 7: Build**

Run: `cd scientific-explorer/web && npm run build`
Expected: build succeeds with no TypeScript errors

- [ ] **Step 8: Commit**

```bash
git add web/src/components/IntegrationSegmentationView.tsx
git commit -m "feat(hatch-merge): wire three new hatch layers into Chromosome tracks"
```

---

### Task 9: Exploratory merge tool panel

**Files:**
- Create: `web/src/components/HatchMergeView.tsx`
- Modify: `web/src/components/IntegrationSegmentationView.tsx` (mount the new panel)

**Interfaces:**
- Consumes: `runHatchMerge`, `MergeResult`, `MergeStep`, `HatchMetric`, `Estimator`, `hatchShapes`, `scoreForMetric` (indirectly, via the returned rows if needed) from `web/src/hatchMerge.ts`; `useIntegration` from `web/src/integrationStore.ts`; `useStore` from `web/src/store.ts`; `Plot` from `./Plot` (same convention as `AlternativeFlankView.tsx`, which does not use `DeferredPlot` since it is not part of the always-rendered tracks list — confirm against Task constraint "keep DeferredPlot" by using `DeferredPlot` here too, since this panel's plots are also off-screen until scrolled to).
- Produces: `export default function HatchMergeView(): JSX.Element` — a self-contained panel with no props, reading `ig.seg.chrom`/`ig.seg.sources` for its default chromosome/source like `AlternativeFlankView` does.

- [ ] **Step 1: Implement the component**

```tsx
// web/src/components/HatchMergeView.tsx
import { useRef, useState } from 'react'
import DeferredPlot from './DeferredPlot'
import { useIntegration } from '../integrationStore'
import { useStore } from '../store'
import { runHatchMerge, hatchShapes, type Estimator, type HatchMetric, type MergeResult } from '../hatchMerge'

const METRIC_NAMES: Record<HatchMetric, string> = { srd: 'SRD', srd_phi: 'SRD/√φ', delta_log2fc: 'delta log2FC' }
function download(text: string, name: string, type: string) {
  const url = URL.createObjectURL(new Blob([text], { type })); const a = document.createElement('a')
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000)
}

export default function HatchMergeView() {
  const ig = useIntegration(), meta = ig.meta!
  const ex = useStore((s) => s.meta)!, presentation = useStore((s) => s.presentation)
  const [source, setSource] = useState('cluster4')
  const [metric, setMetric] = useState<HatchMetric>('srd_phi')
  const [estimator, setEstimator] = useState<Estimator>('mean')
  const [threshold, setThreshold] = useState(3.3)
  const [vetoTransition, setVetoTransition] = useState(false)
  const [vetoAmbiguous, setVetoAmbiguous] = useState(false)
  const [smallMaxBins, setSmallMaxBins] = useState<number | ''>('')
  const [allowGapCrossing, setAllowGapCrossing] = useState(false)
  const [result, setResult] = useState<MergeResult | null>(null)
  const [stepIdx, setStepIdx] = useState(0)
  const [busy, setBusy] = useState(false)
  const [error, setError] = useState('')
  const token = useRef(0)
  const chrom = ig.seg.chrom
  const request = { chrom, source, metric, threshold, estimator, veto_transition: vetoTransition, veto_ambiguous: vetoAmbiguous,
                    small_max_bins: smallMaxBins === '' ? undefined : smallMaxBins, allow_gap_crossing: allowGapCrossing }
  const signature = JSON.stringify(request)
  const [frozen, setFrozen] = useState('')
  const stale = !!result && frozen !== signature

  const run = async () => {
    const id = ++token.current; setBusy(true); setError('')
    try {
      const data = await runHatchMerge(request)
      if (token.current !== id) return
      setResult(data); setFrozen(signature); setStepIdx(data.steps.length - 1)
    } catch (e: any) { if (token.current === id) setError(String(e.message)) }
    finally { if (token.current === id) setBusy(false) }
  }
  const reset = () => setStepIdx(0)
  const undo = () => setStepIdx((i) => Math.max(0, i - 1))
  const step = result?.steps[stepIdx]

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
    <h3>Exploratory adjacent-segment merge</h3>
    <p className="lead">A new, independent exploratory rule: repeatedly merges the eligible adjacent-segment boundary with the smallest absolute score until none remain eligible. This is not the original bootstrap/chain-guard merge prototype (see its results and controls further below) and it does not recompute CN calls on the merged boundaries.</p>
    <div className="controls-row">
      <label>Source<select aria-label="Merge source" value={source} onChange={(e) => setSource(e.target.value)}>{meta.sources.map((s) => <option key={s}>{s}</option>)}</select></label>
      <label>Chromosome<select aria-label="Merge chromosome" value={chrom} onChange={(e) => ig.setSeg({ chrom: e.target.value })}>{ex.chromosomes.map((c) => <option key={c.name}>{c.name}</option>)}</select></label>
      <label>Merge metric<select aria-label="Merge metric" value={metric} onChange={(e) => setMetric(e.target.value as HatchMetric)}><option value="srd">SRD</option><option value="srd_phi">SRD/√φ</option><option value="delta_log2fc">delta log2FC</option></select></label>
      {metric === 'delta_log2fc' && <label>Estimator<select aria-label="Merge estimator" value={estimator} onChange={(e) => setEstimator(e.target.value as Estimator)}><option value="mean">mean</option><option value="median">median</option><option value="iqr_mean">two-sided IQR mean</option></select></label>}
      <label>Threshold<input aria-label="Merge threshold" type="number" step={0.1} value={threshold} onChange={(e) => setThreshold(+e.target.value)} style={{ width: 72 }} /></label>
    </div>
    <div className="controls-row">
      <label className="chk-inline"><input type="checkbox" checked={vetoTransition} onChange={(e) => setVetoTransition(e.target.checked)} />veto boundaries touching a transition-zone segment</label>
      <label className="chk-inline"><input type="checkbox" checked={vetoAmbiguous} onChange={(e) => setVetoAmbiguous(e.target.checked)} />veto boundaries touching an ambiguous-flank segment</label>
      <label className="chk-inline"><input type="checkbox" checked={allowGapCrossing} onChange={(e) => setAllowGapCrossing(e.target.checked)} />allow merging across a genomic gap (off by default)</label>
      <label>Restrict eligibility to boundaries touching a segment ≤ <input aria-label="Small-segment restriction" type="number" min={1} placeholder="off" value={smallMaxBins} onChange={(e) => setSmallMaxBins(e.target.value === '' ? '' : +e.target.value)} style={{ width: 64 }} /> bins</label>
      <button className="btn primary" disabled={busy} onClick={run}>{busy ? 'Merging…' : 'Run merge'}</button>
    </div>
    {error && <p role="alert">{error}</p>}
    {result && step && <>
      <p className={stale ? 'notice' : 'muted'} role="status">{stale ? 'Settings changed; showing the previous run. ' : ''}{result.source} {result.chrom}: {result.original.length} original segments → {result.final.length} after {result.steps.length - 1} merge(s). Metric {METRIC_NAMES[metric]}{metric === 'delta_log2fc' ? ` (${estimator})` : ''}, threshold {threshold}.</p>
      <div className="controls-row">
        <label>Step<input aria-label="Merge step" type="range" min={0} max={result.steps.length - 1} value={stepIdx} onChange={(e) => setStepIdx(+e.target.value)} /></label>
        <span className="muted small">step {stepIdx} of {result.steps.length - 1}{step.removed_boundary !== null ? ` — removed boundary at bin ${step.removed_boundary}, score ${step.score?.toFixed(3)}` : ' — original segmentation'}</span>
        <button className="btn-sm" onClick={reset} disabled={stepIdx === 0}>Reset to original</button>
        <button className="btn-sm" onClick={undo} disabled={stepIdx === 0}>Undo one step</button>
        <button className="btn-sm" onClick={() => setStepIdx(result.steps.length - 1)} disabled={stepIdx === result.steps.length - 1}>Jump to final</button>
        <button className="btn-sm" onClick={() => download(JSON.stringify(result, null, 2), `hatch-merge-${result.chrom}-${result.source}.json`, 'application/json')}>Download JSON</button>
        <button className="btn-sm" onClick={() => {
          const rows = [['step', 'source', 'chromosome', 'merged_start', 'merged_end', 'removed_boundary', 'metric', 'estimator', 'score', 'threshold', 'phi'],
            ...result.steps.filter((s) => s.removed_boundary !== null).map((s) => [s.step, result.source, result.chrom, s.merged_interval?.[0], s.merged_interval?.[1], s.removed_boundary, metric, metric === 'delta_log2fc' ? estimator : '', s.score, threshold, s.phi ?? ''])]
          download(rows.map((r) => r.join(',')).join('\n') + '\n', `hatch-merge-history-${result.chrom}-${result.source}.csv`, 'text/csv')
        }}>Download merge-history CSV</button>
      </div>
      <DeferredPlot data={[segmentTrace(result.original, result.start_bp, result.end_bp, '#94a3b8', 'Original')]}
        layout={{ ...layout, title: { text: 'Original segmentation', font: { size: 12 } }, xaxis: { title: `${result.chrom} position (Mb)` }, yaxis: { visible: false }, shapes: shapesForSegments(result.original, result.start_bp, result.end_bp) }} />
      <DeferredPlot data={[segmentTrace(step.segments, result.start_bp, result.end_bp, '#96491d', `Step ${stepIdx}`)]}
        layout={{ ...layout, title: { text: `Segmentation at step ${stepIdx} (${step.segments.length} segments)`, font: { size: 12 } }, xaxis: { title: `${result.chrom} position (Mb)` }, yaxis: { visible: false }, shapes: [...shapesForSegments(step.segments, result.start_bp, result.end_bp), ...(step.merged_interval ? hatchShapes(result.start_bp[step.merged_interval[0]] / 1e6, result.end_bp[step.merged_interval[1] - 1] / 1e6, '', '/', '#111827') : [])] }} />
      <div className="table-scroll"><table className="tbl"><thead><tr><th>step</th><th>source</th><th>chromosome</th><th>merged interval (bins)</th><th>removed boundary</th><th>metric</th><th>estimator</th><th>score</th><th>threshold</th><th>φ</th></tr></thead><tbody>
        {result.steps.filter((s) => s.removed_boundary !== null).map((s) => <tr key={s.step} className={s.step === stepIdx ? 'sel' : 'clickable'} onClick={() => setStepIdx(s.step)}>
          <td>{s.step}</td><td>{result.source}</td><td>{result.chrom}</td><td>[{s.merged_interval?.[0]}, {s.merged_interval?.[1]})</td><td>{s.removed_boundary}</td>
          <td>{METRIC_NAMES[metric]}</td><td>{metric === 'delta_log2fc' ? estimator : '–'}</td><td>{s.score?.toFixed(3)}</td><td>{threshold}</td><td>{s.phi != null ? s.phi.toFixed(3) : '–'}</td>
        </tr>)}
      </tbody></table></div>
      <p className="small muted">Segment summaries and hatch classifications for this step's segmentation use the "new exploratory hatch layers" controls above once you tick a layer — those recompute from live boundaries and are not tied to this panel's step selection. CN calls are not recomputed on merged boundaries; any CN shown elsewhere in this tab still reflects the original accepted-breakpoint segmentation.</p>
      <details><summary>Merge run provenance and settings</summary><pre style={{ whiteSpace: 'pre-wrap', overflowWrap: 'anywhere' }}>{JSON.stringify({ settings: result.settings, provenance: result.provenance }, null, 2)}</pre></details>
    </>}
  </section>
}
```

- [ ] **Step 2: Mount the panel**

In `web/src/components/IntegrationSegmentationView.tsx`, add the import:
```typescript
import HatchMergeView from './HatchMergeView'
```
And mount it directly after the closing `</div>` of the `<div className="panel" id="cluster-tracks">` block (i.e. as a sibling panel right after Chromosome tracks, before the "Small-segment diagnostics" panel):
```tsx
      <HatchMergeView />
```

- [ ] **Step 3: Build**

Run: `cd scientific-explorer/web && npm run build`
Expected: build succeeds with no TypeScript errors

- [ ] **Step 4: Commit**

```bash
git add web/src/components/HatchMergeView.tsx web/src/components/IntegrationSegmentationView.tsx
git commit -m "feat(hatch-merge): exploratory merge tool panel with history, undo/reset and exports"
```

---

### Task 10: Documentation and handoff

**Files:**
- Modify: `scientific-explorer/START_HERE.md`
- Modify: `scientific-explorer/README.md`
- Create: `scientific-explorer/validation/hatch-merge/HANDOFF.md`

**Interfaces:** none (documentation only).

- [ ] **Step 1: Add a `START_HERE.md` section** (near the existing "Integration · segmentation & CN" paragraph)

```markdown
### New: exploratory hatching and adjacent-segment merging

Three new hatch layers in Chromosome tracks are computed live from the current accepted-breakpoint segmentation and pseudobulk profile, for any chromosome (not only the seven with archived diagnostic tables): **transition zone** (/, opposite-sign flank scores), **merge proposal** (×, at least one flank's |score| below its own threshold), and **ambiguous flank preference** (×, the weaker SRD/√φ flank disagrees with the weaker delta-log2FC flank). Each has its own metric/estimator/threshold controls and is clearly labelled "new" to distinguish it from the archived "flank-consistency hatch"/"SRD-screen proposal" layers above, which still reflect only the original prototype's saved tables.

Below Chromosome tracks, **Exploratory adjacent-segment merge** runs a separate, new rule: it repeatedly merges the eligible adjacent boundary with the smallest absolute score (strictly below your threshold) until none remain eligible, recomputing every score — including the pooled dispersion φ, which is shared across the whole chromosome's segmentation — after each merge. Optional, off-by-default vetoes can exclude boundaries touching a transition or ambiguous segment; an optional, off-by-default size restriction can limit eligibility to boundaries touching a small segment; genomic gaps are never crossed unless you explicitly allow it. Use the step slider (or Reset/Undo) to inspect any intermediate segmentation, and download the full run (JSON) or just the merge-history table (CSV). This is an independent exploratory rule, not the original bootstrap/chain-guard merge prototype, and it does not recompute CN calls on the merged boundaries.
```

- [ ] **Step 2: Add a `README.md` section** (near the existing integration API/route documentation, if any; otherwise near "Integration tabs")

```markdown
### Hatching and merge routes

`POST /api/integration/hatch-scores` `{chrom, source}` returns per-internal-segment SRD, SRD/√φ and delta-log2FC (mean/median/two-sided-IQR-mean) scores against both flanks, plus the chromosome's pooled φ and per-side genomic-gap flags, computed live from `server/hatch_merge.py` — genome-wide, not limited to the archived pilot chromosomes.

`POST /api/integration/hatch-merge` `{chrom, source, metric, threshold, estimator?, veto_transition?, veto_ambiguous?, small_max_bins?, allow_gap_crossing?}` runs the new exploratory greedy merge and returns the full step history (original → final), so the frontend can undo/reset/inspect any intermediate step purely by indexing into the response — no server-side session state. See `server/hatch_merge.py` for the numerical rules and `tests/test_hatch_merge.py` for their test coverage.
```

- [ ] **Step 3: Run the full verification checklist and record it**

Run, from `scientific-explorer/`:
```bash
python3 -m pytest -q tests/test_hatch_merge.py tests/test_experiments.py tests/test_segment_estimators.py tests/test_integration_review.py
npm run build --prefix web
web/node_modules/.bin/esbuild tests/hatch_merge_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/hatch_merge_helpers.mjs && node tests/_bundle/hatch_merge_helpers.mjs
```
Then, in Chrome against a running server (`python3 run.py --port 8766`), walk the task's browser checklist: all hatch toggles and metric selectors work; labels explain both flank scores and why a segment was classified; a merge removes the intended boundary and changes subsequent comparisons; original/merged views, history, reset and exports agree; session restoration and stale-result handling work (toggle a layer, save session, reload, confirm it restores; change a merge setting after running and confirm the stale notice appears); scrolling, chromosome changes and controls remain responsive; no browser console errors.

Write `scientific-explorer/validation/hatch-merge/HANDOFF.md` recording exactly what passed (paste command output/counts), any unresolved limitations (e.g., no CN recomputation on merged boundaries; large chromosomes' merge-run latency if observed; anything found not to reproduce), and a short "how to use the new controls" section addressed to the next person opening the tab. Do not claim the rule is statistically calibrated or biologically validated — state plainly that it is descriptive/exploratory.

- [ ] **Step 4: Commit**

```bash
git add scientific-explorer/START_HERE.md scientific-explorer/README.md scientific-explorer/validation/hatch-merge/HANDOFF.md
git commit -m "docs(hatch-merge): document new controls and record verification handoff"
```

---

## Self-Review Notes (for whoever executes this plan)

- **Spec coverage:** §1 (flank definitions) → Tasks 1-2; §2 (three hatch layers) → Tasks 3-4 (backend) + Task 8 (frontend); §3 (merge rule) → Task 4 (backend) + Task 9 (frontend); §4 (show/preserve results) → Task 9 (history/undo/reset/export) + Task 7 (session persistence of settings); §5 (implementation constraints) → addressed throughout (small backend module, reused formulas, `DeferredPlot` reuse in Task 9, no restriction to pilot chromosomes in Task 4/5); §6 (verification) → Tasks 1-6's tests plus Task 10's checklist and handoff doc.
- **CN recomputation** is explicitly out of scope (Global Constraints and called out again in Task 9's UI copy) since the task says not to silently reuse original CN calls as if recomputed, and implementing/verifying real recomputation is a separably large scope the task doc does not require.
- **Type consistency check:** `HatchScoreRow` field names (`srd_left`, `delta_mean_left`, etc.) match exactly between `server/hatch_merge.py`'s `hatch_scores` row dict (Task 4), the FastAPI response (Task 5, passthrough), and the TypeScript interface (Task 6) — verify this literally at Task 6 Step 3 before moving on, since a mismatch here would silently produce `undefined` scores in the UI rather than a type error (the fetch response is untyped JSON).
