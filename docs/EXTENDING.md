# Extending the explorer

If you (or a coding agent) want to add a new scientific diagnostic,
experiment, scoring method, visualization, or integration panel, this is how
to do it without silently changing what the tool already claims to do.

## The repeatable workflow

1. **Define the scientific question.** What hypothesis does this feature let
   someone test?
2. **Identify what the feature changes:** data, statistic, candidate
   selection, significance testing, segmentation, CN assignment, or
   visualization only. Most new panels in this repository are
   visualization-only over an existing statistic; changing the underlying
   statistic is rare and should be flagged loudly in the PR (see
   [CONTRIBUTING.md](../CONTRIBUTING.md)).
3. **Reuse existing scientific functions where possible.**
   `server/flank_view.py`'s `summarize()` and `server/hatch_merge.py`'s
   ported SRD/pooled-φ functions are both used across otherwise-independent
   modules rather than each reimplementing the statistic.
4. **Keep exploratory alternatives explicitly labelled.** See how the
   hatch-merge layers are labelled "new" in the UI and kept visually and
   textually distinct from the archived consistency/proposal layers in
   `IntegrationSegmentationView.tsx`.
5. **Add the backend contract:** a Pydantic request model and a
   `register_X(app, get_integration)` function following
   `server/hatch_merge_api.py`'s pattern, registered in `server/app.py`
   *before* the SPA catch-all route (see below for why).
6. **Add frontend state + visualization:** a Zustand store slice (see
   `integrationStore.ts`'s `hatchLayers`/`ensureHatchScores` pattern for
   settings + fetch/cache) and a component under `web/src/components/`.
7. **Add unit/parity tests:** a `tests/test_*.py` file mirroring
   `tests/test_hatch_merge.py` (pure-function tests first, route tests with
   a fake store second); a `tests/*_helpers.ts` file if client-side
   classification logic is ported from Python, built with esbuild and run
   under Node the same way `tests/hatch_merge_helpers.ts` is.
8. **Add one browser workflow test:** a `tests/browser_*.cjs` script using
   the repository's existing Playwright harness pattern.
9. **Add validation evidence:** a `validation/<feature>/` directory with a
   `REPORT.md` or `HANDOFF.md` stating what was checked and how, following
   [`validation/hatch-merge/HANDOFF.md`](../validation/hatch-merge/HANDOFF.md)'s
   structure.
10. **Document assumptions and limitations** in the relevant `docs/*.md` file
    and in the UI's own copy. An undocumented limitation is, for this
    repository's purposes, not finished.

## How to add an API route

Create a small module (e.g. `server/my_feature.py`) with pure functions and
no FastAPI/framework imports. Create `server/my_feature_api.py` with Pydantic
request models and a `register_my_feature(app, get_integration)` function.
Add the import and registration call in `server/app.py` immediately after the
existing

```python
from .hatch_merge_api import register_hatch_merge
register_hatch_merge(app, lambda: integration)
```

and before the `DIST`/`StaticFiles` block. Ordering matters: API routes
registered after the SPA catch-all silently receive HTML instead of JSON —
this is a real, previously-hit failure mode, not a hypothetical.

## How to add a frontend workspace/panel

Add a new tab constant to `web/src/store.ts`'s tab state (or the relevant
store, if it's an Integration sub-panel); add a component under
`web/src/components/`; mount it in `web/src/App.tsx`'s conditional-view
block. If it needs settings that should survive a session reload, add them
to the relevant store's `serialize()`/`restore()` functions.

## Session/export conventions

Settings that should persist across a saved/reloaded session belong in a
store's `serialize()`/`restore()`, merged onto sensible defaults on restore
so an older saved session without the new field doesn't crash —
`integrationStore.ts`'s hatch-layer restore does exactly this
(`{ ...defaultHatchLayers, ...(p.hatchLayers || {}) }`). Generated *results*
(not settings) are typically kept local to the component instead —
`AlternativeFlankView.tsx` keeps its generated result in local `useState`
plus a frozen-settings signature comparison, not in a store.

## Avoiding silent mutation of original data

Experiments and edits always start from original X or boundaries and never
accumulate. Every edit path returns both the edited result and enough
provenance (parameters, seed, source hashes) to reconstruct it —
`server/experiments.py`'s `apply_edit` copies its source array
(`np.asarray(source, ...).copy()`) rather than mutating the input in place.

## Stale-result handling

Any component with an expensive "Run"/"Generate" action should compute a
signature from its current settings, compare it to the signature frozen at
the last run, and show an explicit "settings changed, showing the previous
run" notice rather than silently displaying results that no longer match the
visible controls. `HatchMergeView.tsx`'s `stale` check is the pattern to
copy. This is exactly the class of bug a prior whole-branch review in this
repository caught and fixed before merge — not a hypothetical concern.

## Seeds and provenance

Every stochastic operation (bootstrap resampling, empirical class sampling)
takes an explicit, user-visible, editable seed. Every generated result's
JSON export includes the seed, parameters, and a dataset/implementation
identity hash. `flank_view.py`'s `provenance` dict
(`implementation_sha256`, `dataset`, ...) is the template to copy.

## When something belongs in `reference/` rather than runtime code

`reference/` holds copied, byte-for-byte-verified snapshots of external or
prior-session scripts kept for reproducibility auditing. It is read-only
archival material: no file under `server/` or `web/src/` imports from it —
the only references are docstring citations (e.g. `server/hatch_merge.py`
notes its SRD formulas were "ported from
`reference/integration-prototype/src/srd_effect/srd_pruning.py`", meaning
behaviorally reproduced and tested against, never imported at runtime). A
new feature's own code never goes there.

## Worked example: the hatch-merge feature

The feature started from a scientific question (worked example 2 in
[SCIENTIFIC_WORKFLOW.md](SCIENTIFIC_WORKFLOW.md)): is a small segment a real
boundary or an artifact of weak flank support? It touched, in order:
`server/hatch_merge.py` (pure classification + merge functions) →
`server/hatch_merge_api.py` (Pydantic layer models, the `register_hatch_merge`
route) → `web/src/hatchMerge.ts` (types, API client, ported classifiers) →
`web/src/integrationStore.ts` (hatch-layer settings state) →
`web/src/components/IntegrationSegmentationView.tsx` +
`web/src/components/HatchMergeView.tsx` (the hatch-layer controls and the
merge panel) → `tests/test_hatch_merge.py` +
`tests/hatch_merge_helpers.ts` + `tests/browser_hatch_merge_walkthrough.cjs`
(unit, parity, and browser coverage) →
`validation/hatch-merge/HANDOFF.md` (what was checked). The feature went
through an individual-task review process and then a whole-branch review
that caught and fixed six real cross-cutting issues before merge — a
concrete demonstration that step 7 above ("add tests") is not a formality.
Its eligibility model was later redesigned from a single-metric-with-vetoes
scheme to a union-of-independently-configured-layers scheme, on the same
files, which is why several of the file references above describe the
current, not the original, shape of that feature.
