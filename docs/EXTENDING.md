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
   merge-rules panel states its rules are user-defined exploratory
   diagnostics, not validated decision rules, and how the tracks panel in
   `IntegrationSegmentationView.tsx` separates "1. Visualization and Rule 0"
   (colouring is display-only; Rule 0 hatches what the merge starts from)
   from "2. Merge rules", while the remaining archived layer ("retained
   reason") stays labelled as archived.
5. **Add the backend contract:** a Pydantic request model and a
   `register_X(app, get_integration)` function following
   `server/hatch_merge_api.py`'s pattern, registered in `server/app.py`
   *before* the SPA catch-all route (see below for why).
6. **Add frontend state + visualization:** a Zustand store slice (see
   `integrationStore.ts`'s `mergeRules`/`ensureHatchScores` pattern for
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
`integrationStore.ts`'s merge-rule restore does exactly this
(`Array.isArray(p.mergeRules) ? p.mergeRules : migrateHatchLayers(p.hatchLayers)`,
converting pre-merge-rules sessions). Generated *results*
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
`server/hatch_merge.py` (flank scoring, `score_rows`) →
`server/merge_rules.py` (rule validation, evaluation, `run_rule_merge`) →
`server/hatch_merge_api.py` (the `rules` request, `register_hatch_merge`) →
`web/src/hatchMerge.ts` (types, API client) + `web/src/mergeRules.ts`
(TypeScript evaluator port, presets, legacy-session migration) →
`web/src/integrationStore.ts` (`mergeRules`, `allowGapCrossing`) →
`web/src/components/MergeRulesEditor.tsx`, `HatchMergeView.tsx` and
`IntegrationSegmentationView.tsx` (rule editing, merge panel, tracks) →
`tests/test_merge_rules.py` + `tests/merge_rules_helpers.ts` +
`tests/fixtures/merge_rules_parity.json` +
`tests/browser_hatch_merge_walkthrough.cjs` (unit, Python/TypeScript
parity, and browser coverage) → `validation/hatch-merge/HANDOFF.md` (what
was checked). History: the feature first shipped as a
single-metric-with-vetoes scheme, was redesigned to a union of three fixed
hatch layers, and was then replaced by user-defined merge rules. It went
through an individual-task review process and then a whole-branch review
that caught and fixed six real cross-cutting issues before the first
version merged, a concrete demonstration that step 7 above ("add tests")
is not a formality. A later revision (Rule 0 revision 2) dropped the flank-conflict switch, added Rule 0 copies as `flank_threshold` user rules and the `mergeSignsOnly` last phase (a `phase` on every step, ordering by `(phase, rule, score, start)`), and fixed the genomic-gap test; the Python engine and the TypeScript port were changed together and are kept in step by the shared parity fixture.
