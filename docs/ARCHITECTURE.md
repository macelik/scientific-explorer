# Architecture

A FastAPI backend (`server/`) serves both a JSON API and the built React
frontend (`web/dist/`). The React/TypeScript/Vite frontend (`web/src/`) uses
Zustand for state and Plotly for genomic tracks and plots. `reference/` holds
archival, verified-reproducible prior-research tables; `validation/` holds
test evidence and screenshots.

## Backend (`server/`)

| Module | Responsibility |
|---|---|
| `run.py` (repo root) | Entry point: CLI args, env var resolution, launches the FastAPI app. |
| `server/app.py` | Assembles the FastAPI app and registers every route module. **Route registration order matters**: API routes must be registered before the SPA static-file catch-all, or they silently receive HTML instead of JSON. |
| `server/config.py` / `server/data.py` | Resolves the data root and its overrides (`PYEPI_SCIENTIFIC_*` env vars), loads the cohort. |
| `server/science.py` / `server/pyepi_adapter.py` | The adapter around the external `pyEpiAneufinder` package — see [LICENSING.md](LICENSING.md). Registers the package directory without executing its `__init__.py` (which eagerly imports plotting/postprocessing modules) and imports the scientific modules directly. Every statistic here is a direct, unmodified call into `pyEpiAneufinder`. |
| `server/experiments.py` / `server/experiment_api.py` | X experiments: region edits (multiply/add/replace/simulate/extend/duplicate), width series, event-edge diagnostics, recursion-depth inspection (k = 1–10), the four CN comparison rows. |
| `server/segment_estimators.py` | Manual segmentation's alternative segment-mean estimators (arithmetic, upper-tail IQR, two-sided IQR). |
| `server/integration.py` / `server/integration_api.py` | Multimodal clustering/pseudobulk data adapter for the Integration tabs — loads the depth + haplotype integration story's exported artifacts read-only (nothing is re-clustered, re-integrated, or re-segmented here), with background indexing/caching under `cache/`. |
| `server/flank_view.py` | The on-demand "alternative flank view" and paired-cell bootstrap: descriptive flank-comparison recomputation independent of the archived SRD-screen tables. |
| `server/hatch_merge.py` | Flank scoring for the exploratory tools: SRD, SRD/√φ, Δlog2FC, pooled φ, gap flags, and `score_rows()` / `hatch_scores()`. Independent of `flank_view.py` and the archived merge prototype. |
| `server/merge_rules.py` / `server/hatch_merge_api.py` | Rule 0 (criteria checklist, all/any, sign switch, the `mergeSignsOnly` last phase; direction toward the criterion furthest below its threshold) plus user-defined merge rules, including Rule 0 copies (`type: "flank_threshold"`): validation, per-segment evaluation (Rule 0 first, then OR of ANDs), and the one-merge-per-step loop, ordered by (phase, rule, score, start), with φ recomputed every step and a `phase` on every step; the `hatch-scores` and `hatch-merge` routes. See [SCIENTIFIC_WORKFLOW.md](SCIENTIFIC_WORKFLOW.md)'s worked example 2. |
| `server/jobs.py` | Bounded worker-process pool for expensive computation. Workers use multiprocessing `spawn` — worker functions must be importable/picklable at module level, not defined as closures. |

## Frontend (`web/src/`)

| Module | Responsibility |
|---|---|
| `store.ts` | Main cohort/cell/experiment state; session serialize/restore. |
| `integrationStore.ts` | Integration/multimodal state — merge-rule and visualization settings persistence, candidate/score caching, cluster selection. |
| `experimentStore.ts` | Per-cell/chromosome experiment drafts and run history. |
| `integration.ts` | Pure helpers: consistency classes, flank picking, hatch-shape rendering geometry. |
| `hatchMerge.ts` | Score/result types and the API client for `hatch-scores` / `hatch-merge`. |
| `mergeRules.ts` | Rule types (including Rule 0 and its display outcome / `rule0FromSeg`), presets, legacy-session migration, and a 1:1 TypeScript port of `server/merge_rules.py`'s evaluator, so hatches update live while rules are edited. `tests/fixtures/merge_rules_parity.json` is checked against both implementations. |
| `trackMetrics.ts` | Live colour metrics for the tracks (smaller-\|value\| flank), migration of archived colour ids in saved sessions, and derivation of candidate-audit rejection reasons from p-values. |
| `api.ts` | Fetch wrapper. |
| `components/` | Major workspaces: `CohortView`, `ExplorerView`, `ManualView`, `ExperimentView`, `IntegrationClustersView`, `IntegrationSegmentationView`, `HatchMergeView`, `MergeRulesEditor`, `AlternativeFlankView`, `PairedFlankBootstrap`, plus shared `Plot`/`DeferredPlot` Plotly wrappers. `DeferredPlot` exists specifically to defer initializing off-screen charts until scrolled into view — see [`validation/integration-freeze/REPORT.md`](../validation/integration-freeze/REPORT.md) for the responsiveness investigation that motivated it. |

## Data flow: one representative request

Browser action (e.g. clicking "Run merge") → HatchMergeView reads the merge rules from the
integration store → `hatchMerge.ts`'s `runHatchMerge()` posts to
`/api/integration/hatch-merge` → `server/hatch_merge_api.py`'s route resolves
the current integration profile → `server/merge_rules.py`'s `run_rule_merge()`
validates the rules (`validate_rules`; the route maps a `ValueError` to a 422)
and computes the full step history → JSON response → the store
updates → `HatchMergeView.tsx` re-renders the step history and plots.

```mermaid
flowchart TD
    Browser[Browser] --> Components[React components]
    Components --> Stores[Zustand stores]
    Stores --> APIClient["API client (api.ts / hatchMerge.ts)"]
    APIClient --> Routes["FastAPI routes (server/*_api.py)"]
    Routes --> Science["Scientific modules\n(science.py, flank_view.py, hatch_merge.py, merge_rules.py)"]
    Science --> PyEpi["pyEpiAneufinder (external, not vendored)"]
    Science --> Reference["reference/ (archival tables)"]
    Science --> Cache["cache/ (integration index)"]
```

## Tests and validation

`tests/` holds unit/parity tests (pytest), a Python/TypeScript numerical
parity check, and Playwright browser workflow tests. `validation/` holds the
evidence and reports those tests produce. See
[`tests/README.md`](../tests/README.md) and
[`validation/README.md`](../validation/README.md) for detail, and
[VALIDATION.md](VALIDATION.md) for the portfolio-level summary.
