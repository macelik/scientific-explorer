# web/

React 18 + TypeScript + Vite frontend, using Zustand for state and
Plotly.js for genomic tracks and plots.

## Stores (`src/`)

- `store.ts` — main cohort/cell/experiment state; session serialize/restore.
- `integrationStore.ts` — Integration/multimodal state, hatch-layer settings,
  candidate/score caching.
- `experimentStore.ts` — per-cell/chromosome experiment drafts and history.

## Components (`src/components/`)

Major workspaces: `CohortView`, `ExplorerView`, `ManualView`,
`ExperimentView`, `IntegrationClustersView`, `IntegrationSegmentationView`,
`HatchMergeView`, `AlternativeFlankView`, `PairedFlankBootstrap`.

`Plot.tsx` / `DeferredPlot.tsx` are the shared Plotly wrappers.
**`DeferredPlot` defers initialization until scrolled into view** — use it
for any chart that isn't immediately visible on tab load. This exists
because eagerly initializing every chart on tab load caused a real,
previously-investigated responsiveness problem (see
[`../validation/integration-freeze/REPORT.md`](../validation/integration-freeze/REPORT.md)).

## Adding a new panel

See [`../docs/EXTENDING.md`](../docs/EXTENDING.md)'s "How to add a frontend
workspace/panel" and "Session/export conventions" sections rather than
re-deriving the pattern here.
