# server/

The FastAPI backend: serves the JSON API and the built frontend
(`web/dist/`).

- `app.py` — assembles the app and registers every route module. Route
  registration order matters: API routes must be registered before the SPA
  static-file catch-all.
- `config.py` / `data.py` — data root resolution (`PYEPI_SCIENTIFIC_*`
  overrides), cohort loading.
- `science.py` / `pyepi_adapter.py` — the adapter around the external,
  non-vendored `pyEpiAneufinder` package. See
  [`../docs/LICENSING.md`](../docs/LICENSING.md).
- `experiments.py` / `experiment_api.py` — X experiments: region edits,
  width series, event-edge diagnostics, recursion-depth inspection.
- `segment_estimators.py` — manual segmentation's alternative segment-mean
  estimators.
- `integration.py` / `integration_api.py` — multimodal clustering/pseudobulk
  data adapter for the Integration tabs.
- `flank_view.py` — the alternative flank view and paired-cell bootstrap.
- `hatch_merge.py` / `hatch_merge_api.py` — the exploratory hatch-layer
  scoring and adjacent-segment merge module.
- `jobs.py` — bounded worker-process pool (multiprocessing `spawn`; worker
  functions must be importable/picklable).

## Where should new scientific logic live?

A new pure-function module beside the existing ones (`flank_view.py`,
`hatch_merge.py`), never inline in a route handler, with its own
`*_api.py` module for the FastAPI wiring. See
[`../docs/EXTENDING.md`](../docs/EXTENDING.md) for the full workflow.
