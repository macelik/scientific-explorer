# tests/

Three kinds of tests live here:

- **Unit / parity tests** (`test_*.py`), run via `pytest`.
- **Python/TypeScript numerical parity checks** (`parity_reference.py`,
  `parity_browser.mjs`), comparing browser-executed TypeScript statistics
  against the Python reference implementation.
- **Browser workflow tests** (`browser_*.cjs`), driving the actual running
  app with Playwright.

## Files

- `test_experiments.py`, `test_segment_estimators.py`, `test_jobs.py`,
  `test_hatch_merge.py`, `test_integration_review.py` — pytest unit/route
  tests.
- `validate.py`, `validate_integration_review.py` — standalone API
  validation scripts (run against a live server).
- `parity_reference.py` / `parity_browser.mjs` — Python-vs-TypeScript
  numerical parity.
- `experiment_overlays.ts`, `hatch_merge_helpers.ts`, `integration_helpers.ts`
  — TypeScript logic ported from Python, bundled with esbuild and run under
  Node so it can be tested without a browser.
- `browser_experiments.cjs`, `browser_workflows.cjs`,
  `browser_experiment_depth.cjs`, `browser_integration_rendering.cjs`,
  `browser_integration_review.cjs`, `browser_integration_freeze.cjs`,
  `browser_hatch_merge_walkthrough.cjs`, `browser_restore.cjs` — Playwright
  workflow tests, one per feature area.

## Prerequisites

The app must be running on port 8766 for API/browser tests
(`python3 run.py`, from the repo root). Browser tests need Playwright and
Chromium; set `PLAYWRIGHT_MODULE` to an installed Playwright module
directory and `CHROMIUM_PATH` if using a separately installed Chromium,
otherwise normal Node module/browser discovery is used.

## Minimal sanity check after a small change

```bash
python3 -m pytest -q tests/test_hatch_merge.py tests/test_experiments.py tests/test_segment_estimators.py tests/test_integration_review.py
cd web && npm run build
```
