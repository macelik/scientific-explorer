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
  `test_hatch_merge.py`, `test_merge_rules.py`, `test_integration_review.py` — pytest unit/route
  tests.
- `validate.py`, `validate_integration_review.py` — standalone API
  validation scripts (run against a live server).
- `parity_reference.py` / `parity_browser.mjs` — Python-vs-TypeScript
  numerical parity.
- `experiment_overlays.ts`, `hatch_merge_helpers.ts`,
  `merge_rules_helpers.ts`, `cn_recall_helpers.ts` (Changed-mode rows and export TSV writers), `track_helpers.ts` (live track metrics, session migration, audit reasons), `integration_helpers.ts` — TypeScript logic
  ported from Python, bundled with esbuild and run under Node so it can be
  tested without a browser. For example:
  `web/node_modules/.bin/esbuild tests/merge_rules_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rules_helpers.mjs && node tests/_bundle/merge_rules_helpers.mjs`
  Same pattern for the re-call helpers:
  `web/node_modules/.bin/esbuild tests/cn_recall_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/cn_recall_helpers.mjs && node tests/_bundle/cn_recall_helpers.mjs`
- `compare_archived_live.py` — compares the archived flank/screening tables
  with live `hatch-scores` on the pilot chromosomes (needs a running server;
  set `BASE_URL`, e.g. `BASE_URL=http://127.0.0.1:8766 python3 tests/compare_archived_live.py`).
  The recorded result is in `validation/rule0/README.md`.
- `track_helpers.ts` is bundled the same way as the other helpers:
  `web/node_modules/.bin/esbuild tests/track_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/track_helpers.mjs && node tests/_bundle/track_helpers.mjs`
- `fixtures/merge_rules_parity.json` — shared by the Python and TypeScript
  merge-rule evaluator tests.
- `browser_experiments.cjs`, `browser_workflows.cjs`,
  `browser_experiment_depth.cjs`, `browser_integration_rendering.cjs`,
  `browser_integration_review.cjs`, `browser_integration_freeze.cjs`,
  `browser_hatch_merge_walkthrough.cjs` (set `BASE_URL` to target a server other than http://127.0.0.1:8766, e.g. a private `python3 run.py --port 8799`; it covers the sign switch, the last-phase option and phase-2 history on chr7, Rule 0 copies, and the absence of the removed flank-conflict and breakpoint layers), `browser_cn_recall_walkthrough.cjs` (genome merge + CN re-call: button and status stages, Before/After/Changed on both karyograms, impact summary, `merged_breakpoints.tsv` export and the stale notice at τ 0.5 / ref 8 / any; needs a free port such as a private `python3 run.py --port 8799`, takes about a minute on a cold server, and writes `validation/cn-recall/`), `browser_restore.cjs` — Playwright
  workflow tests, one per feature area.

## Prerequisites

The app must be running on port 8766 for API/browser tests
(`python3 run.py`, from the repo root). Browser tests need Playwright and
Chromium; set `PLAYWRIGHT_MODULE` to an installed Playwright module
directory and `CHROMIUM_PATH` if using a separately installed Chromium,
otherwise normal Node module/browser discovery is used.

## Minimal sanity check after a small change

```bash
python3 -m pytest -q tests/test_merge_rules.py tests/test_hatch_merge.py tests/test_experiments.py tests/test_segment_estimators.py tests/test_integration_review.py
cd web && npm run build
```
