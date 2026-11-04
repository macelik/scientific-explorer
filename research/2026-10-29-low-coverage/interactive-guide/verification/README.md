# Guide verification · 2026-11-02

The final guide passed **281 browser checks**. The actual results are in
[`browser_report.json`](browser_report.json). Desktop and mobile screenshots
were also inspected visually.

## What was checked

- All 36 entries render every inference step, their worked example and their
  proposed workflow adaptation. Step changes retain keyboard focus.
- Search, empty results, all six category filters, comparison limits/removal,
  URL-fragment restoration, glossary dismissal, note downloads and print styles.
- All four teaching examples, copy/accessibility ambiguity, reference changes,
  keyboard slider input, absent-allele behavior and model-specific controls.
- Count probabilities against 15 independent SciPy distributions, including
  the tail bucket, at absolute tolerance `1e-11`; six Beta posterior intervals,
  including zero observations and boundary counts, at tolerance `0.0004`.
- All 54 region/stratum selections against the embedded source rows.
- All five views at widths 1440, 1024, 768, 390 and 320 pixels, with no horizontal
  page overflow. Table and step-strip scrolling is intentional.
- No JavaScript runtime errors and no HTTP(S) requests from the local page.
- Source-file hashes, the generated HTML hash, no unexpanded template markers,
  and no external runtime resources.

These checks establish rendering, interaction and teaching-calculation
consistency. They do not validate CNV detection, the biological interpretation
of a region, external website availability, or every browser engine. The
credible interval is numerical quadrature under the stated simplified model,
not calibrated uncertainty for this dataset.

The first launch was blocked by the execution sandbox. Chromium then ran with
the approved local browser-check command. Its first application check caught a
`top` declaration colliding with the browser global. Wrapping the guide in a
local function scope fixed the cause; the subsequent complete runs passed.

## Reproduce

From the repository root, with Python/SciPy and Playwright/Chromium available:

```bash
python3 research/2026-10-29-low-coverage/interactive-guide/build_guide.py
python3 research/2026-10-29-low-coverage/interactive-guide/verification/build_math_reference.py
node research/2026-10-29-low-coverage/interactive-guide/verification/browser_check.cjs
```

Set `PLAYWRIGHT_MODULE` and `CHROMIUM_PATH` if their locations differ from this
workspace. The browser script opens `index.html` through a local file URL and
writes reports, screenshots and an example Markdown download here. It neither
starts an application server nor uses the user's port 8766.

## Repository checks

The repository-prescribed checks were also run successfully:

```bash
python3 -m pytest -q tests/test_merge_rules.py tests/test_hatch_merge.py tests/test_experiments.py tests/test_segment_estimators.py tests/test_integration_review.py
# 166 passed, 6 warnings, 18.48 seconds

cd web
npm run build
# TypeScript checking and Vite production build passed.
```

The warnings concern installed pandas optional dependency versions and existing
FastAPI lifecycle deprecations. No production scientific code was changed.

## Saved visuals

- [Method page, desktop](desktop-method.png)
- [Three-method comparison](desktop-compare.png)
- [Allelic information example](desktop-alleles.png)
- [Saved pilot viewer](desktop-pilot.png)
- [Method page, mobile](mobile-method.png)
- [Count-model example, mobile](mobile-example.png)

Screenshots disable animations and start at the page top so they show the final
layout rather than an intermediate fade or sticky-header position.
