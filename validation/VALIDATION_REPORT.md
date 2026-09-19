# Scientific explorer validation

Validated 2026-09-18. Separate application at http://127.0.0.1:8766; no public deployment. Source, data, caches and sessions belonging to `interactive-explorer` were not modified. All 66 recorded original source hashes still match `original_source_hashes.json`.

## Numerical checks

- 30 Python tests passed across `test_experiments.py` and `test_segment_estimators.py`. Covers immutable regional editing, sampling and duplication, chromosome-edge clipping, total preservation, direct production recursion/permutation/CN parity, undefined CN, and all three manual segment estimators on a real cell.
- One additional worker regression passed: a computation worker does not inherit a parent socket. Uses fresh spawned processes after an old forked worker was found retaining the service port during restart.
- **32/32 inherited scientific checks passed** against the updated server (`scientific-parity.log`). Includes production segment reconstruction and full-genome CN and scale reproduction for all 300 cells, original permutation functions, coordinate mapping, and manual boundary splicing.
- **9,294/9,294 browser/Python statistic comparisons passed** (`browser-python-parity.log`) across four node examples. Statistic implementations were unchanged by the subsequent IQR addition.
- TypeScript checking and the Vite production build passed.
- Narrow-layout inspection caught Plotly containers collapsing to zero height during resize. The shared plot wrapper now reserves the declared chart height; browser geometry assertions check that experiment panels retain their full chart height.

IQR synthetic tests separately verify the original normalization and scale-search functions with full segment lengths as weights. Real-cell API tests verify returned estimates, retained counts, unchanged boundaries, original-X preservation, arithmetic references and changed-call counts. Arithmetic remains the default and calls the original scientific implementation.

## Browser checks

`workflow-results.json`: drawn explorer-window import; arbitrary manual boundary at chromosome-local bin 169; split-test result; undo/redo; full-genome CN; both IQR options; stale-result marking; CSV download; session restoration of the selected estimator and frozen output; drawing a new region on an experiment result. No JavaScript errors.

`browser-results.json`: multiply, empirical loss simulation, leftward duplication, rightward extension, and five additional width-series runs; nine retained results, draft changes, tab navigation, export, reload and upload restore. No JavaScript errors. Latest uncached job walkthrough: initial load 1.35 s; multiply 3.59 s; simulation 3.09 s; duplication 3.30 s; extension 3.76 s. These are single local measurements, not performance guarantees.

Screenshots include `manual-iqr.png`, `experiment-simulate.png`, `width-series.png`, `experiment-narrow.png`, `explorer.png` and `cohort.png`. Demonstration sessions: `manual-iqr-session.json` and `browser-session.json`.

In the saved manual example (TTAGGCTAGGCCGGAA-1, chr9, boundary 169), arithmetic versus two-sided IQR on identical boundaries changes 169/1,031 chromosome bins and 4,544 elsewhere. The shared scale changes from 0.575882 to 0.486180. This is one exploratory result, not evidence that trimming improves accuracy.

## Interpretation and limits

- Interventions operate on normalized, GC-corrected X. Raw preprocessing is not rerun.
- State sampling is seeded, independent empirical resampling conditioned on existing calls. It does not preserve spatial correlation or establish biological ground truth. Duplication deliberately repeats the ordered pattern.
- Extension overwrites existing retained bins without inserting genomic coordinates. Every run begins from original X.
- Automatic experiments use the supplied k=2 recursion, 1,000 permutations and strict local/global p < 0.001. Only the selected chromosome is resegmented; CN fits and global null pools are full-genome.
- Event-edge recovery means the exact bin boundary. Width series use one realization per width; they are not calibrated detection probabilities. The copied production path has no explicit 50/100-bin event-size cutoff.
- Manual segmentation operates on original X, separately from the X experiment tab. The IQR option applies to all full-genome segments, retains the original genome baseline and original bin-count weights, and can alter CN elsewhere through both summaries and scale. Two-sided fences can exclude zeros. Invalid positive-scale fits are reported explicitly.
- History retains at most 12 runs per cell/chromosome. Export retained history or sessions to keep experiments.
- Cursor synchronization across separate parent/child plot components was not added. Existing linked tracks within each node and comparison-window workflows remain available.
- Tests emitted pre-existing optional pandas acceleration-version warnings and FastAPI lifecycle deprecation warnings. These did not affect the passing checks.

## Addendum 2026-09-18 — Integration tabs (multimodal clusters, pseudobulk segmentation & CN)

Inputs resolved: `integrated.h5mu` via the path recorded in `data/prototype/run_manifest.json`
(the symlinks in `integration-story/data/` are broken from this checkout and are reported,
not silently skipped); breakpoints, derived and prototype tables from `integration-story/data/`.
Small prototype tables were copied into `integration-story/data/prototype/` (read-only inputs).

Loader checks (server, at start-up): h5mu depth bins equal the explorer's `count_matrix` bins
(seq/start), cell order identical across mudata obs, both modalities and `count_matrix.h5ad`
(3,460 cells); allelic modality has 25,337 bins and is built from the same var triplets in
`wnn_agg_exp.py` (asserted there); cluster sizes from the h5mu equal the manifest and the
breakpoint table (1286 / 833 / 603 / 556 / 182). `mod/depth/X` is `log1p(X)` (verified), so the
pseudobulk uses the explorer's X.

Reproduced story numbers (computed in the app from the artifacts):

| quantity | app | story |
|---|---|---|
| chr8 score median, cluster 4 / others | 0.947 / −0.03 … −0.08 | ≈ 0.95 / −0.03 … −0.08 |
| resolution-1.0 group containing cluster 4 | group 5 keeps 182, adds 55: 48 from cluster 3 (47 > 0.5), 6 from cluster 0, 1 from cluster 2 | same |
| chr8 tracks (fig 9/10) | cluster 4 small segments at 43–47 (≥1.00), 62–65 (≥1.00), 66–77 Mb (0.50–0.75, o); clusters 0/1/3 at 100–104 Mb (0.20–0.50, x) and 134–138 Mb (≥1.00, o) | same bands and hatches |
| small segments with metrics | 162 (7 pilot chromosomes) | 162 |
| SRD-screen proposals | 15 ambiguous_keep, 109 keep_focal, 39 merge_left, 33 merge_right | same table |

Browser checks (Chromium 1500×1000): both tabs render without console errors after fixes;
3 WebGL contexts only (embeddings), all other plots SVG; cluster chip selection highlights 182
cells in all three embeddings and preselects the cluster in the segmentation tab; clicking a
cohort cell shows its card with the "open in cell explorer" action; a lasso-style selection
(50 cells) is reflected in the controls; clicking cluster-row 4 of the fig 3c karyogram at a chr9
position switches the tracks to chr9 and selects cluster 4; clicking a shaded small segment
opens its flank table; all ten label layers toggle on together (proposal, x/o//, reason, values,
CN strip, candidate audit with accepted / rejected / geometry-blocked markers, allelic overlay);
clicking a diagnostic point jumps the tracks to that chromosome and segment; a saved session
restores tab, chromosome, metric, selected segment and layers; presentation mode disables hover
on every integration plot.

Measured: integration meta 0.9 MB in 0.29 s; pseudobulk block 1.1 MB in 0.16 s; per-cell CN
karyogram 1.3 MB gzip (87 MB decoded in the browser); first-time index 6.0 s, cached reload
0.4 s; clusters tab first render 1.4 s, segmentation tab 0.5 s after the first.

Bugs found and fixed while inspecting: hook-order crash in the clusters view (React #310);
WebGL context exhaustion blanking the embeddings; null AD scores in the candidate audit
crashing the tracks; horizontal segment labels overflowing (now compact and vertical inside
each span).

Limitations: small-segment metrics exist only for the seven pilot chromosomes and for
segments ≤ 100 bins (as in the prototype); the allelic overlay is the mean abs-dBAF where zero
can mean no coverage or balance (coverage fraction is in the hover); the per-cell karyogram
samples rows at "fit" height (switch to 1 px rows to see every cell); no recomputation of
clustering, integration or the modified caller is offered.
