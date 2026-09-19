# pyEpi scientific explorer

Separate React/TypeScript + Plotly and FastAPI application, adapted from the supplied Fable explorer. The original `interactive-explorer` is not modified. See [START_HERE.md](START_HERE.md) for the experimental workflow and scientific interpretation, and [validation/VALIDATION_REPORT.md](validation/VALIDATION_REPORT.md) for verification.

## Launch on Linux

```bash
cd scientific-explorer
python3 -m pip install -r requirements.txt
python3 run.py
```

Open **http://127.0.0.1:8766**. The built frontend is included locally. Optional launcher arguments:

```bash
python3 run.py --port 9001 --data-root /path/to/workspace
```

The data root contains `pyepi_results/`, `default_300cells/` and `pyEpiAneufinder/`. By default it is this application's parent directory. No data are copied or edited by an experiment. `cache/` and `sessions/` belong only to this application.

Overrides use the **PYEPI_SCIENTIFIC_** prefix: `ROOT`, `COHORT`, `H5AD`, `PACKAGE`, `CACHE`, `SESSIONS`, `PORT`, `HOST`. The Fable app uses a different prefix and default port.

## Rebuild the interface

```bash
cd web
npm ci
npm run build
```

For frontend development, `npm run dev` serves port 5174 and proxies API requests to 8766. Run the Python service separately. Python >=3.10 and Node >=18 are required by the inherited stack.

## Workspaces

- **Cohort karyogram:** the real 300-cell categorical matrix, coverage/event filters, cell and chromosome navigation.
- **Explore output:** independently computed parent/child AD, X and log2FC tracks; interactive windows; separate AD/X statistics, distributions and autocorrelation.
- **Manual segmentation:** original-X recursion with user-selected boundaries, exact tests, undo/redo and full-genome CN fitting.
- **X experiments:** original-versus-edited data, empirical state simulations, source-region extensions and duplication, width series, event-edge tests, and four CN comparisons.

The inspector is closed initially to give plots more room. The chromosome overview can be compacted; comparison and experiment actions are directly accessible from Explore. The original presentation mode remains available.

## Tests

Start the service first for API/browser tests:

```bash
python3 -m pytest -q tests/test_experiments.py tests/test_segment_estimators.py tests/test_jobs.py
API=http://127.0.0.1:8766 python3 tests/validate.py
python3 tests/parity_reference.py
web/node_modules/.bin/esbuild tests/_bundle/entry.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/stats_bundle.mjs
API=http://127.0.0.1:8766 node tests/parity_browser.mjs
```

Browser tests require Playwright and Chromium. Set `PLAYWRIGHT_MODULE` to an installed Playwright module directory and `CHROMIUM_PATH` if using a separately installed Chromium; otherwise normal Node module/browser discovery is used:

```bash
node tests/browser_experiments.cjs
node tests/browser_workflows.cjs
```

The tests save screenshots, a demonstration session, and result JSON under `validation/`. `browser_experiments.cjs` runs the transformations and width series, checks stale results, and saves/reloads the session. `browser_workflows.cjs` checks drawn-window import, manual non-argmax splitting, tests, undo/redo, CN fitting, and drawing another region on a completed experiment.

## Scientific scope

Experiments use the existing scientific functions unchanged. They modify normalized X only; they do not simulate raw fragments or rerun GC correction. Automatic experimental segmentation is limited to the selected chromosome at production depth k=2; manual trees operate on original X. The original and edited scenarios each supply their own full-genome row to the global permutation test. Every CN fit uses the full genome.

Simulation samples independently from bins with an exact original same-cell production state (loss 0, base 1, gain 2). At least eight donor bins are required. It is a call-conditioned exploratory model, not independently validated biological truth. Extension resamples the selected source; duplication repeats its ordered pattern. Both overwrite neighbouring retained bins and report chromosome-edge clipping.

The most recent 12 experiments per cell/chromosome are retained in browser state and ordinary saved/downloaded sessions. Export before leaving the page. Individual and history JSON exports contain full arrays, transformation parameters, dataset/scientific-code identities, test parameters, and provenance. CSV exports contain the full experimental event span, including extensions. A changed draft never silently relabels an older result as current.

Manual segmentation also supports arithmetic, upper-tail IQR, and two-sided IQR segment means. See [START_HERE.md](START_HERE.md) for definitions, matched-boundary comparisons and scope. Validation evidence is in [validation/VALIDATION_REPORT.md](validation/VALIDATION_REPORT.md).

## Next design phase

See [the multimodal design handoff](../HANDOFF_FOR_MULTIMODAL_DESIGN_AGENT.md) for current architecture, proposed visualization workspaces, scientific constraints, data contracts and setup lessons.


## Integration tabs (multimodal clustering and pseudobulk segmentation)

Two read-only tabs load the exported artifacts of the depth + haplotype integration story
(`../integration-story/`): nothing is re-clustered, re-integrated or re-segmented.

| Input | Default location | Override |
|---|---|---|
| `integrated.h5mu` (cells, labels at four Leiden resolutions, WNN weights, UMAPs, PCA, graphs, allelic abs-dBAF) | `integration-story/data/integrated.h5mu`, else the path recorded in `data/prototype/run_manifest.json` | `PYEPI_SCIENTIFIC_H5MU` |
| pseudobulk breakpoints (`raw_accepted_breakpoints.tsv`) | `integration-story/data/` | `PYEPI_SCIENTIFIC_CLUSTER_BREAKPOINTS` |
| derived cluster segments / per-cell CN (`data/derived/`) | `integration-story/data/derived` | `PYEPI_SCIENTIFIC_DERIVED` |
| merge-prototype and SRD-screen tables (`data/prototype/`) | `integration-story/data/prototype`, else the verified `scientific-explorer/reference/integration-prototype/regenerated` | `PYEPI_SCIENTIFIC_PROTOTYPE` |
| whole story directory | `<root>/integration-story` | `PYEPI_SCIENTIFIC_INTEGRATION` |

The explorer's own `count_matrix.h5ad` X (not the log1p copy inside the h5mu) is used for the
per-cluster pseudobulk profiles and the per-cell chr8 score. At first start a background index
(one pass over X, allelic pseudobulk, WNN diffusion component, per-cell CN as uint8) takes
~6 s and is cached under `cache/integration_<hash>.npz`; `/api/integration/status` reports it.

- **Integration · clusters** — the three UMAPs (depth, allelic, WNN) linked by cell ID with lasso
  selection, cluster chips, colour by cluster / coverage / weights / chr8 score, a cell card with
  memberships, coverage, weights and kNN neighbour overlap (and "open in cell explorer" for
  cohort cells), plus the fig 1a/1b/3/3b/5 diagnostics, including the resolution comparison for
  the chr8 group.
- **Integration · segmentation & CN** — fig 3c (cluster-level CN) and fig 3d (per-cell CN under the
  cluster segmentation, ordered by cluster) as canvases; chromosome tracks per cluster with
  accepted breakpoints (by recursion depth), small internal segments coloured by a selectable
  metric (log2 effect from median / mean / IQR-trimmed levels, screen log2FC, |SRD|, |SRD/√φ|)
  on a selectable flank (min-effect, min-SRD, max-effect, left, right), binned or continuous
  scale, level basis (chromosome median, genome trimmed mean, raw), visual smoothing, and
  addable label layers (SRD-screen proposal incl. `ambiguous_keep`, flank-consistency x/o//,
  retained reason, length + value, cluster CN strip, candidate audit, allelic pseudobulk);
  small-segment diagnostics (figs 4, 6, 7, 8) linked to the tracks; merge-prototype summary.
- Session JSON carries the integration state (labelling, colour, selection, track settings).

API: `/api/integration/{status,meta,pseudobulk.bin,percell_karyogram.bin,cell/{id},candidates}`.

The integration has been reviewed and now includes an on-demand flank comparison and paired whole-cell bootstrap. See [integration review and verification](validation/integration-review/REVIEW.md).
