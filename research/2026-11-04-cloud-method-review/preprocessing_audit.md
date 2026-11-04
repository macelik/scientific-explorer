# Preprocessing and measurement audit

Scope: the retained 3,460 cells × 25,337 bins of the 2026-11-04 data package, the pinned
upstream source `colomemaria/pyEpiAneufinder@e8305ea` (pyproject 0.3.5), and the workbench at
`allelic-segmentation@6d4f48d`. Code locators below refer to those exact revisions.
"Observed" means measured in an experiment saved under [`experiments/`](experiments/). "Plausible"
means a failure mode the code allows but that we did not observe. All results describe MCF10-A
development data and do not establish copy number.

Following the user's direction during the review, GC-correction internals were audited only
far enough to establish provenance and the size of their effects. The weight of the work
went to diagnostics that carry over to other samples (estimands, segmentation calibration,
CN assignment, region evidence and allelic information).

## Pipeline as executed (order matters)

| # | Stage | Formula / rule | Locator | Inputs → outputs |
|---|---|---|---|---|
| 1 | Windowing + blacklist | 100 kb windows; a window overlapping any blacklist interval is **dropped**; GC, AT, N recorded | `makeWindows.py:19-90` | FASTA, `hg38-blacklist.v2.bed` → `binned_genome.csv` (26,510 windows incl. non-autosomal contigs) |
| 2 | Fragment counting | per-cell fragments per window; `minFrags` cell floor | `render_fragments.py:process_fragments` | fragments → cell × bin counts (not in package) |
| 3 | Bin filter | keep bin if `nnz(bin) ≥ (1 − 0.85) · n_cells` | `pyEpiAneufinder.py:212-214` | removes 1,103 autosomal bins here |
| 4 | Cell filter | keep cell if `nnz(cell) > 0.35 · n_bins` | `:222-224` | — |
| 5 | Outlier filter | IQR(k = 1.5) on total counts **and** on per-cell SD of raw counts; intersection kept | `:29-35, :228-238` | — |
| 6 | GC fit | `skmisc.loess(GC, counts, span = 0.75)`, other settings default; on exception print and **return the input counts** | `render_fragments.py:155-164` | per-cell fitted values |
| 7 | GC correction | `correction = mean(counts) / fit`; `corrected = counts · correction` (no guard for fit ≤ 0 or ≈ 0) | `pyEpiAneufinder.py:256-258` | float |
| 8 | Rounding | `np.rint(corrected).astype(int)` (non-finite values would wrap to INT_MIN) | `:260` | int |
| 9 | Save raw | `layers['raw'] = X` (post-QC, pre-GC counts) | `:263` | `layers/raw` |
| 10 | Clip | `expr[expr < 0] = 0` | `:267` | — |
| 11 | Normalize | each row scaled to sum 1e5 | `:269` | `X` (float64 CSR, explicit zeros dropped) |
| 12 | Pseudobulk (workbench) | unweighted mean of `X` rows per WNN cluster (`wnn_leiden_0.3`) and all cells | `tests/validate_allelic_depth_reproduction.py:29-50`, `server/integration.py:282-319`, `server/cn_recall.py:pseudobulk` | 6 profiles |
| 13 | Segmentation (workbench) | AD argmax over geometry-valid cuts; local permutation (1,000) and genome-pool global permutation, both `p < 0.001`; `k = 6`, `min_child_bins = 20`, `min_terminal_bins = 50` | `server/allelic_scan.py:recursive_getbp`, `tools/segment_depth_genome_pool.py` | `genome_pool/raw_accepted_breakpoints.tsv` |
| 14 | CN assignment | divide by upper-trimmed mean (`lb=False`); segment means; weighted scale search over `s ∈ [0.5·median, min(1.2, 0.75·median)]`, weights `(2, 1.5, 1/2 … 1/6)`; Holmes/Watson/combined | `assign_somy.py:30-55, 165-303` (weights at L206) | `karyogram_cluster_segments.tsv` |

The production per-cell run in `pyepi_results/` used upstream `recursive_getbp_df` with `k = 2` per
cell, which is a different segmentation from the workbench cluster scan in stage 13.

## Reconstruction of `X` (experiment 01)

[report](experiments/01_reconstruct_x/report.html) · [summary](experiments/01_reconstruct_x/summary.json) · [analysis](experiments/01_reconstruct_x/analysis.json)

- With Python 3.10.20, numpy 1.24.4 and scikit-misc 0.5.2, stages 6–11 applied to `layers/raw` and
  `var.GC` reproduce `X` **bit-for-bit for 3,449 of 3,460 cells**. No LOESS fit fell back; no
  non-finite value arose before casting.
- The 11 other cells differ in **17 cell-bins**. In all 17 the rebuild's LOESS prediction is small
  and positive (0.0002–0.013 on a raw-count scale), so the rebuilt value is huge (one fragment
  becomes 2,564 counts), while saved `X` has **0**. Every other bin of those cells matches once
  the historical row total is recovered (residual ≤ 6e-14).
- Five environments (scikit-misc 0.1.4, 0.2.0, 0.3.1, 0.5.2, 0.5.3; numpy 1.24 and 2.5; Python
  3.10 and 3.13) give identical predictions (≤ 8e-14) and none goes ≤ 0 at those bins. The
  historical cause is **unresolved**: the run environment, and any guard or post-processing
  outside the pinned file, were not recorded. This is a provenance finding and gives no licence
  to substitute our rebuild for the historical matrix.
- The recorded handoff environment (`ENVIRONMENT.json`: numpy 1.24.4 with scikit-misc 0.5.2)
  cannot be installed as recorded, because scikit-misc 0.5.2 declares numpy ≥ 1.26.4. We
  installed it with `--no-deps`, and the fit imports and runs.

### Observed GC-step artefacts (in the saved `X`, not just the rebuild)

- 2,969 cell-bins in 327 cells have LOESS prediction ≤ 0. All of them sit in **15 bins below the
  1st GC percentile** (GC 0.22–0.32), i.e. at the LOESS boundary. 1,350 raw-positive cell-bins
  are zeroed by negative-value clipping.
- Small positive predictions amplify single fragments: 861 raw-positive cell-bins in 333 cells are
  multiplied by more than 30 (maximum 717×). In 12 bins, amplified cell-bins make up more than 25%
  of the all-cell mean `X`. Examples are chr1:223.5 Mb (mean `X` 76.8 against a genome mean of
  about 3.9) and chr17:81.7 Mb (79.8). These are isolated low-GC bins. None lies in the six saved
  windows, but each can inflate a segment mean that contains it.
- One cell's corrected total is 3.3× its raw total (`norm_factor` 3.33), driven by a handful
  of amplified bins. Its normalization then shrinks every other bin.
- *Plausible, not observed*: an exactly zero prediction would give `inf`/`NaN`, which
  `astype(int)` wraps to INT_MIN and clipping then turns into 0, hiding the invalid value. A
  fit exception would return the input vector, making every nonzero bin equal the cell mean
  and every zero bin `NaN → 0`.

### GC correction and regional biology are confounded (experiments 04, 07)

- In processed `X`, the deep saved dips appear shallower than in raw counts (chr7 8.3–10.7 Mb:
  0.20× its flanks raw, 0.36× in `X`). Across 300 random windows, the processed-vs-raw log2
  attenuation follows the target-minus-flank GC difference with a slope of −15.6 per unit GC.
  Per-cell LOESS raises low-GC regions.
- The published MCF10-A 19q13-qter gain (GC-rich) is 1.82× genome density in raw counts but
  0.93× in `X`; its pooled allele data show excess imbalance (experiment 06). MYC (2.57× raw,
  2.34× `X`) survives. The correction removes whatever depth structure correlates with GC,
  including real dosage on GC-extreme arms.

## Rounding and clipping (experiments 01, 02)

[report](experiments/02_zero_clip_audit/report.html)

- 544,584 raw-positive cell-bins (0.88% of raw-positive cell-bins) are 0 in `X`. In the rebuild
  this splits into 543,217 by `rint` of a value in [0, 0.5) and 1,350 by negative clipping. No
  raw-zero cell-bin becomes positive. Neither matrix stores explicit zeros.
- The loss depends on coverage: a median of 1.41% of raw-positive bins per cell in the
  low-coverage tercile, 0.76% (mid) and 0.48% (high). Processed `X` therefore carries a
  coverage-dependent distortion that a raw-count model with exposure avoids. In the fixed-region
  contrasts, though, raw-count contrasts are **coverage-invariant** (|log(rr_low/rr_high)| ≤ 0.03
  for the saved windows), so rounding does not explain the saved dips.

## QC filters (experiment 02)

- 1,103 autosomal post-blacklist bins are missing from the retained grid. Bins are removed only by
  stage 3, so these failed the ≥15%-nonzero-cells rule. 818 of them are in runs that are not
  assembly-gap-like (N < 0.5). The longest such runs are chr21:18.9–25.3 Mb (64 bins),
  chr13:53.3–56.2, chr16:60.9–63.8, chr8:82.4–84.6 and chr11:54.4–56.6 Mb.
- **Prespecified check:** the bins covering CDKN2A/CDKN2B (chr9:21.9–22.1 Mb), a published
  homozygous deletion in MCF10-A, are among the removed bins. They have N = 0 and normal GC. A
  clonal homozygous deletion has zero counts in every cell, so it fails the bin filter: **the
  event decides its own exclusion**, and downstream steps see only a gap. The epiAneufinder paper
  describes this filter as removing low-mappability regions, so the effect is general, not
  specific to this dataset. The saved chr9 windows contain 16–18 removed bins (24.0–25.7 Mb).
  Without pre-QC counts, nobody can say whether those are deletions or inaccessible deserts.
- The cell filters (stages 4–5) cannot be audited: removed cells are not in the package. Cell
  filters that use the SD of raw counts could preferentially remove cells with large CN changes.
  That is plausible and untested.

## Normalization and the CN reference (experiment 03)

[report](experiments/03_normalization_aggregation/report.html)

- Per-cell normalization to 1e5 means a broad gain lowers every other bin of that cell. The
  CN step adds a second, genome-wide reference: upper-trimmed mean, then the median segment.
- **All 18 source × estimand scale fits end at the lower grid edge, s = 0.5 × median segment
  mean.** The scale is therefore not estimated: the median-weighted segment is defined as CN 2,
  and calls reduce to fixed thresholds on segment mean / genome median. Holmes calls a loss
  below m/s ≈ 1.37, i.e. ≈ 0.68 × median, and Watson at ≤ 1, i.e. ≤ 0.5 × median. Holmes calls a
  gain above 2.45 (≈ 1.22 × median) and Watson at ≥ 3 (≥ 1.5 × median). The weights make losses
  stricter than gains. The code comment says they "penalize deviations from diploid", but the
  weights fall monotonically above k = 2.
- Removing one chromosome from the denominator changes up to about 5% of the calls on the other
  chromosomes (median 0–1.3%). 40 of the 164 `all` segments lie within 0.1 (m/s units) of a
  cut-point, and moving `best_s` by ±3–5% changes 3–18 segments (299–1,725 bins).

## Aggregation estimand (experiment 03)

- The workbench pseudobulk is the **unweighted mean of per-cell normalized, GC-corrected
  profiles**: every cell has equal weight, however many fragments it has. The alternatives are
  the exposure-weighted ratio Σraw / Σexposure and the equal-weight mean of raw rates. These
  answer different questions, and the latter two have no GC correction. On identical boundaries
  they change 20–37% of bin calls, mostly by roughly doubling gain calls in GC-rich regions. This
  is a sensitivity range, not a ranking.

## Measurement units (experiment 06)

[report](experiments/06_allelic_information/report.html)

- AT = HA + HB holds exactly for every cell and bin. The values are non-negative integers;
  coordinates match the AnnData; the barcode sets are equal but the column order differs, so the
  join is by name.
- AT sums DP over selected phased heterozygous SNP sites; one fragment can count at several sites.
  An overdispersion model fitted to base-called bins (ρ ≈ 0.0038, consistent across AT strata)
  still underpredicts the allelic variance of 31% of uncalled random windows. More importantly,
  at the published MYC amplification the measured effective |dBAF| (0.11) is far below what
  nominal counts imply for any clonal one-copy change. Independent-likelihood or pooled allele
  tests that treat AT as independent observations **overstate power**. SNP-level counts, phase
  blocks and unique-molecule counting are required before an allelic likelihood can be calibrated.

## What would be needed to go further

Pre-QC counts and QC masks; `parameter_configuration.txt` from the historical run; fragment files
(for recounting, bin size and deduplication); SNP-level allele counts with phase-set IDs; a
mappability track matched to the alignment reference; and at least one independently
characterized sample (scDNA/WGS) processed identically. See [input_audit.json](input_audit.json).
