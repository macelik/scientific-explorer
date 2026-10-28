# Allelic segmentation: run the depth AD scan on the allelic signal and compare it with depth

**Status:** agreed in conversation 2026-10-28.
**Builds on:** `2026-10-27-merge-rounds-design.md` (branch `merge-rules`, unmerged). Work happens on branch `allelic-segmentation`, cut from `merge-rules`.

## 1. Goal

The project's aim is a haplotype-aware CNV tool. Today segments and CN labels come from read depth only. This step asks how the read-depth segmentation method behaves on the allelic signal:
- Does AD find segments in the allelic signal?
- Are they similar to the depth segments?
- How do the two correspond?

A new tab runs the same recursive Anderson-Darling scan that produced the depth breakpoints on a per-cluster allelic pseudobulk track. It shows the allelic tracks and segments next to the depth tracks and segments, with a small comparison table.

**Agreed in conversation:**
- **Approach 1:** a background job per signal and settings, cached, covering every source and chromosome.
- **Signal selector:** the user chooses which allelic signal AD runs on. This is exploratory, so the point is to see how each signal behaves.
- **No orientation per cluster.** Orientation against a cluster's own pseudobulk equals |dBAF| at the pseudobulk level, so it is not offered. Per-cell orientation (the `haplo_initial` oriented dBAF) belongs to the later step that scores cells over segments and assigns CN labels with allelic support. It is out of scope here.
- **The comparison is a quick first version:** plain counts and per-segment summaries, with no combined score. Integration of the two signals is designed later, after looking at the segmentation.

**Background, measured on real data:**

Cluster pseudobulks are dense enough to segment. Per cell, only 9% of bin×cell entries have any SNP read, with a median of 1.

| Source | Cells | Bins with AT > 0 | Median AT per bin | Bins with AT ≥ 10 |
|---|---|---|---|---|
| cluster0 | 1,286 | 91% | 51 | 75% |
| cluster1 | 833 | 92% | 62 | 77% |
| cluster2 | 603 | 91% | 41 | 73% |
| cluster3 | 556 | 90% | 36 | 72% |
| cluster4 | 182 | 84% | 10 | 51% |
| all | 3,460 | 94% | 202 | — |

- Plain pyEpi `recursive_getbp` (k = 6, 1000 permutations, alpha 0.001) on cluster2 |dBAF| (AT ≥ 10) takes 1.3 s on chr1 and finds 1 breakpoint, and 0.7 s on chr8, also 1 breakpoint.
- Each cluster's per-bin direction is opposite to the all-cells direction in 21–27% of bins with AT ≥ 10 (33% for cluster4).

## 2. Data and signals

### Input

- **Files:** `HA.parquet` and `HB.parquet` in `<EGA>/haplotype_aware_cnv/hapcnv_outs_from_updated_phasing/aligned_matrices/`.
  - Rows are bins, with `chr`, `start` and `end` columns. Columns are cell barcodes.
  - `AT.parquet` holds HA + HB and is read only for the load check.
- **Config:** a new field `allelic_dir`, overridable with `PYEPI_SCIENTIFIC_ALLELIC`.
  - Its default is that folder, resolved relative to the data root: `<data_root>/../../../haplotype_aware_cnv/hapcnv_outs_from_updated_phasing/aligned_matrices`.
  - When the folder is missing, the tab is disabled with the reason. Nothing else in the app changes.
- **The older `hapcnv_outs` sets use an 18,738-bin grid.** They are not supported.

### Load checks

Each failure disables the tab and lists the problem in its status:
1. `chr`/`start`/`end` equal the app's bins exactly, in order: `count_matrix.h5ad` var `seq`/`start`/`end`, 25,337 bins.
2. The barcode set equals the h5mu `obs` index set. The parquet column order differs, so columns are matched by name.
3. AT = HA + HB, checked on a sample of 200 cells.

### Pseudobulk sums

- For each source (`cluster0`…`cluster4` and `all`), ΣHA and ΣHB per bin over the source's member cells. Membership is the one the integration store already uses for depth sources.
- **Streaming:** the matrices are read in column batches, never whole. About 350 MB each as float32.
- **Kept:** only the sums, 6 × 25,337 × 2 float64.
- **Disk cache:** `<cache_dir>/allelic/pseudobulk_<key>.npz`. The key covers the HA and HB file sizes and modification times, the membership hash and the code version.

### Signals

HA, HB and AT are a source's pseudobulk sums for one bin.

| id | Label | Value |
|---|---|---|
| `dbaf` | signed dBAF | (HB − HA) / (2·AT) |
| `abs_dbaf` | \|dBAF\| | \|HB − HA\| / (2·AT) |
| `shrunk_abs_dbaf` | shrunk \|dBAF\| | \|HB − HA\| / (2·(AT + 2)), the h5mu formula with α = 1 |
| `oriented_all` | oriented to all cells | sign(ΣHB_all − ΣHA_all) · (HB − HA) / (2·AT) |

**Kept bins:** AT ≥ `min_at` (default 10). For `oriented_all`, the all-cells HA and HB must also differ. Every other bin is skipped. Skipped bins never enter the scan and are never imputed.

## 3. The segmentation job

### Scan

- **A port of `recursive_scan` from `SNP_project/beta_revised/src/segmentation/wnn_bp_exp.py`,** behaving the same:
  - At each node: compute `seq_dist_ad`.
  - If the max-AD candidate breaks a size rule, record it as `rejected_constraint` and test the best valid candidate instead. If no valid candidate exists, stop the branch.
  - Accept a candidate when p_local < alpha and p_global < alpha. Recurse into both children until depth k.
- **Size rules:**
  - Both children must have at least `min_child_bins` kept bins.
  - The breakpoint must be at least `min_terminal_bins` kept bins from the chromosome's first and last kept bins.
- **The local permutation test** shuffles within the node. The global test draws from all of the source's kept bins genome-wide, matching the depth scan's use of the whole-genome pseudobulk.
- **Statistics:** pyEpi's `seq_dist_ad`, `dist_ad`, `permutation_test_ad` and `global_permutation_test_ad`, called through `server/pyepi_adapter.py` with the default fixed `random_state`. Nothing statistical is re-implemented.

### Settings and defaults

| Setting | Default | Range |
|---|---|---|
| signal | `abs_dbaf` | the four ids above |
| min_at | 10 | 1–1000 |
| k | 6 | 1–8 |
| n_permutations | 1000 | 100–5000 |
| alpha | 0.001 | 0 < alpha ≤ 0.1 |
| min_child_bins | 20 | 1–500 |
| min_terminal_bins | 50 | 0–500 |

The defaults are the depth run's settings, inferred from `segmentation_candidate_audit.tsv`. Section 7 checks them by reproducing the depth breakpoints.

### Mapping back to the genome

- A breakpoint between kept bins i and i+1 of a chromosome sits at the absolute bin of kept bin i+1. That is the first bin of the right child.
- If skipped bins lie between the two kept bins, the breakpoint also carries `gap_start` (the first skipped bin) and `gap_end` (kept bin i+1). It is drawn as a band.

### Job

- **Request:** `POST /api/allelic/segment` with body `{signal, settings}` returns `{id, status}`. Invalid values return 422 with the field name.
- **Job id:** a hash of the signal, the settings, the pseudobulk key and the sha256 of the scan module.
- **Repeats:** a finished job with the same id returns at once, from memory or from `<cache_dir>/allelic/seg_<id>.json.gz`.
- **Scheduling:** at most one allelic job runs. A new request while one runs is queued; a later different request replaces the queued one. The running job is never killed.
- **Execution:** 6 sources × 22 chromosomes, each (source, chromosome) a task in the existing process pool.
- **Progress:** `{done, total, source, chromosome}`.
- **Failure:** a failed task fails the job with its message. Nothing is cached for a failed job.

**`GET /api/allelic/segment/{id}`** returns:
- `status` (`queued` | `running` | `done` | `error`), `progress` and `error`;
- `settings` and `signal`;
- when done, `result`:
  - `breakpoints`, one row per accepted breakpoint: source, chromosome, absolute_bin, bp_genomic_start, bp_genomic_end, gap_start, gap_end, recursion_depth, ad_score, p_local, p_global;
  - `audit`, every candidate in the columns of `segmentation_candidate_audit.tsv`, with bin fields in absolute bins plus kept-bin indices;
  - `segments`, per source and chromosome: bin_start, bin_end, n_bins, kept_bins, mean, median, at_sum, ha_sum, hb_sum.

**Other routes:**
- `GET /api/allelic/status` returns `{available, problems, notes, allelic_dir, sources}`.
- `GET /api/allelic/tracks?chrom=&sources=&signal=&min_at=` returns per source the per-bin signal (null where skipped), AT, and the bin positions. It reads only the cached sums and needs no job.

## 4. UI: tab "Integration · allelic segmentation"

- **Placement:** after "Integration · segmentation & CN" in `App.tsx`.
- **Component:** `AllelicSegmentationView.tsx`. State lives in the integration store and is saved with the session.

### Controls

- **Signal:** a select with the four labels.
- **min AT:** a number field.
- **AD settings:** a collapsed group with k, permutations, alpha, min child and min terminal.
- **Run segmentation:**
  - Shows the progress line: "running · cluster2 · chr8 · 74/132".
  - If the shown result's signal and settings differ from the current controls, a note reads "Segments shown are from an earlier run (signal …, settings …). Run again to update."
- **Depth segments: current round / original.**
  - **Current round** uses the merge-rounds base: the `/api/integration/base` breakpoints and cluster segments for the session's rounds. It is the default.
  - **Original** uses the on-disk breakpoints and cluster segments.
- **Chromosome select and source checkboxes,** like the segmentation tab.

### Tracks

For each checked source, two panels share one x axis (bin position on the chromosome).

1. **Depth panel.** The depth pseudobulk (`pb_x`) with the depth segments coloured by CN call, drawn as on the segmentation tab.
2. **Allelic panel:**
   - one dot per kept bin, with opacity scaled by AT (clamped to 0.15–1 between min_at and the source's 90th percentile AT);
   - skipped bins as small grey ticks on the x axis;
   - allelic segment means as horizontal lines;
   - accepted allelic breakpoints as solid vertical lines, and gap breakpoints as shaded bands;
   - depth breakpoints as faint dashed vertical lines;
   - y range −0.5 to 0.5 with a zero line for `dbaf` and `oriented_all`, 0 to 0.5 otherwise.

Before any run, the allelic panel shows the dots only, with the hint "Run segmentation to see allelic segments."

## 5. Comparison tables

- **Placement:** under the tracks, in one region titled "Depth vs allelic segments".
- **Scope switch:** this chromosome / genome-wide. It applies to table 1; tables 2 and 3 always cover the chosen chromosome.
- **Inputs:** the matching distance **N** (default 3 bins) is a number field. The tables need a finished run; before one, the region says so.

**Table 1: Breakpoint agreement.** One row per source.

| Column | Definition |
|---|---|
| depth bps | depth breakpoints of the source (current round or original) |
| allelic bps | accepted allelic breakpoints |
| matched | depth–allelic pairs matched within N bins |
| depth only | depth breakpoints without a match |
| allelic only | allelic breakpoints without a match |
| median offset | median \|depth bin − allelic bin\| over matched pairs; for a gap breakpoint, the distance to the gap |

**Table 2: Depth segments seen by the allelic signal.** One row per depth segment of the checked sources:
- source, bin_start, bin_end, CN call;
- kept bins;
- mean and sd of the allelic signal over kept bins;
- allelic bps inside: the count of allelic breakpoints more than N bins from both segment edges.

**Table 3: Allelic segments seen by depth.** One row per allelic segment of the checked sources:
- source, bin_start, bin_end, kept bins, mean signal;
- depth calls spanned: each depth CN call with its share of the segment's bins, largest first, e.g. "1.0: 80% · 1.5: 20%".

**Definitions**, shown as text under the tables:
- **Matching is one-to-one, nearest first.** All depth–allelic pairs within N bins are sorted by distance and accepted greedily, skipping any breakpoint already matched. Ties go to the lower depth bin.
- **A gap breakpoint is at distance 0 from any bin inside its gap.** Outside it, the distance is to the nearer gap edge.
- **Distances count genomic bins,** not kept bins.

**Export:** `allelic_vs_depth.tsv`. It writes all three tables in one file, separated by commented headers, with the run's signal, settings, N and depth-segment choice in the header lines.

## 6. Errors

- **Allelic data missing or a load check failed:** the tab shows the problems from `/api/allelic/status` and disables Run. The other tabs are unaffected.
- **A source with no kept bins on a chromosome:** no breakpoints and one segment with `kept_bins` 0 and a null mean. The tracks show only grey ticks.
- **Job error:** the message is shown next to Run, and the previous result stays on screen with the "earlier run" note.
- **Session with invalid saved settings:** fall back to the defaults with the notice "Saved allelic settings could not be applied: …".

## 7. Testing

- **Reproduction of the depth breakpoints, a validation script:**
  - Run the ported scan on the depth pseudobulk (mean of `X` over each source's cells, as `wnn_bp_exp.py` does), with every bin kept, k = 6, 1000 permutations, alpha 0.001, min child 20 and min terminal 50.
  - It must match `integration-story/data/raw_accepted_breakpoints.tsv` on the same (source, chromosome, absolute_bin) set, and on recursion_depth, p_local and p_global for each.
  - **If it does not match,** the defaults are wrong. Stop and report the differences before building the UI.
- **Unit tests, synthetic, no real data:**
  - the signal formulas, including AT = 0, HA = HB and the oriented sign;
  - kept-bin mapping and gap breakpoints;
  - the size rules, and the fallback to the best valid candidate;
  - the matching rules (one-to-one, ties, gaps);
  - the depth-calls-spanned shares;
  - the job id changes with each setting;
  - repeated requests reuse the job;
  - the cache survives a manager restart;
  - the load checks reject shuffled coordinates, missing barcodes and AT ≠ HA + HB.
- **Route tests:** status, tracks, segment start and poll, and 422 on bad settings.
- **Web helper tests:** settings validation, matching, the shares string, the TSV writer, and the session round trip.
- **Browser walkthrough on the private 8799 server:** open the tab, see the tracks before a run, run |dBAF|, see the breakpoints and tables, switch to signed dBAF and run again, then switch the depth segments to original.
- **Real-data record:** for each signal at the defaults, the allelic breakpoints per source and the table 1 counts genome-wide.

## 8. Out of scope

- per-cell allelic scores, per-cell orientation, and CN labels with allelic support;
- any combined depth + allelic segmentation or score;
- merge rules on allelic segments;
- the older 18,738-bin allelic outputs and the MAF/AI matrices.
