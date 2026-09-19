# Integration story review and generated alternative

Reviewed 2026-09-19. Scope: the supplied `integration-story` figures, explanations and two generator scripts, plus Fable's integration in `scientific-explorer`. Original story artifacts and `interactive-explorer` were preserved. This is a targeted scientific/data-flow review, not independent validation of the modified caller or proof of a biological subclone.

## What is coherent

The overall chain is clear: depth and allelic representations → modality and integrated graphs/embeddings → fixed WNN cluster membership → mean normalized-depth pseudobulk → exported cluster-specific segmentation → cluster and per-cell CN on those boundaries. The existing integration keeps cluster versus cell-level fits distinct and uses the original normalized X rather than the log1p matrix in the depth modality.

The supplied main-workspace and story count matrices have matching X data arrays. Runtime checks validate cell ordering and the depth-bin coordinates used in the integration. The distinction between per-bin profile roughness and whole-cell effect stability is scientifically useful and is retained in the new comparison.

## Confirmed implementation defects addressed

| Finding | Change |
|---|---|
| The user intentionally removed `data/prototype`. Empty flank histograms formatted an undefined median and could crash the segmentation tab. | Show an explicit availability notice; do not render empty source-score or merge-result panels. Offer the independent generator. Missing exports are not zero effects or zero merges. |
| Allelic nonzero fraction was called “covered in …% of cells”. Nonzero imbalance excludes balanced covered bins and cannot identify all covered cells. | Hover labels now say “nonzero imbalance in …% of cells”. |
| Allelic overlay used invalid Plotly axis identifiers such as `y2h`. | Give each overlay a valid numbered axis, mapped to its depth panel. |
| Raw-depth breakpoint/shading shapes used a data-axis height of 1e9. | Use subplot-domain coordinates for vertical extents, preserving the actual data range. |
| Gap logic connected tracks over smaller missing-bin intervals. | Use actual bin ends when supplied and break at omitted genomic intervals. |
| Maximum-effect flank selection mishandled a missing right-side effect. | Select the available finite side. |
| Some cluster/weight/histogram medians selected the upper middle observation for even sample sizes. | Average the two central values. |
| Integration cache identity omitted the depth input and several derived/prototype files. | Fingerprint depth, MuData, boundaries, derived files, prototype files and adapter source; use nanosecond modification times. These are file identity checks, not full content hashes for every large input. |
| Loaded run manifest was overwritten by initialization. | Initialize it before resolving inputs and preserve it. The original manifest was subsequently located upstream and its numerical fields verified by regeneration. |
| Cell-order checking occurred only during cache building. | Validate order before using a cache, reject duplicate MuData IDs, and validate depth bin ends as well as starts/chromosomes. |

## Generated alternative

Open **Integration · segmentation & CN → Alternative flank view → Generate alternative flank view**.

Choose chromosome, source (or selected track rows), and maximum retained-bin width. The generator uses mean normalized X over the original WNN res-0.3 members and the fixed exported accepted breakpoints. It selects internal segments with both neighbors and compares each complete adjacent segment as a flank. This is a separate descriptive comparison; the original merge prototype is now also available below.

- Switch between arithmetic mean, median and two-sided 1.5-IQR mean. Quantiles use NumPy's linear percentile convention and fences are inclusive.
- See left and right log2 effects together, without choosing one flank and discarding the other. Positive/positive means elevation relative to both neighbors; negative/negative means depression relative to both; opposite signs mean the segment is intermediate. These patterns are not themselves CN calls.
- Select a segment to inspect genomic context, actual bin-value distributions, mean/median/IQR levels, original and retained lengths, zero fraction, sample standard deviation and CV.
- Ratios involving a zero level remain undefined; no pseudocount is inserted. Undefined points remain available in the selector and exports.
- The generated profiles come from the integration's float32 pseudobulk cache. The independent numerical check tolerates only its rounding difference from direct float64 matrix means.
- Changing scope marks the frozen result stale. CSV exports include every estimator and both sides; JSON includes profile values, provenance, coordinate conventions and implementation hash. Generated views are local to the tab mount, not added to normal session persistence; download them before leaving if needed.

### Paired whole-cell bootstrap

For a selected segment, **Compute paired cell bootstrap** reads the original member-cell X values and computes each cell's arithmetic-mean depth in the left/segment/right intervals. Each replicate resamples the same cells jointly for all three regions (200 replicates, seed 42). Averaging these cell-level arithmetic means is equivalent to taking the arithmetic mean of the reconstructed mean pseudobulk over that region.

The output gives the observed arithmetic-mean contrast, exploratory 95% percentile interval, bootstrap sample SD (SE), valid/attempted replicates, sign agreement and replicate distribution. Zero-level replicates are excluded and counted. It does **not** bootstrap median or IQR segment estimates; the UI explicitly distinguishes those displayed estimates from this arithmetic-only bootstrap.

Membership and selected boundaries stay fixed. The bootstrap does not account for selecting/clustering/segmenting the same data, does not estimate per-cell call accuracy, and does not establish independent biological replication. Its interval and sign fraction are not calibrated significance or merge decisions. No SRD dispersion estimate is fabricated.

Example saved in `generated-bootstrap.json`: cluster4 chr8 [396,416), 182 cells. The left arithmetic effect is about −1.206 log2 with a percentile interval about [−1.265,−1.135]. It can be stable across member cells while the selected region retains substantial bin-to-bin variation. This is an illustrative result conditional on this dataset and selection.

## Story interpretation that needs care

1. `FIGURES_EXPLAINED.md` writes variance as `mean((x−mean)^2)` while also saying `ddof=1`. Those definitions differ: sample variance divides the squared-deviation sum by n−1. The alternative uses and labels sample SD explicitly.
2. A trimmed mean is not guaranteed to lie between the arithmetic mean and median. It can reduce outlier influence, but the sentence claiming it “sits between them” should not be treated as a mathematical property.
3. A near-zero global Spearman correlation between the chr8 score and coverage does not prove coverage cannot explain or contribute to the pattern. Nonlinear and within-cluster relationships can differ. The visible chr8 pattern is a diagnostic, not validated clonal identity.
4. Negative weight-versus-coverage correlation does not establish that **both** weight tails are low coverage. The README itself reports high median coverage in the bottom depth-weight tail. Inspect the two tails separately before summarizing them together.
5. HA/HB should be called haplotype A/B unless parent-of-origin phasing is established separately. The supplied story calls them maternal/paternal, but that provenance is not demonstrated by these artifacts.
6. The haplotype matrix has only positional feature IDs, whereas the depth matrix has chromosome/start/end. Its equal column count is consistent with alignment but does not independently prove it. The integration now reports that positional alignment is assumed. Recover the construction script or coordinate mapping to verify it.
7. The SRD plots are descriptive effect/noise comparisons. The reference 3.3 is not a calibrated significance cutoff; color correlations do not establish a causal explanation of the discrepancy. Whole-cell stability and per-bin scatter measure different things.
8. Cluster-level CN and member-cell CN share segmentation. Agreement therefore evaluates calling under those boundaries, not independent breakpoint recovery in each cell.
9. The integration's pseudobulk allelic track is the mean of per-cell absolute, coverage-shrunk imbalance. It is not imbalance computed after pooling HA/HB counts. A zero may be balanced or uncovered.
10. Existing original-track smoothing acts on displayed values (including log ratios and their display clipping), so it is not a new smoothed signal for calling. The alternative raw-value track is unsmoothed.

These interpretation notes preserve the original story text. The original diagnostic scripts were subsequently rerun separately as described below.

## Correction: original generators and data found and reproduced

The earlier claim that exact reproduction required the user to supply code/configuration was too broad. Those files were absent from the copied `integration-story` bundle, but the README pointed to the upstream workspace `/home/mami/Desktop/EGA/SNP_project/beta_revised`, where the generators, manifest, candidate audit and reference tables were available. Following those paths resolved the missing inputs; no additional files are needed from the user for these diagnostics.

The original depth matrix, MuData and accepted-breakpoint files are byte-for-byte identical to the story copies (full SHA-256 checked). The three merge-prototype modules match the hashes recorded in the original manifest. Unmodified source snapshots and archived reference tables now live under [`reference/integration-prototype`](../../reference/integration-prototype/README.md), separately from the intentionally deleted story prototype directory.

Both generators were rerun against the local story inputs:

- **Merge prototype:** seven chromosomes × five clusters, 200 whole-cell replicates, seed 20260907, median segment levels, τ=0.2, stability=0.9, original chain guard and validity requirements. Nine tables match the archive byte-for-byte, including 162 flank-score rows and 324 SRD rows. All five replicate-draw hashes and numerical manifest fields match. Results: 0 merges, 299 retained boundaries and 334 segments; 0 merges across the nine-setting sensitivity grid. All 31 prototype PNG figures also match their archived PNGs byte-for-byte.
- **SRD screen:** the separate original generator, 200 replicates, seed 42. Its 196-row score table and 196-row bootstrap table also match byte-for-byte; ten figures were regenerated. The all-cell bootstrap resamples within each cluster and retains cluster-size weights.
- **Candidate audit:** copied unchanged from the original accepted-breakpoint run. The AD segmentation itself was not rerun.

The explorer automatically loads the verified local reference when `integration-story/data/prototype` is absent; `PYEPI_SCIENTIFIC_PROTOTYPE` remains an explicit override. The original diagnostic panels, merge summary and candidate overlays can now use these tables. The new alternative arithmetic bootstrap remains separate: sharing seed 42 with the SRD screen does not imply identical replicate streams or estimators.

See `reference/integration-prototype/regenerated/verification.json` for table/input/code hashes and `reproduce.py` for repeatable commands. The original `make_figures.py` still contains historical paths and also recomputes embeddings; it is not a portable one-command rebuild of the full story. This verification covers SRD, dispersion, screening/bootstrap and merge diagnostics, not a fresh WNN clustering or CN-calling run.

## Verification

- 39 Python checks passed: the previous 30 experiment/estimator checks plus 9 targeted checks for summaries, zero handling, boundaries, fingerprints and paired-cell resampling.
- Frontend helper checks passed for missing-flank selection, even-sized median and breaking a line over one missing genomic bin.
- Live independent matrix validation passed: **10 pseudobulk profiles, 35 internal segments across chr8/chr9, and the 182-cell/200-replicate paired bootstrap**. See `numerical-results.json` and `real-bootstrap-reference.json`.
- TypeScript and production build passed. Browser verification evidence is recorded in `browser-results.json` with screenshots and generated JSON/CSV alongside it.
- All 66 recorded original Fable source hashes still match. Original story files were read-only inputs to this work. Browser checks passed with no JavaScript errors, six generated cluster4 chr8 segments, JSON/CSV exports, bootstrap, stale scope marking, and five valid allelic overlay axes.

Key files: `server/flank_view.py`, integration routes/adapter, `AlternativeFlankView.tsx`, `PairedFlankBootstrap.tsx`, and focused corrections in the existing integration components. Tests: `test_integration_review.py`, `integration_helpers.ts`, `validate_integration_review.py`, `browser_integration_review.cjs`.
