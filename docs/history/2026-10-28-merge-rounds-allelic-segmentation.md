# Merge rounds, allelic segmentation and the depth global test (2026-10-27 – 2026-10-28)

A record of one development session. It covers two branches. `merge-rules` holds the genome merge, the CN re-call and the merge rounds. `allelic-segmentation` was cut from `merge-rules` and holds the allelic tab and the depth global-test change. The format for each scientific change follows [AGENTS.md](../../AGENTS.md): old behaviour, new behaviour, reason, validation.

Everything added here is **exploratory**. It is a set of tools for looking at how segmentation and CN calls behave. Nothing here is a validated pipeline.

## 1. Genome merge, CN re-call and merge rounds (branch `merge-rules`)

- **Spec and plan:**
  - [`docs/superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md`](../superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md)
  - [`docs/superpowers/specs/2026-10-27-merge-rounds-design.md`](../superpowers/specs/2026-10-27-merge-rounds-design.md)
  - [`docs/superpowers/plans/2026-10-27-merge-rounds.md`](../superpowers/plans/2026-10-27-merge-rounds.md)
- **What it does:**
  - "Merge genome & re-call CN" applies the merge rules to every source and chromosome.
  - It then re-calls CN with pyEpi's unmodified `assign_gainloss_new`, first per cluster and then per cell.
  - A re-call can be adopted as a **round**, and later rounds start from the merged segmentation.
  - The karyogram switch has four views: Original · Before · After · Changed.
  - A "Call noise by round" table counts segments, zig-zags and gain-next-to-loss boundaries.
- **Scientific behaviour:**
  - The CN statistic is unchanged.
  - The original breakpoints and calls on disk are never modified.
  - `best_s` is re-fitted per cell on the merged segment means. It is genome-wide, but it changes when the segments change, as the docs explain.
- **Validation:**
  - [`validation/cn-recall/`](../../validation/cn-recall/) and [`validation/hatch-merge/`](../../validation/hatch-merge/).
  - On real data, round 1 at τ 0.5 / ref 8 / any reduced chr8 track breakpoints from 35 to 21, and cluster zig-zags genome-wide from 33 to 7.
  - The same rules applied twice give 0 merges.
- **Deferred:** a spike/dip merge rule.

## 2. Allelic segmentation tab (branch `allelic-segmentation`)

- **Spec and plan:**
  - [`docs/superpowers/specs/2026-10-28-allelic-segmentation-design.md`](../superpowers/specs/2026-10-28-allelic-segmentation-design.md)
  - [`docs/superpowers/plans/2026-10-28-allelic-segmentation.md`](../superpowers/plans/2026-10-28-allelic-segmentation.md)
- **What it does:**
  - **Scan:** runs the depth run's constrained AD scan on a per-cluster allelic pseudobulk track. The data is HA/HB from `haplotype_aware_cnv/hapcnv_outs_from_updated_phasing/aligned_matrices`, on the pyEpi bin grid.
  - **Signals:** signed dBAF, |dBAF|, shrunk |dBAF|, or dBAF oriented to the all-cells direction.
  - **Kept bins:** only bins with AT ≥ min AT go into the scan.
  - **Display:** allelic tracks and segments sit next to the depth tracks and segments.
  - **Comparison tables:** breakpoint agreement, depth segments seen by the allelic signal, and allelic segments seen by depth.
- **Scientific behaviour:**
  - **No new statistic.** The recursion and the geometry check are copied from `beta_revised/src/pyEpiAneufinder/get_breakpoints.py`. The AD statistic and both permutation tests are pyEpi's own, with their fixed seed.
  - **Orientation per cluster is not offered.** At the pseudobulk level it equals |dBAF|. Per-cell orientation is left to the later step that assigns CN with allelic support.
- **Validation:** [`validation/allelic/`](../../validation/allelic/).
  - **Gate:** the copied scan reproduced all 814 original depth breakpoints exactly before anything was built on it. The gate used the chromosome pool, the same way those breakpoints were made.
  - **Result:** at the defaults with the genome pool, AD finds 46–62 allelic breakpoints across all six sources, depending on the signal (4–18 per source). Depth has 117–156 per cluster.

## 3. Depth global permutation test now follows original pyEpi (branch `allelic-segmentation`)

- **Old behaviour:**
  - The depth cluster breakpoints the app loaded (`integration-story/data/raw_accepted_breakpoints.tsv`) were made by `run_bootstrap_diagnostics.py` with `recursive_getbp(seq, seq, …)`.
  - The global permutation test therefore sampled its random segment pairs from **the same chromosome**.
- **New behaviour:**
  - The global test samples from **the cluster pseudobulk's whole genome**, as original pyEpi does. `pyEpiAneufinder.py` passes the cell's whole genome, `counts.X[i, :]`, and `global_permutation_test_ad` draws len_left + len_right random bins from it without replacement. No region is masked.
  - Everything else is unchanged: k 6, 1000 permutations, alpha 0.001, min child 20, min terminal 50, pseudobulk = mean of X per WNN cluster.
- **Reason:** the user asked for the depth module to match how pyEpi originally applies the global test, and chose to replace the app's depth segmentation rather than keep both.
- **Where the data lives:**
  - `tools/segment_depth_genome_pool.py` writes the new breakpoints and audit to `integration-story/data/genome_pool/`. `make_karyograms.py` regenerates that folder's cluster and per-cell CN tables. The config defaults point there.
  - The originals are untouched in `integration-story/data/`. `PYEPI_SCIENTIFIC_CLUSTER_BREAKPOINTS`, `PYEPI_SCIENTIFIC_DERIVED` and `PYEPI_SCIENTIFIC_SEGMENTATION_AUDIT` switch back ([`docs/USAGE.md`](../USAGE.md)).
- **Validation:** [`validation/depth-genome-pool/`](../../validation/depth-genome-pool/).
  - **Breakpoints:** 814 → 823. 767 are shared, 47 appear only in the old set and 56 only in the new one; per source the change is −8 to +6. The run is deterministic.
  - **Tests and gate:** pytest passes, and the gate still reproduces the original 814.
  - **Walkthroughs:** the hatch-merge, cn-recall and allelic walkthroughs pass on the new segmentation. Their genome-pool outputs are kept in `validation/depth-genome-pool/walkthroughs/`, and the archived records in `validation/hatch-merge/` and `validation/cn-recall/` were left as they were.
- **Known consequences:**
  - Saved sessions with merge rounds re-apply their rules to the new breakpoints, and the app now shows a notice when that happens.
  - The alternative flank view and the archived merge-prototype tables still describe the earlier chromosome-pool segmentation, and both the UI and the docs say so.

## Open items

- `tests/browser_hatch_merge_walkthrough.cjs` and `tests/browser_cn_recall_walkthrough.cjs` still write into the archived `validation/hatch-merge/` and `validation/cn-recall/` folders. Re-running them overwrites those records. Point them at a new report folder before the next run.
- Small robustness items remain:
  - A failed allelic job's sibling tasks keep running in the shared worker pool.
  - A corrupt cached allelic result gives an error instead of recomputing.
  - The status check can start two polls while the allelic data is building.
- The allelic tab's loading and error message for merge rounds has no automated test.
- [`docs/notes/2026-10-27-tau-calibration.md`](../notes/2026-10-27-tau-calibration.md) is still deferred: calibrating τ from a known CNV's attenuation.
