# Cloud method review — 2026-11-04

Independent review of the scATAC CNV workflow on `allelic-segmentation` (checkout
`6d4f48d`), using the verified data release `cloud-review-data-2026-11-04` and pinned
pyEpiAneufinder `e8305ea`. MCF10-A is development data; nothing here is CN truth.


> **Start here:** [final_report.html](final_report.html) is the consolidated report for experiments 00–19: findings, what to
> change, what worked, explicit methods with assumptions, and the path toward a haplotype-aware caller.

## Recommendation (details: [recommendation.md](recommendation.md))

Put the next change in **evidence integration and the CN reference**, not in GC/rounding
tuning and not in a new caller:

1. Report a **per-region evidence layer** next to the unchanged pyEpi calls. It covers raw-count
   contrast against a same-width random-window baseline, coverage invariance, compatibility with
   a one-copy change, a per-cell mixture bound, allele observability with an empirically
   calibrated label, and QC-gap adjacency. Its states include *insufficient evidence*.
2. **Keep QC-removed bins as an explicit "unobservable" state**: the bin filter removed the
   published CDKN2A/B homozygous deletion.
3. Get a **cross-sample technical baseline and mappability** before changing formulas. Then
   evaluate a raw-count, exposure-aware model and a joint depth–allele model on independently
   characterized samples, with Numbat-multiome as the published comparator once its inputs exist.

Ranked alternatives: [method_comparison.tsv](method_comparison.tsv). Primary sources and access
limits: [sources.json](sources.json).

## What is established (local, reproducible)

| Finding | Evidence |
|---|---|
| Saved X is the pinned formula on `layers/raw`: 3,449/3,460 cells bit-identical; 17 cell-bins at near-zero LOESS predictions are 0 historically and unexplained in five environments | [exp 01](experiments/01_reconstruct_x/report.html) |
| The GC step amplifies single fragments up to 717× at low-GC boundary bins and zeroes 1,350 raw-positive cell-bins by clipping; rounding removes 0.5–1.4% of raw-positive bins per cell, more in low-coverage cells | [exp 01](experiments/01_reconstruct_x/report.html), [exp 02](experiments/02_zero_clip_audit/report.html) |
| The bin filter removed 1,103 autosomal bins, including the CDKN2A/B bins (published homozygous deletion): an event can remove itself | [exp 02](experiments/02_zero_clip_audit/report.html) |
| Boundaries are reproducible (split-half r > 0.99; window edges in both halves); exchangeable-null acceptance ≈ 1/300 at α = 0.001; significance mostly reflects reproducible non-sampling structure | [exp 05](experiments/05_segmentation_calibration/report.html) |
| Uncalled random 25-bin windows vary 0.44–1.72× their flanks; chr7 saved dips are 0.13–0.34×, below the one-copy floor, coverage-invariant, shared across clusters, with ≤ 2.3% of cells compatible with near-complete loss | [exp 04](experiments/04_fixed_region_contrasts/report.html) |
| chr7 windows keep balanced heterozygous-site observations; 9p21 shows the expected het-site depletion; MYC/5q are allelic positive controls; nominal AT overstates allelic power | [exp 06](experiments/06_allelic_information/report.html) |
| CN scale fits all sit at the grid edge, so calls are thresholds on mean/genome median; the estimand changes 20–37% of calls; one chromosome changes up to 5% elsewhere | [exp 03](experiments/03_normalization_aggregation/report.html) |
| GC correction attenuates low-GC dips and erases the raw 19q13-qter excess (1.82× → 0.93×), whose alleles nonetheless show imbalance | [exp 04](experiments/04_fixed_region_contrasts/report.html), [exp 07](experiments/07_subwindow_profiles/report.html) |

## What remains unknown

The CN state of every saved window; whether the chr7 plateaus are mappability, accessibility or
dosage; the cause of the 17 historical zeros; the contents of QC-removed bins; calibrated allelic
power at SNP level; any tool's performance beyond MCF10-A. Blocked and skipped items, with
reasons: [status.json](status.json).

## Files

- [input_audit.json](input_audit.json): hashes, dimensions, units, alignment, source identity,
  stage availability, missing prerequisites ([report](experiments/00_input_audit/report.html)).
- [preprocessing_audit.md](preprocessing_audit.md): stage-by-stage formulas with code locators,
  failure hypotheses, diagnostics and reconstruction.
- [experiments/](experiments/README.md): scripts, commands, seeds, environment, tables and one
  HTML report per diagnostic (00–19; 08–19 implement the 2026-10-29 method designs: E:1/E:2/E:3 count and
  allele models, B:2/B:3 count changepoint and joint HMM, C:3 baselines, per-modality model choice, major/minor
  haplotype orientation, allele-specific states, cell-level joint mixture, spike-in benchmark, masking diagnosis).
- [final_report.html](final_report.html), built by [final_report/build.py](final_report/build.py) from the experiments' outputs.
- [method_comparison.tsv](method_comparison.tsv), [sources.json](sources.json),
  [recommendation.md](recommendation.md), [status.json](status.json).

Scope notes: GC-correction internals were deliberately not pursued further after the user asked
to prioritize transferable diagnostics. The requested GPT 6.1 Terra delegation was not possible
in this environment, so all work was done directly. Literature hosts (NCBI/PMC/EBI) were blocked
by the session network policy, so claims rest on abstracts and official repositories. No
application code, archived `reference/`/`validation/` file or earlier research file was changed.
