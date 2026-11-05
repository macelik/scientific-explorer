# Recommendation: where the next CNV-method change belongs

Review of `macelik/scientific-explorer@allelic-segmentation` (checkout `6d4f48d`), 2026-11-04.
Evidence labels: **[P]** published, **[L]** local observation on MCF10-A development data
(this review; links point to the saved experiment), **[I]** proposed idea. No result here
establishes copy number at any locus.

## Short answer

Make the next change in **evidence integration and the CN reference**, not in GC/rounding
tuning and not by adopting a new caller. Pair it with one **upstream bookkeeping change**:
keep bins that QC removes as an explicit "unobservable" state.

1. **Add a per-region evidence layer next to the unchanged pyEpi calls** (rank 1 in
   [method_comparison.tsv](method_comparison.tsv)). For each segment it reports: raw-count
   contrast with exposure; where that contrast falls among same-width random windows; whether it
   depends on coverage; whether its size is compatible with a one-copy change; a per-cell mixture
   bound; allelic observability and an empirically calibrated allelic label; and adjacency to
   QC-removed bins. The output is multi-state and includes *insufficient evidence*.
2. **Keep the pre-filter bin grid and QC masks** and carry removed bins as a state instead of a
   silent gap (rank 2).
3. **Obtain a cross-sample technical baseline** (same binning on unrelated samples, plus a
   mappability track) before changing any formula (rank 3). It is the decisive experiment for
   the saved chr7 dips.
4. Only then evaluate a **raw-count, exposure-aware segmentation/assignment** (rank 4) and a
   **joint depth–allele likelihood** (rank 5) against the preserved baseline on independently
   characterized samples. Numbat-multiome is the published comparator (rank 6) once its inputs exist.

## Why

**The saved boundaries are reproducible. Interpreting them is the problem.**
- [L] Cell sampling explains < 0.2% of pseudobulk bin variance; disjoint stratified halves
  correlate at r > 0.99; all 10 saved-window edges that coincide with a breakpoint reappear in
  both halves within ±2 bins (the broad chr9 20–27 Mb window has no breakpoint at its edges) ([exp 05](experiments/05_segmentation_calibration/report.html)).
- [L] On exchangeable nulls, the select-then-test AD procedure accepted 1/300 first splits at
  α = 0.001, so the p-values are not grossly anti-conservative in that setting. This contradicts
  the hypothesis I started with. The real issue is that the null is exchangeable bins, which
  reproducible non-CN structure (accessibility, mappability, GC residue) also rejects. A
  significant boundary therefore says little about CN.
- [L] AD responds to variance as well as mean, and sampling variance follows depth: a pure
  sampling profile still splits at chr7 51.1 and 77.4 Mb, around a broad low-depth region
  ([exp 05](experiments/05_segmentation_calibration/report.html)).

**Depth contrasts at this scale are routinely as large as one-copy changes without being called.**
- [L] Random 25-bin windows that the pipeline leaves entirely uncalled range 0.44–1.72× their
  flanks (5th–95th percentile) in raw counts, and 19 of 226 fall below 0.5×
  ([exp 04](experiments/04_fixed_region_contrasts/report.html)).

**The chr7 dips do not behave like a simple clonal or subclonal one-copy loss.**
- [L] Raw contrasts are 0.13–0.34× the flanks, below the 0.5× floor of a one-copy loss from a
  diploid baseline with unchanged accessibility. They do not depend on coverage, are similar in
  all five clusters (0.11–0.41), and among cells expected to have ≥ 5 target fragments at most
  2.3% have zero. A near-complete loss in a subclone is therefore bounded to about 2% of cells
  ([exp 04](experiments/04_fixed_region_contrasts/report.html)).
- [L] The windows keep heterozygous-site observations (AT 573–3,359 pooled) with no excess
  allelic variance (effective |dBAF| upper bound ≤ 0.05). An imbalance as strong as the one
  measured at the published MYC amplification would be detected (power ≥ 0.85); one half as
  strong usually would not ([exp 06](experiments/06_allelic_information/report.html)).
- [L] chr7 67.4–70.0 Mb is an even plateau at about 0.13× (bin heterogeneity at the 0.3rd
  percentile of random windows) ([exp 07](experiments/07_subwindow_profiles/report.html)).
- Taken together, these are most consistent with a shared **regional baseline** effect
  (accessibility or mappability), possibly combined with a smaller dosage change. They are
  **not** consistent with a clonal one-copy loss explaining the dip by itself. Mappability and a
  cross-sample baseline are the missing data. This remains an interpretation under stated
  assumptions (diploid baseline, AT model), not a demonstration.

**Positive controls show that both evidence types carry information, and where they fail.**
- [P] MCF10-A carries a homozygous CDKN2A/B deletion within a ~4 Mb hemizygous 9p21 loss, MYC
  amplification, and 19q13-qter and 5q13-qter gains ([sources.json](sources.json) S9–S11).
- [L] MYC: raw 2.57× genome density and called gain; strong allelic excess. 5q13-qter: strong
  allelic excess, 44% of bins called gain. The 9p21 zone has the lowest SNP observations per
  fragment (0.20× genome), as expected when clonal LOH removes heterozygous sites. 19q13-qter is
  1.82× in raw counts but 0.93× after GC correction, and is called partly weak loss, while its
  alleles show excess imbalance ([exp 04](experiments/04_fixed_region_contrasts/report.html),
  [exp 06](experiments/06_allelic_information/report.html)).
- [L] The CDKN2A/B bins were **removed by the bin filter** and cannot be called at all
  ([exp 02](experiments/02_zero_clip_audit/report.html)).

**CN assignment is a fixed threshold on the genome median.**
- [L] All 18 scale fits end at the grid's lower edge (s = 0.5 × median), so a loss means
  segment mean < 0.68 × the genome median (Holmes) and ≤ 0.5 × (Watson). Both the reference level
  and the cut-points are set by the same profile. Removing one chromosome changes up to 5% of
  calls elsewhere; changing the estimand changes 20–37%
  ([exp 03](experiments/03_normalization_aggregation/report.html)).

**Preprocessing is reproducible and not the main driver of the target dips.**
- [L] `X` is reproduced bit-for-bit for 3,449 of 3,460 cells from `layers/raw` with the pinned
  source; 17 cell-bins remain unexplained ([exp 01](experiments/01_reconstruct_x/report.html)).
- [L] The GC step does create isolated amplified bins (up to 717×, 12 bins where amplified
  cell-bins make up > 25% of the mean) and coverage-dependent rounding zeros. The saved contrasts,
  however, are coverage-invariant and larger in raw counts than in `X`. GC correction attenuates
  dips in low-GC regions and can remove real gains on GC-rich arms. That is a confounding
  problem a baseline has to address, not something tuning can fix
  ([preprocessing_audit.md](preprocessing_audit.md)).

**Allelic evidence needs empirical calibration.**
- [L] Nominal AT overstates information: MYC's effective |dBAF| is 0.11, far below what nominal
  counts predict for any clonal change, and 31% of uncalled random windows exceed a balanced
  beta-binomial fitted to base-called bins. Independent per-bin allele likelihoods would be
  miscalibrated until SNP-level, phase-block-aware, molecule-aware counts are available.

## Update 2026-11-05: the 2026-10-29 method designs, implemented

The user asked for the saved method records to be tried on raw counts, with an observation model
chosen per modality, then extended to the allelic signal and to major/minor haplotypes. Each item
links its report; all are exploratory and validated by the coordinator, either by independent
recomputation or by a byte-identical rerun.

**Depth (E:1, E:3, B:2, C:3).**
- [L] NB2 with exposure fits cell-level depth; Poisson does not. Quasi-Poisson φ ≈ 3, NB α ≈ 0.01.
  NB reproduces the held-out zero cells of the deep chr7 windows: 30 observed vs 24 NB-predicted vs
  2.5 Poisson-predicted for 8.3–10.7 Mb ([exp 08](experiments/08_exposure_count_models/report.html)).
- [L] Zero inflation is not needed. The bootstrap-"supported" ZINB fits have π ≈ 0.0003 (about one
  cell) and appear just as often in random windows
  ([exp 09](experiments/09_zero_inflation/report.html)).
- [L] Sampling-based z-scores declare almost every region significant. Against the spread of uncalled
  random windows, 27/46 called short losses and 4/21 called short gains remain unusual (BH q < 0.05).
- [L] A raw-count changepoint with sampling-only dispersion (half-split φ = 1.30) produces about 4,200
  segments, and the Pearson-iterated φ is unstable. Any count-likelihood segmentation needs an explicit
  between-bin baseline-variance term ([exp 11](experiments/11_count_changepoint/report.html)). The joint
  HMM confirms this: τ² ≈ 0.26 gives 600 depth runs, versus 17,352 without it
  ([exp 13](experiments/13_joint_hmm/report.html)).
- [L] The baseline choice moves calls. 14/42 called losses and 10/26 called gains change class under at
  least one of seven prespecified baselines. The saved loss windows are stable, except chr9 20–27 under
  a GC-decile baseline ([exp 12](experiments/12_baseline_sensitivity/report.html)).

**Alleles (E:2, observation model, haplotype orientation).**
- [L] Per-modality choice by held-out log score
  ([exp 15](experiments/15_allelic_count_models/report.html)):
  - allele total AT: NB2 with its own, larger α (median ≈ 0.12), conditioned on the region's depth;
    ZINB gains only 0.0004 nats;
  - haplotype split HB | AT: beta-binomial, about 0.030 nats better than binomial; the ρ must be
    region- or cell-level, not the bin-level 0.0038;
  - dependence: AT is not independent of depth (elasticity 0.50).
- [L] Regional beta-binomial on the published controls
  ([exp 10](experiments/10_betabinomial_allele_model/report.html)). Phase-agnostic h: MYC 0.593,
  5q13-qter 0.605. The saved windows are 0.50–0.545 with intervals reaching 0.5. At the empirical ρ,
  24% of uncalled random windows exceed the nominal LR cut, so allelic thresholds must be empirical.
- [L] A phase-orientation HMM recovers most of the GT-oriented loss of signal at gains: 5q |h − 0.5|
  goes from 0.010 to 0.100 and MYC from 0.032 to 0.087. It does not recover 19q or 9p21. 31/164
  segments exceed the random-window 95th percentile. One cluster shows the opposite major haplotype
  (cluster 4, chr8:115–145 Mb), which is a mirrored-event hypothesis
  ([exp 16](experiments/16_major_minor_haplotypes/report.html)).

**Combination: major/minor haplotypes.**
- [L] Explicit allele-specific states (1:0, 0:0, 2:0, 2:1, 2:2, 4:1 at cell fractions 0.1–1) are
  scored against the depth level and the major-haplotype fraction. Each modality has its own
  likelihood and empirical null ([exp 17](experiments/17_allele_specific_states/report.html)).
  - MYC and 5q fit a 2:1 gain with an identified major haplotype.
  - The deep chr7 windows fit **no** modelled state: their depth is below the one-copy floor while
    their alleles are balanced.
  - 22/46 called short losses fit a state. Adding alleles changes the depth-only call in 24 of them.
  - Calls depend on the allelic estimator and its null: 58/67 agree between the phase-agnostic and
    HMM-oriented sources.
- [L] Depth alone carries little information for one-copy changes at this scale. The baseline spread is
  σ_A ≈ 0.39 (log), so a one-copy loss sits about 1.8 σ_A from neutral.

**What this changes.**
1. The evidence layer (stage 1) now has a concrete, tested form:
   - NB2 depth rate ratio, judged against an empirical regional null;
   - NB2 AT conditioned on depth;
   - beta-binomial split oriented by a phase HMM;
   - an allele-specific state score with an explicit "no modelled state fits" outcome.
2. Zero inflation is dropped from the plan for depth.
3. A count-likelihood segmentation (stage 4) must model between-bin baseline variance; sampling-only
   dispersion over-segments by more than an order of magnitude.
4. **Segmentation sensitivity** ([exp 14](experiments/14_spikein_benchmark/report.html),
   [exp 19](experiments/19_masking_two_cut_scan/report.html)) [L]. The benchmark injected events into
   real raw counts (1,800 semi-synthetic realizations) and reran the current AD scan. It recovers clonal
   deep losses (0.1×) in 40–60% of loci at 20–50 bins, but clonal one-copy losses in ≤ 15% and 1.5×
   gains in ≤ 2.5%. Events carried by ≤ 50% of cells are essentially never recovered.
   - Part of the failure is masking: the best single cut in a long segment is elsewhere and gets rejected.
     A two-cut scan rescues deep losses.
   - For one-copy changes, natural windows in the same segment are as extreme as the injected event. A
     better search alone does not fix this; a depth baseline and allele-first segmentation are needed.
5. **Cell-level joint mixture** ([exp 18](experiments/18_cell_level_joint_states/report.html)) [L].
   - Near-clonal 2:1 gains at MYC, 5q13-qter and chr8:115–145 Mb are supported by both depth and alleles.
   - Depth-only cell-level calls flag 7/30 uncalled random windows (23%), and none of those is supported
     by both modalities. Report cell-level clonal/subclonal estimates only where both modalities agree,
     until the depth null includes between-cell accessibility heterogeneity.
6. The consolidated account, with explicit methods and assumptions, is in
   [final_report.html](final_report.html).

## Staged plan with decision criteria

| Stage | Action | Resolves | Success | Failure / stop | Changes the plan if… |
|---|---|---|---|---|---|
| 0 (hygiene, no scientific change) | Persist run parameters, pre-filter grid, QC masks, GC fit diagnostics (predictions ≤ 0, max correction, fallbacks) with every run | Provenance gaps of exp 01–02 | A rerun reproduces X exactly and lists removed bins with reasons | — | — |
| 1 [I] | Implement the region evidence layer as a research tool reading existing outputs (scripts in exp 04/06/07 are the prototype). Labels: depth contrast (percentile vs random windows), dosage compatibility, mixture bound, coverage dependence, allele state (unobserved / insufficient / no excess at observed sites / consistent with imbalance), QC-gap adjacency | Which calls are dosage-compatible, which are not, and which lack evidence | On MCF10-A the labels reproduce the controls (MYC, 5q: depth and allele agree; 9p21: SNP depletion) and flag the chr7 dips as below the one-copy floor | Labels unstable under the seed or to ±20% changes of random-window width | Positive controls do not separate from random windows → allelic or depth layer uninformative at this coverage |
| 2 [I] | Run identical binning/QC on ≥ 2 unrelated scATAC samples, at least one with scDNA/WGS truth and one with an annotated normal population; add a read-length-matched mappability track | Is the chr7 depletion a shared technical/accessibility baseline? | Deep homogeneous dips recurring across unrelated samples at the same loci with normal DNA there → baseline offset justified | Dips sample-specific with DNA-confirmed loss → dosage; then the evidence layer's one-copy floor must accommodate non-diploid baselines | DNA shows deletion in chr7 windows → prioritize sensitivity (raw-count model, rank 4) over specificity |
| 3 [I] | Relax the bin filter on a sample with pre-QC counts; keep removed bins as 'unobservable' | Focal homozygous deletion recall (CDKN2A/B-type events) | Known homozygous deletions become reportable without new false calls in mappability-poor runs (checked against mappability) | Removed runs are dominated by mappability/inaccessible deserts with no truth support | — |
| 4 [I] | Raw-count, exposure-aware segmentation + CN assignment with optional baseline offset; compare with the preserved pyEpi baseline on held-out samples with DNA truth, stratified by event size (100 kb–5 Mb), coverage and reference availability | Measurement and estimand artefacts; calibrated uncertainty | Better detection/false-call trade-off at matched abstention on ≥ 2 held-out samples, including focal events | No gain over baseline, or gains only on MCF10-A | — |
| 5 [I] | Once SNP-level phased counts exist: joint depth–allele model with empirical allelic null; run Numbat-multiome ATAC-bin as comparator where a normal reference exists | CNLoH, allelic support for depth calls | Improved CNLoH/false-loss control on benchmarks without a forced normal cluster | Requires a normal reference that is not available → record limitation | — |

Throughout: do not force a normal cluster from MCF10-A identity or flatness; keep calibration
honest about region selection (the random windows are not neutral truth); do not count caller
agreement as validation.

## What is established, what is not

- **Established (local, reproducible):** provenance of X and its exceptions; QC removal of the
  CDKN2A/B bins; magnitude, coverage invariance and cell homogeneity of the saved contrasts;
  reproducibility of the boundaries; the grid-edge behaviour of the CN scale; the estimand
  sensitivity; allelic observability and attenuation at controls.
- **Not established:** CN state at any saved window; whether the chr7 plateaus are mappability,
  accessibility or dosage; the cause of the 17 historical zeros; any tool's performance beyond
  MCF10-A; calibrated allelic power.
- **Would change this recommendation:** DNA evidence of deletion at the chr7 windows (sensitivity
  would then matter more than specificity); unrelated samples showing no such plateaus with normal
  mappability; or SNP-level data showing the attenuation is a bin-summation artefact (a joint
  model would move up).
