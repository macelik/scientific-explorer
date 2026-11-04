# Proposed next experiments

Status: proposed research protocol, not executed except for the saved
descriptive pilot. These experiments do not replace existing AD statistics,
permutation settings or CN assignment. MCF10-A remains development data;
generalization requires held-out samples and independent DNA evidence.

## 1. Observation and provenance audit

**Question:** which observations can support a count or allelic likelihood?

Inputs: original fragment/BAM preparation records, actual depth-generation
configuration, pileup invocation and candidate SNP list, genotype/phasing inputs,
aligned matrices, H5AD retained-bin grid and available build metadata.

Tasks:

- Reconcile fragment/read/UMI/site-observation definitions. For a few fixed
  candidate and comparison windows, trace observations from source to bin sum
  if original data are available. Record whether molecules can contribute to
  multiple SNP sites or bins. Do not assume AT/depth is a sampling probability.
- Determine what makes a site eligible, called heterozygous, phased and finally
  counted. Separate the numbers of candidate sites, retained sites, observed
  sites and allele observations. Document how clonal LOH could affect this chain.
- Determine what evidence supports orientation across SNPs and phase blocks.
  Absence of PS tags alone does not prove incorrect phasing; neither does a
  chromosome-wide GT orientation prove it is error-free.
- Verify the actual alignment reference and coordinate transformations. Explain
  missing bins without treating every retained-grid gap as biological deletion.

Outputs: `observation_units.json`, `site_audit.tsv`, `phase_audit.tsv`,
`coordinate_audit.tsv`, exact input identities and a reproducible command record
in a new experiment directory. Do not modify the frozen first-pass files.

**Decision:** classify each input as suitable, suitable with a stated
measurement model, or unresolved for the intended likelihood. The raw-count
summary can remain descriptive when provenance is incomplete. An unresolved
allelic unit/phase question prevents calibrated allelic claims, not all further
depth exploration.

## 2. Fixed-region count baseline and evidence sufficiency

**Question:** does coverage-adjusted regional variation require more than an
ordinary sparse-count model, and what allelic alternatives could be distinguished?

Start with the same six candidate windows, explicit subwindow bins and both
flanks. Add comparison regions chosen using predeclared callable-span,
GC/mappability and coverage criteria; none is presumed neutral merely because
it was not called a loss. Use the frozen 300 cells for development and targeted
all-cell checks only when they answer a defined question. The full cohort is
not an independent validation sample.

**Count arm:** compare a Poisson sampling baseline with an overdispersed
NB alternative on raw counts, with per-cell library exposure excluding the
target and documented region/covariate terms. Treat technical/accessibility
baseline and CN as potentially confounded. Compare predicted versus observed
zero fractions, positive-count quantiles and dispersion across coverage strata;
show sensitivity with and without WNN covariates. Assess narrower windows so
that a focal change is not erased by pooling across intact flanks.

Use held-out predictive checks for model adequacy. These do not solve the
selection bias from choosing regions on the full cohort. Inferential claims
require independently specified boundaries or calibration that repeats the
selection process. Do not treat all correlated bins as independent samples.

**Allelic arm:** after experiment 1, retain phase blocks and use unique or
explicitly modeled observations. Report total information, interval widths
and sensitivity to bias, dispersion and phase errors. A proposed sensitivity
grid is haplotype fractions 0.5 (null), 0.6, 0.75 and 0.9; beta-binomial
intraclass correlation 0, 0.02 and 0.10; and documented phase-switch scenarios.
These are illustrative model scenarios, not asserted biological CN states.
Translate a biological hypothesis into an allele fraction only under explicit
purity/accessibility/copy assumptions.

For a prospective simulation, use seed 20261029 with named child streams and
save the generator, versions and every parameter. An initial 1,000 replicates
per small scenario gives a development estimate with Monte Carlo uncertainty;
increase only if that uncertainty changes a decision. A proposed diagnostic
target is power at least 0.8 at size 0.05 for the prespecified alternative,
with uncertainty intervals and calibration reported. These are proposed
research targets, not app defaults or validated operating thresholds. Define
the family of tested regions/alternatives and a multiplicity strategy before
claiming discoveries.

Outputs: raw-count predictive checks, coverage-stratified effect summaries,
phase-block support and power/interval sensitivity tables. Keep evidence labels
separate: allele unobserved, insufficient information, compatible with balance,
evidence against a specified imbalance, or evidence consistent with imbalance.

**Decision:** prioritize further allelic testing only where assumptions and
power support the specified comparison. Consider an extra-zero model only
after simpler model misfit is demonstrated. Even a reproducible extra-zero
component does not distinguish deletion from inaccessible chromatin without
additional evidence.

## 3. Published caller compatibility and independent benchmarks

**Question:** does an existing joint method improve detection without hiding
the no-reference and small-event failure modes?

Pin a Numbat code revision and inspect the actual CLI because its vignette's
mode prose and example command differ. Verify raw-bin counting compatibility,
phase/SNP fields, missing-emission behavior and the normalized reference.
Do not install or rerun on the current sample until these inputs are concrete.
Use the current pyEpi output as the preserved comparison baseline.

Candidate benchmark leads:

- Numbat-multiome's MM/CLL-to-RS evaluation and
  [analysis repository](https://github.com/getzlab/Numbat-multiome_Analysis)
  provide a reproduction lead; sample-level input access and mapping to DNA
  calls need verification. The event-level focal definitions must not be
  mistaken for comprehensive validation of 100 kb–5 Mb events.
  [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12531993/).
- Alleloscope's data statement lists scATAC PRJNA674903, cell-line scDNA
  PRJNA498809 and public ATAC/exome accessions PRJNA532774/PRJNA533341. These
  are discovery leads, not a verified ready-to-download matched benchmark;
  resolve specimen/cell-line relationships and truth granularity first.
  Some patient DNA datasets are listed under dbGaP accessions.
  [Primary study, data availability](https://pmc.ncbi.nlm.nih.gov/articles/PMC9127226/).

Construct a small benchmark manifest containing assay, specimen relationship,
normal/reference availability, event sizes, orthogonal truth provenance,
coverage, build, access conditions and held-out role. Include more than one
biological setting, reference-present and reference-absent conditions, and
known focal events as well as broad alterations. Obtain calibrated truth for
zero-copy events specifically before claiming homozygous-deletion performance.

Evaluate event detection and boundaries, state accuracy where truth supports
it, false calls near technical low-coverage regions, and uncertainty versus
fraction of regions receiving a call. Stratify by event size, coverage, mixture
and reference availability. Use sample-level holdouts; retain ambiguous DNA
truth as ambiguous. Do not use caller agreement as independent truth.

**Decision:** advance a method only if it improves the prespecified detection/
error trade-off over the preserved baseline across held-out settings. When a
method cannot run without a reliable reference, record that limitation rather
than assigning a synthetic normal cluster.

## 4. Add priors or representations only for a demonstrated gap

Apply gap/blacklist/mappability annotations as sensitivity flags first. Compare
unflagged and flagged-region errors on the independent benchmark; quantify
lost real-event sensitivity before considering exclusion. Match build and
read-length assumptions and distinguish code licenses from data-track terms.

The [adult accessibility atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC8664161/)
may inform expected cell-state accessibility after compatible rebinning.
Atlas membership alone does not establish a copy-neutral reference. Test
donor, cell-type and chemistry mismatch before reference adoption.

If annotation or residual confounding remains a specific problem and authentic
peak/fragment inputs exist, compare a simple peak-based representation with
Atacformer. Hold out target CN regions or otherwise assess leakage; evaluate
added value on independent DNA truth and across coverage. PEAKVI/SCALE and
larger foundation models remain optional later alternatives. Do not train
models merely because they were included in the literature screen.

## Deliverables and stopping rules

Each experiment gets a new directory with question, frozen manifest, exact
source/commands, parameters/seed, small output tables, effect/uncertainty
summary and a decision. Stop an individual branch when required information
is absent or the proposed method cannot discriminate the hypotheses under
its assumptions; retain an unresolved result. Move successful branches to a
separate implementation design with the repository's scientific-change
documentation and tests. Those application changes are outside this research pass.
