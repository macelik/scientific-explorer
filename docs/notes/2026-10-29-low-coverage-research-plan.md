# Proposed research plan: CNV evidence under low coverage

Status: all six first-pass reviews, the shared descriptive pilot and coordinator
synthesis completed. See `research/2026-10-29-low-coverage/report.md` for findings
and `synthesis/NEXT_EXPERIMENTS.md` within that workspace for proposed follow-up.
Prepared: 2026-10-29. Repo inspected at `3bb9aec`, branch `allelic-segmentation`.

## Objective and corrected sample context

Develop a general scATAC CNV inference approach that distinguishes evidence for
copy loss, allelic imbalance and technical coverage variation, and explicitly
recognizes cases that cannot be resolved with the available data.

**User clarification:** the current dataset is the non-malignant MCF10-A cell
line. The method must not be tailored to this cell line or to one tumour type.
This supersedes the tumour-sample framing in the
[original handoff](2026-10-29-research-handoff-low-coverage-segments.md).
No labelled normal reference is supplied. Cell-line identity is not a
region-level copy-number truth set. We should not search this dataset for an
assumed immune/stromal compartment or force a normal-cell assignment.

Use MCF10-A for initial feasibility and failure-mode investigation. Reserve
generalization claims for independent datasets and orthogonal truth. Normal-cell
inference is an optional capability for future mixed samples; the core research
must also cover homogeneous samples and samples with no suitable normal cells.

## Context verified locally

- Depth H5AD: 3,460 cells by 25,337 retained 100 kb bins; `X` is processed
  depth, `layers/raw` is a separate sparse count matrix. Matrix metadata does
  not record sample identity or the reference build.
- WNN cluster sizes: 1,286 / 833 / 603 / 556 / 182. The integration story
  documents coverage confounding in the representations. Clusters are not
  independently established clones, cell types or normal/malignant labels.
- HA, HB and AT parquet matrices all match the H5AD coordinates and barcode
  set. Their barcode order differs: join by barcode, never column position.
- The active depth defaults use `integration-story/data/genome_pool/`.
  The allelic scan omits bins below `min_at`; masked values are NaN.
- The H5AD uses 1-based inclusive intervals. The phased VCF's chr1 length is
  248,956,422, supporting GRCh38/hg38, but the exact alignment reference and
  provenance still need confirmation before annotation joins.
- Five integration caches exist. Selecting one by modification time is not a
  sufficient provenance check. Pin the data/configuration identity and explicit
  cache path, or recompute small summaries directly from the pinned inputs.
- The integrated object contains depth and haplotype bin modalities. A
  peak-by-cell matrix, fragment file and gene-activity matrix have not yet been
  located for annotation/representation methods; they must not be assumed.

The repo's current CN values are categorical encodings (`0 = loss`, `1 = base`,
etc.), not absolute copy numbers. In particular, a call encoded as zero is not
a homozygous-deletion label. See `server/data.py` and the existing cohort README.

## Research design

Recommended: shared inputs and pilot cells, independent category reviews, then
central synthesis and a small number of follow-up experiments.

Alternatives considered:

1. Literature only: cheapest initial computation, but input incompatibilities
   and inadequate allelic power can remain hidden.
2. Independent 300-cell samples per agent: easy to dispatch, but sampling
   differences complicate comparison and repeat data preparation.
3. **Shared 300-cell pilot plus existing full-cohort summaries:** comparable,
   inexpensive feasibility checks, with full-cohort confirmation only where
   the pilot identifies a concrete question. This is the recommendation.

### Shared preparation, owned by the coordinating agent

Before any scientific pilot:

1. Freeze a manifest of resolved input paths, file identities, barcode order,
   genomic intervals, cluster labels, repo commit and dirty state. Record
   content hashes for the small exported pilot tables and exact inputs used.
   Validate `AT = HA + HB` for every selected barcode. Audit whether upstream
   counts represent unique fragments, reads or SNP observations; depth and AT
   are coupled measurements, not automatically independent evidence.
2. Use the existing `../default_300cells/cells.csv` membership as the initial
   pilot. It has 100 cells per raw-coverage tercile, 300 unique matched IDs,
   and cluster counts 106 / 72 / 56 / 46 / 20. Preserve its exact membership
   and hash. The historical sampling seed was not verified in this inspection;
   do not invent one or describe this as a fresh random draw.
3. Record the cluster-by-coverage table. Some combinations are absent, so this
   sample does not support every within-coverage cluster comparison. If a
   category needs more cells in a small group, make a separately labelled,
   seeded supplementary sample; do not replace the shared panel silently.
4. Freeze candidate regions from the handoff and the corresponding unmerged
   genome-pool segments: chr7 8.3–10.7, 51.4–53.5, 67.4–70.0 and 70.0–73.3 Mb;
   chr9 approximately 20–27 Mb. Include immediate flanks and matched comparison
   regions. These are hypotheses, not positive/negative ground truth labels.
5. Export raw depth, processed depth, HA/HB/AT, cluster and library-size
   covariates for these regions and controls. Record retained-bin count,
   callable span and genomic gaps. Never rerun cohort filtering on 300 cells:
   that would change the retained grid. Use raw counts for count/zero models.
6. Keep full-cohort pseudobulk and 300-cell outputs distinct. The full-cohort
   pseudobulks provide context; they are not independent validation of a
   300-cell result. Do not discard useful allelic counts merely to enforce a
   300-cell limit on already cached aggregate summaries.

Use seed `20261029` for new pilot sampling or simulations, with named child
streams and recorded RNG/library versions. The seed does not retroactively
describe the legacy 300-cell selection. Save every new selected barcode list.

### Scientific questions shared by all tasks

- Missing allelic observations, evidence compatible with balance and evidence
  against a specified imbalance are different outcomes. Specify the alternative,
  uncertainty and attainable power before claiming absence of LOH.
- Pooling can obscure opposing allelic states. Check phase sets, SNP selection,
  allele-specific mapping bias and possible mixture effects before interpreting
  segment-wide HA/HB sums. Global HA/HB orientation must be justified.
- Observed heterozygous-SNP scarcity can be biological or technical. The source
  and filtering of the SNP list matter, especially without a matched normal.
- Compare AT depletion with **raw** depth depletion as well as processed X.
  Their differing normalization and counting units may affect the handoff's
  ratio comparison. Do not assume that its ratio estimates alone identify a
  technical mechanism.
- Zero counts are observations. A structural-zero component is a model
  hypothesis. Excess zeros or apparent bimodality alone do not identify a
  deletion; low mean, dispersion, cell state and library size need controls.
- A shared dip across WNN groups can reflect a shared event or shared bias.
  Groups and the all-cells aggregate are not independent biological replicates.
- Smoothed/imputed tracks cannot supply independent evidence in unobserved
  regions. Avoid filling missing allelic values with artificial balance.

## Six research-agent briefs

All agents use **GPT 6.1 Terra**, tool identifier `gpt-6.1-terra`, with fresh
context (`fork_turns = none`). Each receives this plan, the handoff, the user
clarification, repository constraints and its own task. No recursive delegation.
Use a compact initial screening pass, followed by deeper review of at most
three promising candidates per task. Stop when the requested decision is
supported, rather than accumulating papers without a decision.

### A. Existing CNV methods: what can be reused?

Screen epiAneufinder/pyEpiAneufinder, copyscAT, Alleloscope, Numbat, CopyKAT,
inferCNV, CaSpER, HoneyBADGER, XClone, CHISEL and SIGNALS. Check exact method
identity where a name is ambiguous.

For each, establish assay, data prerequisites, reference/baseline assumptions,
low-coverage emission or filtering model, allelic pooling, supported copy states,
homozygous-deletion handling, uncertainty output and tested resolution. Separate
direct scATAC applicability from an idea transferred from scRNA or scDNA.
Inspect methods/supplements and relevant source code, including actual licenses.
Do not equate a generic loss state with a homozygous-deletion capability.

Deliver a comparison matrix and top three reusable components, with input gaps
and reasons for rejecting apparent candidates. Optional pilot: check whether a
leading method's input requirements and count ranges fit the shared tables;
do not install and run every complete caller.

### B. Joint segmentation and uncertainty under missing evidence

Research count-based change points, heteroscedastic models, joint depth/allele
HMMs, FACETS/ASCAT-style segmentation, Bayesian change-point approaches, and
selective prediction/abstention. Include masked signal/image methods as a
secondary search, asking what transfers to genomic counts and what does not.

Determine how each approach distinguishes a genuinely low observed count from
missing data, propagates genomic gaps, allows a depth boundary when alleles are
unobserved, and avoids treating depth/AT as independent assays. Check whether
the model requires a reliable normal baseline or diploid anchor.

Deliver two or three candidate model designs with assumptions, identifiability
limits and proposed output semantics. Optional pilot: preserve measured
coverage and compare the evidence available at fixed candidate boundaries;
evaluate a small seeded sensitivity example without replacing the AD pipeline.
Any synthetic result describes its generating assumptions, not biological truth.

### C. Reference strategy with and without suitable normal cells

Research optional internal normal identification, cell-type annotation/atlas
transfer (including Signac/ArchR), robust or iterative baseline estimation,
public panels of normals, and the relevant CopyKAT/Numbat/inferCNV assumptions
identified by task A. Avoid duplicating A's entire caller survey.

Cover mixed tissues, homogeneous cell lines, highly altered samples and samples
containing no suitable normals. Distinguish cell type from malignancy and
relative CN state from absolute ploidy. Identify circular procedures that select
the reference using the same focal CN signal subsequently evaluated.

Deliver a decision tree that can return "no reliable reference identified",
with prerequisites and cross-sample failure modes. Optional pilot: assess
baseline sensitivity using existing depth/allelic summaries. Do not designate
a WNN cluster as normal, or MCF10-A as uniformly diploid, without independent
evidence. Annotation pilots require compatible inputs first.

### D. Representation and foundation models: practical added value

Verify the papers and implementations for atacformer, ChromFound,
EpiFoundation, scBasset, SCALE and PeakVI/scVI-family approaches. Distinguish
preprints, published work and maintained code. Check input feature space,
pretraining data, supported species/builds, pretrained weights, licenses,
compute requirements and out-of-distribution uncertainty.

Ask whether the output is an embedding, cell-type label, accessibility estimate
or an actual calibrated CNV-related prediction. Assess whether CN differences
could drive the representation and whether normalization could remove useful
dosage information. A tumour/normal label is not a substitute for CN inference.

Deliver a ranked feasibility table and an explicit decision on whether any
model earns a later experiment. Optional pilot: input compatibility only, or
one lightweight pretrained inference if data and weights already exist.
No model training, large atlas downloads or forced peak-model input from the
100 kb depth matrix in this first pass.

### E. Per-cell loss evidence, zero counts and allelic power

Research Poisson/NB baselines, hurdle or zero-inflated alternatives, latent
mixtures, beta-binomial allelic models and suitable model-checking/calibration
procedures. Do not assume zero inflation or observable bimodality in advance.

Use raw regional counts with library-size/exposure terms and matched-region
controls; compare remaining signal after coverage and cluster adjustment.
Audit whether a focal event diluted into a wider segment could be missed.
For allelic tests, require adequate phase consistency and account for sampling
noise, overdispersion, contamination and mixtures. Avoid uncalibrated mixture
likelihood-ratio p-values and multiple testing across chosen regions/models.

Deliver an identifiability/power assessment for the candidate regions, proposed
comparisons between technical depletion, hemizygous loss and biallelic loss,
and the evidence needed to separate them. Optional pilot: one shared-panel
region/count summary and a small calibrated model comparison if justified.
Report effects and uncertainty; inability to distinguish alternatives is a
valid finding. Selected depth segments require selection-aware interpretation.

### F. Technical priors, SNP ascertainment and coordinate provenance

Research blacklist/gap masks, mappability, GC, cell-type accessibility priors,
SNP density and phasing quality. Distinguish reference SNP density from SNPs
observed as heterozygous in this sample. Inspect the upstream SNP ascertainment
and counting pipeline before interpreting low AT as independent dropout.

Confirm genome-build and interval conventions before any overlap analysis.
Check what the existing bin filtering has already removed and what cannot be
assessed from the retained grid. Explore soft weights and uncertainty flags as
well as exclusion masks; assess the risk of masking real events systematically.

Deliver an annotation/prior matrix with available sources, build, license,
limitations and required joins. Optional pilot: annotate the fixed regions from
already available local resources and summarize raw-depth versus AT depletion.
Do not download large reference resources until a specific candidate needs them.

## Execution order, output contract and cost controls

Preparation precedes scientific pilots. Wave 1 runs A, E and F concurrently.
The coordinator consolidates input limitations and useful evidence. Wave 2 runs
B, C and D, each receiving that compact consolidated brief. Maximum concurrency
is three research agents plus the coordinator. All six use Terra explicitly;
this is a model selection, not a verified pricing claim.

Each agent saves its results in its assigned directory under
`research/2026-10-29-low-coverage/` and returns links to those files. The workspace
README and `fields.yaml` define the required formats: `findings.md`,
`methods.json`, `sources.json`, `pilot.md` and `status.json`. Agents save progress
as they work; chat-only findings do not count as completion. Shared inputs and
cross-task synthesis have separate coordinator-owned directories.

Each agent records:

1. A short decision summary and screened-method table.
2. Claim-level primary citations: paper/method section and official code/docs
   where relevant. Record version/date and verify licenses separately.
3. Top candidates with inputs, outputs, reference assumptions, failure modes,
   fit to the current data, approximate computational burden and unresolved facts.
4. Clear evidence labels: published/evaluated for this assay; published in
   another assay; untested transfer idea; or exploratory local diagnostic.
5. Either one justified pilot with inputs, seed, exact command, environment,
   outputs and limitations, or a concrete reason to skip computation.
6. One recommended next experiment with a measurable decision criterion.

Papers and public code can be searched openly; local cells/data stay local.
Research reports are saved under the new task-specific directory
`research/2026-10-29-low-coverage/`. Executable scratch pilots, if used, go in
isolated `/tmp` directories; reports preserve their exact commands/source text
for reproducibility. Agents do not edit application code, source datasets,
existing caches, `reference/` or archived `validation/`. They do not touch the
user's server on port 8766. Heavy computation is serialized if needed so that
three simultaneous agents do not each launch unrestricted worker pools.

The 300-cell panel is a feasibility check, not a sufficient benchmark of rare
events, normal-cell recovery or generalization. Agents share its membership;
they do not repeatedly resample until a desired conclusion appears. Targeted
checks on all 3,460 cells may be cheaper and more informative than more pilots;
the coordinator decides this from the first-pass evidence.

## Synthesis and next-experiment selection

The coordinator audits the leading claims and licenses against primary sources,
resolves contradictions across agents and produces a ranked shortlist. Rank
first by relevance and evidence, then input feasibility, honest uncertainty,
computational cost and integration effort. Preserve reasons for exclusions.

Select a small sequence of experiments that can change a decision, such as:

- Do raw depth and SNP ascertainment explain the relative AT depletion?
- Is any specified allelic-loss alternative distinguishable at observed coverage?
- Does a raw-count model require an additional zero component after controls,
  and would such a component distinguish deletion from inaccessible chromatin?
- Can a method express missing allelic evidence without suppressing a plausible
  depth-only event?
- Which reference strategy remains usable when no suitable normal cells exist?

Define cross-dataset evaluation before tuning: independently characterized CN
events from orthogonal DNA measurements where possible; datasets with and
without reference cells; varying coverage, event size, ploidy and mixture;
technical low-coverage regions; held-out samples and parameter settings.
Measure event/boundary accuracy, false calls, allele-specific-state performance,
uncertainty calibration/abstention and runtime, at appropriate cell/region/sample
levels. Do not treat correlated bins as independent replicates.

MCF10-A development results, agreement between computational callers and
simulated perturbations are separate evidence categories. None substitutes for
independent biological validation. Proposed changes to the app or its scientific
calculations require a later implementation plan and the repository's validation
workflow; this document does not change scientific behavior.

## Preliminary source checks informing these briefs

These are starting points, not the completed literature review:

- [Alleloscope paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9127226/): directly
  relevant scATAC work, but its recommended scATAC segmentation uses matched
  DNA evidence. This motivates checking prerequisites rather than assuming a
  caller is directly usable here.
- [Numbat paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10289836/): explains
  dilution by pooling genetically distinct cells and iterative clonal inference.
- [Numbat official workflow](https://kharchenkolab.github.io/numbat/articles/numbat.html):
  documents reduced heterozygous-SNP availability in clonal LOH without matched
  normal cells, motivating the SNP-ascertainment audit.
- [PeakVI paper](https://doi.org/10.1016/j.crmeth.2022.100182): models binary
  region detection with cell/region factors; this is a reason to distinguish
  representation/observation models from raw regional-count CN inference.

## Outstanding inputs

- Confirm exact reference FASTA/build provenance; VCF contig length is supporting
  evidence, not a complete alignment manifest.
- Locate any independent characterization of the particular MCF10-A stock and
  candidate events, without treating canonical cell-line annotations as truth
  for this preparation.
- Establish whether peak/fragment/gene-activity inputs and independent benchmark
  datasets are already available. These affect optional tasks, not whether the
  core method research can proceed.
