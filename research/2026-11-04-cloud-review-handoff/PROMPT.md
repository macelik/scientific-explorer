# Cloud-agent brief: choose the next CNV-method direction

You are reviewing `macelik/scientific-explorer`, branch
`allelic-segmentation`. Work from the latest handoff on this branch and record
the exact checkout commit. The earlier research baseline is commit
`994fb3afcf4caa93f6b52af7267d3668d469eec2`. Do not assume that `main` contains
these research files.

## Objective

Assess the whole scientific workflow and recommend what fits our CNV-detection
problem best. Evaluate whether the most useful next change belongs upstream
(QC, GC correction, rounding, zero clipping, normalization), in segmentation,
in CN assignment, or in how depth and allelic evidence are combined. Support
your recommendation with source inspection, focused reproducible diagnostics,
and primary literature. Do not assume that adding a newer caller is the best
answer, or that our existing preprocessing is correct.

We are developing a general haplotype-aware CNV tool for single-cell ATAC.
MCF10-A, a nonmalignant cell line, is the current development dataset. It is
not a region-level diploid/CN truth set. Recommendations must transfer across
samples, cell states and reference availability.

After QC and preprocessing we average per-cell processed depth within WNN
clusters, segment the pseudobulk with an Anderson–Darling (AD) scan, and assign
relative CN states. Some short depth dips become loss segments while allelic
coverage in those regions is very sparse. A coverage change can arise from
dosage, accessibility, technical effects or sampling. Increases can be
ambiguous too. We need to preserve real focal changes, avoid unsupported calls,
and explicitly represent insufficient evidence.

## Read first, then trace the actual code

1. `AGENTS.md`, `docs/SCIENTIFIC_WORKFLOW.md`, `docs/ARCHITECTURE.md`,
   `docs/EXTENDING.md`, `docs/VALIDATION.md`, `docs/LICENSING.md`.
2. `research/2026-10-29-low-coverage/README.md`, `report.md`,
   `synthesis/NEXT_EXPERIMENTS.md`, `synthesis/uncertainty_register.tsv`,
   and the six agent folders' findings/method/source records.
3. `research/2026-11-04-cloud-review-handoff/README.md`,
   `DATA_MANIFEST.json`, `PACKAGE.json`, and `UPSTREAM_AUDIT_NOTES.md`.
4. Inspect the rest of the repository to map code, scientific assumptions,
   tests, archived evidence and runtime dependencies. Follow relevant call
   paths rather than treating the literature summaries as authoritative.

The original `docs/notes/2026-10-29-research-handoff-low-coverage-segments.md`
is historical context. Its claim that an absent LOH signature rules out a
clonal one-copy loss was revised by the later synthesis. Treat historical
labels such as “real deletion” or “technical artifact” in code/plots as
hypotheses unless independent supporting data are supplied.

## Establish what you can actually reproduce

The Git repository alone does not contain the full AnnData/MuData matrices,
HA/HB/AT parquet files, or the external pyEpi package. Follow the handoff README
to obtain the separately supplied data package and pinned upstream source.
Check checksums and input identities before fitting anything. If the data
asset has not been published or is inaccessible, continue code/literature
review and clearly label numerical work as unavailable. Do not invent results
or replace a missing real dataset with unlabelled synthetic data.

The supplied AnnData has 3,460 cells and 25,337 retained 100 kb bins. `X` is
processed depth; `layers/raw` contains pre-GC fragment counts **on the retained,
post-QC cells and bins**. It does not restore cells/bins discarded by QC.
The MuData depth modality has a different transformation from the explorer's
AnnData `X`; do not interchange them. Join allelic matrices by barcode name
and verify genomic coordinates. Their barcode order differs from AnnData.

Record what is available for each stage: source fragments/BAM, original full
bin grid, pre-QC matrix, filtering masks, run parameters, intermediate GC fits,
pre-clipping values, processed counts, SNP selection and phase information.
Distinguish verified historical execution from the behavior of the currently
available source. Missing pre-QC counts limit experiments with relaxed QC;
they do not prevent audits on the retained matrix.

## Workstream 1: preprocessing and measurement audit

Trace the exact inputs, outputs, formula, order and parameters of each stage.
In the available upstream source, start with `pyEpiAneufinder.py`,
`render_fragments.py`, `assign_somy.py`, and `get_breakpoints.py`. In this repo,
trace `server/config.py`, `server/data.py`, `server/integration.py`,
`server/pyepi_adapter.py`, `server/allelic_data.py`, `server/allelic_scan.py`,
`server/cn_recall.py`, and `tools/segment_depth_genome_pool.py`.

The GC/clipping/normalization implementation is not vendored in this workbench.
Read these exact original-repository locations before proposing changes:

- [Per-cell GC LOESS helper](https://github.com/colomemaria/pyEpiAneufinder/blob/e8305ea2162634d07e92a3bb401df672e48f217b/pyEpiAneufinder/render_fragments.py#L155-L164).
- [GC correction → rounding → clip negatives to zero → normalize each cell to 100,000 → store X](https://github.com/colomemaria/pyEpiAneufinder/blob/e8305ea2162634d07e92a3bb401df672e48f217b/pyEpiAneufinder/pyEpiAneufinder.py#L252-L271).

Both links are pinned to public commit `e8305ea2162634d07e92a3bb401df672e48f217b`,
which matches the local package. Do not substitute today's upstream default
branch when reconstructing this baseline.

Investigate these concrete questions:

- QC: how do bin nonzero-frequency filtering, cell coverage filtering and
  count/dispersion outlier filters affect retained regions? Could a true
  event influence its own inclusion? Which masks or original observations
  are needed to evaluate this? Check sparse explicit-zero handling. Also
  distinguish original read/fragment and cell-quality filters from these
  later matrix filters; audit them where the necessary records exist.
- GC: inspect per-cell LOESS, the correction denominator, boundary behavior,
  very small/zero/negative predictions and fit failures. Separate observed
  problems from plausible failure modes. Check how broad events or sparse
  cells affect the fit. Source defaults are not proof of historical settings.
- Rounding and clipping: quantify raw-positive values that become zero at
  each reproducible step, negative/nonfinite values before clipping, and
  their distribution across GC, coverage and genomic regions. Capture floating
  intermediates before integer casting so invalid values cannot be hidden by
  casting and subsequent clipping. Do not treat
  allowing negative fragment counts as a valid count-model solution.
- Normalization: trace scaling to a fixed per-cell total and any subsequent
  genome-wide scaling during CN assignment. Assess whether broad events or
  preprocessing artifacts change the denominator and therefore contrasts
  elsewhere. Separate normalization effects from real dosage changes.
- Aggregation: compare averaging processed cells with appropriately defined
  aggregation of raw counts plus exposure. These have different estimands;
  do not silently swap them. Explain weighting and cell-mixture effects.
- Measurement: depth counts and selected SNP-site AT are related. Audit what
  each counts, SNP ascertainment, phase blocks and repeated molecule
  contributions before using independent likelihoods or pooled allele tests.

First attempt to reconstruct the existing `X` from retained `layers/raw` and
GC values using the pinned source. Report agreement and any mismatches, with
dependency versions and parameters. An unsuccessful reconstruction is a
provenance finding, not permission to claim the historical pipeline used
your substitute. Do not overwrite the original matrix.

## Workstream 2: focused diagnostic comparisons

Use the saved 300-cell membership for initial diagnostics; it is a fixed
coverage-stratified panel, not an independent validation set. Preserve its
membership and do not silently draw another 300 cells. Use the full retained
matrix for inexpensive coverage/GC summaries when useful. Every stochastic
experiment must record a seed (suggested new-experiment seed: `20261104`),
generator, parameters and environment.

Start with the six saved chr7/chr9 windows, their subwindows and gaps, plus
prespecified comparison regions. Comparisons must include increases as well
as decreases before making a general gain/loss claim. Existing calls and
flanks are not biological truth. Avoid choosing a favorable control after
looking at the result.

Isolate causes: keep regions fixed when evaluating transformed measurements;
keep segmentation settings fixed when evaluating preprocessing; then compare
boundary and CN-assignment consequences. Do not change all stages together
and attribute an improvement to one of them. Select a small set of justified
alternatives after the audit instead of running a large parameter search.

Useful outputs include raw-versus-processed profiles, newly introduced zeros,
finite-value checks, normalization factors, coverage-stratified contrasts,
boundary stability, model predictive checks and allele information. Report
both affected and unaffected regions. Fewer segments, smoother profiles or
agreement with another caller alone do not establish improved CNV detection.
Calibration must account for data-dependent region/boundary selection.

## Workstream 3: choose a method strategy

Critically assess the saved shortlist: the current pyEpi baseline with better
diagnostics/preprocessing; raw-count segmentation with exposure and explicit
missingness; Numbat-multiome ATAC-bin; regional allele models such as
Alleloscope; and technical/accessibility priors. Representation/foundation
models are optional candidates, not presumed CN callers.

Verify leading candidates against current primary papers and official code.
Assess actual assay/input compatibility, reference assumptions, missing-data
behavior, focal-event resolution, compute cost, exact version and license.
Do not describe the original RNA-only Numbat model as the entire current
Numbat-multiome method. Do not claim an unexecuted tool outperforms the baseline.
Do not force a normal cluster from MCF10-A identity, cell-type annotations,
WNN labels or apparent flatness of the same CN signal being evaluated.

Rank strategies by what they resolve, evidence strength, required inputs,
identifiability, implementation effort and cross-sample validation needs.
Explain which uncertainty each proposed next experiment would resolve, what
would count as failure, and what evidence would change your recommendation.

## Deliverables and scope

Write new work under `research/<actual-date>-cloud-method-review/`. Keep the
earlier research, original data, `reference/` and `validation/` intact.

Required outputs:

- `README.md`: concise recommendation, what is established and what remains
  unknown; link every substantive result.
- `input_audit.json`: availability, hashes, dimensions, units, alignment,
  upstream source identity, missing prerequisites and historical uncertainties.
- `preprocessing_audit.md`: stage-by-stage formulas/code locators, failure
  hypotheses, observed diagnostics and historical-reproduction results.
- `method_comparison.tsv`: alternatives, evidence, prerequisites, limitations,
  cost and reason for their rank.
- `experiments/`: runnable scripts, exact commands, parameters/seeds,
  environment, small tables/figures and interpretation for each executed check.
- `sources.json`: primary references with claim locators and access dates.
- `recommendation.md`: whether to change upstream processing, segmentation,
  evidence integration or a combination; a staged plan with decision criteria.
- `status.json`: completed, skipped and blocked work, with concrete reasons.

You may implement isolated research scripts and small reproducible experiments.
This assignment ends in a reviewed recommendation and concrete next steps;
adopting a new production formula or rewriting the application requires a
separate, explicit scientific decision. Follow repository test requirements
for any implementation changes and explain old behavior, proposed behavior,
reason and validation. Do not run archive-writing validation commands against
existing provenance records; direct new checks to the new research folder.

If you delegate independent work, use GPT 6.1 Terra as requested by the user,
give each agent a bounded brief and its own folder, and independently validate
its saved evidence. If that model is unavailable, continue directly and report
the limitation rather than silently selecting another model. Final conclusions
must distinguish published evidence, local observations and proposed ideas.
