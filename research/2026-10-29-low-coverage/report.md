# Low-coverage CNV research: findings and recommended direction

Completed first-pass review, 2026-10-29. All six research tasks ran with
`gpt-6.1-terra`; the coordinator reviewed their saved outputs and prepared this
synthesis. This is a research assessment and descriptive pilot, not a validated
replacement CNV pipeline.

## Recommendation

Retain the current depth pipeline as the comparison baseline. First establish
what the raw and allelic observations actually count, then evaluate regional
depth depletion and allelic information separately. Investigate Numbat-multiome
as the leading published scATAC comparator. A future combined model should
retain depth evidence when alleles are unavailable and explicitly report which
biological explanations remain unresolved.

Do not begin by tuning a zero-inflation threshold, forcing a normal-cell group,
or training a foundation model. The current evidence does not justify those
choices. The development sample is MCF10-A; its non-malignant identity does not
establish region-level copy number. The method must be evaluated across samples
with and without usable normal references.

## Six saved reviews

1. [A — Existing CNV methods](a-cnv-methods/findings.md): caller comparison,
   assay fit, inputs, state definitions and software availability.
2. [B — Joint segmentation](b-joint-segmentation/findings.md): count models,
   missing evidence, uncertainty and adjacent-field transfers.
3. [C — Reference strategies](c-reference-strategies/findings.md): internal,
   external and unanchored analysis, including a concrete public atlas candidate.
4. [D — Representation models](d-representation-models/findings.md): Atacformer,
   ChromFound, EpiFoundation, PEAKVI, SCALE and scBasset compatibility.
5. [E — Per-cell evidence](e-per-cell-evidence/findings.md): zero counts,
   allelic power and the shared 300-cell descriptive pilot.
6. [F — Technical priors](f-technical-priors/findings.md): ascertainment, counting
   units, phasing, gaps, blacklist and mappability.

Each folder includes structured methods, source records, pilot details and
status. The [combined inventory](synthesis/methods_inventory.json) preserves
every method record and its originating task. Unknown facts remain explicit.

## What changes our interpretation of the handoff

**Low AT is not automatically technical dropout.** The inspected builder sums
DP over selected phased heterozygous SNP sites; AT is therefore not established
as unique-fragment coverage. A molecule can potentially contribute at multiple
sites. Its exact historical run and counting filters remain unverified. The
builder also does not preserve phase-set metadata in the exported matrices.
The phasing header supports hg38, while exact alignment-reference provenance
remains unresolved. See [F's audit and source locators](f-technical-priors/findings.md).

**Absence of measured imbalance is not sufficient to exclude loss.** In
addition to sparse observations, ascertainment can remove informative sites:
Numbat's documentation specifically discusses reduced heterozygous-SNP
availability in clonal LOH without matched normal data. This is a plausible
failure mode to investigate here, not a diagnosis of the candidate regions.
[Official workflow](https://kharchenkolab.github.io/numbat/articles/numbat.html).

**A loss label is not an absolute zero-copy call.** The app's categorical CN
encoding uses zero for a loss class. It does not establish homozygous deletion.
Observed zeros, a statistical extra-zero component and biological loss of both
copies are separate concepts. See [E](e-per-cell-evidence/findings.md) and the
repo's `server/data.py` encoding.

**Missing evidence should be an output.** Report depth-change support,
allelic observability/uncertainty, baseline sensitivity and remaining
explanations separately. A model should be able to return a depth-supported
candidate with insufficient allelic evidence. That proposed behavior is not
implemented by this research pass.

## Ranked shortlist

Ranks are research priorities, not measured performance rankings.

| Priority | Candidate | Inputs → outputs | Fit and limitations | Code / license |
|---|---|---|---|---|
| 1 | Existing pyEpi baseline plus separate regional count diagnostics | Existing depth, raw counts, exposure, fixed regions → original calls plus exploratory depletion/zero summaries | Immediately useful; preserve original calls. A regional NB/Poisson diagnostic is a proposed addition, not validated scATAC CN inference. | [pyEpi](https://github.com/colomemaria/pyEpiAneufinder), GPL-3.0-only; new diagnostic not implemented |
| 2 | Numbat-multiome ATAC-bin | Raw binned counts, SNP-level phased counts, bin annotation, normalized ATAC reference → joint CN/clone inference | Direct scATAC comparator; current reference and SNP/count prerequisites unresolved. | [Numbat](https://github.com/kharchenkolab/numbat), MIT; metadata 1.5.2, not a pinned local run |
| 3 | Alleloscope regional allelic model | SNP/cell counts, region definitions, relative depth/reference information → regional haplotype proportions and state estimates | Relevant pooling approach; phase-level inputs and independently supported regions needed. | [Alleloscope](https://github.com/seasoncloud/Alleloscope); current official license/version unverified |
| 4 | Gap-aware count segmentation with conditional allele assessment | Raw counts, exposure, coordinates and callable masks; optional phase-valid alleles → relative boundaries plus evidence states | Research design transferable from count segmentation; no absolute CN anchor or ATAC calibration established. | [B's design](b-joint-segmentation/findings.md); no implementation selected |
| 5 | Technical and accessibility priors | Confirmed build, mappability/gap/blacklist tracks, selected SNPs, suitable atlas fragments → annotations and sensitivity analyses | Prefer soft flags initially. These priors cannot determine whether an individual dip is a deletion. | [F](f-technical-priors/findings.md), [C](c-reference-strategies/findings.md); code and track terms differ |
| 6 | Optional reference annotation / Atacformer representation | Authentic peaks/fragments, compatible vocabulary and labels → embeddings or annotation hypotheses | Useful only if independently evaluated; not a calibrated CN output or a way to create missing allele evidence. | [geniml](https://github.com/databio/geniml), BSD-2-Clause; checkpoint terms unresolved |

Numbat-multiome has published scATAC evaluation against WGS in selected tumour
cohorts. Its examples use normal reference information. Its event-level
arm/focal annotation excludes segments below 10% of an arm or below 200 kb,
so that analysis alone does not establish performance on every small focal
event of interest here. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12531993/).
Its documented input requirements differ from simply supplying HA/HB bin
matrices. [Official vignette](https://kharchenkolab.github.io/numbat/articles/numbat-multiome.html).

Alleloscope's scATAC workflow recommends DNA-supported segmentation. Its
regional pooling is useful to study, but using our candidate boundaries would
be a distinct experiment with its own validation needs.
[Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9127226/).

The Zhang adult human sci-ATAC atlas is a concrete public candidate for
cell-type accessibility comparisons across 30 tissues. It is not a verified
copy-neutral panel. Rebinning, donor/chemistry differences and independent CN
characterization need consideration before use as a reference.
[Atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC8664161/),
[example fragment-data record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5589349).

## What the 300-cell pilot found

The exact legacy panel has 100 cells per raw-coverage tercile and represents all
five WNN clusters. No new random draw or cohort filtering was performed. The
coordinator checked barcode and coordinate alignment and exact AT=HA+HB for
all 300 cells. The [manifest](shared/manifest.json) records hashes and environment;
[preparation](shared/PREPARATION.md) and [pilot](e-per-cell-evidence/pilot.md)
contain exact source and commands.

| Candidate window (Mb) | Retained bins | Cells with raw count zero across window | Cells with AT zero across window | Median per-cell window AT | Mean raw-density ratio, target / left flank; right flank |
|---|---:|---:|---:|---:|---|
| chr7 8.3–10.7 | 24 | 7 / 300 | 217 / 300 | 0 | 0.142; 0.323 |
| chr7 51.4–53.5 | 19 | 22 / 300 | 255 / 300 | 0 | 0.143; 0.199 |
| chr7 67.4–70.0 | 26 | 10 / 300 | 248 / 300 | 0 | 0.146; 0.119 |
| chr7 70.0–73.3 | 33 | 0 / 300 | 144 / 300 | 1 | 0.616; 0.238 |
| chr9 20–27 | 52 | 0 / 300 | 116 / 300 | 1 | 0.632; 0.910 |
| chr9 22.5–26.6 | 25 | 0 / 300 | 204 / 300 | 0 | 0.467; 0.560 |

Values come from [zero-support tables](e-per-cell-evidence/artifacts/zero_support_by_stratum.tsv)
and [flank-density comparisons](e-per-cell-evidence/artifacts/flank_density_comparison.tsv).
Densities use retained span; each flank contains 50 retained bins and can cross
gaps or overlap another candidate. These are not independently known neutral
controls. Coordinates in the export are 1-based inclusive, as documented in
the manifest.

Raw-zero bins comprise approximately 51–82% of these windows, but whole-window
raw-zero cells are much less common. A simple zero-versus-nonzero split is
therefore not sufficient to identify two biological populations. Conversely,
nonzero counts somewhere in a wide window do not exclude a smaller deletion
inside it. No null distribution or excess-zero test was fitted.

Per-cell allelic support is extremely sparse: no selected cell reaches a total
AT of 10 in any listed target window. That is a descriptive regional count,
not the app's per-bin pseudobulk filtering threshold. This observation alone is
not a calibrated power analysis. Pooling, phase and measurement-unit decisions
must precede a numerical allele test. The 300 cells are a feasibility subset,
not independent validation or a rare-event benchmark.

## Next research sequence

The [experiment protocol](synthesis/NEXT_EXPERIMENTS.md) defines inputs, outputs,
comparisons and decision points in detail:

1. Reconcile count units, SNP ascertainment, phase continuity and retained gaps.
2. On fixed regions, compare raw-count baselines and determine where allele
   evidence could distinguish specified alternatives; quantify abstention.
3. Assess Numbat-multiome input compatibility and benchmark on independently
   characterized samples, including a no-normal setting and small events.
4. Only then assess whether technical priors, a joint model or an optional
   representation improves held-out performance over the existing baseline.

The primary recommendation is steps 1 and 2. They resolve the most immediate
uncertainties without changing the app or assuming the candidate events are
true deletions. Model fitting, biological validation and new production
features remain future work.

## Review scope and limitations

This bounded pass screened the requested categories, checked leading claims
against primary papers/documentation, and performed one descriptive pilot.
It did not install/run complete callers, train models, download atlases, fit
mixtures, establish normal cells or determine true CN at any candidate locus.
Some licenses, exact commits, phase/count provenance, reference suitability
and benchmark access remain unresolved; see the
[uncertainty register](synthesis/uncertainty_register.tsv).

All research artifacts are new local files. Source data, application behavior,
existing caches, archived reference/validation files and the user's server
were not modified. [Artifact checks](synthesis/artifact_validation.json) verify
structure and shared-data integrity; they are not scientific validation.
