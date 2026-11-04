# Low-coverage CNV research

Research workspace created 2026-10-29. Scope and task briefs: [research plan](../../docs/notes/2026-10-29-low-coverage-research-plan.md).

## What question are we investigating?

**When a short region has unusually low or high scATAC depth, how can we
determine whether it reflects a copy-number change, accessibility or technical
effects, sampling variation, or insufficient information to make a call?**

This folder started from a concrete problem: the depth segmentation produced
short regions called as losses, especially on chr7, but the allelic track had
very little usable information in those same regions. We wanted to understand
which evidence could distinguish the possible explanations, and how other
methods handle this situation.

### Where this sits in our workflow

The relevant depth workflow is:

```text
QC-filtered cells and genomic bins
    → GC correction and normalization of per-cell depth
    → cluster pseudobulk depth profiles
    → AD segmentation: find supported changes along the chromosome
    → relative CN assignment on the resulting segments
    → assess depth, allelic support and alternative explanations
```

QC is quality control; GC correction addresses dependence on the proportion of
G/C bases. The displayed depth pseudobulk is the per-bin **mean of corrected,
normalized X values** across cells in each cluster, plus an all-cells profile.
It is therefore a processed depth profile; the raw fragment counts remain
separately available in `layers/raw`. See the
[pseudobulk implementation](../../server/integration.py) and
[scientific workflow](../../docs/SCIENTIFIC_WORKFLOW.md).

A stretch of bins with a different depth distribution can produce boundaries
on either side and become a short segment, provided it passes the scan's size
constraints and local/global permutation tests. A dip can lead to a relative
loss call; an increase can lead to a relative gain call. A detected boundary
establishes a change under the segmentation procedure, but does not identify
its biological cause. A short segment is not automatically an erroneous split.

### Why the interpretation is difficult

Several explanations can produce a similar-looking depth dip:

- A genuine copy-number loss, potentially present in only some cells.
- Low chromatin accessibility, which changes what ATAC measures without
  necessarily changing DNA dosage.
- Technical effects, including mapping difficulty, retained-grid gaps or an
  unsuitable normalization/reference baseline.
- Sparse-count fluctuations that produce an unstable apparent contrast.

Pooling cells reduces some sampling variation, but a shared regional effect
can persist in every cluster. The same general ambiguity applies to increases:
more coverage can reflect a gain, greater accessibility or technical effects.
The saved pilot concentrates on low-coverage candidate windows; it is not an
equivalent investigation of gains.

“No signal” needs care here. **A measured dip, an observed count of zero, and
too little evidence to interpret a region are different situations.** In
particular, the depth track may contain a clear dip while the allele track is
unavailable. Allelic counts require observations at selected informative SNPs;
their availability depends on more than total depth. After low-AT masking,
there may be too little information to assess imbalance. That does not show
that alleles are balanced or that a deletion is absent.

The later research synthesis revises an overstrong statement in the original
handoff: lack of an obvious LOH signature does **not**, by itself, rule out
one-copy loss. Similarly, the app's categorical CN value `0` means its loss
class, not a demonstrated absolute zero-copy deletion. These distinctions are
explained in [the synthesis](report.md).

### What this directory contributes

Yes: it investigates **how other methods deal with these ambiguities and what
we could adapt to our workflow**. The six reviews cover existing CNV callers,
joint depth/allele segmentation, reference selection, representation models,
per-cell count/allele evidence, and technical annotations. They consider both
whether a boundary is supported and what can be inferred about the resulting
region, including when the appropriate result is “insufficient evidence.”

The folder also contains a descriptive check of six candidate windows in an
existing 300-cell panel. It found sparse allelic support and showed why many
zero bins are not equivalent to an entirely empty region. It did not determine
which candidate regions are true deletions or run a benchmark of all reviewed
callers. The [next-experiment plan](synthesis/NEXT_EXPERIMENTS.md) describes how
to move from these observations to calibrated comparisons and independent
validation. MCF10-A is the development dataset; the intended method must
generalize beyond this cell line.

## Current state and reading order

Current state: **All six bounded first-pass reviews and the shared descriptive
pilot are complete.** Start with the [research synthesis](report.md) and
[proposed next experiments](synthesis/NEXT_EXPERIMENTS.md). The
[machine-readable status](synthesis/run_status.json) records completion;
[RESUME.md](synthesis/RESUME.md) retains the resolved interruption history.

For a guided walkthrough, open the [interactive method guide](interactive-guide/index.html)
locally. It explains each research entry with steps, examples, workflow
adaptations and citations, and includes interactive teaching models and the
saved 300-cell summaries. See its [usage and provenance notes](interactive-guide/README.md).

For an independent review that also audits QC, GC correction, zero clipping and
normalization, use the [cloud-agent prompt and data handoff](../2026-11-04-cloud-review-handoff/README.md).

The development sample is **MCF10-A, a non-malignant cell line**. The objective is a general scATAC CNV tool; sample identity is not region-level CN truth. No normal population should be assumed to exist.

## Directory ownership

| Directory | Owner | Purpose |
|---|---|---|
| shared/ | Coordinator | Frozen input manifest, shared 300-cell membership and summaries |
| [a-cnv-methods/](a-cnv-methods/findings.md) | Agent A, GPT 6.1 Terra | Existing CNV methods; completed first pass |
| [b-joint-segmentation/](b-joint-segmentation/findings.md) | Agent B, GPT 6.1 Terra | Joint segmentation and uncertainty; completed first pass |
| [c-reference-strategies/](c-reference-strategies/findings.md) | Agent C, GPT 6.1 Terra | Reference strategies; completed first pass |
| [d-representation-models/](d-representation-models/findings.md) | Agent D, GPT 6.1 Terra | Representation and foundation models; completed first pass |
| [e-per-cell-evidence/](e-per-cell-evidence/findings.md) | Agent E, GPT 6.1 Terra | Per-cell loss evidence; completed review and descriptive pilot |
| [f-technical-priors/](f-technical-priors/findings.md) | Agent F, GPT 6.1 Terra | Technical priors and SNP ascertainment; completed first pass |
| synthesis/ | Coordinator | Cross-task review, ranked shortlist and next experiments |

Agents write only in their assigned directory. Shared files are read-only for agents. All six task folders use the same [output schema](fields.yaml). See [execution outline](outline.yaml).

## Required saved outputs in every agent directory

- `findings.md`: decision summary, findings with primary citations, ranked recommendations, rejected options, limitations, next experiment.
- `methods.json`: structured method/prior records using the schema below; unknown facts remain explicit.
- `sources.json`: primary-source register with stable IDs, URLs, access dates and claim locators.
- `pilot.md`: either the pilot question, inputs, exact commands/source, parameters/seed/environment, outputs and limits, or why computation was skipped.
- `status.json`: task status, model, artifact list and unresolved questions; mark complete only after saving and checking all outputs.
- `artifacts/`: optional small diagnostic tables/figures; generated only when needed. No raw input copies.

Save an initial outline early and update it as work progresses. Final chat summaries must link saved files; chat-only findings do not count as completion. Do not manufacture results for an unrun pilot.

## Evidence and provenance

Use claim-level primary references. Distinguish published evidence in scATAC, evidence in another assay, transferred ideas and exploratory local diagnostics. Missing evidence is not balance; a pyEpi loss encoding is not absolute zero copies. Preserve exact barcode membership and bin coordinates, and keep processed depth separate from raw counts. Record uncertainties and conflicting sources.

Application code, input data, existing caches, archived reference/validation files and the user's server on port 8766 are outside the agents' write scope. Pilot scripts execute in isolated /tmp directories; preserve exact source in pilot.md so they can be reproduced. No large downloads, model training, full pipeline reruns or recursive agent spawning in this pass.
