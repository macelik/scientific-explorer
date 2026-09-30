# Scientific workflow

This document explains how this repository is used for method development,
not merely how to operate the user interface. The interactive application is
the visible layer of a longer diagnostic loop — most of the work it supports
happens in the reasoning between one run and the next, not inside a single
screenshot.

## The loop

```text
Observed behaviour
      ↓
Hypothesis
      ↓
Diagnostic / counterfactual experiment
      ↓
Recompute AD / segmentation / CN
      ↓
Inspect ranking + significance + downstream consequences
      ↓
Accept / reject / refine hypothesis
      ↓
Modify method or move to another analysis level
```

No single visualization in this repository claims to be the conclusion of
that loop. Each is a checkpoint — a way to inspect one iteration's
consequences before deciding whether to run another.

## Worked example 1: does the segmentation actually recover a simulated event?

- **Observed behaviour:** the production recursive segmentation accepts or
  rejects candidate breakpoints under a fixed recursion depth and
  significance threshold. It isn't obvious, just from looking at accepted
  calls on real cells, how sensitive that acceptance is to an event's width.
- **Hypothesis:** a simulated copy-number event of a given width and local
  dispersion should be recoverable as an accepted breakpoint, but detection
  should degrade as the event narrows.
- **Diagnostic experiment:** the **X experiments** width-series feature —
  10, 25, 50, 100, and 150-bin scenarios built by sampling from the same
  cell's own genome-wide bins under a matching production call (loss/base/
  gain), never from a separate simulator.
- **Recompute:** original-vs-edited AD curves, both the local and global
  permutation tests (1,000 permutations, seed 42, strict `p < 0.001`), and
  four CN comparison rows (original calls, edited values with production
  boundaries, automatic resegmentation, and diagnostic segmentation with
  event edges imposed).
- **Inspect:** event-edge diagnostics — was each intended boundary accepted,
  what was its rank among candidates in the tested node, and both p-values;
  Holmes/Watson/combined call counts under each boundary treatment.
- **Revise:** the honest finding, already stated in this repository's
  scientific-scope notes, is that "Width-series results describe these
  specific realizations, not calibrated detection probabilities or a
  universal size threshold." A width series is a set of concrete
  observations to reason from, not a calibration curve.

## Worked example 2: is a small segment a real boundary or an artifact of weak flank support?

- **Observed behaviour:** a cluster's pseudobulk accepted-breakpoint
  segmentation sometimes produces adjacent small segments whose flank
  support (SRD, SRD/√φ, delta log2FC) looks weak or internally inconsistent.
- **Hypothesis:** specific adjacent segments should not have been split, or
  the metric used to judge that boundary is itself ambiguous at that
  location.
- **Diagnostic:** the three hatch layers in Chromosome tracks — transition
  zone (opposite-sign flank scores), merge proposal (a flank score below
  threshold), and ambiguous flank preference (the weaker SRD-variant flank
  disagrees with the weaker delta-log2FC flank) — followed by the
  exploratory adjacent-segment merge tool, whose eligibility is the union of
  whichever layers are enabled.
- **Recompute:** one boundary merges per step; pooled dispersion φ is
  recomputed, unconditionally, after every merge, since it is shared across
  the whole chromosome's segmentation.
- **Inspect:** the merge-history table, the step slider, and the
  original-vs-current segmentation plots side by side.
- **Revise:** the tool is explicitly labelled, in its own UI copy, as "an
  independent exploratory rule, not the original bootstrap/chain-guard merge
  prototype," and it "does not recompute CN calls on the merged boundaries" —
  a deliberate scope limit, not an oversight.

## What this repository's diagnostic capabilities support, by loop stage

- **Forming a hypothesis:** Cohort karyogram, Explore output, Integration ·
  clusters (UMAP/cluster inspection).
- **Designing a counterfactual experiment:** X experiments (region edits,
  simulate/extend/duplicate transformations, width series), manual
  breakpoint placement.
- **Recomputing under the perturbation:** automatic resegmentation at
  configurable recursion depth (k = 1–10), full-genome CN fitting, the
  hatch-layer scores and adjacent-segment merge tool, alternative segment
  estimators (arithmetic / upper-tail IQR / two-sided IQR).
- **Inspecting consequences:** event-edge diagnostics, both permutation
  tests, four CN comparison rows, the alternative flank view and paired-cell
  bootstrap, the merge-history table.
- **Revising the hypothesis or method:** provenance-carrying JSON/CSV
  exports (seed, parameters, dataset/implementation identity) so a result
  can be re-examined, disputed, or reproduced later.

## What this demonstrates

The structure of this repository — tests alongside features, validation
reports alongside claims, exploratory work labelled as exploratory — is
evidence of a working style: debugging scientific algorithms by direct
inspection rather than by trusting summary statistics alone, designing
targeted counterfactual experiments instead of guessing from outputs,
tracking provenance so a result can be checked later, and revising a method
only after inspecting what actually happened. It is a description of how the
work was done, not a claim about what the work proved.

## Exploratory vs. validated

See [VALIDATION.md](VALIDATION.md) for what has actually been checked —
byte-for-byte reproduction of archived tables, and browser/parity test
coverage — versus what remains exploratory or descriptive. Some concrete
examples of language this repository uses deliberately: SRD/√φ is "a
dispersion-scaled deviance statistic, not a p-value"; the hatch/merge layers
are "an independent exploratory rule"; the simulation transformation "is not
a validated biological CN simulator." None of these are hedges — they are
precise statements of scope.
