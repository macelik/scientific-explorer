# Saved state and resume instructions

**Resolved:** the user reset usage limits and all three interrupted Terra
reviews subsequently completed. The [final first-pass report](../report.md)
and [current status](run_status.json) supersede the pending-work instructions
below. This document retains the earlier interruption record; there are no
outstanding agent tasks from that interruption.

Saved 2026-10-29. The requested structured workspace is in place and all six
agents were dispatched explicitly with `gpt-6.1-terra` and fresh context.
Three first-pass tasks completed. All three second-wave turns returned a usage
limit error. Their files are preserved; none should be treated as accepted
complete solely because an agent wrote `status: complete` before the failure.

| Task | Saved state | Resume work |
|---|---|---|
| A — callers | Completed bounded review, 11 method records | Read results; do not repeat search. Important remaining input/license/version uncertainties are explicit. |
| E — per-cell evidence | Completed review and descriptive 300-cell pilot | Read saved tables and exact script. No mixture fitting or formal power calculation was performed. |
| F — technical priors | Completed review and local provenance audit | Read count-unit and phasing caveats; historical input provenance remains unresolved. |
| B — joint segmentation | Substantial draft and structured files saved; turn interrupted | Review primary-source claims and completeness, resolve any missing scope, validate, then finalize. |
| C — reference strategies | Substantial draft and structured files saved; turn interrupted | Review decision tree and primary sources, validate, then finalize. |
| D — representations | Method/source records and pilot note saved; findings contains only a progress placeholder | Verify records and complete findings, ranking and limitations before finalizing. |

## Resume without repeating work

1. Wait until Terra capacity is available; do not silently substitute another
   model for the user's selected research model.
2. Read each task's status.json, TASK.md, existing outputs, the workspace README,
   fields.yaml and shared/WAVE1_BRIEF.md. Reuse the frozen shared panel and hashes.
3. Resume B, C and D with explicit `gpt-6.1-terra`, at most three concurrently.
   Ask each to finish and verify existing artifacts, not start from scratch.
4. Audit citations and leading claims, check all JSON fields/source references,
   and remove contradictory progress/completion claims before marking complete.
5. Finish the cross-category shortlist and concrete next-experiment sequence.
   It is not yet a final six-category synthesis.

## Coordinator review already performed

- Confirmed the Numbat-multiome ATAC-bin workflow in the 2025 primary paper and
  current official vignette; corrected the initial RNA-only framing in A/E and
  the shared brief. Its published examples use normal reference information.
  [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12531993/),
  [vignette](https://kharchenkolab.github.io/numbat/articles/numbat-multiome.html).
- The paper's event-level arm/focal annotation excludes segments under 10% of
  an arm or under 200 kb; that particular evaluation does not by itself establish
  performance for all focal windows in this project. Preserve this limitation
  when choosing benchmarks. [Methods: arm/focal annotations](https://pmc.ncbi.nlm.nih.gov/articles/PMC12531993/).
- Shared exports checked all 300 barcodes, matching coordinates, integer
  nonnegative raw counts and exact AT=HA+HB. Full source identities, environment
  and deterministic preparation source are saved under shared/.
- E's pilot is descriptive: numerous empty bins do not imply empty multi-bin
  windows, and neither observation establishes a structural-zero component.
  No p-values, deletion probabilities or calibrated allelic power were produced.
- The proposed model must generalize beyond the non-malignant MCF10-A sample;
  no immune/stromal reference compartment or region-level diploid truth is assumed.

Application code and scientific pipeline behavior were not changed. Existing
input data, caches, reference files, validation archives and port 8766 were not
modified. Research artifacts are new, local files.
