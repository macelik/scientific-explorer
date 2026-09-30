# Validation

Every feature in this repository ships with a companion verification pass:
automated tests, live browser walkthroughs where interactive behavior can't
be unit-tested, and byte-for-byte reproduction checks against archived
reference outputs where a prior result exists to check against.

| Feature area | Authoritative report |
|---|---|
| X experiments (transformations, width series, event-edge diagnostics) | [`validation/VALIDATION_REPORT.md`](../validation/VALIDATION_REPORT.md) |
| Integration review (alternative flank view, paired-cell bootstrap) | [`validation/integration-review/REVIEW.md`](../validation/integration-review/REVIEW.md) |
| Segmentation-tab chart-loading responsiveness | [`validation/integration-freeze/REPORT.md`](../validation/integration-freeze/REPORT.md) |
| Exploratory hatching and adjacent-segment merge | [`validation/hatch-merge/HANDOFF.md`](../validation/hatch-merge/HANDOFF.md) |
| Original prototype byte-for-byte reproduction | [`reference/integration-prototype/README.md`](../reference/integration-prototype/README.md) |
| Reload / recursion-depth session restore | [`validation/reload-recursion/REPORT.md`](../validation/reload-recursion/REPORT.md) |

**Validated** means numerically reproduced, covered by an automated test, or
directly observed in a recorded browser walkthrough. **Exploratory** means a
new descriptive rule or alternative estimator introduced to support
investigation, not independently calibrated against ground truth. See
[SCIENTIFIC_WORKFLOW.md](SCIENTIFIC_WORKFLOW.md)'s closing section for the
general framing, and [`validation/README.md`](../validation/README.md) for
what each file under `validation/` actually is.
