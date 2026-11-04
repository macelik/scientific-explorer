# Shared inputs

Coordinator-owned. Shared inputs are now frozen. Use the exact common cells.tsv membership; do not independently resample. Existing source data remain read-only.

- `manifest.json`: source/output SHA-256 identities, repo state, environment, conventions and completed alignment checks.
- `cells.tsv`: 300 barcodes, cluster, coverage tercile and library-size covariates.
- `regions.tsv`: fixed candidate/comparison windows and 50-retained-bin flanks, including gaps and callable span.
- `region_bins.tsv`: explicit original-grid membership of every exported region.
- `region_counts.tsv`: per-cell raw/processed depth and HA/HB/AT summaries for 21 windows/flanks (6,300 rows).
- `existing_chr7_chr9_segments.tsv`: copied current genome-pool segment records for context; not truth labels.
- `PREPARATION.md`: exact deterministic preparation source and command.
- `WAVE1_BRIEF.md`: cross-task findings and corrections for wave 2.

Regional HA/HB sums are descriptive; global phase continuity is not established. AT counts selected SNP-site observations, not established unique fragments. Flanks can overlap other candidates and are not independent truth-labelled controls.
