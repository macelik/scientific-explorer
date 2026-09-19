# Original integration diagnostics: verified reproduction

Reproduced 2026-09-19 from unmodified upstream scripts and the local story data.
The previous assertion that the original generating code/configuration was unavailable
was incorrect: it was available in `/home/mami/Desktop/EGA/SNP_project/beta_revised`.

All **11 generated TSVs are byte-for-byte identical** to the archived originals:
nine merge-prototype tables, plus the SRD screen's scores and bootstrap summaries.
The depth matrix, MuData and accepted-breakpoint files also match their upstream
counterparts by full SHA-256. Prototype module hashes match the original manifest.

- `analysis/depth_segment_prototype/`: original runner, loader and merge rules.
- `src/`: original screening script and its SRD/bootstrap dependencies.
- `archived/`: original reference tables, manifest and candidate audit.
- `regenerated/`: fresh prototype tables, manifest and 31 figures; verified screen
  tables and original candidate audit are also copied here for the app.
- `screening-regenerated/`: fresh screen tables and ten figures.
- `regenerated/verification.json`: input/code/table hashes, numerical checks and settings.
- `regeneration.log`, `screening-regeneration.log`: generator output.

The prototype's 31 PNGs also match the archived PNGs byte-for-byte. The screen's
ten figures were regenerated; their PNG bytes were not used as a parity criterion.
The observed result is 0 merges, 299 retained boundaries and 334 final segments.
Zero merges is now a verified result, not a fallback for missing tables.

## Repeat

From the workspace root:

```bash
python3 scientific-explorer/reference/integration-prototype/reproduce.py
# Check existing regenerated outputs only:
python3 scientific-explorer/reference/integration-prototype/reproduce.py --verify-only
```

The script uses the local source snapshots and story inputs, verifies exact input
hashes, runs both generators and compares the outputs. It requires the original
scientific dependencies, including anndata, muon and matplotlib (these are not
required merely to launch the app with the saved tables). It removes the success
marker before a new check and writes it only after all table/manifest comparisons
pass. Restart the scientific explorer after regenerating.

The app automatically uses `regenerated/` when the story prototype directory is
absent. An explicit `PYEPI_SCIENTIFIC_PROTOTYPE` takes precedence. The user-deleted
`integration-story/data/prototype` was not recreated, and `interactive-explorer`
was not modified.

## Keep the three bootstrap calculations distinct

| Calculation | Segment estimator | Replicates / seed | Meaning |
|---|---|---|---|
| Original merge prototype | Median of per-bin pseudobulk | 200 / 20260907 | Fraction below the original effect threshold; all valid replicates and chain guard required |
| Original SRD screen | Mean-based rate contrast and SRD | 200 / 42 | Conditional SE, percentile interval and sign agreement; all-source resampling is stratified by cluster |
| Alternative interactive view | Arithmetic mean | 200 / 42 | Paired cell-resampling comparison for the selected three regions; separately implemented draw stream |

Matching seeds across different implementations do not imply matching draws.
All calculations keep membership and accepted boundaries fixed. SRD/√φ and
bootstrap fractions are descriptive diagnostics, not calibrated significance.
The candidate audit is an archived original, not a newly rerun AD search. WNN,
embeddings, accepted breakpoints and CN calls were not recomputed in this step.

The original prototype uses its original 100-kb coordinate convention and gap
plotting policy. Its plots were kept unchanged for reproduction; the interactive
alternative uses actual bin ends. See the [integration review](../../validation/integration-review/REVIEW.md)
for interpretation and alignment caveats.
