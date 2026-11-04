# MCF10-A development data for the cloud method review

Companion data for Scientific Explorer's `allelic-segmentation` branch and
`research/2026-11-04-cloud-review-handoff/PROMPT.md`.

The archive contains 24 original files: the full AnnData with processed depth
and a post-QC raw-count layer, the integration MuData, HA/HB/AT matrices,
existing cohort files, and baseline depth-segmentation/CN outputs. It includes
an input manifest with SHA-256 hashes. Verify the archive with `SHA256SUMS`
before extracting; run the repository's `verify_handoff.py` afterward.

Dimensions: 3,460 cells × 25,337 retained 100 kb bins. MCF10-A is development
data and is not region-level CN truth. Raw counts in AnnData already underwent
cell/bin QC; excluded observations, original reads and independent DNA truth
are not supplied. No new preprocessing or inference was performed when
packaging these files.

The external pyEpi source must be fetched separately at the public commit
specified in `DATA_MANIFEST.json`, then checked against its per-file hashes.
Follow the handoff README for cloud paths, dependencies, limitations and the
research assignment. A normal Git clone does not include these release assets.
