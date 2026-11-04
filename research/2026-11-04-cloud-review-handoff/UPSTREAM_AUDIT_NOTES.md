# Verified starting points for the upstream review

Inspected 2026-11-04. These are source/data observations, not findings that a
particular preprocessing choice caused the candidate segments.

## Source identity

The external package used by the local workbench is the sibling
`../pyEpiAneufinder/`. Its 12 Python files are byte-identical to the package in
the public `colomemaria/pyEpiAneufinder` commit
`e8305ea2162634d07e92a3bb401df672e48f217b`, whose GitHub availability was checked.
The local checkout HEAD, `c77d82c356642ed1b3ace978c74353763a67a28d`, is nine
documentation/analysis commits ahead and was not found on GitHub. There is no
package or `pyproject.toml` difference between these two commits, and no tracked
package modifications were reported locally. Use the public commit for cloud
setup. Per-file hashes are in `DATA_MANIFEST.json`. Fetch the source separately;
the workbench repository's MIT license does not relicense the upstream package.
Historical execution still requires verification. The pinned pyproject declares
version `0.3.5`; this differs from the newer version discussed in the literature
review, so record and compare versions explicitly.

## Preprocessing path visible in that source

In `pyEpiAneufinder/pyEpiAneufinder.py`, approximately lines 201–280:

1. Keep bins meeting a nonzero-frequency criterion across cells.
2. Keep cells meeting a nonzero-bin-frequency criterion on the retained grid.
3. Apply IQR-based filters to cell total counts and raw-count standard
   deviations, intersecting the masks.
4. In the GC-enabled branch, fit LOESS for each cell and compute
   `correction = mean(raw_counts) / fitted_GC_values`.
5. Multiply raw counts by the correction; round with `np.rint` and cast to
   integer.
6. Save post-QC, pre-GC counts as `layers['raw']`.
7. Set negative corrected values to zero.
8. Scale each corrected row to a total of `1e5`; store it in `X`.

`render_fragments.py:get_loess_smoothed` uses a span of `0.75` and, on an
exception, prints the error and returns the input count vector. The reviewer
should investigate the effect of this fallback and denominator behavior;
we have not established that a fit failed in the historical run.

There are older/alternative functions, including `pyEpiAneufinder_v1.py` and
multiple somy-assignment routines. Follow the actual entry points. Additional
scaling during CN assignment is a separate stage from preprocessing `X`.

## Available data and limits

- The two local copies of `count_matrix.h5ad` in `pyepi_results/` and
  `integration-story/data/` are byte-identical (SHA-256
  `b98a4f86cb80b86dc49382ffa90ab44cc6b3f1028f8481c80535cbf39683d5f8`).
- Shape: 3,460 cells × 25,337 retained bins. `obs` contains cell identifiers;
  `var` contains chromosome/start/end/GC; the only layer is `raw`; `uns` is
  empty. The file does not contain a complete preprocessing execution log.
- The current AnnData is already filtered. Relaxing the original QC and
  recovering excluded cells/bins requires a pre-QC matrix or earlier inputs.
  Those are not included in this handoff.
- Cluster labels and integration metadata are in the separate `integrated.h5mu`.
  The explorer's pseudobulk uses the AnnData `X`, not the transformed depth
  modality inside that MuData object.
- Full HA/HB/AT matrices are external parquet files. Their barcode order must
  be joined by name; selected SNP-site DP is not established unique-fragment
  coverage. SNP-level inputs, phase-set continuity and the historical pileup
  invocation remain unresolved by this package.
- The original bin grid and blacklist file are useful audit inputs, but their
  presence does not establish the exact historical settings or recover counts
  for excluded bins.

This is why a cloud review should begin with provenance and reconstruction,
then compare small, justified alternatives. The current research does not yet
show that clipping, normalization or segmentation is the primary cause.
