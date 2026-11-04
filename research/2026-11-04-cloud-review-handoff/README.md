# Cloud review handoff

Give the cloud agent **[PROMPT.md](PROMPT.md)**. It requests a review of the
whole scientific workflow, including QC, GC correction, rounding, zero
clipping, normalization, pseudobulk aggregation, segmentation and CN evidence.
The requested outcome is an evidence-supported recommendation with focused
diagnostics and saved results, before adopting production scientific changes.

Use branch **`allelic-segmentation`**, not an assumed default-branch checkout.
The earlier research baseline is `994fb3a`. Record the actual handoff commit
when starting the cloud task.

## What is available

| Item | In normal Git checkout? | Available locally / in prepared package? |
|---|---|---|
| Workbench code, research, guide, 300-cell regional summaries | Yes | Yes |
| Full AnnData: processed X and post-QC raw layer | No | Yes; 259.49 MiB |
| MuData: integration, labels and graphs | No | Yes; 296.41 MiB |
| HA/HB/AT matrices | No | Yes; about 43.77 MiB combined |
| Existing cohort and segmentation/CN baseline outputs | Only selected summaries | Yes |
| Exact external pyEpi source | No | Fetch the public pinned commit below |
| Pre-QC counts, original fragments/BAM, full historical run provenance | No | Not included |
| Independent region-level DNA truth | No | Not established |

The AnnData has **3,460 cells × 25,337 retained 100 kb bins**. Its `raw` layer
supports revisiting GC correction and normalization on those retained
observations. It cannot recover excluded cells/bins to evaluate relaxed QC.
Read [UPSTREAM_AUDIT_NOTES.md](UPSTREAM_AUDIT_NOTES.md) for the verified source
path and important remaining uncertainties.

## Data delivery

`PACKAGE.json` is the authoritative delivery status and archive identity.
The archive was prepared locally; do not assume that a release exists while
its `download_url` is null. All payload files and hashes are listed in
`DATA_MANIFEST.json`; `PACKAGE_VERIFICATION.json` records archive checks.

The repository is **public**. The user authorized public distribution of this
specific 24-file data package on 2026-11-04. `PACKAGE.json` records whether
publication has completed; a locally prepared archive alone is not a download.

Ordinary GitHub commits block files larger than 100 MiB. GitHub Releases allow
assets under 2 GiB each, making a separate checksum-verified release asset a
suitable delivery route for this prepared package. It keeps the data out of
Git history. See [GitHub file limits](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)
and [release asset limits](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases).
Cloning the repository does not download release assets automatically.

The prepared archive contains 24 original files, approximately 1.26 GiB before
compression. It includes one copy of the AnnData, because the local
`pyepi_results/` and `integration-story/data/` copies have identical hashes.
The package does not contain the pyEpi source; retrieve the exact public
revision separately. It does not contain original sequencing reads.

## Cloud setup

First clone/check out `scientific-explorer` on `allelic-segmentation`. If the
data package has been published, obtain its actual release tag and URL from
`PACKAGE.json`. The following is the proposed download command, to be used only
after that status records publication (or transfer the local archive privately):

```bash
gh release download cloud-review-data-2026-11-04 \
  --repo macelik/scientific-explorer \
  --pattern 'cloud-review-data-2026-11-04.tar.gz' --pattern 'SHA256SUMS' \
  --dir ../cloud-review-download
```

Verify `SHA256SUMS` from within the download directory, then extract the archive
into the directory alongside the workbench checkout. This creates a sibling
`cloud-review-data-2026-11-04/` containing the original data layout. Allow at
least 2 GiB free for the download and extraction, plus room for environments,
caches and experiments. Use backed/sparse loading for initial diagnostics.

Retrieve the exact publicly available upstream source:

```bash
git clone https://github.com/colomemaria/pyEpiAneufinder.git ../pyepi-upstream
git -C ../pyepi-upstream checkout e8305ea2162634d07e92a3bb401df672e48f217b
```

From the workbench root, validate the extracted data and the upstream files:

```bash
python3 research/2026-11-04-cloud-review-handoff/verify_handoff.py \
  --data-root ../cloud-review-data-2026-11-04 \
  --package-dir ../pyepi-upstream/pyEpiAneufinder
```

Configure paths explicitly; the original absolute paths in older manifests
refer to the author's machine and must not be mistaken for cloud paths:

```bash
export PYEPI_SCIENTIFIC_ROOT="$PWD/../cloud-review-data-2026-11-04"
export PYEPI_SCIENTIFIC_PACKAGE="$PWD/../pyepi-upstream/pyEpiAneufinder"
export PYEPI_SCIENTIFIC_H5AD="$PYEPI_SCIENTIFIC_ROOT/pyepi_results/count_matrix.h5ad"
export PYEPI_SCIENTIFIC_H5MU="$PYEPI_SCIENTIFIC_ROOT/integration-story/data/integrated.h5mu"
export PYEPI_SCIENTIFIC_ALLELIC="$PYEPI_SCIENTIFIC_ROOT/aligned_matrices"
export PYEPI_SCIENTIFIC_PROTOTYPE="$PWD/reference/integration-prototype/regenerated"
```

Use an isolated Python environment and read both the workbench requirements
and the pinned upstream `pyproject.toml`. Research tasks additionally need
compatible `anndata`, `pyarrow` and `scikit-misc` versions; record the resolved
environment. `ENVIRONMENT.json` records the currently installed local versions,
not a verified historical run environment. The workbench requirements alone
do not describe every optional research dependency. Run bounded checks before
full scans; the offline
segmentation tool otherwise chooses workers from available CPU count.

## Rebuilding the local package

`build_data_package.py` reads the original external files and writes a new
archive, manifest and checksums. It does not upload anything. Supply a fresh
output directory to preserve an existing archive:

```bash
python3 research/2026-11-04-cloud-review-handoff/build_data_package.py \
  --workspace-root /path/to/original/workspace \
  --allelic-dir /path/to/original/aligned_matrices \
  --output-dir /tmp/new-cloud-handoff
```

It preserves original data bytes; packaging performs no new sampling,
correction, normalization, inference or biological validation.
