"""Prepare a local, checksum-verified data handoff; never upload or alter inputs.

The payload contains original files, with one copy of the identical AnnData.
The external pyEpi source is pinned and hashed, not relicensed or redistributed.
"""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TAG = "cloud-review-data-2026-11-04"


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace-root", type=Path, default=REPO.parent)
    parser.add_argument("--allelic-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("/tmp/cnv-cloud-handoff-2026-11-04"))
    args = parser.parse_args()
    workspace = args.workspace_root.resolve()
    archive = args.output_dir / f"{TAG}.tar.gz"
    if archive.exists():
        raise FileExistsError(f"Preserving existing archive: {archive}; choose a new --output-dir")
    # Rebuild this reviewed handoff only with the same upstream source identity.
    reviewed = json.loads((HERE / "DATA_MANIFEST.json").read_text())
    assert reviewed["upstream"]["commit"] == "e8305ea2162634d07e92a3bb401df672e48f217b"
    for name, digest in reviewed["upstream"]["package_files_sha256"].items():
        assert sha(workspace / "pyEpiAneufinder" / name) == digest, f"Upstream source changed: {name}"
    payload = []
    def add(source, target, role):
        assert source.is_file(), source
        payload.append({"source": source.resolve(), "path": target, "role": role})

    for name in ["count_matrix.h5ad", "binned_genome.csv", "hg38-blacklist.v2.bed", "breakpoints.csv", "clusters.json"]:
        add(workspace / "pyepi_results" / name, "pyepi_results/" + name, "depth inputs or existing per-cell baseline")
    add(workspace / "integration-story/data/integrated.h5mu", "integration-story/data/integrated.h5mu", "cluster labels and integration metadata")
    for path in sorted((workspace / "integration-story/data/genome_pool").rglob("*")):
        if path.is_file():
            add(path, str(path.relative_to(workspace)), "current genome-pool segmentation and derived baseline")
    add(workspace / "integration-story/data/raw_accepted_breakpoints.tsv", "integration-story/data/raw_accepted_breakpoints.tsv", "older chromosome-pool reproduction reference")
    for path in sorted((workspace / "default_300cells").iterdir()):
        if path.is_file() and path.suffix != ".png":
            add(path, "default_300cells/" + path.name, "existing cohort membership and baseline")
    for name in ["HA.parquet", "HB.parquet", "AT.parquet"]:
        add(args.allelic_dir / name, "aligned_matrices/" + name, "selected SNP-site haplotype observations")

    assert len({r["path"] for r in payload}) == len(payload)
    print(f"Hashing {len(payload)} original input files...", flush=True)
    records = []
    for entry in payload:
        source = entry["source"]
        records.append({"path": entry["path"], "size_bytes": source.stat().st_size,
                        "sha256": sha(source), "role": entry["role"]})
    upstream = workspace / "pyEpiAneufinder"
    manifest = {
        "date": "2026-11-04", "repository": "macelik/scientific-explorer",
        "branch": "allelic-segmentation", "repo_commit_at_preparation": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        "sample": "MCF10-A development dataset; not a region-level CN truth set",
        "dimensions": {"cells": 3460, "retained_bins": 25337, "nominal_bin_size_bp": 100000},
        "coordinate_convention": "1-based inclusive",
        "upstream": {"repository": "https://github.com/colomemaria/pyEpiAneufinder.git",
                     "commit": "e8305ea2162634d07e92a3bb401df672e48f217b",
                     "local_checkout_head_at_audit": "c77d82c356642ed1b3ace978c74353763a67a28d",
                     "public_commit_verified": True,
                     "package_files_sha256": {p.name: sha(p) for p in sorted(upstream.glob("*.py"))},
                     "note": "Local HEAD is nine documentation/analysis commits ahead of this public source commit; all 12 package Python files are identical. Fetch the public commit and verify hashes. Source identity does not prove historical execution."},
        "missing_for_complete_upstream_reanalysis": ["pre-QC matrix retaining excluded cells and bins", "original fragments/BAM and definitive alignment/reference/filtering invocation", "historical preprocessing parameters and intermediate GC predictions", "SNP-level selection/phase-set and pileup provenance", "independent region-level DNA truth"],
        "files": records,
    }
    manifest_text = (json.dumps(manifest, indent=2) + "\n").encode()
    (HERE / "DATA_MANIFEST.json").write_bytes(manifest_text)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Creating local archive: {archive}", flush=True)
    with archive.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=1, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w|") as bundle:
                for entry in payload:
                    print("  " + entry["path"], flush=True)
                    bundle.add(entry["source"], arcname=f"{TAG}/{entry['path']}", recursive=False)
                info = tarfile.TarInfo(f"{TAG}/MANIFEST.json")
                info.size = len(manifest_text)
                bundle.addfile(info, io.BytesIO(manifest_text))
    digest = sha(archive)
    assert archive.stat().st_size < 2 * 1024**3, "Archive exceeds GitHub release-asset limit"
    checksum = f"{digest}  {archive.name}\n"
    (args.output_dir / "SHA256SUMS").write_text(checksum)
    package = {"status": "prepared_locally_not_uploaded", "release_tag_proposal": TAG,
               "filename": archive.name, "size_bytes": archive.stat().st_size,
               "sha256": digest, "sha256sums_filename": "SHA256SUMS", "payload_files": len(records),
               "uncompressed_payload_bytes": sum(r["size_bytes"] for r in records),
               "local_output_directory": str(args.output_dir.resolve()),
               "repository_visibility_checked": "public", "public_distribution_confirmation": "This rebuild has not been published; confirm its contents before distribution.",
               "download_url": None}
    (HERE / "PACKAGE.json").write_text(json.dumps(package, indent=2) + "\n")
    print(json.dumps(package, indent=2), flush=True)


if __name__ == "__main__":
    main()
