"""Regenerate original diagnostics with frozen source; verify before app exposure.

Run from any directory. Requires anndata, muon, matplotlib and the scientific
Python dependencies used by the original scripts. --verify-only checks a run
already made by the two original generators.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
INPUTS = {
    "adata": ("count_matrix.h5ad", "b98a4f86cb80b86dc49382ffa90ab44cc6b3f1028f8481c80535cbf39683d5f8"),
    "integrated": ("integrated.h5mu", "b6c9c5bce06b47a3d85d0c26e6b42f2a70dfcb6c72bd009aaab695ae78658bec"),
    "breakpoints": ("raw_accepted_breakpoints.tsv", "d024eab3c62db6dc7d5489e2573ffa08a0fc8a482137b64163db27ab96f61525"),
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def compare_table(actual, expected):
    a, b = pd.read_csv(actual, sep="\t"), pd.read_csv(expected, sep="\t")
    # Export rounding is part of the reference. Only floating-point serialization
    # noise is tolerated; categories, row order, NaNs, booleans and integers agree.
    pd.testing.assert_frame_equal(a, b, check_exact=False, rtol=1e-12, atol=1e-12)
    return {"rows": len(a), "columns": len(a.columns),
            "byte_identical": sha256(actual) == sha256(expected),
            "sha256": sha256(actual), "reference_sha256": sha256(expected)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--verify-only", action="store_true")
    args = p.parse_args()
    archived, regenerated = HERE / "archived", HERE / "regenerated"
    screening = HERE / "screening-regenerated/tables"
    marker = regenerated / "verification.json"
    # Never leave an old success marker if a later regeneration/check fails.
    marker.unlink(missing_ok=True)
    inputs = {}
    for key, (name, expected) in INPUTS.items():
        path = WORKSPACE / "integration-story/data" / name
        actual = sha256(path)
        assert actual == expected, f"Reference input changed: {path}"
        inputs[key] = {"path": str(path), "sha256": actual}
    old = json.loads((archived / "run_manifest.json").read_text())
    for key, name in [("runner", "run_depth_segment_prototype.py"), ("loader", "loader.py"), ("merge_core", "merge_core.py")]:
        assert sha256(HERE / "analysis/depth_segment_prototype" / name).startswith(old["source_modules"][key])
    if not args.verify_only:
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                   NUMBA_CACHE_DIR="/tmp/pyepi-repro-numba", MPLCONFIGDIR="/tmp/pyepi-repro-mpl",
                   PYTHONPATH=str(HERE / "src/srd_effect"))
        common = ["--adata", inputs["adata"]["path"], "--integrated", inputs["integrated"]["path"],
                  "--breakpoints", inputs["breakpoints"]["path"], "--n-bootstrap", "200"]
        commands = [
            ("regeneration.log", [str(HERE / "analysis/depth_segment_prototype/run_depth_segment_prototype.py"),
             *common, "--reference-scores", str(archived / "small_segment_scores.tsv"),
             "--output-dir", str(regenerated), "--seed", "20260907"]),
            ("screening-regeneration.log", [str(HERE / "src/screening/screen_small_segments.py"),
             *common, "--cluster-key", "wnn_leiden_0.3", "--audit", str(archived / "segmentation_candidate_audit.tsv"),
             "--output-dir", str(screening), "--seed", "42"]),
        ]
        for log, cmd in commands:
            print(f"Running {cmd[0]}; log: {HERE / log}", flush=True)
            with (HERE / log).open("w") as f:
                subprocess.run([sys.executable, *cmd], cwd=WORKSPACE, env=env, stdout=f, stderr=subprocess.STDOUT, check=True)
    results = {}
    for original in sorted(archived.glob("*.tsv")):
        if original.name == "segmentation_candidate_audit.tsv":
            continue  # Original audit is preserved, not a new AD segmentation run.
        actual = (screening if original.name in ("small_segment_scores.tsv", "small_segment_bootstrap.tsv") else regenerated) / original.name
        results[original.name] = compare_table(actual, original)
    new = json.loads((regenerated / "run_manifest.json").read_text())
    fields = ["params", "source_modules", "cluster_sizes", "cluster_replicate_hashes",
              "reproduction_check", "totals_default", "n_processed_source_chrom", "sensitivity_grid_max_merges"]
    for field in fields:
        assert new[field] == old[field], f"Manifest mismatch: {field}"
    # Assemble the app's read-only table directory only after all checks succeed.
    for name in ["small_segment_scores.tsv", "small_segment_bootstrap.tsv"]:
        shutil.copy2(screening / name, regenerated / name)
    shutil.copy2(archived / "segmentation_candidate_audit.tsv", regenerated / "segmentation_candidate_audit.tsv")
    report = {"verified": True, "input_sha256": inputs, "tables": results,
              "numeric_tolerance": {"rtol": 1e-12, "atol": 1e-12}, "matched_manifest_fields": fields,
              "prototype": {"seed": 20260907, "replicates": 200, "segment_level": "median"},
              "screening": {"seed": 42, "replicates": 200, "segment_level": "arithmetic mean",
                            "all_source_sampling": "within-cluster resampling with fixed cluster-size weights"},
              "candidate_audit": "archived original; AD caller not rerun",
              "source_sha256": {str(p.relative_to(HERE)): sha256(p) for folder in (HERE / "analysis", HERE / "src") for p in folder.rglob("*.py")},
              "figure_count": len(list((regenerated / "figures").glob("*.png"))) + len(list((HERE / "screening-regenerated/figures").glob("*.png")))}
    marker.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"verified": True, "tables": len(results), "figures": report["figure_count"],
                      "byte_identical_tables": sum(t["byte_identical"] for t in results.values())}, indent=2))


if __name__ == "__main__":
    main()
