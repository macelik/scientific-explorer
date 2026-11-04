"""Experiment 00: input audit -> ../../input_audit.json (+ report.html here).

Recomputes SHA-256 for every data-package file and checks it against the handoff manifest,
verifies the pinned upstream source (git commit and per-file hashes), records dimensions,
units and alignment between AnnData, MuData and HA/HB/AT, and lists what is available or
missing for each preprocessing stage. Read-only; deterministic.
"""
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

HANDOFF = C.REPO / "research" / "2026-11-04-cloud-review-handoff"


def main():
    p = C.default_paths()
    man = json.loads((HANDOFF / "DATA_MANIFEST.json").read_text())
    pkg = json.loads((HANDOFF / "PACKAGE.json").read_text())
    files = []
    for rec in man["files"]:
        fp = p["root"] / rec["path"]
        h = C.sha256(fp) if fp.is_file() else None
        files.append({"path": rec["path"], "size_bytes": fp.stat().st_size if fp.is_file() else None, "sha256": h,
                      "matches_manifest": h == rec["sha256"]})
    up = []
    for name, want in man["upstream"]["package_files_sha256"].items():
        fp = p["package"] / name
        up.append({"file": name, "sha256": C.sha256(fp), "matches_manifest": C.sha256(fp) == want})
    commit = subprocess.run(["git", "-C", str(p["package"].parent), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    d = C.H5AD(p["h5ad"])
    with h5py.File(p["h5mu"], "r") as f:
        mu_cells = f["obs/_index"][...].astype(str)
        mu_obs = list(f["obs"].keys())
        dep = f["mod/depth/X"]
        dep_sample = dep["data"][:100000]
        mods = {m: {"shape": [int(x) for x in f[f"mod/{m}/X"].attrs["shape"]], "layers": list(f[f"mod/{m}"].get("layers", {}).keys()),
                    "uns": list(f[f"mod/{m}/uns"].keys())} for m in f["mod"]}
    lab = C.wnn_labels(p["h5mu"])
    at_cols = [c for c in pq.ParquetFile(p["allelic"] / "AT.parquet").schema_arrow.names if c not in ("chr", "start", "end", "__index_level_0__")]
    e6 = json.loads((HERE.parent / "06_allelic_information" / "summary.json").read_text())["alignment_checks"]
    e1 = json.loads((HERE.parent / "01_reconstruct_x" / "summary.json").read_text())
    e2 = json.loads((HERE.parent / "02_zero_clip_audit" / "summary.json").read_text())
    rb = C.shared_cells()
    audit = {
        "scope": "Inputs actually used by the 2026-11-04 cloud method review. Byte identity and alignment checks; not scientific validation.",
        "workbench": {"repository": "macelik/scientific-explorer", "branch_reviewed": "allelic-segmentation",
                      "checkout_commit": "6d4f48d79e0545fd3fb7ad84b86dbf659b57f107", "earlier_research_baseline": "994fb3afcf4caa93f6b52af7267d3668d469eec2",
                      "work_branch": "claude/eager-thompson-y0nc0t"},
        "data_package": {"release_tag": "cloud-review-data-2026-11-04", "archive_sha256": pkg.get("archive", pkg).get("sha256") if isinstance(pkg.get("archive"), dict) else pkg.get("sha256"),
                         "archive_checksum_verified": True, "download_route": "curl from the GitHub release asset URLs (gh release download is blocked in this session: GraphQL 403)",
                         "files": files, "all_files_match_manifest": all(f["matches_manifest"] for f in files)},
        "upstream_source": {"repository": "https://github.com/colomemaria/pyEpiAneufinder", "pinned_commit": "e8305ea2162634d07e92a3bb401df672e48f217b",
                            "checked_out_commit": commit, "pyproject_version": "0.3.5", "files": up,
                            "all_files_match_manifest": all(u["matches_manifest"] for u in up),
                            "license_note": "pyEpiAneufinder pyproject declares BSD-3-clause; earlier research recorded GPL-3.0-only for the repository. Not resolved here; see sources.json."},
        "anndata": {"path": "pyepi_results/count_matrix.h5ad", "shape": [d.n_cells, d.n_bins],
                    "X": "processed depth: GC-corrected (per-cell LOESS), np.rint to int, negatives set to 0, each row scaled to sum 1e5 (float64, CSR)",
                    "layers/raw": "post-QC, pre-GC integer fragment counts on retained cells and bins (int64, CSR)",
                    "explicit_zeros_stored": e1["explicit_zeros_stored"], "nnz": e1["nnz"],
                    "var": "seq (categorical), start (1-based), end (inclusive), GC", "uns": "empty", "obs": "cellID only",
                    "bins_per_chromosome": pd.Series(d.seq).value_counts(sort=False).to_dict()},
        "mudata": {"path": "integration-story/data/integrated.h5mu", "cells": int(len(mu_cells)), "obs_columns": mu_obs, "modalities": mods,
                   "depth_X_note": "log1p-transformed depth (uns has log1p; sample values non-integer, "
                                   f"min {float(dep_sample.min()):.3g}, max {float(dep_sample.max()):.3g}); not interchangeable with AnnData X",
                   "cluster_key_used": "wnn_leiden_0.3", "cluster_sizes": lab.value_counts().sort_index().to_dict(),
                   "cell_set_equal_anndata": set(mu_cells) == set(d.cells), "cell_order_equal_anndata": bool(np.array_equal(mu_cells, d.cells))},
        "allelic": {"files": ["aligned_matrices/HA.parquet", "aligned_matrices/HB.parquet", "aligned_matrices/AT.parquet"],
                    "unit": "sum over selected phased heterozygous SNP sites of per-site DP split by haplotype; a fragment can contribute at several sites; not unique fragments",
                    "dtype": "float32 stored, integer valued", "checks": e6, "join": "by barcode name (column order differs from AnnData)",
                    "barcodes": len(at_cols)},
        "frozen_panel": {"path": "research/2026-10-29-low-coverage/shared/cells.tsv", "cells": int(len(rb)),
                         "all_in_anndata": bool(set(rb.cell) <= set(d.cells)), "per_tercile": rb.coverage_tercile.value_counts().to_dict()},
        "stage_availability": [
            {"stage": "source fragments / BAM", "available": False, "consequence": "cannot recount, change bin size, deduplicate or audit read filters"},
            {"stage": "full bin grid after blacklist (binned_genome.csv)", "available": True, "consequence": "QC-removed bins can be located (experiment 02) but not their counts"},
            {"stage": "pre-QC cell x bin matrix", "available": False, "consequence": "relaxed bin/cell QC cannot be evaluated; removed regions (e.g. CDKN2A/B) stay unobservable"},
            {"stage": "QC masks / run parameters (parameter_configuration.txt)", "available": False, "consequence": "historical thresholds are inferred from source defaults, not records"},
            {"stage": "intermediate GC fits", "available": "reconstructed", "consequence": "rebuilt with pinned source; 3,449/3,460 cells bit-identical (experiment 01)"},
            {"stage": "pre-clipping floating values", "available": "reconstructed", "consequence": "available only from the rebuild, not historical"},
            {"stage": "processed counts (X)", "available": True, "consequence": "used as the baseline estimand"},
            {"stage": "SNP selection, genotype calls, per-SNP counts, phase blocks", "available": False, "consequence": "allelic tests are bin-level and phase-agnostic; unique-molecule counting unresolved"},
            {"stage": "independent region-level DNA truth", "available": False, "consequence": "published karyotype loci used only as weak, prespecified positive controls"},
        ],
        "historical_uncertainties": [
            "17 cell-bins zero in saved X but large in every tested rebuild environment (experiment 01); historical environment unrecorded.",
            "Historical run parameters (minFrags, thresholds, k for the production per-cell run) not recorded in the package; source defaults assumed.",
            "Recorded handoff environment lists numpy 1.24.4 with scikit-misc 0.5.2, whose metadata require numpy >= 1.26.4: that environment cannot be resolved as recorded (installed here with --no-deps).",
            "Pileup tool/version, UMI/read mode, SNP list and genotyping thresholds behind HA/HB/AT unknown.",
            f"{e2['removed_by_bin_filter']} autosomal bins removed by the bin filter; their contents are unknown without pre-QC counts.",
        ],
        "missing_prerequisites_for_candidates": {
            "Numbat-multiome ATAC-bin": ["fragment file(s) for binning", "SNP-level allele counts from pileup_and_phase", "binned ATAC reference from normal samples"],
            "Alleloscope (scATAC)": ["SNP-level allele counts per cell", "region/segment definitions from DNA (sample-matched bulk DNA recommended)"],
            "raw-count segmentation with exposure": ["none for the retained grid; pre-QC counts needed to include removed bins"],
        },
        "environment": C.environment(),
    }
    C.write_json(C.REVIEW / "input_audit.json", audit)
    frag = [("Result", f"<div class='card'>Data files matching manifest: <b>{sum(f['matches_manifest'] for f in files)}/{len(files)}</b>. "
                       f"Upstream files matching manifest: <b>{sum(u['matches_manifest'] for u in up)}/{len(up)}</b> at commit <code>{commit}</code>. "
                       f"MuData cells equal AnnData cells: {audit['mudata']['cell_set_equal_anndata']} (order equal: {audit['mudata']['cell_order_equal_anndata']}). "
                       f"Allelic checks: {e6}.</div>"),
            ("Stage availability", C.html_table(pd.DataFrame(audit["stage_availability"]))),
            ("Historical uncertainties", "<ul>" + "".join(f"<li>{u}</li>" for u in audit["historical_uncertainties"]) + "</ul>"),
            ("Files", C.html_table(pd.DataFrame(files)))]
    C.write_html_report(HERE / "report.html", "Input audit", "Experiment 00 · identities, dimensions, units and alignment", frag,
                        {"command": "python -I experiments/00_input_audit/run.py", "output": "../../input_audit.json"})
    print(audit["data_package"]["all_files_match_manifest"], audit["upstream_source"]["all_files_match_manifest"], commit,
          audit["mudata"]["cell_set_equal_anndata"], audit["mudata"]["cell_order_equal_anndata"])


if __name__ == "__main__":
    main()
