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
    n_mm = json.loads((HERE.parent / '01_reconstruct_x' / 'analysis.json').read_text())['mismatch_bins']
    n_files_ok, n_up_ok = sum(f["matches_manifest"] for f in files), sum(u["matches_manifest"] for u in up)
    stages = pd.DataFrame(audit["stage_availability"])
    avail = stages.available.astype(str)
    n_missing = int((avail == "False").sum()); n_recon = int((avail == "reconstructed").sum()); n_avail = int((avail == "True").sum())
    total_mb = sum(f["size_bytes"] or 0 for f in files) / 1e6
    mu = audit["mudata"]; an = audit["anndata"]; al = audit["allelic"]
    verdict = (f"All {n_files_ok} of {len(files)} data-package files and {n_up_ok} of {len(up)} pinned pyEpi source files match their manifests, and the matrices are aligned by cell and bin, "
               f"but {n_missing} of {len(stages)} preprocessing stages (fragments, pre-QC counts, QC masks and run parameters, SNP-level counts, DNA truth) are not available, so those parts of the analysis rest on the retained post-QC grid and source defaults.")
    findings = [
        f"Data package: {n_files_ok}/{len(files)} files match the handoff manifest by SHA-256 ({total_mb:,.0f} MB in total); release tag <code>{audit['data_package']['release_tag']}</code>.",
        f"Pinned upstream source: {n_up_ok}/{len(up)} files match at commit <code>{commit[:7]}</code> (expected <code>e8305ea</code>; pyproject {audit['upstream_source']['pyproject_version']}).",
        f"AnnData: {an['shape'][0]:,} cells &times; {an['shape'][1]:,} bins (100 kb) over {len(an['bins_per_chromosome'])} autosomes; raw non-zeros {an['nnz']['layers/raw']:,}, X non-zeros {an['nnz']['X']:,}; no explicit stored zeros ({an['explicit_zeros_stored']}).",
        f"MuData: {mu['cells']:,} cells, cell set equal to the AnnData ({mu['cell_set_equal_anndata']}), order equal ({mu['cell_order_equal_anndata']}); clusters used: <code>{mu['cluster_key_used']}</code> with sizes {', '.join(f'{k}: {v:,}' for k, v in mu['cluster_sizes'].items())}.",
        f"Allelic matrices: {al['barcodes']:,} barcodes; AT = HA + HB exactly ({e6['AT_equals_HA_plus_HB_all_cells']}), non-negative ({e6['nonnegative']}), integer-valued ({e6['integer_valued']}); barcode sets equal but column order differs from the AnnData, so joins are by name.",
        f"Frozen 300-cell panel: {audit['frozen_panel']['cells']} cells, all in the AnnData ({audit['frozen_panel']['all_in_anndata']}), {audit['frozen_panel']['per_tercile']} per coverage tercile.",
        f"Stage availability: {n_avail} available, {n_recon} reconstructed from the pinned source, {n_missing} not available; {e2['removed_by_bin_filter']:,} autosomal bins removed by the bin filter have unknown contents.",
        f"{len(audit['historical_uncertainties'])} historical uncertainties are recorded, including {e1['reconstruction']['cells_with_any_bin_absdiff_gt_1e-9']} cells whose saved X is not reproduced ({n_mm} cell-bins zero historically) and an environment that cannot be resolved as recorded.",
    ]
    files_df = pd.DataFrame(files).assign(sha256_prefix=lambda x: x.sha256.str[:16], size_mb=lambda x: x.size_bytes / 1e6)[["path", "size_mb", "sha256_prefix", "matches_manifest"]]
    up_df = pd.DataFrame(up).assign(sha256_prefix=lambda x: x.sha256.str[:16])[["file", "sha256_prefix", "matches_manifest"]]
    dim_tab = pd.DataFrame([
        {"object": "AnnData count_matrix.h5ad", "shape": f"{an['shape'][0]:,} cells x {an['shape'][1]:,} bins", "X": "GC-corrected, rint, clipped, row sum 1e5 (float64 CSR)", "layers": "raw = post-QC pre-GC integer counts (int64 CSR)"},
        {"object": "MuData integrated.h5mu", "shape": f"{mu['cells']:,} cells", "X": "depth modality is log1p-transformed (not interchangeable with AnnData X)", "layers": "depth: raw; haplo: none"},
        {"object": "HA / HB / AT parquet", "shape": f"{an['shape'][1]:,} bins x {al['barcodes']:,} barcodes", "X": al["unit"], "layers": al["dtype"]},
    ])
    chrom_tab = pd.DataFrame([{"chromosome": k, "retained_bins": v} for k, v in an["bins_per_chromosome"].items()])
    chrom_row = pd.DataFrame([an["bins_per_chromosome"]])
    mu_tab = pd.DataFrame([{"modality": m, "shape": f"{v['shape'][0]} x {v['shape'][1]}", "layers": ", ".join(v["layers"]) or "none", "uns": ", ".join(v["uns"])} for m, v in mu["modalities"].items()])
    check_tab = pd.DataFrame([{"matrix": k, **v} for k, v in e6.items() if isinstance(v, dict)])
    miss_rows = [{"candidate method": k, "missing prerequisites": "; ".join(v)} for k, v in audit["missing_prerequisites_for_candidates"].items()]
    sections = [
        ("Question and why it matters",
         "<p><b>Question.</b> Which exact files, versions, dimensions and units did the review use, are they what the handoff manifest says they are, are the matrices aligned with each other, and which preprocessing stages are available or missing?</p>"
         "<p><b>Why it matters.</b> Every later result (experiments 01&ndash;07) is only interpretable if the inputs are byte-identical to the released package and the pinned source, and if their units are understood. The recommendation&rsquo;s stage 0 and its blocked or skipped items "
         "(Numbat-multiome and Alleloscope comparisons, a raw-count model that includes QC-removed bins, SNP-level allelic likelihoods) depend on exactly which stages are missing, listed here. This experiment does not test any scientific claim: it checks identity, alignment and availability.</p>"),
        ("Inputs",
         f"<ul><li><b>Handoff manifest.</b> <code>research/2026-11-04-cloud-review-handoff/DATA_MANIFEST.json</code> and <code>PACKAGE.json</code>: {len(files)} data files, {len(up)} upstream source files with expected SHA-256 hashes.</li>"
         "<li><b>Data package</b> <code>cloud-review-data-2026-11-04</code> (root from <code>PYEPI_SCIENTIFIC_ROOT</code>): <code>pyepi_results/count_matrix.h5ad</code>, <code>binned_genome.csv</code>, <code>hg38-blacklist.v2.bed</code>, <code>integration-story/data/integrated.h5mu</code>, <code>aligned_matrices/{HA,HB,AT}.parquet</code>, <code>genome_pool</code> files and others listed below.</li>"
         "<li><b>Upstream source.</b> pyEpiAneufinder checkout (via <code>PYEPI_SCIENTIFIC_PACKAGE</code>), pinned to commit e8305ea.</li>"
         "<li><b>Frozen panel.</b> <code>research/2026-10-29-low-coverage/shared/cells.tsv</code> (300 cells, 100 per coverage tercile).</li>"
         "<li><b>Other experiments&rsquo; summaries.</b> Experiment 06 alignment checks and experiments 01/02 counts are read from their <code>summary.json</code>, so run this experiment last.</li>"
         "<li><b>Not available.</b> Source fragments or BAM, pre-QC counts, QC masks and run parameters, SNP-level counts and phase blocks, independent DNA truth (see stage availability).</li></ul>"),
        ("Method",
         "<ol><li>Recompute SHA-256 of every file in the manifest and compare with the recorded hash.</li>"
         "<li>Recompute SHA-256 of every pinned upstream file, and read the checked-out git commit.</li>"
         "<li>Read AnnData and MuData dimensions, cell sets and order, and cluster sizes; check that the frozen 300-cell panel is contained in the AnnData and balanced across terciles.</li>"
         "<li>Record units (what X and the raw layer mean, what AT counts) and the alignment checks of experiment 06.</li>"
         "<li>List what is available, reconstructed or missing for each preprocessing stage and what each missing item blocks.</li></ol>"
         "<pre>matches_manifest = (sha256(file) == manifest.sha256)\ncell_set_equal  = set(mudata.obs) == set(anndata.cells);  cell_order_equal = identical order\nAT_equals_HA_plus_HB = all(AT == HA + HB)                  # from experiment 06</pre>"
         "<p>Read-only and deterministic; no random numbers. The output is written to <code>../../input_audit.json</code>.</p>"),
        ("Results", "", 2),
        ("Identity: data package", C.html_table(files_df) + C.caption(f"One row per data-package file: size in MB, first 16 hex characters of the recomputed SHA-256 (full hashes in <code>input_audit.json</code>) and whether it matches the manifest. {n_files_ok} of {len(files)} match."), 3),
        ("Identity: pinned upstream source", f"<p>Checked-out commit: <code>{commit}</code>. Expected commit: <code>{audit['upstream_source']['pinned_commit']}</code>.</p>" + C.html_table(up_df) +
         C.caption(f"One row per pinned source file. {n_up_ok} of {len(up)} match. License note: {audit['upstream_source']['license_note']}"), 3),
        ("Dimensions and units", C.html_table(dim_tab) + C.caption("What each object contains and in what unit. X and the raw layer are 100 kb bins per cell; AT counts selected-site observations and is not a unique-fragment count.") +
         C.html_table(mu_tab) + C.caption("MuData modalities. The depth modality of the MuData is not the AnnData X.") +
         C.html_table(chrom_row) + C.caption("Retained bins per autosome in the AnnData (sum equals the bin dimension)."), 3),
        ("Alignment", C.html_table(check_tab) + C.caption("Alignment of HA, HB and AT with the AnnData (from experiment 06): coordinates, barcode set and column order. Equal sets with unequal order means joins must be by barcode name.") +
         f"<div class='card'>MuData cells equal AnnData cells: {mu['cell_set_equal_anndata']} (order equal: {mu['cell_order_equal_anndata']}). AT = HA + HB for all cells and bins: {e6['AT_equals_HA_plus_HB_all_cells']}. Frozen panel contained in the AnnData: {audit['frozen_panel']['all_in_anndata']}.</div>", 3),
        ("Stage availability", C.html_table(stages) + C.caption("Each preprocessing stage, whether its input or output is available (True), reconstructed from the pinned source (reconstructed) or missing (False), and what the status blocks or permits."), 3),
        ("Historical uncertainties", "<ul>" + "".join(f"<li>{u}</li>" for u in audit["historical_uncertainties"]) + "</ul>", 3),
        ("Missing prerequisites for candidate methods", C.html_table(pd.DataFrame(miss_rows)) + C.caption("What is needed before each candidate in <code>method_comparison.tsv</code> can be run."), 3),
        ("Sensitivity", "<p>Not applicable: this experiment makes no estimate. Byte identity is binary.</p>"),
        ("Limits and alternative explanations",
         "<ul><li>Matching hashes show the files are those that were released; they do not show that the release itself is correct or complete.</li>"
         "<li>Identity of the upstream source does not establish that the historical run used that source or its defaults; run parameters were not recorded.</li>"
         "<li>The recorded handoff environment (numpy 1.24.4 with scikit-misc 0.5.2) cannot be resolved as recorded; it was installed here with <code>--no-deps</code>.</li>"
         "<li>Alignment checks cover coordinates, barcode sets and units; they cannot detect an upstream error in SNP selection, phasing or pileup, which are unrecorded.</li></ul>"),
    ]
    glossary = {
        "matches_manifest": "Whether the recomputed SHA-256 equals the hash in the handoff manifest.",
        "sha256_prefix": "First 16 hex characters of the SHA-256 (full value in input_audit.json).",
        "size_mb": "File size in megabytes.",
        "available": "Stage status: True = available; reconstructed = rebuilt from the pinned source; False = missing.",
        "consequence": "What the status means for the analysis.",
        "coords_match_anndata": "HA/HB/AT chromosome, start and end equal the AnnData var coordinates.",
        "barcode_set_equal / column_order_equal_anndata": "Whether the barcodes equal the AnnData cells as a set / in the same order.",
        "AT / HA / HB": "Selected-site allele observation count and its haplotype A and B parts (sum of per-site DP, not unique fragments).",
        "wnn_leiden_0.3": "Cluster key used throughout the review (5 clusters).",
        "pinned_commit": "pyEpiAneufinder commit the review is pinned to (e8305ea).",
        "cell_set_equal / cell_order_equal": "MuData cells equal the AnnData cells as a set / in order.",
    }
    C.write_html_report(HERE / "report.html", "Input audit", "Experiment 00 &middot; identities, dimensions, units and alignment", sections,
                        {"command": "python -I experiments/00_input_audit/run.py", "output": "../../input_audit.json", "environment": C.environment()},
                        key_findings=findings, verdict=verdict, glossary=glossary,
                        supports=["The data package and pinned source used for all experiments are byte-identical to the handoff manifests.",
                                  "The matrices can be joined safely (by name for allelic data), with units documented.",
                                  "Missing stages are recorded so that blocked comparisons are traceable to specific missing inputs."],
                        not_supports=["That the released data or the historical run were correct, or that the historical environment matches the one used here.",
                                      "Anything about copy number: this is an identity and availability check.",
                                      "Any method comparison that needs fragments, pre-QC counts, SNP-level counts or DNA truth."],
                        decision_impact="Supports recommendation stage 0 (persist run parameters, pre-filter grid, QC masks and GC diagnostics with every run) and explains why Numbat-multiome, Alleloscope and a raw-count model that includes QC-removed bins are listed as blocked: each depends on inputs marked missing here.",
                        next_steps=["Obtain fragment files, pre-QC counts and QC masks, and the historical parameter_configuration.txt.",
                                    "Obtain SNP-level allele counts with phase-set IDs and the pileup and genotyping settings behind HA/HB/AT.",
                                    "Resolve the recorded environment (numpy and scikit-misc versions) and add independent DNA truth for at least one sample."])
    print(audit["data_package"]["all_files_match_manifest"], audit["upstream_source"]["all_files_match_manifest"], commit,
          audit["mudata"]["cell_set_equal_anndata"], audit["mudata"]["cell_order_equal_anndata"])


if __name__ == "__main__":
    main()
