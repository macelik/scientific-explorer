"""Experiment 02: QC-removed bins (retained-grid gaps) and zeros introduced by preprocessing.

1. Compares the saved full bin grid (pyepi_results/binned_genome.csv, written by
   make_windows *after* blacklist exclusion) with the retained 25,337-bin AnnData grid.
   Autosomal bins present in the grid but absent from the AnnData were removed by the
   cohort nonzero-frequency bin filter (pyEpiAneufinder.py L212-L215, threshold
   (1-0.85) x n_cells) - the only bin-removing step after windowing. Runs of removed bins
   are reported with N content and GC, and checked against prespecified published
   MCF10-A loci (CDKN2A/B homozygous deletion; see sources.json).
2. Uses experiment-01 per-cell outputs to stratify raw-positive -> 0 conversions by raw
   coverage tercile (from the frozen 300-cell panel membership and all cells).
Deterministic. Writes TSV/JSON + report.html here.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

# Prespecified from literature before looking at these bins' data (hg38 gene coordinates, NCBI RefSeq):
# CDKN2A chr9:21,967,752-21,995,301; CDKN2B chr9:22,002,903-22,009,313; MTAP chr9:21,802,636-21,941,115.
LOCI = [("CDKN2A/CDKN2B (published MCF10A homozygous deletion)", "chr9", 21_802_636, 22_009_313)]


def main():
    p = C.default_paths()
    d = C.H5AD(p["h5ad"])
    g = pd.read_csv(p["root"] / "pyepi_results" / "binned_genome.csv", index_col=0)
    g = g[g.chromosome.isin(sorted(set(d.seq)))].copy()
    key = set(zip(d.seq, d.start))
    g["retained"] = [(c, s) in key for c, s in zip(g.chromosome, g.start)]
    assert g.retained.sum() == d.n_bins, "retained bins not all found in the full grid"
    # runs of removed bins
    runs = []
    for ch, sub in g.groupby("chromosome", sort=False):
        sub = sub.sort_values("start")
        rem = (~sub.retained).to_numpy()
        i = 0
        while i < len(sub):
            if rem[i]:
                j = i
                while j + 1 < len(sub) and rem[j + 1]:
                    j += 1
                part = sub.iloc[i:j + 1]
                runs.append({"chrom": ch, "start": int(part.start.iloc[0]), "end": int(part.end.iloc[-1]), "n_bins": j - i + 1,
                             "mean_N": float(part.N.mean()), "mean_GC": float(part.GC.mean()),
                             "assembly_gap_like": bool(part.N.mean() > 0.5)})
                i = j + 1
            else:
                i += 1
    runs = pd.DataFrame(runs)
    runs.to_csv(HERE / "qc_removed_bin_runs.tsv", sep="\t", index=False)
    loci = []
    for name, ch, s, e in LOCI:
        w = g[(g.chromosome == ch) & (g.end >= s) & (g.start <= e)]
        loci.append({"locus": name, "chrom": ch, "start": s, "end": e, "grid_bins": len(w), "retained": int(w.retained.sum()),
                     "removed_by_bin_filter": int((~w.retained).sum()), "N_content": float(w.N.mean())})
    loci = pd.DataFrame(loci)
    reg = C.shared_regions()
    win = []
    for _, r in reg[reg.part == "target"].iterrows():
        w = g[(g.chromosome == r.chrom) & (g.start >= r.requested_start_1based) & (g.end <= r.requested_end_inclusive)]
        win.append({"window": r.parent, "grid_bins": len(w), "retained": int(w.retained.sum()),
                    "removed_non_gap": int(((~w.retained) & (w.N < 0.5)).sum()), "removed_gap_like": int(((~w.retained) & (w.N >= 0.5)).sum())})
    win = pd.DataFrame(win)

    # zeros by coverage
    cells = pd.read_csv(HERE.parent / "01_reconstruct_x" / "per_cell.tsv.gz", sep="\t")
    cells["tercile"] = pd.qcut(cells.raw_total, 3, labels=["low", "mid", "high"])
    cells["frac_rawpos_to_zero"] = cells.n_rawpos_X0_obs / cells.raw_nnz
    zt = cells.groupby("tercile", observed=True).agg(cells=("cell", "size"), median_raw_total=("raw_total", "median"),
                                                     median_frac_rawpos_to_zero=("frac_rawpos_to_zero", "median"),
                                                     max_frac=("frac_rawpos_to_zero", "max"),
                                                     total_negclip=("n_rawpos_negclip_to_0", "sum")).reset_index()
    panel = C.shared_cells()
    pc = cells.merge(panel[["cell", "coverage_tercile"]], on="cell")
    zp = pc.groupby("coverage_tercile").frac_rawpos_to_zero.median().rename("panel_median_frac_rawpos_to_zero").reset_index()
    pb = pd.read_csv(HERE.parent / "01_reconstruct_x" / "per_bin.tsv.gz", sep="\t")
    pb["frac_rawpos_to_zero"] = pb.n_cells_raw_pos_X0_obs / pb.n_cells_raw_pos.clip(lower=1)
    pb["gc_decile"] = pd.qcut(pb.gc, 10, labels=False)
    zg = pb.groupby("gc_decile").agg(gc_lo=("gc", "min"), gc_hi=("gc", "max"), mean_raw=("mean_raw", "mean"),
                                     frac_rawpos_to_zero=("frac_rawpos_to_zero", "mean")).reset_index()

    res = {"experiment": "02_zero_clip_audit", "deterministic": True,
           "autosomal_grid_bins_after_blacklist": int(len(g)), "retained": int(g.retained.sum()), "removed_by_bin_filter": int((~g.retained).sum()),
           "removed_runs": int(len(runs)), "removed_runs_assembly_gap_like": int(runs.assembly_gap_like.sum()),
           "removed_bins_non_gap": int(runs.loc[~runs.assembly_gap_like, "n_bins"].sum()),
           "longest_non_gap_runs": runs[~runs.assembly_gap_like].nlargest(10, "n_bins").to_dict("records"),
           "prespecified_loci": loci.to_dict("records"), "windows": win.to_dict("records"),
           "zeros_by_coverage_tercile_all_cells": zt.astype({"tercile": str}).to_dict("records"),
           "environment": C.environment()}
    C.write_json(HERE / "summary.json", res)

    fig, ax = plt.subplots(1, 2, figsize=(13, 3.6))
    ax[0].scatter(cells.raw_total, 100 * cells.frac_rawpos_to_zero, s=3, alpha=.4, color="#2f5d8a")
    ax[0].set_xscale("log"); ax[0].set_xlabel("raw fragments in retained bins"); ax[0].set_ylabel("% raw-positive bins → 0 in X")
    ax[0].set_title("Rounding/clipping zeros vs cell coverage")
    sub = g[(g.chromosome == "chr9") & (g.start > 15e6) & (g.start < 32e6)]
    ax[1].scatter(sub.start / 1e6, np.where(sub.retained, 1, 0), s=8, c=np.where(sub.retained, "#2f5d8a", "#9b2c2c"))
    ax[1].axvspan(21.80, 22.01, color="#e0b45c", alpha=.4, label="CDKN2A/B (published HD)")
    for _, r in reg[(reg.part == "target") & (reg.chrom == "chr9")].iterrows():
        ax[1].hlines(0.5 + 0.1 * (r.n_bins > 30), r.requested_start_1based / 1e6, r.requested_end_inclusive / 1e6, color="k")
    ax[1].set_yticks([0, 1]); ax[1].set_yticklabels(["removed", "retained"]); ax[1].set_xlabel("chr9 (Mb)")
    ax[1].legend(fontsize=8); ax[1].set_title("chr9p bin retention; black = saved windows")
    b64 = C.fig_to_b64(fig); plt.close(fig)
    frag = [
        ("Question", "<p>Which bins did the cohort nonzero-frequency filter remove, could a real event remove itself, and where do "
         "rounding/clipping zeros concentrate?</p>"),
        ("Published locus check (prespecified)", C.html_table(loci) +
         "<div class='card'>The bins covering CDKN2A/CDKN2B, reported as homozygously deleted in MCF10-A, are <b>absent from the "
         "retained grid</b>. They are not blacklisted (they are in the post-blacklist grid) and contain no N bases, so the "
         "≥15%-nonzero-cells bin filter removed them. The pipeline cannot call this event; it becomes a gap between retained neighbours. "
         "Any clonal homozygous deletion (or a region inaccessible in at least 85% of cells) behaves the same way, in any dataset.</div>"),
        ("Saved windows", C.html_table(win)),
        ("Removed-bin runs (longest non-assembly-gap)", C.html_table(runs[~runs.assembly_gap_like].nlargest(25, "n_bins"))),
        ("Zeros introduced by rounding/clipping", C.html_table(zt) + C.html_table(zp) + C.html_table(zg) + C.img(b64, "zeros and chr9 retention")),
        ("Interpretation", "<ul><li>Removal by the bin filter is outcome-dependent: the stronger and more clonal the loss, the more "
         "likely its bins disappear. A deletion-aware method needs the pre-filter grid with an explicit 'removed: too few nonzero "
         "cells' state, not silent gaps. Without pre-QC counts we cannot say what the other removed runs contain.</li>"
         "<li>Rounding zeros depend on coverage: low-coverage cells lose a larger share of raw-positive bins. That is a "
         "coverage-dependent distortion of the processed measurement, avoided by modelling raw counts with exposure.</li></ul>"),
    ]
    C.write_html_report(HERE / "report.html", "QC gaps and introduced zeros", "Experiment 02 · retained-grid removal and rounding/clipping zeros",
                        frag, {"command": "python -I experiments/02_zero_clip_audit/run.py", "summary": "summary.json", "environment": C.environment()})
    print({k: res[k] for k in ["removed_by_bin_filter", "removed_runs", "removed_runs_assembly_gap_like", "removed_bins_non_gap", "prespecified_loci", "windows"]})
    print(zt)


if __name__ == "__main__":
    main()
