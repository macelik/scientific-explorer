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
    low, mid, high = (zt.set_index("tercile").loc[k] for k in ("low", "mid", "high"))
    nongap = runs[~runs.assembly_gap_like]
    size_bins = pd.cut(nongap.n_bins, [0, 1, 4, 9, 1000], labels=["1 bin", "2-4 bins", "5-9 bins", ">=10 bins"])
    size_tab = nongap.groupby(size_bins, observed=False).agg(runs=("n_bins", "size"), bins=("n_bins", "sum")).reset_index().rename(columns={"n_bins": "run_length"})
    chrom_tab = nongap.groupby("chrom").n_bins.sum().sort_values(ascending=False).head(8).reset_index().rename(columns={"n_bins": "removed_non_gap_bins"})
    cd = loci.iloc[0]
    n_cell_pos_zero = zg.frac_rawpos_to_zero
    win_notes = []
    for _, r in win.iterrows():
        rem = r.removed_non_gap + r.removed_gap_like
        if rem == 0:
            note = f"All {r.grid_bins} post-blacklist bins retained; no QC gap inside the window."
        else:
            note = (f"{rem} of {r.grid_bins} grid bins ({100*rem/r.grid_bins:.0f}%) were removed by the bin filter"
                    f"{' (none assembly-gap-like)' if r.removed_gap_like == 0 else ''}; the window contains unobservable stretches, so a deeper deletion could be hidden here.")
        win_notes.append(note)
    win_show = win.assign(removed_total=win.removed_non_gap + win.removed_gap_like)
    win_show["note"] = win_notes
    win_show = win_show[["window", "grid_bins", "retained", "removed_non_gap", "removed_gap_like", "note"]]
    win_html = C.html_table(win_show).replace("<td>All ", "<td style='text-align:left;white-space:normal'>All ")
    verdict = (f"The bin filter removed {res['removed_by_bin_filter']:,} of {res['autosomal_grid_bins_after_blacklist']:,} autosomal grid bins, including {int(cd.removed_by_bin_filter)} of the {int(cd.grid_bins)} bins overlapping the published CDKN2A/B homozygous deletion, "
               f"so that event cannot be called; separately, rounding/clipping zeros {100*float(low.median_frac_rawpos_to_zero):.2f}% of raw-positive bins in low-coverage cells against {100*float(high.median_frac_rawpos_to_zero):.2f}% in high-coverage cells.")
    findings = [
        f"<b>{res['removed_by_bin_filter']:,}</b> of {res['autosomal_grid_bins_after_blacklist']:,} autosomal post-blacklist bins ({100*res['removed_by_bin_filter']/res['autosomal_grid_bins_after_blacklist']:.1f}%) are missing from the {res['retained']:,}-bin retained grid; they form {res['removed_runs']} runs, of which only {res['removed_runs_assembly_gap_like']} are assembly-gap-like (mean N &ge; 0.5).",
        f"{res['removed_bins_non_gap']:,} removed bins are in runs that are not assembly-gap-like; the longest is {nongap.nlargest(1,'n_bins').iloc[0].chrom}:{nongap.nlargest(1,'n_bins').iloc[0].start/1e6:.1f}&ndash;{nongap.nlargest(1,'n_bins').iloc[0].end/1e6:.1f} Mb ({int(nongap.n_bins.max())} bins, mean GC {nongap.nlargest(1,'n_bins').iloc[0].mean_GC:.2f}).",
        f"Prespecified check: {int(cd.removed_by_bin_filter)} of {int(cd.grid_bins)} grid bins overlapping CDKN2A/CDKN2B (chr9:{int(cd.start):,}&ndash;{int(cd.end):,}) were removed, none with N content ({cd.N_content:.1f}); the published homozygous deletion is not callable in this pipeline.",
        f"Two chr9 saved windows contain removed bins: chr9_20_27 ({int(win.set_index('window').loc['chr9_20_27','removed_non_gap'])} of {int(win.set_index('window').loc['chr9_20_27','grid_bins'])}) and chr9_22.5_26.6 ({int(win.set_index('window').loc['chr9_22.5_26.6','removed_non_gap'])} of {int(win.set_index('window').loc['chr9_22.5_26.6','grid_bins'])}); chr7_51.4_53.5 has {int(win.set_index('window').loc['chr7_51.4_53.5','removed_non_gap'])} of {int(win.set_index('window').loc['chr7_51.4_53.5','grid_bins'])}; the other four windows have none.",
        f"Rounding/clipping zeros depend on coverage: median fraction of raw-positive bins that become 0 is {100*float(low.median_frac_rawpos_to_zero):.2f}% (low tercile, median {float(low.median_raw_total):,.0f} raw fragments), {100*float(mid.median_frac_rawpos_to_zero):.2f}% (mid) and {100*float(high.median_frac_rawpos_to_zero):.2f}% (high, median {float(high.median_raw_total):,.0f}); maximum {100*float(low.max_frac):.2f}% in the low tercile.",
        f"Negative-value clipping contributes {int(zt.total_negclip.sum()):,} zeroed raw-positive cell-bins in total ({int(low.total_negclip)} low, {int(mid.total_negclip)} mid, {int(high.total_negclip)} high tercile), concentrated in the highest-coverage cells.",
        f"By bin GC decile, the mean share of raw-positive cells whose bin becomes 0 is {100*n_cell_pos_zero.iloc[0]:.3f}% in the lowest decile and {100*n_cell_pos_zero.iloc[1:9].max():.3f}% at most in deciles 1-8, but {100*n_cell_pos_zero.iloc[9]:.2f}% in the highest-GC decile (GC {zg.gc_lo.iloc[9]:.2f}&ndash;{zg.gc_hi.iloc[9]:.2f}, mean raw count {zg.mean_raw.iloc[9]:.2f}): the zeros concentrate in GC-rich bins, where the correction multiplier is below 1.",
    ]
    sections = [
        ("Question and why it matters",
         "<p><b>Question.</b> Which bins did the cohort nonzero-frequency filter remove, could a real event remove itself, and where do rounding and clipping zeros concentrate?</p>"
         "<p><b>Why it matters.</b> The recommendation (rank 2) asks to keep QC-removed bins as an explicit &lsquo;unobservable&rsquo; state, because a clonal homozygous deletion has zero counts in every cell and therefore fails the &ge;15%-nonzero-cells bin filter. "
         "It also argues for a raw-count, exposure-aware model, partly because rounding zeros depend on cell coverage in the processed matrix. Both claims need to be measured here.</p>"),
        ("Inputs",
         "<ul><li><b>Full grid.</b> <code>pyepi_results/binned_genome.csv</code> (written by make_windows after blacklist exclusion; columns chromosome, start, end, GC, N), restricted to the chromosomes present in the AnnData (autosomes). "
         f"{res['autosomal_grid_bins_after_blacklist']:,} bins of 100 kb.</li>"
         f"<li><b>Retained grid.</b> The {res['retained']:,} bins of <code>count_matrix.h5ad</code>; every retained bin was matched to the full grid by (chromosome, start).</li>"
         "<li><b>Experiment-01 tables.</b> <code>per_cell.tsv.gz</code> (3,460 cells: raw total, non-zeros, raw-positive-to-zero counts) and <code>per_bin.tsv.gz</code> (per-bin raw-positive and zero counts, GC).</li>"
         "<li><b>Frozen shared inputs.</b> Saved-window coordinates from <code>2026-10-29-low-coverage/shared/regions.tsv</code> (7 target windows) and the 300-cell panel (<code>cells.tsv</code>) for the coverage terciles.</li>"
         "<li><b>Prespecified locus.</b> CDKN2A/CDKN2B/MTAP, chr9:21,802,636&ndash;22,009,313 (hg38 RefSeq gene coordinates, fixed before looking at the bins).</li>"
         "<li><b>Not available.</b> Pre-QC counts: what the removed bins contained (deletion, inaccessible desert, low mappability) cannot be determined; cells removed by cell filters are also absent.</li></ul>"),
        ("Method",
         "<ol><li>Mark each full-grid bin as retained or removed. By pyEpi's code the bin filter (<code>nnz(bin) &ge; (1 - 0.85) &middot; n_cells</code>, L212&ndash;L215) is the only bin-removing step after windowing, so every absent autosomal bin failed it.</li>"
         "<li>Group consecutive removed bins into runs; report each run&rsquo;s mean N fraction and GC. A run with mean N &ge; 0.5 is called assembly-gap-like.</li>"
         "<li>Count removed bins overlapping the prespecified CDKN2A/B interval and inside each saved window.</li>"
         "<li>Stratify raw-positive-to-zero conversions per cell by raw-coverage tercile (<code>pd.qcut</code> on raw total, 3 groups), and per bin by GC decile.</li></ol>"
         "<pre>frac_rawpos_to_zero(cell) = n_rawpos_X0_obs / raw_nnz\nfrac_rawpos_to_zero(bin)  = n_cells_raw_pos_X0_obs / max(n_cells_raw_pos, 1)\nbin filter keeps bin  iff  nnz(bin) &gt;= 0.15 * n_cells        # 0.15 = 1 - 0.85</pre>"
         "<p>Deterministic: no random numbers or thresholds were tuned in this experiment.</p>"),
        ("Results", "", 2),
        ("Which bins were removed by the bin filter",
         C.html_table(size_tab) + C.caption("Non-assembly-gap removed runs grouped by length: <code>runs</code> is the number of runs, <code>bins</code> the number of removed bins. Many single-bin runs mean scattered removal; long runs mean extended low-signal territory.") +
         C.html_table(chrom_tab) + C.caption("The eight chromosomes with the most removed non-assembly-gap bins.") +
         C.html_table(nongap.nlargest(15, "n_bins")) + C.caption("The 15 longest removed runs that are not assembly-gap-like (<code>mean_N</code> and <code>mean_GC</code> are bin means within the run). Full table of all runs: <code>qc_removed_bin_runs.tsv</code>."), 3),
        ("Published locus check (prespecified)",
         C.html_table(loci) +
         C.caption("<code>grid_bins</code> are the full-grid bins overlapping the interval; <code>removed_by_bin_filter</code> were absent from the retained grid. Zero N content means the removal is not an assembly gap.") +
         "<div class='card'>The bins covering CDKN2A/CDKN2B, reported as homozygously deleted in MCF10-A, are <b>mostly absent from the retained grid</b>. They are not blacklisted (they are in the post-blacklist grid) and contain no N bases, so the "
         "&ge;15%-nonzero-cells bin filter removed them. The pipeline cannot call this event; it becomes a gap between retained neighbours. "
         "Any clonal homozygous deletion (or a region inaccessible in at least 85% of cells) behaves the same way, in any dataset. Because there are no pre-QC counts, "
         "this is consistent with, but does not demonstrate, that the removal was caused by the deletion itself.</div>", 3),
        ("Zeros introduced by rounding and clipping",
         C.html_table(zt) + C.caption("Cells split into raw-coverage terciles (all 3,460 cells). <code>median_frac_rawpos_to_zero</code> is the per-cell fraction of raw-positive bins that are 0 in X, medianed within the tercile; <code>total_negclip</code> counts raw-positive cell-bins zeroed by negative clipping. Expect a decreasing median from low to high coverage.") +
         C.html_table(zp) + C.caption("The same median restricted to the frozen 300-cell panel (100 per tercile) used in later experiments, so panel results can be related to the whole matrix.") +
         C.html_table(zg) + C.caption("Bins split into GC deciles: <code>mean_raw</code> is the mean raw count per cell, <code>frac_rawpos_to_zero</code> the mean over bins of the share of raw-positive cells whose bin is 0 in X. Look for the decile where the share is far above the others.") +
         C.img(b64, "zeros and chr9 retention") +
         C.caption("Left: one point per cell, percentage of raw-positive bins that become 0 against cell coverage (log x); a downward trend means low-coverage cells lose more. Right: chr9 15&ndash;32 Mb, each point a grid bin (blue retained, red removed), "
                   "amber band the CDKN2A/B interval, black lines the saved windows; look for red points under the amber band and inside the windows."), 3),
        ("Per-region detail",
         win_html + C.caption("One row per saved target window. <code>removed_non_gap</code> / <code>removed_gap_like</code> are bins missing from the retained grid by N content; the note is derived from these counts. Retention of a window says nothing about its copy number.")),
        ("Sensitivity",
         "<p>No parameters were varied. The bin-filter threshold (15% nonzero cells) and the assembly-gap cut (mean N 0.5) are fixed by the pinned code and by this report&rsquo;s labelling convention; the number of &lsquo;non-gap&rsquo; bins would change with the N cut. "
         "Coverage-tercile results are given both for all cells and the 300-cell panel.</p>"),
        ("Limits and alternative explanations",
         "<ul><li>Removal is outcome-dependent but not diagnostic: a removed run may be a clonal homozygous deletion, a region inaccessible in at least 85% of cells, or low mappability. Without pre-QC counts these cannot be told apart.</li>"
         "<li>The cell filters (stages 4&ndash;5) are not audited; removed cells are not in the package.</li>"
         "<li>The deleted-locus argument rests on published coordinates and a single locus; it shows the mechanism exists, not how often it occurs.</li>"
         "<li>Raw-positive-to-zero fractions are properties of this matrix and its depth; they are not a measure of bias in any called segment.</li></ul>"),
    ]
    glossary = {
        "chrom / start / end": "Chromosome and 1-based window coordinates (hg38, 100 kb bins).",
        "n_bins": "Number of consecutive removed bins in a run.",
        "mean_N / mean_GC": "Mean fraction of N bases / of G+C bases across the bins of a run or locus.",
        "assembly_gap_like": "True if a removed run has mean N &ge; 0.5, i.e. mostly assembly gap.",
        "grid_bins / retained": "Full-grid bins overlapping a locus or window / those present in the retained AnnData grid.",
        "removed_by_bin_filter": "Grid bins overlapping a locus that were removed by the nonzero-frequency bin filter.",
        "removed_non_gap / removed_gap_like": "Removed bins in a window whose own N fraction is below / at or above 0.5.",
        "N_content": "Mean N fraction of the bins overlapping a locus.",
        "tercile": "Third of cells by raw fragment total in retained bins (low / mid / high).",
        "median_raw_total": "Median raw fragments in retained bins within the tercile.",
        "median_frac_rawpos_to_zero / panel_median_frac_rawpos_to_zero": "Median over cells of the share of raw-positive bins that are 0 in X, for all cells / for the 300-cell panel.",
        "max_frac": "Largest per-cell fraction within the tercile.",
        "total_negclip": "Raw-positive cell-bins zeroed because the corrected value was negative.",
        "gc_decile / gc_lo / gc_hi": "Bin GC decile (0 = lowest) and its GC range.",
        "mean_raw": "Mean raw count per cell, averaged over the bins in the decile.",
        "frac_rawpos_to_zero": "Share of raw-positive cell-bins that are 0 in X (per cell or per bin, as labelled).",
        "run_length / runs / bins": "Length class of removed runs, number of runs in the class, and bins they contain.",
    }
    C.write_html_report(HERE / "report.html", "QC gaps and introduced zeros", "Experiment 02 &middot; retained-grid removal and rounding/clipping zeros",
                        sections, {"command": "python -I experiments/02_zero_clip_audit/run.py", "summary": "summary.json", "environment": C.environment()},
                        key_findings=findings, verdict=verdict, glossary=glossary,
                        supports=["The bin filter removes a sizeable set of autosomal bins, mostly in runs that are not assembly gaps, including most bins at the published CDKN2A/B homozygous deletion.",
                                  "Processed X carries a coverage-dependent loss of raw-positive bins through rounding, larger in low-coverage cells.",
                                  "Two chr9 saved windows and one chr7 window contain removed bins, so parts of them are unobservable."],
                        not_supports=["What the other removed runs contain (deletion, inaccessibility, low mappability): pre-QC counts are not available.",
                                      "That rounding zeros explain any saved dip: later experiments show the saved contrasts are coverage-invariant in raw counts.",
                                      "Any statement about copy number at the removed bins."],
                        decision_impact="Supports recommendation ranks 2 and 3: carry bins that QC removes as an explicit &lsquo;unobservable&rsquo; state instead of a silent gap, and relax or log the bin filter on a sample with pre-QC counts. "
                                        "The coverage dependence of zeros is one of the arguments for a raw-count, exposure-aware model (rank 4), but it does not change the conclusion that the next change belongs in evidence integration.",
                        next_steps=["Obtain the pre-filter count matrix and QC masks to see what the removed runs contain.",
                                    "Compare removed runs with a read-length-matched mappability track to separate low-mappability from deletion.",
                                    "Repeat the CDKN2A/B-type check on a sample with an independently known homozygous deletion."])
    print({k: res[k] for k in ["removed_by_bin_filter", "removed_runs", "removed_runs_assembly_gap_like", "removed_bins_non_gap", "prespecified_loci", "windows"]})
    print(zt)


if __name__ == "__main__":
    main()
