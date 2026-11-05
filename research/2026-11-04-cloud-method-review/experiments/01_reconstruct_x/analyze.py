"""Experiment 01b: explain reconstruction mismatches and quantify GC-correction amplification.

Reads experiment-01 outputs (per_cell.tsv.gz, per_bin.tsv.gz, scratch float intermediates)
and the saved AnnData. For each mismatching cell the historical integer vector is
recovered from saved X (X * T_hist / 1e5 with T_hist chosen so that all bins are integers,
estimated from bins that agree), and the differing bins are listed with their LOESS
prediction. Amplification = correction factor mean(raw)/fit for raw-positive bins, read
from the historical X directly where possible (X_int / raw), so it describes the saved
matrix, not only our rebuild. Writes TSV/JSON + report.html. Deterministic.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common as C  # noqa: E402

OUT = Path(__file__).resolve().parent


def env_section():
    files = sorted(f for f in OUT.glob("env_sensitivity_*.tsv") if "comparison" not in f.name)
    if not files:
        return "<p>Not run.</p>"
    mm = pd.read_csv(OUT / "mismatch_bins.tsv", sep="\t")[["row", "bin"]]
    base = pd.read_csv(OUT / "env_sensitivity_py310_np124_skmisc052.tsv", sep="\t")
    rows = []
    for f in files:
        t = pd.read_csv(f, sep="\t")
        m = base.merge(t, on=["row", "bin"], suffixes=("", "_o"))
        h = mm.merge(t, on=["row", "bin"])
        rows.append({"environment": f.stem.replace("env_sensitivity_", ""), "cell_bins_compared": len(m),
                     "max_abs_diff_vs_skmisc_0.5.2": float((m.fit - m.fit_o).abs().max()),
                     "historical_mismatch_bins_with_fit_le0": f"{int((h.fit <= 0).sum())} / {len(h)}"})
    return (C.html_table(pd.DataFrame(rows)) +
            "<p>All tested environments (Python 3.10/3.13, numpy 1.24/2.5, scikit-misc 0.1.4–0.5.3) give the same predictions to "
            "≤1e-13 and none gives a non-positive prediction at the 17 bins that are zero in the saved X. The historical cause of "
            "those zeros is therefore <b>unresolved</b> (a different execution path, guard, version or platform is possible; the "
            "run environment was not recorded). It affects 17 of 87.7 million cell-bins.</p>")


def main():
    p = C.default_paths()
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default=str(p["scratch"]))
    a = ap.parse_args()
    d = C.H5AD(p["h5ad"])
    cells = pd.read_csv(OUT / "per_cell.tsv.gz", sep="\t")
    summ = pd.read_json(OUT / "summary.json", typ="series")
    fit = np.load(Path(a.scratch) / "gc_fitted_float32.npy", mmap_mode="r")
    cor = np.load(Path(a.scratch) / "gc_corrected_prerint_float32.npy", mmap_mode="r")

    # ---- 1. mismatching cells: recover historical integers, list differing bins
    mism = []
    for r in cells.loc[~cells.exact_equal, "row"]:
        raw = d.rows("layers/raw", r, r + 1)[0]
        xo = d.rows("X", r, r + 1)[0]
        rec = np.rint(np.asarray(cor[r], np.float64)); rec[rec < 0] = 0
        ok = (rec > 0) & (xo > 0)
        ratio = np.median(xo[ok] / rec[ok])       # = 1e5 / T_hist if most bins agree
        hist_int = np.rint(xo / ratio)
        resid = np.max(np.abs(xo / ratio - hist_int))
        for b in np.flatnonzero(hist_int != rec):
            mism.append({"row": r, "cell": d.cells[r], "bin": b, "seq": d.seq[b], "start": d.start[b], "gc": d.gc[b],
                         "raw": raw[b], "fit_this_env": float(fit[r, b]), "rebuilt_int": rec[b], "historical_int": hist_int[b],
                         "historical_T": float(1e5 / ratio), "inversion_max_residual": float(resid)})
    mism = pd.DataFrame(mism)
    mism.to_csv(OUT / "mismatch_bins.tsv", sep="\t", index=False)

    # ---- 2. amplification in the saved (historical) X: per raw-positive cell-bin factor = X_int/raw
    T_rec = cells.set_index("row")["int_total"].to_numpy().astype(float)
    edges = np.array([0, 0.5, 1, 2, 3, 5, 10, 30, 100, 1e3, 1e9])
    hist = np.zeros(len(edges) - 1, np.int64)
    big = []      # cell-bins with factor > 30 in the saved X
    pb_share = np.zeros(d.n_bins)   # sum over cells of X from bins with factor > 10
    x_sum = np.zeros(d.n_bins)
    step = 200
    for lo in range(0, d.n_cells, step):
        hi = min(d.n_cells, lo + step)
        raw = d.rows("layers/raw", lo, hi); xo = d.rows("X", lo, hi)
        for r in range(hi - lo):
            gi = lo + r
            rp = raw[r] > 0
            T = T_rec[gi]
            if not cells.exact_equal.iloc[gi]:
                ok = (xo[r] > 0)
                rec = np.rint(np.asarray(cor[gi], np.float64)); rec[rec < 0] = 0
                both = ok & (rec > 0)
                T = 1e5 / np.median(xo[r][both] / rec[both])
            xi = xo[r] * T / 1e5
            fac = np.zeros(d.n_bins); fac[rp] = xi[rp] / raw[r][rp]
            hist += np.histogram(fac[rp], edges)[0]
            m = rp & (fac > 30)
            for b in np.flatnonzero(m):
                big.append({"row": gi, "cell": d.cells[gi], "bin": b, "seq": d.seq[b], "start": d.start[b], "gc": d.gc[b],
                            "raw": raw[r][b], "x_saved": xo[r][b], "factor": fac[b], "fit": float(fit[gi, b])})
            pb_share += np.where(rp & (fac > 10), xo[r], 0)
            x_sum += xo[r]
    big = pd.DataFrame(big)
    big.to_csv(OUT / "amplified_cellbins_factor_gt30.tsv.gz", sep="\t", index=False, compression="gzip")
    amp = pd.DataFrame({"factor_lo": edges[:-1], "factor_hi": edges[1:], "raw_positive_cellbins": hist})
    amp["fraction"] = amp.raw_positive_cellbins / amp.raw_positive_cellbins.sum()
    amp.to_csv(OUT / "amplification_histogram.tsv", sep="\t", index=False)
    share = np.divide(pb_share, x_sum, out=np.zeros_like(x_sum), where=x_sum > 0)
    pb = pd.read_csv(OUT / "per_bin.tsv.gz", sep="\t")
    pb["share_allcell_mean_from_factor_gt10"] = share
    pb.to_csv(OUT / "per_bin.tsv.gz", sep="\t", index=False, compression="gzip")

    # ---- 3. where do near-zero / negative fits occur (GC)?
    gcq = np.quantile(d.gc, [0.001, 0.01, 0.5, 0.99, 0.999])
    low_fit_gc = pb.loc[pb.n_cells_fit_le0 > 0, "gc"]
    res = {
        "mismatch_cells": int((~cells.exact_equal).sum()),
        "mismatch_bins": int(len(mism)),
        "mismatch_fit_this_env_range": [float(mism.fit_this_env.min()), float(mism.fit_this_env.max())] if len(mism) else None,
        "mismatch_pattern_counts": mism.assign(kind=np.where(mism.historical_int == 0, "historical_zero", "historical_nonzero"))
                                       .kind.value_counts().to_dict() if len(mism) else {},
        "inversion_max_residual": float(mism.inversion_max_residual.max()) if len(mism) else None,
        "amplification_histogram": amp.to_dict("records"),
        "saved_X_cellbins_factor_gt30": int(len(big)),
        "saved_X_cells_with_factor_gt30": int(big.cell.nunique()) if len(big) else 0,
        "saved_X_max_factor": float(big.factor.max()) if len(big) else None,
        "bins_where_factor_gt10_contributes_gt10pct_of_allcell_mean": int((share > 0.10).sum()),
        "bins_where_factor_gt10_contributes_gt25pct_of_allcell_mean": int((share > 0.25).sum()),
        "gc_quantiles_0.1_1_50_99_99.9pct": gcq.tolist(),
        "gc_of_bins_with_any_fit_le0": {"n_bins": int(len(low_fit_gc)), "min": float(low_fit_gc.min()), "max": float(low_fit_gc.max()),
                                        "frac_below_gc_q01": float((low_fit_gc < gcq[1]).mean())},
        "norm_factor_extremes": cells.nlargest(5, "norm_factor")[["cell", "raw_total", "int_total", "norm_factor", "n_fit_le0", "fit_min", "corr_max"]].to_dict("records"),
    }
    C.write_json(OUT / "analysis.json", res)

    # ---- figures
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.8))
    sub = np.random.default_rng(C.SEED).choice(d.n_cells, 40, replace=False)
    for r in sub:
        o = np.argsort(d.gc)
        ax[0].plot(d.gc[o], fit[r][o], lw=.6, alpha=.6)
    ax[0].axhline(0, color="k", lw=.8); ax[0].set_xlabel("bin GC"); ax[0].set_ylabel("LOESS prediction (raw-count scale)")
    ax[0].set_title("Per-cell LOESS fits (40 cells, seed 20261104)")
    ax[1].bar(range(len(amp)), amp.raw_positive_cellbins, color="#2f5d8a"); ax[1].set_yscale("log")
    ax[1].set_xticks(range(len(amp))); ax[1].set_xticklabels([f"{lo:g}-{hi:g}" for lo, hi in zip(edges[:-1], edges[1:])], rotation=60, fontsize=7)
    ax[1].set_title("Saved X: integer count / raw count, raw>0 cell-bins"); ax[1].set_ylabel("cell-bins")
    ax[2].scatter(pb.gc, pb.share_allcell_mean_from_factor_gt10, s=2, alpha=.4, color="#9b2c2c")
    ax[2].set_xlabel("bin GC"); ax[2].set_ylabel("share of all-cell mean X"); ax[2].set_title("Share of pseudobulk from cell-bins amplified >10x")
    b64 = C.fig_to_b64(fig); plt.close(fig)

    s = summ
    rp = s["raw_positive_to_zero"]; fit_s = s["fit"]; nf = s["normalization_factor"]; gcb = res["gc_of_bins_with_any_fit_le0"]
    hist_big = amp[amp.factor_lo >= 30].raw_positive_cellbins.sum()
    n_rawpos_total = int(amp.raw_positive_cellbins.sum())
    top_share = pb.loc[pb.share_allcell_mean_from_factor_gt10 > 0.25, ["bin", "seq", "start", "gc", "mean_raw", "mean_x_obs", "share_allcell_mean_from_factor_gt10"]].copy()
    top_share = top_share.sort_values("share_allcell_mean_from_factor_gt10", ascending=False).rename(columns={"share_allcell_mean_from_factor_gt10": "share_mean_X_from_gt10x"})
    extremes = pd.DataFrame(res["norm_factor_extremes"])[["cell", "raw_total", "int_total", "norm_factor", "n_fit_le0", "corr_max"]]
    mism_show = mism[["cell", "seq", "start", "gc", "raw", "fit_this_env", "rebuilt_int", "historical_int", "historical_T"]] if len(mism) else mism
    cells_mm = cells.loc[~cells.exact_equal]
    verdict = (f"The pinned pyEpi GC step reproduces the saved X bit-for-bit for {s['reconstruction']['cells_exact_equal']:,} of {s['cells']:,} cells; "
               f"the {res['mismatch_bins']} remaining cell-bins (in {res['mismatch_cells']} cells) are zero historically but not in any tested environment, so their cause is unresolved, "
               "and the step has no guard for near-zero or negative LOESS predictions.")
    findings = [
        f"<b>{s['reconstruction']['cells_exact_equal']:,} / {s['cells']:,}</b> cells ({100*s['reconstruction']['cells_exact_equal']/s['cells']:.1f}%) reproduce saved <code>X</code> exactly from <code>layers/raw</code> and <code>var.GC</code>; "
        f"{res['mismatch_cells']} cells differ in {res['mismatch_bins']} cell-bins of {s['cells']*s['bins']/1e6:.1f} million.",
        f"In all {res['mismatch_bins']} differing cell-bins the saved value is 0 ({res['mismatch_pattern_counts']}) while our LOESS prediction is small and positive "
        f"({res['mismatch_fit_this_env_range'][0]:.4f}–{res['mismatch_fit_this_env_range'][1]:.4f}), so the rebuild would amplify a single fragment up to {mism.rebuilt_int.max():,.0f} counts. "
        f"Recovered historical integers invert saved X to a residual of {res['inversion_max_residual']:.1e}.",
        f"No LOESS fit fell back to the raw vector ({fit_s['fallback_cells']}) and no non-finite float arose ({fit_s['cell_bins_nonfinite_float']}); "
        f"{fit_s['cell_bins_fit_le0']:,} cell-bins in {fit_s['cells_with_fit_le0']} cells have prediction &le; 0, all in {gcb['n_bins']} bins at GC {gcb['min']:.2f}–{gcb['max']:.2f} ({100*gcb['frac_below_gc_q01']:.0f}% below the 1st GC percentile).",
        f"{rp['total_observed_in_saved_X']:,} raw-positive cell-bins ({100*rp['fraction_of_raw_positive']:.2f}%) are 0 in saved X: {rp['rounding_0_to_0.5']:,} by rounding a value in [0, 0.5) and {rp['negative_clip']:,} by negative clipping in the rebuild; "
        f"{rp['raw_zero_but_X_positive_observed']} raw-zero cell-bins become positive.",
        f"{res['saved_X_cellbins_factor_gt30']:,} raw-positive cell-bins in {res['saved_X_cells_with_factor_gt30']} cells are multiplied by more than 30 in the saved X (maximum {res['saved_X_max_factor']:.0f}x); "
        f"in {res['bins_where_factor_gt10_contributes_gt25pct_of_allcell_mean']} bins, cell-bins amplified more than 10x make up over 25% of the all-cell mean X.",
        f"The most extreme per-cell normalization factor is {nf['max']:.2f} (median {nf['median']:.2f}, 5th percentile {nf['q05']:.2f}); the extreme cell has {res['norm_factor_extremes'][0]['n_fit_le0']} bins with prediction &le; 0.",
        "Five scikit-misc/numpy/Python environments give identical predictions and none yields a non-positive prediction at the differing bins; the historical run environment was not recorded.",
    ]
    sections = [
        ("Question and why it matters",
         "<p><b>Question.</b> Does the pinned pyEpi source (commit e8305ea, <code>pyEpiAneufinder.py</code> L252&ndash;L271 with <code>get_loess_smoothed</code> L155&ndash;L164) "
         "reproduce the saved <code>X</code> from <code>layers/raw</code> and <code>var.GC</code>, and what happens to the floating values before integer casting?</p>"
         "<p><b>Why it matters.</b> The review's recommendation says the next change belongs in evidence integration and the CN reference, not in GC/rounding tuning. That claim only holds if the saved matrix is "
         "the documented formula applied to the retained raw counts (so nothing is hidden in an unrecorded step) and if the GC step's artefacts are quantified rather than assumed. "
         "This experiment establishes provenance of <code>X</code> (preprocessing_audit stages 6&ndash;11) and measures the size of the GC step's effects; whether those effects create called segments is tested in experiments 04&ndash;05.</p>"),
        ("Inputs",
         f"<ul><li><b>Matrix.</b> <code>pyepi_results/count_matrix.h5ad</code>: {s['cells']:,} cells &times; {s['bins']:,} bins (100 kb windows after QC). <code>layers/raw</code> = post-QC, pre-GC integer fragment counts; <code>X</code> = saved processed matrix (float64 CSR, row sums 1e5). "
         f"Stored non-zeros: raw {s['nnz']['layers/raw']:,}, X {s['nnz']['X']:,} (difference {s['nnz']['raw_minus_X']:,}); explicit stored zeros: X {s['explicit_zeros_stored']['X']}, raw {s['explicit_zeros_stored']['layers/raw']}.</li>"
         "<li><b>Bin GC</b> from <code>var.GC</code>. Units: fragments per 100 kb bin per cell (raw), counts per 1e5 (X).</li>"
         f"<li><b>Experiment-01 outputs read by this report:</b> <code>per_cell.tsv.gz</code>, <code>per_bin.tsv.gz</code>, <code>summary.json</code>, scratch float32 arrays <code>gc_fitted_float32.npy</code> and <code>gc_corrected_prerint_float32.npy</code> (cells &times; bins), and <code>env_sensitivity_*.tsv</code>.</li>"
         f"<li><b>Not available.</b> The historical run environment and <code>parameter_configuration.txt</code>; pre-QC counts and removed cells (so stages 3&ndash;5 cannot be reconstructed); fragment files.</li>"
         f"<li><b>Hash.</b> h5ad sha256 <code>{s['inputs']['h5ad_sha256'][:16]}&hellip;</code></li></ul>"),
        ("Method",
         "<ol><li>Load <code>layers/raw</code> and <code>var.GC</code>. For each cell fit <code>skmisc.loess(GC, counts, span = 0.75)</code> with all other settings default (pinned source), giving fitted value <code>fit</code>.</li>"
         "<li>Apply the pinned formula, rounding, clipping and normalization (stages 6&ndash;11 of preprocessing_audit.md).</li>"
         "<li>Compare to saved <code>X</code> cell by cell. Record float intermediates (before <code>rint</code>) in scratch so negative or non-finite values are visible.</li>"
         "<li><b>This script (analyze.py):</b> for each non-matching cell recover the historical integer vector as <code>X*T/1e5</code> with <code>T</code> chosen from bins that agree, list the bins that differ, and compute the amplification factor of every raw-positive cell-bin in the saved matrix.</li></ol>"
         "<pre>correction = mean(raw_cell) / fit(GC)\ncorrected  = raw * correction\ninteger    = rint(corrected)          # astype(int)\ninteger    = max(integer, 0)          # clip negatives\nX          = integer / sum(integer) * 1e5\n\nfactor(cell,bin) = X_int / raw   for raw &gt; 0     # amplification in the SAVED matrix\nT_hist = 1e5 / median(X / rebuilt_int)           # over bins where both are &gt; 0</pre>"
         "<p>Deterministic (no resampling). The only random choice is the 40 cells drawn for the left panel of the figure, with <code>np.random.default_rng(20261104)</code>. "
         "Environment sensitivity (<code>env_sensitivity.py</code>) re-runs only the LOESS fit at the 17 differing cell-bins and at comparison bins in other package versions.</p>"),
        ("Results", "", 2),
        ("Reconstruction of saved X", 
         C.html_table(pd.DataFrame([
             {"quantity": "cells exactly equal", "value": f"{s['reconstruction']['cells_exact_equal']:,} / {s['cells']:,}"},
             {"quantity": "cells with any bin |diff| > 1e-9", "value": s["reconstruction"]["cells_with_any_bin_absdiff_gt_1e-9"]},
             {"quantity": "differing cell-bins (after recovering historical row total)", "value": res["mismatch_bins"]},
             {"quantity": "LOESS fits that fell back to input", "value": fit_s["fallback_cells"]},
             {"quantity": "non-finite float values before casting", "value": fit_s["cell_bins_nonfinite_float"]},
             {"quantity": "historical value at the differing bins", "value": str(res["mismatch_pattern_counts"])},
             {"quantity": "max residual of historical-integer inversion", "value": f"{res['inversion_max_residual']:.2g}"},
         ])) +
         C.caption("Each row is one reconstruction check. Exact equality of 99.7% of cells means the pinned formula is what produced X. The 17 differing bins are the only evidence of an unrecorded step.") +
         "<p>The raw summary reports a large maximum absolute difference "
         f"({s['reconstruction']['max_abs_diff_overall']:,.0f}) because the unrecovered row total shifts every bin of the 11 affected cells; once the historical row total is recovered, every other bin matches (see Method).</p>", 3),
        ("The 17 historical zeros", 
         C.html_table(mism_show) + C.caption(
             "One row per differing cell-bin. <code>raw</code> is the retained fragment count, <code>fit_this_env</code> our LOESS prediction, <code>rebuilt_int</code> what the pinned formula gives here, <code>historical_int</code> what saved X implies (always 0). "
             "Look for: tiny positive fits that would give rebuilt counts in the hundreds to thousands, all at low-GC boundary bins.") +
         f"<p>Affected cells: {len(cells_mm)}. Full list: <code>mismatch_bins.tsv</code>.</p>", 3),
        ("Environment sensitivity", env_section() +
         C.caption("Each row is one Python/numpy/scikit-misc combination. Differences of about 1e-13 and zero non-positive fits at the 17 bins mean the library version does not explain the historical zeros. Full comparison: <code>env_sensitivity_comparison.tsv</code>."), 3),
        ("Non-positive fits and rounding losses",
         f"<p>{fit_s['cell_bins_fit_le0']:,} cell-bins in {fit_s['cells_with_fit_le0']} cells have a LOESS prediction &le; 0. All fall in {gcb['n_bins']} bins at GC {gcb['min']:.3f}&ndash;{gcb['max']:.3f} (GC quantiles 0.1/1/50/99/99.9%: "
         f"{', '.join(f'{q:.3f}' for q in res['gc_quantiles_0.1_1_50_99_99.9pct'])}). Of {rp['total_observed_in_saved_X']:,} raw-positive cell-bins that are 0 in saved X "
         f"({100*rp['fraction_of_raw_positive']:.2f}% of raw-positive cell-bins), the rebuild attributes {rp['rounding_0_to_0.5']:,} to <code>np.rint</code> of a value in [0, 0.5) and {rp['negative_clip']:,} to clipping a negative value "
         f"(total reconstructed {rp['total_recon']:,}; the {rp['total_observed_in_saved_X']-rp['total_recon']} extra observed are the historical zeros above). {rp['raw_zero_but_X_positive_observed']} raw-zero cell-bins become positive.</p>" +
         C.caption("Negative predictions only occur at the lowest-GC edge of the LOESS fit, where the local fit extrapolates. Rounding losses are coverage-dependent and analysed in experiment 02."), 3),
        ("Amplification in the saved X",
         C.html_table(amp) + C.caption("Histogram of factor = integer X count / raw count over raw-positive cell-bins in the saved matrix (columns <code>factor_lo</code>&ndash;<code>factor_hi</code>; <code>fraction</code> of all raw-positive cell-bins). "
                                      "Most bins sit between 0.5 and 3; look at the thin tail above 30.") +
         f"<p>{res['saved_X_cellbins_factor_gt30']:,} cell-bins ({100*res['saved_X_cellbins_factor_gt30']/n_rawpos_total:.4f}% of raw-positive cell-bins) in {res['saved_X_cells_with_factor_gt30']} cells exceed 30x (maximum {res['saved_X_max_factor']:.0f}x). "
         f"{res['bins_where_factor_gt10_contributes_gt10pct_of_allcell_mean']} bins have &gt;10% and {res['bins_where_factor_gt10_contributes_gt25pct_of_allcell_mean']} bins &gt;25% of their all-cell mean X contributed by cell-bins amplified more than 10x.</p>" +
         (C.html_table(top_share) + C.caption("Bins where more than 25% of the all-cell mean X comes from cell-bins amplified more than 10x. <code>mean_raw</code> is the raw mean per cell, <code>mean_x_obs</code> the saved-X mean; a large ratio of the latter to the genome mean (about 3.9 per the preprocessing audit) marks an inflated bin.") if len(top_share) else "") +
         C.img(b64, "LOESS fits and amplification") +
         C.caption("Left: per-cell LOESS fits against bin GC for 40 random cells; the black line is zero, and fits crossing it at low GC give negative or near-zero denominators. Middle: log-scale histogram of the amplification factor. "
                   "Right: per-bin share of the all-cell mean X coming from cell-bins amplified more than 10x, against GC; look for a cluster of bins at the low-GC edge."), 3),
        ("Per-cell normalization", C.html_table(extremes) + C.caption("The five cells with the largest <code>norm_factor</code> = integer total / 1e5-scale divisor. A factor well above 1 means the GC step inflated the cell total (mostly through a few amplified bins), so row normalization then shrinks every other bin of that cell. Columns: <code>raw_total</code> fragments, <code>int_total</code> corrected integer total, <code>n_fit_le0</code> bins with prediction &le; 0, <code>corr_max</code> largest float correction result."), 3),
        ("Per-region detail",
         "<p>This experiment is genome-wide and does not analyse the six saved windows. Per the preprocessing audit, none of the amplified low-GC bins listed above lies in a saved window, but each can inflate a segment mean that contains it; "
         "the saved-window contrasts are examined in experiments 04 and 07.</p>"),
        ("Sensitivity",
         "<p>The only sensitivity run here is across software environments (table above). Sensitivity of results to LOESS span or alternative GC correction was not explored: the review deliberately did not pursue GC internals further.</p>"),
        ("Limits and alternative explanations",
         "<ul><li>The 17 historical zeros could come from a different execution path, an undocumented guard, a different platform or library build; this analysis cannot distinguish them. Our rebuild must not replace the historical matrix.</li>"
         "<li>Only post-QC data are available: the cell and bin filters (stages 3&ndash;5) are taken as given, and cells removed by them cannot be audited.</li>"
         "<li>Exact reproduction shows the formula was applied, not that the formula is appropriate for CN inference.</li>"
         "<li>Amplification factors in the saved X are inferred from X and raw, with the historical row total recovered by a median ratio; for the 11 non-matching cells this adds the inversion assumption (residual shown above).</li>"
         "<li>Not observed but plausible: an exactly zero prediction would give inf/NaN that integer casting wraps to a large negative number and clipping turns into 0; a fit exception would return the raw vector.</li></ul>"),
    ]
    glossary = {
        "layers/raw": "Post-QC, pre-GC-correction fragment counts per cell and 100 kb bin.",
        "X": "Saved processed matrix: GC-corrected, rounded, clipped, each row normalized to sum 1e5.",
        "fit / fit_this_env": "Per-cell LOESS prediction of the raw count as a function of bin GC (span 0.75), on the raw-count scale.",
        "correction": "mean(raw of cell) / fit; the multiplier applied to a bin. Unbounded when fit is near 0, negative when fit is negative.",
        "rebuilt_int / historical_int": "Integer count after rounding and clipping in our rebuild / recovered from saved X.",
        "historical_T": "Recovered integer row total of the historical run for that cell (1e5 / median ratio of X to rebuilt integers).",
        "factor": "Integer X count divided by raw count for a raw-positive cell-bin: how much the GC step scaled it.",
        "factor_lo / factor_hi": "Bin edges of the amplification histogram.",
        "raw_positive_cellbins / fraction": "Number and share of raw-positive cell-bins in the factor bin.",
        "share_mean_X_from_gt10x": "Share of a bin's all-cell mean X contributed by cell-bins whose factor exceeds 10.",
        "norm_factor": "Integer total of a cell divided by 1e5 (as reported by per_cell.tsv.gz); values far from the typical range indicate a distorted cell total.",
        "int_total / raw_total": "Cell sum of integer corrected counts / of raw counts.",
        "n_fit_le0": "Number of bins in the cell where the LOESS prediction is zero or negative.",
        "corr_max": "Largest corrected float value in the cell before rounding.",
        "exact_equal": "True if the rebuilt X equals the saved X to within 1e-9 in every bin.",
        "mean_raw / mean_x_obs": "Mean over cells of raw counts / of saved X for a bin.",
        "max_abs_diff_vs_skmisc_0.5.2": "Largest absolute difference of LOESS predictions between an environment and the reference environment.",
    }
    C.write_html_report(
        OUT / "report.html", "X reconstruction audit", "Experiment 01 &middot; pinned pyEpi GC correction rebuilt on retained raw counts",
        sections, {"command": "python -I experiments/01_reconstruct_x/run.py && python -I experiments/01_reconstruct_x/analyze.py",
                   "summary": "summary.json", "analysis": "analysis.json", "environment": C.environment()},
        key_findings=findings, verdict=verdict, glossary=glossary,
        supports=["Saved X is the pinned GC-correction formula applied to the retained raw layer and GC (bit-identical for most cells).",
                  "The GC step has no guard for near-zero or negative predictions; its artefacts (zeroing, large amplification) are measured at the low-GC edge.",
                  "Rounding to integers removes a small, coverage-dependent fraction of raw-positive cell-bins (see experiment 02)."],
        not_supports=["Why 17 cell-bins are 0 historically (no tested environment reproduces it).",
                      "Whether any of the GC artefacts create or hide a called segment (tested in experiments 04&ndash;05).",
                      "That GC correction is appropriate for CN inference, or that our rebuild may replace the historical matrix."],
        decision_impact="Provenance of X is sound enough to treat the preprocessing as a reproducible measurement layer, which supports the recommendation's point that the next change belongs in evidence integration and the CN reference rather than in GC tuning. "
                        "The recommendation's stage 0 hygiene step (persist run parameters, GC fit diagnostics and QC masks) follows directly from the unexplained zeros and the unguarded denominator.",
        next_steps=["Record the environment, parameters and GC fit diagnostics (predictions &le; 0, maximum correction) with every future run.",
                    "If the historical environment can be recovered, re-run the 17 cell-bins there.",
                    "Test whether the amplified low-GC bins change any saved-window segment mean (experiments 04/05 cover the windows; a segment-level check across the genome is open)."])
    print({k: v for k, v in res.items() if k not in ("amplification_histogram", "norm_factor_extremes")})


if __name__ == "__main__":
    main()
