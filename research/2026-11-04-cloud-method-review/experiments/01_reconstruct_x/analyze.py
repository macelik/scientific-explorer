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
    frag = [
        ("Question", "<p>Does the pinned pyEpi source (e8305ea, <code>pyEpiAneufinder.py</code> L252–L271 with "
         "<code>get_loess_smoothed</code> L155–L164) reproduce the saved <code>X</code> from <code>layers/raw</code> and "
         "<code>var.GC</code>? What happens to the floating values before integer casting?</p>"),
        ("Result", f"<div class='card'><b>{s['reconstruction']['cells_exact_equal']} / {s['cells']}</b> cells reproduce bit-for-bit. "
         f"{res['mismatch_cells']} cells differ, in {res['mismatch_bins']} bins in total; in every differing bin our rebuild's LOESS "
         f"prediction is near zero (range {res['mismatch_fit_this_env_range']}). Recovering the historical integers from saved X "
         f"(residual ≤ {res['inversion_max_residual']:.2g}) shows the historical value was "
         f"{res['mismatch_pattern_counts']} at those bins. This is consistent with platform/library-dependent sign or magnitude "
         "of a near-zero LOESS denominator, not with a different formula. The historical run environment is not recorded, so this "
         "remains a provenance finding.</div>"
         f"<p>No fit fell back to the raw vector; no non-finite float occurred in this environment. "
         f"{s['fit']['cell_bins_fit_le0']} cell-bins in {s['fit']['cells_with_fit_le0']} cells have LOESS prediction ≤ 0 "
         f"(GC range of affected bins {res['gc_of_bins_with_any_fit_le0']['min']:.3f}–{res['gc_of_bins_with_any_fit_le0']['max']:.3f}; "
         f"{100*res['gc_of_bins_with_any_fit_le0']['frac_below_gc_q01']:.0f}% below the 1st GC percentile): the LOESS boundary.</p>"),
        ("Raw-positive values that become zero", f"<p>{s['raw_positive_to_zero']['total_observed_in_saved_X']:,} raw-positive cell-bins are 0 in saved X "
         f"({100*s['raw_positive_to_zero']['fraction_of_raw_positive']:.2f}% of raw-positive cell-bins). Rebuild attribution: "
         f"{s['raw_positive_to_zero']['rounding_0_to_0.5']:,} by <code>np.rint</code> of a value in [0, 0.5) and "
         f"{s['raw_positive_to_zero']['negative_clip']:,} by negative-to-zero clipping (negative LOESS prediction). "
         f"No raw-zero cell-bin becomes positive.</p>"),
        ("Amplification in the saved X", C.html_table(amp) +
         f"<p>{res['saved_X_cellbins_factor_gt30']:,} raw-positive cell-bins in {res['saved_X_cells_with_factor_gt30']} cells were multiplied "
         f"by more than 30 (max {res['saved_X_max_factor']:.0f}×). In {res['bins_where_factor_gt10_contributes_gt10pct_of_allcell_mean']} bins, "
         f"cell-bins amplified more than 10× make up over 10% of the all-cell mean X.</p>" + C.img(b64, "LOESS fits and amplification")),
        ("Mismatching bins", C.html_table(mism) if len(mism) else "<p>none</p>"),
        ("Environment sensitivity (env_sensitivity.py)", env_section()),
        ("Interpretation", "<ul><li>The saved X is the pinned formula applied to the saved raw layer and GC (bit-identical for 99.7% of cells); "
         "the remaining differences sit exactly where <code>correction = mean/fit</code> has a near-zero denominator.</li>"
         "<li>The formula has no guard for small or negative predictions: a negative prediction zeroes a raw-positive bin, a tiny "
         "positive one multiplies it by hundreds to thousands. Both are measurement artefacts of the GC step, observed here at the GC boundary.</li>"
         "<li>Rounding to integers removes about 0.9% of raw-positive cell-bins, concentrated in cells and bins with low expected counts; see experiment 02.</li>"
         "<li>This check does not show whether these artefacts create any called segment; that is tested in experiments 04–05.</li></ul>"),
    ]
    C.write_html_report(OUT / "report.html", "X reconstruction audit", "Experiment 01 · pinned pyEpi GC correction rebuilt on retained raw counts",
                        frag, {"command": "python -I experiments/01_reconstruct_x/run.py && python -I experiments/01_reconstruct_x/analyze.py",
                               "summary": "summary.json", "analysis": "analysis.json", "environment": C.environment()})
    print({k: v for k, v in res.items() if k not in ("amplification_histogram", "norm_factor_extremes")})


if __name__ == "__main__":
    main()
