"""Experiment 01: reconstruct AnnData X from layers/raw + var.GC with the pinned pyEpi source.

Reproduces pyEpiAneufinder.py L252-L271 (commit e8305ea) step by step on the retained
3,460 x 25,337 matrix and keeps every floating intermediate *before* np.rint/astype(int)
so that non-finite or negative values cannot be hidden by casting and clipping.

The LOESS fit is the pinned `render_fragments.get_loess_smoothed` itself (imported from
file, unmodified). A fit failure is detected by object identity (the helper returns its
input on exception). The saved X is never modified.

Deterministic (no random numbers). Writes small tables here and large float32
intermediates to --scratch (not committed).
"""
import argparse
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common as C  # noqa: E402

OUT = Path(__file__).resolve().parent
CHUNK = 40


def work(args):
    lo, hi, h5ad, pkg = args
    import contextlib
    import io
    rf = C.load_pyepi_module("render_fragments", pkg)
    d = C.H5AD(h5ad)
    raw = d.rows("layers/raw", lo, hi)
    xo = d.rows("X", lo, hi)
    gc = d.gc
    nb = d.n_bins
    bins = {k: np.zeros(nb, np.int64) for k in
            ["raw_pos", "raw_pos_X0_obs", "round0", "negclip", "fit_le0", "nonfinite", "neg_float"]}
    bins_sum = {k: np.zeros(nb) for k in ["raw", "x_obs", "corr_float_clipped"]}
    fitted = np.zeros((hi - lo, nb), np.float32)
    corrf = np.zeros((hi - lo, nb), np.float32)
    rows = []
    for r in range(hi - lo):
        counts = raw[r].copy()
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            fit = rf.get_loess_smoothed(counts, gc)
        fallback = fit is counts
        fit = np.asarray(fit, dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            correction = counts.mean() / fit          # pinned L257
            loess_norm_row = counts * correction      # pinned L258
        with np.errstate(invalid="ignore"):
            as_int = np.rint(loess_norm_row).astype(int)  # pinned L260 (may wrap non-finite)
        clipped = as_int.copy()
        clipped[clipped < 0] = 0                      # pinned L267
        tot = clipped.sum()
        xr = clipped / (tot / 1e5)                    # pinned L269
        fin = np.isfinite(loess_norm_row)
        rp = counts > 0
        round0 = rp & fin & (loess_norm_row >= 0) & (loess_norm_row < 0.5)
        negclip = rp & (as_int < 0)
        diff = np.abs(xr - xo[r])
        rows.append({
            "row": lo + r, "cell": d.cells[lo + r], "raw_total": counts.sum(), "raw_nnz": int(rp.sum()),
            "raw_mean": counts.mean(), "fit_fallback": bool(fallback), "fit_stderr": err.getvalue().strip()[:200],
            "fit_min": fit.min(), "fit_max": fit.max(), "n_fit_le0": int((fit <= 0).sum()),
            "n_fit_lt_10pct_mean": int((fit < 0.1 * counts.mean()).sum()),
            "corr_max": float(np.nanmax(np.where(np.isfinite(correction), correction, np.nan))),
            "n_nonfinite_float": int((~fin).sum()), "n_negative_float": int((fin & (loess_norm_row < 0)).sum()),
            "n_rawpos_negative_float": int((rp & fin & (loess_norm_row < 0)).sum()),
            "n_rawpos_round_to_0": int(round0.sum()), "n_rawpos_negclip_to_0": int(negclip.sum()),
            "n_rawpos_X0_recon": int((rp & (clipped == 0)).sum()), "n_rawpos_X0_obs": int((rp & (xo[r] == 0)).sum()),
            "n_raw0_Xpos_obs": int(((~rp) & (xo[r] > 0)).sum()),
            "int_total": int(tot), "norm_factor": tot / 1e5, "raw_over_int_total": counts.sum() / tot,
            "x_obs_rowsum": xo[r].sum(), "max_abs_diff": diff.max(), "n_bins_absdiff_gt_1e-9": int((diff > 1e-9).sum()),
            "exact_equal": bool(np.array_equal(xr, xo[r])),
            "x_obs_times_factor_max_frac": float(np.max(np.abs(xo[r] * (tot / 1e5) - np.rint(xo[r] * (tot / 1e5))))),
        })
        bins["raw_pos"] += rp
        bins["raw_pos_X0_obs"] += rp & (xo[r] == 0)
        bins["round0"] += round0
        bins["negclip"] += negclip
        bins["fit_le0"] += fit <= 0
        bins["nonfinite"] += ~fin
        bins["neg_float"] += fin & (loess_norm_row < 0)
        bins_sum["raw"] += counts
        bins_sum["x_obs"] += xo[r]
        bins_sum["corr_float_clipped"] += np.where(fin, np.maximum(loess_norm_row, 0), 0)
        fitted[r] = fit
        corrf[r] = np.where(fin, loess_norm_row, np.nan)
    return lo, hi, rows, bins, bins_sum, fitted, corrf


def main():
    p = C.default_paths()
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--h5ad", default=str(p["h5ad"]))
    ap.add_argument("--package-dir", default=str(p["package"]))
    ap.add_argument("--scratch", default=str(p["scratch"]))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2)))
    ap.add_argument("--limit", type=int, default=0, help="bounded check: first N cells only")
    a = ap.parse_args()
    t0 = time.time()
    d = C.H5AD(a.h5ad)
    n = a.limit or d.n_cells
    scratch = Path(a.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    fit_mm = np.lib.format.open_memmap(scratch / "gc_fitted_float32.npy", "w+", np.float32, (n, d.n_bins))
    cor_mm = np.lib.format.open_memmap(scratch / "gc_corrected_prerint_float32.npy", "w+", np.float32, (n, d.n_bins))
    tasks = [(lo, min(n, lo + CHUNK), a.h5ad, a.package_dir) for lo in range(0, n, CHUNK)]
    rows, bins, bsum = [], None, None
    with get_context("fork").Pool(a.workers) as pool:
        for lo, hi, rr, b, bs, fit, cor in pool.imap_unordered(work, tasks):
            rows += rr
            bins = b if bins is None else {k: bins[k] + b[k] for k in b}
            bsum = bs if bsum is None else {k: bsum[k] + bs[k] for k in bs}
            fit_mm[lo:hi] = fit
            cor_mm[lo:hi] = cor
    fit_mm.flush(); cor_mm.flush()
    cells = pd.DataFrame(rows).sort_values("row")
    cells.to_csv(OUT / "per_cell.tsv.gz", sep="\t", index=False, compression="gzip")
    pb = pd.DataFrame({"bin": np.arange(d.n_bins), "seq": d.seq, "start": d.start, "end": d.end, "gc": d.gc,
                       **{"n_cells_" + k: v for k, v in bins.items()},
                       "mean_raw": bsum["raw"] / n, "mean_x_obs": bsum["x_obs"] / n,
                       "mean_corr_float_clipped": bsum["corr_float_clipped"] / n})
    pb.to_csv(OUT / "per_bin.tsv.gz", sep="\t", index=False, compression="gzip")
    nnz_raw, nnz_x = int(d.row_nnz("layers/raw")[:n].sum()), int(d.row_nnz("X")[:n].sum())
    summary = {
        "experiment": "01_reconstruct_x", "deterministic": True, "cells": n, "bins": d.n_bins,
        "seconds": round(time.time() - t0, 1),
        "pinned_source": {"render_fragments.get_loess_smoothed": "L155-L164", "pyEpiAneufinder.py": "L252-L271",
                          "commit": "e8305ea2162634d07e92a3bb401df672e48f217b",
                          "loess_settings": "skmisc.loess(gc, counts, span=0.75), all other skmisc defaults"},
        "explicit_zeros_stored": {"X": d.explicit_zeros("X"), "layers/raw": d.explicit_zeros("layers/raw")},
        "nnz": {"layers/raw": nnz_raw, "X": nnz_x, "raw_minus_X": nnz_raw - nnz_x},
        "reconstruction": {
            "cells_exact_equal": int(cells.exact_equal.sum()),
            "cells_with_any_bin_absdiff_gt_1e-9": int((cells["n_bins_absdiff_gt_1e-9"] > 0).sum()),
            "max_abs_diff_overall": float(cells.max_abs_diff.max()),
            "median_max_abs_diff": float(cells.max_abs_diff.median()),
            "total_bins_absdiff_gt_1e-9": int(cells["n_bins_absdiff_gt_1e-9"].sum()),
            "x_obs_rowsum_min_max": [float(cells.x_obs_rowsum.min()), float(cells.x_obs_rowsum.max())],
        },
        "fit": {"fallback_cells": int(cells.fit_fallback.sum()),
                "cells_with_fit_le0": int((cells.n_fit_le0 > 0).sum()),
                "cell_bins_fit_le0": int(cells.n_fit_le0.sum()),
                "cells_with_nonfinite_float": int((cells.n_nonfinite_float > 0).sum()),
                "cell_bins_nonfinite_float": int(cells.n_nonfinite_float.sum()),
                "cell_bins_negative_float": int(cells.n_negative_float.sum())},
        "raw_positive_to_zero": {
            "rounding_0_to_0.5": int(cells.n_rawpos_round_to_0.sum()),
            "negative_clip": int(cells.n_rawpos_negclip_to_0.sum()),
            "total_recon": int(cells.n_rawpos_X0_recon.sum()),
            "total_observed_in_saved_X": int(cells.n_rawpos_X0_obs.sum()),
            "raw_zero_but_X_positive_observed": int(cells.n_raw0_Xpos_obs.sum()),
            "fraction_of_raw_positive": float(cells.n_rawpos_X0_obs.sum() / cells.raw_nnz.sum()),
        },
        "normalization_factor": {q: float(cells.norm_factor.quantile(v)) for q, v in
                                 [("min", 0), ("q05", .05), ("median", .5), ("q95", .95), ("max", 1)]},
        "environment": C.environment(),
        "inputs": {"h5ad": a.h5ad, "h5ad_sha256": C.sha256(a.h5ad)},
        "scratch_outputs": [str(scratch / "gc_fitted_float32.npy"), str(scratch / "gc_corrected_prerint_float32.npy")],
    }
    C.write_json(OUT / "summary.json", summary)
    print(f"done {n} cells in {summary['seconds']} s")


if __name__ == "__main__":
    main()
