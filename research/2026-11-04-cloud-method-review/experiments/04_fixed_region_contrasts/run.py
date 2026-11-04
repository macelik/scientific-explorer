"""Experiment 04: fixed-region raw-count contrasts, per-cell homogeneity and dosage compatibility.

Region families, defined by rule before computing any statistic below:
  saved      - the six saved chr7/chr9 windows and the historical chr9 70-100 Mb comparison,
               with their frozen 50-retained-bin flanks (research/2026-10-29-low-coverage/shared).
  called     - every segment of the `all` source in the saved genome-pool karyogram with
               n_bins <= 100 and cn_call != 1 (46 loss/weak loss, 21 gain/weak gain). Flanks: up to
               50 retained bins on each side within the chromosome. Calls are not truth.
  literature - published MCF10-A loci (sources.json): 9p21 hemizygous zone chr9:20-24 Mb
               (approximate; Cowell 2005 ~4 Mb), MYC chr8:127.0-128.5 Mb, 19q13-qter chr19>=35 Mb,
               5q13-qter chr5>=70 Mb. Compared with the genome-wide density only.
  random     - 300 random windows of 25 consecutive retained bins with no internal gap, seed tuple
               (20261104, 4), numpy PCG64; same flank rule. Not neutral truth: they can contain events.

Statistics (all 3,460 cells, raw counts = layers/raw; processed = X):
  rr_flank   - ratio estimator (sum raw per target bin)/(sum raw per flank bin); exposure cancels.
  rr_genome  - target density over genome-wide density.
  rr by raw-coverage tercile and WNN cluster; coverage_log_slope = log(rr_low/rr_high).
  phi        - Pearson dispersion of y_target | y_target + y_flank around a common binomial
               proportion (phi_cluster: proportion per WNN cluster). ~1 = homogeneous Poisson-like;
               a subclonal dosage change inflates phi relative to the null windows.
  zero test  - cells with y_target = 0 versus the binomial expectation sum_i (1-pi)^n_i, and the
               fraction of cells whose expected target count exceeds 5 yet observe 0 (an upper bound
               on the fraction of cells that could carry a near-complete loss).
  dosage     - rr_flank < 0.5 cannot arise from a one-copy loss of a diploid baseline with unchanged
               accessibility/mappability alone.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

FLANK = 50
LIT = [("lit_9p21_hemizygous_zone_approx", "chr9", 20.0e6, 24.0e6, "loss"),
       ("lit_MYC_8q24", "chr8", 127.0e6, 128.5e6, "gain"),
       ("lit_19q13_qter", "chr19", 35.0e6, 1e9, "gain"),
       ("lit_5q13_qter", "chr5", 70.0e6, 1e9, "gain")]


def flanks(chrom_idx, lo, hi, n=FLANK):
    """chrom_idx: absolute bins of the chromosome (sorted); lo/hi: positions in it of the target [lo, hi)."""
    return chrom_idx[max(0, lo - n):lo], chrom_idx[hi:hi + n]


def build_regions(d, p):
    regs = []
    rb = C.shared_region_bins()
    for rid in rb.region_id.unique():
        if not rid.endswith("__target"):
            continue
        base = rid[: -len("__target")]
        t = rb.loc[rb.region_id == rid, "bin_index"].to_numpy()
        l = rb.loc[rb.region_id == base + "__left_flank", "bin_index"].to_numpy()
        r = rb.loc[rb.region_id == base + "__right_flank", "bin_index"].to_numpy()
        regs.append({"family": "saved", "name": base, "direction": "loss" if "70_100" not in base else "comparison",
                     "target": t, "flank": np.r_[l, r]})
    seg = pd.read_csv(p["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    seg = seg[(seg.group == "all") & (seg.n_bins <= 100) & (seg.cn_call != 1)]
    for _, s in seg.iterrows():
        ci = np.flatnonzero(d.seq == s.chromosome)
        lo, hi = int(np.searchsorted(ci, s.bin_start)), int(np.searchsorted(ci, s.bin_end)) + 1
        l, r = flanks(ci, lo, hi)
        regs.append({"family": "called", "name": f"{s.chromosome}:{d.start[s.bin_start]/1e6:.1f}-{d.end[s.bin_end]/1e6:.1f}_cn{s.cn_call}",
                     "direction": "loss" if s.cn_call < 1 else "gain", "cn_call": s.cn_call,
                     "target": np.arange(s.bin_start, s.bin_end + 1), "flank": np.r_[l, r]})
    for name, ch, a, b, dirn in LIT:
        t = np.flatnonzero((d.seq == ch) & (d.start >= a) & (d.end <= b))
        regs.append({"family": "literature", "name": name, "direction": dirn, "target": t, "flank": np.array([], int)})
    rng = np.random.default_rng((C.SEED, 4))
    ok_start = []
    for ch in dict.fromkeys(d.seq):
        ci = np.flatnonzero(d.seq == ch)
        for k in range(FLANK, len(ci) - 25 - FLANK):
            w = ci[k:k + 25]
            if d.start[w[-1]] - d.start[w[0]] == 24 * 100_000:
                ok_start.append((ch, k))
    pick = rng.choice(len(ok_start), 300, replace=False)
    for j in pick:
        ch, k = ok_start[j]
        ci = np.flatnonzero(d.seq == ch)
        l, r = flanks(ci, k, k + 25)
        regs.append({"family": "random", "name": f"random_{ch}:{d.start[ci[k]]/1e6:.1f}", "direction": "random",
                     "target": ci[k:k + 25], "flank": np.r_[l, r]})
    return regs


def region_sums(d, regs, key, step=250):
    cols, rows = [], []
    for j, r in enumerate(regs):
        for part, off in (("target", 0), ("flank", 1)):
            b = r[part]
            rows += list(b); cols += [2 * j + off] * len(b)
    M = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(d.n_bins, 2 * len(regs)))
    out = np.zeros((d.n_cells, 2 * len(regs)))
    tot = np.zeros(d.n_cells)
    for lo in range(0, d.n_cells, step):
        hi = min(d.n_cells, lo + step)
        blk = d.rows(key, lo, hi)
        out[lo:hi] = (sparse.csr_matrix(blk) @ M).toarray()
        tot[lo:hi] = blk.sum(1)
    return out[:, 0::2], out[:, 1::2], tot


def main():
    p = C.default_paths()
    d = C.H5AD(p["h5ad"])
    lab = C.cell_clusters(d, p)
    regs = build_regions(d, p)
    kseg = pd.read_csv(p["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    call = np.ones(d.n_bins)
    for r in kseg[kseg.group == "all"].itertuples():
        call[r.bin_start:r.bin_end + 1] = r.cn_call
    yt, yf, tot = region_sums(d, regs, "layers/raw")
    xt, xf, _ = region_sums(d, regs, "X")
    terc = pd.qcut(tot, 3, labels=["low", "mid", "high"]).astype(str)
    panel = set(C.shared_cells().cell)
    in_panel = np.array([c in panel for c in d.cells])
    G = tot.sum() / d.n_bins
    Gx = 1e5 / d.n_bins
    rows = []
    for j, r in enumerate(regs):
        nt, nf = len(r["target"]), len(r["flank"])
        T, F = yt[:, j], yf[:, j]
        rec = {"family": r["family"], "region": r["name"], "direction": r["direction"], "cn_call": r.get("cn_call", np.nan),
               "n_target_bins": nt, "n_flank_bins": nf,
               "rr_genome_raw": (T.sum() / nt) / G, "rr_genome_x": (xt[:, j].sum() / nt / d.n_cells) / Gx}
        if nf:
            rr = lambda m: (T[m].sum() / nt) / (F[m].sum() / nf) if F[m].sum() > 0 else np.nan  # noqa: E731
            rec["rr_flank_raw"] = rr(np.ones(d.n_cells, bool))
            rec["rr_flank_x"] = (xt[:, j].sum() / nt) / (xf[:, j].sum() / nf)
            rec["rr_flank_raw_panel300"] = rr(in_panel)
            for t in ["low", "high"]:
                rec[f"rr_flank_raw_{t}cov"] = rr(terc == t)
            rec["coverage_log_slope"] = np.log(rec["rr_flank_raw_lowcov"] / rec["rr_flank_raw_highcov"])
            cl = [rr(lab == c) for c in sorted(set(lab))]
            rec["rr_cluster_min"], rec["rr_cluster_max"] = np.nanmin(cl), np.nanmax(cl)
            n = T + F
            m = n > 0
            pi = T[m].sum() / n[m].sum()
            rec["phi"] = float(np.sum((T[m] - n[m] * pi) ** 2 / (n[m] * pi * (1 - pi))) / (m.sum() - 1))
            pic = np.zeros(d.n_cells)
            for c in set(lab):
                mc = (lab == c) & m
                pic[mc] = T[mc].sum() / n[mc].sum()
            rec["phi_cluster"] = float(np.sum((T[m] - n[m] * pic[m]) ** 2 / (n[m] * pic[m] * (1 - pic[m]))) / (m.sum() - 1 - 5))
            # zeros: expected under a homogeneous binomial split of each cell's T+F
            rec["cells_target_zero"] = int((T == 0).sum())
            rec["cells_target_zero_expected"] = float(np.sum((1 - pi) ** n[m]) + (~m).sum())
            lam = n * pi                                  # expected target count given the cell's local total
            hi_e = lam >= 5
            rec["frac_cells_lam_ge5_with_zero"] = float(((T == 0) & hi_e).sum() / max(1, hi_e.sum()))
            rec["cells_lam_ge5"] = int(hi_e.sum())
            rr0 = rec["rr_flank_raw"]
            rec["dosage_class"] = ("below one-copy-loss floor" if rr0 < 0.5 else
                                   "loss-compatible depth" if rr0 < 0.9 else
                                   "gain-compatible depth" if rr0 > 1.1 else "near flank level")
            rec["implied_loss_fraction_if_1copy"] = float(min(1, 2 * (1 - rr0))) if rr0 < 1 else 0.0
        rec["frac_bins_called_loss_all"] = float((call[r["target"]] < 1).mean()) if nt else np.nan
        rec["frac_bins_called_gain_all"] = float((call[r["target"]] > 1).mean()) if nt else np.nan
        rows.append(rec)
    df = pd.DataFrame(rows)
    nullv = df[df.family == "random"]
    for col in ["rr_flank_raw", "phi", "phi_cluster", "coverage_log_slope"]:
        ref = nullv[col].dropna().to_numpy()
        df[col + "_null_pct"] = df[col].apply(lambda v: float((ref <= v).mean() * 100) if np.isfinite(v) else np.nan)
    df["abs_log_rr_null_pct"] = df.rr_flank_raw.apply(
        lambda v: float((np.abs(np.log(nullv.rr_flank_raw)) <= abs(np.log(v))).mean() * 100) if np.isfinite(v) and v > 0 else np.nan)
    df.to_csv(HERE / "region_contrasts.tsv", sep="\t", index=False)
    pc = pd.DataFrame({"cell": d.cells, "cluster": lab, "coverage_tercile": terc, "raw_total": tot})
    for j, r in enumerate(regs):
        if r["family"] == "saved":
            pc[r["name"] + "__target"] = yt[:, j]; pc[r["name"] + "__flank"] = yf[:, j]
    pc.to_csv(HERE / "per_cell_saved_window_counts.tsv.gz", sep="\t", index=False, compression="gzip")

    # summaries
    q = lambda s: f"{s.quantile(.05):.3g} / {s.median():.3g} / {s.quantile(.95):.3g}"  # noqa: E731
    fam = []
    for (f, dirn), s in df[df.n_flank_bins > 0].groupby(["family", "direction"]):
        fam.append({"family": f, "direction": dirn, "n": len(s), "rr_flank_raw q05/med/q95": q(s.rr_flank_raw),
                    "phi q05/med/q95": q(s.phi), "phi_cluster q05/med/q95": q(s.phi_cluster),
                    "coverage_log_slope q05/med/q95": q(s.coverage_log_slope),
                    "n_below_one_copy_floor": int((s.rr_flank_raw < 0.5).sum()),
                    "n_phi_cluster_above_null_95pct": int((s.phi_cluster_null_pct > 95).sum()),
                    "n_|logrr|_above_null_95pct": int((s.abs_log_rr_null_pct > 95).sum())})
    fam = pd.DataFrame(fam)
    fam.to_csv(HERE / "family_summary.tsv", sep="\t", index=False)
    unc = nullv[(nullv.frac_bins_called_loss_all == 0) & (nullv.frac_bins_called_gain_all == 0)]
    deep = nullv[nullv.rr_flank_raw < 0.5]
    rnd = {"random_windows": len(nullv), "uncalled_random_windows": len(unc),
           "uncalled_rr_q05_q50_q95": unc.rr_flank_raw.quantile([.05, .5, .95]).round(3).tolist(),
           "uncalled_below_one_copy_floor": int((unc.rr_flank_raw < 0.5).sum()),
           "deep_random_windows_rr_lt_0.5": len(deep), "deep_random_with_majority_loss_call": int((deep.frac_bins_called_loss_all > 0.5).sum()),
           "corr_abs_log_rr_vs_log_phi_cluster_random": float(np.corrcoef(np.abs(np.log(nullv.rr_flank_raw)), np.log(nullv.phi_cluster))[0, 1]),
           "phi_cluster_median_deep_random": float(deep.phi_cluster.median()),
           "phi_cluster_median_other_random": float(nullv[nullv.rr_flank_raw >= 0.5].phi_cluster.median())}
    show = ["region", "n_target_bins", "rr_flank_raw", "rr_flank_x", "rr_flank_raw_panel300", "rr_flank_raw_lowcov", "rr_flank_raw_highcov",
            "rr_cluster_min", "rr_cluster_max", "phi", "phi_cluster", "phi_cluster_null_pct", "cells_target_zero", "cells_target_zero_expected",
            "frac_cells_lam_ge5_with_zero", "dosage_class", "abs_log_rr_null_pct"]
    saved = df[df.family == "saved"][show]
    lit = df[df.family == "literature"][["region", "direction", "n_target_bins", "rr_genome_raw", "rr_genome_x"]]
    called = df[df.family == "called"][["region", "direction"] + show[1:]]
    lit = df[df.family == "literature"][["region", "direction", "n_target_bins", "rr_genome_raw", "rr_genome_x",
                                         "frac_bins_called_loss_all", "frac_bins_called_gain_all"]]
    res = {"experiment": "04_fixed_region_contrasts", "seed": [C.SEED, 4], "generator": "numpy.random.default_rng (PCG64)",
           "n_regions": df.family.value_counts().to_dict(), "family_summary": fam.to_dict("records"), "random_window_calibration": rnd,
           "saved": saved.to_dict("records"), "literature": lit.to_dict("records"), "environment": C.environment()}
    C.write_json(HERE / "summary.json", res)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    bins_ = np.linspace(-2.5, 1.5, 60)
    ax[0].hist(np.log2(nullv.rr_flank_raw), bins_, color="#bbb", label="null windows")
    for f, c in [("called", None)]:
        s = df[df.family == f]
        ax[0].hist(np.log2(s[s.direction == "loss"].rr_flank_raw), bins_, color="#9b2c2c", alpha=.6, label="called short losses")
        ax[0].hist(np.log2(s[s.direction == "gain"].rr_flank_raw), bins_, color="#2f5d8a", alpha=.6, label="called short gains")
    for _, r in df[df.family == "saved"].iterrows():
        ax[0].axvline(np.log2(r.rr_flank_raw), color="k", lw=.8)
    ax[0].axvline(-1, color="#8a5a00", ls="--", lw=1, label="one-copy loss of diploid (rr = 0.5)")
    ax[0].set_xlabel("log2 raw rate ratio, target / flanks"); ax[0].legend(fontsize=7); ax[0].set_title("Depth contrast (black = saved windows)")
    for f, c in [("random", "#bbb"), ("called", "#9b2c2c"), ("saved", "k")]:
        s = df[df.family == f]
        ax[1].scatter(np.log2(s.rr_flank_raw), s.phi_cluster, s=10 if f != "saved" else 30, color=c, label=f, alpha=.7)
    ax[1].set_yscale("log"); ax[1].set_xlabel("log2 rr"); ax[1].set_ylabel("phi (cluster-adjusted)"); ax[1].legend(fontsize=7)
    ax[1].set_title("Per-cell heterogeneity vs contrast")
    for f, c in [("random", "#bbb"), ("called", "#9b2c2c"), ("saved", "k")]:
        s = df[df.family == f]
        ax[2].scatter(np.log2(s.rr_flank_raw), s.coverage_log_slope, s=10 if f != "saved" else 30, color=c, alpha=.7)
    ax[2].axhline(0, color="k", lw=.5); ax[2].set_xlabel("log2 rr"); ax[2].set_ylabel("log(rr low-cov / rr high-cov)")
    ax[2].set_title("Coverage dependence of the contrast")
    b64 = C.fig_to_b64(fig); plt.close(fig)
    frag = [
        ("Question", "<p>For fixed regions, how deep is the raw-count depletion or excess, is it shared by all cells (homogeneous) or "
         "concentrated in some (mixture), does it depend on cell coverage or cluster, and which explanations does its size allow? "
         "Losses and gains are both evaluated, alongside prespecified literature loci and random same-width windows.</p>"),
        ("Family summary", C.html_table(fam)),
        ("How deep do uncalled regions go?", C.html_table(pd.DataFrame([rnd]).T.reset_index().rename(columns={"index": "quantity", 0: "value"})) +
         "<div class='card'>Random 25-bin windows that the pipeline leaves entirely uncalled still span roughly "
         f"{rnd['uncalled_rr_q05_q50_q95'][0]}–{rnd['uncalled_rr_q05_q50_q95'][2]}x their flanks in raw counts, and "
         f"{rnd['uncalled_below_one_copy_floor']} of {rnd['uncalled_random_windows']} fall below 0.5x. At this scale, raw depth varies "
         "between neighbouring regions about as much as a one-copy change, so a depth contrast alone cannot be read as dosage. "
         "Per-cell dispersion rises with contrast size in random windows (correlation of |log rr| with log phi: "
         f"{rnd['corr_abs_log_rr_vs_log_phi_cluster_random']:.2f}), so elevated phi is not by itself evidence of a subclonal mixture.</div>"),
        ("Saved windows", C.html_table(saved)),
        ("Prespecified literature loci (relative to genome-wide density)", C.html_table(lit)),
        ("Figure", C.img(b64, "contrasts")),
        ("Rule-selected short called segments (all source)", C.html_table(called, max_rows=80)),
        ("Reading the columns", "<ul><li><b>rr_flank_raw</b>: raw-count ratio estimator; exposure cancels. <b>rr_flank_x</b>: same on "
         "processed X means (equal cell weights after per-cell normalization).</li><li><b>phi_cluster</b>: per-cell dispersion around "
         "cluster-specific proportions; the null-window percentile shows how unusual it is.</li><li><b>frac_cells_lam_ge5_with_zero</b>: "
         "among cells expected to have at least 5 target fragments under homogeneity, the fraction with zero. This bounds the share of cells "
         "that could carry a near-complete loss.</li><li><b>dosage_class</b>: rr below 0.5 cannot come from a one-copy loss on a diploid "
         "baseline alone; it needs a deeper loss in most cells or a lower accessibility/mappability baseline.</li></ul>"),
    ]
    C.write_html_report(HERE / "report.html", "Fixed-region count contrasts", "Experiment 04 · raw-count depth contrasts, homogeneity and dosage compatibility",
                        frag, {"command": "python -I experiments/04_fixed_region_contrasts/run.py", "seed": [C.SEED, 4],
                               "summary": "summary.json", "environment": C.environment()})
    pd.set_option("display.width", 250)
    print(fam.to_string())
    print(saved.to_string())
    print(lit.to_string())


if __name__ == "__main__":
    main()
