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
    sv = df[df.family == "saved"].set_index("region")
    chr7 = sv[sv.index.str.startswith("chr7")]
    chr9t = sv[sv.index.str.startswith("chr9") & ~sv.index.str.contains("70_100")]
    cal = rnd
    lossc = df[(df.family == "called") & (df.direction == "loss")]
    gainc = df[(df.family == "called") & (df.direction == "gain")]
    L = df[df.family == "literature"].set_index("region")
    verdict = (f"The four chr7 saved windows hold {chr7.rr_flank_raw.min():.2f}&ndash;{chr7.rr_flank_raw.max():.2f}x their flanks in raw counts, below the 0.5x floor of a one-copy loss, independent of cell coverage and shared across clusters, "
               f"while uncalled random windows alone span {cal['uncalled_rr_q05_q50_q95'][0]:.2f}&ndash;{cal['uncalled_rr_q05_q50_q95'][2]:.2f}x; this fits a shared regional baseline effect (possibly with a smaller dosage change) better than a simple clonal one-copy loss, and does not establish copy number.")
    findings = [
        f"Region families: {len(sv)} saved windows (6 loss windows + chr9 70&ndash;100 Mb comparison), {int((df.family=='called').sum())} short called segments ({len(lossc)} loss, {len(gainc)} gain; n_bins &le; 100, all-source karyogram), {len(L)} literature loci, {int((df.family=='random').sum())} random 25-bin windows (seed [{C.SEED}, 4]).",
        f"chr7 saved windows: raw rate / flank = {', '.join(f'{v:.2f}' for v in chr7.rr_flank_raw)} (processed X: {', '.join(f'{v:.2f}' for v in chr7.rr_flank_x)}), all &lt; 0.5; the null percentile of |log rr| is {chr7.abs_log_rr_null_pct.min():.0f}&ndash;{chr7.abs_log_rr_null_pct.max():.0f}%.",
        f"Uncalled random windows ({cal['uncalled_random_windows']} of {cal['random_windows']}) range {cal['uncalled_rr_q05_q50_q95'][0]}&ndash;{cal['uncalled_rr_q05_q50_q95'][2]}x (5th&ndash;95th percentile, median {cal['uncalled_rr_q05_q50_q95'][1]}); {cal['uncalled_below_one_copy_floor']} of them fall below 0.5x. {cal['deep_random_windows_rr_lt_0.5']} random windows are below 0.5x, of which {cal['deep_random_with_majority_loss_call']} are majority-called loss.",
        f"Coverage invariance: |log(rr low-coverage / rr high-coverage)| &le; {chr7.coverage_log_slope.abs().max():.3f} for the chr7 windows (&le; {sv.coverage_log_slope.abs().max():.3f} for all saved windows with flanks); cluster-specific ratios stay within {chr7.rr_cluster_min.min():.2f}&ndash;{chr7.rr_cluster_max.max():.2f} for chr7.",
        f"Among cells expected to have &ge; 5 target fragments, at most {100*chr7.frac_cells_lam_ge5_with_zero.max():.1f}% observe zero in the chr7 windows (bound on the share of cells compatible with near-complete loss); raw zero-cell counts exceed the homogeneous expectation (e.g. {int(sv.loc['chr7_67.4_70.0','cells_target_zero'])} observed vs {sv.loc['chr7_67.4_70.0','cells_target_zero_expected']:.0f} expected in chr7_67.4_70.0).",
        f"Dispersion: cluster-adjusted phi is {chr7.phi_cluster.min():.2f}&ndash;{chr7.phi_cluster.max():.2f} for chr7 windows (percentile {chr7.phi_cluster_null_pct.min():.0f}&ndash;{chr7.phi_cluster_null_pct.max():.0f} of random windows), but in random windows phi rises with contrast size (correlation of |log rr| with log phi = {cal['corr_abs_log_rr_vs_log_phi_cluster_random']:.2f}), so elevated phi is not by itself evidence of a subclonal mixture.",
        f"Called short segments: {int((lossc.rr_flank_raw<0.5).sum())} of {len(lossc)} called losses fall below the 0.5x floor (median rr {lossc.rr_flank_raw.median():.2f}); called gains have median rr {gainc.rr_flank_raw.median():.2f}.",
        f"Literature loci relative to genome-wide density (raw / processed X): MYC 8q24 {L.loc['lit_MYC_8q24','rr_genome_raw']:.2f} / {L.loc['lit_MYC_8q24','rr_genome_x']:.2f}; 19q13-qter {L.loc['lit_19q13_qter','rr_genome_raw']:.2f} / {L.loc['lit_19q13_qter','rr_genome_x']:.2f}; 9p21 zone {L.loc['lit_9p21_hemizygous_zone_approx','rr_genome_raw']:.2f} / {L.loc['lit_9p21_hemizygous_zone_approx','rr_genome_x']:.2f}; 5q13-qter {L.loc['lit_5q13_qter','rr_genome_raw']:.2f} / {L.loc['lit_5q13_qter','rr_genome_x']:.2f}.",
    ]

    def saved_note(r):
        s = f"Raw rate {r.rr_flank_raw:.2f}x flanks ({r.dosage_class}); "
        if r.rr_flank_raw < 0.5:
            s += f"deeper than a one-copy loss of a diploid baseline allows; "
        elif r.rr_flank_raw > 1.1:
            s += "above the flanks (comparison window); "
        else:
            s += f"compatible in size with a loss in about {100*r.implied_loss_fraction_if_1copy:.0f}% of cells if one-copy; "
        s += (f"coverage slope {r.coverage_log_slope:+.3f}, clusters {r.rr_cluster_min:.2f}&ndash;{r.rr_cluster_max:.2f}, cluster-adjusted phi at percentile {r.phi_cluster_null_pct:.0f} of random windows; "
              f"at most {100*r.frac_cells_lam_ge5_with_zero:.1f}% of informative cells have zero target counts.")
        return s

    sv_note = sv.reset_index()
    import html as _html
    sv_note["note"] = [_html.unescape(saved_note(r)) for r in sv_note.itertuples()]
    sv_note = sv_note[["region", "n_target_bins", "rr_flank_raw", "implied_loss_fraction_if_1copy", "note"]]
    note_html = C.html_table(sv_note).replace("<td>Raw rate", "<td style='text-align:left;white-space:normal'>Raw rate")

    def lit_note(r):
        s = f"Raw density {r.rr_genome_raw:.2f}x genome, processed X {r.rr_genome_x:.2f}x; "
        s += f"{100*r.frac_bins_called_loss_all:.0f}% of bins called loss, {100*r.frac_bins_called_gain_all:.0f}% gain. "
        if r.direction == "gain" and r.rr_genome_x < 1.0 and r.rr_genome_raw > 1.2:
            s += "The raw excess is not visible after GC correction."
        elif r.direction == "gain" and r.rr_genome_raw > 1.2:
            s += "Raw excess is retained after GC correction."
        elif r.direction == "loss" and r.rr_genome_raw < 0.8:
            s += "Depleted in raw counts, as expected for a hemizygous zone."
        else:
            s += "Modest contrast."
        return s

    lit_note_df = L.reset_index()
    lit_note_df["note"] = [_html.unescape(lit_note(r)) for r in lit_note_df.itertuples()]
    lit_note_df = lit_note_df[["region", "direction", "n_target_bins", "rr_genome_raw", "rr_genome_x", "note"]]
    lit_html = C.html_table(lit_note_df).replace("<td>Raw density", "<td style='text-align:left;white-space:normal'>Raw density")
    savedf = df[df.family == "saved"]
    depth_cols = savedf[["region", "n_target_bins", "rr_flank_raw", "rr_flank_x", "rr_flank_raw_panel300", "abs_log_rr_null_pct", "dosage_class"]]
    het_cols = savedf[["region", "rr_flank_raw_lowcov", "rr_flank_raw_highcov", "rr_cluster_min", "rr_cluster_max", "phi", "phi_cluster", "phi_cluster_null_pct"]]
    zero_cols = savedf[["region", "cells_target_zero", "cells_target_zero_expected", "frac_cells_lam_ge5_with_zero"]]
    called_show = df[df.family == "called"][["region", "direction", "n_target_bins", "rr_flank_raw", "rr_flank_x", "coverage_log_slope", "phi_cluster", "phi_cluster_null_pct", "dosage_class"]]
    sections = [
        ("Question and why it matters",
         "<p><b>Question.</b> For fixed regions, how deep is the raw-count depletion or excess, is it shared by all cells (homogeneous) or concentrated in some (mixture), does it depend on cell coverage or cluster, and which explanations does its size allow?</p>"
         "<p><b>Why it matters.</b> The recommendation&rsquo;s central uncertainty is whether the chr7 dips behave like a simple clonal or subclonal one-copy loss, and whether raw depth contrasts at this scale carry dosage information at all. "
         "This experiment supplies the numbers for the per-region evidence layer (depth contrast against random windows, coverage dependence, compatibility with a one-copy change, per-cell mixture bound) and for the stage-2 argument that a cross-sample baseline is the decisive missing experiment.</p>"),
        ("Inputs",
         f"<ul><li><b>Matrix.</b> <code>count_matrix.h5ad</code> ({d.n_cells:,} cells &times; {d.n_bins:,} bins): <code>layers/raw</code> (pre-GC counts, used for all primary statistics) and <code>X</code> (processed; secondary). Units: fragments per 100 kb bin.</li>"
         "<li><b>Clusters.</b> WNN clusters (<code>wnn_leiden_0.3</code>, 5 clusters) from <code>integrated.h5mu</code>; coverage terciles by raw total per cell.</li>"
         f"<li><b>Region families</b> (defined by rule before computing any statistic): "
         f"<i>saved</i> = the 6 saved chr7/chr9 windows plus the historical chr9 70&ndash;100 Mb comparison, with frozen 50-retained-bin flanks on each side (shared <code>regions.tsv</code>/<code>region_bins.tsv</code>); "
         f"<i>called</i> = every segment of the <code>all</code> source in the saved genome-pool karyogram with n_bins &le; 100 and cn_call &ne; 1 ({len(lossc)} loss/weak loss, {len(gainc)} gain/weak gain; flanks up to 50 retained bins each side within the chromosome); "
         "<i>literature</i> = four published MCF10-A loci (9p21 hemizygous zone chr9:20&ndash;24 Mb approximate, MYC chr8:127.0&ndash;128.5 Mb, 19q13-qter chr19 &ge; 35 Mb, 5q13-qter chr5 &ge; 70 Mb), compared with the genome-wide density only (no flanks); "
         f"<i>random</i> = {int((df.family=='random').sum())} random windows of 25 consecutive retained bins without internal gap (same flank rule). Random windows are not neutral truth: they can contain events.</li>"
         "<li><b>300-cell panel.</b> Frozen shared panel (<code>cells.tsv</code>) used for a consistency column.</li>"
         "<li><b>Not available.</b> Independent CN truth, mappability, cross-sample baselines, SNP-level allele data (experiment 06 uses the allelic matrices).</li></ul>"),
        ("Method",
         "<ol><li>For every region sum raw counts per cell over target bins and over flank bins (all 3,460 cells).</li>"
         "<li>Rate ratio <code>rr_flank</code> is a ratio of sums, so cell exposure cancels. <code>rr_genome</code> compares target density to genome-wide density.</li>"
         "<li>Repeat rr in raw-coverage terciles and WNN clusters; <code>coverage_log_slope = log(rr_low / rr_high)</code>.</li>"
         "<li>Per-cell dispersion <code>phi</code>: Pearson dispersion of the target count given the cell&rsquo;s target+flank total around one common binomial proportion; <code>phi_cluster</code> uses a proportion per cluster (5 degrees of freedom removed).</li>"
         "<li>Zero test: cells with zero target counts versus the homogeneous binomial expectation; among cells expected to have &ge; 5 target fragments, the fraction observing zero (an upper bound on the share of cells that could carry a near-complete loss).</li>"
         "<li>Compare every statistic with the null distribution of the random windows (percentiles).</li></ol>"
         "<pre>rr_flank   = (sum_cells target / n_target_bins) / (sum_cells flank / n_flank_bins)\nrr_genome  = (sum target / n_target_bins) / (total counts / n_bins)\nn_c = T_c + F_c;  pi = sum T / sum n\nphi       = sum_c (T_c - n_c pi)^2 / (n_c pi (1 - pi)) / (n_cells - 1)\nlam_c     = n_c * pi;  frac_zero = #{T_c = 0 and lam_c >= 5} / #{lam_c >= 5}\nimplied_loss_fraction_if_1copy = min(1, 2 (1 - rr_flank))   # share of cells with one-copy loss that would give rr</pre>"
         f"<p>Dosage rule: rr_flank &lt; 0.5 cannot arise from a one-copy loss of a diploid baseline with unchanged accessibility/mappability alone. Random windows: seed tuple ({C.SEED}, 4), numpy PCG64 (<code>default_rng</code>); candidates are 25-bin runs with no internal gap on every chromosome. "
         "Percentiles are of the random-window distribution.</p>"),
        ("Results", "", 2),
        ("Family summary", C.html_table(fam) + C.caption("One row per family and direction. Columns give the 5th / median / 95th percentile of each statistic and counts above the random-window 95th percentile. Compare the random row against the saved and called rows: random windows span a wide range of rr and phi; saved losses have much lower rr and higher phi."), 3),
        ("How deep do uncalled regions go?",
         C.html_table(pd.DataFrame([rnd]).T.reset_index().rename(columns={"index": "quantity", 0: "value"})) +
         C.caption("Calibration against the 300 random windows. <code>uncalled_*</code> restricts to random windows with no bin called loss or gain in the <code>all</code> karyogram; <code>deep_random</code> are windows with rr below 0.5.") +
         "<div class='card'>Random 25-bin windows that the pipeline leaves entirely uncalled still span roughly "
         f"{rnd['uncalled_rr_q05_q50_q95'][0]}&ndash;{rnd['uncalled_rr_q05_q50_q95'][2]}x their flanks in raw counts, and "
         f"{rnd['uncalled_below_one_copy_floor']} of {rnd['uncalled_random_windows']} fall below 0.5x. At this scale, raw depth varies "
         "between neighbouring regions about as much as a one-copy change, so a depth contrast alone cannot be read as dosage. "
         "Per-cell dispersion rises with contrast size in random windows (correlation of |log rr| with log phi: "
         f"{rnd['corr_abs_log_rr_vs_log_phi_cluster_random']:.2f}; median cluster-adjusted phi {rnd['phi_cluster_median_deep_random']:.2f} in deep random windows against {rnd['phi_cluster_median_other_random']:.2f} in the others), so elevated phi is not by itself evidence of a subclonal mixture.</div>", 3),
        ("Saved windows: depth contrast", C.html_table(depth_cols) + C.caption("Raw (<code>rr_flank_raw</code>), processed-X (<code>rr_flank_x</code>) and 300-cell-panel ratios of target to flanks, and where the absolute log ratio sits among the 300 random windows (<code>abs_log_rr_null_pct</code>, 100 = more extreme than all). Look for ratios below 0.5 and null percentiles near 100; X is shallower than raw for the low-GC chr7 windows."), 3),
        ("Saved windows: homogeneity and coverage", C.html_table(het_cols) + C.caption("Ratios in the lowest and highest coverage tercile (equal values mean the contrast does not depend on coverage), the lowest and highest cluster-specific ratios, and per-cell dispersion with its percentile among random windows. Look for narrow cluster ranges (shared by all clusters) and phi elevated relative to random windows."), 3),
        ("Saved windows: zero-count test", C.html_table(zero_cols) + C.caption("<code>cells_target_zero</code> cells have no target fragments, against the number expected if every cell carried the same proportion; <code>frac_cells_lam_ge5_with_zero</code> restricts to cells expected to have at least 5 target fragments, where a zero is informative. Small values bound how many cells could carry a near-complete loss."), 3),
        ("Prespecified literature loci (relative to genome-wide density)", C.html_table(lit) + C.caption("Raw and processed-X density of each published locus divided by the genome-wide density, with the share of its bins called loss or gain. Compare raw with processed values: GC correction can remove a raw excess."), 3),
        ("Rule-selected short called segments (all source)", C.html_table(called_show, max_rows=80) + C.caption("The first 80 short called segments (loss and gain). Full table with all columns: <code>region_contrasts.tsv</code>. Look at <code>dosage_class</code> and <code>phi_cluster_null_pct</code>: called losses frequently sit below the one-copy floor, and called gains are rarely beyond random-window dispersion."), 3),
        ("Figure", C.img(b64, "contrasts") + C.caption("Left: log2 raw rate ratio for random windows (grey), called short losses (red), called short gains (blue); black lines are the saved windows; dashed line is the one-copy-loss level (rr = 0.5). Middle: cluster-adjusted dispersion against contrast for random (grey), called (red) and saved (black) windows. "
                                                   "Right: coverage dependence of the contrast; points near the zero line are coverage-invariant."), 3),
        ("Per-region detail",
         "<h3>Saved windows</h3>" + note_html + C.caption("One row per saved window; the note restates the table values in plain language. <code>implied_loss_fraction_if_1copy</code> = min(1, 2(1 &minus; rr)): the share of cells that would need a one-copy loss to give this ratio (1 means a one-copy loss in all cells is not enough).") +
         "<h3>Literature loci</h3>" + lit_html + C.caption("One row per published locus (approximate extents). Notes restate the density ratios and call shares.")),
        ("Sensitivity",
         "<p>Checks included: 300-cell panel versus all cells (<code>rr_flank_raw_panel300</code>), raw versus processed X, coverage terciles and clusters. Not varied: the 50-bin flank width, the 25-bin random-window width and the 5-fragment threshold of the zero test (a &plusmn;20% width sensitivity is a recommended next check); the percentile of a saved window depends on the random-window distribution, which includes real events.</p>"),
        ("Limits and alternative explanations",
         "<ul><li>Random windows are not neutral truth; they sample the same genome that contains the events, so null percentiles are calibration, not proof.</li>"
         "<li>The one-copy floor assumes a diploid baseline and unchanged accessibility/mappability; a ratio below 0.5 is also expected from a regional baseline offset (accessibility or mappability), possibly combined with a smaller dosage change. These cannot be separated without a cross-sample baseline or mappability track.</li>"
         "<li>Homogeneity (cluster-invariant ratios, bounded zero fractions) argues against a small subclone carrying a deep loss but cannot distinguish a shared baseline from a shared clonal change.</li>"
         "<li>Literature loci are approximate and compared with genome-wide density only; the 9p21 extent is a ~4 Mb approximation.</li>"
         "<li>Called segments are outputs of the same pipeline, so &lsquo;called&rsquo; windows are not independent evidence.</li></ul>"),
    ]
    glossary = {
        "rr_flank_raw": "Ratio of mean raw count per target bin to mean raw count per flank bin, summed over all cells (exposure cancels).",
        "rr_flank_x": "The same ratio on processed X (per-cell normalized, GC-corrected).",
        "rr_flank_raw_panel300": "rr_flank_raw restricted to the frozen 300-cell panel.",
        "rr_flank_raw_lowcov / rr_flank_raw_highcov": "rr_flank_raw among cells in the lowest / highest raw-coverage tercile.",
        "coverage_log_slope": "log(rr in low-coverage cells / rr in high-coverage cells); 0 means coverage-invariant.",
        "rr_cluster_min / rr_cluster_max": "Smallest and largest cluster-specific rr_flank_raw over the 5 WNN clusters.",
        "rr_genome_raw / rr_genome_x": "Target density divided by genome-wide density, in raw counts / processed X.",
        "phi": "Pearson dispersion of per-cell target counts given target+flank totals around one common proportion; 1 = Poisson-like homogeneous.",
        "phi_cluster": "The same with a separate proportion for each WNN cluster.",
        "phi_cluster_null_pct": "Percent of random windows with phi_cluster at or below this value.",
        "abs_log_rr_null_pct": "Percent of random windows whose |log rr| is at or below this region&rsquo;s.",
        "cells_target_zero": "Cells with zero target fragments.",
        "cells_target_zero_expected": "Expected number of such cells under one common binomial proportion.",
        "frac_cells_lam_ge5_with_zero": "Among cells expected to have &ge; 5 target fragments, the fraction observing zero.",
        "dosage_class": "Label from rr_flank_raw: below 0.5 = below one-copy-loss floor; 0.5&ndash;0.9 = loss-compatible depth; above 1.1 = gain-compatible depth; otherwise near flank level.",
        "implied_loss_fraction_if_1copy": "min(1, 2 (1 &minus; rr)): share of cells needing a one-copy loss to explain the ratio.",
        "frac_bins_called_loss_all / frac_bins_called_gain_all": "Share of target bins called loss / gain in the <code>all</code> karyogram.",
        "n_target_bins / n_flank_bins": "Retained bins in the target / flanks.",
        "family / direction": "Region family (saved, called, literature, random) and expected direction.",
        "n_below_one_copy_floor": "Regions in the group with rr_flank_raw &lt; 0.5.",
        "n_phi_cluster_above_null_95pct / n_|logrr|_above_null_95pct": "Regions above the random-window 95th percentile of phi_cluster / of |log rr|.",
        "q05/med/q95": "5th percentile / median / 95th percentile.",
    }
    C.write_html_report(HERE / "report.html", "Fixed-region count contrasts", "Experiment 04 &middot; raw-count depth contrasts, homogeneity and dosage compatibility",
                        sections, {"command": "python -I experiments/04_fixed_region_contrasts/run.py", "seed": [C.SEED, 4],
                                   "summary": "summary.json", "environment": C.environment()},
                        key_findings=findings, verdict=verdict, glossary=glossary,
                        supports=["The chr7 saved windows are coverage-invariant, similar in every cluster and have few zero-count cells among informative cells, which is consistent with a shared regional depth deficit rather than a small deep-loss subclone.",
                                  "Their raw ratios are below the 0.5x floor of a one-copy loss from a diploid baseline with unchanged accessibility, so a one-copy loss alone does not account for the depth.",
                                  "Depth contrasts of one-copy size are common in uncalled random windows, so depth alone is weak evidence of dosage."],
                        not_supports=["The copy-number state of any saved window.",
                                      "Whether the chr7 plateaus are mappability, accessibility or dosage (no cross-sample baseline or mappability track).",
                                      "That elevated dispersion indicates subclonal mixture: random deep windows show the same association."],
                        decision_impact="Supports recommendation ranks 1 and 3: report raw-count contrast with its percentile among same-width random windows, the coverage-invariance and one-copy-floor labels, and obtain a cross-sample baseline and mappability before changing any formula. "
                                        "It argues against reading the chr7 dips as simple clonal one-copy losses, but does not rule out a smaller dosage change on top of a shared baseline.",
                        next_steps=["Run the same binning and QC on unrelated samples, with a read-length-matched mappability track, to test whether the chr7 plateaus recur.",
                                    "Repeat with 25&plusmn;20% window widths and different flank widths to test label stability.",
                                    "Add an independent DNA-based copy-number measurement for MCF10-A (or a sample with truth) to calibrate the one-copy floor."])
    pd.set_option("display.width", 250)
    print(fam.to_string())
    print(saved.to_string())
    print(lit.to_string())


if __name__ == "__main__":
    main()
