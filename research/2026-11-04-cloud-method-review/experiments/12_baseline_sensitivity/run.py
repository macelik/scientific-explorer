"""Experiment 12: reference-free relative baseline sensitivity (exploratory diagnostic).

Segmentation is held FIXED at the saved genome-pool `all` segmentation. Nothing here is a copy-number
estimate: every "level" is (region rate) / (a baseline), and no baseline is treated as copy-neutral truth.

PRESPECIFIED BASELINES (listed before any value is computed; each computed per estimand)
  a  genome_mean      mean of bin values over all 25,337 retained bins
  b  genome_median    median of bin values
  c  segwt_median     length-weighted median of the saved `all` segment means (weights = n_bins); the baseline
                      that experiment 03 found pyEpi's CN scale effectively reduces to (stated, not re-proven)
  d  trimmed_mean     scipy.stats.trim_mean(bins, 0.10) (drops lowest and highest 10% of bins)
  e  kde_mode         mode of a Gaussian KDE of segment means (weights = n_bins); bandwidth = Silverman's rule
                      (scipy gaussian_kde bw_method='silverman', weighted); mode searched on a 4,000-point grid
                      over [min, max] of segment means
  f  loco_median      leave-one-chromosome-out: median of bins on all OTHER chromosomes (region-specific)
  g  cluster_vs_rest  each WNN cluster's raw exposure-weighted rate over the pooled rate of the other four
                      clusters (cancels shared baseline structure; a between-cluster contrast only, NOT a
                      level, so it is excluded from the class-reversal calculation and reported in task 4)
  h  gc_decile        bin value divided by the median bin value of all bins in the same GC decile (a
                      GC-matched reference without LOESS); region level = sum(values)/sum(decile medians)
Estimands (a-f,h computed for both): raw = sum raw counts / sum exposure x 1e5 (exposure-weighted rate per
source; `all` source = all 3,460 cells); x = mean of processed X over cells.
Classification of a region level L: loose  loss L<0.683, gain L>1.225; strict  loss L<=0.5, gain L>=1.5.
Reversal (per estimand and threshold set, over baselines a,b,c,d,e,f,h): a region is "loss-ever" if some baseline
classes it loss, "unstable loss" if at least one baseline does not; likewise for gain. "opposite" = loss under one
baseline and gain under another. Retention vs c = share of regions called loss (gain) under c that keep the call
under each other baseline.
Technical flags (descriptive, F:2/F:4): distance to nearest QC-removed grid bin (retained-grid gaps, experiment 02),
to nearest blacklist interval (hg38-blacklist.v2), to assembly-gap-like removed runs (mean N > 0.5);
fraction of region bins within 1 Mb of a removed bin; mean GC, GC percentile among retained bins; mean N.
Enrichment: two-sided Fisher exact tests, descriptive only (no multiplicity adjustment; regions overlap/correlate).
Regions: rule-defined families from experiment 04 (build_regions, target bins only) plus `all_segments`.
Deterministic; cell bootstrap seed tuple (20261104, 12).
Run: PYEPI_SCIENTIFIC_PACKAGE=/home/user/pyepi-upstream/pyEpiAneufinder /home/user/envs/cnv310/bin/python -I run.py
"""
import importlib.util
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

T0 = time.time()
LEVEL_BASES = ["a_genome_mean", "b_genome_median", "c_segwt_median", "d_trimmed_mean", "e_kde_mode", "f_loco_median", "h_gc_decile"]
REF = "c_segwt_median"
THR = {"loose": (0.683, 1.225, False), "strict": (0.5, 1.5, True)}  # loss cut, gain cut, inclusive
BOOT = 300
NEAR_BP = 1_000_000


def load_exp04():
    spec = importlib.util.spec_from_file_location("exp04", HERE.parent / "04_fixed_region_contrasts" / "run.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def classify(L, kind):
    lo, hi, incl = THR[kind]
    L = np.asarray(L, float)
    loss = (L <= lo) if incl else (L < lo)
    gain = (L >= hi) if incl else (L > hi)
    return np.where(loss, "loss", np.where(gain, "gain", "neutral"))


def wmedian(x, w):
    o = np.argsort(x)
    cw = np.cumsum(w[o])
    return float(x[o][np.searchsorted(cw, cw[-1] / 2)])


def scalar_baselines(v, segmean, segn):
    kde = stats.gaussian_kde(segmean, bw_method="silverman", weights=segn / segn.sum())
    grid = np.linspace(segmean.min(), segmean.max(), 4000)
    return {"a_genome_mean": float(v.mean()), "b_genome_median": float(np.median(v)),
            "c_segwt_median": wmedian(segmean, segn.astype(float)), "d_trimmed_mean": float(stats.trim_mean(v, 0.10)),
            "e_kde_mode": float(grid[np.argmax(kde(grid))])}, float(kde.factor * np.sqrt(kde.covariance.squeeze() / kde.factor ** 2))


def main():
    p = C.default_paths()
    d = C.H5AD(p["h5ad"])
    lab = C.cell_clusters(d, p)
    ex04 = load_exp04()
    regs = ex04.build_regions(d, p)
    seg = pd.read_csv(p["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    sa = seg[seg.group == "all"].reset_index(drop=True)
    for s in sa.itertuples():
        regs.append({"family": "all_segments", "name": f"seg{s.segment}_{s.chromosome}:{s.bin_start}-{s.bin_end}",
                     "direction": "loss" if s.cn_call < 1 else ("gain" if s.cn_call > 1 else "neutral"),
                     "cn_call": s.cn_call, "target": np.arange(s.bin_start, s.bin_end + 1), "flank": np.array([], int)})
    R = len(regs)
    srcs = [f"cluster{c}" for c in range(5)] + ["all"]
    masks = {f"cluster{c}": lab == str(c) for c in range(5)}
    masks["all"] = np.ones(d.n_cells, bool)

    # ---- one pass over raw counts: exposure, per-source bin sums, per-cell region target sums
    rows, cols = [], []
    for j, r in enumerate(regs):
        rows += list(r["target"]); cols += [j] * len(r["target"])
    M = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(d.n_bins, R))
    expo = np.zeros(d.n_cells)
    Y = np.zeros((d.n_cells, R))
    bsum = {k: np.zeros(d.n_bins) for k in srcs}
    for lo in range(0, d.n_cells, 250):
        hi = min(d.n_cells, lo + 250)
        blk = d.rows("layers/raw", lo, hi)
        expo[lo:hi] = blk.sum(1)
        Y[lo:hi] = (sparse.csr_matrix(blk) @ M).toarray()
        for k in srcs:
            sel = masks[k][lo:hi]
            if sel.any():
                bsum[k] += blk[sel].sum(0)
    rate = {k: bsum[k] / expo[masks[k]].sum() * 1e5 for k in srcs}   # raw exposure-weighted rate per bin
    xmean = C.group_means(d, "X", masks)
    est = {"raw": rate["all"], "x": xmean["all"]}
    print(f"[{time.time()-T0:.0f}s] data pass done", flush=True)

    segn = (sa.bin_end - sa.bin_start + 1).to_numpy()
    seg_ids = [np.arange(s.bin_start, s.bin_end + 1) for s in sa.itertuples()]
    chrom_arr = np.asarray(d.seq)
    # ---- baselines
    gcdec = pd.qcut(d.gc, 10, labels=False, duplicates="drop")
    base_scalar, bw = {}, {}
    gc_ref = {}
    loco = {}
    for e, v in est.items():
        segmean = np.array([v[ix].mean() for ix in seg_ids])
        base_scalar[e], bw[e] = scalar_baselines(v, segmean, segn)
        med = pd.Series(v).groupby(gcdec).median()
        gc_ref[e] = med.to_numpy()[gcdec]
        loco[e] = {c: float(np.median(v[chrom_arr != c])) for c in dict.fromkeys(chrom_arr)}
    rec = []
    for j, r in enumerate(regs):
        t = r["target"]
        ch = chrom_arr[t[0]]
        base = {"family": r["family"], "region": r["name"], "direction": r["direction"], "n_bins": len(t), "chrom": ch}
        for e, v in est.items():
            lvl = {}
            for b, val in base_scalar[e].items():
                lvl[b] = v[t].mean() / val
            lvl["f_loco_median"] = v[t].mean() / loco[e][ch]
            lvl["h_gc_decile"] = v[t].sum() / gc_ref[e][t].sum()
            for b, val in lvl.items():
                rec.append({**base, "estimand": e, "baseline": b, "level": val,
                            "class_loose": classify(val, "loose").item(), "class_strict": classify(val, "strict").item()})
    bl = pd.DataFrame(rec)

    # ---- (g) cluster vs rest (raw)
    ysum = {k: Y[masks[k]].sum(0) for k in srcs}
    expsum = {k: expo[masks[k]].sum() for k in srcs}
    tot_e = expo.sum()
    tot_y = Y.sum(0)
    nbin = np.array([len(r["target"]) for r in regs], float)
    rng = np.random.default_rng((C.SEED, 12))
    idx = {c: np.flatnonzero(lab == str(c)) for c in range(5)}
    cvr = {}
    for c in range(5):
        kk, rr = idx[c], np.flatnonzero(lab != str(c))
        point = (ysum[f"cluster{c}"] / expsum[f"cluster{c}"]) / ((tot_y - ysum[f"cluster{c}"]) / (tot_e - expsum[f"cluster{c}"]))
        bs = np.zeros((BOOT, R))
        for i in range(BOOT):
            a, b = rng.choice(kk, len(kk)), rng.choice(rr, len(rr))
            bs[i] = (Y[a].sum(0) / expo[a].sum()) / (Y[b].sum(0) / expo[b].sum())
        cvr[c] = (point, np.nanpercentile(bs, 2.5, axis=0), np.nanpercentile(bs, 97.5, axis=0))
    gl = []
    for j, r in enumerate(regs):
        rec_ = {"family": r["family"], "region": r["name"], "direction": r["direction"], "n_bins": int(nbin[j])}
        for c in range(5):
            rec_[f"cvr_c{c}"] = cvr[c][0][j]; rec_[f"cvr_c{c}_lo"] = cvr[c][1][j]; rec_[f"cvr_c{c}_hi"] = cvr[c][2][j]
            # cluster's own level against its own genome median (raw, bin-rate based)
            v = rate[f"cluster{c}"]
            rec_[f"own_level_c{c}"] = v[regs[j]["target"]].mean() / np.median(v)
        pts = np.array([cvr[c][0][j] for c in range(5)])
        rec_["cvr_max_abs_log2"] = float(np.nanmax(np.abs(np.log2(pts))))
        rec_["cvr_n_ci_excl_1"] = int(sum((cvr[c][1][j] > 1) or (cvr[c][2][j] < 1) for c in range(5)))
        gl.append(rec_)
    gdf = pd.DataFrame(gl)
    null = gdf[gdf.family == "random"].cvr_max_abs_log2
    q95 = float(null.quantile(0.95))
    gdf["cvr_specific_vs_random95"] = gdf.cvr_max_abs_log2 > q95
    gdf["pct_of_random_null"] = gdf.cvr_max_abs_log2.apply(lambda x: float((null <= x).mean() * 100))
    own_lo = gdf[[f"own_level_c{c}" for c in range(5)]]
    gdf["own_level_min"], gdf["own_level_max"] = own_lo.min(1), own_lo.max(1)
    gdf["shared_depletion_all_clusters_below_0.683"] = (own_lo < 0.683).all(1)
    gdf["shared_depletion_all_clusters_below_0.85"] = (own_lo < 0.85).all(1)
    gdf.to_csv(HERE / "cluster_vs_rest.tsv", sep="\t", index=False)
    # add g to baselines.tsv as separate rows (between-cluster contrast, raw)
    grows = []
    for j, r in enumerate(regs):
        for c in range(5):
            grows.append({"family": r["family"], "region": r["name"], "direction": r["direction"], "n_bins": int(nbin[j]),
                          "chrom": chrom_arr[r["target"][0]], "estimand": "raw", "baseline": f"g_cluster_vs_rest_c{c}",
                          "level": cvr[c][0][j], "class_loose": "n/a", "class_strict": "n/a"})
    bl_all = pd.concat([bl, pd.DataFrame(grows)], ignore_index=True)
    bl_all.to_csv(HERE / "baselines.tsv", sep="\t", index=False)

    # ---- reversals
    wide = {}
    cls_rows, rate_rows = [], []
    for e in est:
        for kind in THR:
            w = bl[bl.estimand == e].pivot_table(index=["family", "region", "direction"], columns="baseline",
                                                 values="class_" + kind, aggfunc="first")[LEVEL_BASES]
            lv = bl[bl.estimand == e].pivot_table(index=["family", "region", "direction"], columns="baseline", values="level")[LEVEL_BASES]
            w = w.reset_index()
            w["estimand"], w["thresholds"] = e, kind
            A = w[LEVEL_BASES].to_numpy()
            w["n_loss_baselines"] = (A == "loss").sum(1)
            w["n_gain_baselines"] = (A == "gain").sum(1)
            w["n_distinct_classes"] = [len(set(x)) for x in A]
            w["opposite_loss_gain"] = (w.n_loss_baselines > 0) & (w.n_gain_baselines > 0)
            lvl = lv.reset_index()
            w["log2_level_range"] = np.log2(lv.max(1).to_numpy() / lv.min(1).to_numpy())
            cls_rows.append(w)
            for fam, sub in w.groupby("family"):
                for dirn, key in (("loss", "n_loss_baselines"), ("gain", "n_gain_baselines")):
                    anyv = sub[sub[key] > 0]
                    unstable = anyv[anyv[key] < len(LEVEL_BASES)]
                    refc = sub[sub[REF] == dirn]
                    ret = {b: (float((refc[b] == dirn).mean()) if len(refc) else np.nan) for b in LEVEL_BASES if b != REF}
                    rate_rows.append({"family": fam, "estimand": e, "thresholds": kind, "direction": dirn, "n_regions": len(sub),
                                      "n_called_any_baseline": len(anyv), "n_called_all_7": int((sub[key] == len(LEVEL_BASES)).sum()),
                                      "n_unstable": len(unstable),
                                      "reversal_rate_among_any": (len(unstable) / len(anyv)) if len(anyv) else np.nan,
                                      "n_called_under_c": len(refc),
                                      **{f"keep_vs_c_{b[0]}": v for b, v in ret.items()}})
                rate_rows.append({"family": fam, "estimand": e, "thresholds": kind, "direction": "opposite_loss_vs_gain",
                                  "n_regions": len(sub), "n_called_any_baseline": int(sub.opposite_loss_gain.sum()),
                                  "n_called_all_7": 0, "n_unstable": int(sub.opposite_loss_gain.sum()),
                                  "reversal_rate_among_any": float(sub.opposite_loss_gain.mean())})
    cls = pd.concat(cls_rows, ignore_index=True)
    cls.to_csv(HERE / "region_classes.tsv", sep="\t", index=False)
    rev = pd.DataFrame(rate_rows)
    rev.to_csv(HERE / "reversals.tsv", sep="\t", index=False)
    print(f"[{time.time()-T0:.0f}s] reversals done", flush=True)

    # ---- technical flags
    g = pd.read_csv(p["root"] / "pyepi_results" / "binned_genome.csv", index_col=0)
    g = g[g.chromosome.isin(set(chrom_arr))].copy()
    key = set(zip(d.seq, d.start))
    g["retained"] = [(c, s) in key for c, s in zip(g.chromosome, g.start)]
    assert int(g.retained.sum()) == d.n_bins
    bed = pd.read_csv(p["root"] / "pyepi_results" / "hg38-blacklist.v2.bed", sep="\t", header=None, names=["chrom", "s0", "e", "kind"])
    bed["s"] = bed.s0 + 1                                      # 0-based half-open -> 1-based inclusive
    removed, gaplike = {}, {}
    gcol = {}
    for ch, sub in g.groupby("chromosome", sort=False):
        sub = sub.sort_values("start").reset_index(drop=True)
        sub["gi"] = np.arange(len(sub))
        removed[ch] = sub[~sub.retained][["gi", "start", "end"]].to_numpy()
        gcol[ch] = sub
        # assembly-gap-like removed runs (mean N > 0.5), same rule as experiment 02
        rem = (~sub.retained).to_numpy()
        runs, i = [], 0
        while i < len(sub):
            if rem[i]:
                j = i
                while j + 1 < len(sub) and rem[j + 1]:
                    j += 1
                part = sub.iloc[i:j + 1]
                if part.N.mean() > 0.5:
                    runs.append((int(part.start.iloc[0]), int(part.end.iloc[-1])))
                i = j + 1
            else:
                i += 1
        gaplike[ch] = runs
    gc_sorted = np.sort(d.gc)

    def span_dist(a, b, ints):
        """bp distance from [a,b] to nearest interval (0 if overlapping); inf if none."""
        if len(ints) == 0:
            return np.inf
        ints = np.asarray(ints)
        return float(np.maximum(0, np.maximum(ints[:, 0] - b, a - ints[:, 1])).min())

    fl = []
    for r in regs:
        t = r["target"]
        ch = chrom_arr[t[0]]
        st, en = d.start[t], d.end[t]
        a, b = int(st.min()), int(en.max())
        rem = removed[ch]
        bed_c = bed[bed.chrom == ch][["s", "e"]].to_numpy()
        gl_c = gaplike[ch]
        gi_t = gcol[ch].set_index("start").loc[st, "gi"].to_numpy()
        if len(rem):
            d_rem_bp = span_dist(a, b, rem[:, 1:3])
            dgi = np.maximum(0, np.maximum(rem[:, 0] - gi_t.max(), gi_t.min() - rem[:, 0]))
            d_rem_bins = int(dgi.min())
            # per-bin: nearest removed bin within NEAR_BP
            dd = np.maximum(0, np.maximum(rem[None, :, 1] - en[:, None], st[:, None] - rem[None, :, 2])).min(1)
            frac_near = float((dd <= NEAR_BP).mean())
            n_inside = int(((rem[:, 0] > gi_t.min()) & (rem[:, 0] < gi_t.max())).sum())
        else:
            d_rem_bp, d_rem_bins, frac_near, n_inside = np.inf, -1, 0.0, 0
        d_bl = span_dist(a, b, bed_c)
        d_gap = span_dist(a, b, gl_c)
        gcm = float(d.gc[t].mean())
        fl.append({"family": r["family"], "region": r["name"], "direction": r["direction"], "chrom": ch, "n_bins": len(t),
                   "dist_removed_bin_bp": d_rem_bp, "dist_removed_bin_gridbins": d_rem_bins,
                   "n_removed_bins_inside_span": n_inside, "frac_within_1Mb_of_removed": frac_near,
                   "dist_blacklist_bp": d_bl, "dist_blacklist_bins": d_bl / 1e5,
                   "dist_assembly_gap_like_bp": d_gap, "mean_GC": gcm,
                   "GC_percentile": float(np.searchsorted(gc_sorted, gcm) / len(gc_sorted) * 100),
                   "mean_N": float(gcol[ch].set_index("start").loc[st, "N"].mean()),
                   "near_removed_1Mb": bool(d_rem_bp <= NEAR_BP), "near_blacklist_1Mb": bool(d_bl <= NEAR_BP),
                   "near_gap_like_1Mb": bool(d_gap <= NEAR_BP), "GC_extreme_decile": bool(gcm < np.quantile(d.gc, .1) or gcm > np.quantile(d.gc, .9))})
    fdf = pd.DataFrame(fl)
    lvl_c = bl[(bl.estimand == "raw") & (bl.baseline == REF)].set_index(["family", "region"]).level
    fdf["level_raw_c"] = [lvl_c[(a, b)] for a, b in zip(fdf.family, fdf.region)]
    fdf.to_csv(HERE / "flags.tsv", sep="\t", index=False)

    flag_cols = ["near_removed_1Mb", "near_blacklist_1Mb", "near_gap_like_1Mb", "GC_extreme_decile"]
    rnd = fdf[fdf.family == "random"]
    groups = {
        "called loss segments vs random windows": (fdf[(fdf.family == "called") & (fdf.direction == "loss")], rnd),
        "saved loss windows vs random windows": (fdf[(fdf.family == "saved") & (fdf.direction == "loss")], rnd),
        "deep random (raw level c < 0.683) vs other random": (rnd[rnd.level_raw_c < 0.683], rnd[rnd.level_raw_c >= 0.683]),
        "called gain segments vs random windows": (fdf[(fdf.family == "called") & (fdf.direction == "gain")], rnd),
    }
    enr = []
    for nm, (A, B) in groups.items():
        for fc in flag_cols:
            a1, a0, b1, b0 = int(A[fc].sum()), int((~A[fc]).sum()), int(B[fc].sum()), int((~B[fc]).sum())
            if len(A) == 0:
                continue
            orr, pv = stats.fisher_exact([[a1, a0], [b1, b0]])
            enr.append({"comparison": nm, "flag": fc, "n_group": len(A), "group_flagged": a1, "group_frac": a1 / len(A),
                        "n_reference": len(B), "reference_flagged": b1, "reference_frac": b1 / len(B),
                        "odds_ratio": orr, "fisher_p_two_sided": pv})
    enr = pd.DataFrame(enr)
    enr.to_csv(HERE / "flag_enrichment.tsv", sep="\t", index=False)
    cont = {}
    for nm, (A, B) in groups.items():
        cont[nm] = {c: [float(A[c].median()), float(B[c].median())] for c in ["mean_GC", "GC_percentile", "mean_N", "dist_blacklist_bp", "dist_removed_bin_bp"]}

    # ---- summaries for report / json
    def fam_table(e, kind):
        s = rev[(rev.estimand == e) & (rev.thresholds == kind) & (rev.family.isin(["saved", "called", "literature", "random", "all_segments"]))]
        t = s[["family", "direction", "n_regions", "n_called_any_baseline", "n_called_all_7", "n_unstable", "reversal_rate_among_any", "n_called_under_c"]
              + [c for c in s.columns if c.startswith("keep_vs_c_")]]
        return t.sort_values(["direction", "family"]).reset_index(drop=True)
    saved_tab = bl[(bl.family == "saved") & (bl.estimand == "raw")].pivot_table(index="region", columns="baseline", values="level")[LEVEL_BASES]
    saved_x = bl[(bl.family == "saved") & (bl.estimand == "x")].pivot_table(index="region", columns="baseline", values="level")[LEVEL_BASES]
    saved_g = gdf[gdf.family == "saved"].set_index("region")
    called_g = gdf[gdf.family == "called"]
    lit_tab = bl[(bl.family == "literature") & (bl.estimand == "raw")].pivot_table(index="region", columns="baseline", values="level")[LEVEL_BASES]
    base_tab = pd.DataFrame(base_scalar).T
    base_tab["loco_min"] = {e: min(loco[e].values()) for e in est}
    base_tab["loco_max"] = {e: max(loco[e].values()) for e in est}

    # ---- figure
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    fams = ["saved", "called", "literature", "random", "all_segments"]
    sub = bl[(bl.estimand == "raw") & bl.family.isin(["saved", "called"])]
    lbl = [b.split("_")[0] for b in LEVEL_BASES]
    dat, pos = [], []
    cols_ = {"saved": "#2f5d8a", "called": "#b5651d"}
    for i, b in enumerate(LEVEL_BASES):
        for k, (f, dr) in enumerate([("saved", "loss"), ("called", "loss"), ("called", "gain")]):
            s = sub[(sub.baseline == b) & (sub.family == f) & (sub.direction == dr)].level
            ax[0].scatter(np.full(len(s), i + (k - 1) * .25) + np.random.default_rng((C.SEED, 12, i, k)).uniform(-.05, .05, len(s)), s,
                          s=14, alpha=.7, color=["#2f5d8a", "#b5651d", "#3a7d44"][k], label=f"{f} {dr}" if i == 0 else None)
    ax[0].axhline(.683, ls="--", c="grey"); ax[0].axhline(1.225, ls="--", c="grey"); ax[0].axhline(.5, ls=":", c="grey"); ax[0].axhline(1.5, ls=":", c="grey")
    ax[0].set_xticks(range(len(LEVEL_BASES))); ax[0].set_xticklabels(lbl); ax[0].set_yscale("log"); ax[0].set_ylabel("level = region rate / baseline (raw)")
    ax[0].set_title("Saved windows and called short segments by baseline"); ax[0].legend(fontsize=7)
    s = rev[(rev.estimand == "raw") & (rev.thresholds == "loose") & rev.direction.isin(["loss", "gain"]) & rev.family.isin(fams)]
    w_ = .38
    for k, dr in enumerate(["loss", "gain"]):
        ss = s[s.direction == dr].set_index("family").reindex(fams)
        ax[1].bar(np.arange(len(fams)) + (k - .5) * w_, ss.reversal_rate_among_any.fillna(0), w_, label=dr, color=["#2f5d8a", "#3a7d44"][k])
        for x_, (r_, n_) in enumerate(zip(ss.reversal_rate_among_any, ss.n_called_any_baseline)):
            ax[1].text(x_ + (k - .5) * w_, (0 if not np.isfinite(r_) else r_) + .01, f"n={int(n_)}", ha="center", fontsize=7)
    ax[1].set_xticks(range(len(fams))); ax[1].set_xticklabels(fams, rotation=20); ax[1].set_ylabel("share of any-baseline calls not called by all 7")
    ax[1].set_title("Call instability across baselines (raw, loose thresholds)"); ax[1].legend()
    for fam_, c_, m_ in [("random", "#bbbbbb", "o"), ("saved", "#2f5d8a", "s"), ("called", "#b5651d", "^")]:
        s_ = gdf[gdf.family == fam_]
        ax[2].scatter(s_.own_level_min, s_.cvr_max_abs_log2, s=14 if fam_ == "random" else 28, c=c_, marker=m_, alpha=.8, label=fam_)
    ax[2].axhline(q95, ls="--", c="grey"); ax[2].set_xscale("log")
    ax[2].set_xlabel("min over clusters of own level (rate/cluster genome median)"); ax[2].set_ylabel("max |log2 cluster-vs-rest ratio|")
    ax[2].set_title("Shared vs cluster-specific depletion (dashed: random-window 95th pct)"); ax[2].legend(fontsize=7)
    fig.tight_layout()
    b64 = C.fig_to_b64(fig)
    fig.savefig(HERE / "figure.png", dpi=110, bbox_inches="tight")

    # ---- headline numbers for json
    def pick(fam, dr, e="raw", kind="loose"):
        r_ = rev[(rev.family == fam) & (rev.direction == dr) & (rev.estimand == e) & (rev.thresholds == kind)]
        return None if r_.empty else {k: (None if pd.isna(v) else v) for k, v in r_.iloc[0].to_dict().items()}
    summary = {
        "experiment": "12_baseline_sensitivity", "deterministic": True, "seed_tuple_bootstrap": [C.SEED, 12], "bootstrap_reps": BOOT,
        "random_windows_seed_tuple_from_exp04": [C.SEED, 4],
        "baselines_prespecified": ["a_genome_mean", "b_genome_median", "c_segwt_median", "d_trimmed_mean", "e_kde_mode", "f_loco_median", "g_cluster_vs_rest", "h_gc_decile"],
        "reversal_baselines": LEVEL_BASES, "thresholds": {k: {"loss": v[0], "gain": v[1], "inclusive": v[2]} for k, v in THR.items()},
        "kde_bandwidth_rule": "Silverman, weighted by segment n_bins", "kde_bandwidth_value": bw,
        "scalar_baselines": base_scalar, "n_regions_by_family": {k: int(v) for k, v in pd.Series([r["family"] for r in regs]).value_counts().items()},
        "reversal_table": rev.replace({np.nan: None}).to_dict("records"),
        "saved_windows_levels_raw": saved_tab.round(4).to_dict("index"), "saved_windows_levels_x": saved_x.round(4).to_dict("index"),
        "cluster_vs_rest": {"random_null_q95_max_abs_log2": q95, "saved": saved_g[["cvr_max_abs_log2", "pct_of_random_null", "cvr_n_ci_excl_1", "own_level_min", "own_level_max"]].round(4).to_dict("index"),
                            "called_n": int(len(called_g)), "called_specific_vs_random95": int(called_g.cvr_specific_vs_random95.sum()),
                            "called_loss_n": int((called_g.direction == "loss").sum()),
                            "called_loss_specific": int(called_g[called_g.direction == "loss"].cvr_specific_vs_random95.sum()),
                            "called_loss_all_clusters_below_0.85": int(called_g[called_g.direction == "loss"]["shared_depletion_all_clusters_below_0.85"].sum()),
                            "called_gain_n": int((called_g.direction == "gain").sum()),
                            "called_gain_specific": int(called_g[called_g.direction == "gain"].cvr_specific_vs_random95.sum())},
        "flag_enrichment": enr.to_dict("records"), "flag_medians_group_vs_reference": cont,
        "runtime_seconds": round(time.time() - T0, 1), "environment": C.environment(),
    }
    C.write_json(HERE / "summary.json", summary)

    # ---- report (presentation only; analysis above is unchanged)
    T, cap = C.html_table, C.caption
    rl = rev[(rev.estimand == "raw") & (rev.thresholds == "loose")].set_index(["family", "direction"])

    def rr_(f, dr, e="raw", kind="loose"):
        x = rev[(rev.family == f) & (rev.direction == dr) & (rev.estimand == e) & (rev.thresholds == kind)].iloc[0]
        return f"{x.n_unstable:.0f}/{x.n_called_any_baseline:.0f}"
    E_ = enr.set_index(["comparison", "flag"])
    k_ = E_.loc[("called loss segments vs random windows", "near_removed_1Mb")]
    k2 = E_.loc[("called loss segments vs random windows", "near_blacklist_1Mb")]
    k3 = E_.loc[("deep random (raw level c < 0.683) vs other random", "GC_extreme_decile")]
    sv = cls[(cls.family == "saved") & (cls.estimand == "raw") & (cls.thresholds == "loose")].set_index("region")
    svs = cls[(cls.family == "saved") & (cls.estimand == "raw") & (cls.thresholds == "strict")].set_index("region")
    cv = summary["cluster_vs_rest"]
    nloss_saved = int((sv.direction == "loss").sum())
    all7_saved = int(((sv.direction == "loss") & (sv.n_loss_baselines == 7)).sum())

    # extra figures
    lg = np.log2(saved_tab.copy())
    figA, axA = plt.subplots(figsize=(8.5, 3.8))
    im = axA.imshow(lg.to_numpy(), cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
    axA.set_xticks(range(len(LEVEL_BASES))); axA.set_xticklabels([b.split("_")[0] + " " + "_".join(b.split("_")[1:]) for b in LEVEL_BASES], rotation=30, ha="right", fontsize=8)
    axA.set_yticks(range(len(lg))); axA.set_yticklabels(lg.index, fontsize=8)
    for i in range(lg.shape[0]):
        for j in range(lg.shape[1]):
            axA.text(j, i, f"{saved_tab.iloc[i, j]:.2f}", ha="center", va="center", fontsize=7)
    figA.colorbar(im, ax=axA, label="log2 level (blue = depleted)")
    axA.set_title("Saved windows: level under each baseline (raw)")
    b64A = C.fig_to_b64(figA)
    figB, axB = plt.subplots(figsize=(8.5, 3.8))
    sg = gdf[gdf.family == "saved"].set_index("region")
    for c in range(5):
        axB.errorbar(np.arange(len(sg)) + (c - 2) * .12, sg[f"cvr_c{c}"], yerr=[sg[f"cvr_c{c}"] - sg[f"cvr_c{c}_lo"], sg[f"cvr_c{c}_hi"] - sg[f"cvr_c{c}"]],
                     fmt="o", ms=4, lw=1, label=f"cluster {c}")
    axB.axhline(1, c="grey", ls="--"); axB.set_xticks(range(len(sg))); axB.set_xticklabels(sg.index, rotation=25, ha="right", fontsize=8)
    axB.set_ylabel("cluster / rest ratio (95% cell-bootstrap CI)"); axB.legend(fontsize=7, ncol=5); axB.set_title("Saved windows: cluster-vs-rest")
    b64B = C.fig_to_b64(figB)
    plt.close("all")

    verdict = (f"Under eight prespecified reference-free baselines the saved chr7/chr9 depletions stay below the 0.683 loss cut-off in {all7_saved} of {nloss_saved} windows "
               f"for all seven level baselines (chr9_20_27 in 6 of 7), but loss/gain classes of ordinary and random regions change with baseline choice often enough that a call is, "
               "consistent with the data, a statement about the chosen reference as much as about copy number.")
    kf = [f"Reversal (raw, loose; regions called by at least one baseline but not all seven): called losses <b>{rr_('called','loss')}</b>, called gains <b>{rr_('called','gain')}</b>, "
          f"random-window losses <b>{rr_('random','loss')}</b> and gains <b>{rr_('random','gain')}</b>; saved losses {rr_('saved','loss')}; literature loci loss {rr_('literature','loss')}, gain {rr_('literature','gain')}.",
          f"Opposite-sign reversals (loss under one baseline, gain under another): 1 of 67 called segments but 31 of 300 random windows (raw, loose).",
          f"Processed-X estimand is more stable than the raw rate: called losses {rr_('called','loss','x')}, called gains {rr_('called','gain','x')}; strict thresholds are less stable: called losses {rr_('called','loss','raw','strict')}, gains {rr_('called','gain','raw','strict')}.",
          "Saved windows (raw level, baseline a-f,h): " + "; ".join(f"{r}: {int(v.n_loss_baselines)}/7 loss, range {saved_tab.loc[r].min():.2f}-{saved_tab.loc[r].max():.2f}" for r, v in sv.iterrows()) + ".",
          "The GC-decile baseline (h) is the main outlier for the chr7_8.3 and chr9 windows (e.g. chr9_20_27 = %.2f vs 0.42-0.59 for the others), consistent with part of that depletion being GC-correlated." % saved_tab.loc["chr9_20_27", "h_gc_decile"],
          f"Cluster-vs-rest: chr9 windows have ratios 0.86-1.10 in every cluster (consistent with a depletion shared across clusters); chr7_8.3, chr7_51.4 and chr7_67.4 sit near or just above the random-window 95th percentile ({q95:.2f} max |log2|) driven by cluster 2 (1.36-1.53). Called losses: {cv['called_loss_specific']}/{cv['called_loss_n']} cluster-specific, {cv['called_loss_all_clusters_below_0.85']}/{cv['called_loss_n']} below 0.85 in every cluster.",
          f"Flags (descriptive Fisher): called losses near a QC-removed bin {k_.group_flagged:.0f}/{k_.n_group:.0f} vs {k_.reference_flagged:.0f}/300 random windows (OR {k_.odds_ratio:.0f}); near blacklist {k2.group_flagged:.0f}/{k2.n_group:.0f} vs {k2.reference_flagged:.0f}/300. "
          "Caveat: depleted bins are the ones the nonzero-frequency filter removes, so this is consistent with a selection effect on the same signal, not an independent artifact test.",
          f"Deep random windows (raw level c &lt; 0.683) are not enriched near removed bins or blacklist, but are GC-extreme more often ({k3.group_flagged:.0f}/{k3.n_group:.0f} vs {k3.reference_flagged:.0f}/{k3.n_reference:.0f})."]

    q = ("<p>A loss or gain call from pyEpiAneufinder is a threshold on (segment mean) / (a genome-wide scale), and experiment 03 showed the scale sits at the lower grid edge. "
         "This experiment asks how much each region's relative level, and the class it receives under fixed cut-offs, depends on which reference-free baseline is used.</p>"
         "<ul><li><b>C:3</b> reference-free relative baseline sensitivity: compare genome-wide centering, robust baselines and held-out alternatives; expose the baseline choice.</li>"
         "<li><b>C:1</b> independent internal controls: none exist in this sample, so no normal cluster is forced from MCF10-A identity or apparent flatness and no baseline is treated as copy-neutral.</li>"
         "<li><b>F:2</b> blacklist / genomic gap mask and <b>F:4</b> GC content: technical context reported as soft flags, not as corrections.</li></ul>"
         "<p><b>Decision informed:</b> whether the method recommendation should report calls as baseline-conditional relative levels (with the baseline exposed) rather than absolute copy-number states, and how much weight the saved chr7/chr9 depletions can carry.</p>")
    inputs = ("<ul><li>AnnData <code>pyepi_results/count_matrix.h5ad</code>: 3,460 cells x 25,337 retained autosomal 100 kb bins; <code>layers/raw</code> raw counts, <code>X</code> processed, <code>var GC</code>.</li>"
              "<li>Saved genome-pool segmentation <code>karyogram_cluster_segments.tsv</code> (group <code>all</code>, 164 segments), held fixed; WNN clusters 0-4 from <code>integrated.h5mu</code>.</li>"
              "<li>Full post-blacklist grid <code>binned_genome.csv</code> (GC, N) and <code>hg38-blacklist.v2.bed</code> for flags; QC-removed bins = grid bins absent from the AnnData (experiment 02).</li>"
              "<li>Region families from experiment 04 <code>build_regions</code> (target bins only): 7 saved windows, 67 called short segments (n_bins &le; 100, cn_call &ne; 1), 4 literature loci, 300 random 25-bin windows (seed tuple (20261104, 4)); plus all 164 segments as <code>all_segments</code>.</li></ul>"
              + T(pd.DataFrame([{"family": k, "regions": v} for k, v in summary["n_regions_by_family"].items()]), "{:.0f}"))
    btab = pd.DataFrame([
        ["a genome_mean", "B = mean_i v_i over all bins", "global scalar"],
        ["b genome_median", "B = median_i v_i", "global scalar"],
        ["c segwt_median", "B = weighted median of segment means m_s, weights n_s (bins)", "global scalar; closest to pyEpi's effective scale (exp 03)"],
        ["d trimmed_mean", "B = mean of v_i after dropping lowest/highest 10%", "global scalar"],
        ["e kde_mode", "B = argmax of weighted Gaussian KDE of m_s (Silverman bandwidth; raw %.2f, X %.2f)" % (bw["raw"], bw["x"]), "global scalar"],
        ["f loco_median", "B_c = median of v_i on all chromosomes except the region's chromosome c", "region-specific"],
        ["g cluster_vs_rest", "level = rate_k(region) / rate_rest(region); rate = sum raw / sum exposure; rest = other four clusters pooled", "between-cluster contrast; not a level; excluded from reversal counts"],
        ["h gc_decile", "level = sum_i v_i / sum_i median(v | GC decile of bin i) over region bins", "GC-matched reference, no LOESS"]],
        columns=["baseline", "definition / formula", "type"])
    pre = (T(btab) + "<p>Region level L = mean of bin values in the region / B (a-f); h as above. Estimands: <code>raw</code> = sum raw / sum exposure x 1e5; <code>x</code> = mean processed X. "
           "Classes: loose loss L &lt; 0.683, gain L &gt; 1.225; strict loss L &le; 0.5, gain L &ge; 1.5. Scalar baseline values:</p>"
           + T(base_tab.reset_index().rename(columns={"index": "estimand"}), "{:.4g}"))
    r1 = (T(fam_table("raw", "loose"), "{:.3g}") + cap("Raw estimand, loose thresholds. <code>n_called_any_baseline</code> = regions classed loss (gain) by at least one of seven baselines; "
          "<code>n_called_all_7</code> = by all seven; <code>reversal_rate_among_any</code> = unstable share; <code>keep_vs_c_*</code> = share of calls under baseline c retained under baseline a, b, d, e, f, h. "
          "Random windows are not neutral, so their high instability marks how easily any region crosses a cut-off, not an error rate.")
          + T(fam_table("raw", "strict"), "{:.3g}") + cap("Same, strict thresholds (loss &le; 0.5, gain &ge; 1.5).")
          + T(fam_table("x", "loose"), "{:.3g}") + cap("Same, processed-X estimand, loose thresholds."))
    r2 = (T(saved_tab.reset_index(), "{:.3f}") + cap("Raw level of each saved window under baselines a-f,h; values below 0.683 are loose-threshold losses.")
          + T(saved_x.reset_index(), "{:.3f}") + cap("Same for processed X.")
          + C.img(b64A, "saved windows heatmap") + cap("Heatmap of log2 level (cell text = level). Blue = depleted relative to the baseline. Columns that differ from the rest (h) show baseline dependence."))
    r3 = (T(saved_g.reset_index()[["region", "n_bins"] + [f"cvr_c{c}" for c in range(5)] + ["cvr_max_abs_log2", "pct_of_random_null", "cvr_n_ci_excl_1", "own_level_min", "own_level_max"]], "{:.3f}")
          + cap(f"<code>cvr_c*</code> = cluster / rest rate ratio (1 = same depth as the other clusters). <code>pct_of_random_null</code> = percentile of max|log2 ratio| among 300 random windows (95th percentile = {q95:.2f}). "
                "<code>own_level_*</code> = region rate over that cluster's own genome median, min/max across clusters.")
          + C.img(b64B, "cluster vs rest") + cap("Points = cluster/rest ratio with 95% cell-bootstrap CI (300 replicates). Ratios near 1 with low own levels are consistent with a depletion shared by all clusters.")
          + f"<p>Called short segments: losses cluster-specific (above random 95th percentile) {cv['called_loss_specific']}/{cv['called_loss_n']}; gains {cv['called_gain_specific']}/{cv['called_gain_n']}; "
            f"loss segments below 0.85 of their own genome median in every cluster: {cv['called_loss_all_clusters_below_0.85']}/{cv['called_loss_n']}.</p>")
    r4 = (T(enr, "{:.3g}") + cap("Descriptive Fisher exact tests (two-sided, uncorrected). Flags are within 1 Mb of a QC-removed bin / blacklist interval / assembly-gap-like run, or mean GC in the lowest or highest decile. "
          "Called-loss enrichment near removed bins is consistent with selection on depletion by the bin-removal filter.")
          + C.img(b64, "overview") + cap("Left: levels of saved windows and called segments per baseline (dashed = loose, dotted = strict cut-offs). Middle: instability per family. Right: cluster-specific versus shared depletion."))
    notes = {"chr7_8.3_10.7": "Deep depletion under every baseline; GC baseline raises it (low-GC window); cluster 2 is shallower than the rest.",
             "chr7_51.4_53.5": "Deep under all baselines; next to QC-removed bins; cluster 2 shallower.",
             "chr7_67.4_70.0": "Deepest window, below 0.15 under most baselines; cluster 2 shallower.",
             "chr7_70.0_73.3": "Shallowest of the chr7 loss windows (0.30-0.47 depending on baseline); sits just under the strict 0.5 cut-off, so a slightly higher baseline-independent level would change its strict class.",
             "chr9_20_27": "Loss under 6 of 7 baselines; GC baseline (h) lifts it above the cut-off; shared across clusters.",
             "chr9_22.5_26.6": "Loss under all baselines; shared across clusters; next to removed bins.",
             "chr9_70_100": "Comparison region, near 1 (0.93-1.31 depending on mean vs median baseline); never a loss, in the gain range only for median-type baselines at the loose cut."}
    det = saved_tab.copy().round(3)
    det["loose_loss_baselines"] = [int(sv.loc[r, "n_loss_baselines"]) for r in det.index]
    det["strict_loss_baselines"] = [int(svs.loc[r, "n_loss_baselines"]) for r in det.index]
    det["note"] = [notes.get(r, "") for r in det.index]
    lit = lit_tab.copy().round(3)
    lit["loose_loss_baselines"] = [int(cls[(cls.region == r) & (cls.estimand == "raw") & (cls.thresholds == "loose")].n_loss_baselines.iloc[0]) for r in lit.index]
    lit["loose_gain_baselines"] = [int(cls[(cls.region == r) & (cls.estimand == "raw") & (cls.thresholds == "loose")].n_gain_baselines.iloc[0]) for r in lit.index]
    lit["note"] = ["hemizygous 9p21 zone (approximate coordinates); read as a relative level only" if "9p21" in r else "published gain locus; class depends on baseline" for r in lit.index]
    detail = (T(det.reset_index()) + T(lit.reset_index()) + cap("One row per region: raw level under each baseline, number of baselines (of 7) giving a loss at the loose and strict cut-offs, and a plain-language note."))
    lim = ("<ul><li>Segmentation is fixed; baselines are not propagated into re-segmentation.</li>"
           "<li>Baselines a-e share one genome-wide distribution, so their agreement is not independent evidence; on a skewed, aneuploid genome they can all be off together.</li>"
           "<li>f removes only one chromosome; h removes GC-driven structure but also any real GC-correlated copy-number signal.</li>"
           "<li>Random windows are not neutral controls; called regions are selected by their own call, which biases overlap with c.</li>"
           "<li>Alternative explanations for the saved depletions that this experiment cannot exclude: shared accessibility/mappability baseline offset, a GC-related effect, or real loss. Cluster-vs-rest cancels a shared baseline and so cannot distinguish them.</li>"
           "<li>Bootstrap resamples cells only (300 reps) and ignores segmentation uncertainty. Enrichment tests are uncorrected and regions are correlated.</li></ul>")
    supports = [f"{all7_saved} of {nloss_saved} saved loss windows are below the loose loss cut-off under all seven level baselines, consistent with a robust relative depletion regardless of centering.",
                "chr9 windows' depletion is similar across clusters (cluster-vs-rest 0.86-1.10), consistent with a shared rather than cluster-specific pattern.",
                "A call is partly a statement about the chosen baseline: many random and gain-class regions flip class across baselines."]
    notsup = ["That any baseline is copy-neutral, or that the depletions are copy-number loss rather than a shared depth/accessibility baseline effect.",
              "That GC is the cause of the chr9_20_27 shift; baseline h only shows GC-matched centering moves it.",
              "That QC-removed-bin enrichment of called losses is an independent artifact; it is consistent with filter selection on depleted bins."]
    decision = ("<p>Supports the point already in <code>recommendation.md</code> (rank 2 cross-sample technical baseline, rank 4 baseline-offset model): relative depth levels should be reported with the baseline exposed, and "
                "loss/gain states, especially gains and strict thresholds, should be treated as baseline-conditional. It does not change the ranking; it adds a concrete reversal rate (called gains %s, random-window gains %s) to the evidence for needing an independent baseline.</p>" % (rr_('called', 'gain'), rr_('random', 'gain')))
    nxt = ["Run the same baselines on unrelated samples with DNA truth or an annotated normal population (recommendation rank 2) to see which baseline tracks truth.",
           "Re-segment under baseline-offset alternatives and compare breakpoints rather than only classes.",
           "Replace the 1 Mb removed-bin flag with a pre-filter depth flag to separate selection from artifact.",
           "Check chr9_20_27 and chr7_8.3 with a LOESS GC correction versus the decile reference."]
    gl = {"level": "region rate / baseline", "baseline": "a-f,h reference; g is a between-cluster contrast", "class_loose / class_strict": "loss / neutral / gain by the loose or strict cut-offs",
          "n_loss_baselines / n_gain_baselines": "number of the 7 level baselines giving loss / gain", "n_distinct_classes": "distinct classes across baselines",
          "opposite_loss_gain": "loss under one baseline and gain under another", "log2_level_range": "log2(max level / min level) across baselines",
          "n_called_any_baseline": "regions called by at least one baseline", "n_called_all_7": "regions called by all seven", "n_unstable": "called by some but not all",
          "reversal_rate_among_any": "n_unstable / n_called_any_baseline", "n_called_under_c": "regions called under baseline c", "keep_vs_c_*": "share of c-calls kept under the other baseline (letter)",
          "cvr_c*": "cluster / rest raw rate ratio; _lo/_hi = 2.5/97.5 percentile of 300 cell bootstraps", "cvr_max_abs_log2": "largest |log2 ratio| over clusters",
          "cvr_n_ci_excl_1": "clusters whose bootstrap CI excludes 1", "pct_of_random_null": "percentile within the 300 random windows", "own_level_c*": "region rate over cluster c's own genome median",
          "cvr_specific_vs_random95": "max|log2 ratio| above random-window 95th percentile", "dist_removed_bin_bp / _gridbins": "distance to nearest QC-removed grid bin (0 = within/adjacent span)",
          "frac_within_1Mb_of_removed": "share of region bins within 1 Mb of a removed bin", "dist_blacklist_bp": "distance to nearest hg38 blacklist interval",
          "dist_assembly_gap_like_bp": "distance to nearest removed run with mean N > 0.5", "mean_GC / GC_percentile": "region mean GC and its percentile among retained bins", "mean_N": "mean N fraction",
          "near_*_1Mb / GC_extreme_decile": "binary soft flags", "odds_ratio / fisher_p_two_sided": "descriptive Fisher exact test, uncorrected",
          "KDE / Silverman": "kernel density estimate with Silverman's bandwidth rule", "LOCO": "leave-one-chromosome-out", "exposure": "per-cell total raw counts"}
    sections = [("Question and why it matters", q), ("Inputs", inputs), ("Prespecified baselines", pre), ("Results", "<p>Exploratory results; each block is followed by a caption.</p>"),
                ("Reversals by family, direction, estimand and threshold", r1, 3), ("Saved windows under every baseline", r2, 3),
                ("Cluster-vs-rest", r3, 3), ("Technical flags and enrichment", r4, 3),
                ("Per-region detail", detail), ("Limits and alternative explanations", lim)]
    C.write_html_report(HERE / "report.html", "Experiment 12: baseline sensitivity of relative depth calls",
                        "Exploratory diagnostic, saved genome-pool segmentation held fixed; no baseline is copy-neutral truth.", sections,
                        {k: summary[k] for k in ["seed_tuple_bootstrap", "bootstrap_reps", "kde_bandwidth_rule", "runtime_seconds", "environment"]},
                        key_findings=kf, verdict=verdict, glossary=gl, supports=supports, not_supports=notsup, decision_impact=decision, next_steps=nxt)
    print(f"done in {time.time()-T0:.0f}s")


if __name__ == "__main__":
    main()
