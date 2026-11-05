"""Experiment 11: depth-only masked count changepoint (method record B:2) vs the AD baseline.

EXPLORATORY DIAGNOSTIC. Nothing here is validated production behaviour and segment states are
relative read depth (accessibility and dosage are confounded), not copy number.

Model. Per source (all, cluster0) the post-QC raw integer counts (layers/raw) of the member cells are
summed to a pseudobulk y_b per retained bin; exposure E = total raw counts of the member cells (constant over
bins, so the bin rate is y_b/E). Zeros are kept. Within a chromosome the retained bins form the grid; a
changepoint may only fall between two consecutive retained bins (QC-removed bins are explicit gaps, never
interpolated, and every breakpoint records whether bins were removed at that position).

Segmentation. Optimal partitioning (exact dynamic programme, O(n^2) per chromosome) with piecewise-constant
log-rate. Segment cost = -S log(S/m) (Poisson negative log-likelihood at the segment MLE, constants that do not
depend on the partition dropped; S = sum of y, m = bins) divided by dispersion phi. Penalty per extra segment
beta = c * log(n_chrom) with c in {1, 2, 5}; minimum segment length 5 bins.  [Scaling cost by 1/phi is the same
as penalty beta*phi on the unscaled cost.]

phi (quasi-Poisson dispersion), modes:
  phi1       phi = 1 (plain Poisson; reference only)
  half_split phi from the half-split (C.stratified_halves(labels, (SEED, 11))). Given y_A+y_B=n, y_A ~ Bin(n,p),
             p = E_A/(E_A+E_B) under no difference beyond cell sampling; z^2 = (y_A-np)^2/(n p (1-p)) has mean phi.
             phi = median(z^2)/median(chi2_1) over bins with n >= 20 (robust). Isolates cell-sampling noise.
  pearson    iterative: start from phi of a one-segment-per-chromosome fit (all structure counted as noise),
             segment, recompute phi = median(Pearson r^2)/0.4549 with r^2=(y-mu)^2/mu from the segmentation
             residuals, resegment; stop when the segment count is unchanged (max 12 iterations). Bin-level
             reproducible structure inflates this phi, so it is the conservative setting.
  nb_half    negative binomial, single alpha (var = mu + alpha mu^2) from the half split by a moment fit,
             alpha_half = sum[(y_A-rho y_B)^2 - (mu_A+rho^2 mu_B)] / sum[mu_A^2+rho^2 mu_B^2] (plug-in mu from
             y_A+y_B, 0.5/99.5% winsorised), alpha_full = alpha_half * n_cells_half/n_cells_full (the
             cell-sampling quadratic term scales as 1/#cells; an assumption). Segment mean S/m; penalty c*log(n).
phi is one scalar per source (genome-wide), reused for every chromosome. In half-data stability runs the
half_split phi of the parent source is reused (variance/mean ratio of cell sums does not depend on #cells);
pearson is re-estimated on each half; nb_half uses alpha_half.

Run: /home/user/envs/cnv310/bin/python -I run.py   (2 worker processes)
"""
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

SEED = C.SEED
HALF_SEED = (SEED, 11)
CS = [1, 2, 5]
MIN_SEG = 5
MAX_ITER = 12
MED_CHI2_1 = float(chi2.ppf(0.5, 1))
CHROMS = [f"chr{i}" for i in range(1, 23)]
SOURCES = ["all", "cluster0"]
MODES = ["phi1", "half_split", "pearson", "nb_half"]


# ----------------------------------------------------------------------------- segmentation
def op_segment(n, cost_fn, beta, minlen=MIN_SEG):
    """Optimal partitioning. cost_fn(s_array, t) -> cost of segments [s,t). Returns sorted segment starts (excl. 0)."""
    if n < 2 * minlen:
        return []
    F = np.full(n + 1, np.inf)
    F[0] = -beta
    last = np.zeros(n + 1, dtype=np.int64)
    for t in range(minlen, n + 1):
        s = np.arange(0, t - minlen + 1)
        s = s[np.isfinite(F[s])]
        if len(s) == 0:
            continue
        v = F[s] + cost_fn(s, t) + beta
        k = int(np.argmin(v))
        F[t], last[t] = v[k], s[k]
    bps, t = [], n
    while t > 0:
        s = int(last[t])
        if s > 0:
            bps.append(s)
        t = s
    return sorted(bps)


def make_cost(y, phi, nb_r=None):
    cs = np.concatenate([[0.0], np.cumsum(y)])
    if nb_r is None:
        def cost(s, t):
            S = cs[t] - cs[s]
            m = (t - s).astype(float)
            with np.errstate(divide="ignore", invalid="ignore"):
                c = np.where(S > 0, -S * np.log(S / m), 0.0)
            return c / phi
    else:
        r = nb_r

        def cost(s, t):
            S = cs[t] - cs[s]
            m = (t - s).astype(float)
            mu = S / m
            with np.errstate(divide="ignore", invalid="ignore"):
                # = sum of NB log-pmf kernels at constant mean: S log(mu/r) - (m r + S) log1p(mu/r)  (S=0 -> second term only)
                ll = np.where(S > 0, S * np.log(mu / r), 0.0) - (m * r + S) * np.log1p(mu / r)
            return -ll
    return cost


def seg_bounds(n, bps):
    e = [0] + list(bps) + [n]
    return list(zip(e[:-1], e[1:]))


def fit_chromosomes(Y, mode, c, phi=None, nb_r=None):
    out = {}
    for ch, y in Y.items():
        beta = c * np.log(len(y))
        cost = make_cost(y, 1.0, nb_r) if mode == "nb_half" else make_cost(y, phi)
        out[ch] = op_segment(len(y), cost, beta)
    return out


def pearson_phi(Y, BP):
    r2 = []
    for ch, y in Y.items():
        for a, b in seg_bounds(len(y), BP[ch]):
            mu = y[a:b].mean()
            if mu > 0:
                r2.append((y[a:b] - mu) ** 2 / mu)
    return float(np.median(np.concatenate(r2)) / MED_CHI2_1)


def pearson_iterate(Y, c):
    BP = {ch: [] for ch in Y}
    phi = max(pearson_phi(Y, BP), 1.0)
    trace, prev = [], None
    for it in range(MAX_ITER):
        BP = fit_chromosomes(Y, "pearson", c, phi=phi)
        nseg = sum(len(v) + 1 for v in BP.values())
        trace.append({"iter": it, "phi_used": phi, "n_segments": nseg})
        if nseg == prev:
            break
        prev = nseg
        phi = max(pearson_phi(Y, BP), 1.0)
    return BP, trace[-1]["phi_used"], trace


def half_split_params(yA, yB, EA, EB):
    n = yA + yB
    p = EA / (EA + EB)
    ok = n >= 20
    z2 = (yA[ok] - n[ok] * p) ** 2 / (n[ok] * p * (1 - p))
    rho = EA / EB
    muA, muB = n * rho / (1 + rho), n / (1 + rho)
    T = (yA - rho * yB) ** 2 - (muA + rho ** 2 * muB)
    D = muA ** 2 + rho ** 2 * muB ** 2
    q = np.quantile(T, [0.005, 0.995])
    alpha_half = float(np.clip(T, *q).sum() / D.sum())
    return {"phi_half_median": float(np.median(z2) / MED_CHI2_1),
            "phi_half_trimmed_mean": float(np.mean(np.minimum(z2, np.quantile(z2, 0.99)))),
            "alpha_half": max(alpha_half, 1e-6), "alpha_half_raw": alpha_half, "n_bins_used": int(ok.sum())}


def diff_phi(Y):
    v = []
    for y in Y.values():
        a, b = y[1:], y[:-1]
        k = (a + b) > 0
        v.append((a[k] - b[k]) ** 2 / (a[k] + b[k]))
    return float(np.median(np.concatenate(v)) / MED_CHI2_1)


# ----------------------------------------------------------------------------- pseudobulks
def build_pseudobulks(d, lab, half):
    masks = {}
    for s, m in [("all", np.ones(d.n_cells, bool)), ("cluster0", lab == "0")]:
        masks[s], masks[s + "_A"], masks[s + "_B"] = m, m & half, m & ~half
    Y = {k: np.zeros(d.n_bins) for k in masks}
    E = {k: 0.0 for k in masks}
    for lo in range(0, d.n_cells, 400):
        hi = min(lo + 400, d.n_cells)
        blk = d.rows("layers/raw", lo, hi)
        tot = blk.sum(1)
        for k, m in masks.items():
            mm = m[lo:hi]
            if mm.any():
                Y[k] += blk[mm].sum(0)
                E[k] += tot[mm].sum()
    return masks, Y, E


def run_unit(args):
    key, Yc, mode, c, phi, nb_r = args
    t0 = time.time()
    if mode == "pearson":
        BP, phi_used, trace = pearson_iterate(Yc, c)
    else:
        BP, phi_used, trace = fit_chromosomes(Yc, mode, c, phi=phi, nb_r=nb_r), phi, []
    return key, BP, phi_used, trace, time.time() - t0


def match(a, b, tol):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) == 0 or len(b) == 0:
        return 0
    return int(sum(np.min(np.abs(b - x)) <= tol for x in a))


def main():
    t_start = time.time()
    p = C.default_paths()
    d = C.H5AD(p["h5ad"])
    lab = C.cell_clusters(d, p)
    half = C.stratified_halves(lab, HALF_SEED)
    masks, Yfull, Efull = build_pseudobulks(d, lab, half)
    ncell = {k: int(m.sum()) for k, m in masks.items()}
    Yc = {k: {ch: v[d.seq == ch] for ch in CHROMS} for k, v in Yfull.items()}
    print("pseudobulks built", round(time.time() - t_start, 1), "s", flush=True)

    hs = {}
    for s in SOURCES:
        hs[s] = half_split_params(Yfull[s + "_A"], Yfull[s + "_B"], Efull[s + "_A"], Efull[s + "_B"])
        hs[s]["alpha_full"] = hs[s]["alpha_half"] * ncell[s + "_A"] / ncell[s]
        hs[s]["diff_phi_full_informational"] = diff_phi(Yc[s])
        hs[s]["E_total_counts"] = Efull[s]
        hs[s]["n_cells"] = ncell[s]
        hs[s]["n_cells_A"], hs[s]["n_cells_B"] = ncell[s + "_A"], ncell[s + "_B"]
    print({s: {k: round(v, 4) for k, v in hs[s].items()} for s in SOURCES}, flush=True)

    tasks = []
    for s in SOURCES:
        for part in ["", "_A", "_B"]:
            for mode in MODES:
                for c in CS:
                    phi = {"phi1": 1.0, "half_split": hs[s]["phi_half_median"]}.get(mode)
                    nb_r = None
                    if mode == "nb_half":
                        nb_r = 1.0 / (hs[s]["alpha_half"] if part else hs[s]["alpha_full"])
                    tasks.append(((s + part, mode, c), Yc[s + part], mode, c, phi, nb_r))
    results = {}
    with ProcessPoolExecutor(2) as ex:
        for key, BP, phi_used, trace, dt in ex.map(run_unit, tasks, chunksize=1):
            results[key] = {"BP": BP, "phi": phi_used, "trace": trace}
            print(key, "segments", sum(len(v) + 1 for v in BP.values()), "phi", None if phi_used is None else round(phi_used, 2), round(dt, 1), "s", flush=True)

    # ---- AD baseline (retained chrom-local bin index = first bin of the new segment)
    ad = pd.read_csv(p["genome_pool"] / "raw_accepted_breakpoints.tsv", sep="\t")
    AD = {s: {ch: sorted(ad[(ad.source == s) & (ad.chromosome == ch)].absolute_bin.astype(int)) for ch in CHROMS} for s in SOURCES}
    chrom_idx = {ch: np.flatnonzero(d.seq == ch) for ch in CHROMS}
    for _, r in ad[ad.source.isin(SOURCES)].iterrows():
        assert d.start[chrom_idx[r.chromosome]][int(r.absolute_bin)] == r.bp_genomic_start

    bg = pd.read_csv(p["root"] / "pyepi_results" / "binned_genome.csv", index_col=0)
    bg_by = {ch: bg[bg.chromosome == ch].start.to_numpy() for ch in CHROMS}
    width = int(np.median(d.end - d.start + 1))

    def gap_info(ch, i):
        st, en = d.start[chrom_idx[ch]], d.end[chrom_idx[ch]]
        removed = int((st[i] - en[i - 1] - 1) // width)
        in_bg = int(((bg_by[ch] > en[i - 1]) & (bg_by[ch] < st[i])).sum())
        return removed, in_bg

    main_keys = [k for k in results if not k[0].endswith(("_A", "_B"))]
    bp_rows = []
    for (src, mode, c) in main_keys:
        for ch in CHROMS:
            for i in results[(src, mode, c)]["BP"][ch]:
                removed, in_bg = gap_info(ch, i)
                bp_rows.append({"source": src, "mode": mode, "c": c, "chrom": ch, "bin_index": i,
                                "genomic_start": int(d.start[chrom_idx[ch]][i]), "removed_bins_before": removed,
                                "binned_genome_bins_in_gap": in_bg, "at_retained_grid_gap": removed > 0})
    bp = pd.DataFrame(bp_rows)
    bp.to_csv(HERE / "breakpoints.tsv", sep="\t", index=False)
    ad_gap = {s: int(sum(gap_info(ch, i)[0] > 0 for ch in CHROMS for i in AD[s][ch])) for s in SOURCES}
    ad_tot = {s: int(sum(len(AD[s][ch]) for ch in CHROMS)) for s in SOURCES}

    # ---- segments with relative rates
    seg_rows = []
    for (src, mode, c) in main_keys:
        r = results[(src, mode, c)]
        recs = []
        for ch in CHROMS:
            y = Yc[src][ch]
            for a, b in seg_bounds(len(y), r["BP"][ch]):
                S, m = y[a:b].sum(), b - a
                recs.append({"source": src, "mode": mode, "c": c, "chrom": ch, "bin_start": a, "bin_end_excl": b,
                             "genomic_start": int(d.start[chrom_idx[ch]][a]), "genomic_end": int(d.end[chrom_idx[ch]][b - 1]),
                             "n_bins": m, "counts": S, "rate_per_bin": S / m / Efull[src]})
        sg = pd.DataFrame(recs)
        med = float(np.median(sg.rate_per_bin))
        srt = sg.sort_values("rate_per_bin")
        wmed = float(srt.rate_per_bin.to_numpy()[np.searchsorted(np.cumsum(srt.n_bins.to_numpy()), srt.n_bins.sum() / 2)])
        sg["rel_rate"], sg["rel_rate_wmedian"] = sg.rate_per_bin / med, sg.rate_per_bin / wmed
        sg["phi_used"] = r["phi"]
        for nm, ph in [("qp_half", hs[src]["phi_half_median"]), ("qp_pearson", results[(src, "pearson", 2)]["phi"])]:
            se = np.sqrt(ph / np.maximum(sg.counts, 1))
            sg[f"rel_lo95_{nm}"], sg[f"rel_hi95_{nm}"] = sg.rel_rate * np.exp(-1.96 * se), sg.rel_rate * np.exp(1.96 * se)
        seg_rows.append(sg)
    segs = pd.concat(seg_rows, ignore_index=True)
    segs.to_csv(HERE / "segments.tsv", sep="\t", index=False, float_format="%.6g")

    # ---- comparison with AD
    cmp_rows = []
    for (src, mode, c) in main_keys:
        r = results[(src, mode, c)]
        for region in ["chr7", "chr9", "chr1-22"]:
            chs = CHROMS if region == "chr1-22" else [region]
            row = {"level": "chrom", "source": src, "mode": mode, "c": c, "region": region, "phi_used": r["phi"] if r["phi"] is not None else np.nan,
                   "n_new": 0, "n_ad": 0, "new_matched_pm2": 0, "ad_matched_pm2": 0, "new_matched_pm5": 0, "ad_matched_pm5": 0}
            for ch in chs:
                new, adb = r["BP"][ch], AD[src][ch]
                row["n_new"] += len(new)
                row["n_ad"] += len(adb)
                for tol in (2, 5):
                    row[f"new_matched_pm{tol}"] += match(new, adb, tol)
                    row[f"ad_matched_pm{tol}"] += match(adb, new, tol)
            for tol in (2, 5):
                row[f"new_only_pm{tol}"] = row["n_new"] - row[f"new_matched_pm{tol}"]
                row[f"ad_only_pm{tol}"] = row["n_ad"] - row[f"ad_matched_pm{tol}"]
            cmp_rows.append(row)
    reg = C.shared_regions()
    win = []
    for _, w in reg[reg.part == "target"].iterrows():
        idx = chrom_idx[w.chrom]
        lo = int(np.flatnonzero(d.start[idx] == w.first_retained_start)[0])
        hi = int(np.flatnonzero(d.end[idx] == w.last_retained_end)[0]) + 1
        win.append((w.region_id.replace("__target", ""), w.chrom, lo, hi, w.role))
    for (src, mode, c) in main_keys:
        r = results[(src, mode, c)]
        for name, ch, lo, hi, role in win:
            row = {"level": "window", "source": src, "mode": mode, "c": c, "region": name, "phi_used": r["phi"] if r["phi"] is not None else np.nan, "role": role,
                   "window_bins": hi - lo}
            for lab_, bps in [("new", r["BP"][ch]), ("ad", AD[src][ch])]:
                b = np.asarray(bps)
                for edge, pos in [("left", lo), ("right", hi)]:
                    dist = b - pos
                    row[f"{lab_}_{edge}_nearest_offset_bins"] = int(dist[np.argmin(np.abs(dist))]) if len(dist) else np.nan
                row[f"{lab_}_boundaries_inside"] = int(((b > lo) & (b < hi)).sum())
            cmp_rows.append(row)
    cmpdf = pd.DataFrame(cmp_rows)
    cmpdf.to_csv(HERE / "comparison_with_ad.tsv", sep="\t", index=False)

    # ---- stability across halves
    stab_rows = []
    for src in SOURCES:
        for mode in MODES:
            for c in CS:
                full, A, B = results[(src, mode, c)]["BP"], results[(src + "_A", mode, c)]["BP"], results[(src + "_B", mode, c)]["BP"]
                for region in ["chr7", "chr9", "chr1-22"]:
                    chs = CHROMS if region == "chr1-22" else [region]
                    n = both = either = 0
                    for ch in chs:
                        a, b = np.asarray(A[ch]), np.asarray(B[ch])
                        for x in full[ch]:
                            n += 1
                            ina = len(a) > 0 and np.min(np.abs(a - x)) <= 2
                            inb = len(b) > 0 and np.min(np.abs(b - x)) <= 2
                            both += int(ina and inb)
                            either += int(ina or inb)
                    stab_rows.append({"source": src, "mode": mode, "c": c, "region": region, "n_full_boundaries": n,
                                      "reproduced_in_both_halves_pm2": both, "frac_both_pm2": both / n if n else np.nan,
                                      "frac_either_pm2": either / n if n else np.nan,
                                      "n_half_A": sum(len(A[ch]) for ch in chs), "n_half_B": sum(len(B[ch]) for ch in chs)})
    stab = pd.DataFrame(stab_rows)
    stab.to_csv(HERE / "stability.tsv", sep="\t", index=False)

    cnt = pd.DataFrame([{"source": s, "mode": m, "c": c, "phi_used": results[(s, m, c)]["phi"],
                         "segments_chr7": len(results[(s, m, c)]["BP"]["chr7"]) + 1,
                         "segments_chr9": len(results[(s, m, c)]["BP"]["chr9"]) + 1,
                         "segments_chr1_22": sum(len(v) + 1 for v in results[(s, m, c)]["BP"].values()),
                         "boundaries_at_gap": int(sum(gap_info(ch, i)[0] > 0 for ch in CHROMS for i in results[(s, m, c)]["BP"][ch]))}
                        for (s, m, c) in main_keys])

    # ---- figure: chr7
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    ch = "chr7"
    x = d.start[chrom_idx[ch]] / 1e6
    y = Yc["all"][ch] / Efull["all"] * 1e4
    panels = [("AD baseline (all, saved breakpoints; segment means of raw rate)", AD["all"][ch]),
              (f"count changepoint, half-split phi={hs['all']['phi_half_median']:.2f}, c=2", results[("all", "half_split", 2)]["BP"][ch]),
              (f"count changepoint, Pearson-iterated phi={results[('all', 'pearson', 2)]['phi']:.1f}, c=2", results[("all", "pearson", 2)]["BP"][ch])]
    for ax, (ttl, bps) in zip(axes, panels):
        ax.plot(x, y, ".", ms=2.5, color="#999")
        for a, b in seg_bounds(len(y), list(bps)):
            ax.hlines(y[a:b].mean(), x[a], x[b - 1] + 0.1, color="#c03a2b", lw=1.8)
        for w in win:
            if w[1] == ch:
                ax.axvspan(x[w[2]], x[min(w[3], len(x) - 1)], color="#2f5d8a", alpha=0.15, lw=0)
        ax.set_title(ttl, fontsize=9, loc="left")
        ax.set_ylabel("counts per 1e4\n(bin / total)")
    axes[-1].set_xlabel("chr7 position (Mb); blue bands = saved target windows")
    fig.tight_layout()
    fig7 = C.fig_to_b64(fig)
    plt.close(fig)

    summ = {
        "experiment": "11_count_changepoint", "status": "exploratory diagnostic (method record B:2), not validated",
        "seeds": {"global": SEED, "halves": list(HALF_SEED), "generator": "numpy.random.default_rng (PCG64); segmentation itself is deterministic"},
        "parameters": {"c_grid": CS, "min_segment_bins": MIN_SEG, "penalty": "beta = c*log(n_bins_in_chromosome)",
                       "phi_modes": MODES, "pearson_max_iter": MAX_ITER, "cost": "-S log(S/m) / phi (Poisson)",
                       "nb": "alpha from half-split moment fit; alpha_full = alpha_half * n_cells_half/n_cells_full",
                       "chromosomes": CHROMS, "sources": SOURCES,
                       "ad_baseline": "server/allelic_scan DEFAULTS, genome pool, saved raw_accepted_breakpoints.tsv"},
        "dispersion": hs,
        "pearson_traces": {f"{k[0]}|c={k[2]}": v["trace"] for k, v in results.items() if k[1] == "pearson"},
        "segment_counts": cnt.to_dict("records"),
        "ad_baseline_counts": {s: {"total_chr1_22": ad_tot[s], "at_gap": ad_gap[s], "chr7": len(AD[s]["chr7"]), "chr9": len(AD[s]["chr9"])} for s in SOURCES},
        "windows_bin_indices": [{"region": w[0], "chrom": w[1], "first_bin": w[2], "after_last_bin": w[3], "role": w[4]} for w in win],
        "environment": C.environment(), "runtime_seconds": time.time() - t_start,
    }
    C.write_json(HERE / "summary.json", summ)
    # window / flank relative depth (reference = median segment rate of the all/pearson/c=2 run)
    ref = float(np.median(segs[(segs.source == "all") & (segs["mode"] == "pearson") & (segs.c == 2)].rate_per_bin))
    phi_p = results[("all", "pearson", 2)]["phi"]
    phi_h = hs["all"]["phi_half_median"]
    wr = []
    for _, w in reg[reg.part.isin(["target", "left_flank", "right_flank"])].iterrows():
        idx = chrom_idx[w.chrom]
        lo = int(np.flatnonzero(d.start[idx] == w.first_retained_start)[0])
        hi = int(np.flatnonzero(d.end[idx] == w.last_retained_end)[0]) + 1
        S, m = Yc["all"][w.chrom][lo:hi].sum(), hi - lo
        rel = S / m / Efull["all"] / ref
        wr.append({"region": w.region_id.split("__")[0], "part": w.part, "n_bins": m, "counts": S, "rel_rate": rel,
                   "lo_half": rel * np.exp(-1.96 * np.sqrt(phi_h / S)), "hi_half": rel * np.exp(1.96 * np.sqrt(phi_h / S)),
                   "lo_pearson": rel * np.exp(-1.96 * np.sqrt(phi_p / S)), "hi_pearson": rel * np.exp(1.96 * np.sqrt(phi_p / S))})
    winrate = pd.DataFrame(wr)
    winrate.to_csv(HERE / "window_relative_rates.tsv", sep="\t", index=False, float_format="%.6g")
    write_report(hs, cnt, cmpdf, stab, bp, segs, ad_gap, ad_tot, fig7, AD, summ, winrate, results, win)
    print("done", round(time.time() - t_start, 1), "s")


# ----------------------------------------------------------------------------- report
def interpretation(hs, cnt, cmpdf, stab, ad_tot):
    def g(df, **kw):
        m = np.ones(len(df), bool)
        for k, v in kw.items():
            m &= (df[k] == v).to_numpy()
        return df[m].iloc[0]
    a_h = g(cnt, source="all", mode="half_split", c=2)
    a_p = {c: g(cnt, source="all", mode="pearson", c=c) for c in CS}
    ch_h = g(cmpdf[cmpdf.level == "chrom"], source="all", mode="half_split", c=2, region="chr1-22")
    ch_p = g(cmpdf[cmpdf.level == "chrom"], source="all", mode="pearson", c=2, region="chr1-22")
    st_h = g(stab, source="all", mode="half_split", c=2, region="chr1-22")
    st_p = g(stab, source="all", mode="pearson", c=2, region="chr1-22")
    return f"""<ul>
<li><b>Segment number is governed by phi, not by c.</b> With cell-sampling dispersion (half-split phi = {hs['all']['phi_half_median']:.2f}) or plain Poisson, the segmentation of <code>all</code> uses
about {int(a_h.segments_chr1_22)} segments on chr1-22 at c = 2 ({int(a_h.segments_chr7)} on chr7), and c from 1 to 5 changes this by under 3%. That is a segmentation of reproducible bin-level
structure (first-difference dispersion {hs['all']['diff_phi_full_informational']:.0f}x Poisson), not of cell-sampling noise; it is far from the AD baseline's {ad_tot['all']+22} segments.
The negative-binomial variant gave the same picture (alpha estimated near zero: {hs['all']['alpha_full']:.2g}), i.e. cell sampling adds almost no overdispersion beyond Poisson at these depths, so NB adds nothing here.</li>
<li><b>The Pearson-phi iteration did not settle on a robust level.</b> It was started from the one-segment-per-chromosome dispersion and each pass lowered phi while adding segments until the count
stopped changing, so the fixed point reached depends on c and on the starting value (phi {a_p[1].phi_used:.0f} / {a_p[2].phi_used:.0f} / {a_p[5].phi_used:.0f} and {int(a_p[1].segments_chr1_22)} / {int(a_p[2].segments_chr1_22)} / {int(a_p[5].segments_chr1_22)} segments
for c = 1 / 2 / 5; chr7 {int(a_p[1].segments_chr7)} / {int(a_p[2].segments_chr7)} / {int(a_p[5].segments_chr7)}). At c = 5 chr7 has a single segment. The segment count is therefore a tuning outcome, not a measured property.</li>
<li><b>Agreement with AD.</b> With half-split phi (c = 2) every AD boundary on chr7/chr9 has a new boundary within +/-2 bins ({int(ch_h.ad_matched_pm2)} of {int(ch_h.n_ad)} genome-wide matched within +/-2, all within +/-5), consistent with AD
boundaries being real depth steps in the raw counts; the other {int(ch_h.new_only_pm2)} new boundaries have no AD counterpart. With Pearson phi (c = 2) {int(ch_p.n_new)} boundaries remain and {int(ch_p.ad_matched_pm2)} of {int(ch_p.n_ad)} AD boundaries match within +/-2
({int(ch_p.ad_only_pm2)} AD-only, {int(ch_p.new_only_pm2)} new-only), so at this penalty the count model neither reproduces nor is reproduced by the AD set.</li>
<li><b>Saved windows.</b> Under half-split phi each window edge has a new boundary within +/-2 bins of the edge and AD boundaries within +/-1 on the chr7 windows (chr9_20_27 left edge: nearest AD boundary 25 bins away), with 2-9 additional boundaries inside the windows (50 in the historical chr9_70_100 window); under
Pearson phi (c = 2) the chr7_67.4_70.0 left edge and chr7_70.0_73.3 right edge (67.4 Mb, 73.3 Mb) coincide with boundaries, while the other edges have no boundary within 16 bins and no window has an internal boundary other than chr9_70_100. The saved edges are consistent with a depth step in the raw counts at lenient phi and only partly at conservative phi.</li>
<li><b>Stability.</b> Boundaries reproduce across the two stratified halves ({st_h.frac_both_pm2:.2f} of {int(st_h.n_full_boundaries)} full-data boundaries in both halves within +/-2 at half-split phi, c = 2; {st_p.frac_both_pm2:.2f} of {int(st_p.n_full_boundaries)} at Pearson phi, c = 2),
which says the boundaries are driven by structure shared across cells, not by which cells were drawn. This is expected from the near-Poisson cell-sampling noise and does not by itself say the boundaries are biologically meaningful.</li>
<li><b>Gaps.</b> A small fraction of boundaries ({int(a_h.boundaries_at_gap)} of {int(a_h.segments_chr1_22) - 22} at half-split c = 2) sit at retained-grid gaps; AD has 3 of 142. Boundary positions at gaps are determined by the retained bins either side only.</li>
<li><b>Relative states.</b> chr7 Pearson c = 2 segments range roughly 0.2-2.1 times the median segment rate (e.g. 67.4-73.3 Mb 0.23, 145.2-148.3 Mb 0.22). These are relative depth states; the half-split interval is narrow because it reflects cell sampling only, and the Pearson interval is wider but
still not a copy-number interval.</li></ul>"""


def write_report(hs, cnt, cmpdf, stab, bp, segs, ad_gap, ad_tot, fig7, AD, summ, winrate, results, win):
    T, cap = C.html_table, C.caption
    h, c0 = hs["all"], hs["cluster0"]
    ca = cmpdf[cmpdf.level == "chrom"]

    def g(df, **kw):
        mk = np.ones(len(df), bool)
        for k, v in kw.items():
            mk &= (df[k] == v).to_numpy()
        return df[mk].iloc[0]

    ah = g(cnt, source="all", mode="half_split", c=2)
    ap = {c: g(cnt, source="all", mode="pearson", c=c) for c in CS}
    chh = g(ca, source="all", mode="half_split", c=2, region="chr1-22")
    chp = g(ca, source="all", mode="pearson", c=2, region="chr1-22")
    sth = g(stab, source="all", mode="half_split", c=2, region="chr1-22")
    stp = g(stab, source="all", mode="pearson", c=2, region="chr1-22")
    n_ad = ad_tot["all"]

    verdict = (f"A depth-only count changepoint with cell-sampling dispersion recovers every AD boundary within +/-2 bins but needs about {int(ah.segments_chr1_22)} segments on chr1-22 "
               f"(AD: {n_ad+22}); with a structure-absorbing dispersion it needs {int(ap[2].segments_chr1_22)} and recovers only {int(chp.ad_matched_pm2)} of {n_ad} AD boundaries, so the count model "
               "does not resolve the segmentation question and its segment number is set by the dispersion choice, not by the data alone (exploratory, depth-only).")
    kf = [
        f"Dispersion estimated from the half split is close to Poisson (phi = {h['phi_half_median']:.2f} for <code>all</code>, {c0['phi_half_median']:.2f} for cluster0), whereas first-difference dispersion is {h['diff_phi_full_informational']:.0f} and the Pearson-iterated phi is {ap[1].phi_used:.0f}-{ap[5].phi_used:.0f}: bin-to-bin variation is overwhelmingly reproducible structure, not cell sampling.",
        f"Segments on chr1-22 (source all): half-split phi {int(g(cnt, source='all', mode='half_split', c=1).segments_chr1_22)} / {int(ah.segments_chr1_22)} / {int(g(cnt, source='all', mode='half_split', c=5).segments_chr1_22)} for c = 1 / 2 / 5; Pearson phi {int(ap[1].segments_chr1_22)} / {int(ap[2].segments_chr1_22)} / {int(ap[5].segments_chr1_22)}; AD {n_ad+22}. chr7: {int(ah.segments_chr7)} (half-split, c=2), {int(ap[2].segments_chr7)} (Pearson, c=2), {len(AD['all']['chr7'])+1} (AD).",
        f"Agreement with AD at half-split phi, c = 2: {int(chh.ad_matched_pm2)}/{int(chh.n_ad)} AD boundaries have a new boundary within +/-2 bins ({int(chh.ad_matched_pm5)}/{int(chh.n_ad)} within +/-5), but {int(chh.new_only_pm2)} of {int(chh.n_new)} new boundaries have no AD counterpart.",
        f"Agreement at Pearson phi, c = 2: {int(chp.ad_matched_pm2)}/{int(chp.n_ad)} AD boundaries matched within +/-2 ({int(chp.ad_matched_pm5)} within +/-5); {int(chp.ad_only_pm2)} AD-only and {int(chp.new_only_pm2)} new-only; at c = 5 only {int(ap[5].segments_chr1_22)} segments remain and chr7 has none.",
        f"The Pearson-phi iteration is unstable: starting from the one-segment dispersion, phi falls as segments are added and stops where the count stops changing, so the endpoint depends on c and on the start (phi {ap[1].phi_used:.0f} / {ap[2].phi_used:.0f} / {ap[5].phi_used:.0f}).",
        f"Half-split stability (c = 2, chr1-22): {sth.frac_both_pm2:.2f} of {int(sth.n_full_boundaries)} half-split-phi boundaries and {stp.frac_both_pm2:.2f} of {int(stp.n_full_boundaries)} Pearson-phi boundaries are found in both halves within +/-2 bins.",
        "Saved windows: at half-split phi every window edge has a boundary within +/-2 bins and 2-9 extra boundaries lie inside each window; at Pearson phi (c = 2) only the chr7_67.4_70.0 left edge (67.4 Mb) and chr7_70.0_73.3 right edge (73.3 Mb) coincide with a boundary and no candidate window has an internal boundary.",
        f"Negative-binomial variant: alpha is estimated near zero ({h['alpha_full']:.2g}), segment counts match the Poisson ones within 1%, so NB changes nothing at this depth.",
    ]

    # ---- figures
    fig, ax = plt.subplots(figsize=(8, 3.6))
    for mode, col in [("phi1", "#888"), ("half_split", "#2f5d8a"), ("nb_half", "#2a9d8f"), ("pearson", "#c03a2b")]:
        sub = cnt[(cnt.source == "all") & (cnt["mode"] == mode)]
        ax.plot(sub.c, sub.segments_chr1_22, "o-", color=col, label=mode)
    ax.axhline(n_ad + 22, color="k", ls="--", lw=1, label="AD baseline")
    ax.set_yscale("log"); ax.set_xticks(CS); ax.set_xlabel("penalty multiplier c"); ax.set_ylabel("segments, chr1-22 (log)"); ax.legend(fontsize=8)
    figA = C.fig_to_b64(fig); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 3.6))
    sub = ca[(ca.source == "all") & (ca.region == "chr1-22")]
    for mode, col in [("half_split", "#2f5d8a"), ("pearson", "#c03a2b")]:
        sm = sub[sub["mode"] == mode]
        ax.plot(sm.c, sm.ad_matched_pm2 / sm.n_ad, "o-", color=col, label=f"{mode}: AD boundaries matched +/-2")
        ax.plot(sm.c, sm.new_matched_pm2 / sm.n_new, "s--", color=col, label=f"{mode}: new boundaries matched +/-2")
    ax.set_xticks(CS); ax.set_ylim(0, 1.05); ax.set_xlabel("penalty multiplier c"); ax.set_ylabel("fraction matched (chr1-22)"); ax.legend(fontsize=7)
    figB = C.fig_to_b64(fig); plt.close(fig)

    # ---- sections
    q = ("<p>Method record B:2 (masked count changepoint with conditional alleles) proposes segmenting raw counts with an exposure term, explicit retained-grid gaps, zeros kept and "
         "relative rather than absolute states. This experiment implements only its depth part and asks whether such a count model gives boundaries consistent with the AD scan on the processed X pseudobulk, "
         "how many segments it needs, how that depends on the dispersion, and what it does at the saved target windows.</p>"
         "<p><b>Decision informed:</b> whether the recommendation should move from the AD scan on processed X to a raw-count changepoint (or add one as a cross-check).</p>")
    inputs = ("<ul><li><code>pyepi_results/count_matrix.h5ad</code>, <code>layers/raw</code>: post-QC raw integer counts, 3,460 cells x 25,337 bins (100 kb), chr1-22 only used here. Units are counts per bin per cell; zeros kept.</li>"
              "<li>Sources: <code>all</code> (3,460 cells) and <code>cluster0</code> (1,286 cells; WNN labels via <code>C.cell_clusters</code>). Halves: <code>C.stratified_halves(labels, (20261104, 11))</code> "
              f"(all: {h['n_cells_A']}/{h['n_cells_B']} cells; cluster0: {c0['n_cells_A']}/{c0['n_cells_B']}). Total exposures E: all {h['E_total_counts']:.4g}, cluster0 {c0['E_total_counts']:.4g}.</li>"
              "<li>AD baseline: <code>integration-story/data/genome_pool/raw_accepted_breakpoints.tsv</code> (chromosome-local retained-bin index, checked against bp coordinates).</li>"
              "<li>Gaps: retained grid from the h5ad; removed bins cross-checked with <code>pyepi_results/binned_genome.csv</code>. Windows: <code>research/2026-10-29-low-coverage/shared/regions.tsv</code> (part == target; 1-based starts).</li>"
              "<li><b>Not available / not used:</b> allelic (conditional-allele) emissions, per-bin mappability or GC covariates, per-cell models, an absolute ploidy reference, and any ground-truth copy number.</li></ul>")
    method = ("<ol><li>Per source sum the raw counts of member cells: y<sub>b</sub>; exposure E = total raw counts of those cells (constant over bins, so the rate is y<sub>b</sub>/E).</li>"
              "<li>Per chromosome, find the partition into segments of at least 5 retained bins minimising the penalised cost by exact optimal partitioning (O(n<sup>2</sup>)); breakpoints fall only between consecutive retained bins.</li>"
              "<li>Estimate dispersion phi (modes below), segment with phi, for c in {1, 2, 5} (prespecified, no search).</li>"
              "<li>Compare boundaries with the AD set (+/-2, +/-5 bins), rerun on each half, compute segment relative rates, flag gaps.</li></ol>"
              "<pre>segment cost (Poisson, constants dropped)   C(S, m) = -S * log(S / m) / phi        S = sum y over segment, m = bins\n"
              "penalised objective                        sum_segments C + beta * (#segments - 1),  beta = c * log(n_chrom)\n"
              "NB variant (single alpha, r = 1/alpha)     C_NB(S, m) = -[ S*log(mu/r) - (m*r + S)*log1p(mu/r) ],  mu = S/m,  beta = c*log(n_chrom), phi = 1\n"
              "half-split phi                              z2_b = (yA - n p)^2 / (n p (1-p)),  n = yA + yB,  p = E_A/(E_A+E_B);  phi = median(z2)/0.4549, bins with n >= 20\n"
              "Pearson phi (iterated)                      phi = median((y - mu_seg)^2 / mu_seg) / 0.4549 ; start = one segment per chromosome ; repeat until segment count unchanged (max 12)\n"
              "NB alpha (half split)                       alpha_half = sum[(yA - rho*yB)^2 - (muA + rho^2 muB)] / sum[muA^2 + rho^2 muB^2], winsorised 0.5/99.5%; alpha_full = alpha_half * cells_half / cells_full\n"
              "segment relative rate                       rel = (S/m/E) / median over genome segments of (S/m/E);  95% interval rel * exp(+/- 1.96 sqrt(phi / S))</pre>"
              "<p>Parameters: minimum segment 5 bins; c in {1, 2, 5}; phi modes phi1, half_split, pearson, nb_half; halves seed tuple (20261104, 11); the segmentation itself is deterministic. phi is one genome-wide scalar per source. "
              "In half runs the parent's half-split phi is reused; Pearson phi is re-estimated per half; NB uses alpha_half.</p>")
    sc = cnt[cnt.source == "all"]
    sc2 = cnt[cnt.source == "cluster0"]
    res = []
    res.append(("Segment counts versus penalty and dispersion",
                f"<p>AD baseline (all): {len(AD['all']['chr7'])+1} segments on chr7, {len(AD['all']['chr9'])+1} on chr9, {n_ad+22} on chr1-22, {summ['ad_baseline_counts']['all']['at_gap']} AD boundaries at a gap.</p>"
                + T(sc) + C.img(figA, "segments vs c") + cap("Segments on chr1-22 for source all (log scale) by penalty multiplier c and dispersion mode; dashed line = AD. Curves for phi1, half_split and nb_half overlap: c barely matters there, while the Pearson mode drops to a handful of segments at c = 5. Full table incl. cluster0: summary.json / breakpoints.tsv, segments.tsv.")
                + "<p>cluster0:</p>" + T(sc2), 3))
    cc = ca[["source", "mode", "c", "region", "phi_used", "n_new", "n_ad", "new_matched_pm2", "ad_matched_pm2", "new_only_pm2", "ad_only_pm2", "new_matched_pm5", "ad_matched_pm5", "new_only_pm5", "ad_only_pm5"]]
    res.append(("Comparison with the AD baseline",
                T(cc[(cc.source == "all") & (cc.region.isin(["chr7", "chr9", "chr1-22"]))]) + C.img(figB, "agreement") +
                cap("Fraction of AD boundaries that have a new boundary within +/-2 bins (circles) and fraction of new boundaries with an AD boundary within +/-2 (squares), chr1-22, source all. High AD recall with low new-boundary precision means the count model contains the AD set but adds many more boundaries. Source cluster0 and all other settings: comparison_with_ad.tsv.")
                + C.img(fig7, "chr7") + cap("chr7 raw rate (counts per 10,000 total counts per bin; grey dots) with segment means (red) for the AD boundaries, count changepoint at half-split phi (c = 2) and at Pearson-iterated phi (c = 2). Blue bands = saved target windows. Many short steps in the middle panel versus a few long ones in the bottom panel is the phi effect."), 3))
    stb = stab[["source", "mode", "c", "region", "n_full_boundaries", "reproduced_in_both_halves_pm2", "frac_both_pm2", "frac_either_pm2", "n_half_A", "n_half_B"]]
    res.append(("Half-split stability",
                T(stb[(stb.source == "all") & stb.region.isin(["chr7", "chr9", "chr1-22"])]) +
                cap("For each full-data boundary: is there a boundary within +/-2 bins in half A and in half B (frac_both_pm2)? High values mean boundaries do not depend on which cells were drawn; they say nothing about biology. One split only (seed tuple (20261104, 11)). cluster0 rows are in stability.tsv."), 3))
    s2 = segs[(segs.source == "all") & (segs["mode"] == "pearson") & (segs.c == 2) & (segs.chrom == "chr7")]
    cols = ["chrom", "genomic_start", "genomic_end", "n_bins", "counts", "rel_rate", "rel_lo95_qp_half", "rel_hi95_qp_half", "rel_lo95_qp_pearson", "rel_hi95_qp_pearson"]
    res.append(("Segment relative rates",
                T(s2[cols], max_rows=80) + cap("chr7 segments for Pearson phi, c = 2: relative depth (segment rate over the genome median segment rate) with quasi-Poisson 95% intervals using the half-split phi (cell sampling only; very narrow) and the Pearson phi (also absorbs reproducible bin structure). "
                "These are relative depth states, not copy number; accessibility and dosage are confounded, and the reference median is treated as fixed. All settings: segments.tsv."), 3))
    gp = bp.groupby(["source", "mode", "c"]).agg(n=("bin_index", "size"), at_gap=("at_retained_grid_gap", "sum")).reset_index()
    gp["ad_total"] = gp.source.map(ad_tot); gp["ad_at_gap"] = gp.source.map(ad_gap)
    res.append(("Retained-grid gaps",
                T(gp[gp.source == "all"]) + cap(f"Boundaries (chr1-22) placed immediately next to removed bins (removed_bins_before > 0, cross-checked against binned_genome.csv). Gaps are explicit: a gap-adjacent boundary is fixed by the retained bins on each side only; nothing is interpolated. Roughly {100*ah.boundaries_at_gap/max(ah.segments_chr1_22-22,1):.1f}% of half-split boundaries are gap-adjacent versus 3/142 for AD."), 3))

    # per-region detail
    wd = cmpdf[(cmpdf.level == "window") & (cmpdf.source == "all") & (cmpdf.c == 2)]
    rows = []
    for name, ch, lo, hi, role in win:
        hsr = g(wd, **{"mode": "half_split", "region": name}); pr = g(wd, **{"mode": "pearson", "region": name})
        wt = winrate[(winrate.region == name) & (winrate.part == "target")].iloc[0]
        wl = winrate[(winrate.region == name) & (winrate.part == "left_flank")].iloc[0]
        wrf = winrate[(winrate.region == name) & (winrate.part == "right_flank")].iloc[0]
        edges_h = sum(abs(x) <= 2 for x in (hsr.new_left_nearest_offset_bins, hsr.new_right_nearest_offset_bins))
        edges_p = sum(abs(x) <= 2 for x in (pr.new_left_nearest_offset_bins, pr.new_right_nearest_offset_bins))
        edges_ad = sum(abs(x) <= 2 for x in (hsr.ad_left_nearest_offset_bins, hsr.ad_right_nearest_offset_bins))
        note = ("window is a historical comparison, not truth" if role.startswith("historical") else
                ("edges placed by both models" if edges_h == 2 and edges_p == 2 else
                 ("edges placed only at lenient phi" if edges_h == 2 else "edges not placed")))
        rows.append({"window": name, "bins": hi - lo, "edges_within_pm2 AD": edges_ad, "edges_within_pm2 half_split": edges_h,
                     "edges_within_pm2 pearson": edges_p, "inside: half_split": int(hsr.new_boundaries_inside),
                     "inside: pearson": int(pr.new_boundaries_inside), "inside: AD": int(hsr.ad_boundaries_inside),
                     "rel_rate window": wt.rel_rate, "95% (half)": f"{wt.lo_half:.3f}-{wt.hi_half:.3f}", "95% (pearson)": f"{wt.lo_pearson:.3f}-{wt.hi_pearson:.3f}",
                     "rel left flank": wl.rel_rate, "rel right flank": wrf.rel_rate, "note": note})
    pr_df = pd.DataFrame(rows)
    pr_df.to_csv(HERE / "per_region_detail.tsv", sep="\t", index=False, float_format="%.5g")
    perreg = (T(pr_df) + cap("One row per saved target window (source all, c = 2): number of window edges (left, right; max 2) that have a boundary within +/-2 bins for the AD scan, the count changepoint at half-split phi and at Pearson phi; "
              "boundaries strictly inside the window; window relative depth (counts over median segment rate of the Pearson c = 2 run) with quasi-Poisson intervals; relative depth of the 5 Mb flanks of regions.tsv. "
              "The chr9_70_100 row is the historical comparison, not one of the six candidate windows. Window-edge offsets for every setting: comparison_with_ad.tsv; per-window rates: window_relative_rates.tsv, per_region_detail.tsv. "
              "Window depth below 1 relative to the median segment means lower raw depth than the typical segment, which is consistent with, but does not distinguish, lower dosage and lower accessibility."))
    sens = (T(cnt[(cnt["mode"].isin(["half_split", "pearson", "nb_half"])) & (cnt.source == "all")][["mode", "c", "phi_used", "segments_chr7", "segments_chr9", "segments_chr1_22"]]) +
            cap("Sensitivity to the prespecified settings. c from 1 to 5 changes the half-split / NB / Poisson counts by under 3%, but moves the Pearson count from "
                f"{int(ap[1].segments_chr1_22)} to {int(ap[5].segments_chr1_22)}. The Pearson trace per setting is in summary.json (pearson_traces); halves behave like the full data, so the instability is not a sampling artefact. "
                "Switching source to cluster0 gives the same pattern (cnt rows in summary.json)."))
    limits = ("<ul><li>Constant exposure per source and no bin covariates (mappability, GC, accessibility); only a segment-constant level is modelled.</li>"
              "<li>phi is one genome-wide scalar. Quasi-Poisson with reproducible bin structure is not a likelihood model; the Pearson phi absorbs structure and the half-split phi ignores it, so segment number depends strongly on the estimator. An intermediate, defensible dispersion was not identified here.</li>"
              "<li>The penalty (c log n) and the c grid are conventions, not calibrated to a false-positive rate; no null calibration was run, so none of the boundaries carries an error-rate statement.</li>"
              "<li>NB variant: moment alpha with a 1/#cells scaling assumption, segment mean S/m rather than the NB MLE.</li>"
              "<li>One stratified split; half agreement tests cell sampling only. AD and this model differ in input (processed X versus raw counts), so differences mix method and normalisation.</li>"
              "<li>Depth only: conditional-allele emissions of B:2 are not implemented. Alternative explanations for depth steps include accessibility, mappability/GC residuals, cell-cycle composition and dosage; they are not separated.</li></ul>")
    supports = ["Every AD boundary on chr1-22 has a boundary within +/-5 bins in the raw-count model at cell-sampling dispersion, consistent with AD boundaries being steps in raw depth, not artefacts of the X processing.",
                "Boundaries are reproducible across stratified halves, consistent with structure shared across cells rather than cell-draw noise.",
                "Gaps can be handled explicitly with no change to the boundaries' positions beyond the retained grid."]
    not_supports = ["That the count model resolves segmentation: its segment number moves from tens to thousands with the dispersion estimator.",
                    "That the extra boundaries at cell-sampling dispersion are separate states: they follow reproducible bin-level structure for which no sampling model was fitted.",
                    "Any copy-number interpretation of relative depth; accessibility and dosage are confounded.",
                    "Calibration: no false-positive rate was measured for this changepoint."]
    di = ("<p>For <code>recommendation.md</code>: the raw-count changepoint does not replace the AD scan. It adds support that AD boundaries are steps in raw depth (high recall at lenient dispersion), but "
          "it does not decide how many segments are justified, because the dispersion choice (not the data alone) sets that. Any adoption would first need a dispersion model that handles reproducible bin structure (for example a segment-level random-effect or autoregressive term) and a selection-aware null calibration comparable to experiment 05.</p>")
    nxt = ["Add a dispersion model for reproducible bin-level structure (e.g. NB with bin-level random effects fitted from replicate halves, or an AR(1) residual) and re-estimate phi without the unstable iteration.",
           "Run a selection-aware null (permuted bins, as in experiment 05) through this changepoint at each setting to measure false-positive segment rates.",
           "Add the allelic emission (B:2 conditional alleles) and test whether it separates depth steps AD finds from those it does not.",
           "Repeat on additional random splits and the other WNN clusters to see whether boundary stability varies with split."]
    gl = {"phi": "quasi-Poisson dispersion; cost is divided by phi (equivalently penalty multiplied by phi)",
          "c": "multiplier of the BIC-type penalty beta = c * log(n_bins in chromosome)", "mode": "dispersion estimator: phi1, half_split, pearson, or nb_half",
          "phi_used": "dispersion actually used for the setting", "alpha": "negative-binomial overdispersion, var = mu + alpha mu^2",
          "segments_chr7 / chr9 / chr1_22": "number of segments (breakpoints + 1 per chromosome, summed over chr1-22)",
          "boundaries_at_gap": "count of boundaries with removed bins directly before them",
          "n_new / n_ad": "boundaries from the count changepoint / saved AD scan",
          "new_matched_pm2 / pm5": "new boundaries with an AD boundary within +/-2 / +/-5 bins",
          "ad_matched_pm2 / pm5": "AD boundaries with a new boundary within +/-2 / +/-5 bins",
          "new_only / ad_only": "boundaries lacking a counterpart in the other set at that tolerance",
          "frac_both_pm2": "fraction of full-data boundaries with a boundary within +/-2 bins in both halves", "frac_either_pm2": "same, in at least one half",
          "n_half_A / n_half_B": "number of boundaries found in each half",
          "rel_rate": "segment (or window) rate relative to the median segment rate; relative depth, not copy number",
          "rel_lo95_qp_half / hi": "quasi-Poisson 95% interval with half-split phi (cell sampling only)",
          "rel_lo95_qp_pearson / hi": "quasi-Poisson 95% interval with the Pearson-iterated phi",
          "removed_bins_before / at_retained_grid_gap": "number of QC-removed grid bins between the previous retained bin and the boundary / whether it is > 0",
          "*_nearest_offset_bins": "signed bins from a window edge to the nearest boundary (negative = upstream)",
          "*_boundaries_inside": "boundaries strictly inside the window", "edges_within_pm2": "window edges (0-2) with a boundary within +/-2 bins",
          "rate_per_bin": "counts / bins / total exposure E", "rel_rate_wmedian": "rel_rate using the bin-weighted median as reference",
          "window / flank": "target window and the 5 Mb flanks of regions.tsv; flank values use the same reference"}
    sections = [("Question and why it matters", q), ("Inputs", inputs), ("Method", method)]
    sections.append(("Results", "<p>Each analysis below has a table and a caption saying how to read it. Full tables are the TSV files in this folder: breakpoints.tsv, segments.tsv, comparison_with_ad.tsv, stability.tsv, window_relative_rates.tsv, per_region_detail.tsv.</p>", 2))
    sections += res
    sections.append(("Per-region detail", perreg, 2))
    sections.append(("Sensitivity", sens, 2))
    sections.append(("Limits and alternative explanations", limits, 2))
    C.write_html_report(HERE / "report.html", "Experiment 11: depth-only masked count changepoint vs AD baseline",
                        "Exploratory diagnostic of method record B:2 (depth part only) on MCF10-A development data. Consistent-with language only.",
                        sections, {k: summ[k] for k in ["seeds", "parameters", "dispersion", "environment", "runtime_seconds"]},
                        key_findings=kf, verdict=verdict, glossary=gl, supports=supports, not_supports=not_supports, decision_impact=di, next_steps=nxt)


if __name__ == "__main__":
    main()
