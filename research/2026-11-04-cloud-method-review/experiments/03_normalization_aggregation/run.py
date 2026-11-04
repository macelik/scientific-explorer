"""Experiment 03: aggregation estimand, CN-assignment thresholds and genome-wide denominator sensitivity.

Segmentation is held fixed at the saved genome-pool breakpoints (source-specific). CN states are
assigned with the pinned, unmodified pyEpi assign_gainloss_new (via server/pyepi_adapter), exactly as
server/cn_recall.cluster_rows does.

A. Estimands (same cells, same boundaries):
   x_mean     - mean of processed X over cells (current workbench pseudobulk);
   raw_ratio  - sum of raw counts over cells / sum of cell exposures (raw totals), x 1e5: an
                exposure-weighted rate with no GC correction, rounding or clipping;
   raw_cellmean - mean over cells of raw/exposure x 1e5 (equal cell weights, no GC correction).
   These are different estimands; differences show how much the call depends on the choice.
B. Analytic decision thresholds implied by weighted_scale_search's state weights, expressed as
   segment mean / best_s, and each `all` segment's margin to the nearest threshold.
C. Denominator sensitivity: assign_gainloss_new with one chromosome's bins removed at a time; the
   number of bins on *other* chromosomes whose combined call changes.
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

sys.path.insert(0, str(C.REPO))
from server import cn_recall as cr  # noqa: E402
from server import pyepi_adapter as sci  # noqa: E402

SOURCES = ["cluster0", "cluster1", "cluster2", "cluster3", "cluster4", "all"]


def thresholds():
    """Pairwise decision boundaries of argmin_k w_k (m - s k)^2 in units of m/s (adjacent states)."""
    w = np.concatenate((np.array([2.0, 1.5]), 1 / np.arange(2, 7)))
    out = []
    for k in range(6):
        a, b = np.sqrt(w[k]), np.sqrt(w[k + 1])
        out.append({"between": f"{k}|{k+1}", "w_lo": w[k], "w_hi": w[k + 1], "boundary_m_over_s": (a * k + b * (k + 1)) / (a + b),
                    "unweighted_midpoint": k + 0.5})
    return pd.DataFrame(out)


def main():
    p = C.default_paths()
    sci.load(str(p["package"]))
    d = C.H5AD(p["h5ad"])
    lab = C.cell_clusters(d, p)
    masks = {f"cluster{c}": lab == str(c) for c in range(5)}
    masks["all"] = np.ones(d.n_cells, bool)
    x_mean = C.group_means(d, "X", masks)
    raw_mean = C.group_means(d, "layers/raw", masks)
    tot = np.asarray(d.row_nnz("layers/raw"), float) * 0  # placeholder replaced below
    tot = np.zeros(d.n_cells)
    cm = {k: np.zeros(d.n_bins) for k in masks}
    for lo in range(0, d.n_cells, 250):
        hi = min(d.n_cells, lo + 250)
        blk = d.rows("layers/raw", lo, hi)
        t = blk.sum(1)
        tot[lo:hi] = t
        rate = blk / t[:, None] * 1e5
        for k, m in masks.items():
            sel = m[lo:hi]
            if sel.any():
                cm[k] += rate[sel].sum(0)
    raw_cellmean = {k: cm[k] / masks[k].sum() for k in masks}
    raw_ratio = {k: raw_mean[k] * masks[k].sum() / tot[masks[k]].sum() * 1e5 for k in masks}
    est = {"x_mean": x_mean, "raw_ratio": raw_ratio, "raw_cellmean": raw_cellmean}

    bp = pd.read_csv(p["genome_pool"] / "raw_accepted_breakpoints.tsv", sep="\t")
    chroms = list(dict.fromkeys(d.seq))
    bounds = cr.boundaries_from_bp(bp, SOURCES, chroms)
    saved = pd.read_csv(p["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")

    # harness check + estimand comparison
    calls, rows_all, best = {}, [], {}
    for s in SOURCES:
        seg = cr.seg_labels(d.seq, bounds[s])
        for e, v in est.items():
            rows, comb = cr.cluster_rows(s, int(masks[s].sum()), v[s], seg, d.seq)
            calls[(s, e)] = comb
            best[(s, e)] = rows[0]["best_s"]
            if e == "x_mean":
                ref = saved[saved.group == s].reset_index(drop=True)
                got = pd.DataFrame(rows)
                same = len(ref) == len(got) and np.allclose(ref.cn_call, got.cn_call) and np.allclose(ref.best_s, got.best_s)
                rows_all.append({"source": s, "harness_identical_to_saved_karyogram": bool(same), "segments": len(got)})
    harness = pd.DataFrame(rows_all)
    comp = []
    for s in SOURCES:
        b = calls[(s, "x_mean")]
        for e in ["raw_ratio", "raw_cellmean"]:
            a = calls[(s, e)]
            comp.append({"source": s, "estimand": e, "best_s_x_mean": best[(s, "x_mean")], "best_s": best[(s, e)],
                         "pct_bins_call_changed": 100 * float((a != b).mean()),
                         "bins_loss_x": int((b < 1).sum()), "bins_loss_alt": int((a < 1).sum()),
                         "bins_gain_x": int((b > 1).sum()), "bins_gain_alt": int((a > 1).sum())})
    comp = pd.DataFrame(comp)
    comp.to_csv(HERE / "estimand_call_comparison.tsv", sep="\t", index=False)

    # saved windows / literature loci: call per estimand (all source)
    reg = C.shared_regions()
    loci = [(r.parent, r.chrom, r.first_retained_start, r.last_retained_end) for r in reg[reg.part == "target"].itertuples()]
    loci += [("lit_MYC_8q24", "chr8", 127_000_001, 128_500_000), ("lit_19q13_qter", "chr19", 35_000_001, 10**9),
             ("lit_5q13_qter", "chr5", 70_000_001, 10**9), ("lit_9p21_zone_approx", "chr9", 20_000_001, 24_000_000)]
    lrows = []
    for name, ch, a, b in loci:
        m = (d.seq == ch) & (d.start >= a) & (d.end <= b)
        rec = {"region": name, "bins": int(m.sum())}
        for e in est:
            v = calls[("all", e)][m]
            rec[f"{e}_frac_loss"] = float((v < 1).mean()); rec[f"{e}_frac_gain"] = float((v > 1).mean())
            rec[f"{e}_mean_over_genome_median"] = float(est[e]["all"][m].mean() / np.median(est[e]["all"]))
        lrows.append(rec)
    lt = pd.DataFrame(lrows)
    lt.to_csv(HERE / "region_calls_by_estimand.tsv", sep="\t", index=False)

    # B. thresholds and margins for `all`
    th = thresholds()
    th.to_csv(HERE / "state_thresholds.tsv", sep="\t", index=False)
    seg = cr.seg_labels(d.seq, bounds["all"])
    xs = x_mean["all"]
    norm = xs / sci.trimmed_mean_iqr(xs, lb=False)
    segmeans = pd.Series(norm).groupby(seg).mean()
    seglen = pd.Series(norm).groupby(seg).size()
    med = float(np.median(segmeans[seg].to_numpy()))
    bs, cont, ints = sci.weighted_scale_search(segmeans.values.copy(), seglen.values, med)
    grid = (med * 0.5, min(1.2, med * 0.75))
    ratio = segmeans.values / bs
    bnd = th.boundary_m_over_s.to_numpy()
    margin = np.min(np.abs(ratio[:, None] - np.r_[bnd, 1.0, 3.0][None, :]), axis=1)   # include Watson cut-points 1 and 3
    mt = pd.DataFrame({"segment": segmeans.index, "n_bins": seglen.values, "m_over_s": ratio, "holmes_k": ints, "watson_cont": cont,
                       "min_distance_to_a_cutpoint": margin})
    mt.to_csv(HERE / "all_segment_margins.tsv", sep="\t", index=False)
    s_sens = []
    for f in [0.95, 0.97, 1.03, 1.05]:
        s2 = bs * f
        d2 = ((segmeans.values[:, None] - s2 * np.arange(7)[None, :]) ** 2) * np.concatenate((np.array([2.0, 1.5]), 1 / np.arange(2, 7)))[None, :]
        k2 = np.argmin(d2, 1)
        h1 = np.clip(ints, 1, 3) - 1; h2 = np.clip(k2, 1, 3) - 1
        c2 = np.round(segmeans.values / s2, 2)
        w1 = np.where(cont <= 1, 0, np.where(cont >= 3, 2, 1)); w2 = np.where(c2 <= 1, 0, np.where(c2 >= 3, 2, 1))
        changed = (h1 != h2) | (w1 != w2)
        s_sens.append({"best_s_factor": f, "segments_changed": int(changed.sum()), "bins_changed": int(seglen.values[changed].sum())})
    s_sens = pd.DataFrame(s_sens)

    edge = []
    for src in SOURCES:
        for e in est:
            sg = cr.seg_labels(d.seq, bounds[src])
            v = est[e][src]
            nv = v / sci.trimmed_mean_iqr(v, lb=False)
            sm = pd.Series(nv).groupby(sg).mean()
            md = float(np.median(sm[sg].to_numpy()))
            b_, *_ = sci.weighted_scale_search(sm.values.copy(), pd.Series(nv).groupby(sg).size().values, md)
            lo_, hi_ = md * 0.5, min(1.2, md * 0.75)
            edge.append({"source": src, "estimand": e, "best_s": b_, "grid_lo": lo_, "grid_hi": hi_,
                         "at_lower_edge": bool(abs(b_ - lo_) < 1e-12), "at_upper_edge": bool(abs(b_ - hi_) < 1e-12)})
    edge = pd.DataFrame(edge)
    edge.to_csv(HERE / "scale_grid_edge.tsv", sep="\t", index=False)

    # C. leave-one-chromosome-out denominator sensitivity (x_mean, all sources)
    loo = []
    for s in SOURCES:
        full = calls[(s, "x_mean")]
        segs = cr.seg_labels(d.seq, bounds[s])
        for ch in chroms:
            keep = d.seq != ch
            *_, comb, b2 = sci.assign_gainloss_new(np.asarray(x_mean[s][keep], float), segs[keep])
            ch_changed = np.asarray(comb) != full[keep]
            loo.append({"source": s, "chromosome_removed": ch, "best_s_full": best[(s, "x_mean")], "best_s_loo": float(b2),
                        "bins_changed_elsewhere": int(ch_changed.sum()), "pct_changed_elsewhere": 100 * float(ch_changed.mean())})
    loo = pd.DataFrame(loo)
    loo.to_csv(HERE / "leave_one_chromosome_out.tsv", sep="\t", index=False)
    loo_sum = loo.groupby("source").agg(max_pct_changed=("pct_changed_elsewhere", "max"),
                                        median_pct_changed=("pct_changed_elsewhere", "median"),
                                        worst_chromosome=("pct_changed_elsewhere", lambda x: loo.loc[x.idxmax(), "chromosome_removed"]),
                                        best_s_range=("best_s_loo", lambda x: f"{x.min():.4f}–{x.max():.4f}")).reset_index()

    res = {"experiment": "03_normalization_aggregation", "deterministic": True, "harness": harness.to_dict("records"),
           "estimand_comparison": comp.to_dict("records"), "thresholds": th.to_dict("records"),
           "all_best_s": bs, "all_s_grid": grid, "best_s_at_grid_edge": bool(abs(bs - grid[0]) < 1e-9 or abs(bs - grid[1]) < 1e-9),
           "all_segments_within_0.1_of_cutpoint": int((margin < 0.1).sum()), "all_segments": len(mt),
           "best_s_sensitivity": s_sens.to_dict("records"), "scale_grid_edge": edge.to_dict("records"), "leave_one_chromosome_out": loo_sum.to_dict("records"),
           "environment": C.environment()}
    C.write_json(HERE / "summary.json", res)

    fig, ax = plt.subplots(1, 2, figsize=(14, 4))
    ax[0].hist(ratio, bins=80, weights=seglen.values, color="#2f5d8a")
    for v in bnd[:4]:
        ax[0].axvline(v, color="#9b2c2c", lw=1)
    for v in [1, 3]:
        ax[0].axvline(v, color="#8a5a00", lw=1, ls="--")
    ax[0].set_xlabel("segment mean / best_s (all source)"); ax[0].set_ylabel("bins"); ax[0].set_title("Holmes cut-points (red), Watson (dashed)")
    piv = loo.pivot(index="chromosome_removed", columns="source", values="pct_changed_elsewhere").loc[chroms]
    ax[1].imshow(piv.T.values, aspect="auto", cmap="Reds"); ax[1].set_yticks(range(len(piv.columns))); ax[1].set_yticklabels(piv.columns)
    ax[1].set_xticks(range(len(chroms))); ax[1].set_xticklabels([c[3:] for c in chroms], fontsize=7)
    ax[1].set_title("% of other bins whose call changes when one chromosome is removed")
    b64 = C.fig_to_b64(fig); plt.close(fig)
    frag = [
        ("Question", "<p>Holding segmentation fixed, how much do CN calls depend on the pseudobulk estimand, on the weighted scale search's "
         "thresholds, and on which chromosomes enter the genome-wide denominator?</p>"),
        ("Harness check", C.html_table(harness)),
        ("A. Estimand comparison (calls on identical boundaries)", C.html_table(comp) + C.html_table(lt)),
        ("B. Thresholds implied by the state weights", C.html_table(th) +
         f"<p>For <code>all</code>: best_s = {bs:.4f} on grid [{grid[0]:.4f}, {grid[1]:.4f}]; "
         f"{int((margin < 0.1).sum())} of {len(mt)} segments lie within 0.1 (in m/s units) of a cut-point.</p>" + C.html_table(s_sens) +
         C.html_table(edge) + f"<div class='card'><b>{int(edge.at_lower_edge.sum())} of {len(edge)}</b> source × estimand fits "
         "end at the lower edge of the s grid (s = 0.5 × median segment mean). There the fitted scale is not estimated from the "
         "data: the median-weighted segment is defined as k = 2, and every call becomes a fixed threshold on segment mean / genome "
         "median (Holmes loss below about 0.68, Watson loss at or below 0.5, Holmes gain above about 1.22, Watson gain at or above 1.5). "
         "The reference level is the genome median of the same profile, which broad events and baseline structure can move.</div>" +
         "<div class='card'>The weights (2, 1.5, 0.5, 0.33, 0.25, 0.2, 0.17 for k = 0…6) make the loss/base cut-point fall at "
         "m/s ≈ 1.37 instead of 1.5 and the base/gain cut-point at ≈ 2.45 instead of 2.5. Losses therefore need a deeper "
         "contrast than gains, and the Watson rule (≤ 1 loss, ≥ 3 gain) is stricter still; the combined call takes 0.5/1.5 "
         "when only Holmes moves. The code comment says the weights penalize deviations from diploid, but weights fall "
         "monotonically above k = 2, so higher states are cheaper than k = 2.</div>"),
        ("C. Genome-wide denominator", C.html_table(loo_sum) + C.img(b64, "thresholds and leave-one-out")),
        ("Interpretation", "<ul><li>Calls depend on the estimand: processed-X means, exposure-weighted raw rates and equal-weight raw "
         "rates assign different states on the same boundaries (A). The raw estimands have no GC correction, so this mixes GC with "
         "aggregation effects; read it as a sensitivity range, not as one estimand being right.</li>"
         "<li>The relative CN scale is fitted genome-wide, so a call in one region can change when another chromosome is added or "
         "removed (C). A per-region evidence report should state the reference level it was compared against.</li></ul>"),
    ]
    C.write_html_report(HERE / "report.html", "Normalization and CN assignment", "Experiment 03 · estimand, thresholds and denominator sensitivity with fixed segmentation",
                        frag, {"command": "python -I experiments/03_normalization_aggregation/run.py", "summary": "summary.json", "environment": C.environment()})
    pd.set_option("display.width", 250)
    print(harness.to_string()); print(comp.to_string()); print(lt.to_string()); print(th.to_string()); print(s_sens.to_string()); print(loo_sum.to_string())
    print("best_s", bs, grid, "near cutpoint", int((margin < 0.1).sum()), len(mt))


if __name__ == "__main__":
    main()
