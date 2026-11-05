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
    pc = comp.pct_bins_call_changed
    ca = comp[comp.source == "all"].set_index("estimand")
    gain_ratio = (comp.bins_gain_alt / comp.bins_gain_x)
    thr = th.set_index("between").boundary_m_over_s
    n_edge, n_fits = int(edge.at_lower_edge.sum()), len(edge)
    harness_ok = bool(harness.harness_identical_to_saved_karyogram.all())
    mt_close = int((margin < 0.1).sum())
    lt_i = lt.set_index("region")
    verdict = (f"With segmentation held fixed, changing the pseudobulk estimand changes {pc.min():.0f}&ndash;{pc.max():.0f}% of bin calls, all {n_edge} of {n_fits} CN scale fits sit at the lower grid edge so calls are fixed thresholds on segment mean / genome median, "
               f"and removing one chromosome changes up to {loo_sum.max_pct_changed.max():.1f}% of calls elsewhere.")
    findings = [
        f"The workbench harness reproduces the saved karyogram calls and scales for {int(harness.harness_identical_to_saved_karyogram.sum())} of {len(harness)} sources ({'all identical' if harness_ok else 'NOT all identical'}).",
        f"Estimand choice changes {pc.min():.1f}&ndash;{pc.max():.1f}% of bin calls on identical boundaries (for <code>all</code>: {ca.loc['raw_ratio','pct_bins_call_changed']:.1f}% for raw_ratio, {ca.loc['raw_cellmean','pct_bins_call_changed']:.1f}% for raw_cellmean); gain bins rise {gain_ratio.min():.1f}&ndash;{gain_ratio.max():.1f}x "
        f"(<code>all</code>: {int(ca.loc['raw_ratio','bins_gain_x'])} &rarr; {int(ca.loc['raw_ratio','bins_gain_alt'])}) while loss-bin counts change by smaller amounts in either direction.",
        f"<b>{n_edge} of {n_fits}</b> source &times; estimand scale fits end at the lower grid edge (s = 0.5 &times; median segment mean); for <code>all</code> best_s = {bs:.4f} on [{grid[0]:.4f}, {grid[1]:.4f}], never at the upper edge ({int(edge.at_upper_edge.sum())}).",
        f"Implied cut-points on segment mean / best_s: Holmes loss/base {thr['1|2']:.3f}, base/gain {thr['2|3']:.3f} (unweighted midpoints 1.5 and 2.5); Watson uses 1 and 3. With best_s at 0.5 &times; median this is a loss below about {thr['1|2']*0.5:.2f} &times; the genome median and a gain above about {thr['2|3']*0.5:.2f} &times;.",
        f"{mt_close} of {len(mt)} <code>all</code> segments lie within 0.1 (m/s units) of a cut-point; moving best_s by &minus;5% to +5% changes {int(s_sens.segments_changed.min())}&ndash;{int(s_sens.segments_changed.max())} segments ({int(s_sens.bins_changed.min()):,}&ndash;{int(s_sens.bins_changed.max()):,} bins).",
        f"Leave-one-chromosome-out: removing one chromosome changes at most {loo_sum.max_pct_changed.max():.1f}% of calls elsewhere (largest for {loo_sum.loc[loo_sum.max_pct_changed.idxmax(),'source']} / {loo_sum.loc[loo_sum.max_pct_changed.idxmax(),'worst_chromosome']}); median across removed chromosomes {loo_sum.median_pct_changed.min():.1f}&ndash;{loo_sum.median_pct_changed.max():.1f}% by source.",
        f"The four chr7 saved windows are 100% loss in all three estimands (mean/genome median {lt[lt.region.str.startswith('chr7')].x_mean_mean_over_genome_median.min():.2f}&ndash;{lt[lt.region.str.startswith('chr7')].x_mean_mean_over_genome_median.max():.2f} in processed X); the GC-rich 19q13-qter literature region is {100*lt_i.loc['lit_19q13_qter','x_mean_frac_loss']:.0f}% loss / {100*lt_i.loc['lit_19q13_qter','x_mean_frac_gain']:.0f}% gain with processed X but {100*lt_i.loc['lit_19q13_qter','raw_ratio_frac_gain']:.0f}% gain with raw rates.",
    ]
    region_notes = []
    for r in lt.itertuples():
        fl = [r.x_mean_frac_loss, r.raw_ratio_frac_loss, r.raw_cellmean_frac_loss]
        fg = [r.x_mean_frac_gain, r.raw_ratio_frac_gain, r.raw_cellmean_frac_gain]
        if max(fl) - min(fl) < 0.1 and max(fg) - min(fg) < 0.1:
            note = f"Calls agree across estimands ({100*r.x_mean_frac_loss:.0f}% loss, {100*r.x_mean_frac_gain:.0f}% gain of {r.bins} bins)."
        else:
            note = (f"Calls depend on estimand: loss {100*r.x_mean_frac_loss:.0f}%/{100*r.raw_ratio_frac_loss:.0f}%, gain {100*r.x_mean_frac_gain:.0f}%/{100*r.raw_ratio_frac_gain:.0f}% "
                    f"(processed X / raw rate); mean over genome median {r.x_mean_mean_over_genome_median:.2f} vs {r.raw_ratio_mean_over_genome_median:.2f}.")
        region_notes.append(note)
    lt_show = lt[["region", "bins", "x_mean_frac_loss", "x_mean_frac_gain", "x_mean_mean_over_genome_median", "raw_ratio_frac_loss", "raw_ratio_frac_gain", "raw_ratio_mean_over_genome_median"]].copy()
    lt_show["note"] = region_notes
    lt_html = C.html_table(lt_show).replace("<td>Calls ", "<td style='text-align:left;white-space:normal'>Calls ")
    comp_show = comp[comp.estimand == "raw_ratio"].drop(columns=["estimand"])
    edge_show = edge.groupby("source").agg(x_mean=("best_s", "first"), at_lower_edge_n=("at_lower_edge", "sum"), grid_lo=("grid_lo", "first"), grid_hi=("grid_hi", "first")).reset_index()
    sections = [
        ("Question and why it matters",
         "<p><b>Question.</b> Holding segmentation fixed at the saved genome-pool breakpoints, how much do CN calls depend on the pseudobulk estimand, on the weighted scale search&rsquo;s thresholds, and on which chromosomes enter the genome-wide denominator?</p>"
         "<p><b>Why it matters.</b> The recommendation states that CN assignment is a fixed threshold on the genome median, and that a per-region evidence layer should say which reference level it was compared against. "
         "This experiment quantifies how far the call states move when only the estimand, the scale, or the denominator changes, i.e. how much of a &lsquo;gain&rsquo; or &lsquo;loss&rsquo; label is a property of the assignment rule rather than of the sample.</p>"),
        ("Inputs",
         f"<ul><li><b>Matrix.</b> <code>count_matrix.h5ad</code>: {d.n_cells:,} cells &times; {d.n_bins:,} bins; <code>X</code> (processed) and <code>layers/raw</code>.</li>"
         "<li><b>Groups.</b> WNN clusters (<code>wnn_leiden_0.3</code>, labels 0&ndash;4 from <code>integrated.h5mu</code>) plus <code>all</code> cells: six pseudobulk sources.</li>"
         "<li><b>Boundaries.</b> <code>genome_pool/raw_accepted_breakpoints.tsv</code> (source-specific breakpoints from the workbench segmentation); the saved <code>karyogram_cluster_segments.tsv</code> is used only to check the harness.</li>"
         f"<li><b>Regions.</b> The 7 saved target windows from the frozen shared <code>regions.tsv</code> (retained-bin extents) and 4 literature regions: MYC 8q24 (chr8:127.0&ndash;128.5 Mb), 19q13-qter (chr19 from 35 Mb), 5q13-qter (chr5 from 70 Mb), approximate 9p21 zone (chr9:20&ndash;24 Mb) ({len(lt)} rows). Bin counts are the retained bins in each.</li>"
         "<li><b>Units.</b> Processed X is counts per 1e5 per cell; raw rates are fragments per 1e5 of a cell&rsquo;s retained-bin total. Segment means are divided by the upper-trimmed mean before CN assignment.</li>"
         "<li><b>Not available.</b> Any independent copy-number truth, so no state is scored as right or wrong here; only changes between rules are measured.</li></ul>"),
        ("Method",
         "<ol><li><b>Estimands</b> (same cells, same boundaries): <code>x_mean</code> = mean of processed X over cells (current workbench pseudobulk); <code>raw_ratio</code> = &Sigma; raw counts / &Sigma; cell exposure &times; 1e5 (exposure-weighted, no GC correction, rounding or clipping); "
         "<code>raw_cellmean</code> = mean over cells of raw/exposure &times; 1e5 (equal cell weights, no GC correction).</li>"
         "<li><b>CN assignment</b> uses the pinned, unmodified pyEpi <code>assign_gainloss_new</code> via <code>server/pyepi_adapter</code>, exactly as <code>server/cn_recall.cluster_rows</code> does. Calls are compared bin by bin.</li>"
         "<li><b>Thresholds.</b> For adjacent states k and k+1 the decision boundary of the weighted squared error is solved analytically; each <code>all</code> segment&rsquo;s distance to the nearest cut-point is recorded, including the Watson cut-points 1 and 3.</li>"
         "<li><b>Scale sensitivity.</b> best_s is multiplied by 0.95, 0.97, 1.03, 1.05 and segments whose Holmes or Watson class changes are counted.</li>"
         "<li><b>Leave-one-chromosome-out.</b> <code>assign_gainloss_new</code> is rerun with one chromosome&rsquo;s bins removed from <code>x_mean</code>; bins on the remaining chromosomes whose combined call changes are counted.</li></ol>"
         "<pre>norm       = pseudobulk / upper_trimmed_mean(pseudobulk)\nseg_mean   = mean(norm) per segment;  med = median(seg_mean over bins)\ngrid       s in [0.5*med, min(1.2, 0.75*med)]\nweights    w_k = (2, 1.5, 1/2, 1/3, 1/4, 1/5, 1/6) for k = 0..6\ncost(s)    = sum_seg len * min_k w_k (m - s*k)^2\nboundary   (k|k+1) in m/s = (sqrt(w_k)*k + sqrt(w_k+1)*(k+1)) / (sqrt(w_k) + sqrt(w_k+1))</pre>"
         "<p>Deterministic: no random numbers. Calls are not scored against truth.</p>"),
        ("Results", "", 2),
        ("Harness check", C.html_table(harness) + C.caption("For each source, whether the harness (cluster_rows on x_mean) reproduces the saved karyogram calls and best_s exactly. All True means the differences below come from the rule or estimand, not from a re-implementation."), 3),
        ("A. Estimand comparison on identical boundaries",
         C.html_table(comp_show) + C.caption("Per source, processed-X means against the exposure-weighted raw rate (<code>raw_ratio</code>; the equal-weight <code>raw_cellmean</code> rows are in <code>estimand_call_comparison.tsv</code> and are within "
                                             f"{abs(comp[comp.estimand=='raw_ratio'].pct_bins_call_changed.values - comp[comp.estimand=='raw_cellmean'].pct_bins_call_changed.values).max():.1f} percentage points). "
                                             "<code>pct_bins_call_changed</code> is the share of all bins whose combined call differs; <code>bins_loss_*</code> / <code>bins_gain_*</code> are the numbers of loss / gain bins. Look for a large rise in gain bins with raw rates."), 3),
        ("B. Thresholds implied by the state weights",
         C.html_table(th) + C.caption("Decision boundaries between adjacent states in units of segment mean / best_s. <code>boundary_m_over_s</code> is the weighted cut-point; <code>unweighted_midpoint</code> is where equal weights would put it. The weights shift the loss/base cut to 1.37 and base/gain to 2.45.") +
         f"<p>For <code>all</code>: best_s = {bs:.4f} on grid [{grid[0]:.4f}, {grid[1]:.4f}]; {mt_close} of {len(mt)} segments lie within 0.1 of a cut-point.</p>" +
         C.html_table(s_sens) + C.caption("Sensitivity of the <code>all</code> source to scaling best_s by <code>best_s_factor</code>: segments and bins whose Holmes or Watson class changes.") +
         C.html_table(edge_show) + C.caption("Per source, the fitted scale for x_mean, how many of the three estimands end at the lower grid edge (<code>at_lower_edge_n</code> of 3), and the grid limits. Full table: <code>scale_grid_edge.tsv</code>.") +
         f"<div class='card'><b>{n_edge} of {n_fits}</b> source &times; estimand fits end at the lower edge of the s grid (s = 0.5 &times; median segment mean). There the fitted scale is not estimated from the "
         "data: the median-weighted segment is defined as k = 2, and every call becomes a fixed threshold on segment mean / genome median "
         f"(Holmes loss below about {thr['1|2']*0.5:.2f}, Watson loss at or below 0.5, Holmes gain above about {thr['2|3']*0.5:.2f}, Watson gain at or above 1.5). The reference level is the genome median of the same profile, which broad events and baseline structure can move.</div>"
         "<div class='card'>The weights (2, 1.5, 0.5, 0.33, 0.25, 0.2, 0.17 for k = 0&hellip;6) make losses need a deeper contrast than gains, and the Watson rule (&le; 1 loss, &ge; 3 gain) is stricter still. "
         "The code comment says the weights penalize deviations from diploid, but the weights fall monotonically above k = 2, so higher states are cheaper than k = 2.</div>", 3),
        ("C. Genome-wide denominator (leave one chromosome out)",
         C.html_table(loo_sum) + C.caption("Per source: the largest and median percentage of bins on <em>other</em> chromosomes whose combined call changes when one chromosome is removed, the chromosome responsible for the maximum, and the range of the refitted best_s. Full table: <code>leave_one_chromosome_out.tsv</code>.") +
         C.img(b64, "thresholds and leave-one-out") +
         C.caption("Left: histogram of <code>all</code> segment mean / best_s weighted by bins; red lines are the Holmes cut-points, dashed lines the Watson cut-points 1 and 3. Bars near a line are fragile. "
                   "Right: percentage of other bins whose call changes when the chromosome on the x-axis is removed, per source; look for dark cells."), 3),
        ("Per-region detail", lt_html + C.caption("One row per saved window and literature region, for the <code>all</code> source. <code>*_frac_loss</code> / <code>*_frac_gain</code> are the shares of the region&rsquo;s bins called loss / gain; <code>*_mean_over_genome_median</code> is the region mean divided by the genome median of the same estimand. "
                                                  "The raw_cellmean columns are in <code>region_calls_by_estimand.tsv</code>. The note is derived from the loss and gain shares and is not a CN statement."), 2),
        ("Sensitivity", "<p>Sensitivity analyses are part of the results above: estimand (A), scale (B, best_s &plusmn;3&ndash;5%) and denominator (C). The upper-trimmed-mean normalisation, grid limits and weights are fixed by the pinned code and were not varied.</p>"),
        ("Limits and alternative explanations",
         "<ul><li>The raw estimands lack GC correction, so estimand differences mix GC and aggregation effects; they are a sensitivity range, not a ranking of estimands.</li>"
         "<li>Boundaries are held at the saved segmentation; different estimands would segment differently, so these are conservative about call change.</li>"
         "<li>With best_s at the grid edge the scale is a free parameter of the code, not an estimate; the apparent agreement of cluster fits says little about copy-number levels.</li>"
         "<li>The 9p21 zone and MYC/19q13/5q regions are approximate literature extents, not validated coordinates; &lsquo;loss&rsquo; and &lsquo;gain&rsquo; are labels from the assignment rule, not CN truth.</li></ul>"),
    ]
    glossary = {
        "x_mean": "Mean of processed X over the cells of a group (current workbench pseudobulk).",
        "raw_ratio": "Sum of raw counts over cells divided by the sum of cell raw totals, times 1e5 (exposure-weighted rate).",
        "raw_cellmean": "Mean over cells of (raw / cell raw total) times 1e5; equal cell weights.",
        "best_s / best_s_x_mean": "Fitted scale of the weighted scale search: the segment mean treated as one copy; for x_mean and for the compared estimand.",
        "grid_lo / grid_hi": "Lower and upper limits of the s grid: 0.5 x median and min(1.2, 0.75 x median).",
        "at_lower_edge / at_lower_edge_n": "Whether (number of estimands for which) the fitted scale equals the grid lower limit.",
        "pct_bins_call_changed": "Percentage of all bins whose combined (Holmes/Watson) call differs between estimands.",
        "bins_loss_* / bins_gain_*": "Number of bins called loss / gain with x_mean (_x) or the alternative estimand (_alt).",
        "boundary_m_over_s": "Weighted decision boundary between states k and k+1 in units of segment mean / best_s.",
        "unweighted_midpoint": "Boundary with equal weights (k + 0.5).",
        "w_lo / w_hi": "State weights of the lower and upper adjacent states.",
        "best_s_factor / segments_changed / bins_changed": "Multiplier applied to best_s, and the segments / bins whose class changes.",
        "max_pct_changed / median_pct_changed": "Largest / median over removed chromosomes of the % of bins elsewhere whose call changes.",
        "worst_chromosome": "Removed chromosome producing max_pct_changed.",
        "best_s_range": "Range of best_s refitted after removing each chromosome.",
        "*_frac_loss / *_frac_gain": "Share of a region&rsquo;s bins called loss / gain with the named estimand (<code>all</code> source).",
        "*_mean_over_genome_median": "Region mean of the pseudobulk divided by the genome median of that estimand.",
        "source": "Pseudobulk group: WNN cluster 0&ndash;4 or all cells.",
        "harness_identical_to_saved_karyogram": "Whether the harness reproduces the saved cn_call and best_s for the source.",
    }
    C.write_html_report(HERE / "report.html", "Normalization and CN assignment", "Experiment 03 &middot; estimand, thresholds and denominator sensitivity with fixed segmentation",
                        sections, {"command": "python -I experiments/03_normalization_aggregation/run.py", "summary": "summary.json", "environment": C.environment()},
                        key_findings=findings, verdict=verdict, glossary=glossary,
                        supports=["Call states are sensitive to the pseudobulk estimand, to the position of best_s and to the set of chromosomes in the denominator.",
                                  "The CN scale fit sits at the grid edge for every source and estimand, so calls reduce to thresholds on segment mean / genome median.",
                                  "The chr7 dips are called loss under every estimand tested (consistent with a depth deficit relative to the genome median)."],
                        not_supports=["Which estimand or threshold is right: there is no CN truth here.",
                                      "That any called gain or loss is a copy-number state; the labels follow from the rule and the genome median.",
                                      "Behaviour on other samples or with a different reference."],
                        decision_impact="Supports the recommendation that the per-region evidence layer should report raw-count contrast with an explicit reference level, and that it must not inherit the fixed gain/loss labels as evidence. "
                                        "It also supports deferring any change to the scale search or weights until a cross-sample baseline exists, since the current rule cannot be validated here.",
                        next_steps=["Report the reference level (genome median or an external baseline) with every call.",
                                    "Re-run this comparison with GC-corrected raw-count estimands, to separate GC from aggregation effects.",
                                    "Test the assignment against a sample with independent CN truth before changing weights or the grid."])
    pd.set_option("display.width", 250)
    print(harness.to_string()); print(comp.to_string()); print(lt.to_string()); print(th.to_string()); print(s_sens.to_string()); print(loo_sum.to_string())
    print("best_s", bs, grid, "near cutpoint", int((margin < 0.1).sum()), len(mt))


if __name__ == "__main__":
    main()
