"""Experiment 19: why the AD binary scan misses moderate focal events (masking) and whether a two-cut scan recovers them.

Motivation: experiment 14 (semi-synthetic spike-in) and an independent coordinator check found that clean
0.5x dips of 30-50 bins in the cluster-0 pseudobulk are usually not segmented, while 0.1x dips are; raising
k to 10 did not help.

Design (prespecified):
  Events: clean multiplicative edits of the saved cluster-0 processed-X pseudobulk (no GC recomputation; an
  idealised version of exp 14's thinning), at the first 8 loci of experiment 14's prespecified locus lists
  (same scenario seeds, read from 14's realizations.tsv where present, else drawn with 14's rule),
  types loss_one (x0.5), gain (x1.5), loss_deep (x0.1); widths 20, 30, 50; f = 1 only (pseudobulk factor
  1 - f/2 etc. is equivalent for a mean).
  A. Baseline: workbench scan (server.allelic_scan.scan_chromosome, DEFAULTS, genome pool). For missed edges,
     classify the reason from the audit: the deepest tested node containing the event; was its top
     geometry-valid candidate elsewhere and rejected by the local test, by the global test, or did recursion
     stop for geometry/depth.
  B. Two-cut AD scan ('segment vs rest', the idea behind circular binary segmentation): inside the final
     baseline segment that contains the event, for every window [i, i+w) with 20 <= w <= 120 bins that leaves
     >= 20 bins outside, statistic T(i, w) = AD(window, rest of segment) with the pinned pyEpi dist_ad.
     Detection: max T over the segment exceeds a threshold AND the argmax window's edges are within +-2 bins
     of the injected edges. Threshold: 95th percentile of the same max statistic over unmodified final
     cluster-0 segments of length 150-1,500 bins (all such segments; real background, so it includes any real
     events and is conservative), scaled per segment length by fitting max T ~ a + b log(n).
  C. False positives: fraction of the unmodified segments above the threshold (5% by construction overall;
     reported by length) and, for the injected chromosomes, any new window called away from the event.
Deterministic (locus lists come from exp 14 seeds). Outputs: realizations.tsv, null_segments.tsv,
summary.json, report.html.
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

sys.path.insert(0, str(C.REPO))
from server import allelic_scan as sc  # noqa: E402
from server import pyepi_adapter as sci  # noqa: E402
from server import cn_recall as cr  # noqa: E402

WMIN, WMAX, MINREST, TOL = 20, 120, 20, 2
FACT = {"loss_one": 0.5, "gain": 1.5, "loss_deep": 0.1}
WIDTHS = [20, 30, 50]
NLOC = 8


def two_cut(v, true=None):
    """max over windows of AD(window, rest); returns (Tmax, i, w) or, with true=(i0, w0),
    (Tmax, i, w, T_true, fraction of windows with T >= T_true)."""
    n = len(v)
    best = (-np.inf, None, None)
    allT, t_true = [], np.nan
    for w in range(WMIN, min(WMAX, n - MINREST) + 1):
        for i in range(0, n - w + 1):
            if n - w < MINREST:
                continue
            win = v[i:i + w]
            rest = np.concatenate([v[:i], v[i + w:]])
            t = sci.dist_ad(win, rest)
            allT.append(t)
            if true is not None and i == true[0] and w == true[1]:
                t_true = t
            if t > best[0]:
                best = (t, i, w)
    if true is None:
        return best
    allT = np.asarray(allT)
    return best + (t_true, float(np.mean(allT >= t_true)) if np.isfinite(t_true) else np.nan)


def null_task(args):
    seg_vals, pkg, sid, n = args
    sci.load(pkg)
    t, i, w = two_cut(seg_vals)
    return {"segment": sid, "n_bins": n, "Tmax": t, "argmax_start": i, "argmax_width": w}


def real_task(args):
    (row, v, pool, pkg, ref_cuts) = args
    sci.load(pkg)
    s, e = row["s"], row["e"]
    v2 = v.copy(); v2[s:e] *= FACT[row["type"]]
    g = pool.copy(); g[row["chrom_idx"]] = v2
    out = sc.scan_chromosome(v2, np.ones(len(v2), bool), sc.DEFAULTS, pool=g)
    b = np.array([x["bin"] for x in out["breakpoints"]])
    near = lambda x: int(np.min(np.abs(b - x))) if len(b) else 10**6  # noqa: E731
    le, re_ = near(s), near(e)
    a = pd.DataFrame(out["audit"])
    tested = a[a.selected_for_significance == True]  # noqa: E712
    cont = tested[(tested.node_start_bin <= s) & (tested.node_end_bin >= e)].sort_values("recursion_depth").tail(1)
    reason = "detected" if (le <= TOL and re_ <= TOL) else "unknown"
    if reason != "detected" and len(cont):
        c = cont.iloc[0]
        at_event = min(abs(c.candidate_absolute_bin - s), abs(c.candidate_absolute_bin - e)) <= TOL
        if not c.accepted:
            if c.p_local >= sc.DEFAULTS["alpha"]:
                reason = "top cut " + ("at event" if at_event else "elsewhere") + "; rejected by local test"
            else:
                reason = "top cut " + ("at event" if at_event else "elsewhere") + "; rejected by global test"
        else:
            reason = "accepted cut elsewhere; recursion ended (depth/geometry)"
    # segment of the UNMODIFIED saved segmentation that contains the event (ref_cuts: saved cluster0 cuts)
    bounds = [0] + sorted(int(x) for x in ref_cuts) + [len(v2)]
    k = np.searchsorted(bounds, s, side="right") - 1
    lo, hi = bounds[k], bounds[k + 1]
    while hi < e and k + 2 < len(bounds):   # event spans a saved boundary: merge with the next segment
        k += 1; hi = bounds[k + 1]
    seg = v2[lo:hi]
    if len(seg) >= WMIN + MINREST and (e - s) <= WMAX:
        t, i, w, t_true, frac_ge = two_cut(seg, true=(s - lo, e - s))
    else:
        t, i, w, t_true, frac_ge = np.nan, None, None, np.nan, np.nan
    return {**{k_: row[k_] for k_ in ("type", "width", "locus", "chrom", "s", "e")}, "baseline_left_err": le, "baseline_right_err": re_,
            "baseline_detected": bool(le <= TOL and re_ <= TOL), "baseline_reason": reason,
            "segment_lo": lo, "segment_hi": hi, "segment_len": hi - lo, "two_cut_T": t,
            "two_cut_start": (lo + i) if i is not None else None, "two_cut_end": (lo + i + w) if i is not None else None,
            "two_cut_T_true_window": t_true, "frac_windows_T_ge_true": frac_ge}


VERDICT = DECISION = ""


def main():
    t0 = time.time()
    p = C.default_paths(); pkg = str(p["package"])
    sci.load(pkg)
    d = C.H5AD(p["h5ad"]); lab = C.cell_clusters(d, p)
    pb = C.group_means(d, "X", {"c0": lab == "0"})["c0"]
    bp = pd.read_csv(p["genome_pool"] / "raw_accepted_breakpoints.tsv", sep="\t")
    bp0 = bp[bp.source == "cluster0"]
    # null: unmodified final cluster0 segments 150-1500 bins
    nulls = []
    for ch in dict.fromkeys(d.seq):
        ci = np.flatnonzero(d.seq == ch)
        cuts = [0] + sorted(bp0.loc[bp0.chromosome == ch, "absolute_bin"].astype(int)) + [len(ci)]
        for a_, b_ in zip(cuts[:-1], cuts[1:]):
            if 150 <= b_ - a_ <= 1500:
                nulls.append((pb[ci[a_:b_]], pkg, f"{ch}:{a_}-{b_}", b_ - a_))
    with ProcessPoolExecutor(3) as ex:
        null = pd.DataFrame(list(ex.map(null_task, nulls)))
    null.to_csv(HERE / "null_segments.tsv", sep="\t", index=False)
    coef = np.polyfit(np.log(null.n_bins), null.Tmax, 1)
    resid = null.Tmax - np.polyval(coef, np.log(null.n_bins))
    q95 = float(np.quantile(resid, .95))
    thr = lambda n: float(np.polyval(coef, np.log(n)) + q95)  # noqa: E731
    # loci from exp 14 realizations (same prespecified lists), f = 1
    R14 = pd.read_csv(HERE.parent / "14_spikein_benchmark" / "realizations.tsv", sep="\t")
    R14 = R14[(R14.f == 1.0) & R14.type.isin(list(FACT)) & R14.width.isin(WIDTHS)]
    tasks = []
    for (typ, w), g in R14.groupby(["type", "width"]):
        for _, r in g.sort_values("locus").head(NLOC).iterrows():
            ci = np.flatnonzero(d.seq == r.chrom)
            row = {"type": typ, "width": int(w), "locus": int(r.locus), "chrom": r.chrom, "s": int(r.event_start_local),
                   "e": int(r.event_end_local_excl), "chrom_idx": ci}
            ref = bp0.loc[bp0.chromosome == r.chrom, "absolute_bin"].astype(int).tolist()
            tasks.append((row, pb[ci].copy(), pb, pkg, ref))
    with ProcessPoolExecutor(3) as ex:
        res = pd.DataFrame(list(ex.map(real_task, tasks)))
    res["two_cut_threshold"] = res.segment_len.apply(thr)
    res["two_cut_edges_ok"] = ((res.two_cut_start - res.s).abs() <= TOL) & ((res.two_cut_end - res.e).abs() <= TOL)
    res["two_cut_detected"] = (res.two_cut_T > res.two_cut_threshold) & res.two_cut_edges_ok
    res.to_csv(HERE / "realizations.tsv", sep="\t", index=False)
    null["threshold"] = null.n_bins.apply(thr); null["above"] = null.Tmax > null.threshold
    res["true_window_is_max"] = res.frac_windows_T_ge_true <= (1.0 / 1000)
    summ = res.groupby(["type", "width"]).agg(n=("locus", "size"), baseline_detected=("baseline_detected", "mean"),
                                             two_cut_detected=("two_cut_detected", "mean"),
                                             two_cut_edges_ok=("two_cut_edges_ok", "mean"),
                                             median_frac_windows_T_ge_true=("frac_windows_T_ge_true", "median")).reset_index()
    tt = res.assign(true_above=res.two_cut_T_true_window > res.two_cut_threshold).groupby(["type", "width"]).true_above.mean()
    summ["true_window_above_threshold"] = tt.values
    reasons = res[~res.baseline_detected].groupby(["type", "baseline_reason"]).size().rename("missed").reset_index()
    summ.to_csv(HERE / "detection_summary.tsv", sep="\t", index=False)
    reasons.to_csv(HERE / "miss_reasons.tsv", sep="\t", index=False)
    out = {"experiment": "19_masking_two_cut_scan", "deterministic": True, "n_null_segments": len(null), "threshold_fit": coef.tolist(),
           "threshold_resid_q95": q95, "null_frac_above": float(null.above.mean()), "detection": summ.to_dict("records"),
           "miss_reasons": reasons.to_dict("records"), "seconds": round(time.time() - t0, 1), "environment": C.environment()}
    C.write_json(HERE / "summary.json", out)

    # ---------------- report
    def fig(fn, w=7, h=4):
        f_, ax = plt.subplots(figsize=(w, h)); fn(ax); b = C.fig_to_b64(f_); plt.close(f_); return b

    def f_det(ax):
        for typ, c in [("loss_one", "#9b2c2c"), ("gain", "#2f5d8a"), ("loss_deep", "#000")]:
            g = summ[summ.type == typ]
            ax.plot(g.width, g.baseline_detected, "o--", color=c, label=f"{typ}: AD binary scan")
            ax.plot(g.width, g.two_cut_detected, "s-", color=c, label=f"{typ}: two-cut scan")
        ax.set_xlabel("event width (bins)"); ax.set_ylabel("detected (both edges ±2)"); ax.set_ylim(-0.05, 1.05); ax.legend(fontsize=6)

    def f_null(ax):
        ax.scatter(null.n_bins, null.Tmax, s=10, color="#bbb", label="unmodified segments")
        xs = np.linspace(null.n_bins.min(), null.n_bins.max(), 100); ax.plot(xs, [thr(x) for x in xs], color="#8a5a00", label="threshold")
        for typ, c in [("loss_one", "#9b2c2c"), ("gain", "#2f5d8a"), ("loss_deep", "#000")]:
            g = res[res.type == typ]; ax.scatter(g.segment_len, g.two_cut_T, s=18, color=c, label=typ)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("segment length (bins)"); ax.set_ylabel("max two-cut AD"); ax.legend(fontsize=6)
    b_det, b_null = fig(f_det), fig(f_null)
    lo1 = summ[summ.type == "loss_one"]
    better = float(lo1.two_cut_detected.mean()) > float(lo1.baseline_detected.mean()) + 0.2
    global VERDICT, DECISION
    if better:
        VERDICT = ("Moderate focal changes are missed mainly through masking in long segments; a segment-vs-rest (two-cut) scan inside the "
                   "same segment recovers many of them.")
        DECISION = ("Moves rank 4 (segmentation) forward: a two-cut or seeded scan should be benchmarked with exp 14's realistic injections and "
                    "an FDR-controlled threshold before adoption.")
    else:
        VERDICT = ("Moderate focal changes are missed not only because of the binary search: inside their own segment, natural windows on this real "
                   "background differ from the rest about as much as an injected 0.5× or 1.5× change, so a better search alone (two-cut scan) "
                   "does not isolate them. Detecting one-copy focal changes needs a baseline model of the regional structure. Masking is real for strong events: the two-cut scan recovers deep (0.1×) losses that the binary scan misses, because those windows are the most extreme in their segment.")
        DECISION = ("Shifts weight from segmentation search (rank 9/4) to the baseline (rank 3): without a cross-sample or modelled regional "
                    "baseline, neither the binary AD scan nor a two-cut AD scan separates one-copy focal events from accessibility structure at "
                    "this scale. Deep (0.1×) losses remain detectable by the current scan.")
    sm = summ.set_index(["type", "width"])
    key = [f"Baseline AD binary scan, clean clonal one-copy loss (×0.5): detected {sm.loc[('loss_one',20),'baseline_detected']:.0%} / "
           f"{sm.loc[('loss_one',30),'baseline_detected']:.0%} / {sm.loc[('loss_one',50),'baseline_detected']:.0%} at 20/30/50 bins; "
           f"two-cut scan {sm.loc[('loss_one',20),'two_cut_detected']:.0%} / {sm.loc[('loss_one',30),'two_cut_detected']:.0%} / {sm.loc[('loss_one',50),'two_cut_detected']:.0%}.",
           f"Gain (×1.5): baseline {sm.loc[('gain',20),'baseline_detected']:.0%} / {sm.loc[('gain',30),'baseline_detected']:.0%} / {sm.loc[('gain',50),'baseline_detected']:.0%}; "
           f"two-cut {sm.loc[('gain',20),'two_cut_detected']:.0%} / {sm.loc[('gain',30),'two_cut_detected']:.0%} / {sm.loc[('gain',50),'two_cut_detected']:.0%}.",
           f"Deep loss (×0.1): baseline {sm.loc[('loss_deep',20),'baseline_detected']:.0%} / {sm.loc[('loss_deep',30),'baseline_detected']:.0%} / {sm.loc[('loss_deep',50),'baseline_detected']:.0%}; "
           f"two-cut {sm.loc[('loss_deep',20),'two_cut_detected']:.0%} / {sm.loc[('loss_deep',30),'two_cut_detected']:.0%} / {sm.loc[('loss_deep',50),'two_cut_detected']:.0%}.",
           "Why the baseline misses: " + "; ".join(f"{r.baseline_reason} ({r.missed})" for r in reasons.groupby('baseline_reason').missed.sum().reset_index().itertuples()) + ".",
           f"Two-cut threshold from {len(null)} unmodified cluster-0 segments (150–1,500 bins): {null.above.mean():.0%} of them exceed it by construction (≈5%).",
           "Rank of the true event window among all windows of its unmodified containing segment (median fraction of windows with T ≥ T_true): " +
           "; ".join(f"{t} {w}: {v:.3f}" for (t, w), v in sm.median_frac_windows_T_ge_true.items()) + ". A value well above 0 means natural windows in the same segment differ more from the rest than the injected event does.",
           f"n = {NLOC} loci per cell; intervals are wide (Wilson 95% for 0/8 is 0–32%, for 8/8 is 68–100%)."]
    glossary = {"baseline_detected": "workbench AD scan places boundaries within ±2 bins of both injected edges",
                "baseline_reason": "for misses: what happened in the deepest tested node containing the event",
                "two_cut_T": "max over windows (20–120 bins) of AD(window, rest of the containing segment)",
                "two_cut_threshold": "95th percentile of the same max statistic in unmodified segments of that length (log-length fit)",
                "two_cut_edges_ok": "argmax window's edges within ±2 bins of the injected edges",
                "two_cut_detected": "T above threshold and edges ok", "segment_len": "length of the baseline final segment containing the event"}
    sections = [
        ("Question and why it matters", "<p>Experiment 14 and an independent check showed the workbench AD scan rarely segments a clean "
         "0.5× or 1.5× change of 20–50 bins in a 1,286-cell pseudobulk. If the cause is the top-down search rather than the data, the next "
         "change belongs in segmentation (recommendation stage 4), and it transfers to every sample.</p>"),
        ("Inputs", "<ul><li>Saved cluster-0 processed-X pseudobulk (1,286 cells), saved cluster-0 genome-pool breakpoints.</li>"
         "<li>Event loci: the first 8 prespecified loci per (type, width) of experiment 14 (f = 1).</li>"
         "<li>Events are clean multiplicative edits of the pseudobulk: idealised, no sampling noise added, no GC re-correction.</li></ul>"),
        ("Method", "<pre>" + __doc__.split("Design (prespecified):")[1].split("Deterministic")[0] + "</pre>"),
        ("Results", "", 2),
        ("Detection: binary scan vs two-cut scan", C.img(b_det, "detection") + C.caption("Fraction of loci with both edges recovered within ±2 bins. "
         "Dashed: workbench AD binary scan; solid: two-cut scan inside the containing segment.") + C.html_table(summ), 3),
        ("Why the binary scan misses", C.html_table(reasons) + C.caption("For each missed event, the deepest node that contained the event and "
         "had a candidate tested. 'Top cut elsewhere; rejected' means the node's best single cut was not at the event and failed a test, which "
         "stops the search in that node (masking)."), 3),
        ("Calibration of the two-cut scan", C.img(b_null, "null") + C.caption("Grey: max two-cut statistic of unmodified segments by length; "
         "line: threshold (log-length fit + 95th percentile of residuals); coloured: segments containing injected events."), 3),
        ("Per-realization detail", C.html_table(res.drop(columns=[]), max_rows=80)),
        ("Limits and alternative explanations", "<ul><li>Idealised events (exact multiplicative change of the mean profile, no extra sampling noise, "
         "no GC re-correction); experiment 14 has the realistic thinning version.</li><li>The two-cut threshold comes from real segments that can contain "
         "real events, so it is conservative; it is not a calibrated genome-wide FDR.</li><li>The two-cut scan here is applied inside the baseline's final "
         "segment; a full method would integrate it into the recursion (e.g. circular binary segmentation or seeded/wild binary segmentation) and needs "
         "its own multiplicity control.</li><li>Detectability is not interpretation: a recovered boundary still needs the evidence layer (experiments 08, 17).</li></ul>"),
    ]
    C.write_html_report(HERE / "report.html", "Masking and two-cut scan", "Experiment 19 · why moderate focal events are missed by binary AD segmentation",
                        sections, {"command": "python -I experiments/19_masking_two_cut_scan/run.py", "outputs": ["realizations.tsv", "null_segments.tsv",
                                   "detection_summary.tsv", "miss_reasons.tsv", "summary.json"], "environment": C.environment()},
                        verdict=VERDICT,
                        key_findings=key, glossary=glossary,
                        supports=["A segmentation change (two-cut/CBS-like or seeded binary segmentation) as the next method experiment.",
                                  "Reporting 'not segmented' separately from 'no change' for focal events below ~0.5x contrast."],
                        not_supports=["Any claim about real events in MCF10-A; events here are synthetic.", "A calibrated false-discovery rate for the two-cut scan."],
                        decision_impact=DECISION,
                        next_steps=["Integrate a two-cut statistic into the recursion and rerun experiment 14's realistic injections.",
                                    "Calibrate the two-cut threshold by permuting cells between halves (sampling-only null) plus the background-segment null.",
                                    "Compare with an exposure-aware count changepoint with a baseline-variance term (experiment 11 showed sampling-only dispersion over-segments)."])
    pd.set_option("display.width", 250)
    print(summ.to_string()); print(reasons.to_string()); print(out["null_frac_above"], out["seconds"])


if __name__ == "__main__":
    main()
