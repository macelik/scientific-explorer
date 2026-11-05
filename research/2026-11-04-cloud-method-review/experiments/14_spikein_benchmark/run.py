"""Experiment 14: semi-synthetic spike-in benchmark of the baseline pipeline (known injected events on real background).

PRESPECIFIED DESIGN (written before any run)
============================================
Purpose. NEXT_EXPERIMENTS.md section 3 asks for event-size-stratified evaluation. Here KNOWN events are
injected into the real raw counts of WNN cluster 0 (1,286 cells) and the unchanged baseline pipeline is
run on the result. This is a SIMULATION on a real background. It is NOT biological ground truth, not a
calibrated sensitivity, and accessibility structure is preserved by construction (thinning/adding scales
each bin's own observed counts), so it measures detectability on this background only.

Injection (raw counts, layers/raw, chromosome-local event bins [s, s+w)):
  loss_one   : each chosen cell's counts in the event bins -> Binomial(count, 0.5)   (one-copy loss model)
  loss_deep  : each chosen cell's counts in the event bins -> Binomial(count, 0.1)   (deep loss model)
  gain       : counts + Poisson(0.5 * shape_b * flank_rate_cell), flank_rate_cell = the cell's mean raw count
               over the 50 bins on each side of the event, shape_b = cluster0 raw-pseudobulk of bin b over the
               mean raw-pseudobulk of those flank bins (keeps the bin-specific accessibility profile) (+0.5x model)
  widths w in {5, 10, 20, 30, 50} bins; fraction f in {1.0, 0.5, 0.25} of cluster0 cells (random cells).
  Thinning-as-dosage is an assumption: it models a copy-number change as a proportional change of observed
  counts in the affected cells and nothing else.
Loci. Per scenario (type, w, f; index = enumeration order types x widths x f), a pre-drawn list of 40
  loci (seed tuple (20261104, 14, scenario), numpy PCG64) from autosome stretches with no internal gap
  (all 100 kb start steps inside the event), >= 60 bins from both chromosome ends, and >= 10 bins from every
  cluster0 AD breakpoint in raw_accepted_breakpoints.tsv. Per-realization cells and noise use
  (20261104, 14, scenario, locus). Tasks run locus-major so any time-limited run is balanced across
  scenarios; the first `--n-loci` of the 40 are used (see summary.json for how many were run).
Per realization:
  1. Recompute X for the modified cells with the pinned GC/round/clip/normalise steps
     (pyEpiAneufinder.py L252-L271; render_fragments.get_loess_smoothed span 0.75). X of the unmodified
     cells is the saved X. Because 11 cells do not reconstruct exactly (exp 01), the pseudobulk is updated as
     saved-X mean + (X_new - X_recon_unmodified)/n over the modified cells, so reconstruction error cancels.
  2. Cluster-0 pseudobulk = mean of X over cells. Rescan ONLY the affected chromosome with the workbench AD
     settings (server.allelic_scan.scan_chromosome, DEFAULTS k=6, 1000 perms, alpha 0.001, min_child 20,
     min_terminal 50; global pool = updated genome-wide pseudobulk).
  3. CN: pinned assign_gainloss_new (server.cn_recall.cluster_rows) on the full genome with the new
     boundaries on the affected chromosome and the saved boundaries elsewhere.
Scoring.
  left/right edge detected: some new boundary within +-2 bins of the injected edge (edges are the cuts at
  bin s and bin s+w). detected (primary) = both edges AND CN call in the injected direction (loss: call < 1,
  gain: call > 1; scale checked in summary.json) on >= 50% of event bins. Also reported: edge-wise rates,
  CN-direction-only rate, "net" CN (call changed relative to the unmodified call), boundary error in bins
  (signed, nearest boundary), false new boundaries (new boundaries on the chromosome with no reference
  boundary and no injected edge within +-2 bins), lost reference boundaries (reference boundary with no new
  boundary within +-2), and the measured raw-count rate ratio target/flanks (50 bins each side) of the
  cluster-0 raw pseudobulk after injection, compared with the unmodified ratio.
Harness check. The unmodified chromosome scans (all 22) must reproduce the saved cluster0 breakpoints, and
  unmodified CN calls the saved karyogram_cluster_segments.tsv.
Structural note (not a result): the workbench minimum child size is 20 bins, so events narrower than 20
  bins cannot be a segment by construction.
Reduced grid / determinism: the pinned permutation tests use their own fixed random_state=42, so a scan is a
  deterministic function of its input; the only randomness here is locus choice, cell choice and thinning/Poisson draws.
"""
import argparse
import os
import sys
import time
from multiprocessing import get_context
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
from server import cn_recall as cr  # noqa: E402
from server import pyepi_adapter as sci  # noqa: E402

SEED = C.SEED
TYPES = {"loss_one": ("loss", 0.5), "loss_deep": ("loss", 0.1), "gain": ("gain", 0.5)}
WIDTHS = [5, 10, 20, 30, 50]
FRACS = [1.0, 0.5, 0.25]
N_LOCI_PRESPEC = 40
FLANK = 50
EDGE_TOL = 2
G = {}  # shared (fork-inherited) read-only state


def scenarios():
    out = []
    for t in TYPES:
        for w in WIDTHS:
            for f in FRACS:
                out.append({"scenario": len(out), "type": t, "width": w, "f": f})
    return out


def draw_loci(sc_row, d_seq_chroms, bp_by_chrom):
    """40 prespecified loci for a scenario: list of (chrom, s_local)."""
    w = sc_row["width"]
    cand = []
    for ch, (idx, step_ok) in d_seq_chroms.items():
        n = len(idx)
        bps = np.array(bp_by_chrom.get(ch, []), int)
        for s in range(60, n - 60 - w + 1):
            e = s + w - 1
            if not step_ok[s:e].all():  # steps between consecutive event bins
                continue
            if len(bps):
                dist = np.where(bps <= s, s - bps, np.where(bps >= e + 1, bps - (e + 1), 0))
                if dist.min() < 10:
                    continue
            cand.append((ch, s))
    rng = np.random.default_rng([SEED, 14, sc_row["scenario"]])
    pick = rng.choice(len(cand), N_LOCI_PRESPEC, replace=False)
    return [cand[i] for i in pick], len(cand)


def corrected_x(counts, gc, rf):
    """Pinned pyEpiAneufinder.py L252-L271 for one cell (counts float64)."""
    import contextlib
    import io
    with contextlib.redirect_stderr(io.StringIO()):
        fit = rf.get_loess_smoothed(counts, gc)
    corr = counts.mean() / fit
    as_int = np.rint(counts * corr).astype(int)
    as_int[as_int < 0] = 0
    return as_int / (as_int.sum() / 1e5)


def init_worker():
    G["rf"] = C.load_pyepi_module("render_fragments", G["pkg"])
    sci.load(G["pkg"])


def recon_chunk(rng_):
    lo, hi = rng_
    return lo, np.vstack([corrected_x(G["raw0"][c].astype(np.float64), G["gc"], G["rf"]) for c in range(lo, hi)])


def ref_scan(ch):
    m = G["seq"] == ch
    out = sc.scan_chromosome(G["pb0"][m], np.ones(m.sum(), bool), sc.DEFAULTS, package_dir=G["pkg"], pool=G["pb0"])
    return ch, sorted(b["bin"] for b in out["breakpoints"])


def call_cn(pb, bounds):
    seg = cr.seg_labels(G["seq"], bounds)
    _, comb = cr.cluster_rows("cluster0", G["n"], pb, seg, G["seq"])
    return np.asarray(comb, float)


def realize(task):
    scn, locus, ch, s, deadline = task
    if time.time() > deadline:
        return None
    t0 = time.time()
    typ, w, f = scn["type"], scn["width"], scn["f"]
    kind, par = TYPES[typ]
    n = G["n"]
    cidx = G["chrom_idx"][ch]
    ev = cidx[s:s + w]
    fl = np.concatenate([cidx[s - FLANK:s], cidx[s + w:s + w + FLANK]])
    rng = np.random.default_rng([SEED, 14, scn["scenario"], locus])
    k = int(round(f * n))
    cells = np.sort(rng.choice(n, k, replace=False))
    R = G["raw0"][np.ix_(cells, ev)].astype(np.float64)
    if kind == "loss":
        Rn = rng.binomial(R.astype(np.int64), par).astype(np.float64)
    else:
        cellflank = G["raw0"][np.ix_(cells, fl)].astype(np.float64).mean(1)
        shape = G["pbraw0"][ev] / G["pbraw0"][fl].mean()
        Rn = R + rng.poisson(par * shape[None, :] * cellflank[:, None]).astype(np.float64)
    dx = np.zeros(G["nb"])
    for j, c in enumerate(cells):
        counts = G["raw0"][c].astype(np.float64)
        counts[ev] = Rn[j]
        dx += corrected_x(counts, G["gc"], G["rf"]) - G["xrecon0"][c]
    pb = G["pb0"] + dx / n
    pbraw = G["pbraw0"].copy()
    pbraw[ev] += (Rn - R).sum(0) / n
    ratio_before = G["pbraw0"][ev].mean() / G["pbraw0"][fl].mean()
    ratio_after = pbraw[ev].mean() / G["pbraw0"][fl].mean()
    xratio_before = G["pb0"][ev].mean() / G["pb0"][fl].mean()
    xratio_after = pb[ev].mean() / G["pb0"][fl].mean()
    # rescan the affected chromosome only
    m = G["seq"] == ch
    out = sc.scan_chromosome(pb[m], np.ones(m.sum(), bool), sc.DEFAULTS, package_dir=G["pkg"], pool=pb)
    newb = sorted(b["bin"] for b in out["breakpoints"])
    refb = G["ref_bounds"][ch]
    bounds = {c: list(v) for c, v in G["saved_bounds"].items()}
    bounds[ch] = newb
    comb = call_cn(pb, bounds)
    cn_ev = comb[ev]
    cn0 = G["comb0"][ev]
    direction = (cn_ev < 1) if kind == "loss" else (cn_ev > 1)
    net = (cn_ev < cn0) if kind == "loss" else (cn_ev > cn0)
    nb_arr = np.array(newb, int)

    def nearest(edge):
        return int(nb_arr[np.argmin(np.abs(nb_arr - edge))] - edge) if len(nb_arr) else None
    el, er = nearest(s), nearest(s + w)
    left = el is not None and abs(el) <= EDGE_TOL
    right = er is not None and abs(er) <= EDGE_TOL
    cn_ok = bool(direction.mean() >= 0.5)
    ref_arr = np.array(refb, int)
    false_new = [b for b in newb
                 if (not len(ref_arr) or np.abs(ref_arr - b).min() > EDGE_TOL) and min(abs(b - s), abs(b - (s + w))) > EDGE_TOL]
    lost_ref = [b for b in refb if not len(nb_arr) or np.abs(nb_arr - b).min() > EDGE_TOL]
    return {"scenario": scn["scenario"], "type": typ, "width": w, "f": f, "locus": locus, "chrom": ch,
            "event_start_local": s, "event_end_local_excl": s + w, "genomic_start": int(G["start"][ev[0]]),
            "n_cells_modified": k, "ratio_raw_before": ratio_before, "ratio_raw_after": ratio_after,
            "injected_contrast_raw": ratio_after / ratio_before, "ratio_X_before": xratio_before,
            "ratio_X_after": xratio_after, "injected_contrast_X": xratio_after / xratio_before,
            "n_boundaries_new": len(newb), "n_boundaries_ref": len(refb),
            "left_edge_error_bins": el, "right_edge_error_bins": er,
            "left_detected": bool(left), "right_detected": bool(right), "both_edges": bool(left and right),
            "cn_event_median": float(np.median(cn_ev)), "cn_baseline_median": float(np.median(cn0)),
            "cn_frac_in_direction": float(direction.mean()), "cn_direction": cn_ok,
            "cn_frac_net_change": float(net.mean()), "cn_net": bool(net.mean() >= 0.5),
            "detected": bool(left and right and cn_ok), "detected_net": bool(left and right and net.mean() >= 0.5),
            "n_false_new": len(false_new), "false_new_bins": ";".join(map(str, false_new)),
            "n_lost_ref": len(lost_ref), "seconds": round(time.time() - t0, 2)}


def wilson(k, n):
    from scipy.stats import binomtest
    if n == 0:
        return (np.nan, np.nan)
    ci = binomtest(int(k), int(n)).proportion_ci(0.95, method="wilson")
    return ci.low, ci.high


def summarise(df):
    rows = []
    for (t, w, f), g in df.groupby(["type", "width", "f"], sort=False):
        n = len(g)
        r = {"type": t, "width": w, "f": f, "n": n}
        for col, name in [("detected", "detected"), ("both_edges", "both_edges"), ("left_detected", "left_edge"),
                          ("right_detected", "right_edge"), ("cn_direction", "cn_direction_only"), ("detected_net", "detected_net")]:
            k = int(g[col].sum())
            lo, hi = wilson(k, n)
            r[name] = k / n
            if name in ("detected", "both_edges"):
                r[name + "_lo"], r[name + "_hi"] = lo, hi
        e = np.concatenate([g.left_edge_error_bins.dropna().abs().to_numpy(), g.right_edge_error_bins.dropna().abs().to_numpy()])
        r["median_abs_nearest_edge_err_bins"] = float(np.median(e)) if len(e) else np.nan
        det = g[g.left_detected]["left_edge_error_bins"].abs().tolist() + g[g.right_detected]["right_edge_error_bins"].abs().tolist()
        r["mean_abs_err_when_edge_within_tol"] = float(np.mean(det)) if det else np.nan
        r["mean_false_new_boundaries"] = float(g.n_false_new.mean())
        r["frac_with_false_new"] = float((g.n_false_new > 0).mean())
        r["mean_lost_ref_boundaries"] = float(g.n_lost_ref.mean())
        r["median_injected_contrast_raw"] = float(g.injected_contrast_raw.median())
        r["median_contrast_X"] = float(g.injected_contrast_X.median())
        rows.append(r)
    return pd.DataFrame(rows)


def figures(df, summ):
    figs = {}
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    cols = {1.0: "#2f5d8a", 0.5: "#c27c0e", 0.25: "#7a3b8a"}
    for ax, t in zip(axs, TYPES):
        for f in FRACS:
            s = summ[(summ.type == t) & (summ.f == f)].sort_values("width")
            if not len(s):
                continue
            ax.errorbar(s.width, s.detected, yerr=[s.detected - s.detected_lo, s.detected_hi - s.detected],
                        marker="o", capsize=3, color=cols[f], label=f"f={f}")
        ax.axvline(20, color="grey", ls=":", lw=1)
        ax.set_title(t)
        ax.set_xlabel("event width (100 kb bins)")
        ax.set_xticks(WIDTHS)
        ax.set_ylim(-0.03, 1.03)
    axs[0].set_ylabel("detected: both edges +-2 bins and CN in direction")
    axs[0].legend(fontsize=8)
    figs["width"] = C.fig_to_b64(fig)
    plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    big = df[df.width >= 20]
    rng = np.random.default_rng(0)
    for ax, t in zip(axs, TYPES):
        g = big[big.type == t]
        ax.scatter(g.injected_contrast_raw, g.detected.astype(float) + rng.uniform(-.06, .06, len(g)),
                   c=g.width, cmap="viridis", s=14, alpha=.7)
        ax.axvspan(0.44, 1.72, color="grey", alpha=.12)
        ax.set_xscale("log")
        ax.set_title(t + " (width >= 20)")
        ax.set_xlabel("measured raw ratio target/flank, after / before")
    axs[0].set_ylabel("detected (jittered 0/1)")
    figs["contrast"] = C.fig_to_b64(fig)
    plt.close(fig)
    return figs


def tol_sensitivity(df):
    rows = []
    for tol in (0, 1, 2, 5, 10):
        for (t, w), g in df[df.width >= 20].groupby(["type", "width"], sort=False):
            l = g.left_edge_error_bins.abs() <= tol
            r = g.right_edge_error_bins.abs() <= tol
            rows.append({"tolerance_bins": tol, "type": t, "width": w, "n": len(g), "both_edges_and_cn": float((l & r & g.cn_direction).mean()),
                         "both_edges": float((l & r).mean())})
    return pd.DataFrame(rows)


def report_figures(df, summ):
    figs = figures(df, summ)  # width + contrast
    cols = {"loss_one": "#2f5d8a", "loss_deep": "#c27c0e", "gain": "#7a3b8a"}
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for ax, w in zip(axs, [20, 30, 50]):
        for t in TYPES:
            s = summ[(summ.type == t) & (summ.width == w)].sort_values("f")
            ax.errorbar(s.f, s.detected, yerr=[s.detected - s.detected_lo, s.detected_hi - s.detected], marker="o", capsize=3, color=cols[t], label=t)
        ax.set_title(f"width {w}"); ax.set_xlabel("fraction of cells carrying event"); ax.set_xticks(FRACS); ax.set_ylim(-.03, 1.03)
    axs[0].set_ylabel("detected"); axs[0].legend(fontsize=8)
    figs["frac"] = C.fig_to_b64(fig); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    pooled = df[df.width >= 20].groupby("type").detected.mean()
    ax.bar(list(pooled.index), pooled.values, color=[cols[t] for t in pooled.index]); ax.set_ylim(0, 1)
    ax.set_ylabel("detected, pooled over widths 20-50 and f"); 
    figs["type"] = C.fig_to_b64(fig); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 3.4))
    big = df[df.width >= 20]
    errs = [np.concatenate([big[big.width == w].left_edge_error_bins.dropna().to_numpy(), big[big.width == w].right_edge_error_bins.dropna().to_numpy()]) for w in (20, 30, 50)]
    ax.boxplot(errs, labels=["20", "30", "50"], showfliers=True, flierprops={"markersize": 2})
    ax.axhspan(-EDGE_TOL, EDGE_TOL, color="grey", alpha=.15); ax.set_ylim(-60, 60)
    ax.set_xlabel("event width (bins)"); ax.set_ylabel("signed offset of nearest boundary (bins)")
    figs["err"] = C.fig_to_b64(fig); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 3.4))
    g = df.groupby(["width"]).n_false_new.mean()
    ax.bar([str(x) for x in g.index], g.values, color="#8a5a00"); ax.set_xlabel("event width (bins)"); ax.set_ylabel("mean false new boundaries / realization")
    figs["false"] = C.fig_to_b64(fig); plt.close(fig)
    return figs


def write_report(out, df, summ, summary):
    figs = report_figures(df, summ)
    ts = tol_sensitivity(df)
    ts.to_csv(out / "tolerance_sensitivity.tsv", sep="\t", index=False)
    n = summary["source"]["cells"]
    per = summary["grid"]
    h = summary["harness"]
    cn_values = summary["cn_state_values_baseline"]
    big = df[df.width >= 20]
    small = df[df.width < 20]

    def cell(typ, w, f):
        r = summ[(summ.type == typ) & (summ.width == w) & (summ.f == f)]
        if not len(r):
            return "n/a"
        r = r.iloc[0]
        return f"{r.detected:.2f} (95% CI {r.detected_lo:.2f}-{r.detected_hi:.2f}, n={int(r.n)})"
    errs = pd.concat([big.left_edge_error_bins, big.right_edge_error_bins]).dropna()
    within = errs[errs.abs() <= 10]
    fn_all = df.n_false_new.mean()
    verdict = (f"On this cluster-0 background the baseline pipeline recovered injected events only when they were at least 20 bins wide and carried by most cells "
               f"(pooled detection at widths 20-50: {big.detected.mean():.2f}; widths 5-10: {small.detected.mean():.2f}); this is a semi-synthetic detectability "
               f"measurement, not biological sensitivity.")
    key = [
        f"Width 50, all cells (f=1.0): detection loss_one {cell('loss_one', 50, 1.0)}; loss_deep {cell('loss_deep', 50, 1.0)}; gain {cell('gain', 50, 1.0)}.",
        f"Width 50, 25% of cells (f=0.25): loss_one {cell('loss_one', 50, 0.25)}; loss_deep {cell('loss_deep', 50, 0.25)}; gain {cell('gain', 50, 0.25)}.",
        f"Widths 5 and 10 bins: pooled detection {small.detected.mean():.2f} ({int(small.detected.sum())}/{len(small)}); the workbench min-child rule (20 bins) forbids segments this short, so this is structural.",
        f"Boundary accuracy at widths >=20 (nearest boundary within 10 bins of an edge, n={len(within)} edges of {len(errs)} with any boundary): median |offset| {within.abs().median():.1f} bins, "
        f"fraction within +-2 bins {(errs.abs() <= 2).mean():.2f}.",
        f"False new boundaries on the rescanned chromosome: mean {fn_all:.2f} per realization; {(df.n_false_new > 0).mean():.2f} of realizations had at least one; mean reference boundaries lost {df.n_lost_ref.mean():.2f}.",
        f"Harness check: {h['scans_identical_to_saved']}/{h['chromosomes']} unmodified chromosome scans reproduced the saved cluster0 breakpoints; CN calls identical to saved karyogram: {h['cn_calls_identical_to_saved_karyogram']}; X reconstruction exact for {n - h['x_reconstruction_mismatching_cells']}/{n} cells.",
        f"Measured injected contrast (raw target/flank, after/before) at f=1.0 ranged {df[df.f == 1.0].injected_contrast_raw.min():.2f}-{df[df.f == 1.0].injected_contrast_raw.max():.2f}; "
        f"at f=0.25 {df[df.f == 0.25].injected_contrast_raw.min():.2f}-{df[df.f == 0.25].injected_contrast_raw.max():.2f}, inside the natural 0.44-1.72x range of uncalled windows (experiment 04) for most loss realizations.",
        f"Ran {len(df)} realizations ({per['loci_used_per_scenario_min']}-{per['loci_used_per_scenario_max']} loci of the 40 prespecified per scenario, 45 scenarios), so intervals are wide.",
    ]
    inputs = ("<ul><li><code>pyepi_results/count_matrix.h5ad</code>: <code>layers/raw</code>, <code>X</code>, <code>var.GC</code> (read-only).</li>"
              "<li>WNN cluster 0 labels via <code>common.cell_clusters</code> (integrated.h5mu): " + str(n) + " cells.</li>"
              "<li><code>integration-story/data/genome_pool/raw_accepted_breakpoints.tsv</code> (saved cluster0 AD breakpoints, chromosome-local) and "
              "<code>derived/karyogram_cluster_segments.tsv</code> (harness).</li>"
              "<li>Pinned pyEpi e8305ea: <code>render_fragments.get_loess_smoothed</code>, <code>assign_gainloss_new</code>, AD/permutation tests; workbench <code>server/allelic_scan.py</code>, <code>server/cn_recall.py</code>.</li>"
              "<li>Loci rule: autosomes, gap-free event (all 100 kb start steps), &ge;60 bins from chromosome ends, &ge;10 bins from any saved cluster0 breakpoint.</li>"
              "<li><strong>Not available:</strong> independent copy-number truth, allelic information in this simulation, any real event with known cell fraction; only the raw-count model is perturbed.</li></ul>")
    design = ("<p>Chromosome-local event bins [s, s+w); chosen cells = random fraction f of the cluster (seed tuple (20261104, 14, scenario, locus)). "
              "R = raw counts of a chosen cell in the event bins; flank = 50 bins each side.</p><pre>"
              "loss_one :  R' ~ Binomial(R, 0.5)\nloss_deep:  R' ~ Binomial(R, 0.1)\n"
              "gain     :  R' = R + Poisson(0.5 * shape_b * flank_rate_cell)\n"
              "            flank_rate_cell = mean raw count of the cell over the flank bins\n"
              "            shape_b = pb_raw[b] / mean(pb_raw[flank])   (cluster-0 raw pseudobulk)\n"
              "X' (modified cells) = pinned L252-L271: loess(GC, counts) -> counts*mean/fit -> rint -> clip<0 -> /(sum/1e5)\n"
              "pb' = pb_savedX + sum_cells(X' - X_recon_unmodified)/n\n"
              "detected = boundary within +-2 bins of BOTH edges AND CN call (loss <1, gain >1) on >=50% of event bins</pre>"
              "<p>Grid: 3 types x widths {5,10,20,30,50} x f {1.0,0.5,0.25} = 45 scenarios. 40 loci per scenario are prespecified; "
              f"<strong>{per['loci_used_per_scenario_min']}-{per['loci_used_per_scenario_max']} were run per scenario</strong> (loci taken in order, locus-major) because each realization needs LOESS for up to "
              "1,286 cells plus a chromosome AD scan (~10-25 s), and the budget was 45 minutes on 2 workers. Nothing else in the grid was dropped. "
              "Only the affected chromosome is rescanned; CN is re-called genome-wide using saved boundaries elsewhere.</p>"
              "<p>Seeds: loci (20261104, 14, scenario); cells and noise (20261104, 14, scenario, locus), numpy PCG64. The pinned permutation tests use their own fixed random_state=42, so scans are deterministic given input.</p>")
    t_w = summ[["type", "f", "width", "n", "detected", "detected_lo", "detected_hi", "left_edge", "right_edge", "cn_direction_only"]].sort_values(["type", "f", "width"])
    t_f = summ[summ.width >= 20][["type", "width", "f", "n", "detected", "detected_lo", "detected_hi", "median_injected_contrast_raw"]].sort_values(["type", "width", "f"])
    pt = big.groupby("type").agg(n=("detected", "size"), detected=("detected", "mean"), both_edges=("both_edges", "mean"), cn_direction=("cn_direction", "mean"),
                                 cn_net=("cn_net", "mean"), detected_net=("detected_net", "mean")).reset_index()
    eb = pd.concat([big.assign(edge="left", off=big.left_edge_error_bins), big.assign(edge="right", off=big.right_edge_error_bins)]).dropna(subset=["off"])
    ea = eb.groupby(["type", "width"]).agg(edges=("off", "size"), median_signed=("off", "median"), median_abs=("off", lambda x: x.abs().median()),
                                            frac_within_2=("off", lambda x: (x.abs() <= 2).mean()), frac_within_10=("off", lambda x: (x.abs() <= 10).mean())).reset_index()
    fa = df.groupby(["type", "width"]).agg(n=("n_false_new", "size"), mean_false_new=("n_false_new", "mean"), frac_any_false=("n_false_new", lambda x: (x > 0).mean()),
                                           mean_lost_ref=("n_lost_ref", "mean")).reset_index()
    big2 = big.copy()
    big2["contrast_bin"] = pd.cut(big2.injected_contrast_raw, [0, 0.44, 0.7, 0.9, 1.1, 1.3, 1.72, 10])
    ca = big2.groupby(["type", "contrast_bin"], observed=True).agg(n=("detected", "size"), detected=("detected", "mean")).reset_index()
    ca["contrast_bin"] = ca.contrast_bin.astype(str)
    ts_show = ts.groupby(["tolerance_bins", "type"]).apply(lambda g: pd.Series({"both_edges_and_cn": np.average(g.both_edges_and_cn, weights=g.n), "n": g.n.sum()})).reset_index()
    T = lambda x: C.html_table(x, "{:.3f}")  # noqa: E731
    results = [
        ("Detection by width", T(t_w) + C.img(figs["width"], "detection vs width") +
         C.caption("Points: fraction of realizations meeting the primary detection rule (both edges within +-2 bins and CN in direction); bars: Wilson 95% CI; colours: cell fraction f; dotted line: 20-bin minimum child size. Read left to right: detection requires width above the minimum segment size and rises with width."), 3),
        ("Detection by cell fraction", T(t_f) + C.img(figs["frac"], "detection vs cell fraction") +
         C.caption("Same detection measure against the fraction of cells carrying the event, one panel per width. Steeper decline to the left means the event needs a larger carrier fraction to leave a pseudobulk step."), 3),
        ("Loss versus gain", T(pt) + C.img(figs["type"], "detection by type") +
         C.caption("Pooled over widths 20-50 and all f. detected_net requires the CN call to change versus the unmodified call; cn_net/cn_direction are CN-only rates. Pooled values mix scenarios of very different difficulty, so use them for ordering only; the per-scenario table above is the primary result."), 3),
        ("Boundary accuracy", T(ea) + C.img(figs["err"], "boundary offset") +
         C.caption("Signed offset (bins) of the nearest new boundary from each injected edge, widths >=20, all realizations with any boundary on the chromosome; grey band: +-2 tolerance. Large offsets mean the nearest boundary belongs to something other than this event."), 3),
        ("False boundaries", T(fa) + C.img(figs["false"], "false boundaries") +
         C.caption("New boundaries on the rescanned chromosome with no reference boundary and no injected edge within +-2 bins (mean per realization, by width). Lost reference boundaries are reference boundaries no longer found within +-2 bins."), 3),
        ("Detection versus measured raw contrast", T(ca) + C.img(figs["contrast"], "detection vs measured contrast") +
         C.caption("Each point is one realization (widths >=20, jittered 0/1 detected) against the measured raw-pseudobulk target/flank ratio after injection divided by before. Grey band 0.44-1.72x is the span of target/flank ratios of uncalled random 25-bin windows (experiment 04): contrasts inside it are within natural background variability, so detection there is not expected to be reliable."), 3),
        ("Files", "<p>Full tables: <code>realizations.tsv</code> (one row per realization), <code>detection_summary.tsv</code> (width x f x type), <code>tolerance_sensitivity.tsv</code>, <code>summary.json</code>.</p>", 3),
    ]
    sections = [("Question and why it matters",
                 "<p>NEXT_EXPERIMENTS.md section 3 needs event-size-stratified evaluation: the app's width series describes specific realizations, not detection probabilities. "
                 "This benchmark injects <em>known</em> events into real cluster-0 raw counts and asks how the current baseline pipeline (AD recursive segmentation, k=6, 1000 permutations, alpha 0.001, "
                 "min child 20 / terminal 50 bins; pinned <code>assign_gainloss_new</code>) responds, as a function of width, depth, carrier-cell fraction and local baseline.</p>"
                 "<p class='warn'>Semi-synthetic simulation on a real background; not biological ground truth and not a calibrated sensitivity.</p>"),
                ("Inputs", inputs), ("Design", design), ("Results", "<p>Each subsection gives a focused table and figure; full tables are in the TSV files listed at the end.</p>"),
                ] + results + [
                ("Sensitivity", "<p>Edge tolerance (widths &ge;20; both edges within tolerance and CN in direction), pooled over f and widths per type:</p>" + T(ts_show) +
                 C.caption("A tolerance of 0 bins demands exact edges; 10 bins is lenient. If detection climbs strongly with tolerance, boundaries are found but displaced. Primary analysis uses 2 bins (prespecified). Alternative detection definitions (detected_net, cn_direction_only) are in the detection tables."), 3),
                ("Limits and alternative explanations",
                 "<ul><li>Semi-synthetic: thinning and Poisson addition model a dosage change as a proportional change in observed counts; real events also change accessibility, GC and mappability structure, which is preserved here by construction. This measures detectability on this background only.</li>"
                 "<li>Conditional on random loci, cluster 0 (1,286 cells) and this subset of loci; not a calibrated sensitivity or false-positive rate.</li>"
                 "<li>Loci avoid chromosome ends by 60 bins and existing breakpoints by 10 bins, which favours detection compared with hard cases.</li>"
                 "<li>Detection at widths below 20 is structurally impossible under the min-child rule, an alternative explanation for their zero rate.</li>"
                 "<li>Baseline CN calls already deviate from 1 in some regions; detected_net separates this. Only the affected chromosome is rescanned.</li>"
                 "<li>The gain model's Poisson noise is not overdispersed like real cells.</li></ul>")]
    supports = ["Whether this pseudobulk-AD pipeline can place both edges of injected events as a function of width and carrier fraction, on this background.",
                "That boundary offsets and extra boundaries can be quantified per scenario against a known answer."]
    not_sup = ["Biological sensitivity or specificity, or the true fraction of cells in a real subclone.",
               "Performance for events whose accessibility also changes, or on other clusters or samples.",
               "That a failed detection means the event is absent in real data."]
    decision = ("<p>For research/2026-11-04-cloud-method-review/recommendation.md: the baseline can be described as able to recover events of at least "
                "20 bins when most cells carry them and the pseudobulk contrast leaves the natural range, and as unable by construction to return segments shorter than 20 bins. "
                "Detection at low carrier fractions is given by the tables above (conditional on this simulation); statements about small or subclonal events in the app should cite these rates as exploratory, not as calibrated detection probabilities.</p>")
    nxt = ["Run the remaining prespecified loci (up to 40 per scenario) to narrow the intervals.",
           "Repeat on another cluster and on the all-cell pseudobulk, and with a lower min-child setting as an exploratory alternative.",
           "Add an accessibility-altering injection arm and per-cell overdispersed gain noise.",
           "Rescan the whole genome for false boundaries elsewhere, not only the affected chromosome."]
    gl = {"type": "injection model: loss_one (thin x0.5), loss_deep (thin x0.1), gain (+0.5x)", "width": "event width in 100 kb bins", "f": "fraction of cluster-0 cells carrying the event",
          "n": "number of realizations", "detected": "fraction with both edges within +-2 bins AND CN call in injected direction on >=50% of event bins", "detected_lo / detected_hi": "Wilson 95% interval of detected",
          "both_edges": "fraction with both edges within +-2 bins (CN ignored)", "left_edge / right_edge": "fraction with a boundary within +-2 bins of that edge", "cn_direction_only / cn_direction": "CN call in injected direction on >=50% of event bins (boundaries ignored)",
          "detected_net / cn_net": "as detected / cn_direction but requiring the call to differ from the unmodified call in the injected direction", "left_edge_error_bins / right_edge_error_bins": "signed offset of the nearest new boundary from the injected edge",
          "median_signed / median_abs": "median signed / absolute edge offset in bins", "frac_within_2 / frac_within_10": "fraction of edges whose nearest boundary is within 2 / 10 bins",
          "mean_abs_err_when_edge_within_tol": "mean absolute offset among edges that are within tolerance", "n_false_new / mean_false_new / mean_false_new_boundaries": "new boundaries with no reference boundary or injected edge within +-2 bins",
          "frac_with_false_new / frac_any_false": "fraction of realizations with at least one false new boundary", "n_lost_ref / mean_lost_ref / mean_lost_ref_boundaries": "reference boundaries with no new boundary within +-2 bins",
          "injected_contrast_raw / median_injected_contrast_raw": "raw-pseudobulk target/flank ratio after injection divided by before (50-bin flanks)", "ratio_raw_before / ratio_raw_after": "target/flank raw-pseudobulk ratio before / after injection",
          "ratio_X_before / ratio_X_after / injected_contrast_X": "same ratios on the processed-X pseudobulk", "median_contrast_X": "median injected_contrast_X",
          "median_abs_nearest_edge_err_bins": "median absolute offset of the nearest boundary over both edges", "tolerance_bins": "edge tolerance used in the sensitivity table",
          "both_edges_and_cn": "primary detection rule at the stated tolerance", "contrast_bin": "bin of injected_contrast_raw", "cn_event_median / cn_baseline_median": "median CN call over event bins after injection / unmodified",
          "cn_frac_in_direction / cn_frac_net_change": "fraction of event bins with CN call in direction / changed versus unmodified", "n_cells_modified": "cells carrying the event",
          "locus": "locus index within scenario (0-39)", "event_start_local / event_end_local_excl": "chromosome-local bin range [start, end)", "n_boundaries_new / n_boundaries_ref": "boundaries on the chromosome after injection / unmodified",
          "seconds": "worker time for the realization"}
    C.write_html_report(out / "report.html", "Spike-in benchmark of the baseline pipeline",
                        "Experiment 14 · known injected events on real cluster-0 raw counts · semi-synthetic, not biological ground truth",
                        sections, {k: summary[k] for k in ["seed", "seed_streams", "grid", "harness", "runtime_seconds", "environment", "command"]},
                        key_findings=key, verdict=verdict, glossary=gl, supports=supports, not_supports=not_sup, decision_impact=decision, next_steps=nxt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--n-loci", type=int, default=N_LOCI_PRESPEC)
    ap.add_argument("--limit", type=int, default=0, help="bounded check: run only the first N tasks")
    ap.add_argument("--only-scenarios", default="", help="comma list of scenario indices")
    ap.add_argument("--report-only", action="store_true", help="rebuild tables/figures/report from existing realizations.tsv + summary.json")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--max-minutes", type=float, default=1e9, help="stop starting new realizations after this")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.report_only:
        import json
        df = pd.read_csv(out / "realizations.tsv", sep="\t")
        summ = summarise(df)
        summ.to_csv(out / "detection_summary.tsv", sep="\t", index=False)
        write_report(out, df, summ, json.loads((out / "summary.json").read_text()))
        return
    t_start = time.time()
    p = C.default_paths()
    pkg = str(p["package"])
    d = C.H5AD(p["h5ad"])
    lab = C.cell_clusters(d, p)
    c0 = np.flatnonzero(lab == "0")
    raw0 = np.vstack([d.rows("layers/raw", int(i), int(i) + 1)[0] for i in c0]).astype(np.int32)
    x0 = np.vstack([d.rows("X", int(i), int(i) + 1)[0] for i in c0])
    seq, start, gc = d.seq.copy(), d.start.copy(), d.gc.copy()
    nb = d.n_bins
    d.f.close()
    del d
    chroms = list(dict.fromkeys(seq))
    chrom_idx = {c: np.flatnonzero(seq == c) for c in chroms}
    bp = pd.read_csv(p["genome_pool"] / "raw_accepted_breakpoints.tsv", sep="\t")
    saved_bounds = cr.boundaries_from_bp(bp, ["cluster0"], chroms)["cluster0"]
    n = len(c0)
    G.update(pkg=pkg, raw0=raw0, gc=gc, seq=seq, start=start, nb=nb, n=n, chrom_idx=chrom_idx, saved_bounds=saved_bounds,
             pb0=x0.mean(0), pbraw0=raw0.astype(np.float64).mean(0))
    ctx = get_context("fork")

    # ---- unmodified reconstruction of X for the cluster-0 cells (pinned steps), for exact pseudobulk deltas
    G["rf"] = C.load_pyepi_module("render_fragments", pkg)
    chunks = [(lo, min(n, lo + 40)) for lo in range(0, n, 40)]
    xrecon0 = np.zeros((n, nb))
    with ctx.Pool(a.workers, initializer=init_worker) as pool:
        for lo, blk in pool.imap_unordered(recon_chunk, chunks):
            xrecon0[lo:lo + len(blk)] = blk
    exact = bool(np.array_equal(xrecon0, x0))
    cell_mismatch = int((np.abs(xrecon0 - x0).max(1) > 1e-9).sum())
    G["xrecon0"] = xrecon0
    t_recon = time.time() - t_start

    # ---- reference scans (all chromosomes) + harness check vs saved
    with ctx.Pool(a.workers, initializer=init_worker) as pool:
        ref = dict(pool.map(ref_scan, chroms))
    G["ref_bounds"] = ref
    harness = {c: {"saved": saved_bounds[c], "rerun": ref[c], "identical": saved_bounds[c] == ref[c]} for c in chroms}
    sci.load(pkg)
    comb0 = call_cn(G["pb0"], saved_bounds)
    G["comb0"] = comb0
    saved_seg = pd.read_csv(p["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    sref = saved_seg[saved_seg.group == "cluster0"].reset_index(drop=True)
    rows0, _ = cr.cluster_rows("cluster0", n, G["pb0"], cr.seg_labels(seq, saved_bounds), seq)
    got = pd.DataFrame(rows0)
    cn_ident = bool(len(sref) == len(got) and np.allclose(sref.cn_call, got.cn_call) and np.allclose(sref.best_s, got.best_s))
    cn_values = sorted(set(np.round(comb0[np.isfinite(comb0)], 3).tolist()))
    print("harness: scans identical", sum(v["identical"] for v in harness.values()), "/", len(chroms), "cn identical", cn_ident,
          "x recon exact", exact, "mismatching cells", cell_mismatch, "cn states", cn_values, flush=True)

    # ---- loci
    d_seq_chroms = {c: (chrom_idx[c], np.diff(start[chrom_idx[c]]) == 100000) for c in chroms}
    scs = scenarios()
    if a.only_scenarios:
        keep = {int(x) for x in a.only_scenarios.split(",")}
        scs = [s for s in scs if s["scenario"] in keep]
    loci, ncand = {}, {}
    for s_ in scs:
        loci[s_["scenario"]], ncand[s_["scenario"]] = draw_loci(s_, d_seq_chroms, saved_bounds)
    deadline = t_start + a.max_minutes * 60
    done_path = out / "realizations.tsv"
    done = pd.read_csv(done_path, sep="\t") if done_path.exists() else pd.DataFrame()
    done_keys = set(zip(done.scenario, done.locus)) if len(done) else set()
    tasks = []
    for li in range(a.n_loci):
        for s_ in scs:
            if (s_["scenario"], li) in done_keys:
                continue
            ch, s = loci[s_["scenario"]][li]
            tasks.append((s_, li, ch, s, deadline))
    if a.limit:
        tasks = tasks[:a.limit]
    print(f"{len(tasks)} realizations to run ({len(done_keys)} already done)", flush=True)
    t_run = time.time()
    with ctx.Pool(a.workers, initializer=init_worker) as pool:
        for i, r in enumerate(pool.imap_unordered(realize, tasks, chunksize=1)):
            if r is None:
                continue
            exists = done_path.exists()
            pd.DataFrame([r]).to_csv(done_path, sep="\t", index=False, mode="a" if exists else "w", header=not exists)
            if i % 10 == 0:
                print(f"  {i + 1}/{len(tasks)} {time.time() - t_run:.0f}s", flush=True)
    df = pd.read_csv(done_path, sep="\t").sort_values(["scenario", "locus"]).reset_index(drop=True)
    df.to_csv(done_path, sep="\t", index=False)
    summ = summarise(df)
    summ.to_csv(out / "detection_summary.tsv", sep="\t", index=False)
    runtime = time.time() - t_start
    env = C.environment()
    per_sc = df.groupby("scenario").size()
    summary = {
        "experiment": "14_spikein_benchmark", "semi_synthetic": True, "not_biological_ground_truth": True,
        "seed": SEED, "seed_streams": {"loci": [SEED, 14, "scenario"], "cells_and_noise": [SEED, 14, "scenario", "locus"]},
        "generator": "numpy.random.default_rng (PCG64)",
        "pinned_permutation_random_state": 42,
        "grid": {"types": TYPES, "widths": WIDTHS, "fractions": FRACS, "n_scenarios": 45, "loci_prespecified_per_scenario": N_LOCI_PRESPEC,
                 "loci_used_per_scenario_min": int(per_sc.min()), "loci_used_per_scenario_max": int(per_sc.max()),
                 "n_realizations": len(df), "candidate_starts_per_scenario": ncand},
        "source": {"cluster": "0", "cells": n, "bins": nb}, "flank_bins": FLANK, "edge_tolerance_bins": EDGE_TOL,
        "scan_settings": sc.DEFAULTS, "cn_state_values_baseline": cn_values,
        "harness": {"scans_identical_to_saved": int(sum(v["identical"] for v in harness.values())), "chromosomes": len(chroms),
                    "nonidentical": {c: v for c, v in harness.items() if not v["identical"]},
                    "cn_calls_identical_to_saved_karyogram": cn_ident,
                    "x_reconstruction_exact_all_cluster0_cells": exact, "x_reconstruction_mismatching_cells": cell_mismatch},
        "reference_boundaries": ref,
        "overall_detected": float(df.detected.mean()),
        "runtime_seconds": {"total_this_invocation": runtime, "setup_reconstruct": t_recon, "realizations": time.time() - t_run,
                            "mean_per_realization_worker": float(df.seconds.mean())},
        "workers": a.workers, "environment": env,
        "command": "PYEPI_SCIENTIFIC_PACKAGE=/home/user/pyepi-upstream/pyEpiAneufinder /home/user/envs/cnv310/bin/python -I "
                   + " ".join(sys.argv),
    }
    C.write_json(out / "summary.json", summary)

    write_report(out, df, summ, summary)
    print("done", f"{runtime:.0f}s", "realizations", len(df), flush=True)


if __name__ == "__main__":
    main()
