"""Experiment 05: calibration and stability of the AD scan used for cluster pseudobulks.

A. Selection-aware null: on exchangeable null sequences (bins of the all-cell chr7
   pseudobulk randomly permuted, so no spatial structure remains), run the workbench's
   first-node procedure (argmax AD over geometry-valid cuts, then pinned
   permutation_test_ad at that split and pinned global_permutation_test_ad on a permuted
   genome pool) and count acceptances at alpha=0.001. A calibrated test accepts ~0.1%.
   A max-over-cuts permutation p-value is computed for a subset as the selection-aware
   comparator.
B. Sampling share: split cells into two halves within each WNN cluster (seed tuple
   (20261104, 2)); the half-difference D=(A-B)/2 contains cell-sampling variation only.
   Compare its bin-to-bin variability with the full pseudobulk's; scan D with the
   workbench settings.
C. Boundary stability: scan chr7/chr9 with server.allelic_scan.scan_chromosome (workbench
   DEFAULTS, genome pool) on the full and half pseudobulks of `all` and `cluster0`.
   The full scans are a harness check against the saved genome_pool breakpoints.

Pinned statistics come from pyEpi e8305ea via server/pyepi_adapter. Nothing is written
outside this folder.
"""
import json
import os
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

ALPHA = 0.001
N_PERM = 1000
SEED = C.SEED
CHROMS = ["chr7", "chr9"]


def null_trial(args):
    """One selection-then-test trial on a permuted sequence."""
    i, seq, pool, L, pkg, do_max = args
    sci.load(pkg)
    rng = np.random.default_rng([SEED, 1, L, i])
    s = rng.permutation(seq)[:L]
    g = rng.permutation(pool)
    dv = np.array(sci.seq_dist_ad(s))
    order = np.argsort(-dv, kind="stable") + 1
    cut = None
    for c in order:
        ok, _ = sc.validate_candidate_geometry(int(c), 0, L, 0, L, 20, 50)
        if ok:
            cut = int(c); break
    stat = dv[cut - 1]
    p_loc = sci.permutation_test_ad(s[:cut], s[cut:], stat, N_PERM)
    p_glob = sci.global_permutation_test_ad(g, cut, L - cut, stat, N_PERM)
    p_max = None
    if do_max:
        # selection-aware: permute the sequence, recompute the max over valid cuts (200 perms)
        valid = np.array([sc.validate_candidate_geometry(c, 0, L, 0, L, 20, 50)[0] for c in range(1, L)])
        cnt = 0
        for b in range(200):
            ps = rng.permutation(s)
            pv = np.array(sci.seq_dist_ad(ps))
            if pv[valid].max() >= stat:
                cnt += 1
        p_max = (cnt + 1) / 201
    return {"trial": i, "L": L, "cut": cut, "stat": float(stat), "p_local": p_loc, "p_global": p_glob,
            "accepted": bool(p_loc < ALPHA and p_glob < ALPHA), "p_max_over_cuts_200perm": p_max}


def scan(args):
    name, chrom, values, pool, pkg = args
    sci.load(pkg)
    out = sc.scan_chromosome(values, np.ones(len(values), bool), sc.DEFAULTS, package_dir=pkg, pool=pool)
    return name, chrom, [b["bin"] for b in out["breakpoints"]], out["segments"]


def main():
    t0 = time.time()
    p = C.default_paths()
    pkg = str(p["package"])
    d = C.H5AD(p["h5ad"])
    lab = C.cell_clusters(d, p)
    half = C.stratified_halves(lab, (SEED, 2))
    masks = {"all": np.ones(d.n_cells, bool), "all_A": half, "all_B": ~half,
             "cluster0": lab == "0", "cluster0_A": (lab == "0") & half, "cluster0_B": (lab == "0") & ~half}
    pbs = C.group_means(d, "X", masks)
    np.savez_compressed(HERE / "pseudobulks.npz", **pbs, seq=d.seq, start=d.start)
    workers = os.cpu_count() or 2

    # ---- A. null calibration
    chr7 = pbs["all"][d.seq == "chr7"]
    tasks = [(i, chr7, pbs["all"], L, pkg, i < 20) for L in (120, 400) for i in range(150)]
    with ProcessPoolExecutor(workers) as ex:
        nulls = pd.DataFrame(list(ex.map(null_trial, tasks)))
    nulls.to_csv(HERE / "null_calibration.tsv", sep="\t", index=False)
    calib = nulls.groupby("L").agg(trials=("accepted", "size"), accepted=("accepted", "sum"),
                                   p_local_lt_alpha=("p_local", lambda x: int((x < ALPHA).sum())),
                                   median_p_local=("p_local", "median")).reset_index()
    sub = nulls.dropna(subset=["p_max_over_cuts_200perm"])
    calib_max = sub.groupby("L").agg(n=("p_max_over_cuts_200perm", "size"),
                                     median_p_max=("p_max_over_cuts_200perm", "median"),
                                     frac_p_max_lt_0_05=("p_max_over_cuts_200perm", lambda x: float((x < .05).mean())),
                                     frac_fixed_split_p_lt_alpha=("p_local", lambda x: float((x < ALPHA).mean()))).reset_index()

    # ---- B. sampling share
    rows = []
    for g in ["all", "cluster0"]:
        A, B = pbs[g + "_A"], pbs[g + "_B"]
        D = (A - B) / 2
        M = pbs[g]
        for ch in CHROMS + ["genome"]:
            m = np.ones(d.n_bins, bool) if ch == "genome" else d.seq == ch
            mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))  # noqa: E731
            rows.append({"source": g, "region": ch, "cells": int(masks[g].sum()),
                         "robust_sd_full": mad(M[m]), "robust_sd_half_diff": mad(D[m]),
                         "sampling_share_of_variance": (mad(D[m]) / mad(M[m])) ** 2,
                         "lag1_autocorr_full": float(np.corrcoef(M[m][:-1], M[m][1:])[0, 1]),
                         "lag1_autocorr_half_diff": float(np.corrcoef(D[m][:-1], D[m][1:])[0, 1]),
                         "split_half_corr": float(np.corrcoef(A[m], B[m])[0, 1])})
    share = pd.DataFrame(rows)
    share.to_csv(HERE / "sampling_share.tsv", sep="\t", index=False)

    # ---- C. scans (full, halves, half-difference)
    stasks = []
    for g in ["all", "cluster0"]:
        D = (pbs[g + "_A"] - pbs[g + "_B"]) / 2
        Dn = D + np.mean(pbs[g])   # shift so values are positive like a depth profile; AD is shift-invariant
        for name, vec in [(g, pbs[g]), (g + "_A", pbs[g + "_A"]), (g + "_B", pbs[g + "_B"]), (g + "_halfdiff", Dn)]:
            for ch in CHROMS:
                stasks.append((name, ch, vec[d.seq == ch], vec, pkg))
    with ProcessPoolExecutor(workers) as ex:
        scans = list(ex.map(scan, stasks))
    bps = {(n, ch): b for n, ch, b, _ in scans}
    saved = pd.read_csv(p["genome_pool"] / "raw_accepted_breakpoints.tsv", sep="\t")
    stab = []
    harness = {}
    for g in ["all", "cluster0"]:
        for ch in CHROMS:
            ref = sorted(saved.loc[(saved.source == g) & (saved.chromosome == ch), "absolute_bin"].astype(int))
            full = sorted(bps[(g, ch)])
            harness[f"{g}/{ch}"] = {"saved": len(ref), "rerun": len(full), "identical": ref == full}
            A, B = np.array(bps[(g + "_A", ch)]), np.array(bps[(g + "_B", ch)])
            for b in full:
                dA = int(np.min(np.abs(A - b))) if len(A) else None
                dB = int(np.min(np.abs(B - b))) if len(B) else None
                stab.append({"source": g, "chrom": ch, "bin_in_chrom": b,
                             "genomic_start": int(d.start[np.flatnonzero(d.seq == ch)[b]]),
                             "nearest_A": dA, "nearest_B": dB,
                             "in_both_halves_pm2": bool(dA is not None and dB is not None and dA <= 2 and dB <= 2),
                             "in_both_halves_pm5": bool(dA is not None and dB is not None and dA <= 5 and dB <= 5)})
    stab = pd.DataFrame(stab)
    stab.to_csv(HERE / "boundary_stability.tsv", sep="\t", index=False)
    halfdiff = {f"{g}/{ch}": len(bps[(g + '_halfdiff', ch)]) for g in ["all", "cluster0"] for ch in CHROMS}
    segs = {(n, ch): sg for n, ch, _, sg in scans}
    hd_rows = []
    for g in ["all", "cluster0"]:
        D = (pbs[g + "_A"] - pbs[g + "_B"]) / 2
        for ch in CHROMS:
            m = d.seq == ch
            for sg in segs[(g + "_halfdiff", ch)]:
                sl = slice(sg["bin_start"], sg["bin_end"] + 1)
                hd_rows.append({"source": g, "chrom": ch, "start_mb": d.start[m][sg["bin_start"]] / 1e6,
                                "end_mb": d.start[m][sg["bin_end"]] / 1e6 + 0.1, "n_bins": sg["n_bins"],
                                "sd_half_diff": float(D[m][sl].std()), "mean_full_pseudobulk": float(pbs[g][m][sl].mean())})
    hd = pd.DataFrame(hd_rows)
    hd.to_csv(HERE / "half_difference_segments.tsv", sep="\t", index=False)
    per_half = {f"{n}/{ch}": len(b) for (n, ch), b in bps.items()}
    stab_sum = stab.groupby(["source", "chrom"]).agg(n=("in_both_halves_pm2", "size"), both_pm2=("in_both_halves_pm2", "sum"),
                                                      both_pm5=("in_both_halves_pm5", "sum")).reset_index()

    # candidate-window boundaries (shared regions table, 1-based start coordinates)
    reg = C.shared_regions()
    tgt = reg[reg.part == "target"]
    cand = []
    for _, r in tgt.iterrows():
        sub_s = stab[(stab.chrom == r.chrom)]
        for edge, pos in [("left", r.first_retained_start), ("right", r.last_retained_end + 1)]:
            hit = sub_s.iloc[(sub_s.genomic_start - pos).abs().argsort()[:1]] if len(sub_s) else sub_s
            for _, h in hit.iterrows():
                if abs(h.genomic_start - pos) <= 200_000:
                    cand.append({"region": r.region_id.replace("__target", ""), "edge": edge, "source": h.source,
                                 "bp_start": h.genomic_start, "nearest_A": h.nearest_A, "nearest_B": h.nearest_B,
                                 "in_both_halves_pm2": h.in_both_halves_pm2})
    cand = pd.DataFrame(cand)
    cand.to_csv(HERE / "candidate_edge_stability.tsv", sep="\t", index=False)

    res = {"experiment": "05_segmentation_calibration", "seed": SEED, "generator": "numpy.random.default_rng (PCG64)",
           "seed_streams": {"halves": [SEED, 2], "null_trials": [SEED, 1, "L", "trial"]},
           "settings": sc.DEFAULTS, "alpha": ALPHA, "n_permutations": N_PERM,
           "null_calibration": calib.to_dict("records"), "null_max_over_cuts": calib_max.to_dict("records"),
           "harness_check_vs_saved_genome_pool": harness, "half_difference_breakpoints": halfdiff,
           "breakpoints_per_scan": per_half, "boundary_stability": stab_sum.to_dict("records"),
           "seconds": round(time.time() - t0, 1), "environment": C.environment()}
    C.write_json(HERE / "summary.json", res)

    # ---- figure + report
    fig, ax = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True)
    m = d.seq == "chr7"; x = d.start[m] / 1e6
    ax[0].plot(x, pbs["all_A"][m], lw=.6, label="half A"); ax[0].plot(x, pbs["all_B"][m], lw=.6, label="half B")
    for b in bps[("all", "chr7")]:
        ax[0].axvline(x[b], color="#9b2c2c", lw=.5, alpha=.6)
    ax[0].set_ylabel("mean X (all cells)"); ax[0].legend(fontsize=8); ax[0].set_title("chr7 all-cell half pseudobulks; red = rerun breakpoints")
    ax[1].plot(x, (pbs["all_A"][m] - pbs["all_B"][m]) / 2, lw=.6, color="#555")
    ax[1].set_ylabel("(A − B) / 2"); ax[1].set_xlabel("chr7 position (Mb)")
    ax[1].set_ylim(ax[0].get_ylim()[0] - np.mean(pbs["all"][m]), ax[0].get_ylim()[1] - np.mean(pbs["all"][m]))
    b64 = C.fig_to_b64(fig); plt.close(fig)
    acc = calib.set_index("L")
    from scipy.stats import beta as _beta
    n_tr, n_acc = int(calib.trials.sum()), int(calib.accepted.sum())
    cp_hi = float(_beta.ppf(0.975, n_acc + 1, n_tr - n_acc)) if n_tr > n_acc else 1.0
    sh = share.set_index(["source", "region"])
    harness_txt = '; '.join(k + ' ' + str(v['rerun']) for k, v in harness.items())
    harness_ok = all(v["identical"] for v in harness.values())
    n_edges_hit = len(cand)
    n_edges_both = int(cand.in_both_halves_pm2.sum()) if len(cand) else 0
    n_total_bp = int(stab_sum.n.sum()); n_total_pm2 = int(stab_sum.both_pm2.sum())
    ss = stab_sum.set_index(["source", "chrom"])
    # per-window edge notes
    win_rows = []
    for rid in tgt.region_id.str.replace("__target", "", regex=False):
        c = cand[cand.region == rid] if len(cand) else cand
        if len(c) == 0:
            note = "No breakpoint of the all or cluster0 scan lies within 200 kb of either window edge; the window&rsquo;s edges are not supported by this scan."
        else:
            edges_txt = "; ".join(f"{r.edge} edge: {r.source} breakpoint at {r.bp_start/1e6:.1f} Mb, nearest half-A/half-B breakpoint {r.nearest_A}/{r.nearest_B} bins away" for r in c.itertuples())
            both = bool(c.in_both_halves_pm2.all())
            note = f"{len(c)} of 2 edges coincide with a breakpoint ({edges_txt}); " + ("reproduced in both disjoint halves within &plusmn;2 bins." if both else "not all reproduced in both halves within &plusmn;2 bins.")
        win_rows.append({"window": rid, "edges_with_breakpoint": len(c), "reproduced_both_halves_pm2": int(c.in_both_halves_pm2.sum()) if len(c) else 0, "note": __import__("html").unescape(note)})
    win_df = pd.DataFrame(win_rows)
    win_html = C.html_table(win_df).replace("<td>No breakpoint", "<td style='text-align:left;white-space:normal'>No breakpoint").replace("<td>1 of 2", "<td style='text-align:left;white-space:normal'>1 of 2").replace("<td>2 of 2", "<td style='text-align:left;white-space:normal'>2 of 2")
    verdict = (f"Boundaries are reproducible between disjoint halves of the cells (split-half r = {sh.split_half_corr.min():.3f}&ndash;{sh.split_half_corr.max():.3f}, cell sampling &le; {100*sh.sampling_share_of_variance.max():.2f}% of pseudobulk variance) "
               f"and the select-then-test scan accepted {n_acc} of {n_tr} exchangeable nulls at &alpha; = {ALPHA}, so a significant boundary shows reproducible non-sampling structure, not that the structure is copy number.")
    findings = [
        f"Exchangeable null: {int(acc.loc[120,'accepted'])}/{int(acc.loc[120,'trials'])} (120 bins) and {int(acc.loc[400,'accepted'])}/{int(acc.loc[400,'trials'])} (400 bins) first splits accepted at &alpha; = {ALPHA}; pooled {n_acc}/{n_tr} (upper 95% Clopper&ndash;Pearson bound {100*cp_hi:.1f}%, i.e. not excluding rates above &alpha;).",
        f"Fixed-split p-values ignore the search over cuts: median p_local {acc.loc[120,'median_p_local']:.2f} (120 bins) and {acc.loc[400,'median_p_local']:.3f} (400 bins) versus median max-over-cuts p {calib_max.set_index('L').loc[120,'median_p_max']:.2f} and {calib_max.set_index('L').loc[400,'median_p_max']:.2f} (20 trials each, 200 permutations).",
        f"Cell sampling explains {100*sh.loc[('all','genome'),'sampling_share_of_variance']:.2f}% of the <code>all</code> pseudobulk bin variance (robust SD {sh.loc[('all','genome'),'robust_sd_half_diff']:.3f} of half-difference vs {sh.loc[('all','genome'),'robust_sd_full']:.3f} full) and {100*sh.loc[('cluster0','genome'),'sampling_share_of_variance']:.2f}% for cluster 0; split-half correlation {sh.loc[('all','genome'),'split_half_corr']:.4f} / {sh.loc[('cluster0','genome'),'split_half_corr']:.4f}.",
        f"Harness: the rerun reproduces the saved genome-pool breakpoints exactly for {sum(v['identical'] for v in harness.values())} of {len(harness)} source/chromosome scans ({harness_txt} breakpoints).",
        f"Half-split stability (&plusmn;2 bins in both halves): all/chr7 {int(ss.loc[('all','chr7'),'both_pm2'])}/{int(ss.loc[('all','chr7'),'n'])}, all/chr9 {int(ss.loc[('all','chr9'),'both_pm2'])}/{int(ss.loc[('all','chr9'),'n'])}, cluster0/chr7 {int(ss.loc[('cluster0','chr7'),'both_pm2'])}/{int(ss.loc[('cluster0','chr7'),'n'])}, cluster0/chr9 {int(ss.loc[('cluster0','chr9'),'both_pm2'])}/{int(ss.loc[('cluster0','chr9'),'n'])} ({n_total_pm2} of {n_total_bp} overall).",
        f"{n_edges_hit} of {2*len(tgt)} saved-window edges lie within 200 kb of a rerun breakpoint, and {n_edges_both} of those {n_edges_hit} reappear in both halves within &plusmn;2 bins; windows with no edge breakpoint: {', '.join(win_df[win_df.edges_with_breakpoint==0].window) or 'none'}.",
        f"A sampling-only profile (half-difference) still yields accepted breakpoints in <code>all</code> ({halfdiff['all/chr7']} on chr7 at {', '.join(f'{v:.1f}' for v in hd[(hd.source=='all')&(hd.chrom=='chr7')].start_mb.iloc[1:])} Mb, {halfdiff['all/chr9']} on chr9) and none in cluster 0, with segment SD tracking depth (e.g. chr7 51.1&ndash;77.4 Mb SD {hd[(hd.source=='all')&(hd.chrom=='chr7')].sd_half_diff.iloc[1]:.3f} at pseudobulk mean {hd[(hd.source=='all')&(hd.chrom=='chr7')].mean_full_pseudobulk.iloc[1]:.2f}): AD responds to variance changes as well as means.",
    ]
    sections = [
        ("Question and why it matters",
         "<p><b>Question.</b> Do the pinned local/global permutation tests, as used after choosing the maximum-AD cut, control false acceptance? How much of the cluster pseudobulk&rsquo;s bin-to-bin variation is cell-sampling noise, and are boundaries reproducible between disjoint halves of the same cells?</p>"
         "<p><b>Why it matters.</b> The recommendation states that the saved boundaries are reproducible and that the problem is interpreting them. I started from the hypothesis that the select-then-test procedure is anti-conservative; "
         "this experiment tests that hypothesis and quantifies how much boundary structure is cell sampling versus reproducible non-sampling structure, which determines whether a significant boundary can be used as evidence of copy-number change in the per-region evidence layer.</p>"),
        ("Inputs",
         f"<ul><li><b>Matrix.</b> processed <code>X</code> of <code>count_matrix.h5ad</code> ({d.n_cells:,} cells &times; {d.n_bins:,} bins), averaged per group: <code>all</code> ({int(masks['all'].sum()):,} cells) and WNN <code>cluster0</code> ({int(masks['cluster0'].sum()):,} cells), each also split into halves A and B (<code>pseudobulks.npz</code>).</li>"
         "<li><b>Chromosomes scanned.</b> chr7 and chr9 (the chromosomes of the saved windows). The null uses the all-cell chr7 profile.</li>"
         "<li><b>Saved breakpoints.</b> <code>genome_pool/raw_accepted_breakpoints.tsv</code> (harness check) and the workbench pool of the all-cell pseudobulk for the global permutation.</li>"
         f"<li><b>Windows.</b> The {len(tgt)} saved target windows (shared <code>regions.tsv</code>, 1-based start coordinates).</li>"
         f"<li><b>Settings (workbench DEFAULTS).</b> k = {sc.DEFAULTS['k']}, min_child_bins = {sc.DEFAULTS['min_child_bins']}, min_terminal_bins = {sc.DEFAULTS['min_terminal_bins']}, {N_PERM:,} permutations, &alpha; = {ALPHA}.</li>"
         "<li><b>Not available.</b> Per-cell-level segmentation of the production run (<code>pyepi_results/</code> used k = 2 per cell, a different procedure), and any independent CN truth.</li></ul>"),
        ("Method",
         "<ol><li><b>A. Selection-aware null.</b> Randomly permute the bins of the all-cell chr7 pseudobulk (so no spatial structure remains), take the first L bins (L = 120 and 400; 150 trials each), pick the geometry-valid cut with maximal AD distance, then apply the pinned <code>permutation_test_ad</code> at that split (local) and <code>global_permutation_test_ad</code> against a permuted genome pool (global). Accept if both p &lt; &alpha;. "
         "For the first 20 trials per L a max-over-cuts permutation p-value (200 permutations) is computed as the selection-aware comparator.</li>"
         "<li><b>B. Sampling share.</b> Split cells into two halves within each cluster; the half-difference D = (A &minus; B)/2 contains cell-sampling variation only. Compare robust SD of D with that of the full pseudobulk.</li>"
         "<li><b>C. Boundary stability.</b> Scan chr7 and chr9 with <code>server.allelic_scan.scan_chromosome</code> (workbench DEFAULTS, genome pool) on the full profile, each half and the half-difference, for <code>all</code> and <code>cluster0</code>. A full-scan breakpoint is &lsquo;stable&rsquo; if each half has a breakpoint within &plusmn;2 (or &plusmn;5) bins.</li></ol>"
         "<pre>D = (A - B) / 2                       # pure cell-sampling variation\nsampling_share = (MAD_sd(D) / MAD_sd(full))^2,  MAD_sd(v) = 1.4826 * median(|v - median(v)|)\nsplit_half_corr = corr(A, B)\ncut = argmax over geometry-valid cuts of AD distance;  accept iff p_local &lt; 0.001 and p_global &lt; 0.001\np_max = (1 + #{perms with max valid-cut AD &gt;= observed}) / 201</pre>"
         f"<p>Seeds: halves <code>default_rng(({SEED}, 2))</code>; null trials <code>default_rng([{SEED}, 1, L, trial])</code> (numpy PCG64). Pinned statistics come from pyEpi e8305ea via <code>server/pyepi_adapter</code>. Nothing is written outside this folder.</p>"),
        ("Results", "", 2),
        ("A. Selection-aware null", C.html_table(calib) + C.caption("Per sequence length L: <code>trials</code> run, <code>accepted</code> (both local and global p &lt; &alpha;), trials with local p &lt; &alpha;, and the median local p-value. A calibrated test accepts about 0.1% of trials.") +
         C.html_table(calib_max) + C.caption("Selection-aware comparator on the first 20 trials per L: median max-over-cuts p-value, the fraction below 0.05, and the fraction of those same trials with fixed-split p &lt; &alpha;. Look for max-over-cuts p-values far larger than fixed-split p-values.") +
         f"<div class='card'>On sequences with no spatial structure, the procedure accepted <b>{int(acc.loc[120,'accepted'])}/{int(acc.loc[120,'trials'])}</b> (120 bins) and <b>{int(acc.loc[400,'accepted'])}/{int(acc.loc[400,'trials'])}</b> (400 bins) first splits at &alpha; = {ALPHA}. "
         f"The fixed-split permutation p-value ignores the search over cuts (its median is far below the max-over-cuts p-value), but at &alpha; = {ALPHA} the realised false-acceptance rate on exchangeable data stayed low ({n_acc}/{n_tr}; upper 95% Clopper&ndash;Pearson bound {100*cp_hi:.1f}%). "
         "<b>The prespecified hypothesis that the scan is grossly anti-conservative under an exchangeable null was not supported.</b> This does not address spatially dependent nulls, which part B shows are the relevant case.</div>", 3),
        ("B. Sampling share of variation", C.html_table(share) +
         C.caption("Per source and region: robust SD of the full pseudobulk and of the half-difference, their squared ratio (<code>sampling_share_of_variance</code>), lag-1 autocorrelations (spatial structure) and the split-half correlation. Expect a tiny sampling share, high split-half correlation and near-zero autocorrelation of the half-difference.") +
         "<p>The half-difference contains only cell-sampling variation (same clusters, disjoint cells). A small sampling share and high split-half correlation mean that almost all pseudobulk structure is reproducible between cells; reproducible structure can be dosage, accessibility or a technical baseline, so passing these tests does not by itself indicate CN.</p>", 3),
        ("C. Boundary stability", C.html_table(stab_sum) + C.caption("Number of rerun breakpoints (<code>n</code>) and how many have a breakpoint within &plusmn;2 (<code>both_pm2</code>) or &plusmn;5 (<code>both_pm5</code>) bins in <em>both</em> halves. Full list: <code>boundary_stability.tsv</code>.") +
         "<p>Harness check against the saved genome-pool run:</p>" + C.html_table(pd.DataFrame([{"scan": k, **v} for k, v in harness.items()])) +
         C.caption("Saved vs rerun breakpoint counts and whether the positions are identical; True everywhere means the stability results describe the saved boundaries.") +
         "<p>Breakpoints accepted when the half-difference (sampling-only) profile is scanned: " + ", ".join(f"{k}: {v}" for k, v in halfdiff.items()) + ".</p>", 3),
        ("Sampling-only (half-difference) segments", C.html_table(hd) + C.caption("Segments of the half-difference profile in the <code>all</code> and <code>cluster0</code> sources. <code>sd_half_diff</code> is the SD of D inside the segment and <code>mean_full_pseudobulk</code> the mean of the full pseudobulk there. A segment with lower SD at lower pseudobulk mean means variance, not mean, drives the split.") +
         "<p>A profile that contains only cell-sampling variation still yields accepted breakpoints in <code>all</code>. Its segments differ in noise scale (SD), which tracks the full pseudobulk depth: AD compares whole distributions, so a change in variance alone can be detected. Low-depth regions therefore present a distribution change even before any mean shift.</p>" +
         C.img(b64, "chr7 half pseudobulks") + C.caption("Top: chr7 all-cell pseudobulk for half A and half B (they overlay almost exactly); red lines are the rerun breakpoints. Bottom: (A &minus; B)/2 on the same scale, i.e. the sampling-only component; look for its small amplitude and its smaller spread in the low-depth region around 51&ndash;77 Mb."), 3),
        ("Per-region detail",
         win_html + C.caption("One row per saved window. An edge counts if a rerun breakpoint of the all or cluster0 scan lies within 200 kb of the window&rsquo;s first retained bin or the next bin after its last. Full list: <code>candidate_edge_stability.tsv</code>.") +
         (C.html_table(cand) + C.caption("Edge-level detail: breakpoint source and position, distance in bins to the nearest half-A and half-B breakpoint, and whether both are within &plusmn;2 bins.") if len(cand) else "")),
        ("Sensitivity",
         "<p>Two sequence lengths (120 and 400 bins) for the null; two sources (<code>all</code>, <code>cluster0</code>) and two stability tolerances (&plusmn;2 and &plusmn;5 bins) for stability. Not varied: &alpha;, k, min_child_bins, number of permutations, the single halving seed, and the number of null trials (150 per L; the null acceptance rate is imprecise at this size).</p>"),
        ("Limits and alternative explanations",
         "<ul><li>The exchangeable null removes spatial structure; it shows the permutation tests are not grossly loose on exchangeable data, not that significance of a boundary in structured data means CN. Reproducible accessibility, mappability and GC residue structure also rejects an exchangeable null.</li>"
         "<li>Stability is assessed on halves of the same cells and the same sample: it shows boundaries are not cell-sampling artefacts but not that they are shared across samples or are CN.</li>"
         "<li>150 trials per length cannot distinguish an acceptance rate of 0.1% from several times that (see Clopper&ndash;Pearson bound).</li>"
         "<li>Only chr7 and chr9 were scanned; all half-split results are for those chromosomes.</li>"
         "<li>The production per-cell run used a different segmentation (k = 2 per cell), so these results describe the workbench cluster scan only.</li></ul>"),
    ]
    glossary = {
        "L": "Length in bins of the permuted null sequence.",
        "trials / accepted": "Null trials run / trials in which both local and global permutation p-values were below &alpha;.",
        "p_local_lt_alpha / median_p_local": "Trials with local p &lt; &alpha; / median of the local permutation p-value at the chosen split.",
        "median_p_max / frac_p_max_lt_0_05": "Median of the max-over-cuts p-value (200 permutations) and the fraction below 0.05.",
        "frac_fixed_split_p_lt_alpha": "Fraction of the same trials whose fixed-split local p is below &alpha;.",
        "robust_sd_full / robust_sd_half_diff": "1.4826 &times; median absolute deviation of the full pseudobulk / of the half-difference over bins.",
        "sampling_share_of_variance": "(robust SD of half-difference / robust SD of full)^2: share of pseudobulk bin variance due to cell sampling.",
        "lag1_autocorr_full / lag1_autocorr_half_diff": "Lag-1 correlation between adjacent bins of the full profile / half-difference.",
        "split_half_corr": "Correlation across bins between half A and half B pseudobulks.",
        "n / both_pm2 / both_pm5": "Rerun breakpoints / those with a breakpoint within &plusmn;2 or &plusmn;5 bins in both halves.",
        "bin_in_chrom / genomic_start": "Breakpoint position as bin index in the chromosome / genomic start (bp).",
        "nearest_A / nearest_B": "Distance in bins to the nearest breakpoint of half A / half B scan.",
        "in_both_halves_pm2": "True if both nearest distances are &le; 2 bins.",
        "sd_half_diff / mean_full_pseudobulk": "SD of the half-difference and mean of the full pseudobulk within a segment.",
        "edges_with_breakpoint / reproduced_both_halves_pm2": "Window edges with a rerun breakpoint within 200 kb / of those, reproduced in both halves within &plusmn;2 bins.",
        "saved / rerun / identical": "Saved breakpoint count / rerun count / whether positions are identical.",
        "AD": "Anderson&ndash;Darling distance between the bin-value distributions on either side of a cut.",
        "p_global": "Pinned genome-pool permutation p-value at the chosen split.",
    }
    C.write_html_report(HERE / "report.html", "Segmentation calibration", "Experiment 05 &middot; AD scan null behaviour, sampling share and half-split stability",
                        sections, {"command": "python -I experiments/05_segmentation_calibration/run.py", "summary": "summary.json",
                                   "seed": SEED, "streams": res["seed_streams"], "environment": C.environment()},
                        key_findings=findings, verdict=verdict, glossary=glossary,
                        supports=["Boundaries are not cell-sampling artefacts: halves of the same cells agree almost perfectly and most breakpoints reappear in both halves.",
                                  "On exchangeable nulls the select-then-test AD scan rarely accepts at the stated &alpha;, so the p-values are not grossly anti-conservative in that setting.",
                                  "Window edges that coincide with a breakpoint are reproduced in both halves."],
                        not_supports=["That a significant boundary indicates a copy-number change: reproducible accessibility, mappability and GC structure also reject an exchangeable null.",
                                      "Calibration under spatially dependent nulls or in other samples.",
                                      "Edge support for windows with no nearby breakpoint."],
                        decision_impact="Supports the recommendation that significance of a segmentation boundary should not be used as CN evidence on its own, and that the per-region evidence layer must weigh depth contrast, coverage invariance and allelic evidence instead. "
                                        "It also shows the scan is variance-sensitive in low-depth regions, which argues for evaluating an exposure-aware raw-count segmentation (rank 4) only after a cross-sample baseline exists.",
                        next_steps=["Calibrate against nulls with realistic spatial dependence, e.g. block-permuted or baseline-structured profiles from unrelated samples.",
                                    "Increase null trials to resolve the acceptance rate relative to &alpha;.",
                                    "Repeat the half-split analysis on an independent sample."])
    print(json.dumps({k: res[k] for k in ["null_calibration", "null_max_over_cuts", "harness_check_vs_saved_genome_pool",
                                          "half_difference_breakpoints", "boundary_stability", "seconds"]}, indent=0, default=str))


if __name__ == "__main__":
    main()
