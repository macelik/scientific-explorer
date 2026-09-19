"""run_depth_segment_prototype.py - per-cluster depth-segment merge prototype.

Seven-chromosome pilot.  Reuses the existing mean-over-cells pseudobulks and
accepted breakpoints (no membership change, no AD rerun).  Implements the
effect-based merge rule from docs/plans/2026-09-07-depth-segment-prototype-plan.md:

  - observed median per-bin segment levels on the source mean pseudobulk;
  - whole-cell bootstrap stability (>= stab fraction of replicates with
    |log2(mA/mB)| < tau, all replicates valid required);
  - smallest-observed-|delta| deterministic sequential merge;
  - chain guard: original constituent log2 median range < tau inside any
    proposed merged interval;
  - terminal segments, long-long boundaries, and gap/zero-level cases retained
    and reported, never silently merged.

Writes a full boundary ledger and simplified partitions; preserves all
reference files unchanged.  No CN states are assigned in this pilot.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "analysis/depth_segment_prototype"
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))

import loader  # noqa: E402
import merge_core as mc  # noqa: E402

TAUS = [0.15, 0.20, 0.25]
STABS = [0.80, 0.90, 0.95]
DEFAULT_CHROMOSOMES = ["chr1", "chr2", "chr5", "chr7", "chr8", "chr9", "chr13"]


def log(msg):
    print(f"[depth-prototype] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser(
        description="per-cluster depth-segment merge prototype (7-chr pilot)")
    p.add_argument("--adata", required=True)
    p.add_argument("--integrated", required=True)
    p.add_argument("--cluster-key", default="wnn_leiden_0.3")
    p.add_argument("--breakpoints", required=True)
    p.add_argument("--reference-scores", default=None,
                   help="outputs/screening/small_segment_scores.tsv (coordinate check)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--figures-dir", default=None)
    p.add_argument("--chromosomes", nargs="+", default=DEFAULT_CHROMOSOMES)
    p.add_argument("--max-small-bins", type=int, default=100)
    p.add_argument("--tau", type=float, default=0.20)
    p.add_argument("--stab", type=float, default=0.90)
    p.add_argument("--n-bootstrap", type=int, default=200)
    p.add_argument("--seed", type=int, default=20260907)
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def fingerprint(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def pooled_phi(profile, seg_ids):
    profile = np.asarray(profile, dtype=float)
    seg_ids = np.asarray(seg_ids)
    num = 0.0
    df = 0
    for sid in np.unique(seg_ids):
        p = profile[seg_ids == sid]
        E = len(p)
        if E < 2:
            continue
        m = p.mean()
        if m <= 0:
            continue
        num += float(((p - m) ** 2 / m).sum())
        df += int(E - 1)
    if df == 0:
        return np.nan
    return float(num / df)


def iqr_trimmed_mean(vals):
    """Mean after removing values outside [Q1-1.5IQR, Q3+1.5IQR].

    Robust central level, used as an alternative to the plain mean/median for
    the small-segment flank effect.
    """
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return np.nan
    q1, q3 = np.percentile(v, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    t = v[(v >= lo) & (v <= hi)]
    return float(t.mean()) if len(t) else float(v.mean())


def reproduce_check(adata, members, raw, chrom_of, reference_path, chromosomes,
                    max_small):
    """Reproduce cluster-only small-segment coordinates/medians vs stored
    screening reference (rounding tolerance).  Returns a summary dict."""
    if reference_path is None or not Path(reference_path).exists():
        return {"status": "skipped_no_reference"}
    ref = pd.read_csv(reference_path, sep="\t")
    ref = ref[ref["source"].str.startswith("cluster")]
    X = adata.X
    n_match = 0
    n_rows = 0
    max_median_diff = 0.0
    for ch in chromosomes:
        cols = np.where(chrom_of == ch)[0]
        n_bins = len(cols)
        ref_ch = ref[ref["chromosome"] == ch]
        for src in [f"cluster{i}" for i in range(5)]:
            rows = ref_ch[ref_ch["source"] == src]
            if rows.empty:
                continue
            acc = sorted(int(b) for b in raw[
                (raw["source"] == src) & (raw["chromosome"] == ch)]["absolute_bin"])
            segs, small = loader.small_internal_segments(acc, n_bins, max_small)
            seg_map = {(s, e): k for k, s, e in small}
            pb = np.asarray(X[members[src]][:, cols].mean(axis=0)).ravel()
            exp_rows = {(int(r.bin_start), int(r.bin_end)) for _, r in rows.iterrows()}
            got_rows = set(seg_map.keys())
            if exp_rows != got_rows:
                return {"status": "FAIL_segment_mismatch", "source": src,
                        "chromosome": ch,
                        "extra": sorted(exp_rows - got_rows)[:5],
                        "missing": sorted(got_rows - exp_rows)[:5]}
            for _, r in rows.iterrows():
                n_rows += 1
                key = (int(r.bin_start), int(r.bin_end))
                if key not in seg_map:
                    continue
                s, e = key
                med = float(np.median(pb[s:e]))
                max_median_diff = max(max_median_diff, abs(med - float(r["median_level"])))
                if abs(med - float(r["median_level"])) <= 5e-5:
                    n_match += 1
    status = "OK" if (n_rows and n_match == n_rows) else "FAIL_value_diff"
    return {"status": status, "n_rows": n_rows, "n_match": n_match,
            "max_median_diff": float(max_median_diff)}


def process_source_chr(args, X, members, raw, src, ch, chrom_of, positions_all,
                       cluster_reps):
    """Run the full merge prototype for one (source, chromosome)."""
    cols = np.where(chrom_of == ch)[0]
    n_bins = len(cols)
    positions = positions_all[cols]
    acc = loader.boundary_bins_for(raw, src, ch)
    orig_segs = mc.flat_partition(acc, n_bins)
    small = [(k, s.start, s.end) for k, s in enumerate(orig_segs)
             if k not in (0, len(orig_segs) - 1)
             and s.end - s.start <= args.max_small_bins]
    pb_obs = loader.observed_pseudobulk(X, members[src], cols)
    need_boot = len(small) > 0
    pb_boot = None
    if need_boot:
        pb_boot = loader.replicate_pseudobulks(
            X, members[src], cols, cluster_reps[src]["reps"], args.n_bootstrap)

    n_orig_bounds = max(0, len(orig_segs) - 1)
    if not need_boot and n_orig_bounds == 0:
        return {"skipped": True, "source": src, "chromosome": ch}

    # initial boundary snapshot (default tau/stab)
    osegs, snap = mc.initial_boundary_scores(
        pb_obs, pb_boot, acc, n_bins, args.tau, args.stab,
        max_small=args.max_small_bins, need_all_valid=True)

    # final merge under default parameters (pb_boot None allowed -> no merges)
    res = mc.run_merges(pb_obs, pb_boot, acc, n_bins, args.tau, args.stab,
                        max_small=args.max_small_bins, need_all_valid=True)
    final_segs = res["final_segs"]
    olog = res["orig_log2_medians"]

    # per-source outputs
    init_rows = [{"source": src, "chromosome": ch, **r} for r in snap]
    merge_rows = [{"source": src, "chromosome": ch, **m} for m in res["merges"]]

    retained_rows = []
    for j in range(len(orig_segs) - 1):
        a, b = orig_segs[j], orig_segs[j + 1]
        d = res["decided"].get(j, "retained_not_eligible")
        removed = d.startswith("removed_iter")
        retained_rows.append({
            "source": src, "chromosome": ch,
            "orig_boundary_id": j,
            "boundary_bin": int(b.start),
            "Mb_boundary": round(float(positions[min(int(b.start), len(positions)-1)]) / 1e6, 4),
            "left_orig_seg": int(j),
            "right_orig_seg": int(j + 1),
            "left_bins": int(a.end - a.start),
            "right_bins": int(b.end - b.start),
            "removed": bool(removed),
            "retained": bool(not removed),
            "final_reason": d,
        })

    simplified_rows = []
    for fi, s in enumerate(final_segs):
        med = float(np.median(pb_obs[s.start:s.end])) if s.end > s.start else np.nan
        gr = mc.constituent_log2_range(s.orig, olog) if s.orig else np.nan
        simplified_rows.append({
            "source": src, "chromosome": ch, "seg_id": fi,
            "n_orig_constituents": len(s.orig),
            "orig_seg_ids": ",".join(str(i) for i in s.orig),
            "bin_start": int(s.start), "bin_end": int(s.end),
            "seg_length": int(s.end - s.start),
            "Mb_start": round(float(positions[int(s.start)]) / 1e6, 4),
            "Mb_end": round((float(positions[min(int(s.end)-1, len(positions)-1)]) + 100000) / 1e6, 4),
            "median_obs": round(med, 5) if np.isfinite(med) else np.nan,
            "merged": bool(len(s.orig) > 1),
            "constituent_log2_range": round(gr, 5) if np.isfinite(gr) else np.nan,
        })

    # unresolved: original small internal segments plus their final status
    seg_map = {}
    for fi, s in enumerate(final_segs):
        for o in s.orig:
            seg_map[int(o)] = fi
    unresolved_rows = []
    for k, s, e in small:
        med = float(np.median(pb_obs[s:e]))
        left_reason = res["decided"].get(k - 1, "chrom_start") if k - 1 >= 0 else "chrom_start"
        right_reason = res["decided"].get(k, "chrom_end") if k < len(orig_segs) - 1 else "chrom_end"
        l_rem = isinstance(left_reason, str) and left_reason.startswith("removed_iter")
        r_rem = isinstance(right_reason, str) and right_reason.startswith("removed_iter")
        in_final = seg_map.get(k)
        standalone = (in_final is not None
                      and final_segs[in_final].end - final_segs[in_final].start == e - s)
        unresolved_rows.append({
            "source": src, "chromosome": ch, "orig_seg_id": k,
            "bin_start": int(s), "bin_end": int(e), "seg_length": int(e - s),
            "median_obs": round(med, 5) if np.isfinite(med) else np.nan,
            "standalone_in_simplified": bool(standalone),
            "merged_left": bool(l_rem), "merged_right": bool(r_rem),
            "reason_left": left_reason, "reason_right": right_reason,
        })

    # descriptive SRD rows (never decision inputs) and flank-oriented small
    # segment table (delta_left/right relative to the SEGMENT, per-bin noise)
    srd_rows = []
    flank_rows = []
    if need_boot:
        seg_ids = np.full(len(pb_obs), -1, dtype=int)
        for k, s in enumerate(orig_segs):
            seg_ids[s.start:s.end] = k
        phi = pooled_phi(pb_obs, seg_ids)
        Y = np.array([float(pb_obs[s.start:s.end].sum()) for s in orig_segs])
        E = np.array([float(s.end - s.start) for s in orig_segs])
        snap_by = {int(r["orig_seg_left"]): r for r in snap}
        mobs = mc.segment_medians(pb_obs, orig_segs)

        def _stats(vals):
            vals = np.asarray(vals, dtype=float)
            n = int(len(vals))
            if n == 0:
                return np.nan, np.nan, np.nan, np.nan, 0
            mean = float(vals.mean())
            var = float(vals.var(ddof=1)) if n >= 2 else np.nan
            std = float(np.sqrt(var)) if np.isfinite(var) else np.nan
            cv = float(std / mean) if np.isfinite(std) and np.isfinite(mean) and mean > 0 else np.nan
            return mean, var, std, cv, n

        for k, s, e in small:
            medS = float(mobs[k])
            meanS, var_seg, std_seg, cv_seg, n_seg = _stats(pb_obs[s:e])
            trimS = iqr_trimmed_mean(pb_obs[s:e])
            row = {
                "source": src, "chromosome": ch, "orig_seg_id": k,
                "bin_start": int(s), "bin_end": int(e),
                "seg_length": int(e - s),
                "Mb_start": round(float(positions[s]) / 1e6, 3),
                "Mb_end": round((float(positions[min(e - 1, len(positions) - 1)]) + 100000) / 1e6, 3),
                "median_seg": round(medS, 5) if np.isfinite(medS) else np.nan,
                "mean_seg_perbin": round(meanS, 5) if np.isfinite(meanS) else np.nan,
                "trimmean_seg_perbin": round(trimS, 5) if np.isfinite(trimS) else np.nan,
                "var_seg_perbin": round(var_seg, 6) if np.isfinite(var_seg) else np.nan,
                "std_seg_perbin": round(std_seg, 5) if np.isfinite(std_seg) else np.nan,
                "cv_seg_perbin": round(cv_seg, 5) if np.isfinite(cv_seg) else np.nan,
                "n_bins_seg": n_seg,
            }
            for side, bj in (("left", k - 1), ("right", k)):
                if bj < 0 or bj >= len(orig_segs) - 1:
                    continue
                b_orig = orig_segs[bj]          # seg bj
                b_nxt = orig_segs[bj + 1]       # seg bj+1
                if side == "left":
                    flank_k = b_orig             # neighbour is seg k-1
                    srd = mc.exposure_aware_srd(Y[k], E[k], Y[k - 1], E[k - 1])
                    srd_phi = srd / np.sqrt(phi) if np.isfinite(srd) and phi > 0 else np.nan
                    delta = float(np.log2(medS / mobs[k - 1])) if medS > 0 and mobs[k - 1] > 0 else np.nan
                    fj = k - 1
                else:
                    flank_k = b_nxt
                    srd = mc.exposure_aware_srd(Y[k], E[k], Y[k + 1], E[k + 1])
                    srd_phi = srd / np.sqrt(phi) if np.isfinite(srd) and phi > 0 else np.nan
                    delta = float(np.log2(medS / mobs[k + 1])) if medS > 0 and mobs[k + 1] > 0 else np.nan
                    fj = k + 1
                fmean, fvar, fstd, fcv, fn = _stats(
                    pb_obs[orig_segs[fj].start:orig_segs[fj].end])
                ftrim = iqr_trimmed_mean(pb_obs[orig_segs[fj].start:orig_segs[fj].end])
                snaprow = snap_by.get(bj, {})
                fin_reason = res["decided"].get(bj, "retained_not_eligible")
                merged = isinstance(fin_reason, str) and fin_reason.startswith("removed_iter")
                if side == "left":
                    row["median_left"] = round(float(mobs[k - 1]), 5) if np.isfinite(mobs[k - 1]) else np.nan
                else:
                    row["median_right"] = round(float(mobs[k + 1]), 5) if np.isfinite(mobs[k + 1]) else np.nan
                row[f"delta_{side}"] = round(delta, 4) if np.isfinite(delta) else np.nan
                row[f"srd_phi_{side}"] = round(srd_phi, 4) if np.isfinite(srd_phi) else np.nan
                row[f"frac_under_{side}"] = snaprow.get("frac_under_tau", np.nan)
                row[f"eligible_{side}"] = bool(snaprow.get("eligible", False))
                row[f"merged_{side}"] = bool(merged)
                row[f"reason_{side}"] = fin_reason
                row[f"flank_orig_seg_{side}"] = int(flank_k.orig[0])
                row[f"mean_{side}_flank_perbin"] = round(fmean, 5) if np.isfinite(fmean) else np.nan
                row[f"trimmean_{side}_flank_perbin"] = round(ftrim, 5) if np.isfinite(ftrim) else np.nan
                row[f"var_{side}_flank_perbin"] = round(fvar, 6) if np.isfinite(fvar) else np.nan
                row[f"std_{side}_flank_perbin"] = round(fstd, 5) if np.isfinite(fstd) else np.nan
                row[f"cv_{side}_flank_perbin"] = round(fcv, 5) if np.isfinite(fcv) else np.nan
                row[f"n_bins_{side}_flank"] = int(fn)
                delta_mean = float(np.log2(meanS / fmean)) if meanS > 0 and fmean > 0 else np.nan
                row[f"delta_mean_{side}"] = round(delta_mean, 4) if np.isfinite(delta_mean) else np.nan
                delta_trim = float(np.log2(trimS / ftrim)) if trimS > 0 and ftrim > 0 else np.nan
                row[f"delta_trim_{side}"] = round(delta_trim, 4) if np.isfinite(delta_trim) else np.nan
                srd_rows.append({
                    "source": src, "chromosome": ch, "orig_seg_id": k,
                    "bin_start": int(s), "bin_end": int(e), "side": side,
                    "srd": round(srd, 4) if np.isfinite(srd) else np.nan,
                    "srd_phi": round(srd_phi, 4) if np.isfinite(srd_phi) else np.nan,
                    "pooled_phi": round(phi, 4) if np.isfinite(phi) else np.nan,
                })
            flank_rows.append(row)

    # sensitivity: rerun merges for every (tau, stab) on SAME replicate pb
    sens_rows = []
    for tau in TAUS:
        for stab in STABS:
            r2 = mc.run_merges(pb_obs, pb_boot, acc, n_bins, tau, stab,
                               max_small=args.max_small_bins, need_all_valid=True)
            sens_rows.append({
                "source": src, "chromosome": ch, "tau": tau, "stab": stab,
                "n_orig_boundaries": int(len(orig_segs) - 1),
                "n_merged": int(len(r2["merges"])),
                "n_retained": int(len(orig_segs) - 1 - len(r2["merges"])),
                "n_simplified_segments": int(len(r2["final_segs"])),
            })

    return {
        "skipped": False, "source": src, "chromosome": ch,
        "pb_obs": pb_obs, "positions": positions,
        "orig_segs": orig_segs, "final_segs": final_segs,
        "acc": acc,
        "merges": res["merges"],
        "removed_bins": [int(m["boundary_bin"]) for m in res["merges"]],
        "init_rows": init_rows, "merge_rows": merge_rows,
        "retained_rows": retained_rows, "simplified_rows": simplified_rows,
        "unresolved_rows": unresolved_rows, "srd_rows": srd_rows,
        "flank_rows": flank_rows,
        "sens_rows": sens_rows,
    }


DEPTH_COLORS = ["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725", "#ff7f0e"]
GAP_MB = 0.15  # gap > ~1.5 bins between consecutive retained bins breaks the line


def break_at_gaps(pos, prof):
    """Insert NaN between consecutive retained bins separated by > GAP_MB so the
    profile line does not bridge missing genomic regions."""
    pos = np.asarray(pos, dtype=float)
    prof = np.asarray(prof, dtype=float)
    if len(pos) < 2:
        return pos, prof
    gap = np.diff(pos) > GAP_MB
    if not gap.any():
        return pos, prof
    out_p, out_v = [pos[0]], [prof[0]]
    for i in range(len(pos) - 1):
        if gap[i]:
            out_p.append(np.nan)
            out_v.append(np.nan)
        out_p.append(pos[i + 1])
        out_v.append(prof[i + 1])
    return np.asarray(out_p), np.asarray(out_v)


def figure_chromosome(ch, fig_src, out_fig_dir):
    """Screening-style per-chromosome track (mirrors small_segment_screening).

    fig_src: dict source -> dict(pb, positions_mb, acc_bins, acc_depth,
                                 small_rows, removed_bins).  small_rows carry the
    SRD-screen proposal label per small segment: hatch 'x' = merge_* proposal,
    hatch 'o' = ambiguous_keep; span colour = weaker-flank |SRD/sqrt(phi)|
    (blue < 3.3 reference, red >= 3.3), exactly like the screening figures.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.lines import Line2D
    sources = list(fig_src.keys())
    if not sources:
        return None
    thr = 3.3
    fig, axs = plt.subplots(len(sources), 1, figsize=(15, 2.2 * len(sources)),
                            sharex=True)
    axs = np.atleast_1d(axs)
    # two-hue weaker-flank |SRD/√φ| colour scale (same as screening figures)
    lo_ramp = plt.cm.Blues(np.linspace(0.25, 1.0, 8))
    hi_ramp = plt.cm.Reds(np.linspace(0.35, 1.0, 8))
    cmap = mcolors.ListedColormap(np.vstack([lo_ramp, hi_ramp]))
    bounds = np.concatenate([np.linspace(0, thr, 8), np.linspace(thr, 8.0, 8)[1:]])
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    HATCH = {"merge_left": "x", "merge_right": "x", "ambiguous_keep": "o"}
    for ax, src in zip(axs, sources):
        d = fig_src[src]
        pb, pos = d["pb"], d["positions_mb"]
        med = float(np.median(pb))
        with np.errstate(divide="ignore", invalid="ignore"):
            prof = np.log2(pb / med) if med > 0 else np.full_like(pb, np.nan)
        px, py = break_at_gaps(pos, prof)
        ax.plot(px, py, color="grey", lw=0.7, alpha=0.9, zorder=1)
        ax.axhline(0, color="black", lw=0.5, ls=":")
        # accepted breakpoints coloured by recursion depth (screening style)
        for b, depth in zip(d["acc_bins"], d["acc_depth"]):
            bi = min(int(b), len(pos) - 1)
            col = DEPTH_COLORS[int(depth) % len(DEPTH_COLORS)]
            ax.axvline(pos[bi], color=col, lw=1.1, alpha=0.9, zorder=4)
        # small segments: span colour by weaker-flank |SRD/√φ|, hatch by proposal
        sub = d["small_rows"]
        for _, r in sub.iterrows():
            lo = min(int(r["bin_start"]), len(pos) - 1)
            hi = min(int(r["bin_end"]) - 1, len(pos) - 1)
            vals = [abs(float(v)) for v in (r["srd_phi_left"], r["srd_phi_right"])
                    if pd.notna(v)]
            c = cmap(norm(min(min(vals), 8.0))) if vals else "#bbbbbb"
            h = HATCH.get(str(r["proposal"]), None)
            ax.axvspan(pos[lo], pos[hi], color=c, alpha=0.75,
                       hatch=h, edgecolor="black", lw=0.4, zorder=0)
        # prototype removals (red) when merges do occur
        for b in d.get("removed_bins", []):
            ax.axvline(pos[min(int(b), len(pos) - 1)], color="#d62728",
                       lw=1.8, zorder=5)
        # secondary (top) axis: segment sizes at midpoints, like screening
        segs = mc.flat_partition(sorted(int(b) for b in d["acc_bins"]), len(pb))
        if len(segs) <= 60:
            mids = [pos[min(int((s + e) // 2), len(pos) - 1)] for s, e in
                    [(s.start, s.end) for s in segs]]
            tw = ax.twiny()
            tw.set_xlim(ax.get_xlim())
            tw.set_xticks(mids)
            tw.set_xticklabels([f"{s.end - s.start}" for s in segs], fontsize=5)
            tw.tick_params(axis="x", length=3, labelsize=5, pad=1, top=True,
                           bottom=False, labeltop=True, labelbottom=False)
        ax.set_ylabel(src, fontsize=8)
    axs[-1].set_xlabel(f"{ch} position (Mb)", fontsize=9)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=fig.get_axes(), fraction=0.015, pad=0.01)
    cb.set_label("small-segment weaker-flank |SRD/√φ| (blue < 3.3, red ≥ 3.3)",
                 fontsize=7)
    handles = [
        Line2D([0], [0], color="grey", lw=0.7, label="log2 depth vs chr median"),
        Line2D([0], [0], color="#440154", lw=1.4, label="accepted breakpoint"),
        Line2D([0], [0], color="#2ca02c", lw=0.4, markeredgecolor="black",
               marker="s", markersize=6, markerfacecolor="#bbbbbb",
               label="merge proposal (hatch x)"),
        Line2D([0], [0], color="#2ca02c", lw=0.4, markeredgecolor="black",
               marker="o", markersize=6, markerfacecolor="#bbbbbb",
               label="ambiguous (hatch o)"),
        Line2D([0], [0], color="#d62728", lw=1.8, label="removed (merged) boundary"),
    ]
    fig.legend(handles=handles, loc="upper left", fontsize=8, ncol=3, frameon=True)
    fig.suptitle(f"{ch} — per-cluster depth partition; small-segment spans "
                 "coloured by weaker-flank |SRD/√φ| (SRD-screen proposals hatched)",
                 fontsize=11)
    fig.tight_layout()
    out = out_fig_dir / f"{ch}_small_segments.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out


def figure_chromosome_effect(ch, fig_src, out_fig_dir, delta_col="delta",
                             out_tag="median", tau=0.20, color_metric="delta"):
    """Per-chromosome track with small segments coloured either by the segment's
    min |flank delta| (color_metric='delta') or by the |SRD/√φ| of the flank that
    gives the min delta (color_metric='srd').  Hatches always encode the SRD
    relationship of the min-delta flank.

    Hatch semantics per small segment:
      - pick the flank with min |delta| (that flank's |SRD/√φ| is used);
      - hatch 'x' when that |SRD/√φ| < 3.3 (SRD weak), hatch 'o' when >= 3.3;
      - hatch '/' when the flank giving min |delta| does NOT also give the min
        |SRD/√φ| (delta-min and SRD-min flanks disagree -> inconsistent).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.lines import Line2D
    sources = list(fig_src.keys())
    if not sources:
        return None
    srd_ref = 3.3
    use_srd = color_metric == "srd"

    # collect colour values across sources to set the discrete colourbar top
    cvals = []
    for src in sources:
        sub = fig_src[src].get("small_rows", pd.DataFrame())
        for _, r in sub.iterrows():
            dl = pd.to_numeric(r.get(f"{delta_col}_left"), errors="coerce")
            dr = pd.to_numeric(r.get(f"{delta_col}_right"), errors="coerce")
            if not (np.isfinite(dl) and np.isfinite(dr)):
                continue
            sl = abs(float(pd.to_numeric(r.get("srd_phi_left"), errors="coerce")))
            sr = abs(float(pd.to_numeric(r.get("srd_phi_right"), errors="coerce")))
            if abs(float(dl)) <= abs(float(dr)):
                cvals.append(sl if use_srd else abs(float(dl)))
            else:
                cvals.append(sr if use_srd else abs(float(dr)))
    if use_srd:
        top = max([8.0] + [v for v in cvals if np.isfinite(v)])
        BANDS = [0.0, 1.5, 3.3, 5.0, top]
        LABELS = ["< 1.5", "1.5–3.3", "3.3–5.0", "\u2265 5.0"]
        clabel = "|SRD/√φ| of the min-delta flank"
    else:
        top = max([1.0] + [v for v in cvals if np.isfinite(v)])
        top = max(top, 1.0)
        BANDS = [0.0, 0.20, 0.50, 0.75, 1.00, top]
        LABELS = ["< 0.20", "0.20–0.50", "0.50–0.75", "0.75–1.00", "\u2265 1.00"]
        clabel = f"min(|delta|) of small segment ({delta_col} levels; tau={tau:.2f})"
    band_colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(BANDS) - 1))
    cmap = mcolors.ListedColormap(band_colors)
    bnorm = mcolors.BoundaryNorm(BANDS, cmap.N)
    fig, axs = plt.subplots(len(sources), 1, figsize=(15, 2.2 * len(sources)),
                            sharex=True)
    axs = np.atleast_1d(axs)
    for ax, src in zip(axs, sources):
        d = fig_src[src]
        pb, pos = d["pb"], d["positions_mb"]
        med = float(np.median(pb))
        with np.errstate(divide="ignore", invalid="ignore"):
            prof = np.log2(pb / med) if med > 0 else np.full_like(pb, np.nan)
        px, py = break_at_gaps(pos, prof)
        ax.plot(px, py, color="grey", lw=0.7, alpha=0.9, zorder=1)
        ax.axhline(0, color="black", lw=0.5, ls=":")
        for b, depth in zip(d["acc_bins"], d["acc_depth"]):
            bi = min(int(b), len(pos) - 1)
            ax.axvline(pos[bi], color=DEPTH_COLORS[int(depth) % len(DEPTH_COLORS)],
                       lw=1.1, alpha=0.9, zorder=4)
        sub = d.get("small_rows", pd.DataFrame())
        for _, r in sub.iterrows():
            dl = pd.to_numeric(r.get(f"{delta_col}_left"), errors="coerce")
            dr = pd.to_numeric(r.get(f"{delta_col}_right"), errors="coerce")
            sl = pd.to_numeric(r.get("srd_phi_left"), errors="coerce")
            sr = pd.to_numeric(r.get("srd_phi_right"), errors="coerce")
            if not (np.isfinite(dl) and np.isfinite(dr)):
                continue
            al, ar = abs(float(dl)), abs(float(dr))
            if al <= ar:
                dmin, sval, oth_s = al, abs(float(sl)) if np.isfinite(sl) else np.nan, abs(float(sr)) if np.isfinite(sr) else np.nan
            else:
                dmin, sval, oth_s = ar, abs(float(sr)) if np.isfinite(sr) else np.nan, abs(float(sl)) if np.isfinite(sl) else np.nan
            hatch = None
            if np.isfinite(sval):
                if np.isfinite(oth_s) and oth_s < sval:
                    hatch = "/"  # delta-min flank is not the SRD-min flank
                else:
                    hatch = "x" if sval < srd_ref else "o"
            else:
                hatch = "/"
            cval = sval if use_srd else dmin
            if not np.isfinite(cval):
                c = "#bbbbbb"
            else:
                c = cmap(bnorm(min(cval, top)))
            lo = min(int(r["bin_start"]), len(pos) - 1)
            hi = min(int(r["bin_end"]) - 1, len(pos) - 1)
            ax.axvspan(pos[lo], pos[hi], color=c,
                       alpha=0.85, hatch=hatch, edgecolor="black", lw=0.4, zorder=0)
        for b in d.get("removed_bins", []):
            ax.axvline(pos[min(int(b), len(pos) - 1)], color="#d62728",
                       lw=1.8, zorder=5)
        ax.set_ylabel(src, fontsize=8)
    axs[-1].set_xlabel(f"{ch} position (Mb)", fontsize=9)
    # discrete colourbar with band labels
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=bnorm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=fig.get_axes(), fraction=0.015, pad=0.01)
    centers = [(BANDS[i] + BANDS[i + 1]) / 2 for i in range(len(BANDS) - 1)]
    cb.set_ticks(centers)
    cb.ax.set_yticklabels(LABELS)
    cb.set_label(clabel, fontsize=8)
    from matplotlib.patches import Rectangle
    handles = [
        Line2D([0], [0], color="grey", lw=0.7, label="log2 depth vs chr median"),
        Line2D([0], [0], color="#440154", lw=1.4, label="accepted breakpoint"),
        Rectangle((0, 0), 1, 1, facecolor="#bbbbbb", edgecolor="black",
                  hatch="x", label="hatch x: min-delta flank |SRD/√φ| < 3.3"),
        Rectangle((0, 0), 1, 1, facecolor="#bbbbbb", edgecolor="black",
                  hatch="o", label="hatch o: min-delta flank |SRD/√φ| \u2265 3.3"),
        Rectangle((0, 0), 1, 1, facecolor="#bbbbbb", edgecolor="black",
                  hatch="/", label="hatch /: min-delta flank \u2260 min-SRD flank"),
    ]
    fig.legend(handles=handles, loc="upper left", fontsize=8, ncol=3, frameon=True)
    metric_name = "SRD/√φ" if use_srd else f"min |{delta_col} delta|"
    fig.suptitle(f"{ch} — small segments coloured by {metric_name}, "
                 f"hatched by flank |SRD/√φ| ({out_tag} levels)", fontsize=11)
    fig.tight_layout()
    suffix = f"_{color_metric}"
    out = out_fig_dir / f"{ch}_small_segments_effect_{out_tag}{suffix}.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out


def figure_small_segment_volcano(flank_df, out_fig_dir):
    """Volcano for small segments.

    x = |delta| of the weaker-effect flank; y = |SRD_phi| of that SAME flank
    (same-flank pairing).  Tau and the SRD reference (3.3) are shown as dashed
    lines.  Coloured by the SRD-screen proposal when available.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if flank_df is None or len(flank_df) == 0:
        return None
    df = flank_df.copy()
    dl = pd.to_numeric(df["delta_left"], errors="coerce").abs()
    dr = pd.to_numeric(df["delta_right"], errors="coerce").abs()
    sl = pd.to_numeric(df["srd_phi_left"], errors="coerce").abs()
    sr = pd.to_numeric(df["srd_phi_right"], errors="coerce").abs()
    weaker = dl <= dr
    x = np.where(weaker, dl, dr)
    y = np.where(weaker, sl, sr)
    side = np.where(weaker, "left", "right")
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 2:
        return None
    x, y, side = x[m], y[m], side[m]
    col = {"left": "#1f77b4", "right": "#ff7f0e"}
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for s in ("left", "right"):
        ax.scatter(x[side == s], y[side == s], s=26, alpha=0.75,
                   color=col[s], edgecolor="none", label=f"weaker flank = {s}")
    ax.axvline(0.20, color="black", ls="--", lw=1.0)
    ax.axhline(3.3, color="black", ls="--", lw=1.0)
    ax.set_xlabel("|log2 median effect| of weaker flank")
    ax.set_ylabel("|SRD/√φ| of that same flank")
    ax.set_title("small segments: weaker-flank effect vs its SRD/√φ "
                 "(dashed: tau=0.20, SRD 3.3)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = out_fig_dir / "small_segment_volcano.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out


def figure_srd_log2fc_noise(flank_df, out_fig_dir, delta_col="delta",
                            out_name="srd_log2fc_noise_relationship.png",
                            delta_label="median"):
    """Three panels: SRD_phi vs log2(delta) relationship coloured by
    max(segment, flank-of-that-side) variance / std / CV of per-bin depth.

    delta_col selects whether the log2 effect is computed from median or mean
    per-bin segment levels (flank table columns delta_* vs delta_mean_*).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    if flank_df is None or len(flank_df) == 0:
        return None
    df = flank_df.copy()
    if not {f"{delta_col}_left", f"{delta_col}_right"}.issubset(df.columns):
        return None
    rows = []
    for _, r in df.iterrows():
        for side in ("left", "right"):
            d = pd.to_numeric(r[f"{delta_col}_{side}"], errors="coerce")
            s = pd.to_numeric(r[f"srd_phi_{side}"], errors="coerce")
            if pd.notna(d) and pd.notna(s):
                seg_m = pd.to_numeric(r.get("mean_seg_perbin"), errors="coerce")
                seg_t = pd.to_numeric(r.get("std_seg_perbin"), errors="coerce")
                seg_c = pd.to_numeric(r.get("cv_seg_perbin"), errors="coerce")
                flk_m = pd.to_numeric(r.get(f"mean_{side}_flank_perbin"), errors="coerce")
                flk_t = pd.to_numeric(r.get(f"std_{side}_flank_perbin"), errors="coerce")
                flk_c = pd.to_numeric(r.get(f"cv_{side}_flank_perbin"), errors="coerce")

                def _max2(a, b):
                    if pd.isna(a):
                        return b
                    if pd.isna(b):
                        return a
                    return max(float(a), float(b))

                rows.append({
                    "source": r["source"], "chromosome": r["chromosome"],
                    "orig_seg_id": r["orig_seg_id"], "side": side,
                    "delta": float(d), "srd_phi_abs": abs(float(s)),
                    "mean": _max2(seg_m, flk_m),
                    "std": _max2(seg_t, flk_t),
                    "cv": _max2(seg_c, flk_c),
                })
    if len(rows) < 2:
        return None
    long = pd.DataFrame(rows)

    def _lin_norm(values):
        v = pd.to_numeric(values, errors="coerce")
        v = v[np.isfinite(v)]
        if len(v) == 0:
            return None
        vmin, vmax = float(v.min()), float(v.max())
        if vmin == vmax:
            vmax = vmin + 1e-9
        return vmin, vmax

    import matplotlib.colors as _mc

    CV_BOUNDS = [0.0, 0.6, 0.8, 1.0, 1.2]
    CV_LABELS = ["< 0.6", "0.6–0.8", "0.8–1.0", "1.0–1.2", "≥ 1.2"]
    fig, axs = plt.subplots(1, 3, figsize=(19, 5.8))
    panels = [
        ("std", "std", "max(seg,flank) std (depth, per-bin)"),
        ("mean", "mean", "max(seg,flank) mean (depth, per-bin)"),
        ("cv", "CV", "max(seg,flank) CV = std / mean (per-bin)"),
    ]
    for ax, (metric, key, title) in zip(axs, panels):
        vv = pd.to_numeric(long[metric], errors="coerce")
        mm = vv.notna() & long["srd_phi_abs"].notna() & long["delta"].notna()
        sub = long[mm]
        sv = pd.to_numeric(sub[metric], errors="coerce")
        if metric == "cv":
            top = max(float(sv.max()), 1.2 + 1e-9)
            bounds = [0.0, 0.6, 0.8, 1.0, 1.2, top]
            colors = plt.cm.viridis(np.linspace(0.15, 1.0, len(bounds) - 1))
            cmap = _mc.ListedColormap(colors)
            norm = _mc.BoundaryNorm(bounds, cmap.N)
            sc = ax.scatter(sub["delta"], sub["srd_phi_abs"], c=sv, cmap=cmap,
                            norm=norm, s=24, alpha=0.9, edgecolor="none")
            centers = [(bounds[i] + bounds[i + 1]) / 2
                       for i in range(len(bounds) - 1)]
            cb = fig.colorbar(sc, ax=ax, ticks=centers)
            cb.ax.set_yticklabels(CV_LABELS)
            cb.set_label("CV = std / mean (per-bin, binned)")
        else:
            _res = _lin_norm(sv)
            if _res is None:
                ax.set_title(f"{title}\n(no data)")
                continue
            vmin, vmax = _res
            sc = ax.scatter(sub["delta"], sub["srd_phi_abs"], c=sv, cmap="viridis",
                            s=24, vmin=vmin, vmax=vmax, alpha=0.85, edgecolor="none")
            cb = fig.colorbar(sc, ax=ax)
            cb.set_label(title)
        ax.axhline(3.3, color="grey", ls="--", lw=0.8)
        ax.axvline(0.2, color="grey", ls="--", lw=0.8)
        ax.axvline(-0.2, color="grey", ls="--", lw=0.8)
        ax.set_xlabel(f"log2 {delta_label} effect (delta, signed)")
        ax.set_ylabel("|SRD/√φ|")
        ax.set_title(f"colour = {title}")
    fig.suptitle(f"small segments — |SRD_phi| vs log2 effect ({delta_label} levels)",
                 fontsize=13)
    fig.tight_layout()
    out = out_fig_dir / out_name
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.figures_dir) if args.figures_dir else out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    adata = ad.read_h5ad(args.adata)
    raw = loader.read_breakpoints(args.breakpoints)
    _, members = loader.load_annotations(adata, args.integrated, args.cluster_key)
    chrom_of, positions_all = loader.chrom_bin_positions(adata)
    X = adata.X
    sources = [f"cluster{i}" for i in range(5)]

    repro = reproduce_check(adata, members, raw, chrom_of, args.reference_scores,
                            args.chromosomes, args.max_small_bins)
    log(f"reproduction check: {repro}")

    cluster_reps = {}
    for ci, src in enumerate(sources):
        reps = loader.replicate_draws(len(members[src]), args.n_bootstrap,
                                      args.seed, ci)
        cluster_reps[src] = {"reps": reps, "hash": loader.draw_hash(reps)}

    init_rows, merge_rows, retained_rows = [], [], []
    simplified_rows, unresolved_rows, srd_rows, sens_rows = [], [], [], []
    flank_rows = []
    chr_figs = {}
    eff_figs = {}
    n_processed = 0

    # SRD-screen proposals for the hatches (from the screening reference)
    if args.reference_scores and Path(args.reference_scores).exists():
        screen = pd.read_csv(args.reference_scores, sep="\t")
    else:
        screen = None

    for ch in args.chromosomes:
        if ch not in set(chrom_of):
            log(f"skip {ch}: not present")
            continue
        src_results = {}
        fig_src = {}
        for src in sources:
            out = process_source_chr(args, X, members, raw, src, ch, chrom_of,
                                     positions_all, cluster_reps)
            if out.get("skipped"):
                log(f"{src} {ch}: no segments / no boundaries -> skip")
                continue
            n_processed += 1
            init_rows.extend(out["init_rows"])
            merge_rows.extend(out["merge_rows"])
            retained_rows.extend(out["retained_rows"])
            simplified_rows.extend(out["simplified_rows"])
            unresolved_rows.extend(out["unresolved_rows"])
            srd_rows.extend(out["srd_rows"])
            flank_rows.extend(out["flank_rows"])
            sens_rows.extend(out["sens_rows"])
            src_results[src] = out
            # assemble screening-style figure payload
            acc_raw = raw[(raw["source"] == src) & (raw["chromosome"] == ch)]
            acc_bins = sorted(int(b) for b in acc_raw["absolute_bin"])
            acc_depth = [int(d) for _, d in acc_raw.set_index("absolute_bin")
                         ["recursion_depth"].sort_index().items()]
            small_rows = pd.DataFrame()
            if screen is not None:
                sub = screen[(screen["source"] == src) & (screen["chromosome"] == ch)]
                small_rows = sub.sort_values("bin_start")[
                    ["bin_start", "bin_end", "srd_phi_left", "srd_phi_right",
                     "proposal"]].reset_index(drop=True)
            fig_src[src] = {
                "pb": out["pb_obs"],
                "positions_mb": out["positions"] / 1e6,
                "acc_bins": acc_bins,
                "acc_depth": acc_depth,
                "small_rows": small_rows,
                "removed_bins": out["removed_bins"],
            }
        # whole-sample P_all context row (bottom), like the screening figures
        cols = np.where(chrom_of == ch)[0]
        acc_raw = raw[(raw["source"] == "all") & (raw["chromosome"] == ch)]
        acc_bins = sorted(int(b) for b in acc_raw["absolute_bin"])
        acc_depth = [int(d) for _, d in acc_raw.set_index("absolute_bin")
                     ["recursion_depth"].sort_index().items()]
        small_rows = pd.DataFrame()
        if screen is not None:
            sub = screen[(screen["source"] == "all") & (screen["chromosome"] == ch)]
            small_rows = sub.sort_values("bin_start")[
                ["bin_start", "bin_end", "srd_phi_left", "srd_phi_right",
                 "proposal"]].reset_index(drop=True)
        fig_src["all"] = {
            "pb": loader.observed_pseudobulk(X, list(range(X.shape[0])), cols),
            "positions_mb": positions_all[cols] / 1e6,
            "acc_bins": acc_bins,
            "acc_depth": acc_depth,
            "small_rows": small_rows,
            "removed_bins": [],
        }
        fig_path = figure_chromosome(ch, fig_src, fig_dir)
        if fig_path:
            chr_figs[ch] = str(fig_path)

        # effect figures (median & mean delta): per-cluster rows only, built
        # from the flank rows which carry both delta and delta_mean
        flank_ch = [r for r in flank_rows if r["chromosome"] == ch]
        if flank_ch:
            fch = pd.DataFrame(flank_ch)
            eff_src = {}
            for src in sources:
                sub = fch[fch["source"] == src]
                if not len(sub) or src not in src_results:
                    continue
                acc_raw = raw[(raw["source"] == src) & (raw["chromosome"] == ch)]
                acc_bins = sorted(int(b) for b in acc_raw["absolute_bin"])
                acc_depth = [int(d) for _, d in acc_raw.set_index("absolute_bin")
                             ["recursion_depth"].sort_index().items()]
                eff_src[src] = {
                    "pb": src_results[src]["pb_obs"],
                    "positions_mb": src_results[src]["positions"] / 1e6,
                    "acc_bins": acc_bins,
                    "acc_depth": acc_depth,
                    "small_rows": sub.sort_values("bin_start")[
                        ["bin_start", "bin_end", "delta_left", "delta_right",
                         "delta_mean_left", "delta_mean_right",
                         "srd_phi_left", "srd_phi_right"]].reset_index(drop=True),
                    "removed_bins": src_results[src]["removed_bins"],
                }
            for tag, dcol, cmetric in (("median", "delta", "delta"),
                                       ("mean", "delta_mean", "delta"),
                                       ("mean", "delta_mean", "srd")):
                if not eff_src:
                    continue
                p = figure_chromosome_effect(ch, eff_src, fig_dir,
                                             delta_col=dcol, out_tag=tag,
                                             color_metric=cmetric)
                if p:
                    eff_figs.setdefault(ch, {})[f"{tag}_{cmetric}"] = str(p)
                    log(f"wrote {p}")

    def w(name, rows, columns=None):
        df = pd.DataFrame(rows, columns=columns)
        df.to_csv(out_dir / name, sep="\t", index=False)
        log(f"wrote {out_dir/name} ({len(rows)} rows)")

    MERGE_COLS = ["source", "chromosome", "iteration", "orig_boundary_id",
                  "boundary_bin", "delta_obs", "frac_under_tau", "n_valid",
                  "guard_range", "left_orig", "right_orig"]
    w("initial_boundary_scores.tsv", init_rows)
    w("merge_history.tsv", merge_rows, columns=MERGE_COLS)
    w("retained_boundaries.tsv", retained_rows)
    w("simplified_segments.tsv", simplified_rows)
    w("unresolved_segments.tsv", unresolved_rows)
    w("srd_descriptive.tsv", srd_rows)
    w("small_segment_flank_scores.tsv", flank_rows)
    w("threshold_sensitivity.tsv", sens_rows)

    flank_df = pd.DataFrame(flank_rows) if flank_rows else None
    vfig = figure_small_segment_volcano(flank_df, fig_dir) if flank_df is not None else None
    nfig = figure_srd_log2fc_noise(flank_df, fig_dir, delta_col="delta",
                                   out_name="srd_log2fc_noise_relationship.png",
                                   delta_label="median") if flank_df is not None else None
    nfig_mean = figure_srd_log2fc_noise(
        flank_df, fig_dir, delta_col="delta_mean",
        out_name="srd_log2fc_noise_relationship_mean.png",
        delta_label="mean") if flank_df is not None else None
    if vfig:
        log(f"wrote {vfig}")
    if nfig:
        log(f"wrote {nfig}")
    if nfig_mean:
        log(f"wrote {nfig_mean}")

    mdf = pd.DataFrame(merge_rows, columns=MERGE_COLS)
    log(f"total merges (tau={args.tau}, stab={args.stab}): {len(mdf)}")
    if len(mdf):
        log("merges per source: " + json.dumps(mdf.groupby("source").size().to_dict()))
    if init_rows:
        isdf = pd.DataFrame(init_rows)
        cand = isdf[isdf["candidate"]]
        log("candidate boundaries: %d" % len(cand))
        log("candidate reasons: " + json.dumps(cand["reason"].value_counts().to_dict()))

    # compact review summary per (source, chromosome)
    if init_rows:
        isdf = pd.DataFrame(init_rows)
        mr_df = pd.DataFrame(merge_rows)
        rev = []
        for (src, ch), grp in isdf.groupby(["source", "chromosome"]):
            cand = grp[grp["candidate"]]
            ad = cand["delta_obs"].abs().dropna()
            sub_mr = (mr_df[(mr_df["source"] == src) & (mr_df["chromosome"] == ch)]
                      if len(mr_df) else mr_df)
            rev.append({
                "source": src, "chromosome": ch,
                "n_orig_segments": int(grp["orig_seg_left"].max()) + 2,
                "n_boundaries": int(len(grp)),
                "n_nominated": int(len(cand)),
                "n_eligible_default": int((cand["eligible"]).sum()),
                "n_merged_default": int(len(sub_mr)),
                "min_abs_delta_nominated": round(float(ad.min()), 4) if len(ad) else np.nan,
                "median_abs_delta_nominated": round(float(ad.median()), 4) if len(ad) else np.nan,
                "max_frac_under_default": round(float(cand["frac_under_tau"].max()), 4) if len(cand) else np.nan,
            })
        w("review_summary.tsv", rev)

    cluster_sizes = {s: len(members[s]) for s in sources}
    extra = {
        "reproduction_check": repro,
        "cluster_sizes": cluster_sizes,
        "cluster_replicate_hashes": {s: cluster_reps[s]["hash"] for s in sources},
        "n_processed_source_chrom": n_processed,
        "totals_default": {
            "merged_boundaries": len(merge_rows),
            "retained_boundaries": len(retained_rows),
            "simplified_segments": len(simplified_rows),
        },
        "sensitivity_grid_max_merges": max((r["n_merged"] for r in sens_rows),
                                           default=0),
        "figures": chr_figs,
        "effect_figures": eff_figs,
        "summary_figures": {
            "small_segment_volcano": str(fig_dir / "small_segment_volcano.png")
            if vfig else None,
            "srd_log2fc_noise_relationship": str(fig_dir / "srd_log2fc_noise_relationship.png")
            if nfig else None,
            "srd_log2fc_noise_relationship_mean": str(fig_dir / "srd_log2fc_noise_relationship_mean.png")
            if nfig_mean else None,
        },
        "source_modules": {
            "merge_core": fingerprint(PKG / "merge_core.py"),
            "loader": fingerprint(PKG / "loader.py"),
            "runner": fingerprint(PKG / "run_depth_segment_prototype.py"),
        },
    }
    meta = {
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "description": "per-cluster depth-segment merge prototype, 7-chr pilot",
        "inputs": {
            "adata": args.adata,
            "integrated": args.integrated,
            "breakpoints": args.breakpoints,
        },
        "params": {
            "chromosomes": args.chromosomes,
            "max_small_bins": args.max_small_bins,
            "tau_default": args.tau,
            "stab_default": args.stab,
            "n_bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "sensitivity_taus": TAUS,
            "sensitivity_stabs": STABS,
        },
        "notes": "tau/stab are engineering settings, not calibrated p-values. "
                 "Bootstrap fractions are stability scores, not probabilities.",
    }
    meta.update(extra)
    (out_dir / "run_manifest.json").write_text(json.dumps(meta, indent=2))
    log(f"wrote {out_dir/'run_manifest.json'}")
    log("DONE.")


if __name__ == "__main__":
    main()
