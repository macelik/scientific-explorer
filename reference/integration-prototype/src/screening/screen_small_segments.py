"""screen_small_segments.py — per-source small-segment screening (no union).

For each source (clusters 0-4 + all) and each chromosome, build the source's OWN
flat partition from its ACCEPTED breakpoints (no cross-source union) and score
every small internal segment (<= max-small-bins, two neighbours):

  median_level            per-bin median of the source's mean pseudobulk
  Y / exposure / rate     sum-of-bins level, bin count, rate = Y/E
  delta_left/right_median log2(median_S / median_L|R)   (primary effect, median)
  log2fc_left/right       log2(rate_S / rate_L|R)       (mean-based rate)
  srd_left/right          exposure-aware SRD (sum-based, Poisson deviance)
  srd_phi_left/right      srd / sqrt(pooled_phi(source, chr))

Plots (like the exposure_aware_srd_pruning dry run, but on the constrained
per-source segmentation): one row per source per chromosome with the log2-depth
track, ACCEPTED breakpoints, REJECTED candidates (constraint vs p-value), and
eligible small segments highlighted orange. A score-viz figure shows SRD and
|median-delta| distributions with the SRD 3.3 / log2FC 0.20 reference lines so
a screening threshold can be inspected. No boundary is pruned or merged here.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BETA = HERE.parent.parent
for d in (str(HERE), str(BETA / "analysis/exposure_aware_srd_pruning"),
          str(BETA / "analysis/flank_bootstrap_diagnostics")):
    if d not in sys.path:
        sys.path.insert(0, d)

from srd_pruning import assign_proposal, exposure_aware_srd  # noqa: E402
from flank_dispersion import pooled_phi  # noqa: E402
from flank_bootstrap import bootstrap_indices, replicate_mean_sums  # noqa: E402

warnings.filterwarnings("ignore")

SRD_REF = 3.3
FC_REF = 0.20
DEPTH_COLORS = ["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725", "#ff7f0e"]


def log(msg):
    print(f"[screen] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser(description="per-source small-segment screening (no union, no pruning)")
    p.add_argument("--adata", required=True)
    p.add_argument("--integrated", required=True)
    p.add_argument("--cluster-key", required=True)
    p.add_argument("--breakpoints", required=True)   # raw_accepted_breakpoints.tsv
    p.add_argument("--audit", required=True)          # segmentation_candidate_audit.tsv
    p.add_argument("--output-dir", required=True)
    p.add_argument("--chromosomes", nargs="+",
                   default=["chr1", "chr2", "chr5", "chr7", "chr8", "chr9", "chr13"])
    p.add_argument("--max-small-bins", type=int, default=100)
    p.add_argument("--n-bootstrap", type=int, default=200,
                   help="within-source bootstrap replicates for stability (diagnostic)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def load(args):
    import anndata as ad
    import muon as mu
    adata = ad.read_h5ad(args.adata)
    if "cellID" in adata.obs.columns:
        adata.obs.index = [str(c) for c in adata.obs["cellID"]]
    else:
        adata.obs.index = [str(c) for c in adata.obs.index]
    m = mu.read_h5mu(args.integrated)
    cl = m.obs[args.cluster_key].astype(str)
    cl.index = [str(c) for c in cl.index]
    cell_to_cluster = {str(c): f"cluster{l}" for c, l in cl.items()}
    missing = [c for c in adata.obs.index if c not in cell_to_cluster]
    if missing:
        raise ValueError(f"{len(missing)} cells absent from cluster assignment")
    labels = np.array([cell_to_cluster[c] for c in adata.obs.index])
    members = {f"cluster{i}": list(np.where(labels == f"cluster{i}")[0]) for i in range(5)}
    members["all"] = list(range(adata.shape[0]))
    raw = pd.read_csv(args.breakpoints, sep="\t")
    audit = pd.read_csv(args.audit, sep="\t")
    return adata, members, raw, audit


def flat_segments(accepted_bins, n_bins):
    cuts = [0] + sorted(int(b) for b in accepted_bins) + [n_bins]
    return [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1) if cuts[i + 1] > cuts[i]]


def score_source_chromosome(pb, segs, source, chromosome, positions,
                            max_small, src_cells):
    """pb: per-bin source mean pseudobulk on this chromosome. segs: own flat partition."""
    n = len(segs)
    E = np.array([e - s for s, e in segs])
    Y = np.array([float(pb[s:e].sum()) for s, e in segs])
    rate = Y / E
    med = np.array([float(np.median(pb[s:e])) for s, e in segs])
    # pooled dispersion over this source-chromosome's own segments
    seg_ids = np.full(len(pb), -1, dtype=int)
    for k, (s, e) in enumerate(segs):
        seg_ids[s:e] = k
    phi = pooled_phi(pb, seg_ids)
    rows = []
    for k in range(1, n - 1):  # internal segments only
        s, e = segs[k]
        length = e - s
        if length > max_small:
            continue
        kl, kr = k - 1, k + 1
        mS = med[k]
        with np.errstate(divide="ignore", invalid="ignore"):
            dL = float(np.log2(mS / med[kl])) if mS > 0 and med[kl] > 0 else np.nan
            dR = float(np.log2(mS / med[kr])) if mS > 0 and med[kr] > 0 else np.nan
            fcL = float(np.log2(rate[k] / rate[kl])) if rate[k] > 0 and rate[kl] > 0 else np.nan
            fcR = float(np.log2(rate[k] / rate[kr])) if rate[k] > 0 and rate[kr] > 0 else np.nan
        srdL = exposure_aware_srd(Y[k], E[k], Y[kl], E[kl]) if E[k] and E[kl] else np.nan
        srdR = exposure_aware_srd(Y[k], E[k], Y[kr], E[kr]) if E[k] and E[kr] else np.nan
        # dispersion-scaled SRD (quasi-Poisson reference) is what we score with
        spL = srdL / np.sqrt(phi) if np.isfinite(srdL) and phi > 0 else np.nan
        spR = srdR / np.sqrt(phi) if np.isfinite(srdR) and phi > 0 else np.nan
        # proposal exactly like exposure_aware_srd_pruning, but scored on
        # SRD/√phi and the MEDIAN effect (no union, no merging performed here)
        if np.isfinite(spL) and np.isfinite(spR) and np.isfinite(dL) and np.isfinite(dR):
            prop, pb_bin = assign_proposal(spL, spR, dL, dR, SRD_REF, FC_REF, int(s), int(e))
        else:
            prop, pb_bin = "keep_focal", None
        rows.append({
            "source": source, "chromosome": chromosome, "n_cells": src_cells,
            "seg_id": k,
            "bin_start": int(s), "bin_end": int(e), "seg_length": int(length),
            "Mb_start": round(float(positions[s]), 3),
            "Mb_end": round(float(positions[e - 1] + 0.1), 3),
            "mean_rate": round(float(rate[k]), 4),
            "median_level": round(float(mS), 4),
            "median_level_left": round(float(med[kl]), 4),
            "median_level_right": round(float(med[kr]), 4),
            "delta_left_median": round(dL, 4) if np.isfinite(dL) else np.nan,
            "delta_right_median": round(dR, 4) if np.isfinite(dR) else np.nan,
            "log2fc_left": round(fcL, 4) if np.isfinite(fcL) else np.nan,
            "log2fc_right": round(fcR, 4) if np.isfinite(fcR) else np.nan,
            "srd_left": round(srdL, 4) if np.isfinite(srdL) else np.nan,
            "srd_right": round(srdR, 4) if np.isfinite(srdR) else np.nan,
            "srd_phi_left": round(spL, 4) if np.isfinite(spL) else np.nan,
            "srd_phi_right": round(spR, 4) if np.isfinite(spR) else np.nan,
            "pooled_phi": round(phi, 4) if np.isfinite(phi) else np.nan,
            "proposal": prop,
            "proposed_boundary_bin": int(pb_bin) if pb_bin is not None else np.nan,
        })
    return rows, segs


def _summarize_boot(obs, rep):
    """SE (ddof=1), 2.5/97.5 percentiles, and sign-stability of replicate stats."""
    rep = np.asarray(rep, dtype=float)
    fin = np.isfinite(rep)
    if fin.sum() < 2:
        return np.nan, np.nan, np.nan, np.nan
    se = float(np.std(rep[fin], ddof=1))
    lo, hi = np.percentile(rep[fin], [2.5, 97.5])
    frac = float(np.mean(np.sign(rep[fin]) == np.sign(obs))) if np.isfinite(obs) else np.nan
    return se, float(lo), float(hi), frac


def bootstrap_source_segments(X, chrom_of, members, source, ch, segs, E, n_boot, seed,
                              observed_rows):
    """Per-source stability: resample the source's cells, re-score log2fc and SRD at
    the source's OWN segment coordinates. Returns list of per-segment dicts."""
    cidx = np.where(chrom_of == ch)[0]
    block = X[:, cidx]
    if hasattr(block, "toarray"):
        dens = np.asarray(block.toarray(), dtype=float)
    else:
        dens = np.asarray(block, dtype=float)
    cum = np.cumsum(dens, axis=1)
    n_seg = len(segs)
    q = np.zeros((X.shape[0], n_seg), dtype=float)
    for i, (s, e) in enumerate(segs):
        if e <= 0:
            continue
        q[:, i] = cum[:, e - 1]
        if s > 0:
            q[:, i] -= cum[:, s - 1]
    # constituents: cluster0..4 always (for 'all' resample each and weight)
    member_map = {f"cluster{i}": list(members[f"cluster{i}"]) for i in range(5)}
    indices = bootstrap_indices(member_map, n_boot, seed)
    sizes = {f"cluster{i}": len(members[f"cluster{i}"]) for i in range(5)}
    mats = replicate_mean_sums(q, indices, sizes, n_boot)
    if source == "all":
        Y_obs = q.mean(axis=0)
        Y_rep = mats["all"]
    else:
        src_rows = members[source]
        Y_obs = q[src_rows].mean(axis=0)
        Y_rep = mats[source]
    rows = []
    for k in range(1, n_seg - 1):
        s, e = segs[k]
        if e - s > 100:
            continue
        kl, kr = k - 1, k + 1
        obs = next((r for r in observed_rows if r["bin_start"] == s and r["bin_end"] == e), None)
        if obs is None:
            continue
        E_k = E[k]
        rate_obs = Y_obs / E
        # replicate rates/log2fc
        with np.errstate(divide="ignore", invalid="ignore"):
            fcl_rep = np.where((Y_rep[:, k] > 0) & (Y_rep[:, kl] > 0),
                               np.log2((Y_rep[:, k] / E[k]) / (Y_rep[:, kl] / E[kl])), np.nan)
            fcr_rep = np.where((Y_rep[:, k] > 0) & (Y_rep[:, kr] > 0),
                               np.log2((Y_rep[:, k] / E[k]) / (Y_rep[:, kr] / E[kr])), np.nan)
        srds = []
        for side, kr2 in (("L", kl), ("R", kr)):
            v = np.empty(n_boot)
            for b in range(n_boot):
                try:
                    v[b] = exposure_aware_srd(Y_rep[b, k], E[k], Y_rep[b, kr2], E[kr2])
                except ValueError:
                    v[b] = np.nan
            srds.append(v)
        se_fcl, lo_fcl, hi_fcl, sf_fcl = _summarize_boot(obs["log2fc_left"], fcl_rep)
        se_fcr, lo_fcr, hi_fcr, sf_fcr = _summarize_boot(obs["log2fc_right"], fcr_rep)
        se_sl, lo_sl, hi_sl, sf_sl = _summarize_boot(obs["srd_left"], srds[0])
        se_sr, lo_sr, hi_sr, sf_sr = _summarize_boot(obs["srd_right"], srds[1])
        rows.append({
            "source": source, "chromosome": ch, "seg_id": obs["seg_id"],
            "bin_start": s, "bin_end": e, "seg_length": e - s,
            "se_log2fc_left": se_fcl, "ci_lo_log2fc_left": lo_fcl, "ci_hi_log2fc_left": hi_fcl,
            "sign_agree_log2fc_left": sf_fcl,
            "se_log2fc_right": se_fcr, "ci_lo_log2fc_right": lo_fcr, "ci_hi_log2fc_right": hi_fcr,
            "sign_agree_log2fc_right": sf_fcr,
            "se_srd_left": se_sl, "ci_lo_srd_left": lo_sl, "ci_hi_srd_left": hi_sl,
            "sign_agree_srd_left": sf_sl,
            "se_srd_right": se_sr, "ci_lo_srd_right": lo_sr, "ci_hi_srd_right": hi_sr,
            "sign_agree_srd_right": sf_sr,
        })
    return rows


def make_chromosome_figure(fig_dir, args, ch, chrom_of, var, X, members, raw, scores):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cidx = np.where(chrom_of == ch)[0]
    pos = var.iloc[cidx]["start"].to_numpy() / 1e6
    sources = [f"cluster{i}" for i in range(5)] + ["all"]
    n_cells = {s: len(members[s]) for s in sources}
    fig, axs = plt.subplots(len(sources), 1, figsize=(15, 2.2 * len(sources)), sharex=True)
    axs = np.atleast_1d(axs)
    # colour small-segment spans by their WEAKER flank |SRD/√φ| (min of the two),
    # with a two-hue scale that clearly separates below vs above the 3.3 reference
    import matplotlib.colors as mcolors
    thr = SRD_REF
    wv = scores[["srd_phi_left", "srd_phi_right"]].abs().min(axis=1)
    vmax = max(8.0, float(wv.max())) if len(scores) else 8.0
    n_lo, n_hi = 8, 8
    lo_ramp = plt.cm.Blues(np.linspace(0.25, 1.0, n_lo))    # below 3.3 -> blue
    hi_ramp = plt.cm.Reds(np.linspace(0.35, 1.0, n_hi))     # at/above 3.3 -> red
    cmap = mcolors.ListedColormap(np.vstack([lo_ramp, hi_ramp]))
    bounds = np.concatenate([np.linspace(0, thr, n_lo),
                             np.linspace(thr, vmax, n_hi)[1:]])
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    for ax, src in zip(axs, sources):
        if src == "all":
            pb = np.asarray(X.mean(axis=0)).ravel()
        else:
            pb = np.asarray(X[members[src]].mean(axis=0)).ravel()
        pb = pb[cidx]
        chrm = float(np.median(pb))
        with np.errstate(divide="ignore", invalid="ignore"):
            prof = np.log2(pb / chrm)
        ax.plot(pos, prof, color="grey", lw=0.7, alpha=0.9, zorder=1)
        ax.axhline(0, color="black", lw=0.5, ls=":")
        # accepted breakpoints (this source): vertical lines coloured by recursion depth
        acc = raw[(raw["source"] == src) & (raw["chromosome"] == ch)]
        for _, r in acc.iterrows():
            b = min(int(r["absolute_bin"]), len(pos) - 1)
            col = DEPTH_COLORS[int(r["recursion_depth"]) % len(DEPTH_COLORS)]
            ax.axvline(pos[b], color=col, lw=1.1, alpha=0.9, zorder=4)
        # small segments highlighted: span colour by weaker-flank |SRD/√φ|
        sub = scores[(scores["source"] == src) & (scores["chromosome"] == ch)]
        for _, r in sub.iterrows():
            lo = min(int(r["bin_start"]), len(pos) - 1)
            hi = min(int(r["bin_end"]) - 1, len(pos) - 1)
            vals = [abs(float(v)) for v in (r["srd_phi_left"], r["srd_phi_right"])
                    if pd.notna(v)]
            c = cmap(norm(min(min(vals), vmax))) if vals else "#bbbbbb"
            ax.axvspan(pos[lo], pos[hi], color=c, alpha=0.75, zorder=0)
        # secondary (top) axis: segment size in bins, ticked at segment midpoints
        segs = flat_segments(sorted(acc["absolute_bin"].astype(int)), len(cidx))
        if len(segs) <= 60:
            mids = [pos[min(int((s + e) // 2), len(pos) - 1)] for s, e in segs]
            tw = ax.twiny()
            tw.set_xlim(ax.get_xlim())
            tw.set_xticks(mids)
            tw.set_xticklabels([f"{e - s}" for s, e in segs], fontsize=5)
            tw.tick_params(axis="x", length=3, labelsize=5, pad=1, top=True,
                           bottom=False, labeltop=True, labelbottom=False)
        ax.set_ylabel(f"{src}\n(n={n_cells[src]})", fontsize=8)
    axs[-1].set_xlabel(f"{ch} position (Mb)", fontsize=9)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=fig.get_axes(), fraction=0.015, pad=0.01)
    cb.set_label("small-segment weaker-flank |SRD/√φ| (blue < 3.3 = not isolated, "
                 "red \u2265 3.3)", fontsize=7)
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color="grey", lw=0.7, label="log2 depth vs chr median"),
        Line2D([0], [0], color="#440154", lw=1.4, label="accepted breakpoint"),
        Line2D([0], [0], color="#ff7f0e", lw=6, alpha=0.75,
               label="small segment (\u2264 100 bins)"),
    ]
    fig.legend(handles=handles, loc="upper left", fontsize=8, ncol=3, frameon=True)
    fig.suptitle(f"{ch} — per-source small segments (no union); span colour = weaker-"
                 "flank |SRD/√φ| (red = both flanks \u2265 3.3, blue = at least one weak)",
                 fontsize=11)
    fig.tight_layout()
    out = fig_dir / f"{ch}_small_segment_screening.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    log(f"saved {out}")


def score_viz(fig_dir, scores):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    ax = axs[0]
    ax.scatter(scores["srd_left"], scores["srd_right"], s=10, alpha=0.5, color="#1f77b4")
    ax.axhline(SRD_REF, color="grey", ls="--", lw=0.8)
    ax.axhline(-SRD_REF, color="grey", ls="--", lw=0.8)
    ax.axvline(SRD_REF, color="grey", ls="--", lw=0.8)
    ax.axvline(-SRD_REF, color="grey", ls="--", lw=0.8)
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("SRD_left")
    ax.set_ylabel("SRD_right")
    ax.set_title("small-segment SRD (both flanks) — ref ±3.3")
    ax = axs[1]
    ax.scatter(scores["delta_left_median"].abs(), scores["delta_right_median"].abs(),
               s=10, alpha=0.5, color="#2ca02c")
    ax.axhline(FC_REF, color="grey", ls="--", lw=0.8)
    ax.axvline(FC_REF, color="grey", ls="--", lw=0.8)
    ax.set_xlabel("|delta_left_median| (log2)")
    ax.set_ylabel("|delta_right_median| (log2)")
    ax.set_title("small-segment median-effect magnitude — ref ±0.20")
    fig.tight_layout()
    out = fig_dir / "small_segment_srd_fc_scatter.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    log(f"saved {out}")


def make_srd_scatter(fig_dir, scores):
    """Descriptive SRD and SRD/√φ, left vs right, coloured by source."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    src_col = {f"cluster{i}": plt.get_cmap("tab10")(i) for i in range(5)}
    src_col["all"] = "black"
    order = [f"cluster{i}" for i in range(5)] + ["all"]
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    for ax, (colx, coly, title) in zip(
            axs, [("srd_left", "srd_right", "descriptive SRD"),
                  ("srd_phi_left", "srd_phi_right", "SRD/√φ (dispersion-scaled)")]):
        for src in order:
            sub = scores[scores["source"] == src]
            ax.scatter(sub[colx], sub[coly], s=11, alpha=0.6, color=src_col[src],
                       label=src)
        ax.axhline(SRD_REF, color="grey", ls="--", lw=0.8)
        ax.axhline(-SRD_REF, color="grey", ls="--", lw=0.8)
        ax.axvline(SRD_REF, color="grey", ls="--", lw=0.8)
        ax.axvline(-SRD_REF, color="grey", ls="--", lw=0.8)
        ax.axhline(0, color="black", lw=0.5)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_xlabel(f"{colx}")
        ax.set_ylabel(f"{coly}")
        ax.set_title(title)
        ax.legend(fontsize=6, ncol=2)
    fig.suptitle("small-segment left vs right flank SRD (coloured by source)", fontsize=11)
    fig.tight_layout()
    out = fig_dir / "small_segment_srd_scatter.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    log(f"saved {out}")


def make_size_dependence(fig_dir, scores):
    """|SRD| and |SRD/√φ| vs segment length (does the statistic inflate with size?)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    L = scores["seg_length"]
    for ax, col, color in zip(axs,
                              [("srd_left", "srd_right"), ("srd_phi_left", "srd_phi_right")],
                              ["#1f77b4", "#ff7f0e"]):
        m = scores[[col[0], col[1]]].abs().max(axis=1)
        ax.scatter(L, m, s=12, alpha=0.5, color=color)
        ax.set_xlabel("segment length (bins)")
        ax.set_ylabel("max |SRD| flank (desc / √φ)")
        lab = col[0].replace("_left", "")
        rho, p = spearmanr(L, m)
        ax.set_title(f"size dependence — {lab} (Spearman ρ={rho:.2f})")
    fig.suptitle("Does SRD inflate/shrink with segment length?", fontsize=11)
    fig.tight_layout()
    out = fig_dir / "small_segment_size_dependence.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    log(f"saved {out}")
    for name in ["desc SRD", "SRD/√φ"]:
        pass
    # report correlations
    for col in [("srd_left", "srd_right"), ("srd_phi_left", "srd_phi_right")]:
        m = scores[[col[0], col[1]]].abs().max(axis=1)
        rho, p = spearmanr(scores["seg_length"], m)
        print(f"[screen] |{col[0].replace('_left','')}| vs seg_length: Spearman rho={rho:.3f} (p={p:.3g})",
              flush=True)


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    if not out_dir.is_absolute():
        out_dir = BETA / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir.parent / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    adata, members, raw, audit = load(args)
    chrom_of = adata.var["seq"].astype(str).to_numpy()
    var = adata.var
    X = adata.X
    sources = [f"cluster{i}" for i in range(5)] + ["all"]

    all_rows = []
    boot_rows = []
    for ch in args.chromosomes:
        if ch not in set(chrom_of):
            log(f"skip {ch}: not present")
            continue
        for src in sources:
            acc_bins = sorted(raw[(raw["source"] == src) & (raw["chromosome"] == ch)]["absolute_bin"].astype(int))
            n_bins = int((chrom_of == ch).sum())
            segs = flat_segments(acc_bins, n_bins)
            if src == "all":
                pb = np.asarray(X.mean(axis=0)).ravel()
            else:
                pb = np.asarray(X[members[src]].mean(axis=0)).ravel()
            cidx = np.where(chrom_of == ch)[0]
            pb_chr = pb[cidx]
            positions = var.iloc[cidx]["start"].to_numpy()
            rows, _ = score_source_chromosome(pb_chr, segs, src, ch, positions,
                                              args.max_small_bins, len(members[src]))
            all_rows.extend(rows)
            if args.n_bootstrap > 0 and rows:
                E = np.array([e - s for s, e in segs])
                boot_rows.extend(bootstrap_source_segments(
                    X, chrom_of, members, src, ch, segs, E, args.n_bootstrap,
                    args.seed, rows))
    scores = pd.DataFrame(all_rows)
    ch_rank = {c: i for i, c in enumerate(args.chromosomes)}
    scores["_ch"] = scores["chromosome"].map(ch_rank)
    scores = scores.sort_values(["_ch", "source", "bin_start"]).drop(columns="_ch").reset_index(drop=True)
    scores.to_csv(out_dir / "small_segment_scores.tsv", sep="\t", index=False)
    log(f"wrote small_segment_scores.tsv ({len(scores)} small segments)")
    if boot_rows:
        boots = pd.DataFrame(boot_rows)
        boots.to_csv(out_dir / "small_segment_bootstrap.tsv", sep="\t", index=False)
        log(f"wrote small_segment_bootstrap.tsv ({len(boots)} rows, n_bootstrap={args.n_bootstrap})")

    for ch in args.chromosomes:
        make_chromosome_figure(fig_dir, args, ch, chrom_of, var, X, members, raw, scores)
    score_viz(fig_dir, scores)
    make_srd_scatter(fig_dir, scores)
    make_size_dependence(fig_dir, scores)
    log(f"proposal counts: {scores['proposal'].value_counts().to_dict()}")
    log("DONE.")


if __name__ == "__main__":
    main()
