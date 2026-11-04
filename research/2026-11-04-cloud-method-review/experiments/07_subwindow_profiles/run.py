"""Experiment 07: per-bin subwindow profiles of the saved windows, with gaps, GC and allele observations.

For each saved window (frozen shared regions) and the 300 random windows of experiment 04
(same rule and seed), per-bin all-cell raw sums are expressed relative to the mean per-bin
raw sum of the region's own flanks. Within-window heterogeneity = SD of log2 per-bin ratios;
an even dosage step predicts a flat plateau (heterogeneity near the random-window level),
while patchy depletion suggests bin-level baseline (accessibility/mappability) structure.
This is a descriptive contrast, not a test that separates the two.
Also reports mean GC and processed-X / raw attenuation per window, and draws per-window
profiles (raw, processed X, AT per raw fragment) with QC-removed grid positions marked.
Deterministic apart from the reused random-window rule (seed tuple (20261104, 4)).
"""
import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

spec = importlib.util.spec_from_file_location("exp04", HERE.parent / "04_fixed_region_contrasts" / "run.py")
E4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(E4)


def main():
    p = C.default_paths()
    d = C.H5AD(p["h5ad"])
    allm = {"all": np.ones(d.n_cells, bool)}
    raw = C.group_means(d, "layers/raw", allm)["all"]
    xm = C.group_means(d, "X", allm)["all"]
    t = pq.read_table(p["allelic"] / "AT.parquet")
    names = [c for c in t.column_names if c not in ("chr", "start", "end", "__index_level_0__")]
    at = np.zeros(d.n_bins)
    for c in names:
        at += t.column(c).to_numpy(zero_copy_only=False)
    del t
    g = pd.read_csv(p["root"] / "pyepi_results" / "binned_genome.csv", index_col=0)
    regs = E4.build_regions(d, p)
    rows = []
    for r in regs:
        if r["family"] not in ("saved", "random"):
            continue
        tb, fb = r["target"], r["flank"]
        ratio = raw[tb] / raw[fb].mean()
        lr = np.log2(np.clip(ratio, 1e-3, None))
        third = max(1, len(tb) // 3)
        rows.append({"family": r["family"], "region": r["name"], "n_bins": len(tb), "rr_flank_raw": raw[tb].mean() / raw[fb].mean(),
                     "rr_flank_x": xm[tb].mean() / xm[fb].mean(), "sd_log2_bin_ratio": float(np.std(lr)),
                     "frac_bins_ratio_lt_0.5": float((ratio < 0.5).mean()), "frac_bins_ratio_gt_0.9": float((ratio > 0.9).mean()),
                     "rr_first_third": raw[tb[:third]].mean() / raw[fb].mean(), "rr_middle": raw[tb[third:-third]].mean() / raw[fb].mean() if len(tb) > 2 else np.nan,
                     "rr_last_third": raw[tb[-third:]].mean() / raw[fb].mean(),
                     "gc_target": float(d.gc[tb].mean()), "gc_flank": float(d.gc[fb].mean()),
                     "x_over_raw_attenuation_log2": float(np.log2((xm[tb].mean() / xm[fb].mean()) / (raw[tb].mean() / raw[fb].mean())))})
    df = pd.DataFrame(rows)
    rnd = df[df.family == "random"]
    # compare heterogeneity at matched depth: random windows with rr < 0.6
    deep = rnd[rnd.rr_flank_raw < 0.6]
    df["sd_pct_vs_random"] = df.sd_log2_bin_ratio.apply(lambda v: float((rnd.sd_log2_bin_ratio <= v).mean() * 100))
    df["sd_pct_vs_deep_random"] = df.sd_log2_bin_ratio.apply(lambda v: float((deep.sd_log2_bin_ratio <= v).mean() * 100))
    df.to_csv(HERE / "subwindow_summary.tsv", sep="\t", index=False)
    gc_slope = np.polyfit(rnd.gc_target - rnd.gc_flank, rnd.x_over_raw_attenuation_log2, 1)
    res = {"experiment": "07_subwindow_profiles", "random_sd_q05_q50_q95": rnd.sd_log2_bin_ratio.quantile([.05, .5, .95]).tolist(),
           "deep_random_n": len(deep), "deep_random_sd_median": float(deep.sd_log2_bin_ratio.median()),
           "random_attenuation_vs_gc_difference_slope_intercept": gc_slope.tolist(),
           "saved": df[df.family == "saved"].to_dict("records"), "environment": C.environment()}
    C.write_json(HERE / "summary.json", res)

    saved = [r for r in regs if r["family"] == "saved"]
    fig, axes = plt.subplots(len(saved), 1, figsize=(13, 2.3 * len(saved)))
    for ax, r in zip(axes, saved):
        b = np.r_[r["flank"][: len(r["flank"]) // 2], r["target"], r["flank"][len(r["flank"]) // 2:]]
        b = np.unique(b)
        ch = d.seq[b[0]]
        x = d.start[b] / 1e6
        fm = raw[r["flank"]].mean(); fx = xm[r["flank"]].mean(); fa = (at[r["flank"]] / raw[r["flank"]].clip(min=1e-9)).mean()
        ax.plot(x, raw[b] / fm, ".-", ms=2, lw=.6, color="#2f5d8a", label="raw / flank")
        ax.plot(x, xm[b] / fx, ".-", ms=2, lw=.6, color="#9b2c2c", label="processed X / flank")
        ax.plot(x, (at[b] / raw[b].clip(min=1e-9)) / fa, ".", ms=2, color="#8a5a00", label="AT per raw / flank")
        ax.axvspan(d.start[r["target"][0]] / 1e6, d.end[r["target"][-1]] / 1e6, color="#ddd", alpha=.5, zorder=0)
        gg = g[(g.chromosome == ch) & (g.start >= d.start[b[0]]) & (g.end <= d.end[b[-1]])]
        rem = gg[~gg.set_index(["chromosome", "start"]).index.isin(list(zip(d.seq[b], d.start[b])))]
        for s in rem.start:
            ax.axvline(s / 1e6, color="k", lw=.4, alpha=.5)
        ax.axhline(0.5, color="#888", ls=":", lw=.8); ax.axhline(1, color="#888", lw=.5)
        ax.set_ylim(0, 2.5); ax.set_title(r["name"] + " (grey = target, black ticks = QC-removed bins)", fontsize=9)
    axes[0].legend(fontsize=7, ncol=3)
    axes[-1].set_xlabel("Mb")
    b64 = C.fig_to_b64(fig); plt.close(fig)
    sv = df[df.family == "saved"]
    frag = [
        ("Question", "<p>Within each saved window, is the depth change an even plateau or patchy? Do processed values attenuate it "
         "relative to raw counts, and does that track GC? Where are QC-removed bins?</p>"),
        ("Saved windows", C.html_table(sv)),
        ("Reference: random 25-bin windows", f"<p>SD of log2 per-bin ratio, 5/50/95th percentiles: "
         f"{', '.join(f'{v:.3f}' for v in res['random_sd_q05_q50_q95'])}; random windows with rr &lt; 0.6 (n = {len(deep)}): median "
         f"{res['deep_random_sd_median']:.3f}. Processed-vs-raw log2 attenuation regressed on target−flank GC difference across "
         f"random windows: slope {gc_slope[0]:.2f}, intercept {gc_slope[1]:.3f}.</p>"),
        ("Profiles", C.img(b64, "window profiles")),
        ("Interpretation", "<ul><li>Per-bin heterogeneity is a descriptive shape cue. Wide windows can contain a narrower event, so a "
         "patchy profile does not exclude a focal change, and an even plateau does not prove dosage.</li>"
         "<li>Processed X attenuates or exaggerates contrasts depending on the GC difference between target and flanks: per-cell GC "
         "LOESS rescales low-GC regions, so GC correction and regional biology are confounded at this scale.</li></ul>"),
    ]
    C.write_html_report(HERE / "report.html", "Subwindow profiles", "Experiment 07 · per-bin shape, GC attenuation and gaps in the saved windows",
                        frag, {"command": "python -I experiments/07_subwindow_profiles/run.py", "summary": "summary.json", "environment": C.environment()})
    pd.set_option("display.width", 250)
    print(sv.to_string()); print(res["random_sd_q05_q50_q95"], len(deep), res["deep_random_sd_median"], gc_slope)


if __name__ == "__main__":
    main()
