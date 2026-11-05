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
    import json as _json
    e2w = pd.DataFrame(_json.loads((HERE.parent / "02_zero_clip_audit" / "summary.json").read_text())["windows"]).set_index("window")
    sv = df[df.family == "saved"].copy()
    sv["removed_bins_in_window"] = sv.region.map(lambda w: int(e2w.loc[w, "removed_non_gap"] + e2w.loc[w, "removed_gap_like"]) if w in e2w.index else np.nan)
    sv["gc_diff"] = sv.gc_target - sv.gc_flank
    sv["x_over_raw_pred_from_gc_slope"] = gc_slope[1] + gc_slope[0] * sv.gc_diff
    chr7 = sv[sv.region.str.startswith("chr7")]
    sdq = res["random_sd_q05_q50_q95"]
    flat = sv[sv.sd_pct_vs_random <= 5]
    verdict = (f"chr7 67.4&ndash;70.0 Mb is an even plateau at {sv.set_index('region').loc['chr7_67.4_70.0','rr_flank_raw']:.2f}x its flanks (per-bin heterogeneity at percentile {sv.set_index('region').loc['chr7_67.4_70.0','sd_pct_vs_random']:.1f} of random windows), "
               f"the other saved windows are uneven or shallow, and processed X attenuates the deepest low-GC dips (log2 attenuation up to {sv.x_over_raw_attenuation_log2.max():.2f}) in line with the GC slope across random windows ({gc_slope[0]:.1f} per unit GC difference).")
    findings = [
        f"Random 25-bin windows: SD of the log2 per-bin ratio is {sdq[0]:.2f}/{sdq[1]:.2f}/{sdq[2]:.2f} (5th/50th/95th percentile); the {len(deep)} random windows with rr &lt; 0.6 have median SD {res['deep_random_sd_median']:.2f}.",
        f"chr7_67.4_70.0: SD {sv.set_index('region').loc['chr7_67.4_70.0','sd_log2_bin_ratio']:.2f}, {100*sv.set_index('region').loc['chr7_67.4_70.0','frac_bins_ratio_lt_0.5']:.0f}% of bins below 0.5x flank mean; thirds {sv.set_index('region').loc['chr7_67.4_70.0','rr_first_third']:.2f}/{sv.set_index('region').loc['chr7_67.4_70.0','rr_middle']:.2f}/{sv.set_index('region').loc['chr7_67.4_70.0','rr_last_third']:.2f} &mdash; the flattest saved window (percentile {sv.set_index('region').loc['chr7_67.4_70.0','sd_pct_vs_deep_random']:.0f} among deep random windows).",
        f"chr7_70.0_73.3: SD {sv.set_index('region').loc['chr7_70.0_73.3','sd_log2_bin_ratio']:.2f} (percentile {sv.set_index('region').loc['chr7_70.0_73.3','sd_pct_vs_random']:.0f} of random windows); thirds {sv.set_index('region').loc['chr7_70.0_73.3','rr_first_third']:.2f}/{sv.set_index('region').loc['chr7_70.0_73.3','rr_middle']:.2f}/{sv.set_index('region').loc['chr7_70.0_73.3','rr_last_third']:.2f} &mdash; a deep middle with shallower shoulders, not a plateau.",
        f"chr7_8.3_10.7 and chr7_51.4_53.5 are fully depleted (all bins below 0.5x) with moderate heterogeneity (SD {sv.set_index('region').loc['chr7_8.3_10.7','sd_log2_bin_ratio']:.2f} and {sv.set_index('region').loc['chr7_51.4_53.5','sd_log2_bin_ratio']:.2f}; percentiles {sv.set_index('region').loc['chr7_8.3_10.7','sd_pct_vs_random']:.0f} and {sv.set_index('region').loc['chr7_51.4_53.5','sd_pct_vs_random']:.0f} of random windows).",
        f"chr9 windows are shallower and uneven: rr {sv.set_index('region').loc['chr9_20_27','rr_flank_raw']:.2f} (SD {sv.set_index('region').loc['chr9_20_27','sd_log2_bin_ratio']:.2f}) and {sv.set_index('region').loc['chr9_22.5_26.6','rr_flank_raw']:.2f} (SD {sv.set_index('region').loc['chr9_22.5_26.6','sd_log2_bin_ratio']:.2f}); only {100*sv.set_index('region').loc['chr9_20_27','frac_bins_ratio_lt_0.5']:.0f}% and {100*sv.set_index('region').loc['chr9_22.5_26.6','frac_bins_ratio_lt_0.5']:.0f}% of bins fall below 0.5x; {int(sv.set_index('region').loc['chr9_20_27','removed_bins_in_window'])} and {int(sv.set_index('region').loc['chr9_22.5_26.6','removed_bins_in_window'])} grid bins in them were removed by QC (experiment 02).",
        f"Processed X attenuates the contrast: rr in X vs raw = {', '.join(f'{r.rr_flank_x:.2f} vs {r.rr_flank_raw:.2f}' for r in chr7.itertuples())} for the four chr7 windows (log2 attenuation {chr7.x_over_raw_attenuation_log2.min():.2f}&ndash;{chr7.x_over_raw_attenuation_log2.max():.2f}); their target GC is {chr7.gc_diff.min():+.3f} to {chr7.gc_diff.max():+.3f} relative to flanks.",
        f"Across the {len(rnd)} random windows, log2 attenuation regresses on target&minus;flank GC with slope {gc_slope[0]:.1f} and intercept {gc_slope[1]:.3f}; the regression-predicted attenuation for the chr7 windows is {chr7.x_over_raw_pred_from_gc_slope.min():.2f}&ndash;{chr7.x_over_raw_pred_from_gc_slope.max():.2f} against observed {chr7.x_over_raw_attenuation_log2.min():.2f}&ndash;{chr7.x_over_raw_attenuation_log2.max():.2f}.",
    ]

    def note(r):
        s = f"{int(r.n_bins)} bins, raw rate {r.rr_flank_raw:.2f}x flanks; "
        s += f"thirds {r.rr_first_third:.2f}/{r.rr_middle:.2f}/{r.rr_last_third:.2f}; SD {r.sd_log2_bin_ratio:.2f} (percentile {r.sd_pct_vs_random:.0f} of random windows). "
        if r.family == "saved" and r.sd_pct_vs_random <= 5 and r.frac_bins_ratio_lt_0_5 == 1.0:
            s += "Even plateau: every bin below 0.5x with little spread. "
        elif r.frac_bins_ratio_lt_0_5 == 1.0:
            s += "Every bin is below 0.5x, with moderate spread. "
        elif r.frac_bins_ratio_lt_0_5 >= 0.5:
            s += f"Uneven: {100*r.frac_bins_ratio_lt_0_5:.0f}% of bins below 0.5x and {100*r.frac_bins_ratio_gt_0_9:.0f}% above 0.9x. "
        elif r.rr_flank_raw > 1.1:
            s += "Above the flanks (comparison window, wide). "
        else:
            s += f"Shallow and uneven: {100*r.frac_bins_ratio_lt_0_5:.0f}% of bins below 0.5x and {100*r.frac_bins_ratio_gt_0_9:.0f}% above 0.9x. "
        s += f"GC target&minus;flank {r.gc_diff:+.3f}; log2 X/raw attenuation {r.x_over_raw_attenuation_log2:+.2f} (GC-slope prediction {r.x_over_raw_pred_from_gc_slope:+.2f}); {int(r.removed_bins_in_window)} QC-removed grid bins in the window."
        return s

    sv = sv.rename(columns={"frac_bins_ratio_lt_0.5": "frac_bins_ratio_lt_0_5", "frac_bins_ratio_gt_0.9": "frac_bins_ratio_gt_0_9"})
    reg_tab = sv[["region"]].copy()
    reg_tab["note"] = [__import__("html").unescape(note(r)) for r in sv.itertuples()]
    reg_html = C.html_table(reg_tab).replace("<td>chr", "<td>chr").replace("<td>", "<td style='text-align:left;white-space:normal'>")
    shape_cols = sv[["region", "n_bins", "rr_flank_raw", "rr_first_third", "rr_middle", "rr_last_third"]]
    het_cols = sv[["region", "sd_log2_bin_ratio", "sd_pct_vs_random", "sd_pct_vs_deep_random", "frac_bins_ratio_lt_0_5", "frac_bins_ratio_gt_0_9"]]
    gc_cols = sv[["region", "gc_target", "gc_flank", "gc_diff", "rr_flank_raw", "rr_flank_x", "x_over_raw_attenuation_log2", "x_over_raw_pred_from_gc_slope"]]
    sections = [
        ("Question and why it matters",
         "<p><b>Question.</b> Within each saved window, is the depth change an even plateau or patchy? Do processed values attenuate it relative to raw counts, and does that track GC? Where are QC-removed bins?</p>"
         "<p><b>Why it matters.</b> The recommendation reads the chr7 dips as most consistent with a shared regional baseline possibly combined with a smaller dosage change. "
         "An even plateau is what a dosage step would predict, while patchy depletion points to bin-level accessibility or mappability structure; the GC attenuation tests the claim that GC correction and regional biology are confounded at this scale. "
         "This is a descriptive contrast, not a test that separates the two.</p>"),
        ("Inputs",
         f"<ul><li><b>Depth.</b> All-cell mean per-bin raw counts (<code>layers/raw</code>) and processed <code>X</code> from <code>count_matrix.h5ad</code> ({d.n_cells:,} cells &times; {d.n_bins:,} bins); per-bin GC from <code>var.GC</code>.</li>"
         "<li><b>Allele observations.</b> Per-bin pooled AT from <code>AT.parquet</code> divided by raw counts, drawn in the figure (site-observation count, not unique fragments).</li>"
         "<li><b>QC-removed grid positions.</b> From <code>pyepi_results/binned_genome.csv</code> versus the retained grid (counts per window from experiment 02).</li>"
         f"<li><b>Regions.</b> The {len(sv)} saved windows (6 windows + chr9 70&ndash;100 Mb comparison) with frozen 50-bin flanks on each side, and the {len(rnd)} random 25-bin windows of experiment 04 (same rule and seed).</li>"
         "<li><b>Not available.</b> Per-fragment positions (no sub-100-kb resolution), mappability, and pre-QC counts for removed bins.</li></ul>"),
        ("Method",
         "<ol><li>Express each target bin&rsquo;s all-cell raw mean relative to the mean of the region&rsquo;s flank bins (ratio per bin).</li>"
         "<li>Heterogeneity = SD of log2 per-bin ratios (ratios clipped at 1e-3 before log); also the share of bins below 0.5 and above 0.9, and the ratio in the first, middle and last third of the window.</li>"
         "<li>Place each saved window&rsquo;s SD among the 300 random windows, and among the random windows with rr &lt; 0.6 (matched depth).</li>"
         "<li>Processed-vs-raw attenuation: log2 of (rr in X / rr in raw). Regress it on the target&minus;flank GC difference across the random windows and compare the saved windows with the regression line.</li></ol>"
         "<pre>ratio_j        = raw_j / mean(raw over flank bins)             # per target bin j\nsd_log2        = SD_j log2(max(ratio_j, 1e-3))\nthirds         = mean(raw) over first / middle / last len//3 bins / mean(flank)\natten_log2     = log2( (mean X_target / mean X_flank) / (mean raw_target / mean raw_flank) )\natten ~ slope * (GC_target - GC_flank) + intercept      # least squares over random windows</pre>"
         f"<p>Deterministic apart from the reused random-window rule (seed tuple ({C.SEED}, 4), numpy PCG64, shared with experiment 04). Percentiles are of the random-window distribution.</p>"),
        ("Results", "", 2),
        ("Reference: random 25-bin windows",
         C.html_table(pd.DataFrame([{"quantity": "SD of log2 per-bin ratio, 5th percentile", "value": sdq[0]}, {"quantity": "50th percentile", "value": sdq[1]}, {"quantity": "95th percentile", "value": sdq[2]},
                                    {"quantity": f"random windows with rr &lt; 0.6 (n)", "value": len(deep)}, {"quantity": "their median SD", "value": res["deep_random_sd_median"]},
                                    {"quantity": "attenuation vs GC-difference slope", "value": gc_slope[0]}, {"quantity": "intercept", "value": gc_slope[1]}])) +
         C.caption("Reference values from the 300 random windows. An even plateau should sit at the low end of the SD distribution (near or below the 5th percentile); patchy depletion at or above the median. Deep random windows are slightly less heterogeneous than random windows overall."), 3),
        ("Window shape: thirds", C.html_table(shape_cols) + C.caption("Raw rate relative to the flank mean in the whole window and in its first, middle and last third. A plateau has equal thirds; a shoulder shows larger values in the outer thirds. Windows differ in width, so thirds are not directly comparable between windows."), 3),
        ("Per-bin heterogeneity", C.html_table(het_cols) + C.caption("<code>sd_log2_bin_ratio</code> and its percentile among all random windows (<code>sd_pct_vs_random</code>) and among deep random windows (<code>sd_pct_vs_deep_random</code>); shares of bins below 0.5x and above 0.9x of the flank mean. Low percentiles with all bins below 0.5x indicate an even plateau. Note that window widths (19&ndash;295 bins) differ from the 25-bin random windows, which affects SD comparisons for the wider windows."), 3),
        ("GC and processed-vs-raw attenuation", C.html_table(gc_cols) + C.caption("Mean GC of target and flanks, rr in raw and processed X, the observed log2 attenuation (positive = X shallower than raw) and the attenuation predicted from the random-window regression on the GC difference. Look for observed values near the prediction, and for large positive values in low-GC windows."), 3),
        ("Profiles", C.img(b64, "window profiles") + C.caption("One panel per saved window. Blue: raw per-bin all-cell mean over the flank mean; red: processed X over its flank mean; dots: AT per raw fragment over its flank mean. The grey band is the target, vertical black ticks mark QC-removed grid bins, dotted line 0.5 and solid line 1. Look for a flat blue band under 0.5 (plateau), for red lying above blue (attenuation), and for gaps of black ticks."), 3),
        ("Per-region detail", reg_html + C.caption("One row per saved window; the note restates the table values. Windows differ in width, so percentiles for the wide windows are indicative only."), 2),
        ("Sensitivity",
         "<p>Heterogeneity is reported against all random windows and against matched-depth random windows (rr &lt; 0.6); the attenuation prediction uses a linear regression over 300 windows. Not varied: the third-based split, the log2 clipping at 1e-3, the 50-bin flank width and the 25-bin random width.</p>"),
        ("Limits and alternative explanations",
         "<ul><li>Per-bin heterogeneity is a descriptive shape cue. Wide windows can contain a narrower event, so a patchy profile does not exclude a focal change, and an even plateau does not prove dosage.</li>"
         "<li>Bins are 100 kb pooled over all cells: sub-bin structure and cell-level mixture are not visible here.</li>"
         "<li>Processed X attenuates or exaggerates contrasts depending on the GC difference between target and flanks: per-cell GC LOESS rescales low-GC regions, so GC correction and regional biology are confounded at this scale; the regression shows association, not mechanism.</li>"
         "<li>SD depends on the number of bins and on depth; saved windows (19&ndash;295 bins) are compared with 25-bin random windows.</li>"
         "<li>QC-removed bins are positions, not measurements; their content is unknown.</li></ul>"),
    ]
    glossary = {
        "n_bins": "Retained 100 kb bins in the target.",
        "rr_flank_raw / rr_flank_x": "Mean raw / processed-X value in the target divided by the same in the flanks.",
        "rr_first_third / rr_middle / rr_last_third": "Raw rate relative to the flank mean in the first, middle and last third of the window.",
        "sd_log2_bin_ratio": "SD of log2 of per-bin raw ratio to the flank mean; low = even, high = patchy.",
        "sd_pct_vs_random / sd_pct_vs_deep_random": "Percent of all random windows / random windows with rr &lt; 0.6 with SD at or below this window&rsquo;s.",
        "frac_bins_ratio_lt_0_5 / frac_bins_ratio_gt_0_9": "Share of target bins whose raw ratio to the flank mean is below 0.5 / above 0.9.",
        "gc_target / gc_flank / gc_diff": "Mean bin GC in the target / flanks, and their difference.",
        "x_over_raw_attenuation_log2": "log2 of (rr in X) / (rr in raw); positive means processed X is shallower than raw (for losses).",
        "x_over_raw_pred_from_gc_slope": "Attenuation predicted by the random-window regression on gc_diff.",
        "removed_bins_in_window": "Grid bins in the window missing from the retained grid (experiment 02).",
        "AT per raw": "Pooled selected-site allele observations per raw fragment in a bin.",
    }
    C.write_html_report(HERE / "report.html", "Subwindow profiles", "Experiment 07 &middot; per-bin shape, GC attenuation and gaps in the saved windows",
                        sections, {"command": "python -I experiments/07_subwindow_profiles/run.py", "summary": "summary.json", "seed": [C.SEED, 4], "environment": C.environment()},
                        key_findings=findings, verdict=verdict, glossary=glossary,
                        supports=["chr7 67.4&ndash;70.0 Mb is an even plateau, consistent with a uniform change over the window (dosage or a regional baseline).",
                                  "The other saved windows are less even; chr7_70.0_73.3 and the chr9 windows have wide bin-to-bin spread.",
                                  "GC correction attenuates low-GC dips, with attenuation tracking the target&minus;flank GC difference across random windows."],
                        not_supports=["That an even plateau is dosage rather than a regional accessibility or mappability baseline.",
                                      "That patchy windows lack a focal change.",
                                      "Any statement about the contents of QC-removed bins."],
                        decision_impact="Supports the recommendation that raw-count contrast, not processed X, should carry the depth evidence, and that GC correction can hide real or baseline structure; "
                                        "the plateau at chr7 67.4&ndash;70.0 Mb makes a cross-sample baseline and mappability track the decisive missing experiment (rank 3).",
                        next_steps=["Compare the same windows in unrelated samples with a mappability track, to see whether the plateau recurs.",
                                    "Repeat the heterogeneity comparison with width-matched random windows for each saved window.",
                                    "Test whether the 19q13-qter raw excess lost by GC correction has the same shape in an independent sample."])
    pd.set_option("display.width", 250)
    print(sv.to_string()); print(res["random_sd_q05_q50_q95"], len(deep), res["deep_random_sd_median"], gc_slope)


if __name__ == "__main__":
    main()
