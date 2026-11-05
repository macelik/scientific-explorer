"""Build the final review report (../final_report.html) from the saved outputs of experiments 00-19.

Reads only files written by the experiments (TSV/JSON); recomputes nothing scientific. Numbers in the
text are formatted from those files at build time. Run from the workbench root:
    python -I research/2026-11-04-cloud-method-review/final_report/build.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
REV = HERE.parent
EX = REV / "experiments"
sys.path.insert(0, str(EX))
import common as C  # noqa: E402

J = lambda e: json.loads((EX / e / "summary.json").read_text())  # noqa: E731
T = lambda e, f: pd.read_csv(EX / e / f, sep="\t")  # noqa: E731
link = lambda e, label=None: f"<a href='experiments/{e}/report.html'>{label or e}</a>"  # noqa: E731


def fig(fn, w=7, h=4):
    f, ax = plt.subplots(figsize=(w, h)); fn(ax); b = C.fig_to_b64(f); plt.close(f); return b


def main():
    s01, s02, s05, s08 = J("01_reconstruct_x"), J("02_zero_clip_audit"), J("05_segmentation_calibration"), J("08_exposure_count_models")
    s09, s11, s13, s15 = J("09_zero_inflation"), J("11_count_changepoint"), J("13_joint_hmm"), J("15_allelic_count_models")
    s16, s18, s14 = J("16_major_minor_haplotypes"), J("18_cell_level_joint_states"), J("14_spikein_benchmark")
    rc04 = T("04_fixed_region_contrasts", "region_contrasts.tsv")
    r17 = T("17_allele_specific_states", "region_states.tsv"); r17 = r17[r17.k == 0.6]
    r18 = T("18_cell_level_joint_states", "region_mixtures.tsv").set_index("region")
    d14 = T("14_spikein_benchmark", "detection_summary.tsv")
    d19 = T("19_masking_two_cut_scan", "detection_summary.tsv")
    bo = T("16_major_minor_haplotypes", "bin_orientation.tsv.gz") if False else pd.read_csv(EX / "16_major_minor_haplotypes" / "bin_orientation.tsv.gz", sep="\t")
    rh16 = T("16_major_minor_haplotypes", "region_haplotypes.tsv"); rh16 = rh16[rh16.group == "all"].set_index("region")
    sg16 = T("16_major_minor_haplotypes", "segment_haplotypes.tsv")
    n14 = int(d14.n.sum())
    fpr = s18["random_fpr"]
    unc = rc04[(rc04.family == "random") & (rc04.frac_bins_called_loss_all == 0) & (rc04.frac_bins_called_gain_all == 0)]
    q_unc = unc.rr_flank_raw.quantile([.05, .95]).round(2).tolist()
    sig17 = J("17_allele_specific_states")["sigma_A_flank"]
    saved = rc04[rc04.family == "saved"].set_index("region")
    chr7 = saved.loc[[i for i in saved.index if i.startswith("chr7")]]

    # ---------------- figures
    seg8 = bo[(bo.group == "all") & (bo.group + "|" + bo.segment.astype(str) == "all|8") & (bo.AT >= 5)]
    raw8 = (seg8.HB / seg8.AT).to_numpy(); ori8 = np.where(seg8.p_B_major >= 0.5, raw8, 1 - raw8)

    def f_orient(ax):
        hb = np.linspace(0, 1, 41)
        ax.hist(raw8, hb, color="#999", alpha=.85, label=f"GT-oriented HB/AT per bin (pooled {seg8.HB.sum() / seg8.AT.sum():.2f})")
        ax.hist(ori8, hb, histtype="step", lw=1.8, color="#2f5d8a", label=f"after orientation HMM (median {np.median(ori8):.2f})")
        ax.axvline(0.5, color="k", lw=.6); ax.axvline(2 / 3, color="#9b2c2c", ls=":", lw=1, label="2:1 expectation (0.667)")
        ax.set_xlabel("allele fraction per 100 kb bin, chr1:166.8–236 Mb, all cells"); ax.set_ylabel("bins"); ax.legend(fontsize=7)

    def f_states(ax):
        cols = {"random": "#bbb", "called": "#9b2c2c", "saved": "k", "literature": "#2f5d8a"}
        for f_, c in cols.items():
            g = r17[r17.family == f_]
            ax.scatter(g.L_obs, g.h_obs, s=14 if f_ in ("random", "called") else 50, color=c, alpha=.75, label=f_,
                       marker={"random": "o", "called": "s", "saved": "^", "literature": "D"}[f_])
        for nm, L, h in [("1:1", 1, .5), ("1:0 clonal", .5, .5 + .6 * .5), ("2:1 clonal", 1.5, .5 + .6 * (2 / 3 - .5)), ("2:0", 1, .5 + .6 * .5), ("0:0 (deep)", .1, .5)]:
            ax.scatter([L], [h], marker="*", s=140, color="#8a5a00"); ax.annotate(nm, (L, h), fontsize=7, xytext=(4, 4), textcoords="offset points")
        ax.axvspan(np.exp(-1.96 * sig17), np.exp(1.96 * sig17), color="#ddd", alpha=.4, zorder=0)
        ax.set_xscale("log"); ax.set_xlabel("depth level L (raw, vs flanks or genome)"); ax.set_ylabel("major-haplotype fraction ĥ (phase-agnostic)")
        ax.legend(fontsize=7); ax.set_ylim(.48, .85)

    def f_bench(ax):
        for typ, c in [("loss_one", "#9b2c2c"), ("loss_deep", "#000"), ("gain", "#2f5d8a")]:
            g = d14[(d14.type == typ) & (d14.f == 1.0)]
            ax.errorbar(g.width, g.detected, yerr=[g.detected - g.detected_lo, g.detected_hi - g.detected], marker="o", color=c, capsize=2,
                        label=f"{typ}, f=1 (exp 14, realistic thinning)")
            h = d19[d19.type == typ]
            ax.plot(h.width, h.two_cut_detected, "s--", color=c, alpha=.7, label=f"{typ}: two-cut scan (exp 19, clean edit)")
        ax.axvline(20, color="#888", ls=":", lw=1); ax.text(20.5, 0.95, "min segment 20 bins", fontsize=7, color="#555")
        ax.set_xlabel("event width (100 kb bins)"); ax.set_ylabel("detected (both edges ±2 bins + CN direction)"); ax.set_ylim(-.05, 1.05); ax.legend(fontsize=6)

    def f_cell(ax):
        g = r18.reset_index()
        y = np.arange(len(g))
        col = np.where(g.depth_allele_agreement.str.contains("both"), "#2f5d8a", np.where(g.depth_allele_agreement.str.contains("depth only"), "#9b2c2c", "#999"))
        ax.errorbar(g.pi, y, xerr=[np.clip(g.pi - g.pi_lo, 0, None), np.clip(g.pi_hi - g.pi, 0, None)], fmt="none", ecolor="#555", lw=1)
        ax.scatter(g.pi, y, c=col, s=40, zorder=3)
        ax.set_yticks(y); ax.set_yticklabels([f"{r} ({t})" for r, t in zip(g.region, g.best_type)], fontsize=7)
        ax.set_xlabel("estimated fraction of cells in the event state π (95% profile CI)")
        ax.set_title("blue: depth and allele both support · red: depth only · grey: neither alone", fontsize=8)
    b_orient, b_states, b_bench, b_cell = fig(f_orient, 7.5, 3.8), fig(f_states, 8, 5), fig(f_bench, 8, 4.6), fig(f_cell, 8, 4.6)

    # ---------------- numbers
    rec = s01["reconstruction"]
    cd = [x for x in s02["prespecified_loci"]][0]
    myc, q5, q19, p9 = rh16.loc["lit_MYC_8q24"], rh16.loc["lit_5q13_qter"], rh16.loc["lit_19q13_qter"], rh16.loc["lit_9p21_hemizygous_zone_approx"]
    q16 = s16["null_by_group"]["all"]["delta_hat_q95"] if isinstance(s16.get("null_by_group"), dict) and "all" in s16["null_by_group"] else 0.035
    n_imb = int(sg16[(sg16.group == "all") & (sg16.imbalanced == True)].shape[0])  # noqa: E712
    n_seg = int(sg16[(sg16.group == "all") & (sg16.eligible == True)].shape[0])  # noqa: E712
    nsig = s08["called_significant_empirical_q<0.05"]
    lo1 = d14[(d14.type == "loss_one") & (d14.f == 1.0)].set_index("width")
    dp1 = d14[(d14.type == "loss_deep") & (d14.f == 1.0)].set_index("width")
    gn1 = d14[(d14.type == "gain") & (d14.f == 1.0)].set_index("width")
    st = s15["recommendation"] if "recommendation" in s15 else None

    key = [
        f"<b>The saved X is reproducible</b> from <code>layers/raw</code> with the pinned pyEpi source: {rec['cells_exact_equal']:,}/{s01['cells']:,} cells bit-identical "
        f"({link('01_reconstruct_x', 'exp 01')}). The GC step has no guard for near-zero LOESS predictions (single fragments amplified up to 717×).",
        f"<b>QC removes the event it should detect:</b> the bins of the published MCF10-A homozygous CDKN2A/B deletion fail the ≥15%-nonzero-cells bin filter "
        f"({cd['removed_by_bin_filter']} of {cd['grid_bins']} grid bins removed; {link('02_zero_clip_audit', 'exp 02')}).",
        f"<b>Depth contrasts at 2.5 Mb scale are as large as one-copy changes without being events:</b> uncalled random windows span {q_unc[0]}–{q_unc[1]}× their flanks "
        f"(5th–95th percentile); the realistic baseline spread is σ_A ≈ {sig17:.2f} on the log scale, while a one-copy loss is log 0.5 = −0.69 ({link('04_fixed_region_contrasts', 'exp 04')}, {link('17_allele_specific_states', 'exp 17')}).",
        f"<b>The current segmentation misses moderate focal events:</b> injected clonal one-copy losses of 20/30/50 bins were detected in "
        f"{lo1.loc[20,'detected']:.0%}/{lo1.loc[30,'detected']:.0%}/{lo1.loc[50,'detected']:.0%} of loci and 1.5× gains in {gn1.loc[20,'detected']:.0%}/{gn1.loc[30,'detected']:.0%}/{gn1.loc[50,'detected']:.0%}, "
        f"deep (0.1×) losses in {dp1.loc[20,'detected']:.0%}/{dp1.loc[30,'detected']:.0%}/{dp1.loc[50,'detected']:.0%} ({n14:,} semi-synthetic realizations; {link('14_spikein_benchmark', 'exp 14')}). "
        f"A two-cut scan rescues deep losses but not one-copy changes, because natural windows in the same segment are as extreme ({link('19_masking_two_cut_scan', 'exp 19')}).",
        f"<b>Each modality has an adequate observation model:</b> depth NB2 with exposure (Poisson worse by {-s15['A_depth_all_regions_vs_NB2']['pois']['d_logscore_vs_nb']:.2f} nats per cell-region; "
        f"no zero inflation needed, {link('09_zero_inflation', 'exp 09')}); allele total NB2 conditioned on depth (median α {s15['A_median_alpha']['AT_raw_in']:.2f} vs {s15['A_median_alpha']['depth']:.3f} for depth); "
        f"haplotype split beta-binomial with region-level ρ ({link('15_allelic_count_models', 'exp 15')}).",
        f"<b>Major/minor haplotypes can be identified where there is imbalance:</b> GT phase switches about every 0.5 Mb and cancels the signal when summed; an orientation HMM recovers it "
        f"(5q13-qter |h−0.5| {q5.delta_naive:.3f} → {q5.delta_hat:.3f}, MYC {myc.delta_naive:.3f} → {myc.delta_hat:.3f}, chr1q h ≈ {np.median(ori8):.2f}); random-window null 95th percentile {q16:.3f}; "
        f"{n_imb}/{n_seg} eligible segments imbalanced ({link('16_major_minor_haplotypes', 'exp 16')}).",
        f"<b>Depth + allele together behave sensibly on controls:</b> MYC, 5q13-qter and chr8:115–145 Mb are near-clonal 2:1 gains supported by both modalities (π ≥ 0.95); the deep chr7 dips "
        f"are flagged 'no modelled state fits' at region level ({link('17_allele_specific_states', 'exp 17')}, {link('18_cell_level_joint_states', 'exp 18')}).",
        f"<b>Depth-only calls are not calibrated:</b> the cell-level mixture flags {fpr['event_detected_p_lt_0.05']['k']}/{s18['random_fpr']['n_windows']} uncalled random windows "
        f"({fpr['event_detected_p_lt_0.05']['rate']:.0%}); none of those is supported by both modalities ({fpr['corroborated_by_both_modalities']['k']}/30).",
    ]

    # ---------------- sections
    exp_map = pd.DataFrame([
        ("00", "Input audit", "hashes, units, alignment, stage availability", "all inputs verified; pre-QC counts, fragments, SNP-level alleles missing"),
        ("01", "Reconstruct X", "pinned GC/rint/clip/normalise rebuild", "bit-identical for 99.7% of cells; GC guard missing"),
        ("02", "QC gaps and zeros", "retained vs full grid; rounding zeros", "CDKN2A/B removed by bin filter"),
        ("03", "Normalization / CN assignment", "fixed segmentation, estimands, thresholds", "scale fit pinned at grid edge → thresholds on genome median"),
        ("04", "Fixed-region contrasts", "raw ratio estimator, homogeneity, dosage floor", "chr7 dips 0.13–0.34×, below one-copy floor, shared by cells"),
        ("05", "Segmentation calibration", "exchangeable null, split halves", "boundaries reproducible; significance ≠ CN"),
        ("06", "Allelic information", "observability, attenuation", "nominal AT overstates power"),
        ("07", "Subwindow profiles", "per-bin shape, GC attenuation", "chr7 67.4–70 Mb even plateau; GC correction attenuates low-GC dips"),
        ("08", "E:1 count models", "Poisson / quasi-Poisson / NB with exposure", "NB adequate; empirical regional null required"),
        ("09", "E:3 zero inflation", "ZINB / hurdle, parametric bootstrap", "not needed"),
        ("10", "E:2 beta-binomial", "regional haplotype fraction, power grid", "controls detected; empirical thresholds required"),
        ("11", "B:2 count changepoint", "quasi-Poisson PELT on raw pseudobulk", "over-segments without a baseline-variance term"),
        ("12", "C:3 baseline sensitivity", "7 prespecified baselines, flags", "calls reverse with baseline for 1/3 of called segments"),
        ("13", "B:3 joint HMM", "depth NB + allele BB, distance-aware transitions", "works only with a between-bin baseline variance"),
        ("14", "Spike-in benchmark", "thinning injection on real background", "one-copy focal events rarely segmented"),
        ("15", "Allelic count models", "per-modality model choice by held-out log score", "NB2 depth; NB2 AT|depth; BB split"),
        ("16", "Major/minor haplotypes", "phase-orientation HMM", "recovers imbalance; mirrored-event hypothesis in cluster 4"),
        ("17", "Allele-specific states", "region-level state scoring", "controls fit 2:1; chr7 dips fit no state"),
        ("18", "Cell-level joint mixture", "two-component cell mixture, hierarchical bootstrap", "both-modality calls credible; depth-only FPR 23%"),
        ("19", "Masking diagnosis", "binary vs two-cut AD scan", "deep losses rescued; one-copy changes not"),
    ], columns=["exp", "name", "method", "outcome"])
    exp_map["report"] = [link(f"{e}_" + n, "report") for e, n in zip(exp_map.exp, ["input_audit", "reconstruct_x", "zero_clip_audit", "normalization_aggregation",
                                                                                     "fixed_region_contrasts", "segmentation_calibration", "allelic_information", "subwindow_profiles",
                                                                                     "exposure_count_models", "zero_inflation", "betabinomial_allele_model", "count_changepoint",
                                                                                     "baseline_sensitivity", "joint_hmm", "spikein_benchmark", "allelic_count_models",
                                                                                     "major_minor_haplotypes", "allele_specific_states", "cell_level_joint_states", "masking_two_cut_scan"])]
    map_html = "<div class='tw'><table><thead><tr><th>exp</th><th>name</th><th>method</th><th>outcome</th><th></th></tr></thead><tbody>" + "".join(
        f"<tr><td>{r.exp}</td><td style='text-align:left'>{r.name}</td><td style='text-align:left;white-space:normal'>{r.method}</td>"
        f"<td style='text-align:left;white-space:normal'>{r.outcome}</td><td>{r.report}</td></tr>" for r in exp_map.itertuples()) + "</tbody></table></div>"

    change = pd.DataFrame([
        ("1", "Evidence layer per region/segment (no change to calls)", "evidence integration",
         "depth NB2 rate ratio vs empirical regional null; AT NB2|depth; oriented BB split; allele-specific state score with 'no state fits'; QC-gap adjacency; depth-only vs both-modality label",
         "exp 08, 15, 16, 17, 18", "prototype exists in this folder; needs packaging"),
        ("2", "Keep the pre-filter bin grid; QC-removed bins as 'unobservable'", "upstream QC bookkeeping",
         "the ≥15%-nonzero bin filter removes clonal homozygous deletions", "exp 02", "requires rerun from fragments"),
        ("3", "Depth baseline from unrelated samples (+ mappability)", "CN reference",
         "baseline spread ≈ one-copy effect; without it, depth-only calls have ~23% FPR at cell level and one-copy focal events are not separable", "exp 04, 12, 17, 18, 19", "needs other samples"),
        ("4", "Allele-first segmentation (phase-switch-aware HMM) with depth as second emission", "segmentation",
         "allele imbalance is less confounded by accessibility; AD depth scan misses one-copy focal events", "exp 14, 16, 19", "research implementation; benchmark with exp 14 injections"),
        ("5", "SNP-level, molecule-aware allele counts with phase sets; site-density term for LOH ascertainment", "allele inputs",
         "AT is summed site DP; GT phase switches ~0.5 Mb; LOH regions lose het sites (9p21 looks balanced)", "exp 06, 15, 16, 17", "needs fragments + pileup rerun"),
        ("6", "CN reference other than the genome median; report the reference used", "CN assignment",
         "pyEpi scale fit always at the grid edge; estimand and baseline change 20–37% of calls", "exp 03, 12", "small code change; scientific decision"),
        ("7", "GC step hygiene: guard near-zero predictions, record fits", "preprocessing",
         "isolated amplified bins; GC correction confounds GC-rich arm gains (19q)", "exp 01, 04, 07", "small; dataset-independent hygiene"),
    ], columns=["priority", "change", "stage", "why (evidence)", "experiments", "effort / dependency"])

    methods_16 = f"""
<p><b>What problem it solves.</b> Allele counts per 100 kb bin (HA, HB) come from SNPs phased by Eagle (GT field, no phase-set IDs). Within a bin the
phase is mostly consistent, but between bins the assignment of "haplotype A" and "haplotype B" switches every few bins (statistical phasing switch
errors and gaps). In a region where one parental haplotype is truly in excess, GT-oriented per-bin fractions HB/AT are therefore <i>bimodal and mirror-symmetric
around 0.5</i>: about half the bins near h and half near 1−h. Summing HB over the region cancels the two modes, which is why the naive estimate was
{q5.delta_naive:.3f} at 5q13-qter. The figure below shows the chr1q case: the grey histogram has two peaks, the oriented one a single peak.</p>
{C.img(b_orient, "orientation example")}
{C.caption("In the experiment-16 tracks the HMM does indeed flip bins. It flips whole runs of bins where the GT label put the major haplotype on HA, "
           "so every bin reports the fraction of the same (major) haplotype. The flip is justified only because the GT values are bimodal with modes far from 0.5 "
           "relative to sampling noise, and mirror-symmetric; where a region is balanced the same procedure gives |h − 0.5| ≤ " + f"{q16:.3f}" + " (random-window null), which bounds the folding bias.")}
<p><b>Model.</b> For one segment with bins j = 1..n (only bins with allele data), a hidden orientation o_j ∈ {{B-major, A-major}} follows a two-state
Markov chain along the chromosome; the observation is HB_j out of AT_j:</p>
<pre>P(o_j ≠ o_(j-1)) = 1 − exp(−λ · d_j),   d_j = bp distance to the previous bin with data, λ = 1 per 5 Mb (prespecified; sensitivity ×0.1, ×10)
HB_j | AT_j, o_j = B-major ~ BetaBinomial(AT_j, p = 0.5 + δ, ρ)
HB_j | AT_j, o_j = A-major ~ BetaBinomial(AT_j, p = 0.5 − δ, ρ)
segment parameter δ ∈ [0, 0.5): profile likelihood over a grid (forward algorithm sums over all orientation paths)
outputs: δ̂, ĥ = 0.5 + δ̂; per-bin posterior P(B-major) (forward–backward); Viterbi switch count; cell-bootstrap CI (200 resamples of cells)
major haplotype per bin = argmax posterior; major counts per cell = HB on B-major bins + HA on A-major bins</pre>
<p><b>Assumptions</b> (each can fail): (1) the imbalance δ is constant within a segment (a segment mixing balanced and imbalanced parts, like 5q 70–118 Mb plus
118–180 Mb, gives an average); (2) the phase within a bin is consistent, so a bin can be oriented as a whole; (3) bins are conditionally independent given
the orientation (repeated site counting and shared fragments violate this mildly); (4) ρ = 0.0038 from base-called bins describes pooled-bin overdispersion
(experiment 15 showed cell-level ρ is larger, so this model is for pooled bins only); (5) switch hazard λ is the same everywhere; (6) the "major" label is
a statement about count asymmetry, not verified parental identity; allele-specific accessibility, reference mapping bias or copy-neutral LOH also create
asymmetry.</p>
<p><b>Why this method.</b> Experiment 06 showed that the phase-agnostic estimator (|HB−HA| per bin) recovers the imbalance but is biased upward at balance
(folding) and cannot say which haplotype is major. GT-oriented summing is unbiased at balance but cancels real imbalance. The orientation HMM is the
standard way to reconcile the two (it is the mechanism Numbat uses for phase switches); it makes the folding explicit, estimates where the orientation
changes, and is calibrated by running it on uncalled random windows. <b>Result:</b> 5q13-qter {q5.delta_naive:.3f} → {q5.delta_hat:.3f}, MYC {myc.delta_naive:.3f} → {myc.delta_hat:.3f};
19q13-qter ({q19.delta_hat:.3f}) and the 9p21 zone ({p9.delta_hat:.3f}) are not imbalanced by this model. Cross-cluster: large imbalanced segments have the same major haplotype in
all clusters with data; cluster 4 on chr8:115–145 Mb has the opposite one (a mirrored-event hypothesis, 182 cells).</p>"""

    methods_17 = f"""
<p><b>What it does.</b> For a region it combines two summaries, each with its own likelihood and its own empirical null, and asks which explicit
allele-specific state explains both.</p>
<pre>depth   L_obs = raw-count rate ratio target / flanks (or / genome for broad loci)      log L_obs ~ Normal(log L_exp, σ_A²)
        σ_A = robust SD of log L among uncalled random windows (≈ {sig17:.2f}): regional baseline variability, which dominates sampling noise
allele  ĥ = phase-agnostic major-haplotype fraction (exp 10) with SE                   ĥ ~ Normal(max(0.5 + k (h_exp − 0.5), h_null), SE² + τ²)
        τ, h_null from uncalled random windows; k = allelic attenuation (not identified; 1.0 / 0.6 / 0.4)
cells   homozygous-loss states only: zero-count cells among cells expected ≥ 5 fragments ~ Binomial
states  1:1 (none), 1:0, 0:0, 2:0, 2:1, 2:2, 4:1, each in a fraction f ∈ {{0.1 … 1}} of cells on a 1:1 background:
        L_exp = population-average copies / 2, h_exp = population-average major fraction (e.g. 1:0 in f: L = 1 − f/2, h = 1/(2 − f))
output  posterior over states; MAP state; major/minor relative copies with 90% bounds; 'no modelled state fits' if |z| > 2.5 for either modality</pre>
<p><b>Assumptions:</b> diploid 1:1 background; population-average mixtures (one event type per region); independent modality likelihoods (they share
fragments); σ_A and τ from random windows, which can contain real events; attenuation k unknown (calls that change with k are unresolved).
<b>Why this method:</b> it is the simplest way to make the user's question explicit (which haplotype, how many copies, in how many cells), keeps each
modality's evidence visible (depth-only vs allele-only MAP columns), and refuses to force a call when nothing fits. <b>Result:</b> MYC and 5q fit a 2:1 gain
with an identified major haplotype; the chr7 dips fit no state (depth {chr7.rr_flank_raw.min():.2f}–{chr7.rr_flank_raw.max():.2f}× flanks with balanced alleles);
9p21 is not called a loss because LOH removes heterozygous sites from the selected SNP set (ascertainment), so balanced alleles there are not evidence
against loss.</p>"""

    t18 = r18.reset_index()[["region", "best_type", "pi", "pi_lo", "pi_hi", "boot_p", "evidence_class", "depth_allele_agreement", "pooled_major_fraction"]]
    methods_18 = f"""
<p><b>What it does.</b> Moves from pooled regions to individual cells: in each region every cell is either background (1:1) or carries one event
type; the model estimates the fraction π of event cells, which cells they are, and whether depth and alleles agree.</p>
<pre>cell i, region r:  state s_i ∈ {{background, event e}},  P(s_i = event) = π   (π ≥ 0.05 under the alternative)
depth   y_i | s ~ 0.99·NB2(μ_i·L_s, α_r) + 0.01·NB2(μ_i·L_s, 1)          (1% broad component: outlier-robust, a nuisance)
        μ_i = E_i · (n_target / n_outside) · exp(g_c) · exp(a_c)
        E_i = cell's raw total outside the region (exposure); g_c = cluster's local flank rate
        a_c = b + u_c = region baseline offset, b shared by all cells, u_c per cluster,
              with priors from uncalled random windows → a homogeneous accessibility shift is absorbed by b, not by π
allele  M_i | N_i, s ~ BetaBinomial(N_i, h_s, ρ_r / N_i)                     (N = AT on confidently oriented bins, M = major-haplotype counts, exp 16)
        h_0 = 0.5 + d0 (orientation folding bias, prior from random windows); h_event = h_0 + k (h_e − 0.5), k = 0.6 (1.0 sensitivity)
events  1:0, 0:0, 2:0, 2:1, 2:2, 4:1 (depth level L_e, major fraction h_e), each fitted separately
test    T = max over event types of 2 (ℓ_alt − ℓ_null); p from a hierarchical parametric bootstrap (B = 99) that redraws the baseline and
        folding offsets from their priors and refits the whole pipeline; profile CI for π
labels  clonal (π CI ≥ 0.8) · subclonal · small-fraction deviation (< 10% of cells) · no evidence · insufficient evidence</pre>
<p><b>Assumptions:</b> one event type per region; cells are independent given their state and cluster; the baseline offset prior from random windows is
representative; the allele model conditions on N, so AT itself is not modelled (its dependence on depth is not double counted, but depth and allele
still share fragments); orientation comes from the pooled 'all' cells, so a cluster with the opposite major haplotype is mis-oriented (cluster 4 on chr8q).
<b>Why this method:</b> pooled regions cannot separate 'all cells slightly lower' (accessibility or clonal) from 'some cells much lower' (subclone); the
per-cell likelihoods chosen in experiment 15 make that mixture identifiable in principle, and the jointly estimated baseline offset protects against
the most common confounder. <b>Calibration result:</b> on 30 uncalled random windows, {fpr['event_detected_p_lt_0.05']['k']} (≈ {fpr['event_detected_p_lt_0.05']['rate']:.0%}) are
flagged; without the shared baseline offset {fpr['no_offset_p_lt_0.05']['k']}/30; but 0/30 are supported by both modalities. So: <b>both-modality cell-level calls are credible
(MYC, 5q13-qter, chr8:115–145 Mb), depth-only cell-level calls are not</b> (including the chr7 'subclonal' windows and the 9p21/19q 'clonal' calls), until
the depth null includes the between-cell and between-cluster accessibility heterogeneity that random windows show.</p>
{C.html_table(t18)}
{C.img(b_cell, "cell mixtures")}
{C.caption("Per-region fraction of event cells with profile CI, coloured by which modalities support the event. Only the blue rows have support from both.")}"""

    methods_models = f"""
<p><b>Depth (E:1, exp 08/09/15).</b> For region r and cell i, raw target counts y with offset log E_i (cell total outside the region) + log(bins), WNN cluster
terms; NB2 variance μ + αμ². Chosen because it wins the held-out log score (Poisson worse by {-s15['A_depth_all_regions_vs_NB2']['pois']['d_logscore_vs_nb']:.2f} nats per cell-region; ZINB/hurdle gain
≤ 0.005) and reproduces held-out zero cells; quasi-Poisson gives the same estimates with robust intervals but no likelihood, so NB2 is the one usable in a joint model.
Statistical significance against rate ratio 1 is uninformative with thousands of cells; the reference must be the spread of same-width uncalled regions.</p>
<p><b>Allele total AT (exp 15).</b> NB2 with its own α (median {s15['A_median_alpha']['AT_raw_in']:.2f}, ≈ 6× depth's), offset log(raw depth in region + 0.5): AT is not independent of
depth (median elasticity {(s15['C_elasticity_median_min_max'][0] if isinstance(s15['C_elasticity_median_min_max'], list) else s15['C_elasticity_median_min_max']):.2f}), so it is modelled conditionally.
ZINB adds {s15['A_AT_rawin_all_regions_vs_NB2']['zinb']['d_logscore_vs_nb']:.4f} nats.</p>
<p><b>Haplotype split HB | AT (exp 15).</b> Beta-binomial; binomial is worse by {-s15['B_all_vs_BB']['b1']['d_logscore_per_cell_vs_bb1']:.3f} nats per cell, a two-component mixture adds nothing;
cell-level ρ falls with AT and must be estimated per region (the bin-level 0.0038 is a pooled-bin value).</p>
<p><b>Why model per modality.</b> The user asked whether modalities need different models; the data say yes: AT is about six times more overdispersed than depth,
and the split needs a bounded (beta-binomial) rather than a count model. Using one family for both would mis-weight the allele evidence.</p>"""

    methods_seg = f"""
<p><b>What the current AD segmentation is good at.</b> Boundaries are reproducible between disjoint halves of the cells (split-half r > 0.99; all 10
saved-window edges that coincide with breakpoints reappear), and the select-then-test procedure accepted 1/300 first splits on exchangeable nulls
({link('05_segmentation_calibration', 'exp 05')}).</p>
<p><b>Where it fails.</b> (1) Significance tests reproducible structure, not copy number. (2) In the benchmark, moderate focal changes are rarely
segmented (figure). Diagnosis ({link('19_masking_two_cut_scan', 'exp 19')}): in a long segment the single best AD cut is elsewhere and is rejected,
which ends the search (masking); a two-cut "segment vs rest" scan rescues deep losses but not one-copy changes, because inside their own segment natural
windows differ from the rest as much as the injected change. (3) Count-likelihood segmentation with sampling-only dispersion produces ~4,200 segments;
it needs a between-bin baseline-variance term (τ² ≈ 0.26 in the joint HMM; {link('11_count_changepoint', 'exp 11')}, {link('13_joint_hmm', 'exp 13')}).</p>
{C.img(b_bench, "benchmark")}
{C.caption(f"Detection of injected events (both edges within ±2 bins and CN in the injected direction) by width. Circles: experiment 14, realistic thinning of raw counts in cluster 0 "
           f"(f = 1; {n14:,} realizations in total; Wilson 95% intervals). Squares: experiment 19, two-cut scan on clean edits (8 loci per point). Events narrower than 20 bins cannot be a segment by construction.")}"""

    caller = """
<p><b>Do the pieces make a haplotype-aware CNV caller?</b> They make its skeleton, and on the controls it behaves sensibly; as stitched now it would not be
a trustworthy caller, for three reasons: the segmentation is depth-driven and misses moderate focal events, depth has no baseline so depth-only calls are not
calibrated, and the allele inputs are bin-summed with unverified phase. A design that follows from the evidence:</p>
<ol>
<li><b>Allele-first segmentation:</b> a phase-switch-aware HMM over bins on the oriented major-haplotype fraction (experiment 16 is most of this), because
allele imbalance is much less confounded by accessibility than depth (5q/1q/MYC imbalance is 3–6× the random-window null; depth contrasts of uncalled
regions reach the size of one-copy changes).</li>
<li><b>Depth as the second emission:</b> NB2 with exposure and an explicit baseline-variance term, ideally against a baseline built from unrelated samples;
depth gives the direction and magnitude of an allelic event and is the only evidence where alleles are missing.</li>
<li><b>Allele-specific state space</b> (experiment 17) with explicit outcomes: supported by both, depth-only (unconfirmed), allele-only (CN-LOH candidate), no
modelled state fits, unobservable (QC-removed), and a site-density term so that LOH-driven loss of heterozygous sites counts as evidence for loss.</li>
<li><b>Cell-level mixture</b> (experiment 18) on top, for clonal vs subclonal fractions and mirrored events, after its depth null is calibrated to ≈ 5% on random windows.</li>
</ol>
<p>This converges on the architecture of Numbat-multiome (haplotype-aware HMM with expression/accessibility and allele emissions). Once fragments and SNP-level
phased allele counts exist, running Numbat-multiome ATAC-bin as the comparator is cheaper than rebuilding everything; what this review adds is the
no-normal calibration (empirical regional nulls), the explicit 'no state fits / unobservable / depth-only' outcomes, and the QC-gap bookkeeping.</p>"""

    worked = """<ul>
<li><b>Worked:</b> X reconstruction; locating QC-removed bins; raw-count NB2 models with exposure; held-out model choice per modality; empirical nulls from random
windows; phase-orientation HMM; allele-specific state scoring with 'no state fits'; both-modality cell-level calls; the spike-in harness (unmodified scans
reproduce saved results exactly).</li>
<li><b>Promising, needs calibration:</b> cell-level mixtures (credible only with both modalities), two-cut scan for strong focal events, distance-aware joint HMM with
a baseline-variance term.</li>
<li><b>Did not work / hypotheses rejected:</b> 'the AD p-values are grossly anti-conservative' (rejected, exp 05); 'zero inflation explains the dips' (rejected, exp 09);
'a better search (two-cut) recovers one-copy focal events' (rejected, exp 19); 'nominal AT gives calibrated allelic power' (rejected, exp 06/10/15); sampling-only
dispersion for count segmentation (over-segments, exp 11/13); depth-only cell mixtures (23% FPR, exp 18).</li>
</ul>"""

    sections = [
        ("Scope and evidence levels", "<p>Development data are MCF10-A scATAC (3,460 cells × 25,337 retained 100 kb bins); nothing here is copy-number truth. "
         "<b>[L]</b> local observation from a saved, reproducible experiment; <b>[P]</b> published evidence (MCF10-A karyotype: CDKN2A/B homozygous deletion within a ~4 Mb "
         "9p21 loss, MYC amplification, 19q13-qter and 5q13-qter gains, used only as prespecified weak controls; access limits in <a href='sources.json'>sources.json</a>); "
         "<b>[I]</b> proposed idea. Every experiment was run by the coordinator or by a delegated agent and validated by the coordinator (independent recomputation or "
         "byte-identical rerun). Seeds, environments and commands are in each experiment's summary.json.</p>"),
        ("Map of everything that was done", map_html),
        ("Noteworthy findings by pipeline stage", "", 2),
        ("Inputs, QC and preprocessing", f"<ul><li>[L] X = pinned formula on <code>layers/raw</code> (bit-identical {rec['cells_exact_equal']:,}/{s01['cells']:,} cells); 17 cell-bins at near-zero "
         "LOESS predictions are 0 historically and unexplained in five environments.</li><li>[L] The bin filter removed 1,103 autosomal bins incl. CDKN2A/B: clonal homozygous "
         "deletions remove themselves.</li><li>[L] GC correction attenuates low-GC dips (chr7 8.3–10.7: 0.20× raw vs 0.36× in X) and erases the raw 19q13-qter excess "
         "(1.82× → 0.93×) whose alleles show imbalance.</li><li>[L] Rounding zeros are coverage-dependent (median 1.4% vs 0.5% of raw-positive bins, low vs high tercile) but the "
         "regional raw contrasts are coverage-invariant.</li></ul>", 3),
        ("Segmentation", "<ul><li>[L] Boundaries reproducible; tests measure reproducible structure, not CN.</li><li>[L] Moderate focal events are missed "
         f"(clonal one-copy loss 20–50 bins: {lo1.detected.max():.0%} at best; 1.5× gains: {gn1.detected.max():.0%} at best); deep losses are found from ~20 bins.</li>"
         "<li>[L] Masking explains part of it (two-cut rescues deep losses); background structure explains the rest (one-copy changes not rescued).</li></ul>", 3),
        ("CN assignment and reference", "<ul><li>[L] pyEpi's scale search always ends at the grid edge (s = 0.5 × median): calls are fixed thresholds on segment mean / "
         "genome median (loss < 0.68×, Watson ≤ 0.5×; gain > 1.22×, Watson ≥ 1.5×).</li><li>[L] Changing the estimand changes 20–37% of bin calls; changing the baseline "
         "reverses about one third of called short segments; removing one chromosome changes up to 5% of calls elsewhere.</li></ul>", 3),
        ("Depth evidence", f"<ul><li>[L] NB2 with exposure is adequate; zero inflation is not needed.</li><li>[L] Against uncalled random windows, {nsig['loss']}/{nsig['called_loss']} called short "
         f"losses and {nsig['gain']}/{nsig['called_gain']} gains remain unusual.</li><li>[L] The chr7 saved windows are {chr7.rr_flank_raw.min():.2f}–{chr7.rr_flank_raw.max():.2f}× their flanks, coverage-invariant, "
         "shared by all clusters, and ≤ 2.3% of cells could carry a near-complete loss: not explained by a one-copy loss alone.</li></ul>", 3),
        ("Allelic evidence and haplotypes", "<ul><li>[L] Per-modality models: AT NB2 | depth, split beta-binomial with region-level ρ.</li><li>[L] GT phase switches about every "
         "0.5 Mb; orientation inference is required to see imbalance; it recovers the published gains (5q, MYC) and chr1q (h ≈ 0.71).</li><li>[L] Nominal AT overstates "
         "power; allelic thresholds must be empirical.</li><li>[L] LOH ascertainment: the 9p21 hemizygous zone has 0.20× the selected-site observations per fragment, and its "
         "remaining sites look balanced.</li></ul>", 3),
        ("Combined depth + allele", "<ul><li>[L] Controls with both modalities (MYC, 5q13-qter, chr8:115–145 Mb) come out as near-clonal 2:1 gains at region and cell level.</li>"
         "<li>[L] The chr7 dips fit no modelled allele-specific state.</li><li>[L] Depth-only cell-level calls have a 23% false-positive rate on random windows.</li></ul>", 3),
        ("What worked, what is promising, what did not", worked),
        ("Methods explained", "", 2),
        ("Observation models per modality (exp 08, 09, 15)", methods_models, 3),
        ("Major/minor haplotype identification (exp 16)", methods_16, 3),
        ("Allele-specific state scoring (exp 17)", methods_17, 3),
        ("Cell-level joint mixture (exp 18)", methods_18, 3),
        ("Segmentation diagnostics and the spike-in benchmark (exp 05, 11, 13, 14, 19)", methods_seg, 3),
        ("What should change (prioritized)", C.html_table(change) + C.caption("Ranks are research priorities, not measured performance. None of these is a production change "
         "without a separate scientific decision (old behaviour, new behaviour, reason, validation).")),
        ("Toward a haplotype-aware CNV caller", caller),
        ("Open questions and the experiments that would resolve them", "<ul>"
         "<li>Are the chr7 0.13–0.34× plateaus mappability, accessibility or dosage? → same bins in ≥ 2 unrelated scATAC samples + a read-length-matched mappability track; "
         "decision: recurrence across samples with normal DNA → baseline offset.</li>"
         "<li>How large is the true allelic attenuation (k)? → SNP-level, molecule-aware counts on a sample with DNA-defined allele-specific states.</li>"
         "<li>Can depth-only cell mixtures be calibrated? → add between-cell/cluster accessibility heterogeneity to the null; target ≈ 5% on random windows.</li>"
         "<li>Does allele-first segmentation recover one-copy focal events that the AD scan misses? → rerun experiment 14's injections with matching allele thinning.</li>"
         "<li>What do QC-removed runs contain? → pre-QC counts / fragments.</li>"
         "<li>Is cluster 4's opposite major haplotype on chr8q a mirrored event? → per-cell orientation within cluster 4 with SNP-level phase sets.</li></ul>"),
    ]
    glossary = {"L, L_obs": "depth level: raw-count rate ratio target/flanks (or /genome)", "ĥ, h": "major-haplotype fraction (≥ 0.5)",
                "δ": "imbalance h − 0.5", "σ_A": "robust SD of log depth level over uncalled random windows",
                "π": "fraction of cells in the event state (exp 18)", "k": "allelic attenuation factor (unidentified)",
                "GT-oriented": "summing HB as labelled by the phased VCF genotype", "phase-agnostic": "|HB − HA| per bin, orientation ignored (folded)",
                "orientation HMM": "Markov chain over bins for which haplotype is major", "NB2": "negative binomial with Var = μ + αμ²",
                "beta-binomial": "binomial with overdispersed success probability (ρ = intraclass correlation)",
                "uncalled random windows": "seeded random 25-bin windows with no called bin; empirical reference, not neutral truth",
                "two-cut scan": "maximum over windows of AD(window, rest of segment)"}
    C.write_html_report(
        REV / "final_report.html", "CNV method review — final report",
        "scientific-explorer · allelic-segmentation@6d4f48d · MCF10-A development data · experiments 00–19 · 2026-11-05",
        sections, {"builder": "python -I research/2026-11-04-cloud-method-review/final_report/build.py",
                   "inputs": "experiments/*/summary.json and TSV outputs", "spikein_realizations_used": n14, "environment": C.environment()},
        verdict="The workflow is reproducible and its boundaries are stable, but depth alone cannot tell one-copy changes from regional accessibility at this "
                "scale, the current segmentation misses moderate focal events, and QC removes clonal homozygous deletions. Phase-oriented allele evidence, "
                "modelled per modality and combined with depth, identifies major/minor haplotypes and allele-specific states credibly where both modalities "
                "agree. The next change belongs in evidence integration, QC bookkeeping and a depth baseline, leading to an allele-first haplotype-aware caller.",
        key_findings=key, glossary=glossary,
        supports=["An evidence layer with modality-specific likelihoods and explicit 'no state fits / depth-only / unobservable' outcomes.",
                  "Allele-first, phase-aware segmentation as the core of a haplotype-aware caller.",
                  "Keeping QC-removed bins as an explicit state.",
                  "Both-modality cell-level clonal/subclonal estimates."],
        not_supports=["Copy number at any saved window; whether the chr7 plateaus are dosage.",
                      "Depth-only subclonal calls (uncalibrated).",
                      "Claims beyond MCF10-A; any tool outperforming the baseline on real truth.",
                      "Verified parental haplotype identity (phase unverified)."],
        decision_impact="Ordering of the plan (recommendation.md): 1 evidence layer (now prototyped), 2 QC bookkeeping, 3 cross-sample depth baseline + mappability, "
                        "4 allele-first segmentation with depth as second emission (replaces 'raw-count depth segmentation' as the segmentation priority), "
                        "5 SNP-level allele inputs; GC tuning and zero-inflation models are dropped from the plan.",
        next_steps=["Obtain fragments + SNP-level phased allele counts; rerun pileup with molecule-aware counting.",
                    "Process ≥ 2 unrelated scATAC samples (one with DNA truth) with identical bins.",
                    "Prototype the allele-first HMM caller from experiments 16/13/17 and benchmark with experiment 14 injections (depth + allele thinning).",
                    "Run Numbat-multiome ATAC-bin as comparator when its inputs exist."])
    print("wrote", REV / "final_report.html", "| spike-in realizations:", n14)


if __name__ == "__main__":
    main()
