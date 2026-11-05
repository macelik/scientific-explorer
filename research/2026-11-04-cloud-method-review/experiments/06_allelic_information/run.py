"""Experiment 06: allelic observability, information and evidence labels for fixed regions.

Inputs: HA/HB/AT parquet (selected phased heterozygous SNP-site DP summed per bin; AT is a
site-observation count, not established unique fragments), AnnData layers/raw, WNN labels.
Joins by barcode name; checks coordinates and AT == HA + HB for every cell.

Per region (the same rule-defined families as experiment 04; pooled over all 3,460 cells):
  at_per_raw      - pooled AT / pooled raw fragments; relative to flanks and genome. Low values mean
                    fewer selected-site observations per fragment (fewer heterozygous sites called,
                    LOH-driven ascertainment, homozygosity or technical; not distinguishable here).
  S               - mean of z_j^2 = (HB_j - HA_j)^2 / AT_j over bins with AT_j >= 5. Phase-agnostic:
                    the sign of each bin's imbalance is not used, so inconsistent orientation
                    between bins does not cancel signal.
  thresholds      - model null: beta-binomial(p = 0.5, rho) with rho estimated from bins of `all`
                    segments called base (cn_call == 1; not truth), simulated with each region's own
                    AT_j; empirical null: 95th percentile of S over random windows.
  power           - same simulation with p = 0.5 +/- delta (random sign per bin), for a one-copy loss in
                    a fraction f of cells, delta = f / (2 (2 - f)), f in {1, 0.5, 0.25}, and a one-copy
                    gain, delta = f / (2 (2 + f)), f in {1, 0.5}. Assumes the region's sites stay
                    observable; under clonal LOH without a matched normal, sites are not ascertained.
  label           - unobserved (< 3 bins with AT >= 5); insufficient (power for clonal one-copy loss
                    < 0.8); consistent with imbalance (S above both null 95th percentiles);
                    against clonal one-copy imbalance at observed sites (S below the 5th
                    percentile under that alternative); otherwise compatible with balance, not excluded.
Seeds: numpy PCG64 default_rng((20261104, 6, region_index)); 2,000 simulations per scenario.
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

NSIM = 2000
MIN_AT = 5
SCEN = {"loss_f1": 1 / 2, "loss_f0.5": 0.5 / (2 * 1.5), "loss_f0.25": 0.25 / (2 * 1.75),
        "gain_f1": 1 / (2 * 3), "gain_f0.5": 0.5 / (2 * 2.5)}


def bb_sample(rng, n, p, rho, size):
    p = np.clip(p, 1e-3, 1 - 1e-3)   # delta = 0.5 (complete LOH) would otherwise give p = 0 or 1
    if rho <= 0:
        return rng.binomial(n, p, size=(size, len(n)))
    a = p * (1 - rho) / rho
    b = (1 - p) * (1 - rho) / rho
    q = rng.beta(a, b, size=(size, len(n)))
    return rng.binomial(n, q)


def main():
    pth = C.default_paths()
    d = C.H5AD(pth["h5ad"])
    cells = list(d.cells)
    mats = {}
    checks = {}
    for nm in ("HA", "HB", "AT"):
        t = pq.read_table(pth["allelic"] / f"{nm}.parquet")
        df_idx = t.select([c for c in t.column_names if c in ("chr", "start", "end")]).to_pandas()
        if "chr" not in df_idx.columns:   # coordinates stored in the pandas index
            df_idx = t.to_pandas()[[]].reset_index()
        ok = (df_idx["chr"].astype(str).to_numpy() == d.seq).all() and (df_idx["start"].to_numpy() == d.start).all() \
            and (df_idx["end"].to_numpy() == d.end).all()
        names = [c for c in t.column_names if c not in ("chr", "start", "end", "__index_level_0__")]
        checks[nm] = {"coords_match_anndata": bool(ok), "n_barcodes": len(names), "barcode_set_equal": set(names) == set(cells),
                      "column_order_equal_anndata": names == cells}
        if not ok or set(names) != set(cells):
            raise SystemExit(f"alignment failure in {nm}: {checks[nm]}")
        mats[nm] = np.column_stack([t.column(c).to_numpy(zero_copy_only=False) for c in cells]).astype(np.float64)  # bins x cells
        del t
    HA, HB, AT = mats["HA"], mats["HB"], mats["AT"]
    checks["AT_equals_HA_plus_HB_all_cells"] = bool(np.array_equal(AT, HA + HB))
    checks["nonnegative"] = bool((HA >= 0).all() and (HB >= 0).all())
    checks["integer_valued"] = bool(np.all(np.mod(AT, 1) == 0))
    ha, hb, at = HA.sum(1), HB.sum(1), AT.sum(1)
    del mats, HA, HB, AT
    raw = C.group_means(d, "layers/raw", {"all": np.ones(d.n_cells, bool)})["all"] * d.n_cells

    # rho from bins inside `all` base-called segments
    kseg = pd.read_csv(pth["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    call = np.ones(d.n_bins)
    for r in kseg[kseg.group == "all"].itertuples():
        call[r.bin_start:r.bin_end + 1] = r.cn_call
    z2 = np.where(at > 0, (hb - ha) ** 2 / np.where(at > 0, at, 1), np.nan)
    base = (call == 1) & (at >= 20)
    strata = pd.qcut(at[base], 5)
    rho_rows = []
    for s, g in pd.DataFrame({"z2": z2[base], "n": at[base], "s": strata}).groupby("s", observed=True):
        var = np.median(g.z2) / 0.4549   # median of chi2_1 = 0.4549
        rho_rows.append({"AT_stratum": str(s), "bins": len(g), "median_AT": g.n.median(), "variance_inflation": var,
                         "rho_hat": (var - 1) / (g.n.median() - 1)})
    rho_tab = pd.DataFrame(rho_rows)
    rho = float(np.median(rho_tab.rho_hat.clip(lower=0)))

    regs = E4.build_regions(d, pth)
    gen_atr = at.sum() / raw.sum()
    out = []
    for i, r in enumerate(regs):
        t, f = r["target"], r["flank"]
        rng = np.random.default_rng((C.SEED, 6, i))
        rec = {"family": r["family"], "region": r["name"], "direction": r["direction"], "n_target_bins": len(t),
               "AT_sum": at[t].sum(), "raw_sum": raw[t].sum(), "at_per_raw_vs_genome": (at[t].sum() / raw[t].sum()) / gen_atr,
               "bins_AT_ge5": int((at[t] >= MIN_AT).sum()), "bins_AT_0": int((at[t] == 0).sum())}
        if len(f):
            rec["at_per_raw_vs_flank"] = (at[t].sum() / raw[t].sum()) / (at[f].sum() / raw[f].sum())
        m = t[at[t] >= MIN_AT]
        rec["S_obs"] = float(np.mean(z2[m])) if len(m) else np.nan
        if len(m) >= 3:
            n = at[m]
            v = 1 + (n - 1) * rho
            e = (z2[m] - v) / (4 * n)                    # per-bin unbiased estimate of delta^2 under the BB model
            rec["delta_eff"] = float(np.sqrt(max(0.0, e.mean())))
            bs = rng.choice(len(e), size=(500, len(e)))
            rec["delta_eff_ci_lo"], rec["delta_eff_ci_hi"] = [float(np.sqrt(max(0.0, q))) for q in np.quantile(e[bs].mean(1), [.025, .975])]
        rec["frac_bins_called_loss_all"] = float((call[t] < 1).mean()) if len(t) else np.nan
        rec["frac_bins_called_gain_all"] = float((call[t] > 1).mean()) if len(t) else np.nan
        rec["_m"] = m
        out.append(rec)
    # empirically anchored alternative: the effective |dBAF| measured at the published MYC amplification
    ctrl = {r["region"]: r.get("delta_eff") for r in out if r["family"] == "literature"}
    SCEN["anchored_MYC_delta_eff"] = ctrl["lit_MYC_8q24"]
    SCEN["anchored_half_MYC"] = ctrl["lit_MYC_8q24"] / 2
    for i, rec in enumerate(out):
        m = rec.pop("_m")
        if len(m) < 3:
            continue
        rng = np.random.default_rng((C.SEED, 6, i, 1))
        n = at[m].astype(int)

        def S(p):
            k = bb_sample(rng, n, p, rho, NSIM)
            return np.mean((2 * k - n) ** 2 / n, axis=1)
        rec["S_null_model_q95"] = float(np.quantile(S(np.full(len(n), 0.5)), .95))
        for nm, dl in SCEN.items():
            sign = rng.choice([-1, 1], size=len(n))
            rec[f"power_{nm}"] = float((S(0.5 + sign * dl) > rec["S_null_model_q95"]).mean())
    df = pd.DataFrame(out)
    rnd = df[df.family == "random"]
    unc = rnd[(rnd.frac_bins_called_loss_all == 0) & (rnd.frac_bins_called_gain_all == 0)]
    emp95 = float(np.nanquantile(unc.S_obs, .95))
    model_calib = {"random_uncalled": len(unc), "frac_uncalled_random_S_above_model_q95": float((unc.S_obs > unc.S_null_model_q95).mean()),
                   "frac_all_random_S_above_model_q95": float((rnd.S_obs > rnd.S_null_model_q95).mean()),
                   "control_delta_eff": ctrl}

    def label(r):
        """Evidence labels. No 'against' label: the alternative's attenuation is not calibrated (see report)."""
        if r.bins_AT_ge5 < 3:
            return "allele unobserved"
        if r.S_obs > r.S_null_model_q95 and r.delta_eff_ci_lo > 0:
            return "consistent with imbalance"
        if r.power_anchored_MYC_delta_eff < 0.8:
            return "insufficient information"
        return "no excess imbalance at observed sites"
    df["S_null_empirical_random_q95"] = emp95
    df["label"] = df.apply(label, axis=1)
    df.to_csv(HERE / "region_allelic_information.tsv", sep="\t", index=False)
    rho_tab.to_csv(HERE / "overdispersion_estimate.tsv", sep="\t", index=False)
    fam = df.groupby(["family", "direction", "label"]).size().rename("regions").reset_index()
    fam.to_csv(HERE / "label_counts.tsv", sep="\t", index=False)
    atr = df.groupby(["family", "direction"]).agg(n=("region", "size"), median_at_per_raw_vs_genome=("at_per_raw_vs_genome", "median"),
                                                  median_at_per_raw_vs_flank=("at_per_raw_vs_flank", "median"),
                                                  median_power_loss_f1=("power_loss_f1", "median"),
                                                  median_power_loss_f0_5=("power_loss_f0.5", "median"),
                                                  median_power_gain_f1=("power_gain_f1", "median"),
                                                  median_delta_eff=("delta_eff", "median"),
                                                  median_power_anchored_MYC=("power_anchored_MYC_delta_eff", "median"),
                                                  median_power_anchored_half_MYC=("power_anchored_half_MYC", "median")).reset_index()
    res = {"experiment": "06_allelic_information", "seed": "default_rng((20261104, 6, region_index))", "n_sim": NSIM, "min_AT": MIN_AT,
           "alignment_checks": checks, "rho_used": rho, "scenarios_delta": SCEN, "empirical_uncalled_random_S_q95": emp95, "model_null_calibration": model_calib,
           "label_counts": fam.to_dict("records"), "family_medians": atr.to_dict("records"), "environment": C.environment()}
    C.write_json(HERE / "summary.json", res)

    show = ["region", "n_target_bins", "AT_sum", "bins_AT_ge5", "bins_AT_0", "at_per_raw_vs_genome", "at_per_raw_vs_flank", "S_obs",
            "S_null_model_q95", "power_loss_f1", "power_loss_f0.5", "power_loss_f0.25", "power_gain_f1", "delta_eff", "delta_eff_ci_lo", "delta_eff_ci_hi", "power_anchored_MYC_delta_eff", "power_anchored_half_MYC", "label"]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4))
    for f, c in [("random", "#bbb"), ("called", "#9b2c2c"), ("saved", "k"), ("literature", "#2f5d8a")]:
        s = df[df.family == f]
        ax[0].scatter(s.at_per_raw_vs_genome, s.S_obs, s=12 if f != "saved" else 35, color=c, label=f, alpha=.75)
    ax[0].axhline(emp95, color="#8a5a00", ls="--", lw=1, label="uncalled random-window S 95th pct")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("AT per raw fragment, relative to genome")
    ax[0].set_ylabel("S = mean (HB−HA)²/AT"); ax[0].legend(fontsize=7); ax[0].set_title("Allelic observability vs excess allelic variance")
    s = df[df.family == "random"]
    ax[1].scatter(s.AT_sum, s.power_loss_f1, s=8, color="#bbb", label="random: clonal 1-copy loss")
    ax[1].scatter(s.AT_sum, s["power_loss_f0.25"], s=8, color="#9b2c2c", alpha=.5, label="random: 1-copy loss in 25% of cells")
    sv = df[df.family == "saved"]
    ax[1].scatter(sv.AT_sum, sv.power_loss_f1, s=35, color="k", label="saved windows (clonal)")
    ax[1].set_xscale("log"); ax[1].set_xlabel("pooled AT in region"); ax[1].set_ylabel("power"); ax[1].legend(fontsize=7)
    ax[1].set_title("Power at the region's own AT profile")
    b64 = C.fig_to_b64(fig); plt.close(fig)
    dfi = df.set_index("region")
    sv = df[df.family == "saved"]
    sv7 = sv[sv.region.str.startswith("chr7")]
    lit = df[df.family == "literature"].set_index("region")
    myc = lit.loc["lit_MYC_8q24"]
    n_mean_myc = myc.AT_sum / myc.n_target_bins
    pred_myc = 1 + (n_mean_myc - 1) * rho + 4 * SCEN["gain_f1"] ** 2 * n_mean_myc
    n_unc_above = int(round(model_calib["frac_uncalled_random_S_above_model_q95"] * model_calib["random_uncalled"]))
    lab_tab = df.groupby(["family", "direction", "label"]).size().rename("regions").reset_index()
    lab_rand = rnd.assign(label=df.loc[rnd.index, "label"]).label.value_counts()
    verdict = (f"The chr7 saved windows keep allele observations ({int(sv7.bins_AT_ge5.min())}&ndash;{int(sv7.bins_AT_ge5.max())} bins with AT &ge; {MIN_AT}) and show no excess allelic variance (effective |dBAF| upper bound &le; {sv7.delta_eff_ci_hi.max():.2f}), "
               f"which would have been detectable if as strong as at MYC (power {sv7.power_anchored_MYC_delta_eff.min():.2f}&ndash;{sv7.power_anchored_MYC_delta_eff.max():.2f}) but usually not at half that strength; "
               f"nominal AT overstates allelic information (MYC effective |dBAF| {myc.delta_eff:.2f} against a nominal prediction of S &asymp; {pred_myc:.0f} versus observed {myc.S_obs:.1f}).")
    findings = [
        f"Alignment: AT = HA + HB exactly for every cell and bin ({checks['AT_equals_HA_plus_HB_all_cells']}); values are non-negative integers ({checks['nonnegative'] and checks['integer_valued']}); coordinates match the AnnData; barcode sets are equal but column order differs, so the join is by name.",
        f"Overdispersion: &rho; estimated on base-called bins is {rho:.4f} (median over {len(rho_tab)} AT strata, range {rho_tab.rho_hat.min():.4f}&ndash;{rho_tab.rho_hat.max():.4f}), while the variance inflation rises from {rho_tab.variance_inflation.iloc[0]:.2f} to {rho_tab.variance_inflation.iloc[-1]:.2f} with AT.",
        f"The balanced beta-binomial model underpredicts allelic variance in uncalled random windows: {n_unc_above} of {model_calib['random_uncalled']} ({100*model_calib['frac_uncalled_random_S_above_model_q95']:.0f}%) exceed their own model 95th percentile ({100*model_calib['frac_all_random_S_above_model_q95']:.0f}% of all {len(rnd)} random windows); the empirical 95th percentile of S in uncalled random windows is {emp95:.1f}.",
        f"Attenuation at positive controls: MYC S = {myc.S_obs:.1f} (model 95th percentile {myc.S_null_model_q95:.1f}), effective |dBAF| {myc.delta_eff:.3f} (bootstrap 95% {myc.delta_eff_ci_lo:.3f}&ndash;{myc.delta_eff_ci_hi:.3f}), far below the nominal clonal-gain value {SCEN['gain_f1']:.3f} that would give S &asymp; {pred_myc:.0f}; 5q13-qter {lit.loc['lit_5q13_qter','delta_eff']:.3f}; 19q13-qter {lit.loc['lit_19q13_qter','delta_eff']:.3f}.",
        f"chr7 saved windows: pooled AT {int(sv7.AT_sum.min()):,}&ndash;{int(sv7.AT_sum.max()):,}; S = {sv7.S_obs.min():.2f}&ndash;{sv7.S_obs.max():.2f} against model 95th percentiles {sv7.S_null_model_q95.min():.2f}&ndash;{sv7.S_null_model_q95.max():.2f}; effective |dBAF| = {sv7.delta_eff.max():.2f} in all four; selected-site observations per fragment {sv7.at_per_raw_vs_flank.min():.2f}&ndash;{sv7.at_per_raw_vs_flank.max():.2f}x the flanks.",
        f"Power anchored to MYC&rsquo;s effective imbalance ({SCEN['anchored_MYC_delta_eff']:.3f}) is {sv7.power_anchored_MYC_delta_eff.min():.2f}&ndash;{sv7.power_anchored_MYC_delta_eff.max():.2f} for the chr7 windows; at half that ({SCEN['anchored_half_MYC']:.3f}) it is {sv7.power_anchored_half_MYC.min():.2f}&ndash;{sv7.power_anchored_half_MYC.max():.2f}. Theoretical powers for nominal one-copy changes are 1.0 and are upper bounds not used for labels.",
        f"9p21 hemizygous zone: {lit.loc['lit_9p21_hemizygous_zone_approx','at_per_raw_vs_genome']:.2f}x genome selected-site observations per fragment ({int(lit.loc['lit_9p21_hemizygous_zone_approx','bins_AT_0'])} of {int(lit.loc['lit_9p21_hemizygous_zone_approx','n_target_bins'])} bins have AT = 0), the lowest of the literature loci, as expected when clonal LOH removes heterozygous sites.",
        f"Labels over {len(df)} regions: random windows {int(lab_rand.get('consistent with imbalance',0))} of {len(rnd)} 'consistent with imbalance'; short called losses {int(((df.family=='called')&(df.direction=='loss')&(df.label=='consistent with imbalance')).sum())} of {int(((df.family=='called')&(df.direction=='loss')).sum())}; "
        f"saved windows {int((sv.label=='consistent with imbalance').sum())} of {len(sv)} (chr9_70_100 comparison). The MYC control is labelled &lsquo;{myc.label}&rsquo; although its S exceeds the model 95th percentile, because the label rule also requires a bootstrap lower bound of effective |dBAF| above 0 and 15 informative bins give a lower bound of {myc.delta_eff_ci_lo:.3f}.",
    ]

    def region_note(r):
        s = f"{int(r.bins_AT_ge5)} of {int(r.n_target_bins)} bins have AT &ge; {MIN_AT}"
        if r.bins_AT_0:
            s += f" ({int(r.bins_AT_0)} with AT = 0)"
        s += f"; {r.at_per_raw_vs_genome:.2f}x genome selected-site observations per fragment"
        if np.isfinite(getattr(r, 'at_per_raw_vs_flank', np.nan)):
            s += f" ({r.at_per_raw_vs_flank:.2f}x flanks)"
        s += f". S = {r.S_obs:.2f} vs model 95th percentile {r.S_null_model_q95:.2f}"
        s += " (above)" if r.S_obs > r.S_null_model_q95 else " (not above)"
        s += f"; effective |dBAF| {r.delta_eff:.3f} (CI {r.delta_eff_ci_lo:.3f}&ndash;{r.delta_eff_ci_hi:.3f}). "
        s += f"Power for the MYC-sized imbalance {r.power_anchored_MYC_delta_eff:.2f}, for half of it {r.power_anchored_half_MYC:.2f}. "
        if r.label == "no excess imbalance at observed sites" and r.S_obs > r.S_null_model_q95 and r.delta_eff_ci_lo <= 0:
            s += "S exceeds the model threshold, but the bootstrap lower bound of effective |dBAF| is 0, so the preset label rule does not call it imbalance; read as unresolved at this bin count."
        elif r.label == "no excess imbalance at observed sites":
            s += "No excess allelic variance at observed sites" + (" and enough power to see a MYC-sized imbalance." if r.power_anchored_MYC_delta_eff >= 0.8 else "; too little power to call.")
        elif r.label == "consistent with imbalance":
            s += "Allelic variance in excess of the balanced model" + (", but below the empirical random-window threshold." if r.S_obs < emp95 else ", also above the empirical random-window threshold.")
        return s

    reg_show = pd.concat([sv, df[df.family == "literature"]])
    reg_tab = reg_show[["region", "label"]].copy()
    reg_tab["note"] = [__import__("html").unescape(region_note(r)) for r in reg_show.itertuples()]
    reg_html = C.html_table(reg_tab).replace("<td>", "<td style='text-align:left;white-space:normal'>")
    sv_cols = ["region", "n_target_bins", "AT_sum", "bins_AT_ge5", "bins_AT_0", "at_per_raw_vs_genome", "at_per_raw_vs_flank"]
    st_cols = ["region", "S_obs", "S_null_model_q95", "delta_eff", "delta_eff_ci_lo", "delta_eff_ci_hi", "label"]
    pw_cols = ["region", "power_loss_f1", "power_loss_f0.5", "power_loss_f0.25", "power_gain_f1", "power_anchored_MYC_delta_eff", "power_anchored_half_MYC"]
    region_df = reg_show
    lab_pivot = lab_tab.pivot_table(index=["family", "direction"], columns="label", values="regions", fill_value=0).reset_index()
    lab_pivot.columns.name = None
    atr_show = atr[["family", "direction", "n", "median_at_per_raw_vs_genome", "median_at_per_raw_vs_flank", "median_delta_eff", "median_power_anchored_MYC", "median_power_anchored_half_MYC"]]
    scen_tab = pd.DataFrame([{"scenario": k, "delta (|p - 0.5|)": v} for k, v in SCEN.items()])
    checks_tab = pd.DataFrame([{"matrix": k, **v} for k, v in checks.items() if isinstance(v, dict)])
    sections = [
        ("Question and why it matters",
         "<p><b>Question.</b> For each fixed region, are allele observations present, how many selected-site observations come per fragment, could a specified imbalance be detected, and which evidence label is supported? Missing alleles are kept distinct from balance.</p>"
         "<p><b>Why it matters.</b> The recommendation proposes an allelic observability and calibrated-label component in the per-region evidence layer, and a joint depth&ndash;allele model later. Both rest on whether allelic data are informative at this coverage and whether nominal counts can be used as independent observations. "
         "The chr7 dips are the key case: a clonal one-copy loss would normally leave allelic imbalance at heterozygous sites, while a shared baseline effect would not.</p>"),
        ("Inputs",
         f"<ul><li><b>Allelic matrices.</b> <code>aligned_matrices/HA.parquet</code>, <code>HB.parquet</code>, <code>AT.parquet</code>: bins &times; cells ({d.n_bins:,} &times; {d.n_cells:,}); HA and HB are the two haplotype counts and AT their sum, taken from selected phased heterozygous SNP-site depth (DP) summed per bin. AT is a site-observation count, not established unique fragments. Pooled over all cells.</li>"
         "<li><b>Depth.</b> <code>layers/raw</code> of <code>count_matrix.h5ad</code> (fragments), WNN labels not used here.</li>"
         f"<li><b>Regions.</b> The same rule-defined families as experiment 04 (imported from its code): {int((df.family=='saved').sum())} saved (6 windows + chr9 70&ndash;100 Mb comparison), {int((df.family=='called').sum())} short called segments of the <code>all</code> source, {int((df.family=='literature').sum())} literature loci (9p21 zone approx., MYC 8q24, 19q13-qter, 5q13-qter), {len(rnd)} random 25-bin windows.</li>"
         "<li><b>Base-called bins</b> (cn_call = 1 in the <code>all</code> karyogram, AT &ge; 20) are used only to estimate &rho;; calls are not truth.</li>"
         "<li><b>Not available.</b> SNP-level counts, phase-set IDs, unique-molecule counts and a matched normal; clonal LOH regions therefore lose ascertainment and imbalance there cannot be measured at sites that were never selected.</li></ul>"),
        ("Method",
         "<ol><li>Check alignment (coordinates, barcode set, AT = HA + HB, integer values) and join by barcode name.</li>"
         "<li>Estimate beta-binomial overdispersion &rho; from the variance of z<sup>2</sup> in base-called bins stratified by AT quintile; use the median over strata.</li>"
         f"<li>For each region compute S, the mean z<sup>2</sup> over bins with AT &ge; {MIN_AT}; phase-agnostic, so inconsistent orientation between bins does not cancel signal.</li>"
         "<li>Model null: simulate beta-binomial (p = 0.5, &rho;) with each region&rsquo;s own AT per bin; threshold = 95th percentile of S. Empirical null: 95th percentile of S over uncalled random windows.</li>"
         "<li>Power: same simulation with p = 0.5 &plusmn; &delta; (random sign per bin). Scenarios: nominal one-copy loss in a fraction f of cells, nominal one-copy gain, and imbalance anchored to the effective |dBAF| measured at the published MYC amplification and half of it.</li>"
         "<li>Effective |dBAF| (&delta;<sub>eff</sub>) from the bias-corrected mean of z<sup>2</sup>; bootstrap 95% interval over bins (500 resamples).</li>"
         "<li>Label from preset rules (below).</li></ol>"
         "<pre>z_j^2      = (HB_j - HA_j)^2 / AT_j           # per bin, AT_j &gt;= 5\nS          = mean_j z_j^2\nVar(z^2)   under BB:  v_j = 1 + (AT_j - 1) * rho\ndelta_eff  = sqrt(max(0, mean_j (z_j^2 - v_j) / (4 AT_j)))\nrho_stratum = (median(z^2)/0.4549 - 1) / (median_AT - 1)    # 0.4549 = median of chi2_1\nloss (fraction f):  delta = f / (2 (2 - f));  gain:  delta = f / (2 (2 + f))\nlabel:  unobserved if bins_AT&gt;=5 &lt; 3;  consistent with imbalance if S &gt; S_null_model_q95 and CI_lo(delta_eff) &gt; 0;\n        insufficient information if power(MYC-anchored) &lt; 0.8;  else no excess imbalance at observed sites</pre>"
         f"<p>Seeds: <code>numpy.random.default_rng(({C.SEED}, 6, region_index))</code> for the bootstrap and the model null; <code>({C.SEED}, 6, region_index, 1)</code> for power simulations; {NSIM:,} simulations per scenario. "
         "There is no &lsquo;against loss&rsquo; label: the attenuation of the alternative is not calibrated (see below).</p>"),
        ("Results", "", 2),
        ("Alignment and units", C.html_table(checks_tab) + C.caption("Per matrix: whether coordinates match the AnnData, the barcode count, set equality and whether column order equals the AnnData order. All sets match and the join is by name.") +
         f"<p>AT = HA + HB for every cell and bin: {checks['AT_equals_HA_plus_HB_all_cells']}; non-negative: {checks['nonnegative']}; integer-valued: {checks['integer_valued']}.</p>", 3),
        ("Overdispersion (bins in base-called segments; not truth)", C.html_table(rho_tab) + C.caption(f"Per AT quintile of base-called bins: number of bins, median AT, variance inflation (median z<sup>2</sup> / 0.4549) and the implied &rho;. A stable &rho; across strata supports a single beta-binomial parameter; rho used = {rho:.4g}. Variance inflation rises with AT because the overdispersion is proportional to AT."), 3),
        ("Null calibration", C.html_table(pd.DataFrame([{"quantity": "uncalled random windows", "value": model_calib["random_uncalled"]},
                                                         {"quantity": "share of uncalled random windows with S above their model 95th percentile", "value": round(model_calib["frac_uncalled_random_S_above_model_q95"], 3)},
                                                         {"quantity": "share of all random windows with S above their model 95th percentile", "value": round(model_calib["frac_all_random_S_above_model_q95"], 3)},
                                                         {"quantity": "empirical 95th percentile of S, uncalled random windows", "value": round(emp95, 2)}])) +
         C.caption("If the balanced model were adequate, about 5% of uncalled random windows would exceed their threshold. A much larger share means allelic structure (or model misfit) is widespread outside depth calls, so the empirical null is more conservative."), 3),
        ("Effective imbalance at the positive controls",
         C.html_table(lit.reset_index()[["region", "bins_AT_ge5", "S_obs", "S_null_model_q95", "delta_eff", "delta_eff_ci_lo", "delta_eff_ci_hi"]]) +
         C.caption("S and effective |dBAF| for the four published loci (pooled, all cells). The nominal clonal-gain value is &delta; = 0.167, which would give a much larger S; the measured effective &delta; is a fraction of it.") +
         f"<div class='card'>Under the beta-binomial model with nominal AT as the number of observations, a clonal one-copy gain at MYC&rsquo;s depth predicts S of about 1 + (n&minus;1)&rho; + 4&delta;<sup>2</sup>n &asymp; {pred_myc:.0f} (n &asymp; {n_mean_myc:,.0f} per bin), but the published MYC amplification shows S = {myc.S_obs:.1f} "
         f"(effective |dBAF| {myc.delta_eff:.3f}). Measured imbalance is therefore strongly attenuated relative to nominal counts: repeated site observations per molecule, within-bin phase inconsistency, genotype errors and mixed cell states are all candidates. "
         "Theoretical power columns are upper bounds and are not used for labels; power is anchored to the effective imbalance observed at MYC and half of it.</div>", 3),
        ("Observability and information in the saved windows", C.html_table(region_df[region_df.family == "saved"][sv_cols]) +
         C.caption("Pooled allele observations per window: total AT, bins with AT &ge; 5 and with AT = 0, and selected-site observations per raw fragment relative to the genome and to the flanks. Low ratios mean fewer selected heterozygous-site observations per fragment (consistent with, but not proof of, LOH-driven loss of ascertainment, homozygosity or technical effects)."), 3),
        ("Statistics and labels", C.html_table(region_df[st_cols]) + C.caption("S, the model 95th percentile, effective |dBAF| with bootstrap 95% interval, and the resulting label for saved windows and literature loci. Look for S above its threshold together with a positive interval lower bound."), 3),
        ("Power", C.html_table(region_df[pw_cols]) + C.caption("Simulated power at each region&rsquo;s own AT profile. <code>power_loss_f*</code> and <code>power_gain_*</code> use nominal one-copy changes (upper bounds); the last two columns use the effective imbalance measured at MYC and half of it, which is the realistic scale. Scenario deltas below.") +
         C.html_table(scen_tab) + C.caption("Scenario values of &delta; = |p &minus; 0.5|."), 3),
        ("Family medians and label counts", C.html_table(atr_show) + C.caption("Median observability, effective imbalance and anchored power by region family (random windows are the calibration set).") +
         C.html_table(lab_pivot) + C.caption("Number of regions per family and direction by evidence label. Random windows show how often the label arises without any called event."), 3),
        ("Figure", C.img(b64, "allelic information") +
         C.caption("Left: observations per fragment relative to the genome (x) against S (y) for random (grey), called (red), saved (black) and literature (blue) regions; dashed line is the 95th percentile of S in uncalled random windows. "
                   "Look for the controls (MYC, 5q13-qter) far above the cloud and for the 9p21 zone far left. Right: simulated power against pooled AT; points at the top are well powered for nominal changes, while a 25%-of-cells loss (red) is less well powered at low AT."), 3),
        ("Called short segments", C.html_table(df[df.family == "called"][["direction", "region", "bins_AT_ge5", "at_per_raw_vs_genome", "S_obs", "S_null_model_q95", "delta_eff", "power_anchored_MYC_delta_eff", "label"]], max_rows=80) +
         C.caption("Short called segments (first 80). Full columns: <code>region_allelic_information.tsv</code>."), 3),
        ("Per-region detail", reg_html + C.caption("One row per saved window and literature locus; the note restates the numbers in plain language and flags where the preset label rule and the S comparison disagree."), 2),
        ("Sensitivity",
         f"<p>The alternatives are anchored to MYC&rsquo;s effective imbalance and half of it ({SCEN['anchored_MYC_delta_eff']:.3f} and {SCEN['anchored_half_MYC']:.3f}); nominal scenarios give an upper bound. &rho; was estimated per AT stratum ({rho_tab.rho_hat.min():.4f}&ndash;{rho_tab.rho_hat.max():.4f}) and the median used. "
         f"Not varied: the AT &ge; {MIN_AT} cut, the {NSIM:,} simulations, the 500-resample bootstrap and the 0.8 power cut in the label rule.</p>"),
        ("Limits and alternative explanations",
         "<ul><li>AT sums SNP-site DP; one fragment can count at several sites, which &rho; partly absorbs. Unique-molecule counting and SNP-level phase blocks are not available in this package.</li>"
         "<li>Power assumes the region&rsquo;s sites remain ascertained. In clonal LOH without a matched normal, heterozygous sites tend not to be called, so low AT per fragment is the expected signature and imbalance cannot be measured at sites that were never selected.</li>"
         "<li>Labels describe pooled all-cell evidence; a subclonal event in one cluster is diluted.</li>"
         "<li>The label rule depends on a bootstrap lower bound over few bins; MYC (15 informative bins) is labelled &lsquo;no excess imbalance&rsquo; despite S well above the model threshold, so labels near the bin-count limit are fragile.</li>"
         "<li>&lsquo;No excess imbalance&rsquo; is not &lsquo;balanced&rsquo;: absence of detectable imbalance with power only for MYC-sized effects.</li>"
         f"<li>{n_unc_above} of {model_calib['random_uncalled']} uncalled random windows exceed the balanced model, so the null is not calibrated for per-window decisions.</li></ul>"),
    ]
    glossary = {
        "HA / HB / AT": "Haplotype A and B counts per bin and cell; AT = HA + HB, selected-site observation count (sum of DP over selected phased heterozygous SNP sites).",
        "AT_sum / raw_sum": "Pooled AT / pooled raw fragments in the region over all cells.",
        "bins_AT_ge5 / bins_AT_0": "Bins in the target with pooled AT &ge; 5 / equal to 0.",
        "at_per_raw_vs_genome / at_per_raw_vs_flank": "Pooled AT per raw fragment in the region relative to the genome / to its flanks; low values mean fewer selected-site observations per fragment.",
        "S / S_obs": "Mean over bins with AT &ge; 5 of (HB &minus; HA)<sup>2</sup> / AT; sign-agnostic allelic excess variance (1 under a balanced binomial).",
        "S_null_model_q95": "95th percentile of S in simulations from a balanced beta-binomial with &rho; and the region&rsquo;s own AT profile.",
        "S_null_empirical_random_q95": "95th percentile of S over uncalled random windows.",
        "rho / rho_hat": "Beta-binomial overdispersion; per-stratum estimate from base-called bins.",
        "variance_inflation": "Median z<sup>2</sup> divided by 0.4549, in base-called bins of an AT stratum.",
        "delta_eff (|dBAF|)": "Effective absolute deviation of the B-allele fraction from 0.5 after removing the overdispersion contribution.",
        "delta_eff_ci_lo / delta_eff_ci_hi": "Bootstrap 95% interval over bins of delta_eff.",
        "power_loss_f* / power_gain_f*": "Simulated power to exceed the model threshold for a nominal one-copy loss/gain in a fraction f of cells (upper bound).",
        "power_anchored_MYC_delta_eff / power_anchored_half_MYC": "Power for an imbalance equal to / half of the effective |dBAF| measured at MYC.",
        "label": "Evidence label: allele unobserved; consistent with imbalance; insufficient information; no excess imbalance at observed sites.",
        "frac_bins_called_loss_all / frac_bins_called_gain_all": "Share of the region&rsquo;s bins called loss / gain in the <code>all</code> karyogram (not truth).",
        "median_*": "Median of the named column over regions of the family.",
        "AT_stratum / bins / median_AT": "AT quintile of base-called bins, their number and median AT.",
        "regions": "Number of regions with the label.",
        "n / family / direction": "Number of regions; region family; expected direction.",
    }
    C.write_html_report(HERE / "report.html", "Allelic information audit", "Experiment 06 &middot; observability, power and evidence labels for fixed regions",
                        sections, {"command": "python -I experiments/06_allelic_information/run.py", "summary": "summary.json", "seed": f"default_rng(({C.SEED}, 6, region_index))", "n_sim": NSIM, "environment": C.environment()},
                        key_findings=findings, verdict=verdict, glossary=glossary,
                        supports=["Allele observations exist in the chr7 saved windows and show no excess allelic variance, with enough power to see a MYC-sized imbalance but usually not half of it.",
                                  "The 9p21 zone has the lowest selected-site observations per fragment, as expected when clonal LOH removes heterozygous sites.",
                                  "MYC and 5q13-qter carry strong allelic excess; nominal counts overstate allelic information, so likelihoods that treat AT as independent observations would be miscalibrated."],
                        not_supports=["That the chr7 windows are balanced or diploid: only that no excess imbalance is detectable at observed sites, at limited power for weaker effects.",
                                      "Calibrated per-region allelic calls: a large share of uncalled random windows exceed the balanced model.",
                                      "An &lsquo;against clonal loss&rsquo; statement: the attenuation of the alternative is not calibrated."],
                        decision_impact="Supports recommendation ranks 1 and 5: carry an empirically calibrated allelic label (unobserved / insufficient / no excess / consistent with imbalance) in the evidence layer, and defer any joint depth&ndash;allele likelihood until SNP-level, phase-block-aware, molecule-aware counts exist. "
                                        "The chr7 result does not support a clonal one-copy loss as the sole explanation of the dips but cannot exclude weaker or subclonal changes.",
                        next_steps=["Obtain SNP-level phased counts with phase-set IDs and unique-molecule counting to calibrate an allelic likelihood.",
                                    "Revisit the label rule for regions with few informative bins (e.g. MYC) once an empirical null for the bootstrap bound exists.",
                                    "Add a matched normal or panel of unrelated samples so LOH-driven ascertainment can be distinguished from technical loss of sites."])
    pd.set_option("display.width", 250)
    print(checks, rho, model_calib)
    print(rho_tab.to_string())
    print(atr.to_string())
    print(fam.to_string())
    print(df[df.family.isin(["saved", "literature"])][show].to_string())


if __name__ == "__main__":
    main()
