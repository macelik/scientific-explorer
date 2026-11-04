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
    frag = [
        ("Question", "<p>For each fixed region, are allele observations present, how many selected-site observations come per "
         "fragment, could a specified imbalance be detected, and which evidence label is supported? Missing alleles are kept "
         "distinct from balance.</p>"),
        ("Alignment checks", "<pre>" + str(checks) + "</pre>"),
        ("Overdispersion (bins in base-called segments; not truth)", C.html_table(rho_tab) + f"<p>rho used: {rho:.4g}</p>"),
        ("Null calibration", "<pre>" + str(model_calib) + f"</pre><p>Empirical null: 95th percentile of S over uncalled random windows = {emp95:.3g}.</p>"),
        ("Family medians", C.html_table(atr)),
        ("Label counts", C.html_table(fam, max_rows=60)),
        ("Saved windows", C.html_table(df[df.family == "saved"][show])),
        ("Literature loci", C.html_table(df[df.family == "literature"][show])),
        ("Figure", C.img(b64, "allelic information")),
        ("Called short segments", C.html_table(df[df.family == "called"][["direction"] + show], max_rows=80)),
        ("Why there is no 'against loss' label", "<div class='card'>Under the beta-binomial model with nominal AT as the number of "
         "observations, a clonal one-copy gain at MYC's depth predicts S of about 1 + (n−1)ρ + 4δ²n ≈ 146, but the published MYC "
         f"amplification shows S = {df.loc[df.region == 'lit_MYC_8q24', 'S_obs'].iloc[0]:.1f} (effective |dBAF| "
         f"{ctrl['lit_MYC_8q24']:.3f}; 5q13-qter {ctrl['lit_5q13_qter']:.3f}). Measured imbalance is therefore strongly attenuated "
         "relative to nominal counts: repeated site observations per molecule, within-bin phase inconsistency, genotype errors and "
         "mixed cell states are all candidates. Theoretical power columns (loss_f*, gain_f*) are upper bounds and are not used for "
         "labels. Power is instead anchored to the effective imbalance observed at MYC and half of it.</div>"),
        ("Interpretation", "<ul><li><b>chr7 saved windows</b>: allele observations exist (13–23 bins with AT ≥ 5) and show no excess "
         "allelic variance (effective |dBAF| upper 95% bound ≤ about 0.05). An imbalance as strong as the one measured at MYC would "
         "be detected (power ≥ 0.85); one half as strong would usually not (power 0.23–0.77). Selected-site observations per fragment "
         "are 0.44–0.76× the flanks.</li>"
         "<li><b>9p21</b>: the published hemizygous zone has the lowest observations per fragment (0.20× genome), as expected when "
         "heterozygous sites in a clonal LOH region are not ascertained, and the sites that remain are labelled consistent with "
         "imbalance. The saved chr9 windows have a positive effective |dBAF| (bootstrap CI excludes 0) but do not pass the model threshold.</li>"
         "<li><b>Gain controls</b>: MYC and 5q13-qter show strong excess allelic variance. 19q13-qter, which processed X does not call as a "
         "gain (experiment 04), is consistent with imbalance; that is allelic support for a change the depth track loses.</li>"
         "<li>About 31% of uncalled random windows exceed the balanced-model threshold, so allelic structure (or model misfit) is widespread "
         "outside depth calls. A joint method has to calibrate its allelic null empirically, not from nominal counts.</li></ul>"),
        ("Limits", "<ul><li>AT sums SNP-site DP; one fragment can count at several sites, which rho partly absorbs. Unique-molecule "
         "counting and SNP-level phase blocks are not available in this package.</li><li>Power assumes the region's sites remain "
         "ascertained. In clonal LOH without a matched normal, heterozygous sites tend not to be called, so low AT per fragment is the "
         "expected signature and imbalance cannot be measured at sites that were never selected.</li><li>Labels describe pooled "
         "all-cell evidence; a subclonal event in one cluster is diluted.</li></ul>"),
    ]
    C.write_html_report(HERE / "report.html", "Allelic information audit", "Experiment 06 · observability, power and evidence labels for fixed regions",
                        frag, {"command": "python -I experiments/06_allelic_information/run.py", "summary": "summary.json", "environment": C.environment()})
    pd.set_option("display.width", 250)
    print(checks, rho, model_calib)
    print(rho_tab.to_string())
    print(atr.to_string())
    print(fam.to_string())
    print(df[df.family.isin(["saved", "literature"])][show].to_string())


if __name__ == "__main__":
    main()
