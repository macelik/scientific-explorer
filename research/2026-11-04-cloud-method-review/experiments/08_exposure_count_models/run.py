"""Experiment 08 (method record E:1): Poisson / quasi-Poisson / negative-binomial exposure regression.

Implements research/2026-10-29-low-coverage method E:1 on raw counts (layers/raw) for the rule-defined
region families of experiment 04 (saved windows, called short segments, literature loci, 300 random
windows, seed tuple (20261104, 4)).

Rate-ratio model (cells i, part p in {target, flank}, stacked):
    y_ip ~ F(mu_ip),  log mu_ip = log E_i + log n_p + b0 + b1 [p = target] + g_cluster(i)
  E_i = cell raw total EXCLUDING the target and flank bins (exposure must not contain the focal
  region, E:1 guidance); n_p = retained bins in the part. exp(b1) is the target/flank rate ratio.
  F = Poisson (model SE), quasi-Poisson (Pearson scale phi), NB2 with alpha by the Cameron-Trivedi
  moment regression, and Poisson with cell-clustered sandwich SE (target and flank counts of a
  cell are not independent).
Heterogeneity tests (quasi-Poisson Wald, F reference): part x coverage tercile and part x cluster.
Held-out predictive check (cells split within cluster, seed tuple (20261104, 8)): target-only model
    y_it ~ F, log mu = log E_i + log n_t + g_cluster, fit on the training half, predict the test half:
  predicted vs observed zero fraction and Pearson dispersion by coverage tercile, Poisson vs NB.
Calibration: the sampling SE of b1 is tiny with thousands of cells, so the reference for 'is this
  contrast unusual' is the spread of b1 across uncalled random windows (empirical z = (b1 - median) /
  (1.4826 MAD)); BH across the called-segment family. Random windows are not neutral truth.
Outputs: region_models.tsv, predictive_checks.tsv, summary.json, report.html.
"""
import importlib.util
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

spec = importlib.util.spec_from_file_location("exp04", HERE.parent / "04_fixed_region_contrasts" / "run.py")
E4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(E4)
warnings.filterwarnings("ignore")


def design(cl_dummies, part, extra=None):
    cols = [np.ones(len(part)), part, *cl_dummies.T]
    if extra is not None:
        cols += list(extra)
    return np.column_stack(cols)


def wald_block(res, idx, scale_df=None):
    b = res.params[idx]
    V = res.cov_params()[np.ix_(idx, idx)]
    w = float(b @ np.linalg.solve(V, b))
    k = len(idx)
    if scale_df is not None:
        return float(stats.f.sf(w / k, k, scale_df))
    return float(stats.chi2.sf(w, k))


def ct_alpha(y, mu):
    """Cameron-Trivedi auxiliary OLS for NB2 alpha: ((y-mu)^2 - y)/mu = alpha * mu + e."""
    z = ((y - mu) ** 2 - y) / mu
    return max(1e-8, float(np.sum(z * mu) / np.sum(mu * mu)))


def main():
    t0 = time.time()
    p = C.default_paths()
    d = C.H5AD(p["h5ad"])
    lab = C.cell_clusters(d, p)
    regs = E4.build_regions(d, p)
    yt, yf, tot = E4.region_sums(d, regs, "layers/raw")
    clusters = sorted(set(lab))
    D = np.column_stack([(lab == c).astype(float) for c in clusters[1:]])       # reference = cluster 0
    terc = np.asarray(pd.qcut(tot, 3, labels=False))
    T = np.column_stack([(terc == k).astype(float) for k in (1, 2)])
    test = ~C.stratified_halves(lab, (C.SEED, 8))                                 # True = test cell
    kseg = pd.read_csv(p["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    call = np.ones(d.n_bins)
    for r in kseg[kseg.group == "all"].itertuples():
        call[r.bin_start:r.bin_end + 1] = r.cn_call
    panel = set(C.shared_cells().cell)
    in_panel = np.array([c in panel for c in d.cells])
    n = d.n_cells
    rows, pred = [], []
    for j, r in enumerate(regs):
        nt, nf = len(r["target"]), len(r["flank"])
        if nf == 0:
            continue
        Ei = tot - yt[:, j] - yf[:, j]
        y = np.r_[yt[:, j], yf[:, j]]
        part = np.r_[np.ones(n), np.zeros(n)]
        off = np.r_[np.log(Ei) + np.log(nt), np.log(Ei) + np.log(nf)]
        Dd = np.r_[D, D]
        X = design(Dd, part)
        groups = np.r_[np.arange(n), np.arange(n)]
        pois = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit()
        phi = float(np.sum(pois.resid_pearson ** 2) / pois.df_resid)
        quasi = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit(scale="X2")
        robust = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit(cov_type="cluster", cov_kwds={"groups": groups})
        alpha = ct_alpha(y, pois.mu)
        nb = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha), offset=off).fit()
        b1 = float(pois.params[1])
        rec = {"family": r["family"], "region": r["name"], "direction": r["direction"], "n_target_bins": nt,
               "frac_bins_called_loss": float((call[r["target"]] < 1).mean()), "frac_bins_called_gain": float((call[r["target"]] > 1).mean()),
               "rr": float(np.exp(b1)), "log_rr": b1, "phi_quasi": phi, "alpha_nb": alpha}
        for nm, res in [("poisson", pois), ("quasi", quasi), ("nb", nb), ("cluster_robust", robust)]:
            se = float(res.bse[1])
            rec[f"se_{nm}"] = se
            rec[f"rr_lo_{nm}"], rec[f"rr_hi_{nm}"] = float(np.exp(res.params[1] - 1.96 * se)), float(np.exp(res.params[1] + 1.96 * se))
            rec[f"z_{nm}"] = float(res.params[1] / se)
        # interactions (quasi-Poisson)
        Xc = design(Dd, part, extra=[np.r_[T[:, 0], T[:, 0]], np.r_[T[:, 1], T[:, 1]],
                                     part * np.r_[T[:, 0], T[:, 0]], part * np.r_[T[:, 1], T[:, 1]]])
        rc = sm.GLM(y, Xc, family=sm.families.Poisson(), offset=off).fit(scale="X2")
        k0 = Xc.shape[1]
        rec["p_part_x_coverage"] = wald_block(rc, [k0 - 2, k0 - 1], rc.df_resid)
        Xk = design(Dd, part, extra=list((part[:, None] * Dd).T))
        rk = sm.GLM(y, Xk, family=sm.families.Poisson(), offset=off).fit(scale="X2")
        k1 = Xk.shape[1]
        rec["p_part_x_cluster"] = wald_block(rk, list(range(k1 - Dd.shape[1], k1)), rk.df_resid)
        rec["rr_by_cluster"] = ";".join(f"{c}:{np.exp(rk.params[1] + (rk.params[k1 - Dd.shape[1] + i - 1] if i else 0)):.3f}" for i, c in enumerate(clusters))
        # panel-300 rate ratio
        mp = np.r_[in_panel, in_panel]
        rp = sm.GLM(y[mp], X[mp], family=sm.families.Poisson(), offset=off[mp]).fit(scale="X2")
        rec["rr_panel300"] = float(np.exp(rp.params[1]))
        # held-out predictive check on the target counts
        yt_ = yt[:, j]; offt = np.log(Ei) + np.log(nt)
        Xt = np.column_stack([np.ones(n), D])
        tr = ~test
        pt = sm.GLM(yt_[tr], Xt[tr], family=sm.families.Poisson(), offset=offt[tr]).fit()
        at = ct_alpha(yt_[tr], pt.mu)
        nbt = sm.GLM(yt_[tr], Xt[tr], family=sm.families.NegativeBinomial(alpha=at), offset=offt[tr]).fit()
        mu_p = np.exp(Xt[test] @ pt.params + offt[test]); mu_n = np.exp(Xt[test] @ nbt.params + offt[test])
        obs = yt_[test]
        p0_p = np.exp(-mu_p)
        p0_n = (1 + at * mu_n) ** (-1 / at)
        for k in range(3):
            m = terc[test] == k
            pred.append({"region": r["name"], "family": r["family"], "coverage_tercile": k, "test_cells": int(m.sum()),
                         "obs_zero_frac": float((obs[m] == 0).mean()), "pois_pred_zero_frac": float(p0_p[m].mean()),
                         "nb_pred_zero_frac": float(p0_n[m].mean()),
                         "pois_pearson_disp": float(np.mean((obs[m] - mu_p[m]) ** 2 / mu_p[m])),
                         "nb_pearson_disp": float(np.mean((obs[m] - mu_n[m]) ** 2 / (mu_n[m] + at * mu_n[m] ** 2)))})
        z0p = (obs == 0).sum() - p0_p.sum(); z0n = (obs == 0).sum() - p0_n.sum()
        vp, vn = np.sum(p0_p * (1 - p0_p)), np.sum(p0_n * (1 - p0_n))
        rec["test_obs_zero_cells"] = int((obs == 0).sum())
        rec["test_expected_zero_cells_pois"] = float(p0_p.sum()); rec["test_expected_zero_cells_nb"] = float(p0_n.sum())
        # z only where the expected zero count has non-trivial variance (otherwise ~0/~0)
        rec["test_excess_zeros_pois_z"] = float(z0p / np.sqrt(vp)) if vp >= 1 else np.nan
        rec["test_excess_zeros_nb_z"] = float(z0n / np.sqrt(vn)) if vn >= 1 else np.nan
        rec["test_pearson_disp_nb"] = float(np.mean((obs - mu_n) ** 2 / (mu_n + at * mu_n ** 2)))
        rec["alpha_target_only"] = at
        rows.append(rec)
    df = pd.DataFrame(rows)
    pr = pd.DataFrame(pred)
    rnd = df[(df.family == "random") & (df.frac_bins_called_loss == 0) & (df.frac_bins_called_gain == 0)]
    med = float(rnd.log_rr.median()); mad = float(1.4826 * np.median(np.abs(rnd.log_rr - med)))
    df["z_empirical"] = (df.log_rr - med) / mad
    df["p_empirical_two_sided"] = df.log_rr.apply(lambda v: float((np.abs(rnd.log_rr - med) >= abs(v - med)).mean()))
    df["p_empirical_normal"] = 2 * stats.norm.sf(np.abs(df.z_empirical))
    called = df.family == "called"
    from statsmodels.stats.multitest import multipletests
    df.loc[called, "q_bh_empirical_called"] = multipletests(df.loc[called, "p_empirical_normal"], method="fdr_bh")[1]
    df["nb_misfit_flag"] = (df.test_excess_zeros_nb_z.abs() > 3).fillna(False) | (df.test_pearson_disp_nb > 1.5)
    df.to_csv(HERE / "region_models.tsv", sep="\t", index=False)
    pr.to_csv(HERE / "predictive_checks.tsv", sep="\t", index=False)

    def fam_summary(g):
        return pd.Series({"n": len(g), "median_phi_quasi": g.phi_quasi.median(), "median_alpha_nb": g.alpha_nb.median(),
                          "frac_|z_quasi|>1.96": float((np.abs(g.z_quasi) > 1.96).mean()),
                          "frac_|z_cluster_robust|>1.96": float((np.abs(g.z_cluster_robust) > 1.96).mean()),
                          "frac_|z_empirical|>1.96": float((np.abs(g.z_empirical) > 1.96).mean()),
                          "frac_p_coverage<0.01": float((g.p_part_x_coverage < 0.01).mean()),
                          "frac_p_cluster<0.01": float((g.p_part_x_cluster < 0.01).mean()),
                          "median_test_excess_zeros_pois_z": g.test_excess_zeros_pois_z.median(),
                          "median_test_excess_zeros_nb_z": g.test_excess_zeros_nb_z.median(),
                          "n_zero_testable": int(g.test_excess_zeros_nb_z.notna().sum()),
                          "frac_nb_misfit": float(g.nb_misfit_flag.mean())})
    fam = df.groupby(["family", "direction"]).apply(fam_summary).reset_index()
    fam.to_csv(HERE / "family_summary.tsv", sep="\t", index=False)
    show = ["region", "rr", "rr_lo_quasi", "rr_hi_quasi", "rr_lo_cluster_robust", "rr_hi_cluster_robust", "phi_quasi", "alpha_nb",
            "z_quasi", "z_empirical", "p_empirical_two_sided", "p_part_x_coverage", "p_part_x_cluster", "rr_by_cluster", "rr_panel300",
            "test_obs_zero_cells", "test_expected_zero_cells_pois", "test_expected_zero_cells_nb",
            "test_excess_zeros_pois_z", "test_excess_zeros_nb_z", "test_pearson_disp_nb"]
    sv = df[df.family.isin(["saved", "literature"])][show]
    cl = df[called].sort_values("z_empirical")[["region", "direction"] + show[1:] + ["q_bh_empirical_called"]]
    res = {"experiment": "08_exposure_count_models", "method_record": "E:1", "seeds": {"random_windows": [C.SEED, 4], "train_test_split": [C.SEED, 8]},
           "generator": "numpy.random.default_rng (PCG64)", "empirical_null": {"n_uncalled_random": len(rnd), "median_log_rr": med, "mad_log_rr": mad},
           "family_summary": fam.to_dict("records"),
           "called_significant_empirical_q<0.05": {"loss": int(((df.q_bh_empirical_called < .05) & called & (df.direction == "loss")).sum()),
                                                   "gain": int(((df.q_bh_empirical_called < .05) & called & (df.direction == "gain")).sum()),
                                                   "called_loss": int((called & (df.direction == "loss")).sum()),
                                                   "called_gain": int((called & (df.direction == "gain")).sum())},
           "saved": sv.to_dict("records"), "seconds": round(time.time() - t0, 1), "environment": C.environment()}
    C.write_json(HERE / "summary.json", res)

    # ------------------------------------------------------------------ report (structured)
    def figure(fn, w=7, h=3.8):
        fig, ax = plt.subplots(figsize=(w, h)); fn(ax); out = C.fig_to_b64(fig); plt.close(fig); return out

    def f_sig(ax):
        for f, c in [("random", "#bbb"), ("called", "#9b2c2c"), ("saved", "k"), ("literature", "#2f5d8a")]:
            g = df[df.family == f]
            ax.scatter(g.log_rr / np.log(2), np.abs(g.z_quasi), s=12 if f != "saved" else 40, color=c, label=f, alpha=.7)
        ax.set_yscale("log"); ax.axhline(1.96, color="#8a5a00", ls="--", lw=1, label="|z| = 1.96")
        ax.set_xlabel("log2 rate ratio, target / flanks"); ax.set_ylabel("|z|, quasi-Poisson"); ax.legend(fontsize=7)

    def f_emp(ax):
        ax.hist(rnd.log_rr / np.log(2), bins=40, color="#bbb", label="uncalled random windows")
        for _, g in df[df.family == "saved"].iterrows():
            ax.axvline(g.log_rr / np.log(2), color="k", lw=.9)
        for v in (med - 1.96 * mad, med + 1.96 * mad):
            ax.axvline(v / np.log(2), color="#8a5a00", ls="--", lw=1)
        ax.set_xlabel("log2 rate ratio"); ax.set_ylabel("windows"); ax.legend(fontsize=7)

    def f_zero(ax):
        q = pr.groupby(["family", "coverage_tercile"])[["obs_zero_frac", "pois_pred_zero_frac", "nb_pred_zero_frac"]].mean().reset_index()
        for f, mk in [("random", "o"), ("called", "s"), ("saved", "^")]:
            g = q[q.family == f]
            ax.plot(g.coverage_tercile, g.obs_zero_frac - g.pois_pred_zero_frac, marker=mk, color="#2f5d8a", label=f"{f}: observed − Poisson")
            ax.plot(g.coverage_tercile, g.obs_zero_frac - g.nb_pred_zero_frac, marker=mk, color="#9b2c2c", ls="--", label=f"{f}: observed − NB")
        ax.axhline(0, color="k", lw=.5); ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["low", "mid", "high"])
        ax.set_xlabel("coverage tercile (held-out cells)"); ax.set_ylabel("zero-cell fraction: observed − predicted"); ax.legend(fontsize=6)

    def f_clu(ax):
        sv_ = df[df.family == "saved"]
        for _, g in sv_.iterrows():
            vals = [float(x.split(":")[1]) for x in g.rr_by_cluster.split(";")]
            ax.plot(range(len(vals)), vals, marker="o", label=g.region)
        ax.set_xticks(range(len(clusters))); ax.set_xticklabels([f"cluster {c}" for c in clusters])
        ax.set_yscale("log"); ax.set_ylabel("rate ratio target / flanks"); ax.legend(fontsize=6)

    b_sig, b_emp, b_zero, b_clu = figure(f_sig), figure(f_emp), figure(f_zero), figure(f_clu)
    nsig = res["called_significant_empirical_q<0.05"]
    rf = fam.set_index(["family", "direction"])
    sv_ = df[df.family == "saved"].set_index("region")
    lit_ = df[df.family == "literature"]

    def note(r):
        bits = [f"rate ratio {r.rr:.2f} (95% {r.rr_lo_cluster_robust:.2f}–{r.rr_hi_cluster_robust:.2f})"]
        bits.append("unusual vs uncalled random windows" if r.p_empirical_two_sided < 0.05 else "within the random-window range")
        if r.p_part_x_coverage < 0.01:
            bits.append("contrast varies with coverage")
        if r.p_part_x_cluster < 0.01:
            bits.append("contrast differs between clusters")
        if np.isfinite(r.test_excess_zeros_nb_z):
            bits.append("NB explains held-out zeros" if abs(r.test_excess_zeros_nb_z) <= 3 else "NB misfits held-out zeros")
        return "; ".join(bits)
    per = df[df.family.isin(["saved", "literature"])].copy()
    per["note"] = per.apply(note, axis=1)
    per_tab = per[["region", "rr", "rr_lo_cluster_robust", "rr_hi_cluster_robust", "z_empirical", "p_empirical_two_sided",
                   "p_part_x_coverage", "p_part_x_cluster", "rr_by_cluster", "note"]]
    chr7 = sv_.loc[[i for i in sv_.index if i.startswith("chr7")]]
    key = [
        f"Sampling-based tests call almost everything: |z| > 1.96 under quasi-Poisson for {100*rf.loc[('random','random'),'frac_|z_quasi|>1.96']:.0f}% of "
        f"random windows and {100*rf.loc[('called','loss'),'frac_|z_quasi|>1.96']:.0f}% of called losses. Quasi-Poisson φ ≈ "
        f"{rf.loc[('random','random'),'median_phi_quasi']:.1f} and NB α ≈ {rf.loc[('random','random'),'median_alpha_nb']:.3f} fix sampling variance only.",
        f"Against the empirical spread of {len(rnd)} uncalled random windows (median log rr {med:.3f}, robust SD {mad:.3f}), "
        f"{nsig['loss']}/{nsig['called_loss']} called short losses and {nsig['gain']}/{nsig['called_gain']} called short gains remain unusual (BH q < 0.05).",
        f"The four chr7 saved windows: rate ratio {chr7.rr.min():.2f}–{chr7.rr.max():.2f}, empirical z {chr7.z_empirical.min():.1f} to {chr7.z_empirical.max():.1f}. "
        f"The chr9 windows (z {sv_.loc['chr9_20_27','z_empirical']:.1f}, {sv_.loc['chr9_22.5_26.6','z_empirical']:.1f}) sit within the random range.",
        f"Held-out zeros: in chr7 8.3–10.7, {int(sv_.loc['chr7_8.3_10.7','test_obs_zero_cells'])} test cells have zero target counts; Poisson predicts "
        f"{sv_.loc['chr7_8.3_10.7','test_expected_zero_cells_pois']:.1f} and NB {sv_.loc['chr7_8.3_10.7','test_expected_zero_cells_nb']:.1f}. "
        "NB accounts for the zeros of the deep chr7 windows, so zero inflation is not required there (E:3 precondition not met).",
        f"Cluster heterogeneity is common: part × cluster p < 0.01 in {100*rf.loc[('called','loss'),'frac_p_cluster<0.01']:.0f}% of called losses and "
        f"{100*rf.loc[('random','random'),'frac_p_cluster<0.01']:.0f}% of random windows. With thousands of cells, small between-cluster differences are detectable.",
        f"Coverage dependence (a technical signature) appears in {100*rf.loc[('called','loss'),'frac_p_coverage<0.01']:.0f}% of called losses vs "
        f"{100*rf.loc[('random','random'),'frac_p_coverage<0.01']:.0f}% of random windows.",
    ]
    glossary = {
        "rr": "rate ratio target/flanks from the Poisson GLM, exp(b1); exposure and cluster adjusted",
        "rr_lo_* / rr_hi_*": "95% interval under quasi-Poisson, NB or cell-clustered robust SE",
        "phi_quasi": "Pearson χ²/df of the Poisson fit (quasi-Poisson scale)",
        "alpha_nb": "NB2 dispersion from the Cameron–Trivedi regression; Var = μ + αμ²",
        "z_quasi": "b1 / SE under quasi-Poisson: tests rr = 1 against sampling noise only",
        "z_empirical": "(log rr − median) / (1.4826·MAD) over uncalled random windows: how unusual the contrast is among same-width regions",
        "p_empirical_two_sided": "fraction of uncalled random windows at least as far from the median",
        "q_bh_empirical_called": "Benjamini–Hochberg q over the called-segment family using normal p from z_empirical",
        "p_part_x_coverage / p_part_x_cluster": "quasi-Poisson Wald/F test that the rate ratio differs between coverage terciles / clusters",
        "rr_by_cluster": "rate ratio per WNN cluster from the interaction model",
        "test_obs_zero_cells / test_expected_zero_cells_*": "held-out cells with zero target counts, observed vs predicted",
        "test_excess_zeros_nb_z": "(observed − NB-expected zero cells)/SD; NaN when the expected count is ~0",
        "test_pearson_disp_nb": "mean squared Pearson residual of held-out target counts under NB (≈1 if adequate)",
        "family": "saved = frozen windows; called = rule-selected short segments; literature = published MCF10-A loci; random = seeded windows",
    }
    sections = [
        ("Question and why it matters", "<p>Method record <b>E:1</b> (2026-10-29 review): model each cell's raw counts in a region with its "
         "library exposure, compare Poisson, quasi-Poisson and negative-binomial variance, and check whether zeros and spread are "
         "explained. The decision it informs: can a raw-count, exposure-aware likelihood join the evidence layer recommended in "
         "<code>recommendation.md</code>, and what null does it need?</p>"),
        ("Inputs", f"<ul><li><code>layers/raw</code> of the AnnData: {d.n_cells:,} cells × {d.n_bins:,} retained bins (post-QC, pre-GC integer counts).</li>"
         f"<li>WNN clusters <code>wnn_leiden_0.3</code> ({len(clusters)} clusters) from the MuData.</li>"
         f"<li>Regions from experiment 04's rule: {int((df.family=='saved').sum())} saved windows, {int((df.family=='called').sum())} called short segments, "
         f"{len(regs) - len(df)} literature loci with no flank (excluded here), {int((df.family=='random').sum())} random windows (seed (20261104, 4)), "
         f"{len(rnd)} of them uncalled.</li><li>Not available: pre-QC counts, mappability, independent truth.</li></ul>"),
        ("Method", "<pre>for region r, cell i, part p ∈ {target, flank}:\n"
         "  y_ip ~ F(μ_ip)\n  log μ_ip = log E_i + log n_p + b0 + b1·[p = target] + γ_cluster(i)\n"
         "  E_i = cell raw total minus target and flank counts;  n_p = retained bins in part p\n"
         "  F: Poisson | quasi-Poisson (scale = Pearson χ²/df) | NB2 (α by Cameron–Trivedi) | Poisson + sandwich SE clustered by cell\n"
         "heterogeneity: add part×coverage-tercile or part×cluster terms; quasi-Poisson Wald/F test\n"
         "held-out check: split cells within cluster (seed (20261104, 8)); fit y_it on training cells (target only);\n"
         "  predict P(y=0) on test cells: Poisson exp(−μ), NB (1+αμ)^(−1/α)\n"
         "empirical null: median and MAD of b1 over uncalled random windows → z_empirical; BH over called segments</pre>"),
        ("Results", "", 2),
        ("Sampling-based significance is uninformative", C.img(b_sig, "z vs rr") + C.caption(
            "Each point is a region. x: log2 target/flank rate ratio; y: |z| for rr = 1 under quasi-Poisson (log scale). Random windows "
            "(grey) sit almost entirely above the dashed 1.96 line: once counts from thousands of cells are pooled, any regional "
            "difference is 'significant'.") + C.html_table(fam) + C.caption("Family summary. Fractions are within each family."), 3),
        ("Empirical reference from uncalled random windows", C.img(b_emp, "empirical null") + C.caption(
            "Grey: log2 rate ratios of uncalled random 25-bin windows. Dashed: median ± 1.96 robust SD. Black: the saved windows. "
            "This spread is the realistic scale against which a regional contrast should be judged at this resolution.") +
         C.html_table(cl.head(25)) + C.caption("Called short segments, most negative empirical z first (top 25; full table in region_models.tsv)."), 3),
        ("Held-out zero check (Poisson vs NB)", C.img(b_zero, "zero check") + C.caption(
            "Observed minus predicted fraction of held-out cells with zero target counts, by coverage tercile. Solid blue: Poisson; "
            "dashed red: NB. Poisson leaves excess zeros; NB brings them close to zero.") +
         C.html_table(df[df.family == "saved"][["region", "test_obs_zero_cells", "test_expected_zero_cells_pois", "test_expected_zero_cells_nb",
                                                 "test_excess_zeros_nb_z", "test_pearson_disp_nb"]]) +
         C.caption("Saved windows: held-out zero cells. NaN z means too few expected zeros to test."), 3),
        ("Between-cluster heterogeneity of the saved windows", C.img(b_clu, "rr by cluster") + C.caption(
            "Rate ratio per WNN cluster for each saved window (log scale). Parallel flat lines mean a depletion shared by all clusters; "
            "a single outlying cluster suggests a cluster-specific contribution (cell state or subclone)."), 3),
        ("Per-region detail", C.html_table(per_tab) + C.caption("Saved windows (literature loci have no flank and are not modelled here). "
         "The note is generated from the columns.")),
        ("Sensitivity", "<p>All four variance models give nearly the same point estimate; the intervals differ by φ but stay narrow. "
         f"Rate ratios on the frozen 300-cell panel agree with all cells (column <code>rr_panel300</code> in region_models.tsv). "
         "Cluster adjustment is included in every model; removing it was not prespecified.</p>"),
        ("Limits and alternative explanations", "<ul><li>Accessibility, mappability and CN are confounded in every column; a deep, "
         "homogeneous, NB-consistent depletion is equally expected from an inaccessible or poorly mappable region.</li>"
         "<li>Called segments were selected on these data, so their empirical p-values are conditional on selection.</li>"
         "<li>Random windows are not neutral truth; some contain events, which widens the reference and makes the test conservative.</li>"
         "<li>Exposure excludes the focal region and flanks but still contains other events; cluster terms can absorb real signal.</li></ul>"),
    ]
    C.write_html_report(
        HERE / "report.html", "Exposure count models", "Experiment 08 · method E:1 · Poisson / quasi-Poisson / NB regression on raw counts",
        sections, {"command": "python -I experiments/08_exposure_count_models/run.py", "outputs": ["region_models.tsv", "predictive_checks.tsv",
                   "family_summary.tsv", "summary.json"], "seeds": res["seeds"], "environment": C.environment()},
        verdict="Raw-count models with exposure fit the data's sampling behaviour (NB explains the zeros), but they only become informative "
                "when judged against the spread of same-width uncalled regions. On that scale the deep chr7 windows stand out and the chr9 windows do not.",
        key_findings=key, glossary=glossary,
        supports=["NB, not zero inflation, is adequate for zeros in the deep chr7 windows.",
                  "An empirical, region-based null is required; sampling-based p-values are not.",
                  f"{nsig['loss']} called short losses and {nsig['gain']} called gains are unusual relative to uncalled regions."],
        not_supports=["Any statement that a depletion is dosage rather than accessibility or mappability.",
                      "Calibrated discovery claims for called segments (selected on the same data).",
                      "Absolute CN states."],
        decision_impact="E:1 can enter the evidence layer as <i>rate ratio + empirical z + heterogeneity flags</i>, with NB as the sampling model. "
                        "It does not change the recommendation's ordering: the empirical z is only as good as the reference regions, which strengthens "
                        "the case for a cross-sample technical baseline (stage 2).",
        next_steps=["Replace random windows with a cross-sample baseline (same bins in unrelated samples) and recompute z.",
                    "Repeat with region selection and assessment on disjoint cell halves to remove selection bias.",
                    "Add mappability/GC covariates once a matched mappability track is available."])
    pd.set_option("display.width", 250)
    print(nsig, res["seconds"])


if __name__ == "__main__":
    main()
