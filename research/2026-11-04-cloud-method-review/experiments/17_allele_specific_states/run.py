"""Experiment 17: combine depth and allele evidence into allele-specific relative copy-number states.

Modality-specific likelihoods, combined per region (rule-defined families of experiment 04):
  depth   L_obs  = raw-count rate ratio target/flanks (rr_flank_raw, exp 04); for broad literature
                   loci without flanks, target/genome density (rr_genome_raw).
                   log L_obs ~ Normal(log L_exp, sigma_A^2), sigma_A = robust SD (1.4826 MAD) of
                   log L over UNCALLED random windows (exp 04/08): the realistic regional
                   baseline (accessibility/mappability) variability; sampling SE is negligible.
  allele  h_obs  = phase-agnostic major-haplotype fraction (beta-binomial MLE, group 'all',
                   empirical rho; exp 10), folded to >= 0.5, with SE from its profile interval.
                   h_obs ~ Normal(m, se^2 + tau^2),  m = max(0.5 + k (h_exp - 0.5), h_null_median)
                   tau = robust SD of h_obs over uncalled random windows; h_null_median absorbs
                   the folding bias at balance; k = allelic attenuation (not identified here; grid
                   1.0 / 0.6 / 0.4; exp 06/10/13 measured strong attenuation at published gains).
                   Regions with < 3 bins carrying AT contribute no allele likelihood (missing, not balanced).
  cells   (homozygous-loss states only) among cells expected to have >= 5 target fragments
                   (exp 04), zeros ~ Binomial(n, f + (1 - f) p0), p0 = 0.005.
States: event type x cell fraction f in {0.1, 0.2, ..., 1.0}, diploid baseline 1:1 = level 1.
  none (1:1)                    L = 1,             h = 0.5
  loss1  1:0 in f              L = 1 - f/2,       h = 1 / (2 - f)
  homdel 0:0 in f              L = 1 - f,         h = 0.5
  cnloh  2:0 in f              L = 1,             h = (1 + f) / 2
  gain1  2:1 in f              L = 1 + f/2,       h = (1 + f) / (2 + f)
  gain2  2:2 in f              L = 1 + f,         h = 0.5
  amp    4:1 in f              L = 1 + 1.5 f,     h = (1 + 3 f) / (2 + 3 f)
Prior: uniform over event types (0.5 on none to reflect that most regions are not events is a
sensitivity), uniform over f. Posterior per type; MAP state; posterior-predictive fit check:
depth z = (log L_obs - log L_exp)/sigma_A and allele z under the MAP state; |z| > 2.5 for the MAP
flags 'no modelled state fits' (unmodelled baseline effect or state). Also reported: the best
fit if the depth term is dropped (allele-only) and if alleles are dropped (depth-only).
Major/minor relative copies: posterior means of a = major, b = minor copies averaged over cells
(population average under each state), with 90% credible intervals.
Deterministic; no random numbers.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

E = HERE.parent
FR = np.round(np.arange(0.1, 1.01, 0.1), 2)
P0 = 0.005


def states():
    out = [("none", 0.0, 1.0, 0.5, 1.0, 1.0)]
    for f in FR:
        out += [("loss1", f, 1 - f / 2, 1 / (2 - f), 1.0, 1 - f),
                ("homdel", f, 1 - f, 0.5, 1 - f, 1 - f),
                ("cnloh", f, 1.0, (1 + f) / 2, 1 + f, 1 - f),
                ("gain1", f, 1 + f / 2, (1 + f) / (2 + f), 1 + f, 1.0),
                ("gain2", f, 1 + f, 0.5, 1 + f, 1 + f),
                ("amp", f, 1 + 1.5 * f, (1 + 3 * f) / (2 + 3 * f), 1 + 3 * f, 1.0)]
    return pd.DataFrame(out, columns=["type", "f", "L_exp", "h_exp", "major", "minor"])


def score(row, S, sig_A, tau, hmed, k, use_depth=True, use_allele=True, prior_none=None):
    ll = np.zeros(len(S))
    if use_depth and np.isfinite(row.L_obs) and row.L_obs > 0:
        ll += stats.norm.logpdf(np.log(row.L_obs), np.log(np.clip(S.L_exp, 1e-3, None)), sig_A)
    has_allele = use_allele and np.isfinite(row.h_obs) and row.allele_bins >= 3
    if has_allele:
        m = np.maximum(0.5 + k * (S.h_exp - 0.5), hmed)
        ll += stats.norm.logpdf(row.h_obs, m, np.sqrt(row.h_se ** 2 + tau ** 2))
    hd = (S.type == "homdel").to_numpy()
    if np.isfinite(row.cells_lam5) and row.cells_lam5 > 0:
        z = int(round(row.frac_lam5_zero * row.cells_lam5))
        pz = np.where(hd, S.f + (1 - S.f) * P0, P0)
        ll += stats.binom.logpmf(z, int(row.cells_lam5), pz)
    # prior: uniform over types, uniform over f within type
    ntypes = S.type.nunique()
    lp = np.where(S.type == "none", np.log(1 / ntypes), np.log(1 / ntypes / len(FR)))
    if prior_none is not None:
        lp = np.where(S.type == "none", np.log(prior_none), np.log((1 - prior_none) / (ntypes - 1) / len(FR)))
    post = ll + lp
    post = np.exp(post - post.max()); post /= post.sum()
    return post, ll, has_allele


def main():
    p = C.default_paths()
    rc = pd.read_csv(E / "04_fixed_region_contrasts" / "region_contrasts.tsv", sep="\t")
    rc = rc[rc.family != "null"]
    al = pd.read_csv(E / "10_betabinomial_allele_model" / "region_estimates.tsv", sep="\t")
    al = al[(al.group == "all") & (al.rho_setting == "empirical")][["family", "region", "h_agn", "h_agn_lo", "h_agn_hi", "n_bins_AT_ge1", "any_bin_called_all"]]
    df = rc.merge(al, on=["family", "region"], how="left")
    df["L_obs"] = np.where(df.n_flank_bins > 0, df.rr_flank_raw, df.rr_genome_raw)
    df["depth_reference"] = np.where(df.n_flank_bins > 0, "flanks", "genome")
    df["h_obs"] = df.h_agn
    df["h_se"] = (df.h_agn_hi - df.h_agn_lo) / (2 * 1.96)
    df["allele_bins"] = df.n_bins_AT_ge1.fillna(0)
    df["cells_lam5"] = df.cells_lam_ge5
    df["frac_lam5_zero"] = df.frac_cells_lam_ge5_with_zero
    unc = df[(df.family == "random") & (df.frac_bins_called_loss_all == 0) & (df.frac_bins_called_gain_all == 0)]
    mad = lambda v: float(1.4826 * np.median(np.abs(v - np.median(v))))  # noqa: E731
    sig_flank = mad(np.log(unc.rr_flank_raw))
    sig_genome = mad(np.log(unc.rr_genome_raw))
    hmed = float(unc.h_obs.median())
    tau = mad(unc.h_obs)
    S = states()
    rows, typ_rows = [], []
    for k in (1.0, 0.6, 0.4):
        for _, r in df.iterrows():
            sig = sig_flank if r.depth_reference == "flanks" else sig_genome
            post, ll, has_allele = score(r, S, sig, tau, hmed, k)
            i = int(np.argmax(post))
            st = S.iloc[i]
            dz = (np.log(r.L_obs) - np.log(st.L_exp)) / sig if r.L_obs > 0 else np.nan
            az = np.nan
            if has_allele:
                m = max(0.5 + k * (st.h_exp - 0.5), hmed)
                az = (r.h_obs - m) / np.sqrt(r.h_se ** 2 + tau ** 2)
            ptype = pd.Series(post).groupby(S.type).sum()
            a_mean = float(np.sum(post * S.major)); b_mean = float(np.sum(post * S.minor))
            order = np.argsort(S.major.to_numpy())
            cdf = np.cumsum(post[order]); maj = S.major.to_numpy()[order]
            a_lo, a_hi = maj[np.searchsorted(cdf, .05)], maj[min(len(cdf) - 1, np.searchsorted(cdf, .95))]
            order = np.argsort(S.minor.to_numpy()); cdf = np.cumsum(post[order]); mn = S.minor.to_numpy()[order]
            b_lo, b_hi = mn[np.searchsorted(cdf, .05)], mn[min(len(cdf) - 1, np.searchsorted(cdf, .95))]
            p_d, _, _ = score(r, S, sig, tau, hmed, k, use_allele=False)
            p_a, _, ha = score(r, S, sig, tau, hmed, k, use_depth=False)
            p_n, _, _ = score(r, S, sig, tau, hmed, k, prior_none=0.5)
            fit = "fits" if (abs(dz) <= 2.5 and (not np.isfinite(az) or abs(az) <= 2.5)) else "no modelled state fits"
            rec = {"k": k, "family": r.family, "region": r.region, "direction": r.direction, "depth_reference": r.depth_reference,
                   "L_obs": r.L_obs, "h_obs": r.h_obs, "h_se": r.h_se, "allele_used": has_allele,
                   "map_type": st.type, "map_f": st.f, "map_post": float(post[i]), "map_depth_z": dz, "map_allele_z": az, "fit": fit,
                   "major_mean": a_mean, "major_lo90": a_lo, "major_hi90": a_hi, "minor_mean": b_mean, "minor_lo90": b_lo, "minor_hi90": b_hi,
                   "map_type_depth_only": S.type.iloc[int(np.argmax(p_d))], "map_type_allele_only": S.type.iloc[int(np.argmax(p_a))] if ha else "no alleles",
                   "map_type_prior_none_0.5": S.type.iloc[int(np.argmax(p_n))],
                   "frac_bins_called_loss": r.frac_bins_called_loss_all, "frac_bins_called_gain": r.frac_bins_called_gain_all}
            for t, v in ptype.items():
                rec[f"P_{t}"] = float(v)
            rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(HERE / "region_states.tsv", sep="\t", index=False)

    # ---- second allele source: phase-orientation HMM of experiment 16 (k = 0.6), own null
    E16 = E / "16_major_minor_haplotypes"
    hmm_rows, hmm_meta = [], {}
    if (E16 / "region_haplotypes.tsv").exists():
        rh = pd.read_csv(E16 / "region_haplotypes.tsv", sep="\t")
        rh = rh[rh.group == "all"].set_index("region")
        sh = pd.read_csv(E16 / "segment_haplotypes.tsv", sep="\t")
        sh = sh[sh.group == "all"]
        rn = pd.read_csv(E16 / "random_null.tsv", sep="\t")
        rn = rn[(rn.group == "all") & rn.eligible]
        hmed16 = float(rn.h_hat.median()); tau16 = mad(rn.h_hat)
        hmm_meta = {"h_null_median_exp16": hmed16, "tau_h_exp16": tau16, "n_random_exp16": int(len(rn))}
        d16 = df.copy()
        d16["h_obs"] = np.nan; d16["h_se"] = np.nan; d16["allele_bins"] = 0
        for i, r in d16.iterrows():
            src = None
            if r.region in rh.index:
                src = rh.loc[r.region]
            elif r.family == "called":
                ch = r.region.split(":")[0]; a = float(r.region.split(":")[1].split("-")[0]) * 1e6
                m = sh[(sh.chromosome == ch) & ((sh.start_bp - 1 - a).abs() < 60_000)]
                if len(m):
                    src = m.iloc[0]
            if src is not None and bool(src.get("eligible", True)) and np.isfinite(src.h_hat):
                hh = max(float(src.h_hat), 1 - float(src.h_hat))
                lo, hi = src.get("h_ci_lo", np.nan), src.get("h_ci_hi", np.nan)
                se = (hi - lo) / (2 * 1.96) if np.isfinite(lo) and np.isfinite(hi) and hi > lo else 0.01
                d16.loc[i, ["h_obs", "h_se", "allele_bins"]] = [hh, se, float(src.bins_AT_ge5)]
        for _, r in d16[d16.family.isin(["saved", "literature", "called"])].iterrows():
            sig = sig_flank if r.depth_reference == "flanks" else sig_genome
            post, _, has_allele = score(r, S, sig, tau16, hmed16, 0.6)
            i = int(np.argmax(post)); st = S.iloc[i]
            hmm_rows.append({"family": r.family, "region": r.region, "direction": r.direction, "L_obs": r.L_obs, "h_hmm": r.h_obs,
                             "allele_used": has_allele, "map_type_hmm": st.type, "map_f_hmm": st.f, "map_post_hmm": float(post[i]),
                             "major_mean_hmm": float(np.sum(post * S.major)), "minor_mean_hmm": float(np.sum(post * S.minor))})
    hmm = pd.DataFrame(hmm_rows)
    if len(hmm):
        hmm = hmm.merge(out[out.k == 0.6][["region", "h_obs", "map_type", "map_f"]].rename(columns={"h_obs": "h_agnostic", "map_type": "map_type_agnostic",
                                                                                                  "map_f": "map_f_agnostic"}), on="region", how="left")
        hmm["agree"] = hmm.map_type_hmm == hmm.map_type_agnostic
        hmm.to_csv(HERE / "region_states_hmm_allele.tsv", sep="\t", index=False)
    prim = out[out.k == 0.6]

    def fam_tab(g):
        return pd.Series({"n": len(g), "fits": int((g.fit == "fits").sum()),
                          **{f"MAP_{t}": int((g.map_type == t).sum()) for t in ["none", "loss1", "homdel", "cnloh", "gain1", "gain2", "amp"]},
                          "allele_changes_depth_only_call": int((g.map_type != g.map_type_depth_only).sum())})
    fam = prim.groupby(["family", "direction"])[list(prim.columns)].apply(fam_tab).reset_index()
    fam.to_csv(HERE / "family_summary.tsv", sep="\t", index=False)
    ksens = out[out.family.isin(["saved", "literature"])].pivot_table(index="region", columns="k", values="map_type", aggfunc="first")
    ksens.columns = [f"MAP_k{c}" for c in ksens.columns]
    ksens = ksens.reset_index()
    ksens.to_csv(HERE / "k_sensitivity.tsv", sep="\t", index=False)
    res = {"experiment": "17_allele_specific_states", "deterministic": True, "inputs": ["04 region_contrasts.tsv", "10 region_estimates.tsv"],
           "sigma_A_flank": sig_flank, "sigma_A_genome": sig_genome, "h_null_median": hmed, "tau_h": tau, "k_grid": [1.0, 0.6, 0.4],
           "primary_k": 0.6, "cell_fractions": FR.tolist(), "p0_zero": P0,
           "family_summary_k0.6": fam.to_dict("records"), "environment": C.environment()}
    res["hmm_allele_source"] = hmm_meta
    if len(hmm):
        res["hmm_vs_agnostic_agreement_called"] = [int(hmm[hmm.family == "called"].agree.sum()), int((hmm.family == "called").sum())]
    C.write_json(HERE / "summary.json", res)

    # ------------------------------------------------------------------ figures
    def figure(fn, w=7, h=4.2):
        fig, ax = plt.subplots(figsize=(w, h)); fn(ax); b = C.fig_to_b64(fig); plt.close(fig); return b

    def f_plane(ax):
        L = np.exp(np.linspace(np.log(0.08), np.log(3.5), 200))
        for t, c in [("loss1", "#9b2c2c"), ("homdel", "#000"), ("cnloh", "#8a5a00"), ("gain1", "#2f5d8a"), ("gain2", "#5a8a2f"), ("amp", "#7a3fa0")]:
            g = S[S.type == t]
            ax.plot(g.L_exp, 0.5 + 0.6 * (g.h_exp - 0.5), "-", color=c, label=f"{t} (k=0.6)", lw=1.5)
        ax.scatter([1], [0.5], color="k", marker="*", s=80, label="none 1:1")
        for f, mk, c in [("random", "o", "#bbb"), ("called", "s", "#9b2c2c"), ("saved", "^", "k"), ("literature", "D", "#2f5d8a")]:
            g = prim[prim.family == f]
            ax.scatter(g.L_obs, g.h_obs, s=14 if f in ("random", "called") else 45, marker=mk, color=c, alpha=.7, edgecolor="none")
        ax.axvspan(np.exp(-1.96 * sig_flank), np.exp(1.96 * sig_flank), color="#ddd", alpha=.35, zorder=0)
        ax.set_xscale("log"); ax.set_xlabel("depth level L (raw, relative to flanks / genome)"); ax.set_ylabel("major-haplotype fraction ĥ")
        ax.legend(fontsize=6, ncol=2); ax.set_ylim(0.48, 0.85)
    b_plane = figure(f_plane, 8, 5)

    def f_post(ax):
        g = prim[prim.family.isin(["saved", "literature"])].set_index("region")
        cols = [c for c in g.columns if c.startswith("P_")]
        bottom = np.zeros(len(g))
        colors = {"P_none": "#bbb", "P_loss1": "#9b2c2c", "P_homdel": "#000", "P_cnloh": "#8a5a00", "P_gain1": "#2f5d8a", "P_gain2": "#5a8a2f", "P_amp": "#7a3fa0"}
        for c in cols:
            ax.barh(range(len(g)), g[c], left=bottom, color=colors.get(c, "#999"), label=c[2:]); bottom += g[c].to_numpy()
        ax.set_yticks(range(len(g))); ax.set_yticklabels(g.index, fontsize=7); ax.set_xlabel("posterior probability (k = 0.6)"); ax.legend(fontsize=6, ncol=4)
    b_post = figure(f_post, 8, 4.5)

    # ------------------------------------------------------------------ report
    sv = prim[prim.family.isin(["saved", "literature"])]

    def note(r):
        s = f"L = {r.L_obs:.2f}, ĥ = {r.h_obs:.3f}" if np.isfinite(r.h_obs) else f"L = {r.L_obs:.2f}, no alleles"
        s += f"; MAP {r.map_type}" + (f" in {r.map_f:.0%} of cells" if r.map_type != "none" else "")
        s += f" (P = {r.map_post:.2f}); depth z {r.map_depth_z:+.1f}"
        if np.isfinite(r.map_allele_z):
            s += f", allele z {r.map_allele_z:+.1f}"
        if r.fit != "fits":
            s += "; no modelled state fits → regional baseline effect or unmodelled state"
        if r.map_type != r.map_type_depth_only:
            s += f"; alleles change the depth-only call ({r.map_type_depth_only} → {r.map_type})"
        return s
    per = sv.copy(); per["note"] = per.apply(note, axis=1)
    per_tab = per[["region", "L_obs", "h_obs", "map_type", "map_f", "map_post", "major_mean", "major_lo90", "major_hi90", "minor_mean",
                   "minor_lo90", "minor_hi90", "fit", "note"]]
    fs = fam.set_index(["family", "direction"])
    sw = prim[prim.family == "saved"]
    chr7 = sw[sw.region.str.startswith("chr7")]
    lit = prim[prim.family == "literature"].set_index("region")
    def lab(r):
        return r.map_type + ("" if r.map_type == "none" else " f=%.1f" % r.map_f)
    key = [
        f"Baseline variability used as depth uncertainty: σ_A = {sig_flank:.3f} (log scale, flank-relative) and {sig_genome:.3f} (genome-relative), "
        f"from {len(unc)} uncalled random windows. Allelic null: median ĥ {hmed:.3f}, spread τ = {tau:.3f}.",
        f"chr7 saved windows (L {chr7.L_obs.min():.2f}–{chr7.L_obs.max():.2f}, ĥ {chr7.h_obs.min():.3f}–{chr7.h_obs.max():.3f}): MAP states "
        + ", ".join(r.region + " → " + lab(r) for _, r in chr7.iterrows()) + "; " +
        f"{int((chr7.fit != 'fits').sum())} of {len(chr7)} flagged 'no modelled state fits'.",
        f"Called short losses: {int(fs.loc[('called','loss'),'fits'])}/{int(fs.loc[('called','loss'),'n'])} fit a modelled state; "
        f"MAP one-copy loss {int(fs.loc[('called','loss'),'MAP_loss1'])}, homozygous loss {int(fs.loc[('called','loss'),'MAP_homdel'])}, none {int(fs.loc[('called','loss'),'MAP_none'])}.",
        f"Called short gains: {int(fs.loc[('called','gain'),'fits'])}/{int(fs.loc[('called','gain'),'n'])} fit; MAP gain1 {int(fs.loc[('called','gain'),'MAP_gain1'])}, "
        f"gain2 {int(fs.loc[('called','gain'),'MAP_gain2'])}, amp {int(fs.loc[('called','gain'),'MAP_amp'])}, none {int(fs.loc[('called','gain'),'MAP_none'])}.",
        "Literature controls (k = 0.6): " + "; ".join(i.replace("lit_", "") + ": " + lab(r) + " (" + r.fit + ")" for i, r in lit.iterrows()) + ".",
        f"Adding alleles changes the depth-only MAP type in {int(prim.map_type.ne(prim.map_type_depth_only).sum())} of {len(prim)} regions.",
    ]
    glossary = {
        "L_obs": "depth level: raw-count rate ratio target/flanks (or target/genome for literature loci)",
        "h_obs / h_se": "phase-agnostic major-haplotype fraction from exp 10 (≥ 0.5) and its SE",
        "k": "allelic attenuation: observed imbalance = k × expected imbalance (not identified; 1.0/0.6/0.4)",
        "map_type / map_f": "most probable event type and cell fraction (none = 1:1)",
        "P_<type>": "posterior probability of the event type, summed over cell fractions",
        "map_depth_z / map_allele_z": "standardized residual of depth / allele under the MAP state",
        "fit": "'no modelled state fits' when |z| > 2.5 for either modality under the MAP state",
        "major_mean / minor_mean (lo90, hi90)": "posterior mean relative copies of the major / minor haplotype (population average), 90% credible bounds",
        "map_type_depth_only / allele_only": "MAP type when the other modality is dropped",
        "σ_A": "robust SD of log depth level among uncalled random windows (baseline variability)",
        "τ": "robust SD of ĥ among uncalled random windows (allelic null spread)",
        "loss1, homdel, cnloh, gain1, gain2, amp": "1:0, 0:0, 2:0, 2:1, 2:2, 4:1 in a fraction f of cells on a 1:1 background",
    }
    sections = [
        ("Question and why it matters", "<p>The user asked to move from depth to the allelic signal, use the model that suits each modality, "
         "and combine them to identify major and minor haplotypes. This experiment scores explicit allele-specific states against a region's "
         "depth level and major-haplotype fraction, each with its own likelihood and its own empirical null, and reports when no simple state "
         "explains both. It is the prototype for the evidence layer in <code>recommendation.md</code>.</p>"),
        ("Inputs", f"<ul><li>Depth: experiment 04 <code>region_contrasts.tsv</code> (raw counts, 3,460 cells): {len(df)} regions "
         f"({', '.join(f'{k}: {v}' for k, v in df.family.value_counts().items())}).</li>"
         "<li>Alleles: experiment 10 <code>region_estimates.tsv</code> (group all, empirical ρ = 0.0038): phase-agnostic ĥ with profile interval.</li>"
         "<li>Per-cell zeros: experiment 04 (cells expected to have ≥ 5 target fragments that observed 0).</li>"
         "<li>Not available: SNP-level phase sets (ĥ is phase-agnostic), mappability, independent truth; MCF10-A baseline assumed 1:1.</li></ul>"),
        ("Method", "<pre>" + __doc__.split("Modality-specific likelihoods, combined per region (rule-defined families of experiment 04):")[1].split("Deterministic")[0] + "</pre>"),
        ("Results", "", 2),
        ("Where regions fall relative to the modelled states", C.img(b_plane, "L-h plane") + C.caption(
            "x: depth level (log scale); y: major-haplotype fraction. Coloured curves: expected (L, ĥ) for each event type as the cell fraction "
            "goes from 0.1 to 1 (attenuation k = 0.6). Grey band: ±1.96 σ_A of baseline variability. Points: random (grey), called short "
            "segments (red squares), saved windows (black triangles), literature loci (blue diamonds). A region far from every curve is not "
            "explained by these dosage states at this baseline variability."), 3),
        ("Family summary (k = 0.6)", C.html_table(fam) + C.caption("Counts of MAP types and of regions where a modelled state fits."), 3),
        ("Posterior over event types for saved windows and controls", C.img(b_post, "posteriors") + C.caption(
            "Stacked posterior probability of each event type (k = 0.6)."), 3),
        ("Sensitivity to the attenuation k", C.html_table(ksens) + C.caption("MAP type at k = 1.0, 0.6, 0.4. A call that changes with k depends on how "
         "much the allelic signal is attenuated, which these data do not identify."), 3),
        ("Per-region detail", C.html_table(per_tab) + C.caption("k = 0.6. Notes are generated from the columns.")),
        ("Allele source sensitivity: phase-orientation HMM (experiment 16)", (C.html_table(hmm[hmm.family.isin(["saved", "literature"])]) +
         C.caption(f"Same depth term; alleles from experiment 16's HMM-oriented major-haplotype fraction with its own null (median "
                   f"{hmm_meta.get('h_null_median_exp16', float('nan')):.3f}, spread {hmm_meta.get('tau_h_exp16', float('nan')):.3f}, "
                   f"{hmm_meta.get('n_random_exp16', 0)} random windows). Called segments: MAP type agrees with the phase-agnostic source in "
                   f"{int(hmm[hmm.family == 'called'].agree.sum()) if len(hmm) else 0} of {int((hmm.family == 'called').sum()) if len(hmm) else 0}; full table in region_states_hmm_allele.tsv."))
         if len(hmm) else "<p>Experiment 16 outputs not found.</p>"),
        ("Called short segments", C.html_table(prim[prim.family == "called"][["region", "direction", "L_obs", "h_obs", "map_type", "map_f", "map_post",
                                                                                 "major_mean", "minor_mean", "fit", "map_type_depth_only"]], max_rows=80)),
        ("Limits and alternative explanations", "<ul><li>The state grid assumes a 1:1 diploid background and population-average mixtures; real "
         "samples can have other baselines, several subclones or allele-specific accessibility.</li>"
         "<li>k is not identified. Published MCF10-A gains (MYC, 5q) show strong attenuation, but their true allele-specific states are not "
         "known here, so k cannot be calibrated.</li><li>σ_A comes from random windows, which can contain events and are not neutral truth.</li>"
         "<li>Depth and allele observations share fragments; treating their likelihoods as independent overstates certainty somewhat.</li>"
         "<li>ĥ is phase-agnostic and pooled over bins; mirrored events in different cells cancel.</li>"
         "<li><b>Ascertainment:</b> heterozygous sites are selected in the same sample. In a clonal hemizygous region (e.g. the published "
         "9p21 loss) truly lost-heterozygosity sites are not selected, so the remaining sites look balanced and the allele term pulls "
         "toward 'none'. Experiment 06 found the expected signature there (0.20× selected-site observations per fragment). A state model "
         "needs an explicit site-density term, or a matched normal, before balanced alleles can count against a loss.</li></ul>"),
    ]
    C.write_html_report(
        HERE / "report.html", "Allele-specific states", "Experiment 17 · depth + allele combined into major/minor haplotype states",
        sections, {"command": "python -I experiments/17_allele_specific_states/run.py", "outputs": ["region_states.tsv", "family_summary.tsv",
                   "k_sensitivity.tsv", "summary.json"], "environment": C.environment()},
        verdict="Combining depth with alleles separates regions that behave like simple allele-specific dosage changes from regions whose depth "
                "change has no matching allelic pattern; the deep chr7 dips belong to the second group.",
        key_findings=key, glossary=glossary,
        supports=["Reporting major/minor relative copies with credible bounds per region, with the modality contributions shown separately.",
                  "Flagging depth changes that no modelled dosage state explains given the observed alleles.",
                  "Using allele evidence to change or confirm a depth-only call."],
        not_supports=["Absolute copy number or verified haplotype identity (phase is not verified).",
                      "A calibrated attenuation k; calls that change with k remain unresolved.",
                      "Distinguishing accessibility from dosage when alleles are balanced and depth is moderately reduced."],
        decision_impact="This is the concrete form of rank 1 in <code>method_comparison.tsv</code> (evidence layer) extended with rank 5 (joint "
                        "depth–allele evidence): modality-specific likelihoods, empirical nulls for each, and an explicit 'no state fits' outcome. "
                        "To become a method it needs a calibrated k (SNP-level, phase-aware counts) and an external baseline for σ_A.",
        next_steps=["Calibrate which allele source (phase-agnostic or HMM-oriented) to use on samples with DNA-defined allele-specific states.",
                    "Add a selected-site-density likelihood (AT per fragment vs genome) so that LOH-driven ascertainment loss counts as evidence for loss.",
                    "Replace phase-agnostic ĥ with the HMM-oriented major-haplotype fraction from experiment 16 when available.",
                    "Estimate k on samples with DNA-defined allele-specific states.",
                    "Run per cluster to expose subclonal or mirrored events."])
    pd.set_option("display.width", 250)
    print(sig_flank, sig_genome, hmed, tau)
    print(fam.to_string()); print(per_tab.drop(columns="note").to_string()); print(ksens.to_string())


if __name__ == "__main__":
    main()
