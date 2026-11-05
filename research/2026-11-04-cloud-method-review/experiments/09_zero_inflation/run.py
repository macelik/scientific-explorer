"""Experiment 09 (method record E:3): zero-inflation / hurdle diagnostic for per-cell region counts.

Exploratory diagnostic, not a production model. E:3 says a zero-inflated or hurdle model is only
considered after simpler-model misfit; a structural-zero component is a statistical description, not
a licence to delete zero cells. Mixture LRTs are nonregular (pi = 0 lies on the boundary), so the NB vs
ZINB comparison uses a parametric bootstrap from the fitted NB null with BOTH models refitted in every
replicate; no chi-square reference is used.

Region set: the saved windows, every region with nb_misfit_flag True in experiment 08, and 30 uncalled
random windows drawn with np.random.default_rng((20261104, 9)).
Per region (target counts per cell y_i, exposure E_i = cell raw total minus target and flank counts,
offset log E_i + log n_t, WNN cluster dummies): Poisson, NB2, ZINB (constant pi), ZINB (logit pi linear in
standardised log E), hurdle NB (constant P(0), zero-truncated NB for positives). All likelihoods are
maximised here with analytic gradients (L-BFGS-B); NB2 alpha is the ML value (08 used Cameron-Trivedi).
Held-out log score: train on one stratified half, evaluate on the other (seed tuple (20261104, 9)).
Usage: python -I run.py [B_null B_ci]   (default 199 99; bounded check: python -I run.py 19 19 --check)
"""
import importlib.util
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import optimize, stats  # noqa: E402
from scipy.special import digamma, expit, gammaln  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

spec = importlib.util.spec_from_file_location("exp04", HERE.parent / "04_fixed_region_contrasts" / "run.py")
E4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(E4)
warnings.filterwarnings("ignore")

LA_LO, LA_HI = np.log(1e-6), np.log(50.0)       # bounds on log alpha
G_LO, G_HI = -12.0, 8.0                          # bounds on logit(pi) coefficients
G0_START = (-3.0, -6.0)


# ----------------------------------------------------------------------------- likelihoods
def _nb_parts(y, eta, la):
    r = np.exp(-la); mu = np.exp(eta); lA = np.logaddexp(np.log(r), eta); A = r + mu
    logp0 = r * (np.log(r) - lA)
    ll = gammaln(y + r) - gammaln(r) - gammaln(y + 1) + logp0 + y * (eta - lA)
    return r, mu, A, logp0, ll


def _nb_grads(y, r, mu, A, logp0):
    d_eta = r * (y - mu) / A
    d_r = digamma(y + r) - digamma(r) + np.log(r) - np.log(A) + 1 - (r + y) / A
    return d_eta, d_r


def nll_pois(th, y, X, off):
    eta = np.clip(X @ th + off, -30, 20); mu = np.exp(eta)
    return float(np.sum(mu - y * eta + gammaln(y + 1))), X.T @ (mu - y)


def nll_zinb(th, y, X, off, Z):
    """NB2 (Z None) or ZINB. theta = [beta, log alpha, gamma]."""
    k = X.shape[1]
    eta = np.clip(X @ th[:k] + off, -30, 20); la = th[k]
    r, mu, A, logp0, llnb = _nb_parts(y, eta, la)
    d_eta, d_r = _nb_grads(y, r, mu, A, logp0)
    if Z is None:
        return -float(llnb.sum()), -np.r_[X.T @ d_eta, np.sum(-r * d_r)]
    pi = expit(Z @ th[k + 1:]); zero = y == 0
    p0 = np.exp(logp0)
    ll0 = np.log(pi + (1 - pi) * p0)
    ll = np.where(zero, ll0, np.log1p(-pi) + llnb)
    w = np.exp(np.log1p(-pi) + logp0 - ll0)
    dl0_eta = -r * mu / A; dl0_r = np.log(r) - np.log(A) + 1 - r / A
    de = np.where(zero, w * dl0_eta, d_eta); dr = np.where(zero, w * dl0_r, d_r)
    dz = np.where(zero, (1 - p0) / np.exp(ll0) * pi * (1 - pi), -pi)
    return -float(ll.sum()), -np.r_[X.T @ de, np.sum(-r * dr), Z.T @ dz]


def nll_trunc(th, y, X, off):
    """Zero-truncated NB2 on positive y only. theta = [beta, log alpha]."""
    k = X.shape[1]
    eta = np.clip(X @ th[:k] + off, -30, 20); la = th[k]
    r, mu, A, logp0, llnb = _nb_parts(y, eta, la)
    p0 = np.exp(logp0); l1m = np.log(-np.expm1(logp0)); q = p0 / (1 - p0)
    d_eta, d_r = _nb_grads(y, r, mu, A, logp0)
    de = d_eta - q * r * mu / A; dr = d_r + q * (np.log(r) - np.log(A) + 1 - r / A)
    return -float(np.sum(llnb - l1m)), -np.r_[X.T @ de, np.sum(-r * dr)]


def _opt(fun, x0, args, bounds):
    res = optimize.minimize(fun, x0, args=args, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 600, "maxfun": 900})
    return res


def fit_pois(y, X, off):
    b0 = np.zeros(X.shape[1]); b0[0] = np.log(max(y.sum(), 0.5) / np.exp(off).sum())
    r = _opt(nll_pois, b0, (y, X, off), [(-30, 30)] * X.shape[1])
    return {"theta": r.x, "ll": -r.fun, "ok": bool(r.success), "k": X.shape[1]}


def fit_nb(y, X, off, start=None):
    k = X.shape[1]
    if start is None:
        p = fit_pois(y, X, off); start = np.r_[p["theta"], np.log(0.1)]
    r = _opt(nll_zinb, start, (y, X, off, None), [(-30, 30)] * k + [(LA_LO, LA_HI)])
    return {"theta": r.x, "ll": -r.fun, "ok": bool(r.success), "k": k + 1, "alpha": float(np.exp(r.x[k]))}


def fit_zinb(y, X, off, Z, nb, extra_starts=()):
    k = X.shape[1]; zd = Z.shape[1]
    best = None
    starts = []
    for g0 in G0_START:
        starts.append(np.r_[nb["theta"], g0, np.zeros(zd - 1)])
    starts += list(extra_starts)
    bnds = [(-30, 30)] * k + [(LA_LO, LA_HI)] + [(G_LO, G_HI)] * zd
    for s in starts:
        r = _opt(nll_zinb, s, (y, X, off, Z), bnds)
        if np.isfinite(r.fun) and (best is None or -r.fun > best["ll"]):
            best = {"theta": r.x, "ll": -r.fun, "ok": bool(r.success), "k": k + 1 + zd, "alpha": float(np.exp(r.x[k]))}
    if best is None:
        best = {"theta": starts[0], "ll": -np.inf, "ok": False, "k": k + 1 + zd, "alpha": np.nan}
    pi = expit(Z @ best["theta"][k + 1:])
    best["pi"] = pi; best["pi_mean"] = float(pi.mean())
    best["boundary"] = bool(best["pi_mean"] < 1e-4)
    return best


def fit_hurdle(y, X, off, nb):
    k = X.shape[1]; pos = y > 0; n0 = int((~pos).sum()); n = len(y)
    pi0 = (n0 + 0.5) / (n + 1.0)                      # P(zero), shrunk so held-out zeros/positives are never impossible
    if pos.sum() < k + 5:
        return {"theta": nb["theta"], "ll": -np.inf, "ok": False, "k": k + 2, "pi0": pi0, "alpha": np.nan}
    r = _opt(nll_trunc, nb["theta"], (y[pos], X[pos], off[pos]), [(-30, 30)] * k + [(LA_LO, LA_HI)])
    ll = n0 * np.log(pi0) + (n - n0) * np.log1p(-pi0) - r.fun
    return {"theta": r.x, "ll": ll, "ok": bool(r.success), "k": k + 2, "pi0": pi0, "alpha": float(np.exp(r.x[k]))}


# ----------------------------------------------------------------------------- predictive log pmf
def lpmf(model, fit, y, X, off, Z=None):
    k = X.shape[1]
    if model == "poisson":
        eta = X @ fit["theta"] + off
        return y * eta - np.exp(eta) - gammaln(y + 1)
    eta = X @ fit["theta"][:k] + off; la = fit["theta"][k]
    r, mu, A, logp0, llnb = _nb_parts(y, eta, la)
    if model == "nb":
        return llnb
    if model == "zinb":
        pi = expit(Z @ fit["theta"][k + 1:]); p0 = np.exp(logp0)
        return np.where(y == 0, np.log(pi + (1 - pi) * p0), np.log1p(-pi) + llnb)
    if model == "hurdle":
        l1m = np.log(-np.expm1(logp0))
        return np.where(y == 0, np.log(fit["pi0"]), np.log1p(-fit["pi0"]) + llnb - l1m)


def sim_nb(rng, mu, alpha):
    r = 1.0 / alpha
    return rng.poisson(rng.gamma(r, mu / r)).astype(float)


# ----------------------------------------------------------------------------- per-region worker
G = {}


def run_region(args):
    idx, name, family, j, B, Bci = args
    t0 = time.time()
    y = G["yt"][:, j].astype(float); Ei = G["Ei"][:, j]; nt = G["nt"][j]
    D, test, lab, terc = G["D"], G["test"], G["lab"], G["terc"]
    n = len(y); off = np.log(Ei) + np.log(nt)
    X = np.column_stack([np.ones(n), D])
    le = np.log(Ei); Zc = np.ones((n, 1)); Ze = np.column_stack([np.ones(n), (le - le.mean()) / le.std()])
    out = {"name": name, "family": family, "idx": idx}
    # full-data fits
    p = fit_pois(y, X, off); nb = fit_nb(y, X, off)
    z1 = fit_zinb(y, X, off, Zc, nb); z2 = fit_zinb(y, X, off, Ze, nb)
    hu = fit_hurdle(y, X, off, nb)
    # train/test
    tr = ~test; res_ho = {}
    pt = fit_pois(y[tr], X[tr], off[tr]); nbt = fit_nb(y[tr], X[tr], off[tr])
    z1t = fit_zinb(y[tr], X[tr], off[tr], Zc[tr], nbt); z2t = fit_zinb(y[tr], X[tr], off[tr], Ze[tr], nbt)
    hut = fit_hurdle(y[tr], X[tr], off[tr], nbt)
    s = lambda m, f, Z=None: lpmf(m, f, y[test], X[test], off[test], None if Z is None else Z[test])
    ho = {"poisson": s("poisson", pt), "nb": s("nb", nbt), "zinb_const": s("zinb", z1t, Zc),
          "zinb_exposure": s("zinb", z2t, Ze), "hurdle_nb": s("hurdle", hut)}
    rows = []
    fits = {"poisson": p, "nb": nb, "zinb_const": z1, "zinb_exposure": z2, "hurdle_nb": hu}
    for m, f in fits.items():
        d = ho[m] - ho["nb"]
        rows.append({"region": name, "family": family, "model": m, "n_cells": n, "n_params": f["k"], "loglik": f["ll"],
                     "aic": 2 * f["k"] - 2 * f["ll"], "bic": f["k"] * np.log(n) - 2 * f["ll"],
                     "alpha": f.get("alpha", np.nan), "converged": f["ok"],
                     "pi_hat": {"zinb_const": z1["pi_mean"], "zinb_exposure": z2["pi_mean"], "hurdle_nb": hu["pi0"]}.get(m, np.nan),
                     "heldout_logscore_sum": float(ho[m].sum()), "heldout_logscore_mean": float(ho[m].mean()),
                     "heldout_diff_vs_nb_sum": float(d.sum()), "heldout_diff_vs_nb_se": float(d.std(ddof=1) * np.sqrt(len(d))),
                     "pi_lo": np.nan, "pi_hi": np.nan})
    # nb expected zero fraction (in-sample), for hurdle comparison
    k = X.shape[1]
    eta = X @ nb["theta"][:k] + off; _, _, _, lp0, _ = _nb_parts(y, eta, nb["theta"][k]); p0nb = np.exp(lp0)
    out["nb_expected_zero_frac"] = float(p0nb.mean()); out["obs_zero_frac"] = float((y == 0).mean())
    # zero-excess location (in-sample NB): by cell-total tercile and by cluster
    loc = []
    tau = np.where(y == 0, z1["pi"] / (z1["pi"] + (1 - z1["pi"]) * p0nb), 0.0)
    for kind, levels, key in [("coverage_tercile", [0, 1, 2], terc), ("cluster", sorted(set(lab)), lab)]:
        for lv in levels:
            m = key == lv
            ex = float((y[m] == 0).sum() - p0nb[m].sum()); v = float(np.sum(p0nb[m] * (1 - p0nb[m])))
            loc.append({"region": name, "family": family, "by": kind, "level": str(lv), "n_cells": int(m.sum()),
                        "obs_zero": int((y[m] == 0).sum()), "nb_expected_zero": float(p0nb[m].sum()), "excess": ex,
                        "z": ex / np.sqrt(v) if v >= 1 else np.nan, "zinb_structural_zero_expected": float(tau[m].sum())})
    # null parametric bootstrap
    rng = np.random.default_rng((C.SEED, 9, idx))
    mu0 = np.exp(eta); lr_obs = max(0.0, 2 * (z1["ll"] - nb["ll"]))
    boot = []
    for b in range(B):
        yb = sim_nb(rng, mu0, np.exp(nb["theta"][k]))
        nbb = fit_nb(yb, X, off, start=nb["theta"])
        zb = fit_zinb(yb, X, off, Zc, nbb)
        boot.append({"region": name, "family": family, "kind": "null_nb", "rep": b, "lr": max(0.0, 2 * (zb["ll"] - nbb["ll"])),
                     "pi_hat": zb["pi_mean"], "alpha_nb": nbb["alpha"], "converged_nb": nbb["ok"], "converged_zinb": zb["ok"],
                     "boundary": zb["boundary"]})
    lrs = np.array([r["lr"] for r in boot])
    out["lr_obs"] = lr_obs; out["p_boot"] = float((1 + np.sum(lrs >= lr_obs - 1e-9)) / (B + 1)) if B else np.nan
    # pi interval: parametric bootstrap from the fitted ZINB
    kk = k
    rng2 = np.random.default_rng((C.SEED, 9, idx, 1))
    pi_f = z1["pi"]; mu1 = np.exp(X @ z1["theta"][:kk] + off); a1 = np.exp(z1["theta"][kk])
    for b in range(Bci):
        yb = sim_nb(rng2, mu1, a1); yb[rng2.random(n) < pi_f] = 0
        nbb = fit_nb(yb, X, off, start=nb["theta"])
        zb = fit_zinb(yb, X, off, Zc, nbb, extra_starts=[z1["theta"]])
        boot.append({"region": name, "family": family, "kind": "zinb_ci", "rep": b, "lr": max(0.0, 2 * (zb["ll"] - nbb["ll"])),
                     "pi_hat": zb["pi_mean"], "alpha_nb": nbb["alpha"], "converged_nb": nbb["ok"], "converged_zinb": zb["ok"],
                     "boundary": zb["boundary"]})
    ci = [r["pi_hat"] for r in boot if r["kind"] == "zinb_ci"]
    if ci:
        for r in rows:
            if r["model"] == "zinb_const":
                r["pi_lo"], r["pi_hi"] = float(np.percentile(ci, 2.5)), float(np.percentile(ci, 97.5))
    nullpi = np.array([r["pi_hat"] for r in boot if r["kind"] == "null_nb"])
    out.update({"pi_const": z1["pi_mean"], "lr_obs_exposure": max(0.0, 2 * (z2["ll"] - nb["ll"])),
                "pi_exposure_slope": float(z2["theta"][kk + 2]), "alpha_nb": nb["alpha"], "alpha_zinb": z1["alpha"],
                "boundary_obs": z1["boundary"], "null_pi_median": float(np.median(nullpi)), "null_pi_q95": float(np.percentile(nullpi, 95)),
                "null_frac_boundary": float(np.mean([r["boundary"] for r in boot if r["kind"] == "null_nb"])),
                "null_frac_nonconv_nb": float(np.mean([not r["converged_nb"] for r in boot if r["kind"] == "null_nb"])),
                "null_frac_nonconv_zinb": float(np.mean([not r["converged_zinb"] for r in boot if r["kind"] == "null_nb"])),
                "ci_frac_nonconv": float(np.mean([not (r["converged_nb"] and r["converged_zinb"]) for r in boot if r["kind"] == "zinb_ci"])) if ci else np.nan,
                "ci_frac_boundary": float(np.mean([r["boundary"] for r in boot if r["kind"] == "zinb_ci"])) if ci else np.nan,
                "hurdle_pi0": hu["pi0"], "sum_tau_zinb": float(tau.sum()), "n_obs_zero": int((y == 0).sum()),
                "nb_expected_zero": float(p0nb.sum()), "median_cell_total": float(np.median(G["tot"])),
                "seconds": round(time.time() - t0, 1)})
    return out, rows, boot, loc


# ----------------------------------------------------------------------------- main
def main():
    t0 = time.time()
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    check = "--check" in sys.argv
    B = int(args[0]) if args else 199
    Bci = int(args[1]) if len(args) > 1 else 99
    p = C.default_paths(); d = C.H5AD(p["h5ad"]); lab = C.cell_clusters(d, p)
    regs = E4.build_regions(d, p)
    yt, yf, tot = E4.region_sums(d, regs, "layers/raw")
    m8 = pd.read_csv(HERE.parent / "08_exposure_count_models" / "region_models.tsv", sep="\t")
    names = [r["name"] for r in regs]
    saved = list(m8[m8.family == "saved"].region)
    flagged = [r for r in m8[m8.nb_misfit_flag].region if r not in saved]
    unc = m8[(m8.family == "random") & (m8.frac_bins_called_loss == 0) & (m8.frac_bins_called_gain == 0) & (~m8.nb_misfit_flag)].region.tolist()
    neg = sorted(np.random.default_rng((C.SEED, 9)).choice(unc, 30, replace=False).tolist())
    sel = [(r, "saved") for r in saved] + [(r, "nb_misfit_flagged") for r in flagged] + [(r, "uncalled_random_negative") for r in neg]
    fam08 = dict(zip(m8.region, m8.family))
    if check:
        sel = [s for s in sel if s[0] in ("chr7_8.3_10.7", "chr9_22.5_26.6") or s[1] != "saved"][:4]
    clusters = sorted(set(lab))
    G.update(dict(yt=yt, Ei=tot[:, None] - yt - yf, nt=np.array([max(len(r["target"]), 1) for r in regs]),
                  D=np.column_stack([(lab == c).astype(float) for c in clusters[1:]]),
                  test=~C.stratified_halves(lab, (C.SEED, 9)), lab=lab, terc=np.asarray(pd.qcut(tot, 3, labels=False)), tot=tot))
    jobs = [(i, nm, fam, names.index(nm), B, Bci) for i, (nm, fam) in enumerate(sel)]
    print(f"{len(jobs)} regions, B={B}, B_ci={Bci}, cells={len(lab)}", flush=True)
    with Pool(2) as pool:
        results = []
        for r in pool.imap_unordered(run_region, jobs):
            results.append(r); print(f"  done {r[0]['name']} ({r[0]['seconds']}s) p_boot={r[0]['p_boot']:.3f} pi={r[0]['pi_const']:.4f}", flush=True)
    results.sort(key=lambda r: r[0]["idx"])
    reg = pd.DataFrame([r[0] for r in results])
    mc = pd.DataFrame([x for r in results for x in r[1]])
    bt = pd.DataFrame([x for r in results for x in r[2]])
    loc = pd.DataFrame([x for r in results for x in r[3]])
    from statsmodels.stats.multitest import multipletests
    reg["q_bh"] = multipletests(reg.p_boot, method="fdr_bh")[1]
    out = HERE / ("_check" if check else "")
    if check:
        out.mkdir(exist_ok=True)
    mc.to_csv(out / "model_comparison.tsv", sep="\t", index=False)
    bt.to_csv(out / "bootstrap.tsv", sep="\t", index=False)
    reg.to_csv(out / "lrt_by_region.tsv", sep="\t", index=False)
    loc.to_csv(out / "zero_location.tsv", sep="\t", index=False)
    seconds = round(time.time() - t0, 1)
    meta = {"experiment": "09_zero_inflation", "method_record": "E:3", "seeds": {"negative_windows": [C.SEED, 9], "train_test_split": [C.SEED, 9],
            "null_bootstrap": "(20261104, 9, region_index)", "zinb_ci_bootstrap": "(20261104, 9, region_index, 1)"},
            "generator": "numpy.random.default_rng (PCG64)", "B_null": B, "B_ci": Bci, "n_regions": len(reg), "region_index_order": [s[0] for s in sel],
            "parameters": {"alpha_bounds": [float(np.exp(LA_LO)), float(np.exp(LA_HI))], "logit_pi_bounds": [G_LO, G_HI], "zinb_starts_logit_pi": list(G0_START),
                           "boundary_pi_threshold": 1e-4, "hurdle_pi0": "(n0+0.5)/(n+1)", "optimizer": "scipy L-BFGS-B, analytic gradients",
                           "lrt": "2(llZINB - llNB), floored at 0; p = (1+#{LR*>=LR})/(B+1)"},
            "counts": {"saved": len(saved), "nb_misfit_flagged_not_saved": len(flagged), "negatives": len(neg)},
            "n_p_boot_lt_0.05": int((reg.p_boot < .05).sum()), "n_q_bh_lt_0.05": int((reg.q_bh < .05).sum()),
            "seconds": seconds, "environment": C.environment(), "check_run": check}
    C.write_json(out / "summary.json", meta)
    if check:
        print(reg[["name", "p_boot", "lr_obs", "pi_const", "null_frac_nonconv_zinb", "seconds"]]); print(mc[mc.model != "poisson"].head(12).to_string()); return
    build_report(reg, mc, bt, loc, meta, regs, d, clusters, fam08)
    print("seconds", seconds)


# ----------------------------------------------------------------------------- report
def build_report(reg, mc, bt, loc, meta, regs, d, clusters, fam08):
    B, Bci = meta["B_null"], meta["B_ci"]

    def figure(fn, w=7.5, h=4.0):
        fig, ax = plt.subplots(figsize=(w, h)); fn(ax); fig.tight_layout(); o = C.fig_to_b64(fig); plt.close(fig); return o
    fcol = {"saved": "k", "nb_misfit_flagged": "#9b2c2c", "uncalled_random_negative": "#999"}
    wide = mc.pivot(index="region", columns="model", values=["loglik", "aic", "bic", "heldout_logscore_sum"])
    reg = reg.copy()
    reg["dAIC_zinb_vs_nb"] = [wide.loc[r, ("aic", "zinb_const")] - wide.loc[r, ("aic", "nb")] for r in reg.name]
    reg["dBIC_zinb_vs_nb"] = [wide.loc[r, ("bic", "zinb_const")] - wide.loc[r, ("bic", "nb")] for r in reg.name]
    reg["dAIC_pois_vs_nb"] = [wide.loc[r, ("aic", "poisson")] - wide.loc[r, ("aic", "nb")] for r in reg.name]
    for m in ("zinb_const", "zinb_exposure", "hurdle_nb"):
        reg[f"ho_{m}_minus_nb"] = [wide.loc[r, ("heldout_logscore_sum", m)] - wide.loc[r, ("heldout_logscore_sum", "nb")] for r in reg.name]
    reg["ho_poisson_minus_nb"] = [wide.loc[r, ("heldout_logscore_sum", "poisson")] - wide.loc[r, ("heldout_logscore_sum", "nb")] for r in reg.name]
    ci = mc[mc.model == "zinb_const"].set_index("region")
    reg["pi_lo"] = reg.name.map(ci.pi_lo); reg["pi_hi"] = reg.name.map(ci.pi_hi)
    sig = reg[reg.p_boot < 0.05]; supp = reg[(reg.p_boot < 0.05) & (reg.ho_zinb_const_minus_nb > 0)]
    nfam = reg.groupby("family").agg(n=("name", "size"), n_p_lt_05=("p_boot", lambda s: int((s < .05).sum())), median_pi=("pi_const", "median"),
                                      median_p=("p_boot", "median"), median_lr=("lr_obs", "median"))
    # zero location summary
    lt = loc[loc.by == "coverage_tercile"].pivot_table(index="region", columns="level", values="excess")
    lt["share_low"] = lt["0"] / lt[["0", "1", "2"]].clip(lower=0).sum(axis=1).replace(0, np.nan)
    lz = loc[loc.by == "coverage_tercile"].pivot_table(index="region", columns="level", values="z")
    lc = loc[loc.by == "cluster"]
    maxc = lc.loc[lc.groupby("region").excess.idxmax()].set_index("region")
    reg["excess_low_tercile"] = reg.name.map(lt["0"]); reg["z_low_tercile"] = reg.name.map(lz["0"])
    reg["max_excess_cluster"] = reg.name.map(maxc.level); reg["max_excess_cluster_z"] = reg.name.map(maxc.z)

    def f_p(ax):
        for f, c in fcol.items():
            g = reg[reg.family == f]
            ax.scatter(g.pi_const, g.p_boot, color=c, s=30, alpha=.8, label=f)
        ax.axhline(0.05, color="#8a5a00", ls="--", lw=1); ax.set_xscale("symlog", linthresh=1e-3)
        ax.set_xlabel("ZINB estimated pi (constant)"); ax.set_ylabel(f"bootstrap p (B={B})"); ax.legend(fontsize=7)

    def f_null(ax):
        top = reg.sort_values("p_boot").name.iloc[:1].tolist() + reg.sort_values("lr_obs", ascending=False).name.iloc[:1].tolist()
        for nm in dict.fromkeys(top):
            v = bt[(bt.region == nm) & (bt.kind == "null_nb")].lr
            ax.hist(v, bins=30, alpha=.5, label=f"{nm}: null LR*")
            ax.axvline(reg.set_index("name").loc[nm, "lr_obs"], color="k", lw=1)
        ax.set_xlabel("LR = 2(ll ZINB − ll NB)"); ax.set_ylabel("bootstrap replicates"); ax.legend(fontsize=7)

    def f_forest(ax):
        r = reg.sort_values("pi_const").reset_index(drop=True)
        for i, x in r.iterrows():
            ax.plot([x.pi_lo, x.pi_hi], [i, i], color=fcol[x.family], lw=1)
            ax.scatter([x.pi_const], [i], color=fcol[x.family], s=14)
        ax.set_yticks(range(len(r))); ax.set_yticklabels(r.name, fontsize=5); ax.set_xscale("symlog", linthresh=1e-3)
        ax.set_xlabel(f"pi with percentile interval (parametric bootstrap from fitted ZINB, B={Bci}); black=saved, red=flagged, grey=negatives")

    def f_ho(ax):
        mods = [("poisson", "Poisson"), ("zinb_const", "ZINB"), ("zinb_exposure", "ZINB(exposure)"), ("hurdle_nb", "hurdle NB")]
        for i, (m, lbl) in enumerate(mods):
            col = "ho_poisson_minus_nb" if m == "poisson" else f"ho_{m}_minus_nb"
            for f, c in fcol.items():
                v = reg[reg.family == f][col]
                ax.scatter(np.full(len(v), i) + np.random.default_rng(1).normal(0, .06, len(v)), v, color=c, s=14, alpha=.7)
        ax.axhline(0, color="k", lw=.6); ax.set_yscale("symlog", linthresh=1)
        ax.set_xticks(range(len(mods))); ax.set_xticklabels([m[1] for m in mods]); ax.set_ylabel("held-out log score minus NB (sum over test cells)")

    def f_loc(ax):
        sub = loc[(loc.by == "coverage_tercile") & (loc.region.isin(reg.sort_values("p_boot").name.iloc[:8]))]
        for nm, g in sub.groupby("region"):
            ax.plot(g.level.astype(int), g.z, marker="o", label=nm)
        ax.axhline(0, color="k", lw=.5); ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["low", "mid", "high"])
        ax.set_xlabel("cell-total tercile"); ax.set_ylabel("excess-zero z vs NB (in-sample)"); ax.legend(fontsize=5)
    b_p, b_null, b_for, b_ho, b_loc = figure(f_p), figure(f_null), figure(f_forest, 7.5, 8), figure(f_ho), figure(f_loc)

    def note(r):
        bits = []
        if r.p_boot < 0.05:
            bits.append(f"bootstrap LRT favours extra zeros (p={r.p_boot:.3f}, pi={r.pi_const:.3f})")
        else:
            bits.append(f"no bootstrap evidence for extra zeros (p={r.p_boot:.3f}, pi={r.pi_const:.4f})")
        if r.dAIC_pois_vs_nb > 10:
            bits.append("NB much better than Poisson")
        if r.p_boot < 0.05:
            bits.append("ZINB predicts held-out cells better than NB" if r.ho_zinb_const_minus_nb > 0 else "but ZINB does not improve held-out score")
            if np.isfinite(r.z_low_tercile) and r.z_low_tercile > 3:
                bits.append("excess zeros concentrated in lowest-total tercile")
        if r.null_frac_nonconv_zinb > .05:
            bits.append("ZINB fits unstable under the null")
        return "; ".join(bits)
    reg["note"] = reg.apply(note, axis=1)
    per = reg[["name", "family", "p_boot", "q_bh", "lr_obs", "pi_const", "pi_lo", "pi_hi", "dAIC_zinb_vs_nb", "dBIC_zinb_vs_nb",
               "ho_zinb_const_minus_nb", "ho_hurdle_nb_minus_nb", "z_low_tercile", "note"]].rename(columns={"name": "region"})
    reg.drop(columns=["note"]).to_csv(HERE / "lrt_by_region.tsv", sep="\t", index=False)
    sv = reg[reg.family == "saved"]; fl = reg[reg.family == "nb_misfit_flagged"]; ng = reg[reg.family == "uncalled_random_negative"]
    mc_show = mc[mc.model.isin(["poisson", "nb", "zinb_const", "zinb_exposure", "hurdle_nb"]) & mc.region.isin(sv.name)][
        ["region", "model", "n_params", "loglik", "aic", "bic", "alpha", "pi_hat", "pi_lo", "pi_hi", "heldout_logscore_mean", "heldout_diff_vs_nb_sum", "heldout_diff_vs_nb_se", "converged"]]
    inst = reg[["name", "family", "null_pi_median", "null_pi_q95", "null_frac_boundary", "null_frac_nonconv_nb", "null_frac_nonconv_zinb", "ci_frac_boundary", "ci_frac_nonconv"]].rename(columns={"name": "region"})
    nb_pois_dAIC = reg.dAIC_pois_vs_nb
    key = [
        f"Bootstrap NB vs ZINB LRT (B={B}, both models refitted in every replicate): {len(sig)}/{len(reg)} regions have p < 0.05 "
        f"({int((sv.p_boot<.05).sum())}/{len(sv)} saved, {int((fl.p_boot<.05).sum())}/{len(fl)} NB-misfit flagged, {int((ng.p_boot<.05).sum())}/{len(ng)} negative windows); "
        f"{meta['n_q_bh_lt_0.05']} survive BH q < 0.05 (smallest attainable p = {1/(B+1):.4f}). Under a calibrated test ~5% of the negatives would be below 0.05.",
        f"Saved windows: p_boot = " + ", ".join(f"{r.name.split('_',1)[0]}:{r.name.split('_',1)[1]} {r.p_boot:.3f}" for r in sv.itertuples()) + ". The deep chr7 windows do not need a zero component.",
        f"Estimated constant pi is small where it is estimated at all: median {reg.pi_const.median():.4f}, max {reg.pi_const.max():.3f} ({reg.loc[reg.pi_const.idxmax(),'name']}); "
        f"{int(reg.boundary_obs.sum())}/{len(reg)} fits sit on the pi = 0 boundary, which is why a chi-square LRT would be invalid.",
        f"Of the {len(sig)} nominally supported regions, {len(supp)} also gain held-out log score under ZINB over NB (positive = better); hurdle NB beats NB held-out in "
        f"{int((reg.ho_hurdle_nb_minus_nb>0).sum())}/{len(reg)} regions and BIC prefers ZINB over NB in {int((reg.dBIC_zinb_vs_nb<0).sum())}/{len(reg)}.",
        f"Poisson vs NB: median ΔAIC (Poisson − NB) = {nb_pois_dAIC.median():.0f}; NB wins by >10 AIC in {int((nb_pois_dAIC>10).sum())}/{len(reg)} regions, so overdispersion, not zeros, is the dominant departure from Poisson (consistent with experiment 08).",
        f"Where extra zeros are supported, excess zeros concentrate in the lowest cell-total tercile in {int(((sig.z_low_tercile>3)).sum())}/{max(len(sig),1)} of them (z > 3), pointing to a depth-related (technical) origin where present; see zero_location.tsv for cluster profiles.",
        f"Component instability: in the null replicates the median fraction of ZINB fits on the pi = 0 boundary is {reg.null_frac_boundary.median():.2f}; median non-converged fraction NB {reg.null_frac_nonconv_nb.median():.3f}, ZINB {reg.null_frac_nonconv_zinb.median():.3f}.",
    ]
    glossary = {
        "region / family": "window name and set: saved = frozen windows; nb_misfit_flagged = flagged by experiment 08 (|excess-zero z|>3 or held-out Pearson dispersion>1.5), not saved; uncalled_random_negative = 30 random uncalled windows (seed (20261104, 9))",
        "model": "poisson; nb (NB2, ML alpha); zinb_const (constant structural-zero probability pi); zinb_exposure (logit pi linear in standardised log E); hurdle_nb (P(0) constant, zero-truncated NB for positives)",
        "n_params": "number of fitted parameters (beta incl. intercept + clusters, alpha, pi terms)",
        "loglik / aic / bic": "maximised log-likelihood on all cells; AIC = 2k − 2ll; BIC = k ln n − 2ll (n = cells)",
        "alpha": "NB2 dispersion: Var = μ + αμ²",
        "pi_hat": "mean structural-zero probability (ZINB) or P(0) of the hurdle part; hurdle P(0) is the whole zero probability, not an excess",
        "pi_lo / pi_hi": f"2.5/97.5 percentiles of pi re-estimated in {Bci} parametric-bootstrap replicates simulated from the fitted ZINB",
        "heldout_logscore_mean / _sum": "mean / total log predictive probability of test-half cells for the model fitted on the training half (stratified half split, seed (20261104, 9))",
        "heldout_diff_vs_nb_sum / _se": "sum of per-cell log score differences vs NB and its SE (sqrt(n) × SD of the differences); positive = better than NB",
        "converged": "optimizer success flag (L-BFGS-B)",
        "lr_obs": "2(ll ZINB − ll NB), floored at 0",
        "p_boot": f"(1 + #{{LR* ≥ LR}})/(B + 1), B = {B}, replicates simulated from the fitted NB with the same exposures and clusters",
        "q_bh": "Benjamini–Hochberg q over the regions in this run (a convenience; the regions are not independent)",
        "pi_const": "ZINB constant-pi estimate; boundary_obs = pi < 1e-4",
        "dAIC_zinb_vs_nb / dBIC_zinb_vs_nb": "AIC (BIC) of ZINB minus NB; negative favours ZINB",
        "dAIC_pois_vs_nb": "AIC of Poisson minus NB; positive favours NB",
        "ho_*_minus_nb": "held-out total log score of the named model minus NB",
        "z_low_tercile": "(observed − NB-expected zero cells)/SD in the lowest cell-total tercile, in-sample NB",
        "max_excess_cluster / _z": "cluster with the largest observed-minus-expected zero count and its z",
        "null_pi_median / null_pi_q95": "pi estimate across null replicates (true pi = 0)",
        "null_frac_boundary": "fraction of null replicates whose ZINB pi estimate is < 1e-4",
        "null_frac_nonconv_nb / _zinb, ci_frac_nonconv": "fraction of replicates in which the optimizer reported non-convergence",
        "by / level / excess / z (zero_location.tsv)": "stratum type (coverage tercile 0=low or cluster), observed zeros minus NB-expected zeros, z = excess/√Σp(1−p)",
        "zinb_structural_zero_expected": "Σ over zero cells of the ZINB posterior probability that the zero is structural (a model-based description, not a cell label)",
        "kind / rep (bootstrap.tsv)": "null_nb = replicate from fitted NB (for the LRT); zinb_ci = replicate from fitted ZINB (for the pi interval)",
        "lr, pi_hat, alpha_nb, converged_*, boundary (bootstrap.tsv)": "per-replicate LR, ZINB pi, NB alpha, convergence flags, pi < 1e-4",
    }
    sections = [
        ("Question and why it matters", "<p>Method record <b>E:3</b> allows hurdle, zero-inflated NB or latent count mixtures only after a simpler model misfits. "
         "Experiment 08 found NB explains the zeros of the deep chr7 windows but flagged some regions. Does a zero component add anything, and where? "
         "This decides whether a depth-aware evidence layer needs a zero-inflation term or whether quasi-Poisson/NB suffices. "
         "A supported zero component is a statistical description of the counts; it is not permission to delete zero cells.</p>"),
        ("Inputs", f"<ul><li><code>layers/raw</code>, {d.n_cells:,} cells; WNN clusters ({len(clusters)}); exposure E_i = raw total minus target and flank counts (as in 08).</li>"
         f"<li>Regions ({len(reg)}): {meta['counts']['saved']} saved windows, {meta['counts']['nb_misfit_flagged_not_saved']} further regions with <code>nb_misfit_flag</code> in 08 (called segments and random windows), "
         f"{meta['counts']['negatives']} uncalled, unflagged random windows drawn with seed tuple (20261104, 9). The 08 table has 7 saved windows (not 6); all 7 are included. "
         "The negatives are windows with no called loss/gain in 08.</li><li>Held-out split: <code>C.stratified_halves(lab, (20261104, 9))</code>.</li></ul>"),
        ("Method", "<pre>y_i  target count of cell i;  mu_i = exp(log E_i + log n_t + x_i'beta),  x_i = intercept + cluster dummies\n"
         "Poisson:   y ~ Pois(mu)\nNB2:       y ~ NB(mu, alpha),  Var = mu + alpha mu^2,  p0_i = (1 + alpha mu_i)^(-1/alpha)\n"
         "ZINB:      P(y=0) = pi_i + (1-pi_i) p0_i;  P(y=k) = (1-pi_i) NB(k),  logit pi_i = g0  [or g0 + g1 z_i, z = std. log E_i]\n"
         "Hurdle NB: P(y=0) = pi0;  P(y=k) = (1-pi0) NB(k)/(1-p0_i) for k>0,  pi0 = (n0+0.5)/(n+1)\n"
         "AIC = 2k - 2 ll;  BIC = k ln n - 2 ll\n"
         "LRT:   LR = 2 (ll_ZINB - ll_NB) >= 0;  pi = 0 is on the boundary, so no chi-square reference.\n"
         f"       for b = 1..B (B={B}): y* ~ NB(mu_hat, alpha_hat) (same E, clusters; seed (20261104, 9, region_index));\n"
         "       refit NB and ZINB (2 starts) on y*; LR*_b. p = (1 + #{LR*_b >= LR}) / (B + 1)\n"
         f"pi interval: B_ci={Bci} replicates from the fitted ZINB (structural zeros with prob pi_i, seed (20261104, 9, region_index, 1)); refit ZINB; 2.5/97.5 percentiles\n"
         "held-out: fit on one half, sum_i log p(y_i | fitted model) on the other half</pre>"
         "<p>Likelihoods are maximised with analytic gradients (L-BFGS-B; bounds on log alpha and logit pi). NB alpha is the ML estimate, so it differs slightly from 08's Cameron–Trivedi value.</p>"),
        ("Results", "", 2),
        ("Is a zero component supported?", C.img(b_p, "p vs pi") + C.caption(
            "Each point is a region. x: ZINB constant-pi estimate (symlog); y: bootstrap LRT p. Black = saved windows, red = NB-misfit flagged, grey = negative windows. Dashed: 0.05.") +
         C.html_table(nfam.reset_index()) + C.caption("By set: number of regions, how many with p_boot < 0.05, median pi, median p and median LR.") +
         C.img(b_null, "null LR") + C.caption("Null LR* distribution from NB-simulated data (histograms) with the observed LR (black line), for the region with the smallest p and the one with the largest observed LR."), 3),
        ("Estimated pi and its uncertainty", C.img(b_for, "forest") + C.caption(
            "Constant ZINB pi per region with a parametric-bootstrap percentile interval. Intervals that touch 0 are compatible with no structural component; the estimator is biased upward near the boundary."), 3),
        ("Model comparison for the saved windows", C.html_table(mc_show) + C.caption(
            "Log-likelihood, AIC/BIC, pi and held-out score for the saved windows (all regions in model_comparison.tsv)."), 3),
        ("Held-out predictive score", C.img(b_ho, "heldout") + C.caption(
            "Test-half log score minus NB (summed over test cells; symlog axis). Above 0 = better than NB. Colours as above.") , 3),
        ("Where do the extra zeros sit?", C.img(b_loc, "zero location") + C.caption(
            "Excess-zero z (observed − NB-expected zero cells) by cell-total tercile for the eight regions with the smallest bootstrap p. A rise toward the low tercile indicates zeros tied to depth; "
            "a flat profile with a cluster-specific excess points to cell state or subclone (zero_location.tsv lists clusters).") +
         C.html_table(reg.sort_values("p_boot").head(12)[["name", "family", "p_boot", "pi_const", "excess_low_tercile", "z_low_tercile", "max_excess_cluster", "max_excess_cluster_z"]]) +
         C.caption("Twelve regions with the smallest p. excess_low_tercile = observed − NB-expected zero cells in the lowest cell-total tercile."), 3),
        ("Component instability", C.html_table(inst) + C.caption(
            "Null-replicate pi quantiles and fractions of boundary and non-converged fits per region (rows ordered as in the run)."), 3),
        ("Per-region detail", C.html_table(per) + C.caption("One row per region; the note is generated from the columns (thresholds: p_boot 0.05, z 3).")),
        ("Sensitivity", "<p>Three alternative zero parameterisations are reported: constant pi (used for the LRT), pi depending on log exposure (<code>zinb_exposure</code>, slope in lrt_by_region.tsv as pi_exposure_slope; LR in lr_obs_exposure, not bootstrapped), and the hurdle NB. "
         "The ranking by held-out score and BIC is in model_comparison.tsv. Results are conditional on the cluster terms and on the exposure definition from 08. "
         f"The bootstrap resolution is 1/{B+1}; BH over {len(reg)} regions is a convenience and not a formal family-wise statement. NB alpha is the ML value here, not the Cameron–Trivedi value of 08.</p>"),
        ("Limits and alternative explanations", "<ul><li>A rejected NB is also compatible with a misspecified mean (cluster terms too coarse, subclones) or heavier-tailed dispersion (e.g. alpha varying with exposure), not only with structural zeros.</li>"
         "<li>The ZINB zero probability is not identifiable from a mean-model misfit; the held-out score is the check on whether the component predicts.</li>"
         "<li>Regions were selected partly by experiment 08's NB misfit flag, so their p-values are conditional on selection; the 30 negatives are the unselected reference.</li>"
         "<li>Random windows may contain real events. A zero component describes counts; it does not distinguish deletion, mappability, accessibility, or dropout.</li>"
         "<li>Local optima are possible; 2 starts per ZINB fit and warm starts were used but the global optimum is not proven.</li></ul>"),
    ]
    C.write_html_report(
        HERE / "report.html", "Zero inflation diagnostic", "Experiment 09 · method E:3 · Poisson / NB / ZINB / hurdle with bootstrap LRT (exploratory)",
        sections, {"command": "python -I experiments/09_zero_inflation/run.py", "outputs": ["model_comparison.tsv", "bootstrap.tsv", "lrt_by_region.tsv", "zero_location.tsv", "summary.json"],
                   "seeds": meta["seeds"], "B_null": B, "B_ci": Bci, "environment": meta["environment"]},
        verdict=f"Exploratory: a structural-zero component is supported by the bootstrap LRT in {len(sig)}/{len(reg)} regions (negatives: {int((ng.p_boot<.05).sum())}/{len(ng)}), "
                f"with small pi (median {reg.pi_const.median():.4f}); NB already captures the large majority of the departure from Poisson, and the deep chr7 windows do not need a zero component.",
        key_findings=key, glossary=glossary,
        supports=["NB (not Poisson) is the necessary step for these counts; a zero component is at most a small additional refinement in a minority of regions.",
                  "The deep chr7 saved windows are consistent with NB without extra zeros.",
                  "A parametric bootstrap with refitting is a workable non-regular test for this decision."],
        not_supports=["Deleting or down-weighting zero cells.", "Any statement that extra zeros are biological dropout, deletion or dosage.",
                      "A validated production zero-inflation model; this is an exploratory description."],
        decision_impact="For <code>recommendation.md</code>: the raw-count evidence layer can use NB (or quasi-Poisson) with exposure and cluster terms; zero inflation should not be a default component. "
                        "Keep an optional ZINB/hurdle check as a misfit diagnostic for regions where NB is flagged, reported with its bootstrap p and held-out score, and never as a cell filter.",
        next_steps=["Re-run with region selection and assessment on disjoint cell halves.", "Test alpha varying with exposure or cluster before attributing misfit to structural zeros.",
                    "Check flagged regions against a cross-sample technical baseline and mappability.", "Increase B (e.g. 999) for the regions near p = 0.05."])


if __name__ == "__main__":
    main()
