"""Experiment 15: observation models for the allelic modality at the CELL level (exploratory diagnostic).

Question: a later combined depth + allele model needs one adequate likelihood per modality. Depth is
well described by NB2 with exposure (exp 08). Is the same true for the allele total AT and for the
haplotype split HB | AT, per cell x region?

Units. y_ir = count for cell i in region r (sum over the region's retained target bins). AT = summed
per-SNP-site DP over selected phased het SNPs (NOT unique fragments; AT == HA + HB exactly, exp 06).
Regions: 7 saved windows, 4 literature loci, 60 uncalled random 25-bin windows (seed tuple
(20261104, 15), windows whose bins all have cn_call == 1 in the `all` karyogram), 67 called short
segments (exp 04 rule).

A. Count models for y (modality 'AT'; and 'depth' for context):
     log mu_i = offset_i + b0 + g_cluster(i)            (offset coefficient fixed at 1)
     offsets: 'raw_in' = log(cell raw depth in region + 0.5)                     [AT only]
              'out'    = log(cell AT total outside region + 0.5) for AT; log(cell raw total outside region + 0.5) for depth
     Poisson; NB2 (var = mu + a mu^2); ZIP / ZINB (P(excess zero) = expit(g0 + g1 z), z = standardised offset);
     hurdle-NB (P(0) = expit(g0 + g1 z), zero-truncated NB for y > 0); sensitivity: NB2 with free offset slope.
   Maximum likelihood on the TRAIN half of stratified_halves(lab, (20261104, 15)); score = mean per-cell log
   predictive probability on the TEST half (primary: train = True half; reverse direction = sensitivity).
   Held-out zero fraction vs predicted P(Y=0) by cell coverage tercile (tercile of cell raw total).
B. Haplotype split HB | AT, cells with AT >= 1: Bin(0.5), Bin(p), BetaBin(0.5, rho), BetaBin(p, rho),
   2-component mixture w Bin(p) + (1-w) Bin(1-p) (orientation unknown); rho per region (MLE) and pooled per
   family (shared rho, per-region p; grid MLE, bootstrap over cells). Joint comparison on log P(HA, HB):
   independent NB(HA) x NB(HB) (and ZINB) with offset log(raw_in + 0.5) vs NB(AT) x conditional split.
C. Dependence: Pearson-residual correlation of AT and raw depth after exposure adjustment (permutation
   within cluster x coverage tercile), elasticity of AT on region depth (robust Poisson GLM, offset log AT
   outside), AT / depth ratio by coverage tercile and cluster.
D. Recommendations; CIs are bootstrap over test cells, seed tuple (20261104, 15, 1).
Exploratory diagnostic of observation models; nothing here is copy-number truth.

Usage: python -I run.py [--check] [--report-only]
"""
import argparse
import html
import importlib.util
import json
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy import optimize, sparse, stats  # noqa: E402
from scipy.special import expit, gammaln  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

spec = importlib.util.spec_from_file_location("exp04", HERE.parent / "04_fixed_region_contrasts" / "run.py")
E4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(E4)
warnings.filterwarnings("ignore")

SEED = (C.SEED, 15)
SEED_BOOT = (C.SEED, 15, 1)
N_RANDOM = 60
N_BOOT = 1000
N_PERM = 300
MODELS = ["pois", "nb", "zip", "zinb", "hnb"]
ALL_MODELS = MODELS + ["nb_fs"]
MNAME = {"pois": "Poisson", "nb": "NB2", "zip": "ZIP", "zinb": "ZINB", "hnb": "hurdle-NB", "nb_fs": "NB2 free offset slope"}
RHO_GRID = np.r_[0.0, np.geomspace(1e-5, 0.95, 50)]
SM_ = ["b05", "b1", "bb05", "bb1", "mix"]
SNAME = {"b05": "Binomial p=0.5", "b1": "Binomial p", "bb05": "Beta-binomial p=0.5", "bb1": "Beta-binomial p,rho",
         "mix": "Mixture Bin(p)/Bin(1-p)", "bb_pooled": "Beta-binomial pooled family rho"}
CONFIGS = [("AT", "raw_in"), ("AT", "out"), ("depth", "out")]


# ------------------------------------------------------------------ count-model likelihoods
def softplus(x):
    return np.logaddexp(0, x)


def lp_pois(y, eta):
    return y * eta - np.exp(eta) - gammaln(y + 1)


def lp_nb(y, eta, la):
    r = np.exp(-la)                       # size = 1/alpha
    lr = -la
    den = np.logaddexp(lr, eta)
    return gammaln(y + r) - gammaln(r) - gammaln(y + 1) + r * (lr - den) + y * (eta - den)


def lp0_pois(eta):
    return -np.exp(eta)


def lp0_nb(eta, la):
    return -np.log1p(np.exp(la + eta)) / np.exp(la)


def _split(kind, par, k):
    la = par[k] if kind != "zip" else None
    g = par[k + 1:k + 3] if kind in ("zinb", "hnb") else par[k:k + 2]
    return la, g


def logp(kind, par, y, X, off, Z):
    """Per-cell log predictive probability. par = [beta (k) | log alpha | gamma (2)] by kind."""
    k = X.shape[1]
    eta = np.clip(X @ par[:k] + off, -30, 20)
    if kind == "pois":
        return lp_pois(y, eta)
    if kind in ("nb", "nb_fs"):
        return lp_nb(y, eta, par[k])
    la, g = _split(kind, par, k)
    z = Z @ g
    lpi, l1m = -softplus(-z), -softplus(z)
    if kind == "zip":
        base, l0 = lp_pois(y, eta), lp0_pois(eta)
    elif kind == "zinb":
        base, l0 = lp_nb(y, eta, la), lp0_nb(eta, la)
    else:   # hurdle NB
        l0 = np.minimum(lp0_nb(eta, la), -1e-12)
        pos = l1m + lp_nb(y, eta, la) - np.log(-np.expm1(l0))
        return np.where(y == 0, lpi, pos)
    return np.where(y == 0, np.logaddexp(lpi, l1m + l0), l1m + base)


def p_zero(kind, par, X, off, Z):
    k = X.shape[1]
    eta = np.clip(X @ par[:k] + off, -30, 20)
    if kind == "pois":
        return np.exp(lp0_pois(eta))
    if kind in ("nb", "nb_fs"):
        return np.exp(lp0_nb(eta, par[k]))
    la, g = _split(kind, par, k)
    pi = expit(Z @ g)
    if kind == "hnb":
        return pi
    l0 = lp0_pois(eta) if kind == "zip" else lp0_nb(eta, la)
    return pi + (1 - pi) * np.exp(l0)


def irls_pois(y, X, off, it=40):
    k = X.shape[1]
    b = np.zeros(k)
    b[0] = np.log((y.sum() + 0.5) / np.exp(off).sum())
    for _ in range(it):
        mu = np.exp(np.clip(X @ b + off, -30, 20))
        zz = X @ b + (y - mu) / np.maximum(mu, 1e-12)
        XW = X * mu[:, None]
        nb_ = np.linalg.solve(XW.T @ X + 1e-8 * np.eye(k), XW.T @ zz)
        done = np.max(np.abs(nb_ - b)) < 1e-8
        b = nb_
        if done:
            break
    return b


def fit(kind, y, X, off, Z):
    k = X.shape[1]
    b0 = irls_pois(y, X, off)
    if kind == "pois":
        return b0
    start, bnds = list(b0), [(-30, 30)] * k
    if kind != "zip":
        start.append(np.log(0.1)); bnds.append((-13, 4))
    if kind in ("zip", "zinb", "hnb"):
        start += [-2.5, 0.0] if kind != "hnb" else [-0.5, 0.0]
        bnds += [(-12, 12), (-12, 12)]
    kk = "nb" if kind == "nb_fs" else kind
    r = optimize.minimize(lambda p: -logp(kk, p, y, X, off, Z).sum(), np.array(start, float), method="L-BFGS-B",
                          bounds=bnds, options={"maxiter": 300})
    return r.x


def design(dm, extra=None):
    cols = [np.ones(len(dm))] + list(dm.T)
    if extra is not None:
        cols.append(extra)
    return np.column_stack(cols)


def score_models(y, off, D, tr, te, models=ALL_MODELS):
    mu_o, sd_o = off[tr].mean(), max(off[tr].std(), 1e-6)
    z = (off - mu_o) / sd_o
    X, Xfs = design(D), design(D, extra=z)
    Z = np.column_stack([np.ones(len(y)), z])
    out = {}
    for m in models:
        if m == "nb_fs":
            zo = np.zeros(len(y))
            p = fit("nb", y[tr], Xfs[tr], zo[tr], Z[tr])
            out[m] = (logp("nb", p, y[te], Xfs[te], zo[te], Z[te]), p, p_zero("nb", p, Xfs[te], zo[te], Z[te]))
        else:
            p = fit(m, y[tr], X[tr], off[tr], Z[tr])
            out[m] = (logp(m, p, y[te], X[te], off[te], Z[te]), p, p_zero(m, p, X[te], off[te], Z[te]))
    out["_offset_sd"] = sd_o
    return out


# ------------------------------------------------------------------ haplotype-split likelihoods
def lbb(k, n, p, rho):
    """beta-binomial log pmf, mean p, intra-class correlation rho (rho -> 0: binomial)."""
    p = np.clip(p, 1e-6, 1 - 1e-6)
    lc = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    if rho <= 1e-9:
        return lc + k * np.log(p) + (n - k) * np.log1p(-p)
    s = (1 - rho) / rho
    a, b = p * s, (1 - p) * s
    return lc + gammaln(k + a) + gammaln(n - k + b) - gammaln(n + a + b) + gammaln(a + b) - gammaln(a) - gammaln(b)


def lmix(k, n, p, w):
    return np.logaddexp(np.log(max(w, 1e-12)) + lbb(k, n, p, 0.0), np.log(max(1 - w, 1e-12)) + lbb(k, n, 1 - p, 0.0))


def fit_split(k, n):
    P = {"b05": ()}
    ph = float(np.clip(k.sum() / max(n.sum(), 1), 1e-4, 1 - 1e-4))
    P["b1"] = (ph,)
    r = optimize.minimize_scalar(lambda lr: -lbb(k, n, 0.5, np.exp(lr)).sum(), bounds=(np.log(1e-7), np.log(0.95)), method="bounded")
    P["bb05"] = (float(np.exp(r.x)),)
    r = optimize.minimize(lambda q: -lbb(k, n, expit(q[0]), np.exp(q[1])).sum(), [np.log(ph / (1 - ph)), np.log(1e-3)],
                          method="L-BFGS-B", bounds=[(-8, 8), (np.log(1e-7), np.log(0.95))])
    P["bb1"] = (float(expit(r.x[0])), float(np.exp(r.x[1])))
    best = None
    for p0, w0 in ((0.6, 0.5), (0.8, 0.5), (0.7, 0.3)):
        r = optimize.minimize(lambda q: -lmix(k, n, q[0], q[1]).sum(), [p0, w0], method="L-BFGS-B", bounds=[(0.5, 0.999), (0.0, 1.0)])
        if best is None or r.fun < best.fun:
            best = r
    P["mix"] = (float(best.x[0]), float(best.x[1]))
    return P


def split_logp(P, k, n):
    return {"b05": lbb(k, n, 0.5, 0.0), "b1": lbb(k, n, P["b1"][0], 0.0), "bb05": lbb(k, n, 0.5, P["bb05"][0]),
            "bb1": lbb(k, n, P["bb1"][0], P["bb1"][1]), "mix": lmix(k, n, *P["mix"])}


# ------------------------------------------------------------------ data
def load_allelic(pth, d, regs):
    rows, cols = [], []
    for j, r in enumerate(regs):
        rows += list(r["target"]); cols += [j] * len(r["target"])
    M = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(d.n_bins, len(regs)))
    cells = list(d.cells)
    res, tot, checks = {}, {}, {}
    for nm in ("HA", "HB", "AT"):
        t = pq.read_table(pth["allelic"] / f"{nm}.parquet")
        di = t.select(["chr", "start", "end"]).to_pandas()
        if "chr" not in di.columns:           # coordinates restored as the pandas index
            di = di.index.to_frame(index=False)
        ok = (di["chr"].astype(str).to_numpy() == d.seq).all() and (di["start"].to_numpy() == d.start).all() and (di["end"].to_numpy() == d.end).all()
        names = [c for c in t.column_names if c not in ("chr", "start", "end", "__index_level_0__")]
        checks[nm] = {"coords_match_anndata": bool(ok), "barcode_set_equal": set(names) == set(cells), "column_order_equal_anndata": names == cells}
        if not ok or set(names) != set(cells):
            raise SystemExit(f"alignment failure {nm}: {checks[nm]}")
        A = np.column_stack([t.column(c).to_numpy(zero_copy_only=False) for c in cells]).astype(np.float32)
        del t
        res[nm] = np.asarray((M.T @ A).T, dtype=np.float64)
        tot[nm] = A.sum(0, dtype=np.float64)
        checks[nm]["integer_valued"] = bool(np.all(np.mod(res[nm], 1) == 0))
        del A
    checks["AT_equals_HA_plus_HB_in_regions"] = bool(np.array_equal(res["AT"], res["HA"] + res["HB"]))
    return res, tot, checks


def random_uncalled(d, call, rng, n):
    ok = []
    for ch in dict.fromkeys(d.seq):
        ci_ = np.flatnonzero(d.seq == ch)
        for kk in range(E4.FLANK, len(ci_) - 25 - E4.FLANK):
            w = ci_[kk:kk + 25]
            if d.start[w[-1]] - d.start[w[0]] == 24 * 100_000 and np.all(call[w] == 1):
                ok.append((ch, kk))
    out = []
    for j in rng.choice(len(ok), n, replace=False):
        ch, kk = ok[j]
        ci_ = np.flatnonzero(d.seq == ch)
        out.append({"family": "random", "name": f"random15_{ch}:{d.start[ci_[kk]]/1e6:.1f}", "direction": "random",
                    "target": ci_[kk:kk + 25], "flank": np.array([], int)})
    return out, len(ok)


def boot_weights(rng, n, B):
    return rng.multinomial(n, np.full(n, 1.0 / n), size=B).astype(np.float32) / n


def ci3(v):
    return float(np.mean(v)), float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))


# ------------------------------------------------------------------ compute
def compute(args, out):
    t0 = time.time()
    pth = C.default_paths()
    rho06 = json.load(open(HERE.parent / "06_allelic_information" / "summary.json"))["rho_used"]
    d = C.H5AD(pth["h5ad"])
    lab = C.cell_clusters(d, pth)
    kseg = pd.read_csv(pth["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    call = np.ones(d.n_bins)
    for r in kseg[kseg.group == "all"].itertuples():
        call[r.bin_start:r.bin_end + 1] = r.cn_call
    regs0 = [r for r in E4.build_regions(d, pth) if r["family"] in ("saved", "called", "literature")]
    rr, n_eligible = random_uncalled(d, call, np.random.default_rng(SEED), N_RANDOM)
    regs = regs0 + rr
    if args.check:
        keep = []
        for f_, k in (("saved", 2), ("literature", 2), ("called", 3), ("random", 3)):
            keep += [r for r in regs if r["family"] == f_][:k]
        regs = keep
    R = len(regs)
    fam = np.array([r["family"] for r in regs])
    fams = list(dict.fromkeys(fam))
    print(f"regions {R}; eligible uncalled random windows {n_eligible}", flush=True)
    RAW, _, rawtot = E4.region_sums(d, regs, "layers/raw")
    al, altot, checks = load_allelic(pth, d, regs)
    AT, HA, HB, ATtot = al["AT"], al["HA"], al["HB"], altot["AT"]
    print("loaded", round(time.time() - t0), "s", checks, flush=True)

    clusters = sorted(set(lab))
    D = np.column_stack([(lab == c).astype(float) for c in clusters[1:]])
    terc = np.asarray(pd.qcut(rawtot, 3, labels=False))
    half = C.stratified_halves(lab, SEED)            # True = train (primary)
    te = ~half
    n = d.n_cells
    rng_b = np.random.default_rng(SEED_BOOT)
    NB_ = 60 if args.check else N_BOOT
    dirs = {"primary": (half, ~half), "reverse": (~half, half)}

    # ------------------------------------------------ A
    LS, ZP, PAR = {}, {}, {}
    for (mod, ofs) in CONFIGS:
        for dn, (tr, tt) in dirs.items():
            LS[(mod, ofs, dn)] = {m: np.zeros((tt.sum(), R), np.float32) for m in ALL_MODELS}
        ZP[(mod, ofs)] = {m: np.zeros((te.sum(), R), np.float32) for m in ALL_MODELS}
        for j in range(R):
            y = AT[:, j] if mod == "AT" else RAW[:, j]
            if mod == "AT" and ofs == "raw_in":
                off = np.log(RAW[:, j] + 0.5)
            elif mod == "AT":
                off = np.log(ATtot - AT[:, j] + 0.5)
            else:
                off = np.log(rawtot - RAW[:, j] + 0.5)
            for dn, (tr, tt) in dirs.items():
                res = score_models(y, off, D, tr, tt)
                for m in ALL_MODELS:
                    LS[(mod, ofs, dn)][m][:, j] = res[m][0]
                    if dn == "primary":
                        ZP[(mod, ofs)][m][:, j] = res[m][2]
                if dn == "primary":
                    k = D.shape[1] + 1
                    zt = np.column_stack([np.ones(tt.sum()), ((off - off[tr].mean()) / max(off[tr].std(), 1e-6))[tt]])
                    PAR[(mod, ofs, j)] = {"alpha_nb": float(np.exp(res["nb"][1][k])), "alpha_zinb": float(np.exp(res["zinb"][1][k])),
                                          "pi_zinb_mean": float(expit(zt @ res["zinb"][1][k + 1:k + 3]).mean()),
                                          "pi_zip_mean": float(expit(zt @ res["zip"][1][k:k + 2]).mean()),
                                          "pi_hurdle_zero_mean": float(expit(zt @ res["hnb"][1][k + 1:k + 3]).mean()),
                                          "nb_free_offset_slope": float(res["nb_fs"][1][k] / res["_offset_sd"]), "y_test": y[tt]}
            if (j + 1) % 25 == 0:
                print(f"  A {mod}/{ofs} region {j+1}/{R} {round(time.time()-t0)}s", flush=True)
        print(f"A done {mod}/{ofs} {round(time.time()-t0)}s", flush=True)

    rows = []
    for (mod, ofs) in CONFIGS:
        for j, r in enumerate(regs):
            P = PAR[(mod, ofs, j)]
            yte = P["y_test"]
            rec = {"modality": mod, "offset": ofs, "family": r["family"], "region": r["name"], "n_target_bins": len(r["target"]),
                   "mean_count_test": float(yte.mean()), "zero_frac_test": float((yte == 0).mean()),
                   "alpha_nb": P["alpha_nb"], "alpha_zinb": P["alpha_zinb"], "pi_zinb_mean": P["pi_zinb_mean"], "pi_zip_mean": P["pi_zip_mean"],
                   "pi_hurdle_zero_mean": P["pi_hurdle_zero_mean"], "nb_free_offset_slope": P["nb_free_offset_slope"]}
            for m in ALL_MODELS:
                rec[f"ls_{m}"] = float(LS[(mod, ofs, "primary")][m][:, j].mean())
                rec[f"ls_rev_{m}"] = float(LS[(mod, ofs, "reverse")][m][:, j].mean())
            for m in MODELS:
                rec[f"pz_{m}"] = float(ZP[(mod, ofs)][m][:, j].mean())
            for t in range(3):
                mt = terc[te] == t
                rec[f"obs_zero_t{t}"] = float((yte[mt] == 0).mean())
                for m in ("pois", "nb", "zinb", "hnb"):
                    rec[f"pz_{m}_t{t}"] = float(ZP[(mod, ofs)][m][mt, j].mean())
            n0 = (yte == 0).sum()
            for m in ("pois", "nb", "zinb", "hnb"):
                pz_ = ZP[(mod, ofs)][m][:, j].astype(np.float64)
                rec[f"zz_{m}"] = float((n0 - pz_.sum()) / np.sqrt(max((pz_ * (1 - pz_)).sum(), 1e-12)))
                for t in range(3):
                    mt = terc[te] == t
                    rec[f"zz_{m}_t{t}"] = float(((yte[mt] == 0).sum() - pz_[mt].sum()) / np.sqrt(max((pz_[mt] * (1 - pz_[mt])).sum(), 1e-12)))
            rec["best_model"] = max(MODELS, key=lambda m: rec[f"ls_{m}"])
            rows.append(rec)
    at_tab = pd.DataFrame(rows)
    at_tab.to_csv(out / "at_models.tsv", sep="\t", index=False)

    W = boot_weights(rng_b, te.sum(), NB_)
    fr = []
    for (mod, ofs) in CONFIGS:
        for fname in fams + ["all"]:
            sel = np.ones(R, bool) if fname == "all" else fam == fname
            base = LS[(mod, ofs, "primary")]["nb"][:, sel]
            baser = LS[(mod, ofs, "reverse")]["nb"][:, sel]
            for m in ALL_MODELS:
                if m == "nb":
                    continue
                dl = LS[(mod, ofs, "primary")][m][:, sel] - base
                mean, lo, hi = ci3((W @ dl).mean(1))
                fr.append({"modality": mod, "offset": ofs, "family": fname, "model": m, "regions": int(sel.sum()),
                           "d_logscore_vs_nb": float(dl.mean()), "ci_lo": lo, "ci_hi": hi,
                           "d_logscore_reverse_split": float((LS[(mod, ofs, "reverse")][m][:, sel] - baser).mean()),
                           "frac_regions_better_than_nb": float(np.mean(dl.mean(0) > 0))})
    fam_tab = pd.DataFrame(fr)
    fam_tab.to_csv(out / "at_family_summary.tsv", sep="\t", index=False)
    zr = []
    for (mod, ofs) in CONFIGS:
        for fname in fams:
            sub = at_tab[(at_tab.modality == mod) & (at_tab.offset == ofs) & (at_tab.family == fname)]
            for t in range(3):
                z = {"modality": mod, "offset": ofs, "family": fname, "coverage_tercile": t, "regions": len(sub),
                     "obs_zero_frac": sub[f"obs_zero_t{t}"].mean()}
                for m in ("pois", "nb", "zinb", "hnb"):
                    z[f"pred_zero_{m}"] = sub[f"pz_{m}_t{t}"].mean()
                zr.append(z)
    pd.DataFrame(zr).to_csv(out / "zero_by_tercile.tsv", sep="\t", index=False)
    print("A tables", round(time.time() - t0), "s", flush=True)

    # ------------------------------------------------ B
    SPL = {dn: {m: np.zeros((tt.sum(), R), np.float32) for m in SM_} for dn, (tr, tt) in dirs.items()}
    sp_rows = []
    for j, r in enumerate(regs):
        for dn, (tr, tt) in dirs.items():
            m1 = AT[tr, j] >= 1
            P = fit_split(HB[tr, j][m1], AT[tr, j][m1])
            lp = split_logp(P, HB[tt, j], AT[tt, j])
            m2 = AT[tt, j] >= 1
            for m in SM_:
                SPL[dn][m][:, j] = np.where(m2, lp[m], 0.0)
        mall = AT[:, j] >= 1
        Pa = fit_split(HB[mall, j], AT[mall, j])
        a2 = mall & (AT[:, j] >= 2)
        sp = {"family": r["family"], "region": r["name"], "cells_AT_ge1": int(mall.sum()), "cells_AT_ge2": int(a2.sum()),
              "mean_AT_given_ge1": float(AT[mall, j].mean()), "max_AT": int(AT[:, j].max()),
              "p_hat_HB_over_AT": Pa["b1"][0], "rho_mle_p05": Pa["bb05"][0], "rho_mle": Pa["bb1"][1], "p_bb": Pa["bb1"][0],
              "mix_p": Pa["mix"][0], "mix_w": Pa["mix"][1], "rho_over_exp06": Pa["bb1"][1] / rho06}
        llb = lbb(HB[mall, j], AT[mall, j], Pa["b1"][0], 0.0).sum()
        llbb = lbb(HB[mall, j], AT[mall, j], Pa["bb1"][0], Pa["bb1"][1]).sum()
        sp["lrt_bb_vs_binom"] = float(max(0.0, 2 * (llbb - llb)))
        sp["p_lrt"] = float(0.5 * stats.chi2.sf(sp["lrt_bb_vs_binom"], 1)) if sp["lrt_bb_vs_binom"] > 0 else 1.0
        if a2.any():
            ph = Pa["b1"][0]
            kk, nn = HB[a2, j], AT[a2, j]
            sp["pearson_disp_AT_ge2"] = float(np.mean((kk - nn * ph) ** 2 / (nn * ph * (1 - ph))))
        for dn, pre in (("primary", "ls_"), ("reverse", "ls_rev_")):
            tt = dirs[dn][1]
            inf = (AT[tt, j] >= 1).sum()
            for m in SM_:
                sp[f"{pre}{m}"] = float(SPL[dn][m][:, j].sum() / max(inf, 1))   # mean over held-out cells with AT >= 1
        sp["test_cells_AT_ge1"] = int((AT[te, j] >= 1).sum())
        sp_rows.append(sp)

    # cell-level Pearson dispersion of HB | AT by AT stratum (comparable to the bin-level variance inflation of exp 06)
    drows = []
    strata_AT = [("2", 2, 2), ("3-4", 3, 4), ("5-9", 5, 9), ("10-19", 10, 19), (">=20", 20, 10 ** 9)]
    for fname in fams:
        idx = np.flatnonzero(fam == fname)
        for lab_, lo_, hi_ in strata_AT:
            num, cnt_, nsum = 0.0, 0, 0.0
            for j in idx:
                ph_ = sp_rows[j]["p_hat_HB_over_AT"]
                m_ = (AT[:, j] >= lo_) & (AT[:, j] <= hi_)
                if m_.any():
                    kk, nn = HB[m_, j], AT[m_, j]
                    num += float(np.sum((kk - nn * ph_) ** 2 / (nn * ph_ * (1 - ph_)))); cnt_ += int(m_.sum()); nsum += float(nn.sum())
            if cnt_:
                phi_ = num / cnt_
                nbar = nsum / cnt_
                drows.append({"family": fname, "AT_stratum": lab_, "cell_region_pairs": cnt_, "mean_AT": nbar, "pearson_dispersion": phi_,
                              "implied_rho": (phi_ - 1) / (nbar - 1) if nbar > 1 else np.nan})
    pd.DataFrame(drows).to_csv(out / "split_dispersion_by_AT.tsv", sep="\t", index=False)

    pooled, SPL_pooled = [], np.zeros((te.sum(), R), np.float32)
    for fname in fams:
        idx = np.flatnonzero(fam == fname)
        cell_i, Lg = [], []
        for j in idx:
            mall = AT[:, j] >= 1
            ph = float(np.clip(HB[mall, j].sum() / max(AT[mall, j].sum(), 1), 1e-4, 1 - 1e-4))
            Lg.append(np.column_stack([lbb(HB[mall, j], AT[mall, j], ph, g) for g in RHO_GRID]))
            cell_i.append(np.flatnonzero(mall))
        cell_i = np.concatenate(cell_i); Lg = np.vstack(Lg)
        tot_ll = Lg.sum(0)
        gb = int(np.argmax(tot_ll))
        Wb = rng_b.multinomial(n, np.full(n, 1.0 / n), size=60 if args.check else 400).astype(np.float64)
        bsr = np.array([RHO_GRID[int(np.argmax(wb[cell_i] @ Lg))] for wb in Wb])
        lrt = 2 * (tot_ll[gb] - tot_ll[0])
        rho_f = float(RHO_GRID[gb])
        pooled.append({"family": fname, "regions": len(idx), "region_cell_rows": len(cell_i), "rho_pooled": rho_f,
                       "rho_ci_lo": float(np.quantile(bsr, 0.025)), "rho_ci_hi": float(np.quantile(bsr, 0.975)),
                       "median_region_rho_mle": float(np.median([sp_rows[j]["rho_mle"] for j in idx])),
                       "lrt_vs_binomial": float(lrt), "p_lrt_half_chi2": float(0.5 * stats.chi2.sf(lrt, 1)) if lrt > 0 else 1.0,
                       "exp06_bin_level_rho": rho06, "ratio_to_exp06": rho_f / rho06})
        for j in idx:      # held-out score with family rho (train-half region p)
            m1 = AT[half, j] >= 1
            ph = float(np.clip(HB[half, j][m1].sum() / max(AT[half, j][m1].sum(), 1), 1e-4, 1 - 1e-4))
            m2 = AT[te, j] >= 1
            v = np.where(m2, lbb(HB[te, j], AT[te, j], ph, rho_f), 0.0)
            SPL_pooled[:, j] = v
            sp_rows[j]["ls_bb_pooled"] = float(v.sum() / max(m2.sum(), 1))
    pooled_tab = pd.DataFrame(pooled)
    pd.DataFrame(sp_rows).to_csv(out / "split_models.tsv", sep="\t", index=False)
    pooled_tab.to_csv(out / "split_pooled_rho.tsv", sep="\t", index=False)

    W2 = boot_weights(rng_b, te.sum(), NB_)
    spf = []
    for fname in fams + ["all"]:
        sel = np.ones(R, bool) if fname == "all" else fam == fname
        ninf = max(int((AT[te][:, sel] >= 1).sum()), 1)
        base = SPL["primary"]["bb1"][:, sel]
        for m in SM_ + ["bb_pooled"]:
            if m == "bb1":
                continue
            arr = SPL_pooled[:, sel] if m == "bb_pooled" else SPL["primary"][m][:, sel]
            dl = arr - base
            mean, lo, hi = ci3((W2 @ dl).mean(1))
            spf.append({"family": fname, "model": m, "regions": int(sel.sum()), "informative_test_cell_region_pairs": ninf,
                        "d_logscore_per_cell_vs_bb1": float(dl.mean()), "ci_lo": lo, "ci_hi": hi,
                        "d_logscore_per_informative_pair": float(dl.sum() / ninf),
                        "d_logscore_reverse_split_per_cell": float(((SPL["reverse"][m][:, sel] if m != "bb_pooled" else SPL["reverse"]["bb1"][:, sel])
                                                                   - SPL["reverse"]["bb1"][:, sel]).mean())})
    pd.DataFrame(spf).to_csv(out / "split_family_summary.tsv", sep="\t", index=False)
    print("B split", round(time.time() - t0), "s", flush=True)

    JM = {k: np.zeros((te.sum(), R), np.float32) for k in
          ["NB_HA+NB_HB", "ZINB_HA+ZINB_HB", "NB_AT+Binom", "NB_AT+BB", "ZINB_AT+BB", "ZINB_AT+mixture"]}
    for j in range(R):
        off = np.log(RAW[:, j] + 0.5)
        ra = score_models(HA[:, j], off, D, half, te, models=["nb", "zinb"])
        rb = score_models(HB[:, j], off, D, half, te, models=["nb", "zinb"])
        JM["NB_HA+NB_HB"][:, j] = ra["nb"][0] + rb["nb"][0]
        JM["ZINB_HA+ZINB_HB"][:, j] = ra["zinb"][0] + rb["zinb"][0]
        nbA, zbA = LS[("AT", "raw_in", "primary")]["nb"][:, j], LS[("AT", "raw_in", "primary")]["zinb"][:, j]
        JM["NB_AT+Binom"][:, j] = nbA + SPL["primary"]["b1"][:, j]
        JM["NB_AT+BB"][:, j] = nbA + SPL["primary"]["bb1"][:, j]
        JM["ZINB_AT+BB"][:, j] = zbA + SPL["primary"]["bb1"][:, j]
        JM["ZINB_AT+mixture"][:, j] = zbA + SPL["primary"]["mix"][:, j]
    jr = []
    for fname in fams + ["all"]:
        sel = np.ones(R, bool) if fname == "all" else fam == fname
        base = JM["NB_AT+BB"][:, sel]
        for k_, v in JM.items():
            if k_ == "NB_AT+BB":
                continue
            dl = v[:, sel] - base
            mean, lo, hi = ci3((W2 @ dl).mean(1))
            jr.append({"family": fname, "joint_model": k_, "regions": int(sel.sum()), "d_joint_logscore_vs_NB_AT_plus_BB": float(dl.mean()),
                       "ci_lo": lo, "ci_hi": hi, "mean_joint_logscore_NB_AT_plus_BB": float(base.mean())})
    pd.DataFrame(jr).to_csv(out / "split_joint_HA_HB.tsv", sep="\t", index=False)
    print("B joint", round(time.time() - t0), "s", flush=True)

    # ------------------------------------------------ C
    X0 = design(D)
    Z0 = np.column_stack([np.ones(n), np.zeros(n)])
    sid = pd.factorize(pd.Series(lab.astype(str)) + "|" + pd.Series(terc.astype(str)))[0]
    strata_idx = [np.flatnonzero(sid == s) for s in range(sid.max() + 1)]
    rng_p = np.random.default_rng(SEED_BOOT)
    npm = 50 if args.check else N_PERM
    T2 = np.column_stack([(terc == t).astype(float) for t in (1, 2)])
    dep = []
    for j, r in enumerate(regs):
        offA, offR = np.log(ATtot - AT[:, j] + 0.5), np.log(rawtot - RAW[:, j] + 0.5)
        k = X0.shape[1]
        pa, pr = fit("nb", AT[:, j], X0, offA, Z0), fit("nb", RAW[:, j], X0, offR, Z0)
        muA, aA = np.exp(np.clip(X0 @ pa[:k] + offA, -30, 20)), np.exp(pa[k])
        muR, aR = np.exp(np.clip(X0 @ pr[:k] + offR, -30, 20)), np.exp(pr[k])
        eA = (AT[:, j] - muA) / np.sqrt(muA + aA * muA ** 2)
        eR = (RAW[:, j] - muR) / np.sqrt(muR + aR * muR ** 2)
        rp = float(np.corrcoef(eA, eR)[0, 1])
        cnt = 0
        for _ in range(npm):
            e2 = eR.copy()
            for ii in strata_idx:
                e2[ii] = e2[rng_p.permutation(ii)]
            cnt += abs(np.corrcoef(eA, e2)[0, 1]) >= abs(rp)
        xl = np.log(RAW[:, j] + 0.5)
        try:
            ge = sm.GLM(AT[:, j], np.column_stack([X0, xl - xl.mean()]), family=sm.families.Poisson(), offset=offA).fit(cov_type="HC0")
            el, el_se = float(ge.params[-1]), float(ge.bse[-1])
        except Exception:   # noqa: BLE001
            el, el_se = np.nan, np.nan
        rec = {"family": r["family"], "region": r["name"], "AT_total": float(AT[:, j].sum()), "raw_total": float(RAW[:, j].sum()),
               "AT_per_raw": float(AT[:, j].sum() / max(RAW[:, j].sum(), 1)),
               "resid_corr_pearson": rp, "resid_corr_spearman": float(stats.spearmanr(eA, eR)[0]), "p_perm_resid_corr": (cnt + 1) / (npm + 1),
               "elasticity_AT_on_depth": el, "elasticity_se_robust": el_se}
        for t in range(3):
            mt = terc == t
            rec[f"AT_per_raw_t{t}"] = float(AT[mt, j].sum() / max(RAW[mt, j].sum(), 1))
        rec["log_ratio_low_over_high"] = float(np.log(max(rec["AT_per_raw_t0"], 1e-9) / max(rec["AT_per_raw_t2"], 1e-9)))
        cr = [AT[lab == c, j].sum() / max(RAW[lab == c, j].sum(), 1) for c in clusters]
        rec["AT_per_raw_cluster_min"], rec["AT_per_raw_cluster_max"] = float(min(cr)), float(max(cr))
        offr = np.log(RAW[:, j] + 0.5)
        for nm, cols, base_c in (("tercile", T2, D), ("cluster", D, T2)):
            Xh = np.column_stack([np.ones(n), base_c, cols])
            try:
                gh = sm.GLM(AT[:, j], Xh, family=sm.families.Poisson(), offset=offr).fit(scale="X2")
                idx = list(range(Xh.shape[1] - cols.shape[1], Xh.shape[1]))
                b, V = gh.params[idx], gh.cov_params()[np.ix_(idx, idx)]
                rec[f"p_ratio_{nm}"] = float(stats.f.sf(float(b @ np.linalg.solve(V, b)) / len(idx), len(idx), gh.df_resid))
            except Exception:   # noqa: BLE001
                rec[f"p_ratio_{nm}"] = np.nan
        dep.append(rec)
    pd.DataFrame(dep).to_csv(out / "dependence.tsv", sep="\t", index=False)
    print("C", round(time.time() - t0), "s", flush=True)

    meta = {"regions": [{"family": r["family"], "name": r["name"], "n_target_bins": int(len(r["target"]))} for r in regs],
            "alignment_checks": checks, "n_eligible_uncalled_random_windows": n_eligible, "rho_exp06": rho06,
            "n_cells": n, "n_train_cells": int(half.sum()), "n_test_cells": int(te.sum()), "clusters": [str(c) for c in clusters],
            "bootstrap_reps_models": NB_, "bootstrap_reps_rho": 60 if args.check else 400, "permutations": npm,
            "compute_seconds": time.time() - t0, "check_mode": bool(args.check),
            "frac_cells_AT_zero_median_region": float(np.median((AT == 0).mean(0)))}
    C.write_json(out / "_compute_meta.json", meta)
    return meta


# ------------------------------------------------------------------ report
def report(out, args):
    # implemented in the report stage (reads only the TSV/JSON outputs written above)
    raise NotImplementedError


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--out", default=str(HERE))
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if not args.report_only:
        compute(args, out)
    report(out, args)
    print("total", round(time.time() - t0), "s")


if __name__ == "__main__":
    main()
