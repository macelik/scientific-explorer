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
SM_ = ["b05", "b1", "bb05", "bb1", "bbphi", "mix"]
SNAME = {"b05": "Binomial p=0.5", "b1": "Binomial p", "bb05": "Beta-binomial p=0.5", "bb1": "Beta-binomial p,rho",
         "mix": "Mixture Bin(p)/Bin(1-p)", "bbphi": "Beta-binomial constant variance inflation", "bb_pooled": "Beta-binomial pooled family rho"}
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


def lbb_phi(k, n, p, phi):
    """beta-binomial whose rho depends on the cell's n so that the variance inflation 1 + (n-1) rho = phi is constant."""
    p = np.clip(p, 1e-6, 1 - 1e-6)
    rho = np.where(n > 1, np.clip((phi - 1.0) / np.maximum(n - 1.0, 1.0), 0.0, 0.95), 0.0)
    lc = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    binom = lc + k * np.log(p) + (n - k) * np.log1p(-p)
    rs = np.maximum(rho, 1e-6)
    s = (1 - rs) / rs
    a, b = p * s, (1 - p) * s
    bb = lc + gammaln(k + a) + gammaln(n - k + b) - gammaln(n + a + b) + gammaln(a + b) - gammaln(a) - gammaln(b)
    return np.where(rho > 1e-6, bb, binom)


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
    r = optimize.minimize(lambda q: -lbb_phi(k, n, expit(q[0]), 1.0 + np.exp(q[1])).sum(), [np.log(ph / (1 - ph)), np.log(0.3)],
                          method="L-BFGS-B", bounds=[(-8, 8), (np.log(1e-4), np.log(30.0))])
    P["bbphi"] = (float(expit(r.x[0])), float(1.0 + np.exp(r.x[1])))
    best = None
    for p0, w0 in ((0.6, 0.5), (0.8, 0.5), (0.7, 0.3)):
        r = optimize.minimize(lambda q: -lmix(k, n, q[0], q[1]).sum(), [p0, w0], method="L-BFGS-B", bounds=[(0.5, 0.999), (0.0, 1.0)])
        if best is None or r.fun < best.fun:
            best = r
    P["mix"] = (float(best.x[0]), float(best.x[1]))
    return P


def split_logp(P, k, n):
    return {"b05": lbb(k, n, 0.5, 0.0), "b1": lbb(k, n, P["b1"][0], 0.0), "bb05": lbb(k, n, 0.5, P["bb05"][0]),
            "bb1": lbb(k, n, P["bb1"][0], P["bb1"][1]),
            "bbphi": lbb_phi(k, n, P["bbphi"][0], P["bbphi"][1]), "mix": lmix(k, n, *P["mix"])}


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
    orow = []      # offset choice for AT: NB2 with raw_in offset minus NB2 with 'out' offset (same train/test split)
    for fname in fams + ["all"]:
        sel = np.ones(R, bool) if fname == "all" else fam == fname
        for m in ("nb", "zinb"):
            dl = LS[("AT", "raw_in", "primary")][m][:, sel] - LS[("AT", "out", "primary")][m][:, sel]
            mean, lo, hi = ci3((W @ dl).mean(1))
            orow.append({"family": fname, "model": m, "regions": int(sel.sum()), "d_logscore_raw_in_minus_out": float(dl.mean()),
                         "ci_lo": lo, "ci_hi": hi, "frac_regions_raw_in_better": float(np.mean(dl.mean(0) > 0))})
    pd.DataFrame(orow).to_csv(out / "offset_comparison.tsv", sep="\t", index=False)
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
              "phi_mle": Pa["bbphi"][1], "mix_p": Pa["mix"][0], "mix_w": Pa["mix"][1], "rho_over_exp06": Pa["bb1"][1] / rho06}
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
FAM_ORDER = ["saved", "literature", "called", "random", "all"]
CFG_LABEL = {("AT", "raw_in"): "AT, offset log(raw depth in region + 0.5)",
             ("AT", "out"): "AT, offset log(AT outside region + 0.5)",
             ("depth", "out"): "depth, offset log(raw outside region + 0.5)"}
EQUIV = 0.002   # judgement margin (nats per cell-region): |difference| below this is treated as practically equivalent


def fci(d_, lo, hi, nd=4):
    return f"{d_:+.{nd}f} [{lo:+.{nd}f}, {hi:+.{nd}f}]"


def report(out, args):
    rd = lambda n_: pd.read_csv(out / n_, sep="\t")  # noqa: E731
    meta = json.load(open(out / "_compute_meta.json"))
    at, fam_t, zt = rd("at_models.tsv"), rd("at_family_summary.tsv"), rd("zero_by_tercile.tsv")
    offc, sp, pooled = rd("offset_comparison.tsv"), rd("split_models.tsv"), rd("split_pooled_rho.tsv")
    spf, disp, joint, dep = rd("split_family_summary.tsv"), rd("split_dispersion_by_AT.tsv"), rd("split_joint_HA_HB.tsv"), rd("dependence.tsv")
    rho06 = meta["rho_exp06"]
    ov06 = pd.read_csv(HERE.parent / "06_allelic_information" / "overdispersion_estimate.tsv", sep="\t")
    paths = C.default_paths()

    def G(mod, ofs, fam_, model):
        return fam_t[(fam_t.modality == mod) & (fam_t.offset == ofs) & (fam_t.family == fam_) & (fam_t.model == model)].iloc[0]

    def S(fam_, model):
        return spf[(spf.family == fam_) & (spf.model == model)].iloc[0]

    def J(fam_, model):
        return joint[(joint.family == fam_) & (joint.joint_model == model)].iloc[0]

    # ---------------- A tables
    rows = []
    for (mod, ofs) in CONFIGS:
        for f_ in FAM_ORDER:
            r = {"model class": CFG_LABEL[(mod, ofs)], "family": f_, "regions": int(G(mod, ofs, f_, "pois").regions)}
            for m in ["pois", "zip", "zinb", "hnb", "nb_fs"]:
                g = G(mod, ofs, f_, m)
                r[f"{MNAME[m]} minus NB2"] = fci(g.d_logscore_vs_nb, g.ci_lo, g.ci_hi)
            r["% regions ZINB > NB2"] = round(100 * float(G(mod, ofs, f_, "zinb").frac_regions_better_than_nb))
            rows.append(r)
    tabA1 = pd.DataFrame(rows)

    zc = []
    for (mod, ofs) in CONFIGS:
        for f_ in FAM_ORDER[:-1]:
            s_ = at[(at.modality == mod) & (at.offset == ofs) & (at.family == f_)]
            zc.append({"model class": CFG_LABEL[(mod, ofs)], "family": f_, "regions": len(s_),
                       "median test zero fraction": s_.zero_frac_test.median(), "median NB2 predicted zero fraction": s_.pz_nb.median(),
                       "median Poisson predicted": s_.pz_pois.median(), "median ZINB predicted": s_.pz_zinb.median(),
                       "median z (NB2)": s_.zz_nb.median(), "% |z|>2 (NB2)": 100 * (s_.zz_nb.abs() > 2).mean(),
                       "median z (ZINB)": s_.zz_zinb.median(), "% |z|>2 (ZINB)": 100 * (s_.zz_zinb.abs() > 2).mean(),
                       "median NB2 alpha": s_.alpha_nb.median(), "median ZINB mean pi": s_.pi_zinb_mean.median()})
    tabA2 = pd.DataFrame(zc)

    zz = zt[(zt.modality == "AT")].copy()
    for m in ("pois", "nb", "zinb", "hnb"):
        zz[f"{MNAME[m]} pred - obs (pts)"] = 100 * (zz[f"pred_zero_{m}"] - zz.obs_zero_frac)
    zz["observed zero (%)"] = 100 * zz.obs_zero_frac
    tabA3 = zz[["offset", "family", "coverage_tercile", "regions", "observed zero (%)"] + [f"{MNAME[m]} pred - obs (pts)" for m in ("pois", "nb", "zinb", "hnb")]]
    tabA3 = tabA3.rename(columns={"coverage_tercile": "coverage tercile (0 = low)"})

    off_t = offc.copy()
    off_t["raw_in minus out"] = [fci(a, b, c) for a, b, c in zip(off_t.d_logscore_raw_in_minus_out, off_t.ci_lo, off_t.ci_hi)]
    off_t = off_t[["family", "model", "regions", "raw_in minus out", "frac_regions_raw_in_better"]]

    # ---------------- A figures
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), sharex=False)
    cols = {"zinb": "#2f5d8a", "hnb": "#b5651d", "nb_fs": "#6b8e23"}
    for ax, (mod, ofs) in zip(axs, CONFIGS):
        for i, f_ in enumerate(FAM_ORDER):
            for k, m in enumerate(["zinb", "hnb", "nb_fs"]):
                g = G(mod, ofs, f_, m)
                y = i + (k - 1) * 0.22
                ax.errorbar(g.d_logscore_vs_nb, y, xerr=[[g.d_logscore_vs_nb - g.ci_lo], [g.ci_hi - g.d_logscore_vs_nb]], fmt="o", ms=4,
                            color=cols[m], label=MNAME[m] if i == 0 else None, capsize=2)
        ax.axvline(0, color="grey", lw=.8); ax.axvspan(-EQUIV, EQUIV, color="grey", alpha=.15)
        ax.set_yticks(range(len(FAM_ORDER))); ax.set_yticklabels(FAM_ORDER); ax.set_title(CFG_LABEL[(mod, ofs)], fontsize=8)
        ax.set_xlabel("held-out log score minus NB2 (nats per cell-region)")
    axs[0].legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    figA1 = C.fig_to_b64(fig); plt.close(fig)

    fig, axs = plt.subplots(1, 4, figsize=(13, 3.4), sharey=True)
    for ax, f_ in zip(axs, FAM_ORDER[:-1]):
        s_ = zt[(zt.modality == "AT") & (zt.offset == "raw_in") & (zt.family == f_)]
        for m, c in (("pois", "#999"), ("nb", "#2f5d8a"), ("zinb", "#b5651d"), ("hnb", "#6b8e23")):
            ax.plot(s_.coverage_tercile, 100 * (s_[f"pred_zero_{m}"] - s_.obs_zero_frac), "o-", color=c, label=MNAME[m])
        ax.axhline(0, color="k", lw=.7); ax.set_title(f"{f_} (observed zero {100*s_.obs_zero_frac.mean():.0f}% avg)", fontsize=9)
        ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["low", "mid", "high"]); ax.set_xlabel("cell coverage tercile")
    axs[0].set_ylabel("predicted - observed zero fraction (points)"); axs[0].legend(fontsize=7)
    fig.tight_layout()
    figA2 = C.fig_to_b64(fig); plt.close(fig)

    # ---------------- B tables
    brow = []
    for f_ in FAM_ORDER:
        r = {"family": f_, "regions": int(S(f_, "b1").regions), "informative test cell-region pairs": int(S(f_, "b1").informative_test_cell_region_pairs)}
        for m in ["b05", "b1", "bb05", "bbphi", "mix", "bb_pooled"]:
            g = S(f_, m)
            r[f"{SNAME[m]} minus BB(p,rho)"] = fci(g.d_logscore_per_cell_vs_bb1, g.ci_lo, g.ci_hi)
        brow.append(r)
    tabB1 = pd.DataFrame(brow)

    reg_sp = sp.copy()
    srow = []
    for f_ in FAM_ORDER[:-1]:
        s_ = reg_sp[reg_sp.family == f_]
        srow.append({"family": f_, "regions": len(s_), "median cells with AT>=1": s_.cells_AT_ge1.median(), "median mean AT | AT>=1": s_.mean_AT_given_ge1.median(),
                     "median region rho (BB MLE)": s_.rho_mle.median(), "IQR rho lo": s_.rho_mle.quantile(.25), "IQR rho hi": s_.rho_mle.quantile(.75),
                     "median region phi (constant-inflation MLE)": s_.phi_mle.median(),
                     "median Pearson dispersion, AT>=2": s_.pearson_disp_AT_ge2.median(),
                     "% regions LRT BB vs Bin p<0.01": 100 * (s_.p_lrt < 0.01).mean()})
    tabB2 = pd.DataFrame(srow)
    pt = pooled.copy()
    pt["rho pooled [bootstrap 95% CI]"] = [f"{a:.4f} [{b:.4f}, {c:.4f}]" for a, b, c in zip(pt.rho_pooled, pt.rho_ci_lo, pt.rho_ci_hi)]
    tabB3 = pt[["family", "regions", "region_cell_rows", "rho pooled [bootstrap 95% CI]", "median_region_rho_mle", "lrt_vs_binomial", "exp06_bin_level_rho", "ratio_to_exp06"]]
    spear = stats.spearmanr(np.log(reg_sp.mean_AT_given_ge1), np.log(reg_sp.rho_mle.clip(lower=1e-6)))
    ov = ov06.rename(columns={"AT_stratum": "exp06 bin-level AT stratum", "bins": "bins", "median_AT": "median AT", "variance_inflation": "bin-level variance inflation", "rho_hat": "bin-level rho"})
    dsp = disp.copy()
    allfam = disp.groupby("AT_stratum", sort=False).apply(lambda g: pd.Series({"family": "pooled", "cell_region_pairs": g.cell_region_pairs.sum(),
                                                                           "mean_AT": np.average(g.mean_AT, weights=g.cell_region_pairs),
                                                                           "pearson_dispersion": np.average(g.pearson_dispersion, weights=g.cell_region_pairs),
                                                                           "implied_rho": np.nan})).reset_index()
    allfam["implied_rho"] = (allfam.pearson_dispersion - 1) / (allfam.mean_AT - 1)
    dsp = pd.concat([dsp, allfam[dsp.columns]], ignore_index=True)

    jt = joint.copy()
    jt["minus NB(AT) x BB"] = [fci(a, b, c) for a, b, c in zip(jt.d_joint_logscore_vs_NB_AT_plus_BB, jt.ci_lo, jt.ci_hi)]
    tabB5 = jt.pivot(index="family", columns="joint_model", values="minus NB(AT) x BB").reindex(FAM_ORDER).reset_index()

    # ---------------- B figure
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4))
    cm = {"saved": "#b5651d", "literature": "#6b8e23", "called": "#2f5d8a", "random": "#888"}
    for f_, c in cm.items():
        s_ = reg_sp[reg_sp.family == f_]
        axs[0].scatter(s_.mean_AT_given_ge1, s_.rho_mle.clip(lower=1e-5), s=14, color=c, label=f_, alpha=.8)
    nn = np.geomspace(1.2, 60, 50)
    axs[0].plot(nn, np.clip((1.4 - 1) / (nn - 1), 0, 1), "k--", lw=1, label="constant inflation phi = 1.4")
    axs[0].axhline(rho06, color="r", lw=1); axs[0].text(1.3, rho06 * 1.15, "exp 06 bin-level rho = 0.0038", color="r", fontsize=8)
    axs[0].set_xscale("log"); axs[0].set_yscale("log"); axs[0].set_xlabel("mean AT per cell (cells with AT >= 1)"); axs[0].set_ylabel("cell-level beta-binomial rho (region MLE)")
    axs[0].legend(fontsize=7)
    for f_, c in cm.items():
        s_ = disp[disp.family == f_]
        axs[1].plot(s_.mean_AT, s_.pearson_dispersion, "o-", color=c, label=f_)
    axs[1].scatter(ov06.median_AT, ov06.variance_inflation, marker="s", color="r", label="exp 06 bins (summed over cells)")
    axs[1].axhline(1, color="k", lw=.7); axs[1].set_xscale("log"); axs[1].set_xlabel("mean AT of the unit (cell-region or bin)")
    axs[1].set_ylabel("variance inflation of HB | AT (Pearson dispersion)"); axs[1].legend(fontsize=7)
    fig.tight_layout()
    figB = C.fig_to_b64(fig); plt.close(fig)

    # ---------------- C
    crow = []
    for f_ in FAM_ORDER[:-1]:
        s_ = dep[dep.family == f_]
        crow.append({"family": f_, "regions": len(s_), "median residual corr (Pearson)": s_.resid_corr_pearson.median(), "min": s_.resid_corr_pearson.min(),
                     "max": s_.resid_corr_pearson.max(), "% permutation p<0.01": 100 * (s_.p_perm_resid_corr < 0.01).mean(),
                     "median elasticity of AT on depth": s_.elasticity_AT_on_depth.median(), "elasticity min": s_.elasticity_AT_on_depth.min(),
                     "elasticity max": s_.elasticity_AT_on_depth.max(), "median AT per raw fragment": s_.AT_per_raw.median(),
                     "median log(ratio low/high tercile)": s_.log_ratio_low_over_high.median(),
                     "% ratio varies by tercile (p<0.01)": 100 * (s_.p_ratio_tercile < 0.01).mean(),
                     "% ratio varies by cluster (p<0.01)": 100 * (s_.p_ratio_cluster < 0.01).mean()})
    tabC = pd.DataFrame(crow)
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
    for k, (col, ttl) in enumerate([("resid_corr_pearson", "Pearson-residual correlation, AT vs depth"), ("elasticity_AT_on_depth", "elasticity of AT on region depth (1 = proportional)"),
                                     ("log_ratio_low_over_high", "log(AT/depth, low tercile / high tercile)")]):
        for i, f_ in enumerate(FAM_ORDER[:-1]):
            v = dep[dep.family == f_][col].to_numpy()
            axs[k].scatter(v, np.full(len(v), i) + np.random.default_rng(SEED_BOOT).uniform(-.18, .18, len(v)), s=10, color=cm[f_], alpha=.8)
        axs[k].set_yticks(range(4)); axs[k].set_yticklabels(FAM_ORDER[:-1]); axs[k].set_title(ttl, fontsize=8)
        axs[k].axvline(1 if col == "elasticity_AT_on_depth" else 0, color="k", lw=.7)
    fig.tight_layout()
    figC = C.fig_to_b64(fig); plt.close(fig)

    # ---------------- numbers for text
    pA = G("AT", "raw_in", "all", "pois"); zA = G("AT", "raw_in", "all", "zinb"); hA = G("AT", "raw_in", "all", "hnb"); fsA = G("AT", "raw_in", "all", "nb_fs")
    zO = G("AT", "out", "all", "zinb"); pO = G("AT", "out", "all", "pois"); zipA = G("AT", "raw_in", "all", "zip")
    pD = G("depth", "out", "all", "pois"); zD = G("depth", "out", "all", "zinb"); hD = G("depth", "out", "all", "hnb")
    oN = offc[(offc.family == "all") & (offc.model == "nb")].iloc[0]
    b1, bb05, mixx, bbphi, bbp = S("all", "b1"), S("all", "bb05"), S("all", "mix"), S("all", "bbphi"), S("all", "bb_pooled")
    jN, jZ, jZB = J("all", "NB_HA+NB_HB"), J("all", "ZINB_HA+ZINB_HB"), J("all", "ZINB_AT+BB")
    a_at = at[(at.modality == "AT") & (at.offset == "raw_in")]
    a_dp = at[(at.modality == "depth")]
    low_ex = zt[(zt.modality == "AT") & (zt.offset == "raw_in") & (zt.coverage_tercile == 0) & (zt.family.isin(["called", "random"]))]
    ex_pts = float((100 * (low_ex.obs_zero_frac - low_ex.pred_zero_nb)).mean())
    ex_z_pts = float((100 * (low_ex.obs_zero_frac - low_ex.pred_zero_zinb)).mean())
    phi_pool = allfam.set_index("AT_stratum").pearson_dispersion
    zero_all_pct = 100 * float(meta["frac_cells_AT_zero_median_region"])
    zex_small = abs(zA.d_logscore_vs_nb) < EQUIV and abs(zO.d_logscore_vs_nb) < EQUIV
    rho_fam = pooled.set_index("family")
    dep_all = dep
    el_med, el_lo, el_hi = dep_all.elasticity_AT_on_depth.median(), dep_all.elasticity_AT_on_depth.min(), dep_all.elasticity_AT_on_depth.max()
    rc_med, rc_lo, rc_hi = dep_all.resid_corr_pearson.median(), dep_all.resid_corr_pearson.min(), dep_all.resid_corr_pearson.max()
    ter_frac = 100 * (dep_all.p_ratio_tercile < 0.01).mean(); clu_frac = 100 * (dep_all.p_ratio_cluster < 0.01).mean()
    off_better = "log depth in the region" if oN.d_logscore_raw_in_minus_out > 0 else "log AT total outside the region"

    verdict = (f"Per cell and region, the allele total AT is adequately described by a negative binomial (NB2) with an exposure offset and cluster terms: Poisson is clearly worse "
               f"({fci(pA.d_logscore_vs_nb, pA.ci_lo, pA.ci_hi)} nats per cell-region vs NB2) and adding a zero-inflation or hurdle component changes the held-out score by "
               f"{zA.d_logscore_vs_nb:+.4f} (ZINB), i.e. by less than the {EQUIV} nats margin used here, although low-coverage cells show a small excess of zeros "
               f"(about {ex_pts:.1f} points above NB2 expectation). The haplotype split HB | AT needs a beta-binomial (or equivalent variance-inflated binomial) conditional on AT, "
               f"not a plain binomial ({fci(b1.d_logscore_per_cell_vs_bb1, b1.ci_lo, b1.ci_hi)} nats per cell-region vs BB), and its cell-level rho is not one number: its region median falls from "
               f"{reg_sp[reg_sp.family == 'saved'].rho_mle.median():.2f} (saved windows, ~1.5 observations per informative cell) to {reg_sp[reg_sp.family == 'literature'].rho_mle.median():.3f} (literature loci, tens of observations), while the variance inflation stays near {phi_pool.mean():.1f}.")
    kf = [
        f"<b>Poisson is rejected for both count modalities.</b> NB2 beats Poisson by {-pA.d_logscore_vs_nb:.3f} nats per cell-region for AT (95% CI {pA.ci_lo:+.3f} to {pA.ci_hi:+.3f} for Poisson minus NB2, "
        f"{meta['n_test_cells']} held-out cells) and by {-pD.d_logscore_vs_nb:.2f} for raw depth ({pD.ci_lo:+.2f} to {pD.ci_hi:+.2f}). ZIP is also worse than NB2 ({zipA.d_logscore_vs_nb:+.4f}), so a zero term on a Poisson base does not help.",
        f"<b>Excess zero mass in AT is small.</b> ZINB minus NB2 = {fci(zA.d_logscore_vs_nb, zA.ci_lo, zA.ci_hi)} (raw-depth offset) and {fci(zO.d_logscore_vs_nb, zO.ci_lo, zO.ci_hi)} (AT-outside offset); "
        f"hurdle-NB {fci(hA.d_logscore_vs_nb, hA.ci_lo, hA.ci_hi)}. ZINB is preferred in {100*zA.frac_regions_better_than_nb:.0f}% of regions but by amounts below {EQUIV}. "
        f"In the lowest coverage tercile of called and random regions the observed zero fraction exceeds the NB2 prediction by {ex_pts:.1f} points on average (ZINB: {ex_z_pts:+.1f}).",
        f"<b>AT is more overdispersed than depth.</b> Median NB2 alpha is {a_at.alpha_nb.median():.2f} for AT (saved windows {a_at[a_at.family=='saved'].alpha_nb.median():.2f}) versus {a_dp.alpha_nb.median():.3f} for depth; "
        f"the dispersion parameter cannot be shared between modalities. Median {zero_all_pct:.0f}% of cells have AT = 0 in a region; the saved chr7/chr9 windows have {100*a_at[a_at.family=='saved'].zero_frac_test.median():.0f}% zeros, which NB2 given depth predicts "
        f"({100*a_at[a_at.family=='saved'].pz_nb.median():.0f}%).",
        f"<b>Offset choice:</b> NB2 with {off_better} as exposure scores better by {fci(oN.d_logscore_raw_in_minus_out, oN.ci_lo, oN.ci_hi)} nats (raw-depth offset minus AT-outside offset), but the raw-depth offset puts the region's own fragments in the exposure and builds in the depth-AT dependence (C).",
        f"<b>HB | AT is overdispersed relative to a binomial but needs no mixture.</b> Binomial(p) minus BB(p,rho): {fci(b1.d_logscore_per_cell_vs_bb1, b1.ci_lo, b1.ci_hi)}; the orientation-unknown two-component mixture is {fci(mixx.d_logscore_per_cell_vs_bb1, mixx.ci_lo, mixx.ci_hi)}; "
        f"a beta-binomial with constant variance inflation phi is {fci(bbphi.d_logscore_per_cell_vs_bb1, bbphi.ci_lo, bbphi.ci_hi)} (all relative to a region-specific BB).",
        f"<b>Cell-level rho is far above the bin-level 0.0038 but the inflation is similar.</b> Pooled family rho: saved {rho_fam.loc['saved','rho_pooled']:.3f}, called {rho_fam.loc['called','rho_pooled']:.3f}, random {rho_fam.loc['random','rho_pooled']:.3f}, literature {rho_fam.loc['literature','rho_pooled']:.4f} "
        f"({rho_fam.loc['literature','ratio_to_exp06']:.1f} to {rho_fam.loc['random','ratio_to_exp06']:.0f} times exp 06). Region rho falls with mean AT per cell (Spearman {spear[0]:.2f} on log-log); "
        f"cell-level Pearson dispersion of HB | AT is {phi_pool.min():.2f} to {phi_pool.max():.2f} across AT strata, versus 1.20 to 1.53 for exp 06 bins. rho is therefore not transferable between units; the variance inflation is.",
        f"<b>HA and HB as separate counts predict almost as well as AT x conditional split.</b> Joint log P(HA,HB) relative to NB(AT) x BB: independent NB {fci(jN.d_joint_logscore_vs_NB_AT_plus_BB, jN.ci_lo, jN.ci_hi)}, independent ZINB {fci(jZ.d_joint_logscore_vs_NB_AT_plus_BB, jZ.ci_lo, jZ.ci_hi)}, ZINB(AT) x BB {fci(jZB.d_joint_logscore_vs_NB_AT_plus_BB, jZB.ci_lo, jZB.ci_hi)}. All within {EQUIV} nats.",
        f"<b>AT and depth are dependent.</b> After cluster and exposure adjustment the Pearson-residual correlation is {rc_med:.2f} (regions {rc_lo:.2f} to {rc_hi:.2f}; permutation p &lt; 0.01 in {100*(dep_all.p_perm_resid_corr<0.01).mean():.0f}% of regions) "
        f"and AT scales sub-proportionally with region depth (median elasticity {el_med:.2f}, range {el_lo:.2f} to {el_hi:.2f}). AT per raw fragment does not vary with coverage tercile (p&lt;0.01 in {ter_frac:.0f}% of regions) but does vary by cluster in {clu_frac:.0f}%.",
    ]

    # recommendation
    rec_rows = [
        {"modality": "(1) read depth", "recommended observation model": "NB2 (var = mu + alpha mu^2), log mu = log(exposure) + cluster terms; exposure = cell raw total outside the region",
         "evidence (held-out log score, nats per cell-region, bootstrap over cells)": f"Poisson minus NB2 {fci(pD.d_logscore_vs_nb, pD.ci_lo, pD.ci_hi, 3)}; ZINB minus NB2 {fci(zD.d_logscore_vs_nb, zD.ci_lo, zD.ci_hi)}; hurdle-NB {fci(hD.d_logscore_vs_nb, hD.ci_lo, hD.ci_hi)}",
         "not needed / rejected": "Poisson; zero-inflation (raw depth has almost no zeros except in the saved windows)",
         "caveat": f"median alpha {a_dp.alpha_nb.median():.3f}; quasi-Poisson (exp 08) is a variance correction without a likelihood and cannot be scored or used in a joint model"},
        {"modality": "(2) allele total AT", "recommended observation model": "NB2 with its own alpha; condition on the region's own depth (offset log(raw depth in region + 0.5)); ZINB only as a sensitivity",
         "evidence (held-out log score, nats per cell-region, bootstrap over cells)": f"Poisson minus NB2 {fci(pA.d_logscore_vs_nb, pA.ci_lo, pA.ci_hi, 3)}; ZINB minus NB2 {fci(zA.d_logscore_vs_nb, zA.ci_lo, zA.ci_hi)} / {fci(zO.d_logscore_vs_nb, zO.ci_lo, zO.ci_hi)} (two offsets); hurdle {fci(hA.d_logscore_vs_nb, hA.ci_lo, hA.ci_hi)}",
         "not needed / rejected": "Poisson and ZIP; hurdle; free offset slope adds nothing (" + fci(fsA.d_logscore_vs_nb, fsA.ci_lo, fsA.ci_hi) + ")",
         "caveat": f"depth-conditioning beats the AT-outside exposure by {fci(oN.d_logscore_raw_in_minus_out, oN.ci_lo, oN.ci_hi)} nats, so AT is not conditionally independent of depth (see C); median alpha {a_at.alpha_nb.median():.2f} = {a_at.alpha_nb.median()/a_dp.alpha_nb.median():.0f} times depth alpha; a small excess of zeros remains in low-coverage cells ({ex_pts:.1f} points); AT is correlated with depth (see C)"},
        {"modality": "(3) haplotype split HB | AT", "recommended observation model": "Conditional beta-binomial given AT; dispersion as a variance inflation phi (rho = (phi-1)/(AT-1)) or a region-specific rho, never a single rho carried across units; no mixture",
         "evidence (held-out log score, nats per cell-region, bootstrap over cells)": f"Binomial(p) minus BB {fci(b1.d_logscore_per_cell_vs_bb1, b1.ci_lo, b1.ci_hi)}; constant-phi BB {fci(bbphi.d_logscore_per_cell_vs_bb1, bbphi.ci_lo, bbphi.ci_hi)}; mixture {fci(mixx.d_logscore_per_cell_vs_bb1, mixx.ci_lo, mixx.ci_hi)}; family-pooled constant-rho BB {fci(bbp.d_logscore_per_cell_vs_bb1, bbp.ci_lo, bbp.ci_hi)}",
         "not needed / rejected": "plain binomial; one fixed rho per family (loses by the pooled-rho row; constant phi is within the margin of region-specific rho); two-component mixture; separate NB for HA and HB is equivalent in prediction (within " + str(EQUIV) + ") but does not isolate the allele-sensitive part",
         "caveat": f"cell-level rho is {rho_fam.loc['literature','rho_pooled']:.3f} to {rho_fam.loc['saved','rho_pooled']:.2f} by family vs the bin-level 0.0038; do not reuse the exp 06 rho at cell level"},
    ]
    tabD = pd.DataFrame(rec_rows)

    # per-region detail (saved + literature)
    det = []
    for _, r in sp[sp.family.isin(["saved", "literature"])].iterrows():
        a_ = at[(at.modality == "AT") & (at.offset == "raw_in") & (at.region == r.region)].iloc[0]
        d_ = dep[dep.region == r.region].iloc[0]
        dz = a_.ls_zinb - a_.ls_nb
        note = []
        note.append(f"AT is zero in {100*a_.zero_frac_test:.0f}% of held-out cells; NB2 given depth predicts {100*a_.pz_nb:.0f}%" +
                    (" (zeros explained by low depth and overdispersion)" if abs(a_.pz_nb - a_.zero_frac_test) < 0.02 else " (a zero excess remains)"))
        note.append("an extra zero term is not worth it" if abs(dz) < EQUIV else ("a zero term helps" if dz > 0 else "a zero term overfits"))
        note.append(f"in cells with &ge;2 observations the split is {'strongly' if r.rho_mle > 0.1 else 'mildly'} overdispersed (rho {r.rho_mle:.3f})")
        note.append(f"AT follows depth with elasticity {d_.elasticity_AT_on_depth:.2f}")
        det.append({"family": r.family, "region": r.region, "cells with AT>=1": int(r.cells_AT_ge1), "mean AT | AT>=1": r.mean_AT_given_ge1,
                    "zero fraction (test)": a_.zero_frac_test, "best AT model": MNAME[a_.best_model], "NB2 alpha": a_.alpha_nb, "ZINB minus NB2 (nats)": dz,
                    "rho (BB MLE)": r.rho_mle, "phi (const. inflation)": r.phi_mle, "resid corr AT-depth": d_.resid_corr_pearson, "note": "; ".join(note)})
    tabDet = pd.DataFrame(det)

    # sensitivity table
    sens = []
    for lab_, g_ in (("Poisson", pA), ("ZINB", zA), ("hurdle-NB", hA), ("NB2 free offset slope", fsA)):
        sens.append({"AT raw_in offset, all regions": lab_, "primary split (train = half True)": fci(g_.d_logscore_vs_nb, g_.ci_lo, g_.ci_hi),
                     "reverse split (train and test swapped), point estimate": g_.d_logscore_reverse_split,
                     "AT-outside offset, primary": fci(G("AT", "out", "all", {"Poisson": "pois", "ZINB": "zinb", "hurdle-NB": "hnb", "NB2 free offset slope": "nb_fs"}[lab_]).d_logscore_vs_nb,
                                                      G("AT", "out", "all", {"Poisson": "pois", "ZINB": "zinb", "hurdle-NB": "hnb", "NB2 free offset slope": "nb_fs"}[lab_]).ci_lo,
                                                      G("AT", "out", "all", {"Poisson": "pois", "ZINB": "zinb", "hurdle-NB": "hnb", "NB2 free offset slope": "nb_fs"}[lab_]).ci_hi)})
    tabS1 = pd.DataFrame(sens)
    sp_rev = spf[spf.family == "all"][["model", "d_logscore_per_cell_vs_bb1", "d_logscore_reverse_split_per_cell"]].copy()
    sp_rev["model"] = sp_rev.model.map(SNAME)
    sp_rev = sp_rev.rename(columns={"d_logscore_per_cell_vs_bb1": "primary (minus BB p,rho)", "d_logscore_reverse_split_per_cell": "reverse split"})
    slope = at[(at.modality == "AT")].groupby("offset").nb_free_offset_slope.median()

    glossary = {
        "modality": "depth = raw fragment count in the region (layers/raw); AT = summed per-SNP-site allele depth over selected phased heterozygous SNPs (not unique fragments); HA/HB = AT split by haplotype, AT == HA + HB.",
        "offset / exposure": "known per-cell multiplier of the mean. raw_in: log(cell raw depth in region + 0.5). out: log(cell AT total outside the region + 0.5) for AT, log(cell raw total outside the region + 0.5) for depth. Coefficient fixed at 1.",
        "family": "saved (7 saved windows), literature (4 loci), called (67 short segments of the `all` karyogram, calls are not truth), random (60 uncalled 25-bin windows, seed tuple (20261104, 15)), all (the 138 regions pooled).",
        "region": "target bins only; AT/HA/HB/depth are summed over the region's retained bins per cell.",
        "Poisson, NB2, ZIP, ZINB, hurdle-NB": "count likelihoods. NB2 variance mu + alpha mu^2. ZIP/ZINB: with probability pi the count is a structural zero, pi = expit(g0 + g1 z), z the standardised offset. Hurdle-NB: P(0) = pi, positives from a zero-truncated NB2.",
        "NB2 free offset slope": "NB2 where the offset enters with a free coefficient instead of 1 (sensitivity).",
        "ls_<model>": "mean held-out log predictive probability per test cell (nats), train = stratified half True, test = other half; ls_rev_ is the reverse direction.",
        "d_logscore_vs_nb / X minus NB2": "mean over regions and held-out cells of log score(model) - log score(NB2); positive = model predicts better than NB2. Brackets: 95% bootstrap over test cells (1000 resamples, seed tuple (20261104, 15, 1)).",
        "frac_regions_better_than_nb": "fraction of regions where the model's mean held-out score exceeds NB2's.",
        "alpha_nb / alpha_zinb": "NB2 dispersion alpha fitted on the training half (var = mu + alpha mu^2).",
        "pi_zinb_mean, pi_zip_mean, pi_hurdle_zero_mean": "mean (over test cells) fitted structural-zero (ZINB/ZIP) or zero-hurdle probability.",
        "nb_free_offset_slope": "fitted coefficient on the offset (in offset units) when the exposure coefficient is free; 1 = proportional to exposure.",
        "zero_frac_test, obs_zero_t*": "observed fraction of held-out cells with count 0 (overall; in coverage tercile 0 = low, 1, 2 = high; tercile of cell raw total).",
        "pz_<model>, pz_<model>_t*": "predicted P(Y = 0) averaged over held-out cells (overall, per tercile).",
        "zz_<model>, zz_<model>_t*": "z = (observed zero count - sum of predicted P0) / sqrt(sum P0 (1 - P0)); |z| > 2 flags a zero count the model does not reproduce.",
        "best_model": "model with highest mean held-out score among Poisson, NB2, ZIP, ZINB, hurdle-NB for that region; differences are often below the practical margin.",
        "practical margin": f"{EQUIV} nats per cell-region, a judgement for this review, not a validated threshold.",
        "HB | AT models": "Binomial p=0.5; Binomial p (region p fitted on training cells); BB p=0.5 (rho only); BB(p,rho); BB constant inflation phi (rho = (phi-1)/(AT-1) per cell, phi fitted); mixture w Bin(p) + (1-w) Bin(1-p), p >= 0.5 (orientation unknown); bb_pooled = BB with one rho per family and region p.",
        "rho, rho_mle": "beta-binomial intra-class correlation (variance inflation of HB given AT = 1 + (AT - 1) rho); rho_mle is the maximum-likelihood value per region over cells with AT >= 1 (all cells).",
        "phi_mle / pearson_dispersion": "constant variance inflation phi (1 = binomial); pearson_dispersion = mean of (HB - AT p)^2 / (AT p (1-p)) over cells (or cell-region pairs) with AT in the stratum.",
        "rho_pooled / rho_ci": "family-level rho shared by regions (per-region p fixed at its binomial MLE), grid maximum likelihood; CI = bootstrap over cells (400 resamples) on the grid, so it is as coarse as the grid.",
        "ratio_to_exp06": f"rho_pooled / {rho06:.4f}, the bin-level beta-binomial rho from experiment 06 (bins pooled over cells).",
        "lrt_vs_binomial": "2 x log-likelihood gain of BB over binomial, shared p per region; p_lrt uses 0.5 chi-square(1) (boundary).",
        "d_logscore_per_cell_vs_bb1": "mean over regions and held-out cells of log score(model) - log score(BB(p,rho)); cells with AT = 0 contribute 0 to every split model. d_logscore_per_informative_pair divides by informative pairs only.",
        "joint models": "log P(HA, HB) on held-out cells: independent NB (or ZINB) for HA and HB with offset log(raw in region + 0.5), versus NB (or ZINB) for AT x conditional split.",
        "resid_corr_pearson / spearman": "correlation across cells of Pearson residuals of AT (NB2, offset log AT outside) and raw depth (NB2, offset log raw outside), both with cluster terms. p_perm: two-sided permutation p, residuals of depth shuffled within cluster x coverage tercile (300 permutations; minimum possible p = 1/301).",
        "elasticity_AT_on_depth": "coefficient b in AT ~ exposure x exp(cluster) x depth_in_region^b (Poisson GLM, robust HC0 SE in elasticity_se_robust). b = 1: AT proportional to depth in the region; b = 0: unrelated beyond exposure and cluster.",
        "AT_per_raw, AT_per_raw_t*": "sum AT / sum raw depth over all cells (and in coverage tercile t). log_ratio_low_over_high: log of tercile 0 over tercile 2. p_ratio_tercile / p_ratio_cluster: quasi-Poisson Wald F test of AT/depth heterogeneity (offset log(raw depth + 0.5)).",
        "cells_AT_ge1, cells_AT_ge2, mean_AT_given_ge1, max_AT": "cells (of 3,460) with AT at least 1 / 2; mean AT among cells with AT >= 1; largest per-cell AT in the region.",
        "p_hat_HB_over_AT, p_bb, mix_p, mix_w": "region HB/AT, BB mean, mixture p (>= 0.5) and weight of the Bin(p) component. These describe fit, not copy number.",
    }

    caveat_pairs = [
        "These are observation-model fits per region; they say nothing about copy number. A model that predicts AT better is not evidence about dosage.",
        "Held-out scores differ by a few thousandths of a nat for zero-inflation. The margin of 0.002 is a judgement; with 1,730 test cells a difference of 0.0004 is statistically resolvable but not practically important.",
        "Bootstrap over cells treats regions as fixed and cells as exchangeable within the half; cells in the same cluster share structure, and the two halves share the same donors/library, so intervals are optimistic for transfer to a new sample.",
        "AT counts selected-site depth, not independent fragments. Reads from one haplotype molecule at several SNPs, or duplicate fragments, create within-cell correlation (high rho at small AT) that has nothing to do with biology.",
        "The raw-depth offset includes the region's own fragments. That is useful for prediction but makes AT and depth mechanically dependent; the AT-outside offset avoids it but ignores region-specific accessibility.",
        "Region choice: saved and literature families have 7 and 4 regions only; called segments were selected because the pipeline called them; random windows are restricted to bins with cn_call == 1 in the `all` karyogram and are not neutral truth.",
        "ZINB and hurdle fits use a constant-plus-offset zero model. A zero process tied to the number of heterozygous sites in reach (not observable per cell here) could be misrepresented, so 'no excess zeros' means no excess that this simple zero model can find.",
        "Phase labels HA/HB come from the pipeline's SNP phasing; mislabelled phase would inflate rho. The orientation-unknown mixture did not outperform, but phasing quality is not assessed here.",
        "Only one dataset (MCF10-A development data) and one 100 kb bin grid; regions are 25 bins or more.",
    ]

    sections = [
        ("Question and why it matters", f"""<div class='card'><p>A combined depth + allele likelihood multiplies a depth term and an allele term. If the allele term uses a wrong
count model (for example Poisson for AT, or a binomial split without overdispersion), the allele evidence is over-weighted and overdispersed noise is read as signal. Depth is already known to be
NB2 with exposure (experiment 08). This experiment asks the same question of the allelic modality at the <b>cell</b> level, region by region: how should AT be modelled (many cells have AT = 0), how
should the haplotype split HB | AT be modelled, can HA and HB simply be modelled as two counts, and how dependent are AT and depth? The aim is to pick an observation model per modality by held-out evidence, not by convention.</p></div>"""),
        ("Inputs", f"""<ul>
<li>AnnData <code>count_matrix.h5ad</code> (layers/raw; {meta['n_cells']} cells) and WNN cluster labels ({len(meta['clusters'])} clusters).</li>
<li>HA, HB, AT parquet (bins x barcodes, joined by barcode name; coordinates equal to AnnData; AT == HA + HB in all regions: {meta['alignment_checks']['AT_equals_HA_plus_HB_in_regions']}; integer valued). Column order differs from AnnData, set equality holds; matrices were re-ordered by name.</li>
<li>Regions ({len(meta['regions'])}): 7 saved windows, 4 literature loci, 67 called short segments (experiment 04 rule), and 60 random 25-bin windows drawn without replacement from {meta['n_eligible_uncalled_random_windows']:,} eligible windows (no gap, every bin cn_call == 1 in the `all` karyogram), seed tuple (20261104, 15).</li>
<li>Train/test: <code>stratified_halves(labels, (20261104, 15))</code>, {meta['n_train_cells']} train, {meta['n_test_cells']} test cells; the reverse direction is a sensitivity. Experiment 06 pooled bin-level rho = {rho06:.4f} and experiment 08 depth results are read for comparison.</li></ul>"""),
        ("Method", """<pre>
A. count y_ir (cell i, region r), cluster dummies g_c
   log mu_ir = offset_ir + b0 + g_c(i)                       offset coefficient fixed at 1
   Poisson : P(y) = mu^y e^-mu / y!
   NB2     : var = mu + alpha mu^2                           alpha fitted by ML
   ZIP/ZINB: P(0) = pi + (1-pi) f(0), P(y>0) = (1-pi) f(y)   pi = expit(g0 + g1 z), z = standardised offset
   hurdle  : P(0) = pi,  P(y>0) = (1-pi) f_NB(y) / (1 - f_NB(0))
   offsets : raw_in = log(depth_in_region + 0.5)             [AT only]
             out    = log(AT outside region + 0.5)  (AT);  log(raw outside region + 0.5)  (depth)
   score   : mean over held-out cells of log P(y_ir | train fit); difference to NB2 averaged over regions of a family,
             bootstrap over test cells (1000 resamples, seed (20261104,15,1)).
B. HB | AT = n, cells with n >= 1
   Bin(0.5); Bin(p); BB(0.5,rho); BB(p,rho); BB(p, rho_n = (phi-1)/(n-1)) (constant variance inflation);
   w Bin(p) + (1-w) Bin(1-p) with p >= 0.5.   BB variance = n p (1-p) [1 + (n-1) rho]
   pooled rho per family: grid MLE of one rho shared by regions (p per region fixed), bootstrap over cells (400).
   joint log P(HA,HB): NB(HA) NB(HB) | ZINB(HA) ZINB(HB) | NB(AT) Bin(HB|AT) | NB(AT) BB | ZINB(AT) BB | ZINB(AT) mixture   (offset log(raw_in + 0.5))
C. e_AT, e_depth = Pearson residuals of NB2 fits (cluster terms, outside exposures); corr(e_AT, e_depth), permutation within cluster x tercile;
   log E[AT] = log AT_out + g_c + b (log(depth_in + 0.5) - mean);  AT/depth by tercile and cluster (quasi-Poisson Wald F).
All optimisation: L-BFGS-B, numerical gradients; Poisson IRLS start. No scientific code from the pipeline is changed.
</pre>"""),
        ("Results", "<p>Each sub-section has its table, figure and caption. Positive differences favour the candidate model.</p>"),
        ("A. Allele total AT as a count (and depth for context)",
         "<h4>Held-out log score relative to NB2</h4>" + C.html_table(tabA1) +
         C.caption("Mean held-out log predictive probability of each model minus that of NB2, nats per cell-region, over regions of the family; brackets = 95% bootstrap over held-out cells. Negative = worse than NB2. The last column is the share of regions where ZINB has the higher mean score.") +
         C.img(figA1, "Held-out score differences of ZINB, hurdle and free-slope NB relative to NB2") +
         C.caption(f"Zoomed view of the three alternatives to NB2 for each model class and family (points: mean difference; bars: 95% bootstrap over cells; grey band: the +/-{EQUIV} nats practical margin). Poisson and ZIP are off scale and are in the table.") +
         "<h4>Zero calibration</h4>" + C.html_table(tabA2) +
         C.caption("Held-out zero count against model prediction, per region, summarised by family. z compares the observed number of zero cells with the sum of predicted zero probabilities; |z| > 2 means the model does not reproduce the zeros. alpha = NB2 dispersion; pi = fitted structural-zero probability in the ZINB.") +
         C.img(figA2, "Predicted minus observed zero fraction by coverage tercile") +
         C.caption("For AT with the raw-depth offset: predicted minus observed zero fraction (percentage points) by cell coverage tercile, averaged over regions of each family. Zero = model reproduces the zeros; below zero = the model predicts too few zeros.") +
         C.html_table(tabA3, floatfmt="{:.2f}") +
         C.caption("The same quantities as numbers (AT only, both offsets). Predicted minus observed zero fraction in percentage points.") +
         "<h4>Exposure choice for AT</h4>" + C.html_table(off_t) +
         C.caption("Held-out score of the raw-depth offset minus that of the AT-outside offset (same model, same split); positive favours conditioning on the region's own depth. Brackets: bootstrap over cells."), 3),
        ("B. Haplotype split HB | AT",
         C.html_table(tabB1) +
         C.caption("Held-out log score of each split model minus the region-specific beta-binomial BB(p,rho), nats per cell-region (cells with AT = 0 contribute 0 to all). Negative = worse than BB(p,rho). Brackets: bootstrap over test cells.") +
         "<h4>Cell-level overdispersion versus the bin-level estimate</h4>" + C.html_table(tabB2) +
         C.caption(f"Per-region beta-binomial fits (all cells with AT >= 1) summarised by family. rho is the intra-class correlation; phi the variance inflation 1 + (n - 1) rho; Pearson dispersion uses cells with AT >= 2. Spearman correlation of region rho with mean AT per cell: {spear[0]:.2f} (log-log).") +
         C.html_table(tabB3) + C.caption(f"Family-pooled rho with its bootstrap interval (grid resolution is coarse) against the experiment 06 bin-level value {rho06:.4f}.") +
         C.img(figB, "Cell-level rho and variance inflation versus AT") +
         C.caption("Left: per-region cell-level beta-binomial rho against mean AT per cell; dashed line is rho implied by a constant variance inflation of 1.4; red line is the experiment 06 bin-level rho. Right: variance inflation of HB | AT in cell-region pairs by AT stratum (circles, per family) and in experiment 06 bins pooled over cells (red squares).") +
         C.html_table(dsp, floatfmt="{:.3f}") + C.caption("Pearson dispersion of HB | AT by AT stratum in cell-region pairs (pooled row: pair-weighted mean over families) and the rho it implies, (phi - 1) / (mean AT - 1).") +
         C.html_table(ov, floatfmt="{:.4f}") + C.caption("Experiment 06 bin-level estimate for reference (bins in base-called `all` segments with AT >= 20, read from its overdispersion_estimate.tsv).") +
         "<h4>HA and HB as separate counts</h4>" + C.html_table(tabB5) +
         C.caption("Joint held-out log P(HA, HB) of each candidate minus NB(AT) x BB(HB | AT), nats per cell-region, brackets 95% bootstrap over cells. NB_HA+NB_HB treats the two haplotype counts as independent NB with the raw-depth offset."), 3),
        ("C. Dependence between AT and depth",
         C.html_table(tabC, floatfmt="{:.3g}") +
         C.caption("Per family: correlation of AT and depth Pearson residuals (after exposure and cluster adjustment), permutation test, elasticity of AT on region depth, and heterogeneity of AT per raw fragment across coverage tercile and cluster (quasi-Poisson Wald F).") +
         C.img(figC, "AT-depth dependence by region") +
         C.caption("One point per region. Left: residual correlation (0 = independent after adjustment). Middle: elasticity of AT on region depth (1 = proportional, dashed line at 1). Right: log ratio of AT per fragment in the lowest to highest coverage tercile (0 = no coverage dependence)."), 3),
        ("D. Recommendation per modality",
         C.html_table(tabD) +
         C.caption("Recommended observation model per modality with the held-out evidence (nats per cell-region relative to the stated reference; brackets are 95% bootstrap over cells). Choices rest on predictive adequacy of the observation model, not on copy-number truth."), 3),
        ("Per-region detail", C.html_table(tabDet, floatfmt="{:.3g}") +
         C.caption("One row per saved window and literature locus (AT with the raw-depth offset). The note is generated from the numbers in the row, not an interpretation of biology.")),
        ("Sensitivity", f"<ul><li>Reverse train/test direction: the point estimates agree in sign and size with the primary split (table).</li>"
         f"<li>Two exposure definitions for AT give the same ranking of models; the free offset slope fits at median {slope.get('raw_in', float('nan')):.2f} (raw-depth offset) and {slope.get('out', float('nan')):.2f} (AT-outside offset) per unit of offset (1 = proportional exposure, as assumed); freeing it changes the held-out score by less than {EQUIV}.</li>"
         f"<li>Zero models use a constant-plus-offset logit; this was not varied.</li></ul>" + C.html_table(tabS1, floatfmt="{:+.4f}") +
         C.caption("AT, all 138 regions: alternatives minus NB2. Left columns raw-depth offset; the reverse split column gives the point estimate with train and test swapped.") +
         C.html_table(sp_rev, floatfmt="{:+.4f}") + C.caption("Split models minus BB(p,rho), all regions: primary direction and reverse direction (point estimates).")),
        ("Limits and alternative explanations", "<ul>" + "".join(f"<li>{x}</li>" for x in caveat_pairs) + "</ul>"),
    ]
    supports = [
        "NB2 with exposure is an adequate count likelihood for AT; Poisson, ZIP are not.",
        "Zeros in AT are mostly explained by depth and overdispersion; the leftover excess is a small low-coverage effect.",
        "A conditional beta-binomial with variance-inflation dispersion is adequate for HB | AT; a plain binomial and a fixed rho carried over from the bin level are not.",
        "AT and depth are positively dependent per cell beyond exposure, so a combined model cannot treat them as independent given copy number.",
    ]
    not_supports = [
        "Anything about copy number or allelic imbalance at any locus; these are observation-model comparisons.",
        "That no structural-zero process exists: a zero model linked to heterozygous-site density per cell could not be fitted from these data.",
        "That rho = 0.0038 is wrong at the bin level; only that it does not transfer to cell level (the variance inflation does).",
        "Independence of the halves: both are from the same cells' library, so generalisation to another sample is untested.",
    ]
    decision_impact = ("Supports recommendation.md item 4 (a joint depth-allele likelihood evaluated later) only with specific components: NB2 for depth, NB2 with its own alpha for AT, and a beta-binomial split with "
                       "constant-phi dispersion, plus an explicit AT-depth dependence term (shared cell effect or AT conditioned on depth) because residual correlation is "
                       f"{rc_med:.2f}. It argues against reusing the bin-level rho = {rho06:.4f} at cell level and against Poisson/binomial shortcuts. It does not change the rank-1 evidence layer.")
    next_steps = [
        "Score the constant-phi beta-binomial against a phi that depends on SNP count per bin (needs the per-site table, not available here).",
        "Fit a joint depth-AT-split likelihood with a shared cell random effect and check whether the residual AT-depth correlation falls to zero.",
        "Repeat on an independent sample and with a held-out donor/library to test transfer of alpha and phi.",
        "Obtain heterozygous-site counts per cell and region to test whether the small low-coverage zero excess is 'no informative SNP in reach'.",
        "Add a phasing-quality check (random re-labelling of phase) to separate technical from biological split overdispersion.",
    ]
    prov = {"command": "/home/user/envs/cnv310/bin/python -I research/2026-11-04-cloud-method-review/experiments/15_allelic_count_models/run.py",
            "seeds": {"regions_random": list(SEED), "train_test_halves": list(SEED), "bootstrap_and_permutation": list(SEED_BOOT)},
            "parameters": {"n_random_windows": N_RANDOM, "bootstrap_reps_models": meta["bootstrap_reps_models"], "bootstrap_reps_rho": meta["bootstrap_reps_rho"],
                           "permutations": meta["permutations"], "practical_margin_nats": EQUIV, "rho_grid_points": len(RHO_GRID)},
            "inputs_sha256": {"count_matrix.h5ad": C.sha256(paths["h5ad"]), **{f"{n_}.parquet": C.sha256(paths["allelic"] / f"{n_}.parquet") for n_ in ("AT", "HA", "HB")}},
            "compute_seconds": meta["compute_seconds"], "check_mode": meta["check_mode"], "environment": C.environment()}
    C.write_html_report(out / "report.html", "Experiment 15: observation models for allelic counts per cell",
                        "Which likelihood fits depth, the allele total AT and the haplotype split at the cell level? Exploratory diagnostic on MCF10-A development data; model comparison, not copy-number truth.",
                        sections, prov, key_findings=kf, verdict=verdict, glossary=glossary, supports=supports, not_supports=not_supports,
                        decision_impact=decision_impact, next_steps=next_steps)
    summ = {"experiment": "15_allelic_count_models", "seeds": prov["seeds"], "parameters": prov["parameters"], "environment": prov["environment"],
            "alignment_checks": meta["alignment_checks"], "regions": len(meta["regions"]), "n_train_cells": meta["n_train_cells"], "n_test_cells": meta["n_test_cells"],
            "A_AT_rawin_all_regions_vs_NB2": {m: dict(G("AT", "raw_in", "all", m)[["d_logscore_vs_nb", "ci_lo", "ci_hi", "frac_regions_better_than_nb"]]) for m in ["pois", "zip", "zinb", "hnb", "nb_fs"]},
            "A_AT_out_all_regions_vs_NB2": {m: dict(G("AT", "out", "all", m)[["d_logscore_vs_nb", "ci_lo", "ci_hi"]]) for m in ["pois", "zip", "zinb", "hnb", "nb_fs"]},
            "A_depth_all_regions_vs_NB2": {m: dict(G("depth", "out", "all", m)[["d_logscore_vs_nb", "ci_lo", "ci_hi"]]) for m in ["pois", "zip", "zinb", "hnb", "nb_fs"]},
            "A_offset_rawin_minus_out_NB2": dict(oN[["d_logscore_raw_in_minus_out", "ci_lo", "ci_hi"]]),
            "A_low_tercile_zero_excess_vs_NB2_points": ex_pts, "A_low_tercile_zero_excess_vs_ZINB_points": ex_z_pts,
            "A_median_alpha": {"AT_raw_in": float(a_at.alpha_nb.median()), "depth": float(a_dp.alpha_nb.median())},
            "B_all_vs_BB": {m: dict(S("all", m)[["d_logscore_per_cell_vs_bb1", "ci_lo", "ci_hi"]]) for m in ["b05", "b1", "bb05", "bbphi", "mix", "bb_pooled"]},
            "B_pooled_rho": {r.family: {"rho": r.rho_pooled, "ci": [r.rho_ci_lo, r.rho_ci_hi], "ratio_to_exp06": r.ratio_to_exp06} for r in pooled.itertuples()},
            "B_rho_vs_meanAT_spearman_loglog": float(spear[0]), "B_pearson_dispersion_pooled_by_AT_stratum": {k: float(v) for k, v in phi_pool.items()},
            "B_joint_all": {r.joint_model: {"d": r.d_joint_logscore_vs_NB_AT_plus_BB, "ci": [r.ci_lo, r.ci_hi]} for r in joint[joint.family == "all"].itertuples()},
            "C_residual_corr_median_min_max": [rc_med, rc_lo, rc_hi], "C_elasticity_median_min_max": [el_med, el_lo, el_hi],
            "C_pct_regions_ratio_varies_by_tercile_p01": ter_frac, "C_pct_regions_ratio_varies_by_cluster_p01": clu_frac,
            "recommendation": {"depth": "NB2 with exposure", "AT": "NB2 conditioned on region depth (offset log raw depth in region + 0.5) with its own alpha (ZINB only as sensitivity)", "split": "conditional beta-binomial with constant variance inflation phi"},
            "practical_margin_nats": EQUIV, "compute_seconds": meta["compute_seconds"], "check_mode": meta["check_mode"]}
    C.write_json(out / "summary.json", summ)


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
