"""Experiment 18: cell-level joint depth + allele mixture for allele-specific relative states (exploratory diagnostic).

QUESTION. What fraction pi of cells carries an allele-specific relative state e in a region (overall and per WNN cluster),
which haplotype is major, and is the evidence clonal, subclonal, absent or insufficient?  Observation models follow exp 15:
depth ~ NB2 with exposure = the cell's raw total OUTSIDE the region; allele ~ beta-binomial conditional on AT (N) with a
region-specific, N-dependent rho.  Orientation of alleles to major/minor follows exp 16 (`all` group).  The state grid is exp 17's.
Nothing here is copy-number truth; "event" = a relative state on a 1:1 background (MCF10-A baseline assumed, unverified).

PRESPECIFIED MODEL (per region r, cell i; y = raw target count, N = AT on oriented bins, M = counts on the inferred major haplotype)
  cell state s_i in {0 = background 1:1, 1 = event e}, P(s_i = 1) = pi  (overall fit: one pi; cluster fits: pi_c).
  depth     y_i | s ~ NB2(mu_i * L_s, alpha),  var = m + alpha m^2,  alpha = 0.02 (exp 15 median 0.0199), L_0 = 1
            mu_i = E_i * (n_t / n_out) * exp(g_c) * exp(a_c),     E_i = raw total of cell i outside the region,
            n_t / n_out = target bins / bins outside the region.
            g_c  = cluster-c FLANK term (regions <= 100 bins, or the frozen saved flanks):  log[ sum_i in c F_i / (n_f/n_out * sum_i in c E_i) ]
                   i.e. mu uses the cluster's local flank rate;  g_c = 0 (genome-relative) for regions > 100 bins without frozen flanks.
            a_c = region baseline log ratio target/(flank or genome) in cluster c: a nuisance, a_c = b + u_c with
                   b ~ N(m, sigma_b^2) SHARED by all cells (the region's baseline offset, estimated jointly with pi),
                   u_c ~ N(0, sigma_u^2) per cluster; m, sigma_b, sigma_u are robust (median / 1.4826 MAD) summaries over
                   UNCALLED random 25-bin windows (calibration set, disjoint from the 30-window FPR panel).
            Because b is estimated jointly and penalised by the random-window spread, a homogeneous accessibility shift can be
            absorbed by b instead of masquerading as pi > 0; a shift larger than ~sigma_b still looks like an event.
  allele    M_i | N_i, s ~ BetaBinomial(N_i, h_s, rho_i),   rho_i = rho_r / N_i  (so overdispersion factor 1 + rho_r (N-1)/N saturates,
            consistent with exp 15's rho falling with AT and Pearson dispersion 1.27 -> 1.58 with AT; bin-level rho 0.0038 NOT reused).
            rho_r: maximum likelihood per region at the pooled h (grid), held fixed afterwards.
            h_0 = 0.5 + d0 ;  h_1 = clip(h_0 + k (h_e - 0.5)),  h_e = major-haplotype fraction of state e (exp 17 grid), k = attenuation.
            d0 = background major-fraction offset (orientation was chosen from the same counts, so h_0 > 0.5 even when balanced):
            nuisance, d0 ~ N(d_med, tau^2) with d_med, tau from the same uncalled random windows run through the same orientation rule.
            Cells with N_i = 0 contribute no allele term.  AT total itself is NOT modelled: conditioning on N_i removes the
            AT-depth dependence from the allele term; the residual dependence (depth and allele share fragments) is quantified by
            dropping the allele term (sensitivity).
  events    state e : (L_e, h_e)  : loss1 (0.5, 1.0)  homdel (0.02, 0.5)  cnloh (1.0, 1.0)  gain1 (1.5, 2/3)  gain2 (2.0, 0.5)  amp (2.5, 0.8)
            (L_homdel = 0.02, not 0, to allow background leakage; each type fitted separately; best type = highest penalised log-likelihood).
  orientation  bin oriented by exp 16 `all` posterior (B-major if p_B_major >= 0.5); only bins with certainty max(p, 1-p) >= 0.8 are used;
            M_i = sum of HB on B-major bins + HA on A-major bins.  Own-cluster orientation (exp 16 cluster groups) is only a descriptive column.
INFERENCE.  Penalised log-likelihood  l = sum_i log[(1-pi) f_0i + pi f_1i] - 1/2 (a-m)' S^-1 (a-m) - 1/2 ((d0-d_med)/tau)^2, maximised by L-BFGS-B with
  analytic gradients (a, d0, pi free; pi in [0,1]).  Per-cell P(event) = pi f_1 / mix.  95% CI for pi: profile likelihood (grid 0.02, 2*dl <= 3.84).
  LRT_type = 2 (l_alt - l_null); T = max over the 6 types (selection included in the null).  p = (1 + #{T* >= T}) / (B + 1) from a PARAMETRIC
  BOOTSTRAP: y*, M* simulated from the fitted null (MAP a, d0; N, E fixed), refit with the identical pipeline; B = 99, PCG64 seed
  (20261104, 18, region_index) (random windows: region_index = 1000 + j).  Cluster fits (K = 1) use prior mean = mean of the overall null a-hat,
  sd sigma_u, and their own d_med / tau per cluster; their bootstrap p uses the same replicates (max over types per cluster).
CLASSIFICATION (prespecified).  p < 0.05: 'clonal' if profile CI lower bound >= 0.8; 'subclonal' if CI upper bound < 0.95 and lower > 0;
  otherwise 'event, clonal vs subclonal unresolved'.  p >= 0.05: 'no evidence' if CI upper <= 0.25, else 'insufficient evidence'.
  Clonal vs subclonal depends on the nuisance baseline: a clonal event whose depth change is within sigma_b is carried by the allele term only.
FPR.  30 uncalled random 25-bin windows (rule of exp 04 'random': consecutive retained bins, no gap; no called bin in the `all` karyogram),
  drawn with seed (20261104, 18), pushed through the same pipeline (own calibration windows excluded); empirical rate of p < 0.05 and of
  'subclonal event detected' / 'clonal'.
SENSITIVITY.  k in {0.6 primary, 1.0}; with vs without allele; depth removed (allele only); with vs without the shared baseline offset
  (offset removed = b fixed at m, so a baseline shift cannot be absorbed).
Run:  /home/user/envs/cnv310/bin/python -I research/2026-11-04-cloud-method-review/experiments/18_cell_level_joint_states/run.py [--check]
"""
import importlib.util
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
from scipy import optimize, sparse  # noqa: E402
from scipy.special import digamma, gammaln  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

spec = importlib.util.spec_from_file_location("exp04", HERE.parent / "04_fixed_region_contrasts" / "run.py")
E4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(E4)

CHECK = "--check" in sys.argv
SCR = Path("/tmp/claude-0/-home-user-scientific-explorer/b1079aa7-80fd-5756-88cc-d9253976037e/scratchpad/exp18_check")
OUT = SCR if CHECK else HERE
OUT.mkdir(parents=True, exist_ok=True)

ALPHA, EPS_OUT, PI_MIN = 0.02, 0.01, 0.02
B_MAIN = 9 if CHECK else 99
B_VAR = 9 if CHECK else 49
K_PRIMARY, K_ALT = 0.6, 1.0
CERT = 0.8
N_PANEL, N_CALIB = (6, 30) if CHECK else (30, 120)
TYPES = {"loss1": (0.5, 1.0), "homdel": (0.02, 0.5), "cnloh": (1.0, 1.0), "gain1": (1.5, 2 / 3), "gain2": (2.0, 0.5), "amp": (2.5, 0.8)}
PI_GRID = np.round(np.arange(0.0, 1.0001, 0.02), 2)
SEED18 = (C.SEED, 18)
CLUST = ["0", "1", "2", "3", "4"]
mad = lambda v: float(1.4826 * np.median(np.abs(np.asarray(v) - np.median(v))))  # noqa: E731


# ============================================================================= data preparation
def make_regions(d, p):
    regs = [r for r in E4.build_regions(d, p) if r["family"] in ("saved", "literature")]
    for r in regs:
        r["mode"] = "flank" if len(r["flank"]) else "genome"
    # literature loci short enough to have a local flank get the exp-04 flank rule (50 bins each side)
    for r in regs:
        if r["family"] == "literature" and len(r["target"]) <= 100:
            ci = np.flatnonzero(d.seq == d.seq[r["target"][0]])
            lo = int(np.searchsorted(ci, r["target"][0])); hi = int(np.searchsorted(ci, r["target"][-1])) + 1
            l, rr = E4.flanks(ci, lo, hi)
            r["flank"] = np.r_[l, rr]; r["mode"] = "flank"
    t = np.flatnonzero((d.seq == "chr8") & (d.start >= 115e6) & (d.end <= 145.2e6))
    regs.append({"family": "hypothesis", "name": "chr8_115_145_cluster4_hypothesis", "direction": "gain", "target": t, "flank": np.array([], int), "mode": "genome"})
    return regs


def make_windows(d, call, taken, n_total, rng_seed):
    """Random 25-bin windows (exp-04 rule), uncalled, not touching analysed regions; seed (20261104, 18)."""
    rng = np.random.default_rng(rng_seed)
    cand = []
    for ch in dict.fromkeys(d.seq):
        ci = np.flatnonzero(d.seq == ch)
        for k in range(E4.FLANK, len(ci) - 25 - E4.FLANK):
            w = ci[k:k + 25]
            if d.start[w[-1]] - d.start[w[0]] == 24 * 100_000 and (call[w] == 1).all() and not taken[w].any():
                cand.append((ch, k))
    order = rng.permutation(len(cand))
    wins = []
    used = np.zeros(d.n_bins, bool)
    for j in order:
        ch, k = cand[j]
        ci = np.flatnonzero(d.seq == ch)
        w = ci[k:k + 25]
        if used[w].any():
            continue
        l, r = E4.flanks(ci, k, k + 25)
        wins.append({"family": "random", "name": f"random_{ch}:{d.start[w[0]]/1e6:.1f}", "direction": "random", "target": w, "flank": np.r_[l, r], "mode": "flank"})
        used[w] = True
        if len(wins) == n_total:
            break
    return wins, len(cand)


def orientation_tables(d):
    bo = pd.read_csv(HERE.parent / "16_major_minor_haplotypes" / "bin_orientation.tsv.gz", sep="\t")
    tabs = {}
    for g, s in bo.groupby("group"):
        usedB = np.zeros(d.n_bins); usedA = np.zeros(d.n_bins)
        b = s["bin"].to_numpy(); pB = s["p_B_major"].to_numpy(); cert = np.maximum(pB, 1 - pB)
        ok = cert >= CERT
        usedB[b[ok & (pB >= 0.5)]] = 1.0
        usedA[b[ok & (pB < 0.5)]] = 1.0
        tabs[g] = (usedB, usedA)
    return tabs


def indicator(bins_list, n_bins, weight=None):
    rows, cols, vals = [], [], []
    for j, b in enumerate(bins_list):
        b = np.asarray(b, int)
        rows += list(b); cols += [j] * len(b)
        vals += list(np.ones(len(b)) if weight is None else weight[b])
    return sparse.csr_matrix((vals, (rows, cols)), shape=(n_bins, len(bins_list)))


def prepare(d, p, lab, units, tabs, own_units):
    """Per-cell sums for every unit: Y (target raw), F (flank raw), M, N (oriented, `all`), plus own-orientation M, N for `own_units`."""
    U = len(units)
    T = indicator([u["target"] for u in units], d.n_bins)
    F = indicator([u["flank"] for u in units], d.n_bins)
    uB, uA = tabs["all"]
    WB = indicator([u["target"] for u in units], d.n_bins, uB)
    WA = indicator([u["target"] for u in units], d.n_bins, uA)
    WN = indicator([u["target"] for u in units], d.n_bins, uB + uA)
    own = {}
    for c in CLUST:
        oB, oA = tabs["cluster" + c]
        own[c] = (indicator([units[j]["target"] for j in own_units], d.n_bins, oB), indicator([units[j]["target"] for j in own_units], d.n_bins, oA),
                  indicator([units[j]["target"] for j in own_units], d.n_bins, oB + oA))
    cells = list(d.cells)
    n = d.n_cells
    Y = np.zeros((n, U)); Fl = np.zeros((n, U)); tot = np.zeros(n); M = np.zeros((n, U)); N = np.zeros((n, U))
    Mo = np.zeros((n, len(own_units))); No = np.zeros((n, len(own_units)))
    pfs = {nm: pq.ParquetFile(p["allelic"] / f"{nm}.parquet") for nm in ("HA", "HB", "AT")}
    names = pfs["AT"].schema_arrow.names
    if not all(set(pfs[k].schema_arrow.names) == set(names) for k in pfs) or set(names) - {"chr", "start", "end"} != set(cells):
        raise SystemExit("barcode alignment failure")
    ix = pq.read_table(p["allelic"] / "AT.parquet", columns=["chr", "start", "end"]).to_pandas().index
    assert (ix.get_level_values("chr").astype(str).to_numpy() == d.seq).all() and (ix.get_level_values("start").to_numpy() == d.start).all() \
        and (ix.get_level_values("end").to_numpy() == d.end).all(), "coord mismatch"
    step = 250
    for lo in range(0, n, step):
        hi = min(n, lo + step)
        blk = d.rows("layers/raw", lo, hi)
        Y[lo:hi] = T.T.dot(blk.T).T; Fl[lo:hi] = F.T.dot(blk.T).T; tot[lo:hi] = blk.sum(1)
        cols = cells[lo:hi]
        HA = pfs["HA"].read(columns=cols).to_pandas()[cols].to_numpy()      # bins x cells
        HB = pfs["HB"].read(columns=cols).to_pandas()[cols].to_numpy()
        AT = pfs["AT"].read(columns=cols).to_pandas()[cols].to_numpy()
        assert np.allclose(HA + HB, AT)
        M[lo:hi] = (WB.T.dot(HB) + WA.T.dot(HA)).T; N[lo:hi] = WN.T.dot(AT).T
        lc = lab[lo:hi]
        for c in CLUST:
            sel = lc == c
            if sel.any():
                wb, wa, wn = own[c]
                Mo[lo:hi][sel] = (wb.T.dot(HB[:, sel]) + wa.T.dot(HA[:, sel])).T; No[lo:hi][sel] = wn.T.dot(AT[:, sel]).T
        del HA, HB, AT
    return dict(Y=Y, F=Fl, tot=tot, M=M, N=N, Mo=Mo, No=No)


# ============================================================================= likelihood machinery
class Dat:
    """Cell-level arrays for one fit. y, mu0 length n; cl in 0..K-1; allele cells = N > 0."""

    def __init__(self, y, mu0, cl, K, N, M, rho_r, alpha=ALPHA):
        self.y, self.mu0, self.cl, self.K = y, mu0, cl, K
        self.alpha, self.R = alpha, 1.0 / alpha
        self.cy = gammaln(y + self.R) - gammaln(self.R) - gammaln(y + 1)
        am = N > 0
        self.am = am
        self.Na, self.Ma = N[am], M[am]
        self.rho_r = rho_r
        self.ca = 1.0 / np.minimum(rho_r / self.Na, 0.99) - 1.0
        self.lc = gammaln(self.Na + 1) - gammaln(self.Ma + 1) - gammaln(self.Na - self.Ma + 1)
        self.cA = gammaln(self.ca) - gammaln(self.Na + self.ca)
        self.cl_a = cl[am]


def depth_ll_score(dat, mu):
    """Log-likelihood of y under (1-eps) NB2(mu, alpha) + eps NB2(mu, 1) (outlier-robust background) and its derivative wrt log mu."""
    y = dat.y
    la = dat.cy + dat.R * np.log(dat.R / (dat.R + mu)) + y * (np.log(mu) - np.log(dat.R + mu))
    lw = -np.log1p(mu) + y * (np.log(mu) - np.log1p(mu))              # NB2 with alpha = 1 (gammaln terms cancel)
    lf = np.logaddexp(np.log1p(-EPS_OUT) + la, np.log(EPS_OUT) + lw)
    rw = np.exp(np.log(EPS_OUT) + lw - lf)
    return lf, (1 - rw) * (y - mu) / (1 + dat.alpha * mu) + rw * (y - mu) / (1 + mu)


def bb_ll_dh(dat, h, grad=True):
    a, b = h * dat.ca, (1 - h) * dat.ca
    ll = dat.lc + dat.cA + gammaln(dat.Ma + a) + gammaln(dat.Na - dat.Ma + b) - gammaln(a) - gammaln(b)
    if not grad:
        return ll, None
    dh = dat.ca * (digamma(dat.Ma + a) - digamma(a) - digamma(dat.Na - dat.Ma + b) + digamma(b))
    return ll, dh


class Prior:
    def __init__(self, m_vec, sig_u, sig_b, d_med, tau):
        K = len(m_vec)
        self.m = np.asarray(m_vec, float)
        s = sig_b ** 2 / (sig_u ** 2 + K * sig_b ** 2) if sig_b > 0 else 0.0
        self.Sinv = (np.eye(K) - s * np.ones((K, K))) / sig_u ** 2
        self.d_med, self.tau = d_med, tau


def objective(x, dat, pr, e, k, use_d, use_a, pi_fixed):
    """Return (-penalised loglik, gradient)."""
    K = dat.K
    a = x[:K]; z = x[K]
    pi = pi_fixed if pi_fixed is not None else x[K + 1]
    L, he = TYPES[e] if e else (1.0, 0.5)
    mu = dat.mu0 * np.exp(a[dat.cl])
    n = len(dat.y)
    l0 = np.zeros(n); l1 = np.zeros(n)
    s0 = s1 = None
    if use_d:
        l0, s0 = depth_ll_score(dat, mu)
        l1, s1 = depth_ll_score(dat, mu * L)
    d0 = pr.d_med + pr.tau * z
    h0 = np.clip(0.5 + d0, 0.02, 0.98)
    h1 = np.clip(h0 + k * (he - 0.5), 0.02, 0.98)
    g0 = g1 = None
    if use_a:
        ll0, dh0 = bb_ll_dh(dat, h0); ll1, dh1 = bb_ll_dh(dat, h1)
        l0 = l0.copy(); l1 = l1.copy()
        l0[dat.am] += ll0; l1[dat.am] += ll1
        un0 = 1.0 if 0.02 < 0.5 + d0 < 0.98 else 0.0
        un1 = 1.0 if (0.02 < 0.5 + d0 + k * (he - 0.5) < 0.98) else 0.0
    pen = 0.5 * (a - pr.m) @ pr.Sinv @ (a - pr.m) + 0.5 * z ** 2
    grad = np.zeros_like(x)
    if pi_fixed is not None and pi_fixed == 0.0:
        J = l0.sum(); w = np.zeros(n)
    else:
        pc = min(max(pi, 1e-9), 1 - 1e-9)
        lm = np.logaddexp(np.log1p(-pc) + l0, np.log(pc) + l1)
        J = lm.sum()
        w = np.exp(np.log(pc) + l1 - lm)
        if pi_fixed is None:
            grad[K + 1] = -(w / pc - (1 - w) / (1 - pc)).sum()
    if use_d:
        sc = (1 - w) * s0 + w * s1
        grad[:K] = -np.bincount(dat.cl, weights=sc, minlength=K)
    if use_a:
        wa = w[dat.am]
        grad[K] = -pr.tau * (((1 - wa) * dh0 * un0 + wa * dh1 * un1).sum())
    grad[:K] += pr.Sinv @ (a - pr.m)
    grad[K] += z
    return -(J - pen), grad


def fit(dat, pr, e, k, use_d=True, use_a=True, pi_fixed=None, x0=None, pi0=0.3):
    K = dat.K
    if x0 is None:
        a0 = np.array([np.log((dat.y[dat.cl == c].sum() + 0.5) / ((dat.mu0[dat.cl == c]).sum() + 0.5)) for c in range(K)])
        a0 = 0.5 * a0 + 0.5 * pr.m
        x0 = np.r_[a0, 0.0]
        if pi_fixed is None:
            x0 = np.r_[x0, pi0]
    x0 = np.array(x0, float)
    bnd = [(-4, 4)] * K + [(-6, 6)]
    if not use_d:
        x0[:K] = pr.m; bnd = [(m, m) for m in pr.m] + [(-6, 6)]
    if not use_a:
        x0[K] = 0.0; bnd[K] = (0.0, 0.0)
    if pi_fixed is None:
        bnd = bnd + [(PI_MIN, 1.0)]
    res = optimize.minimize(objective, x0, args=(dat, pr, e, k, use_d, use_a, pi_fixed), jac=True, method="L-BFGS-B", bounds=bnd,
                            options={"maxiter": 200, "ftol": 1e-10, "gtol": 1e-6})
    x = res.x
    return {"ll": -res.fun, "x": x, "pi": (pi_fixed if pi_fixed is not None else float(x[K + 1])), "type": e}


def posterior_w(dat, pr, fitres, k, use_d=True, use_a=True):
    K = dat.K; x = fitres["x"]; e = fitres["type"]; pi = fitres["pi"]
    L, he = TYPES[e]
    mu = dat.mu0 * np.exp(x[:K][dat.cl])
    l0 = np.zeros(len(dat.y)); l1 = np.zeros(len(dat.y))
    if use_d:
        l0 += depth_ll_score(dat, mu)[0]; l1 += depth_ll_score(dat, mu * L)[0]
    if use_a:
        d0 = pr.d_med + pr.tau * x[K]
        h0 = np.clip(0.5 + d0, 0.02, 0.98); h1 = np.clip(h0 + k * (he - 0.5), 0.02, 0.98)
        l0[dat.am] += bb_ll_dh(dat, h0, False)[0]; l1[dat.am] += bb_ll_dh(dat, h1, False)[0]
    pc = min(max(pi, 1e-9), 1 - 1e-9)
    lm = np.logaddexp(np.log1p(-pc) + l0, np.log(pc) + l1)
    return np.exp(np.log(pc) + l1 - lm)


def scan_types(dat, pr, k, use_d=True, use_a=True, null=None):
    """Null fit + all 6 types; returns null, list of alt fits (sorted best first)."""
    if null is None:
        null = fit(dat, pr, None, k, use_d, use_a, pi_fixed=0.0)
    alts = []
    K = dat.K
    for e in TYPES:
        best = None
        for pi0 in (0.3, 0.95):                      # two starts: mostly-background and mostly-event
            x0 = np.r_[null["x"], pi0]
            x0[:K] = x0[:K] - np.log(1 + pi0 * (TYPES[e][0] - 1)) if use_d else x0[:K]
            f = fit(dat, pr, e, k, use_d, use_a, x0=x0)
            if best is None or f["ll"] > best["ll"]:
                best = f
        alts.append(best)
    alts.sort(key=lambda f: -f["ll"])
    return null, alts


def profile_ci(dat, pr, best, k, use_d=True, use_a=True):
    lls = []
    x = best["x"][:-1]
    for pi in PI_GRID:
        f = fit(dat, pr, best["type"], k, use_d, use_a, pi_fixed=float(pi), x0=x)
        lls.append(f["ll"])
    lls = np.array(lls)
    top = max(lls.max(), best["ll"])
    pts = np.r_[PI_GRID, best["pi"]]; vals = np.r_[lls, best["ll"]]
    o = np.argsort(pts); pts, vals = pts[o], vals[o]
    dev = 2 * (top - vals)                     # deviance along pi
    ip = int(np.argmin(np.abs(pts - best["pi"])))
    lo, hi = float(pts[ip]), float(pts[ip])
    j = ip
    while j > 0 and dev[j - 1] <= 3.84:
        j -= 1
    lo = float(pts[j]) if j == 0 else float(pts[j - 1] + (pts[j] - pts[j - 1]) * (dev[j - 1] - 3.84) / (dev[j - 1] - dev[j]))
    j = ip
    while j < len(pts) - 1 and dev[j + 1] <= 3.84:
        j += 1
    hi = float(pts[j]) if j == len(pts) - 1 else float(pts[j] + (pts[j + 1] - pts[j]) * (3.84 - dev[j]) / (dev[j + 1] - dev[j]))
    return lo, hi, lls


def estimate_rho(Nv, Mv):
    am = Nv > 0
    if am.sum() < 5:
        return 0.5
    class Z: pass
    h = Mv[am].sum() / Nv[am].sum()
    h = float(np.clip(h, 0.05, 0.95))
    best, br = -np.inf, 0.3
    for r in np.exp(np.linspace(np.log(0.003), np.log(3.0), 50)):
        z = Z(); z.Na, z.Ma = Nv[am], Mv[am]; z.ca = 1.0 / np.minimum(r / z.Na, 0.99) - 1
        z.lc = gammaln(z.Na + 1) - gammaln(z.Ma + 1) - gammaln(z.Na - z.Ma + 1); z.cA = gammaln(z.ca) - gammaln(z.Na + z.ca)
        ll = bb_ll_dh(z, h, False)[0].sum()
        if ll > best:
            best, br = ll, r
    return float(br)


ALPHA_GRID = np.array([0.003, 0.006, 0.01, 0.015, 0.02, 0.03, 0.045, 0.07, 0.1, 0.15])


def estimate_alpha(Y, mu0, lab_i, K, Nv, Mv, rho_r, pr):
    """Region-specific NB2 dispersion: grid ML of the depth-only null (baseline a_c refitted at each alpha). Held fixed afterwards."""
    best = None
    for al in ALPHA_GRID:
        dat = Dat(Y, mu0, lab_i, K, Nv, Mv, rho_r, al)
        f = fit(dat, pr, None, K_PRIMARY, True, False, pi_fixed=0.0)
        if best is None or f["ll"] > best[0]:
            best = (f["ll"], al, dat)
    return float(best[1]), best[2]


def classify(p, lo, hi):
    if p < 0.05:
        if lo >= 0.8:
            return "clonal"
        if hi < 0.95 and lo > 0:
            return "subclonal"
        return "event, clonal vs subclonal unresolved"
    return "no evidence" if hi <= 0.25 else "insufficient evidence"


def simulate(rng, dat, nullfit, pr, mu0_full=None):
    K = dat.K; x = nullfit["x"]
    mu = dat.mu0 * np.exp(x[:K][dat.cl])
    y = rng.poisson(rng.gamma(dat.R, mu / dat.R))
    d0 = pr.d_med + pr.tau * x[K]
    h0 = float(np.clip(0.5 + d0, 0.02, 0.98))
    Mn = np.zeros(len(dat.y)); Nn = np.zeros(len(dat.y))
    Nn[dat.am] = dat.Na
    q = rng.beta(h0 * dat.ca, (1 - h0) * dat.ca)
    Mn[dat.am] = rng.binomial(dat.Na.astype(int), q)
    return y.astype(float), Mn, Nn


# ============================================================================= per-region analysis
def calib_stats(prep, units_idx, units, lab_i):
    """Per-window calibration summaries (see docstring)."""
    out = []
    for j in units_idx:
        u = units[j]
        nt, nf, nb = len(u["target"]), len(u["flank"]), prep["nb"]
        Y = prep["Y"][:, j]; Fl = prep["F"][:, j]; E = prep["tot"] - Y
        rec = {}
        a_f, a_g = np.zeros(5), np.zeros(5)
        for c in range(5):
            s = lab_i == c
            a_g[c] = np.log((Y[s].sum() + .5) / (E[s].sum() * nt / (nb - nt) + .5))
            a_f[c] = np.log((Y[s].sum() + .5) / (Fl[s].sum() * nt / nf + .5)) if nf else np.nan
        rec["a_g"], rec["a_f"] = a_g, a_f
        N = prep["N"][:, j]; M = prep["M"][:, j]
        rec["d_all"] = M.sum() / N.sum() - 0.5 if N.sum() > 0 else np.nan
        rec["d_c"] = np.array([M[lab_i == c].sum() / N[lab_i == c].sum() - 0.5 if N[lab_i == c].sum() > 0 else np.nan for c in range(5)])
        out.append(rec)
    return out


def priors_from_calib(cs, mode):
    key = "a_f" if mode == "flank" else "a_g"
    A = np.array([r[key] for r in cs])                      # windows x 5
    mean_w = A.mean(1)
    m = float(np.median(mean_w)); sig_b = mad(mean_w)
    dev = (A - mean_w[:, None]).ravel()
    sig_u = mad(dev) * np.sqrt(5 / 4)                       # deviations from a 5-mean lose 1/5 of the variance
    dall = np.array([r["d_all"] for r in cs]); dc = np.array([r["d_c"] for r in cs])
    return {"m": m, "sig_b": max(sig_b, 0.02), "sig_u": max(sig_u, 0.02), "d_med": float(np.nanmedian(dall)), "tau": max(mad(dall[np.isfinite(dall)]), 0.003),
            "d_med_c": np.nanmedian(dc, axis=0), "tau_c": np.maximum([mad(dc[:, c][np.isfinite(dc[:, c])]) for c in range(5)], 0.003)}


def build_dat(u, prep, j, lab_i, mode, y=None, M=None, N=None, rho_r=None, cells=None):
    nb = prep["nb"]; nt, nf = len(u["target"]), len(u["flank"])
    Y = prep["Y"][:, j] if y is None else y
    E = prep["tot"] - prep["Y"][:, j]
    Fl = prep["F"][:, j]
    g = np.zeros(5)
    if mode == "flank":
        for c in range(5):
            s = lab_i == c
            g[c] = np.log((Fl[s].sum() + .5) / (E[s].sum() * nf / (nb - nt) + .5))
    mu0 = E * (nt / (nb - nt)) * np.exp(g[lab_i])
    Nv = prep["N"][:, j] if N is None else N
    Mv = prep["M"][:, j] if M is None else M
    return Y, mu0, Nv, Mv


def analyse_region(task):
    """Full analysis of one region (or random window): overall + (optionally) cluster fits, bootstrap, variants."""
    (u, j, prep_slice, lab_i, mode, cal, ridx, opts) = task
    t0 = time.time()
    K = 5
    prep = prep_slice
    Y, mu0, Nv, Mv = build_dat(u, prep, j, lab_i, mode)
    rho_r = estimate_rho(Nv, Mv)
    k = K_PRIMARY
    m_vec = np.full(K, cal["m"])
    pr = Prior(m_vec, cal["sig_u"], cal["sig_b"], cal["d_med"], cal["tau"])
    pr_noff = Prior(m_vec, cal["sig_u"], 1e-4, cal["d_med"], cal["tau"])
    alpha_r, dat = estimate_alpha(Y, mu0, lab_i, K, Nv, Mv, rho_r, pr)
    res = {"name": u["name"], "family": u["family"], "mode": mode, "rho_r": rho_r, "alpha_r": alpha_r, "n_cells": len(Y), "n_cells_N_gt0": int((Nv > 0).sum()), "sum_N": float(Nv.sum()),
           "pooled_h": float(Mv.sum() / Nv.sum()) if Nv.sum() > 0 else np.nan}
    null, alts = scan_types(dat, pr, k)
    best = alts[0]
    T = max(0.0, 2 * (best["ll"] - null["ll"]))
    lo, hi, _ = profile_ci(dat, pr, best, k)
    w = posterior_w(dat, pr, best, k)
    comp = {}
    for nm, (ud, ua) in {"depth_only": (True, False), "allele_only": (False, True)}.items():
        n_, a_ = scan_types(dat, pr, k, ud, ua)
        comp[nm] = {"type_best": a_[0]["type"], "pi": a_[0]["pi"], "LRT": max(0.0, 2 * (a_[0]["ll"] - n_["ll"])),
                    "LRT_for_joint_best_type": max(0.0, 2 * (next(f for f in a_ if f["type"] == best["type"])["ll"] - n_["ll"])),
                    "pi_for_joint_best_type": next(f for f in a_ if f["type"] == best["type"])["pi"]}
    res.update(best_type=best["type"], pi=best["pi"], pi_lo=lo, pi_hi=hi, LRT=T, frac_cells_P90=float((w > 0.9).mean()), mean_P_event=float(w.mean()),
               a_hat_null=null["x"][:K].tolist(), d0_null=float(cal["d_med"] + cal["tau"] * null["x"][K]), d0_alt=float(cal["d_med"] + cal["tau"] * best["x"][K]),
               b_hat_alt=float(best["x"][:K].mean()), all_types={f["type"]: {"pi": f["pi"], "LRT": max(0.0, 2 * (f["ll"] - null["ll"]))} for f in alts},
               components=comp, prior=cal_summary(cal))
    # ---- bootstrap (overall + cluster)
    rng = np.random.default_rng(SEED18 + (ridx,))
    clus_data = None
    do_cl = opts["clusters"]
    abar_null = float(null["x"][:K].mean())
    cl_prs, cl_dats = [], []
    if do_cl:
        for c in range(K):
            s = lab_i == c
            prc = Prior([abar_null], cal["sig_u"], 1e-4, cal["d_med_c"][c], cal["tau_c"][c])
            cl_prs.append(prc)
            cl_dats.append((s, Dat(Y[s], mu0[s], np.zeros(s.sum(), int), 1, Nv[s], Mv[s], rho_r, alpha_r), prc))
        crow = []
        for c in range(K):
            s, dc, prc = cl_dats[c]
            nl, al = scan_types(dc, prc, k)
            b = al[0]
            lo_c, hi_c, _ = profile_ci(dc, prc, b, k)
            wc = posterior_w(dc, prc, b, k)
            sub = {}
            for nm, (ud, ua) in {"depth_only": (True, False), "allele_only": (False, True)}.items():
                n_, a_ = scan_types(dc, prc, k, ud, ua)
                f = next(f for f in a_ if f["type"] == b["type"])
                sub[nm] = max(0.0, 2 * (f["ll"] - n_["ll"]))
            crow.append({"cluster": CLUST[c], "n_cells": int(s.sum()), "n_cells_N_gt0": int((Nv[s] > 0).sum()), "best_type": b["type"], "pi": b["pi"], "pi_lo": lo_c,
                         "pi_hi": hi_c, "LRT": max(0.0, 2 * (b["ll"] - nl["ll"])), "frac_cells_P90": float((wc > 0.9).mean()), "LRT_depth_only": sub["depth_only"],
                         "LRT_allele_only": sub["allele_only"], "pooled_h_all_orient": float(Mv[s].sum() / Nv[s].sum()) if Nv[s].sum() > 0 else np.nan,
                         "pooled_h_own_orient": float(prep["Mo"][s, opts["own_index"]].sum() / prep["No"][s, opts["own_index"]].sum()) if opts["own_index"] is not None and prep["No"][s, opts["own_index"]].sum() > 0 else np.nan})
        res["clusters"] = crow
    boots = {"main": [], "k1": None}
    Tb = []; Tc = [[] for _ in range(K)]
    for b_ in range(opts["B"]):
        yb, Mb, Nb = simulate(rng, dat, null, pr)
        db = Dat(yb, mu0, lab_i, K, Nb, Mb, rho_r, alpha_r)
        nb_, ab_ = scan_types(db, pr, k)
        Tb.append(max(0.0, 2 * (ab_[0]["ll"] - nb_["ll"])))
        if do_cl:
            for c in range(K):
                s, dc, prc = cl_dats[c]
                dcb = Dat(yb[s], mu0[s], np.zeros(s.sum(), int), 1, Nb[s], Mb[s], rho_r, alpha_r)
                nl, al = scan_types(dcb, prc, k)
                Tc[c].append(max(0.0, 2 * (al[0]["ll"] - nl["ll"])))
    Tb = np.array(Tb)
    res["boot_p"] = float((1 + (Tb >= T).sum()) / (len(Tb) + 1)); res["boot_B"] = len(Tb); res["boot_T_q95"] = float(np.quantile(Tb, .95)) if len(Tb) else np.nan
    res["class"] = classify(res["boot_p"], lo, hi)
    if do_cl:
        for c in range(K):
            Tcc = np.array(Tc[c]); r = res["clusters"][c]
            r["boot_p"] = float((1 + (Tcc >= r["LRT"]).sum()) / (len(Tcc) + 1)); r["class"] = classify(r["boot_p"], r["pi_lo"], r["pi_hi"])
    # ---- depth/allele agreement for the joint best type
    dsup = comp["depth_only"]["LRT_for_joint_best_type"] > 3.84; asup = comp["allele_only"]["LRT_for_joint_best_type"] > 3.84
    L_e, h_e = TYPES[best["type"]]
    d_info = abs(np.log(L_e)) > 0.05; a_info = abs(h_e - 0.5) > 0.05
    res["agreement"] = ("depth and allele both support" if dsup and asup else "depth only" if dsup and not asup else "allele only" if asup and not dsup else "neither alone")
    res["agreement_note"] = ("" if d_info else "type has no depth change; ") + ("" if a_info else "type has no allele change; ")
    # ---- variants: k=1.0, no offset (LRT, pi, type; bootstrap p for 'no offset' and 'depth only' where requested)
    var = []
    for vname, (prv, kv, ud, ua) in {"k=1.0": (pr, K_ALT, True, True), "no shared baseline offset": (pr_noff, k, True, True),
                                       "depth only": (pr, k, True, False), "allele only": (pr, k, False, True)}.items():
        nv, av = scan_types(dat, prv, kv, ud, ua)
        bv = av[0]
        Tv = max(0.0, 2 * (bv["ll"] - nv["ll"]))
        pv = np.nan
        if opts["var_boot"] and vname in opts["var_boot"]:
            rngv = np.random.default_rng(SEED18 + (ridx, 7))
            Tvb = []
            for _ in range(opts["B_var"]):
                yb, Mb, Nb = simulate(rngv, dat, nv, prv)
                db = Dat(yb, mu0, lab_i, K, Nb, Mb, rho_r, alpha_r)
                n2, a2 = scan_types(db, prv, kv, ud, ua)
                Tvb.append(max(0.0, 2 * (a2[0]["ll"] - n2["ll"])))
            pv = float((1 + (np.array(Tvb) >= Tv).sum()) / (len(Tvb) + 1))
        lo_v, hi_v, _ = profile_ci(dat, prv, bv, kv, ud, ua) if opts["var_boot"] and vname in opts["var_boot"] else (np.nan, np.nan, None)
        var.append({"variant": vname, "best_type": bv["type"], "pi": bv["pi"], "LRT": Tv, "boot_p": pv, "pi_lo": lo_v, "pi_hi": hi_v,
                    "class": classify(pv, lo_v, hi_v) if np.isfinite(pv) else ""})
    res["variants"] = var
    res["seconds"] = time.time() - t0
    return res


def cal_summary(cal):
    return {k: (v if np.isscalar(v) else list(map(float, v))) for k, v in cal.items() if k in ("m", "sig_b", "sig_u", "d_med", "tau", "d_med_c", "tau_c", "n_windows")}


# ============================================================================= driver
def main():
    t0 = time.time()
    p = C.default_paths()
    d = C.H5AD(p["h5ad"])
    lab = C.cell_clusters(d, p)
    lab_i = np.array([CLUST.index(str(x)) for x in lab])
    kseg = pd.read_csv(p["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    call = np.ones(d.n_bins)
    for r in kseg[kseg.group == "all"].itertuples():
        call[r.bin_start:r.bin_end + 1] = r.cn_call
    regs = make_regions(d, p)
    taken = np.zeros(d.n_bins, bool)
    for r in regs:
        taken[r["target"]] = True; taken[r["flank"]] = True
    wins, n_cand = make_windows(d, call, taken, N_PANEL + N_CALIB, SEED18)
    units = regs + wins
    tabs = orientation_tables(d)
    own_units = list(range(len(regs)))
    print(f"[{time.time()-t0:.0f}s] {len(regs)} regions, {len(wins)} windows (of {n_cand} candidates); preparing per-cell sums", flush=True)
    prep = prepare(d, p, lab, units, tabs, own_units)
    prep["nb"] = d.n_bins
    print(f"[{time.time()-t0:.0f}s] prepared", flush=True)
    panel = list(range(len(regs), len(regs) + N_PANEL)); calib = list(range(len(regs) + N_PANEL, len(units)))
    cs = calib_stats(prep, calib, units, lab_i)
    cal_f = priors_from_calib(cs, "flank"); cal_g = priors_from_calib(cs, "genome")
    cal_f["n_windows"] = cal_g["n_windows"] = len(calib)
    print("calibration flank:", {k: np.round(v, 4) for k, v in cal_f.items()}, flush=True)
    print("calibration genome:", {k: np.round(v, 4) for k, v in cal_g.items()}, flush=True)
    tasks = []
    for j, u in enumerate(regs):
        cal = cal_f if u["mode"] == "flank" else cal_g
        tasks.append((u, 0, prep_for(prep, j, own=True), lab_i, u["mode"], cal, j, {"clusters": True, "B": B_MAIN, "B_var": B_VAR, "var_boot": [], "own_index": 0}))
    for jj, j in enumerate(panel):
        u = units[j]
        tasks.append((u, 0, prep_for(prep, j), lab_i, "flank", cal_f, 1000 + jj,
                      {"clusters": False, "B": B_MAIN, "B_var": B_VAR, "var_boot": ["no shared baseline offset", "depth only"], "own_index": None}))
    nproc = 1 if CHECK else 2
    with Pool(nproc) as pool:
        results = pool.map(analyse_region, tasks, chunksize=1)
    print(f"[{time.time()-t0:.0f}s] analysis done", flush=True)
    return finish(d, regs, wins, panel, cal_f, cal_g, results, lab_i, call, n_cand, t0)


def prep_for(prep, j, own=False):
    """Slim copy of the per-cell sums needed for unit j (keeps task pickles small)."""
    s = {"nb": prep["nb"], "tot": prep["tot"], "Y": prep["Y"][:, [j]], "F": prep["F"][:, [j]], "M": prep["M"][:, [j]], "N": prep["N"][:, [j]]}
    s["Mo"] = prep["Mo"][:, [j]] if own else None
    s["No"] = prep["No"][:, [j]] if own else None
    return s


def finish(d, regs, wins, panel, cal_f, cal_g, results, lab_i, call, n_cand, t0):
    import json
    C.write_json(OUT / "raw_results.json", results)
    return results


if __name__ == "__main__":
    main()
