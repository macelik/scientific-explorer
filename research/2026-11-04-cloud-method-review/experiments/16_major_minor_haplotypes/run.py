"""Experiment 16: major / minor haplotype inference with a phase-orientation HMM (exploratory diagnostic).

Goal: in regions with allelic imbalance, estimate which haplotype is MAJOR in each bin when the orientation of the
GT-based phase between bins is not verified, so that depth and allele evidence can later be combined into
allele-specific relative copy number. Major / minor labels are INFERRED from allele counts; they are not verified phase.

Data: HA / HB / AT parquet (bins x barcodes; AT = summed SNP-site DP over selected phased heterozygous SNPs, HA + HB = AT;
Eagle GT phasing, no phase-set IDs). Joined to AnnData by barcode name. WNN clusters via common.cell_clusters.

Model (per segment, pooled over the cells of a group; see report "Method"):
  hidden orientation o_j in {B-major, A-major}; HB_j | AT_j, o_j ~ BetaBinomial(AT_j, 0.5 +/- delta, rho)
  (+ for B-major, - for A-major); consecutive informative bins switch orientation with probability
  s_j = 0.5 (1 - exp(-2 lambda d_j)), d_j = bp distance between the bin starts, lambda = 1 switch per 5 Mb (prespecified,
  not fitted); sensitivity lambda x 0.1 and x 10. delta in [0, 0.49] by profile likelihood on a grid; rho = 0.0038 from exp 06.
Seeds: numpy PCG64 default_rng((20261104, 16, stream, ...)); streams 1 segments, 2 regions, 3 random null windows
(random windows themselves: E4.build_regions, seed (20261104, 4)).

Run:  /home/user/envs/cnv310/bin/python -I research/2026-11-04-cloud-method-review/experiments/16_major_minor_haplotypes/run.py [--check]
"""
import importlib.util
import json
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
from scipy import sparse  # noqa: E402
from scipy.special import gammaln  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

spec = importlib.util.spec_from_file_location("exp04", HERE.parent / "04_fixed_region_contrasts" / "run.py")
E4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(E4)

CHECK = "--check" in sys.argv
OUT = Path("/tmp/claude-0/-home-user-scientific-explorer/b1079aa7-80fd-5756-88cc-d9253976037e/scratchpad/exp16_check") if CHECK else HERE
OUT.mkdir(parents=True, exist_ok=True)

SEED16 = (C.SEED, 16)
LAMBDA0 = 1.0 / 5e6            # prespecified switch hazard: one orientation switch per 5 Mb
MIN_AT_BIN = 5                 # bin counts as "well covered"
MIN_BINS_AT5 = 8               # segment eligibility: bins with AT >= 5
D_FINE = np.round(np.arange(0.0, 0.4901, 0.0025), 4)
D_COARSE = np.round(np.arange(0.0, 0.4901, 0.01), 4)
B_SEG = 20 if CHECK else 200
B_REG = 50 if CHECK else 200
N_NULL = 12 if CHECK else 60
CONFIGS = {"base": {"rho": None, "lam": LAMBDA0},          # rho filled from exp 06
           "hazard_x0.1": {"rho": None, "lam": LAMBDA0 * 0.1},
           "hazard_x10": {"rho": None, "lam": LAMBDA0 * 10},
           "rho_0": {"rho": 0.0, "lam": LAMBDA0},
           "rho_0.02": {"rho": 0.02, "lam": LAMBDA0}}
GROUPS = ["all", "cluster0", "cluster1", "cluster2", "cluster3", "cluster4"]
G = {}          # group -> dict(HB, AT csr bins x cells, nc)
BINS = {}       # seq, start, end


# ----------------------------------------------------------------------------- HMM
def log_emissions(k, n, deltas, rho):
    """log P(HB=k | AT=n) for B-major and A-major, up to terms that do not depend on delta or orientation.
    k, n: (..., 1) arrays; deltas: (G,). Returns two arrays (..., G)."""
    if rho <= 0:
        p = 0.5 + deltas
        lp, lq = np.log(p), np.log(1 - p)
        return k * lp + (n - k) * lq, k * lq + (n - k) * lp
    c = (1 - rho) / rho
    a, b = (0.5 + deltas) * c, (0.5 - deltas) * c
    const = gammaln(a) + gammaln(b) - gammaln(c)
    return (gammaln(k + a) + gammaln(n - k + b) - const, gammaln(k + b) + gammaln(n - k + a) - const)


def switch_prob(dist, lam):
    return 0.5 * (1 - np.exp(-2 * lam * dist))


def forward_ll(K, N, dist, deltas, rho, lam):
    """Marginal log-likelihood (up to constants) for every (column of K, delta). K, N: (nbins, ncol)."""
    nb, ncol = K.shape
    ll = np.zeros((ncol, len(deltas)))
    aB = np.full((ncol, len(deltas)), 0.5)
    aA = aB.copy()
    s = switch_prob(dist, lam)
    for j in range(nb):
        if j > 0:
            pB, pA = aB * (1 - s[j - 1]) + aA * s[j - 1], aA * (1 - s[j - 1]) + aB * s[j - 1]
        else:
            pB, pA = aB, aA
        eB, eA = log_emissions(K[j][:, None], N[j][:, None], deltas, rho)
        m = np.maximum(eB, eA)
        wB, wA = pB * np.exp(eB - m), pA * np.exp(eA - m)
        tot = wB + wA
        ll += np.log(tot) + m
        aB, aA = wB / tot, wA / tot
    return ll


def posterior(k, n, dist, delta, rho, lam):
    """Forward-backward and Viterbi at a fixed delta. Returns P(B-major) per bin, expected switches, Viterbi switches."""
    nb = len(k)
    eB, eA = log_emissions(k[:, None].astype(float), n[:, None].astype(float), np.array([delta]), rho)
    eB, eA = eB[:, 0], eA[:, 0]
    m = np.maximum(eB, eA)
    E = np.stack([np.exp(eB - m), np.exp(eA - m)], 1)
    s = switch_prob(dist, lam)
    al = np.zeros((nb, 2))
    c = np.zeros(nb)
    prev = np.array([0.5, 0.5])
    for j in range(nb):
        pr = prev if j == 0 else np.array([prev[0] * (1 - s[j - 1]) + prev[1] * s[j - 1], prev[1] * (1 - s[j - 1]) + prev[0] * s[j - 1]])
        w = pr * E[j]
        c[j] = w.sum()
        al[j] = w / c[j]
        prev = al[j]
    be = np.ones((nb, 2))
    for j in range(nb - 1, 0, -1):
        sj = s[j - 1]
        v = E[j] * be[j]
        be[j - 1] = np.array([(1 - sj) * v[0] + sj * v[1], sj * v[0] + (1 - sj) * v[1]]) / c[j]
    gam = al * be
    gam /= gam.sum(1, keepdims=True)
    exp_sw = 0.0
    for j in range(1, nb):
        sj = s[j - 1]
        x01 = al[j - 1, 0] * sj * E[j, 1] * be[j, 1] / c[j]
        x10 = al[j - 1, 1] * sj * E[j, 0] * be[j, 0] / c[j]
        exp_sw += x01 + x10
    # Viterbi
    le = np.stack([eB, eA], 1)
    dp = np.zeros((nb, 2))
    bp = np.zeros((nb, 2), int)
    dp[0] = np.log(0.5) + le[0]
    for j in range(1, nb):
        sj = np.clip(s[j - 1], 1e-300, 1 - 1e-12)
        for st in (0, 1):
            stay = dp[j - 1, st] + np.log(1 - sj)
            sw = dp[j - 1, 1 - st] + np.log(sj)
            bp[j, st] = st if stay >= sw else 1 - st
            dp[j, st] = max(stay, sw) + le[j, st]
    path = np.zeros(nb, int)
    path[-1] = int(np.argmax(dp[-1]))
    for j in range(nb - 1, 0, -1):
        path[j - 1] = bp[j, path[j]]
    return gam[:, 0], float(exp_sw), int((np.diff(path) != 0).sum())


def argmax_delta(ll, deltas, parab=False):
    """ll: (ncol, G). Returns delta-hat per column (optionally with a parabolic vertex refinement) ."""
    i = ll.argmax(1)
    d = deltas[i].astype(float)
    if parab:
        ok = (i > 0) & (i < len(deltas) - 1)
        ii = np.where(ok, i, 1)
        y0, y1, y2 = ll[np.arange(len(i)), ii - 1], ll[np.arange(len(i)), ii], ll[np.arange(len(i)), ii + 1]
        den = y0 - 2 * y1 + y2
        off = np.where((den < 0) & ok, 0.5 * (y0 - y2) / np.where(den < 0, den, -1), 0.0)
        d = np.clip(d + off * (deltas[1] - deltas[0]), 0, 0.49)
    return d


def delta_eff_cols(K, N, rho):
    """Phase-agnostic effective |dBAF| (exp 06 formula) for every column, over bins with N >= MIN_AT_BIN."""
    m = N >= MIN_AT_BIN
    Ns = np.where(m, N, 1.0)
    z2 = (2 * K - N) ** 2 / Ns
    e = (z2 - (1 + (Ns - 1) * rho)) / (4 * Ns)
    cnt = m.sum(0)
    mean = np.where(cnt > 0, (e * m).sum(0) / np.maximum(cnt, 1), np.nan)
    return np.sqrt(np.maximum(0.0, mean))


# ----------------------------------------------------------------------------- tasks
def run_task(t):
    g = G[t["group"]]
    cfg = t["cfg"]
    rho, lam = cfg["rho"], cfg["lam"]
    idx = np.asarray(t["bins"])
    shb, sat = g["HB"][idx], g["AT"][idx]
    hb = np.asarray(shb.sum(1)).ravel()
    at = np.asarray(sat.sum(1)).ravel()
    inf = at > 0
    rec = {k: t[k] for k in ("kind", "group", "name", "cfg_name")}
    rec.update({"n_bins": len(idx), "n_cells": g["nc"], "bins_AT_gt0": int(inf.sum()), "bins_AT_ge5": int((at >= MIN_AT_BIN).sum()),
                "AT_sum": float(at.sum()), "HB_sum": float(hb.sum())})
    rec["eligible"] = bool(rec["bins_AT_ge5"] >= MIN_BINS_AT5)
    ii = np.flatnonzero(inf)
    if rec["AT_sum"] > 0:
        f = hb.sum() / at.sum()
        rec.update(f_naive=f, h_naive=max(f, 1 - f), delta_naive=abs(f - 0.5))
    if not rec["eligible"]:
        return rec, None
    k, n = hb[ii], at[ii]
    dist = np.diff(BINS["start"][idx[ii]]).astype(float)
    ll = forward_ll(k[:, None], n[:, None], dist, D_FINE, rho, lam)[0]
    i0 = int(ll.argmax())
    dh = float(D_FINE[i0])
    pB, esw, vsw = posterior(k, n, dist, dh, rho, lam)
    rec.update(delta_hat=dh, h_hat=0.5 + dh, llr_vs_balanced=float(2 * (ll[i0] - ll[0])),
               span_bp=float(BINS["end"][idx[-1]] - BINS["start"][idx[0]]),
               delta_eff=float(delta_eff_cols(k[:, None], n[:, None], rho)[0]))
    rec["expected_switches"] = esw if dh > 0 else np.nan
    rec["viterbi_switches"] = vsw if dh > 0 else np.nan
    rec["mean_orientation_certainty"] = float(np.mean(np.maximum(pB, 1 - pB))) if dh > 0 else np.nan
    # cell bootstrap
    B = t["B"]
    if B > 0:
        rng = np.random.default_rng(t["seed"])
        nc = g["nc"]
        W = rng.multinomial(nc, np.full(nc, 1.0 / nc), size=B).T.astype(float)
        Kb = shb[ii] @ W
        Nb = sat[ii] @ W
        llb = forward_ll(Kb, Nb, dist, D_COARSE, rho, lam)
        db = argmax_delta(llb, D_COARSE, parab=True)
        eb = delta_eff_cols(Kb, Nb, rho)
        fb = Kb.sum(0) / np.maximum(Nb.sum(0), 1)
        nb_ = np.abs(fb - 0.5)
        lo, hi = np.quantile(db, [0.025, 0.975])
        rec.update(h_ci_lo=0.5 + lo, h_ci_hi=0.5 + hi, delta_eff_ci_lo=float(np.quantile(eb, 0.025)), delta_eff_ci_hi=float(np.quantile(eb, 0.975)),
                   delta_naive_ci_lo=float(np.quantile(nb_, 0.025)), delta_naive_ci_hi=float(np.quantile(nb_, 0.975)), boot_n=B)
    # per-bin table
    bt = pd.DataFrame({"bin": idx[ii], "HB": k, "AT": n, "p_B_major": pB})
    return rec, bt


def init_groups(paths, d, lab):
    cells = list(d.cells)
    mats = {}
    for nm in ("HA", "HB", "AT"):
        t = pq.read_table(paths["allelic"] / f"{nm}.parquet")
        df_idx = t.select([c for c in t.column_names if c in ("chr", "start", "end")]).to_pandas()
        if "chr" not in df_idx.columns:
            df_idx = t.to_pandas()[[]].reset_index()
        ok = (df_idx["chr"].astype(str).to_numpy() == d.seq).all() and (df_idx["start"].to_numpy() == d.start).all() and (df_idx["end"].to_numpy() == d.end).all()
        names = [c for c in t.column_names if c not in ("chr", "start", "end", "__index_level_0__")]
        if not ok or set(names) != set(cells):
            raise SystemExit(f"alignment failure in {nm}")
        mats[nm] = sparse.csr_matrix(np.column_stack([t.column(c).to_numpy(zero_copy_only=False) for c in cells]).astype(np.float64))
        del t
    checks = {"coords_match_anndata": True, "barcode_set_equal": True, "AT_equals_HA_plus_HB": bool((mats["AT"] - mats["HA"] - mats["HB"]).count_nonzero() == 0),
              "nonnegative": bool(mats["HA"].min() >= 0 and mats["HB"].min() >= 0), "n_nonzero_AT_entries": int(mats["AT"].nnz)}
    for gname in GROUPS:
        sel = np.ones(len(cells), bool) if gname == "all" else (lab == gname.replace("cluster", ""))
        ix = np.flatnonzero(sel)
        G[gname] = {"HB": mats["HB"][:, ix].tocsr(), "AT": mats["AT"][:, ix].tocsr(), "nc": len(ix), "cells": ix}
    return checks


def build_tasks(d, paths, call, kseg, regs):
    tasks = []
    # segments (all groups), config base, bootstrap
    for gi, gname in enumerate(GROUPS):
        s = kseg[kseg.group == gname].reset_index(drop=True)
        if CHECK:
            s = s.head(100)
        for r in s.itertuples():
            tasks.append(dict(kind="segment", group=gname, name=f"{gname}|{r.segment}", cfg_name="base", cfg=None, bins=np.arange(r.bin_start, r.bin_end + 1),
                              B=B_SEG, seed=SEED16 + (1, gi, int(r.segment))))
    # regions: saved + literature in all groups with bootstrap
    for ri, r in enumerate(regs):
        if r["family"] not in ("saved", "literature"):
            continue
        for gi, gname in enumerate(GROUPS):
            tasks.append(dict(kind="region", group=gname, name=r["name"], cfg_name="base", cfg=None, bins=r["target"], B=B_REG, seed=SEED16 + (2, gi, ri)))
    # null windows
    rnd = [r for r in regs if r["family"] == "random"]
    unc = [r for r in rnd if (call[r["target"]] == 1).all()]
    nulls = unc[:N_NULL]
    for wi, r in enumerate(nulls):
        for gi, gname in enumerate(GROUPS):
            tasks.append(dict(kind="null", group=gname, name=r["name"], cfg_name="base", cfg=None, bins=r["target"], B=0, seed=SEED16 + (3, gi, wi)))
    return tasks, len(rnd), len(unc), nulls


def sensitivity_tasks(base_tasks):
    out = []
    for cn in ("hazard_x0.1", "hazard_x10", "rho_0", "rho_0.02"):
        for t in base_tasks:
            if t["group"] != "all":
                continue
            u = dict(t)
            u["cfg_name"], u["B"] = cn, 0
            out.append(u)
    return out


def execute(tasks, nproc=2):
    for t in tasks:
        t["cfg"] = CONFIGS[t["cfg_name"]]
    with Pool(nproc) as p:
        res = p.map(run_task, tasks, chunksize=4)
    return res


# ----------------------------------------------------------------------------- analysis
def q95(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.quantile(x, 0.95)) if len(x) else np.nan


def cross_cluster(seg_all, bin_tabs, imb_flag, ks):
    """Per-bin posterior orientation agreement between groups over imbalanced, confidently oriented bins."""
    nb = BINS["n_bins"]
    P = {g: np.full(nb, np.nan) for g in GROUPS}
    I = {g: np.zeros(nb, bool) for g in GROUPS}
    for (g, name), bt in bin_tabs.items():
        if bt is None or not imb_flag.get((g, name), False):
            continue
        P[g][bt["bin"].to_numpy()] = bt["p_B_major"].to_numpy()
        I[g][bt["bin"].to_numpy()] = True
    conf = {g: I[g] & ((P[g] >= 0.8) | (P[g] <= 0.2)) for g in GROUPS}
    side = {g: P[g] > 0.5 for g in GROUPS}
    pair = []
    cl = GROUPS[1:]
    for a in range(len(GROUPS)):
        for b in range(a + 1, len(GROUPS)):
            ga, gb = GROUPS[a], GROUPS[b]
            m = conf[ga] & conf[gb]
            ag = int((side[ga][m] == side[gb][m]).sum())
            pair.append({"group_a": ga, "group_b": gb, "bins_compared": int(m.sum()), "bins_same_major": ag, "bins_opposite_major": int(m.sum()) - ag,
                         "agree_fraction": ag / m.sum() if m.sum() else np.nan})
    rows = []
    for r in seg_all.itertuples():
        if not r.imbalanced:
            continue
        bt = bin_tabs[("all", r.name)]
        sel = np.zeros(nb, bool)
        sel[bt["bin"].to_numpy()] = True
        for c in cl:
            m = sel & conf["all"] & conf[c]
            n = int(m.sum())
            ag = float((side["all"][m] == side[c][m]).mean()) if n else np.nan
            cls = "insufficient (<10 bins)" if n < 10 else ("same major haplotype" if ag >= 0.8 else ("opposite major haplotype" if ag <= 0.2 else "mixed"))
            rows.append({"all_segment": r.name, "chromosome": r.chromosome, "start_bp": r.start_bp, "end_bp": r.end_bp, "cn_call_all": r.cn_call,
                         "all_h_hat": r.h_hat, "cluster": c, "cluster_imbalanced_bins_in_segment": int((sel & I[c]).sum()),
                         "bins_compared": n, "agree_fraction": ag, "class": cls})
    return pd.DataFrame(pair), pd.DataFrame(rows)


def percell(segs_imb, bin_tabs, lab):
    """Per-cell major-haplotype fraction in imbalanced `all` segments (orientation from the `all` posterior)."""
    rng = np.random.default_rng(SEED16 + (5,))
    g = G["all"]
    rows, store = [], {}
    for r in segs_imb.itertuples():
        bt = bin_tabs[("all", r.name)]
        idx = bt["bin"].to_numpy()
        mB = (bt["p_B_major"].to_numpy() >= 0.5).astype(float)
        hbs, ats = g["HB"][idx], g["AT"][idx]
        D = sparse.diags(mB)
        I = sparse.identity(len(idx), format="csr")
        kmaj = np.asarray((D @ hbs).sum(0) + ((I - D) @ (ats - hbs)).sum(0)).ravel()
        n = np.asarray(ats.sum(0)).ravel()
        for cname in ["all"] + [str(i) for i in range(5)]:
            sel = (n >= 5) & (np.ones(len(n), bool) if cname == "all" else lab == cname)
            if sel.sum() < 5:
                continue
            k, nn = kmaj[sel], n[sel]
            h = k.sum() / nn.sum()
            fr = k / nn
            phi = float((((k - nn * h) ** 2) / (nn * h * (1 - h))).sum() / (sel.sum() - 1))
            bs = rng.integers(0, sel.sum(), size=(200, sel.sum()))
            hb_ = k[bs].sum(1) / nn[bs].sum(1)
            # approximate parametric expectation of phi: beta-binomial at the pooled h with rho
            c_ = (1 - CONFIGS["base"]["rho"]) / CONFIGS["base"]["rho"]
            ph = []
            for _ in range(100):
                q = rng.beta(h * c_, (1 - h) * c_, size=len(nn))
                ks = rng.binomial(nn.astype(int), q)
                hs = ks.sum() / nn.sum()
                ph.append((((ks - nn * hs) ** 2) / (nn * hs * (1 - hs))).sum() / (len(nn) - 1))
            rows.append({"segment": r.name, "chromosome": r.chromosome, "start_bp": r.start_bp, "end_bp": r.end_bp, "cluster": cname,
                         "cells_AT_ge5": int(sel.sum()), "median_AT_per_cell": float(np.median(nn)), "h_pooled_oriented": h,
                         "h_ci_lo": float(np.quantile(hb_, .025)), "h_ci_hi": float(np.quantile(hb_, .975)),
                         "cell_fraction_median": float(np.median(fr)), "cell_fraction_q25": float(np.quantile(fr, .25)),
                         "cell_fraction_q75": float(np.quantile(fr, .75)), "frac_cells_major_below_half": float((fr < 0.5).mean()),
                         "dispersion_phi": phi, "dispersion_phi_expected_bb": float(np.mean(ph)), "dispersion_phi_bb_q95": float(np.quantile(ph, .95))})
            store[(r.name, cname)] = fr
    return pd.DataFrame(rows), store


def main():
    t0 = time.time()
    paths = C.default_paths()
    d = C.H5AD(paths["h5ad"])
    lab = C.cell_clusters(d, paths).astype(str)
    checks = init_groups(paths, d, lab)
    BINS.update(seq=d.seq, start=d.start, end=d.end, n_bins=d.n_bins)
    rho = json.loads((HERE.parent / "06_allelic_information" / "summary.json").read_text())["rho_used"]
    for k in ("base", "hazard_x0.1", "hazard_x10"):
        CONFIGS[k]["rho"] = rho
    kseg = pd.read_csv(paths["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    call = np.ones(d.n_bins)
    for r in kseg[kseg.group == "all"].itertuples():
        call[r.bin_start:r.bin_end + 1] = r.cn_call
    regs = E4.build_regions(d, paths)
    tasks, n_rand, n_unc, nulls = build_tasks(d, paths, call, kseg, regs)
    print(f"tasks {len(tasks)}; random windows {n_rand}, uncalled {n_unc}, null used {len(nulls)}", flush=True)
    res = execute(tasks)
    sens = sensitivity_tasks(tasks)
    sres = execute(sens)
    print(f"computed in {time.time() - t0:.0f}s", flush=True)

    recs = pd.DataFrame([r for r, _ in res])
    bin_tabs = {(r["group"], r["name"]): bt for (r, bt), t in zip(res, tasks) if t["kind"] == "segment"}
    reg_tabs = {(r["group"], r["name"]): bt for (r, bt), t in zip(res, tasks) if t["kind"] == "region"}
    srecs = pd.DataFrame([r for r, _ in sres])
    meta = {"rho_base": rho, "lambda_per_bp": LAMBDA0, "hazard_text": "1 orientation-switch event per 5 Mb (prespecified)"}

    # null thresholds per group (base config)
    nul = recs[(recs.kind == "null") & recs.eligible].copy()
    nul["span_bp"] = nul["span_bp"].astype(float)
    nul["switch_rate_per_10Mb"] = nul["viterbi_switches"] / (nul["span_bp"] / 1e7)
    null_q = {}
    for g in GROUPS:
        s = nul[nul.group == g]
        null_q[g] = {"n_windows": int(len(s)), "delta_hat_q95": q95(s.delta_hat), "delta_hat_median": float(s.delta_hat.median()), "delta_hat_max": float(s.delta_hat.max()),
                     "viterbi_switches_q95": q95(s.viterbi_switches), "expected_switches_q95": q95(s.expected_switches),
                     "switch_rate_per_10Mb_q95": q95(s.switch_rate_per_10Mb), "delta_eff_q95": q95(s.delta_eff),
                     "frac_windows_delta_hat_gt0": float((s.delta_hat > 0).mean())}
    # segments
    seg = recs[recs.kind == "segment"].copy()
    sk = kseg.copy()
    sk["name"] = sk.group + "|" + sk.segment.astype(str)
    seg = seg.merge(sk[["name", "chromosome", "bin_start", "bin_end", "cn_call", "continuous_cn"]], on="name", how="left")
    seg["start_bp"] = d.start[seg.bin_start.to_numpy()]
    seg["end_bp"] = d.end[seg.bin_end.to_numpy()]
    seg["segment"] = seg.name.str.split("|").str[1].astype(int)
    seg["null_delta_q95_group"] = seg.group.map(lambda g: null_q[g]["delta_hat_q95"])
    seg["imbalanced"] = seg.eligible & (seg.delta_hat > seg.null_delta_q95_group)
    seg["span_bp"] = seg["end_bp"] - seg["start_bp"]
    seg["switch_rate_per_10Mb"] = seg.viterbi_switches / (seg.span_bp / 1e7)
    seg["recovered_share"] = np.where(seg.delta_eff > seg.delta_naive + 0.005, (seg.delta_hat - seg.delta_naive) / (seg.delta_eff - seg.delta_naive), np.nan)
    seg["switch_rate_gt_null_q95"] = seg.imbalanced & (seg.switch_rate_per_10Mb > seg.group.map(lambda g: null_q[g]["switch_rate_per_10Mb_q95"]))
    seg = seg.sort_values(["group", "segment"]).reset_index(drop=True)
    segcols = ["group", "segment", "name", "chromosome", "start_bp", "end_bp", "n_bins", "cn_call", "continuous_cn", "n_cells", "bins_AT_gt0", "bins_AT_ge5", "AT_sum", "eligible",
               "delta_hat", "h_hat", "h_ci_lo", "h_ci_hi", "llr_vs_balanced", "null_delta_q95_group", "imbalanced", "h_naive", "delta_naive", "delta_naive_ci_lo", "delta_naive_ci_hi",
               "delta_eff", "delta_eff_ci_lo", "delta_eff_ci_hi", "recovered_share", "expected_switches", "viterbi_switches", "switch_rate_per_10Mb", "switch_rate_gt_null_q95",
               "mean_orientation_certainty", "boot_n"]
    seg_out = seg[segcols].rename(columns={"name": "segment_id"})
    seg_out.to_csv(OUT / "segment_haplotypes.tsv", sep="\t", index=False, float_format="%.5g")

    # bin orientation
    imb_flag = {(r.group, r.name): bool(r.imbalanced) for r in seg.itertuples()}
    dh = {(r.group, r.name): r.delta_hat for r in seg.itertuples()}
    parts = []
    for (g, name), bt in bin_tabs.items():
        if bt is None:
            continue
        b = bt.copy()
        b.insert(0, "segment", int(name.split("|")[1]))
        b.insert(0, "group", g)
        b["chromosome"] = d.seq[b.bin.to_numpy()]
        b["start_bp"] = d.start[b.bin.to_numpy()]
        b["p_A_major"] = 1 - b.p_B_major
        b["orientation_call"] = np.where(b.p_B_major >= 0.8, "B", np.where(b.p_B_major <= 0.2, "A", "uncertain"))
        b["delta_hat_segment"] = dh[(g, name)]
        b["segment_imbalanced"] = imb_flag[(g, name)]
        parts.append(b)
    bo = pd.concat(parts, ignore_index=True)
    bo = bo[["group", "segment", "bin", "chromosome", "start_bp", "HB", "AT", "p_B_major", "p_A_major", "orientation_call", "delta_hat_segment", "segment_imbalanced"]]
    bo.to_csv(OUT / "bin_orientation.tsv.gz", sep="\t", index=False, float_format="%.5g", compression="gzip")

    # random null table
    nl = recs[recs.kind == "null"].copy()
    nl["switch_rate_per_10Mb"] = nl["viterbi_switches"] / (nl["span_bp"] / 1e7)
    nl = nl.rename(columns={"name": "window"})
    nl["h_hat"] = 0.5 + nl.delta_hat
    nl[["group", "window", "n_bins", "n_cells", "bins_AT_ge5", "AT_sum", "eligible", "delta_hat", "h_hat", "llr_vs_balanced", "delta_naive", "delta_eff", "expected_switches",
        "viterbi_switches", "switch_rate_per_10Mb"]].to_csv(OUT / "random_null.tsv", sep="\t", index=False, float_format="%.5g")

    # regions
    rg = recs[recs.kind == "region"].copy().rename(columns={"name": "region"})
    fam = {r["name"]: r["family"] for r in regs}
    rg["family"] = rg.region.map(fam)
    rg["recovered_share"] = np.where(rg.delta_eff > rg.delta_naive + 0.005, (rg.delta_hat - rg.delta_naive) / (rg.delta_eff - rg.delta_naive), np.nan)
    rg["hazard_note"] = ""
    rgcols = ["group", "family", "region", "n_bins", "n_cells", "bins_AT_ge5", "AT_sum", "eligible", "h_naive", "delta_naive", "delta_naive_ci_lo", "delta_naive_ci_hi",
              "delta_eff", "delta_eff_ci_lo", "delta_eff_ci_hi", "delta_hat", "h_hat", "h_ci_lo", "h_ci_hi", "llr_vs_balanced", "recovered_share", "expected_switches", "viterbi_switches",
              "mean_orientation_certainty"]
    rg[rgcols].to_csv(OUT / "region_haplotypes.tsv", sep="\t", index=False, float_format="%.5g")

    # cross-cluster
    pair, cc = cross_cluster(seg[seg.group == "all"], bin_tabs, imb_flag, kseg)
    cc.to_csv(OUT / "cluster_consistency.tsv", sep="\t", index=False, float_format="%.4g")

    # per cell
    imb_all = seg[(seg.group == "all") & seg.imbalanced]
    pc, pcstore = percell(imb_all, bin_tabs, lab)
    pc.to_csv(OUT / "percell_support.tsv", sep="\t", index=False, float_format="%.5g")

    # sensitivity
    sn = srecs.copy()
    cmp_rows, ctrl_rows = [], []
    base_all = seg[seg.group == "all"].set_index("name")
    for cn in ["base", "hazard_x0.1", "hazard_x10", "rho_0", "rho_0.02"]:
        if cn == "base":
            ss = base_all.reset_index()
            nn_ = nul[nul.group == "all"]
            rr = rg[rg.group == "all"]
        else:
            x = sn[sn.cfg_name == cn]
            ss = x[x.kind == "segment"].copy()
            nn_ = x[(x.kind == "null") & x.eligible]
            rr = x[x.kind == "region"].rename(columns={"name": "region"})
        qq = q95(nn_.delta_hat)
        ss = ss[ss.eligible].copy()
        ss["span_bp"] = ss["span_bp"].astype(float)
        b = base_all.loc[ss.name]
        ok = ss.name.isin(base_all.index)
        cmp_rows.append({"config": cn, "rho": CONFIGS[cn]["rho"], "lambda_per_Mb": CONFIGS[cn]["lam"] * 1e6, "eligible_all_segments": len(ss), "null_delta_q95": qq,
                         "null_median_delta": float(nn_.delta_hat.median()),
                         "imbalanced_segments": int((ss.delta_hat > qq).sum()),
                         "median_abs_diff_delta_vs_base": float(np.median(np.abs(ss.delta_hat.to_numpy() - b.delta_hat.to_numpy()))),
                         "spearman_delta_vs_base": float(pd.Series(ss.delta_hat.to_numpy()).corr(pd.Series(b.delta_hat.to_numpy()), method="spearman")),
                         "median_viterbi_switches_imbalanced": float(ss[ss.delta_hat > qq].viterbi_switches.median()),
                         "total_viterbi_switches_imbalanced": float(ss[ss.delta_hat > qq].viterbi_switches.sum())})
        for r in rr.itertuples():
            if r.region in ("lit_MYC_8q24", "lit_5q13_qter", "lit_19q13_qter", "lit_9p21_hemizygous_zone_approx"):
                ctrl_rows.append({"config": cn, "region": r.region, "delta_hat": r.delta_hat, "h_hat": 0.5 + r.delta_hat, "viterbi_switches": r.viterbi_switches})
    sens_df = pd.DataFrame(cmp_rows)
    sens_ctrl = pd.DataFrame(ctrl_rows)
    sens_df.to_csv(OUT / "sensitivity.tsv", sep="\t", index=False, float_format="%.5g")
    sens_ctrl.to_csv(OUT / "sensitivity_controls.tsv", sep="\t", index=False, float_format="%.5g")

    ctx = dict(d=d, paths=paths, checks=checks, meta=meta, seg=seg, bo=bo, nul=nul, nl=nl, null_q=null_q, rg=rg, pair=pair, cc=cc, pc=pc, pcstore=pcstore,
               sens_df=sens_df, sens_ctrl=sens_ctrl, reg_tabs=reg_tabs, bin_tabs=bin_tabs, regs=regs, n_rand=n_rand, n_unc=n_unc, nulls=nulls, t0=t0, lab=lab)
    report(ctx)


def report(ctx):
    seg, rg, nl, nul, null_q, pair, cc, pc = ctx["seg"], ctx["rg"], ctx["nl"], ctx["nul"], ctx["null_q"], ctx["pair"], ctx["cc"], ctx["pc"]
    sens_df, sens_ctrl, reg_tabs, bin_tabs = ctx["sens_df"], ctx["sens_ctrl"], ctx["reg_tabs"], ctx["bin_tabs"]
    d, rho = ctx["d"], ctx["meta"]["rho_base"]
    runtime = time.time() - ctx["t0"]
    ra = rg[rg.group == "all"].set_index("region")
    sa = seg[seg.group == "all"]
    imb = sa[sa.imbalanced]
    imb_counts = seg.groupby("group").agg(segments=("name", "size"), eligible=("eligible", "sum"), imbalanced=("imbalanced", "sum"))
    n_sw_gt = int(imb.switch_rate_gt_null_q95.sum())
    qa = null_q["all"]
    ctrl = ["lit_MYC_8q24", "lit_5q13_qter", "lit_19q13_qter", "lit_9p21_hemizygous_zone_approx"]
    cname = {"lit_MYC_8q24": "MYC chr8:127.0-128.5 Mb", "lit_5q13_qter": "5q13-qter (chr5 >= 70 Mb)", "lit_19q13_qter": "19q13-qter (chr19 >= 35 Mb)",
             "lit_9p21_hemizygous_zone_approx": "9p21 zone chr9:20-24 Mb (approx.)"}

    def f3(x):
        return "n/a" if x is None or not np.isfinite(x) else f"{x:.3f}"

    # ---- tables
    rec_tab = []
    for rname in ctrl + [r for r in ra.index if ra.loc[r, "family"] == "saved"]:
        r = ra.loc[rname]
        rec_tab.append({"region": rname, "bins": r.n_bins, "AT_sum": r.AT_sum, "delta_naive": r.delta_naive, "delta_hmm": r.delta_hat,
                        "delta_hmm_ci_lo": r.h_ci_lo - 0.5, "delta_hmm_ci_hi": r.h_ci_hi - 0.5, "delta_phase_agnostic": r.delta_eff,
                        "h_naive": r.h_naive, "h_hat": r.h_hat, "recovered_share": r.recovered_share, "viterbi_switches": r.viterbi_switches,
                        "orientation_certainty": r.mean_orientation_certainty})
    rec_tab = pd.DataFrame(rec_tab)
    seg_tab = imb[["segment", "chromosome", "start_bp", "end_bp", "n_bins", "cn_call", "AT_sum", "h_naive", "h_hat", "h_ci_lo", "h_ci_hi", "delta_eff", "viterbi_switches",
                   "switch_rate_per_10Mb"]].copy()
    seg_tab["start_Mb"] = seg_tab.start_bp / 1e6
    seg_tab["end_Mb"] = seg_tab.end_bp / 1e6
    seg_tab = seg_tab[["segment", "chromosome", "start_Mb", "end_Mb", "n_bins", "cn_call", "AT_sum", "h_naive", "h_hat", "h_ci_lo", "h_ci_hi", "delta_eff", "viterbi_switches", "switch_rate_per_10Mb"]]
    null_tab = pd.DataFrame(null_q).T.reset_index().rename(columns={"index": "group"})
    nonsw = nul[(nul.group == "all")]
    pw = pair.copy()
    ccl = cc.groupby(["cluster", "class"]).size().unstack(fill_value=0).reset_index()
    opp = cc[cc["class"] == "opposite major haplotype"].copy()
    opp["start_Mb"] = opp.start_bp / 1e6
    opp["end_Mb"] = opp.end_bp / 1e6
    opp = opp[["all_segment", "chromosome", "start_Mb", "end_Mb", "cn_call_all", "all_h_hat", "cluster", "bins_compared", "agree_fraction"]]
    mixed = cc[cc["class"] == "mixed"][["all_segment", "chromosome", "all_h_hat", "cluster", "bins_compared", "agree_fraction"]]
    per_seg = cc[cc.bins_compared >= 10].groupby("all_segment")["class"].agg(lambda s: pd.Series({"same": (s == "same major haplotype").sum(), "opposite": (s == "opposite major haplotype").sum(), "mixed": (s == "mixed").sum()}).to_dict())
    n_seg_cmp = len(per_seg)
    n_seg_opp = int(sum(1 for v in per_seg if v["opposite"] > 0))
    n_seg_allsame = int(sum(1 for v in per_seg if v["opposite"] == 0 and v["mixed"] == 0))
    pcs = pc[pc.cluster != "all"]
    pc_view = pc[pc.segment.isin(["all|8", "all|39", "all|76"])][["segment", "chromosome", "cluster", "cells_AT_ge5", "median_AT_per_cell", "h_pooled_oriented", "h_ci_lo", "h_ci_hi",
                                                                  "cell_fraction_median", "cell_fraction_q25", "cell_fraction_q75", "frac_cells_major_below_half", "dispersion_phi",
                                                                  "dispersion_phi_expected_bb", "dispersion_phi_bb_q95"]]
    pc_sum = pcs.groupby("cluster").agg(segments=("segment", "nunique"), median_h_pooled_oriented=("h_pooled_oriented", "median"),
                                        median_frac_cells_below_half=("frac_cells_major_below_half", "median"), median_phi=("dispersion_phi", "median"),
                                        median_phi_expected_bb=("dispersion_phi_expected_bb", "median"),
                                        segments_phi_above_bb_q95=("dispersion_phi", lambda s: int((s.to_numpy() > pcs.loc[s.index, "dispersion_phi_bb_q95"].to_numpy()).sum()))).reset_index()
    reg_pivot = rg.pivot_table(index=["family", "region"], columns="group", values="h_hat").reset_index()
    reg_pivot_naive = rg.pivot_table(index=["family", "region"], columns="group", values="h_naive").reset_index()
    sens_piv = sens_ctrl.pivot_table(index="region", columns="config", values="h_hat").reset_index()

    # ---- figures
    figs = {}
    fig, ax = plt.subplots(4, 3, figsize=(14, 12), gridspec_kw={"width_ratios": [2.2, 2.2, 1]})
    panels = [("segment all|8 (chr1:166.8-236 Mb)", bin_tabs[("all", "all|8")]), ("5q13-qter", reg_tabs[("all", "lit_5q13_qter")]),
              ("19q13-qter", reg_tabs[("all", "lit_19q13_qter")]), ("MYC chr8:127-128.5 Mb", reg_tabs[("all", "lit_MYC_8q24")])]
    hb = np.linspace(0, 1, 41)
    for row, (ttl, bt) in zip(ax, panels):
        bt = bt[bt.AT >= 5]
        x = d.start[bt.bin.to_numpy()] / 1e6
        raw = (bt.HB / bt.AT).to_numpy()
        ori = np.where(bt.p_B_major >= 0.5, raw, 1 - raw)
        row[0].scatter(x, raw, s=6, c=np.where(raw >= 0.5, "#9b2c2c", "#2f5d8a"), alpha=.8)
        row[0].axhline(0.5, color="k", lw=.5); row[0].set_ylim(0, 1)
        row[0].set_title(f"{ttl}: GT-oriented HB/AT per bin (pooled {bt.HB.sum() / bt.AT.sum():.3f})", fontsize=8)
        row[0].set_ylabel("HB / AT")
        row[1].scatter(x, ori, s=6, color="#2f5d8a")
        row[1].axhline(0.5, color="k", lw=.5); row[1].set_ylim(0, 1)
        r2 = row[1].twinx(); r2.plot(x, bt.p_B_major, color="#9b2c2c", lw=.6, alpha=.6); r2.set_ylim(-0.05, 1.05)
        r2.set_ylabel("P(B-major)", color="#9b2c2c", fontsize=7)
        row[1].set_title(f"after HMM orientation: major-haplotype fraction (median {np.median(ori):.3f})", fontsize=8)
        row[2].hist(raw, hb, color="#999", alpha=.8, label="GT-oriented")
        row[2].hist(ori, hb, histtype="step", color="#2f5d8a", lw=1.5, label="HMM-oriented")
        row[2].axvline(0.5, color="k", lw=.5); row[2].set_title("distribution", fontsize=8)
        for a_ in row[:2]:
            a_.set_xlabel("Mb")
    ax[0, 2].legend(fontsize=7)
    fig.tight_layout()
    figs["orient"] = C.fig_to_b64(fig)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 3.8))
    names = rec_tab.region.tolist()
    y = np.arange(len(names))
    ax.barh(y - 0.27, rec_tab.delta_naive, 0.27, color="#999", label="naive GT-oriented |HB/AT - 0.5|")
    ax.barh(y, rec_tab.delta_hmm, 0.27, xerr=[np.clip(rec_tab.delta_hmm - rec_tab.delta_hmm_ci_lo, 0, None), np.clip(rec_tab.delta_hmm_ci_hi - rec_tab.delta_hmm, 0, None)], color="#2f5d8a", label="HMM-oriented delta_hat (cell-bootstrap 95% CI)")
    ax.barh(y + 0.27, rec_tab.delta_phase_agnostic, 0.27, color="#8a5a00", label="phase-agnostic |dBAF| (exp 06 formula)")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("|haplotype fraction - 0.5|")
    ax.legend(fontsize=7)
    fig.tight_layout()
    figs["recover"] = C.fig_to_b64(fig)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    bins_ = np.linspace(0, 0.3, 31)
    ax[0].hist(nonsw.delta_hat, bins=bins_, color="#bbb", label=f"{len(nonsw)} uncalled random 25-bin windows")
    ax[0].axvline(qa["delta_hat_q95"], color="#9b2c2c", ls="--", label=f"null 95th pct {qa['delta_hat_q95']:.3f}")
    ax[0].set_xlabel("delta_hat"); ax[0].legend(fontsize=7); ax[0].set_title("Null (all cells)", fontsize=9)
    ax[1].hist(sa.delta_hat, bins=bins_, color="#2f5d8a", label=f"{len(sa)} saved `all` segments")
    ax[1].axvline(qa["delta_hat_q95"], color="#9b2c2c", ls="--")
    ax[1].set_xlabel("delta_hat"); ax[1].legend(fontsize=7); ax[1].set_title("Saved segments", fontsize=9)
    fig.tight_layout()
    figs["null"] = C.fig_to_b64(fig)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    M = np.full((6, 6), np.nan)
    for r in pair.itertuples():
        a_, b_ = GROUPS.index(r.group_a), GROUPS.index(r.group_b)
        M[a_, b_] = M[b_, a_] = r.agree_fraction
    im = ax[0].imshow(M, vmin=0, vmax=1, cmap="RdBu")
    ax[0].set_xticks(range(6)); ax[0].set_xticklabels([g.replace("cluster", "c") for g in GROUPS], fontsize=8)
    ax[0].set_yticks(range(6)); ax[0].set_yticklabels([g.replace("cluster", "c") for g in GROUPS], fontsize=8)
    for a_ in range(6):
        for b_ in range(6):
            if np.isfinite(M[a_, b_]):
                ax[0].text(b_, a_, f"{M[a_, b_]:.2f}", ha="center", va="center", fontsize=7)
    ax[0].set_title("Fraction of confident bins with the same major haplotype", fontsize=9)
    fig.colorbar(im, ax=ax[0], fraction=.046)
    ccl_ = cc.groupby(["cluster", "class"]).size().unstack(fill_value=0)
    ccl_.plot.bar(stacked=True, ax=ax[1], color=["#9b2c2c", "#bbb", "#d9a441", "#2f5d8a"][:ccl_.shape[1]])
    ax[1].set_title("Imbalanced `all` segments: class per cluster", fontsize=9)
    ax[1].legend(fontsize=6)
    fig.tight_layout()
    figs["xcl"] = C.fig_to_b64(fig)
    plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    for a, sg in zip(ax, ["all|8", "all|39", "all|76"]):
        data = [ctx_fr for ctx_fr in (ctx["pcstore"].get((sg, c)) for c in ["0", "1", "2", "3", "4"])]
        a.boxplot([x if x is not None else [np.nan] for x in data], showfliers=False)
        a.axhline(0.5, color="k", lw=.5)
        a.set_xticklabels(["c0", "c1", "c2", "c3", "c4"])
        t_ = seg[(seg.group == "all") & (seg.name == sg)].iloc[0]
        a.set_title(f"{sg}: {t_.chromosome}:{t_.start_bp/1e6:.0f}-{t_.end_bp/1e6:.0f} Mb", fontsize=9)
    ax[0].set_ylabel("per-cell major-haplotype fraction (cells with AT >= 5)")
    fig.tight_layout()
    figs["cell"] = C.fig_to_b64(fig)
    plt.close(fig)

    # ---- numbers for text
    r5, rM, r19, r9 = ra.loc["lit_5q13_qter"], ra.loc["lit_MYC_8q24"], ra.loc["lit_19q13_qter"], ra.loc["lit_9p21_hemizygous_zone_approx"]
    chr7 = ra.loc[[i for i in ra.index if i.startswith("chr7_")]]
    sc = sens_df.set_index("config")
    n_opp = int((cc["class"] == "opposite major haplotype").sum())
    n_same = int((cc["class"] == "same major haplotype").sum())
    n_mix = int((cc["class"] == "mixed").sum())
    n_ins = int(cc["class"].str.startswith("insufficient").sum())
    c4 = pc[(pc.segment == "all|76") & (pc.cluster == "4")].iloc[0]
    med_phi, med_phi_bb = float(pcs.dispersion_phi.median()), float(pcs.dispersion_phi_expected_bb.median())
    n_phi_hi = int((pcs.dispersion_phi.to_numpy() > pcs.dispersion_phi_bb_q95.to_numpy()).sum())
    sw5 = r5.viterbi_switches
    n_small = int((imb.delta_hat < 0.06).sum())
    span5 = float(d.end[np.flatnonzero((d.seq == "chr5") & (d.start >= 70e6))[-1]] - 70e6) / 1e6

    cl_imb = ", ".join(g[-1] + ": " + str(int(imb_counts.loc[g, "imbalanced"])) for g in GROUPS[1:])
    chr7_h = ", ".join(f"{ra.loc[i, 'h_hat']:.3f}" for i in chr7.index)
    verdict = (f"Exploratory: inferring the phase orientation of each bin with an HMM recovers most of the allelic imbalance that GT-oriented summing loses at "
               f"5q13-qter (|h-0.5| {r5.delta_naive:.3f} naive to {r5.delta_hat:.3f} oriented; phase-agnostic {r5.delta_eff:.3f}) and at MYC ({rM.delta_naive:.3f} to {rM.delta_hat:.3f}; "
               f"phase-agnostic {rM.delta_eff:.3f}), little at 19q13-qter ({r19.delta_naive:.3f} to {r19.delta_hat:.3f}; phase-agnostic {r19.delta_eff:.3f}) and none at the 9p21 zone. "
               f"Major/minor labels are inferred, not verified phase.")
    kf = [
        f"<b>Orientation inference recovers the attenuated imbalance at strongly imbalanced loci.</b> 5q13-qter: naive h {r5.h_naive:.3f} vs HMM h_hat {r5.h_hat:.3f} (95% CI {r5.h_ci_lo:.3f}-{r5.h_ci_hi:.3f}), phase-agnostic delta {r5.delta_eff:.3f} (recovered share {r5.recovered_share:.2f}); MYC: {rM.h_naive:.3f} vs {rM.h_hat:.3f} (CI {rM.h_ci_lo:.3f}-{rM.h_ci_hi:.3f}), share {rM.recovered_share:.2f}.",
        f"<b>Weak or sparse loci are not rescued.</b> 19q13-qter: {r19.h_naive:.3f} vs {r19.h_hat:.3f} (phase-agnostic delta {r19.delta_eff:.3f}, share {r19.recovered_share:.2f}; hazard x0.1 gives {sens_ctrl[(sens_ctrl.config=='hazard_x0.1')&(sens_ctrl.region=='lit_19q13_qter')].h_hat.iloc[0]:.3f}, rho=0 gives {sens_ctrl[(sens_ctrl.config=='rho_0')&(sens_ctrl.region=='lit_19q13_qter')].h_hat.iloc[0]:.3f}). The 9p21 zone: {r9.h_naive:.3f} vs {r9.h_hat:.3f} although the phase-agnostic delta is {r9.delta_eff:.3f} (CI {r9.delta_eff_ci_lo:.3f}-{r9.delta_eff_ci_hi:.3f}); {int(r9.bins_AT_ge5)} bins with AT>=5 and 1,486 AT in total.",
        f"<b>GT orientation flips often between bins.</b> 5q13-qter needs {int(sw5)} Viterbi orientation switches over about {span5:.0f} Mb ({sw5/span5:.2f}/Mb) against a prior of {LAMBDA0*1e6:.2f}/Mb; {n_sw_gt} of {len(imb)} imbalanced `all` segments exceed the null 95th percentile of switch rate ({qa['switch_rate_per_10Mb_q95']:.1f} per 10 Mb). Switch counts depend on the prespecified hazard (5q: {sens_ctrl[(sens_ctrl.config=='hazard_x0.1')&(sens_ctrl.region=='lit_5q13_qter')].viterbi_switches.iloc[0]:.0f} at x0.1, {sens_ctrl[(sens_ctrl.config=='hazard_x10')&(sens_ctrl.region=='lit_5q13_qter')].viterbi_switches.iloc[0]:.0f} at x10).",
        f"<b>Null calibration.</b> In {qa['n_windows']} uncalled random 25-bin windows (all cells) delta_hat has median {qa['delta_hat_median']:.3f}, 95th percentile {qa['delta_hat_q95']:.3f} (max {qa['delta_hat_max']:.3f}); {int(imb_counts.loc['all','imbalanced'])} of {int(imb_counts.loc['all','eligible'])} eligible `all` segments exceed it (clusters: {cl_imb}). Null windows are 25 bins, longer segments are less noisy, so the threshold is conservative for them.",
        f"<b>Cross-cluster orientation is mostly the same haplotype.</b> Across {n_same + n_opp + n_mix} (imbalanced all-segment x cluster) comparisons with >=10 confident bins, {n_same} are the same major haplotype, {n_mix} mixed and {n_opp} opposite; the {n_opp} opposite cases are all cluster 4 on chr8:115-145 Mb (agreement 0.06-0.15), a candidate mirrored event in 182 cells (hypothesis). {n_ins} comparisons had too few bins.",
        f"<b>Per-cell fractions are not binomially homogeneous.</b> Across imbalanced segments the median dispersion phi is {med_phi:.2f} vs {med_phi_bb:.2f} expected for beta-binomial cells with rho = {rho:.4f}; {n_phi_hi}/{len(pcs)} segment-cluster rows exceed the simulated 95th percentile. At chr8:115-134 Mb {c4.frac_cells_major_below_half:.0%} of cluster-4 cells have major fraction < 0.5 (oriented by the pooled call) vs about 17-20% in clusters 0-3.",
        f"<b>Saved chr7 windows.</b> chr7:8.3-10.7 and 70.0-73.3 Mb show a small oriented imbalance (delta_hat {ra.loc['chr7_8.3_10.7','delta_hat']:.3f}, {ra.loc['chr7_70.0_73.3','delta_hat']:.3f}; null 95th pct {qa['delta_hat_q95']:.3f}) though their phase-agnostic delta is 0.000; 51.4-53.5 and 67.4-70.0 Mb are 0.000. This is weak, near the null edge, and far below the 5q/MYC values.",
        f"<b>Sensitivity.</b> Number of imbalanced `all` segments (each setting with its own null): base {int(sc.loc['base','imbalanced_segments'])}, hazard x0.1 {int(sc.loc['hazard_x0.1','imbalanced_segments'])}, x10 {int(sc.loc['hazard_x10','imbalanced_segments'])}, rho=0 {int(sc.loc['rho_0','imbalanced_segments'])}, rho=0.02 {int(sc.loc['rho_0.02','imbalanced_segments'])}; Spearman of delta_hat vs base {sc.drop('base').spearman_delta_vs_base.min():.2f}-{sc.drop('base').spearman_delta_vs_base.max():.2f}. Large effects (5q, MYC, chr1q) are stable; small effects are not.",
    ]

    glossary = {
        "group": "cell set: all (3,460 cells) or WNN cluster 0-4 (cluster0..4 is that cluster's own saved segmentation).",
        "segment / segment_id / all_segment": "segment index in karyogram_cluster_segments.tsv; segment_id = group|index.",
        "chromosome, start_bp, end_bp": "genomic extent of the segment (absolute bins 100 kb; bp from AnnData var).",
        "n_bins, n_cells": "bins in the segment; cells pooled for the estimate.",
        "cn_call": "pyEpi call of the saved segment (0, 0.5, 1, 1.5, 2); a call, not truth.",
        "bins_AT_gt0 / bins_AT_ge5": "bins with any / at least 5 pooled site observations; eligibility needs bins_AT_ge5 >= 8.",
        "AT_sum": "summed SNP-site DP over selected phased heterozygous SNPs (not unique fragments).",
        "eligible": "enough allele data to fit (bins_AT_ge5 >= 8).",
        "delta_hat, h_hat": "profile-likelihood estimate of the major-haplotype excess delta (grid 0-0.49, step 0.0025) under the orientation HMM; h_hat = 0.5 + delta_hat is the inferred major-haplotype fraction. Boundary 0 means no imbalance fitted.",
        "h_ci_lo, h_ci_hi": "percentile cell-bootstrap 95% interval of h_hat (200 resamples per segment/region with 0.01 grid and parabolic refinement; boot_n); seeds (20261104,16,stream,group,item). Lower end is bounded at 0.5 so the interval is optimistic near zero. The bootstrap uses a coarser grid than the point estimate and can sit slightly above it (5q13-qter: point 0.600, interval 0.603-0.607), so the interval is not guaranteed to contain the point estimate.",
        "llr_vs_balanced": "2 (log-lik at delta_hat - log-lik at delta = 0), same HMM.",
        "null_delta_q95_group": "95th percentile of delta_hat over uncalled random 25-bin windows for the same cell set.",
        "imbalanced": "eligible and delta_hat > null_delta_q95_group (a screening rule, not a test).",
        "h_naive, delta_naive": "GT-oriented estimate: f = sum HB / sum AT; h_naive = max(f, 1 - f); delta_naive = |f - 0.5|. The max makes it a lower-bounded fold; the sign is the GT label only.",
        "delta_naive_ci_lo/hi, delta_eff_ci_lo/hi": "cell-bootstrap percentile intervals of the naive and phase-agnostic quantities.",
        "delta_eff": "phase-agnostic effective |dBAF| of exp 06 (sqrt(max(0, mean[((HB-HA)^2/AT - (1+(AT-1) rho)) / (4 AT)])) over bins with AT >= 5).",
        "recovered_share": "(delta_hat - delta_naive) / (delta_eff - delta_naive); only defined when delta_eff exceeds delta_naive by > 0.005. Can exceed 1 or be negative; not a probability.",
        "expected_switches": "posterior expected number of orientation switches between consecutive informative bins at delta_hat (NaN when delta_hat = 0).",
        "viterbi_switches": "switches on the most probable orientation path at delta_hat (NaN when delta_hat = 0).",
        "switch_rate_per_10Mb": "viterbi_switches per 10 Mb of segment span.",
        "switch_rate_gt_null_q95": "imbalanced segment whose switch rate exceeds the null 95th percentile of switch rate.",
        "mean_orientation_certainty": "mean over informative bins of max(P(B-major), P(A-major)).",
        "bin, HB, AT, p_B_major, p_A_major": "(bin_orientation.tsv.gz) absolute bin index; summed HB and AT of the group; posterior probability that haplotype B is major in the bin; 1 - that.",
        "orientation_call": "B if p_B_major >= 0.8, A if <= 0.2, else uncertain.",
        "delta_hat_segment, segment_imbalanced": "segment-level values attached to each bin.",
        "window, switch counts (random_null.tsv)": "random 25-bin window name from exp 04 build_regions (family random, no called bin) and the same statistics as segments.",
        "bins_compared, agree_fraction, class (cluster_consistency.tsv)": "bins of an imbalanced `all` segment where `all` and that cluster's own segment are both imbalanced and confident (P >= 0.8 or <= 0.2); fraction with the same major haplotype; class: same >= 0.8, opposite <= 0.2, mixed otherwise, insufficient < 10 bins.",
        "cluster_imbalanced_bins_in_segment": "bins of the `all` segment that fall in an imbalanced own-segment of that cluster.",
        "cells_AT_ge5, h_pooled_oriented, cell_fraction_*, frac_cells_major_below_half": "(percell_support.tsv) cells with >= 5 site observations in the segment; pooled oriented major fraction; per-cell quantiles; share of cells below 0.5. Orientation comes from the pooled `all` posterior, so values are biased towards the major side.",
        "dispersion_phi, dispersion_phi_expected_bb, dispersion_phi_bb_q95": "Pearson dispersion of per-cell major counts around the pooled fraction; mean and 95th percentile under 100 simulated beta-binomial cell sets with the same AT per cell and rho (approximate, since rho is bin-level).",
        "config, null_delta_q95, imbalanced_segments, spearman_delta_vs_base": "(sensitivity.tsv) setting; recomputed null 95th percentile; count of imbalanced `all` segments under it; rank correlation of delta_hat with the base run.",
        "family": "region family from exp 04: saved window, literature locus.",
    }

    method_pre = """Observed per segment (cells pooled): for each informative bin j (AT_j > 0): k_j = HB_j, n_j = AT_j = HA_j + HB_j.
Hidden orientation o_j in {B, A}; B-major means haplotype B is the major haplotype in bin j.
Emission:  k_j | n_j, o_j = B ~ BetaBinomial(n_j, p = 0.5 + delta, rho)      o_j = A:  p = 0.5 - delta
           BetaBinomial(k; n, a, b) with a = p c, b = (1-p) c, c = (1 - rho)/rho,  rho = 0.0038 (exp 06; rho = 0 means Binomial)
Transition between consecutive informative bins at base-pair distance d_j (bin start to bin start):
           s_j = P(o_{j+1} != o_j) = 0.5 (1 - exp(-2 * lambda * d_j)),   lambda = 1 / 5 Mb  (prespecified)
           sensitivity: lambda x 0.1 and x 10.   Initial state 1/2, 1/2.
Likelihood: L(delta) = sum over orientation paths of  prod_j P(k_j | n_j, o_j; delta) prod_j T(o_j, o_{j+1}; d_j)   (scaled forward algorithm)
Estimate:   delta_hat = argmax_{delta in [0, 0.49]} L(delta)  (profile likelihood on a grid, step 0.0025); h_hat = 0.5 + delta_hat
            no EM is needed: delta is one parameter and the grid is exhaustive.
Posterior:  P(o_j = B | data; delta_hat) by forward-backward; expected switches = sum_j P(o_j != o_{j+1} | data);
            Viterbi path for the switch count.
Cell bootstrap: resample cells with replacement (multinomial weights, seed (20261104,16,stream,group,item)), recompute k_j, n_j, refit
            delta on the 0.01 grid with parabolic refinement; 95% percentile interval.
Naive:      f = sum HB / sum AT;  h_naive = max(f, 1-f).
Phase-agnostic (exp 06):  delta_eff = sqrt(max(0, mean_j [ ((2k_j - n_j)^2 / n_j - (1 + (n_j - 1) rho)) / (4 n_j) ]))  over bins with n_j >= 5.
Imbalanced:  delta_hat > 95th percentile of delta_hat over 60 uncalled random 25-bin windows (same cell set).
Orientation labels are anchored to the HB column (GT phase), so P(B-major) can be compared between cell groups on the same bin."""

    sections = [
        ("Question and why it matters",
         "<p>Allele-specific relative copy number needs, for each region with allelic imbalance, which haplotype carries the extra (or the retained) copy. "
         "Experiment 06 showed that summing HB over bins in GT orientation attenuates imbalance strongly (the published MYC amplification gives an effective |dBAF| of about 0.106 phase-agnostically, but the GT-oriented sum shows far less) "
         "because the GT phase between bins is not verified and the package has no phase-set IDs. Numbat (Gao et al. 2023) handles this with a haplotype-orientation HMM; this experiment applies the same idea at 100 kb-bin resolution to these data "
         "and asks (i) how well the orientation and the major-haplotype fraction are recovered, (ii) how much of the attenuation at the published MCF10-A gains is recovered, "
         "(iii) whether the same haplotype is major in every cluster, and (iv) what the null distribution looks like. It is an exploratory diagnostic, not a validated phasing method.</p>"),
        ("Inputs",
         "<ul><li>HA / HB / AT parquet (bins x barcodes) joined by barcode name to the AnnData var / obs; checks: " + json.dumps(ctx["checks"]) + ".</li>"
         "<li>Saved segmentations: <code>karyogram_cluster_segments.tsv</code>, groups all and cluster0-4 (955 segments).</li>"
         "<li>WNN clusters via <code>common.cell_clusters</code> (sizes 1286 / 833 / 603 / 556 / 182).</li>"
         f"<li>Regions: the seven <code>part == 'target'</code> windows of <code>research/2026-10-29-low-coverage/shared/regions.tsv</code> (six candidate windows plus the historical chr9 70-100 Mb comparison, as in exp 04) and the four literature loci of exp 04.</li>"
         f"<li>Null: first {len(ctx['nulls'])} of the {ctx['n_unc']} random windows (of {ctx['n_rand']}, exp 04 rule, seed (20261104, 4)) with no called bin in the `all` segmentation.</li>"
         f"<li>rho = {rho:.5f} from exp 06 summary.json; hazard lambda = 1 switch / 5 Mb (prespecified, not fitted).</li></ul>"),
        ("Method", "<pre>" + method_pre + "</pre>"
         "<p>Segments need bins_AT_ge5 >= 8 to be fitted; the others are reported as not eligible. The same pipeline runs on regions (all cells and each cluster) and on random windows. "
         "The hazard is a prior chosen before looking at results; because the observed switch counts are far above it, the x10 sensitivity is informative.</p>"),
        ("Results",
         "<p>All numbers below come from the saved tables in this folder. Haplotype labels (A/B major) are inferred from allele counts.</p>"),
        ("Orientation and switches",
         C.img(figs["orient"], "orientation tracks") + C.caption("How to read: left, the GT-oriented HB/AT of each bin (AT >= 5), coloured by side of 0.5; middle, the same bins after HMM orientation (fraction on the inferred major haplotype), with the posterior P(B-major) in red; right, the two distributions. In imbalanced segments the GT-oriented values are bimodal and mirror-symmetric around 0.5 (for chr1:166.8-236 Mb about half the bins near 0.29 and half near 0.71), in runs of a few bins: the GT phase label switches between runs, so the pooled GT sum cancels towards 0.5. The HMM flips the runs where B is the minor allele, which moves every bin to the same side. That is only meaningful where the per-bin deviation is far above sampling noise and mirror-symmetric; in uncalled random windows the same procedure gives |h - 0.5| <= 0.035 (null 95th percentile), which bounds the folding bias.") +
         "<p>Imbalanced `all` segments:</p>" + C.html_table(seg_tab, floatfmt="{:.3f}") +
         C.caption(f"The {len(imb)} `all` segments with delta_hat above the null 95th percentile ({qa['delta_hat_q95']:.3f}). h_naive = GT-oriented fold; h_hat = HMM; delta_eff = phase-agnostic. "
                   "chr1q (cn_call 2) and chr5q/chr8q segments dominate; many small cn_call 0 or 0.5 segments are near the threshold."), 3),
        ("Recovered imbalance at the controls",
         C.html_table(rec_tab, floatfmt="{:.3f}") + C.img(figs["recover"], "recovered imbalance") +
         C.caption("|h - 0.5| for the literature loci and saved windows, all cells. Grey: GT-oriented sum. Blue: HMM-oriented (error bars: cell-bootstrap 95% CI). Brown: phase-agnostic exp 06 estimate. "
                   "recovered_share = (HMM - naive) / (phase-agnostic - naive) where defined.") +
         f"<p>A one-copy gain on one haplotype in a clonal diploid population would give h = 2/3 and a clonal 2+1 to 3+1 pair would give 0.67-0.75; the observed h_hat of {r5.h_hat:.2f} (5q) and {rM.h_hat:.2f} (MYC) are lower, which is compatible with a subclonal fraction, "
         "attenuation by repeated site counting or mapping, or a different copy-number state. This is a reference point, not a validation.</p>", 3),
        ("Cross-cluster consistency",
         C.img(figs["xcl"], "cross-cluster") + C.caption("Left: for every pair of cell groups, the fraction of bins (both imbalanced, both confident) with the same major haplotype; 0.5 would be unrelated orientation. Right: classification of each imbalanced `all` segment per cluster.") +
         "<h4>Pairwise (bins)</h4>" + C.html_table(pw, floatfmt="{:.3f}") + "<h4>Class counts per cluster</h4>" + C.html_table(ccl) +
         f"<p>{n_seg_cmp} imbalanced `all` segments have at least one cluster with >= 10 comparable bins; {n_seg_allsame} have only same-haplotype clusters, {n_seg_opp} have at least one cluster with the opposite major haplotype.</p>"
         "<h4>Opposite major haplotype (candidate mirrored events; hypotheses)</h4>" + C.html_table(opp, floatfmt="{:.3f}") +
         "<h4>Mixed</h4>" + C.html_table(mixed, floatfmt="{:.3f}") +
         C.caption("An opposite call means that in the bins where both groups are confident, the cluster's major haplotype is the `all` call's minor haplotype. It is consistent with independent events on different parental haplotypes in that cluster, but also with cluster-specific allelic accessibility; it does not show either."), 3),
        ("Null calibration",
         C.img(figs["null"], "null") + C.html_table(null_tab, floatfmt="{:.3f}") +
         C.caption("Uncalled random 25-bin windows run through the same pipeline: delta_hat distribution (all cells) and the same quantiles for each cell group. "
                   "viterbi_switches_q95 and switch_rate_per_10Mb_q95 are 95th percentiles over windows.") +
         f"<p>Imbalanced segments per group (delta_hat above that group's null 95th percentile):</p>" + C.html_table(imb_counts.reset_index()) +
         f"<p>With 60 windows the 95th percentile rests on three windows above it; its sampling error is large. Because delta_hat is bounded at 0 it is biased upwards in noise ({qa['frac_windows_delta_hat_gt0']:.0%} of null windows have delta_hat > 0). "
         f"{n_small} of the {len(imb)} imbalanced `all` segments have delta_hat below 0.06, i.e. within about 1.7x the null 95th percentile; they are the least secure calls.</p>", 3),
        ("Per-cell support",
         C.img(figs["cell"], "per cell") + C.caption("Per-cell major-haplotype fraction (cells with AT >= 5 in the segment), orientation from the pooled `all` posterior, by cluster for three imbalanced segments. Cluster 4 on chr8:115-134 Mb falls below 0.5.") +
         C.html_table(pc_view, floatfmt="{:.3f}") + "<h4>Summary over all imbalanced `all` segments</h4>" + C.html_table(pc_sum, floatfmt="{:.3f}") +
         C.caption("dispersion_phi near 1 means cell fractions vary as expected from sampling at a common fraction. Values above the simulated 95th percentile mean extra cell-to-cell or within-cell variation (subclonal imbalance, repeated site counting, or the approximation of the simulation). "
                   "Orientation is taken from the pooled call, which inflates the major fraction."), 3),
        ("Per-region detail",
         "<h4>Haplotype fraction h_hat (HMM) by group</h4>" + C.html_table(reg_pivot, floatfmt="{:.3f}") + "<h4>Naive h by group</h4>" + C.html_table(reg_pivot_naive, floatfmt="{:.3f}") +
         "<h4>All cells, full detail</h4>" + C.html_table(rg[rg.group == "all"].drop(columns=["hazard_note"]), floatfmt="{:.3f}") +
         C.caption("Saved windows (six candidates plus the historical chr9 70-100 Mb comparison) and four literature loci, per cell group. Small regions (< 50 AT-bearing bins) have wide intervals; values of 0.500 mean the fit reached the boundary.") +
         f"<p>chr7 saved windows: h_hat {chr7_h} (windows {', '.join(chr7.index)}). "
         "The 9p21 zone and chr9 saved windows have phase-agnostic delta of 0.13-0.19 in few bins, but no consistent oriented imbalance is fitted; imbalance concentrated in a few bins does not fit a constant-delta model.</p>"),
        ("Sensitivity",
         C.html_table(sens_df, floatfmt="{:.4f}") + "<h4>Control loci h_hat</h4>" + C.html_table(sens_piv, floatfmt="{:.3f}") +
         C.caption("Each setting recomputes the null (60 windows) and the imbalanced count on the 164 `all` segments. hazard x0.1 / x10 change lambda; rho 0 and 0.02 change the beta-binomial overdispersion (base 0.0038)."), 3),
        ("Limits and alternative explanations",
         "<ul><li><b>Allele-specific accessibility / mapping bias.</b> HB / HA are chromatin accessibility reads; allele-specific accessibility (imprinting, regulatory variants) or reference mapping bias gives a haplotype excess without copy-number change. A consistent orientation would look the same.</li>"
         "<li><b>Copy-number-neutral LOH and loss.</b> Imbalance does not say gain or loss; depth is needed for that. In loss or LOH regions heterozygous sites are poorly ascertained (AT per fragment is lowest at 9p21 in exp 06), so the allele fraction there is a biased sample.</li>"
         "<li><b>No SNP-level phase sets.</b> Phase is the sum over selected SNPs per bin; the HMM works on bin-level orientation. Observed switch rates (about 1 per Mb at 5q) are much higher than the prior; a model with SNP-level phasing would be more accurate.</li>"
         "<li><b>Repeated site counting.</b> AT sums SNP-site DP, not unique fragments; the 0.0038 rho only partly absorbs it. h_hat could be attenuated or intervals too narrow; cell-bootstrap CIs do not include bin-level model error.</li>"
         "<li><b>Constant-delta assumption.</b> One delta per segment; heterogeneous imbalance inside a segment, subclones, or segments mixing different events are averaged.</li>"
         "<li><b>Boundary and threshold.</b> delta_hat >= 0 and the 25-bin null make small-segment calls near the threshold uncertain; the null has only 60 windows.</li>"
         "<li><b>Orientation labels are inferred.</b> Major/minor is a statement about count asymmetry; haplotype identity across segments is not verified (GT phase, no phase sets).</li>"
         "<li>Per-cell orientation uses the pooled call and so is biased towards the major side.</li></ul>"),
    ]
    supports = [
        "Orientation inference raises the estimated imbalance at 5q13-qter and MYC by roughly a factor 10 and 3 relative to the GT-oriented sum, close to the phase-agnostic estimate.",
        "Large imbalanced segments (chr1q, chr5q, chr8q) are oriented the same way in all clusters with enough data.",
        "A cluster-4-specific opposite orientation at chr8:115-145 Mb is a candidate for a mirrored (hypothesis) or cluster-specific allelic event.",
        "GT bin orientation is not stable over Mb distances in this data (switch counts far above the prior).",
    ]
    not_supports = [
        "Verified haplotype identity, phase blocks, or allele-specific copy number; depth was not used.",
        "Any statement that the imbalance is a copy-number change (allele-specific accessibility and mapping bias cannot be excluded).",
        "Recovery at weakly imbalanced or sparse loci: 19q13-qter is only partly recovered and 9p21 is not, and the chr7 window values are within about 1.3x the null threshold.",
        "A calibrated significance test: the null has 60 windows of 25 bins and delta_hat is boundary-biased.",
        "Mirrored events as fact: cluster 4 is 182 cells and the opposite call involves three adjacent chr8 segments.",
    ]
    decision = ("<p>recommendation.md (rank 5, joint depth-allele likelihood) assumed that allelic evidence is available as a signed BAF. In these data the GT-oriented sum is attenuated "
                "(5q13-qter naive |h-0.5| " + f"{r5.delta_naive:.3f}" + "), so a joint method needs an explicit latent phase-orientation state, as in Numbat, before depth and allele evidence are combined. "
                "This result supports adding that state to the evidence layer (and to experiment 13's joint HMM) but not adopting any particular major/minor call as a CN statement. "
                "Cluster-specific orientation matters: pooled `all` orientation hides opposite haplotypes in a subcluster (chr8q, cluster 4), so allele-specific CN should be estimated per cluster and compared, not only for `all`. "
                "No scientific behaviour of the existing pipeline is changed by this experiment.</p>")
    nxt = [
        "Run the orientation HMM jointly with the depth track (experiment 13) so that depth can break ties between gain and loss at imbalanced segments.",
        "Replace bin-level HB by SNP-level counts with local phase blocks (or read-backed phasing) and compare switch counts with the prior.",
        "Validate on a sample with known phase (or a matched normal / bulk WGS BAF) before calling any haplotype major in a published locus.",
        "Test the cluster-4 chr8q opposite orientation by splitting cluster 4 and by an independent pseudobulk of allele counts; check allele-specific accessibility at heterozygous sites in that interval.",
        "Replace the constant-delta assumption by a per-bin state (e.g. 0, delta1, delta2) and use a larger, size-matched null (>= 300 windows per segment length class).",
    ]
    prov = {"command": "/home/user/envs/cnv310/bin/python -I research/2026-11-04-cloud-method-review/experiments/16_major_minor_haplotypes/run.py",
            "seed_tuple": list(SEED16), "streams": {"segments": 1, "regions": 2, "null": 3, "percell": 5}, "runtime_seconds": runtime,
            "hazard_per_bp": LAMBDA0, "rho": rho, "grid_fine": [0, 0.49, 0.0025], "grid_coarse": [0, 0.49, 0.01], "B_segments": B_SEG, "B_regions": B_REG, "null_windows": len(ctx["nulls"]),
            "environment": C.environment()}
    C.write_html_report(OUT / "report.html", "Major and minor haplotype inference",
                        "Experiment 16 · phase-orientation HMM for major / minor haplotypes in imbalanced regions (exploratory)", sections, prov,
                        key_findings=kf, verdict=verdict, glossary=glossary, supports=supports, not_supports=not_supports, decision_impact=decision, next_steps=nxt)
    summ = {"experiment": "16_major_minor_haplotypes", "status": "exploratory diagnostic; major/minor labels inferred, not verified phase", "seeds": {"tuple": list(SEED16), "streams": prov["streams"]},
            "hazard": {"lambda_per_bp": LAMBDA0, "text": "1 orientation switch per 5 Mb; s = 0.5(1-exp(-2 lambda d)); prespecified; sensitivity x0.1, x10"},
            "parameters": {"rho": rho, "delta_grid_fine": prov["grid_fine"], "delta_grid_bootstrap": prov["grid_coarse"], "min_bins_AT_ge5": MIN_BINS_AT5, "bootstrap_segments": B_SEG, "bootstrap_regions": B_REG,
                           "null_windows": len(ctx["nulls"]), "random_windows_total": ctx["n_rand"], "random_windows_uncalled": ctx["n_unc"]},
            "alignment_checks": ctx["checks"], "null_by_group": null_q, "imbalanced_by_group": imb_counts.reset_index().to_dict("records"),
            "controls_all_cells": rec_tab.to_dict("records"), "cross_cluster_pairs": pair.to_dict("records"),
            "cross_cluster_counts": {"same": n_same, "opposite": n_opp, "mixed": n_mix, "insufficient": n_ins},
            "percell_median_phi": med_phi, "percell_median_phi_expected": med_phi_bb, "percell_rows_phi_above_bb_q95": n_phi_hi, "percell_rows": len(pcs),
            "sensitivity": sens_df.to_dict("records"), "sensitivity_controls": sens_ctrl.to_dict("records"),
            "imbalanced_all_segments_switch_rate_above_null_q95": n_sw_gt, "runtime_seconds": runtime, "environment": C.environment(),
            "files": ["run.py", "segment_haplotypes.tsv", "bin_orientation.tsv.gz", "cluster_consistency.tsv", "random_null.tsv", "region_haplotypes.tsv", "percell_support.tsv", "sensitivity.tsv", "sensitivity_controls.tsv", "report.html"]}
    C.write_json(OUT / "summary.json", summ)
    print(json.dumps({k: summ[k] for k in ("null_by_group", "cross_cluster_counts", "percell_median_phi", "percell_rows_phi_above_bb_q95", "percell_rows")}, indent=1, default=str)[:3000])


if __name__ == "__main__":
    main()
