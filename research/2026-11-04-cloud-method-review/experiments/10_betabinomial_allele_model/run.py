"""Experiment 10: phase-aware beta-binomial allele model (method record E:2), exploratory diagnostic.

Data: HA/HB/AT parquet (bins x barcodes; AT = summed per-SNP-site DP, NOT unique fragments; phase orientation across
bins is GT-based and unverified; no phase-set IDs). Joined to AnnData by barcode name. Per region and per group
(pooled over all cells, and each WNN cluster) bins with AT >= 1 are used.

Model: HB_j | AT_j ~ BetaBinomial(AT_j, p_j, rho)   (mean p_j, intraclass correlation rho; rho = 0 is binomial)
  phase-agnostic : p_j = 0.5 + delta or 0.5 - delta with probability 1/2 each, independently per bin, delta in [0, 0.49]
  oriented       : p_j = 0.5 + delta for every bin (trusts GT phase), delta in [-0.49, 0.49]
Haplotype fraction h = 0.5 + delta. ML on a delta grid; profile-likelihood 95% interval = {delta : 2(llmax - ll) <= 3.841};
LR = 2(llmax - ll(delta = 0)). rho fixed at the empirical estimate (bins in base-called `all` segments, as in exp 06) and
at the sensitivity values 0, 0.02, 0.10.

Seeds: numpy PCG64 default_rng((20261029, stream, task_index[, scenario_index])), stream 1 = power grid, 2 = label
simulations. Replicates: 1,000 per scenario. Size is calibrated within each (profile, data ICC, model, analysis rho) cell
by using the 95th percentile of the h = 0.5 LR distribution as the critical value.
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
from scipy.special import gammaln  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

SEED = 20261029
REPS = 1000
LR_CUT = 3.841458820694124           # chi2_1 0.95
LR_CUT_ONESIDED = 2.705543454095404  # chi2_1 0.90 (boundary mixture for delta >= 0)
HALF = LR_CUT / 2
RHO_SENS = [0.0, 0.02, 0.10]
DMAX = 0.49
STEP_EST, STEP_SIM = 0.0025, 0.01
H_ALT = [0.6, 0.75, 0.9]
PHASES = ["none", "random", "switch"]
MAX_ELEMS = 6_000_000
MIN_BINS = 3


# ----------------------------------------------------------------------------- likelihood
def dgrid(step):
    k = int(round(DMAX / step))
    return np.arange(-k, k + 1) * step


def ll_grids(k, n, rho, dg):
    """k: (R, m) HB counts; n: (m,) AT. Returns oriented ll (R, G) over signed dg and agnostic ll (R, G//2+1) over dg >= 0.
    Terms constant in delta (binomial coefficient, normalizers depending on n and rho only) are dropped."""
    R, m = k.shape
    G = len(dg)
    p = 0.5 + dg
    ori = np.empty((R, G))
    agn = np.empty((R, G // 2 + 1))
    chunk = max(1, MAX_ELEMS // (G * m))
    for s in range(0, R, chunk):
        kk = k[s:s + chunk, None, :]
        if rho > 0:
            S = (1 - rho) / rho
            a, b = (p * S)[None, :, None], ((1 - p) * S)[None, :, None]
            logf = gammaln(kk + a) + gammaln(n - kk + b) - (gammaln(a) + gammaln(b))
        else:
            logf = kk * np.log(p)[None, :, None] + (n - kk) * np.log1p(-p)[None, :, None]
        ori[s:s + chunk] = logf.sum(-1)
        c = G // 2
        agn[s:s + chunk] = np.logaddexp(logf[:, c:, :], logf[:, c::-1, :]).sum(-1) - m * np.log(2)
    return ori, agn


def summarize(ll, d, interp=False):
    """ll: (R, G') over increasing d containing 0. Returns dict of arrays."""
    i0 = int(np.argmin(np.abs(d)))
    im = ll.argmax(1)
    mx = ll[np.arange(len(ll)), im]
    ok = ll >= (mx - HALF)[:, None]
    lo_i = ok.argmax(1)
    hi_i = ok.shape[1] - 1 - ok[:, ::-1].argmax(1)
    lo, hi = d[lo_i].astype(float), d[hi_i].astype(float)
    if interp:
        for r in range(len(ll)):
            li, hj = lo_i[r], hi_i[r]
            if li > 0:
                y0, y1 = ll[r, li - 1], ll[r, li]
                lo[r] = d[li - 1] + (mx[r] - HALF - y0) / (y1 - y0) * (d[li] - d[li - 1])
            if hj < len(d) - 1:
                y0, y1 = ll[r, hj], ll[r, hj + 1]
                hi[r] = d[hj] + (y0 - (mx[r] - HALF)) / (y0 - y1) * (d[hj + 1] - d[hj])
    return {"delta_hat": d[im], "lo": lo, "hi": hi, "LR": np.maximum(2 * (mx - ll[:, i0]), 0), "llmax": mx, "ll": ll}


def fit(k, n, rho, step=STEP_EST):
    dg = dgrid(step)
    ori, agn = ll_grids(np.asarray(k, float)[None, :], np.asarray(n, float), rho, dg)
    c = len(dg) // 2
    so, sa = summarize(ori, dg, True), summarize(agn, dg[c:], True)
    i75 = int(round(0.25 / step))
    return {"agn": {kk: (float(v[0]) if kk != "ll" else v[0]) for kk, v in sa.items()},
            "ori": {kk: (float(v[0]) if kk != "ll" else v[0]) for kk, v in so.items()},
            "LR75_agn": float(2 * (sa["llmax"][0] - agn[0, i75])), "LR75_ori": float(2 * (so["llmax"][0] - ori[0, c + i75]))}


# ----------------------------------------------------------------------------- simulation
def simulate_k(rng, n, delta, rho, phase, R):
    m = len(n)
    if phase == "none":
        s = np.ones((R, m))
    elif phase == "random":
        s = rng.choice([-1.0, 1.0], size=(R, m))
    else:   # one switch at a random interior bin: bins >= sw flipped
        sw = rng.integers(1, m, size=R) if m > 1 else np.ones(R, int)
        s = np.where(np.arange(m)[None, :] >= sw[:, None], -1.0, 1.0)
    p = 0.5 + s * delta
    if rho > 0:
        S = (1 - rho) / rho
        q = rng.beta(p * S, (1 - p) * S)
        return rng.binomial(n.astype(np.int64), q)
    return rng.binomial(n.astype(np.int64), p)


def sim_stats(rng, n, delta, rho_data, phase, rho_an_list, R=REPS):
    """Returns {rho_an: {'agn': summary, 'ori': summary}} on the simulation grid."""
    k = simulate_k(rng, n, delta, rho_data, phase, R)
    dg = dgrid(STEP_SIM)
    c = len(dg) // 2
    out = {}
    for ra in rho_an_list:
        ori, agn = ll_grids(k.astype(float), n.astype(float), ra, dg)
        out[ra] = {"agn": summarize(agn, dg[c:]), "ori": summarize(ori, dg)}
    return out


def power_task(args):
    task_idx, name, profile, n, rho_emp = args
    n = np.asarray(n, float)
    rows = []
    sc = 0

    def key(icc, model, ra):
        return (icc, model, ra)
    icc_list = [("0", 0.0), ("0.02", 0.02), ("0.10", 0.10), ("emp", rho_emp)]
    null = {}
    for icc_lab, icc in icc_list:
        rng = np.random.default_rng((SEED, 1, task_idx, sc)); sc += 1
        res = sim_stats(rng, n, 0.0, icc, "random", sorted({rho_emp, icc}))
        for ra_lab, ra in (("empirical", rho_emp), ("matched", icc)):
            if ra_lab == "matched" and icc_lab == "emp":
                continue
            for model in ("agn", "ori"):
                s = res[ra][model]
                crit = float(np.quantile(s["LR"], 0.95))
                null[(icc_lab, model, ra_lab)] = crit
                rows.append(dict(region=name, profile=profile, model=model, analysis_rho=ra_lab, data_icc=icc_lab, phase_scenario="n/a (h=0.5)",
                                 scenario="null", h_true=0.5, n_bins=len(n), AT_sum=n.sum(), reps=REPS, crit_LR_q95=crit,
                                 size_nominal_chi2_95=float((s["LR"] > LR_CUT).mean()), size_nominal_onesided=float((s["LR"] > LR_CUT_ONESIDED).mean()),
                                 power=np.nan, mc_se=np.nan, median_ci_width_h=float(np.median(s["hi"] - s["lo"])),
                                 mean_h_hat=float(np.mean(0.5 + s["delta_hat"])), ci_cover_truth=np.nan))
    alts = [("grid", h - 0.5, icc_lab, icc) for h in H_ALT for icc_lab, icc in icc_list[:3]]
    alts += [("anchored_MYC_delta_eff", None, "emp", rho_emp), ("anchored_half_MYC", None, "emp", rho_emp)]
    for scen, dlt, icc_lab, icc in alts:
        for phase in PHASES:
            rng = np.random.default_rng((SEED, 1, task_idx, sc)); sc += 1
            d_use = dlt if dlt is not None else (MYC_DELTA if scen == "anchored_MYC_delta_eff" else MYC_DELTA / 2)
            res = sim_stats(rng, n, d_use, icc, phase, sorted({rho_emp, icc}))
            for ra_lab, ra in (("empirical", rho_emp), ("matched", icc)):
                if ra_lab == "matched" and icc_lab == "emp":
                    continue
                for model in ("agn", "ori"):
                    s = res[ra][model]
                    crit = null[(icc_lab, model, ra_lab)]
                    pw = float((s["LR"] > crit).mean())
                    dtrue = d_use if model == "ori" else abs(d_use)
                    rows.append(dict(region=name, profile=profile, model=model, analysis_rho=ra_lab, data_icc=icc_lab, phase_scenario=phase,
                                     scenario=scen, h_true=0.5 + d_use, n_bins=len(n), AT_sum=n.sum(), reps=REPS, crit_LR_q95=crit,
                                     size_nominal_chi2_95=np.nan, size_nominal_onesided=np.nan, power=pw, mc_se=float(np.sqrt(pw * (1 - pw) / REPS)),
                                     median_ci_width_h=float(np.median(s["hi"] - s["lo"])), mean_h_hat=float(np.mean(0.5 + s["delta_hat"])),
                                     ci_cover_truth=float(((s["lo"] <= dtrue) & (s["hi"] >= dtrue)).mean()),
                                     power_nominal_chi2_95=float((s["LR"] > LR_CUT).mean())))
    return rows


def label_task(args):
    idx, n, rho = args
    n = np.asarray(n, float)
    out = {}
    for j, delta in enumerate((0.0, 0.25)):
        rng = np.random.default_rng((SEED, 2, idx, j))
        s = sim_stats(rng, n, delta, rho, "random", [rho])[rho]["agn"]
        out[j] = s["LR"]
    crit = float(np.quantile(out[0], 0.95))
    pw = float((out[1] > crit).mean())
    return idx, crit, pw


# ----------------------------------------------------------------------------- data
def load_sums(d, pth):
    cl = C.cell_clusters(d, pth)
    groups = ["all"] + sorted(set(cl))
    sums = {nm: np.zeros((len(groups), d.n_bins)) for nm in ("HA", "HB", "AT")}
    checks = {}
    cells = list(d.cells)
    for nm in ("HA", "HB", "AT"):
        t = pq.read_table(pth["allelic"] / f"{nm}.parquet")
        idx = t.to_pandas()[[]].reset_index()   # (chr, start, end) are stored in the pandas index
        ok = (idx["chr"].astype(str).to_numpy() == d.seq).all() and (idx["start"].to_numpy() == d.start).all() and (idx["end"].to_numpy() == d.end).all()
        names = [c for c in t.column_names if c not in ("chr", "start", "end", "__index_level_0__")]
        checks[nm] = {"coords_match_anndata": bool(ok), "barcode_set_equal": set(names) == set(cells), "column_order_equal_anndata": names == cells}
        if not ok or set(names) != set(cells):
            raise SystemExit(f"alignment failure {nm}: {checks[nm]}")
        for ci, c in enumerate(cells):
            col = t.column(c).to_numpy(zero_copy_only=False).astype(np.float64)
            sums[nm][0] += col
            g = groups.index(cl[ci])
            sums[nm][g] += col
        del t
    checks["AT_equals_HA_plus_HB_all_groups"] = bool(np.array_equal(sums["AT"], sums["HA"] + sums["HB"]))
    return groups, sums, checks


def estimate_rho(ha, hb, at, call):
    z2 = np.where(at > 0, (hb - ha) ** 2 / np.where(at > 0, at, 1), np.nan)
    base = (call == 1) & (at >= 20)
    strata = pd.qcut(at[base], 5)
    rows = []
    for s, g in pd.DataFrame({"z2": z2[base], "n": at[base], "s": strata}).groupby("s", observed=True):
        var = np.median(g.z2) / 0.4549
        rows.append({"AT_stratum": str(s), "bins": len(g), "median_AT": g.n.median(), "variance_inflation": var, "rho_hat": (var - 1) / (g.n.median() - 1)})
    tab = pd.DataFrame(rows)
    return float(np.median(tab.rho_hat.clip(lower=0))), tab


MYC_DELTA = None


def main():
    global MYC_DELTA
    t0 = time.time()
    pth = C.default_paths()
    d = C.H5AD(pth["h5ad"])
    spec = importlib.util.spec_from_file_location("exp04", HERE.parent / "04_fixed_region_contrasts" / "run.py")
    E4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(E4)
    groups, sums, checks = load_sums(d, pth)
    ha, hb, at = sums["HA"], sums["HB"], sums["AT"]
    kseg = pd.read_csv(pth["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    call = np.ones(d.n_bins)
    for r in kseg[kseg.group == "all"].itertuples():
        call[r.bin_start:r.bin_end + 1] = r.cn_call
    rho_emp, rho_tab = estimate_rho(ha[0], hb[0], at[0], call)
    s06 = json.loads((HERE.parent / "06_allelic_information" / "summary.json").read_text())
    MYC_DELTA = float(s06["model_null_calibration"]["control_delta_eff"]["lit_MYC_8q24"])
    rho_settings = [("empirical", rho_emp)] + [(f"{r:g}", r) for r in RHO_SENS]
    regs = E4.build_regions(d, pth)
    print(f"loaded {time.time()-t0:.0f}s rho_emp={rho_emp:.5f} (exp06: {s06['rho_used']:.5f}) MYC delta_eff={MYC_DELTA:.4f}", flush=True)

    # ---- estimates
    rows = []
    for ri, r in enumerate(regs):
        t = r["target"]
        called_any = bool((call[t] != 1).any())
        for gi, g in enumerate(groups):
            if r["family"] == "random" and g != "all":
                continue
            m = t[at[gi, t] >= 1]
            base = dict(family=r["family"], region=r["name"], region_index=ri, direction=r["direction"], group=g, n_target_bins=len(t),
                        n_bins_AT_ge1=len(m), AT_sum=at[gi, t].sum(), any_bin_called_all=called_any)
            if len(m) < MIN_BINS:
                rows.append(base | dict(rho_setting="empirical", rho=rho_emp))
                continue
            k, n = hb[gi, m], at[gi, m]
            for rl, rv in rho_settings:
                f = fit(k, n, rv)
                a, o = f["agn"], f["ori"]
                rows.append(base | dict(rho_setting=rl, rho=rv,
                                        h_agn=0.5 + a["delta_hat"], h_agn_lo=0.5 + a["lo"], h_agn_hi=0.5 + a["hi"], LR_agn=a["LR"], LR75_agn=f["LR75_agn"],
                                        h_ori=0.5 + o["delta_hat"], h_ori_lo=0.5 + o["lo"], h_ori_hi=0.5 + o["hi"], LR_ori=o["LR"], LR75_ori=f["LR75_ori"],
                                        absdelta_ori_minus_agn=abs(o["delta_hat"]) - a["delta_hat"]))
    est = pd.DataFrame(rows)
    print(f"estimates {time.time()-t0:.0f}s rows={len(est)}", flush=True)

    # ---- random-window null (pooled, all rho settings)
    pool_all = est[(est.group == "all") & (est.family == "random") & est.h_agn.notna()]
    rn = pool_all.pivot_table(index=["region", "any_bin_called_all", "n_bins_AT_ge1", "AT_sum"], columns="rho_setting",
                              values=["h_agn", "LR_agn", "h_ori", "LR_ori"]).reset_index()
    rn.columns = ["_".join(c for c in col if c) for col in rn.columns]
    unc = pool_all[(~pool_all.any_bin_called_all) & (pool_all.rho_setting == "empirical")]
    calib = {"n_random_windows_with_estimate": int(pool_all.region.nunique()), "n_uncalled": int(len(unc)), "rho_settings": {}}
    for rl, _ in rho_settings:
        u = pool_all[(~pool_all.any_bin_called_all) & (pool_all.rho_setting == rl)]
        a = pool_all[pool_all.rho_setting == rl]
        calib["rho_settings"][rl] = {
            "uncalled": {mdl: {"q95_LR": float(u[f"LR_{mdl}"].quantile(.95)), "median_LR": float(u[f"LR_{mdl}"].median()),
                               "frac_gt_3.841": float((u[f"LR_{mdl}"] > LR_CUT).mean()), "frac_gt_2.706": float((u[f"LR_{mdl}"] > LR_CUT_ONESIDED).mean()),
                               "frac_CI_excludes_0.5": float((u[f"h_{mdl}_lo"] > 0.5).mean()) if mdl == "agn" else float(((u.h_ori_lo > 0.5) | (u.h_ori_hi < 0.5)).mean())}
                         for mdl in ("agn", "ori")},
            "all_random": {mdl: {"frac_gt_3.841": float((a[f"LR_{mdl}"] > LR_CUT).mean())} for mdl in ("agn", "ori")}}
    emp95 = {"agn": float(unc.LR_agn.quantile(.95)), "ori": float(unc.LR_ori.quantile(.95))}
    calib["empirical_q95_uncalled_empirical_rho"] = emp95

    # ---- label simulations (pooled, agnostic, rho = empirical, random orientation)
    pe = est[(est.group == "all") & (est.rho_setting == "empirical")].copy()
    lab_idx = pe.index[pe.h_agn.notna()].tolist()
    jobs = []
    for ix in lab_idx:
        ri = int(pe.at[ix, "region_index"]); m = regs[ri]["target"]; m = m[at[0, m] >= 1]
        jobs.append((ri, at[0, m], rho_emp))
    with Pool(2) as pool:
        lab_res = {ri: (c, p) for ri, c, p in pool.imap_unordered(label_task, jobs, chunksize=4)}
    print(f"label sims {time.time()-t0:.0f}s", flush=True)
    pe["sim_crit_LR_q95"] = pe.region_index.map(lambda i: lab_res.get(i, (np.nan, np.nan))[0])
    pe["power_h075_random_orient"] = pe.region_index.map(lambda i: lab_res.get(i, (np.nan, np.nan))[1])

    def label(r):
        if np.isnan(r.h_agn):
            return "allele unobserved"
        if r.h_agn_lo > 0.5 and r.LR_agn > emp95["agn"]:
            return "evidence consistent with imbalance"
        if r.power_h075_random_orient < 0.8:
            return "insufficient information"
        if r.LR75_agn > LR_CUT:
            return "evidence against specified imbalance (h=0.75; conditional on AT-as-independent-observations)"
        return "compatible with balance"
    pe["label"] = pe.apply(label, axis=1)
    est = est.merge(pe[["family", "region", "group", "rho_setting", "sim_crit_LR_q95", "power_h075_random_orient", "label"]],
                    on=["family", "region", "group", "rho_setting"], how="left")

    # ---- power grid tasks: saved windows pooled + cluster0
    tasks = []
    ti = 0
    for ri, r in enumerate(regs):
        if r["family"] != "saved":
            continue
        for gname in ("all", "0"):
            gi = groups.index(gname)
            m = r["target"][at[gi, r["target"]] >= 1]
            if len(m) >= MIN_BINS:
                tasks.append((ti, r["name"], "pooled" if gname == "all" else "cluster0", at[gi, m], rho_emp))
            ti += 1
    with Pool(2) as pool:
        prow = [x for part in pool.imap_unordered(power_task, tasks) for x in part]
    pg = pd.DataFrame(prow)
    print(f"power grid {time.time()-t0:.0f}s rows={len(pg)}", flush=True)

    # ---- outputs
    est.drop(columns=["region_index"]).to_csv(HERE / "region_estimates.tsv", sep="\t", index=False)
    pg.to_csv(HERE / "power_grid.tsv", sep="\t", index=False)
    rn.to_csv(HERE / "random_null.tsv", sep="\t", index=False)
    lab_counts = pe.groupby(["family", "label"]).size().rename("regions").reset_index()
    saved_pool = est[(est.family.isin(["saved", "literature"])) & (est.group.isin(["all", "0"]))]
    summ = {"experiment": "10_betabinomial_allele_model", "method_record": "E:2",
            "seed": SEED, "generator": "numpy.random.default_rng (PCG64), seed tuples (20261029, stream, task_index[, scenario_index]); stream 1 = power grid, 2 = label sims",
            "replicates_per_scenario": REPS, "grid_step_estimation": STEP_EST, "grid_step_simulation": STEP_SIM, "delta_max": DMAX,
            "lr_cutoffs": {"chi2_1_95": LR_CUT, "boundary_onesided_chi2_1_90": LR_CUT_ONESIDED}, "min_bins": MIN_BINS, "min_AT_per_bin": 1,
            "rho_empirical": rho_emp, "rho_exp06": s06["rho_used"], "rho_sensitivity": RHO_SENS, "MYC_delta_eff_exp06": MYC_DELTA,
            "haplotype_fractions": [0.5] + H_ALT, "data_icc": [0, 0.02, 0.10, "empirical"], "phase_scenarios": PHASES,
            "label_rule": "unobserved (<3 bins AT>=1); consistent with imbalance (agnostic h CI excludes 0.5 and LR > uncalled-random-window q95); insufficient (simulated power at h=0.75, rho=empirical, random orientation < 0.8); against (LR for h=0.75 > 3.841 and power >= 0.8; conditional on AT as independent observations); else compatible with balance",
            "alignment_checks": checks, "random_window_calibration": calib, "label_counts": lab_counts.to_dict("records"),
            "n_regions": len(regs), "runtime_seconds": time.time() - t0, "environment": C.environment()}
    C.write_json(HERE / "summary.json", summ)
    make_report(est, pg, calib, emp95, lab_counts, checks, summ, rho_tab, rho_emp)
    print(f"done {time.time()-t0:.0f}s")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(saved_pool[saved_pool.rho_setting == "empirical"][["region", "group", "n_bins_AT_ge1", "AT_sum", "h_agn", "h_agn_lo", "h_agn_hi", "LR_agn", "h_ori", "h_ori_lo", "h_ori_hi", "LR_ori", "label"]].to_string())
    print(json.dumps(calib, indent=1))
    print(lab_counts.to_string())


GLOSSARY = {
    "h_agn / h_agn_lo / h_agn_hi": "Haplotype fraction h = 0.5 + delta under the phase-agnostic model (each bin's orientation unknown), ML estimate and 95% profile-likelihood interval; h is by construction >= 0.5, so an interval whose lower end is 0.5 is compatible with balance.",
    "h_ori / h_ori_lo / h_ori_hi": "Same, oriented model: every bin shares the GT-based HB orientation, so h can fall below 0.5 (HB fraction lower than HA).",
    "LR_agn / LR_ori": "Likelihood-ratio statistic 2(max ll - ll at delta = 0) for each model. Nominal chi-square(1) 95% cutoff is 3.84; the one-sided boundary cutoff is 2.71.",
    "LR75_agn": "2(max ll - ll at h = 0.75): evidence against the specified imbalance h = 0.75 (cutoff 3.84).",
    "absdelta_ori_minus_agn": "|delta_oriented| - delta_agnostic. Strongly negative means the agnostic model finds imbalance that a fixed orientation does not, i.e. bins disagree in orientation (or the state is mixed).",
    "rho / rho_setting": "Beta-binomial intraclass correlation used for the fit: 'empirical' (from base-called bins) or the sensitivity values 0, 0.02, 0.1.",
    "n_bins_AT_ge1 / AT_sum": "Bins in the region with at least one allele observation, and the sum of their AT (summed per-SNP-site DP, not unique fragments).",
    "any_bin_called_all": "True if any bin lies in a segment with cn_call != 1 in the group 'all' karyogram segments (random windows with False are 'uncalled').",
    "crit_LR_q95": "95th percentile of the LR in the simulated h = 0.5 data for the same profile, ICC, model and analysis rho: the size-calibrated critical value.",
    "sim_crit_LR_q95 / power_h075_random_orient": "Per-region label simulation (1,000 reps): critical value under the model null (rho empirical) and power to reject balance when h = 0.75 with random per-bin orientation.",
    "power / mc_se": "Fraction of replicates with LR above crit_LR_q95; binomial Monte Carlo standard error sqrt(p(1-p)/1000).",
    "size_nominal_chi2_95": "Rejection rate at h = 0.5 using the nominal cutoff 3.84 (should be about 0.05 if the model is calibrated; lower is conservative).",
    "median_ci_width_h": "Median width of the 95% profile interval in h across replicates (grid step 0.01).",
    "ci_cover_truth": "Fraction of replicate intervals containing the true delta (agnostic: |delta|).",
    "data_icc / analysis_rho": "Intraclass correlation used to generate the data ('emp' = empirical) / rho assumed in the fit ('empirical' fixed, or 'matched' to the generating ICC).",
    "phase_scenario": "none: all bins share one orientation; random: independent random orientation per bin; switch: one orientation switch at a random interior bin.",
    "scenario": "grid: prespecified h x ICC x phase grid; anchored_MYC_delta_eff / anchored_half_MYC: delta fixed at the exp 06 MYC effective |dBAF| (0.106) or half; null: h = 0.5.",
    "label": "Evidence label (see Method). 'against' is conditional on AT being independent observations.",
    "frac_CI_excludes_0.5": "Fraction of random windows whose interval excludes balance.",
}


def make_report(est, pg, calib, emp95, lab_counts, checks, summ, rho_tab, rho_emp):
    pe = est[(est.group == "all") & (est.rho_setting == "empirical")]
    # ---- figure 1: calibration, estimates, power
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    rnd = pe[(pe.family == "random") & pe.LR_agn.notna()]
    u = rnd[~rnd.any_bin_called_all]
    ax[0].hist(np.log10(u.LR_agn + 0.01), bins=40, color="#999")
    ax[0].axvline(np.log10(LR_CUT + .01), color="#9b2c2c", ls="--", label="nominal 3.84")
    ax[0].axvline(np.log10(emp95["agn"] + .01), color="#8a5a00", ls=":", label=f"empirical q95 = {emp95['agn']:.1f}")
    ax[0].set_xlabel("log10(LR + 0.01), uncalled random windows"); ax[0].legend(fontsize=7); ax[0].set_title("Empirical null (agnostic, rho empirical)")
    sp = pe[pe.family.isin(["saved", "literature"]) & pe.h_agn.notna()].reset_index(drop=True)
    y = np.arange(len(sp))
    ax[1].errorbar(sp.h_agn, y, xerr=[sp.h_agn - sp.h_agn_lo, sp.h_agn_hi - sp.h_agn], fmt="o", color="k", label="phase-agnostic", ms=4)
    ax[1].errorbar(sp.h_ori, y + .25, xerr=[sp.h_ori - sp.h_ori_lo, sp.h_ori_hi - sp.h_ori], fmt="s", color="#2f5d8a", label="oriented", ms=4)
    ax[1].set_yticks(y); ax[1].set_yticklabels(sp.region, fontsize=6); ax[1].axvline(.5, color="#bbb")
    ax[1].set_xlabel("haplotype fraction h (95% profile interval)"); ax[1].legend(fontsize=7); ax[1].set_title("Pooled estimates")
    q = pg[(pg.scenario == "grid") & (pg.profile == "pooled") & (pg.analysis_rho == "empirical") & (pg.data_icc == "0.02") & (pg.model == "agn") & (pg.phase_scenario == "random")]
    for nm, g in q.groupby("region"):
        ax[2].plot(g.h_true, g.power, "o-", label=nm, ms=3)
    ax[2].axhline(.8, color="#9b2c2c", ls="--", lw=1); ax[2].set_xlabel("true haplotype fraction"); ax[2].set_ylabel("power (size-calibrated 0.05)")
    ax[2].set_title("Agnostic power, ICC 0.02"); ax[2].legend(fontsize=6)
    fig.tight_layout(); b1 = C.fig_to_b64(fig); plt.close(fig)
    # ---- figure 2: oriented vs agnostic under phase scenarios, and ICC dependence at h = 0.6
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    q = pg[(pg.scenario == "grid") & (pg.profile == "pooled") & (pg.analysis_rho == "empirical") & (pg.h_true == 0.75) & (pg.data_icc == "0.02")]
    ph = ["none", "random", "switch"]
    for j, mdl in enumerate(["agn", "ori"]):
        v = q[q.model == mdl].groupby("phase_scenario").power.median().reindex(ph)
        ax[0].bar(np.arange(3) + j * .38, v.values, width=.36, label={"agn": "phase-agnostic", "ori": "oriented"}[mdl])
    ax[0].set_xticks(np.arange(3) + .19); ax[0].set_xticklabels(ph); ax[0].set_ylabel("median power across saved windows (h=0.75)"); ax[0].legend(fontsize=7)
    ax[0].set_title("Phase-switch scenarios")
    q = pg[(pg.scenario == "grid") & (pg.profile == "pooled") & (pg.analysis_rho == "empirical") & (pg.h_true == 0.6) & (pg.model == "agn") & (pg.phase_scenario == "random")]
    for icc, g in q.groupby("data_icc"):
        ax[1].plot(range(len(g)), g.sort_values("region").power, "o-", label=f"ICC {icc}")
    ax[1].set_xticks(range(7)); ax[1].set_xticklabels(sorted(q.region.unique()), rotation=60, fontsize=6, ha="right"); ax[1].axhline(.8, color="#9b2c2c", ls="--", lw=1)
    ax[1].set_ylabel("power at h = 0.6"); ax[1].legend(fontsize=7); ax[1].set_title("Dependence on data ICC (agnostic)")
    fig.tight_layout(); b2 = C.fig_to_b64(fig); plt.close(fig)

    cr = calib["rho_settings"]
    crs = pd.DataFrame([{"rho_setting": rl, "model": m, "n_uncalled": calib["n_uncalled"], "q95_LR": v["uncalled"][m]["q95_LR"], "frac_LR_gt_3.84": v["uncalled"][m]["frac_gt_3.841"],
                         "frac_LR_gt_2.71": v["uncalled"][m]["frac_gt_2.706"], "all_300_gt_3.84": v["all_random"][m]["frac_gt_3.841"]} for rl, v in cr.items() for m in ("agn", "ori")])
    cols = ["region", "group", "n_bins_AT_ge1", "AT_sum", "h_agn", "h_agn_lo", "h_agn_hi", "LR_agn", "h_ori", "h_ori_lo", "h_ori_hi", "LR_ori", "absdelta_ori_minus_agn"]
    sv = est[est.family.isin(["saved", "literature"]) & est.group.isin(["all", "0"]) & (est.rho_setting == "empirical") & est.h_agn.notna()][cols]
    sens = est[est.family.isin(["saved", "literature"]) & (est.group == "all") & est.h_agn.notna()][["region", "rho_setting", "h_agn", "h_agn_lo", "h_agn_hi", "LR_agn", "h_ori"]]
    sens = sens.assign(rho_setting=pd.Categorical(sens.rho_setting, ["empirical", "0", "0.02", "0.1"])).sort_values(["region", "rho_setting"])
    hd = pg[(pg.scenario == "grid") & (pg.analysis_rho == "empirical") & (pg.h_true == 0.75) & (pg.data_icc == "0.02") & (pg.profile == "pooled")][
        ["region", "model", "phase_scenario", "n_bins", "AT_sum", "crit_LR_q95", "power", "mc_se", "median_ci_width_h"]]
    pw_icc = pg[(pg.scenario == "grid") & (pg.analysis_rho == "empirical") & (pg.profile == "pooled") & (pg.model == "agn") & (pg.phase_scenario == "random")].pivot_table(
        index=["region", "data_icc"], columns="h_true", values="power").round(3).reset_index()
    pw_icc.columns = [str(c) if not isinstance(c, float) else f"power_h{c}" for c in pw_icc.columns]
    c0 = pg[(pg.scenario == "grid") & (pg.analysis_rho == "empirical") & (pg.profile == "cluster0") & (pg.model == "agn") & (pg.phase_scenario == "random") & (pg.h_true == 0.75)].pivot_table(
        index="region", columns="data_icc", values="power").round(3).reset_index()
    c0.columns = ["region"] + [f"power_h0.75_ICC{c}" for c in c0.columns[1:]]
    anc = pg[pg.scenario.str.startswith("anchored") & (pg.profile == "pooled")].pivot_table(index=["region", "model", "scenario"], columns="phase_scenario", values="power").round(3).reset_index()
    nul = pg[(pg.scenario == "null") & (pg.profile == "pooled") & (pg.analysis_rho == "empirical")].pivot_table(index=["region", "model"], columns="data_icc", values="size_nominal_chi2_95").round(3).reset_index()
    nul.columns = ["region", "model"] + [f"size_nominal_ICC{c}" for c in nul.columns[2:]]
    cov = pg[(pg.scenario == "grid") & (pg.analysis_rho == "matched") & (pg.profile == "pooled") & (pg.model == "agn") & (pg.phase_scenario == "random") & (pg.h_true == 0.6)][
        ["region", "data_icc", "power", "mc_se", "median_ci_width_h", "ci_cover_truth"]]

    why = {"evidence consistent with imbalance": "interval excludes 0.5 and LR exceeds the uncalled random-window 95th percentile",
           "compatible with balance": "interval includes 0.5; h = 0.75 not excluded or power adequate",
           "insufficient information": "power for h = 0.75 below 0.8", "allele unobserved": "fewer than 3 bins with AT >= 1"}
    det = []
    for r in pe[pe.family.isin(["saved", "literature"])].itertuples():
        lab = r.label if isinstance(r.label, str) else "allele unobserved"
        short = "evidence against specified imbalance" if lab.startswith("evidence against") else lab
        if short.startswith("evidence against"):
            reason = f"h = 0.75 excluded (LR75 = {r.LR75_agn:.0f}), simulated power {r.power_h075_random_orient:.2f}; estimate h = {r.h_agn:.3f} [{r.h_agn_lo:.3f}, {r.h_agn_hi:.3f}], LR = {r.LR_agn:.1f} below the empirical q95 {emp95['agn']:.1f}; conditional on AT as independent observations"
        elif short.startswith("evidence consistent"):
            reason = f"h = {r.h_agn:.3f} [{r.h_agn_lo:.3f}, {r.h_agn_hi:.3f}], LR = {r.LR_agn:.1f} > empirical q95 {emp95['agn']:.1f}; oriented h = {r.h_ori:.3f} (orientation partly inconsistent between bins if it differs)"
        else:
            reason = why.get(short, "")
        det.append({"region": r.region, "family": r.family, "bins_AT_ge1": r.n_bins_AT_ge1, "label": short, "reason": reason})
    det = pd.DataFrame(det)
    sv_pool = pe[pe.family == "saved"]
    n_cons_rand = int(((pe.family == "random") & pe.label.astype(str).str.startswith("evidence consistent")).sum())
    myc = pe[pe.region == "lit_MYC_8q24"].iloc[0]; q5 = pe[pe.region == "lit_5q13_qter"].iloc[0]; q19 = pe[pe.region == "lit_19q13_qter"].iloc[0]; p9 = pe[pe.region == "lit_9p21_hemizygous_zone_approx"].iloc[0]
    e_ = cr["empirical"]["uncalled"]["agn"]; e2 = cr["0.02"]["uncalled"]["agn"]; e0 = cr["0"]["uncalled"]["agn"]
    ori_sw = pg[(pg.scenario == "grid") & (pg.analysis_rho == "empirical") & (pg.profile == "pooled") & (pg.h_true == 0.75) & (pg.data_icc == "0.02") & (pg.model == "ori")].groupby("phase_scenario").power.median()
    anc_agn = pg[(pg.scenario == "anchored_MYC_delta_eff") & (pg.profile == "pooled") & (pg.model == "agn") & (pg.phase_scenario == "random")].power
    anc_half = pg[(pg.scenario == "anchored_half_MYC") & (pg.profile == "pooled") & (pg.model == "agn") & (pg.phase_scenario == "random")].power
    p06 = pg[(pg.scenario == "grid") & (pg.analysis_rho == "empirical") & (pg.profile == "pooled") & (pg.model == "agn") & (pg.phase_scenario == "random") & (pg.h_true == 0.6)]
    p06_0, p06_2, p06_10 = [p06[p06.data_icc == k].power for k in ("0", "0.02", "0.10")]
    verdict = ("The beta-binomial haplotype-fraction model recovers the strong controls (MYC h = %.3f, 5q13-qter h = %.3f) and finds no allelic imbalance beyond the empirical null at the six saved chr7/chr9 windows, but its nominal "
               "chi-square calibration fails in uncalled random windows (%.0f%% exceed 3.84), so it is consistent with serving only as an empirical-null-calibrated diagnostic, not yet as an evidence layer."
               % (myc.h_agn, q5.h_agn, 100 * e_["frac_gt_3.841"]))
    kf = [
        f"MYC control: phase-agnostic h = {myc.h_agn:.3f} (95% interval {myc.h_agn_lo:.3f}-{myc.h_agn_hi:.3f}, LR = {myc.LR_agn:.0f}); 5q13-qter h = {q5.h_agn:.3f} ({q5.h_agn_lo:.3f}-{q5.h_agn_hi:.3f}); 19q13-qter h = {q19.h_agn:.3f} ({q19.h_agn_lo:.3f}-{q19.h_agn_hi:.3f}). The 9p21 zone is not distinguishable from balance (h = {p9.h_agn:.3f}, {p9.h_agn_lo:.3f}-{p9.h_agn_hi:.3f}).",
        "Saved windows (pooled, rho empirical): " + "; ".join(f"{r.region} h = {r.h_agn:.3f} ({r.h_agn_lo:.3f}-{r.h_agn_hi:.3f})" for r in sv_pool.itertuples() if r.region != 'chr9_70_100') + f". chr9_70_100 (comparison window, 282 bins) h = {sv_pool[sv_pool.region=='chr9_70_100'].h_agn.iloc[0]:.3f} with LR = {sv_pool[sv_pool.region=='chr9_70_100'].LR_agn.iloc[0]:.0f}, small but beyond the empirical null given AT = 131,184.",
        f"Calibration is poor at the empirical rho {rho_emp:.4f}: {100*e_['frac_gt_3.841']:.1f}% of {calib['n_uncalled']} uncalled random windows exceed LR 3.84 (agnostic; oriented {100*cr['empirical']['uncalled']['ori']['frac_gt_3.841']:.1f}%), versus 5% nominal; the empirical 95th percentile is LR = {emp95['agn']:.1f}. With rho = 0.02 the exceedance falls to {100*e2['frac_gt_3.841']:.1f}% and with rho = 0 it rises to {100*e0['frac_gt_3.841']:.1f}%.",
        f"Power at h = 0.75 (size-calibrated 0.05, ICC 0.02, agnostic): 1.00 in all 7 saved windows (1,000 reps, MC SE at most 0.01). At h = 0.6 power depends strongly on ICC: median {p06_0.median():.2f} (ICC 0), {p06_2.median():.2f} (ICC 0.02), {p06_10.median():.2f} (ICC 0.10); range {p06_10.min():.2f}-{p06_10.max():.2f} at ICC 0.10.",
        f"Phase switches: the agnostic model's power is the same across orientation scenarios by construction, but the oriented model at h = 0.75 falls from 1.00 (none) to {ori_sw['switch']:.2f} (one switch) and {ori_sw['random']:.2f} (random orientation), median across saved windows.",
        f"Attenuated scenario (MYC effective |dBAF| = {MYC_DELTA:.3f}): agnostic power {anc_agn.min():.2f}-{anc_agn.max():.2f} across saved windows; at half that imbalance {anc_half.min():.2f}-{anc_half.max():.2f}, so weaker real-world imbalance would often be missed in the shortest windows.",
        f"Labels (pooled, rho empirical): saved windows 6 'against h = 0.75' and 1 'consistent with imbalance' (chr9_70_100); literature 3 consistent / 1 against (9p21); 300 random windows: {n_cons_rand} 'consistent with imbalance' ({100*n_cons_rand/300:.0f}%) even after using the empirical threshold, and no region is labelled 'insufficient' because the h = 0.75 alternative is large.",
        f"Oriented vs agnostic: for MYC the oriented estimate is h = {myc.h_ori:.3f} (below 0.5) against agnostic {myc.h_agn:.3f}, and 5q13-qter oriented {q5.h_ori:.3f} versus agnostic {q5.h_agn:.3f}; the large drop in oriented LR (e.g. 5q13-qter {q5.LR_ori:.0f} vs {q5.LR_agn:.0f}) is consistent with GT-based orientation not being consistent between bins.",
    ]
    sections = [
        ("Question and why it matters",
         "<p>Method record <b>E:2 &lsquo;Phase-aware beta-binomial allele model&rsquo;</b> (research/2026-10-29-low-coverage/synthesis/methods_inventory.json; e-per-cell-evidence/findings.md). "
         "It asks whether modelling the allele counts HB | AT as beta-binomial with a haplotype fraction h, with and without trusting the GT-based phase orientation, gives (i) interpretable estimates with intervals, "
         "(ii) a null that is calibrated on regions that were not called, and (iii) enough power at the prespecified alternatives (NEXT_EXPERIMENTS.md section 2, Allelic arm). "
         "The decision it informs is whether E:2 can enter an allelic evidence layer next to depth-based calls, and with which evidence labels. Exp 06 used a phase-agnostic bin-level statistic S and found strong attenuation relative to nominal counts; this experiment replaces S by a likelihood with an explicit orientation model.</p>"),
        ("Inputs",
         "<ul><li>HA / HB / AT parquet (aligned_matrices; 25,337 bins x 3,460 barcodes, float32 integer-valued), joined to AnnData by barcode name; coordinates and barcode sets checked equal; AT = HA + HB verified for the pooled and all five cluster sums.</li>"
         "<li>Groups: pooled (3,460 cells) and WNN clusters (wnn_leiden_0.3) 0-4 (1,286 / 833 / 603 / 556 / 182 cells).</li>"
         "<li>Regions (exp 04 families): 7 saved windows (6 chr7/chr9 windows plus chr9_70_100 as comparison), 67 called short segments, 4 literature loci, 300 random 25-bin windows (226 uncalled in the <code>all</code> karyogram).</li>"
         "<li>karyogram_cluster_segments.tsv (group all) for the called/uncalled status and for base-called bins used for rho; exp 06 summary.json for the MYC effective |dBAF|.</li>"
         "<li>Units: counts per 100 kb bin. <b>Not available</b>: unique fragments or molecules (AT sums per-SNP-site depth), phase-set IDs, a verified cross-bin phase orientation, a matched normal, copy-number truth.</li></ul>" + f"<pre>{checks}</pre>"),
        ("Method",
         "<p>Step 1. Per region and group take bins with AT &ge; 1. Step 2. Fit two models by grid maximum likelihood (delta step 0.0025, |delta| &le; 0.49).</p>"
         "<pre>HB_j | AT_j ~ BetaBinomial(AT_j, p_j, rho),  a = p(1-rho)/rho,  b = (1-p)(1-rho)/rho  (rho = 0: Binomial)\n"
         "phase-agnostic: p_j = 0.5 + delta or 0.5 - delta, prob 1/2 each, independent per bin; delta in [0, 0.49]\n"
         "     log L(delta) = sum_j log[ 0.5 f(HB_j; p=0.5+delta) + 0.5 f(HB_j; p=0.5-delta) ]\n"
         "oriented:       p_j = 0.5 + delta for all j; delta in [-0.49, 0.49]\n"
         "h = 0.5 + delta;  LR = 2 (max log L - log L(0));  95% interval = {delta: 2 (max log L - log L(delta)) <= 3.841}\n"
         "evidence against h = 0.75: LR75 = 2 (max log L - log L(0.25)) > 3.841</pre>"
         f"<p>Step 3. rho is fixed at the empirical value {rho_emp:.5f} (median over five AT strata of bins in base-called <code>all</code> segments, cn_call = 1, AT &ge; 20; identical to exp 06) and, as sensitivity, at 0, 0.02, 0.10.</p>"
         + C.html_table(rho_tab) +
         f"<p>Step 4 (calibration). The same fit is applied to the 300 random windows pooled; uncalled windows (no bin with cn_call &ne; 1) give an empirical null for the LR; its 95th percentile (agnostic {emp95['agn']:.2f}, oriented {emp95['ori']:.2f}) is the threshold for the &lsquo;consistent with imbalance&rsquo; label.</p>"
         "<p>Step 5 (prospective power). Generator numpy PCG64 <code>default_rng((20261029, stream, task, scenario))</code>, stream 1 = power grid, 2 = label simulations; 1,000 replicates per scenario. Data: each saved window's own pooled AT profile (and cluster 0's); "
         "h in {0.5, 0.6, 0.75, 0.9} x data ICC {0, 0.02, 0.10} x phase {none, random per bin, one switch at a random interior bin}; plus the attenuated scenario with delta = MYC effective |dBAF| and half of it (ICC = empirical). "
         "Analysis rho either the fixed empirical value (&lsquo;empirical&rsquo;) or the generating ICC (&lsquo;matched&rsquo;). Size is calibrated per cell by the 95th percentile of the h = 0.5 LR (shared across phase scenarios since the sign is irrelevant at h = 0.5); "
         "power = fraction of alternative replicates above it, with binomial MC SE. Simulation grid step 0.01.</p>"
         "<p>Step 6 (labels, pooled, agnostic, rho empirical): <b>allele unobserved</b> (&lt; 3 bins with AT &ge; 1); <b>evidence consistent with imbalance</b> (interval excludes 0.5 and LR above the empirical random-window q95); "
         "<b>insufficient information</b> (simulated power &lt; 0.8 for h = 0.75 with random orientation); <b>evidence against specified imbalance</b> (LR75 &gt; 3.84 and power &ge; 0.8); otherwise <b>compatible with balance</b>. "
         "The &lsquo;against&rsquo; label is conditional on treating AT as independent observations, which exp 06 showed is optimistic.</p>"),
        ("Results", "<p>Full tables: region_estimates.tsv, power_grid.tsv, random_null.tsv, summary.json.</p>", 2),
        ("Region estimates: saved windows and literature loci", C.html_table(sv, max_rows=40) + C.img(b1, "estimates, calibration, power") +
         C.caption("Left: LR in uncalled random windows (log scale) with the nominal cutoff and the empirical 95th percentile. Middle: pooled haplotype fraction with 95% profile interval, agnostic (black) and oriented (blue). "
                   "Right: size-calibrated power of the agnostic model by true h. Look for: MYC and 5q13-qter far from 0.5; chr7/chr9 windows overlapping 0.5; agnostic and oriented diverging where bins disagree on orientation. Table rows with group 0 are cluster 0 only."), 3),
        ("Calibration in random windows", C.html_table(crs) +
         C.caption("Fraction of uncalled random windows with LR above the nominal cutoff, by rho and model. A calibrated test gives about 0.05. Values far above 0.05 mean the fixed-rho beta-binomial null is too narrow for real windows; values near 0 mean it is too wide. The empirical rho sits on the anti-conservative side."), 3),
        ("Power at the prespecified alternative", C.html_table(hd, max_rows=40) + C.img(b2, "phase scenarios and ICC dependence") +
         C.caption("Table: pooled saved windows at h = 0.75, ICC 0.02, analysis rho empirical. Left panel: oriented vs agnostic power by phase scenario. Right panel: agnostic power at h = 0.6 by data ICC. Look for the oriented model losing power under random or switching orientation while the agnostic model does not, and for the collapse at ICC 0.10."), 3),
        ("Power across the h x ICC grid (agnostic, random orientation)", C.html_table(pw_icc, max_rows=40) +
         C.caption("Power at size-calibrated 0.05 (1,000 replicates; MC SE up to 0.016) by data ICC and true h. The analysis rho is the fixed empirical value, so ICC 0.10 rows show the cost of misspecified overdispersion. Cluster 0 (smaller depth) at h = 0.75 is shown next."), 3),
        ("Cluster 0 profile", C.html_table(c0) + C.caption("Power at h = 0.75 using cluster 0's own AT profile by data ICC; lower than pooled where depth is lower."), 3),
        ("Attenuated (MYC-anchored) scenario", C.html_table(anc, max_rows=60) +
         C.caption("Power when delta equals the exp 06 MYC effective |dBAF| (0.106) or half of it. Columns are phase scenarios. Look for half-strength imbalance falling below 0.8 in the short windows."), 3),
        ("Size and interval behaviour", C.html_table(nul, max_rows=40) + C.html_table(cov) +
         C.caption("Top: rejection rate at h = 0.5 with the nominal 3.84 cutoff when data are generated from the stated ICC and fitted with the empirical rho (about 0.05 expected if calibrated). Bottom: matched-rho analysis at h = 0.6, with interval width and coverage of the true |delta|."), 3),
        ("Per-region detail", C.html_table(det, max_rows=20) +
         C.caption(f"One row per saved window and literature control. Label counts over all families: {lab_counts.to_dict('records')}."), 2),
        ("Sensitivity", C.html_table(sens, max_rows=80) +
         C.caption("Pooled estimates for the empirical rho and the prespecified rho values. Larger assumed rho widens intervals and pulls weak estimates to 0.5, so weak apparent imbalance (chr7_8.3_10.7, chr7_70.0_73.3, chr9_70_100) is sensitive to rho, whereas MYC and 5q13-qter are not.") +
         "<p>Across the prespecified alternatives: h = 0.75 and 0.9 are detected with power &ge; 0.83 in all windows at all ICCs (lowest: chr9_22.5_26.6, ICC 0.10, h = 0.75); h = 0.6 is detected reliably only at low ICC and in deep windows.</p>"),
        ("Limits and alternative explanations",
         "<ul><li>AT is summed site depth: one fragment can contribute at several SNP sites, so counts are not independent. This inflates the effective n and is why the &lsquo;against&rsquo; label (and every power value) is conditional on that assumption; exp 06 showed the assumption is optimistic.</li>"
         "<li>Calibration failure: with the empirical rho the null is anti-conservative in real uncalled windows; the empirical q95 is therefore used for &lsquo;consistent with imbalance&rsquo;. Even so, "
         f"{n_cons_rand} of 300 random windows receive that label, which indicates widespread allelic structure or model misfit outside depth calls.</li>"
         "<li>A significant LR is also consistent with technical causes (genotype or phasing errors, reference bias, ascertainment, mixed cell states), not only a copy-number change.</li>"
         "<li>Phase orientation across bins is GT-based (Eagle) and unverified; agreement between oriented and agnostic fits would not verify phase, and there are no phase-set IDs.</li>"
         "<li>Fixed rho does not capture region-specific overdispersion; grid resolution 0.01 in simulations; labels are for pooled all-cell data and dilute subclonal events.</li>"
         "<li>The phase-agnostic power is independent of the orientation scenario by construction. MCF10-A is development data, not CN truth; the power values describe the stated model scenarios, not a validated operating threshold.</li></ul>"),
    ]
    supports = ["Strong allelic signal at MYC and 5q13-qter is recovered with intervals well away from 0.5 (consistent with a copy-number gain in cells in which sites remain observable).",
                "No imbalance of the h = 0.75 type in the six saved chr7/chr9 windows; these windows have power ~1 for that alternative under the model.",
                "Oriented and phase-agnostic fits diverge in the gain controls, consistent with orientation inconsistency between bins."]
    not_supports = ["Calibrated nominal p-values: uncalled random windows exceed the nominal cutoff far more often than 5% at the empirical rho.",
                    "Any claim that the saved windows are balanced at weaker imbalance (h around 0.55-0.6): power there is limited and depends on ICC.",
                    "Verified phase orientation, or an absolute copy-number state, from these counts.",
                    "The &lsquo;against&rsquo; label as a reliable exclusion: it relies on AT counting independent observations."]
    decision = ("E:2 can enter an evidence layer now only as an <b>exploratory, empirically thresholded</b> diagnostic (phase-agnostic LR compared with an uncalled-window null; oriented fit reported as a phase-consistency check). "
                "The nominal-chi-square version should not be used, the &lsquo;against specified imbalance&rsquo; label should stay flagged as conditional, and labels need a calibrated rho (a value between the empirical 0.0038 and 0.02 would be needed for near-nominal size; not tested here) or a unique-molecule count before any threshold is treated as operating guidance. "
                "This adds a likelihood-based allelic arm to the recommendation in research/2026-11-04-cloud-method-review/recommendation.md but does not change the depth-based recommendation.")
    nxt = ["Re-run calibration with rho chosen so that uncalled random windows give about 5% exceedance, and report the resulting power loss (rho between 0.0038 and 0.02).",
           "Obtain unique-fragment or molecule counts (or phase blocks) and repeat the fit with them in place of AT.",
           "Add a region-level random-effect or bin-block bootstrap to replace fixed rho.",
           "Fit chr9_70_100 and the 14 called short segments labelled 'consistent' on cluster subsets to see whether imbalance is clonal or cluster-specific.",
           "Check orientation consistency with an external phased reference (phase sets) before trusting the oriented model."]
    C.write_html_report(HERE / "report.html", "Beta-binomial allele model", "Experiment 10 &middot; method record E:2 &middot; phase-aware haplotype-fraction estimation, calibration and power",
                        sections, {"command": "PYEPI_SCIENTIFIC_PACKAGE=/home/user/pyepi-upstream/pyEpiAneufinder /home/user/envs/cnv310/bin/python -I experiments/10_betabinomial_allele_model/run.py [--report-only]",
                                   "summary": "summary.json", "environment": C.environment()},
                        key_findings=kf, verdict=verdict, glossary=GLOSSARY, supports=supports, not_supports=not_supports, decision_impact=decision, next_steps=nxt)


def report_only():
    global MYC_DELTA
    summ = json.loads((HERE / "summary.json").read_text())
    MYC_DELTA = summ["MYC_delta_eff_exp06"]
    est = pd.read_csv(HERE / "region_estimates.tsv", sep="\t")
    est["rho_setting"] = est["rho_setting"].astype(str)
    pg = pd.read_csv(HERE / "power_grid.tsv", sep="\t", keep_default_na=False, na_values=[""], dtype={"data_icc": str})
    pth = C.default_paths()
    d = C.H5AD(pth["h5ad"])
    groups, sums, checks = load_sums(d, pth)
    kseg = pd.read_csv(pth["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    call = np.ones(d.n_bins)
    for r in kseg[kseg.group == "all"].itertuples():
        call[r.bin_start:r.bin_end + 1] = r.cn_call
    rho_emp, rho_tab = estimate_rho(sums["HA"][0], sums["HB"][0], sums["AT"][0], call)
    lab = est[(est.group == "all") & (est.rho_setting == "empirical")].groupby(["family", "label"]).size().rename("regions").reset_index()
    make_report(est, pg, summ["random_window_calibration"], summ["random_window_calibration"]["empirical_q95_uncalled_empirical_rho"], lab, checks, summ, rho_tab, rho_emp)


if __name__ == "__main__":
    report_only() if "--report-only" in sys.argv else main()
