"""Experiment 13: transparent joint depth + allele HMM on the `all` and `cluster0` raw pseudobulks.

EXPLORATORY DIAGNOSTIC (method record B:3, "covariance-aware joint segmentation / Bayesian HMM
(proposed)"). Hidden states are RELATIVE depth levels crossed with an allelic state. They are not
copy number, accessibility and dosage are confounded, and AT (summed SNP-site DP) is not independent of
depth. Nothing here is validated; it isolates what the allele layer and the between-bin baseline
variance term contribute.

Command:  /home/user/envs/cnv310/bin/python -I research/2026-11-04-cloud-method-review/experiments/13_joint_hmm/run.py

Seeds: halves = np.random.default_rng((20261104, 13)) via C.stratified_halves on WNN labels. The HMM
itself is deterministic (no random numbers).
"""
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
from scipy.special import gammaln  # noqa: E402
from scipy.stats import betabinom, trim_mean  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import common as C  # noqa: E402

# ----------------------------------------------------------------------------- prespecified settings
R_GRID = np.array([0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0])      # relative depth levels (x GC-decile baseline)
HAZARD = 0.01          # switch hazard per 100 kb of bp distance to the previous retained bin
P_SW_CAP = 0.5         # upper bound on the switch probability (large QC-removed gaps)
PI_BAL_R1, PI_IMB_R1 = 0.55, 0.05   # restart/switch distribution; remaining 0.40 shared by the other 12 states
CN_NEUTRAL_DELTA_NOMINAL = 0.25     # nominal |delta| for the copy-neutral-imbalance state (arbitrary, stated)
MIN_AT_RHO = 20        # exp 06 convention for rho estimation
POSTERIOR_FLAG = 0.8
SEED_HALVES = (C.SEED, 13)
GROUPS = ["all", "cluster0"]
VARIANTS = {"full": dict(tau=True, allele=True, base="cn1"), "no_baseline_var": dict(tau=False, allele=True, base="cn1"),
            "no_allele": dict(tau=True, allele=False, base="cn1"),
            "baseline_all_bins": dict(tau=True, allele=True, base="all"),       # extra: decile medians over all bins instead of cn_call==1 bins
            "baseline_median_only": dict(tau=True, allele=True, base="cn1_median")}  # extra: no mean-matching scalar (pure decile median)
LIT = [("lit_9p21_hemizygous_zone_approx", "chr9", 20.0e6, 24.0e6), ("lit_MYC_8q24", "chr8", 127.0e6, 128.5e6),
       ("lit_19q13_qter", "chr19", 35.0e6, 1e9), ("lit_5q13_qter", "chr5", 70.0e6, 1e9)]

# joint state s = 2*depth_index + allele (0 balanced, 1 imbalanced)
S_DEPTH = np.repeat(np.arange(len(R_GRID)), 2)
S_R = R_GRID[S_DEPTH]
S_IMB = np.tile([0, 1], len(R_GRID))
NS = len(S_R)
S_NAME = [f"r{S_R[s]:g}_{'imb' if S_IMB[s] else 'bal'}" for s in range(NS)]
LOSS, GAIN = S_R < 1, S_R > 1


def nominal_delta(r):
    """|delta| expected for a clonal one-allele change giving relative depth r = copies/2: |r-1|/(2r), capped at 0.5.
    r=1 (copy-neutral imbalance) uses CN_NEUTRAL_DELTA_NOMINAL."""
    r = np.asarray(r, float)
    return np.where(np.isclose(r, 1.0), CN_NEUTRAL_DELTA_NOMINAL, np.minimum(np.abs(r - 1) / (2 * r), 0.5))


# ----------------------------------------------------------------------------- data
def load_pseudobulks(d, lab):
    half = C.stratified_halves(lab, SEED_HALVES)
    masks = {"all": np.ones(d.n_cells, bool), "all_A": half, "all_B": ~half,
             "cluster0": lab == "0", "cluster0_A": (lab == "0") & half, "cluster0_B": (lab == "0") & ~half}
    names = list(masks)
    sums = np.zeros((len(names), d.n_bins))
    tot = np.zeros(d.n_cells)
    for lo in range(0, d.n_cells, 250):
        hi = min(d.n_cells, lo + 250)
        blk = d.rows("layers/raw", lo, hi)
        tot[lo:hi] = blk.sum(1)
        for k, nm in enumerate(names):
            sel = masks[nm][lo:hi]
            if sel.any():
                sums[k] += blk[sel].sum(0)
    y = {nm: sums[k] for k, nm in enumerate(names)}
    e2 = {nm: float((tot[masks[nm]] ** 2).sum()) for nm in names}       # sum of squared cell totals
    return masks, y, tot, e2


def load_allelic(d, pth, masks):
    out = {}
    cells = list(d.cells)
    checks = {}
    for nm in ("HB", "AT"):
        t = pq.read_table(pth["allelic"] / f"{nm}.parquet")
        idx = t.to_pandas()[[]].reset_index()
        cols = [c for c in t.column_names if c not in ("chr", "start", "end", "__index_level_0__")]
        ok = (idx["chr"].astype(str).to_numpy() == d.seq).all() and (idx["start"].to_numpy() == d.start).all() \
            and (idx["end"].to_numpy() == d.end).all()
        checks[nm] = {"coords_match_anndata": bool(ok), "barcode_set_equal": set(cols) == set(cells)}
        if not ok or set(cols) != set(cells):
            raise SystemExit(f"alignment failure {nm}: {checks[nm]}")
        M = np.column_stack([t.column(c).to_numpy(zero_copy_only=False) for c in cells]).astype(np.float32)
        del t
        for g in GROUPS:
            out[(nm, g)] = M[:, masks[g]].sum(1, dtype=np.float64)
        del M
    return out, checks


# ----------------------------------------------------------------------------- model pieces
def nb_logpmf(y, mu, phi):
    n = 1.0 / phi
    return y * np.log(mu / (mu + n)) + n * np.log(n / (n + mu)) + gammaln(y + n) - gammaln(n) - gammaln(y + 1)


def estimate_rho(hb, at, call):
    ha = at - hb
    z2 = np.where(at > 0, (hb - ha) ** 2 / np.where(at > 0, at, 1), np.nan)
    base = (call == 1) & (at >= MIN_AT_RHO)
    strata = pd.qcut(at[base], 5)
    rows = []
    for s, g in pd.DataFrame({"z2": z2[base], "n": at[base], "s": strata}).groupby("s", observed=True):
        var = np.median(g.z2) / 0.4549
        rows.append({"AT_stratum": str(s), "bins": len(g), "median_AT": g.n.median(), "variance_inflation": var,
                     "rho_hat": (var - 1) / (g.n.median() - 1)})
    tab = pd.DataFrame(rows)
    return float(np.median(tab.rho_hat.clip(lower=0))), tab


def delta_eff(hb, at, idx, rho, min_at=5):
    ha = at - hb
    m = idx[at[idx] >= min_at]
    if len(m) < 3:
        return np.nan
    n = at[m]
    z2 = (hb[m] - ha[m]) ** 2 / n
    e = (z2 - (1 + (n - 1) * rho)) / (4 * n)
    return float(np.sqrt(max(0.0, e.mean())))


def allele_loglik(hb, at, rho, deltas_by_state):
    """(T, NS) log emission. AT<1 -> 0 (marginalised: no evidence either way)."""
    obs = at >= 1
    k, n = hb[obs], at[obs]
    out = np.zeros((len(hb), NS))
    cache = {}

    def bb(p):
        p = float(np.clip(p, 1e-3, 1 - 1e-3))
        if p not in cache:
            cache[p] = betabinom.logpmf(k, n, p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho)
        return cache[p]
    for s in range(NS):
        if not S_IMB[s]:
            ll = bb(0.5)
        else:
            dl = deltas_by_state[s]
            ll = np.logaddexp(bb(0.5 + dl), bb(0.5 - dl)) - np.log(2)
        out[obs, s] = ll
    return out


def pi_vec():
    pi = np.full(NS, 0.40 / (NS - 2))
    pi[np.flatnonzero((S_R == 1) & (S_IMB == 0))] = PI_BAL_R1
    pi[np.flatnonzero((S_R == 1) & (S_IMB == 1))] = PI_IMB_R1
    return pi / pi.sum()


def lse(a):
    m = a.max()
    return m + np.log(np.exp(a - m).sum())


def fwd_bwd_viterbi(le, psw, logpi):
    """Forward-backward and Viterbi in log space for A_ij = (1-p)1[i=j] + p pi_j (p depends on the gap before t)."""
    T = le.shape[0]
    lp, lq = np.log(psw), np.log1p(-psw)
    la = np.empty((T, NS))
    la[0] = logpi + le[0]
    for t in range(1, T):
        L = lse(la[t - 1])
        la[t] = le[t] + np.logaddexp(lq[t] + la[t - 1], lp[t] + logpi + L)
    lb = np.zeros((T, NS))
    for t in range(T - 2, -1, -1):
        w = le[t + 1] + lb[t + 1]
        lb[t] = np.logaddexp(lq[t + 1] + w, lp[t + 1] + lse(logpi + w))
    lg = la + lb
    lg -= np.apply_along_axis(lse, 1, lg)[:, None]
    post = np.exp(lg)
    # Viterbi
    dl = np.empty((T, NS))
    ptr = np.zeros((T, NS), dtype=np.int64)
    dl[0] = logpi + le[0]
    stay_c = None
    for t in range(1, T):
        stay = np.log((1 - psw[t]) + psw[t] * np.exp(logpi)) + dl[t - 1]
        im = int(np.argmax(dl[t - 1]))
        sw = lp[t] + logpi + dl[t - 1][im]
        use_stay = stay >= sw
        dl[t] = le[t] + np.where(use_stay, stay, sw)
        ptr[t] = np.where(use_stay, np.arange(NS), im)
    path = np.empty(T, dtype=np.int64)
    path[-1] = int(np.argmax(dl[-1]))
    for t in range(T - 1, 0, -1):
        path[t - 1] = ptr[t, path[t]]
    return post, path, float(lse(la[-1]))


def run_hmm(d, spec_, rho, deltas_by_state, phi_total, base_rate, E, allele_on, logpi):
    """Run per chromosome. Returns posterior (n_bins, NS), path (n_bins,), log-likelihood."""
    n = d.n_bins
    post = np.zeros((n, NS))
    path = np.zeros(n, dtype=np.int64)
    ll_tot = 0.0
    y, hb, at = spec_["y"], spec_["hb"], spec_["at"]
    for ch in pd.unique(d.seq):
        ix = np.flatnonzero(d.seq == ch)
        mu = E * base_rate[ix][:, None] * R_GRID[S_DEPTH][None, :]
        le = nb_logpmf(y[ix][:, None], mu, phi_total)
        if allele_on:
            le = le + allele_loglik(hb[ix], at[ix], rho, deltas_by_state)
        dist = np.r_[1e5, np.diff(d.start[ix])].astype(float)
        psw = np.minimum(1 - np.exp(-HAZARD * dist / 1e5), P_SW_CAP)
        psw[0] = 0.5
        p_, v_, ll_ = fwd_bwd_viterbi(le, psw, logpi)
        post[ix], path[ix] = p_, v_
        ll_tot += ll_
    return post, path, ll_tot


# ----------------------------------------------------------------------------- summaries
def region_masks(d, regs):
    out = {}
    for r in regs.itertuples():
        out[r.region_id.replace("__target", "")] = (r.region_id.replace("__target", ""), "saved_window",
                                                     (d.seq == r.chrom) & (d.start >= r.first_retained_start) & (d.end <= r.last_retained_end))
    for name, ch, lo, hi in LIT:
        out[name] = (name, "literature", (d.seq == ch) & (d.start >= lo + 1) & (d.end <= hi))
    return out


def summarise(post, path, mask, at, allele_on):
    pm = post[mask].mean(0)
    n_at = int((at[mask] >= 1).sum())
    mx = float(pm.max())
    dep = np.array([pm[S_DEPTH == i].sum() for i in range(len(R_GRID))])
    reasons = []
    if mx < POSTERIOR_FLAG:
        reasons.append("max_state_prob<0.8")
    if n_at == 0:
        reasons.append("no_allele_obs")
    if not allele_on:
        reasons.append("allele_layer_off")
    vr = S_R[path[mask]]
    vals, cnt = np.unique(vr, return_counts=True)
    row = {"n_bins": int(mask.sum()), "bins_AT_ge1": n_at, "mass_depth_loss": float(pm[LOSS].sum()),
           "mass_depth_gain": float(pm[GAIN].sum()), "mass_imbalance": float(pm[S_IMB == 1].sum()) if allele_on else np.nan,
           "max_joint_state_prob": mx, "max_joint_state": S_NAME[int(pm.argmax())], "max_depth_prob": float(dep.max()),
           "map_depth_r": float(R_GRID[dep.argmax()]), "viterbi_modal_r": float(vals[cnt.argmax()]),
           "viterbi_modal_r_frac": float(cnt.max() / cnt.sum()), "insufficient_evidence": bool(("max_state_prob<0.8" in reasons) or ("no_allele_obs" in reasons)),
           "flag_reasons": ";".join(reasons)}
    for i, r in enumerate(R_GRID):
        row[f"p_r{r:g}"] = float(dep[i])
    return row


def viterbi_runs(path, d):
    runs = []
    for ch in pd.unique(d.seq):
        ix = np.flatnonzero(d.seq == ch)
        p = path[ix]
        cut = np.r_[0, np.flatnonzero(np.diff(p) != 0) + 1, len(p)]
        for a, b in zip(cut[:-1], cut[1:]):
            runs.append((ch, ix[a], ix[b - 1], a, b - 1, int(p[a])))
    return runs


def class_of(r):
    return np.where(r < 1, "loss", np.where(r > 1, "gain", "neutral"))


# ----------------------------------------------------------------------------- main
def main():
    t0 = time.time()
    pth = C.default_paths()
    d = C.H5AD(pth["h5ad"])
    lab = C.cell_clusters(d, pth)
    masks, ys, tot, e2 = load_pseudobulks(d, lab)
    allele, align_checks = load_allelic(d, pth, masks)
    kseg = pd.read_csv(pth["genome_pool"] / "derived" / "karyogram_cluster_segments.tsv", sep="\t")
    adbp = pd.read_csv(pth["genome_pool"] / "raw_accepted_breakpoints.tsv", sep="\t")
    regs = C.shared_regions()
    regs = regs[regs.part == "target"]
    rmasks = region_masks(d, regs)
    decile = pd.qcut(d.gc, 10, labels=False, duplicates="drop").astype(int)
    logpi = np.log(pi_vec())
    chr_pos = {ch: np.flatnonzero(d.seq == ch) for ch in pd.unique(d.seq)}
    local_idx = np.zeros(d.n_bins, int)
    for ch, ix in chr_pos.items():
        local_idx[ix] = np.arange(len(ix))

    params, post_rows, seg_rows, cmp_rows, sens_rows, tabs = {}, [], [], [], [], {}
    results = {}
    for g in GROUPS:
        y = ys[g]
        E = float(y.sum())
        rate = y / E
        call = np.ones(d.n_bins)
        for r in kseg[kseg.group == g].itertuples():
            call[r.bin_start:r.bin_end + 1] = r.cn_call
        # GC-decile baseline (estimated once per pseudobulk). Primary: median rate per GC decile over bins pyEpi calls base
        # (cn_call == 1; circular, stated). Variant baseline_all_bins: median over all bins.
        base_bins = call == 1
        base_dec_all = np.array([np.median(rate[decile == k]) for k in range(10)])
        base_dec_cn1 = np.array([np.median(rate[(decile == k) & base_bins]) for k in range(10)])
        # Mean-matching scalar c: the NB mean parameterisation needs g to track the MEAN level of base bins; with large
        # between-bin dispersion the median of rate/g sits below its mean (right skew), so the pure decile median would
        # make r=1 sit at the median and push most base bins to r=1.25. g = c * decile median, c = 5%-trimmed mean of
        # rate/median over cn_call==1 bins (one global scalar per pseudobulk). baseline_median_only uses c = 1.
        cs = {}
        base_rates = {"cn1_median": base_dec_cn1[decile]}
        for bk, bd in (("cn1", base_dec_cn1), ("all", base_dec_all)):
            cs[bk] = float(trim_mean(rate[base_bins] / bd[decile][base_bins], 0.05))
            base_rates[bk] = bd[decile] * cs[bk]
        cs["cn1_median"] = 1.0
        base_dec = base_dec_cn1 * cs["cn1"]
        # sampling dispersion from the stratified half split
        A, B = ys[g + "_A"], ys[g + "_B"]
        EA, EB = A.sum(), B.sum()
        v = (A / EA - B / EB) ** 2 - rate * (1 / EA + 1 / EB)
        per_bin = v / rate ** 2
        phi_pair = float(trim_mean(per_bin[rate > 0], 0.05))       # = phi_A + phi_B for the two half pools
        phi_A = max(phi_pair / 2, 0.0)
        # phi of a pool of cells = phi_cell * sum e_i^2 / (sum e_i)^2 -> rescale half -> full
        fac = (e2[g] / E ** 2) / (e2[g + "_A"] / EA ** 2)
        phi_s = max(phi_A * fac, 1e-8)   # floor: negative half-split estimate -> sampling-only NB is effectively Poisson
        # between-bin baseline variance from bins in uncalled (cn_call == 1) regions, per baseline choice
        tau2s, rob_vars = {}, {}
        for bk, br in base_rates.items():
            q = rate[base_bins] / br[base_bins]
            rob_vars[bk] = (1.4826 * np.median(np.abs(q - np.median(q)))) ** 2
            tau2s[bk] = max(0.0, float(rob_vars[bk] - (np.median(1 / (E * br[base_bins])) + phi_s)))
        base_rate = base_rates["cn1"]
        rob_var, tau2 = rob_vars["cn1"], tau2s["cn1"]
        samp_part = float(np.median(1 / (E * base_rate[base_bins])) + phi_s)
        # allele parameters
        hb, at = allele[("HB", g)], allele[("AT", g)]
        rho, rho_tab = estimate_rho(hb, at, call)
        tabs[f"rho_{g}"] = rho_tab
        # attenuation factor from the two published gain loci (MYC, 5q13-qter), per source
        att = {}
        for name, ch, lo, hi in LIT:
            if name in ("lit_MYC_8q24", "lit_5q13_qter"):
                m = np.flatnonzero(rmasks[name][2])
                r_obs = float(y[m].sum() / (E * base_rate[m]).sum())
                de = delta_eff(hb, at, m, rho)
                att[name] = {"r_obs": r_obs, "delta_eff": de, "delta_nominal_from_r_obs": float(nominal_delta(r_obs)),
                             "ratio": de / float(nominal_delta(r_obs))}
        a_fac = float(np.clip(np.mean([v_["ratio"] for v_ in att.values()]), 0.05, 1.0))
        deltas = {s: float(a_fac * nominal_delta(S_R[s])) for s in range(NS) if S_IMB[s]}
        lit9 = np.flatnonzero(rmasks["lit_9p21_hemizygous_zone_approx"][2])
        params[g] = {"cells": int(masks[g].sum()), "E_total_fragments": E, "gc_decile_baseline_rate": base_dec.tolist(),
                     "phi_pair_trimmed_mean": phi_pair, "phi_sampling_full_pool": phi_s, "cell_total_scaling_half_to_full": float(fac),
                     "robust_var_ratio_in_cn1_bins": float(rob_var), "sampling_part_of_that": samp_part, "tau2_between_bin": tau2,
                     "n_cn1_bins": int(base_bins.sum()), "tau2_by_baseline": tau2s, "mean_matching_scalar_c": cs, "gc_decile_median_rate_cn1_bins": base_dec_cn1.tolist(), "gc_decile_median_rate_all_bins": base_dec_all.tolist(), "phi_total_full": phi_s + tau2, "phi_total_no_baseline_var": phi_s,
                     "rho": rho, "attenuation_loci": att, "attenuation_factor_a": a_fac,
                     "delta_per_imbalanced_state": {S_NAME[s]: dv for s, dv in deltas.items()},
                     "unused_context_delta_eff_9p21_zone": delta_eff(hb, at, lit9, rho)}
        print(g, {k: v_ for k, v_ in params[g].items() if k not in ("gc_decile_baseline_rate", "delta_per_imbalanced_state", "attenuation_loci")}, flush=True)
        print(g, "att", att, flush=True)

        spec_ = {"y": y, "hb": hb, "at": at}
        for vn, vs in VARIANTS.items():
            phi_tot = (phi_s + tau2s[vs["base"]]) if vs["tau"] else phi_s
            post, path, ll = run_hmm(d, spec_, rho, deltas, phi_tot, base_rates[vs["base"]], E, vs["allele"], logpi)
            results[(g, vn)] = (post, path)
            params[g].setdefault("loglik", {})[vn] = ll
            for key, (nm, fam, m) in rmasks.items():
                row = {"source": g, "variant": vn, "region": nm, "family": fam,
                       "chrom": d.seq[m][0], "start_bp": int(d.start[m].min()), "end_bp": int(d.end[m].max())}
                row.update(summarise(post, path, m, at, vs["allele"]))
                row["pyepi_cn_call_modal"] = float(pd.Series(call[m]).mode().iloc[0])
                post_rows.append(row)
            # segments (Viterbi runs) for full variant only
            if vn == "full":
                for ch, i0, i1, l0, l1, s in viterbi_runs(path, d):
                    mm = np.zeros(d.n_bins, bool); mm[i0:i1 + 1] = True
                    sm = summarise(post, path, mm, at, True)
                    seg_rows.append({"source": g, "chrom": ch, "bin_start_in_chrom": l0, "bin_end_in_chrom": l1,
                                     "start_bp": int(d.start[i0]), "end_bp": int(d.end[i1]), "n_bins": i1 - i0 + 1,
                                     "state": S_NAME[s], "rel_depth_r": float(S_R[s]), "allele_state": "imbalanced" if S_IMB[s] else "balanced",
                                     "delta_if_imbalanced": deltas.get(s, np.nan),
                                     "mean_post_of_state": float(post[mm, s].mean()), "mass_depth_loss": sm["mass_depth_loss"],
                                     "mass_depth_gain": sm["mass_depth_gain"], "mass_imbalance": sm["mass_imbalance"],
                                     "max_joint_state_prob": sm["max_joint_state_prob"], "bins_AT_ge1": sm["bins_AT_ge1"],
                                     "insufficient_evidence": sm["insufficient_evidence"],
                                     "modal_pyepi_cn_call": float(pd.Series(call[mm]).mode().iloc[0])})
            # comparison: boundaries vs AD and confusion vs cn_call
            hcls = class_of(S_R[path])
            ccls = np.where(call < 1, "loss", np.where(call > 1, "gain", "neutral"))
            for hc in ("loss", "neutral", "gain"):
                for cc in ("loss", "neutral", "gain"):
                    cmp_rows.append({"table": "confusion_bins_hmm_class_vs_pyepi_cn_call", "source": g, "variant": vn,
                                     "key1": f"hmm_{hc}", "key2": f"pyepi_{cc}", "n": int(((hcls == hc) & (ccls == cc)).sum()), "value": np.nan, "note": ""})
            if vn == "full":
                for rv in R_GRID:
                    for cv in np.unique(call):
                        n_ = int(((S_R[path] == rv) & (call == cv)).sum())
                        if n_:
                            cmp_rows.append({"table": "confusion_bins_hmm_r_vs_pyepi_cn_call_value", "source": g, "variant": vn,
                                             "key1": f"hmm_r={rv:g}", "key2": f"pyepi_cn_call={cv:g}", "n": n_, "value": np.nan, "note": ""})
            # boundary agreement (depth-level boundaries only)
            hb_b, ad_b, L_tot = {}, {}, {}
            for ch, ix in chr_pos.items():
                dp = S_DEPTH[path[ix]]
                hb_b[ch] = np.flatnonzero(np.diff(dp) != 0) + 1
                ad_b[ch] = np.sort(adbp[(adbp.source == g) & (adbp.chromosome == ch)].absolute_bin.astype(int).to_numpy())
                L_tot[ch] = len(ix)
            for k in (2, 5):
                n_ad = sum(len(v_) for v_ in ad_b.values())
                n_h = sum(len(v_) for v_ in hb_b.values())
                ad_hit = sum(int(any(abs(h - a) <= k for h in hb_b[ch])) for ch in ad_b for a in ad_b[ch])
                h_hit = sum(int(any(abs(h - a) <= k for a in ad_b[ch])) for ch in hb_b for h in hb_b[ch])
                chance = float(np.mean([1 - (1 - (2 * k + 1) / L_tot[ch]) ** len(hb_b[ch]) for ch in ad_b for _ in ad_b[ch]])) if n_ad else np.nan
                cmp_rows.append({"table": f"boundary_match_pm{k}bins", "source": g, "variant": vn, "key1": "AD_breakpoints_with_hmm_boundary",
                                 "key2": f"{ad_hit}/{n_ad}", "n": n_ad, "value": ad_hit / max(n_ad, 1),
                                 "note": f"chance-level reference if hmm boundaries were placed at random per chromosome: {chance:.3f}"})
                cmp_rows.append({"table": f"boundary_match_pm{k}bins", "source": g, "variant": vn, "key1": "hmm_depth_boundaries_with_AD_breakpoint",
                                 "key2": f"{h_hit}/{n_h}", "n": n_h, "value": h_hit / max(n_h, 1), "note": "AD = raw_accepted_breakpoints.tsv; neither is truth"})
            # diagnostic: bins in pyEpi cn_call==1 called imbalanced / non-neutral by HMM
            cmp_rows.append({"table": "diagnostic_cn1_bins", "source": g, "variant": vn, "key1": "frac_cn1_bins_viterbi_imbalanced",
                             "key2": "", "n": int((call == 1).sum()), "value": float(S_IMB[path[call == 1]].mean()), "note": "pyEpi cn_call==1 is not truth"})
            cmp_rows.append({"table": "diagnostic_cn1_bins", "source": g, "variant": vn, "key1": "frac_cn1_bins_viterbi_depth_not_r1",
                             "key2": "", "n": int((call == 1).sum()), "value": float((S_R[path[call == 1]] != 1).mean()), "note": ""})
            print(g, vn, "done", round(time.time() - t0, 1), "s", flush=True)

    posts = pd.DataFrame(post_rows)
    segs = pd.DataFrame(seg_rows)
    cmp_df = pd.DataFrame(cmp_rows)
    # ---- sensitivity: compare each variant with full
    key = ["source", "region"]
    full = posts[posts.variant == "full"].set_index(key)
    for vn in ("no_baseline_var", "no_allele", "baseline_all_bins", "baseline_median_only"):
        alt = posts[posts.variant == vn].set_index(key)
        for idx in full.index:
            f, a = full.loc[idx], alt.loc[idx]
            sens_rows.append({"source": idx[0], "region": idx[1], "family": f.family, "variant": vn,
                              "full_mass_loss": f.mass_depth_loss, "alt_mass_loss": a.mass_depth_loss,
                              "full_mass_gain": f.mass_depth_gain, "alt_mass_gain": a.mass_depth_gain,
                              "full_mass_imb": f.mass_imbalance, "alt_mass_imb": a.mass_imbalance,
                              "full_map_depth_r": f.map_depth_r, "alt_map_depth_r": a.map_depth_r,
                              "full_max_joint_prob": f.max_joint_state_prob, "alt_max_joint_prob": a.max_joint_state_prob,
                              "full_insufficient": f.insufficient_evidence, "alt_insufficient": a.insufficient_evidence,
                              "map_depth_changed": bool(f.map_depth_r != a.map_depth_r),
                              "loss_gain_class_changed": bool((f.mass_depth_loss > .5) != (a.mass_depth_loss > .5) or (f.mass_depth_gain > .5) != (a.mass_depth_gain > .5))})
    for g in GROUPS:
        pf = results[(g, "full")][1]
        for vn in ("no_baseline_var", "no_allele", "baseline_all_bins", "baseline_median_only"):
            pv = results[(g, vn)][1]
            sens_rows.append({"source": g, "region": "GENOME_ALL_AUTOSOME_BINS", "family": "genome", "variant": vn,
                              "full_mass_loss": np.nan, "alt_mass_loss": np.nan,
                              "full_map_depth_r": np.nan, "alt_map_depth_r": np.nan,
                              "bins_depth_class_changed_frac": float((class_of(S_R[pf]) != class_of(S_R[pv])).mean()),
                              "bins_depth_r_changed_frac": float((S_R[pf] != S_R[pv]).mean()),
                              "depth_runs_full": len({(c, i0) for c, i0, i1, l0, l1, s in viterbi_runs(S_DEPTH[pf], d)}),
                              "depth_runs_alt": len({(c, i0) for c, i0, i1, l0, l1, s in viterbi_runs(S_DEPTH[pv], d)})})
    sens = pd.DataFrame(sens_rows)

    posts.to_csv(HERE / "posteriors_regions.tsv", sep="\t", index=False)
    segs.to_csv(HERE / "segments.tsv", sep="\t", index=False)
    cmp_df.to_csv(HERE / "comparison.tsv", sep="\t", index=False)
    sens.to_csv(HERE / "sensitivity.tsv", sep="\t", index=False)

    # ---- figures
    figs = {}
    sv = posts[(posts.variant == "full") & (posts.family == "saved_window")]
    for g in GROUPS:
        post, path = results[(g, "full")]
        y = ys[g]; E = y.sum(); rate = y / E
        base_rate = np.array(params[g]["gc_decile_baseline_rate"])[decile]
        fig, axs = plt.subplots(3, 2, figsize=(13, 7.5), sharex="col")
        for j, ch in enumerate(["chr7", "chr9"]):
            ix = chr_pos[ch]; x = d.start[ix] / 1e6
            ax = axs[0, j]
            ax.scatter(x, np.log2(rate[ix] / base_rate[ix]), s=3, color="#999999", label="raw rate / GC-decile baseline")
            ax.step(x, np.log2(S_R[path[ix]]), where="mid", color="#0072B2", lw=1.6, label="Viterbi relative depth r")
            ax.set_ylim(-2.3, 1.6); ax.set_ylabel("log2 ratio"); ax.set_title(f"{g}: {ch}"); ax.legend(fontsize=7, loc="lower right")
            axs[1, j].plot(x, post[ix][:, LOSS].sum(1), color="#D55E00", lw=1.2, label="posterior mass, depth loss (r<1)")
            axs[1, j].plot(x, post[ix][:, GAIN].sum(1), color="#009E73", lw=1.2, label="posterior mass, depth gain (r>1)")
            axs[1, j].set_ylim(-0.03, 1.03); axs[1, j].legend(fontsize=7, loc="center right"); axs[1, j].set_ylabel("posterior")
            axs[2, j].plot(x, post[ix][:, S_IMB == 1].sum(1), color="#CC79A7", lw=1.2, label="posterior mass, allelic imbalance")
            axs[2, j].plot(x, (at_ := allele[("AT", g)][ix]) / max(at_.max(), 1), color="#bbbbbb", lw=.7, label="AT / max AT (observability)")
            axs[2, j].set_ylim(-0.03, 1.03); axs[2, j].legend(fontsize=7, loc="center right"); axs[2, j].set_xlabel(f"{ch} position (Mb)")
            for _, r in sv[(sv.source == g) & (sv.chrom == ch)].iterrows():
                for a in axs[:, j]:
                    a.axvspan(r.start_bp / 1e6, r.end_bp / 1e6, color="#E69F00", alpha=.15, lw=0)
        fig.tight_layout()
        figs[g] = C.fig_to_b64(fig); plt.close(fig)

    fig, axs = plt.subplots(2, 1, figsize=(13, 5.5))
    for j, ch in enumerate(["chr7", "chr9"]):
        ix = chr_pos[ch]; x = d.start[ix] / 1e6
        for vn, col, lw in (("full", "#0072B2", 1.8), ("no_baseline_var", "#D55E00", 0.8), ("no_allele", "#009E73", 1.2)):
            axs[j].step(x, np.log2(S_R[results[("all", vn)][1][ix]]), where="mid", color=col, lw=lw, label=vn, alpha=.9)
        axs[j].set_ylabel("log2 Viterbi r"); axs[j].set_title(f"all cells: {ch}"); axs[j].legend(fontsize=7, loc="lower right")
    axs[1].set_xlabel("position (Mb)")
    fig.tight_layout(); figs["sens"] = C.fig_to_b64(fig); plt.close(fig)

    summary = {"experiment": "13_joint_hmm", "status": "exploratory diagnostic; relative depth states are not copy number",
               "command": "/home/user/envs/cnv310/bin/python -I research/2026-11-04-cloud-method-review/experiments/13_joint_hmm/run.py",
               "seeds": {"stratified_halves": list(SEED_HALVES), "hmm": "deterministic"},
               "settings": {"R_GRID": R_GRID.tolist(), "hazard_per_100kb": HAZARD, "switch_prob_cap": P_SW_CAP,
                            "pi_balanced_r1": PI_BAL_R1, "pi_imbalanced_r1": PI_IMB_R1, "pi_other_each": 0.40 / (NS - 2),
                            "copy_neutral_imbalance_nominal_delta": CN_NEUTRAL_DELTA_NOMINAL, "posterior_flag_threshold": POSTERIOR_FLAG,
                            "min_AT_for_rho": MIN_AT_RHO, "variants": VARIANTS, "state_names": S_NAME},
               "alignment_checks": align_checks, "parameters": params,
               "seconds": round(time.time() - t0, 1), "environment": C.environment()}
    tabsd = {k: v.to_dict("records") for k, v in tabs.items()}
    summary["rho_tables"] = tabsd
    C.write_json(HERE / "summary.json", summary)
    write_report(posts, segs, cmp_df, sens, params, figs, summary)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    cols = ["source", "variant", "region", "n_bins", "bins_AT_ge1", "mass_depth_loss", "mass_depth_gain", "mass_imbalance",
            "max_joint_state_prob", "max_joint_state", "map_depth_r", "insufficient_evidence", "pyepi_cn_call_modal"]
    print(posts[cols].to_string())
    print(cmp_df[~cmp_df.table.str.startswith("confusion_bins_hmm_r")].to_string())
    print(sens[sens.family == "genome"].to_string())
    print("seconds", round(time.time() - t0, 1))


def reason(r):
    """Plain-language reason for a region row."""
    bits = []
    if r.bins_AT_ge1 == 0:
        bits.append("no allele observations")
    top = f"largest single joint state {r.max_joint_state} has posterior {r.max_joint_state_prob:.2f}"
    if r.max_joint_state_prob < POSTERIOR_FLAG:
        bits.append(f"{top} (<0.8): posterior is spread over several depth levels or allele states")
    else:
        bits.append(top)
    if r.mass_depth_loss > 0.8:
        bits.append("depth evidence favours a loss-level state")
    elif r.mass_depth_gain > 0.8:
        bits.append("depth evidence favours a gain-level state")
    elif r.mass_depth_loss < 0.2 and r.mass_depth_gain < 0.2:
        bits.append("depth evidence favours the baseline level")
    else:
        bits.append("depth direction is not settled")
    if r.mass_imbalance > 0.8:
        bits.append("allele layer favours imbalance")
    elif r.mass_imbalance < 0.2:
        bits.append("allele layer favours balance (where observed)")
    return "; ".join(bits)


def write_report(posts, segs, cmp_df, sens, params, figs, summary):
    f = posts[posts.variant == "full"].copy()
    f["reason"] = f.apply(reason, axis=1)
    pa, pc = params["all"], params["cluster0"]
    sv = f[f.family == "saved_window"]
    lit = f[f.family == "literature"]

    def g(src, reg, col, var="full"):
        return float(posts[(posts.source == src) & (posts.region == reg) & (posts.variant == var)][col].iloc[0])
    chr7 = ["chr7_8.3_10.7", "chr7_51.4_53.5", "chr7_67.4_70.0", "chr7_70.0_73.3"]
    n_ins = int(f.insufficient_evidence.sum())

    def cmpv(tbl, src, key1, var="full"):
        r = cmp_df[(cmp_df.table == tbl) & (cmp_df.source == src) & (cmp_df.key1 == key1) & (cmp_df.variant == var)]
        return r.iloc[0]
    ad2 = cmpv("boundary_match_pm2bins", "all", "AD_breakpoints_with_hmm_boundary")
    ad5 = cmpv("boundary_match_pm5bins", "all", "AD_breakpoints_with_hmm_boundary")
    h2 = cmpv("boundary_match_pm2bins", "all", "hmm_depth_boundaries_with_AD_breakpoint")
    conf = cmp_df[(cmp_df.table == "confusion_bins_hmm_class_vs_pyepi_cn_call") & (cmp_df.variant == "full")]

    def cf(src, a, b):
        return int(conf[(conf.source == src) & (conf.key1 == f"hmm_{a}") & (conf.key2 == f"pyepi_{b}")].n.iloc[0])
    nb = int(conf[conf.source == "all"].n.sum())
    agree_all = (cf("all", "loss", "loss") + cf("all", "neutral", "neutral") + cf("all", "gain", "gain")) / nb
    genome = sens[sens.family == "genome"].set_index(["source", "variant"])
    nseg_all = int(((segs.source == "all")).sum())
    ndepth_all = int(genome.loc[("all", "no_baseline_var")].depth_runs_full)
    ndepth_nb = int(genome.loc[("all", "no_baseline_var")].depth_runs_alt)
    cn1_imb = float(cmpv("diagnostic_cn1_bins", "all", "frac_cn1_bins_viterbi_imbalanced").value)
    cn1_notr1 = float(cmpv("diagnostic_cn1_bins", "all", "frac_cn1_bins_viterbi_depth_not_r1").value)
    s_reg = sens[sens.family != "genome"]
    na = s_reg[(s_reg.variant == "no_allele")]
    nb_ = s_reg[(s_reg.variant == "no_baseline_var")]
    ins_gain_noallele = int((na.alt_insufficient & ~na.full_insufficient).sum())
    ins_lost_noallele = int((~na.alt_insufficient & na.full_insufficient).sum())

    loss7 = ", ".join("%.2f" % g("all", r, "mass_depth_loss") for r in chr7)
    imb7 = ", ".join("%.2f" % g("all", r, "mass_imbalance") for r in chr7)
    verdict = ("A transparent joint depth + allele HMM reproduces the large chr7 loss-level windows and the MYC gain-level locus as relative states, "
               "but only because a large between-bin baseline-variance term is included; the allele layer mainly changes confidence and the imbalance mass, "
               "and 9p21 and the weak-gain loci remain insufficient evidence. Exploratory; relative states are not copy number.")
    key_findings = [
        f"<b>Baseline variance dominates the depth emission.</b> Cell-sampling dispersion is tiny (phi_s = {pa['phi_sampling_full_pool']:.1e} all; the cluster0 half-split estimate was "
        f"negative and floored at {pc['phi_sampling_full_pool']:.0e}) but the between-bin baseline variance in cn_call==1 bins is tau2 = {pa['tau2_between_bin']:.3f} "
        f"(about {np.sqrt(pa['tau2_between_bin']):.2f} relative SD per bin). Without it the HMM fragments the genome: {ndepth_nb} depth runs vs {ndepth_all} in the full model (all-cell pseudobulk), "
        f"{genome.loc[('all','no_baseline_var')].bins_depth_class_changed_frac:.0%} of bins change loss/neutral/gain class, and {int(nb_.groupby('source').loss_gain_class_changed.sum().sum())} of {len(nb_)} region summaries flip class.",
        f"<b>chr7 saved windows (all cells):</b> posterior mass on depth-loss states {loss7} "
        f"(8.3-10.7, 51.4-53.5, 67.4-70.0, 70.0-73.3 Mb); imbalance mass {imb7}. "
        f"The loss-level depth call is consistent with pyEpi cn_call 0 there; the allele layer shows no imbalance excess at the observed sites.",
        f"<b>chr9 windows are unresolved:</b> loss mass 9p 20-27 Mb {g('all','chr9_20_27','mass_depth_loss'):.2f}, 22.5-26.6 Mb {g('all','chr9_22.5_26.6','mass_depth_loss'):.2f}, "
        f"9p21 literature zone {g('all','lit_9p21_hemizygous_zone_approx','mass_depth_loss'):.2f}, all flagged insufficient evidence (max joint-state probability "
        f"{g('all','chr9_20_27','max_joint_state_prob'):.2f}/{g('all','chr9_22.5_26.6','max_joint_state_prob'):.2f}/{g('all','lit_9p21_hemizygous_zone_approx','max_joint_state_prob'):.2f}); "
        f"only {int(g('all','lit_9p21_hemizygous_zone_approx','bins_AT_ge1'))} of {int(g('all','lit_9p21_hemizygous_zone_approx','n_bins'))} bins in the 9p21 zone carry any allele observation.",
        f"<b>Gain controls:</b> MYC gain mass {g('all','lit_MYC_8q24','mass_depth_gain'):.2f} (imbalance mass {g('all','lit_MYC_8q24','mass_imbalance'):.2f}, max joint-state {g('all','lit_MYC_8q24','max_joint_state_prob'):.2f}); "
        f"5q13-qter gain mass {g('all','lit_5q13_qter','mass_depth_gain'):.2f} (imbalance {g('all','lit_5q13_qter','mass_imbalance'):.2f}, flagged insufficient); "
        f"19q13-qter gain mass {g('all','lit_19q13_qter','mass_depth_gain'):.2f} vs loss {g('all','lit_19q13_qter','mass_depth_loss'):.2f} (unresolved). "
        f"pyEpi cn_call is 2 at MYC and 1 at 5q and 19q.",
        f"<b>Agreement with AD and pyEpi (all cells; neither is truth):</b> {int(ad2.key2.split('/')[0])}/{int(ad2.n)} AD breakpoints ({ad2.value:.0%}) have an HMM depth boundary within 2 bins "
        f"(random-placement reference about {float(ad2.note.split(': ')[1]):.0%}), {ad5.value:.0%} within 5 bins; but only {h2.value:.0%} of the {int(h2.n)} HMM depth boundaries lie within 2 bins of an AD breakpoint. "
        f"Bin-level class agreement with pyEpi cn_call is {agree_all:.0%}; the main disagreement is {cf('all','gain','neutral')} bins the HMM places at gain level (r > 1) where pyEpi says neutral.",
        f"<b>Allele layer in sensitivity:</b> switching it off changes {genome.loc[('all','no_allele')].bins_depth_class_changed_frac:.0%} (all) and {genome.loc[('cluster0','no_allele')].bins_depth_class_changed_frac:.0%} (cluster0) of Viterbi depth classes, "
        f"changes the modal depth level in {int(na.map_depth_changed.sum())} of {len(na)} region summaries, and removes the ability to flag imbalance; it lowers the max joint-state probability "
        f"(insufficient-evidence flags rise in {ins_gain_noallele} region summaries, fall in {ins_lost_noallele}).",
        f"<b>Calibration warning:</b> {cn1_imb:.1%} of pyEpi cn_call==1 bins are Viterbi-imbalanced and {cn1_notr1:.0%} have a depth state other than r = 1 (all cells), so state calls in single bins are not a calibrated error rate. "
        f"{n_ins} of {len(f)} region summaries (full model, both sources) are flagged insufficient evidence.",
    ]
    sec = []
    sec.append(("Question and why it matters",
                "<p>Method record <b>B:3</b> (covariance-aware joint segmentation / Bayesian HMM, proposed) argues that depth and allele counts should be modelled together, with missing allele data treated as "
                "no evidence, distance-aware transitions and abstention. This experiment asks, on the actual MCF10-A development data: (i) what a simple, fully inspectable version of that model says in the six saved windows and four "
                "literature loci; (ii) how its relative states and boundaries compare with the AD breakpoints and pyEpi <code>cn_call</code>; and (iii) what the allele layer and the between-bin baseline-variance "
                "term each contribute.</p><p><b>Decision informed:</b> whether the review's recommendation should add a joint model (and which component would carry its value) or keep depth-only segmentation with a separate allelic evidence label.</p>"))
    sec.append(("Inputs",
                "<ul><li><b>Depth:</b> <code>pyepi_results/count_matrix.h5ad</code> layer <code>layers/raw</code> (3,460 cells x 25,337 retained 100 kb autosomal bins, chr1-22; raw fragment counts). Pseudobulks: all 3,460 cells and WNN cluster 0 "
                "(1,286 cells, <code>wnn_leiden_0.3</code> from <code>integrated.h5mu</code>). Exposure E = pooled raw counts of the pseudobulk (2.01e8 all; 5.12e7 cluster0). GC from <code>var/GC</code>.</li>"
                "<li><b>Allele:</b> <code>aligned_matrices/HB.parquet</code> and <code>AT.parquet</code> (bins x barcodes, joined by barcode name; coordinates checked against the h5ad; float32 integer-valued). AT = summed SNP-site depth, "
                "not unique fragments. HA = AT - HB. Counts are pooled over the same cells as the depth track. Bins with AT &lt; 1 are marginalised (no evidence either way).</li>"
                "<li><b>Calls for comparison / parameter estimation:</b> <code>karyogram_cluster_segments.tsv</code> (pyEpi cn_call per bin, groups all and cluster0) and <code>raw_accepted_breakpoints.tsv</code> (AD breakpoints; chromosome-local bin index). "
                "Regions: <code>research/2026-10-29-low-coverage/shared/regions.tsv</code> rows with part == target (6 windows, bins within first_retained_start..last_retained_end) and four literature loci (chr9 20-24 Mb, chr8 127.0-128.5 Mb, chr19 &gt;= 35 Mb, chr5 &gt;= 70 Mb).</li>"
                "<li><b>Not available:</b> phase-set / SNP-level phase information (imbalance is phase-agnostic), unique-molecule allele counts, a matched normal, per-cell state calls, ground-truth copy number.</li></ul>"))
    model = """<pre>
Observations per retained bin b, chromosome by chromosome (pseudobulk of a set of cells):
  y_b = pooled raw fragment count, E = sum_b y_b;   k_b = HB_b, n_b = AT_b pooled over the same cells.
Hidden joint state s = (r, a):  r in {0.25, 0.5, 0.75, 1, 1.25, 1.5, 2} (relative depth x baseline),  a in {balanced, imbalanced}   -> 14 states
  loss states: r &lt; 1;  gain states: r &gt; 1;  (r = 1, imbalanced) = copy-neutral-imbalance state.

Depth emission     y_b | s ~ NB( mean mu_bs = E * g(GC_b) * r_s ,  var = mu + phi * mu^2 ),   phi = phi_s + tau2
  g(GC decile)  = c * median over cn_call==1 bins of (y/E) within the bin's GC decile   (estimated once per pseudobulk)
  c             = 5%-trimmed mean of (y/E)/median_dec over cn_call==1 bins (one scalar; makes r = 1 the MEAN level of base bins,
                  because NB means, not medians, are parameterised and between-bin dispersion is right-skewed)
  phi_s         = sampling dispersion: with halves A,B from C.stratified_halves(WNN labels, (20261104, 13)),
                  phi_A + phi_B = trimmed mean over bins of [ (y_A/E_A - y_B/E_B)^2 - (y/E)(1/E_A + 1/E_B) ] / (y/E)^2,
                  phi_s = phi_A * [sum_i e_i^2 / E^2]_full / [sum_i e_i^2 / E_A^2]_half   (e_i = cell totals), floored at 1e-8
  tau2          = max(0, (1.4826 MAD of (y/E)/g over cn_call==1 bins)^2 - median(1/mu) - phi_s)   (extra between-bin baseline variance)

Allele emission    k_b | balanced   ~ BetaBin(n_b, p = 0.5, rho)
                   k_b | imbalanced ~ 0.5 BetaBin(n_b, 0.5 + delta_s, rho) + 0.5 BetaBin(n_b, 0.5 - delta_s, rho)      (phase-agnostic)
                   BetaBin parameters a = p(1-rho)/rho, b = (1-p)(1-rho)/rho ; p clipped to [0.001, 0.999]
                   n_b &lt; 1: log emission = 0 for every state (marginalised; never read as balanced)
  rho     = exp-06 estimator: median over 5 AT strata of (median[(HB-HA)^2/AT]/0.4549 - 1)/(median AT - 1), cn_call==1 bins with AT &gt;= 20
  delta_s = a * |r-1| / (2r) (cap 0.5) for r != 1  [a clonal one-allele change giving relative depth r = copies/2]; a * 0.25 for r = 1 (copy-neutral imbalance, arbitrary)
  a       = mean over MYC and 5q13-qter of  delta_eff / ( |r_obs-1| / (2 r_obs) ),  r_obs = observed depth ratio in the locus, delta_eff as in exp 06

Transitions   P(s_t = j | s_{t-1} = i) = (1 - p_t) 1[i = j] + p_t * pi_j ,   p_t = min( 1 - exp(-0.01 * gap_bp / 1e5), 0.5 ),  gap_bp = start_t - start_{t-1}
  (adjacent 100 kb bins: p = 0.010; a 1 Mb QC gap: 0.095; a 7 Mb gap: 0.5)    first bin of a chromosome: p = 0.5 from pi
  pi = 0.55 on (r=1, balanced), 0.05 on (r=1, imbalanced), 0.40 shared equally by the other 12 states.   Hazard and pi are prespecified, not tuned.
Inference   log-space forward-backward (posterior per bin) and Viterbi (maximum path), independently per chromosome. Deterministic.
Region / segment summary: mean of bin posteriors; insufficient evidence = max joint-state probability of that mean &lt; 0.8, or no bin with AT &gt;= 1.
</pre>"""
    sec.append(("Method", model + "<p>Seeds: only the stratified half split is random, <code>np.random.default_rng((20261104, 13))</code>. The HMM is deterministic. "
                "Variants: <code>full</code>; <code>no_baseline_var</code> (tau2 = 0, sampling-only NB); <code>no_allele</code> (allele emission removed); extra checks <code>baseline_all_bins</code> "
                "(decile medians over all bins instead of cn_call==1 bins) and <code>baseline_median_only</code> (c = 1).</p>"))
    ptab = []
    for src, p in params.items():
        ptab.append({"source": src, "cells": p["cells"], "E": p["E_total_fragments"], "phi_s": p["phi_sampling_full_pool"], "tau2": p["tau2_between_bin"],
                     "c": p["mean_matching_scalar_c"]["cn1"], "rho": p["rho"], "a": p["attenuation_factor_a"],
                     "delta r=0.5 imb": p["delta_per_imbalanced_state"]["r0.5_imb"], "delta r=1 imb": p["delta_per_imbalanced_state"]["r1_imb"],
                     "delta r=1.5 imb": p["delta_per_imbalanced_state"]["r1.5_imb"], "delta r=2 imb": p["delta_per_imbalanced_state"]["r2_imb"]})
    att = [{"source": src, "locus": k, **v} for src, p in params.items() for k, v in p["attenuation_loci"].items()]
    res = []
    res.append(("Estimated parameters", C.html_table(pd.DataFrame(ptab), floatfmt="{:.4g}") + C.html_table(pd.DataFrame(att), floatfmt="{:.4g}") +
                C.caption("One row per pseudobulk. phi_s: cell-sampling dispersion (half split; cluster0 estimate was negative and floored). tau2: between-bin baseline variance. "
                          "a: empirical attenuation of allelic imbalance relative to nominal copy-based expectation, from MYC and 5q13-qter (second table); applied to losses and copy-neutral imbalance by extrapolation. "
                          "Full numbers, GC-decile baselines and log-likelihoods are in summary.json."), 3))
    show = ["source", "region", "n_bins", "bins_AT_ge1", "mass_depth_loss", "mass_depth_gain", "mass_imbalance", "max_joint_state_prob", "max_joint_state", "map_depth_r", "insufficient_evidence", "pyepi_cn_call_modal"]
    res.append(("Posterior mass in the saved windows and literature loci", C.html_table(f[f.family == "saved_window"][show]) + C.html_table(f[f.family == "literature"][show]) +
                C.caption("Mean posterior over the bins of each region (full model). Loss/gain mass = total posterior on r &lt; 1 / r &gt; 1; imbalance mass = total posterior on imbalanced states. "
                          "max_joint_state_prob below 0.8 triggers the insufficient-evidence flag. pyepi_cn_call_modal is the most common pyEpi call in the region (0.5 loss, 1 base, 2 gain levels; not truth). Table: posteriors_regions.tsv."), 3))
    res.append(("Chromosome profiles: chr7 and chr9", "".join(f"<p><b>{src}</b></p>" + C.img(b, f"HMM {src}") for src, b in figs.items() if src in ("all", "cluster0")) +
                C.caption("Top: grey points = raw fragment rate over the GC-decile baseline (log2); blue = Viterbi relative depth. Middle: posterior mass on depth-loss (red) and depth-gain (green) states. "
                          "Bottom: posterior mass on imbalanced states (pink) and AT relative to its chromosome maximum (grey; where it is near zero the allele layer is silent). Orange bands = the saved windows. "
                          "Relative depth, not copy number."), 3))
    segq = segs[segs.source == "all"]
    seg_summ = segq.groupby("state").agg(runs=("n_bins", "size"), bins=("n_bins", "sum"), median_run_bins=("n_bins", "median"), mean_max_joint_prob=("max_joint_state_prob", "mean")).reset_index()
    segs_by = segs.groupby("source").agg(runs=("n_bins", "size"), median_run_bins=("n_bins", "median"), insufficient_runs=("insufficient_evidence", "sum")).reset_index()
    res.append(("Segments (Viterbi runs)", C.html_table(segs_by) + C.html_table(seg_summ) +
                C.caption("Viterbi runs of identical joint state across all autosomes. Upper table per pseudobulk; lower table by state for the all-cell pseudobulk. Many short runs have max joint-state probability under 0.8, "
                          "i.e. the segment boundary or state is not well determined. Full list with per-segment posterior masses and flags: segments.tsv."), 3))
    cm = conf[conf.source == "all"][["key1", "key2", "n"]].pivot(index="key1", columns="key2", values="n").reset_index()
    cm2 = conf[conf.source == "cluster0"][["key1", "key2", "n"]].pivot(index="key1", columns="key2", values="n").reset_index()
    bm = cmp_df[(cmp_df.table.str.startswith("boundary")) & (cmp_df.variant == "full")][["table", "source", "key1", "key2", "value", "note"]]
    res.append(("Comparison with AD breakpoints and pyEpi cn_call", "<p>Bin-level class table (rows HMM, columns pyEpi), all cells then cluster0:</p>" + C.html_table(cm) + C.html_table(cm2) +
                "<p>Boundary agreement (depth-level boundaries):</p>" + C.html_table(bm) +
                C.caption("Class: loss = r &lt; 1 or cn_call &lt; 1, gain = r &gt; 1 or cn_call &gt; 1. Both are relative states from the same pseudobulk and neither is truth. The fraction of AD breakpoints with an HMM boundary nearby is far above a random-placement reference, "
                          "but most HMM depth boundaries have no AD breakpoint nearby, so the HMM is finer-grained than the AD scan. Finer table (r vs cn_call value) and all variants: comparison.tsv."), 3))
    sec.append(("Results", "", 2))
    for h, fr, lv in res:
        sec.append((h, fr, lv))
    det = f[["source", "region", "family", "mass_depth_loss", "mass_depth_gain", "mass_imbalance", "insufficient_evidence", "reason"]].copy()
    sec.append(("Per-region detail", C.html_table(det) + C.caption("One row per region and pseudobulk (full model). reason is generated from the numbers: it states the largest joint-state probability and which layer favours what. "
                                                                  "Insufficient evidence means the model could not commit; it is not evidence of balance or neutrality.")))
    sv_t = sens[sens.family != "genome"]
    cols = ["source", "region", "variant", "full_mass_loss", "alt_mass_loss", "full_mass_gain", "alt_mass_gain", "full_mass_imb", "alt_mass_imb", "full_map_depth_r", "alt_map_depth_r",
            "full_max_joint_prob", "alt_max_joint_prob", "full_insufficient", "alt_insufficient"]
    gen = sens[sens.family == "genome"][["source", "variant", "bins_depth_class_changed_frac", "bins_depth_r_changed_frac", "depth_runs_full", "depth_runs_alt"]]
    sens_html = ("<h3>Genome-wide effect of each change</h3>" + C.html_table(gen) +
                 C.caption("Fractions are of all autosomal bins whose Viterbi depth class (or level) differs from the full model; depth_runs counts Viterbi depth-level runs. "
                           "Removing the baseline-variance term multiplies the number of runs about 30-fold.") +
                 "<h3>Saved windows and literature loci, no baseline-variance term</h3>" + C.html_table(sv_t[sv_t.variant == "no_baseline_var"][cols]) +
                 "<h3>Saved windows and literature loci, allele emission switched off</h3>" + C.html_table(sv_t[sv_t.variant == "no_allele"][cols]) +
                 C.caption("full_* = full model; alt_* = variant. With the allele layer off, imbalance mass is undefined and the posterior cannot distinguish balanced from imbalanced states, so max joint-state probability is mechanically lower "
                           "(the 14 states collapse to 7 identifiable ones); compare depth masses and map_depth_r instead. Extra baseline checks (baseline_all_bins, baseline_median_only) are in sensitivity.tsv and posteriors_regions.tsv.") +
                 "<h3>Sensitivity figure</h3>" + C.img(figs["sens"], "sensitivity paths") +
                 C.caption("Viterbi log2 relative depth on chr7 and chr9 (all cells) for the full model, without the baseline-variance term and without the allele layer. The sampling-only model chases bin-to-bin noise; "
                           "dropping the allele layer changes a minority of Viterbi depth levels (see the genome-wide table)."))
    sec.append(("Sensitivity", sens_html))
    sec.append(("Limits and alternative explanations",
                "<ul><li><b>Relative states, not CN.</b> The r grid is prespecified, accessibility and dosage are confounded, and the baseline (GC-decile level of pyEpi cn_call==1 bins) is itself circular. A region the HMM calls r = 1 can sit at any copy number the pseudobulk baseline happens to share.</li>"
                "<li><b>Baseline variance is the main driver and is poorly understood.</b> A relative SD near 0.5 per bin after GC adjustment indicates mappability/accessibility structure that is not modelled. It was treated as independent per bin although it is regional, which makes posteriors overconfident where baseline structure is shared across neighbouring bins; "
                "a state that fits a regional baseline offset (for example a mappability-driven offset in 19q or 9p) is indistinguishable from copy-level change in this model.</li>"
                "<li><b>AT is not independent of depth.</b> It sums SNP-site depth, so allele evidence rises and falls with the same accessibility; it is not an independent confirmation, and an AT-poor region (9p21) gives no allele evidence at all.</li>"
                "<li><b>Allelic null.</b> exp 06 found about 31% of uncalled random windows exceed the balanced beta-binomial; here that appears as imbalance mass in some base-level regions. rho and the attenuation factor a come from the same data (two gain loci), a is applied to losses by extrapolation, and the copy-neutral-imbalance delta (0.25 nominal) is arbitrary.</li>"
                "<li><b>Not tuned or calibrated.</b> Hazard, pi, grid and the 0.8 flag threshold are fixed by hand. Posterior probabilities describe this model only. No simulation recovery study was run, so false-positive rates of the HMM are unknown (the cn_call==1 diagnostics above are not error rates because cn_call==1 is not truth).</li>"
                "<li>Pooled all-cell and cluster0 pseudobulks dilute subclonal events; cluster0 half-split dispersion was negative (floored), so its sampling-only variant is effectively Poisson.</li></ul>"))
    supports = ["Large chr7 windows have most posterior mass on loss-level states in both pseudobulks (all: " + ", ".join(f"{g('all', r, 'mass_depth_loss'):.2f}" for r in chr7) + ") with little imbalance mass, consistent with pyEpi and with the AD structure.",
                "MYC is called gain-level with imbalance by both layers (gain mass 1.00, imbalance mass 1.00), consistent with exp 06's strong allelic excess there.",
                f"AD breakpoints are mostly (about {ad2.value:.0%} within 2 bins) matched by an HMM depth boundary, far above a random-placement reference.",
                "A between-bin baseline-variance term is necessary: sampling-only NB fragments the genome into thousands of runs."]
    not_supports = ["That the relative depth levels (0.25, 0.5, ...) are copy numbers; r = 0.25 states in chr7 mean only 'lower than the baseline by a large factor'.",
                    "Any claim about 9p21, 19q13-qter or 5q13-qter at confident posterior: all are flagged insufficient evidence in the full model (5q has moderate gain mass, about 0.6, split between r = 1.5 and 2 and between balanced and imbalanced states).",
                    "That the allele layer adds independent confirmation: it is built from the same SNP-site depth and is silent where AT is low.",
                    "That HMM boundaries finer than the AD breakpoints are real; most HMM depth boundaries have no AD breakpoint within 5 bins and may follow baseline structure.",
                    "Calibrated probabilities or error rates of any kind."]
    decision = ("<p>This supports keeping the recommendation's position that depth segmentation should be paired with a <b>separately reported allelic evidence label and an abstention rule</b>, and adds a quantified caveat: "
                "a joint model's value here comes mostly from (a) treating missing alleles as no evidence, (b) abstaining when the posterior is split, and (c) an explicit between-bin baseline-variance term. "
                "It does not justify replacing the AD scan with an HMM as production behaviour: the HMM cannot be validated without a recovery simulation, and its calls depend on a baseline whose regional structure is unmodelled. "
                "Suggested wording for recommendation.md: 'a joint depth+allele HMM is a useful exploratory abstention layer; adopting it requires baseline-variance modelling and a simulation-based calibration first.'</p>")
    next_steps = ["Simulation recovery: plant relative depth steps and allele imbalance into the empirical baseline (resampled uncalled blocks) and measure state recovery, boundary error and false-state rate under this exact model.",
                  "Replace the independent per-bin baseline variance with an autocorrelated (regional) baseline or an external reference (matched normal / mappability track) and rerun; check whether the r = 1.25 'gain' excess over pyEpi neutral bins shrinks.",
                  "Replace the fixed r grid by a hierarchical level prior and compare the posterior at 5q and 19q.",
                  "Test sensitivity to the transition hazard (0.003, 0.03) and to pi, which were fixed by hand.",
                  "Use SNP-level data with phase blocks and unique-molecule counts, if obtainable, to replace the phase-agnostic mixture."]
    glossary = {
        "source": "Pseudobulk: all (3,460 cells) or cluster0 (WNN cluster 0, 1,286 cells).",
        "variant": "full; no_baseline_var (tau2 = 0); no_allele (allele emission off); baseline_all_bins; baseline_median_only (see Method).",
        "region / family": "Saved window (regions.tsv part == target) or literature locus.",
        "n_bins / bins_AT_ge1": "Retained bins in the region / bins with pooled AT of at least 1 (bins with AT = 0 contribute no allele evidence).",
        "mass_depth_loss / mass_depth_gain": "Mean over the region's bins of the posterior probability of states with r &lt; 1 / r &gt; 1 (relative depth, not CN).",
        "mass_imbalance": "Mean posterior on imbalanced allele states; blank/NaN when the allele layer is off.",
        "max_joint_state_prob / max_joint_state": "Largest probability and identity of a joint (r, allele) state in the region-mean posterior; below 0.8 flags insufficient evidence.",
        "max_depth_prob / map_depth_r": "Largest depth-level probability (allele summed out) and that level.",
        "viterbi_modal_r / viterbi_modal_r_frac": "Most frequent Viterbi relative depth in the region and its share of bins.",
        "p_r0.25 ... p_r2": "Posterior mass on each relative depth level, allele summed out.",
        "insufficient_evidence / flag_reasons": "True if max_joint_state_prob &lt; 0.8 or no allele observation; reasons listed (allele_layer_off is informational).",
        "pyepi_cn_call_modal": "Most common pyEpi cn_call in the region (relative level; not truth).",
        "state (segments.tsv)": "Joint state name, e.g. r0.5_bal = relative depth 0.5, balanced allele.",
        "rel_depth_r / allele_state / delta_if_imbalanced": "Relative depth level of the run, its allele state and the |delta| assumed for imbalanced states.",
        "mean_post_of_state": "Mean posterior of the run's own state over its bins.",
        "bin_start_in_chrom / bin_end_in_chrom": "Chromosome-local bin indices (same convention as the AD breakpoints' absolute_bin).",
        "table / key1 / key2 / n / value (comparison.tsv)": "Comparison table type, row and column keys, count and fraction.",
        "boundary_match_pm2bins / pm5bins": "Fraction of AD breakpoints with an HMM depth boundary within 2 / 5 bins (and vice versa); chance reference in note.",
        "full_* / alt_* (sensitivity.tsv)": "Statistic in the full model / in the variant named in variant.",
        "map_depth_changed / loss_gain_class_changed": "Modal depth level changed / loss or gain majority (mass &gt; 0.5) changed relative to the full model.",
        "bins_depth_class_changed_frac / bins_depth_r_changed_frac": "Genome-wide fraction of bins whose Viterbi class / level differs from the full model.",
        "phi_s / tau2 / phi": "NB dispersion terms: cell sampling / extra between-bin baseline / total (var = mu + phi mu^2).",
        "rho": "Beta-binomial overdispersion of allele counts (exp 06 estimator).",
        "a / delta": "Empirical attenuation factor of allelic imbalance relative to nominal copy-based expectation / |p - 0.5| assigned to an imbalanced state.",
        "c": "Mean-matching scalar applied to GC-decile medians.",
        "AD / cn_call": "Anderson-Darling breakpoints (raw_accepted_breakpoints.tsv) / pyEpi relative level per bin; neither is truth."}
    C.write_html_report(HERE / "report.html", "Joint depth + allele HMM",
                        "Experiment 13 · method record B:3 · exploratory relative-state model on the all and cluster0 raw pseudobulks", sec,
                        {"command": summary["command"], "seeds": summary["seeds"], "settings": summary["settings"], "environment": summary["environment"],
                         "files": ["posteriors_regions.tsv", "segments.tsv", "comparison.tsv", "sensitivity.tsv", "summary.json"]},
                        key_findings=key_findings, verdict=verdict, glossary=glossary, supports=supports, not_supports=not_supports,
                        decision_impact=decision, next_steps=next_steps)


if __name__ == "__main__":
    main()
