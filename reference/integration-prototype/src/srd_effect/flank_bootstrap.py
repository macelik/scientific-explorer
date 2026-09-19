"""bootstrap.py — paired within-cluster resampling and summaries (plan section 2).

For each replicate we draw n_g cell indices WITH replacement within each cluster
and use the SAME indices across every segment and both flanks, preserving
within-cell dependence (we resample whole cells, never bins or flanks
independently). The whole-dataset replicate is the cell-count-weighted mean of
the cluster mean-sum vectors (``core.weighted_all``), not a sum and not an
equal-weight average.

Per flank, SE is the sample SD of the replicate effects (ddof=1) — do NOT divide
by sqrt(n_replicates). ``Z_boot = observed / SE`` is a bootstrap-standardized
diagnostic, not a calibrated Z-statistic. CIs are 95% percentile intervals.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import flank_core as core

CI_LO, CI_HI = 2.5, 97.5


def bootstrap_indices(cluster_members, n_replicates, seed):
    """Per cluster, n_replicates arrays of its own member cell indices (with replacement).

    Deterministic via ``numpy.random.SeedSequence(seed)``: one child generator
    per (cluster, replicate) pair.
    """
    base = np.random.SeedSequence(seed)
    labels = sorted(cluster_members.keys())
    n_clusters = len(labels)
    children = base.spawn(n_clusters * n_replicates)   # spawn ONCE, not per pair
    out = {}
    for ci, lab in enumerate(labels):
        member = np.asarray(cluster_members[lab], dtype=int)
        n_g = len(member)
        if n_g == 0:
            raise ValueError(f"cluster {lab} has no cells")
        reps = []
        for r in range(n_replicates):
            rng = np.random.default_rng(children[ci * n_replicates + r])
            reps.append(member[rng.integers(0, n_g, size=n_g)])
        out[lab] = reps
    return out


def replicate_mean_sums(q, indices, cluster_sizes, n_replicates):
    """Mean-sum matrices per source: dict source -> (n_replicates, n_segments).

    ``cluster_sizes`` may be a dict keyed by cluster label or an array in the
    sorted-label order.
    """
    labels = sorted(indices.keys())
    if isinstance(cluster_sizes, dict):
        sizes_arr = np.asarray([int(cluster_sizes[lab]) for lab in labels], dtype=float)
    else:
        sizes_arr = np.asarray(cluster_sizes, dtype=float)
    out = {}
    for lab in labels:
        rows = np.empty((n_replicates, q.shape[1]), dtype=float)
        for r in range(n_replicates):
            rows[r] = core.mean_sums_from_indices(q, indices[lab][r])
        out[lab] = rows
    all_rows = np.empty((n_replicates, q.shape[1]), dtype=float)
    for r in range(n_replicates):
        all_rows[r] = core.weighted_all(
            np.stack([out[lab][r] for lab in labels]), sizes_arr)
    out["all"] = all_rows
    return out


def build_segment_neighbours(segments):
    """Left/right neighbour segment indices per segment (None at chromosome ends)."""
    n = len(segments)
    left, right = [None] * n, [None] * n
    by_chr = {}
    for i, seg in enumerate(segments):
        by_chr.setdefault(seg["chromosome"], []).append(i)
    for ch, idxs in by_chr.items():
        idxs = sorted(idxs, key=lambda i: (segments[i]["seg_bin_start"],
                                           segments[i]["seg_bin_end"]))
        for k, i in enumerate(idxs):
            if k > 0:
                left[i] = idxs[k - 1]
            if k < len(idxs) - 1:
                right[i] = idxs[k + 1]
    return left, right


def flank_delta_matrices(m, left, right):
    """(dL, dR), each (n_reps, n_segments), from a mean-level matrix ``m``.

    A log ratio is NaN when either mean is zero (unmeasured zero, not a real
    level); the ``zero_mean`` mask is returned too. Chromosome-end segments get
    NaN on the missing side.
    """
    n_reps, n_seg = m.shape
    dL = np.full((n_reps, n_seg), np.nan)
    dR = np.full((n_reps, n_seg), np.nan)
    zero_mask = np.zeros((n_reps, n_seg), dtype=bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        for j in range(n_seg):
            mj = m[:, j]
            z = mj <= 0
            if left[j] is not None:
                ml = m[:, left[j]]
                r = mj / ml
                bad = z | (ml <= 0)
                dL[:, j] = np.where(bad, np.nan, np.log2(r))
                zero_mask[:, j] = bad
            if right[j] is not None:
                mr = m[:, right[j]]
                r = mj / mr
                bad = z | (mr <= 0)
                dR[:, j] = np.where(bad, np.nan, np.log2(r))
                zero_mask[:, j] |= bad
    return dL, dR, zero_mask


def summarize_flank(delta_obs, delta_star):
    """Bootstrap summary for one flank (see module docstring for semantics).

    Returns a dict with se, z_boot, ci_lo, ci_hi, median, finite_frac,
    nonfinite_replicates, zero_bootstrap_variance, and a withheld reason when
    applicable.
    """
    star = np.asarray(delta_star, dtype=float)
    finite = np.isfinite(star)
    frac = float(finite.mean())
    if frac < 1.0:
        return {"se": np.nan, "z_boot": np.nan, "ci_lo": np.nan, "ci_hi": np.nan,
                "median": np.nan, "finite_frac": frac,
                "nonfinite_replicates": int((~finite).sum()),
                "zero_bootstrap_variance": False,
                "withheld": "nonfinite_replicates"}
    se = float(np.std(star, ddof=1))
    lo, hi = float(np.percentile(star, CI_LO)), float(np.percentile(star, CI_HI))
    med = float(np.median(star))
    if se == 0.0:
        return {"se": 0.0, "z_boot": np.nan, "ci_lo": lo, "ci_hi": hi, "median": med,
                "finite_frac": 1.0, "nonfinite_replicates": 0,
                "zero_bootstrap_variance": True, "withheld": "zero_bootstrap_variance"}
    return {"se": se, "z_boot": float(delta_obs / se), "ci_lo": lo, "ci_hi": hi,
            "median": med, "finite_frac": 1.0, "nonfinite_replicates": 0,
            "zero_bootstrap_variance": False, "withheld": ""}


def summarize_source(source, n_cells, Y_obs, Y_star, exposures, segments,
                     left, right, n_replicates):
    """Per-segment summary rows for one source.

    ``Y_obs``: observed mean-sum vector (n_segments). ``Y_star``: replicate
    mean-sums (n_replicates, n_segments). Returns a list of dicts, one per
    segment, with observed + bootstrap left/right deltas, SEs, Z_boot, CIs,
    covariance, SRDs, and status flags.
    """
    E = exposures
    m_obs = Y_obs / E
    m_star = Y_star / E
    dL_obs, dR_obs, zm_obs = flank_delta_matrices(m_obs[None, :], left, right)
    dL_star, dR_star, zm_star = flank_delta_matrices(m_star, left, right)

    too_few = n_cells < 2

    def withheld():
        return {"se": np.nan, "z_boot": np.nan, "ci_lo": np.nan, "ci_hi": np.nan,
                "median": np.nan, "finite_frac": np.nan, "nonfinite_replicates": 0,
                "zero_bootstrap_variance": False, "withheld": "too_few_cells"}

    def _missing_flank():
        return {"se": np.nan, "z_boot": np.nan, "ci_lo": np.nan, "ci_hi": np.nan,
                "median": np.nan, "finite_frac": np.nan, "nonfinite_replicates": 0,
                "zero_bootstrap_variance": False, "withheld": "terminal_flank_missing"}

    rows = []
    for j, seg in enumerate(segments):
        has_l = left[j] is not None
        has_r = right[j] is not None
        if not has_l and not has_r:
            continue  # single-segment chromosome: no comparison possible
        # SRD uses the exposure-aware formula on Y sums (not means)
        if has_l:
            srd_l, st_l = core.srd_target_with_status(Y_obs[j], E[j],
                                                      Y_obs[left[j]], E[left[j]])
        else:
            srd_l, st_l = np.nan, {"uninformative_zero_pair": False}
        if has_r:
            srd_r, st_r = core.srd_target_with_status(Y_obs[j], E[j],
                                                      Y_obs[right[j]], E[right[j]])
        else:
            srd_r, st_r = np.nan, {"uninformative_zero_pair": False}
        obs_l = float(dL_obs[0, j]) if (has_l and np.isfinite(dL_obs[0, j])) else np.nan
        obs_r = float(dR_obs[0, j]) if (has_r and np.isfinite(dR_obs[0, j])) else np.nan
        if not has_l:
            sum_l = _missing_flank()
        elif too_few:
            sum_l = withheld()
        else:
            sum_l = summarize_flank(obs_l, dL_star[:, j])
        if not has_r:
            sum_r = _missing_flank()
        elif too_few:
            sum_r = withheld()
        else:
            sum_r = summarize_flank(obs_r, dR_star[:, j])
        cov = np.nan
        both = np.isfinite(dL_star[:, j]) & np.isfinite(dR_star[:, j])
        if both.sum() >= 2:
            cov = float(np.cov(dL_star[both, j], dR_star[both, j], ddof=1)[0, 1])
        rows.append({
            "source": source,
            "n_cells": n_cells,
            "chromosome": seg["chromosome"],
            "segment_id": seg["segment_id"],
            "seg_bin_start": int(seg["seg_bin_start"]),
            "seg_bin_end": int(seg["seg_bin_end"]),
            "seg_length": int(seg["seg_length"]),
            "is_terminal": bool(seg.get("is_terminal", not (has_l and has_r))),
            "exposure": int(E[j]),
            "mean_level": float(m_obs[j]),
            "mean_left": float(m_obs[left[j]]) if has_l else np.nan,
            "mean_right": float(m_obs[right[j]]) if has_r else np.nan,
            "srd_left": srd_l,
            "srd_left_zero_pair": st_l["uninformative_zero_pair"],
            "srd_right": srd_r,
            "srd_right_zero_pair": st_r["uninformative_zero_pair"],
            "zero_mean": bool(zm_obs[0, j]),
            "delta_left_obs": obs_l,
            "delta_right_obs": obs_r,
            "delta_left_median": sum_l["median"],
            "delta_right_median": sum_r["median"],
            "delta_left_ci_lo": sum_l["ci_lo"],
            "delta_left_ci_hi": sum_l["ci_hi"],
            "delta_right_ci_lo": sum_r["ci_lo"],
            "delta_right_ci_hi": sum_r["ci_hi"],
            "se_left": sum_l["se"],
            "se_right": sum_r["se"],
            "z_boot_left": sum_l["z_boot"],
            "z_boot_right": sum_r["z_boot"],
            "cov_left_right": cov,
            "finite_frac_left": sum_l["finite_frac"],
            "finite_frac_right": sum_r["finite_frac"],
            "nonfinite_replicates_left": sum_l["nonfinite_replicates"],
            "nonfinite_replicates_right": sum_r["nonfinite_replicates"],
            "withheld_left": sum_l["withheld"],
            "withheld_right": sum_r["withheld"],
            "calibration_status": "not_validated",
        })
    return rows
