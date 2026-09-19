"""Scientific computations for the explorer. Every statistic here is a direct
call into the supplied pyEpiAneufinder functions (via pyepi_adapter); nothing
is re-implemented. Functions that run in worker processes take plain arrays.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

from . import pyepi_adapter as sci


def _ensure(pkg_dir: Optional[str]) -> None:
    if not sci._loaded and pkg_dir:
        sci.load(pkg_dir)


# ------------------------------------------------------------------ nodes
def compute_node(x_node: np.ndarray) -> dict:
    """AD curve over every split of a node plus its two baselines.

    Returns d (length n-1; d[b-1] is the AD statistic of boundary b, i.e.
    left = x[:b], right = x[b:]), the argmax boundary, and baselines.
    """
    n = int(len(x_node))
    if n >= 2:
        d = np.asarray(sci.seq_dist_ad(x_node), dtype=np.float64)
        b = int(np.argmax(d)) + 1
        ad = float(d[b - 1])
    else:
        d = np.zeros(0)
        b = None
        ad = None
    nonfinite = int(np.sum(~np.isfinite(x_node)))
    return {
        "n": n,
        "ad": d.tolist(),
        "argmax_b": b,
        "argmax_ad": ad,
        "baseline_trimmed": float(sci.trimmed_mean_iqr(x_node, lb=False)) if n else None,
        "baseline_median": float(np.median(x_node)) if n else None,
        "x_nonfinite": nonfinite,
    }


def genome_baselines(xfull: np.ndarray) -> dict:
    return {
        "baseline_trimmed": float(sci.trimmed_mean_iqr(xfull, lb=False)),
        "baseline_median": float(np.median(xfull)),
        "zero_fraction": float(np.mean(xfull == 0)),
        "x_nonfinite": int(np.sum(~np.isfinite(xfull))),
    }


# ------------------------------------------------------------------ tests
def run_split_tests(x_node: np.ndarray, b: int, xfull: np.ndarray, n_permutations: int, seed: int,
                    pkg_dir: Optional[str] = None) -> dict:
    """Exact pyEpi local + global permutation tests for boundary b of a node.
    Runs in a worker process."""
    _ensure(pkg_dir)
    left, right = x_node[:b], x_node[b:]
    if len(left) == 0 or len(right) == 0:
        raise ValueError("boundary must leave both children non-empty")
    obs = float(sci.dist_ad(left, right))
    p_local = float(sci.permutation_test_ad(left, right, obs, n_permutations, random_state=seed))
    p_global = float(sci.global_permutation_test_ad(xfull, len(left), len(right), obs, n_permutations, random_state=seed))
    return {
        "observed_ad": obs,
        "p_local": p_local,
        "p_global": p_global,
        "p_local_exact": f"{int(round(p_local * (n_permutations + 1)))}/{n_permutations + 1}",
        "p_global_exact": f"{int(round(p_global * (n_permutations + 1)))}/{n_permutations + 1}",
        "n_permutations": n_permutations,
        "seed": seed,
        "n_left": int(len(left)),
        "n_right": int(len(right)),
    }


# ------------------------------------------------------------- CN calling
def assign_cn(xfull: np.ndarray, clusters: List[int]) -> dict:
    cn_int, cn_cont, holmes, watson, cn5, best_s = sci.assign_gainloss_new(xfull, list(clusters))
    return {
        "cn_int": np.asarray(cn_int),
        "cn_cont": np.asarray(cn_cont, dtype=float),
        "cn5": np.asarray(cn5, dtype=float),
        "best_s": float(best_s),
    }


def segment_table(xfull: np.ndarray, clusters: np.ndarray, res: dict, offset: int, n: int, starts: np.ndarray,
                  ends: np.ndarray) -> List[dict]:
    """Per-segment table for one chromosome (chromosome-local bin bounds)."""
    sl = slice(offset, offset + n)
    cl = clusters[sl]
    x = xfull[sl]
    base = float(sci.trimmed_mean_iqr(xfull, lb=False))
    rows = []
    for seg_id in np.unique(cl):
        idx = np.where(cl == seg_id)[0]
        s, e = int(idx[0]), int(idx[-1]) + 1
        seg = x[s:e]
        info = res.get('segment_info', {}).get(int(seg_id))
        rows.append({
            "segment_id": int(seg_id),
            "start_bin": s,
            "end_bin": e,
            "n_bins": int(e - s),
            "start_bp": int(starts[offset + s]),
            "end_bp": int(ends[offset + e - 1]),
            "mean_x": float(np.mean(seg)),
            "median_x": float(np.median(seg)),
            "mean_norm": float(np.mean(seg) / base) if base > 0 else None,
            "summary_x": float(info['value'] * base) if info else float(np.mean(seg)),
            "summary_norm": float(info['value']) if info else (float(np.mean(seg)/base) if base > 0 else None),
            "retained_bins": int(info['retained']) if info else len(seg),
            "cn5": float(res["cn5"][offset + s]),
            "cn_int": int(res["cn_int"][offset + s]),
            "cn_cont": float(res["cn_cont"][offset + s]),
        })
    return rows


# ---------------------------------------------------------- event indexing
def index_cell(cell: str, xfull: np.ndarray, chroms: List[dict], prod_bps: Dict[str, List[int]],
               pkg_dir: Optional[str] = None) -> List[dict]:
    """Depth-0 candidate and both depth-1 candidates for every chromosome of a
    cell, computed from actual node AD curves. `accepted` means the candidate
    is a production breakpoint AND its parent was accepted (production only
    recurses into accepted splits), so it respects production's gating."""
    _ensure(pkg_dir)
    rows = []
    for c in chroms:
        x = xfull[c["offset"]:c["offset"] + c["n"]]
        prod = set(prod_bps.get(c["name"], []))
        n = len(x)
        if n < 2:
            continue
        d0 = np.asarray(sci.seq_dist_ad(x))
        b0 = int(np.argmax(d0)) + 1
        acc0 = b0 in prod
        rows.append({"cell": cell, "chrom": c["name"], "node": "root", "bp": b0, "ad": float(d0[b0 - 1]),
                     "accepted": bool(acc0), "parent_accepted": True, "n_left": b0, "n_right": n - b0})
        left, right = x[:b0], x[b0:]
        if len(left) >= 2:
            dL = np.asarray(sci.seq_dist_ad(left))
            bL = int(np.argmax(dL)) + 1
            rows.append({"cell": cell, "chrom": c["name"], "node": "left", "bp": bL, "ad": float(dL[bL - 1]),
                         "accepted": bool(acc0 and bL in prod), "parent_accepted": bool(acc0),
                         "n_left": bL, "n_right": len(left) - bL})
        if len(right) >= 2:
            dR = np.asarray(sci.seq_dist_ad(right))
            bR = int(np.argmax(dR)) + 1
            rows.append({"cell": cell, "chrom": c["name"], "node": "right", "bp": b0 + bR, "ad": float(dR[bR - 1]),
                         "accepted": bool(acc0 and (b0 + bR) in prod), "parent_accepted": bool(acc0),
                         "n_left": bR, "n_right": len(right) - bR})
    return rows
