"""merge_core.py - pure per-cluster depth-segment merge prototype.

Implements the 2026-09-07 depth-segment prototype rule WITHOUT any IO:

  * segment level m = per-bin MEDIAN of a mean-over-cells depth pseudobulk;
  * a boundary between adjacent current segments A,B is NOMINATED only when it
    does not touch a chromosome-end (terminal) segment and at least one of the
    two flanking segments is an internal segment <= ``max_small`` bins
    ("size decides where we look, never whether we remove");
  * observed effect delta = log2(m_A / m_B);
  * a nominated boundary is ELIGIBLE iff observed |delta| < tau AND at least
    ``stab``-fraction of whole-cell bootstrap replicates satisfy
    |delta_rep| < tau, AND every replicate effect is valid (zero/undefined
    replicate medians make the boundary unresolved, i.e. retained);
  * a merge additionally requires the CHAIN GUARD: inside the proposed merged
    interval the range of the original constituent segments' log2 median
    levels must stay < tau;
  * among eligible + guard-passing candidates we merge the one with the
    smallest observed |delta| (exact ties -> lower genomic coordinate),
    recompute from the SAME saved replicate per-bin pseudobulks, and repeat
    until none remains eligible.

Terminal segments and long-long boundaries are never merged.  Nothing here
prunes a segment because it is small, and SRD/phi never decides a merge.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np


@dataclass
class Seg:
    """One current segment on one chromosome.

    start/end are half-open chromosome-local bin coordinates.  orig is the
    sorted tuple of ORIGINAL (pre-merge) flat-partition segment ids contained.
    """
    start: int
    end: int
    orig: tuple = field(default_factory=tuple)


def flat_partition(boundary_bins, n_bins):
    """Half-open flat segments from accepted boundary bin coords (end-exclusive)."""
    cuts = [0] + sorted(int(b) for b in boundary_bins) + [int(n_bins)]
    segs = []
    for i in range(len(cuts) - 1):
        s, e = cuts[i], cuts[i + 1]
        if e > s:
            segs.append(Seg(s, e, (len(segs),)))
    return segs


def segment_medians(pb, segs):
    """Per-segment median of the (observed) per-bin pseudobulk."""
    pb = np.asarray(pb, dtype=float)
    return np.array([float(np.median(pb[s.start:s.end])) for s in segs])


def segment_replicate_medians(pb_boot, segs):
    """(n_boot, n_seg) medians of each replicate pseudobulk per segment."""
    pb_boot = np.asarray(pb_boot, dtype=float)
    n_boot = pb_boot.shape[0]
    out = np.empty((n_boot, len(segs)), dtype=float)
    for j, s in enumerate(segs):
        out[:, j] = np.median(pb_boot[:, s.start:s.end], axis=1)
    return out


def _log2ratio(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log2(a / b)


def original_log2_medians(pb_obs, orig_segs):
    """log2 of the observed median of each ORIGINAL segment (NaN if median<=0)."""
    med = segment_medians(pb_obs, orig_segs)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(med > 0, np.log2(med), np.nan)


def constituent_log2_range(orig_ids, orig_log2med):
    """Range of original-constituent log2 medians; NaN if undefined."""
    if len(orig_ids) == 0:
        return np.nan
    vals = np.asarray(orig_log2med, dtype=float)[list(orig_ids)]
    if np.isnan(vals).any():
        return np.nan
    return float(np.nanmax(vals) - np.nanmin(vals))


def merged_seg(left: Seg, right: Seg) -> Seg:
    return Seg(left.start, right.end, left.orig + right.orig)


def boundary_is_terminal_touching(i, n_seg):
    """Boundary i separates current segments i, i+1.  True if either touches a
    chromosome-end (terminal) segment (current index 0 or n_seg-1)."""
    return i <= 0 or i + 1 >= n_seg - 1


def boundary_candidate(i, segs, max_small):
    """Nominated boundary: flanks both non-terminal AND at least one flank is an
    internal segment <= max_small bins."""
    if boundary_is_terminal_touching(i, len(segs)):
        return False
    len_l = segs[i].end - segs[i].start
    len_r = segs[i + 1].end - segs[i + 1].start
    return min(len_l, len_r) <= max_small


def evaluate_boundary(mobs_i, mobs_j, mboot_i, mboot_j, tau, stab, n_boot,
                      need_all_valid=True):
    """Eligibility components for the boundary between segments i and j.

    Returns dict with delta_obs, frac_under_tau, n_valid, eligible, reason.
    reason is None when eligible (subject to chain guard elsewhere).
    """
    d_obs = np.nan
    frac_under = np.nan
    n_valid = n_boot
    eligible = False
    reason = None
    if mobs_i > 0 and mobs_j > 0:
        d_obs = float(_log2ratio(mobs_i, mobs_j))
        if mboot_i is not None and mboot_j is not None:
            d_rep = _log2ratio(mboot_i, mboot_j)
            ok = np.isfinite(d_rep)
            n_valid = int(ok.sum())
            if need_all_valid and n_valid < n_boot:
                reason = "invalid_replicates"
            else:
                under = np.abs(d_rep) < tau
                frac_under = float(under[ok].mean()) if ok.any() else np.nan
                if abs(d_obs) < tau and (frac_under if np.isfinite(frac_under) else 0) >= stab:
                    eligible = True
                elif abs(d_obs) >= tau:
                    reason = "effect_ge_tau"
                else:
                    reason = "stability_below_cutoff"
        else:
            # no bootstrap available: unresolved unless observed delta < tau
            # with no stability evidence would be overreach -> retain
            reason = "no_bootstrap_available"
    else:
        reason = "observed_median_zero"
    return {
        "delta_obs": d_obs,
        "frac_under_tau": frac_under,
        "n_valid": int(n_valid),
        "eligible": bool(eligible),
        "reason": reason,
    }


def initial_boundary_scores(pb_obs, pb_boot, boundary_bins, n_bins, tau, stab,
                            max_small=100, need_all_valid=True):
    """Snapshot scores for every ORIGINAL boundary before any merge.

    Returns (orig_segs, rows); row j corresponds to the boundary between
    original segments j and j+1.
    """
    orig_segs = flat_partition(boundary_bins, n_bins)
    n_seg = len(orig_segs)
    mobs = segment_medians(pb_obs, orig_segs)
    mboot = None
    if pb_boot is not None:
        mboot = segment_replicate_medians(pb_boot, orig_segs)
    olog = original_log2_medians(pb_obs, orig_segs)
    n_boot = pb_boot.shape[0] if pb_boot is not None else 0
    rows = []
    for j in range(n_seg - 1):
        a, b = orig_segs[j], orig_segs[j + 1]
        len_l = a.end - a.start
        len_r = b.end - b.start
        cand = boundary_candidate(j, orig_segs, max_small)
        ev = evaluate_boundary(mobs[j], mobs[j + 1],
                               mboot[:, j] if mboot is not None else None,
                               mboot[:, j + 1] if mboot is not None else None,
                               tau, stab, n_boot, need_all_valid=need_all_valid)
        if not cand:
            ev["eligible"] = False
            if boundary_is_terminal_touching(j, n_seg):
                ev["reason"] = "terminal_unchanged"
            elif len_l > max_small and len_r > max_small:
                ev["reason"] = "long_long_not_candidate"
            else:
                ev["reason"] = "not_candidate"
        gr = constituent_log2_range(a.orig + b.orig, olog)
        rows.append({
            "orig_seg_left": int(a.orig[0]),
            "orig_seg_right": int(b.orig[0]),
            "boundary_bin": int(b.start),
            "len_left": int(len_l),
            "len_right": int(len_r),
            "median_left": float(mobs[j]) if np.isfinite(mobs[j]) else np.nan,
            "median_right": float(mobs[j + 1]) if np.isfinite(mobs[j + 1]) else np.nan,
            "delta_obs": ev["delta_obs"],
            "frac_under_tau": ev["frac_under_tau"],
            "n_valid": ev["n_valid"],
            "eligible": ev["eligible"],
            "candidate": bool(cand),
            "guard_range": gr,
            "reason": ev["reason"] or ("eligible" if ev["eligible"] else "unresolved"),
        })
    return orig_segs, rows


def has_nominated_boundary(orig_segs, max_small):
    n_seg = len(orig_segs)
    return any(boundary_candidate(i, orig_segs, max_small) for i in range(n_seg - 1))


def run_merges(pb_obs, pb_boot, boundary_bins, n_bins, tau, stab,
               max_small=100, need_all_valid=True):
    """Deterministic sequential merge over one chromosome.

    Returns dict with orig_segs, seg_medians_obs, orig_log2_medians,
    final_segs, merges, decided (every original boundary -> final reason),
    final_boundary_order, n_boot.
    """
    n_boot = pb_boot.shape[0] if pb_boot is not None else 0
    orig_segs = flat_partition(boundary_bins, n_bins)
    olog = original_log2_medians(pb_obs, orig_segs)
    mobs_orig = segment_medians(pb_obs, orig_segs)

    segs = list(orig_segs)
    merges = []
    outcome = [None] * max(0, len(orig_segs) - 1)
    last_reason = {}  # orig boundary id -> last reason it failed as a candidate

    def boundary_orig_id(segs_, i):
        """Original boundary id removed when merging current segments i and i+1."""
        return int(segs_[i].orig[-1])

    iteration = 0
    while True:
        n_seg = len(segs)
        if n_seg < 3:
            break
        mobs = segment_medians(pb_obs, segs)
        mboot = None
        if pb_boot is not None:
            mboot = segment_replicate_medians(pb_boot, segs)
        scored = []
        for i in range(n_seg - 1):
            oid = boundary_orig_id(segs, i)
            if not boundary_candidate(i, segs, max_small):
                continue
            ev = evaluate_boundary(
                mobs[i], mobs[i + 1],
                mboot[:, i] if mboot is not None else None,
                mboot[:, i + 1] if mboot is not None else None,
                tau, stab, n_boot, need_all_valid=need_all_valid)
            merged = merged_seg(segs[i], segs[i + 1])
            gr = constituent_log2_range(merged.orig, olog)
            guard_ok = bool(np.isfinite(gr) and gr < tau)
            eligible = ev["eligible"] and guard_ok
            reason = ev["reason"]
            if ev["eligible"] and not guard_ok:
                reason = "chain_guard"
            if reason is not None and outcome[oid] is None:
                last_reason[oid] = reason  # last time it was NOT mergeable
            scored.append({
                "oid": oid,
                "i": i,
                "delta_obs": ev["delta_obs"],
                "frac_under_tau": ev["frac_under_tau"],
                "n_valid": ev["n_valid"],
                "eligible": bool(eligible),
                "guard_range": gr,
                "guard_ok": bool(guard_ok),
                "reason": reason,
                "abs_delta": abs(ev["delta_obs"]) if np.isfinite(ev["delta_obs"]) else np.inf,
                "boundary_bin": int(segs[i + 1].start),
            })
        eligible_candidates = [s for s in scored if s["eligible"] and s["guard_ok"]]
        if not eligible_candidates:
            break
        # smallest observed |delta| first; "exact ties" -> lower genomic
        # coordinate.  abs-delta is quantized to 10 decimals so float noise on
        # genuinely equal effects cannot change the deterministic order.
        eligible_candidates.sort(key=lambda s: (round(s["abs_delta"], 10),
                                                s["boundary_bin"]))
        pick = eligible_candidates[0]
        i = pick["i"]
        iteration += 1
        oid = boundary_orig_id(segs, i)
        outcome[oid] = {"removed": True, "iteration": iteration}
        merges.append({
            "iteration": iteration,
            "orig_boundary_id": oid,
            "boundary_bin": int(segs[i + 1].start),
            "delta_obs": pick["delta_obs"],
            "frac_under_tau": pick["frac_under_tau"],
            "n_valid": pick["n_valid"],
            "guard_range": pick["guard_range"],
            "left_orig": list(segs[i].orig),
            "right_orig": list(segs[i + 1].orig),
        })
        segs = segs[:i] + [merged_seg(segs[i], segs[i + 1])] + segs[i + 2:]

    # decided maps EVERY original boundary id -> final reason (removed/kept)
    decided = {}
    for oid in range(len(orig_segs) - 1):
        if outcome[oid] is not None and outcome[oid]["removed"]:
            decided[oid] = "removed_iter%02d" % outcome[oid]["iteration"]
    n_seg = len(segs)
    for i in range(n_seg - 1):
        oid = boundary_orig_id(segs, i)
        if oid in decided:
            continue
        if boundary_is_terminal_touching(i, n_seg):
            decided[oid] = "terminal_unchanged"
        elif oid in last_reason:
            # last time this boundary was a merge candidate it did not qualify
            decided[oid] = last_reason[oid]
        elif min(segs[i].end - segs[i].start, segs[i + 1].end - segs[i + 1].start) > max_small:
            decided[oid] = "long_long_not_candidate"
        else:
            decided[oid] = "never_candidate"

    final_oid_order = []
    for i in range(n_seg - 1):
        final_oid_order.append(boundary_orig_id(segs, i))

    return {
        "orig_segs": orig_segs,
        "seg_medians_obs": mobs_orig,
        "orig_log2_medians": olog,
        "final_segs": segs,
        "merges": merges,
        "decided": decided,
        "final_boundary_order": final_oid_order,
        "n_boot": n_boot,
    }


def exposure_aware_srd(y_target, exposure_target, y_reference, exposure_reference):
    """Signed root deviance between two exposures (srd_pruning formula).

    Descriptive only; never used for merge decisions in this prototype.
    """
    if exposure_target <= 0 or exposure_reference <= 0:
        return np.nan
    if y_target < 0 or y_reference < 0:
        return np.nan
    y_sum = y_target + y_reference
    if y_sum == 0:
        return 0.0
    rate_t = y_target / exposure_target
    rate_r = y_reference / exposure_reference
    lambda_hat = y_sum / (exposure_target + exposure_reference)
    mu_t = exposure_target * lambda_hat
    mu_r = exposure_reference * lambda_hat
    from scipy.special import xlogy
    D = 2.0 * (xlogy(y_target, y_target / mu_t) + xlogy(y_reference, y_reference / mu_r))
    return float(np.sign(rate_t - rate_r) * math.sqrt(max(float(D), 0.0)))
