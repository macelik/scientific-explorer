"""
srd_pruning.py — Exposure-aware signed root deviance (SRD) scoring

Pure functions for:
  - the exposure-aware SRD statistic (SciPy xlogy Poisson deviance)
  - flat final-segment reconstruction from accepted breakpoints
  - candidate scoring and non-destructive merge-proposal assignment

Design follows `2026-08-24-exposure-aware-srd-pruning-spike.md`. The depth
vector is the NONNEGATIVE, NON-LOG-TRANSFORMED pseudobulk depth used for
segmentation. log2 vs chromosome-median values are for plotting only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import xlogy

MAX_SMALL_SEGMENT_BINS = 100
DEFAULT_MIN_ABS_SRD = 3.3
DEFAULT_MIN_ABS_LOG2FC = 0.20


def exposure_aware_srd(y_target, exposure_target, y_reference, exposure_reference) -> float:
    """Signed root deviance between two Poisson-like exposures.

    Returns sign(rate_t - rate_r) * sqrt(max(D, 0)) with
        D = 2 * [xlogy(y_t, y_t/mu_t) + xlogy(y_r, y_r/mu_r)]
        lambda_hat = (y_t + y_r) / (e_t + e_r);  mu_* = e_* * lambda_hat
    Returns 0.0 when y_t + y_r == 0.
    """
    if exposure_target <= 0 or exposure_reference <= 0:
        raise ValueError("exposure must be strictly positive")
    if y_target < 0 or y_reference < 0:
        raise ValueError("depth values must be nonnegative")

    y_sum = y_target + y_reference
    if y_sum == 0:
        return 0.0

    e_t, e_r = float(exposure_target), float(exposure_reference)
    rate_t = y_target / e_t
    rate_r = y_reference / e_r
    lambda_hat = y_sum / (e_t + e_r)
    mu_t = e_t * lambda_hat
    mu_r = e_r * lambda_hat
    D = 2.0 * (xlogy(y_target, y_target / mu_t) + xlogy(y_reference, y_reference / mu_r))
    return float(np.sign(rate_t - rate_r) * np.sqrt(max(float(D), 0.0)))


def build_final_segments(accepted_bins, n_chrom_bins):
    """Flat ordered final segments from accepted breakpoint bin coordinates.

    Returns a list of (start_bin, end_bin) with end exclusive.
    Chromosome start/end are implicit (0, n_chrom_bins).
    """
    cuts = [0] + sorted(int(b) for b in accepted_bins) + [n_chrom_bins]
    return [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1) if cuts[i + 1] > cuts[i]]


def build_final_segments_with_depths(accepted, n_chrom_bins):
    """Like build_final_segments but retains the boundary recursion depths.

    accepted: iterable of (breakpoint_bin, depth). Chromosome edges get depth 0.
    Returns a list of (start, end, left_boundary_depth, right_boundary_depth).
    """
    ordered = sorted(accepted, key=lambda x: int(x[0]))
    segs = []
    prev_end = 0
    prev_depth = 0
    for b, d in ordered:
        b = int(b)
        if b > prev_end:
            segs.append((prev_end, b, prev_depth, int(d)))
        prev_end = b
        prev_depth = int(d)
    if prev_end < n_chrom_bins:
        segs.append((prev_end, n_chrom_bins, prev_depth, 0))
    return segs


def _segment_stats(depth, start, end):
    """Return (sum_of_depth, n_valid_bins) over [start, end)."""
    seg = depth[start:end]
    valid = np.isfinite(seg) & (seg >= 0)
    return float(seg[valid].sum()), int(valid.sum())


def assign_proposal(srd_left, srd_right, log2fc_left, log2fc_right,
                    min_abs_srd, min_abs_log2fc, seg_start, seg_end):
    """Dry-run merge proposal for one small segment.

    Returns (proposal, proposed_boundary_bin). proposed_boundary_bin is the bin
    coordinate of the boundary to remove for merge_* proposals, else None.
    """
    strong_l = abs(srd_left) >= min_abs_srd and abs(log2fc_left) >= min_abs_log2fc
    strong_r = abs(srd_right) >= min_abs_srd and abs(log2fc_right) >= min_abs_log2fc

    if strong_l and strong_r:
        if np.sign(srd_left) == np.sign(srd_right) and np.sign(srd_left) != 0:
            return "keep_focal", None
        return "ambiguous_keep", None
    if strong_r:            # weak left, strong right -> merge into left
        return "merge_left", int(seg_start)
    if strong_l:            # strong left, weak right -> merge into right
        return "merge_right", int(seg_end)
    # both weak -> merge toward the neighbour with smaller |log2FC|
    if abs(log2fc_left) <= abs(log2fc_right):
        return "merge_left", int(seg_start)
    return "merge_right", int(seg_end)


def score_candidates(depth, segments, group, chromosome,
                     max_small_bins=MAX_SMALL_SEGMENT_BINS):
    """Compute per-candidate stats (SRDs, log2 rate differences) — no proposal.

    segments: list of (start, end) final segments (end exclusive). depth: the
    nonnegative pseudobulk depth vector for the chromosome.
    Terminal segments (>1 neighbours missing) and segments longer than
    max_small_bins are skipped. No rate_* columns: the two log2FC columns ARE
    log2(rate_S/rate_L) and log2(rate_S/rate_R).
    """
    stats = [_segment_stats(depth, s, e) for s, e in segments]

    rows = []
    for idx in range(1, len(segments) - 1):
        s_start, s_end = segments[idx]
        if s_end - s_start > max_small_bins:
            continue
        l_start, l_end = segments[idx - 1]
        r_start, r_end = segments[idx + 1]
        y_s, e_s = stats[idx]
        y_l, e_l = stats[idx - 1]
        y_r, e_r = stats[idx + 1]

        srd_left = exposure_aware_srd(y_s, e_s, y_l, e_l)
        srd_right = exposure_aware_srd(y_s, e_s, y_r, e_r)
        rate_s = y_s / e_s if e_s else 0.0
        rate_l = y_l / e_l if e_l else 0.0
        rate_r = y_r / e_r if e_r else 0.0
        log2fc_left = float(np.log2((rate_s + 1e-12) / (rate_l + 1e-12)))
        log2fc_right = float(np.log2((rate_s + 1e-12) / (rate_r + 1e-12)))

        rows.append({
            "group": group,
            "chromosome": chromosome,
            "segment_id": idx,
            "seg_bin_start": s_start,
            "seg_bin_end": s_end,
            "seg_length": s_end - s_start,
            "left_bin_start": l_start,
            "left_bin_end": l_end,
            "right_bin_start": r_start,
            "right_bin_end": r_end,
            "srd_left": round(srd_left, 4),
            "srd_right": round(srd_right, 4),
            "log2fc_left": round(log2fc_left, 4),
            "log2fc_right": round(log2fc_right, 4),
        })
    return rows


def apply_thresholds(rows, min_abs_srd, min_abs_log2fc):
    """Add strong flags, proposal, and proposed boundary for a threshold pair.

    `rows` are the stats from score_candidates; the list is not mutated.
    """
    out = []
    for r in rows:
        srd_l = r["srd_left"]
        srd_r = r["srd_right"]
        fc_l = r["log2fc_left"]
        fc_r = r["log2fc_right"]
        strong_l = bool(abs(srd_l) >= min_abs_srd and abs(fc_l) >= min_abs_log2fc)
        strong_r = bool(abs(srd_r) >= min_abs_srd and abs(fc_r) >= min_abs_log2fc)
        proposal, boundary = assign_proposal(srd_l, srd_r, fc_l, fc_r,
                                             min_abs_srd, min_abs_log2fc,
                                             r["seg_bin_start"], r["seg_bin_end"])
        row = dict(r)
        row["strong_left"] = strong_l
        row["strong_right"] = strong_r
        row["proposal"] = proposal
        row["proposed_boundary_bin"] = boundary
        out.append(row)
    return out


def score_small_segments(depth, segments, group, chromosome,
                         max_small_bins=MAX_SMALL_SEGMENT_BINS,
                         min_abs_srd=DEFAULT_MIN_ABS_SRD,
                         min_abs_log2fc=DEFAULT_MIN_ABS_LOG2FC):
    """score_candidates + apply_thresholds (primary parameters)."""
    return apply_thresholds(
        score_candidates(depth, segments, group, chromosome, max_small_bins=max_small_bins),
        min_abs_srd, min_abs_log2fc)


def mark_conflicts(rows):
    """Flag candidate segments that nominate the same boundary for removal.

    A conflict is when two or more candidates within the SAME group and
    chromosome propose removing the same boundary bin. Accepts a list of dicts
    or a pandas DataFrame and returns the same type.
    """
    if isinstance(rows, pd.DataFrame):
        was_df = True
        rows_list = rows.to_dict("records")
    else:
        was_df = False
        rows_list = rows
    by_boundary = {}
    for i, r in enumerate(rows_list):
        b = r.get("proposed_boundary_bin")
        if b is None or pd.isna(b):
            continue
        key = (r.get("group"), r.get("chromosome"), b)
        by_boundary.setdefault(key, []).append(i)
    conflict_set = set()
    for key, idxs in by_boundary.items():
        if len(idxs) > 1:
            conflict_set.update(idxs)
    for i, r in enumerate(rows_list):
        r["conflict"] = i in conflict_set
    return pd.DataFrame(rows_list) if was_df else rows_list


def merge_until_converged(depth, segments,
                          max_small_bins=MAX_SMALL_SEGMENT_BINS,
                          min_abs_srd=DEFAULT_MIN_ABS_SRD,
                          min_abs_log2fc=DEFAULT_MIN_ABS_LOG2FC,
                          max_iter=100):
    """Iteratively merge non-focal small segments with rescoring.

    Rules:
      - merge_left / merge_right: remove the proposed boundary (S merges into that neighbour).
      - both-weak: merge toward the smaller |log2FC| side (existing assign_proposal rule).
      - ambiguous_keep (transition): merge toward the smaller-|SRD| side.
      - conflicts (two adjacent segments both nominating the same boundary):
        the boundary appears once in the removal set -> both merge into one segment.
      - keep_focal: kept (boundary retained).

    All proposed boundaries in one pass are removed simultaneously, then scores
    are recomputed (segments only grow, so the loop terminates; max_iter caps it).

    Returns (final_segments, removed_boundaries, history).
      final_segments: list of (start, end)
      removed_boundaries: set of bin positions removed
      history: list of (iteration, sorted removed-bins in that pass)
    """
    segs = [tuple(s) for s in segments]
    removed = set()
    history = []
    for it in range(max_iter):
        stats = score_candidates(depth, segs, "g", "chr", max_small_bins=max_small_bins)
        proposals = apply_thresholds(stats, min_abs_srd, min_abs_log2fc)
        to_remove = set()
        for r in proposals:
            p = r["proposal"]
            if p == "merge_left":
                to_remove.add(int(r["seg_bin_start"]))
            elif p == "merge_right":
                to_remove.add(int(r["seg_bin_end"]))
            elif p == "ambiguous_keep":
                if abs(r["srd_left"]) <= abs(r["srd_right"]):
                    to_remove.add(int(r["seg_bin_start"]))
                else:
                    to_remove.add(int(r["seg_bin_end"]))
        if not to_remove:
            break
        new_segs = []
        for s, e in segs:
            if new_segs and s in to_remove:
                new_segs[-1] = (new_segs[-1][0], e)
            else:
                new_segs.append((s, e))
        removed.update(to_remove)
        history.append((it, sorted(to_remove)))
        segs = new_segs
    return segs, removed, history
