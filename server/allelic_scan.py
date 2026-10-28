"""Constrained AD breakpoint scan of one chromosome track (spec 2026-10-28-allelic-segmentation-design.md §3).

`validate_candidate_geometry` and `recursive_getbp` are copied from
SNP_project/beta_revised/src/pyEpiAneufinder/get_breakpoints.py, the function
run_bootstrap_diagnostics.py used to write raw_accepted_breakpoints.tsv. Only
the recursion and the geometry check are copied; the statistics are pyEpi's
own, called through pyepi_adapter.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from . import pyepi_adapter as sci
from .config import Config

DEFAULTS = {'k': 6, 'n_permutations': 1000, 'alpha': 0.001, 'min_child_bins': 20, 'min_terminal_bins': 50}


def _ensure(package_dir: Optional[str] = None) -> None:
    if not sci._loaded:
        sci.load(package_dir or Config.from_env().package_dir)


def validate_candidate_geometry(local_cut, node_start, node_end, chromosome_start, chromosome_end,
                                min_child_bins, min_terminal_bins):
    """Copied verbatim (logic) from beta_revised get_breakpoints.py; terminal rules apply only to
    children touching the true chromosome start/end."""
    absolute_cut = node_start + local_cut
    left_child_bins = local_cut
    right_child_bins = (node_end - node_start) - local_cut
    reasons = []
    if left_child_bins < min_child_bins:
        reasons.append("left_child_too_short")
    if right_child_bins < min_child_bins:
        reasons.append("right_child_too_short")
    if node_start == chromosome_start and (absolute_cut - chromosome_start) < min_terminal_bins:
        reasons.append("left_terminal_too_short")
    if node_end == chromosome_end and (chromosome_end - absolute_cut) < min_terminal_bins:
        reasons.append("right_terminal_too_short")
    return (len(reasons) == 0, tuple(reasons))


def recursive_getbp(seq_data, seq_data_cell, k=3, depth=0, node_start=0, n_permutations=500,
                    alpha=0.05, chromosome_start=None, chromosome_end=None,
                    min_child_bins=20, min_terminal_bins=50, audit=None):
    """Copied from beta_revised get_breakpoints.py. Candidates in descending AD order; the first
    geometry-valid one is tested; if it fails, the node stops. Every examined candidate is audited."""
    if depth >= k or len(seq_data) <= 1:
        return []
    if chromosome_start is None:
        chromosome_start = 0
    if chromosome_end is None:
        chromosome_end = len(seq_data)
    if audit is None:
        audit = []
    dist_vector = np.array(sci.seq_dist_ad(seq_data))
    if dist_vector.size == 0:
        return []
    n = len(seq_data)
    node_end = node_start + n
    order = np.argsort(-dist_vector, kind="stable") + 1
    for rank, i in enumerate(order, start=1):
        i = int(i)
        abs_cut = node_start + i
        ad_score = float(dist_vector[i - 1])
        geometry_valid, reasons = validate_candidate_geometry(
            i, node_start, node_end, chromosome_start, chromosome_end, min_child_bins, min_terminal_bins)
        row = {
            "source": None, "chromosome": None, "recursion_depth": depth,
            "node_start_bin": node_start, "node_end_bin": node_end,
            "candidate_rank": rank, "candidate_local_cut": i, "candidate_absolute_bin": abs_cut,
            "ad_score": ad_score, "left_child_bins": i, "right_child_bins": n - i,
            "geometry_valid": geometry_valid, "rejection_reasons": reasons,
            "selected_for_significance": False, "p_local": np.nan, "p_global": np.nan, "accepted": False,
        }
        audit.append(row)
        if not geometry_valid:
            continue
        row["selected_for_significance"] = True
        bp_dist = dist_vector[i - 1]
        left_seg = seq_data[:i]
        right_seg = seq_data[i:]
        pval = sci.permutation_test_ad(left_seg, right_seg, bp_dist, n_permutations)
        row["p_local"] = pval
        if pval >= alpha:
            return []
        pval_global = sci.global_permutation_test_ad(seq_data_cell, len(left_seg), len(right_seg), bp_dist, n_permutations)
        row["p_global"] = pval_global
        accepted = (pval < alpha) and (pval_global < alpha)
        row["accepted"] = accepted
        if not accepted:
            return []
        bp_global = node_start + i
        left_bps = recursive_getbp(left_seg, seq_data_cell, k=k, depth=depth + 1, node_start=node_start,
                                   n_permutations=n_permutations, alpha=alpha,
                                   chromosome_start=chromosome_start, chromosome_end=chromosome_end,
                                   min_child_bins=min_child_bins, min_terminal_bins=min_terminal_bins, audit=audit)
        right_bps = recursive_getbp(right_seg, seq_data_cell, k=k, depth=depth + 1, node_start=node_start + i,
                                    n_permutations=n_permutations, alpha=alpha,
                                    chromosome_start=chromosome_start, chromosome_end=chromosome_end,
                                    min_child_bins=min_child_bins, min_terminal_bins=min_terminal_bins, audit=audit)
        return left_bps + [(bp_global, bp_dist, pval, pval_global)] + right_bps
    audit.append({
        "source": None, "chromosome": None, "recursion_depth": depth,
        "node_start_bin": node_start, "node_end_bin": node_end,
        "candidate_rank": None, "candidate_local_cut": None, "candidate_absolute_bin": None,
        "ad_score": np.nan, "left_child_bins": None, "right_child_bins": None,
        "geometry_valid": False, "rejection_reasons": ("no_valid_candidate",),
        "selected_for_significance": False, "p_local": np.nan, "p_global": np.nan, "accepted": False,
    })
    return []


def _plain(v):
    """JSON-safe scalar: numpy → Python, NaN → None."""
    if v is None:
        return None
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return None if not np.isfinite(v) else float(v)
    return v


def scan_chromosome(values, kept, settings, ha=None, hb=None, package_dir=None, pool=None) -> dict:
    """Scan one chromosome on its kept bins; map cuts back to chromosome-local bins.

    `pool` is what the global permutation test samples from: original pyEpi passes the
    whole cell genome-wide (pyEpiAneufinder.py, `counts.X[i, :]`); None uses this
    chromosome's kept values, as run_bootstrap_diagnostics.py did for the depth breakpoints.

    A cut between kept bins i and i+1 sits at the bin of kept bin i+1; skipped bins
    between them are reported as [gap_start, gap_end]. Segment sums and stats cover kept bins.
    """
    _ensure(package_dir)
    values = np.asarray(values, dtype=np.float64)
    kept = np.asarray(kept, dtype=bool)
    n = len(values)
    idx = np.flatnonzero(kept)
    seq = values[idx]
    audit: list = []
    accepted = []
    if len(seq) > 1:
        cell = seq if pool is None else np.asarray(pool, dtype=np.float64)
        accepted = recursive_getbp(seq, cell, k=int(settings['k']), depth=0, node_start=0,
                                   n_permutations=int(settings['n_permutations']), alpha=float(settings['alpha']),
                                   chromosome_start=0, chromosome_end=len(seq),
                                   min_child_bins=int(settings['min_child_bins']),
                                   min_terminal_bins=int(settings['min_terminal_bins']), audit=audit)
    rows = []
    for r in audit:
        cut = r['candidate_absolute_bin']
        row = {key: _plain(v) for key, v in r.items() if key != 'rejection_reasons'}
        row['rejection_reasons'] = ';'.join(r['rejection_reasons'])
        row['bin'] = int(idx[cut]) if cut is not None and 0 <= cut < len(idx) else None
        rows.append(row)
    depth_of = {r['candidate_absolute_bin']: r['recursion_depth'] for r in audit if r['accepted']}
    bps = []
    for cut, ad, pl, pg in accepted:
        b, prev = int(idx[cut]), int(idx[cut - 1])
        gap = (prev + 1, b) if b - prev > 1 else (None, None)
        bps.append({'bin': b, 'gap_start': gap[0], 'gap_end': gap[1], 'recursion_depth': int(depth_of[cut]),
                    'ad_score': float(ad), 'p_local': float(pl), 'p_global': float(pg)})
    bounds = [0] + [c for c, *_ in accepted] + [len(idx)]
    segs = []
    for j in range(len(bounds) - 1):
        k0, k1 = bounds[j], bounds[j + 1]
        b0 = 0 if j == 0 else int(idx[k0])
        b1 = n - 1 if j == len(bounds) - 2 else int(idx[k1]) - 1
        part, kb = seq[k0:k1], idx[k0:k1]
        row = {'bin_start': b0, 'bin_end': b1, 'n_bins': b1 - b0 + 1, 'kept_bins': int(k1 - k0),
               'mean': float(part.mean()) if len(part) else None,
               'median': float(np.median(part)) if len(part) else None}
        if ha is not None and hb is not None:
            hs, bs = float(np.asarray(ha)[kb].sum()), float(np.asarray(hb)[kb].sum())
            row.update({'ha_sum': hs, 'hb_sum': bs, 'at_sum': hs + bs})
        segs.append(row)
    return {'audit': rows, 'breakpoints': bps, 'segments': segs}


def scan_task(task):
    """Process-pool worker: (source, chrom, values, kept, ha, hb, settings, package_dir, pool)."""
    source, chrom, values, kept, ha, hb, settings, package_dir, pool = task
    return source, chrom, scan_chromosome(values, kept, settings, ha=ha, hb=hb, package_dir=package_dir, pool=pool)
