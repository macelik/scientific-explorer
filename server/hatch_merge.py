"""New exploratory hatching and adjacent-segment merging for the segmentation
tab's Chromosome tracks. Deliberately separate from the original bootstrap/
chain-guard merge prototype (analysis/depth_segment_prototype/merge_core.py)
and from the archived SRD screen: this is a new, independent, exploratory
rule, not a replacement or recomputation of either.

Level estimators (mean/median/two-sided-IQR-mean) are reused unmodified from
flank_view.summarize. The signed-root-deviance and pooled-Pearson-dispersion
formulas are ported from reference/integration-prototype/src/srd_effect
(srd_pruning.exposure_aware_srd, flank_dispersion.pooled_phi); tests in
test_hatch_merge.py assert numerical agreement with those reference functions
on shared examples so the two never silently drift apart.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.special import xlogy

from .flank_view import summarize

METRICS = ('srd', 'srd_phi', 'delta_log2fc')
ESTIMATORS = ('mean', 'median', 'iqr_mean')


def flat_partition(boundaries, n_bins):
    """Half-open chromosome-local segments from internal boundary bin offsets."""
    cuts = [0, *sorted(int(b) for b in boundaries), int(n_bins)]
    return [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1) if cuts[i + 1] > cuts[i]]


def segment_summary(profile, s, e):
    """summarize() over profile[s:e], plus start/end/sum (sum needed for SRD)."""
    x = np.asarray(profile, dtype=float)[s:e]
    return {**summarize(x), 'start': int(s), 'end': int(e), 'sum': float(x.sum())}


def signed_srd(y_t, e_t, y_r, e_r):
    """Signed root deviance between two Poisson-like exposures.

    Ported from reference/integration-prototype/src/srd_effect/srd_pruning.py
    (exposure_aware_srd): D = 2*[xlogy(y_t,y_t/mu_t) + xlogy(y_r,y_r/mu_r)],
    lambda_hat = (y_t+y_r)/(e_t+e_r), mu_* = e_* * lambda_hat. Both-zero pair
    returns 0.0. Raises on nonpositive exposure or negative depth.
    """
    if e_t <= 0 or e_r <= 0:
        raise ValueError('exposure must be strictly positive')
    if y_t < 0 or y_r < 0:
        raise ValueError('depth values must be nonnegative')
    y_sum = y_t + y_r
    if y_sum == 0:
        return 0.0
    lam = y_sum / (e_t + e_r)
    mu_t, mu_r = e_t * lam, e_r * lam
    d = 2.0 * (xlogy(y_t, y_t / mu_t) + xlogy(y_r, y_r / mu_r))
    return float(np.sign(y_t / e_t - y_r / e_r) * np.sqrt(max(float(d), 0.0)))


def pooled_phi(profile, seg_bounds):
    """Pooled Pearson dispersion over segments with >=2 bins and mean>0.

    Ported from reference/integration-prototype/src/srd_effect/flank_dispersion.py
    (pooled_phi): sum_S sum_b (p-m_S)^2/m_S / sum_S (E_S-1), taking (start,end)
    bounds instead of a seg_ids array. This is pooled across EVERY segment of
    the current partition for one source/chromosome, so it must be recomputed
    whenever that partition changes (see run_merge below): SRD/sqrt(phi) is
    not a per-segment-pair quantity.
    """
    profile = np.asarray(profile, dtype=float)
    num, df = 0.0, 0
    for s, e in seg_bounds:
        p = profile[s:e]
        if len(p) < 2:
            continue
        m = p.mean()
        if m <= 0:
            continue
        num += float(((p - m) ** 2 / m).sum())
        df += len(p) - 1
    return float(num / df) if df > 0 else None


def estimator_level(summary, estimator):
    if estimator not in ESTIMATORS:
        raise ValueError(f'Unknown estimator {estimator!r}')
    return summary[estimator if estimator != 'iqr_mean' else 'iqr_mean']


def delta_log2fc(seg_summary, flank_summary, estimator):
    """Signed log2(segment level / flank level) for the selected estimator.

    Reuses flank_view.summarize's estimators via estimator_level. No
    pseudocount: undefined (None) whenever either level is not > 0.
    """
    a, b = estimator_level(seg_summary, estimator), estimator_level(flank_summary, estimator)
    return float(np.log2(a / b)) if a > 0 and b > 0 else None


def flank_score(metric, seg_summary, flank_summary, phi, estimator='mean'):
    """Signed score of seg relative to flank for the selected metric.

    metric='delta_log2fc' ignores phi. metric='srd' returns the signed root
    deviance directly. metric='srd_phi' divides by sqrt(phi) (never by phi
    itself) and is None when phi is missing, nonfinite or non-positive.
    """
    if metric not in METRICS:
        raise ValueError(f'Unknown metric {metric!r}')
    if metric == 'delta_log2fc':
        return delta_log2fc(seg_summary, flank_summary, estimator)
    srd = signed_srd(seg_summary['sum'], seg_summary['end'] - seg_summary['start'],
                      flank_summary['sum'], flank_summary['end'] - flank_summary['start'])
    if metric == 'srd':
        return srd
    if phi is None or not np.isfinite(phi) or phi <= 0:
        return None
    return float(srd / np.sqrt(phi))


def is_gap(var_start, var_end, left_end_bin, right_start_bin):
    """True if the two retained bins straddling a boundary are not genomically
    contiguous (a filtered/omitted interval sits between them)."""
    return float(var_start[right_start_bin]) != float(var_end[left_end_bin])
