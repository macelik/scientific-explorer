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


def _finite(v):
    return v is not None and np.isfinite(v)


def classify_transition(left_score, right_score):
    """True iff both scores are finite and have strictly opposite signs.

    Zero never counts as an opposite sign; a missing/nonfinite score never
    counts as a transition. Descriptive only, not a merge decision.
    """
    if not (_finite(left_score) and _finite(right_score)):
        return False
    return (left_score * right_score) < 0


def weaker_side(left_score, right_score):
    """Side with the smaller absolute score. 'tie' on exact equality (not
    treated as evidence either way); None if either side is missing."""
    if left_score is None or right_score is None:
        return None
    al, ar = abs(left_score), abs(right_score)
    if al == ar:
        return 'tie'
    return 'left' if al < ar else 'right'


def classify_merge_proposal(left_score, right_score, threshold):
    """Eligible sides under abs(score) < threshold (strict), and the weaker
    eligible side when both qualify. A terminal segment's missing side is
    never invented as eligible."""
    el = _finite(left_score) and abs(left_score) < threshold
    er = _finite(right_score) and abs(right_score) < threshold
    side = None
    if el and er:
        w = weaker_side(left_score, right_score)
        side = 'both' if w == 'tie' else w
    elif el:
        side = 'left'
    elif er:
        side = 'right'
    return {'eligible_left': bool(el), 'eligible_right': bool(er), 'weaker_side': side}


def classify_ambiguous(srd_left, srd_right, fc_left, fc_right):
    """'ambiguous' iff the weaker SRD-variant flank differs from the weaker
    delta-log2FC flank; 'consistent' if they agree; 'no_unique_preference' on
    an exact tie in either metric; 'undefined' if any input is missing or
    nonfinite. Requires an internal segment (both sides present by construction)."""
    vals = (srd_left, srd_right, fc_left, fc_right)
    if any(v is None or not np.isfinite(v) for v in vals):
        return 'undefined'
    srd_side, fc_side = weaker_side(srd_left, srd_right), weaker_side(fc_left, fc_right)
    if srd_side == 'tie' or fc_side == 'tie':
        return 'no_unique_preference'
    return 'consistent' if srd_side == fc_side else 'ambiguous'


def hatch_scores(profile, boundaries, n_bins, var_start, var_end, chrom_offset):
    """Per internal segment: SRD, SRD/sqrt(phi) and delta-log2FC (all three
    estimators) against both flanks, plus the shared pooled phi and per-side
    genomic-gap flags. Computed genome-wide from live boundaries/profile, not
    restricted to chromosomes with archived diagnostic tables."""
    segs = flat_partition(boundaries, n_bins)
    summaries = [segment_summary(profile, s, e) for s, e in segs]
    phi = pooled_phi(profile, segs)
    rows = []
    for i in range(1, len(segs) - 1):
        seg, left, right = summaries[i], summaries[i - 1], summaries[i + 1]
        row = {
            'start': seg['start'], 'end': seg['end'], 'phi': phi,
            'gap_left': is_gap(var_start, var_end, chrom_offset + segs[i - 1][1] - 1, chrom_offset + segs[i][0]),
            'gap_right': is_gap(var_start, var_end, chrom_offset + segs[i][1] - 1, chrom_offset + segs[i + 1][0]),
        }
        for m in ('srd', 'srd_phi'):
            row[f'{m}_left'] = flank_score(m, seg, left, phi)
            row[f'{m}_right'] = flank_score(m, seg, right, phi)
        for est in ESTIMATORS:
            row[f'delta_{est}_left'] = delta_log2fc(seg, left, est)
            row[f'delta_{est}_right'] = delta_log2fc(seg, right, est)
        rows.append(row)
    return {'segments': [{'start': s, 'end': e} for s, e in segs], 'phi': phi, 'rows': rows}


def _segment_status(idx, segs, summaries, phi, metric, estimator):
    """Transition/ambiguous status of the internal segment at idx, or None for
    a terminal segment (fewer than two neighbors)."""
    if idx <= 0 or idx >= len(segs) - 1:
        return None
    seg, left, right = summaries[idx], summaries[idx - 1], summaries[idx + 1]
    l_score = flank_score(metric, seg, left, phi, estimator)
    r_score = flank_score(metric, seg, right, phi, estimator)
    srd_l = flank_score('srd_phi', seg, left, phi)
    srd_r = flank_score('srd_phi', seg, right, phi)
    fc_l = delta_log2fc(seg, left, estimator)
    fc_r = delta_log2fc(seg, right, estimator)
    return {
        'transition': classify_transition(l_score, r_score),
        'ambiguous': classify_ambiguous(srd_l, srd_r, fc_l, fc_r) == 'ambiguous',
    }


def _vetoed(i, segs, summaries, phi, metric, estimator, veto_transition, veto_ambiguous):
    if not (veto_transition or veto_ambiguous):
        return False
    for idx in (i, i + 1):
        status = _segment_status(idx, segs, summaries, phi, metric, estimator)
        if status is None:
            continue
        if veto_transition and status['transition']:
            return True
        if veto_ambiguous and status['ambiguous']:
            return True
    return False


def _eligible_boundaries(segs, summaries, phi, metric, estimator, threshold,
                          veto_transition, veto_ambiguous, small_max_bins,
                          allow_gap_crossing, var_start, var_end, chrom_offset):
    out = []
    for i in range(len(segs) - 1):
        left_s, left_e = segs[i]
        right_s, right_e = segs[i + 1]
        if not allow_gap_crossing and is_gap(var_start, var_end, chrom_offset + left_e - 1, chrom_offset + right_s):
            continue
        if small_max_bins is not None and (left_e - left_s) > small_max_bins and (right_e - right_s) > small_max_bins:
            continue
        score = flank_score(metric, summaries[i], summaries[i + 1], phi, estimator)
        if score is None or not np.isfinite(score) or abs(score) >= threshold:
            continue
        if _vetoed(i, segs, summaries, phi, metric, estimator, veto_transition, veto_ambiguous):
            continue
        out.append((i, score))
    return out


def run_merge(profile, boundaries, n_bins, var_start, var_end, chrom_offset, metric, threshold,
              estimator='mean', veto_transition=False, veto_ambiguous=False, small_max_bins=None,
              allow_gap_crossing=False, max_steps=10000):
    """Greedy exploratory adjacent-segment merge, operating on one source's
    one-chromosome profile/boundaries. Repeatedly merges the eligible boundary
    (abs(score) < threshold, strict) with the smallest absolute score,
    tie-breaking on the leftmost boundary bin, recomputing every score
    (including the shared pooled phi) after each merge, until none remain
    eligible. Transition/ambiguous vetoes are off by default and explicit when
    set; a genomic-gap boundary is never crossed unless allow_gap_crossing is
    set. This is a new exploratory rule, independent of the original
    bootstrap/chain-guard merge prototype.
    """
    if metric not in METRICS:
        raise ValueError(f'Unknown merge metric {metric!r}')
    if estimator not in ESTIMATORS:
        raise ValueError(f'Unknown log2FC estimator {estimator!r}')
    if not (isinstance(threshold, (int, float)) and np.isfinite(threshold) and threshold > 0):
        raise ValueError('Threshold must be a positive finite number')
    segs = flat_partition(boundaries, n_bins)
    original = [{'start': s, 'end': e} for s, e in segs]
    steps = [{'step': 0, 'segments': list(original), 'removed_boundary': None,
              'merged_interval': None, 'score': None, 'phi': None}]
    guard = 0
    while True:
        guard += 1
        if guard > max_steps:
            raise RuntimeError('Merge did not terminate within max_steps')
        summaries = [segment_summary(profile, s, e) for s, e in segs]
        phi = pooled_phi(profile, segs) if metric != 'delta_log2fc' else None
        candidates = _eligible_boundaries(segs, summaries, phi, metric, estimator, threshold,
                                           veto_transition, veto_ambiguous, small_max_bins,
                                           allow_gap_crossing, var_start, var_end, chrom_offset)
        if not candidates:
            break
        i, score = min(candidates, key=lambda t: (abs(t[1]), segs[t[0]][1]))
        removed_boundary = segs[i][1]
        merged = (segs[i][0], segs[i + 1][1])
        segs = segs[:i] + [merged] + segs[i + 2:]
        steps.append({
            'step': len(steps), 'segments': [{'start': s, 'end': e} for s, e in segs],
            'removed_boundary': removed_boundary, 'merged_interval': [merged[0], merged[1]],
            'score': score, 'phi': phi,
        })
    return {
        'original': original, 'final': steps[-1]['segments'], 'steps': steps,
        'settings': {'metric': metric, 'estimator': estimator, 'threshold': threshold,
                     'veto_transition': veto_transition, 'veto_ambiguous': veto_ambiguous,
                     'small_max_bins': small_max_bins, 'allow_gap_crossing': allow_gap_crossing},
    }
