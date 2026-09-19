"""dispersion.py — SRD variance diagnostics separate from cell-bootstrap uncertainty.

Plan section 3:
- per-segment Fano factors (mean, variance, usable-bin count, status) from the
  scored mean pseudobulk, per source;
- pooled Pearson dispersion phi_g over segments with positive mean and >=2 usable
  bins:  sum_S sum_b (P-m_S)^2/m_S  /  sum_S (E_S-1);
- residual correlation at physical-bin lags 1..10 using only pairs within the
  same segment at the requested genomic separation (adjacent array rows after
  bin filtering are NOT treated as physically adjacent);
- srd_phi_diagnostic = srd / sqrt(phi_g) for comparison only (not a calibration).

All outputs are descriptive; there is no real-data p-value column.
"""

from __future__ import annotations

import numpy as np


def per_segment_fano(profile, seg_ids):
    """Per-segment mean, variance, Fano, usable-bin count, and status.

    ``profile`` is the per-bin mean pseudobulk; ``seg_ids`` gives each bin's
    segment id (length == len(profile)). Segments with nonpositive mean or fewer
    than 2 usable bins get status and NaN Fano.
    """
    profile = np.asarray(profile, dtype=float)
    seg_ids = np.asarray(seg_ids)
    rows = []
    for sid in np.unique(seg_ids):
        p = profile[seg_ids == sid]
        E = int(len(p))
        m = float(p.mean())
        var = float(p.var())
        if E < 2 or m <= 0:
            rows.append({"segment_id": sid, "E": E, "mean": m, "variance": var,
                         "fano": np.nan, "status": "excluded"})
        else:
            rows.append({"segment_id": sid, "E": E, "mean": m, "variance": var,
                         "fano": var / m, "status": "ok"})
    return rows


def pooled_phi(profile, seg_ids):
    """Pooled Pearson dispersion across segments with positive mean and E>=2."""
    profile = np.asarray(profile, dtype=float)
    seg_ids = np.asarray(seg_ids)
    num = 0.0
    df = 0
    for sid in np.unique(seg_ids):
        p = profile[seg_ids == sid]
        E = len(p)
        if E < 2:
            continue
        m = p.mean()
        if m <= 0:
            continue
        num += float(((p - m) ** 2 / m).sum())
        df += int(E - 1)
    if df == 0:
        return np.nan
    return float(num / df)


def _bin_width(positions):
    d = np.diff(np.sort(np.asarray(positions, dtype=float)))
    d = d[d > 0]
    return float(np.median(d)) if len(d) else 1.0


def residual_lag_correlations(profile, seg_ids, positions, usable_mask=None,
                              max_lag=10):
    """Pearson correlation of within-segment residuals at physical-bin lags.

    Returns ``{lag: (corr, n_pairs)}``. A pair of bins contributes to a lag only
    if both are usable, belong to the same segment, and their genomic separation
    (rounded to bin widths, from ``positions``) equals the lag. Missing/masked
    bins never bridge a lag.
    """
    profile = np.asarray(profile, dtype=float)
    seg_ids = np.asarray(seg_ids)
    positions = np.asarray(positions, dtype=float)
    if usable_mask is None:
        usable_mask = np.ones(len(profile), dtype=bool)
    else:
        usable_mask = np.asarray(usable_mask, dtype=bool)
    bw = _bin_width(positions[usable_mask])
    if not np.isfinite(bw) or bw <= 0:
        return {k: (np.nan, 0) for k in range(1, max_lag + 1)}

    # residuals per segment
    resid = np.full(len(profile), np.nan, dtype=float)
    for sid in np.unique(seg_ids):
        idx = np.where((seg_ids == sid) & usable_mask)[0]
        if len(idx) >= 2:
            m = profile[idx].mean()
            resid[idx] = profile[idx] - m

    pairs = {k: ([], []) for k in range(1, max_lag + 1)}
    for sid in np.unique(seg_ids):
        idx = np.where((seg_ids == sid) & usable_mask)[0]
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                lag = int(round((positions[idx[b]] - positions[idx[a]]) / bw))
                if 1 <= lag <= max_lag:
                    pairs[lag][0].append(resid[idx[a]])
                    pairs[lag][1].append(resid[idx[b]])
    out = {}
    for k in range(1, max_lag + 1):
        xs, ys = pairs[k]
        n = len(xs)
        if n >= 3 and np.std(xs) > 0 and np.std(ys) > 0:
            out[k] = (float(np.corrcoef(xs, ys)[0, 1]), n)
        else:
            out[k] = (np.nan, n)
    return out


def srd_phi_diagnostic(srd, phi):
    """srd / sqrt(phi) when phi is positive and finite; else NaN (comparison only)."""
    phi = float(phi)
    if not np.isfinite(phi) or phi <= 0:
        return np.nan
    return float(srd / np.sqrt(phi))
