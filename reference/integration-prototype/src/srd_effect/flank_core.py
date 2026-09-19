"""core.py — interval sums, rates, flank effects, and SRD.

Plan section-1 definitions (flank-effects-bootstrap-diagnostics):

    q_{c,S}   = sum_{b in S} X_{c,b}                     per-cell per-segment sum
    Y_{g,S}   = (1/n_g) sum_{c in g} q_{c,S}             sum over bins of the mean pseudobulk
    m_{g,S}   = Y_{g,S} / E_S                            corrected mean depth per usable bin
    Delta_L   = log2(m_S / m_L),  Delta_R = log2(m_S / m_R)   (retained SEPARATELY)
    flank_asymmetry = log2(m_L / m_R) = Delta_R - Delta_L      (annotation only)

There is deliberately NO combined log2(Delta_L / Delta_R) effect: it is undefined
when a contrast is zero or signs differ, and it vanishes for a symmetric focal
event (L=10, S=5, R=10 gives Delta_L = Delta_R = -1 yet the combined value is 0).
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import issparse
from scipy.special import xlogy

from srd_pruning import exposure_aware_srd  # reuse the exposure-aware SRD formula


def flank_effects(left, segment, right):
    """Delta_L = log2(segment/left), Delta_R = log2(segment/right)."""
    dl = float(np.log2(segment / left))
    dr = float(np.log2(segment / right))
    return dl, dr


def flank_asymmetry(m_left, m_right):
    """log2(m_left / m_right) = Delta_R - Delta_L; positive when left is higher.

    Annotation only: the target segment cancels, so this is NOT evidence for or
    against a focal segment and must not be used as a merge criterion.
    """
    return float(np.log2(m_left / m_right))


def srd_target(y_target, e_target, y_reference, e_reference):
    """Target-minus-reference SRD on corrected mean profiles.

    Reuses the tested exposure-aware SRD from the SRD spike. Raises on invalid
    exposures and negative/nonfinite inputs; a both-zero pair returns 0.0.
    """
    for val, name in ((y_target, "y_target"), (e_target, "e_target"),
                      (y_reference, "y_reference"), (e_reference, "e_reference")):
        if not np.isfinite(val):
            raise ValueError(f"{name} must be finite, got {val}")
    if e_target <= 0 or e_reference <= 0:
        raise ValueError("exposures must be strictly positive")
    if y_target < 0 or y_reference < 0:
        raise ValueError("depth values must be nonnegative")
    return float(exposure_aware_srd(y_target, e_target, y_reference, e_reference))


def srd_target_with_status(y_target, e_target, y_reference, e_reference):
    """srd_target plus the ``uninformative_zero_pair`` status.

    Both totals zero -> R=0 and uninformative_zero_pair=True; one zero total can
    give a finite SRD. Same input validation as srd_target.
    """
    val = srd_target(y_target, e_target, y_reference, e_reference)
    return val, {"uninformative_zero_pair": bool(y_target == 0 and y_reference == 0)}


def build_cell_segment_sums(X, chrom_of, segments):
    """cells x n_segments matrix of per-cell per-segment sums of ``X``.

    ``X`` is the sparse/dense cells-by-bins matrix (GC-corrected normalized
    depth). Each segment dict needs ``chromosome``, ``seg_bin_start``,
    ``seg_bin_end`` (half-open local bin coordinates on that chromosome).
    Preserves cell order and segment order. Per chromosome the block is
    densified once and a cumsum gives every segment sum in O(n_bins).
    """
    n_cells = X.shape[0]
    q = np.zeros((n_cells, len(segments)), dtype=float)
    by_chr = {}
    for i, seg in enumerate(segments):
        by_chr.setdefault(seg["chromosome"], []).append(i)
    for ch, idxs in by_chr.items():
        cols = np.where(chrom_of == ch)[0]
        block = X[:, cols]
        if issparse(block):
            dens = np.asarray(block.toarray(), dtype=float)
        else:
            dens = np.asarray(block, dtype=float)
        cum = np.cumsum(dens, axis=1)
        n_chr = dens.shape[1]
        for i in idxs:
            lo, hi = int(segments[i]["seg_bin_start"]), int(segments[i]["seg_bin_end"])
            hi = min(hi, n_chr)
            if hi <= 0:
                q[:, i] = 0.0
            else:
                s = cum[:, hi - 1]
                if lo > 0:
                    s = s - cum[:, lo - 1]
                q[:, i] = s
    return q


def segment_exposures(segments, usable_mask):
    """Number of usable bins per segment, from a frozen common bin mask."""
    out = []
    for seg in segments:
        # caller passes segments with local bin coordinates; usable_mask is the
        # per-bin mask aligned to the segment's chromosome slice
        lo, hi = int(seg["seg_bin_start"]), int(seg["seg_bin_end"])
        out.append(int(np.asarray(usable_mask[lo:hi], dtype=bool).sum()))
    return np.asarray(out, dtype=int)


def mean_sums_from_indices(cell_sums, indices):
    """Mean of the sampled cell rows of a cell-by-segment-sums matrix.

    ``cell_sums`` is (n_cells, n_segments); ``indices`` are cell row indices
    (with replacement). Returns the segment-length mean vector. Sampling a whole
    row at once preserves the within-cell dependence across all segments and
    both flanks.
    """
    return np.asarray(cell_sums[np.asarray(indices)], dtype=float).mean(axis=0)


def weighted_all(cluster_mean_sums, sizes):
    """Cell-count-weighted mean of cluster mean-sum vectors (plan section 2).

    ``cluster_mean_sums`` is (n_clusters, n_segments); ``sizes`` are the cluster
    cell counts. The whole-dataset pseudobulk is the size-weighted average, not
    an equal-weight average and not a plain sum.
    """
    sizes = np.asarray(sizes, dtype=float)
    if sizes.sum() <= 0:
        raise ValueError("total cluster size must be positive")
    return np.average(np.asarray(cluster_mean_sums, dtype=float), axis=0, weights=sizes)


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def deviance_terms(y, mu):
    """2 * xlogy(y, y/mu); returns 0 where y==0 (xlogy handles it)."""
    return 2.0 * xlogy(np.asarray(y, dtype=float), np.asarray(y, dtype=float) / np.asarray(mu, dtype=float))
