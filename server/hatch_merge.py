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
