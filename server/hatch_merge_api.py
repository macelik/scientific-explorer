"""New exploratory hatching/merge routes. Kept separate from
integration_api.py's flank-view/bootstrap routes, which back the existing
(unmodified) Alternative flank view."""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Optional

import numpy as np
from fastapi import HTTPException
from pydantic import BaseModel, Field

from . import hatch_merge as hm
from .integration import IntegrationStore


class HatchScoresRequest(BaseModel):
    chrom: str
    source: str


class HatchMergeRequest(BaseModel):
    chrom: str
    source: str
    metric: str
    threshold: float
    estimator: str = 'mean'
    veto_transition: bool = False
    veto_ambiguous: bool = False
    small_max_bins: Optional[int] = Field(default=None, ge=1)
    allow_gap_crossing: bool = False


def _resolve(ig: IntegrationStore, chrom: str, source: str):
    if ig.arrays is None:
        raise HTTPException(503, 'Pseudobulk index is not ready')
    if source not in ig.members:
        raise HTTPException(422, f'Unknown source {source!r}')
    idx = np.flatnonzero(ig.var_seq == chrom)
    if not len(idx):
        raise HTTPException(422, f'Unknown chromosome {chrom!r}')
    if not np.array_equal(idx, np.arange(idx[0], idx[-1] + 1)):
        raise HTTPException(422, 'Chromosome bins are not contiguous')
    si = list(ig.members).index(source)  # ig.members is a dict {source: cell-index-array}; list() gives its keys in insertion order
    profile = ig.arrays['pb_x'][si, idx]
    bp = ig.tables['breakpoints']
    boundaries = bp.loc[(bp.source == source) & (bp.chromosome == chrom), 'absolute_bin'].tolist()
    return profile, boundaries, idx


def register_hatch_merge(app, get_integration):
    def _ig() -> IntegrationStore:
        ig = get_integration()
        if ig is None or not getattr(ig, 'available', False):
            raise HTTPException(503, 'integration data not loaded')
        return ig

    @app.post('/api/integration/hatch-scores')
    def hatch_scores(req: HatchScoresRequest):
        ig = _ig()
        try:
            profile, boundaries, idx = _resolve(ig, req.chrom, req.source)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e
        result = hm.hatch_scores(profile, boundaries, len(idx), ig.var_start, ig.var_end, chrom_offset=int(idx[0]))
        return {**result, 'source': req.source, 'chrom': req.chrom,
                'start_bp': ig.var_start[idx].tolist(), 'end_bp': ig.var_end[idx].tolist()}

    @app.post('/api/integration/hatch-merge')
    def hatch_merge(req: HatchMergeRequest):
        ig = _ig()
        try:
            profile, boundaries, idx = _resolve(ig, req.chrom, req.source)
            result = hm.run_merge(profile, boundaries, len(idx), ig.var_start, ig.var_end, chrom_offset=int(idx[0]),
                                   metric=req.metric, threshold=req.threshold, estimator=req.estimator,
                                   veto_transition=req.veto_transition, veto_ambiguous=req.veto_ambiguous,
                                   small_max_bins=req.small_max_bins, allow_gap_crossing=req.allow_gap_crossing)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e
        return {**result, 'source': req.source, 'chrom': req.chrom,
                'start_bp': ig.var_start[idx].tolist(), 'end_bp': ig.var_end[idx].tolist(),
                'provenance': {'dataset': ig.identity,
                               'implementation_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                               'scope': 'New exploratory merge rule; independent of the original bootstrap/chain-guard '
                                        'merge prototype and the archived SRD screen. No CN recomputation on merged boundaries.'}}
