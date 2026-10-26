from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from server import hatch_merge as hm


def test_flat_partition_half_open_and_drops_empty():
    assert hm.flat_partition([3, 3, 7], 10) == [(0, 3), (3, 7), (7, 10)]
    assert hm.flat_partition([], 5) == [(0, 5)]


def test_segment_summary_adds_start_end_sum_to_summarize():
    profile = np.array([2., 2., 4., 4.])
    s = hm.segment_summary(profile, 1, 3)
    assert s['start'] == 1 and s['end'] == 3
    assert s['sum'] == pytest.approx(6.0)
    assert s['mean'] == pytest.approx(3.0)


def test_signed_srd_matches_reference_examples():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'reference' / 'integration-prototype' / 'src' / 'srd_effect'))
    from srd_pruning import exposure_aware_srd  # reference-only import, test comparison only
    cases = [(10., 5., 10., 5.), (20., 5., 5., 5.), (0., 4., 0., 4.), (7., 3., 2., 9.), (5., 5., 20., 5.)]
    for y_t, e_t, y_r, e_r in cases:
        assert hm.signed_srd(y_t, e_t, y_r, e_r) == pytest.approx(exposure_aware_srd(y_t, e_t, y_r, e_r))


def test_signed_srd_zero_pair_and_validation():
    assert hm.signed_srd(0., 4., 0., 4.) == 0.0
    with pytest.raises(ValueError):
        hm.signed_srd(1., 0., 1., 4.)
    with pytest.raises(ValueError):
        hm.signed_srd(-1., 4., 1., 4.)


def test_pooled_phi_matches_reference_and_handles_short_segments():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'reference' / 'integration-prototype' / 'src' / 'srd_effect'))
    from flank_dispersion import pooled_phi as ref_pooled_phi
    rng = np.random.default_rng(0)
    profile = rng.poisson(5, size=40).astype(float)
    bounds = [(0, 5), (5, 6), (6, 20), (20, 40)]  # one segment with only 1 bin
    seg_ids = np.zeros(40, dtype=int)
    for i, (s, e) in enumerate(bounds):
        seg_ids[s:e] = i
    assert hm.pooled_phi(profile, bounds) == pytest.approx(ref_pooled_phi(profile, seg_ids))


def test_pooled_phi_none_when_no_segment_qualifies():
    assert hm.pooled_phi(np.array([0., 0., 0.]), [(0, 1), (1, 3)]) is None


def test_delta_log2fc_undefined_on_zero_level_no_pseudocode():
    seg = hm.segment_summary(np.array([0., 0.]), 0, 2)
    flank = hm.segment_summary(np.array([2., 2.]), 0, 2)
    assert hm.delta_log2fc(seg, flank, 'mean') is None
    flank2 = hm.segment_summary(np.array([4., 4.]), 0, 2)
    assert hm.delta_log2fc(hm.segment_summary(np.array([2., 2.]), 0, 2), flank2, 'mean') == pytest.approx(-1.0)


def test_flank_score_srd_phi_none_when_phi_undefined():
    seg = hm.segment_summary(np.array([5., 5.]), 0, 2)
    flank = hm.segment_summary(np.array([2., 2.]), 0, 2)
    assert hm.flank_score('srd_phi', seg, flank, None) is None
    assert hm.flank_score('srd_phi', seg, flank, 0.0) is None
    assert hm.flank_score('srd', seg, flank, None) is not None


def test_is_gap_detects_bp_discontinuity():
    var_start = np.array([0, 100, 250])
    var_end = np.array([100, 200, 350])
    assert hm.is_gap(var_start, var_end, 1, 2) is True   # 200 != 250
    assert hm.is_gap(var_start, var_end, 0, 1) is False  # 100 == 100


def test_is_gap_one_based_inclusive_coordinates():
    # pyEpiAneufinder bins: [800001, 900000], [900001, 1000000], then a 100 kb hole
    var_start = np.array([800001, 900001, 1100001], dtype=float)
    var_end = np.array([900000, 1000000, 1200000], dtype=float)
    assert hm.is_gap(var_start, var_end, 0, 1) is False
    assert hm.is_gap(var_start, var_end, 1, 2) is True


def test_classify_transition_opposite_same_zero_missing():
    assert hm.classify_transition(1.0, -1.0) is True
    assert hm.classify_transition(1.0, 1.0) is False
    assert hm.classify_transition(0.0, -1.0) is False    # zero is not an opposite sign
    assert hm.classify_transition(None, -1.0) is False
    assert hm.classify_transition(float('nan'), 1.0) is False


def test_weaker_side_uses_absolute_magnitude_when_both_negative():
    assert hm.weaker_side(-1.0, -3.0) == 'left'
    assert hm.weaker_side(-3.0, -1.0) == 'right'
    assert hm.weaker_side(-2.0, 2.0) == 'tie'
    assert hm.weaker_side(None, 1.0) is None


def test_classify_merge_proposal_strict_threshold_and_sides():
    assert hm.classify_merge_proposal(0.5, 5.0, 1.0) == {'eligible_left': True, 'eligible_right': False, 'weaker_side': 'left'}
    assert hm.classify_merge_proposal(1.0, 5.0, 1.0)['eligible_left'] is False  # strict <, not <=
    assert hm.classify_merge_proposal(0.5, 0.5, 1.0)['weaker_side'] == 'both'
    assert hm.classify_merge_proposal(None, 5.0, 1.0) == {'eligible_left': False, 'eligible_right': False, 'weaker_side': None}


def test_classify_merge_proposal_terminal_segment_one_neighbor():
    # a terminal segment only ever supplies one side; the other stays ineligible/None, never invented
    assert hm.classify_merge_proposal(0.2, None, 1.0) == {'eligible_left': True, 'eligible_right': False, 'weaker_side': 'left'}


def test_classify_ambiguous_differs_ties_and_undefined():
    assert hm.classify_ambiguous(1.0, 3.0, 3.0, 1.0) == 'ambiguous'      # weaker SRD=left, weaker log2fc=right
    assert hm.classify_ambiguous(1.0, 3.0, 1.0, 3.0) == 'consistent'
    assert hm.classify_ambiguous(2.0, 2.0, 1.0, 3.0) == 'no_unique_preference'  # exact SRD tie
    assert hm.classify_ambiguous(1.0, 3.0, None, 3.0) == 'undefined'
    assert hm.classify_ambiguous(1.0, float('nan'), 1.0, 3.0) == 'undefined'


def _toy_profile_and_bounds():
    # chr-local bins 0..30; boundaries at 10 and 20 give three segments
    profile = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)])
    var_start = (np.arange(30) * 100).astype(float)
    var_end = var_start + 100
    return profile, [10, 20], 30, var_start, var_end


def test_hatch_scores_reports_every_internal_segment_with_phi_and_gap_flags():
    profile, boundaries, n_bins, var_start, var_end = _toy_profile_and_bounds()
    out = hm.hatch_scores(profile, boundaries, n_bins, var_start, var_end, chrom_offset=0)
    assert [tuple(s.values()) for s in out['segments']] == [(0, 10), (10, 20), (20, 30)]
    assert len(out['rows']) == 1  # only the middle segment has both flanks
    row = out['rows'][0]
    assert row['start'] == 10 and row['end'] == 20
    assert row['gap_left'] is False and row['gap_right'] is False
    assert row['phi'] == out['phi']
    assert row['delta_mean_left'] == pytest.approx(np.log2(4.2 / 4.0))
    assert row['srd_left'] is not None and row['srd_phi_right'] is not None


def test_hatch_scores_flags_a_genomic_gap():
    profile, boundaries, n_bins, var_start, var_end = _toy_profile_and_bounds()
    var_end = var_end.copy()
    var_end[9] += 500  # open a gap between bin 9 and bin 10 (the first boundary)
    out = hm.hatch_scores(profile, boundaries, n_bins, var_start, var_end, chrom_offset=0)
    assert out['rows'][0]['gap_left'] is True


from fastapi import FastAPI
from fastapi.testclient import TestClient


class _FakeStore:
    """Minimal stand-in for IntegrationStore covering what the routes touch."""
    def __init__(self):
        self.identity = {'hash': 'fake', 'files': []}
        self.var_seq = np.array(['chr1'] * 30)
        self.var_start = (np.arange(30) * 100).astype(float)
        self.var_end = self.var_start + 100
        self.members = {'cluster0': np.arange(5), 'all': np.arange(5)}
        pb = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)])
        self.arrays = {'pb_x': np.tile(pb, (2, 1))}
        import pandas as pd
        self.tables = {'breakpoints': pd.DataFrame({'source': ['cluster0', 'cluster0'],
                                                      'chromosome': ['chr1', 'chr1'],
                                                      'absolute_bin': [10, 20]})}
        self.available = True


def _client():
    from server.hatch_merge_api import register_hatch_merge
    app = FastAPI()
    store = _FakeStore()
    register_hatch_merge(app, lambda: store)
    return TestClient(app)


def test_hatch_scores_route_returns_rows_and_bp_arrays():
    resp = _client().post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'cluster0'})
    assert resp.status_code == 200
    body = resp.json()
    assert body['source'] == 'cluster0' and body['chrom'] == 'chr1'
    assert len(body['start_bp']) == 30
    assert len(body['rows']) == 1


def test_hatch_scores_route_rejects_unknown_source():
    resp = _client().post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'nope'})
    assert resp.status_code == 422
