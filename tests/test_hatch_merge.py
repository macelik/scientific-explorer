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
    cases = [(10., 5., 10., 5.), (20., 5., 5., 5.), (0., 4., 0., 4.), (7., 3., 2., 9.)]
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
