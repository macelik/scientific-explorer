from pathlib import Path
import json
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from server import hatch_merge as hm
from server import merge_rules as mr

FIXTURE = Path(__file__).resolve().parent / 'fixtures' / 'merge_rules_parity.json'


def _row(start=10, end=20, gap_left=False, gap_right=False, **scores):
    row = {'start': start, 'end': end, 'phi': 1.0, 'gap_left': gap_left, 'gap_right': gap_right}
    for key in ('srd', 'srd_phi', 'delta_mean', 'delta_median', 'delta_iqr_mean'):
        row[f'{key}_left'] = None
        row[f'{key}_right'] = None
    row.update(scores)
    return row


def _rule(rid, conditions, merge_metric='srd_phi', merge_estimator='mean', enabled=True, pattern='/'):
    return {'id': rid, 'name': rid, 'enabled': enabled, 'conditions': conditions,
            'mergeInto': {'metric': merge_metric, 'estimator': merge_estimator},
            'hatch': {'pattern': pattern, 'color': '#000000'}}


WEAK_SRD_PHI = {'kind': 'weak_flank', 'metric': 'srd_phi', 'estimator': 'mean', 'threshold': 1.0}
GOOD = _rule('ok', [WEAK_SRD_PHI])


def test_score_rows_matches_hatch_scores_rows():
    profile = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)])
    var_start = (np.arange(30) * 100).astype(float)
    var_end = var_start + 100
    phi, rows = hm.score_rows(profile, hm.flat_partition([10, 20], 30), var_start, var_end, 0)
    full = hm.hatch_scores(profile, [10, 20], 30, var_start, var_end, 0)
    assert rows == full['rows']
    assert phi == full['phi']


def test_weak_flank_is_strict_and_needs_both_flanks():
    row = _row(srd_phi_left=0.5, srd_phi_right=-3.0)
    assert mr.condition_holds(WEAK_SRD_PHI, row) is True
    assert mr.condition_holds({**WEAK_SRD_PHI, 'threshold': 0.5}, row) is False
    assert mr.condition_holds(WEAK_SRD_PHI, _row(srd_phi_left=0.5, srd_phi_right=None)) is False
    assert mr.condition_holds(WEAK_SRD_PHI, _row(srd_phi_left=float('nan'), srd_phi_right=0.1)) is False


def test_sign_mismatch_zero_same_sign_and_estimator_columns():
    c = {'kind': 'sign_mismatch', 'metric': 'delta_log2fc', 'estimator': 'median'}
    assert mr.condition_holds(c, _row(delta_median_left=0.2, delta_median_right=-0.1)) is True
    assert mr.condition_holds(c, _row(delta_median_left=0.0, delta_median_right=-0.1)) is False
    assert mr.condition_holds(c, _row(delta_median_left=0.2, delta_median_right=0.1)) is False
    assert mr.condition_holds(c, _row(delta_mean_left=0.2, delta_mean_right=-0.1)) is False


def test_ambiguous_condition_uses_selected_srd_variant_and_ties_are_not_ambiguous():
    c = {'kind': 'ambiguous', 'srdMetric': 'srd_phi', 'estimator': 'mean'}
    row = _row(srd_phi_left=1.0, srd_phi_right=3.0, srd_left=3.0, srd_right=1.0,
               delta_mean_left=3.0, delta_mean_right=1.0)
    assert mr.condition_holds(c, row) is True
    assert mr.condition_holds({**c, 'srdMetric': 'srd'}, row) is False
    assert mr.condition_holds(c, _row(srd_phi_left=2.0, srd_phi_right=2.0,
                                       delta_mean_left=3.0, delta_mean_right=1.0)) is False


def test_ambiguous_condition_on_real_scores():
    seg_c = np.array([3.0 + (0.1 if i % 2 == 0 else -0.1) for i in range(30)])
    profile = np.concatenate([np.full(2, 1.0), np.full(2, 2.0), seg_c])
    var_start = (np.arange(34) * 100).astype(float)
    _, rows = hm.score_rows(profile, hm.flat_partition([2, 4], 34), var_start, var_start + 100, 0)
    assert mr.condition_holds({'kind': 'ambiguous', 'srdMetric': 'srd_phi', 'estimator': 'mean'}, rows[0]) is True


def test_max_size_is_inclusive():
    c = {'kind': 'max_size', 'maxBins': 10}
    assert mr.condition_holds(c, _row(start=10, end=20)) is True
    assert mr.condition_holds(c, _row(start=10, end=21)) is False


def test_evaluate_rules_and_within_or_across():
    rows = [_row(start=10, end=20, srd_phi_left=0.5, srd_phi_right=-3.0),
            _row(start=20, end=25, srd_phi_left=2.0, srd_phi_right=-1.5)]
    both = _rule('and', [WEAK_SRD_PHI, {'kind': 'max_size', 'maxBins': 5}])
    assert [e['matched'] for e in mr.evaluate_rules(rows, [both])] == [[], []]
    size = _rule('size', [{'kind': 'max_size', 'maxBins': 5}])
    weak = _rule('weak', [WEAK_SRD_PHI])
    assert [e['matched'] for e in mr.evaluate_rules(rows, [size, weak])] == [['weak'], ['size']]


def test_claiming_rule_is_topmost_and_its_merge_into_sets_direction():
    row = _row(start=20, end=25, srd_phi_left=0.5, srd_phi_right=-3.0,
               delta_mean_left=-2.0, delta_mean_right=0.1)
    first = _rule('first', [{'kind': 'max_size', 'maxBins': 5}], merge_metric='delta_log2fc')
    second = _rule('second', [WEAK_SRD_PHI])
    ev = mr.evaluate_rules([row], [first, second])[0]
    assert ev['matched'] == ['first', 'second'] and ev['claimed_by'] == 'first'
    assert ev['side'] == 'right' and ev['score'] == pytest.approx(0.1) and ev['blocked'] is None
    ev = mr.evaluate_rules([row], [second, first])[0]
    assert ev['claimed_by'] == 'second' and ev['side'] == 'left' and ev['score'] == pytest.approx(0.5)


def test_merge_direction_tie_goes_left():
    ev = mr.evaluate_rules([_row(srd_phi_left=-2.0, srd_phi_right=2.0)],
                           [_rule('r', [{'kind': 'max_size', 'maxBins': 100}])])[0]
    assert ev['side'] == 'left' and ev['score'] == pytest.approx(2.0) and ev['blocked'] is None


def test_nonfinite_merge_into_flags_but_blocks():
    ev = mr.evaluate_rules([_row(srd_phi_left=None, srd_phi_right=1.0)],
                           [_rule('r', [{'kind': 'max_size', 'maxBins': 100}])])[0]
    assert ev['matched'] == ['r'] and ev['claimed_by'] == 'r'
    assert ev['blocked'] == 'nonfinite' and ev['side'] is None and ev['score'] is None


def test_gap_blocks_without_redirect_and_can_be_allowed():
    row = _row(srd_phi_left=0.2, srd_phi_right=5.0, gap_left=True)
    rules = [_rule('r', [WEAK_SRD_PHI])]
    ev = mr.evaluate_rules([row], rules)[0]
    assert ev['side'] == 'left' and ev['blocked'] == 'gap'
    ev = mr.evaluate_rules([row], rules, allow_gap_crossing=True)[0]
    assert ev['side'] == 'left' and ev['blocked'] is None


@pytest.mark.parametrize('rules', [
    [],
    [_rule('off', [WEAK_SRD_PHI], enabled=False)],
    [_rule('empty', [])],
    [_rule('k', [{'kind': 'bogus'}])],
    [_rule('m', [{'kind': 'weak_flank', 'metric': 'bogus', 'estimator': 'mean', 'threshold': 1.0}])],
    [_rule('e', [{'kind': 'sign_mismatch', 'metric': 'delta_log2fc', 'estimator': 'bogus'}])],
    [_rule('s', [{'kind': 'ambiguous', 'srdMetric': 'delta_log2fc', 'estimator': 'mean'}])],
    [_rule('t0', [{**WEAK_SRD_PHI, 'threshold': 0}])],
    [_rule('tn', [{**WEAK_SRD_PHI, 'threshold': float('nan')}])],
    [_rule('ti', [{**WEAK_SRD_PHI, 'threshold': float('inf')}])],
    [_rule('tb', [{**WEAK_SRD_PHI, 'threshold': True}])],
    [_rule('b0', [{'kind': 'max_size', 'maxBins': 0}])],
    [_rule('bf', [{'kind': 'max_size', 'maxBins': 2.5}])],
    [_rule('mi', [WEAK_SRD_PHI], merge_metric='bogus')],
    [_rule('p', [WEAK_SRD_PHI], pattern='#')],
    [GOOD, _rule('ok', [WEAK_SRD_PHI])],
    [{'name': 'no id', 'enabled': True, 'conditions': [WEAK_SRD_PHI]}],
    'not a list',
])
def test_validate_rules_rejects(rules):
    with pytest.raises(ValueError):
        mr.validate_rules(rules)


def test_validate_rules_ignores_disabled_and_returns_enabled_in_order():
    broken_but_off = _rule('off', [], enabled=False)
    second = _rule('two', [{'kind': 'max_size', 'maxBins': 3}])
    assert [r['id'] for r in mr.validate_rules([GOOD, broken_but_off, second])] == ['ok', 'two']


def _equivalence_rows():
    rows = []
    for profile, bounds in [
        (np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)]), [10, 20]),
        (np.concatenate([np.full(10, 8.0), np.full(10, 4.0), np.full(10, 2.0)]), [10, 20]),
        (np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(5, 8.0), np.full(15, 9.0)])
         + np.where(np.arange(40) % 2 == 0, 0.05, -0.05), [10, 20, 25]),
    ]:
        n = len(profile)
        vs = (np.arange(n) * 100).astype(float)
        rows += hm.score_rows(profile, hm.flat_partition(bounds, n), vs, vs + 100, 0)[1]
    return rows


@pytest.mark.parametrize('metric,estimator', [('srd', 'mean'), ('srd_phi', 'mean'), ('delta_log2fc', 'median')])
def test_presets_reproduce_old_layer_classifications(metric, estimator):
    for row in _equivalence_rows():
        l = mr.row_score(row, 'left', metric, estimator)
        r = mr.row_score(row, 'right', metric, estimator)
        assert mr.condition_holds({'kind': 'sign_mismatch', 'metric': metric, 'estimator': estimator}, row) \
            == hm.classify_transition(l, r)
        assert mr.condition_holds({'kind': 'weak_flank', 'metric': metric, 'estimator': estimator, 'threshold': 0.5}, row) \
            == (hm.classify_merge_proposal(l, r, 0.5)['weaker_side'] is not None)
        srd_l, srd_r = mr.row_score(row, 'left', 'srd_phi'), mr.row_score(row, 'right', 'srd_phi')
        fc_l = mr.row_score(row, 'left', 'delta_log2fc', estimator)
        fc_r = mr.row_score(row, 'right', 'delta_log2fc', estimator)
        assert mr.condition_holds({'kind': 'ambiguous', 'srdMetric': 'srd_phi', 'estimator': estimator}, row) \
            == (hm.classify_ambiguous(srd_l, srd_r, fc_l, fc_r) == 'ambiguous')


def _staircase(gap_at=None):
    levels = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(5, 8.0), np.full(15, 9.0)])
    profile = levels + np.where(np.arange(40) % 2 == 0, 0.05, -0.05)
    var_start = (np.arange(40) * 100).astype(float)
    var_end = var_start + 100
    if gap_at is not None:
        var_end[gap_at - 1] += 500
    return profile, [10, 20, 25], 40, var_start, var_end


def _delta_rule(rid, conditions):
    return _rule(rid, conditions, merge_metric='delta_log2fc', merge_estimator='mean')


WEAK_DELTA = {'kind': 'weak_flank', 'metric': 'delta_log2fc', 'estimator': 'mean', 'threshold': 0.5}


def test_run_rule_merge_within_a_rule_smallest_score_first():
    profile, bounds, n, vs, ve = _staircase()
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [_delta_rule('r', [WEAK_DELTA])])
    assert [s['removed_boundary'] for s in out['steps'][1:]] == [10, 25]
    assert out['steps'][1]['segment'] == [10, 20] and out['steps'][1]['side'] == 'left'
    assert out['steps'][1]['score'] == pytest.approx(np.log2(4.2 / 4.0))
    assert out['steps'][1]['rule_id'] == 'r' and out['steps'][1]['rule_name'] == 'r'
    assert out['final'] == [{'start': 0, 'end': 20}, {'start': 20, 'end': 40}]


def test_run_rule_merge_rule_order_beats_score():
    profile, bounds, n, vs, ve = _staircase()
    first = _delta_rule('small-and-weak', [{'kind': 'max_size', 'maxBins': 5}, WEAK_DELTA])
    second = _delta_rule('weak', [WEAK_DELTA])
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [first, second])
    s1, s2 = out['steps'][1], out['steps'][2]
    assert (s1['removed_boundary'], s1['rule_id'], s1['side']) == (25, 'small-and-weak', 'right')
    assert (s2['removed_boundary'], s2['rule_id'], s2['side']) == (10, 'weak', 'left')


def test_run_rule_merge_equal_scores_merge_leftmost_first():
    profile = np.concatenate([np.full(15, 4.0), np.full(5, 20.0)])
    vs = (np.arange(20) * 100).astype(float)
    out = mr.run_rule_merge(profile, [5, 10, 15], 20, vs, vs + 100, 0,
                            [_delta_rule('r', [{**WEAK_DELTA, 'threshold': 1.0}])])
    assert out['steps'][1]['removed_boundary'] == 5


def test_run_rule_merge_recomputes_phi_every_step():
    profile, bounds, n, vs, ve = _staircase()
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [_delta_rule('r', [WEAK_DELTA])])
    for k in (1, 2):
        before = [(s['start'], s['end']) for s in out['steps'][k - 1]['segments']]
        assert out['steps'][k]['phi'] == pytest.approx(hm.pooled_phi(profile, before))
    assert out['steps'][1]['phi'] != pytest.approx(out['steps'][2]['phi'])


def test_run_rule_merge_gap_blocks_without_redirect():
    profile, bounds, n, vs, ve = _staircase(gap_at=10)
    rules = [_delta_rule('r', [{**WEAK_DELTA, 'threshold': 0.1}])]
    blocked = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, rules)
    assert len(blocked['steps']) == 1 and blocked['final'] == blocked['original']
    allowed = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, rules, allow_gap_crossing=True)
    assert [s['removed_boundary'] for s in allowed['steps'][1:]] == [10]


def test_run_rule_merge_no_internal_segments():
    profile, _, n, vs, ve = _staircase()
    out = mr.run_rule_merge(profile, [20], n, vs, ve, 0, [_delta_rule('r', [WEAK_DELTA])])
    assert len(out['steps']) == 1 and out['final'] == out['original']


def test_run_rule_merge_terminates_when_every_flag_is_nonfinite():
    vs = (np.arange(40) * 100).astype(float)
    out = mr.run_rule_merge(np.zeros(40), [10, 20, 30], 40, vs, vs + 100, 0,
                            [_rule('r', [{'kind': 'max_size', 'maxBins': 100}])])
    assert len(out['steps']) == 1


def test_run_rule_merge_echoes_settings_and_step_zero_shape():
    profile, bounds, n, vs, ve = _staircase()
    rules = [_delta_rule('r', [WEAK_DELTA]), _rule('off', [], enabled=False)]
    out = mr.run_rule_merge(profile, bounds, n, vs, ve, 0, rules)
    assert out['settings'] == {'rules': rules, 'allow_gap_crossing': False}
    assert out['steps'][0] == {'step': 0, 'segments': out['original'], 'removed_boundary': None,
                               'merged_interval': None, 'segment': None, 'rule_id': None,
                               'rule_name': None, 'side': None, 'score': None, 'phi': None}
    with pytest.raises(ValueError):
        mr.run_rule_merge(profile, bounds, n, vs, ve, 0, [])


from fastapi import FastAPI
from fastapi.testclient import TestClient


class _FakeStore:
    def __init__(self):
        import pandas as pd
        self.identity = {'hash': 'fake', 'files': []}
        self.var_seq = np.array(['chr1'] * 30)
        self.var_start = (np.arange(30) * 100).astype(float)
        self.var_end = self.var_start + 100
        self.members = {'cluster0': np.arange(5), 'all': np.arange(5)}
        pb = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)])
        self.arrays = {'pb_x': np.tile(pb, (2, 1))}
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


def test_hatch_merge_route_runs_rules_and_reports_rule_per_step():
    resp = _client().post('/api/integration/hatch-merge', json={
        'chrom': 'chr1', 'source': 'cluster0', 'rules': [_delta_rule('weak', [WEAK_DELTA])]})
    assert resp.status_code == 200
    body = resp.json()
    assert body['final'] == [{'start': 0, 'end': 20}, {'start': 20, 'end': 30}]
    assert body['steps'][1]['rule_id'] == 'weak' and body['steps'][1]['side'] == 'left'
    assert body['settings']['allow_gap_crossing'] is False
    assert body['provenance']['dataset'] == {'hash': 'fake', 'files': []}


@pytest.mark.parametrize('rules', [
    [],
    [_rule('off', [WEAK_DELTA], enabled=False)],
    [_rule('empty', [])],
    [_rule('t', [{**WEAK_DELTA, 'threshold': -1.0}])],
    [_rule('k', [{'kind': 'bogus'}])],
    [_rule('p', [WEAK_DELTA], pattern='#')],
    [_rule('dup', [WEAK_DELTA]), _rule('dup', [WEAK_DELTA])],
])
def test_hatch_merge_route_rejects_bad_rules(rules):
    resp = _client().post('/api/integration/hatch-merge', json={'chrom': 'chr1', 'source': 'cluster0', 'rules': rules})
    assert resp.status_code == 422


def test_hatch_merge_route_requires_rules_field():
    resp = _client().post('/api/integration/hatch-merge', json={'chrom': 'chr1', 'source': 'cluster0'})
    assert resp.status_code == 422
