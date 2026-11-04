import numpy as np
import pandas as pd
import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from server import cn_recall as cr
from server import merge_rounds as mrd
from test_cn_recall import LOOSE, NOOP, synthetic_inputs

MID = [{**NOOP[0], 'combine': 'all', 'tau': 0.3, 'srdRef': 1000.0}]
R = lambda rules, gap=False: {'rules': rules, 'allow_gap_crossing': gap}


def _rows(df):
    return sorted(map(tuple, df[['source', 'chromosome', 'absolute_bin']].astype(str).values.tolist()))


def test_fold_empty_is_original():
    inp = synthetic_inputs(percell=False)
    bp, per = mrd.fold(inp, [])
    assert _rows(bp) == _rows(inp.breakpoints) and per == []


def test_fold_one_round_equals_merge_genome():
    inp = synthetic_inputs(percell=False)
    bp, per = mrd.fold(inp, [R(LOOSE)])
    m = cr.merge_genome(inp, LOOSE)
    assert _rows(bp) == _rows(m['merged_breakpoints'])
    assert per[0]['merges'] == sum(p['merges'] for p in m['pairs']) > 0


def test_same_rules_twice_is_stable():
    inp = synthetic_inputs(percell=False)
    one, _ = mrd.fold(inp, [R(LOOSE)])
    two, per = mrd.fold(inp, [R(LOOSE), R(LOOSE)])
    assert _rows(one) == _rows(two) and per[1]['merges'] == 0


def test_round_two_is_subset_with_consistent_counts():
    inp = synthetic_inputs(percell=False)
    one, p1 = mrd.fold(inp, [R(MID)])
    two, p2 = mrd.fold(inp, [R(MID), R(LOOSE)])
    assert set(_rows(two)) <= set(_rows(one)) <= set(_rows(inp.breakpoints))
    first = {(p['source'], p['chrom']): p for p in p1[0]['pairs']}
    for p in p2[1]['pairs']:
        assert p['original'] == first[(p['source'], p['chrom'])]['merged']
        assert p['merged'] == p['original'] - p['merges']


def test_normalize_rounds_names_the_round():
    assert mrd.normalize_rounds(None) == [] and mrd.normalize_rounds([]) == []
    assert mrd.normalize_rounds([{'rules': LOOSE}]) == [R(LOOSE)]
    with pytest.raises(ValueError, match=r'^round 2: '):
        mrd.normalize_rounds([R(LOOSE), R([])])
    with pytest.raises(ValueError, match=r'^round 1: rules must be a list'):
        mrd.normalize_rounds([{'rules': 'x'}])


def test_rounds_hash_is_canonical():
    assert mrd.rounds_hash([R(LOOSE)]) == mrd.rounds_hash([{'allow_gap_crossing': False, 'rules': LOOSE}])
    assert mrd.rounds_hash([R(LOOSE)]) != mrd.rounds_hash([R(LOOSE, True)])
    assert len(mrd.rounds_hash([])) == 64


def test_cache_reuses_entry_and_keys_on_breakpoints_object():
    c = mrd.BaseCache(2)
    inp = synthetic_inputs(percell=False)
    a = c.get(inp, [R(LOOSE)])
    assert c.get(inp, [R(LOOSE)]) is a
    other = synthetic_inputs(percell=False)
    assert c.get(other, [R(LOOSE)]) is not a


def test_boundaries_from_bp():
    bp = pd.DataFrame({'source': ['a', 'a', 'b'], 'chromosome': ['chr1', 'chr1', 'chr2'], 'absolute_bin': [9, 3, 4]})
    assert cr.boundaries_from_bp(bp, ['a', 'b'], ['chr1', 'chr2']) == {'a': {'chr1': [3, 9], 'chr2': []}, 'b': {'chr1': [], 'chr2': [4]}}


class _FakeStore:
    def __init__(self):
        self.identity = {'hash': 'fake', 'files': []}
        self.var_seq = np.array(['chr1'] * 30)
        self.var_start = (np.arange(30) * 100).astype(float)
        self.var_end = self.var_start + 100
        self.members = {'cluster0': np.arange(5), 'all': np.arange(5)}
        pb = np.concatenate([np.full(10, 4.0), np.full(10, 4.2), np.full(10, 20.0)])
        self.arrays = {'pb_x': np.tile(pb, (2, 1))}
        self.tables = {'breakpoints': pd.DataFrame({'source': ['cluster0', 'cluster0'], 'chromosome': ['chr1', 'chr1'],
                                                    'absolute_bin': [10, 20]})}
        self.available = True


def _client():
    from server.hatch_merge_api import register_hatch_merge
    app = FastAPI()
    store = _FakeStore()
    register_hatch_merge(app, lambda: store)
    return TestClient(app)


def test_hatch_scores_without_rounds_equals_empty_rounds():
    c = _client()
    a = c.post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'cluster0'}).json()
    b = c.post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'cluster0', 'rounds': []}).json()
    assert a == b and len(a['rows']) == 1


def test_hatch_routes_read_the_base():
    c = _client()
    rounds = [R(LOOSE)]
    s = c.post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'cluster0', 'rounds': rounds}).json()
    assert len(s['rows']) < 1 and len(s['segments']) == 2
    m = c.post('/api/integration/hatch-merge', json={'chrom': 'chr1', 'source': 'cluster0', 'rules': LOOSE, 'rounds': rounds}).json()
    assert len(m['original']) == len(s['segments'])


def test_hatch_routes_reject_invalid_round():
    c = _client()
    r = c.post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'cluster0', 'rounds': [R([])]})
    assert r.status_code == 422 and r.json()['detail'].startswith('round 1: ')
    r = c.post('/api/integration/hatch-merge', json={'chrom': 'chr1', 'source': 'cluster0', 'rules': LOOSE, 'rounds': [R([])]})
    assert r.status_code == 422 and r.json()['detail'].startswith('round 1: ')


def test_hatch_merge_records_rounds_in_settings_only_when_present():
    c = _client()
    body = {'chrom': 'chr1', 'source': 'cluster0', 'rules': LOOSE}
    assert 'rounds' not in c.post('/api/integration/hatch-merge', json=body).json()['settings']
    assert 'rounds' not in c.post('/api/integration/hatch-merge', json={**body, 'rounds': []}).json()['settings']
    m = c.post('/api/integration/hatch-merge', json={**body, 'rounds': [R(LOOSE)]}).json()
    assert m['settings']['rounds'] == [R(LOOSE)]
