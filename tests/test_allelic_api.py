import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from fastapi import FastAPI
from fastapi.testclient import TestClient

from server import allelic_api as aa
from server import allelic_data as ald
from server import allelic_scan as sc
from allelic_fixtures import make_allelic

FAST = {'min_at': 5, 'k': 2, 'n_permutations': 100, 'alpha': 0.05, 'min_child_bins': 10, 'min_terminal_bins': 10, 'global_pool': 'genome'}


def _client(tmp_path, f=None):
    f = f or make_allelic(tmp_path)
    d = ald.AllelicData(f['dir'], str(tmp_path / 'cache'), f['var_seq'], f['var_start'], f['var_end'], f['cell_ids'], f['members'], f['sources'])
    ex = ThreadPoolExecutor(2)
    app = FastAPI()
    mgr = aa.register_allelic(app, lambda: d, lambda: ex, None, str(tmp_path / 'cache'))
    c = TestClient(app)
    c.get('/api/allelic/status'); d.wait(30)
    return c, d, mgr, f


def _wait(c, jid, timeout=120):
    t0 = time.time()
    while time.time() - t0 < timeout:
        v = c.get(f'/api/allelic/segment/{jid}').json()
        if v['status'] in ('done', 'error'):
            return v
        time.sleep(0.05)
    raise AssertionError('job did not finish')


def test_normalize_defaults_and_errors():
    s = aa.normalize('abs_dbaf', {})
    assert s == {'min_at': 10, 'k': 6, 'n_permutations': 1000, 'alpha': 0.001, 'min_child_bins': 20, 'min_terminal_bins': 50, 'global_pool': 'genome'}
    for bad, field in [({'alpha': 0}, 'alpha'), ({'alpha': 0.2}, 'alpha'), ({'k': 2.5}, 'k'), ({'min_at': 0}, 'min_at'),
                       ({'n_permutations': 'x'}, 'n_permutations'), ({'zzz': 1}, 'zzz'), ({'global_pool': 'cell'}, 'global_pool')]:
        try:
            aa.normalize('abs_dbaf', bad)
        except ValueError as e:
            assert str(e).startswith(f'{field}:')
        else:
            raise AssertionError(bad)
    try:
        aa.normalize('baf', {})
    except ValueError as e:
        assert str(e).startswith('signal:')


def test_status_and_tracks(tmp_path):
    c, d, _, f = _client(tmp_path)
    s = c.get('/api/allelic/status').json()
    assert s['state'] == 'ready' and s['available'] and s['sources'] == f['sources']
    t = c.get('/api/allelic/tracks', params={'chrom': 'chr1', 'sources': 'cluster1,all', 'signal': 'dbaf', 'min_at': 5}).json()
    assert t['n'] == 120 and [x['source'] for x in t['tracks']] == ['cluster1', 'all']
    v = t['tracks'][0]['values']
    assert len(v) == 120 and all(x is None or -0.5 <= x <= 0.5 for x in v)
    assert c.get('/api/allelic/tracks', params={'chrom': 'chrZ', 'sources': 'all'}).status_code == 404
    assert c.get('/api/allelic/tracks', params={'chrom': 'chr1', 'sources': 'cluster9'}).status_code == 404
    assert c.get('/api/allelic/tracks', params={'chrom': 'chr1', 'sources': 'all', 'signal': 'baf'}).status_code == 422


def test_segment_runs_and_reuses(tmp_path):
    c, d, _, f = _client(tmp_path)
    v = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': FAST}).json()
    done = _wait(c, v['id'])
    assert done['status'] == 'done' and done['settings'] == FAST and done['signal'] == 'abs_dbaf'
    bps = done['result']['breakpoints']
    c1 = [b for b in bps if b['source'] == 'cluster1' and b['chromosome'] == 'chr1']
    assert any(abs(b['absolute_bin'] - 60) <= 2 for b in c1)
    b = c1[0]
    assert b['bp_genomic_start'] == int(f['var_start'][b['absolute_bin']])
    segs = [s for s in done['result']['segments'] if s['source'] == 'cluster0' and s['chromosome'] == 'chr2']
    assert segs[0]['bin_start'] == 0 and segs[-1]['bin_end'] == 79
    assert {r['source'] for r in done['result']['audit']} <= set(f['sources'])
    again = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': FAST}).json()
    assert again['id'] == v['id'] and again['status'] == 'done'
    other = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': {**FAST, 'alpha': 0.01}}).json()
    assert other['id'] != v['id']
    _wait(c, other['id'])


def test_global_pool_reaches_the_scan(tmp_path, monkeypatch):
    c, d, _, f = _client(tmp_path)
    pools = {}
    real = sc.scan_task

    def spy(task):
        pools[(task[-1] is None, task[0], task[1])] = None if task[-1] is None else len(task[-1])
        return real(task)
    monkeypatch.setattr(sc, 'scan_task', spy)
    _wait(c, c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': FAST}).json()['id'])
    kept_total = sum(int(d.track('cluster0', ch, 'abs_dbaf', 5)['kept'].sum()) for ch in d.chromosomes)
    assert pools[(False, 'cluster0', 'chr1')] == kept_total == pools[(False, 'cluster0', 'chr2')]
    _wait(c, c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': {**FAST, 'global_pool': 'chromosome'}}).json()['id'])
    assert (True, 'cluster0', 'chr1') in pools


def test_bad_requests(tmp_path):
    c, *_ = _client(tmp_path)
    r = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': {'alpha': 0}})
    assert r.status_code == 422 and r.json()['detail'].startswith('alpha:')
    assert c.post('/api/allelic/segment', json={'signal': 'nope', 'settings': {}}).status_code == 422
    assert c.get('/api/allelic/segment/abc').status_code == 404


def test_cache_survives_a_new_manager(tmp_path):
    f = make_allelic(tmp_path)
    c, *_ = _client(tmp_path, f)
    jid = c.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()['id']
    first = _wait(c, jid)
    c2, *_ = _client(tmp_path, f)
    v = c2.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()
    assert v['id'] == jid and v['status'] == 'done' and v['result'] == first['result']


def test_queue_replaces_the_waiting_request(tmp_path, monkeypatch):
    c, *_ = _client(tmp_path)
    gate = threading.Event()
    real = sc.scan_task

    def slow(task):
        gate.wait(30)
        return real(task)
    monkeypatch.setattr(sc, 'scan_task', slow)
    a = c.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()
    b = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': FAST}).json()
    cc = c.post('/api/allelic/segment', json={'signal': 'shrunk_abs_dbaf', 'settings': FAST}).json()
    assert a['status'] == 'running' and b['status'] == 'queued' and cc['status'] == 'queued'
    vb = c.get(f"/api/allelic/segment/{b['id']}").json()
    assert vb['status'] == 'error' and vb['error'] == 'replaced by a newer request'
    gate.set()
    assert _wait(c, a['id'])['status'] == 'done'
    assert _wait(c, cc['id'])['status'] == 'done'


def test_failed_task_fails_the_job_and_is_not_cached(tmp_path, monkeypatch):
    c, *_ = _client(tmp_path)

    def boom(task):
        raise RuntimeError('kaput')
    monkeypatch.setattr(sc, 'scan_task', boom)
    v = _wait(c, c.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()['id'])
    assert v['status'] == 'error' and 'kaput' in v['error'] and v['result'] is None
    monkeypatch.undo()
    again = c.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()
    assert again['status'] in ('running', 'queued')
    assert _wait(c, again['id'])['status'] == 'done'
