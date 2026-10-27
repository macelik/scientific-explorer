import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from fastapi import FastAPI
from fastapi.testclient import TestClient

from server import call_noise as cn
from test_cn_recall import LOOSE, NOOP, synthetic_inputs

LOOSE2 = [{**LOOSE[0], 'tau': 4.0}]
R = lambda rules: {'rules': rules, 'allow_gap_crossing': False}


def _api(inp, keep=3):
    from server.cn_recall_api import register_cn_recall
    app = FastAPI()
    ex = ThreadPoolExecutor(2)
    register_cn_recall(app, lambda: (inp, {'hash': 'fake'}, ex), keep=keep)
    return TestClient(app)


def _wait(c, jid, timeout=120):
    t0 = time.time()
    while time.time() - t0 < timeout:
        v = c.get(f'/api/integration/recall/{jid}').json()
        if v['status'] != 'running':
            return v
        time.sleep(0.05)
    raise AssertionError('job did not finish')


def _recall(c, rules, rounds=()):
    return _wait(c, c.post('/api/integration/recall', json={'rules': rules, 'rounds': list(rounds)}).json()['id'])


def test_base_original():
    inp = synthetic_inputs()
    c = _api(inp)
    b = c.post('/api/integration/base', json={'rounds': []}).json()
    assert b['n_rounds'] == 0 and len(b['breakpoints']) == len(inp.breakpoints)
    assert b['percell'] == {'status': 'original', 'job': None}
    assert b['noise']['cluster0']['genome']['cluster'] == cn.cluster_noise(inp.before_cluster_segments)['cluster0']['genome']
    assert b['noise']['cluster0']['genome']['cells']['n_cells'] == 6 and b['noise_original'] == b['noise']


def test_adopt_then_base_uses_the_pinned_recall():
    inp = synthetic_inputs()
    c = _api(inp)
    v = _recall(c, LOOSE)
    assert c.post(f"/api/integration/recall/{v['id']}/pin", json={'rounds': [R(LOOSE)]}).status_code == 200
    b = c.post('/api/integration/base', json={'rounds': [R(LOOSE)]}).json()
    assert b['percell'] == {'status': 'ready', 'job': v['id']} and b['n_rounds'] == 1
    assert len(b['breakpoints']) == len(v['merge']['merged_breakpoints'])
    assert b['cluster_segments'] == v['cluster_segments']
    assert b['noise']['cluster0']['genome'] == v['diff']['noise']['cluster0']['genome']
    assert b['noise_original']['cluster0']['genome']['cluster']['segments'] > b['noise']['cluster0']['genome']['cluster']['segments']


def test_pin_refuses_wrong_recipe_unknown_and_invalid():
    inp = synthetic_inputs()
    c = _api(inp)
    v = _recall(c, LOOSE)
    assert c.post(f"/api/integration/recall/{v['id']}/pin", json={'rounds': [R(LOOSE2)]}).status_code == 409
    assert c.post('/api/integration/recall/nope/pin', json={'rounds': [R(LOOSE)]}).status_code == 404
    assert c.post(f"/api/integration/recall/{v['id']}/pin", json={'rounds': [R([])]}).status_code == 422


def test_restart_path_recomputes_base_percell():
    inp = synthetic_inputs()
    ref = _api(inp)
    v = _recall(ref, LOOSE)
    want = ref.get(f"/api/integration/recall/{v['id']}/percell.bin").content
    c = _api(inp)   # a fresh manager: nothing pinned, as after a restart
    b = c.post('/api/integration/base', json={'rounds': [R(LOOSE)]}).json()
    assert b['percell']['status'] == 'running' and b['percell']['job']
    done = _wait(c, b['percell']['job'])
    assert done['status'] == 'done' and done['percell_ready']
    b2 = c.post('/api/integration/base', json={'rounds': [R(LOOSE)]}).json()
    assert b2['percell'] == {'status': 'ready', 'job': b['percell']['job']}
    assert b2['noise']['cluster0']['genome']['cells'] is not None
    assert c.get(f"/api/integration/recall/{b['percell']['job']}/percell.bin").content == want


def test_recall_on_a_base_reports_against_the_base():
    inp = synthetic_inputs()
    c = _api(inp)
    first = _recall(c, LOOSE)
    c.post(f"/api/integration/recall/{first['id']}/pin", json={'rounds': [R(LOOSE)]})
    v = _recall(c, LOOSE, [R(LOOSE)])
    assert v['status'] == 'done' and v['id'] != first['id']
    assert v['settings'] == {'rules': LOOSE, 'allow_gap_crossing': False, 'rounds': [R(LOOSE)]}
    assert v['provenance']['scope'] == ('Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on a '
                                        'rule-merged segmentation (1 earlier rounds); originals unchanged.')
    s0 = {d['source']: d for d in first['diff']['sources']}
    for d in v['diff']['sources']:
        assert d['merges'] == 0 and d['bins_changed'] == 0 and d['segments_before'] == s0[d['source']]['segments_after']
    assert v['diff']['cells']['n_changed_cells'] == 0


def test_recall_on_unpinned_base_waits_for_base_calls():
    inp = synthetic_inputs()
    c = _api(inp)
    v = _recall(c, LOOSE, [R(LOOSE)])
    assert v['status'] == 'done' and v['diff']['cells']['n_changed_cells'] == 0


def test_rounds_change_the_job():
    inp = synthetic_inputs(percell=False)
    c = _api(inp)
    a = c.post('/api/integration/recall', json={'rules': LOOSE}).json()
    b = c.post('/api/integration/recall', json={'rules': LOOSE, 'rounds': []}).json()
    d = c.post('/api/integration/recall', json={'rules': LOOSE, 'rounds': [R(NOOP)]}).json()
    assert a['id'] == b['id'] and d['id'] != a['id']
    _wait(c, a['id']); _wait(c, d['id'])


def test_invalid_round_is_422_everywhere():
    inp = synthetic_inputs(percell=False)
    c = _api(inp)
    for url, body in (('/api/integration/base', {'rounds': [R([])]}),
                      ('/api/integration/recall', {'rules': LOOSE, 'rounds': [R([])]})):
        r = c.post(url, json=body)
        assert r.status_code == 422 and r.json()['detail'].startswith('round 1: ')


def test_pinned_job_survives_eviction_and_repin_releases_it():
    inp = synthetic_inputs()
    c = _api(inp, keep=1)
    a = _recall(c, LOOSE)
    c.post(f"/api/integration/recall/{a['id']}/pin", json={'rounds': [R(LOOSE)]})
    for tau in (1e-9, 2e-9, 3e-9):
        _recall(c, [{**NOOP[0], 'tau': tau}])
    assert c.get(f"/api/integration/recall/{a['id']}").status_code == 200
    b = _recall(c, LOOSE2, [R(LOOSE)])
    assert c.post(f"/api/integration/recall/{b['id']}/pin", json={'rounds': [R(LOOSE), R(LOOSE2)]}).status_code == 200
    _recall(c, [{**NOOP[0], 'tau': 4e-9}])
    assert c.get(f"/api/integration/recall/{a['id']}").status_code == 404
    assert c.get(f"/api/integration/recall/{b['id']}").status_code == 200


def test_wait_percell_survives_eviction_from_jobs():
    import threading
    from server.cn_recall_api import RecallManager
    m = RecallManager(keep=1)
    jid = m._new('h', {}, 's', {}, 'k', 'percell')
    job = m.jobs[jid]
    del m.jobs[jid]                      # evicted while a re-call still holds the job object
    u8 = np.zeros((1, 1), dtype=np.uint8)

    def finish():
        job.update(_u8=u8, _best=None, status='done')
        job['_event'].set()
    threading.Timer(0.05, finish).start()
    got, best = m.wait_percell(job)
    assert got is u8 and best is None
