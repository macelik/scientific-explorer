"""CN re-call pipeline (spec docs/superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md)."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import importlib.util

import numpy as np
import pandas as pd
import pytest

from server import cn_recall as cr
from server import merge_rules as mr
from server.config import Config

PKG = Config.from_env().package_dir
NOOP = [{"id": "rule0", "type": "flank_threshold", "name": "Rule 0", "criteria": ["delta_mean", "srd_phi"], "combine": "all",
         "tau": 1e-9, "srdRef": 1e-9, "maxBins": 100, "mergeSignMismatch": False, "mergeSignsOnly": False,
         "hatch": {"pattern": "x", "color": "#b91c1c"}}]
LOOSE = [{**NOOP[0], "combine": "any", "tau": 5.0, "srdRef": 1000.0, "mergeSignMismatch": True}]
STORY = Path(__file__).resolve().parents[2] / 'integration-story'


def synthetic_inputs(percell=True, cell_cluster=None):
    rng = np.random.default_rng(7)
    var_seq = np.array(['chr1'] * 40 + ['chr2'] * 30)
    local = np.concatenate([np.arange(40), np.arange(30)])
    var_start = (1 + local * 100000).astype(float)
    var_end = var_start + 99999
    cl = np.array([0] * 6 + [1] * 6) if cell_cluster is None else np.asarray(cell_cluster)
    mean = np.full((12, 70), 10.0)
    mean[6:, 10:20] = 15.0          # cluster1 gain on chr1
    mean[:6, 40:55] = 5.0           # cluster0 loss on chr2
    X = rng.poisson(mean).astype(float)
    sources = ['cluster0', 'cluster1', 'all']
    members = {'cluster0': np.arange(6), 'cluster1': np.arange(6, 12), 'all': np.arange(12)}
    bp = pd.DataFrame([(s, 'chr1', b) for s in sources for b in (10, 20, 23)] +
                      [(s, 'chr2', b) for s in ('cluster0', 'all') for b in (15, 17)],
                      columns=['source', 'chromosome', 'absolute_bin'])
    bp['recursion_depth'] = 1
    pb = np.stack([X[members[s]].mean(axis=0) for s in sources])
    inp = cr.RecallInputs(var_seq=var_seq, var_start=var_start, var_end=var_end, pb_x=pb.astype(np.float32),
                          breakpoints=bp, members=members, cell_ids=[f'c{i}' for i in range(12)], cell_cluster=cl,
                          read_rows=lambda lo, hi: X[lo:hi].copy(), before_cluster_segments=pd.DataFrame(),
                          before_percell_u8=None, before_percell_best_s=np.full(12, np.nan), package_dir=PKG,
                          sources=sources)
    # "before" = the same pipeline on the original boundaries
    orig = {s: {c: sorted(bp.loc[(bp.source == s) & (bp.chromosome == c), 'absolute_bin'].tolist()) for c in ('chr1', 'chr2')} for s in sources}
    segs = {s: cr.seg_labels(var_seq, orig[s]) for s in sources}
    rows = []
    for k, s in enumerate(sources):
        rows += cr.cluster_rows(s, len(members[s]), pb[k], segs[s], var_seq)[0]
    inp.before_cluster_segments = pd.DataFrame(rows)
    if percell:
        u8 = np.zeros((12, 70), np.uint8); bs = np.zeros(12)
        for i in range(12):
            u8[i:i + 1], bs[i:i + 1] = cr.call_cells(X[i:i + 1], segs.get(f'cluster{cl[i]}', segs['all']), PKG)
        inp.before_percell_u8, inp.before_percell_best_s = u8, bs
    return inp


def _run(inp, rules):
    log = []
    def update(**kw):
        log.append(kw)
    with ThreadPoolExecutor(2) as ex:
        out = cr.run_recall(inp, rules, False, ex, update, chunk=5, max_inflight=2)
    state = {}
    for kw in log:
        state.update(kw)
    return out, state, [kw['stage'] for kw in log if 'stage' in kw]


def test_seg_labels_matches_make_karyograms():
    script = STORY / 'scripts' / 'make_karyograms.py'
    if not script.exists():
        pytest.skip('make_karyograms.py not available')
    spec = importlib.util.spec_from_file_location('mk', script); mk = importlib.util.module_from_spec(spec); spec.loader.exec_module(mk)
    var_seq = np.array(['chr1'] * 40 + ['chr2'] * 30 + ['chr3'] * 5)
    raw = pd.DataFrame({'source': ['g'] * 3, 'chromosome': ['chr1', 'chr1', 'chr2'], 'absolute_bin': [10, 23, 15]})
    want = mk.build_seg_labels(raw, 'g', ['chr1', 'chr2', 'chr3'], var_seq, len(var_seq))
    got = cr.seg_labels(var_seq, {'chr1': [23, 10], 'chr2': [15]})
    assert np.array_equal(got, want)


def test_merge_genome_noop_keeps_every_boundary():
    inp = synthetic_inputs(percell=False)
    m = cr.merge_genome(inp, NOOP)
    assert m['merged_breakpoints'].equals(inp.breakpoints)
    assert all(p['merges'] == 0 and p['original'] == p['merged'] for p in m['pairs'])
    assert len(m['pairs']) == 3 * 2


def test_merge_genome_loose_merges_subset_and_counts():
    inp = synthetic_inputs(percell=False)
    m = cr.merge_genome(inp, LOOSE)
    assert sum(p['merges'] for p in m['pairs']) > 0
    for p in m['pairs']:
        assert p['merges'] == p['original'] - p['merged']
        orig = set(inp.breakpoints.loc[(inp.breakpoints.source == p['source']) & (inp.breakpoints.chromosome == p['chrom']), 'absolute_bin'])
        assert set(m['boundaries'][p['source']][p['chrom']]) <= orig
    c1chr2 = next(p for p in m['pairs'] if p['source'] == 'cluster1' and p['chrom'] == 'chr2')
    assert (c1chr2['original'], c1chr2['merges']) == (1, 0)
    idx = np.flatnonzero(inp.var_seq == 'chr1')
    direct = mr.run_rule_merge(inp.pb_x[0, idx], [10, 20, 23], len(idx), inp.var_start, inp.var_end, 0, LOOSE)
    assert m['boundaries']['cluster0']['chr1'] == [s['start'] for s in direct['final'][1:]]
    kept = m['merged_breakpoints']
    assert list(kept.columns) == list(inp.breakpoints.columns)


def test_merge_genome_names_failing_pair():
    inp = synthetic_inputs(percell=False)
    inp.pb_x = inp.pb_x[:1]              # source index 1 out of range
    with pytest.raises(cr.RecallError, match=r'^merge cluster1/chr1: '):
        cr.merge_genome(inp, NOOP)


def test_cluster_rows_columns_and_quirk():
    inp = synthetic_inputs(percell=False)
    df = inp.before_cluster_segments
    assert list(df.columns) == ['group', 'n_cells', 'best_s', 'segment', 'chromosome', 'bin_start', 'bin_end', 'n_bins',
                                'continuous_cn', 'integer_cn', 'watson', 'cn_call']
    assert (df['continuous_cn'] == df['cn_call'].round(2)).all()
    assert df.groupby('group')['n_bins'].sum().tolist() == [70, 70, 70]
    assert (df.groupby('group')['best_s'].nunique() == 1).all()


def test_cluster_u8_matches_encoded_combined():
    inp = synthetic_inputs(percell=False)
    seg = cr.seg_labels(inp.var_seq, {'chr1': [10, 20, 23], 'chr2': [15, 17]})
    rows, combined = cr.cluster_rows('cluster0', 6, inp.pb_x[0].astype(float), seg, inp.var_seq)
    assert np.array_equal(cr.cluster_u8_from_segments(pd.DataFrame(rows), 'cluster0', 70), cr.encode_states(combined))
    assert np.array_equal(cr.encode_states(np.array([0.5, np.nan, 2.0])), np.array([1, 255, 4], np.uint8))


def test_transitions_counts():
    acc = np.zeros(65536, np.int64)
    cr.add_transitions(acc, np.array([2, 2, 2, 3, 255], np.uint8), np.array([2, 3, 3, 3, 2], np.uint8))
    assert cr.transitions_list(acc) == [{'before': 1.0, 'after': 1.5, 'bins': 2}, {'before': None, 'after': 1.0, 'bins': 1}]


def test_run_recall_noop_reproduces_before():
    inp = synthetic_inputs()
    out, st, stages = _run(inp, NOOP)
    assert stages == ['merge', 'cluster', 'percell']
    assert all(d['bins_changed'] == 0 for d in st['diff']['sources'])
    assert st['diff']['transitions'] == {'clusters': [], 'cells': []}
    assert np.array_equal(out['percell_u8'], inp.before_percell_u8)
    assert st['diff']['cells']['n_changed_cells'] == 0 and st['percell_ready'] is True
    assert pd.DataFrame(st['cluster_segments']).equals(inp.before_cluster_segments)


def test_run_recall_loose_changes_segments_and_reports():
    inp = synthetic_inputs()
    out, st, _ = _run(inp, LOOSE)
    src = {d['source']: d for d in st['diff']['sources']}
    assert src['cluster0']['segments_after'] < src['cluster0']['segments_before']
    assert st['merge']['pairs'] and isinstance(st['merge']['merged_breakpoints'], list)
    assert out['percell_u8'].shape == (12, 70)
    pc = st['diff']['cells']['per_cell']
    assert [r['cell'] for r in pc] == inp.cell_ids
    assert all(r['bins_changed'] == int((inp.before_percell_u8[i] != out['percell_u8'][i]).sum()) for i, r in enumerate(pc))


def test_run_recall_unknown_cluster_fails_with_stage_message():
    inp = synthetic_inputs(cell_cluster=[0] * 6 + [1] * 5 + [7])
    with pytest.raises(cr.RecallError, match=r'^percell: cluster 7 has no source'):
        _run(inp, NOOP)


# ---- real data: re-calling on the unmerged segmentation reproduces the original files exactly
DATA = STORY / 'data'
_REAL = all((DATA / p).exists() for p in ('count_matrix.h5ad', 'raw_accepted_breakpoints.tsv',
                                          'derived/karyogram_cluster_segments.tsv', 'derived/percell_cn_states.npy',
                                          'derived/percell_best_s.tsv'))


def real_inputs():
    import h5py
    from server.integration import _csr_rows, _read_cat
    h5 = h5py.File(DATA / 'count_matrix.h5ad', 'r')
    var_seq = _read_cat(h5['var'], 'seq').astype(str)
    n_bins = len(var_seq)
    best = pd.read_csv(DATA / 'derived/percell_best_s.tsv', sep='\t')
    assert np.array_equal(h5['obs/cellID'][...].astype(str), best['cell'].to_numpy().astype(str))
    cl = best['wnn_cluster'].to_numpy().astype(int)
    members = {f'cluster{k}': np.flatnonzero(cl == k) for k in range(5)}
    members['all'] = np.arange(len(cl))
    inp = cr.RecallInputs(var_seq=var_seq, var_start=h5['var/start'][...], var_end=h5['var/end'][...],
                          pb_x=np.zeros((6, n_bins), np.float32),
                          breakpoints=pd.read_csv(DATA / 'raw_accepted_breakpoints.tsv', sep='\t'),
                          members=members, cell_ids=best['cell'].astype(str).tolist(), cell_cluster=cl,
                          read_rows=lambda lo, hi: _csr_rows(h5['X'], lo, hi, n_bins),
                          before_cluster_segments=pd.read_csv(DATA / 'derived/karyogram_cluster_segments.tsv', sep='\t'),
                          before_percell_u8=None, before_percell_best_s=best['best_s'].to_numpy(), package_dir=PKG)
    return inp


@pytest.mark.skipif(not _REAL, reason='integration-story data not available')
def test_real_noop_cluster_identity():
    inp = real_inputs()
    pb = cr.pseudobulk(inp)
    inp.pb_x = pb.astype(np.float32)
    m = cr.merge_genome(inp, NOOP)
    assert len(m['merged_breakpoints']) == len(inp.breakpoints)
    rows = []
    for k, s in enumerate(inp.sources):
        rows += cr.cluster_rows(s, len(inp.members[s]), pb[k], cr.seg_labels(inp.var_seq, m['boundaries'][s]), inp.var_seq)[0]
    got = pd.DataFrame(rows)
    want = inp.before_cluster_segments
    assert len(got) == len(want)
    for col in want.columns:
        assert got[col].tolist() == want[col].tolist(), col


@pytest.mark.skipif(not _REAL, reason='integration-story data not available')
def test_real_noop_percell_identity_sample():
    inp = real_inputs()
    m = cr.merge_genome(inp, NOOP)
    cn = np.load(DATA / 'derived/percell_cn_states.npy', mmap_mode='r')
    for i in np.random.default_rng(0).choice(len(inp.cell_ids), 20, replace=False):
        seg = cr.seg_labels(inp.var_seq, m['boundaries'][f'cluster{inp.cell_cluster[i]}'])
        u8, bs = cr.call_cells(inp.read_rows(int(i), int(i) + 1), seg, PKG)
        assert np.array_equal(u8[0], cr.encode_states(np.asarray(cn[i]))), int(i)
        assert np.isclose(bs[0], inp.before_percell_best_s[i], rtol=0, atol=1e-9), int(i)


import gzip
import time
from fastapi import FastAPI
from fastapi.testclient import TestClient


def _api(inp_factory, keep=3):
    from server.cn_recall_api import register_cn_recall
    app = FastAPI()
    ex = ThreadPoolExecutor(2)
    register_cn_recall(app, lambda: (inp_factory(), {'hash': 'fake'}, ex), keep=keep)
    return TestClient(app)


def _wait(client, jid, timeout=120):
    t0 = time.time()
    while time.time() - t0 < timeout:
        v = client.get(f'/api/integration/recall/{jid}').json()
        if v['status'] != 'running':
            return v
        time.sleep(0.05)
    raise AssertionError('recall did not finish')


def test_recall_route_rejects_bad_rules():
    c = _api(lambda: synthetic_inputs(percell=False))
    assert c.post('/api/integration/recall', json={'rules': []}).status_code == 422
    assert c.post('/api/integration/recall', json={}).status_code == 422


def test_recall_lifecycle_view_and_percell_bin():
    inp = synthetic_inputs()
    c = _api(lambda: inp)
    r = c.post('/api/integration/recall', json={'rules': LOOSE, 'allow_gap_crossing': False})
    assert r.status_code == 200 and set(r.json()) == {'id', 'hash', 'status'}
    v = _wait(c, r.json()['id'])
    assert set(v) == {'id', 'hash', 'status', 'stage', 'progress', 'error', 'settings', 'provenance', 'merge',
                      'cluster_segments', 'diff', 'percell_ready'}
    assert v['status'] == 'done' and v['stage'] == 'done' and v['percell_ready'] is True
    assert v['settings'] == {'rules': LOOSE, 'allow_gap_crossing': False}
    assert v['provenance']['scope'] == ('Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on the '
                                        'rule-merged segmentation; originals unchanged.')
    assert v['provenance']['dataset'] == {'hash': 'fake'} and len(v['provenance']['implementation_sha256']) == 64
    b = c.get(f"/api/integration/recall/{v['id']}/percell.bin")
    assert b.status_code == 200 and len(b.content) == 12 * 70


def test_recall_reuses_job_by_hash():
    inp = synthetic_inputs()
    c = _api(lambda: inp)
    a = c.post('/api/integration/recall', json={'rules': LOOSE}).json()
    b = c.post('/api/integration/recall', json={'rules': LOOSE}).json()
    d = c.post('/api/integration/recall', json={'rules': LOOSE, 'allow_gap_crossing': True}).json()
    assert a['id'] == b['id'] and a['hash'] == b['hash'] and d['id'] != a['id']
    _wait(c, a['id']); _wait(c, d['id'])


def test_recall_unknown_id_404():
    c = _api(lambda: synthetic_inputs(percell=False))
    assert c.get('/api/integration/recall/nope').status_code == 404
    assert c.get('/api/integration/recall/nope/percell.bin').status_code == 404


def test_recall_missing_percell_publishes_cluster_then_errors():
    c = _api(lambda: synthetic_inputs(percell=False))
    jid = c.post('/api/integration/recall', json={'rules': NOOP}).json()['id']
    v = _wait(c, jid)
    assert v['status'] == 'error' and v['error'].startswith('percell')
    assert v['cluster_segments'] and v['diff']['sources'] and v['percell_ready'] is False
    assert c.get(f'/api/integration/recall/{jid}/percell.bin').status_code == 409


def test_recall_retention_evicts_oldest():
    inp = synthetic_inputs(percell=False)
    c = _api(lambda: inp, keep=2)
    ids = []
    for tau in (1e-9, 2e-9, 3e-9):
        ids.append(c.post('/api/integration/recall', json={'rules': [{**NOOP[0], 'tau': tau}]}).json()['id'])
        _wait(c, ids[-1])
    assert c.get(f'/api/integration/recall/{ids[0]}').status_code == 404
    assert c.get(f'/api/integration/recall/{ids[2]}').status_code == 200


def test_recall_failures_do_not_evict_done_job():
    good = synthetic_inputs()
    bad = synthetic_inputs(percell=False)
    state = {'inp': good}
    c = _api(lambda: state['inp'], keep=2)
    gid = c.post('/api/integration/recall', json={'rules': NOOP}).json()['id']
    assert _wait(c, gid)['status'] == 'done'
    state['inp'] = bad
    for tau in (2e-9, 3e-9, 4e-9, 5e-9):
        jid = c.post('/api/integration/recall', json={'rules': [{**NOOP[0], 'tau': tau}]}).json()['id']
        assert _wait(c, jid)['status'] == 'error'
    assert c.get(f'/api/integration/recall/{gid}').status_code == 200


def test_recall_unwrapped_exception_names_stage():
    inp = synthetic_inputs()
    inp.before_percell_u8 = np.zeros((12, 3), np.uint8)   # wrong width: raises a plain error in the unwrapped per-cell diff
    c = _api(lambda: inp)
    v = _wait(c, c.post('/api/integration/recall', json={'rules': NOOP}).json()['id'])
    assert v['status'] == 'error' and v['error'].startswith('percell: ')


# ---- merge rounds (2026-10-27 merge-rounds spec)
from server import call_noise as cn


def test_run_recall_publishes_noise():
    inp = synthetic_inputs()
    out, st, _ = _run(inp, NOOP)
    noise = st['diff']['noise']
    assert noise['cluster0']['genome']['cluster'] == cn.cluster_noise(inp.before_cluster_segments)['cluster0']['genome']
    segs = {s: cr.seg_labels(inp.var_seq, cr.boundaries_from_bp(inp.breakpoints, inp.sources, inp.chromosomes())[s]) for s in inp.sources}
    want = cn.cells_noise(inp.before_percell_u8, inp.cell_cluster, segs, inp.var_seq, inp.sources)
    assert noise['cluster0']['genome']['cells'] == want['cluster0']['genome']
    assert out['noise_cells'] == want and np.allclose(out['percell_best_s'], inp.before_percell_best_s)


def test_run_calls_reproduces_percell_on_the_given_breakpoints():
    inp = synthetic_inputs()
    log = []
    with ThreadPoolExecutor(2) as ex:
        out = cr.run_calls(inp, ex, lambda **kw: log.append(kw), chunk=5, max_inflight=2)
    assert np.array_equal(out['percell_u8'], inp.before_percell_u8)
    assert np.allclose(out['percell_best_s'], inp.before_percell_best_s)
    assert [kw['stage'] for kw in log if 'stage' in kw] == ['percell']
    last = {}
    for kw in log:
        last.update(kw)
    assert last['percell_ready'] is True and last['diff']['noise']['cluster0']['genome']['cells']['n_cells'] == 6


def test_run_recall_waits_for_base_percell():
    inp = synthetic_inputs()
    u8, bs = inp.before_percell_u8, inp.before_percell_best_s
    calls = []
    inp.before_percell_u8 = None
    inp.before_percell_wait = lambda: (calls.append(1), (u8, bs))[1]
    out, st, _ = _run(inp, NOOP)
    assert calls == [1] and st['diff']['cells']['n_changed_cells'] == 0


def test_run_recall_wait_failure_names_stage():
    inp = synthetic_inputs()
    inp.before_percell_u8 = None
    def boom():
        raise ValueError('nope')
    inp.before_percell_wait = boom
    with pytest.raises(cr.RecallError, match=r'^percell base: ValueError: nope'):
        _run(inp, NOOP)
