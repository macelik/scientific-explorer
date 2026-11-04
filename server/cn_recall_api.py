"""Genome-wide merge + CN re-call jobs (in memory, latest `keep` kept).
Spec: docs/superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md §3."""
from __future__ import annotations

import gzip
import hashlib
import json
import threading
import uuid
from collections import deque
from pathlib import Path
from dataclasses import replace
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from fastapi import HTTPException, Response
from pydantic import BaseModel

from . import call_noise as cn
from . import cn_recall as cr
from . import merge_rounds as mrd
from . import merge_rules as mr

_HERE = Path(__file__).resolve().parent
IMPLEMENTATION_SHA256 = hashlib.sha256(b''.join((_HERE / f).read_bytes() for f in ('merge_rules.py', 'hatch_merge.py', 'cn_recall.py', 'merge_rounds.py', 'call_noise.py'))).hexdigest()
SCOPE = ('Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on the rule-merged segmentation; '
         'originals unchanged.')
SCOPE_ROUNDS = ('Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on a rule-merged segmentation '
                '({n} earlier rounds); originals unchanged.')
VIEW_KEYS = ('id', 'hash', 'status', 'stage', 'progress', 'error', 'settings', 'provenance', 'merge',
             'cluster_segments', 'diff', 'percell_ready')


class RecallRequest(BaseModel):
    rules: List[Dict[str, Any]]
    allow_gap_crossing: bool = False
    rounds: List[Dict[str, Any]] = []


class BaseRequest(BaseModel):
    rounds: List[Dict[str, Any]] = []


class PinRequest(BaseModel):
    rounds: List[Dict[str, Any]]


def recall_hash(rules, allow_gap_crossing, identity, rounds=None) -> str:
    payload = {'rules': rules, 'allow_gap_crossing': allow_gap_crossing}
    if rounds:
        payload['rounds'] = rounds
    body = json.dumps(payload, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256((body + '|' + str(identity.get('hash')) + '|' + IMPLEMENTATION_SHA256).encode()).hexdigest()


def calls_hash(rounds, identity) -> str:
    return hashlib.sha256(('calls|' + mrd.rounds_hash(rounds) + '|' + str(identity.get('hash')) + '|' + IMPLEMENTATION_SHA256).encode()).hexdigest()


def _records(df) -> list:
    if df is None or not len(df):
        return []
    return df.astype(object).where(pd.notna(df), None).to_dict(orient='records')


class RecallManager:
    def __init__(self, keep: int = 3):
        self.keep = keep
        self.jobs: Dict[str, dict] = {}
        self.by_hash: Dict[str, str] = {}
        self.finished: deque = deque()      # done jobs
        self.failed: deque = deque()        # error jobs, budgeted separately so failures never evict a good result
        self.pinned: Optional[str] = None   # the current base's per-cell job; never evicted
        self.lock = threading.RLock()

    def _new(self, h, settings, scope, identity, result_key, stage):
        jid = uuid.uuid4().hex[:12]
        self.jobs[jid] = {'id': jid, 'hash': h, 'status': 'running', 'stage': stage, 'progress': None, 'error': None,
                          'settings': settings,
                          'provenance': {'dataset': identity, 'implementation_sha256': IMPLEMENTATION_SHA256, 'scope': scope},
                          'merge': None, 'cluster_segments': None, 'diff': None, 'percell_ready': False,
                          '_u8': None, '_gz': None, '_best': None, '_noise_cells': None,
                          '_result_key': result_key, '_event': threading.Event()}
        self.by_hash[h] = jid
        return jid

    def _reuse(self, h):
        jid = self.by_hash.get(h)
        return jid if jid and self.jobs[jid]['status'] != 'error' else None

    def start(self, inp, rules, allow_gap_crossing, identity, executor, rounds=()) -> dict:
        rounds = list(rounds)
        h = recall_hash(rules, allow_gap_crossing, identity, rounds)
        with self.lock:
            jid = self._reuse(h)
            if jid:
                return self.view(jid)
            settings = {'rules': rules, 'allow_gap_crossing': allow_gap_crossing}
            if rounds:
                settings['rounds'] = rounds
            scope = SCOPE_ROUNDS.format(n=len(rounds)) if rounds else SCOPE
            key = mrd.rounds_hash(rounds + [{'rules': rules, 'allow_gap_crossing': allow_gap_crossing}])
            jid = self._new(h, settings, scope, identity, key, 'merge')
        threading.Thread(target=self._run, args=(jid, lambda update: cr.run_recall(inp, rules, allow_gap_crossing, executor, update)),
                         daemon=True).start()
        return self.view(jid)

    def start_calls(self, inp, rounds, identity, executor) -> dict:
        """Calls-only per-cell job on inp.breakpoints (the base of `rounds`); pinned at once."""
        h = calls_hash(rounds, identity)
        with self.lock:
            jid = self._reuse(h)
            if not jid:
                jid = self._new(h, {'rules': [], 'allow_gap_crossing': False, 'rounds': list(rounds)},
                                SCOPE_ROUNDS.format(n=len(rounds)), identity, mrd.rounds_hash(list(rounds)), 'percell')
                threading.Thread(target=self._run, args=(jid, lambda update: cr.run_calls(inp, executor, update)), daemon=True).start()
            self.pin(jid)
        return self.view(jid)

    def _run(self, jid, fn):
        job = self.jobs[jid]

        def update(**kw):
            with self.lock:
                job.update(kw)
        try:
            out = fn(update)
            update(_u8=out['percell_u8'], _best=out.get('percell_best_s'), _noise_cells=out.get('noise_cells'),
                   status='done', stage='done')
        except Exception as e:  # noqa: BLE001
            if isinstance(e, cr.RecallError):
                msg = str(e)
            else:
                with self.lock:
                    stage = job['stage']
                msg = f'{stage}: {type(e).__name__}: {e}'
            update(status='error', error=msg)
        finally:
            with self.lock:
                if jid != self.pinned:
                    self._enqueue(jid)
            job['_event'].set()

    def _enqueue(self, jid):
        q = self.finished if self.jobs[jid]['status'] == 'done' else self.failed
        q.append(jid)
        while len(q) > self.keep:
            old = self.jobs.pop(q.popleft(), None)
            if old and self.by_hash.get(old['hash']) == old['id']:
                del self.by_hash[old['hash']]

    def pin(self, jid):
        with self.lock:
            old = self.pinned
            if old == jid:
                return
            self.pinned = jid
            for q in (self.finished, self.failed):
                if jid in q:
                    q.remove(jid)
            if old and old in self.jobs and self.jobs[old]['status'] != 'running':
                self._enqueue(old)

    def base_percell(self, key):
        """('ready' | 'running' | 'missing', job id, job) for the per-cell calls of the base with rounds_hash `key`.
        The job object is captured under the lock so later eviction from `jobs` cannot make it vanish."""
        with self.lock:
            p = self.jobs.get(self.pinned) if self.pinned else None
            if p and p['_result_key'] == key and p['status'] != 'error':
                return ('ready' if p['status'] == 'done' else 'running'), p['id'], p
            for j in self.jobs.values():
                if j['_result_key'] == key and j['status'] == 'done' and j['_u8'] is not None:
                    return 'ready', j['id'], j
            return 'missing', None, None

    def wait_percell(self, job):
        """Wait on the calls job object itself (not a `jobs` lookup), so eviction cannot make it 'gone'."""
        job['_event'].wait()
        if job['status'] != 'done' or job['_u8'] is None:
            raise cr.RecallError(f"percell base: {job['error']}")
        return job['_u8'], job['_best']

    def view(self, jid):
        with self.lock:
            job = self.jobs.get(jid)
            return None if job is None else {k: job[k] for k in VIEW_KEYS}

    def percell_gz(self, jid):
        with self.lock:
            job = self.jobs.get(jid)
            if job is None:
                raise KeyError(jid)
            if not job['percell_ready'] or job['_u8'] is None:
                return None
            if job['_gz'] is None:
                job['_gz'] = gzip.compress(np.ascontiguousarray(job['_u8']).tobytes(), compresslevel=4)
            return job['_gz']


def inputs_from_integration(ig, package_dir: str) -> cr.RecallInputs:
    import h5py
    from .integration import _csr_rows
    seg = ig.tables.get('karyogram_cluster_segments')
    if seg is None:
        raise HTTPException(503, 'karyogram_cluster_segments.tsv is missing; nothing to compare against')
    h5 = h5py.File(ig.h5ad_path, 'r')
    n_bins = ig.n_bins
    best = ig.cells['percell_best_s'].to_numpy(dtype=float) if 'percell_best_s' in ig.cells else np.full(len(ig.cell_ids), np.nan)
    return cr.RecallInputs(var_seq=np.asarray(ig.var_seq), var_start=ig.var_start, var_end=ig.var_end,
                           pb_x=ig.arrays['pb_x'], breakpoints=ig.tables['breakpoints'], members=ig.members,
                           cell_ids=[str(c) for c in ig.cell_ids],
                           cell_cluster=pd.to_numeric(ig.cells['wnn_leiden_0.3']).to_numpy().astype(int),
                           read_rows=lambda lo, hi: _csr_rows(h5['X'], lo, hi, n_bins),
                           before_cluster_segments=seg, before_percell_u8=ig.arrays.get('percell_cn_u8'),
                           before_percell_best_s=best, package_dir=package_dir, sources=list(ig.members))


def integration_context(get_integration, get_jobs, package_dir):
    def ctx():
        ig = get_integration()
        if ig is None or not getattr(ig, 'available', False):
            raise HTTPException(503, 'integration data not loaded')
        if ig.arrays is None:
            raise HTTPException(503, 'Pseudobulk index is not ready')
        jobs = get_jobs()
        if jobs is None:
            raise HTTPException(503, 'job pool not ready')
        return inputs_from_integration(ig, package_dir), ig.identity, jobs.background
    return ctx


def register_cn_recall(app, get_context, keep: int = 3):
    manager = RecallManager(keep=keep)

    caches = {'pb64': {}, 'original_noise': {}}

    def pb64(inp, identity):
        k = str(identity.get('hash'))
        if k not in caches['pb64']:
            caches['pb64'][k] = cr.pseudobulk(inp)
        return caches['pb64'][k]

    def original_noise(inp, identity):
        k = (str(identity.get('hash')), id(inp.breakpoints))
        if k not in caches['original_noise']:
            cells = None
            if inp.before_percell_u8 is not None:
                b = cr.boundaries_from_bp(inp.breakpoints, inp.sources, inp.chromosomes())
                segs = {s: cr.seg_labels(inp.var_seq, b[s]) for s in inp.sources}
                cells = cn.cells_noise(inp.before_percell_u8, inp.cell_cluster, segs, inp.var_seq, inp.sources)
            caches['original_noise'][k] = cn.combine(cn.cluster_noise(inp.before_cluster_segments), cells)
        return caches['original_noise'][k]

    def base_state(inp, rounds, identity, executor):
        """Cache entry for the base, plus its per-cell status; starts the calls-only recompute when missing."""
        try:
            e = mrd.CACHE.get(inp, rounds)
            if 'cluster_segments' not in e['extra']:
                pb = pb64(inp, identity)
                b = cr.boundaries_from_bp(e['breakpoints'], inp.sources, inp.chromosomes())
                rows = []
                for k, s in enumerate(inp.sources):
                    rows += cr.cluster_rows(s, len(inp.members[s]), pb[k], cr.seg_labels(inp.var_seq, b[s]), inp.var_seq)[0]
                e['extra']['cluster_segments'] = rows
                e['extra']['cluster_noise'] = cn.cluster_noise(pd.DataFrame(rows))
        except Exception as ex:  # noqa: BLE001
            raise HTTPException(500, f'base: {type(ex).__name__}: {ex}') from ex
        status, jid, job = manager.base_percell(e['hash'])
        if status == 'missing':
            jid = manager.start_calls(replace(inp, breakpoints=e['breakpoints']), rounds, identity, executor)['id']
            with manager.lock:
                job = manager.jobs[jid]
            status = 'running'
        elif status == 'ready':
            manager.pin(jid)
        return e, status, jid, job

    def _rounds_or_422(raw):
        try:
            return mrd.normalize_rounds(raw)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e

    @app.post('/api/integration/recall')
    def start_recall(req: RecallRequest):
        try:
            mr.validate_rules(req.rules)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e
        rounds = _rounds_or_422(req.rounds)
        inp, identity, executor = get_context()
        if rounds:
            e, status, jid, job = base_state(inp, rounds, identity, executor)
            u8 = best = None
            if status == 'ready':
                with manager.lock:
                    u8, best = job['_u8'], job['_best']
            inp = replace(inp, breakpoints=e['breakpoints'],
                          before_cluster_segments=pd.DataFrame(e['extra']['cluster_segments']),
                          before_percell_u8=u8,
                          before_percell_best_s=best if best is not None else np.full(len(inp.cell_ids), np.nan),
                          before_percell_wait=(lambda j=job: manager.wait_percell(j)) if status == 'running' else None)
        v = manager.start(inp, req.rules, req.allow_gap_crossing, identity, executor, rounds)
        return {'id': v['id'], 'hash': v['hash'], 'status': v['status']}

    @app.post('/api/integration/recall/{jid}/pin')
    def pin_recall(jid: str, req: PinRequest):
        rounds = _rounds_or_422(req.rounds)
        with manager.lock:
            job = manager.jobs.get(jid)
            if job is None:
                raise HTTPException(404, f'unknown re-call {jid!r}')
            if job['status'] != 'done':
                raise HTTPException(409, 'only a finished re-call can be adopted')
            if job['_result_key'] != mrd.rounds_hash(rounds):
                raise HTTPException(409, 'this re-call does not produce the given rounds (computed on a different base)')
            manager.pin(jid)
        return {'pinned': jid}

    @app.post('/api/integration/base')
    def base(req: BaseRequest):
        rounds = _rounds_or_422(req.rounds)
        inp, identity, executor = get_context()
        orig = original_noise(inp, identity)
        if not rounds:
            return {'hash': mrd.rounds_hash([]), 'n_rounds': 0, 'breakpoints': _records(inp.breakpoints), 'rounds': [],
                    'cluster_segments': _records(inp.before_cluster_segments), 'noise': orig, 'noise_original': orig,
                    'percell': {'status': 'original' if inp.before_percell_u8 is not None else 'missing', 'job': None}}
        e, status, jid, job = base_state(inp, rounds, identity, executor)
        cells = job['_noise_cells'] if status == 'ready' else None
        return {'hash': e['hash'], 'n_rounds': len(rounds), 'breakpoints': _records(e['breakpoints']), 'rounds': e['rounds'],
                'cluster_segments': e['extra']['cluster_segments'],
                'noise': cn.combine(e['extra']['cluster_noise'], cells), 'noise_original': orig,
                'percell': {'status': status, 'job': jid}}

    @app.get('/api/integration/recall/{jid}')
    def recall_status(jid: str):
        v = manager.view(jid)
        if v is None:
            raise HTTPException(404, f'unknown re-call {jid!r}')
        return v

    @app.get('/api/integration/recall/{jid}/percell.bin')
    def recall_percell(jid: str):
        try:
            gz = manager.percell_gz(jid)
        except KeyError:
            raise HTTPException(404, f'unknown re-call {jid!r}') from None
        if gz is None:
            raise HTTPException(409, 'per-cell re-call not ready')
        return Response(content=gz, media_type='application/octet-stream',
                        headers={'Content-Encoding': 'gzip', 'X-Encoding': 'uint8 = 2*CN state, 255 = missing; rows = h5mu obs order'})

    return manager
