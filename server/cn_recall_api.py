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
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from fastapi import HTTPException, Response
from pydantic import BaseModel

from . import cn_recall as cr
from . import merge_rules as mr

_HERE = Path(__file__).resolve().parent
IMPLEMENTATION_SHA256 = hashlib.sha256(b''.join((_HERE / f).read_bytes() for f in ('merge_rules.py', 'hatch_merge.py', 'cn_recall.py'))).hexdigest()
SCOPE = ('Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on the rule-merged segmentation; '
         'originals unchanged.')
VIEW_KEYS = ('id', 'hash', 'status', 'stage', 'progress', 'error', 'settings', 'provenance', 'merge',
             'cluster_segments', 'diff', 'percell_ready')


class RecallRequest(BaseModel):
    rules: List[Dict[str, Any]]
    allow_gap_crossing: bool = False


def recall_hash(rules, allow_gap_crossing, identity) -> str:
    body = json.dumps({'rules': rules, 'allow_gap_crossing': allow_gap_crossing}, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256((body + '|' + str(identity.get('hash')) + '|' + IMPLEMENTATION_SHA256).encode()).hexdigest()


class RecallManager:
    def __init__(self, keep: int = 3):
        self.keep = keep
        self.jobs: Dict[str, dict] = {}
        self.by_hash: Dict[str, str] = {}
        self.finished: deque = deque()
        self.lock = threading.RLock()

    def start(self, inp, rules, allow_gap_crossing, identity, executor) -> dict:
        h = recall_hash(rules, allow_gap_crossing, identity)
        with self.lock:
            jid = self.by_hash.get(h)
            if jid and self.jobs[jid]['status'] != 'error':
                return self.view(jid)
            jid = uuid.uuid4().hex[:12]
            self.jobs[jid] = {'id': jid, 'hash': h, 'status': 'running', 'stage': 'merge', 'progress': None, 'error': None,
                              'settings': {'rules': rules, 'allow_gap_crossing': allow_gap_crossing},
                              'provenance': {'dataset': identity, 'implementation_sha256': IMPLEMENTATION_SHA256, 'scope': SCOPE},
                              'merge': None, 'cluster_segments': None, 'diff': None, 'percell_ready': False,
                              '_u8': None, '_gz': None}
            self.by_hash[h] = jid
        threading.Thread(target=self._run, args=(jid, inp, rules, allow_gap_crossing, executor), daemon=True).start()
        return self.view(jid)

    def _run(self, jid, inp, rules, allow_gap_crossing, executor):
        job = self.jobs[jid]

        def update(**kw):
            with self.lock:
                job.update(kw)
        try:
            out = cr.run_recall(inp, rules, allow_gap_crossing, executor, update)
            update(_u8=out['percell_u8'], status='done', stage='done')
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
                self.finished.append(jid)
                while len(self.finished) > self.keep:
                    old = self.jobs.pop(self.finished.popleft(), None)
                    if old and self.by_hash.get(old['hash']) == old['id']:
                        del self.by_hash[old['hash']]

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

    @app.post('/api/integration/recall')
    def start_recall(req: RecallRequest):
        try:
            mr.validate_rules(req.rules)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e
        inp, identity, executor = get_context()
        v = manager.start(inp, req.rules, req.allow_gap_crossing, identity, executor)
        return {'id': v['id'], 'hash': v['hash'], 'status': v['status']}

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
