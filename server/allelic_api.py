"""Allelic segmentation routes and jobs (spec 2026-10-28-allelic-segmentation-design.md §3).
One job runs at a time; one request waits; finished results are cached in memory and on disk."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import threading
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from fastapi import HTTPException
from pydantic import BaseModel

from . import allelic_data as ald
from . import allelic_scan as sc

_HERE = Path(__file__).resolve().parent
IMPLEMENTATION_SHA256 = hashlib.sha256(b''.join((_HERE / f).read_bytes() for f in ('allelic_scan.py', 'allelic_data.py', 'allelic_api.py'))).hexdigest()
DEFAULT_SETTINGS = {'min_at': 10, **sc.DEFAULTS, 'global_pool': 'genome'}
POOLS = ('genome', 'chromosome')
LIMITS = {'min_at': (1, 1000, True), 'k': (1, 8, True), 'n_permutations': (100, 5000, True),
          'alpha': (0.0, 0.1, False), 'min_child_bins': (1, 500, True), 'min_terminal_bins': (0, 500, True)}
VIEW_KEYS = ('id', 'status', 'progress', 'error', 'signal', 'settings', 'result')


class SegmentRequest(BaseModel):
    signal: str = 'abs_dbaf'
    settings: Dict[str, Any] = {}


def normalize(signal, settings) -> dict:
    if signal not in ald.SIGNALS:
        raise ValueError(f"signal: must be one of {', '.join(ald.SIGNALS)}")
    for key in settings:
        if key not in LIMITS and key != 'global_pool':
            raise ValueError(f'{key}: unknown setting')
    out = {}
    for key, (lo, hi, integer) in LIMITS.items():
        v = settings.get(key, DEFAULT_SETTINGS[key])
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            raise ValueError(f'{key}: not a number')
        if integer and float(v) != int(v):
            raise ValueError(f'{key}: must be a whole number')
        if key == 'alpha':
            if not (lo < v <= hi):
                raise ValueError(f'{key}: must be > {lo} and ≤ {hi}')
        elif not (lo <= v <= hi):
            raise ValueError(f'{key}: must be ≥ {lo} and ≤ {hi}')
        out[key] = int(v) if integer else float(v)
    pool = settings.get('global_pool', DEFAULT_SETTINGS['global_pool'])
    if pool not in POOLS:
        raise ValueError(f"global_pool: must be {' or '.join(POOLS)}")
    out['global_pool'] = pool
    return out


def job_id(signal, settings, data_key) -> str:
    body = json.dumps({'signal': signal, 'settings': settings}, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(f'{body}|{data_key}|{IMPLEMENTATION_SHA256}'.encode()).hexdigest()[:16]


class AllelicManager:
    def __init__(self, cache_dir):
        self.cache_dir = os.path.join(cache_dir, 'allelic')
        self.jobs: Dict[str, dict] = {}
        self.running: Optional[str] = None
        self.queued: Optional[str] = None
        self.lock = threading.RLock()

    def _path(self, jid):
        return os.path.join(self.cache_dir, f'seg_{jid}.json.gz')

    def view(self, jid):
        with self.lock:
            j = self.jobs.get(jid)
            return None if j is None else {k: j[k] for k in VIEW_KEYS}

    def start(self, signal, settings, data, executor, package_dir) -> dict:
        jid = job_id(signal, settings, data.key)
        with self.lock:
            j = self.jobs.get(jid)
            if j and j['status'] != 'error':
                return self.view(jid)
            if os.path.isfile(self._path(jid)):
                with gzip.open(self._path(jid), 'rt') as fh:
                    result = json.load(fh)
                self.jobs[jid] = {'id': jid, 'status': 'done', 'progress': None, 'error': None,
                                  'signal': signal, 'settings': settings, 'result': result}
                return self.view(jid)
            self.jobs[jid] = {'id': jid, 'status': 'queued', 'progress': None, 'error': None, 'signal': signal,
                              'settings': settings, 'result': None, '_args': (data, executor, package_dir)}
            if self.running is None:
                self._launch(jid)
            else:
                if self.queued and self.queued != jid:
                    old = self.jobs[self.queued]
                    old.update(status='error', error='replaced by a newer request')
                    old.pop('_args', None)
                self.queued = jid
            return self.view(jid)

    def _launch(self, jid):
        self.running = jid
        self.jobs[jid]['status'] = 'running'
        threading.Thread(target=self._run, args=(jid,), daemon=True).start()

    def _run(self, jid):
        job = self.jobs[jid]
        data, executor, package_dir = job['_args']
        try:
            result = self._compute(job, data, executor, package_dir)
            os.makedirs(self.cache_dir, exist_ok=True)
            tmp = self._path(jid) + '.tmp'
            with gzip.open(tmp, 'wt') as fh:
                json.dump(result, fh)
            os.replace(tmp, self._path(jid))
            with self.lock:
                job.update(status='done', result=result)
        except Exception as e:  # noqa: BLE001
            with self.lock:
                job.update(status='error', error=str(e))
        finally:
            with self.lock:
                job.pop('_args', None)
                self.running = None
                nxt, self.queued = self.queued, None
                if nxt and self.jobs[nxt]['status'] == 'queued':
                    self._launch(nxt)

    def _compute(self, job, data, executor, package_dir) -> dict:
        signal, settings = job['signal'], job['settings']
        scan_settings = {k: settings[k] for k in sc.DEFAULTS}
        tasks = []
        for source in data.sources:
            tracks = {chrom: data.track(source, chrom, signal, settings['min_at']) for chrom in data.chromosomes}
            # original pyEpi: the global test samples from the whole genome (kept bins, genome order)
            pool = (np.concatenate([t['values'][t['kept']] for t in tracks.values()])
                    if settings['global_pool'] == 'genome' else None)
            for chrom, t in tracks.items():
                tasks.append((source, chrom, t['values'], t['kept'], t['ha'], t['hb'], scan_settings, package_dir, pool))
        futures = [executor.submit(sc.scan_task, t) for t in tasks]
        out = {'breakpoints': [], 'audit': [], 'segments': []}
        for i, (fut, t) in enumerate(zip(futures, tasks)):
            try:
                source, chrom, res = fut.result()
            except Exception as e:  # noqa: BLE001
                raise RuntimeError(f'{t[0]} {t[1]}: {type(e).__name__}: {e}') from e
            off = data.offset(chrom)
            for b in res['breakpoints']:
                g = off + b['bin']
                out['breakpoints'].append({'source': source, 'chromosome': chrom, 'absolute_bin': b['bin'],
                                           'bp_genomic_start': int(data.var_start[g]), 'bp_genomic_end': int(data.var_end[g]),
                                           'gap_start': b['gap_start'], 'gap_end': b['gap_end'],
                                           'recursion_depth': b['recursion_depth'], 'ad_score': b['ad_score'],
                                           'p_local': b['p_local'], 'p_global': b['p_global']})
            for r in res['audit']:
                out['audit'].append({**r, 'source': source, 'chromosome': chrom})
            for s in res['segments']:
                out['segments'].append({'source': source, 'chromosome': chrom, **s})
            with self.lock:
                job['progress'] = {'done': i + 1, 'total': len(tasks), 'source': source, 'chromosome': chrom}
        return out


def register_allelic(app, get_data, get_executor, package_dir, cache_dir) -> AllelicManager:
    manager = AllelicManager(cache_dir)

    def data_or_503():
        d = get_data()
        if d is None:
            raise HTTPException(503, 'integration data not loaded')
        d.start()
        return d

    def ready():
        d = data_or_503()
        if d.state != 'ready':
            raise HTTPException(503, f"allelic data {d.state}" + (': ' + '; '.join(d.problems) if d.problems else ''))
        return d

    @app.get('/api/allelic/status')
    def allelic_status():
        d = get_data()
        if d is None:
            return {'available': False, 'state': 'error', 'problems': ['integration data not loaded'], 'notes': [],
                    'allelic_dir': '', 'sources': []}
        d.start()
        return d.status()

    @app.get('/api/allelic/tracks')
    def allelic_tracks(chrom: str, sources: str, signal: str = 'abs_dbaf', min_at: int = 10):
        d = ready()
        try:
            s = normalize(signal, {'min_at': min_at})
        except ValueError as e:
            raise HTTPException(422, str(e))
        if chrom not in d.chromosomes:
            raise HTTPException(404, f'unknown chromosome {chrom}')
        names = [x for x in sources.split(',') if x]
        bad = [x for x in names if x not in d.sources]
        if bad:
            raise HTTPException(404, f"unknown source {', '.join(bad)}")
        tracks = []
        for src in names:
            t = d.track(src, chrom, signal, s['min_at'])
            tracks.append({'source': src, 'values': [None if v != v else round(float(v), 5) for v in t['values']],
                           'at': [int(a) for a in t['at']]})
        return {'chrom': chrom, 'signal': signal, 'min_at': s['min_at'], 'n': len(tracks[0]['values']) if tracks else 0, 'tracks': tracks}

    @app.post('/api/allelic/segment')
    def allelic_segment(req: SegmentRequest):
        try:
            settings = normalize(req.signal, req.settings)
        except ValueError as e:
            raise HTTPException(422, str(e))
        d = ready()
        ex = get_executor()
        if ex is None:
            raise HTTPException(503, 'job pool not ready')
        return manager.start(req.signal, settings, d, ex, package_dir)

    @app.get('/api/allelic/segment/{jid}')
    def allelic_job(jid: str):
        v = manager.view(jid)
        if v is None:
            raise HTTPException(404, f'unknown allelic job {jid}')
        return v

    return manager
