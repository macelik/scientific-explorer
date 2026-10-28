"""Allelic pseudobulks (spec 2026-10-28-allelic-segmentation-design.md §2): load checks,
ΣHA/ΣHB per source streamed from the aligned parquet matrices, a disk cache, and the four signals."""
from __future__ import annotations

import hashlib
import os
import threading
from typing import Dict, List

import numpy as np
import pyarrow.parquet as pq

SIGNALS = ('dbaf', 'abs_dbaf', 'shrunk_abs_dbaf', 'oriented_all')
CODE_VERSION = 'allelic-data-1'
COORD_COLS = ('chr', 'start', 'end')
_NON_CELL = set(COORD_COLS) | {'__index_level_0__'}


class AllelicLoadError(Exception):
    pass


def _path(allelic_dir, name):
    return os.path.join(allelic_dir, f'{name}.parquet')


def matrix_columns(path) -> List[str]:
    return [c for c in pq.ParquetFile(path).schema_arrow.names if c not in _NON_CELL]


def read_coords(path):
    t = pq.read_table(path, columns=list(COORD_COLS))
    return [t.column(c).to_numpy() for c in COORD_COLS]


def signal_values(ha, hb, signal, min_at, ha_all=None, hb_all=None):
    """Per-bin signal (NaN where skipped) and the kept mask. Kept: AT >= min_at (and, for
    oriented_all, the all-cells pseudobulk has a direction)."""
    if signal not in SIGNALS:
        raise ValueError(f'signal: unknown {signal!r}')
    ha = np.asarray(ha, dtype=np.float64); hb = np.asarray(hb, dtype=np.float64)
    at = ha + hb
    kept = at >= min_at
    with np.errstate(divide='ignore', invalid='ignore'):
        d = (hb - ha) / (2 * at)
        if signal == 'dbaf':
            v = d
        elif signal == 'abs_dbaf':
            v = np.abs(d)
        elif signal == 'shrunk_abs_dbaf':
            v = np.abs(hb - ha) / (2 * (at + 2))
        else:
            sign = np.sign(np.asarray(hb_all, dtype=np.float64) - np.asarray(ha_all, dtype=np.float64))
            kept = kept & (sign != 0)
            v = sign * d
    return np.where(kept, v, np.nan), kept


def check_inputs(allelic_dir, var_seq, var_start, var_end, cell_ids) -> List[str]:
    missing = [n for n in ('HA', 'HB', 'AT') if not os.path.isfile(_path(allelic_dir, n))]
    if missing:
        return [f"allelic matrices not found in {allelic_dir}: {', '.join(m + '.parquet' for m in missing)} (set PYEPI_SCIENTIFIC_ALLELIC)"]
    problems = []
    want = [np.asarray(var_seq).astype(str), np.asarray(var_start).astype(np.int64), np.asarray(var_end).astype(np.int64)]
    for name in ('HA', 'HB'):
        got = read_coords(_path(allelic_dir, name))
        if len(got[0]) != len(want[0]) or not (np.array_equal(got[0].astype(str), want[0])
                                               and np.array_equal(got[1].astype(np.int64), want[1])
                                               and np.array_equal(got[2].astype(np.int64), want[2])):
            problems.append(f'{name}.parquet coordinates differ from the app bins ({len(got[0])} rows vs {len(want[0])})')
        cols = set(matrix_columns(_path(allelic_dir, name)))
        ids = set(map(str, cell_ids))
        if cols != ids:
            problems.append(f'{name}.parquet barcodes differ from the h5mu cells: {len(ids - cols)} missing, {len(cols - ids)} extra')
    return problems


def sum_by_source(allelic_dir, cell_ids, members: Dict[str, np.ndarray], sources, batch=200, at_check=200):
    """ΣHA, ΣHB per source (rows in `sources` order), streamed in barcode batches; columns are
    matched by barcode name, never by position. Checks AT == HA + HB on the first `at_check` cells."""
    cell_ids = [str(c) for c in cell_ids]
    n_bins = pq.ParquetFile(_path(allelic_dir, 'HA')).metadata.num_rows
    ha = np.zeros((len(sources), n_bins)); hb = np.zeros((len(sources), n_bins))
    src_rows = {s: i for i, s in enumerate(sources)}
    cell_src = {}
    for s, idx in members.items():
        if s != 'all' and s in src_rows:
            for i in idx:
                cell_src[int(i)] = src_rows[s]
    all_row = src_rows.get('all')
    for lo in range(0, len(cell_ids), batch):
        names = cell_ids[lo:lo + batch]
        a = pq.read_table(_path(allelic_dir, 'HA'), columns=names).to_pandas()[names].to_numpy(np.float64)
        b = pq.read_table(_path(allelic_dir, 'HB'), columns=names).to_pandas()[names].to_numpy(np.float64)
        if lo < at_check:
            m = names[:max(0, at_check - lo)]
            t = pq.read_table(_path(allelic_dir, 'AT'), columns=m).to_pandas()[m].to_numpy(np.float64)
            if not np.allclose(t, a[:, :len(m)] + b[:, :len(m)]):
                raise AllelicLoadError('AT.parquet differs from HA + HB')
        for j in range(len(names)):
            r = cell_src.get(lo + j)
            if r is not None:
                ha[r] += a[:, j]; hb[r] += b[:, j]
        if all_row is not None:
            ha[all_row] += a.sum(axis=1); hb[all_row] += b.sum(axis=1)
    return ha, hb


class AllelicData:
    def __init__(self, allelic_dir, cache_dir, var_seq, var_start, var_end, cell_ids, members, sources):
        self.allelic_dir, self.cache_dir = allelic_dir, cache_dir
        self.var_seq = np.asarray(var_seq).astype(str)
        self.var_start, self.var_end = np.asarray(var_start), np.asarray(var_end)
        self.cell_ids = [str(c) for c in cell_ids]
        self.members = {s: np.asarray(members[s]) for s in sources}
        self.sources = list(sources)
        self.chromosomes = list(dict.fromkeys(self.var_seq.tolist()))
        self._offset = {c: int(np.flatnonzero(self.var_seq == c)[0]) for c in self.chromosomes}
        self.state, self.problems, self.notes = 'idle', [], []
        self.ha = self.hb = None
        self.key = None
        self._lock = threading.Lock()
        self._done = threading.Event()

    def start(self):
        with self._lock:
            if self.state != 'idle':
                return
            self.state = 'building'
        threading.Thread(target=self._build, daemon=True).start()

    def wait(self, timeout=None):
        return self._done.wait(timeout)

    def _cache_key(self):
        h = hashlib.sha256(CODE_VERSION.encode())
        for name in ('HA', 'HB'):
            st = os.stat(_path(self.allelic_dir, name))
            h.update(f'{name}:{st.st_size}:{st.st_mtime_ns}'.encode())
        for s in self.sources:
            h.update(s.encode()); h.update(np.ascontiguousarray(self.members[s], dtype=np.int64).tobytes())
        h.update(','.join(self.cell_ids).encode())
        return h.hexdigest()[:20]

    def _build(self):
        try:
            problems = check_inputs(self.allelic_dir, self.var_seq, self.var_start, self.var_end, self.cell_ids)
            if problems:
                raise AllelicLoadError('; '.join(problems))
            key = self._cache_key()
            path = os.path.join(self.cache_dir, 'allelic', f'pseudobulk_{key}.npz')
            if os.path.isfile(path):
                z = np.load(path)
                ha, hb = z['ha'], z['hb']
                self.notes.append(f'pseudobulk sums loaded from {path}')
            else:
                ha, hb = sum_by_source(self.allelic_dir, self.cell_ids, self.members, self.sources)
                os.makedirs(os.path.dirname(path), exist_ok=True)
                tmp = path + '.tmp.npz'
                np.savez(tmp, ha=ha, hb=hb)
                os.replace(tmp, path)
            self.ha, self.hb, self.key = ha, hb, key
            self.state = 'ready'
        except Exception as e:  # noqa: BLE001
            self.problems.append(str(e) if isinstance(e, AllelicLoadError) else f'{type(e).__name__}: {e}')
            self.state = 'error'
        finally:
            self._done.set()

    def status(self) -> dict:
        return {'available': self.state == 'ready', 'state': self.state, 'problems': list(self.problems),
                'notes': list(self.notes), 'allelic_dir': self.allelic_dir, 'sources': self.sources}

    def offset(self, chrom) -> int:
        return self._offset[chrom]

    def track(self, source, chrom, signal, min_at) -> dict:
        si = self.sources.index(source)
        sl = np.flatnonzero(self.var_seq == chrom)
        ai = self.sources.index('all')
        ha, hb = self.ha[si, sl], self.hb[si, sl]
        v, kept = signal_values(ha, hb, signal, min_at, self.ha[ai, sl], self.hb[ai, sl])
        return {'values': v, 'kept': kept, 'at': ha + hb, 'ha': ha, 'hb': hb}
