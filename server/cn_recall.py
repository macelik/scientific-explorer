"""Genome-wide rule merge and CN re-call on the merged segmentation.

Spec: docs/superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md.
Stage 1 runs merge_rules.run_rule_merge on every (source, chromosome) with the
same inputs as /api/integration/hatch-merge. Stages 2-3 re-call CN with pyEpi
assign_gainloss_new (unmodified) exactly as integration-story/scripts/
make_karyograms.py does: the cluster pseudobulk is the float64 mean of member
rows, and cells of cluster k use cluster k's segmentation. Nothing here writes
the original files.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from . import call_noise as cn
from . import merge_rules as mr
from . import pyepi_adapter as sci
from .config import Config

SOURCES = ['cluster0', 'cluster1', 'cluster2', 'cluster3', 'cluster4', 'all']
MISSING = 255


class RecallError(RuntimeError):
    """A stage failed; the message starts with the stage (and source/chromosome)."""


@dataclass
class RecallInputs:
    var_seq: np.ndarray
    var_start: np.ndarray
    var_end: np.ndarray
    pb_x: np.ndarray                      # (len(sources), n_bins) float32, merge profile (as hatch-merge)
    breakpoints: pd.DataFrame             # raw accepted breakpoints; absolute_bin is chromosome-local
    members: Dict[str, np.ndarray]        # source -> cell indices (row order of read_rows)
    cell_ids: List[str]
    cell_cluster: np.ndarray              # int cluster per cell
    read_rows: Callable[[int, int], np.ndarray]   # dense float64 count rows [lo, hi)
    before_cluster_segments: pd.DataFrame
    before_percell_u8: Optional[np.ndarray]
    before_percell_best_s: np.ndarray
    package_dir: str
    sources: List[str] = field(default_factory=lambda: list(SOURCES))
    before_percell_wait: Optional[Callable[[], tuple]] = None   # base per-cell (u8, best_s) still being computed

    def chromosomes(self) -> List[str]:
        return list(dict.fromkeys(self.var_seq.tolist()))


def _ensure(package_dir: Optional[str] = None) -> None:
    if not sci._loaded:
        sci.load(package_dir or Config.from_env().package_dir)


def merge_genome(inp: RecallInputs, rules, allow_gap_crossing=False, progress=None) -> dict:
    """Stage 1: merge every (source, chromosome) to completion. Kept boundaries
    are the starts of the final segments after the first (chromosome-local)."""
    chroms = inp.chromosomes()
    bp = inp.breakpoints
    pairs, kept = [], {}
    total, done = len(inp.sources) * len(chroms), 0
    for si, source in enumerate(inp.sources):
        kept[source] = {}
        for chrom in chroms:
            try:
                idx = np.flatnonzero(inp.var_seq == chrom)
                bounds = sorted(int(b) for b in bp.loc[(bp.source == source) & (bp.chromosome == chrom), 'absolute_bin'])
                res = mr.run_rule_merge(inp.pb_x[si, idx], bounds, len(idx), inp.var_start, inp.var_end,
                                        chrom_offset=int(idx[0]), rules=rules, allow_gap_crossing=allow_gap_crossing)
            except Exception as e:  # noqa: BLE001
                raise RecallError(f'merge {source}/{chrom}: {type(e).__name__}: {e}') from e
            kept[source][chrom] = [s['start'] for s in res['final'][1:]]
            steps = res['steps'][1:]
            pairs.append({'source': source, 'chrom': chrom, 'original': len(res['original']), 'merged': len(res['final']),
                          'merges': len(steps), 'phase2': sum(1 for s in steps if s['phase'] == 2)})
            done += 1
            if progress:
                progress(done, total)
    keep = {(s, c, b) for s, d in kept.items() for c, bs in d.items() for b in bs}
    mask = np.array([(r.source, r.chromosome, int(r.absolute_bin)) in keep for r in bp.itertuples(index=False)], dtype=bool)
    return {'pairs': pairs, 'boundaries': kept, 'merged_breakpoints': bp[mask] if len(bp) else bp}


def seg_labels(var_seq: np.ndarray, boundaries_by_chrom: Dict[str, List[int]]) -> np.ndarray:
    """Segment id per bin from chromosome-local cuts (make_karyograms.build_seg_labels)."""
    seg = np.full(len(var_seq), -1, dtype=int)
    sid = 0
    for ch in dict.fromkeys(var_seq.tolist()):
        cidx = np.flatnonzero(var_seq == ch)
        cuts = [0] + sorted(int(b) for b in boundaries_by_chrom.get(ch, [])) + [len(cidx)]
        for s, e in zip(cuts[:-1], cuts[1:]):
            if e > s:
                seg[cidx[s]:cidx[e - 1] + 1] = sid
                sid += 1
    return seg


def boundaries_from_bp(bp: pd.DataFrame, sources, chroms) -> Dict[str, Dict[str, List[int]]]:
    """Chromosome-local boundaries per source/chromosome from a breakpoints table."""
    out = {s: {c: [] for c in chroms} for s in sources}
    for r in bp.itertuples(index=False):
        if r.source in out and r.chromosome in out[r.source]:
            out[r.source][r.chromosome].append(int(r.absolute_bin))
    return {s: {c: sorted(v) for c, v in d.items()} for s, d in out.items()}


def pseudobulk(inp: RecallInputs, chunk=250, progress=None) -> np.ndarray:
    """Float64 mean of member rows per source, one pass over the count rows."""
    n = len(inp.cell_ids)
    sums = np.zeros((len(inp.sources), len(inp.var_seq)))
    masks = []
    for s in inp.sources:
        m = np.zeros(n, dtype=bool); m[inp.members[s]] = True; masks.append(m)
    for lo in range(0, n, chunk):
        hi = min(n, lo + chunk)
        rows = inp.read_rows(lo, hi)
        for k, m in enumerate(masks):
            sel = m[lo:hi]
            if sel.any():
                sums[k] += rows[sel].sum(axis=0)
        if progress:
            progress(hi, n)
    counts = np.array([len(inp.members[s]) for s in inp.sources], dtype=float)
    return sums / counts[:, None]


def encode_states(v) -> np.ndarray:
    """uint8 = rint(2*state) for finite states, 255 otherwise (IntegrationStore._build_arrays)."""
    v = np.asarray(v, dtype=np.float32)
    out = np.full(v.shape, MISSING, dtype=np.uint8)
    ok = np.isfinite(v)
    out[ok] = np.rint(v[ok] * 2).astype(np.uint8)
    return out


def cluster_rows(group: str, n_cells: int, x: np.ndarray, seg: np.ndarray, var_seq: np.ndarray):
    """One karyogram_cluster_segments.tsv row per segment, as make_karyograms.py
    writes them (continuous_cn is the rounded combined state there, kept as is)."""
    _ensure()
    _, _, holmes, watson, combined, best_s = sci.assign_gainloss_new(np.asarray(x, dtype=float), seg)
    combined = np.asarray(combined, dtype=float)
    rows = []
    for s in np.unique(seg):
        idx = np.flatnonzero(seg == s)
        rows.append({'group': group, 'n_cells': int(n_cells), 'best_s': round(float(best_s), 4), 'segment': int(s),
                     'chromosome': str(var_seq[idx[0]]), 'bin_start': int(idx[0]), 'bin_end': int(idx[-1]),
                     'n_bins': int(len(idx)), 'continuous_cn': round(float(combined[idx[0]]), 2),
                     'integer_cn': int(holmes[idx[0]]), 'watson': int(watson[idx[0]]), 'cn_call': float(combined[idx[0]])})
    return rows, combined


def cluster_u8_from_segments(df: pd.DataFrame, group: str, n_bins: int) -> np.ndarray:
    out = np.full(n_bins, MISSING, dtype=np.uint8)
    if len(df):
        for r in df[df['group'] == group].itertuples(index=False):
            out[int(r.bin_start):int(r.bin_end) + 1] = int(round(float(r.cn_call) * 2))
    return out


def call_cells(rows: np.ndarray, seg: np.ndarray, package_dir: str):
    """Worker: per-cell assign_gainloss_new on one cluster's segmentation."""
    _ensure(package_dir)
    u8 = np.full((len(rows), len(seg)), MISSING, dtype=np.uint8)
    best = np.full(len(rows), np.nan)
    for r in range(len(rows)):
        *_, combined, best_s = sci.assign_gainloss_new(np.asarray(rows[r], dtype=float), seg)
        u8[r] = encode_states(combined)
        best[r] = float(best_s)
    return u8, best


def add_transitions(acc: np.ndarray, before: np.ndarray, after: np.ndarray) -> None:
    ch = before != after
    if ch.any():
        acc += np.bincount(before[ch].astype(np.int64) * 256 + after[ch].astype(np.int64), minlength=65536)


def _state(code: int):
    return None if code == MISSING else code / 2


def transitions_list(acc: np.ndarray) -> List[dict]:
    nz = np.flatnonzero(acc)
    out = [(int(acc[c]), int(c)) for c in nz]
    out.sort(key=lambda t: (-t[0], t[1]))
    return [{'before': _state(c // 256), 'after': _state(c % 256), 'bins': n} for n, c in out]


def _f(v):
    return None if v is None or not np.isfinite(v) else float(v)


def _source_diff(source, pairs, before_df, best_after, before, after) -> dict:
    mine = [p for p in pairs if p['source'] == source]
    b = before_df[before_df['group'] == source]['best_s'] if len(before_df) else pd.Series(dtype=float)
    called = int((before != MISSING).sum())
    changed = int((before != after).sum())
    return {'source': source, 'segments_before': sum(p['original'] for p in mine), 'segments_after': sum(p['merged'] for p in mine),
            'merges': sum(p['merges'] for p in mine), 'phase2': sum(p['phase2'] for p in mine),
            'best_s_before': _f(b.iloc[0]) if len(b) else None, 'best_s_after': best_after,
            'bins_changed': changed, 'bins_called': called, 'pct_changed': round(100 * changed / called, 3) if called else 0.0}


def call_all_cells(inp: RecallInputs, segs, executor, update, chunk=100, max_inflight=8):
    """Stage 3: per-cell assign_gainloss_new, each cell under its cluster's segmentation."""
    n_bins, n = len(inp.var_seq), len(inp.cell_ids)
    out_u8 = np.full((n, n_bins), MISSING, dtype=np.uint8)
    best = np.full(n, np.nan)
    pending: deque = deque()
    done = 0

    def collect(item):
        nonlocal done
        idx, fut = item
        try:
            u8, bs = fut.result()
        except Exception as e:  # noqa: BLE001
            raise RecallError(f'percell {inp.cell_ids[int(idx[0])]}: {type(e).__name__}: {e}') from e
        out_u8[idx] = u8
        best[idx] = bs
        done += len(idx)
        update(progress={'done': done, 'total': n})

    for lo in range(0, n, chunk):
        hi = min(n, lo + chunk)
        cl = np.asarray(inp.cell_cluster[lo:hi])
        for k in np.unique(cl):
            if f'cluster{int(k)}' not in segs:
                raise RecallError(f'percell: cluster {int(k)} has no source')
        x = inp.read_rows(lo, hi)
        for k in np.unique(cl):
            sel = np.flatnonzero(cl == k)
            while len(pending) >= max_inflight:
                collect(pending.popleft())
            pending.append((lo + sel, executor.submit(call_cells, x[sel], segs[f'cluster{int(k)}'], inp.package_dir)))
    while pending:
        collect(pending.popleft())
    return out_u8, best


def _segs_from_bp(inp: RecallInputs):
    b = boundaries_from_bp(inp.breakpoints, inp.sources, inp.chromosomes())
    return {s: seg_labels(inp.var_seq, b[s]) for s in inp.sources}


def run_calls(inp: RecallInputs, executor, update, chunk=100, max_inflight=8) -> dict:
    """Calls-only job: per-cell calls on inp.breakpoints as given (a base), no merge, no diff."""
    segs = _segs_from_bp(inp)
    update(stage='percell', progress={'done': 0, 'total': len(inp.cell_ids)})
    u8, best = call_all_cells(inp, segs, executor, update, chunk, max_inflight)
    cells = cn.cells_noise(u8, inp.cell_cluster, segs, inp.var_seq, inp.sources)
    update(diff={'sources': [], 'cells': None, 'transitions': {'clusters': [], 'cells': None}, 'noise': cn.combine({}, cells)},
           percell_ready=True)
    return {'percell_u8': u8, 'percell_best_s': best, 'noise_cells': cells}


def run_recall(inp: RecallInputs, rules, allow_gap_crossing, executor, update, chunk=100, max_inflight=8) -> dict:
    """Stages 1-3 plus the diff. Publishes through update(**fields); returns
    {'percell_u8': ...}. Raises RecallError naming the failing stage."""
    n_bins, n = len(inp.var_seq), len(inp.cell_ids)
    update(stage='merge', progress={'done': 0, 'total': len(inp.sources) * len(inp.chromosomes())})
    merge = merge_genome(inp, rules, allow_gap_crossing, progress=lambda d, t: update(progress={'done': d, 'total': t}))
    segs = {s: seg_labels(inp.var_seq, merge['boundaries'][s]) for s in inp.sources}

    update(stage='cluster', progress={'done': 0, 'total': n})
    try:
        pb = pseudobulk(inp, progress=lambda d, t: update(progress={'done': d, 'total': t}))
        _ensure(inp.package_dir)
    except Exception as e:  # noqa: BLE001
        raise RecallError(f'cluster pseudobulk: {type(e).__name__}: {e}') from e
    rows, src_diffs, acc = [], [], np.zeros(65536, np.int64)
    for k, source in enumerate(inp.sources):
        try:
            r, combined = cluster_rows(source, len(inp.members[source]), pb[k], segs[source], inp.var_seq)
        except Exception as e:  # noqa: BLE001
            raise RecallError(f'cluster {source}: {type(e).__name__}: {e}') from e
        rows += r
        before = cluster_u8_from_segments(inp.before_cluster_segments, source, n_bins)
        after = encode_states(combined)
        if source != 'all':
            add_transitions(acc, before, after)
        src_diffs.append(_source_diff(source, merge['pairs'], inp.before_cluster_segments, r[0]['best_s'] if r else None, before, after))
    cluster_nz = cn.cluster_noise(pd.DataFrame(rows))
    diff = {'sources': src_diffs, 'cells': None, 'transitions': {'clusters': transitions_list(acc), 'cells': None},
            'noise': cn.combine(cluster_nz, None)}
    bp_records = merge['merged_breakpoints'].astype(object).where(pd.notna(merge['merged_breakpoints']), None).to_dict(orient='records')
    update(merge={'pairs': merge['pairs'], 'merged_breakpoints': bp_records}, cluster_segments=rows, diff=diff)

    before_u8, before_best = inp.before_percell_u8, inp.before_percell_best_s
    if before_u8 is None and inp.before_percell_wait is not None:
        update(stage='percell', progress=None)
        try:
            before_u8, before_best = inp.before_percell_wait()
        except RecallError:
            raise
        except Exception as e:  # noqa: BLE001
            raise RecallError(f'percell base: {type(e).__name__}: {e}') from e
    if before_u8 is None:
        raise RecallError('percell: the original per-cell CN matrix is missing')
    update(stage='percell', progress={'done': 0, 'total': n})
    out_u8, best = call_all_cells(inp, segs, executor, update, chunk, max_inflight)

    cacc = np.zeros(65536, np.int64)
    per_cell, pcts, dbs = [], [], []
    for i in range(n):
        b, a = np.asarray(before_u8[i]), out_u8[i]
        add_transitions(cacc, b, a)
        changed, called = int((b != a).sum()), int((b != MISSING).sum())
        pcts.append(100 * changed / called if called else 0.0)
        bb, ba = _f(before_best[i]), _f(best[i])
        if bb is not None and ba is not None:
            dbs.append(ba - bb)
        per_cell.append({'cell': inp.cell_ids[i], 'wnn_cluster': int(inp.cell_cluster[i]), 'best_s_before': bb,
                         'best_s_after': ba, 'bins_changed': changed})
    cells_diff = {'n_cells': n, 'n_changed_cells': sum(1 for r in per_cell if r['bins_changed']),
                     'median_pct_changed': round(float(np.median(pcts)), 3) if pcts else 0.0,
                     'best_s_change': {'median': _f(np.median(dbs)) if dbs else None, 'min': _f(min(dbs)) if dbs else None,
                                       'max': _f(max(dbs)) if dbs else None},
                     'per_cell': per_cell}
    cells_nz = cn.cells_noise(out_u8, inp.cell_cluster, segs, inp.var_seq, inp.sources)
    diff = {**diff, 'cells': cells_diff, 'transitions': {**diff['transitions'], 'cells': transitions_list(cacc)},
            'noise': cn.combine(cluster_nz, cells_nz)}
    update(diff=diff, percell_ready=True)
    return {'percell_u8': out_u8, 'percell_best_s': best, 'noise_cells': cells_nz}
