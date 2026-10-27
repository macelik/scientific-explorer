"""Call-noise statistics: zig-zags and gain-next-to-loss boundaries, per
chromosome and genome-wide, for cluster rows and per-cell calls.
Spec: docs/superpowers/specs/2026-10-27-merge-rounds-design.md §4."""
from __future__ import annotations

import math
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

MISSING = 255


def _finite(v) -> bool:
    return v is not None and math.isfinite(float(v))


def seq_noise(states, lengths) -> dict:
    """One chromosome's segment-state sequence; uncalled segments are dropped first."""
    kept = [(float(s), int(n)) for s, n in zip(states, lengths) if _finite(s)]
    v = [s for s, _ in kept]
    ln = [n for _, n in kept]
    zl = [ln[i] for i in range(1, len(v) - 1) if v[i] != v[i - 1] and v[i] != v[i + 1] and v[i - 1] == v[i + 1]]
    gla = sum(1 for a, b in zip(v[:-1], v[1:]) if (a > 1 and b < 1) or (a < 1 and b > 1))
    return {'segments': len(v), 'zigzags': len(zl), 'zigzag_lengths': zl, 'gain_loss_adjacent': gla}


def _median(xs) -> Optional[float]:
    return None if not xs else float(np.median(xs))


def _cluster_entry(segments, zigzags, lengths, gla) -> dict:
    return {'segments': segments, 'zigzags': zigzags, 'zigzag_median_bins': _median(lengths), 'gain_loss_adjacent': gla}


def cluster_noise(df) -> dict:
    """Per source (karyogram_cluster_segments columns): every chromosome, then 'genome'."""
    out: Dict[str, dict] = {}
    if df is None or not len(df):
        return out
    for g, d in df.groupby('group', sort=False):
        per, lens, tot = {}, [], [0, 0, 0]
        for c, e in d.sort_values('bin_start', kind='stable').groupby('chromosome', sort=False):
            r = seq_noise(e['cn_call'].tolist(), e['n_bins'].tolist())
            per[str(c)] = _cluster_entry(r['segments'], r['zigzags'], r['zigzag_lengths'], r['gain_loss_adjacent'])
            lens += r['zigzag_lengths']
            tot[0] += r['segments']; tot[1] += r['zigzags']; tot[2] += r['gain_loss_adjacent']
        per['genome'] = _cluster_entry(tot[0], tot[1], lens, tot[2])
        out[str(g)] = per
    return out


def _layout(seg: np.ndarray, var_seq: np.ndarray) -> Dict[str, tuple]:
    """Per chromosome (genome order): segment first-bin indices and lengths."""
    starts = np.r_[0, np.flatnonzero(np.diff(seg)) + 1]
    ends = np.r_[starts[1:], len(seg)]
    out: Dict[str, tuple] = {}
    for s, e in zip(starts.tolist(), ends.tolist()):
        if seg[s] < 0:
            continue
        f, ln = out.setdefault(str(var_seq[s]), ([], []))
        f.append(s); ln.append(e - s)
    return {c: (np.asarray(f), ln) for c, (f, ln) in out.items()}


def cells_noise(u8, cell_cluster, segs, var_seq, sources) -> dict:
    """Per-cell calls (uint8 = 2*state, 255 = not called), each cell under its own
    cluster's segmentation; 'all' pools every cell."""
    layouts = {s: _layout(np.asarray(segs[s]), np.asarray(var_seq)) for s in segs}
    acc: Dict[str, Dict[str, dict]] = {s: {} for s in sources}

    def add(src, key, segments, zigzags, lengths, gla):
        a = acc[src].setdefault(key, {'n': 0, 'segments': 0, 'zigzags': 0, 'gla': 0, 'lens': []})
        a['n'] += 1; a['segments'] += segments; a['zigzags'] += zigzags; a['gla'] += gla; a['lens'] += lengths

    for i in range(len(u8)):
        own = f'cluster{int(cell_cluster[i])}'
        if own not in layouts:
            continue
        targets = [t for t in (own, 'all') if t in acc]
        row = np.asarray(u8[i])
        tot = [0, 0, [], 0]
        for chrom, (first, lens) in layouts[own].items():
            codes = row[first].astype(float)
            r = seq_noise(np.where(codes == MISSING, np.nan, codes / 2).tolist(), lens)
            for t in targets:
                add(t, chrom, r['segments'], r['zigzags'], r['zigzag_lengths'], r['gain_loss_adjacent'])
            tot[0] += r['segments']; tot[1] += r['zigzags']; tot[2] += r['zigzag_lengths']; tot[3] += r['gain_loss_adjacent']
        for t in targets:
            add(t, 'genome', tot[0], tot[1], tot[2], tot[3])

    out = {}
    for src, d in acc.items():
        if not d:
            continue
        keys = [k for k in d if k != 'genome'] + (['genome'] if 'genome' in d else [])
        out[src] = {k: {'n_cells': d[k]['n'], 'segments_mean': round(d[k]['segments'] / d[k]['n'], 2),
                        'zigzags_mean': round(d[k]['zigzags'] / d[k]['n'], 2), 'zigzag_median_bins': _median(d[k]['lens']),
                        'gain_loss_adjacent_mean': round(d[k]['gla'] / d[k]['n'], 2)} for k in keys}
    return out


def combine(cluster, cells) -> dict:
    cells = cells or {}
    out = {}
    for src in dict.fromkeys([*cluster, *cells]):
        keys = list(dict.fromkeys([*cluster.get(src, {}), *cells.get(src, {})]))
        if 'genome' in keys:
            keys = [k for k in keys if k != 'genome'] + ['genome']
        out[src] = {k: {'cluster': cluster.get(src, {}).get(k), 'cells': cells.get(src, {}).get(k)} for k in keys}
    return out
