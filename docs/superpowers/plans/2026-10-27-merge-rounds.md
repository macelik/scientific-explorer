# Merge Rounds Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let the user adopt a genome re-call's merged segmentation as the new working base (repeatable, undoable), have all of section 2 read from that base, and show a "Call noise by round" table under the chromosome tracks.

**Architecture:**
- The base is a recipe `rounds = [{rules, allow_gap_crossing}]`, folded on the server over the original breakpoints with the existing `cn_recall.merge_genome`, and cached.
- Existing routes take an optional `rounds`, and a new `/api/integration/base` returns the base's breakpoints, cluster calls and noise.
- A base's per-cell calls come from a pinned job: either the adopted re-call, or a calls-only recompute after a restart.
- The client keeps `rounds` in the integration store and in sessions.

**Tech Stack:** FastAPI + numpy/pandas (server), pyEpi `assign_gainloss_new` via `server/pyepi_adapter.py` (unmodified), React 18 + TypeScript + Zustand + Vite (web), pytest, esbuild-bundled TS helper suites under Node, Playwright walkthroughs.

**Spec:** `docs/superpowers/specs/2026-10-27-merge-rounds-design.md`. Read it before starting any task.

## Global Constraints

- **Originals are never written.** pyEpi calling is unchanged.
- **Without `rounds` nothing changes.** With `rounds` left out or `[]`:
  - `/hatch-scores`, `/hatch-merge` and `/recall` return exactly what they return today;
  - the re-call hash formula, `settings` and the scope sentence are unchanged;
  - every existing test passes unmodified.
- **Merge signal:** every round merges on `pb_x` (the original pseudobulk). Only the boundaries carry over.
- **Invalid round:** gives HTTP 422 with message `round <i>: <validate_rules message>`, where i counts from 1.
- **Pinning:** exactly one pinned job at a time. A pinned job is never evicted.
- **Noise definitions,** shown as visible text under the table, verbatim:
  - "Calls are read per chromosome as the sequence of segment states. Segments with no call are skipped. Gain = state > 1 (weak gain included), loss = state < 1."
  - "Zig-zag: an internal segment whose state differs from both neighbours, while the two neighbours have the same state as each other."
  - "Gain next to loss: a boundary with a gain on one side and a loss on the other. This can be real (for example an unbalanced rearrangement), so it is counted but never merged automatically."
  - "Cell values: the mean per cell over all cells of the cluster. The zig-zag median length is taken over all of those cells' zig-zags. For all, the cells are every cell, each read under its own cluster's segmentation."
- **Scope sentence with rounds:** `Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on a rule-merged segmentation ({n} earlier rounds); originals unchanged.`
- **Commits:**
  - Stage files by name only.
  - Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - Never merge or push.
- **Servers:** never touch the user's server on port 8766. Browser checks use a private server, `python3 run.py --port 8799`, and you stop only that one.
- **Hands off:** leave `docs/superpowers/plans/2026-09-21-typesafe-jev-judgment-layers.md` and `docs/superpowers/plans/2026-09-30-portfolio-repository-polish.md` untouched.

## Review Focus

1. **Stale base data after undo or adopt.** A hatch-scores or `/base` response for the old recipe that arrives late must never be shown. Expected: every async write is guarded by the recipe signature it was requested for.
2. **Session restore with rounds that no longer validate** (for example older rule shapes). Expected: the rounds reset to `[]`, with the notice "Saved merge rounds could not be applied: …", and the app keeps working.
3. **Re-call on a base whose per-cell calls are missing** (server restarted). Expected: the job waits for the calls-only recompute and finishes with a real diff. It never compares against the original matrix.
4. **Adopting a re-call computed on a different base.** Expected: the button is disabled and the server refuses with 409.
5. **Eviction pressure.** Many re-calls after adopting must not evict the pinned job, or the base's per-cell matrix would disappear.

Tests for 1 are in Task 5 (signature guards, by inspection in review) and Task 7 (walkthrough undo). Task 5 covers 2. Task 4 covers 3, 4 and 5.

---

### Task 1: Call-noise statistics (`server/call_noise.py`)

**Files:**
- Create: `server/call_noise.py`
- Test: `tests/test_call_noise.py`

**Interfaces:**
- Produces:
  - `seq_noise(states, lengths) -> {'segments': int, 'zigzags': int, 'zigzag_lengths': list[int], 'gain_loss_adjacent': int}`
  - `cluster_noise(df) -> {source: {chrom | 'genome': ClusterNoise}}`, where `ClusterNoise = {'segments', 'zigzags', 'zigzag_median_bins', 'gain_loss_adjacent'}`
  - `cells_noise(u8, cell_cluster, segs, var_seq, sources) -> {source: {chrom | 'genome': CellsNoise}}`, where `CellsNoise = {'n_cells', 'segments_mean', 'zigzags_mean', 'zigzag_median_bins', 'gain_loss_adjacent_mean'}`
  - `combine(cluster, cells) -> {source: {key: {'cluster': ClusterNoise | None, 'cells': CellsNoise | None}}}`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_call_noise.py
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from server import call_noise as cn
from server import cn_recall as cr

NAN = float('nan')


def test_zigzag_counted_with_its_length():
    r = cn.seq_noise([1, 2, 1], [10, 3, 10])
    assert r == {'segments': 3, 'zigzags': 1, 'zigzag_lengths': [3], 'gain_loss_adjacent': 0}


def test_staircase_is_not_a_zigzag():
    assert cn.seq_noise([0, 1, 2], [5, 5, 5])['zigzags'] == 0
    assert cn.seq_noise([0, 2, 1], [5, 5, 5])['zigzags'] == 0   # neighbours differ


def test_weak_gain_zigzag():
    assert cn.seq_noise([1, 1.5, 1], [9, 2, 9])['zigzags'] == 1


def test_gain_next_to_loss():
    assert cn.seq_noise([2, 0], [5, 5])['gain_loss_adjacent'] == 1
    assert cn.seq_noise([1.5, 0.5, 1], [5, 5, 5])['gain_loss_adjacent'] == 1
    assert cn.seq_noise([1, 2, 1, 0], [5, 5, 5, 5])['gain_loss_adjacent'] == 0


def test_uncalled_segments_are_skipped():
    r = cn.seq_noise([1, NAN, 2, None, 1], [5, 1, 4, 1, 5])
    assert r == {'segments': 3, 'zigzags': 1, 'zigzag_lengths': [4], 'gain_loss_adjacent': 0}


def test_edges_and_empty():
    assert cn.seq_noise([2, 1], [3, 3])['zigzags'] == 0
    assert cn.seq_noise([2], [3])['zigzags'] == 0
    assert cn.seq_noise([], []) == {'segments': 0, 'zigzags': 0, 'zigzag_lengths': [], 'gain_loss_adjacent': 0}


def _df():
    rows = []
    def seg(g, c, s, n, call):
        rows.append({'group': g, 'chromosome': c, 'bin_start': s, 'bin_end': s + n - 1, 'n_bins': n, 'cn_call': call})
    seg('cluster0', 'chr1', 0, 4, 1); seg('cluster0', 'chr1', 4, 2, 2); seg('cluster0', 'chr1', 6, 4, 1)
    seg('cluster0', 'chr2', 10, 3, 2); seg('cluster0', 'chr2', 13, 3, 0)
    seg('all', 'chr1', 0, 10, 1); seg('all', 'chr2', 10, 6, 1)
    return pd.DataFrame(rows)


def test_cluster_noise_per_chromosome_and_genome():
    out = cn.cluster_noise(_df())
    assert list(out) == ['cluster0', 'all']
    assert out['cluster0']['chr1'] == {'segments': 3, 'zigzags': 1, 'zigzag_median_bins': 2.0, 'gain_loss_adjacent': 0}
    assert out['cluster0']['chr2'] == {'segments': 2, 'zigzags': 0, 'zigzag_median_bins': None, 'gain_loss_adjacent': 1}
    assert out['cluster0']['genome'] == {'segments': 5, 'zigzags': 1, 'zigzag_median_bins': 2.0, 'gain_loss_adjacent': 1}
    assert list(out['cluster0']) == ['chr1', 'chr2', 'genome']
    assert cn.cluster_noise(pd.DataFrame()) == {}


def test_cells_noise_equals_seq_noise_per_cell():
    var_seq = np.array(['chr1'] * 6 + ['chr2'] * 4)
    segs = {'cluster0': cr.seg_labels(var_seq, {'chr1': [2, 4], 'chr2': [2]}),
            'cluster1': cr.seg_labels(var_seq, {'chr1': [3], 'chr2': []}),
            'all': cr.seg_labels(var_seq, {})}
    # u8 = 2*state; 255 = not called
    u8 = np.array([[2, 2, 4, 4, 2, 2, 4, 4, 0, 0],      # cluster0: chr1 1,2,1 (zig-zag, 2 bins); chr2 2,0 (gain|loss)
                   [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],      # cluster0: flat
                   [0, 0, 0, 255, 255, 255, 2, 2, 2, 2]], dtype=np.uint8)   # cluster1: chr1 0, missing; chr2 1
    out = cn.cells_noise(u8, np.array([0, 0, 1]), segs, var_seq, ['cluster0', 'cluster1', 'all'])
    assert out['cluster0']['chr1'] == {'n_cells': 2, 'segments_mean': 3.0, 'zigzags_mean': 0.5, 'zigzag_median_bins': 2.0, 'gain_loss_adjacent_mean': 0.0}
    assert out['cluster0']['genome'] == {'n_cells': 2, 'segments_mean': 5.0, 'zigzags_mean': 0.5, 'zigzag_median_bins': 2.0, 'gain_loss_adjacent_mean': 0.5}
    assert out['cluster1']['genome'] == {'n_cells': 1, 'segments_mean': 2.0, 'zigzags_mean': 0.0, 'zigzag_median_bins': None, 'gain_loss_adjacent_mean': 0.0}
    assert out['all']['genome']['n_cells'] == 3 and out['all']['genome']['segments_mean'] == 4.0


def test_combine_pairs_levels():
    c = {'cluster0': {'chr1': {'segments': 1}, 'genome': {'segments': 1}}}
    k = {'cluster0': {'chr1': {'n_cells': 2}, 'genome': {'n_cells': 2}}}
    assert cn.combine(c, k) == {'cluster0': {'chr1': {'cluster': {'segments': 1}, 'cells': {'n_cells': 2}},
                                             'genome': {'cluster': {'segments': 1}, 'cells': {'n_cells': 2}}}}
    assert cn.combine(c, None)['cluster0']['chr1'] == {'cluster': {'segments': 1}, 'cells': None}


SEGS = Path(__file__).resolve().parents[2] / 'integration-story' / 'data' / 'derived' / 'karyogram_cluster_segments.tsv'


@pytest.mark.skipif(not SEGS.exists(), reason='integration-story data not available')
def test_real_original_cluster_noise_matches_measured():
    out = cn.cluster_noise(pd.read_csv(SEGS, sep='\t'))
    g = out['cluster2']['genome']
    assert (g['segments'], g['zigzags'], g['gain_loss_adjacent']) == (144, 33, 6)
    assert [out[s]['genome']['zigzags'] for s in ('cluster0', 'cluster1', 'cluster3', 'cluster4', 'all')] == [33, 33, 40, 33, 37]
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest tests/test_call_noise.py -q`
Expected: FAIL with `ImportError: cannot import name 'call_noise'`.

- [ ] **Step 3: Write the implementation**

```python
# server/call_noise.py
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python3 -m pytest tests/test_call_noise.py -q`
Expected: PASS, with the real-data test skipped only when the data is absent. On this machine the data is present and the test must pass.

- [ ] **Step 5: Commit**

```bash
git add server/call_noise.py tests/test_call_noise.py
git commit -m "feat: call-noise statistics (zig-zags, gain next to loss)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Rounds folding and the hatch routes (`server/merge_rounds.py`)

**Files:**
- Create: `server/merge_rounds.py`
- Modify: `server/cn_recall.py` (add `boundaries_from_bp`)
- Modify: `server/hatch_merge_api.py` (requests gain `rounds`; `_resolve` reads the base)
- Test: `tests/test_merge_rounds.py`

**Interfaces:**
- Consumes: `cn_recall.merge_genome(inp, rules, allow_gap_crossing)` → `{pairs, boundaries, merged_breakpoints}`, where `inp` needs `var_seq, var_start, var_end, pb_x, breakpoints, sources, chromosomes()`.
- Produces:
  - `cn_recall.boundaries_from_bp(bp, sources, chroms) -> {source: {chrom: sorted list[int]}}`
  - `merge_rounds.normalize_rounds(rounds) -> list[{'rules', 'allow_gap_crossing'}]`, raising `ValueError('round <i>: …')`
  - `merge_rounds.rounds_hash(rounds) -> str`, the sha256 hex of canonical JSON
  - `merge_rounds.BaseInputs` (dataclass: `var_seq, var_start, var_end, pb_x, breakpoints, sources`; `chromosomes()`)
  - `merge_rounds.fold(inp, rounds) -> (breakpoints_df, [{'pairs', 'merges'}])`
  - `merge_rounds.BaseCache.get(inp, rounds) -> entry`, where `entry = {'source', 'hash', 'breakpoints', 'rounds', 'extra': {}}`. Entries are keyed by the identity of `inp.breakpoints` plus `rounds_hash`. `extra` is free for Task 4.
  - `merge_rounds.CACHE`, a module-level `BaseCache(8)`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_merge_rounds.py
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
    assert a == b and len(a['rows']) == 3


def test_hatch_routes_read_the_base():
    c = _client()
    rounds = [R(LOOSE)]
    s = c.post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'cluster0', 'rounds': rounds}).json()
    assert len(s['rows']) < 3
    m = c.post('/api/integration/hatch-merge', json={'chrom': 'chr1', 'source': 'cluster0', 'rules': LOOSE, 'rounds': rounds}).json()
    assert len(m['original']) == len(s['rows'])


def test_hatch_routes_reject_invalid_round():
    c = _client()
    r = c.post('/api/integration/hatch-scores', json={'chrom': 'chr1', 'source': 'cluster0', 'rounds': [R([])]})
    assert r.status_code == 422 and r.json()['detail'].startswith('round 1: ')
    r = c.post('/api/integration/hatch-merge', json={'chrom': 'chr1', 'source': 'cluster0', 'rules': LOOSE, 'rounds': [R([])]})
    assert r.status_code == 422 and r.json()['detail'].startswith('round 1: ')
```

Before relying on `len(a['rows']) == 3`, check what the hatch-scores response's `rows` counts for 2 boundaries, using `tests/test_hatch_merge.py::test_hatch_scores_route_returns_rows_and_bp_arrays`. If rows are per internal segment rather than per segment, change both assertions to the matching numbers. Keep the relation the same: fewer rows with the base than without, and the hatch-merge `original` count equal to the base's segment count. Record the adjustment in the report.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest tests/test_merge_rounds.py -q`
Expected: FAIL with `ImportError: cannot import name 'merge_rounds'`.

- [ ] **Step 3: Add `boundaries_from_bp` to `server/cn_recall.py`**

Place it directly after `seg_labels`:

```python
def boundaries_from_bp(bp: pd.DataFrame, sources, chroms) -> Dict[str, Dict[str, List[int]]]:
    """Chromosome-local boundaries per source/chromosome from a breakpoints table."""
    out = {s: {c: [] for c in chroms} for s in sources}
    for r in bp.itertuples(index=False):
        if r.source in out and r.chromosome in out[r.source]:
            out[r.source][r.chromosome].append(int(r.absolute_bin))
    return {s: {c: sorted(v) for c, v in d.items()} for s, d in out.items()}
```

- [ ] **Step 4: Write `server/merge_rounds.py`**

```python
# server/merge_rounds.py
"""Merge rounds: the working segmentation is the original accepted breakpoints
with each round's genome-wide rule merge applied in turn. Every round merges on
the original pb_x; only the boundaries carry over.
Spec: docs/superpowers/specs/2026-10-27-merge-rounds-design.md §2."""
from __future__ import annotations

import hashlib
import json
import threading
from collections import OrderedDict
from dataclasses import dataclass, field, replace
from typing import List

import numpy as np
import pandas as pd

from . import cn_recall as cr
from . import merge_rules as mr


def normalize_rounds(rounds) -> List[dict]:
    if rounds is None:
        return []
    if not isinstance(rounds, list):
        raise ValueError('rounds must be a list')
    out = []
    for i, r in enumerate(rounds, 1):
        if not isinstance(r, dict) or not isinstance(r.get('rules'), list):
            raise ValueError(f'round {i}: rules must be a list')
        try:
            mr.validate_rules(r['rules'])
        except ValueError as e:
            raise ValueError(f'round {i}: {e}') from e
        out.append({'rules': r['rules'], 'allow_gap_crossing': bool(r.get('allow_gap_crossing', False))})
    return out


def rounds_hash(rounds) -> str:
    return hashlib.sha256(json.dumps(rounds, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


@dataclass
class BaseInputs:
    """What merge_genome needs; RecallInputs also satisfies this shape."""
    var_seq: np.ndarray
    var_start: np.ndarray
    var_end: np.ndarray
    pb_x: np.ndarray
    breakpoints: pd.DataFrame
    sources: List[str] = field(default_factory=list)

    def chromosomes(self) -> List[str]:
        return list(dict.fromkeys(self.var_seq.tolist()))


def fold(inp, rounds):
    """Apply each (normalized) round in turn; returns (breakpoints, per-round summaries)."""
    bp, per = inp.breakpoints, []
    for i, r in enumerate(rounds, 1):
        try:
            m = cr.merge_genome(replace(inp, breakpoints=bp), r['rules'], r['allow_gap_crossing'])
        except cr.RecallError as e:
            raise cr.RecallError(f'round {i}: {e}') from e
        bp = m['merged_breakpoints'].reset_index(drop=True)
        per.append({'pairs': m['pairs'], 'merges': sum(p['merges'] for p in m['pairs'])})
    return bp, per


class BaseCache:
    def __init__(self, size: int = 8):
        self.size = size
        self._d: OrderedDict = OrderedDict()
        self.lock = threading.Lock()

    def get(self, inp, rounds) -> dict:
        h = rounds_hash(rounds)
        key = (id(inp.breakpoints), h)
        with self.lock:
            e = self._d.get(key)
            if e is not None and e['source'] is inp.breakpoints:
                self._d.move_to_end(key)
                return e
        bp, per = fold(inp, rounds)
        e = {'source': inp.breakpoints, 'hash': h, 'breakpoints': bp, 'rounds': per, 'extra': {}}
        with self.lock:
            self._d[key] = e
            self._d.move_to_end(key)
            while len(self._d) > self.size:
                self._d.popitem(last=False)
        return e


CACHE = BaseCache(8)
```

- [ ] **Step 5: Teach the hatch routes about `rounds`, in `server/hatch_merge_api.py`**

Add the import `from . import merge_rounds as mrd` and `from . import cn_recall as cr` next to the existing imports. Add `rounds: List[Dict[str, Any]] = []` to both `HatchScoresRequest` and `HatchMergeRequest`. Replace `_resolve` with:

```python
def _resolve(ig: IntegrationStore, chrom: str, source: str, rounds=None):
    if ig.arrays is None:
        raise HTTPException(503, 'Pseudobulk index is not ready')
    if source not in ig.members:
        raise HTTPException(422, f'Unknown source {source!r}')
    idx = np.flatnonzero(ig.var_seq == chrom)
    if not len(idx):
        raise HTTPException(422, f'Unknown chromosome {chrom!r}')
    if not np.array_equal(idx, np.arange(idx[0], idx[-1] + 1)):
        raise HTTPException(422, 'Chromosome bins are not contiguous')
    si = list(ig.members).index(source)  # ig.members is a dict {source: cell-index-array}; list() gives its keys in insertion order
    profile = ig.arrays['pb_x'][si, idx]
    bp = ig.tables['breakpoints']
    norm = mrd.normalize_rounds(rounds)   # ValueError -> 422 in the routes
    if norm:
        base = mrd.BaseInputs(var_seq=np.asarray(ig.var_seq), var_start=ig.var_start, var_end=ig.var_end,
                              pb_x=ig.arrays['pb_x'], breakpoints=ig.tables['breakpoints'], sources=list(ig.members))
        try:
            bp = mrd.CACHE.get(base, norm)['breakpoints']
        except cr.RecallError as e:
            raise HTTPException(500, str(e)) from e
    boundaries = bp.loc[(bp.source == source) & (bp.chromosome == chrom), 'absolute_bin'].tolist()
    return profile, boundaries, idx
```

In both routes, pass `req.rounds` as the fourth argument: `_resolve(ig, req.chrom, req.source, req.rounds)`. Both routes already turn `ValueError` into 422. Change nothing else in the responses.

- [ ] **Step 6: Run the tests to verify they pass**

Run: `python3 -m pytest tests/test_merge_rounds.py tests/test_hatch_merge.py tests/test_merge_rules.py tests/test_cn_recall.py -q`
Expected: all PASS.

- [ ] **Step 7: Commit**

```bash
git add server/merge_rounds.py server/cn_recall.py server/hatch_merge_api.py tests/test_merge_rounds.py
git commit -m "feat: merge rounds folded on the server; hatch routes accept rounds

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Re-call pipeline: calls-only job, waiting for the base, noise in the diff

**Files:**
- Modify: `server/cn_recall.py`
- Test: `tests/test_cn_recall.py` (append tests only; existing tests unchanged)

**Interfaces:**
- Consumes:
  - `call_noise.cluster_noise`, `cells_noise`, `combine` (Task 1)
  - `boundaries_from_bp` (Task 2)
- Produces:
  - `RecallInputs.before_percell_wait: Optional[Callable[[], tuple]] = None`, which returns `(u8, best_s)`
  - `call_all_cells(inp, segs, executor, update, chunk=100, max_inflight=8) -> (u8, best)`
  - `run_calls(inp, executor, update, chunk=100, max_inflight=8) -> {'percell_u8', 'percell_best_s', 'noise_cells'}`
  - `run_recall(...)` now also returns `'percell_best_s'` and `'noise_cells'`, and publishes `diff['noise']`:
    - after stage 2: cluster level only;
    - after stage 3: cluster level and cells.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_cn_recall.py`)

```python
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest tests/test_cn_recall.py -q -k "noise or run_calls or waits or wait_failure"`
Expected: FAIL (`KeyError: 'noise'` / `AttributeError: run_calls`).

- [ ] **Step 3: Implement in `server/cn_recall.py`**

1. Add the import `from . import call_noise as cn` next to the other imports.

2. Add a field to `RecallInputs` after `sources`:

```python
    before_percell_wait: Optional[Callable[[], tuple]] = None   # base per-cell (u8, best_s) still being computed
```

3. Move the per-cell loop out of `run_recall` into a function that returns the matrix. It is the same code as today's `collect` / submit loop, with the result arrays as locals:

```python
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
```

4. Changes in `run_recall`:
   - **After building `rows` and `diff` (stage 2):** add `cluster_nz = cn.cluster_noise(pd.DataFrame(rows))` and put `'noise': cn.combine(cluster_nz, None)` into `diff` before the `update(merge=…, cluster_segments=rows, diff=diff)` call.
   - **Replace the stage-3 preamble.** The current lines from `if inp.before_percell_u8 is None:` through the end of the `while pending:` loop become:

```python
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
```

   - **In the per-cell diff loop that follows:** read `before_u8[i]` and `before_best[i]` instead of `inp.before_percell_u8[i]` and `inp.before_percell_best_s[i]`.
   - **At the end:**

```python
    cells_nz = cn.cells_noise(out_u8, inp.cell_cluster, segs, inp.var_seq, inp.sources)
    diff = {**diff, 'cells': cells_diff, 'transitions': {**diff['transitions'], 'cells': transitions_list(cacc)},
            'noise': cn.combine(cluster_nz, cells_nz)}
    update(diff=diff, percell_ready=True)
    return {'percell_u8': out_u8, 'percell_best_s': best, 'noise_cells': cells_nz}
```

   `segs` in `run_recall` stays built from `merge['boundaries']`, as today.

- [ ] **Step 4: Run the whole re-call suite**

Run: `python3 -m pytest tests/test_cn_recall.py -q`
Expected: all PASS, old and new.

- [ ] **Step 5: Commit**

```bash
git add server/cn_recall.py tests/test_cn_recall.py
git commit -m "feat: calls-only per-cell job, base wait, and noise in the re-call diff

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Re-call API: rounds, pinning, calls jobs and `/api/integration/base`

**Files:**
- Modify: `server/cn_recall_api.py`
- Test: `tests/test_merge_rounds_api.py`

**Interfaces:**
- Consumes:
  - `merge_rounds.normalize_rounds`, `rounds_hash`, `CACHE` (Task 2)
  - `cn_recall.run_recall`, `run_calls`, `boundaries_from_bp`, `seg_labels`, `cluster_rows`, `pseudobulk` (Task 3)
  - `call_noise.*` (Task 1)
- Produces these routes:

| Route | Body | Response |
|---|---|---|
| `POST /api/integration/recall` | `{rules, allow_gap_crossing, rounds=[]}` | |
| `POST /api/integration/recall/{id}/pin` | `{rounds}` | `{pinned: id}`; 404 / 409 / 422 |
| `POST /api/integration/base` | `{rounds}` | `{hash, n_rounds, breakpoints, rounds, cluster_segments, noise, noise_original, percell: {status, job}}` |

  `percell.status` is one of `original`, `ready`, `running`, `missing`. `/base` never returns `missing` for a non-empty recipe, because it starts the recompute.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_merge_rounds_api.py
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest tests/test_merge_rounds_api.py -q`
Expected: FAIL (404 for `/api/integration/base`, and similar).

- [ ] **Step 3: Implement in `server/cn_recall_api.py`**

1. Imports and constants:

```python
from dataclasses import replace
from typing import Optional

from . import call_noise as cn
from . import merge_rounds as mrd

IMPLEMENTATION_SHA256 = hashlib.sha256(b''.join((_HERE / f).read_bytes() for f in ('merge_rules.py', 'hatch_merge.py', 'cn_recall.py', 'merge_rounds.py', 'call_noise.py'))).hexdigest()
SCOPE_ROUNDS = ('Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on a rule-merged segmentation '
                '({n} earlier rounds); originals unchanged.')
```

2. Request models:
   - Add `rounds: List[Dict[str, Any]] = []` to `RecallRequest`.
   - Add these two models:

```python
class BaseRequest(BaseModel):
    rounds: List[Dict[str, Any]] = []


class PinRequest(BaseModel):
    rounds: List[Dict[str, Any]]
```

3. Hash. Keep the old formula exactly when there are no rounds:

```python
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
```

4. Replace `RecallManager` with the version below. It keeps the existing behaviour and adds `pinned`, calls jobs, `_result_key`, `_best`, `_noise_cells` and a done event:

```python
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
        """('ready' | 'running' | 'missing', job id) for the per-cell calls of the base with rounds_hash `key`."""
        with self.lock:
            p = self.jobs.get(self.pinned) if self.pinned else None
            if p and p['_result_key'] == key and p['status'] != 'error':
                return ('ready' if p['status'] == 'done' else 'running'), p['id']
            for j in self.jobs.values():
                if j['_result_key'] == key and j['status'] == 'done' and j['_u8'] is not None:
                    return 'ready', j['id']
            return 'missing', None

    def wait_percell(self, jid):
        job = self.jobs.get(jid)
        if job is None:
            raise cr.RecallError(f'percell base: job {jid} is gone')
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
```

5. In `register_cn_recall`, keep the three existing routes, and add the helpers and routes below inside it, after `manager = RecallManager(keep=keep)`. `start_recall` changes as shown.

```python
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
        status, jid = manager.base_percell(e['hash'])
        if status == 'missing':
            jid = manager.start_calls(replace(inp, breakpoints=e['breakpoints']), rounds, identity, executor)['id']
            status = 'running'
        elif status == 'ready':
            manager.pin(jid)
        return e, status, jid

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
            e, status, jid = base_state(inp, rounds, identity, executor)
            u8 = best = None
            if status == 'ready':
                with manager.lock:
                    job = manager.jobs[jid]
                    u8, best = job['_u8'], job['_best']
            inp = replace(inp, breakpoints=e['breakpoints'],
                          before_cluster_segments=pd.DataFrame(e['extra']['cluster_segments']),
                          before_percell_u8=u8,
                          before_percell_best_s=best if best is not None else np.full(len(inp.cell_ids), np.nan),
                          before_percell_wait=(lambda j=jid: manager.wait_percell(j)) if status == 'running' else None)
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
        e, status, jid = base_state(inp, rounds, identity, executor)
        cells = manager.jobs[jid]['_noise_cells'] if status == 'ready' else None
        return {'hash': e['hash'], 'n_rounds': len(rounds), 'breakpoints': _records(e['breakpoints']), 'rounds': e['rounds'],
                'cluster_segments': e['extra']['cluster_segments'],
                'noise': cn.combine(e['extra']['cluster_noise'], cells), 'noise_original': orig,
                'percell': {'status': status, 'job': jid}}
```

   - **Cluster segments on a base:** `test_adopt_then_base_uses_the_pinned_recall` asserts that `/base` returns the adopted re-call's `cluster_segments`. Both are computed from the float64 member mean (`cr.pseudobulk`) with `cluster_rows`, so they match. If a float difference shows up, it means `pb64` is not `cr.pseudobulk(inp)`. Fix the cause; do not loosen the test.
   - **Numpy values in JSON:** `_records` must return JSON-safe values, as `run_recall`'s `bp_records` does today. If FastAPI fails to encode numpy scalars, convert with `.item()` inside `_records`.

- [ ] **Step 4: Run the server test suites**

Run: `python3 -m pytest tests/ -q`
Expected: all PASS, including the unchanged `tests/test_cn_recall.py` API tests (key sets, settings, scope sentence with no rounds).

- [ ] **Step 5: Real-data smoke test on a private server** (not 8766)

```bash
python3 run.py --port 8799 &   # note the PID; stop only this one when finished
sleep 20
R='{"rules":[{"id":"rule0","type":"flank_threshold","name":"Rule 0","criteria":["delta_mean","srd_phi"],"combine":"any","tau":0.5,"srdRef":8,"maxBins":100,"mergeSignMismatch":false,"mergeSignsOnly":false,"hatch":{"pattern":"x","color":"#b91c1c"}}],"allow_gap_crossing":false}'
curl -s -XPOST localhost:8799/api/integration/base -H 'content-type: application/json' -d "{\"rounds\":[$R]}" | python3 -c "import json,sys; b=json.load(sys.stdin); print(b['n_rounds'], len(b['breakpoints']), b['percell'], b['noise']['cluster2']['genome']['cluster'], b['noise_original']['cluster2']['genome']['cluster'])"
```

Expected:
- `1`, then fewer breakpoints than the original (946 rows minus 358 = 588), then `{'status': 'running', …}`;
- cluster2 genome zig-zags 7 against original 33.

Poll the job until it is done (about 60 s), and query `/base` again: the status is `ready`, and cells noise is present, with about 9 zig-zags per cell. Put the numbers in the report, then stop the 8799 server.

- [ ] **Step 6: Commit**

```bash
git add server/cn_recall_api.py tests/test_merge_rounds_api.py
git commit -m "feat: re-calls on a base, pinning, calls-only recompute and /api/integration/base

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Web helpers, API and store (`web/src/mergeRounds.ts`, `cnRecall.ts`, `hatchMerge.ts`, `integrationStore.ts`)

**Files:**
- Create: `web/src/mergeRounds.ts`
- Modify: `web/src/cnRecall.ts`, `web/src/hatchMerge.ts`, `web/src/integrationStore.ts`
- Test: `tests/merge_rounds_helpers.ts`; add `tests/_bundle/merge_rounds_helpers.mjs` to `.gitignore`, next to the other bundles.

**Interfaces:**
- Consumes the Task 4 routes and response shapes.
- Produces:
  - **From `mergeRounds.ts`:**
    - types: `WireRound`, `Round`, `Noise`, `ClusterNoise`, `CellsNoise`, `BaseResponse`, `NoiseColumn`;
    - helpers: `wireRounds`, `sameRounds`, `roundLabel`, `roundStripText`, `baseHeading`, `karyoModes`, `noiseRows`, `noiseCellsN`, `noiseTsv`, `roundsHeader`, `roundsFromSession`, `NOISE_DEFINITIONS`;
    - API calls: `fetchBase`, `pinRecallApi`.
  - **New store fields and actions:**
    - fields: `rounds`, `base`, `baseSig`, `baseError`, `roundsNotice`, `baseKaryo`, `baseCalls`;
    - actions: `setRounds`, `refreshBase`, `adoptRecall`, `undoRound`, `resetRounds`.
  - **Type change:** `KaryoMode` gains `'original'`.

- [ ] **Step 1: Write the failing TS helper suite**

```ts
// tests/merge_rounds_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/merge_rounds_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rounds_helpers.mjs
//           node tests/_bundle/merge_rounds_helpers.mjs
import { wireRounds, sameRounds, roundLabel, roundStripText, baseHeading, karyoModes, noiseRows, noiseCellsN, noiseTsv, roundsHeader, roundsFromSession, NOISE_DEFINITIONS, type Noise, type Round } from '../web/src/mergeRounds'
import { recallHeader, type RecallStatus } from '../web/src/cnRecall'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}

const r0 = { id: 'rule0', type: 'flank_threshold', tau: 0.5, srdRef: 8, combine: 'any' }
const noise: Noise = { cluster2: {
  chr1: { cluster: { segments: 10, zigzags: 2, zigzag_median_bins: 12.5, gain_loss_adjacent: 1 }, cells: null },
  genome: { cluster: { segments: 144, zigzags: 33, zigzag_median_bins: 40, gain_loss_adjacent: 6 }, cells: { n_cells: 612, segments_mean: 160.5, zigzags_mean: 31, zigzag_median_bins: null, gain_loss_adjacent_mean: 13.42 } } } }
const round: Round = { rules: [r0, { id: 'u1', enabled: true }, { id: 'u2', enabled: false }], allow_gap_crossing: false, label: 'x', merges: 358, noise }

eq(wireRounds([round]), [{ rules: round.rules, allow_gap_crossing: false }], 'wire drops display fields')
eq(sameRounds([round], [{ allow_gap_crossing: false, rules: round.rules } as any]), true, 'same recipe regardless of extra fields and key order')
eq(sameRounds([round], []), false, 'different recipe')
eq(roundLabel([r0]), 'τ 0.5 · ref 8 · any', 'label: Rule 0 only')
eq(roundLabel(round.rules as any[]), 'τ 0.5 · ref 8 · any + 1 user rule', 'label: enabled user rules counted')
eq(roundStripText([]), 'Original', 'strip: original')
eq(roundStripText([{ ...round, label: 'τ 0.5 · ref 8 · any' }]), 'Original → Round 1 (τ 0.5 · ref 8 · any, 358 merges)', 'strip: one round')
eq(baseHeading(0), 'Base: original', 'heading 0')
eq(baseHeading(1), 'Base: round 1 (1 round of merges on the original)', 'heading 1')
eq(baseHeading(2), 'Base: round 2 (2 rounds of merges on the original)', 'heading 2')
eq(karyoModes(0, false), [{ mode: 'before', label: 'Before', enabled: true }, { mode: 'after', label: 'After', enabled: false }, { mode: 'changed', label: 'Changed', enabled: false }], 'modes: no rounds, no re-call')
eq(karyoModes(1, true).map((m) => [m.mode, m.enabled]), [['original', true], ['before', true], ['after', true], ['changed', true]], 'modes: rounds and a re-call')

const cols = [{ name: 'Original', noise }, { name: 'After', noise: { cluster2: { genome: { cluster: { segments: 97, zigzags: 7, zigzag_median_bins: 107, gain_loss_adjacent: 1 }, cells: null } } }, cellsPending: true }]
const rows = noiseRows(cols, 'cluster2', 'genome')
eq(rows.map((r) => [r.level, r.measure, ...r.values]), [
  ['cluster', 'segments', '144', '97'], ['cluster', 'zig-zags', '33', '7'], ['cluster', 'zig-zag median length (bins)', '40', '107'], ['cluster', 'gain next to loss', '6', '1'],
  ['cells', 'segments per cell', '160.5', 're-calling…'], ['cells', 'zig-zags per cell', '31.0', 're-calling…'], ['cells', 'zig-zag median length (bins)', '—', 're-calling…'], ['cells', 'gain next to loss per cell', '13.4', 're-calling…'],
], 'noise rows: formatting, missing median, pending cells')
eq(noiseRows(cols, 'cluster2', 'chr1')[2].values, ['12.5', '—'], 'noise rows: half-bin median and absent key')
eq(noiseRows([{ name: 'Original', noise: null }], 'cluster2', 'genome')[0].values, ['—'], 'noise rows: no data')
eq(noiseCellsN(cols, 'cluster2', 'genome'), 612, 'cells n from first column with cells')

const hdr = roundsHeader([round])
eq(hdr, [`round 1 rules: ${JSON.stringify(round.rules)} · allow_gap_crossing: false`], 'rounds header')
const t = noiseTsv([{ name: 'original', noise }, { name: 'round1', noise }], ['cluster2'], hdr).split('\n')
eq(t[0], hdr.map((h) => `# ${h}`)[0], 'tsv header line')
eq(t[1], 'source\tchromosome\tlevel\tmeasure\toriginal\tround1', 'tsv columns')
eq(t[2], 'cluster2\tchr1\tcluster\tsegments\t10\t10', 'tsv first row')
eq(t.filter((l) => l.startsWith('cluster2\tgenome\tcells\tzig-zag median')).length, 1, 'tsv cells rows present')
eq(t.find((l) => l.startsWith('cluster2\tchr1\tcells\tsegments')), 'cluster2\tchr1\tcells\tsegments per cell\t\t', 'tsv empty cells values')

eq(roundsFromSession({}), [], 'session without rounds')
eq(roundsFromSession({ rounds: [round] }), [round], 'session with rounds')
eq(roundsFromSession({ rounds: 'x' }), [], 'session with garbage rounds')
eq(NOISE_DEFINITIONS.length, 4, 'four definitions')
eq(NOISE_DEFINITIONS[1].startsWith('Zig-zag: an internal segment whose state differs from both neighbours'), true, 'zig-zag definition text')

const st = { settings: { rules: [r0], allow_gap_crossing: false, rounds: [{ rules: [r0], allow_gap_crossing: false }] }, provenance: {} } as unknown as RecallStatus
eq(recallHeader(st)[0], `base rounds: ${JSON.stringify([{ rules: [r0], allow_gap_crossing: false }])}`, 'recall header lists base rounds first')
eq(recallHeader({ ...st, settings: { rules: [r0], allow_gap_crossing: false } } as RecallStatus)[0].startsWith('rules: '), true, 'no base rounds: header unchanged')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) } else console.log('all merge_rounds helper checks passed')
```

- [ ] **Step 2: Run it to verify it fails**

Run: `web/node_modules/.bin/esbuild tests/merge_rounds_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/merge_rounds_helpers.mjs && node tests/_bundle/merge_rounds_helpers.mjs`
Expected: build error, `Could not resolve "../web/src/mergeRounds"`.

- [ ] **Step 3: Write `web/src/mergeRounds.ts`**

```ts
// Merge rounds: the working segmentation recipe, the call-noise table and their exports.
// Spec: docs/superpowers/specs/2026-10-27-merge-rounds-design.md
import type { ClusterBreakpoint, ClusterSegment } from './types'
import { tsv, type KaryoMode, type RecallPair } from './cnRecall'

export interface WireRound { rules: unknown[]; allow_gap_crossing: boolean }
export interface ClusterNoise { segments: number; zigzags: number; zigzag_median_bins: number | null; gain_loss_adjacent: number }
export interface CellsNoise { n_cells: number; segments_mean: number; zigzags_mean: number; zigzag_median_bins: number | null; gain_loss_adjacent_mean: number }
export type Noise = Record<string, Record<string, { cluster: ClusterNoise | null; cells: CellsNoise | null }>>
export interface Round extends WireRound { label: string; merges: number; noise: Noise | null }
export interface BaseResponse {
  hash: string; n_rounds: number; breakpoints: ClusterBreakpoint[]; rounds: { pairs: RecallPair[]; merges: number }[]
  cluster_segments: ClusterSegment[]; noise: Noise; noise_original: Noise
  percell: { status: 'original' | 'ready' | 'running' | 'missing'; job: string | null }
}
export interface NoiseColumn { name: string; noise: Noise | null; cellsPending?: boolean }

export const wireRounds = (rs: WireRound[]): WireRound[] => rs.map((r) => ({ rules: r.rules, allow_gap_crossing: r.allow_gap_crossing }))
export const sameRounds = (a: WireRound[], b: WireRound[]) => JSON.stringify(wireRounds(a)) === JSON.stringify(wireRounds(b))

export function roundLabel(rules: any[]): string {
  const r0 = rules[0] ?? {}
  const users = rules.slice(1).filter((r) => r && r.enabled).length
  return `τ ${r0.tau} · ref ${r0.srdRef} · ${r0.combine}${users ? ` + ${users} user rule${users > 1 ? 's' : ''}` : ''}`
}
export const roundStripText = (rounds: Round[]) => ['Original', ...rounds.map((r, i) => `Round ${i + 1} (${r.label}, ${r.merges} merges)`)].join(' → ')
export const baseHeading = (n: number) => (n ? `Base: round ${n} (${n} round${n > 1 ? 's' : ''} of merges on the original)` : 'Base: original')

const MODE_LABELS: Record<KaryoMode, string> = { original: 'Original', before: 'Before', after: 'After', changed: 'Changed' }
export function karyoModes(nRounds: number, hasAfter: boolean): { mode: KaryoMode; label: string; enabled: boolean }[] {
  const modes: KaryoMode[] = nRounds ? ['original', 'before', 'after', 'changed'] : ['before', 'after', 'changed']
  return modes.map((mode) => ({ mode, label: MODE_LABELS[mode], enabled: mode === 'original' || mode === 'before' || hasAfter }))
}

const MEASURES: { level: 'cluster' | 'cells'; key: string; label: string; kind: 'count' | 'median' | 'mean' }[] = [
  { level: 'cluster', key: 'segments', label: 'segments', kind: 'count' },
  { level: 'cluster', key: 'zigzags', label: 'zig-zags', kind: 'count' },
  { level: 'cluster', key: 'zigzag_median_bins', label: 'zig-zag median length (bins)', kind: 'median' },
  { level: 'cluster', key: 'gain_loss_adjacent', label: 'gain next to loss', kind: 'count' },
  { level: 'cells', key: 'segments_mean', label: 'segments per cell', kind: 'mean' },
  { level: 'cells', key: 'zigzags_mean', label: 'zig-zags per cell', kind: 'mean' },
  { level: 'cells', key: 'zigzag_median_bins', label: 'zig-zag median length (bins)', kind: 'median' },
  { level: 'cells', key: 'gain_loss_adjacent_mean', label: 'gain next to loss per cell', kind: 'mean' },
]
const fmtValue = (v: number | null | undefined, kind: string) =>
  v === null || v === undefined ? '—' : kind === 'count' ? String(v) : kind === 'median' ? (Number.isInteger(v) ? String(v) : v.toFixed(1)) : v.toFixed(1)
const entry = (n: Noise | null, source: string, key: string) => n?.[source]?.[key]

export function noiseRows(cols: NoiseColumn[], source: string, key: string) {
  return MEASURES.map((m) => ({
    level: m.level, measure: m.label,
    values: cols.map((c) => {
      const e = entry(c.noise, source, key)
      const lv = m.level === 'cluster' ? e?.cluster : e?.cells
      if (!lv && m.level === 'cells' && c.cellsPending) return 're-calling…'
      return fmtValue(lv ? (lv as any)[m.key] : null, m.kind)
    }),
  }))
}
export function noiseCellsN(cols: NoiseColumn[], source: string, key: string): number | null {
  for (const c of cols) { const n = entry(c.noise, source, key)?.cells?.n_cells; if (n !== undefined) return n }
  return null
}
export const roundsHeader = (rounds: WireRound[]) => rounds.map((r, i) => `round ${i + 1} rules: ${JSON.stringify(r.rules)} · allow_gap_crossing: ${r.allow_gap_crossing}`)
export function noiseTsv(cols: NoiseColumn[], sources: string[], header: string[]): string {
  const rows: Record<string, unknown>[] = []
  for (const source of sources) {
    const keys: string[] = []
    for (const c of cols) for (const k of Object.keys(c.noise?.[source] ?? {})) if (!keys.includes(k)) keys.push(k)
    for (const key of keys) for (const m of MEASURES) {
      const r: Record<string, unknown> = { source, chromosome: key, level: m.level, measure: m.label }
      for (const c of cols) { const e = entry(c.noise, source, key); const lv = m.level === 'cluster' ? e?.cluster : e?.cells; r[c.name] = lv ? (lv as any)[m.key] : null }
      rows.push(r)
    }
  }
  return tsv(header, ['source', 'chromosome', 'level', 'measure', ...cols.map((c) => c.name)], rows)
}
export const roundsFromSession = (p: any): Round[] => (Array.isArray(p?.rounds) ? p.rounds : [])

export const NOISE_DEFINITIONS = [
  'Calls are read per chromosome as the sequence of segment states. Segments with no call are skipped. Gain = state > 1 (weak gain included), loss = state < 1.',
  'Zig-zag: an internal segment whose state differs from both neighbours, while the two neighbours have the same state as each other.',
  'Gain next to loss: a boundary with a gain on one side and a loss on the other. This can be real (for example an unbalanced rearrangement), so it is counted but never merged automatically.',
  "Cell values: the mean per cell over all cells of the cluster. The zig-zag median length is taken over all of those cells' zig-zags. For all, the cells are every cell, each read under its own cluster's segmentation.",
]

async function postJson<T>(url: string, body: unknown): Promise<T> {
  const r = await fetch(url, { method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify(body) })
  if (!r.ok) { let msg = `${r.status}`; try { const b = await r.json(); msg = typeof b.detail === 'string' ? b.detail : JSON.stringify(b.detail ?? b) } catch { /* keep status */ } throw new Error(msg) }
  return r.json() as Promise<T>
}
export const fetchBase = (rounds: WireRound[]) => postJson<BaseResponse>('/api/integration/base', { rounds })
export const pinRecallApi = (id: string, rounds: WireRound[]) => postJson<{ pinned: string }>(`/api/integration/recall/${encodeURIComponent(id)}/pin`, { rounds })
```

- [ ] **Step 4: Update `web/src/cnRecall.ts`**

Change these lines:

```ts
export type KaryoMode = 'original' | 'before' | 'after' | 'changed'
export interface RecallBody { rules: unknown[]; allow_gap_crossing: boolean; rounds?: { rules: unknown[]; allow_gap_crossing: boolean }[] }
```

In `RecallStatus`:
- `settings` becomes `{ rules: unknown[]; allow_gap_crossing: boolean; rounds?: { rules: unknown[]; allow_gap_crossing: boolean }[] }`.
- `diff` gains an optional `noise?: import('./mergeRounds').Noise`. Use `import type { Noise } from './mergeRounds'` at the top: type-only, no runtime cycle.

`recallHeader` lists the base rounds first, when there are any:

```ts
export const recallHeader = (st: RecallStatus) => [
  ...(st.settings.rounds?.length ? [`base rounds: ${JSON.stringify(st.settings.rounds)}`] : []),
  `rules: ${JSON.stringify(st.settings.rules)}`, `allow_gap_crossing: ${st.settings.allow_gap_crossing}`, `provenance: ${JSON.stringify(st.provenance)}`]
```

In `RecallImpact.tsx`, change `const LABELS: Record<KaryoMode, string>` to include `original: 'Original'`, so it still type-checks. Task 6 replaces this switch anyway.

- [ ] **Step 5: Update `web/src/hatchMerge.ts`**

```ts
export interface MergeRequest { chrom: string; source: string; rules: AnyRule[]; allow_gap_crossing: boolean; rounds?: { rules: unknown[]; allow_gap_crossing: boolean }[] }
export const fetchHatchScores = (chrom: string, source: string, rounds: { rules: unknown[]; allow_gap_crossing: boolean }[] = []) =>
  postJson<HatchScoresResult>('/api/integration/hatch-scores', { chrom, source, rounds })
```

- [ ] **Step 6: Update `web/src/integrationStore.ts`**

1. Imports: add `import { fetchBase, pinRecallApi, roundLabel, roundsFromSession, sameRounds, wireRounds, type BaseResponse, type Round } from './mergeRounds'`.

2. Add to `IntegrationState`:

```ts
  rounds: Round[]; base: BaseResponse | null; baseSig: string | null; baseError: string | null; roundsNotice: string | null
  baseKaryo: Uint8Array | null; baseCalls: RecallStatus | null
  setRounds: (rounds: Round[]) => void
  refreshBase: () => Promise<void>
  adoptRecall: () => Promise<void>
  undoRound: () => void
  resetRounds: () => void
```

3. Initial values: `rounds: [], base: null, baseSig: null, baseError: null, roundsNotice: null, baseKaryo: null, baseCalls: null,`.

4. Actions. Add them next to `startRecall`. `setRounds` clears everything derived from the base. The view's effect (Task 6) then calls `refreshBase`; `setRounds` never fetches itself.

```ts
  setRounds: (rounds) => set({ rounds, base: null, baseSig: null, baseError: null, baseKaryo: null, baseCalls: null, hatchScores: {}, hatchScoresError: null,
    recall: null, recallKaryo: null, recallSignature: null, recallError: null, karyoMode: 'before' }),
  undoRound: () => { set({ roundsNotice: null }); get().setRounds(get().rounds.slice(0, -1)) },
  resetRounds: () => { set({ roundsNotice: null }); get().setRounds([]) },
  refreshBase: async () => {
    const wire = wireRounds(get().rounds), sig = JSON.stringify(wire)
    set({ baseSig: sig, baseError: null })
    let b: BaseResponse
    try { b = await fetchBase(wire) } catch (e: any) {
      if (get().baseSig !== sig) return
      const msg = String(e.message || e)
      if (wire.length && /^round \d+: /.test(msg)) { get().setRounds([]); set({ roundsNotice: `Saved merge rounds could not be applied: ${msg}` }); return }
      set({ baseError: `segmentation base: ${msg}` }); return
    }
    if (get().baseSig !== sig) return
    set({ base: b })
    const job = b.percell.job
    if (b.percell.status === 'ready' && job && !get().baseKaryo) {
      try { const k = await fetchRecallPercell(job); if (get().baseSig === sig) set({ baseKaryo: k }) }
      catch (e: any) { if (get().baseSig === sig) set({ baseError: `base per-cell: ${e.message || e}` }) }
    } else if (b.percell.status === 'running' && job) {
      const poll = async () => {
        if (get().baseSig !== sig) return
        let st: RecallStatus
        try { st = await fetchRecall(job) } catch (e: any) { if (get().baseSig === sig) set({ baseError: `base per-cell ${job}: ${e.message || e}` }); return }
        if (get().baseSig !== sig) return
        set({ baseCalls: st })
        if (st.status === 'running') { setTimeout(poll, 1000); return }
        if (st.status === 'done') await get().refreshBase()   // now 'ready': loads the matrix and the cells noise
      }
      await poll()
    }
  },
  adoptRecall: async () => {
    const st = get().recall, cur = get().rounds
    if (!st || st.status !== 'done' || !sameRounds(st.settings.rounds ?? [], cur)) return
    const round: Round = { rules: st.settings.rules, allow_gap_crossing: st.settings.allow_gap_crossing, label: roundLabel(st.settings.rules as any[]),
      merges: (st.merge?.pairs ?? []).reduce((a, p) => a + p.merges, 0), noise: st.diff?.noise ?? null }
    const next = [...cur, round]
    try { await pinRecallApi(st.id, wireRounds(next)) } catch (e: any) { set({ recallError: `adopt: ${e.message || e}` }); return }
    const k = get().recallKaryo
    set({ roundsNotice: null })
    get().setRounds(next)
    if (k) set({ baseKaryo: k })   // the adopted re-call's per-cell calls are the new base's
  },
```

5. `ensureHatchScores` sends the rounds, and keeps a response only if the recipe has not changed since the request:

```ts
  ensureHatchScores: async (source, chrom) => {
    const k = `${source}|${chrom}`
    if (get().hatchScores[k]) return
    const wire = wireRounds(get().rounds), sig = JSON.stringify(wire)
    try { const r = await fetchHatchScores(chrom, source, wire); if (JSON.stringify(wireRounds(get().rounds)) === sig) set({ hatchScores: { ...get().hatchScores, [k]: r }, hatchScoresError: null }) }
    catch (e: any) { if (JSON.stringify(wireRounds(get().rounds)) === sig) set({ hatchScoresError: `hatch-scores: ${e.message || e}` }) }
  },
```

6. Sessions:
   - `serialize` adds `rounds: s.rounds`.
   - `restore` adds `rounds: roundsFromSession(p)` to its `set({...})`, and also resets `base: null, baseSig: null, baseKaryo: null, baseCalls: null, hatchScores: {}`.

- [ ] **Step 7: Run the helper suites and the build**

Run:

```bash
for s in merge_rounds_helpers cn_recall_helpers merge_rules_helpers hatch_merge_helpers track_helpers integration_helpers experiment_overlays; do web/node_modules/.bin/esbuild tests/$s.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/$s.mjs >/dev/null && node tests/_bundle/$s.mjs | tail -1 || echo "FAILED $s"; done
cd web && npm run build
```

Expected: every suite prints its "all … passed" line, and the build succeeds. Fix type errors where `KaryoMode` is used exhaustively.

- [ ] **Step 8: Commit**

```bash
git add web/src/mergeRounds.ts web/src/cnRecall.ts web/src/hatchMerge.ts web/src/integrationStore.ts web/src/components/RecallImpact.tsx tests/merge_rounds_helpers.ts .gitignore
git commit -m "feat: merge-rounds helpers, base API client and store (rounds saved with the session)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: UI: round strip, adopt, four-way switch, base-driven tracks, noise table

**Files:**
- Create: `web/src/components/RoundStrip.tsx`, `web/src/components/NoiseTable.tsx`
- Modify: `web/src/components/RecallImpact.tsx`, `web/src/components/HatchMergeView.tsx`, `web/src/components/IntegrationSegmentationView.tsx`, `web/src/components/integration.css`

**Interfaces:**
- Consumes the Task 5 store fields and helpers.
- Produces these UI hooks, which the walkthrough in Task 7 relies on:

| Element | Role / attribute | Accessible name |
|---|---|---|
| Round strip | `<p aria-label="Segmentation base">` | |
| Undo button | button | `Undo last round` |
| Back button | button | `Back to original` |
| Adopt button | button | `Adopt as round N` |
| Noise table | `role="region"` | `Call noise by round` |
| Scope toggle | `role="group"` | `Noise scope`, with buttons `this chromosome` / `genome-wide` |
| Noise export button | button | `noise_by_round.tsv` |
| Definitions | `<ul aria-label="Noise definitions">` | |

- [ ] **Step 1: `web/src/components/RoundStrip.tsx`**

```tsx
import { useIntegration } from '../integrationStore'
import { roundStripText } from '../mergeRounds'

/** Section 2: the segmentation base (original → rounds) with undo and back-to-original. */
export default function RoundStrip() {
  const ig = useIntegration()
  return <div className="round-strip">
    <p className="small" aria-label="Segmentation base"><b>Segmentation base:</b>{' '}
      {['Original', ...ig.rounds.map((r, i) => `Round ${i + 1} (${r.label}, ${r.merges} merges)`)].map((t, i, a) =>
        <span key={i} title={i ? JSON.stringify(ig.rounds[i - 1].rules, null, 1) : 'original accepted breakpoints'} className={i === a.length - 1 ? 'round-current' : ''}>{i ? ' → ' : ''}{t}</span>)}
    </p>
    <span className="sr-only">{roundStripText(ig.rounds)}</span>
    <button className="btn-sm" disabled={!ig.rounds.length} onClick={() => ig.undoRound()}>Undo last round</button>
    <button className="btn-sm" disabled={!ig.rounds.length} onClick={() => ig.resetRounds()}>Back to original</button>
    {ig.roundsNotice && <p className="notice small" role="status">{ig.roundsNotice}</p>}
    {ig.baseError && <p className="notice small" role="status">{ig.baseError}</p>}
  </div>
}
```

If `.sr-only` does not already exist in `integration.css`, add it (`position:absolute; width:1px; height:1px; overflow:hidden; clip:rect(0 0 0 0); white-space:nowrap`). Also add `.round-strip { display:flex; flex-wrap:wrap; gap:.5rem; align-items:center }` and `.round-current { font-weight:600 }`.

- [ ] **Step 2: `HatchMergeView.tsx`**

- Import `RoundStrip`, plus `baseHeading` and `wireRounds` from `../mergeRounds`.
- Add `const wire = wireRounds(ig.rounds)`.
- Make these replacements:
  - `const request = { chrom, source, rules: [rule0, ...ig.mergeRules], allow_gap_crossing: ig.allowGapCrossing, rounds: wire }`
  - `const recallBody = { rules: [rule0, ...ig.mergeRules], allow_gap_crossing: ig.allowGapCrossing, rounds: wire }`
  - h3: `<h3>2. Merge rules — exploratory adjacent-segment merge <span className="muted small">· {baseHeading(ig.rounds.length)}</span></h3>`
- Render `<RoundStrip />` directly after the lead `<p>`, before `<MergeRulesEditor />`.
- Status line: `' · compare on the karyograms above (Before · After · Changed)'` stays.
- Below the provenance `details`, change the text "The rule hatches on the Chromosome tracks are computed from the original accepted-breakpoint segmentation" to "…computed from the current segmentation base".

- [ ] **Step 3: `RecallImpact.tsx`**

Replace the mode buttons with `karyoModes`, and add the adopt button after the four exports:

```tsx
import { karyoModes, sameRounds } from '../mergeRounds'
// …
  const modes = karyoModes(ig.rounds.length, hasAfter)
  const onBase = !!st && sameRounds(st.settings.rounds ?? [], ig.rounds)
// in the segmented group:
        {modes.map((m) => <button key={m.mode} className={`btn-sm${ig.karyoMode === m.mode ? ' on' : ''}`} aria-pressed={ig.karyoMode === m.mode} disabled={!m.enabled} onClick={() => ig.set({ karyoMode: m.mode })}>{m.label}</button>)}
// after the cn_change_summary.tsv button:
        {st?.status === 'done' && (onBase
          ? <button className="btn-sm primary" onClick={() => ig.adoptRecall()}>{`Adopt as round ${ig.rounds.length + 1}`}</button>
          : <span className="muted small">computed on a different base</span>)}
```

Remove the now unused `LABELS` constant. Add one line under the switch when rounds exist: `<span className="muted small">Before = {ig.rounds.length ? `round ${ig.rounds.length}` : 'original'}; Original = the on-disk calls.</span>`, shown only when `ig.rounds.length > 0`.

- [ ] **Step 4: `IntegrationSegmentationView.tsx`, base-driven data**

At the top of the component, after `const meta = …`, add:

```tsx
  const roundsSig = JSON.stringify(wireRounds(ig.rounds))
  const nRounds = ig.rounds.length
  const baseBps = nRounds ? (ig.base?.breakpoints ?? null) : (meta?.breakpoints ?? null)
  const baseClusterSegs = nRounds ? (ig.base?.cluster_segments ?? null) : (meta?.cluster_segments ?? null)
  const baseCells = nRounds ? ig.baseKaryo : ig.karyo
  const beforeName = nRounds ? `round ${nRounds}` : 'before'
  useEffect(() => { if (meta) ig.refreshBase() }, [meta, roundsSig])
```

Import `wireRounds` from `../mergeRounds` and `NoiseTable` from `./NoiseTable`. Then:

1. **Hatch-scores effect:** deps become `[seg.chrom, seg.sources, meta, roundsSig]`.
2. **`clSegs`:** `useMemo(() => clusterSegmentsLocal(baseClusterSegs, exMeta.chromosomes), [baseClusterSegs, exMeta.chromosomes])`. This drives the CN strip on the tracks, `beforeHover` and `segTable`.
3. **Tracks:** in `trackFig`, replace `meta.breakpoints` with `baseBps` (the guard is `!baseBps`, the loop is `baseBps.filter(…)`), and add `baseBps` to the deps.
4. **Heading paragraph:** replace `{meta.breakpoints?.length} rows` with `{baseBps?.length ?? '…'} rows{nRounds ? ` after ${nRounds} merge round${nRounds > 1 ? 's' : ''}` : ''}`. When `nRounds` is non-zero, append the sentence " The alternative flank view and the archived merge-prototype tables still use the original accepted breakpoints."
5. **Mode resolution:**

```tsx
  const hasAfter = !!ig.recall?.cluster_segments
  const mode = (ig.karyoMode === 'after' || ig.karyoMode === 'changed') && !hasAfter ? 'before' : ig.karyoMode === 'original' && !nRounds ? 'before' : ig.karyoMode
```

6. **`clusterRows`:**
   - `original` draws `meta.cluster_segments` with the sub-label `original: best_s … · N segments`.
   - `before` draws `baseClusterSegs` with `${nRounds ? `round ${nRounds}: ` : ''}best_s … · N segments`.
   - `after` / `changed` draw the re-call, with `before` = the `baseClusterSegs` state row.
   - If `baseClusterSegs` is null, return `[]`.
   - Deps: `[meta, nBins, mode, recallClusterSegs, baseClusterSegs, nRounds]`.
7. **Cluster hover** in after/changed: the before state comes from `clusterStateRow(baseClusterSegs!, src, nBins)` and the base's `best_s`. The first line is `${beforeName} ${stateLabel(b)} → after ${stateLabel(a)}`, and `best_s ${bsB} → ${bsA}` uses the base's best_s.
8. **`cellRows`:**

```tsx
  const percellAfter = mode === 'after' || mode === 'changed' ? ig.recallKaryo : null
  const percellBefore = mode === 'original' ? ig.karyo : baseCells
  const cellRows: KaryoRow[] = useMemo(() => {
    if (!percellBefore || ((mode === 'after' || mode === 'changed') && !percellAfter)) return []
    return cellOrder.map((i) => {
      const before = percellBefore.subarray(i * nBins, (i + 1) * nBins)
      const pc = recallPerCell?.[i]
      const base = { label: cells[i].cell, band: clusterColor(clusterOf(cells[i], 'wnn_leiden_0.3'), levels) }
      if (!percellAfter) return { ...base, sub: `cluster ${clusterOf(cells[i], 'wnn_leiden_0.3')}${mode === 'original' || !nRounds ? ` · best_s ${fmt(cells[i].percell_best_s, 3)}` : ` · ${beforeName}`}`, data: before }
      return { ...base, sub: `cluster ${clusterOf(cells[i], 'wnn_leiden_0.3')} · best_s ${fmt(pc?.best_s_before ?? null, 3)} → ${fmt(pc?.best_s_after ?? null, 3)} · ${pc?.bins_changed ?? 0} bins changed`, data: percellAfter.subarray(i * nBins, (i + 1) * nBins), before: mode === 'changed' ? before : undefined }
    })
  }, [percellBefore, percellAfter, mode, recallPerCell, cellOrder, cells, levels, nBins, nRounds, beforeName])
```

9. **Per-cell hover:** the before value reads `percellBefore`, with the prefix `${beforeName} … → after …`.
10. **Per-cell placeholder:**
    - In `before` mode with rounds and no `baseKaryo`, show `re-calling base per-cell CN… ${done} / ${total}`, using `ig.baseCalls?.progress`, or `base per-cell CN not available: ${ig.baseCalls?.error ?? ig.baseError ?? ''}` when the job errored.
    - Keep the existing after/changed text.
    - Otherwise show `per-cell CN matrix not loaded`.
11. **Noise table:** render `<NoiseTable />` directly after the "Segments on {chrom}" `</details>`.
12. **Loading notice:** when `nRounds && !ig.base`, show `<p className="muted small">loading segmentation base…</p>` at the top of the tracks panel.

- [ ] **Step 5: `web/src/components/NoiseTable.tsx`**

```tsx
import { useState } from 'react'
import { useIntegration } from '../integrationStore'
import { downloadText } from '../cnRecall'
import { NOISE_DEFINITIONS, noiseCellsN, noiseRows, noiseTsv, roundsHeader, sameRounds, type NoiseColumn } from '../mergeRounds'
import { SOURCE_COLORS } from '../integration'

/** Call noise by round, under the chromosome tracks (spec §4). */
export default function NoiseTable() {
  const ig = useIntegration()
  const [scope, setScope] = useState<'chrom' | 'genome'>('genome')
  const key = scope === 'genome' ? 'genome' : ig.seg.chrom
  const st = ig.recall
  const pending = !!st?.diff?.noise && sameRounds(st.settings.rounds ?? [], ig.rounds)
  const cols: NoiseColumn[] = [
    { name: 'Original', noise: ig.base?.noise_original ?? null },
    ...ig.rounds.map((r, i) => ({ name: `Round ${i + 1}`, noise: r.noise })),
    ...(pending ? [{ name: 'After', noise: st!.diff!.noise!, cellsPending: !st!.percell_ready }] : []),
  ]
  const tsvCols = cols.map((c, i) => ({ ...c, name: i === 0 ? 'original' : c.name === 'After' ? 'after' : `round${i}` }))
  const header = [...roundsHeader(ig.rounds), ...(pending ? [`after rules: ${JSON.stringify(st!.settings.rules)} · allow_gap_crossing: ${st!.settings.allow_gap_crossing}`] : [])]
  return <div className="noise-table" role="region" aria-label="Call noise by round">
    <h4>Call noise by round</h4>
    <div className="controls-row">
      <div className="segmented" role="group" aria-label="Noise scope">
        <button className={`btn-sm${scope === 'chrom' ? ' on' : ''}`} aria-pressed={scope === 'chrom'} onClick={() => setScope('chrom')}>this chromosome</button>
        <button className={`btn-sm${scope === 'genome' ? ' on' : ''}`} aria-pressed={scope === 'genome'} onClick={() => setScope('genome')}>genome-wide</button>
      </div>
      <span className="muted small">{scope === 'genome' ? 'all chromosomes' : ig.seg.chrom} · sources shown on the tracks</span>
    </div>
    {!ig.base ? <p className="muted small">loading…</p> : <div className="table-scroll"><table className="tbl small">
      <thead><tr><th>source</th><th>level</th><th>measure</th>{cols.map((c) => <th key={c.name}>{c.name}</th>)}</tr></thead>
      <tbody>{ig.seg.sources.flatMap((src) => noiseRows(cols, src, key).map((r, i) => <tr key={`${src}-${i}`}>
        {i === 0 && <td rowSpan={8} style={{ color: SOURCE_COLORS[src] }}>{src}</td>}
        {i % 4 === 0 && <td rowSpan={4}>{r.level === 'cluster' ? 'cluster row' : `cells (n = ${noiseCellsN(cols, src, key) ?? '…'})`}</td>}
        <td>{r.measure}</td>{r.values.map((v, j) => <td key={j}>{v}</td>)}
      </tr>))}</tbody>
    </table></div>}
    <ul className="small muted" aria-label="Noise definitions">{NOISE_DEFINITIONS.map((d) => <li key={d}>{d}</li>)}</ul>
    <button className="btn-sm" disabled={!ig.base} onClick={() => downloadText(noiseTsv(tsvCols, ig.seg.sources, header), 'noise_by_round.tsv')}>noise_by_round.tsv</button>
  </div>
}
```

Check that `SOURCE_COLORS` is exported from `../integration` (it is imported from there in `IntegrationSegmentationView.tsx`). Add `.noise-table { margin-top: .75rem }` to `integration.css`.

- [ ] **Step 6: Build and run the helper suites**

Run the Task 5 Step 7 loop and `cd web && npm run build`.
Expected: all suites pass and the build succeeds.

- [ ] **Step 7: Manual check on a private server (8799)**

1. Start `python3 run.py --port 8799` and open it with the Playwright harness, or check its HTML and API with `curl`.
2. Confirm the strip shows `Original` and the noise table shows an Original column.
3. Leave the end-to-end flow to Task 7.
4. Stop the 8799 server.

- [ ] **Step 8: Commit**

```bash
git add web/src/components/RoundStrip.tsx web/src/components/NoiseTable.tsx web/src/components/RecallImpact.tsx web/src/components/HatchMergeView.tsx web/src/components/IntegrationSegmentationView.tsx web/src/components/integration.css
git commit -m "feat: round strip, adopt, Original/Before/After/Changed, base-driven tracks and call-noise table

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Walkthrough, validation record and docs

**Files:**
- Modify: `tests/browser_cn_recall_walkthrough.cjs`
- Regenerate: `validation/cn-recall/cn_recall_check_results.json`, `validation/cn-recall/cn-recall.png`
- Modify: `README.md`, `docs/USAGE.md`, `docs/SCIENTIFIC_WORKFLOW.md`, `docs/ARCHITECTURE.md`, `tests/README.md`

**Interfaces:**
- Consumes the Task 6 UI hooks. The existing walkthrough steps stay as they are: with no rounds, Before equals the original.

- [ ] **Step 1: Extend the walkthrough**

Insert the following before the block that resets τ to 0.2. At that point the re-call at τ 0.5 / ref 8 / any is done and Before has been re-selected.

- **Count the track breakpoints before adopting.** In the tracks Plotly figure, count `layout.shapes` of `type === 'line'` whose `yref` ends with ` domain`. Read it through `page.evaluate` on the tracks plot element (`#cluster-tracks .js-plotly-plot`, via `el.layout.shapes`). Store it as `bpBefore`.
- **Read the Original genome-wide zig-zags.** On the noise region (`getByLabel('Call noise by round')`), click `genome-wide`. Read the first source block's `zig-zags` row, Original column, and store it as `zzOriginal`.
- **Adopt.** Restore τ 0.5 / ref 8 / any if the stale-notice step changed it, so the re-call is on the current rules. Click `Adopt as round 1`. Then:
  - wait until `getByLabel('Segmentation base')` contains `Round 1 (`;
  - wait until the noise table header has a `Round 1` column;
  - read `bpAfter` the same way as `bpBefore`, and the Round 1 zig-zags as `zzRound1`;
  - assert `bpAfter < bpBefore` and `zzRound1 < zzOriginal`;
  - assert `getByLabel('Noise definitions')` is visible and contains `Zig-zag: an internal segment`.
- **Re-call on the base with the same rules.** Click `Merge genome and re-call CN` and wait for `per-cell ✓ · done` (up to 600 s). Sum the impact table's merges column and assert it is `0`.
- **Undo.** Click `Undo last round`. Wait until the strip text is exactly `Segmentation base: Original`. Read the breakpoints again and assert the count equals `bpBefore`.
- **Record the results.** Add these fields to the result JSON: `rounds: {bpBefore, bpAfter, bpUndo, zzOriginal, zzRound1, round2Merges, definitionsVisible}`.

Keep the console-error assertion at the end covering the new steps.

- [ ] **Step 2: Run the walkthroughs on a private server**

```bash
python3 run.py --port 8799 &   # note the PID
PLAYWRIGHT_MODULE=/home/mami/Desktop/EGA/job_search/career-ops/node_modules/playwright \
CHROMIUM_PATH=/home/mami/.cache/ms-playwright/chromium_headless_shell-1234/chrome-headless-shell-linux64/chrome-headless-shell \
BASE_URL=http://127.0.0.1:8799 node tests/browser_cn_recall_walkthrough.cjs
PLAYWRIGHT_MODULE=… CHROMIUM_PATH=… BASE_URL=http://127.0.0.1:8799 node tests/browser_hatch_merge_walkthrough.cjs
kill <that PID>
```

Expected: both print PASS. The validation JSON contains the `rounds` block with `bpAfter < bpBefore`, `zzRound1 < zzOriginal`, `round2Merges === 0` and `bpUndo === bpBefore`.

- [ ] **Step 3: Docs**

- **`README.md`:** one bullet under the re-call feature: "Adopt a re-call as a merge round; section 2 then works on the merged segmentation; undo or go back to the original; call-noise table by round."
- **`docs/USAGE.md`:** a "Merge rounds" subsection:
  - adopt, then new rules, then re-call, then undo;
  - what Original/Before/After/Changed mean;
  - the noise table, with the four definitions verbatim from Global Constraints;
  - the export file name;
  - that rounds are saved with the session;
  - that after a server restart the base's per-cell calls are recomputed (about 60 s).
- **`docs/SCIENTIFIC_WORKFLOW.md`:**
  - every round merges on the original pseudobulk;
  - the same rules twice give 0 merges;
  - zig-zags with both steps pointing the same way are kept by Rule 0 by design, and the spike/dip rule is parked;
  - gain next to loss can be real.
- **`docs/ARCHITECTURE.md`:**
  - `merge_rounds.py` and `call_noise.py`;
  - `rounds` on the routes, `/api/integration/base`, pinning and calls-only jobs;
  - the store fields.
- **`tests/README.md`:** the new pytest files and the `merge_rounds_helpers.ts` command.

- [ ] **Step 4: Full regression**

```bash
python3 -m pytest tests/ -q
for s in merge_rounds_helpers cn_recall_helpers merge_rules_helpers hatch_merge_helpers track_helpers integration_helpers experiment_overlays; do web/node_modules/.bin/esbuild tests/$s.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/$s.mjs >/dev/null && node tests/_bundle/$s.mjs | tail -1 || echo "FAILED $s"; done
cd web && npm run build
```

Expected: everything passes.

- [ ] **Step 5: Commit**

```bash
git add tests/browser_cn_recall_walkthrough.cjs validation/cn-recall/cn_recall_check_results.json validation/cn-recall/cn-recall.png README.md docs/USAGE.md docs/SCIENTIFIC_WORKFLOW.md docs/ARCHITECTURE.md tests/README.md
git commit -m "docs, test: merge-rounds walkthrough, validation record and docs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
