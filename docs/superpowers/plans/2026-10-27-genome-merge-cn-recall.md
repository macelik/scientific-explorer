# Genome-wide merge and CN re-call: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task by task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One action merges every chromosome of every source with the current merge rules. It then re-calls CN per cluster and per cell with pyEpi `assign_gainloss_new`, unmodified. The cluster and per-cell karyograms switch between Before, After and Changed. An impact summary and TSV exports accompany them.

**Architecture:**
- **Pipeline.** A pure module `server/cn_recall.py` holds stages 1–3 and the diff. It takes plain inputs (`RecallInputs`), so tests need no `IntegrationStore`.
- **Server state.** `server/cn_recall_api.py` holds an in-memory job manager. Each job runs on its own thread, keeps the latest 3 results, and reuses a job by hash. It also serves the routes. Per-cell CN calls run in chunks on the existing background process pool.
- **Frontend.** `web/src/cnRecall.ts` has the types, fetchers, pixel and TSV helpers. A `RecallImpact` component holds the switch, the summary and the exports. `IntegrationSegmentationView` and `HatchMergeView` are wired to the integration store.

**Tech stack:**
- Python 3.10, FastAPI, numpy, pandas, h5py, pytest.
- React 18, TypeScript, Zustand, Vite.
- esbuild + node for the TS helper tests; Playwright for the browser walkthrough.

**Spec:** `docs/superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md`

## Global constraints

- **pyEpi unmodified.** It is called only through `server/pyepi_adapter.py`, which takes `package_dir` from `Config.from_env().package_dir`.
- **Originals never written:** `karyogram_cluster_segments.tsv`, `percell_cn_states.npy`, `percell_best_s.tsv`, `raw_accepted_breakpoints.tsv` and the cache npz.
- **Source order:** `['cluster0','cluster1','cluster2','cluster3','cluster4','all']`, the same as `server/integration.SOURCES` and the rows of `pb_x`.
- **Merging** uses the float32 `pb_x` row and the same `run_rule_merge` arguments as `/api/integration/hatch-merge`.
- **CN pseudobulk** is the float64 mean of the member rows of `count_matrix.h5ad`.
- **Cells:** a cell in cluster k uses source `f'cluster{k}'`.
- **Per-cell encoding:** uint8 `rint(2·state)` for finite states and 255 otherwise, as in `IntegrationStore._build_arrays`. In Changed mode, code **252** (`UNCHANGED`) is drawn `#f1f1f1`.
- **Cluster rows:** columns are exactly `group, n_cells, best_s, segment, chromosome, bin_start, bin_end, n_bins, continuous_cn, integer_cn, watson, cn_call`.
  - `best_s` is rounded to 4 decimals.
  - `continuous_cn = round(combined[first bin], 2)`. This copies a quirk of `make_karyograms.py`: it stores the combined state there, not the Watson continuous score.
- **Provenance scope string, verbatim:** `Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on the rule-merged segmentation; originals unchanged.`
- **Results:** kept in memory only, with the 3 most recent kept.
- **Commits:**
  - Stage files by name only, never `git add -A` or `git add .`.
  - Never touch `docs/superpowers/plans/2026-09-21-typesafe-jev-judgment-layers.md` or `docs/superpowers/plans/2026-09-30-portfolio-repository-polish.md`; both are untracked.
  - Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Servers:** never stop or restart the user's server on port 8766. Browser runs use a private `python3 run.py --port 8799`, stopped afterwards.
- **UI strings, verbatim:**

  | Element | String |
  |---|---|
  | Button | `Merge genome & re-call CN` (aria-label `Merge genome and re-call CN`) |
  | Status element | aria-label `Genome re-call status` |
  | Switch group | aria-label `Karyogram mode`, buttons `Before`, `After`, `Changed` |
  | Impact block | aria-label `CN re-call impact` |
  | Stale notice | `Rules changed since this re-call; showing the previous one.` |
  | Cross-chromosome note | `A merge on one chromosome can change calls elsewhere: best_s is fitted over the whole genome.` |

## Review focus

1. **A cell whose cluster has no source:** a label outside 0–4 must fail stage 3 with a message naming it, not KeyError. The test is in Task 1, `test_run_recall_unknown_cluster_fails_with_stage_message`.
2. **A source with no breakpoints on a chromosome:** this is one segment, with no error and 0 merges. Task 1's synthetic data has none for `cluster1` on `chr2`, so this is covered by `test_merge_genome_loose_merges_subset_and_counts`.
3. **The original per-cell matrix is missing:** the cluster results must still be published, and the job ends in `error` with the message starting `percell`. The test is in Task 2, `test_recall_missing_percell_publishes_cluster_then_errors`.
4. **An evicted or unknown job id:** a 404, which the UI shows as a re-call error rather than spinning forever. The test is in Task 2, `test_recall_retention_evicts_oldest`. The UI path is in Task 4's poll error handling.
5. **Editing rules after a run:** the stale notice appears, and a new run replaces the old result without the polls mixing. Checked in the Task 5 walkthrough (`staleAfterEdit`).

---

### Task 1: The pipeline module `server/cn_recall.py`

**Files:**
- Create: `server/cn_recall.py`
- Test: `tests/test_cn_recall.py` (new)

**Interfaces:**
- Consumes:
  - `merge_rules.run_rule_merge(profile, boundaries, n_bins, var_start, var_end, chrom_offset, rules, allow_gap_crossing)`, which returns `{'original','final','steps'}`. `final` holds `{'start','end'}` in chromosome-local bins. `steps[1:]` carry `phase`.
  - `pyepi_adapter.load(pkg_dir)` / `pyepi_adapter.assign_gainloss_new(x, seg)` → `(bin_cn_int, bin_cn_cont, holmes, watson, combined, best_s)`.
- Produces:
  - `SOURCES`, `MISSING = 255`, `RecallError`, `RecallInputs`.
  - `merge_genome(inp, rules, allow_gap_crossing=False, progress=None) -> {'pairs': [...], 'boundaries': {source: {chrom: [int]}}, 'merged_breakpoints': pd.DataFrame}`.
  - `seg_labels(var_seq, boundaries_by_chrom) -> np.ndarray[int]`.
  - `pseudobulk(inp, chunk=250, progress=None) -> np.ndarray (len(sources), n_bins) float64`.
  - `encode_states(v) -> uint8 array`.
  - `cluster_rows(group, n_cells, x, seg, var_seq) -> (rows: list[dict], combined: np.ndarray)`.
  - `cluster_u8_from_segments(df, group, n_bins) -> uint8 array`.
  - `call_cells(rows, seg, package_dir) -> (u8 (k, n_bins), best_s (k,))`. This is top-level and picklable.
  - `add_transitions(acc, before, after)`, `transitions_list(acc) -> [{'before','after','bins'}]`.
  - `run_recall(inp, rules, allow_gap_crossing, executor, update, chunk=100, max_inflight=8) -> {'percell_u8': np.ndarray}`. It reports through `update(**fields)` with the fields `stage`, `progress`, `merge`, `cluster_segments`, `diff` and `percell_ready`.

- [ ] **Step 1: Write the failing tests** in `tests/test_cn_recall.py`:

```python
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
            u8[i:i + 1], bs[i:i + 1] = cr.call_cells(X[i:i + 1], segs[f'cluster{cl[i]}'], PKG)
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest tests/test_cn_recall.py -q`
Expected: collection error, `ModuleNotFoundError: No module named 'server.cn_recall'`.

- [ ] **Step 3: Write `server/cn_recall.py`**

```python
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

from . import merge_rules as mr
from . import pyepi_adapter as sci

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

    def chromosomes(self) -> List[str]:
        return list(dict.fromkeys(self.var_seq.tolist()))


def _ensure(package_dir: str) -> None:
    if not sci._loaded:
        sci.load(package_dir)


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
    diff = {'sources': src_diffs, 'cells': None, 'transitions': {'clusters': transitions_list(acc), 'cells': None}}
    bp_records = merge['merged_breakpoints'].astype(object).where(pd.notna(merge['merged_breakpoints']), None).to_dict(orient='records')
    update(merge={'pairs': merge['pairs'], 'merged_breakpoints': bp_records}, cluster_segments=rows, diff=diff)

    if inp.before_percell_u8 is None:
        raise RecallError('percell: the original per-cell CN matrix is missing')
    update(stage='percell', progress={'done': 0, 'total': n})
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
            source = f'cluster{int(k)}'
            if source not in segs:
                raise RecallError(f'percell: cluster {int(k)} has no source')
        x = inp.read_rows(lo, hi)
        for k in np.unique(cl):
            sel = np.flatnonzero(cl == k)
            while len(pending) >= max_inflight:
                collect(pending.popleft())
            pending.append((lo + sel, executor.submit(call_cells, x[sel], segs[f'cluster{int(k)}'], inp.package_dir)))
    while pending:
        collect(pending.popleft())

    cacc = np.zeros(65536, np.int64)
    per_cell, pcts, dbs = [], [], []
    for i in range(n):
        b, a = np.asarray(inp.before_percell_u8[i]), out_u8[i]
        add_transitions(cacc, b, a)
        changed, called = int((b != a).sum()), int((b != MISSING).sum())
        pcts.append(100 * changed / called if called else 0.0)
        bb, ba = _f(inp.before_percell_best_s[i]), _f(best[i])
        if bb is not None and ba is not None:
            dbs.append(ba - bb)
        per_cell.append({'cell': inp.cell_ids[i], 'wnn_cluster': int(inp.cell_cluster[i]), 'best_s_before': bb,
                         'best_s_after': ba, 'bins_changed': changed})
    diff['cells'] = {'n_cells': n, 'n_changed_cells': sum(1 for r in per_cell if r['bins_changed']),
                     'median_pct_changed': round(float(np.median(pcts)), 3) if pcts else 0.0,
                     'best_s_change': {'median': _f(np.median(dbs)) if dbs else None, 'min': _f(min(dbs)) if dbs else None,
                                       'max': _f(max(dbs)) if dbs else None},
                     'per_cell': per_cell}
    diff['transitions']['cells'] = transitions_list(cacc)
    update(diff=diff, percell_ready=True)
    return {'percell_u8': out_u8}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python3 -m pytest tests/test_cn_recall.py -q`
Expected: all pass. The two real-data tests run when `../integration-story/data` is present, which it is on this machine.

If pyEpi raises on the 70-bin synthetic genome, lengthen the synthetic chromosomes (for example 400 and 300 bins) and scale the breakpoints. Do not change the pipeline.

- [ ] **Step 5: Run the full suite:** `python3 -m pytest -q`. Expected: every test passes.

- [ ] **Step 6: Commit**

```bash
git add server/cn_recall.py tests/test_cn_recall.py
git commit -m "feat: genome-wide merge and CN re-call pipeline (cn_recall)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: The job manager and API, `server/cn_recall_api.py`

**Files:**
- Create: `server/cn_recall_api.py`
- Modify: `server/app.py`. Add the registration after `register_hatch_merge(app, lambda: integration)`.
- Test: append to `tests/test_cn_recall.py`

**Interfaces:**
- Consumes: from Task 1, `cr.run_recall`, `cr.RecallInputs`, `cr.RecallError`, `cr.SOURCES`. From the existing code, `server.integration._csr_rows` and `IntegrationStore` fields: `arrays['pb_x']`, `arrays.get('percell_cn_u8')`, `tables['breakpoints']`, `tables.get('karyogram_cluster_segments')`, `members`, `cell_ids`, `cells['wnn_leiden_0.3']`, `cells['percell_best_s']`, `h5ad_path`, `var_seq`, `var_start`, `var_end` and `identity`.
- Produces:
  - `RecallManager(keep=3)`, with `.start(inp, rules, allow_gap_crossing, identity, executor) -> view`, `.view(jid) -> dict | None` and `.percell_gz(jid) -> bytes | None`.
  - `register_cn_recall(app, get_context, keep=3)`, where `get_context() -> (RecallInputs, identity: dict, executor)` or raises `HTTPException`.
  - `integration_context(get_integration, get_jobs, package_dir)`.
  - Routes: `POST /api/integration/recall`, `GET /api/integration/recall/{jid}` and `GET /api/integration/recall/{jid}/percell.bin`.
  - The view keys are exactly `id, hash, status, stage, progress, error, settings, provenance, merge, cluster_segments, diff, percell_ready`.

- [ ] **Step 1: Write the failing tests** by appending to `tests/test_cn_recall.py`:

```python
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
```

- [ ] **Step 2: Run the tests to verify they fail.** Run `python3 -m pytest tests/test_cn_recall.py -q -k recall_`. Expected: `ModuleNotFoundError: No module named 'server.cn_recall_api'`.

- [ ] **Step 3: Write `server/cn_recall_api.py`**

```python
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
            update(status='error', error=str(e) if isinstance(e, cr.RecallError) else f'{type(e).__name__}: {e}')
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
```

In `server/app.py`, directly after `register_hatch_merge(app, lambda: integration)`:

```python
from .cn_recall_api import integration_context, register_cn_recall
register_cn_recall(app, integration_context(lambda: integration, lambda: jobs, cfg.package_dir))
```

- [ ] **Step 4: Run the tests.** Run `python3 -m pytest tests/test_cn_recall.py -q`, then `python3 -m pytest -q`. Expected: all pass.

- [ ] **Step 5: Smoke test on a private server.** Start `python3 run.py --port 8799` in the background and wait for `/api/health` to report the integration index `done`. Then:

```bash
curl -s -X POST localhost:8799/api/integration/recall -H 'content-type: application/json' \
  -d '{"rules":[{"id":"rule0","type":"flank_threshold","name":"Rule 0","criteria":["delta_mean","srd_phi"],"combine":"any","tau":0.5,"srdRef":8,"maxBins":100,"mergeSignMismatch":true,"mergeSignsOnly":false,"hatch":{"pattern":"x","color":"#b91c1c"}}]}'
```

Poll `GET /api/integration/recall/<id>` until `status` is `done`. That takes about 1–3 min with the process pool. Record:
- the elapsed time;
- the `diff.sources` merges, which should sum to about 383;
- `diff.cells.n_changed_cells`.

Stop the 8799 server. Never touch 8766.

- [ ] **Step 6: Commit**

```bash
git add server/cn_recall_api.py server/app.py tests/test_cn_recall.py
git commit -m "feat: /api/integration/recall jobs for genome merge + CN re-call

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Frontend logic in `web/src/cnRecall.ts`, plus the karyogram pixel rule

**Files:**
- Create: `web/src/cnRecall.ts`
- Modify: `web/src/components/ClusterKaryogram.tsx` (the `KaryoRow` type and the pixel loop)
- Test: `tests/cn_recall_helpers.ts` (new). Add its run line to `tests/README.md` in Task 5.

**Interfaces:**
- Consumes: the Task 2 view JSON, and `ClusterSegment` from `web/src/types.ts`.
- Produces:
  - `UNCHANGED = 252`.
  - Types: `KaryoMode`, `RecallStatus`, `RecallBody`, `RecallSourceDiff`, `RecallTransition`, `RecallCellRow`.
  - Pure helpers: `karyoCode(after, before?)`, `clusterStateRow(segs, group, nBins)`, `stateLabel(code)`, `transitionText(t)`, `recallProgressText(st)`, `recallHeader(st)`, `tsv(header, columns, rows)`.
  - Writers: `mergedBreakpointsTsv(st)`, `recalledClusterSegmentsTsv(st)`, `recalledPercellTsv(st)`, `changeSummaryTsv(st)`.
  - Fetchers: `startRecallApi(body)`, `fetchRecall(id)`, `fetchRecallPercell(id)`. Also `downloadText(text, name)`.
  - `KaryoRow.before?: Uint8Array`.

- [ ] **Step 1: Write the failing test** `tests/cn_recall_helpers.ts`:

```ts
// tests/cn_recall_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/cn_recall_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/cn_recall_helpers.mjs
//           node tests/_bundle/cn_recall_helpers.mjs
import { UNCHANGED, karyoCode, clusterStateRow, stateLabel, transitionText, recallProgressText, tsv, recallHeader, mergedBreakpointsTsv, recalledClusterSegmentsTsv, recalledPercellTsv, changeSummaryTsv, type RecallStatus } from '../web/src/cnRecall'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}

eq(karyoCode(2), 2, 'no before: after code')
eq(karyoCode(2, 2), UNCHANGED, 'same state: unchanged')
eq(karyoCode(3, 2), 3, 'changed: after code')
eq(karyoCode(255, 255), 255, 'missing in both stays missing')
eq(karyoCode(255, 2), 255, 'called -> missing shows missing')

const segs = [
  { group: 'cluster0', n_cells: 5, best_s: 0.5, segment: 0, chromosome: 'chr1', bin_start: 0, bin_end: 2, n_bins: 3, continuous_cn: 1, integer_cn: 2, watson: 1, cn_call: 1 },
  { group: 'cluster0', n_cells: 5, best_s: 0.5, segment: 1, chromosome: 'chr1', bin_start: 3, bin_end: 4, n_bins: 2, continuous_cn: 1.5, integer_cn: 3, watson: 1, cn_call: 1.5 },
  { group: 'all', n_cells: 9, best_s: 0.6, segment: 0, chromosome: 'chr1', bin_start: 0, bin_end: 5, n_bins: 6, continuous_cn: 0, integer_cn: 1, watson: 0, cn_call: 0 },
]
eq([...clusterStateRow(segs, 'cluster0', 6)], [2, 2, 2, 3, 3, 255], 'cluster state row')
eq(stateLabel(3), '1.5 (weak gain)', 'state label')
eq(stateLabel(255), 'not called', 'missing label')
eq(transitionText({ before: 1, after: 1.5, bins: 12 }), '1 (baseline) → 1.5 (weak gain): 12 bins', 'transition text')
eq(transitionText({ before: null, after: 1, bins: 1 }), 'not called → 1 (baseline): 1 bins', 'transition from missing')

const st: RecallStatus = {
  id: 'j1', hash: 'h', status: 'running', stage: 'percell', progress: { done: 1240, total: 3460 }, error: null,
  settings: { rules: [{ id: 'rule0' }], allow_gap_crossing: false }, provenance: { scope: 's' },
  merge: { pairs: [], merged_breakpoints: [{ source: 'cluster0', chromosome: 'chr1', absolute_bin: 3, recursion_depth: 1 }] },
  cluster_segments: segs,
  diff: {
    sources: [{ source: 'cluster0', segments_before: 3, segments_after: 2, merges: 1, phase2: 0, best_s_before: 0.51, best_s_after: 0.5, bins_changed: 2, bins_called: 6, pct_changed: 33.333 }],
    cells: { n_cells: 2, n_changed_cells: 1, median_pct_changed: 5, best_s_change: { median: 0, min: -0.01, max: 0.01 },
      per_cell: [{ cell: 'a', wnn_cluster: 0, best_s_before: 0.5, best_s_after: 0.49, bins_changed: 3 }, { cell: 'b', wnn_cluster: 1, best_s_before: null, best_s_after: 0.6, bins_changed: 0 }] },
    transitions: { clusters: [{ before: 1, after: 1.5, bins: 2 }], cells: [{ before: 0, after: 1, bins: 3 }] },
  },
  percell_ready: false,
}
eq(recallProgressText(st), 'merging ✓ → cluster CN ✓ → per-cell 1,240 / 3,460', 'progress percell')
eq(recallProgressText({ ...st, stage: 'merge', progress: { done: 3, total: 144 } }), 'merging… 3 / 144', 'progress merge')
eq(recallProgressText({ ...st, stage: 'cluster', progress: { done: 250, total: 3460 } }), 'merging ✓ → cluster CN… 250 / 3,460', 'progress cluster')
eq(recallProgressText({ ...st, status: 'done', stage: 'done' }), 'merging ✓ → cluster CN ✓ → per-cell ✓ · done', 'progress done')
eq(recallProgressText({ ...st, status: 'error', error: 'percell: boom' }), 'error: percell: boom', 'progress error')
eq(recallProgressText(null), '', 'progress null')

eq(tsv(['a: 1'], ['x', 'y'], [{ x: 1, y: null }, { x: 'q', y: 2 }]), '# a: 1\nx\ty\n1\t\nq\t2\n', 'tsv writer')
eq(recallHeader(st), ['rules: [{"id":"rule0"}]', 'allow_gap_crossing: false', 'provenance: {"scope":"s"}'], 'header lines')
eq(mergedBreakpointsTsv(st).split('\n').slice(3, 5), ['source\tchromosome\tabsolute_bin\trecursion_depth', 'cluster0\tchr1\t3\t1'], 'merged breakpoints tsv')
eq(recalledClusterSegmentsTsv(st).split('\n')[3], 'group\tn_cells\tbest_s\tsegment\tchromosome\tbin_start\tbin_end\tn_bins\tcontinuous_cn\tinteger_cn\twatson\tcn_call', 'cluster segments columns')
eq(recalledPercellTsv(st).split('\n').slice(3, 6), ['cell\twnn_cluster\tbest_s_before\tbest_s_after\tbins_changed', 'a\t0\t0.5\t0.49\t3', 'b\t1\t\t0.6\t0'], 'percell tsv')
const summary = changeSummaryTsv(st).split('\n')
eq(summary[3], 'source\tsegments_before\tsegments_after\tmerges\tphase2\tbest_s_before\tbest_s_after\tbins_changed\tbins_called\tpct_changed', 'summary columns')
eq(summary.slice(summary.indexOf('scope\tbefore\tafter\tbins')), ['scope\tbefore\tafter\tbins', 'clusters\t1\t1.5\t2', 'cells\t0\t1\t3', ''], 'summary transitions')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) } else console.log('all cn_recall helper checks passed')
```

- [ ] **Step 2: Run it to verify it fails.** Run `web/node_modules/.bin/esbuild tests/cn_recall_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/cn_recall_helpers.mjs`. Expected: the build fails with `Could not resolve "../web/src/cnRecall"`.

- [ ] **Step 3: Write `web/src/cnRecall.ts`**

```ts
// Genome-wide merge + CN re-call: types, fetchers, karyogram pixel rule and TSV exports.
// Spec: docs/superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md
import type { ClusterSegment } from './types'
import { STORY_CN_LABELS } from './integration'

export const UNCHANGED = 252
export type KaryoMode = 'before' | 'after' | 'changed'
export interface RecallBody { rules: unknown[]; allow_gap_crossing: boolean }
export interface RecallPair { source: string; chrom: string; original: number; merged: number; merges: number; phase2: number }
export interface RecallSourceDiff { source: string; segments_before: number; segments_after: number; merges: number; phase2: number; best_s_before: number | null; best_s_after: number | null; bins_changed: number; bins_called: number; pct_changed: number }
export interface RecallTransition { before: number | null; after: number | null; bins: number }
export interface RecallCellRow { cell: string; wnn_cluster: number; best_s_before: number | null; best_s_after: number | null; bins_changed: number }
export interface RecallCells { n_cells: number; n_changed_cells: number; median_pct_changed: number; best_s_change: { median: number | null; min: number | null; max: number | null }; per_cell: RecallCellRow[] }
export interface RecallStatus {
  id: string; hash: string; status: 'running' | 'done' | 'error'; stage: 'merge' | 'cluster' | 'percell' | 'done'
  progress: { done: number; total: number } | null; error: string | null
  settings: { rules: unknown[]; allow_gap_crossing: boolean }; provenance: Record<string, unknown>
  merge: { pairs: RecallPair[]; merged_breakpoints: Record<string, unknown>[] } | null
  cluster_segments: ClusterSegment[] | null
  diff: { sources: RecallSourceDiff[]; cells: RecallCells | null; transitions: { clusters: RecallTransition[]; cells: RecallTransition[] | null } } | null
  percell_ready: boolean
}

async function asJson<T>(r: Response): Promise<T> {
  if (!r.ok) { let msg = `${r.status}`; try { const b = await r.json(); msg = typeof b.detail === 'string' ? b.detail : JSON.stringify(b.detail ?? b) } catch { /* keep status */ } throw new Error(msg) }
  return r.json() as Promise<T>
}
export const startRecallApi = async (body: RecallBody) =>
  asJson<{ id: string; hash: string; status: string }>(await fetch('/api/integration/recall', { method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify(body) }))
export const fetchRecall = async (id: string) => asJson<RecallStatus>(await fetch(`/api/integration/recall/${encodeURIComponent(id)}`))
export async function fetchRecallPercell(id: string): Promise<Uint8Array> {
  const r = await fetch(`/api/integration/recall/${encodeURIComponent(id)}/percell.bin`)
  if (!r.ok) throw new Error(`${r.status}`)
  return new Uint8Array(await r.arrayBuffer())
}

/** Pixel code: with a `before` value (Changed mode), a called bin whose state did not change is drawn as UNCHANGED. */
export function karyoCode(after: number, before?: number): number {
  return before !== undefined && before === after && after !== 255 ? UNCHANGED : after
}
export function clusterStateRow(segs: ClusterSegment[], group: string, nBins: number): Uint8Array {
  const a = new Uint8Array(nBins).fill(255)
  for (const g of segs) if (g.group === group) a.fill(Math.round(g.cn_call * 2), g.bin_start, g.bin_end + 1)
  return a
}
const stateText = (s: number | null) => (s === null ? 'not called' : `${s} (${STORY_CN_LABELS[s]})`)
export const stateLabel = (code: number) => stateText(code === 255 ? null : code / 2)
export const transitionText = (t: RecallTransition) => `${stateText(t.before)} → ${stateText(t.after)}: ${t.bins} bins`

const n = (v: number) => v.toLocaleString('en-US')
export function recallProgressText(st: RecallStatus | null): string {
  if (!st) return ''
  if (st.status === 'error') return `error: ${st.error}`
  const p = st.progress ? ` ${n(st.progress.done)} / ${n(st.progress.total)}` : ''
  if (st.status === 'done') return 'merging ✓ → cluster CN ✓ → per-cell ✓ · done'
  if (st.stage === 'merge') return `merging…${p}`
  if (st.stage === 'cluster') return `merging ✓ → cluster CN…${p}`
  return `merging ✓ → cluster CN ✓ → per-cell${p}`
}

export function tsv(header: string[], columns: string[], rows: Record<string, unknown>[]): string {
  const cell = (v: unknown) => (v === null || v === undefined ? '' : String(v))
  return [...header.map((h) => `# ${h}`), columns.join('\t'), ...rows.map((r) => columns.map((c) => cell(r[c])).join('\t'))].join('\n') + '\n'
}
export const recallHeader = (st: RecallStatus) => [`rules: ${JSON.stringify(st.settings.rules)}`, `allow_gap_crossing: ${st.settings.allow_gap_crossing}`, `provenance: ${JSON.stringify(st.provenance)}`]
const BREAKPOINT_COLUMNS = ['source', 'chromosome', 'absolute_bin', 'bp_genomic_start', 'bp_genomic_end', 'recursion_depth', 'ad_score', 'p_local', 'p_global', 'n_cells']
export const CLUSTER_SEGMENT_COLUMNS = ['group', 'n_cells', 'best_s', 'segment', 'chromosome', 'bin_start', 'bin_end', 'n_bins', 'continuous_cn', 'integer_cn', 'watson', 'cn_call']
export function mergedBreakpointsTsv(st: RecallStatus): string {
  const rows = st.merge?.merged_breakpoints ?? []
  return tsv(recallHeader(st), rows.length ? Object.keys(rows[0]) : BREAKPOINT_COLUMNS, rows)
}
export const recalledClusterSegmentsTsv = (st: RecallStatus) => tsv(recallHeader(st), CLUSTER_SEGMENT_COLUMNS, (st.cluster_segments ?? []) as unknown as Record<string, unknown>[])
export const recalledPercellTsv = (st: RecallStatus) => tsv(recallHeader(st), ['cell', 'wnn_cluster', 'best_s_before', 'best_s_after', 'bins_changed'], (st.diff?.cells?.per_cell ?? []) as unknown as Record<string, unknown>[])
export function changeSummaryTsv(st: RecallStatus): string {
  const cols = ['source', 'segments_before', 'segments_after', 'merges', 'phase2', 'best_s_before', 'best_s_after', 'bins_changed', 'bins_called', 'pct_changed']
  const tr = [...(st.diff?.transitions.clusters ?? []).map((t) => ({ scope: 'clusters', ...t })), ...(st.diff?.transitions.cells ?? []).map((t) => ({ scope: 'cells', ...t }))]
  return tsv(recallHeader(st), cols, (st.diff?.sources ?? []) as unknown as Record<string, unknown>[]) + '\n' + tsv([], ['scope', 'before', 'after', 'bins'], tr)
}
export function downloadText(text: string, name: string) {
  const url = URL.createObjectURL(new Blob([text], { type: 'text/tab-separated-values' })); const a = document.createElement('a')
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000)
}
```

Note: `changeSummaryTsv` ends with `…\n` + `\n` + `scope…`. So after the source rows there is one empty line, then the transitions block, and the test's slice ends with `''`.

- [ ] **Step 4: Change `web/src/components/ClusterKaryogram.tsx`**
  - **Type and import.** Change `export interface KaryoRow { label: string; band?: string; data: Uint8Array; sub?: string }` to `export interface KaryoRow { label: string; band?: string; data: Uint8Array; sub?: string; before?: Uint8Array }`. Add `import { UNCHANGED, karyoCode } from '../cnRecall'`.
  - **Pixel loop.** Replace the two lines

    ```ts
    const v = r.data[colBin[px]]; const o = off + px * 4
    if (v === 255 || v === undefined) { d[o] = 230; d[o + 1] = 230; d[o + 2] = 230 } else { const c = pal[v]; d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2] }
    ```

    with

    ```ts
    const b = colBin[px]; const v = karyoCode(r.data[b], r.before ? r.before[b] : undefined); const o = off + px * 4
    if (v === 255 || v === undefined) { d[o] = 230; d[o + 1] = 230; d[o + 2] = 230 } else if (v === UNCHANGED) { d[o] = 241; d[o + 1] = 241; d[o + 2] = 241 } else { const c = pal[v]; d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2] }
    ```

  - Rows without `before` render exactly as before.

- [ ] **Step 5: Run the helper test and the build**

Run: `web/node_modules/.bin/esbuild tests/cn_recall_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/cn_recall_helpers.mjs && node tests/_bundle/cn_recall_helpers.mjs`
Expected: `all cn_recall helper checks passed`.

Also run the three existing helper suites (commands in `tests/README.md`) and `cd web && npm run build`. Expected: all pass, and the build succeeds.

- [ ] **Step 6: Commit**

```bash
git add web/src/cnRecall.ts web/src/components/ClusterKaryogram.tsx tests/cn_recall_helpers.ts
git commit -m "feat: CN re-call client helpers and Changed-mode karyogram pixels

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: The UI (store, button and status, switch, summary, exports, hovers)

**Files:**
- Modify: `web/src/integrationStore.ts`
- Modify: `web/src/components/HatchMergeView.tsx`
- Create: `web/src/components/RecallImpact.tsx`
- Modify: `web/src/components/IntegrationSegmentationView.tsx`
- Modify: `web/src/components/integration.css` (small styles)

**Interfaces:**
- Consumes: everything from Task 3.
- Produces store fields: `recall: RecallStatus | null`, `recallKaryo: Uint8Array | null`, `recallSignature: string | null`, `recallError: string | null` and `karyoMode: KaryoMode`. It also adds the action `startRecall(body: RecallBody): Promise<void>`. None of these are serialized.

- [ ] **Step 1: Store.** In `web/src/integrationStore.ts`:
  - **Import:** `import { fetchRecall, fetchRecallPercell, startRecallApi, type KaryoMode, type RecallBody, type RecallStatus } from './cnRecall'`.
  - **`IntegrationState`:** add `recall: RecallStatus | null; recallKaryo: Uint8Array | null; recallSignature: string | null; recallError: string | null; karyoMode: KaryoMode; startRecall: (body: RecallBody) => Promise<void>`.
  - **Initial state:** add `recall: null, recallKaryo: null, recallSignature: null, recallError: null, karyoMode: 'before',`.
  - **Action:**

```ts
  startRecall: async (body) => {
    const signature = JSON.stringify(body)
    set({ recallError: null })
    let id: string
    try { id = (await startRecallApi(body)).id } catch (e: any) { set({ recallError: String(e.message || e) }); return }
    if (get().recall?.id !== id) set({ recall: null, recallKaryo: null, karyoMode: 'before' })
    set({ recallSignature: signature })
    const poll = async () => {
      if (get().recallSignature !== signature) return
      let st: RecallStatus
      try { st = await fetchRecall(id) } catch (e: any) { if (get().recallSignature === signature) set({ recallError: `re-call ${id}: ${e.message || e}` }); return }
      if (get().recallSignature !== signature) return
      set({ recall: st })
      if (st.percell_ready && !get().recallKaryo) {
        try { const k = await fetchRecallPercell(id); if (get().recallSignature === signature) set({ recallKaryo: k }) }
        catch (e: any) { set({ recallError: `per-cell re-call: ${e.message || e}` }) }
      }
      if (st.status === 'running') setTimeout(poll, 1000)
    }
    await poll()
  },
```

  `serialize`/`restore` stay unchanged, so the re-call state is not saved.

- [ ] **Step 2: Button and status in `HatchMergeView.tsx`.**
  - **Import:** `import { recallProgressText } from '../cnRecall'`.
  - **Values**, after `const stale = …`:

```ts
  const recallBody = { rules: [rule0, ...ig.mergeRules], allow_gap_crossing: ig.allowGapCrossing }
  const recallSig = JSON.stringify(recallBody)
  const recallRunning = ig.recall?.status === 'running' && ig.recallSignature === recallSig
```

  - **Button**, right after the `Run merge` button inside the same `controls-row`:

```tsx
      <button className="btn" aria-label="Merge genome and re-call CN" disabled={!runnable || recallRunning} onClick={() => ig.startRecall(recallBody)}>Merge genome &amp; re-call CN</button>
```

  - **Status line**, right after that `controls-row` `</div>`:

```tsx
    {(ig.recall || ig.recallError) && <p className={ig.recallError || ig.recall?.status === 'error' ? 'notice small' : 'muted small'} role="status" aria-label="Genome re-call status">
      {ig.recallError ? `error: ${ig.recallError}` : recallProgressText(ig.recall)}
      {ig.recall && ig.recallSignature !== recallSig ? ' · Rules changed since this re-call; showing the previous one.' : ''}
      {ig.recall?.cluster_segments ? ' · compare on the karyograms above (Before · After · Changed)' : ''}
    </p>}
```

  - **Lead paragraph.** Replace the sentence `and CN calls are not recomputed on merged boundaries.` with `and "Run merge" does not recompute CN; "Merge genome & re-call CN" applies these rules to every chromosome of every cluster and re-calls CN with pyEpi on the merged segmentation (compare on the karyograms above).`

- [ ] **Step 3: Create `web/src/components/RecallImpact.tsx`**

```tsx
import { useIntegration } from '../integrationStore'
import { changeSummaryTsv, downloadText, mergedBreakpointsTsv, recallProgressText, recalledClusterSegmentsTsv, recalledPercellTsv, transitionText, type KaryoMode } from '../cnRecall'

const LABELS: Record<KaryoMode, string> = { before: 'Before', after: 'After', changed: 'Changed' }
const f4 = (v: number | null | undefined) => (v === null || v === undefined ? '—' : v.toFixed(4))

/** Before · After · Changed switch for both karyograms, plus the re-call impact summary and exports. */
export default function RecallImpact() {
  const ig = useIntegration(); const st = ig.recall
  const hasAfter = !!st?.cluster_segments
  const d = st?.diff
  return <div className="recall-impact" aria-label="CN re-call impact">
    <div className="controls-row">
      <span className="muted small">CN calls</span>
      <div className="segmented" role="group" aria-label="Karyogram mode">
        {(['before', 'after', 'changed'] as const).map((m) => <button key={m} className={`btn-sm${ig.karyoMode === m ? ' on' : ''}`} aria-pressed={ig.karyoMode === m} disabled={m !== 'before' && !hasAfter} onClick={() => ig.set({ karyoMode: m })}>{LABELS[m]}</button>)}
      </div>
      {st ? <span className="muted small">{recallProgressText(st)}</span> : <span className="muted small">Run "Merge genome &amp; re-call CN" in section 2 (Merge rules) to compare.</span>}
      {ig.karyoMode === 'changed' && <span className="seg-legend"><span><i style={{ background: '#f1f1f1' }} />unchanged</span><span>other colours: the after state of changed bins</span></span>}
    </div>
    {d && <>
      <div className="table-wrap"><table className="tbl small">
        <thead><tr><th>source</th><th>segments before → after</th><th>merges (phase 2)</th><th>best_s before → after</th><th>bins changed</th></tr></thead>
        <tbody>{d.sources.map((s) => <tr key={s.source}><td>{s.source}</td><td>{s.segments_before} → {s.segments_after}</td><td>{s.merges} ({s.phase2})</td><td>{f4(s.best_s_before)} → {f4(s.best_s_after)}</td><td>{s.bins_changed.toLocaleString('en-US')} ({s.pct_changed.toFixed(2)}%)</td></tr>)}</tbody>
      </table></div>
      {d.cells && <p className="small">cells with any change: {d.cells.n_changed_cells.toLocaleString('en-US')} of {d.cells.n_cells.toLocaleString('en-US')} · median {d.cells.median_pct_changed.toFixed(2)}% bins changed · best_s change median {f4(d.cells.best_s_change.median)} [{f4(d.cells.best_s_change.min)}, {f4(d.cells.best_s_change.max)}]</p>}
      <p className="small">Top transitions, clusters: {d.transitions.clusters.slice(0, 5).map(transitionText).join(' · ') || 'none'}{d.transitions.cells ? <> · cells: {d.transitions.cells.slice(0, 5).map(transitionText).join(' · ') || 'none'}</> : null}</p>
      <p className="muted small">A merge on one chromosome can change calls elsewhere: best_s is fitted over the whole genome.</p>
      <div className="controls-row">
        <button className="btn-sm" onClick={() => downloadText(mergedBreakpointsTsv(st!), 'merged_breakpoints.tsv')}>merged_breakpoints.tsv</button>
        <button className="btn-sm" onClick={() => downloadText(recalledClusterSegmentsTsv(st!), 'recalled_cluster_segments.tsv')}>recalled_cluster_segments.tsv</button>
        <button className="btn-sm" disabled={!d.cells} onClick={() => downloadText(recalledPercellTsv(st!), 'recalled_percell_best_s.tsv')}>recalled_percell_best_s.tsv</button>
        <button className="btn-sm" onClick={() => downloadText(changeSummaryTsv(st!), 'cn_change_summary.tsv')}>cn_change_summary.tsv</button>
      </div>
    </>}
  </div>
}
```

  Check `web/src/styles.css` for existing table classes (`tbl`, `table-wrap`) and use whatever the history table in `HatchMergeView.tsx` uses. Add to `integration.css`:

```css
.recall-impact { margin: 6px 0 10px; }
.recall-impact .segmented { display: inline-flex; gap: 0; }
.recall-impact .segmented .btn-sm.on { background: #1f2937; color: #fff; }
.recall-impact .table-wrap { overflow-x: auto; }
```

- [ ] **Step 4: Wire the karyograms in `IntegrationSegmentationView.tsx`.**
  - **Imports:** `import RecallImpact from './RecallImpact'` and `import { clusterStateRow, stateLabel } from '../cnRecall'`.
  - **Effective mode**, near the other derived values:

```ts
  const mode = ig.karyoMode !== 'before' && ig.recall?.cluster_segments ? ig.karyoMode : 'before'
  const clSegsAfter = useMemo(() => clusterSegmentsLocal(ig.recall?.cluster_segments ?? null, exMeta.chromosomes), [ig.recall, exMeta.chromosomes])
```

  - **`clusterRows`:** replace the memo with:

```ts
  const clusterRows: KaryoRow[] = useMemo(() => {
    if (!meta || !meta.cluster_segments) return []
    const after = mode !== 'before' ? ig.recall!.cluster_segments! : null
    return meta.sources.map((src) => {
      const before = clusterStateRow(meta.cluster_segments!, src, nBins)
      const label = `${src} (n=${meta.source_sizes[src]})`, band = SOURCE_COLORS[src]
      if (!after) {
        const own = meta.cluster_segments!.filter((g) => g.group === src)
        return { label, sub: `best_s ${own[0]?.best_s.toFixed(4)} · ${own.length} segments`, data: before, band }
      }
      const own = after.filter((g) => g.group === src)
      return { label, sub: `after re-call: best_s ${own[0]?.best_s.toFixed(4)} · ${own.length} segments`, data: clusterStateRow(after, src, nBins), before: mode === 'changed' ? before : undefined, band }
    })
  }, [meta, nBins, mode, ig.recall])
```

  - **`cellRows`:**

```ts
  const percellAfter = mode !== 'before' ? ig.recallKaryo : null
  const cellRows: KaryoRow[] = useMemo(() => {
    if (!ig.karyo || (mode !== 'before' && !percellAfter)) return []
    return cellOrder.map((i) => {
      const before = ig.karyo!.subarray(i * nBins, (i + 1) * nBins)
      const pc = ig.recall?.diff?.cells?.per_cell[i]
      const base = { label: cells[i].cell, band: clusterColor(clusterOf(cells[i], 'wnn_leiden_0.3'), levels) }
      if (!percellAfter) return { ...base, sub: `cluster ${clusterOf(cells[i], 'wnn_leiden_0.3')} · best_s ${fmt(cells[i].percell_best_s, 3)}`, data: before }
      return { ...base, sub: `cluster ${clusterOf(cells[i], 'wnn_leiden_0.3')} · best_s ${fmt(pc?.best_s_before ?? null, 3)} → ${fmt(pc?.best_s_after ?? null, 3)} · ${pc?.bins_changed ?? 0} bins changed`, data: percellAfter.subarray(i * nBins, (i + 1) * nBins), before: mode === 'changed' ? before : undefined }
    })
  }, [ig.karyo, percellAfter, mode, ig.recall, cellOrder, cells, levels, nBins])
```

  - **Placement:** `<RecallImpact />` goes in the cluster-level panel, right after its `controls-row` `</div>` and before `<ClusterKaryogram rows={clusterRows} …`.
  - **Cluster hover.** Move the current `hoverText` arrow body into a local `const beforeHover = (ri: number, bin: number, chrom: string): string[] => { … }`, unchanged. Then pass:

```tsx
hoverText={(ri, bin, chrom) => {
  if (mode === 'before') return beforeHover(ri, bin, chrom)
  const src = meta.sources[ri]; const c = exMeta.chromosomes.find((x) => x.name === chrom)!; const local = bin - c.offset
  const b = clusterStateRow(meta.cluster_segments!, src, nBins)[bin], a = clusterRows[ri].data[bin]
  const g = clSegsAfter.find((x) => x.group === src && x.chromosome === chrom && local >= x.s && local < x.e)
  const bsB = meta.cluster_segments!.find((x) => x.group === src)?.best_s, bsA = g?.best_s
  return [`before ${stateLabel(b)} → after ${stateLabel(a)}`, ...(g ? [`after segment ${g.segment}: bins [${g.s}, ${g.e}) · ${g.n_bins} bins · integer ${g.integer_cn} · watson ${g.watson}`] : []), `best_s ${bsB?.toFixed(4)} → ${bsA?.toFixed(4)}`]
}}
```

  - **Per-cell hover.** Inside the existing per-cell `hoverText`, prepend this line when `percellAfter` is set: `` `before ${stateLabel(ig.karyo![cellOrder[ri] * nBins + bin])} → after ${stateLabel(percellAfter[cellOrder[ri] * nBins + bin])}` ``. The per-cell `hoverText` then needs its `bin` argument: `(ri, bin) => …`.
  - **Pending per-cell results.** Replace the per-cell fallback `<div className="muted small">per-cell CN matrix not loaded</div>` with:

```tsx
: <div className="muted small">{mode !== 'before' ? `re-calling per-cell CN… ${ig.recall?.stage === 'percell' && ig.recall.progress ? `${ig.recall.progress.done.toLocaleString('en-US')} / ${ig.recall.progress.total.toLocaleString('en-US')}` : ''}` : 'per-cell CN matrix not loaded'}</div>}
```

  - Keep every existing prop on both `ClusterKaryogram` elements: `onClick`, `selectedRow`, `palette`, `rowHeight` and so on.

- [ ] **Step 5: Build and run the helper suites.** Run `cd web && npm run build`, and run all four helper suites. Expected: they pass and the build succeeds.

- [ ] **Step 6: Manual check on a private server.** Start `python3 run.py --port 8799`. Open the Integration tab, set Rule 0 to τ 0.5, ref 8, `any`, and click "Merge genome & re-call CN". Watch the status line advance. Check that After and Changed enable and that both karyograms switch. Stop the 8799 server.

- [ ] **Step 7: Commit**

```bash
git add web/src/integrationStore.ts web/src/components/HatchMergeView.tsx web/src/components/RecallImpact.tsx web/src/components/IntegrationSegmentationView.tsx web/src/components/integration.css
git commit -m "feat: Before/After/Changed karyograms, genome re-call button, impact summary and exports

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Browser walkthrough, validation record and docs

**Files:**
- Create: `tests/browser_cn_recall_walkthrough.cjs`
- Create, generated by the walkthrough: `validation/cn-recall/cn_recall_check_results.json` and `validation/cn-recall/cn-recall.png`
- Modify: `tests/README.md`, `README.md`, `docs/USAGE.md`, `docs/SCIENTIFIC_WORKFLOW.md`, `docs/ARCHITECTURE.md`
- Modify: `docs/superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md` §5. Change "extend `tests/browser_hatch_merge_walkthrough.cjs`" to "a separate `tests/browser_cn_recall_walkthrough.cjs` with the same harness". The ruling: the merge walkthrough is already 267 dense lines, and a separate file keeps both maintainable.

**Interfaces:** consumes the UI strings in Global constraints.

- [ ] **Step 1: Write `tests/browser_cn_recall_walkthrough.cjs`**

```js
// Genome merge + CN re-call walkthrough (2026-10-27): button, status, Before/After/Changed on both karyograms,
// impact summary, export, stale notice.
// Requires a private server (python3 run.py --port 8799; BASE_URL=http://127.0.0.1:8799) and Playwright + Chromium:
//   PLAYWRIGHT_MODULE=/path/to/node_modules/playwright CHROMIUM_PATH=/path/to/chrome BASE_URL=http://127.0.0.1:8799 \
//   node tests/browser_cn_recall_walkthrough.cjs
// Writes validation/cn-recall/cn_recall_check_results.json and validation/cn-recall/cn-recall.png.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const BASE_URL=process.env.BASE_URL||'http://127.0.0.1:8766';
const fs=require('fs'),path=require('path'),assert=require('assert/strict');
const PALETTE=[[123,31,162],[206,147,216],[76,175,80],[255,138,101],[171,0,0]];
(async()=>{
 const browser=await chromium.launch({headless:process.env.HEADED!=='1',executablePath:process.env.CHROMIUM_PATH});
 try {
 const page=await browser.newPage({viewport:{width:1400,height:1100},acceptDownloads:true}),errors=[];
 page.on('pageerror',e=>errors.push(String(e)));
 page.on('console',msg=>{if(msg.type()==='error')errors.push('console: '+msg.text())});
 await page.goto(BASE_URL);await page.getByText('300 matching cells',{exact:true}).waitFor();
 await page.getByRole('button',{name:'Integration · segmentation & CN',exact:true}).click();
 await page.getByText(/Showing 760 of 3460 rows/).waitFor({timeout:60000});
 const clusterPanel=page.locator('.panel',{has:page.getByRole('heading',{name:/Cluster-level CN/})});
 const clusterCanvas=clusterPanel.locator('canvas').first();
 const cellPanel=page.locator('.panel',{has:page.getByRole('heading',{name:/Per-cell CN under the cluster segmentation/})});
 // pixels exactly equal to a palette colour inside the plot area (labels/bands excluded: x >= 164, y >= 20)
 const palettePixels=()=>clusterCanvas.evaluate((cv,P)=>{const d=cv.getContext('2d').getImageData(164,20,cv.width-164,cv.height-20).data;let n=0;for(let i=0;i<d.length;i+=4)if(P.some(p=>p[0]===d[i]&&p[1]===d[i+1]&&p[2]===d[i+2]))n++;return n},PALETTE);
 const digest=()=>clusterCanvas.evaluate(cv=>{const d=cv.getContext('2d').getImageData(0,0,cv.width,cv.height).data;let h=0;for(let i=0;i<d.length;i+=7)h=(h*31+d[i])>>>0;return h});
 const modeBtn=n=>page.getByRole('group',{name:'Karyogram mode'}).getByRole('button',{name:n,exact:true});
 const afterDisabledInitially=await modeBtn('After').isDisabled(),changedDisabledInitially=await modeBtn('Changed').isDisabled();
 const beforeDigest=await digest(),beforePalette=await palettePixels();

 await page.getByLabel('Rule 0 tau').fill('0.5');await page.getByLabel('Rule 0 SRD/√φ ref').fill('8');await page.getByLabel('Rule 0 combine').selectOption('any');
 const panel=page.locator('.hatch-merge');await panel.scrollIntoViewIfNeeded();
 const t0=Date.now();
 await panel.getByRole('button',{name:'Merge genome and re-call CN'}).click();
 const status=page.getByLabel('Genome re-call status');
 await page.waitForFunction(()=>/cluster CN ✓|error:/.test(document.querySelector('[aria-label="Genome re-call status"]')?.textContent||''),{timeout:300000});
 const midText=await status.innerText();
 const afterEnabledMid=!(await modeBtn('After').isDisabled());
 await modeBtn('After').click();
 const perCellPendingShown=(await cellPanel.getByText(/re-calling per-cell CN/).count())>0||/per-cell ✓/.test(await status.innerText());
 await page.waitForFunction(()=>/per-cell ✓ · done|error:/.test(document.querySelector('[aria-label="Genome re-call status"]')?.textContent||''),{timeout:600000});
 const doneText=await status.innerText(),elapsedS=Math.round((Date.now()-t0)/1000);

 const impact=page.getByLabel('CN re-call impact');
 const rows=await impact.locator('tbody tr').evaluateAll(trs=>trs.map(tr=>[...tr.children].map(td=>td.textContent)));
 const merges=rows.reduce((s,r)=>s+parseInt(r[2],10),0),binsChanged=rows.reduce((s,r)=>s+parseInt(r[4].replace(/,/g,''),10),0);
 const cellLine=await impact.getByText(/cells with any change:/).innerText();
 const perCellCanvasAfter=await cellPanel.locator('canvas').count();
 await modeBtn('Changed').click();await page.waitForTimeout(400);
 const changedPalette=await palettePixels();
 await clusterPanel.screenshot({path:path.join('validation','cn-recall','cn-recall.png')}).catch(async()=>{fs.mkdirSync(path.join('validation','cn-recall'),{recursive:true});await clusterPanel.screenshot({path:path.join('validation','cn-recall','cn-recall.png')})});
 const [dl]=await Promise.all([page.waitForEvent('download'),impact.getByRole('button',{name:'merged_breakpoints.tsv'}).click()]);
 const exported=fs.readFileSync(await dl.path(),'utf8').split('\n');
 await modeBtn('Before').click();await page.waitForTimeout(400);
 const backDigest=await digest();
 await page.getByLabel('Rule 0 tau').fill('0.4');await page.waitForTimeout(300);
 const staleAfterEdit=(await status.innerText()).includes('Rules changed since this re-call; showing the previous one.');
 await page.getByLabel('Rule 0 tau').fill('0.2');await page.getByLabel('Rule 0 SRD/√φ ref').fill('3.3');await page.getByLabel('Rule 0 combine').selectOption('all');

 const result={date:new Date().toISOString(),settings:{tau:0.5,srdRef:8,combine:'any'},elapsedS,midText,doneText,merges,binsChanged,cellLine,
   afterDisabledInitially,changedDisabledInitially,afterEnabledMid,perCellPendingShown,perCellCanvasAfter,beforePalette,changedPalette,
   backToBeforeIdentical:backDigest===beforeDigest,exportHeader:exported.slice(0,4),staleAfterEdit,errors};
 fs.mkdirSync(path.join('validation','cn-recall'),{recursive:true});
 fs.writeFileSync(path.join('validation','cn-recall','cn_recall_check_results.json'),JSON.stringify(result,null,2));
 assert.equal(afterDisabledInitially,true);assert.equal(changedDisabledInitially,true);
 assert.ok(!doneText.includes('error:'),doneText);
 assert.equal(afterEnabledMid,true,'After enables once cluster CN is ready');
 assert.equal(perCellPendingShown,true);
 assert.ok(merges>0,'a loose Rule 0 should merge something genome-wide');
 assert.ok(binsChanged>0,'re-calling should change at least one cluster bin');
 assert.ok(perCellCanvasAfter>=1,'per-cell karyogram should draw in After mode');
 assert.ok(changedPalette>0&&changedPalette<beforePalette,'Changed draws only changed bins in colour');
 assert.equal(result.backToBeforeIdentical,true,'Before must equal the original rendering');
 assert.ok(exported[0].startsWith('# rules: ')&&exported[3].includes('absolute_bin'),'export header + columns');
 assert.equal(staleAfterEdit,true,'editing a rule after a re-call shows the stale notice');
 assert.deepEqual(errors,[],'console/page errors present');
 console.log('PASS cn re-call: button, status, Before/After/Changed, summary, export, stale notice');
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1)});
```

- [ ] **Step 2: Run it against a private server**

```bash
python3 run.py --port 8799   # in the background; wait for /api/health integration index "done"
PLAYWRIGHT_MODULE=/home/mami/Desktop/EGA/job_search/career-ops/node_modules/playwright \
CHROMIUM_PATH=/home/mami/.cache/ms-playwright/chromium_headless_shell-1234/chrome-headless-shell-linux64/chrome-headless-shell \
BASE_URL=http://127.0.0.1:8799 node tests/browser_cn_recall_walkthrough.cjs
```

Expected: `PASS cn re-call: …`. Then run `tests/browser_hatch_merge_walkthrough.cjs` against the same server with the same env. Expected: PASS, with no regression. Stop the 8799 server, and remove any stray `.playwright-mcp/` directory.

- [ ] **Step 3: Docs**
  - **`tests/README.md`:** add the `cn_recall_helpers.ts` run line, next to the other helper suites, and the walkthrough run line, next to the other browser scripts.
  - **`docs/USAGE.md`:** a short "Merge genome & re-call CN" subsection under the merge-rules section. Cover:
    - the button and the status stages;
    - the Before · After · Changed switch, the summary and the four exports;
    - "in memory only; the latest 3 kept; lost on restart".
  - **`docs/SCIENTIFIC_WORKFLOW.md`:** a paragraph. It says that the re-call uses pyEpi `assign_gainloss_new` unmodified, the same as `make_karyograms.py`, and that the no-op identity test reproduces the original files. It says that `best_s` is fitted genome-wide, so one merge can change calls elsewhere. And it says these are exploratory diagnostics.
  - **`docs/ARCHITECTURE.md`:** add `server/cn_recall.py`, `server/cn_recall_api.py` and the three routes, plus `web/src/cnRecall.ts` and `RecallImpact.tsx`.
  - **`README.md`:** add one feature bullet.

- [ ] **Step 4: Final regression.** Run `python3 -m pytest -q`, all four TS helper suites and `cd web && npm run build`. Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add tests/browser_cn_recall_walkthrough.cjs validation/cn-recall/cn_recall_check_results.json validation/cn-recall/cn-recall.png tests/README.md README.md docs/USAGE.md docs/SCIENTIFIC_WORKFLOW.md docs/ARCHITECTURE.md docs/superpowers/specs/2026-10-27-genome-merge-cn-recall-design.md
git commit -m "docs, test: genome merge + CN re-call walkthrough, validation record and docs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
