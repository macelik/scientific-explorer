# Allelic Segmentation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A new tab, "Integration · allelic segmentation". It runs the depth run's constrained AD scan on a chosen per-cluster allelic pseudobulk signal and shows the allelic tracks and segments next to the depth tracks and segments, with three comparison tables.

**Architecture:**
- **Server:** three new modules.
  - `allelic_scan.py`: a copy of the constrained `recursive_getbp`, with statistics from `pyepi_adapter`.
  - `allelic_data.py`: load checks, streamed ΣHA/ΣHB per source with a disk cache, and the four signals.
  - `allelic_api.py`: routes plus a cached one-at-a-time job manager running on the existing background process pool.
- **Web:**
  - `allelic.ts`: types, API calls, settings checks, matching and the comparison tables.
  - `allelicStore.ts`: a zustand store, saved with the session under `allelic`.
  - `AllelicSegmentationView.tsx`: the tab.
- **Gate:** before any API or UI work, a validation script must show that the copied scan reproduces all depth breakpoints.

**Tech Stack:** Python 3.10, FastAPI, numpy, pandas, pyarrow, h5py, pytest; React 18, zustand 4, Plotly 2.35, TypeScript 5.6, esbuild (tests), Playwright (walkthrough).

**Spec:** `docs/superpowers/specs/2026-10-28-allelic-segmentation-design.md`

## Global Constraints

**Process and repository:**
- **Branch:** `allelic-segmentation`, cut from `merge-rules`. Do not merge or push.
- **Staging:** stage files by name only. Never `git add -A` or `git add .`.
- **Untracked files:** never stage or touch `docs/superpowers/plans/2026-09-21-typesafe-jev-judgment-layers.md` or `docs/superpowers/plans/2026-09-30-portfolio-repository-polish.md`.
- **Commit trailer:** every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Repo metadata:** do not change it (remotes, git config, CI files).
- **Servers:**
  - The server on port 8766 belongs to the user. Never stop, restart or signal it.
  - For live checks, start a private server with `python3 run.py --port 8799` from the repo root, and stop only that process when done.
- **Subagents:** implementers never dispatch subagents.

**Data and statistics:**
- **Allelic data:** `HA.parquet`, `HB.parquet` and `AT.parquet` in `<EGA>/haplotype_aware_cnv/hapcnv_outs_from_updated_phasing/aligned_matrices/`, 25,337 bins × 3,460 barcodes.
  - The default path is `<data_root>/../../../haplotype_aware_cnv/hapcnv_outs_from_updated_phasing/aligned_matrices`. `PYEPI_SCIENTIFIC_ALLELIC` overrides it.
- **Sources:** `["cluster0", "cluster1", "cluster2", "cluster3", "cluster4", "all"]`. Membership is the integration store's `members`.
- **Signals:** ids `dbaf`, `abs_dbaf`, `shrunk_abs_dbaf`, `oriented_all`, with the formulas in spec §2. A bin is kept when AT ≥ `min_at`; for `oriented_all` it is also kept only when ΣHB_all ≠ ΣHA_all.
- **Settings defaults:**

  | Setting | Default | Range |
  |---|---|---|
  | min_at | 10 | integer 1–1000 |
  | k | 6 | integer 1–8 |
  | n_permutations | 1000 | integer 100–5000 |
  | alpha | 0.001 | 0 < alpha ≤ 0.1 |
  | min_child_bins | 20 | integer 1–500 |
  | min_terminal_bins | 50 | integer 0–500 |
  | global_pool | `genome` | `genome` or `chromosome` |

- **Scan:** a copy of `recursive_getbp` and `validate_candidate_geometry` from `SNP_project/beta_revised/src/pyEpiAneufinder/get_breakpoints.py`.
  - The statistics are `sci.seq_dist_ad`, `sci.permutation_test_ad` and `sci.global_permutation_test_ad` from `server/pyepi_adapter.py`, with their default `random_state`.
  - It is called per chromosome as `recursive_getbp(seq, pool, …)`, with `chromosome_start=0` and `chromosome_end=len(seq)` in kept-bin space.
  - **Global pool:** with `global_pool = genome` (the default, original pyEpi), `pool` is the source's kept-bin values over all chromosomes in genome order. Original `pyEpiAneufinder.py` passes the whole cell, `counts.X[i, :]`. With `chromosome`, `pool = seq`, as `run_bootstrap_diagnostics.py` made the depth breakpoints.
- **Bins:** positions are chromosome-local, like `absolute_bin` in `raw_accepted_breakpoints.tsv`. Integration `cluster_segments` use genome-absolute inclusive `bin_start`/`bin_end`; convert them with `clusterSegmentsLocal`.

**User-facing text (verbatim):**
- tab label `Integration · allelic segmentation`
- "Segments shown are from an earlier run (signal …, settings …). Run again to update."
- "Run segmentation to see allelic segments."
- "Saved allelic settings could not be applied: …"
- export file name `allelic_vs_depth.tsv`

## Review Focus

1. **Cluster 4 and low-coverage chromosomes.** A source with fewer than 2 kept bins on a chromosome must give one segment with `kept_bins` 0 or 1, no breakpoints and no crash. Pinned in Task 1 (`test_no_or_one_kept_bin`).
2. **Barcode order.** The parquet column order differs from the h5mu order. Sums must follow membership by name, never by position. Pinned in Task 2 (`test_sums_follow_barcodes_not_column_order`).
3. **A gap breakpoint next to a depth breakpoint.** A depth breakpoint inside an allelic gap must match at distance 0, and one N bins outside the gap must still match. Pinned in Task 4 helper tests.
4. **Switching signal or settings while a job runs.** The old poll must stop and the shown result must not jump to the superseded job. Pinned in Task 4 store test (`run supersedes`).
5. **A restored session with junk allelic settings.** Fall back to the defaults with the notice. Pinned in Task 4 (`allelicFromSession`).

---

### Task 1: Constrained scan module and the depth-reproduction gate

**Files:**
- Create: `server/allelic_scan.py`
- Create: `tests/test_allelic_scan.py`
- Create: `tests/validate_allelic_depth_reproduction.py`
- Create (by running the script): `validation/allelic/depth_reproduction.json`

**Interfaces:**
- Consumes: `server/pyepi_adapter.py` (`load`, `_loaded`, `seq_dist_ad`, `permutation_test_ad`, `global_permutation_test_ad`) and `server/config.Config.from_env().package_dir`.
- Produces:
  - `allelic_scan.DEFAULTS: dict` with keys k, n_permutations, alpha, min_child_bins, min_terminal_bins.
  - `allelic_scan.scan_chromosome(values, kept, settings, ha=None, hb=None, package_dir=None, pool=None) -> {'audit': list[dict], 'breakpoints': list[dict], 'segments': list[dict]}`. `pool` is the global-test pool; `None` means the chromosome's own kept values.
    - Breakpoint keys: `bin, gap_start, gap_end, recursion_depth, ad_score, p_local, p_global`.
    - Segment keys: `bin_start, bin_end, n_bins, kept_bins, mean, median`, plus `ha_sum, hb_sum, at_sum` when `ha`/`hb` are given.
    - All bins are chromosome-local; `bin_end` is inclusive.
  - `allelic_scan.scan_task(task: tuple) -> (source, chrom, result)`, a picklable worker. `task = (source, chrom, values, kept, ha, hb, settings, package_dir, pool)`.

- [ ] **Step 1: Write the failing tests**

`tests/test_allelic_scan.py`:

```python
import numpy as np

from server import allelic_scan as sc

FAST = {'k': 1, 'n_permutations': 200, 'alpha': 0.01, 'min_child_bins': 20, 'min_terminal_bins': 20}


def _step(n_left=80, n_right=80, lo=0.05, hi=0.30, sd=0.01, seed=1):
    rng = np.random.default_rng(seed)
    return np.r_[rng.normal(lo, sd, n_left), rng.normal(hi, sd, n_right)]


def test_geometry_rules_match_the_depth_run():
    ok, why = sc.validate_candidate_geometry(10, 0, 100, 0, 100, 20, 50)
    assert not ok and why == ('left_child_too_short', 'left_terminal_too_short')
    ok, why = sc.validate_candidate_geometry(30, 20, 100, 0, 100, 20, 50)
    assert ok and why == ()          # internal left edge: no terminal rule
    ok, why = sc.validate_candidate_geometry(60, 0, 100, 0, 100, 20, 50)
    assert not ok and why == ('right_terminal_too_short',)


def test_step_gives_one_breakpoint_and_two_segments():
    x = _step()
    out = sc.scan_chromosome(x, np.ones(len(x), bool), FAST)
    assert [b['bin'] for b in out['breakpoints']] == [80]
    b = out['breakpoints'][0]
    assert b['gap_start'] is None and b['gap_end'] is None and b['recursion_depth'] == 0
    assert b['p_local'] < 0.01 and b['p_global'] < 0.01
    segs = out['segments']
    assert [(s['bin_start'], s['bin_end'], s['kept_bins']) for s in segs] == [(0, 79, 80), (80, 159, 80)]
    assert abs(segs[0]['mean'] - 0.05) < 0.01 and abs(segs[1]['mean'] - 0.30) < 0.01


def test_skipped_bins_map_back_with_a_gap():
    x = _step()
    kept = np.ones(len(x), bool)
    kept[75:85] = False                       # the change at 80 falls inside skipped bins
    out = sc.scan_chromosome(x, kept, FAST)
    assert len(out['breakpoints']) == 1
    b = out['breakpoints'][0]
    assert (b['bin'], b['gap_start'], b['gap_end']) == (85, 75, 85)
    segs = out['segments']
    assert [(s['bin_start'], s['bin_end'], s['n_bins'], s['kept_bins']) for s in segs] == [(0, 84, 85, 75), (85, 159, 75, 75)]


def test_terminal_rule_blocks_a_step_near_the_end():
    x = _step(n_left=10, n_right=150)
    out = sc.scan_chromosome(x, np.ones(len(x), bool), {**FAST, 'min_child_bins': 5, 'min_terminal_bins': 20})
    assert 10 not in [b['bin'] for b in out['breakpoints']]
    top = [r for r in out['audit'] if r['candidate_rank'] == 1][0]
    assert top['candidate_absolute_bin'] == 10 and 'left_terminal_too_short' in top['rejection_reasons']


def test_noise_stops_at_the_first_tested_candidate():
    x = np.random.default_rng(3).normal(0.1, 0.01, 160)
    out = sc.scan_chromosome(x, np.ones(len(x), bool), {**FAST, 'k': 3})
    assert out['breakpoints'] == []
    tested = [r for r in out['audit'] if r['selected_for_significance']]
    assert len(tested) == 1 and tested[0]['recursion_depth'] == 0 and tested[0]['p_local'] >= 0.01


def test_no_or_one_kept_bin():
    x = _step()
    none = sc.scan_chromosome(x, np.zeros(len(x), bool), FAST)
    assert none['breakpoints'] == [] and none['audit'] == []
    assert none['segments'] == [{'bin_start': 0, 'bin_end': 159, 'n_bins': 160, 'kept_bins': 0, 'mean': None, 'median': None}]
    one = np.zeros(len(x), bool); one[5] = True
    s = sc.scan_chromosome(x, one, FAST)['segments']
    assert len(s) == 1 and s[0]['kept_bins'] == 1 and abs(s[0]['mean'] - x[5]) < 1e-12


def test_sums_cover_kept_bins_only():
    x = _step()
    kept = np.ones(len(x), bool); kept[:5] = False
    ha = np.ones(len(x)); hb = np.full(len(x), 2.0)
    s = sc.scan_chromosome(x, kept, FAST, ha=ha, hb=hb)['segments'][0]
    assert (s['ha_sum'], s['hb_sum'], s['at_sum']) == (75.0, 150.0, 225.0)


def test_segments_partition_the_chromosome_and_run_is_deterministic():
    x = np.r_[_step(seed=4), _step(seed=5) + 0.1]
    settings = {**FAST, 'k': 3}
    a = sc.scan_chromosome(x, np.ones(len(x), bool), settings)
    b = sc.scan_chromosome(x, np.ones(len(x), bool), settings)
    assert a == b
    segs = a['segments']
    assert segs[0]['bin_start'] == 0 and segs[-1]['bin_end'] == len(x) - 1
    assert all(segs[i]['bin_end'] + 1 == segs[i + 1]['bin_start'] for i in range(len(segs) - 1))


def test_scan_task_round_trip():
    x = _step()
    src, chrom, out = sc.scan_task(('cluster1', 'chr2', x, np.ones(len(x), bool), None, None, FAST, None, None))
    assert (src, chrom) == ('cluster1', 'chr2') and [b['bin'] for b in out['breakpoints']] == [80]


def test_global_test_uses_the_given_pool(monkeypatch):
    sc._ensure()
    seen = []
    real = sc.sci.global_permutation_test_ad

    def spy(seq_data, *a, **k):
        seen.append(len(seq_data))
        return real(seq_data, *a, **k)
    monkeypatch.setattr(sc.sci, 'global_permutation_test_ad', spy)
    x = _step()
    kept = np.ones(len(x), bool); kept[:4] = False
    sc.scan_chromosome(x, kept, FAST)
    pool = np.random.default_rng(9).normal(0.1, 0.05, 500)
    sc.scan_chromosome(x, kept, FAST, pool=pool)
    assert seen == [156, 500]                 # chromosome kept bins, then the given genome pool
```

- [ ] **Step 2: Run them to check they fail**

Run: `python3 -m pytest -q tests/test_allelic_scan.py`
Expected: FAIL, `ImportError: cannot import name 'allelic_scan'`.

- [ ] **Step 3: Write the module**

`server/allelic_scan.py`:

```python
"""Constrained AD breakpoint scan of one chromosome track (spec 2026-10-28-allelic-segmentation-design.md §3).

`validate_candidate_geometry` and `recursive_getbp` are copied from
SNP_project/beta_revised/src/pyEpiAneufinder/get_breakpoints.py, the function
run_bootstrap_diagnostics.py used to write raw_accepted_breakpoints.tsv. Only
the recursion and the geometry check are copied; the statistics are pyEpi's
own, called through pyepi_adapter.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from . import pyepi_adapter as sci
from .config import Config

DEFAULTS = {'k': 6, 'n_permutations': 1000, 'alpha': 0.001, 'min_child_bins': 20, 'min_terminal_bins': 50}


def _ensure(package_dir: Optional[str] = None) -> None:
    if not sci._loaded:
        sci.load(package_dir or Config.from_env().package_dir)


def validate_candidate_geometry(local_cut, node_start, node_end, chromosome_start, chromosome_end,
                                min_child_bins, min_terminal_bins):
    """Copied verbatim (logic) from beta_revised get_breakpoints.py; terminal rules apply only to
    children touching the true chromosome start/end."""
    absolute_cut = node_start + local_cut
    left_child_bins = local_cut
    right_child_bins = (node_end - node_start) - local_cut
    reasons = []
    if left_child_bins < min_child_bins:
        reasons.append("left_child_too_short")
    if right_child_bins < min_child_bins:
        reasons.append("right_child_too_short")
    if node_start == chromosome_start and (absolute_cut - chromosome_start) < min_terminal_bins:
        reasons.append("left_terminal_too_short")
    if node_end == chromosome_end and (chromosome_end - absolute_cut) < min_terminal_bins:
        reasons.append("right_terminal_too_short")
    return (len(reasons) == 0, tuple(reasons))


def recursive_getbp(seq_data, seq_data_cell, k=3, depth=0, node_start=0, n_permutations=500,
                    alpha=0.05, chromosome_start=None, chromosome_end=None,
                    min_child_bins=20, min_terminal_bins=50, audit=None):
    """Copied from beta_revised get_breakpoints.py. Candidates in descending AD order; the first
    geometry-valid one is tested; if it fails, the node stops. Every examined candidate is audited."""
    if depth >= k or len(seq_data) <= 1:
        return []
    if chromosome_start is None:
        chromosome_start = 0
    if chromosome_end is None:
        chromosome_end = len(seq_data)
    if audit is None:
        audit = []
    dist_vector = np.array(sci.seq_dist_ad(seq_data))
    if dist_vector.size == 0:
        return []
    n = len(seq_data)
    node_end = node_start + n
    order = np.argsort(-dist_vector, kind="stable") + 1
    for rank, i in enumerate(order, start=1):
        i = int(i)
        abs_cut = node_start + i
        ad_score = float(dist_vector[i - 1])
        geometry_valid, reasons = validate_candidate_geometry(
            i, node_start, node_end, chromosome_start, chromosome_end, min_child_bins, min_terminal_bins)
        row = {
            "source": None, "chromosome": None, "recursion_depth": depth,
            "node_start_bin": node_start, "node_end_bin": node_end,
            "candidate_rank": rank, "candidate_local_cut": i, "candidate_absolute_bin": abs_cut,
            "ad_score": ad_score, "left_child_bins": i, "right_child_bins": n - i,
            "geometry_valid": geometry_valid, "rejection_reasons": reasons,
            "selected_for_significance": False, "p_local": np.nan, "p_global": np.nan, "accepted": False,
        }
        audit.append(row)
        if not geometry_valid:
            continue
        row["selected_for_significance"] = True
        bp_dist = dist_vector[i - 1]
        left_seg = seq_data[:i]
        right_seg = seq_data[i:]
        pval = sci.permutation_test_ad(left_seg, right_seg, bp_dist, n_permutations)
        row["p_local"] = pval
        if pval >= alpha:
            return []
        pval_global = sci.global_permutation_test_ad(seq_data_cell, len(left_seg), len(right_seg), bp_dist, n_permutations)
        row["p_global"] = pval_global
        accepted = (pval < alpha) and (pval_global < alpha)
        row["accepted"] = accepted
        if not accepted:
            return []
        bp_global = node_start + i
        left_bps = recursive_getbp(left_seg, seq_data_cell, k=k, depth=depth + 1, node_start=node_start,
                                   n_permutations=n_permutations, alpha=alpha,
                                   chromosome_start=chromosome_start, chromosome_end=chromosome_end,
                                   min_child_bins=min_child_bins, min_terminal_bins=min_terminal_bins, audit=audit)
        right_bps = recursive_getbp(right_seg, seq_data_cell, k=k, depth=depth + 1, node_start=node_start + i,
                                    n_permutations=n_permutations, alpha=alpha,
                                    chromosome_start=chromosome_start, chromosome_end=chromosome_end,
                                    min_child_bins=min_child_bins, min_terminal_bins=min_terminal_bins, audit=audit)
        return left_bps + [(bp_global, bp_dist, pval, pval_global)] + right_bps
    audit.append({
        "source": None, "chromosome": None, "recursion_depth": depth,
        "node_start_bin": node_start, "node_end_bin": node_end,
        "candidate_rank": None, "candidate_local_cut": None, "candidate_absolute_bin": None,
        "ad_score": np.nan, "left_child_bins": None, "right_child_bins": None,
        "geometry_valid": False, "rejection_reasons": ("no_valid_candidate",),
        "selected_for_significance": False, "p_local": np.nan, "p_global": np.nan, "accepted": False,
    })
    return []


def _plain(v):
    """JSON-safe scalar: numpy → Python, NaN → None."""
    if v is None:
        return None
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return None if not np.isfinite(v) else float(v)
    return v


def scan_chromosome(values, kept, settings, ha=None, hb=None, package_dir=None, pool=None) -> dict:
    """Scan one chromosome on its kept bins; map cuts back to chromosome-local bins.

    `pool` is what the global permutation test samples from: original pyEpi passes the
    whole cell genome-wide (pyEpiAneufinder.py, `counts.X[i, :]`); None uses this
    chromosome's kept values, as run_bootstrap_diagnostics.py did for the depth breakpoints.

    A cut between kept bins i and i+1 sits at the bin of kept bin i+1; skipped bins
    between them are reported as [gap_start, gap_end]. Segment sums and stats cover kept bins.
    """
    _ensure(package_dir)
    values = np.asarray(values, dtype=np.float64)
    kept = np.asarray(kept, dtype=bool)
    n = len(values)
    idx = np.flatnonzero(kept)
    seq = values[idx]
    audit: list = []
    accepted = []
    if len(seq) > 1:
        cell = seq if pool is None else np.asarray(pool, dtype=np.float64)
        accepted = recursive_getbp(seq, cell, k=int(settings['k']), depth=0, node_start=0,
                                   n_permutations=int(settings['n_permutations']), alpha=float(settings['alpha']),
                                   chromosome_start=0, chromosome_end=len(seq),
                                   min_child_bins=int(settings['min_child_bins']),
                                   min_terminal_bins=int(settings['min_terminal_bins']), audit=audit)
    rows = []
    for r in audit:
        cut = r['candidate_absolute_bin']
        row = {key: _plain(v) for key, v in r.items() if key != 'rejection_reasons'}
        row['rejection_reasons'] = ';'.join(r['rejection_reasons'])
        row['bin'] = int(idx[cut]) if cut is not None and 0 <= cut < len(idx) else None
        rows.append(row)
    depth_of = {r['candidate_absolute_bin']: r['recursion_depth'] for r in audit if r['accepted']}
    bps = []
    for cut, ad, pl, pg in accepted:
        b, prev = int(idx[cut]), int(idx[cut - 1])
        gap = (prev + 1, b) if b - prev > 1 else (None, None)
        bps.append({'bin': b, 'gap_start': gap[0], 'gap_end': gap[1], 'recursion_depth': int(depth_of[cut]),
                    'ad_score': float(ad), 'p_local': float(pl), 'p_global': float(pg)})
    bounds = [0] + [c for c, *_ in accepted] + [len(idx)]
    segs = []
    for j in range(len(bounds) - 1):
        k0, k1 = bounds[j], bounds[j + 1]
        b0 = 0 if j == 0 else int(idx[k0])
        b1 = n - 1 if j == len(bounds) - 2 else int(idx[k1]) - 1
        part, kb = seq[k0:k1], idx[k0:k1]
        row = {'bin_start': b0, 'bin_end': b1, 'n_bins': b1 - b0 + 1, 'kept_bins': int(k1 - k0),
               'mean': float(part.mean()) if len(part) else None,
               'median': float(np.median(part)) if len(part) else None}
        if ha is not None and hb is not None:
            hs, bs = float(np.asarray(ha)[kb].sum()), float(np.asarray(hb)[kb].sum())
            row.update({'ha_sum': hs, 'hb_sum': bs, 'at_sum': hs + bs})
        segs.append(row)
    return {'audit': rows, 'breakpoints': bps, 'segments': segs}


def scan_task(task):
    """Process-pool worker: (source, chrom, values, kept, ha, hb, settings, package_dir, pool)."""
    source, chrom, values, kept, ha, hb, settings, package_dir, pool = task
    return source, chrom, scan_chromosome(values, kept, settings, ha=ha, hb=hb, package_dir=package_dir, pool=pool)
```

- [ ] **Step 4: Run the tests to check they pass**

Run: `python3 -m pytest -q tests/test_allelic_scan.py`
Expected: 10 passed.

If `test_step_gives_one_breakpoint_and_two_segments` finds a cut at 79 or 81 rather than 80, do not loosen the expected values. Print `seq_dist_ad` around the step and report the off-by-one: the 1-based cut convention must match `recursive_getbp`.

- [ ] **Step 5: Write the depth-reproduction gate**

`tests/validate_allelic_depth_reproduction.py`:

```python
"""Gate for the allelic scan (spec §7). Runs the copied constrained scan on the depth
pseudobulks exactly as run_bootstrap_diagnostics.py built them (mean of count_matrix X
over each wnn_leiden_0.3 cluster, plus all cells), every bin kept, at the defaults and with
the chromosome global pool (pool=None; that run called recursive_getbp(seq, seq, ...)),
and compares with raw_accepted_breakpoints.tsv. This checks the copied recursion; the
allelic tab defaults to original pyEpi's genome-wide pool.

Run from the repo root:  python3 tests/validate_allelic_depth_reproduction.py
Writes validation/allelic/depth_reproduction.json; exits 1 on any difference.
"""
import json
import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import h5py
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from server import allelic_scan as sc  # noqa: E402
from server.config import Config  # noqa: E402
from server.integration import _read_cat  # noqa: E402


def main() -> int:
    import anndata as ad
    t0 = time.time()
    cfg = Config.from_env()
    h5mu = cfg.h5mu or os.path.join(cfg.integration_dir, 'data', 'integrated.h5mu')
    depth = ad.read_h5ad(cfg.h5ad)
    cells = depth.obs['cellID'].astype(str).to_numpy()
    with h5py.File(h5mu, 'r') as f:
        obs = f['obs/_index'][...].astype(str)
        lab = _read_cat(f['obs'], 'wnn_leiden_0.3').astype(str)
    of = dict(zip(obs, lab))
    cell_lab = np.array([of[c] for c in cells])
    groups = {f'cluster{c}': cell_lab == c for c in sorted(set(lab), key=int)}
    groups['all'] = np.ones(len(cells), bool)
    chrom_of = depth.var['seq'].astype(str).to_numpy()
    chroms = list(dict.fromkeys(chrom_of))
    tasks = []
    for g, mask in groups.items():
        x_pb = np.asarray(depth.X[mask].mean(axis=0)).ravel()
        for ch in chroms:
            seq = x_pb[chrom_of == ch]
            tasks.append((g, ch, seq, np.ones(len(seq), bool), None, None, sc.DEFAULTS, cfg.package_dir, None))
    workers = max(1, (os.cpu_count() or 4) - 2)
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as ex:
        results = list(ex.map(sc.scan_task, tasks))
    got = {}
    for g, ch, out in results:
        for b in out['breakpoints']:
            got[(g, ch, b['bin'])] = b
    ref = pd.read_csv(cfg.cluster_breakpoints, sep='\t')
    want = {(r.source, r.chromosome, int(r.absolute_bin)): r for r in ref.itertuples(index=False)}
    missing = sorted(set(want) - set(got))
    extra = sorted(set(got) - set(want))
    mismatched = []
    for key in sorted(set(want) & set(got)):
        w, b = want[key], got[key]
        if (int(w.recursion_depth) != b['recursion_depth'] or abs(float(w.p_local) - b['p_local']) > 1e-12
                or abs(float(w.p_global) - b['p_global']) > 1e-12 or abs(float(w.ad_score) - b['ad_score']) > 1e-6 * max(1.0, abs(b['ad_score']))):
            mismatched.append({'key': list(key), 'want': {'depth': int(w.recursion_depth), 'p_local': float(w.p_local), 'p_global': float(w.p_global), 'ad': float(w.ad_score)},
                               'got': {'depth': b['recursion_depth'], 'p_local': b['p_local'], 'p_global': b['p_global'], 'ad': b['ad_score']}})
    report = {'settings': sc.DEFAULTS, 'global_pool': 'chromosome', 'n_reference': len(want), 'n_reproduced': len(got),
              'missing': [list(k) for k in missing[:50]], 'n_missing': len(missing),
              'extra': [list(k) for k in extra[:50]], 'n_extra': len(extra),
              'mismatched': mismatched[:50], 'n_mismatched': len(mismatched),
              'seconds': round(time.time() - t0, 1), 'workers': workers,
              'reference': cfg.cluster_breakpoints, 'h5ad': cfg.h5ad, 'package_dir': cfg.package_dir}
    out_dir = os.path.join(ROOT, 'validation', 'allelic')
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'depth_reproduction.json'), 'w') as f:
        json.dump(report, f, indent=1)
    ok = not (missing or extra or mismatched)
    print(f"reference {len(want)} · reproduced {len(got)} · missing {len(missing)} · extra {len(extra)} · "
          f"mismatched {len(mismatched)} · {report['seconds']} s → {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
```

- [ ] **Step 6: Run the gate on real data**

Run: `python3 tests/validate_allelic_depth_reproduction.py`
Expected: `reference 814 · reproduced 814 · missing 0 · extra 0 · mismatched 0 · … → PASS`. It takes several minutes on the real data.

**If the gate fails, stop.** Report `BLOCKED` with the summary line and the first entries of `missing`, `extra` and `mismatched` from `validation/allelic/depth_reproduction.json`. Do not change the scan, the defaults or the tolerances to make it pass. The controller takes the result to the user.

- [ ] **Step 7: Commit**

```bash
git add server/allelic_scan.py tests/test_allelic_scan.py tests/validate_allelic_depth_reproduction.py validation/allelic/depth_reproduction.json
git commit -m "feat: constrained AD scan for allelic tracks, with depth-breakpoint reproduction gate

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Allelic data — load checks, pseudobulk sums, signals

**Files:**
- Create: `server/allelic_data.py`
- Create: `tests/allelic_fixtures.py`
- Create: `tests/test_allelic_data.py`
- Modify: `server/config.py` (add `allelic_dir`)

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces:
  - `allelic_data.SIGNALS = ('dbaf', 'abs_dbaf', 'shrunk_abs_dbaf', 'oriented_all')`.
  - `allelic_data.signal_values(ha, hb, signal, min_at, ha_all=None, hb_all=None) -> (values: float64 with NaN where skipped, kept: bool)`.
  - `allelic_data.check_inputs(allelic_dir, var_seq, var_start, var_end, cell_ids) -> list[str]` (problems).
  - `allelic_data.sum_by_source(allelic_dir, cell_ids, members, sources, batch=200, at_check=200) -> (ha (S, n) float64, hb (S, n) float64)`. Raises `AllelicLoadError` when AT ≠ HA + HB.
  - `allelic_data.AllelicData(allelic_dir, cache_dir, var_seq, var_start, var_end, cell_ids, members, sources)`, with:
    - attributes `.state` (`'idle' | 'building' | 'ready' | 'error'`), `.problems`, `.notes`, `.key`, `.sources`, `.chromosomes`;
    - methods `.start()`, `.wait(timeout)`, `.status() -> dict`, `.offset(chrom) -> int`, `.track(source, chrom, signal, min_at) -> dict(values, kept, at, ha, hb)`.
  - `Config.allelic_dir: str`.
  - `tests/allelic_fixtures.make_allelic(tmp_path, ...)`: a synthetic parquet folder and the matching var/cell/membership arrays (used by Task 3).

- [ ] **Step 1: Write the fixture helper**

`tests/allelic_fixtures.py`:

```python
"""Synthetic aligned HA/HB/AT parquet folders for allelic tests (same layout as
hapcnv_outs_from_updated_phasing/aligned_matrices: bins × barcodes, chr/start/end index)."""
import numpy as np
import pandas as pd

CHROMS = (('chr1', 120), ('chr2', 80))


def make_allelic(tmp_path, n_cells=8, seed=0, shuffle_columns=True, at_offset=0.0, coords_shift=0, drop_barcode=False):
    rng = np.random.default_rng(seed)
    var_seq = np.concatenate([[c] * n for c, n in CHROMS])
    var_start = np.concatenate([np.arange(n) * 100000 + 1 for _, n in CHROMS])
    var_end = var_start + 99999
    n_bins = len(var_seq)
    cell_ids = [f'BC{i:03d}-1' for i in range(n_cells)]
    lab = np.array([i % 2 for i in range(n_cells)])                   # cluster0 / cluster1 alternate
    members = {'cluster0': np.where(lab == 0)[0], 'cluster1': np.where(lab == 1)[0], 'all': np.arange(n_cells)}
    ha = rng.poisson(6, size=(n_bins, n_cells)).astype(np.float32)
    hb = rng.poisson(6, size=(n_bins, n_cells)).astype(np.float32)
    # cluster1, chr1 bins 60..119: haplotype B lost (HB small) → strong imbalance step at bin 60
    c1 = members['cluster1']
    hb[60:120][:, c1] = rng.poisson(1, size=(60, len(c1)))
    order = list(rng.permutation(n_cells)) if shuffle_columns else list(range(n_cells))
    cols = [cell_ids[i] for i in order]
    if drop_barcode:
        cols = cols[:-1]
        order = order[:-1]
    idx = pd.MultiIndex.from_arrays([var_seq, var_start + coords_shift, var_end + coords_shift], names=['chr', 'start', 'end'])
    d = tmp_path / 'aligned_matrices'
    d.mkdir(parents=True, exist_ok=True)
    for name, m in (('HA', ha), ('HB', hb), ('AT', ha + hb + at_offset)):
        pd.DataFrame(m[:, order], index=idx, columns=cols).to_parquet(d / f'{name}.parquet')
    return {'dir': str(d), 'var_seq': var_seq, 'var_start': var_start, 'var_end': var_end, 'cell_ids': cell_ids,
            'members': members, 'sources': ['cluster0', 'cluster1', 'all'], 'ha': ha.astype(np.float64), 'hb': hb.astype(np.float64)}
```

- [ ] **Step 2: Write the failing tests**

`tests/test_allelic_data.py`:

```python
import numpy as np
import pytest

from server import allelic_data as ald
from allelic_fixtures import make_allelic


def test_signal_formulas_and_kept_bins():
    ha = np.array([30., 10., 5., 0., 4.]); hb = np.array([10., 10., 1., 0., 8.])
    v, k = ald.signal_values(ha, hb, 'dbaf', 10)
    assert k.tolist() == [True, True, False, False, True]
    assert v[0] == pytest.approx(-0.25) and v[1] == 0.0 and np.isnan(v[2]) and np.isnan(v[3]) and v[4] == pytest.approx(1 / 6)
    v, _ = ald.signal_values(ha, hb, 'abs_dbaf', 10)
    assert v[0] == pytest.approx(0.25)
    v, _ = ald.signal_values(ha, hb, 'shrunk_abs_dbaf', 10)
    assert v[0] == pytest.approx(20 / (2 * 42))
    ha_all = np.array([50., 20., 9., 1., 7.]); hb_all = np.array([20., 20., 3., 1., 2.])
    v, k = ald.signal_values(ha, hb, 'oriented_all', 10, ha_all, hb_all)
    assert k.tolist() == [True, False, False, False, True]          # bin 1: all-cells HA == HB → skipped
    assert v[0] == pytest.approx(0.25) and v[4] == pytest.approx(-1 / 6)
    with pytest.raises(ValueError):
        ald.signal_values(ha, hb, 'baf', 10)


def test_check_inputs_good(tmp_path):
    f = make_allelic(tmp_path)
    assert ald.check_inputs(f['dir'], f['var_seq'], f['var_start'], f['var_end'], f['cell_ids']) == []


def test_check_inputs_rejects_coordinates_and_barcodes(tmp_path):
    f = make_allelic(tmp_path / 'a', coords_shift=1)
    p = ald.check_inputs(f['dir'], f['var_seq'], f['var_start'], f['var_end'], f['cell_ids'])
    assert any('coordinates' in x for x in p)
    g = make_allelic(tmp_path / 'b', drop_barcode=True)
    p = ald.check_inputs(g['dir'], g['var_seq'], g['var_start'], g['var_end'], g['cell_ids'])
    assert any('barcodes' in x for x in p)
    p = ald.check_inputs(str(tmp_path / 'nope'), g['var_seq'], g['var_start'], g['var_end'], g['cell_ids'])
    assert any('not found' in x for x in p)


def test_sums_follow_barcodes_not_column_order(tmp_path):
    f = make_allelic(tmp_path, shuffle_columns=True)
    ha, hb = ald.sum_by_source(f['dir'], f['cell_ids'], f['members'], f['sources'], batch=3, at_check=5)
    for i, s in enumerate(f['sources']):
        m = f['members'][s]
        assert np.allclose(ha[i], f['ha'][:, m].sum(1)) and np.allclose(hb[i], f['hb'][:, m].sum(1))


def test_sums_reject_at_mismatch(tmp_path):
    f = make_allelic(tmp_path, at_offset=1.0)
    with pytest.raises(ald.AllelicLoadError, match='AT'):
        ald.sum_by_source(f['dir'], f['cell_ids'], f['members'], f['sources'])


def _data(f, cache):
    return ald.AllelicData(f['dir'], str(cache), f['var_seq'], f['var_start'], f['var_end'], f['cell_ids'], f['members'], f['sources'])


def test_allelic_data_builds_then_loads_from_cache(tmp_path, monkeypatch):
    f = make_allelic(tmp_path)
    d = _data(f, tmp_path / 'cache')
    assert d.status()['state'] == 'idle'
    d.start(); d.wait(30)
    assert d.state == 'ready', d.problems
    t = d.track('cluster1', 'chr1', 'abs_dbaf', 10)
    assert len(t['values']) == 120 and t['kept'].dtype == bool
    assert np.nanmean(t['values'][60:]) > 0.3 > np.nanmean(t['values'][:60])
    assert d.offset('chr2') == 120 and d.chromosomes == ['chr1', 'chr2']
    monkeypatch.setattr(ald, 'sum_by_source', lambda *a, **k: (_ for _ in ()).throw(AssertionError('cache miss')))
    d2 = _data(f, tmp_path / 'cache')
    d2.start(); d2.wait(30)
    assert d2.state == 'ready' and d2.key == d.key


def test_allelic_data_reports_missing_folder(tmp_path):
    f = make_allelic(tmp_path)
    d = ald.AllelicData(str(tmp_path / 'missing'), str(tmp_path / 'cache'), f['var_seq'], f['var_start'], f['var_end'], f['cell_ids'], f['members'], f['sources'])
    d.start(); d.wait(30)
    s = d.status()
    assert s['state'] == 'error' and not s['available'] and 'PYEPI_SCIENTIFIC_ALLELIC' in ' '.join(s['problems'])
```

- [ ] **Step 3: Run them to check they fail**

Run: `python3 -m pytest -q tests/test_allelic_data.py`
Expected: FAIL, `ImportError: cannot import name 'allelic_data'`.

- [ ] **Step 4: Write the module**

`server/allelic_data.py`:

```python
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
```

- [ ] **Step 5: Add the config field**

In `server/config.py`:
- Add `allelic_dir: str = ""` after `cluster_breakpoints: str = ""` in the dataclass.
- In `from_env`, pass the argument after `cluster_breakpoints=...`:

```python
            allelic_dir=os.environ.get("PYEPI_SCIENTIFIC_ALLELIC", os.path.normpath(os.path.join(
                root, "..", "..", "..", "haplotype_aware_cnv", "hapcnv_outs_from_updated_phasing", "aligned_matrices"))),
```

- [ ] **Step 6: Run the tests to check they pass**

Run: `python3 -m pytest -q tests/test_allelic_data.py`
Expected: 7 passed.

Also run: `python3 -c "from server.config import Config; import os; print(os.path.isfile(os.path.join(Config.from_env().allelic_dir, 'HA.parquet')))"`
Expected: `True`.

- [ ] **Step 7: Commit**

```bash
git add server/allelic_data.py server/config.py tests/allelic_fixtures.py tests/test_allelic_data.py
git commit -m "feat: allelic pseudobulk sums, load checks and the four allelic signals

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Allelic routes and the segmentation job

**Files:**
- Create: `server/allelic_api.py`
- Create: `tests/test_allelic_api.py`
- Modify: `server/app.py` (register the routes)

**Interfaces:**
- Consumes:
  - Task 1: `allelic_scan.DEFAULTS` and `scan_task`.
  - Task 2: `allelic_data.SIGNALS`, `AllelicData` and `signal_values`.
- Produces (HTTP):
  - **`GET /api/allelic/status`** → `{available, state, problems, notes, allelic_dir, sources}`. The first call starts the build.
  - **`GET /api/allelic/tracks?chrom=&sources=a,b&signal=&min_at=`** → `{chrom, signal, min_at, n, tracks: [{source, values: (float|null)[], at: int[]}]}`.
    - 503 while the data is not ready.
    - 404 for an unknown chromosome or source.
    - 422 for a bad signal or min_at.
  - **`POST /api/allelic/segment`** with body `{signal, settings}` → job view. 422 detail names the field: `"<field>: <message>"`.
  - **`GET /api/allelic/segment/{id}`** → job view, 404 if unknown.
  - **Job view:** `{id, status: 'queued'|'running'|'done'|'error', progress: {done, total, source, chromosome}|null, error, signal, settings, result|null}`.
    - `result = {breakpoints, audit, segments}`.
    - Breakpoint rows: `source, chromosome, absolute_bin, bp_genomic_start, bp_genomic_end, gap_start, gap_end, recursion_depth, ad_score, p_local, p_global`.
    - Segment rows: `source, chromosome, bin_start, bin_end, n_bins, kept_bins, mean, median, at_sum, ha_sum, hb_sum`.
    - Audit rows: the scan's audit columns with source and chromosome filled in.
  - **Python:** `allelic_api.normalize(signal, settings) -> dict` (full settings including `global_pool`, raises `ValueError('<field>: …')`), `allelic_api.DEFAULT_SETTINGS` and `allelic_api.register_allelic(app, get_data, get_executor, package_dir, cache_dir) -> AllelicManager`.

- [ ] **Step 1: Write the failing tests**

`tests/test_allelic_api.py`:

```python
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from fastapi import FastAPI
from fastapi.testclient import TestClient

from server import allelic_api as aa
from server import allelic_data as ald
from server import allelic_scan as sc
from allelic_fixtures import make_allelic

FAST = {'min_at': 5, 'k': 2, 'n_permutations': 100, 'alpha': 0.05, 'min_child_bins': 10, 'min_terminal_bins': 10, 'global_pool': 'genome'}


def _client(tmp_path, f=None):
    f = f or make_allelic(tmp_path)
    d = ald.AllelicData(f['dir'], str(tmp_path / 'cache'), f['var_seq'], f['var_start'], f['var_end'], f['cell_ids'], f['members'], f['sources'])
    ex = ThreadPoolExecutor(2)
    app = FastAPI()
    mgr = aa.register_allelic(app, lambda: d, lambda: ex, None, str(tmp_path / 'cache'))
    c = TestClient(app)
    c.get('/api/allelic/status'); d.wait(30)
    return c, d, mgr, f


def _wait(c, jid, timeout=120):
    t0 = time.time()
    while time.time() - t0 < timeout:
        v = c.get(f'/api/allelic/segment/{jid}').json()
        if v['status'] in ('done', 'error'):
            return v
        time.sleep(0.05)
    raise AssertionError('job did not finish')


def test_normalize_defaults_and_errors():
    s = aa.normalize('abs_dbaf', {})
    assert s == {'min_at': 10, 'k': 6, 'n_permutations': 1000, 'alpha': 0.001, 'min_child_bins': 20, 'min_terminal_bins': 50, 'global_pool': 'genome'}
    for bad, field in [({'alpha': 0}, 'alpha'), ({'alpha': 0.2}, 'alpha'), ({'k': 2.5}, 'k'), ({'min_at': 0}, 'min_at'),
                       ({'n_permutations': 'x'}, 'n_permutations'), ({'zzz': 1}, 'zzz'), ({'global_pool': 'cell'}, 'global_pool')]:
        try:
            aa.normalize('abs_dbaf', bad)
        except ValueError as e:
            assert str(e).startswith(f'{field}:')
        else:
            raise AssertionError(bad)
    try:
        aa.normalize('baf', {})
    except ValueError as e:
        assert str(e).startswith('signal:')


def test_status_and_tracks(tmp_path):
    c, d, _, f = _client(tmp_path)
    s = c.get('/api/allelic/status').json()
    assert s['state'] == 'ready' and s['available'] and s['sources'] == f['sources']
    t = c.get('/api/allelic/tracks', params={'chrom': 'chr1', 'sources': 'cluster1,all', 'signal': 'dbaf', 'min_at': 5}).json()
    assert t['n'] == 120 and [x['source'] for x in t['tracks']] == ['cluster1', 'all']
    v = t['tracks'][0]['values']
    assert len(v) == 120 and all(x is None or -0.5 <= x <= 0.5 for x in v)
    assert c.get('/api/allelic/tracks', params={'chrom': 'chrZ', 'sources': 'all'}).status_code == 404
    assert c.get('/api/allelic/tracks', params={'chrom': 'chr1', 'sources': 'cluster9'}).status_code == 404
    assert c.get('/api/allelic/tracks', params={'chrom': 'chr1', 'sources': 'all', 'signal': 'baf'}).status_code == 422


def test_segment_runs_and_reuses(tmp_path):
    c, d, _, f = _client(tmp_path)
    v = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': FAST}).json()
    done = _wait(c, v['id'])
    assert done['status'] == 'done' and done['settings'] == FAST and done['signal'] == 'abs_dbaf'
    bps = done['result']['breakpoints']
    c1 = [b for b in bps if b['source'] == 'cluster1' and b['chromosome'] == 'chr1']
    assert any(abs(b['absolute_bin'] - 60) <= 2 for b in c1)
    b = c1[0]
    assert b['bp_genomic_start'] == int(f['var_start'][b['absolute_bin']])
    segs = [s for s in done['result']['segments'] if s['source'] == 'cluster0' and s['chromosome'] == 'chr2']
    assert segs[0]['bin_start'] == 0 and segs[-1]['bin_end'] == 79
    assert {r['source'] for r in done['result']['audit']} <= set(f['sources'])
    again = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': FAST}).json()
    assert again['id'] == v['id'] and again['status'] == 'done'
    other = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': {**FAST, 'alpha': 0.01}}).json()
    assert other['id'] != v['id']
    _wait(c, other['id'])


def test_global_pool_reaches_the_scan(tmp_path, monkeypatch):
    c, d, _, f = _client(tmp_path)
    pools = {}
    real = sc.scan_task

    def spy(task):
        pools[(task[-1] is None, task[0], task[1])] = None if task[-1] is None else len(task[-1])
        return real(task)
    monkeypatch.setattr(sc, 'scan_task', spy)
    _wait(c, c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': FAST}).json()['id'])
    kept_total = sum(int(d.track('cluster0', ch, 'abs_dbaf', 5)['kept'].sum()) for ch in d.chromosomes)
    assert pools[(False, 'cluster0', 'chr1')] == kept_total == pools[(False, 'cluster0', 'chr2')]
    _wait(c, c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': {**FAST, 'global_pool': 'chromosome'}}).json()['id'])
    assert (True, 'cluster0', 'chr1') in pools


def test_bad_requests(tmp_path):
    c, *_ = _client(tmp_path)
    r = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': {'alpha': 0}})
    assert r.status_code == 422 and r.json()['detail'].startswith('alpha:')
    assert c.post('/api/allelic/segment', json={'signal': 'nope', 'settings': {}}).status_code == 422
    assert c.get('/api/allelic/segment/abc').status_code == 404


def test_cache_survives_a_new_manager(tmp_path):
    f = make_allelic(tmp_path)
    c, *_ = _client(tmp_path, f)
    jid = c.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()['id']
    first = _wait(c, jid)
    c2, *_ = _client(tmp_path, f)
    v = c2.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()
    assert v['id'] == jid and v['status'] == 'done' and v['result'] == first['result']


def test_queue_replaces_the_waiting_request(tmp_path, monkeypatch):
    c, *_ = _client(tmp_path)
    gate = threading.Event()
    real = sc.scan_task

    def slow(task):
        gate.wait(30)
        return real(task)
    monkeypatch.setattr(sc, 'scan_task', slow)
    a = c.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()
    b = c.post('/api/allelic/segment', json={'signal': 'abs_dbaf', 'settings': FAST}).json()
    cc = c.post('/api/allelic/segment', json={'signal': 'shrunk_abs_dbaf', 'settings': FAST}).json()
    assert a['status'] == 'running' and b['status'] == 'queued' and cc['status'] == 'queued'
    vb = c.get(f"/api/allelic/segment/{b['id']}").json()
    assert vb['status'] == 'error' and vb['error'] == 'replaced by a newer request'
    gate.set()
    assert _wait(c, a['id'])['status'] == 'done'
    assert _wait(c, cc['id'])['status'] == 'done'


def test_failed_task_fails_the_job_and_is_not_cached(tmp_path, monkeypatch):
    c, *_ = _client(tmp_path)

    def boom(task):
        raise RuntimeError('kaput')
    monkeypatch.setattr(sc, 'scan_task', boom)
    v = _wait(c, c.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()['id'])
    assert v['status'] == 'error' and 'kaput' in v['error'] and v['result'] is None
    monkeypatch.undo()
    again = c.post('/api/allelic/segment', json={'signal': 'dbaf', 'settings': FAST}).json()
    assert again['status'] in ('running', 'queued')
    assert _wait(c, again['id'])['status'] == 'done'
```

- [ ] **Step 2: Run them to check they fail**

Run: `python3 -m pytest -q tests/test_allelic_api.py`
Expected: FAIL, `ImportError: cannot import name 'allelic_api'`.

- [ ] **Step 3: Write the module**

`server/allelic_api.py`:

```python
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
```

- [ ] **Step 4: Run the tests to check they pass**

Run: `python3 -m pytest -q tests/test_allelic_api.py`
Expected: 8 passed.

- [ ] **Step 5: Register the routes in the app**

In `server/app.py`, after the `register_cn_recall(...)` line at the end of the file, add:

```python
from .allelic_api import register_allelic
from .allelic_data import AllelicData
from .integration import SOURCES as _INTEGRATION_SOURCES

allelic_data: Optional[AllelicData] = None


def _allelic() -> Optional[AllelicData]:
    global allelic_data
    if allelic_data is None and integration is not None and integration.available and store is not None:
        allelic_data = AllelicData(cfg.allelic_dir, cfg.cache_dir, store.var_seq, store.var_start, store.var_end,
                                   list(integration.cell_ids), integration.members, _INTEGRATION_SOURCES)
    return allelic_data


register_allelic(app, _allelic, lambda: jobs.background if jobs else None, cfg.package_dir, cfg.cache_dir)
```

- [ ] **Step 6: Check the routes on a private server**

```bash
python3 run.py --port 8799 > /tmp/allelic_8799.log 2>&1 &
echo $! > /tmp/allelic_8799.pid
for i in $(seq 1 120); do curl -sf http://127.0.0.1:8799/api/allelic/status | grep -q '"state":"ready"' && break; sleep 2; done
curl -s http://127.0.0.1:8799/api/allelic/status
curl -s 'http://127.0.0.1:8799/api/allelic/tracks?chrom=chr8&sources=cluster2&signal=abs_dbaf&min_at=10' | head -c 300; echo
```

Expected:
- the status shows `"state":"ready"` and `"problems":[]`;
- the tracks response has `"n":` equal to chr8's bin count, with numbers and nulls in `values`.

Then start a run and poll it until done:

```bash
JID=$(curl -s -X POST -H 'content-type: application/json' -d '{"signal":"abs_dbaf","settings":{}}' http://127.0.0.1:8799/api/allelic/segment | python3 -c 'import sys,json;print(json.load(sys.stdin)["id"])')
for i in $(seq 1 300); do S=$(curl -s http://127.0.0.1:8799/api/allelic/segment/$JID | python3 -c 'import sys,json;d=json.load(sys.stdin);print(d["status"], d["progress"])'); echo "$S"; case "$S" in done*|error*) break;; esac; sleep 2; done
kill $(cat /tmp/allelic_8799.pid)
```

Expected: the job ends `done`. Record the elapsed time and the number of breakpoints per source in the report. Stop only the 8799 process.

- [ ] **Step 7: Run the whole Python suite**

Run: `python3 -m pytest -q`
Expected: all pass (220 earlier + 25 new).

- [ ] **Step 8: Commit**

```bash
git add server/allelic_api.py server/app.py tests/test_allelic_api.py
git commit -m "feat: allelic status, tracks and cached segmentation job routes

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Web logic, API client, store and session

**Files:**
- Create: `web/src/allelic.ts`
- Create: `web/src/allelicStore.ts`
- Create: `tests/allelic_helpers.ts`
- Create: `tests/allelic_store.ts`
- Modify: `web/src/store.ts` (session serialize/restore)

**Interfaces:**
- Consumes: Task 3 HTTP routes and their shapes; `web/src/types.ts` (`ClusterBreakpoint`, `ClusterSegment`); `web/src/cnRecall.ts` `downloadText` (used in Task 5).
- Produces:
  - **From `allelic.ts`:**
    - types: `SignalId`, `AdSettings`, `AllelicStatus`, `AllelicTracks`, `AllelicBreakpoint`, `AllelicSegment`, `AllelicJob`, `AllelicView`, `LocalDepthSeg`, `AgreementRow`, `DepthSegRow`, `AllelicSegRow`;
    - constants: `SIGNALS`, `DEFAULT_SETTINGS`, `SETTING_FIELDS`, `DEFAULT_VIEW`, `COMPARISON_DEFINITIONS`;
    - functions: `signalLabel`, `settingsProblems`, `sameRun`, `bpDistance`, `matchBreakpoints`, `agreementRows`, `depthSegmentRows`, `callShares`, `sharesText`, `allelicSegmentRows`, `runHeader`, `comparisonTsv`, `allelicFromSession`;
    - the `allelicApi` object.
  - **From `allelicStore.ts`:** `useAllelic` (state `AllelicView & {status, statusError, tracks, tracksKey, tracksError, shown, active, runError, notice}`, actions `loadStatus`, `setView`, `ensureTracks`, `run`, `serialize`, `restore`) and `timing.pollMs`.
  - **Session:** saved under the key `allelic`.

- [ ] **Step 1: Write the failing helper tests**

`tests/allelic_helpers.ts`:

```ts
// tests/allelic_helpers.ts
// Run with: web/node_modules/.bin/esbuild tests/allelic_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/allelic_helpers.mjs
//           node tests/_bundle/allelic_helpers.mjs
import { DEFAULT_SETTINGS, DEFAULT_VIEW, agreementRows, allelicFromSession, allelicSegmentRows, bpDistance, callShares, comparisonTsv, depthSegmentRows, matchBreakpoints, sameRun, settingsProblems, sharesText, type AllelicBreakpoint, type AllelicSegment } from '../web/src/allelic'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}
const abp = (source: string, chromosome: string, bin: number, gap: [number, number] | null = null): AllelicBreakpoint => ({ source, chromosome, absolute_bin: bin, bp_genomic_start: 0, bp_genomic_end: 0, gap_start: gap ? gap[0] : null, gap_end: gap ? gap[1] : null, recursion_depth: 0, ad_score: 1, p_local: 0.001, p_global: 0.001 })
const dbp = (source: string, chromosome: string, bin: number) => ({ source, chromosome, absolute_bin: bin, bp_genomic_start: 0, bp_genomic_end: 0, recursion_depth: 0, ad_score: 1, p_local: 0.001, p_global: 0.001, n_cells: 1 })

eq(settingsProblems(DEFAULT_SETTINGS), [], 'defaults are valid')
eq(settingsProblems({ ...DEFAULT_SETTINGS, alpha: 0, k: 2.5 }), ['k: must be a whole number', 'alpha: must be > 0 and ≤ 0.1'], 'settings problems name the field')
eq(sameRun({ signal: 'abs_dbaf', settings: DEFAULT_SETTINGS }, 'abs_dbaf', { ...DEFAULT_SETTINGS }), true, 'same run')
eq(sameRun({ signal: 'abs_dbaf', settings: DEFAULT_SETTINGS }, 'dbaf', DEFAULT_SETTINGS), false, 'different signal')
eq(sameRun(null, 'dbaf', DEFAULT_SETTINGS), false, 'no run')
eq(sameRun({ signal: 'dbaf', settings: DEFAULT_SETTINGS }, 'dbaf', { ...DEFAULT_SETTINGS, global_pool: 'chromosome' }), false, 'pool is part of the run')
eq(settingsProblems({ ...DEFAULT_SETTINGS, global_pool: 'cell' as any }), ['global_pool: must be genome or chromosome'], 'pool checked')

eq(bpDistance(50, { bin: 53, gap_start: null, gap_end: null }), 3, 'plain distance')
eq(bpDistance(47, { bin: 50, gap_start: 45, gap_end: 50 }), 0, 'inside a gap is 0')
eq(bpDistance(41, { bin: 50, gap_start: 45, gap_end: 50 }), 4, 'left of a gap: distance to its start')
eq(bpDistance(53, { bin: 50, gap_start: 45, gap_end: 50 }), 3, 'right of a gap: distance to its end')

const m = matchBreakpoints([10, 12, 40], [{ bin: 11, gap_start: null, gap_end: null }, { bin: 44, gap_start: null, gap_end: null }, { bin: 90, gap_start: null, gap_end: null }], 3)
eq(m.pairs.map((p) => [p.d, p.a, p.dist]), [[0, 0, 1]], 'one-to-one: tie on distance goes to the lower depth bin')
eq([m.depthOnly, m.allelicOnly], [[1, 2], [1, 2]], 'unmatched on both sides')
const g = matchBreakpoints([47, 60], [{ bin: 50, gap_start: 45, gap_end: 50 }, { bin: 58, gap_start: null, gap_end: null }], 3)
eq(g.pairs.map((p) => [p.d, p.a, p.dist]), [[0, 0, 0], [1, 1, 2]], 'gap match at distance 0, then nearest')

const depthBps = [dbp('cluster1', 'chr1', 10), dbp('cluster1', 'chr1', 40), dbp('cluster1', 'chr2', 5), dbp('cluster2', 'chr1', 7)]
const allBps = [abp('cluster1', 'chr1', 11), abp('cluster1', 'chr2', 30), abp('cluster1', 'chr2', 6)]
eq(agreementRows(['cluster1', 'cluster2'], ['chr1', 'chr2'], depthBps, allBps, 3), [
  { source: 'cluster1', depth: 3, allelic: 3, matched: 2, depthOnly: 1, allelicOnly: 1, medianOffset: 1 },
  { source: 'cluster2', depth: 1, allelic: 0, matched: 0, depthOnly: 1, allelicOnly: 0, medianOffset: null }], 'agreement rows over chromosomes')
eq(agreementRows(['cluster1'], ['chr2'], depthBps, allBps, 3)[0].matched, 1, 'agreement rows for one chromosome')

const segs = [{ group: 'cluster1', s: 0, e: 20, cn_call: 1 }, { group: 'cluster1', s: 20, e: 60, cn_call: 1.5 }]
const values: Record<string, (number | null)[]> = { cluster1: Array.from({ length: 60 }, (_, i) => (i < 20 ? (i % 2 ? 0.1 : null) : 0.3)) }
eq(depthSegmentRows(segs, values, [abp('cluster1', 'chr1', 40), abp('cluster1', 'chr1', 22)], 'chr1', 3), [
  { source: 'cluster1', bin_start: 0, bin_end: 19, cn_call: 1, kept: 10, mean: 0.1, sd: 0, allelicInside: 0 },
  { source: 'cluster1', bin_start: 20, bin_end: 59, cn_call: 1.5, kept: 40, mean: 0.3, sd: 0, allelicInside: 1 }], 'depth segments: kept, mean, sd, allelic bps more than N from both edges')

eq(callShares({ bin_start: 10, bin_end: 29 }, segs), [{ call: 1, frac: 0.5 }, { call: 1.5, frac: 0.5 }], 'shares: equal split, smaller call first')
eq(sharesText([{ call: 1, frac: 0.8 }, { call: 1.5, frac: 0.2 }]), '1.0: 80% · 1.5: 20%', 'shares text')
eq(sharesText([]), 'no depth call', 'no overlap')
const aseg: AllelicSegment = { source: 'cluster1', chromosome: 'chr1', bin_start: 0, bin_end: 29, n_bins: 30, kept_bins: 25, mean: 0.12, median: 0.1, at_sum: 1, ha_sum: 0, hb_sum: 1 }
eq(allelicSegmentRows([aseg, { ...aseg, chromosome: 'chr2' }], segs, ['cluster1'], 'chr1'), [{ source: 'cluster1', bin_start: 0, bin_end: 29, kept: 25, mean: 0.12, spanned: '1.0: 67% · 1.5: 33%' }], 'allelic segment rows for the chosen chromosome')

const tsv = comparisonTsv(['signal: abs_dbaf'], agreementRows(['cluster2'], ['chr1'], depthBps, allBps, 3), [], [])
eq(tsv.split('\n').slice(0, 4), ['# signal: abs_dbaf', '# table 1: breakpoint agreement', 'source\tdepth_bps\tallelic_bps\tmatched\tdepth_only\tallelic_only\tmedian_offset', 'cluster2\t1\t0\t0\t1\t0\tNA'], 'TSV: header, table 1')
eq(tsv.includes('# table 2: depth segments seen by the allelic signal') && tsv.includes('# table 3: allelic segments seen by depth'), true, 'TSV: all three tables')

eq(allelicFromSession(undefined), { view: DEFAULT_VIEW, notice: null }, 'no saved state: defaults')
const saved = { ...DEFAULT_VIEW, signal: 'dbaf', chrom: 'chr3', sources: ['cluster4'], depthSource: 'original', matchN: 5, scope: 'genome', settings: { ...DEFAULT_SETTINGS, alpha: 0.01 } }
eq(allelicFromSession(JSON.parse(JSON.stringify(saved))), { view: saved, notice: null }, 'round trip')
const junk = allelicFromSession({ signal: 'baf', settings: { alpha: 7 }, chrom: 3, matchN: -1 })
eq(junk.view, DEFAULT_VIEW, 'junk: defaults')
eq(junk.notice?.startsWith('Saved allelic settings could not be applied: '), true, 'junk: notice')

if (failures) { console.error(`${failures} failure(s)`); process.exit(1) } else console.log('all allelic helper checks passed')
```

- [ ] **Step 2: Run them to check they fail**

Run: `web/node_modules/.bin/esbuild tests/allelic_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/allelic_helpers.mjs && node tests/_bundle/allelic_helpers.mjs`
Expected: the esbuild error `Could not resolve "../web/src/allelic"`.

- [ ] **Step 3: Write `web/src/allelic.ts`**

```ts
// Allelic segmentation: types, API, settings, comparison tables.
// Spec: docs/superpowers/specs/2026-10-28-allelic-segmentation-design.md
import type { ClusterBreakpoint } from './types'

export type SignalId = 'dbaf' | 'abs_dbaf' | 'shrunk_abs_dbaf' | 'oriented_all'
export const SIGNALS: { id: SignalId; label: string; signed: boolean; formula: string }[] = [
  { id: 'dbaf', label: 'signed dBAF', signed: true, formula: '(HB − HA) / (2·AT)' },
  { id: 'abs_dbaf', label: '|dBAF|', signed: false, formula: '|HB − HA| / (2·AT)' },
  { id: 'shrunk_abs_dbaf', label: 'shrunk |dBAF|', signed: false, formula: '|HB − HA| / (2·(AT + 2))' },
  { id: 'oriented_all', label: 'oriented to all cells', signed: true, formula: 'sign(ΣHB_all − ΣHA_all) · (HB − HA) / (2·AT)' },
]
export const signalLabel = (id: string) => SIGNALS.find((s) => s.id === id)?.label ?? id

export type GlobalPool = 'genome' | 'chromosome'
export interface AdSettings { min_at: number; k: number; n_permutations: number; alpha: number; min_child_bins: number; min_terminal_bins: number; global_pool: GlobalPool }
export const DEFAULT_SETTINGS: AdSettings = { min_at: 10, k: 6, n_permutations: 1000, alpha: 0.001, min_child_bins: 20, min_terminal_bins: 50, global_pool: 'genome' }
export const POOL_NOTE = 'The depth breakpoints were made with the chromosome pool; original pyEpi uses the genome pool.'
export const SETTING_FIELDS: { key: Exclude<keyof AdSettings, 'global_pool'>; label: string; min: number; max: number; integer: boolean; minExclusive?: boolean }[] = [
  { key: 'min_at', label: 'min AT', min: 1, max: 1000, integer: true },
  { key: 'k', label: 'k (max depth)', min: 1, max: 8, integer: true },
  { key: 'n_permutations', label: 'permutations', min: 100, max: 5000, integer: true },
  { key: 'alpha', label: 'alpha', min: 0, max: 0.1, integer: false, minExclusive: true },
  { key: 'min_child_bins', label: 'min child (kept bins)', min: 1, max: 500, integer: true },
  { key: 'min_terminal_bins', label: 'min from chromosome end (kept bins)', min: 0, max: 500, integer: true },
]
export function settingsProblems(s: AdSettings): string[] {
  const out: string[] = []
  for (const f of SETTING_FIELDS) {
    const v = (s as any)[f.key]
    if (typeof v !== 'number' || !Number.isFinite(v)) { out.push(`${f.key}: not a number`); continue }
    if (f.integer && !Number.isInteger(v)) out.push(`${f.key}: must be a whole number`)
    if ((f.minExclusive ? v <= f.min : v < f.min) || v > f.max) out.push(`${f.key}: must be ${f.minExclusive ? '>' : '≥'} ${f.min} and ≤ ${f.max}`)
  }
  if (s.global_pool !== 'genome' && s.global_pool !== 'chromosome') out.push('global_pool: must be genome or chromosome')
  return out
}
export const sameRun = (run: { signal: string; settings: AdSettings } | null, signal: string, s: AdSettings) =>
  !!run && run.signal === signal && run.settings.global_pool === s.global_pool && SETTING_FIELDS.every((f) => run.settings[f.key] === s[f.key])

export interface AllelicStatus { available: boolean; state: 'idle' | 'building' | 'ready' | 'error'; problems: string[]; notes: string[]; allelic_dir: string; sources: string[] }
export interface AllelicTracks { chrom: string; signal: SignalId; min_at: number; n: number; tracks: { source: string; values: (number | null)[]; at: number[] }[] }
export interface AllelicBreakpoint { source: string; chromosome: string; absolute_bin: number; bp_genomic_start: number; bp_genomic_end: number; gap_start: number | null; gap_end: number | null; recursion_depth: number; ad_score: number; p_local: number; p_global: number }
export interface AllelicSegment { source: string; chromosome: string; bin_start: number; bin_end: number; n_bins: number; kept_bins: number; mean: number | null; median: number | null; at_sum: number; ha_sum: number; hb_sum: number }
export interface AllelicJob {
  id: string; status: 'queued' | 'running' | 'done' | 'error'; error: string | null; signal: SignalId; settings: AdSettings
  progress: { done: number; total: number; source: string | null; chromosome: string | null } | null
  result: { breakpoints: AllelicBreakpoint[]; audit: any[]; segments: AllelicSegment[] } | null
}

async function json<T>(r: Response): Promise<T> {
  if (!r.ok) { let msg = `${r.status}`; try { const b = await r.json(); msg = typeof b.detail === 'string' ? b.detail : JSON.stringify(b.detail ?? b) } catch { /* keep status */ } throw Object.assign(new Error(msg), { status: r.status }) }
  return r.json() as Promise<T>
}
const enc = encodeURIComponent
export const allelicApi = {
  status: () => fetch('/api/allelic/status').then((r) => json<AllelicStatus>(r)),
  tracks: (chrom: string, sources: string[], signal: SignalId, minAt: number) =>
    fetch(`/api/allelic/tracks?chrom=${enc(chrom)}&sources=${enc(sources.join(','))}&signal=${signal}&min_at=${minAt}`).then((r) => json<AllelicTracks>(r)),
  start: (signal: SignalId, settings: AdSettings) =>
    fetch('/api/allelic/segment', { method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify({ signal, settings }) }).then((r) => json<AllelicJob>(r)),
  job: (id: string) => fetch(`/api/allelic/segment/${enc(id)}`).then((r) => json<AllelicJob>(r)),
}

// ---------------- comparison (spec §5); all bins chromosome-local
export interface AllelicBp { bin: number; gap_start: number | null; gap_end: number | null }
export function bpDistance(depthBin: number, a: AllelicBp): number {
  if (a.gap_start !== null && a.gap_end !== null) {
    if (depthBin >= a.gap_start && depthBin <= a.gap_end) return 0
    return depthBin < a.gap_start ? a.gap_start - depthBin : depthBin - a.gap_end
  }
  return Math.abs(depthBin - a.bin)
}
export function matchBreakpoints(depth: number[], allelic: AllelicBp[], N: number) {
  const cand: { d: number; a: number; dist: number }[] = []
  depth.forEach((db, i) => allelic.forEach((ab, j) => { const dist = bpDistance(db, ab); if (dist <= N) cand.push({ d: i, a: j, dist }) }))
  cand.sort((x, y) => x.dist - y.dist || depth[x.d] - depth[y.d] || allelic[x.a].bin - allelic[y.a].bin)
  const usedD = new Set<number>(), usedA = new Set<number>(), pairs: typeof cand = []
  for (const c of cand) if (!usedD.has(c.d) && !usedA.has(c.a)) { usedD.add(c.d); usedA.add(c.a); pairs.push(c) }
  return { pairs, depthOnly: depth.map((_, i) => i).filter((i) => !usedD.has(i)), allelicOnly: allelic.map((_, j) => j).filter((j) => !usedA.has(j)) }
}
const median = (v: number[]) => { const s = [...v].sort((a, b) => a - b), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2 }

export interface AgreementRow { source: string; depth: number; allelic: number; matched: number; depthOnly: number; allelicOnly: number; medianOffset: number | null }
export function agreementRows(sources: string[], chroms: string[], depthBps: ClusterBreakpoint[], allelicBps: AllelicBreakpoint[], N: number): AgreementRow[] {
  return sources.map((source) => {
    const row: AgreementRow = { source, depth: 0, allelic: 0, matched: 0, depthOnly: 0, allelicOnly: 0, medianOffset: null }
    const offs: number[] = []
    for (const chrom of chroms) {
      const d = depthBps.filter((b) => b.source === source && b.chromosome === chrom).map((b) => b.absolute_bin)
      const a = allelicBps.filter((b) => b.source === source && b.chromosome === chrom).map((b) => ({ bin: b.absolute_bin, gap_start: b.gap_start, gap_end: b.gap_end }))
      const m = matchBreakpoints(d, a, N)
      row.depth += d.length; row.allelic += a.length; row.matched += m.pairs.length; row.depthOnly += m.depthOnly.length; row.allelicOnly += m.allelicOnly.length
      offs.push(...m.pairs.map((p) => p.dist))
    }
    row.medianOffset = offs.length ? median(offs) : null
    return row
  })
}

/** Depth cluster segment, chromosome-local half-open [s, e) as returned by clusterSegmentsLocal. */
export interface LocalDepthSeg { group: string; s: number; e: number; cn_call: number }
export interface DepthSegRow { source: string; bin_start: number; bin_end: number; cn_call: number; kept: number; mean: number | null; sd: number | null; allelicInside: number }
const round6 = (x: number) => Math.round(x * 1e6) / 1e6
export function depthSegmentRows(segs: LocalDepthSeg[], values: Record<string, (number | null)[]>, allelicBps: AllelicBreakpoint[], chrom: string, N: number): DepthSegRow[] {
  return segs.filter((g) => values[g.group]).map((g) => {
    const v = values[g.group].slice(g.s, g.e).filter((x): x is number => x !== null && Number.isFinite(x))
    const mean = v.length ? v.reduce((a, b) => a + b, 0) / v.length : null
    const sd = mean === null ? null : Math.sqrt(v.reduce((a, b) => a + (b - mean) ** 2, 0) / v.length)
    const inside = allelicBps.filter((b) => b.source === g.group && b.chromosome === chrom && b.absolute_bin - g.s > N && g.e - b.absolute_bin > N).length
    return { source: g.group, bin_start: g.s, bin_end: g.e - 1, cn_call: g.cn_call, kept: v.length, mean: mean === null ? null : round6(mean), sd: sd === null ? null : round6(sd), allelicInside: inside }
  })
}
export function callShares(seg: { bin_start: number; bin_end: number }, depthSegs: LocalDepthSeg[]): { call: number; frac: number }[] {
  const n = seg.bin_end - seg.bin_start + 1, acc = new Map<number, number>()
  for (const d of depthSegs) { const lo = Math.max(d.s, seg.bin_start), hi = Math.min(d.e - 1, seg.bin_end); if (hi >= lo) acc.set(d.cn_call, (acc.get(d.cn_call) ?? 0) + hi - lo + 1) }
  return [...acc].map(([call, k]) => ({ call, frac: k / n })).sort((a, b) => b.frac - a.frac || a.call - b.call)
}
export const sharesText = (sh: { call: number; frac: number }[]) => (sh.length ? sh.map((x) => `${x.call.toFixed(1)}: ${Math.round(x.frac * 100)}%`).join(' · ') : 'no depth call')
export interface AllelicSegRow { source: string; bin_start: number; bin_end: number; kept: number; mean: number | null; spanned: string }
export function allelicSegmentRows(segs: AllelicSegment[], depthSegs: LocalDepthSeg[], sources: string[], chrom: string): AllelicSegRow[] {
  return segs.filter((s) => s.chromosome === chrom && sources.includes(s.source)).map((s) => ({
    source: s.source, bin_start: s.bin_start, bin_end: s.bin_end, kept: s.kept_bins, mean: s.mean,
    spanned: sharesText(callShares(s, depthSegs.filter((d) => d.group === s.source))),
  }))
}

export const COMPARISON_DEFINITIONS = [
  'Matching is one-to-one, nearest first: all depth–allelic pairs within N bins are sorted by distance and accepted greedily, skipping any breakpoint already matched; ties go to the lower depth bin.',
  'An allelic breakpoint inside a gap of skipped bins is at distance 0 from any bin of that gap; outside it, the distance is to the nearer gap edge.',
  'Distances count genomic bins, not kept bins. Median offset is taken over matched pairs.',
  'Depth segments: kept bins have AT ≥ min AT; mean and sd (population) of the allelic signal over kept bins; "allelic bps inside" counts allelic breakpoints more than N bins from both segment edges.',
  'Allelic segments: each depth CN call with its share of the segment\'s bins, largest first.',
]
const num = (v: number | null, d = 4) => (v === null || !Number.isFinite(v) ? 'NA' : Number.isInteger(v) ? String(v) : v.toFixed(d))
export function runHeader(job: AllelicJob, N: number, depthSource: 'current' | 'original', nRounds: number, scope: string): string[] {
  return [`signal: ${job.signal} (${signalLabel(job.signal)})`, `settings: ${JSON.stringify(job.settings)}`, `match N: ${N} bins`,
    `depth segments: ${depthSource === 'current' ? `current round (${nRounds} merge rounds)` : 'original'}`, `table 1 scope: ${scope}`]
}
export function comparisonTsv(header: string[], t1: AgreementRow[], t2: DepthSegRow[], t3: AllelicSegRow[]): string {
  const lines = header.map((h) => `# ${h}`)
  lines.push('# table 1: breakpoint agreement', ['source', 'depth_bps', 'allelic_bps', 'matched', 'depth_only', 'allelic_only', 'median_offset'].join('\t'))
  for (const r of t1) lines.push([r.source, r.depth, r.allelic, r.matched, r.depthOnly, r.allelicOnly, num(r.medianOffset, 1)].join('\t'))
  lines.push('# table 2: depth segments seen by the allelic signal', ['source', 'bin_start', 'bin_end', 'cn_call', 'kept_bins', 'mean', 'sd', 'allelic_bps_inside'].join('\t'))
  for (const r of t2) lines.push([r.source, r.bin_start, r.bin_end, r.cn_call, r.kept, num(r.mean), num(r.sd), r.allelicInside].join('\t'))
  lines.push('# table 3: allelic segments seen by depth', ['source', 'bin_start', 'bin_end', 'kept_bins', 'mean', 'depth_calls_spanned'].join('\t'))
  for (const r of t3) lines.push([r.source, r.bin_start, r.bin_end, r.kept, num(r.mean), r.spanned].join('\t'))
  return lines.join('\n') + '\n'
}

// ---------------- session (spec §4)
export interface AllelicView { signal: SignalId; settings: AdSettings; chrom: string; sources: string[]; depthSource: 'current' | 'original'; matchN: number; scope: 'chrom' | 'genome' }
export const DEFAULT_VIEW: AllelicView = { signal: 'abs_dbaf', settings: DEFAULT_SETTINGS, chrom: 'chr8', sources: ['cluster0', 'cluster1', 'cluster2', 'cluster3', 'cluster4'], depthSource: 'current', matchN: 3, scope: 'chrom' }
export function allelicFromSession(p: any): { view: AllelicView; notice: string | null } {
  if (!p || typeof p !== 'object') return { view: DEFAULT_VIEW, notice: null }
  const bad: string[] = []
  const signal = SIGNALS.some((s) => s.id === p.signal) ? p.signal : (bad.push(`signal ${JSON.stringify(p.signal)}`), DEFAULT_VIEW.signal)
  const settings = { ...DEFAULT_SETTINGS, ...(p.settings && typeof p.settings === 'object' ? p.settings : {}) }
  const sp = settingsProblems(settings); if (sp.length) bad.push(...sp)
  const chrom = typeof p.chrom === 'string' ? p.chrom : (bad.push('chromosome'), DEFAULT_VIEW.chrom)
  const sources = Array.isArray(p.sources) && p.sources.every((s: unknown) => typeof s === 'string') ? p.sources : (bad.push('sources'), DEFAULT_VIEW.sources)
  const depthSource = p.depthSource === 'original' || p.depthSource === 'current' ? p.depthSource : (bad.push('depth segments'), DEFAULT_VIEW.depthSource)
  const matchN = Number.isInteger(p.matchN) && p.matchN >= 0 && p.matchN <= 50 ? p.matchN : (bad.push('match N'), DEFAULT_VIEW.matchN)
  const scope = p.scope === 'chrom' || p.scope === 'genome' ? p.scope : (bad.push('scope'), DEFAULT_VIEW.scope)
  if (bad.length) return { view: DEFAULT_VIEW, notice: `Saved allelic settings could not be applied: ${bad.join('; ')}` }
  return { view: { signal, settings: sp.length ? DEFAULT_SETTINGS : settings, chrom, sources, depthSource, matchN, scope }, notice: null }
}
```

- [ ] **Step 4: Run the helper tests to check they pass**

Run: `web/node_modules/.bin/esbuild tests/allelic_helpers.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/allelic_helpers.mjs && node tests/_bundle/allelic_helpers.mjs`
Expected: `all allelic helper checks passed`.

- [ ] **Step 5: Write the failing store test**

`tests/allelic_store.ts`:

```ts
// tests/allelic_store.ts
// Run with: web/node_modules/.bin/esbuild tests/allelic_store.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/allelic_store.mjs
//           node tests/_bundle/allelic_store.mjs
import { useAllelic, timing } from '../web/src/allelicStore'
import { DEFAULT_SETTINGS, DEFAULT_VIEW } from '../web/src/allelic'

let failures = 0
function eq(actual: unknown, expected: unknown, label: string) {
  const a = JSON.stringify(actual), e = JSON.stringify(expected)
  if (a !== e) { failures++; console.error(`FAIL ${label}: got ${a}, expected ${e}`) } else console.log(`ok ${label}`)
}
timing.pollMs = 1
const calls: string[] = []
const jobs: Record<string, any[]> = {}
const reply = (body: unknown, status = 200) => ({ ok: status < 400, status, json: async () => body }) as any
;(globalThis as any).fetch = async (url: string, init?: any) => {
  calls.push(`${init?.method ?? 'GET'} ${url}`)
  if (url.startsWith('/api/allelic/status')) return reply({ available: true, state: 'ready', problems: [], notes: [], allelic_dir: '/x', sources: ['cluster0', 'all'] })
  if (url.startsWith('/api/allelic/tracks')) return reply({ chrom: 'chr8', signal: 'abs_dbaf', min_at: 10, n: 2, tracks: [] })
  if (url === '/api/allelic/segment') { const sig = JSON.parse(init.body).signal; return reply(jobs[sig].shift()) }
  const id = url.split('/').pop()!
  const sig = Object.keys(jobs).find((k) => k.startsWith(id.slice(0, 1)))!
  return reply(jobs[sig].shift())
}
const job = (id: string, signal: string, status: string) => ({ id, status, signal, settings: DEFAULT_SETTINGS, progress: null, error: null, result: status === 'done' ? { breakpoints: [], audit: [], segments: [] } : null })

;(async () => {
  const s = useAllelic.getState()
  await s.loadStatus()
  eq(useAllelic.getState().status?.state, 'ready', 'status loaded')
  await useAllelic.getState().ensureTracks()
  await useAllelic.getState().ensureTracks()
  eq(calls.filter((c) => c.includes('/tracks')).length, 1, 'tracks fetched once per key')

  jobs.abs_dbaf = [job('a1', 'abs_dbaf', 'queued'), job('a1', 'abs_dbaf', 'running'), job('a1', 'abs_dbaf', 'done')]
  await useAllelic.getState().run()
  eq([useAllelic.getState().shown?.id, useAllelic.getState().active], ['a1', null], 'run polls until done')

  // run supersedes: a second run started while the first polls; the first must not overwrite shown
  jobs.dbaf = [job('d1', 'dbaf', 'running'), job('d1', 'dbaf', 'running'), job('d1', 'dbaf', 'done')]
  jobs.shrunk_abs_dbaf = [job('s1', 'shrunk_abs_dbaf', 'done')]
  useAllelic.getState().setView({ signal: 'dbaf' })
  const first = useAllelic.getState().run()
  useAllelic.getState().setView({ signal: 'shrunk_abs_dbaf' })
  await useAllelic.getState().run()
  await first
  eq(useAllelic.getState().shown?.id, 's1', 'run supersedes')

  jobs.abs_dbaf = [{ ...job('a2', 'abs_dbaf', 'error'), error: 'kaput' }]
  useAllelic.getState().setView({ signal: 'abs_dbaf', settings: { ...DEFAULT_SETTINGS, alpha: 0.01 } })
  await useAllelic.getState().run()
  eq([useAllelic.getState().runError, useAllelic.getState().shown?.id], ['kaput', 's1'], 'error keeps the previous result')

  useAllelic.getState().setView({ settings: { ...DEFAULT_SETTINGS, alpha: 0 } })
  await useAllelic.getState().run()
  eq(useAllelic.getState().runError, 'alpha: must be > 0 and ≤ 0.1', 'invalid settings never reach the server')

  const saved = useAllelic.getState().serialize()
  eq(Object.keys(saved).sort(), Object.keys(DEFAULT_VIEW).sort(), 'serialize writes the view only')
  useAllelic.getState().restore({ signal: 'nope' })
  eq([useAllelic.getState().signal, useAllelic.getState().shown, useAllelic.getState().notice?.startsWith('Saved allelic settings could not be applied')], ['abs_dbaf', null, true], 'restore: defaults, result cleared, notice')

  if (failures) { console.error(`${failures} failure(s)`); process.exit(1) } else console.log('all allelic store checks passed')
})()
```

- [ ] **Step 6: Write `web/src/allelicStore.ts`**

```ts
// Allelic segmentation tab state (spec §4). Saved with the session under `allelic`.
import { create } from 'zustand'
import { DEFAULT_VIEW, allelicApi, allelicFromSession, settingsProblems, type AllelicJob, type AllelicStatus, type AllelicTracks, type AllelicView } from './allelic'

export const timing = { pollMs: 1000 }
const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms))
const tracksKeyOf = (v: AllelicView) => `${v.chrom}|${v.sources.join(',')}|${v.signal}|${v.settings.min_at}`
const pickView = (s: AllelicView): AllelicView => ({ signal: s.signal, settings: s.settings, chrom: s.chrom, sources: s.sources, depthSource: s.depthSource, matchN: s.matchN, scope: s.scope })
let runToken = 0

interface AllelicState extends AllelicView {
  status: AllelicStatus | null; statusError: string | null
  tracks: AllelicTracks | null; tracksKey: string | null; tracksError: string | null
  shown: AllelicJob | null; active: AllelicJob | null; runError: string | null; notice: string | null
  loadStatus: () => Promise<void>
  setView: (p: Partial<AllelicView>) => void
  ensureTracks: () => Promise<void>
  run: () => Promise<void>
  serialize: () => AllelicView
  restore: (p: any) => void
}

export const useAllelic = create<AllelicState>((set, get) => ({
  ...DEFAULT_VIEW,
  status: null, statusError: null, tracks: null, tracksKey: null, tracksError: null, shown: null, active: null, runError: null, notice: null,
  loadStatus: async () => {
    try {
      const status = await allelicApi.status()
      set({ status, statusError: null })
      if (status.state === 'idle' || status.state === 'building') setTimeout(() => { get().loadStatus() }, 2000)
    } catch (e: any) { set({ statusError: e.message || String(e) }) }
  },
  setView: (p) => set(p),
  ensureTracks: async () => {
    const s = get()
    if (s.status?.state !== 'ready' || !s.sources.length) return
    const key = tracksKeyOf(s)
    if (key === s.tracksKey) return
    set({ tracksKey: key, tracksError: null })
    try {
      const tracks = await allelicApi.tracks(s.chrom, s.sources, s.signal, s.settings.min_at)
      if (get().tracksKey === key) set({ tracks })
    } catch (e: any) { if (get().tracksKey === key) set({ tracksError: e.message || String(e), tracksKey: null }) }
  },
  run: async () => {
    const s = get()
    const problems = settingsProblems(s.settings)
    if (problems.length) { set({ runError: problems.join('; ') }); return }
    const token = ++runToken
    set({ runError: null })
    try {
      let job = await allelicApi.start(s.signal, s.settings)
      if (token !== runToken) return
      set({ active: job })
      while (job.status === 'queued' || job.status === 'running') {
        await sleep(timing.pollMs)
        if (token !== runToken) return
        job = await allelicApi.job(job.id)
        if (token !== runToken) return
        set({ active: job })
      }
      if (job.status === 'done') set({ shown: job, active: null })
      else set({ active: null, runError: job.error || 'segmentation failed' })
    } catch (e: any) { if (token === runToken) set({ active: null, runError: e.message || String(e) }) }
  },
  serialize: () => pickView(get()),
  restore: (p) => {
    runToken++
    const { view, notice } = allelicFromSession(p)
    set({ ...view, notice, shown: null, active: null, runError: null, tracks: null, tracksKey: null, tracksError: null })
  },
}))
```

- [ ] **Step 7: Hook the session**

In `web/src/store.ts`:
- Add the import near the other store imports: `import { useAllelic } from './allelicStore'`.
- In `serialize`, add `allelic: useAllelic.getState().serialize(),` after the `integration: useIntegration.getState().serialize(),` line.
- In `restore`, add `useAllelic.getState().restore(p.allelic)` after `useIntegration.getState().restore(p.integration)`.

- [ ] **Step 8: Run all web tests and the build**

```bash
for t in allelic_helpers allelic_store; do web/node_modules/.bin/esbuild tests/$t.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/$t.mjs && node tests/_bundle/$t.mjs || exit 1; done
for t in merge_rules_helpers merge_rounds_helpers merge_rounds_store cn_recall_helpers hatch_merge_helpers track_helpers integration_helpers experiment_overlays; do web/node_modules/.bin/esbuild tests/$t.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/$t.mjs && node tests/_bundle/$t.mjs > /dev/null || { echo "FAIL $t"; exit 1; }; done
cd web && npm run build
```

Expected:
- `all allelic helper checks passed` and `all allelic store checks passed`;
- the 8 earlier suites run silently, with no `FAIL` line;
- the build succeeds.

- [ ] **Step 9: Commit**

```bash
git add web/src/allelic.ts web/src/allelicStore.ts web/src/store.ts tests/allelic_helpers.ts tests/allelic_store.ts
git commit -m "feat: allelic segmentation web logic, store and session state

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: The tab

**Files:**
- Create: `web/src/components/AllelicSegmentationView.tsx`
- Modify: `web/src/store.ts` (`Tab` type)
- Modify: `web/src/App.tsx` (nav button and panel)
- Modify: `web/src/components/integration.css` (allelic styles)

**Interfaces:**
- Consumes:
  - Task 4's `useAllelic` and `allelic.ts` exports.
  - `useIntegration` (`meta`, `pb`, `rounds`, `base`, `load`, `ensurePseudobulk`, `refreshBase`).
  - `useStore` (`meta`, `bins`, `presentation`).
  - From `integration.ts`: `SOURCE_COLORS`, `STORY_CN`, `DEPTH_COLORS`, `levelTrack`, `withGaps`, `clusterSegmentsLocal`.
  - `mergeRounds.wireRounds`, `cnRecall.downloadText`, `Plot`/`DeferredPlot`, `MB`.
- Produces: the tab `allelicseg`.

- [ ] **Step 1: Add the tab type and the nav entry**

In `web/src/store.ts`, change the `Tab` type to:

```ts
export type Tab = 'cohort' | 'explore' | 'manual' | 'experiment' | 'clusters' | 'clusterseg' | 'allelicseg'
```

In `web/src/App.tsx`:
- Add the import: `import AllelicSegmentationView from './components/AllelicSegmentationView'`.
- After the `Integration · segmentation &amp; CN` nav button, add:

```tsx
          <button className={tab === 'allelicseg' ? 'on' : ''} onClick={() => setTab('allelicseg')}>Integration · allelic segmentation</button>
```

- After the `clusterseg` panel `div`, add:

```tsx
          <div style={{ display: tab === 'allelicseg' ? 'block' : 'none' }}>{tab === 'allelicseg' && <ErrorBoundary label="Integration · allelic segmentation"><AllelicSegmentationView /></ErrorBoundary>}</div>
```

- [ ] **Step 2: Write the component**

`web/src/components/AllelicSegmentationView.tsx`:

```tsx
// Integration · allelic segmentation (spec docs/superpowers/specs/2026-10-28-allelic-segmentation-design.md §4–5)
import { useEffect, useMemo, useRef } from 'react'
import Plot from './DeferredPlot'
import { downloadPlot } from './Plot'
import { useAllelic } from '../allelicStore'
import { useIntegration } from '../integrationStore'
import { useStore } from '../store'
import { MB } from '../coords'
import { DEPTH_COLORS, SOURCE_COLORS, STORY_CN, clusterSegmentsLocal, levelTrack, withGaps } from '../integration'
import { wireRounds } from '../mergeRounds'
import { downloadText } from '../cnRecall'
import { COMPARISON_DEFINITIONS, POOL_NOTE, SETTING_FIELDS, SIGNALS, agreementRows, allelicSegmentRows, comparisonTsv, depthSegmentRows, runHeader, sameRun, signalLabel, type AdSettings, type LocalDepthSeg } from '../allelic'

const fmt = (v: number | null | undefined, d = 3) => (v === null || v === undefined || !Number.isFinite(v) ? 'N/A' : v.toFixed(d))

export default function AllelicSegmentationView() {
  const al = useAllelic()
  const ig = useIntegration()
  const exMeta = useStore((s) => s.meta)!
  const bins = useStore((s) => s.bins)!
  const presentation = useStore((s) => s.presentation)
  const plotEl = useRef<any>(null)
  const meta = ig.meta
  const roundsSig = JSON.stringify(wireRounds(ig.rounds))
  const nRounds = ig.rounds.length
  useEffect(() => { ig.load(); al.loadStatus() }, [])
  useEffect(() => { if (meta) ig.ensurePseudobulk() }, [meta])
  useEffect(() => { if (meta && nRounds) ig.refreshBase() }, [meta, roundsSig])
  useEffect(() => { al.ensureTracks() }, [al.status?.state, al.chrom, al.sources.join(','), al.signal, al.settings.min_at])

  const useOriginal = al.depthSource === 'original' || nRounds === 0
  const depthBps = useOriginal ? (meta?.breakpoints ?? null) : (ig.base?.breakpoints ?? null)
  const depthSegsAll = useOriginal ? (meta?.cluster_segments ?? null) : (ig.base?.cluster_segments ?? null)
  const chromMeta = exMeta.chromosomes.find((c) => c.name === al.chrom) ?? exMeta.chromosomes[0]
  const localSegs: LocalDepthSeg[] = useMemo(() => clusterSegmentsLocal(depthSegsAll, exMeta.chromosomes).filter((g) => g.chromosome === chromMeta.name).map((g) => ({ group: g.group, s: g.s, e: g.e, cn_call: g.cn_call })), [depthSegsAll, chromMeta])
  const shown = al.shown
  const stale = !!shown && !sameRun(shown, al.signal, al.settings)
  const signed = SIGNALS.find((s) => s.id === al.signal)!.signed
  const sources = (meta?.sources ?? []).filter((s) => al.sources.includes(s))
  const trackValues = useMemo(() => Object.fromEntries((al.tracks?.tracks ?? []).map((t) => [t.source, t.values])), [al.tracks])

  // ---------------- tracks
  const fig = useMemo(() => {
    if (!meta || !ig.pb || !al.tracks || al.tracks.chrom !== chromMeta.name || !sources.length) return null
    const nBins = exMeta.dataset.identity.n_bins
    const mb: number[] = []; for (let i = 0; i < chromMeta.n; i++) mb.push(bins[2 * (chromMeta.offset + i)] / MB)
    const endsMb = Array.from({ length: chromMeta.n }, (_, i) => bins[2 * (chromMeta.offset + i) + 1] / MB)
    const xEnd = endsMb[endsMb.length - 1]
    const xAt = (b: number) => (b < chromMeta.n ? mb[b] : xEnd)
    const data: any[] = [], shapes: any[] = [], layout: any = {}
    const nP = sources.length * 2, rowH = 1 / nP, gap = 0.05 / nP
    sources.forEach((src, k) => {
      const si = meta.sources.indexOf(src)
      const dAx = 2 * k + 1, aAx = 2 * k + 2
      const dRef = `y${dAx === 1 ? '' : dAx}`, aRef = `y${aAx}`
      const dDom: [number, number] = [1 - (2 * k + 1) * rowH + gap, 1 - 2 * k * rowH]
      const aDom: [number, number] = [1 - (2 * k + 2) * rowH + gap, 1 - (2 * k + 1) * rowH]
      // depth panel: pseudobulk, CN strip, depth breakpoints
      const lt = levelTrack(ig.pb!.subarray(si * nBins, (si + 1) * nBins), chromMeta, 'chrmedian', ig.pb!.subarray(si * nBins, (si + 1) * nBins))
      const g = withGaps(mb, lt.y, endsMb)
      data.push({ type: 'scatter', mode: 'lines', x: g.x, y: g.y, line: { color: '#6b7280', width: 0.8 }, connectgaps: false, xaxis: 'x', yaxis: dRef, hovertemplate: `${src} depth · %{x:.2f} Mb · %{y:.2f}<extra></extra>`, showlegend: false })
      layout[`yaxis${dAx === 1 ? '' : dAx}`] = { domain: dDom, range: [-4.6, 3.2], title: { text: `${src}<br><span style="font-size:8px">depth log2</span>`, font: { size: 10, color: SOURCE_COLORS[src] } }, gridcolor: '#f1f5f9', zeroline: false, tickfont: { size: 8 } }
      for (const cs of localSegs.filter((c) => c.group === src)) shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: xAt(cs.s), x1: xAt(cs.e), y0: dDom[1] - 0.04 * rowH, y1: dDom[1], fillcolor: STORY_CN[cs.cn_call], line: { width: 0 }, layer: 'below' })
      const dB = (depthBps ?? []).filter((b) => b.source === src && b.chromosome === chromMeta.name)
      for (const b of dB) {
        const x = xAt(b.absolute_bin)
        shapes.push({ type: 'line', xref: 'x', yref: `${dRef} domain`, x0: x, x1: x, y0: 0, y1: 1, line: { color: DEPTH_COLORS[Math.min(b.recursion_depth, DEPTH_COLORS.length - 1)], width: 1.2 }, name: 'depth breakpoint' })
        shapes.push({ type: 'line', xref: 'x', yref: `${aRef} domain`, x0: x, x1: x, y0: 0, y1: 1, line: { color: '#9ca3af', width: 1, dash: 'dash' }, opacity: 0.6, name: 'depth breakpoint (echo)' })
      }
      // allelic panel: dots by AT, skipped ticks, segment means, breakpoints and gap bands
      const t = al.tracks!.tracks.find((x) => x.source === src)
      if (t) {
        const ats = t.at.filter((a, i) => t.values[i] !== null).sort((a, b) => a - b)
        const p90 = ats.length ? ats[Math.floor(0.9 * (ats.length - 1))] : al.settings.min_at
        const opa = (a: number) => Math.max(0.15, Math.min(1, p90 > al.settings.min_at ? 0.15 + 0.85 * (a - al.settings.min_at) / (p90 - al.settings.min_at) : 1))
        const xs: number[] = [], ys: number[] = [], os: number[] = [], txt: string[] = [], skipped: number[] = []
        t.values.forEach((v, i) => { if (v === null) skipped.push(mb[i]); else { xs.push(mb[i]); ys.push(v); os.push(opa(t.at[i])); txt.push(`AT ${t.at[i]}`) } })
        data.push({ type: 'scatter', mode: 'markers', x: xs, y: ys, text: txt, marker: { size: 3, color: SOURCE_COLORS[src], opacity: os }, xaxis: 'x', yaxis: aRef, hovertemplate: `${src} ${signalLabel(al.signal)} %{y:.3f} · %{text}<extra></extra>`, showlegend: false, name: 'allelic bins' })
        const lo = signed ? -0.5 : 0
        if (skipped.length) data.push({ type: 'scatter', mode: 'markers', x: skipped, y: skipped.map(() => lo), marker: { symbol: 'line-ns-open', size: 4, color: '#9ca3af' }, xaxis: 'x', yaxis: aRef, hoverinfo: 'skip', showlegend: false, name: 'skipped bins' })
      }
      layout[`yaxis${aAx}`] = { domain: aDom, range: signed ? [-0.5, 0.5] : [0, 0.5], title: { text: `${signalLabel(al.signal)}`, font: { size: 9, color: '#b5651d' } }, gridcolor: '#f8fafc', zeroline: signed, zerolinecolor: '#9ca3af', tickfont: { size: 8 } }
      if (shown && !stale) {
        for (const s of shown.result!.segments.filter((x) => x.source === src && x.chromosome === chromMeta.name && x.mean !== null))
          shapes.push({ type: 'line', xref: 'x', yref: aRef, x0: xAt(s.bin_start), x1: xAt(s.bin_end + 1), y0: s.mean, y1: s.mean, line: { color: '#111827', width: 2 }, name: 'allelic segment mean' })
        for (const b of shown.result!.breakpoints.filter((x) => x.source === src && x.chromosome === chromMeta.name)) {
          if (b.gap_start !== null && b.gap_end !== null) shapes.push({ type: 'rect', xref: 'x', yref: `${aRef} domain`, x0: xAt(b.gap_start), x1: xAt(b.gap_end), y0: 0, y1: 1, fillcolor: '#b5651d', opacity: 0.18, line: { width: 0 }, name: 'allelic gap breakpoint' })
          else shapes.push({ type: 'line', xref: 'x', yref: `${aRef} domain`, x0: xAt(b.absolute_bin), x1: xAt(b.absolute_bin), y0: 0, y1: 1, line: { color: '#b5651d', width: 1.6 }, name: 'allelic breakpoint' })
        }
      }
    })
    Object.assign(layout, { height: Math.max(300, 120 * nP + 60), margin: { l: 70, r: 14, t: 20, b: 40 }, paper_bgcolor: 'white', plot_bgcolor: 'white', font: { size: 11 }, hovermode: presentation ? false : 'closest', shapes, showlegend: false, uirevision: chromMeta.name, xaxis: { title: { text: `${chromMeta.name} position (Mb)` }, range: [mb[0], xEnd], anchor: `y${nP}`, gridcolor: '#f1f5f9' } })
    return { data, layout }
  }, [meta, ig.pb, al.tracks, chromMeta, sources.join(','), localSegs, depthBps, shown, stale, signed, al.signal, al.settings.min_at, presentation, bins])

  // ---------------- tables
  const chroms = al.scope === 'genome' ? exMeta.chromosomes.map((c) => c.name) : [chromMeta.name]
  const t1 = useMemo(() => (shown && depthBps ? agreementRows(sources, chroms, depthBps, shown.result!.breakpoints, al.matchN) : []), [shown, depthBps, sources.join(','), chroms.join(','), al.matchN])
  const t2 = useMemo(() => (shown ? depthSegmentRows(localSegs.filter((g) => sources.includes(g.group)), trackValues, shown.result!.breakpoints, chromMeta.name, al.matchN) : []), [shown, localSegs, trackValues, sources.join(','), chromMeta, al.matchN])
  const t3 = useMemo(() => (shown ? allelicSegmentRows(shown.result!.segments, localSegs, sources, chromMeta.name) : []), [shown, localSegs, sources.join(','), chromMeta])

  const setSetting = (k: Exclude<keyof AdSettings, 'global_pool'>, v: string) => al.setView({ settings: { ...al.settings, [k]: v === '' ? NaN : +v } })
  const st = al.status
  const disabled = !st || st.state !== 'ready'
  const progress = al.active?.progress
  return (
    <div className="allelic-view">
      <div className="section-heading"><div><div className="eyebrow">INTEGRATION · ALLELIC</div><h2>Allelic segmentation</h2>
        <p>The depth run's constrained AD scan, run on a per-cluster allelic pseudobulk signal. Depth tracks and segments are shown above each allelic track for comparison.</p></div></div>
      {al.notice && <p className="notice small">{al.notice}</p>}
      {st && st.state === 'error' && <div className="notice small" role="alert">Allelic data unavailable: {st.problems.join('; ')}</div>}
      {st && st.state === 'building' && <p className="muted small">Building allelic pseudobulks from {st.allelic_dir}…</p>}
      {al.statusError && <div className="notice small" role="alert">{al.statusError}</div>}

      <div className="panel">
        <div className="controls-row">
          <label>signal<select aria-label="Allelic signal" value={al.signal} onChange={(e) => al.setView({ signal: e.target.value as any })}>{SIGNALS.map((s) => <option key={s.id} value={s.id}>{s.label}</option>)}</select></label>
          <span className="muted small">{SIGNALS.find((s) => s.id === al.signal)!.formula}</span>
          <label>min AT<input aria-label="min AT" type="number" value={Number.isFinite(al.settings.min_at) ? al.settings.min_at : ''} onChange={(e) => setSetting('min_at', e.target.value)} style={{ width: 70 }} /></label>
          <button className="btn-sm" disabled={disabled || !!al.active} onClick={() => al.run()}>Run segmentation</button>
          {al.active && <span className="small" aria-live="polite">{al.active.status} · {progress ? `${progress.source} · ${progress.chromosome} · ${progress.done}/${progress.total}` : 'starting'}</span>}
          {al.runError && <span className="notice small" role="alert">{al.runError}</span>}
        </div>
        <details className="allelic-ad-settings"><summary>AD settings</summary>
          <div className="controls-row">{SETTING_FIELDS.filter((f) => f.key !== 'min_at').map((f) => <label key={f.key}>{f.label}<input aria-label={f.label} type="number" step={f.integer ? 1 : 0.0001} value={Number.isFinite(al.settings[f.key]) ? al.settings[f.key] : ''} onChange={(e) => setSetting(f.key, e.target.value)} style={{ width: 80 }} /></label>)}
            <label>global test pool<select aria-label="Global test pool" value={al.settings.global_pool} onChange={(e) => al.setView({ settings: { ...al.settings, global_pool: e.target.value as any } })}><option value="genome">genome (original pyEpi)</option><option value="chromosome">chromosome (depth run)</option></select></label></div>
          <p className="muted small">{POOL_NOTE}</p>
        </details>
        {stale && <p className="notice small">Segments shown are from an earlier run (signal {signalLabel(shown!.signal)}, settings {JSON.stringify(shown!.settings)}). Run again to update.</p>}
        <div className="controls-row">
          <label>chromosome<select aria-label="Allelic chromosome" value={chromMeta.name} onChange={(e) => al.setView({ chrom: e.target.value })}>{exMeta.chromosomes.map((c) => <option key={c.name}>{c.name}</option>)}</select></label>
          {(meta?.sources ?? []).map((s) => <label key={s} className="chk"><input type="checkbox" checked={al.sources.includes(s)} onChange={(e) => al.setView({ sources: e.target.checked ? [...al.sources, s] : al.sources.filter((x) => x !== s) })} /><span style={{ color: SOURCE_COLORS[s] }}>{s}</span></label>)}
          <div className="segmented" role="group" aria-label="Depth segments">
            <button className={`btn-sm${al.depthSource === 'current' ? ' on' : ''}`} aria-pressed={al.depthSource === 'current'} onClick={() => al.setView({ depthSource: 'current' })}>current round{nRounds ? ` (${nRounds})` : ''}</button>
            <button className={`btn-sm${al.depthSource === 'original' ? ' on' : ''}`} aria-pressed={al.depthSource === 'original'} onClick={() => al.setView({ depthSource: 'original' })}>original</button>
          </div>
        </div>
        {!shown && !disabled && <p className="muted small">Run segmentation to see allelic segments.</p>}
        {al.tracksError && <div className="notice small" role="alert">{al.tracksError}</div>}
        {fig ? <><Plot data={fig.data} layout={fig.layout} config={{ displayModeBar: presentation ? false : 'hover' }} onReady={(el) => { plotEl.current = el }} /><button className="btn-xs" onClick={() => plotEl.current && downloadPlot(plotEl.current, `allelic_${chromMeta.name}`)}>PNG</button></>
          : <div className="muted small">{sources.length ? 'loading tracks…' : 'select at least one cluster row'}</div>}
      </div>

      <div className="panel allelic-compare" role="region" aria-label="Depth vs allelic segments">
        <h3>Depth vs allelic segments</h3>
        {!shown ? <p className="muted small">The comparison needs a finished run.</p> : <>
          <div className="controls-row">
            <div className="segmented" role="group" aria-label="Comparison scope">
              <button className={`btn-sm${al.scope === 'chrom' ? ' on' : ''}`} aria-pressed={al.scope === 'chrom'} onClick={() => al.setView({ scope: 'chrom' })}>this chromosome</button>
              <button className={`btn-sm${al.scope === 'genome' ? ' on' : ''}`} aria-pressed={al.scope === 'genome'} onClick={() => al.setView({ scope: 'genome' })}>genome-wide</button>
            </div>
            <label>match N (bins)<input aria-label="Match N" type="number" min={0} max={50} value={al.matchN} onChange={(e) => al.setView({ matchN: Math.max(0, Math.min(50, Math.round(+e.target.value || 0))) })} style={{ width: 60 }} /></label>
            <button className="btn-sm" onClick={() => downloadText(comparisonTsv(runHeader(shown, al.matchN, al.depthSource, nRounds, al.scope === 'genome' ? 'genome-wide' : chromMeta.name), t1, t2, t3), 'allelic_vs_depth.tsv')}>allelic_vs_depth.tsv</button>
          </div>
          <h4>Breakpoint agreement · {al.scope === 'genome' ? 'genome-wide' : chromMeta.name}</h4>
          <div className="table-scroll"><table className="tbl small" aria-label="Breakpoint agreement"><thead><tr><th>source</th><th>depth bps</th><th>allelic bps</th><th>matched</th><th>depth only</th><th>allelic only</th><th>median offset</th></tr></thead>
            <tbody>{t1.map((r) => <tr key={r.source}><td style={{ color: SOURCE_COLORS[r.source] }}>{r.source}</td><td>{r.depth}</td><td>{r.allelic}</td><td>{r.matched}</td><td>{r.depthOnly}</td><td>{r.allelicOnly}</td><td>{fmt(r.medianOffset, 1)}</td></tr>)}</tbody></table></div>
          <h4>Depth segments seen by the allelic signal · {chromMeta.name}</h4>
          <div className="table-scroll"><table className="tbl small" aria-label="Depth segments seen by the allelic signal"><thead><tr><th>source</th><th>bins</th><th>CN call</th><th>kept bins</th><th>mean</th><th>sd</th><th>allelic bps inside</th></tr></thead>
            <tbody>{t2.map((r) => <tr key={`${r.source}${r.bin_start}`}><td style={{ color: SOURCE_COLORS[r.source] }}>{r.source}</td><td>{r.bin_start}–{r.bin_end}</td><td>{r.cn_call.toFixed(1)}</td><td>{r.kept}</td><td>{fmt(r.mean)}</td><td>{fmt(r.sd)}</td><td>{r.allelicInside || ''}</td></tr>)}</tbody></table></div>
          <h4>Allelic segments seen by depth · {chromMeta.name}</h4>
          <div className="table-scroll"><table className="tbl small" aria-label="Allelic segments seen by depth"><thead><tr><th>source</th><th>bins</th><th>kept bins</th><th>mean</th><th>depth calls spanned</th></tr></thead>
            <tbody>{t3.map((r) => <tr key={`${r.source}${r.bin_start}`}><td style={{ color: SOURCE_COLORS[r.source] }}>{r.source}</td><td>{r.bin_start}–{r.bin_end}</td><td>{r.kept}</td><td>{fmt(r.mean)}</td><td>{r.spanned}</td></tr>)}</tbody></table></div>
          <ul className="small muted" aria-label="Comparison definitions">{COMPARISON_DEFINITIONS.map((d) => <li key={d}>{d}</li>)}</ul>
        </>}
      </div>
    </div>
  )
}
```

- [ ] **Step 3: Add the styles**

Append to `web/src/components/integration.css`:

```css
/* allelic segmentation tab */
.allelic-view .allelic-ad-settings { margin: 6px 0; }
.allelic-view .allelic-ad-settings summary { cursor: pointer; font-size: 12px; color: #374151; }
.allelic-view .chk { display: inline-flex; align-items: center; gap: 4px; }
.allelic-compare h4 { margin: 12px 0 4px; font-size: 12px; }
```

- [ ] **Step 4: Build**

Run: `cd web && npm run build`
Expected: success, with no TypeScript errors.

If `onReady`, `config` or `downloadPlot` signatures differ from `IntegrationSegmentationView.tsx`, follow that file's usage exactly. Check `ClusterSegment` field names (`group`, `chromosome`, `cn_call`) against `web/src/types.ts`.

- [ ] **Step 5: Check it on a private server**

```bash
python3 run.py --port 8799 > /tmp/allelic_8799.log 2>&1 & echo $! > /tmp/allelic_8799.pid
```

Open `http://127.0.0.1:8799` in a headless Playwright session (or the Claude-in-Chrome tools, in a new tab only) and check:
1. The tab "Integration · allelic segmentation" opens. On chr8 with all five clusters, there are 10 panels: depth, then allelic, for each cluster.
2. Before a run, the hint "Run segmentation to see allelic segments." shows and the dots are visible.
3. Run segmentation shows the progress line, then the segment means, breakpoints and three tables.
4. Changing the signal to "signed dBAF" shows the earlier-run note, and the y range becomes −0.5 to 0.5.
5. "original" changes the depth breakpoints only when merge rounds exist. With no rounds the two are the same.

Then stop only the 8799 process: `kill $(cat /tmp/allelic_8799.pid)`.

- [ ] **Step 6: Commit**

```bash
git add web/src/components/AllelicSegmentationView.tsx web/src/App.tsx web/src/store.ts web/src/components/integration.css
git commit -m "feat: Integration · allelic segmentation tab with tracks and comparison tables

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Walkthrough, real-data record and docs

**Files:**
- Create: `tests/browser_allelic_walkthrough.cjs`
- Create (by running): `validation/allelic/allelic_check_results.json`, `validation/allelic/allelic.png`, `validation/allelic/README.md`
- Modify: `README.md`, `docs/USAGE.md`, `docs/ARCHITECTURE.md`, `tests/README.md`

**Interfaces:**
- Consumes: everything above, through the browser and HTTP.
- Produces: the validation record and the docs.

- [ ] **Step 1: Write the walkthrough**

`tests/browser_allelic_walkthrough.cjs`:

```js
// Allelic segmentation walkthrough (2026-10-28): tab, tracks before a run, |dBAF| run, tables, signal switch, depth original.
// Requires a private server (python3 run.py --port 8799; BASE_URL=http://127.0.0.1:8799) and Playwright + Chromium:
//   PLAYWRIGHT_MODULE=/path/to/node_modules/playwright CHROMIUM_PATH=/path/to/chrome BASE_URL=http://127.0.0.1:8799 \
//   node tests/browser_allelic_walkthrough.cjs
// Writes validation/allelic/allelic_check_results.json and validation/allelic/allelic.png.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const BASE_URL=process.env.BASE_URL||'http://127.0.0.1:8766';
const fs=require('fs'),path=require('path'),assert=require('assert/strict');
(async()=>{
 const browser=await chromium.launch({headless:process.env.HEADED!=='1',executablePath:process.env.CHROMIUM_PATH});
 try {
  const page=await browser.newPage({viewport:{width:1400,height:1300},acceptDownloads:true}),errors=[];
  page.on('pageerror',e=>errors.push(String(e)));
  page.on('console',m=>{if(m.type()==='error')errors.push('console: '+m.text())});
  await page.goto(BASE_URL);await page.getByText('300 matching cells',{exact:true}).waitFor();
  await page.getByRole('button',{name:'Integration · allelic segmentation',exact:true}).click();
  await page.getByText('Run segmentation to see allelic segments.').waitFor({timeout:180000});
  await page.waitForFunction(()=>document.querySelectorAll('.allelic-view .js-plotly-plot').length>0,null,{timeout:120000});
  const panels=await page.evaluate(()=>Object.keys(document.querySelector('.allelic-view .js-plotly-plot').layout).filter(k=>/^yaxis\d*$/.test(k)).length);
  const t0=Date.now();
  await page.getByRole('button',{name:'Run segmentation'}).click();
  await page.getByRole('region',{name:'Depth vs allelic segments'}).getByRole('table',{name:'Breakpoint agreement'}).waitFor({timeout:900000});
  const runSeconds=Math.round((Date.now()-t0)/1000);
  const readTable=name=>page.getByRole('table',{name}).evaluate(t=>[...t.querySelectorAll('tbody tr')].map(r=>[...r.cells].map(c=>c.textContent.trim())));
  const chr8=await readTable('Breakpoint agreement');
  await page.getByRole('group',{name:'Comparison scope'}).getByRole('button',{name:'genome-wide'}).click();
  const genome=await readTable('Breakpoint agreement');
  const shapes=()=>page.evaluate(()=>{const s=document.querySelector('.allelic-view .js-plotly-plot').layout.shapes||[];const c={};for(const x of s)c[x.name||'other']=(c[x.name||'other']||0)+1;return c});
  const shapesAbs=await shapes();
  const [dl]=await Promise.all([page.waitForEvent('download'),page.getByRole('button',{name:'allelic_vs_depth.tsv'}).click()]);
  const tsv=fs.readFileSync(await dl.path(),'utf8');
  await page.getByLabel('Allelic signal').selectOption('dbaf');
  await page.getByText(/Segments shown are from an earlier run/).waitFor();
  const yRange=await page.evaluate(()=>document.querySelector('.allelic-view .js-plotly-plot').layout.yaxis2.range);
  await page.getByRole('group',{name:'Depth segments'}).getByRole('button',{name:'original'}).click();
  const out=path.join(__dirname,'..','validation','allelic');fs.mkdirSync(out,{recursive:true});
  await page.screenshot({path:path.join(out,'allelic.png'),fullPage:true});
  const result={panels,runSeconds,chr8,genome,shapesAbs,tsvHead:tsv.split('\n').slice(0,8),yRangeSigned:yRange,errors};
  fs.writeFileSync(path.join(out,'allelic_check_results.json'),JSON.stringify(result,null,1));
  assert.equal(panels,10,'5 clusters × (depth + allelic)');
  assert.ok(genome.length===5,'one agreement row per cluster');
  assert.ok(tsv.startsWith('# signal: abs_dbaf'),'export header');
  assert.deepEqual(yRange,[-0.5,0.5],'signed y range');
  assert.deepEqual(errors,[],'no page errors');
  console.log('PASS',JSON.stringify({runSeconds,genome}));
 } finally { await browser.close() }
})().catch(e=>{console.error(e);process.exit(1)});
```

- [ ] **Step 2: Run it on a private server**

```bash
python3 run.py --port 8799 > /tmp/allelic_8799.log 2>&1 & echo $! > /tmp/allelic_8799.pid
BASE_URL=http://127.0.0.1:8799 node tests/browser_allelic_walkthrough.cjs
```

Expected: a `PASS` line. Use the same `PLAYWRIGHT_MODULE`/`CHROMIUM_PATH` setup the existing walkthroughs use (see `tests/README.md`).

- [ ] **Step 3: Record each signal on real data**

With the 8799 server still running:

```bash
python3 - <<'EOF'
import json, time, urllib.request
B = 'http://127.0.0.1:8799'
def post(u, b): return json.load(urllib.request.urlopen(urllib.request.Request(B + u, json.dumps(b).encode(), {'content-type': 'application/json'})))
def get(u): return json.load(urllib.request.urlopen(B + u))
rec = {}
for sig in ('abs_dbaf', 'dbaf', 'shrunk_abs_dbaf', 'oriented_all'):
    t0 = time.time(); v = post('/api/allelic/segment', {'signal': sig, 'settings': {}})
    while v['status'] in ('queued', 'running'):
        time.sleep(2); v = get(f"/api/allelic/segment/{v['id']}")
    bps = v['result']['breakpoints'] if v['result'] else []
    rec[sig] = {'status': v['status'], 'seconds': round(time.time() - t0), 'breakpoints_per_source': {s: sum(b['source'] == s for b in bps) for s in ('cluster0','cluster1','cluster2','cluster3','cluster4','all')}, 'gap_breakpoints': sum(b['gap_start'] is not None for b in bps)}
    print(sig, rec[sig])
json.dump(rec, open('validation/allelic/signals.json', 'w'), indent=1)
EOF
kill $(cat /tmp/allelic_8799.pid)
```

Expected: four `done` lines. Stop only the 8799 process.

- [ ] **Step 4: Write `validation/allelic/README.md`**

The record, in this order:
1. **Depth reproduction gate:** the result from `depth_reproduction.json`.
2. **Walkthrough result.**
3. **Per-signal table:** breakpoints per source and gap breakpoints, from `signals.json`.
4. **Genome-wide table 1 for |dBAF|,** from `allelic_check_results.json`.
5. **The settings used.**

Write only what was measured. If a signal finds almost no breakpoints, say so plainly; that is a finding, not a failure.

- [ ] **Step 5: Update the docs**

Add these sections, following each file's existing style. Quote the comparison definitions verbatim from `COMPARISON_DEFINITIONS`.

- **`README.md`:** one line in the feature list for the allelic segmentation tab.
- **`docs/USAGE.md`:** an "Allelic segmentation" section covering:
  - the data path and `PYEPI_SCIENTIFIC_ALLELIC`;
  - the four signals with their formulas;
  - min AT and the kept bins;
  - the AD settings and their defaults (the depth run's);
  - the global test pool (genome, original pyEpi, by default; chromosome, as the depth run used) and why it matters for the comparison;
  - Run and the earlier-run note;
  - the depth-segments switch;
  - the three tables with the comparison definitions verbatim, and `allelic_vs_depth.tsv`.
- **`docs/ARCHITECTURE.md`:**
  - the three server modules and their responsibilities;
  - the routes;
  - the job cache in `cache/allelic/`;
  - `allelicStore` and the session key `allelic`.
- **`tests/README.md`:**
  - `test_allelic_scan.py`, `test_allelic_data.py` and `test_allelic_api.py` in the pytest list;
  - `allelic_helpers.ts` and `allelic_store.ts` with their esbuild commands;
  - `validate_allelic_depth_reproduction.py`;
  - `browser_allelic_walkthrough.cjs` with its private-server note.

- [ ] **Step 6: Run everything once more**

```bash
python3 -m pytest -q
for t in allelic_helpers allelic_store merge_rules_helpers merge_rounds_helpers merge_rounds_store cn_recall_helpers hatch_merge_helpers track_helpers integration_helpers experiment_overlays; do web/node_modules/.bin/esbuild tests/$t.ts --bundle --format=esm --platform=node --outfile=tests/_bundle/$t.mjs && node tests/_bundle/$t.mjs > /dev/null || { echo "FAIL $t"; exit 1; }; done
cd web && npm run build
```

Expected: everything passes.

- [ ] **Step 7: Commit**

```bash
git add tests/browser_allelic_walkthrough.cjs validation/allelic/allelic_check_results.json validation/allelic/allelic.png validation/allelic/signals.json validation/allelic/README.md README.md docs/USAGE.md docs/ARCHITECTURE.md tests/README.md
git commit -m "docs, test: allelic segmentation walkthrough, real-data record and docs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
