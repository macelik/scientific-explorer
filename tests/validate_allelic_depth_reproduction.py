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
