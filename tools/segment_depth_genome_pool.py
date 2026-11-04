"""Depth cluster segmentation with original pyEpi's genome-wide global permutation test.

run_bootstrap_diagnostics.py (which made raw_accepted_breakpoints.tsv) sampled the global
test from each chromosome's own bins. Original pyEpi passes the whole cell (counts.X[i, :]),
so here each cluster pseudobulk's full genome is the pool. Output goes to a new folder; the
original files are never touched, and writing into <integration_dir>/data itself is refused.

Run from the repo root:  python3 tools/segment_depth_genome_pool.py [--out-dir DIR]
Writes raw_accepted_breakpoints.tsv, segmentation_candidate_audit.tsv and run.json.
"""
import argparse
import json
import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tests'))
from server import allelic_scan as sc  # noqa: E402
from server.config import Config  # noqa: E402
from validate_allelic_depth_reproduction import load_depth_pseudobulks  # noqa: E402

BP_COLS = ['source', 'chromosome', 'absolute_bin', 'bp_genomic_start', 'bp_genomic_end',
           'recursion_depth', 'ad_score', 'p_local', 'p_global', 'n_cells']
AUDIT_COLS = ['source', 'chromosome', 'recursion_depth', 'node_start_bin', 'node_end_bin', 'candidate_rank',
              'candidate_local_cut', 'candidate_absolute_bin', 'bp_genomic_start', 'bp_genomic_end', 'ad_score',
              'left_child_bins', 'right_child_bins', 'geometry_valid', 'rejection_reasons',
              'selected_for_significance', 'p_local', 'p_global', 'accepted']


def main() -> int:
    cfg = Config.from_env()
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out-dir', default=os.path.join(cfg.integration_dir, 'data', 'genome_pool'))
    args = ap.parse_args()
    out_dir = os.path.abspath(args.out_dir)
    if out_dir == os.path.abspath(os.path.join(cfg.integration_dir, 'data')):
        print(f'refusing to write into {out_dir}: it holds the original depth files', file=sys.stderr)
        return 2
    t0 = time.time()
    depth, pseudo, chrom_of, chroms, sizes = load_depth_pseudobulks(cfg)
    starts = depth.var['start'].to_numpy()
    ends = depth.var['end'].to_numpy()
    tasks, local_idx = [], {}
    for g, x_pb in pseudo.items():
        for ch in chroms:
            idx = np.flatnonzero(chrom_of == ch)
            local_idx[ch] = idx
            tasks.append((g, ch, x_pb[idx], np.ones(len(idx), bool), None, None, sc.DEFAULTS, cfg.package_dir, x_pb))
    workers = max(1, (os.cpu_count() or 4) - 2)
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as ex:
        results = list(ex.map(sc.scan_task, tasks))
    bp_rows, audit_rows = [], []
    for g, ch, out in results:  # task order is source order, then genome order
        idx = local_idx[ch]
        for b in sorted(out['breakpoints'], key=lambda r: r['bin']):
            gi = idx[b['bin']]
            bp_rows.append({'source': g, 'chromosome': ch, 'absolute_bin': b['bin'], 'bp_genomic_start': int(starts[gi]),
                            'bp_genomic_end': int(ends[gi]), 'recursion_depth': b['recursion_depth'], 'ad_score': b['ad_score'],
                            'p_local': b['p_local'], 'p_global': b['p_global'], 'n_cells': sizes[g]})
        for r in out['audit']:
            r = dict(r, source=g, chromosome=ch)
            gi = idx[r['bin']] if r.get('bin') is not None else None
            r['bp_genomic_start'] = int(starts[gi]) if gi is not None else None
            r['bp_genomic_end'] = int(ends[gi]) if gi is not None else None
            audit_rows.append(r)
    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame(bp_rows, columns=BP_COLS).to_csv(os.path.join(out_dir, 'raw_accepted_breakpoints.tsv'), sep='\t', index=False)
    audit = pd.DataFrame(audit_rows, columns=AUDIT_COLS)
    for c in ['candidate_rank', 'candidate_local_cut', 'candidate_absolute_bin', 'bp_genomic_start', 'bp_genomic_end',
              'left_child_bins', 'right_child_bins']:
        audit[c] = audit[c].astype('Int64')  # keep integers integer despite empty cells
    audit.to_csv(os.path.join(out_dir, 'segmentation_candidate_audit.tsv'), sep='\t', index=False)
    per_source = pd.Series([r['source'] for r in bp_rows]).value_counts().reindex(list(pseudo)).fillna(0).astype(int).to_dict()
    rec = {'settings': sc.DEFAULTS, 'pool': 'genome (original pyEpi)', 'seconds': round(time.time() - t0, 1), 'workers': workers,
           'n_breakpoints': len(bp_rows), 'breakpoints_per_source': per_source, 'n_audit_rows': len(audit_rows),
           'n_bins': int(len(chrom_of)), 'h5ad': cfg.h5ad, 'package_dir': cfg.package_dir}
    with open(os.path.join(out_dir, 'run.json'), 'w') as f:
        json.dump(rec, f, indent=1)
    print(json.dumps(rec))
    return 0


if __name__ == '__main__':
    sys.exit(main())
