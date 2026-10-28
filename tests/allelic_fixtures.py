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
