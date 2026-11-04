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
