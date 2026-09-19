from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from server import experiments as exp
from server.config import Config
from server.data import DataStore


def test_edit_does_not_mutate_source_and_uses_half_open_bounds():
    source = np.array([2., 0., 4., 6., 8.])
    edited, info = exp.apply_edit(source, 1, 4, 'multiply', 0.5, False)
    np.testing.assert_array_equal(source, [2, 0, 4, 6, 8])
    np.testing.assert_array_equal(edited, [2, 0, 2, 3, 8])
    assert info['clipped_bins'] == 0


def test_addition_reports_clipping_and_total_preservation_is_explicit():
    edited, info = exp.apply_edit(np.array([2., 0., 4., 6.]), 0, 2, 'add', -3, True)
    np.testing.assert_allclose(edited, [0, 0, 4.8, 7.2])
    assert info['clipped_bins'] == 2
    assert info['normalization_factor'] == pytest.approx(1.2)
    assert edited.sum() == pytest.approx(12)


def test_replace_changes_zero_mass():
    edited, _ = exp.apply_edit(np.array([0., 2., 4.]), 0, 2, 'set', 3, False)
    np.testing.assert_array_equal(edited, [3, 3, 4])


@pytest.mark.parametrize('start,end,op,value,normalize', [
    (-1, 2, 'multiply', 1, False), (1, 1, 'multiply', 1, False),
    (0, 4, 'multiply', 1, False), (0, 2, 'multiply', -1, False),
    (0, 2, 'set', -1, False), (0, 2, 'add', float('nan'), False),
    (0, 3, 'set', 0, True), (0, 2, 'unknown', 1, False),
])
def test_invalid_edits_are_rejected(start, end, op, value, normalize):
    with pytest.raises(ValueError):
        exp.apply_edit(np.array([1., 2., 3.]), start, end, op, value, normalize)


def test_invalid_source_is_rejected():
    with pytest.raises(ValueError):
        exp.apply_edit(np.array([1., np.inf, 3.]), 0, 1, 'multiply', 2, False)


@pytest.fixture(scope='module')
def data():
    return DataStore(Config.from_env())


@pytest.fixture(scope='module')
def result(data):
    cell = 'TTAGGCTAGGCCGGAA-1'
    c = data.chrom_by_name['chr9']
    return exp.run_experiment(
        data.x_row(cell), data.production_clusters(cell),
        [{'name': v.name, 'offset': v.offset, 'n': v.n} for v in data.chroms],
        data.var_start, data.var_end,
        dict(cell=cell, chrom='chr9', start=169, end=210, operation='multiply', value=0.5, preserve_total=False),
        Config.from_env().package_dir,
    )


def test_original_pipeline_parity_and_fixed_boundary_control(data, result):
    from pyEpiAneufinder.get_breakpoints import recursive_getbp
    cell = 'TTAGGCTAGGCCGGAA-1'
    c = data.chrom_by_name['chr9']; full = data.x_row(cell).copy()
    edited, _ = exp.apply_edit(full, c.offset + 169, c.offset + 210, 'multiply', 0.5, False)
    expected = recursive_getbp(edited[c.offset:c.offset+c.n], edited, k=2, n_permutations=1000, alpha=0.001)
    assert result['edited']['boundaries'] == [v[0] for v in expected]
    prod = data.cn[data.cell_row[cell], c.offset:c.offset+c.n]
    np.testing.assert_array_equal(result['cn']['original']['cn5'], prod)
    fixed = data.sci.assign_gainloss_new(edited, data.production_clusters(cell).tolist())
    np.testing.assert_array_equal(result['cn']['fixed']['cn5'], fixed[4][c.offset:c.offset+c.n])
    assert result['cn']['fixed']['best_s'] == pytest.approx(fixed[5])
    np.testing.assert_array_equal(data.x_row(cell), full)


def test_original_boundary_is_retested_against_edited_full_genome(data, result):
    c = data.chrom_by_name['chr9']; full = data.x_row('TTAGGCTAGGCCGGAA-1')
    edited, _ = exp.apply_edit(full, c.offset+169, c.offset+210, 'multiply', 0.5, False)
    b = result['original']['nodes'][0]['argmax_b']; x = edited[c.offset:c.offset+c.n]
    obs = data.sci.dist_ad(x[:b], x[b:])
    pg = data.sci.global_permutation_test_ad(edited, b, len(x)-b, obs, 1000, random_state=42)
    assert result['edited_at_original_boundary']['p_global'] == pg
    assert result['edited_at_original_boundary']['observed_ad'] == pytest.approx(obs)


def test_splicing_keeps_other_chromosome_boundaries():
    prod = np.array([1, 1, 2, 2, 3, 3, 4, 4])
    edited = exp.splice_labels(prod, 2, 4, [1])
    np.testing.assert_array_equal(np.diff(edited[:2]) != 0, np.diff(prod[:2]) != 0)
    np.testing.assert_array_equal(np.diff(edited[6:]) != 0, np.diff(prod[6:]) != 0)
    assert len(set(edited[:2]) & set(edited[2:6])) == 0
    assert len(set(edited[6:]) & set(edited[2:6])) == 0
    np.testing.assert_array_equal(np.diff(edited[2:6]) != 0, [True, False, False])


def test_zero_genome_has_explicit_undefined_cn(data):
    result = exp.fit_cn(np.zeros(4), np.array([1, 1, 2, 2]), 0, 4, [])
    assert result['status'] == 'undefined'
    assert result['cn5'] is None


def test_duplicate_left_preserves_order_and_clips_at_chromosome_edge():
    source = np.array([9., 9., 9., 1., 2., 3., 8., 8.])
    request = dict(start=3, end=6, operation='duplicate', direction='left', length=5, seed=42, preserve_total=False)
    edited, info = exp.apply_region_edit(source, 0, 8, request)
    np.testing.assert_array_equal(edited, [1, 2, 3, 1, 2, 3, 8, 8])
    assert info['actual_extension'] == 3
    assert info['requested_extension'] == 5
    assert info['event_start'] == 0 and info['event_end'] == 6
    np.testing.assert_array_equal(source, [9, 9, 9, 1, 2, 3, 8, 8])


def test_extend_resamples_only_source_and_is_seeded():
    source = np.array([8., 1., 4., 0., 7., 7., 7., 7.])
    req = dict(start=1, end=4, operation='extend', direction='right', length=4, seed=71, preserve_total=False)
    one, info = exp.apply_region_edit(source, 0, 8, req)
    two, _ = exp.apply_region_edit(source, 0, 8, req)
    np.testing.assert_array_equal(one, two)
    np.testing.assert_array_equal(one[:4], source[:4])
    assert set(one[4:]) <= {0, 1, 4}
    assert info['event_start'] == 1 and info['event_end'] == 8


def test_simulate_uses_only_requested_production_state_pool():
    source = np.arange(1., 25.)
    calls = np.repeat([0., 1., 2.], 8)
    req = dict(start=3, end=7, operation='simulate', cn_state='loss', seed=42, preserve_total=False)
    edited, info = exp.apply_region_edit(source, 0, 24, req, calls)
    assert set(edited[3:7]) <= set(source[:8])
    np.testing.assert_array_equal(edited[:3], source[:3])
    assert info['donor_count'] == 8
    assert info['sampling'] == 'Independent empirical resampling; spatial correlation is not preserved'


def test_simulate_rejects_missing_donor_class():
    with pytest.raises(ValueError, match='eight'):
        exp.apply_region_edit(np.ones(10), 0, 10,
            dict(start=1,end=5,operation='simulate',cn_state='loss',seed=42,preserve_total=False), np.ones(10))


def test_completed_event_index_can_return_status_without_deadlock(monkeypatch):
    import threading
    from server import app as service
    monkeypatch.setattr(service, 'events_state', {'status': 'done', 'rows': [], 'done': 300, 'total': 300})
    completed = threading.Event()
    worker = threading.Thread(target=lambda: (service.start_event_index(), completed.set()), daemon=True)
    worker.start()
    assert completed.wait(1), 'Repeated index request deadlocks while reading status under the same lock'
