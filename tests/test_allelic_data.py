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
