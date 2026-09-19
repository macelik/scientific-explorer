from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from server import segment_estimators as est
from server import pyepi_adapter as sci
from server.config import Config


@pytest.fixture(autouse=True)
def scientific_functions():
    sci.load(Config.from_env().package_dir)


def test_upper_only_iqr_retains_zero_and_two_sided_can_remove_it():
    x = np.array([0., 2., 2., 2., 2., 20.])
    a = est.segment_summary(x, 'iqr_upper')
    b = est.segment_summary(x, 'iqr_two_sided')
    assert a['value'] == pytest.approx(1.6)
    assert a['retained'] == 5
    assert b['value'] == 2
    assert b['retained'] == 4
    np.testing.assert_array_equal(x, [0,2,2,2,2,20])


def test_default_matches_original_caller_exactly():
    x = np.array([0.,2.,2.,2.,2.,20.,4.,4.,4.,4.,4.,40.])
    labels = np.repeat([1,2],6)
    reference = sci.assign_gainloss_new(x, labels.tolist())
    result = est.assign_with_estimator(x, labels, 'arithmetic')
    for name,index in [('cn_int',0),('cn_cont',1),('cn5',4)]:
        np.testing.assert_array_equal(result[name], reference[index])
    assert result['best_s'] == reference[5]


@pytest.mark.parametrize('method,lower', [('iqr_upper',False),('iqr_two_sided',True)])
def test_iqr_uses_original_normalization_and_full_lengths(method, lower):
    x = np.array([0.,2.,2.,2.,2.,20.,4.,4.,4.,4.,4.,40.])
    labels = np.repeat([1,2],6)
    original = x.copy()
    baseline = sci.trimmed_mean_iqr(x, lb=False)
    means = np.array([sci.trimmed_mean_iqr(x[labels==i]/baseline,lb=lower) for i in [1,2]])
    median = np.median(np.repeat(means,6))
    scale,continuous,integer = sci.weighted_scale_search(means.copy(), np.array([6,6]), median)
    result = est.assign_with_estimator(x,labels,method)
    assert result['best_s'] == scale
    np.testing.assert_array_equal(result['cn_int'], np.repeat(integer,6))
    np.testing.assert_array_equal(result['cn_cont'], np.repeat(continuous,6))
    assert result['baseline'] == baseline
    np.testing.assert_array_equal(x,original)


def test_iqr_zero_scale_is_explicitly_rejected():
    with pytest.raises(ValueError,match='positive'):
        est.assign_with_estimator(np.zeros(8),np.repeat([1,2],4),'iqr_upper')


def test_invalid_method_does_not_silently_use_arithmetic():
    with pytest.raises(ValueError):
        est.segment_summary(np.ones(8),'unknown')


@pytest.mark.parametrize('method', ['arithmetic', 'iqr_upper', 'iqr_two_sided'])
def test_real_cell_api_retains_boundaries_and_reports_estimator_effect(monkeypatch, method):
    from server import app as service
    from server.data import DataStore
    data = DataStore(Config.from_env())
    monkeypatch.setattr(service, 'store', data)
    cell, chrom = 'TTAGGCTAGGCCGGAA-1', 'chr9'
    x = data.x_row(cell).copy()
    labels = data.spliced_clusters(cell, chrom, [169])
    c = data.chrom_by_name[chrom]
    sl = slice(c.offset,c.offset+c.n)
    result = service.post_assign(service.AssignReq(cell=cell,chrom=chrom,boundaries=[169],segment_estimator=method))
    expected = est.assign_with_estimator(x,labels,method)
    arithmetic = sci.assign_gainloss_new(x,labels.tolist())
    np.testing.assert_array_equal(result['manual']['cn5'],expected['cn5'][sl])
    np.testing.assert_array_equal(result['arithmetic_manual']['cn5'],arithmetic[4][sl])
    assert result['manual']['best_s'] == expected['best_s']
    assert result['genome_baseline_trimmed'] == expected['baseline']
    assert result['estimator_changed_bins'] == np.count_nonzero(expected['cn5'][sl]!=arithmetic[4][sl])
    assert result['production']['matches_export']
    assert [(s['start_bin'],s['end_bin']) for s in result['manual']['segments']] == [(0,169),(169,c.n)]
    for seg in result['manual']['segments']:
        local = x[c.offset+seg['start_bin']:c.offset+seg['end_bin']]
        summary = est.segment_summary(local,method)
        assert seg['summary_x'] == pytest.approx(summary['value'])
        assert seg['retained_bins'] == summary['retained']
    np.testing.assert_array_equal(data.x_row(cell),x)
