from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from server.flank_view import summarize, compare_profile


def test_summaries_keep_zeros_and_report_retained_lengths():
    s=summarize(np.array([0.,2,2,2,2,20]))
    assert s['mean']==pytest.approx(28/6)
    assert s['median']==2
    assert s['iqr_mean']==2
    assert s['retained']==4
    assert s['n']==6
    assert s['zero_fraction']==pytest.approx(1/6)
    assert s['std']==pytest.approx(np.std([0,2,2,2,2,20],ddof=1))


def test_internal_segments_only_both_flanks_and_no_pseudocount():
    x=np.array([2,2,4,4,0,0,1,1],dtype=float)
    before=x.copy()
    rows=compare_profile(x,[2,4,6],max_bins=2)
    assert [(r['start'],r['end']) for r in rows]==[(2,4),(4,6)]
    assert rows[0]['left']['effects']['mean']==1
    assert rows[0]['right']['effects']['mean'] is None
    assert rows[1]['left']['effects']['mean'] is None
    assert rows[0]['left']['start']==0
    assert rows[0]['right']['end']==6
    assert compare_profile(x,[2,4,6],max_bins=1)==[]
    np.testing.assert_array_equal(x,before)


@pytest.mark.parametrize('cuts',[[0,3],[2,8],[-1],[2.5]])
def test_invalid_boundary_is_rejected(cuts):
    with pytest.raises(ValueError):compare_profile(np.ones(8),cuts,10)


def test_cache_identity_changes_with_depth_and_manifest(tmp_path):
    from types import SimpleNamespace
    from server.integration import IntegrationStore
    s=IntegrationStore.__new__(IntegrationStore)
    s.cfg=SimpleNamespace(derived_dir=str(tmp_path),prototype_dir=str(tmp_path),cluster_breakpoints=str(tmp_path/'bp.tsv'))
    s.h5mu_path=str(tmp_path/'m.h5mu');s.h5ad_path=str(tmp_path/'x.h5ad')
    for p in [s.h5mu_path,s.h5ad_path,s.cfg.cluster_breakpoints]:Path(p).write_text('first')
    before=s._identity()['hash']
    Path(s.h5ad_path).write_text('different input depth')
    assert s._identity()['hash']!=before
    before=s._identity()['hash']
    (tmp_path/'run_manifest.json').write_text('{"revision":2}')
    assert s._identity()['hash']!=before


def test_paired_cell_bootstrap_preserves_shared_cell_effect():
    from server.flank_view import paired_bootstrap
    # Every cell has a 2x middle/left contrast; shared cell depth varies greatly.
    cells=np.array([[1,2,4],[10,20,40],[100,200,400]],dtype=float)
    r=paired_bootstrap(cells,replicates=200,seed=42)
    assert r['left']['observed']==1
    assert r['right']['observed']==-1
    assert r['left']['ci95']==[1,1]
    assert r['left']['valid']==200
    assert r['left']['se']==0
    assert r==paired_bootstrap(cells,replicates=200,seed=42)


def test_bootstrap_zero_flank_is_undefined_not_regularized():
    from server.flank_view import paired_bootstrap
    r=paired_bootstrap(np.array([[0,1,2],[0,2,4]],dtype=float),20,42)
    assert r['left']['observed'] is None
    assert r['left']['valid']==0
    assert r['left']['ci95'] is None
