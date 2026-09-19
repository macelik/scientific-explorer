"""Live API checks against independent source-matrix calculations."""
from pathlib import Path
import json,urllib.request
import h5py,numpy as np
from scipy.sparse import csr_matrix
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'scientific-explorer/validation/integration-review'

def post(path,body):
    req=urllib.request.Request('http://127.0.0.1:8766/api/integration/'+path,data=json.dumps(body).encode(),headers={'Content-Type':'application/json'})
    with urllib.request.urlopen(req,timeout=90) as r:return json.load(r)

with h5py.File(ROOT/'integration-story/data/count_matrix.h5ad','r') as f:
    g=f['X']; X=csr_matrix((g['data'][...],g['indices'][...],g['indptr'][...]),shape=tuple(g.attrs['shape']))
    v=f['var'];seq=v['seq/categories'][...].astype(str)[v['seq/codes'][...]]
with h5py.File(ROOT/'integration-story/data/integrated.h5mu','r') as f:
    g=f['obs/wnn_leiden_0.3'];labels=g['categories'][...].astype(str)[g['codes'][...]]
count=0
for chrom in ['chr8','chr9']:
    columns=np.flatnonzero(seq==chrom)
    result=post('flank-view',dict(chrom=chrom,sources=['cluster0','cluster1','cluster2','cluster3','cluster4'],max_bins=100))
    for k in range(5):
        pb=np.asarray(X[labels==str(k)][:,columns].mean(axis=0)).ravel()
        np.testing.assert_allclose(result['profiles'][f'cluster{k}'],pb,rtol=1e-6,atol=1e-6)
        for row in [r for r in result['rows'] if r['source']==f'cluster{k}']:
            sample=pb[row['start']:row['end']]
            assert abs(row['summary']['mean']-sample.mean())<1e-5
            for side in ['left','right']:
                flank=pb[row[side]['start']:row[side]['end']]
                expected=np.log2(sample.mean()/flank.mean()) if sample.mean()>0 and flank.mean()>0 else None
                if expected is None:assert row[side]['effects']['mean'] is None
                else:assert abs(row[side]['effects']['mean']-expected)<1e-5
            count+=1
    if chrom=='chr8':
        row=next(r for r in result['rows'] if r['source']=='cluster4')
        b=post('flank-bootstrap',dict(source='cluster4',chrom=chrom,start=row['start'],end=row['end']))
        levels=np.column_stack([np.asarray(X[labels=='4'][:,columns[a:z]].mean(axis=1)).ravel() for a,z in b['region_bounds']])
        assert b['n_cells']==182
        sampled=levels[np.random.default_rng(42).integers(0,182,(200,182))].mean(axis=1)
        for side,j in [('left',0),('right',2)]:
            effect=np.log2(sampled[:,1]/sampled[:,j]);np.testing.assert_allclose(b[side]['effects'],effect,rtol=1e-10,atol=1e-10)
        (OUT/'real-bootstrap-reference.json').write_text(json.dumps(b,indent=2))
(OUT/'numerical-results.json').write_text(json.dumps(dict(chromosomes=['chr8','chr9'],profiles_checked=10,internal_segments_checked=count,bootstrap_cells=182,paired_replicates=200,passed=True),indent=2))
print(f'PASS: 10 pseudobulk profiles, {count} internal segments, independent 182-cell paired bootstrap')
