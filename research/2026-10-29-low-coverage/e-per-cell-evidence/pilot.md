# Shared-panel descriptive pilot

Question: how often are whole-region and individual bins zero, how sparse is AT, and how do raw and processed target densities compare with both flanks? No inferential model, p-value, random sampling, allele fraction or power calculation was used.

Inputs: coordinator-frozen `shared/region_counts.tsv` (SHA-256 ee9113acea13baa3c9ffd8da176b2a09c11b0521437e5fa460751febbc9ea718), `regions.tsv` (e38857779ca98bb9d775e9d404942102fa7b727d995784de1b4a18fbf4eb239c), `cells.tsv` (10a37543b38f79782436f434e8844eae294929d66815e0cebc90fa6b33cc23c4); see shared manifest for original provenance and fixed legacy membership. Seed: none, deterministic summary. Environment: Python 3.10.12, pandas 2.2.0, numpy 1.24.4. Pandas reported older optional numexpr and bottleneck versions; operations completed.

From repository root: `python3 /tmp/e_pilot_summary.py`. Outputs: `artifacts/zero_support_by_stratum.tsv`, `artifacts/flank_density_comparison.tsv`. Raw and processed densities are per retained megabase; target/flank is the ratio of 300-cell mean densities. AT is site-observation count, potentially repeated molecules, and cannot establish independent allelic trials; phase continuity unknown. Flanks may overlap other candidate regions, and findings do not control selection, accessibility, or coverage confounding.

Exact source:

```python
import pandas as pd, numpy as np, platform, hashlib
from pathlib import Path
root=Path('research/2026-10-29-low-coverage')
s=root/'shared'; out=root/'e-per-cell-evidence/artifacts'
d=pd.read_csv(s/'region_counts.tsv',sep='\t'); c=pd.read_csv(s/'cells.tsv',sep='\t'); reg=pd.read_csv(s/'regions.tsv',sep='\t')
d=d.merge(c[['cell','coverage_tercile']],on='cell',validate='many_to_one').merge(reg[['region_id','parent','part','role','retained_bp']],on='region_id',validate='many_to_one')
d=d[(d.role.isin(['candidate','candidate_broad_window'])) & (d.part=='target')].copy()
rows=[]
for (region,kind,level),g in pd.concat([d.assign(kind='all',level='all'),d.assign(kind='coverage_tercile',level=d.coverage_tercile.astype(str)),d.assign(kind='cluster',level=d.cluster.astype(str))]).groupby(['parent','kind','level'],sort=True):
 a=g.at_sum
 rows.append(dict(region=region,stratum=kind,level=level,n_cells=len(g),n_bins=int(g.n_bins.iloc[0]),raw_whole_zero_cells=int((g.raw_sum==0).sum()),raw_zero_bin_fraction=g.raw_zero_bins.sum()/g.n_bins.sum(),at_whole_zero_cells=int((a==0).sum()),at_zero_bin_fraction=g.at_zero_bins.sum()/g.n_bins.sum(),at_sum_q25=a.quantile(.25),at_sum_median=a.median(),at_sum_q75=a.quantile(.75),at_sum_max=a.max(),cells_at_ge10=int((a>=10).sum())))
pd.DataFrame(rows).to_csv(out/'zero_support_by_stratum.tsv',sep='\t',index=False,float_format='%.4f')
x=pd.read_csv(s/'region_counts.tsv',sep='\t').merge(reg[['region_id','parent','part','role','retained_bp']],on='region_id',validate='many_to_one'); x=x[x.role.isin(['candidate','candidate_broad_window'])].copy()
# Sum across cells divided by retained bp; ratio of densities, not mean of per-cell ratios.
x['raw_per_mb']=x.raw_sum/x.retained_bp*1e6; x['processed_per_mb']=x.processed_sum/x.retained_bp*1e6
z=x.groupby(['parent','part'])[['raw_per_mb','processed_per_mb']].mean().unstack('part'); rows=[]
for region in z.index:
 row={'region':region}
 for metric in ['raw_per_mb','processed_per_mb']:
  for part in ['target','left_flank','right_flank']: row[metric+'_'+part]=z.loc[region,(metric,part)]
  for part in ['left_flank','right_flank']: row[metric+'_target_over_'+part]=z.loc[region,(metric,'target')]/z.loc[region,(metric,part)]
 rows.append(row)
pd.DataFrame(rows).to_csv(out/'flank_density_comparison.tsv',sep='\t',index=False,float_format='%.4f')
print('Python',platform.python_version(),'pandas',pd.__version__,'numpy',np.__version__)
for f in ['region_counts.tsv','regions.tsv','cells.tsv']:
 print(f,hashlib.sha256((s/f).read_bytes()).hexdigest())
print(pd.DataFrame(rows).to_string(index=False))
print(pd.DataFrame([r for r in pd.read_csv(out/'zero_support_by_stratum.tsv',sep='\t').to_dict('records') if r['stratum']=='all']).to_string(index=False))
```
