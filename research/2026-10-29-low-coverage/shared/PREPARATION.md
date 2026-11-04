# Deterministic panel preparation

Command from repository root: `python3 /tmp/prepare_low_coverage_research_panel.py`

No random sampling or statistical inference was performed. All 300 selected cells were checked for exact AT = HA + HB and matching bin coordinates. Exported regional sums do not establish phase coherence or common depth/AT counting units.

```python
from pathlib import Path
import hashlib, json, platform, subprocess
import h5py, numpy as np, pandas as pd, pyarrow, pyarrow.parquet as pq

root = Path.cwd()
out = root / 'research/2026-10-29-low-coverage/shared'
out.mkdir(parents=True, exist_ok=True)
depth = (root / '../integration-story/data/count_matrix.h5ad').resolve()
mu = (root / '../integration-story/data/integrated.h5mu').resolve()
cohort = (root / '../default_300cells/cells.csv').resolve()
allelic = Path('/home/mami/Desktop/EGA/haplotype_aware_cnv/hapcnv_outs_from_updated_phasing/aligned_matrices')
segments = (root / '../integration-story/data/genome_pool/derived/karyogram_cluster_segments.tsv').resolve()

def identity(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    st = path.stat()
    return {'path': str(path), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns, 'sha256': h.hexdigest()}

def cat(group):
    if isinstance(group, h5py.Dataset):
        return group.asstr()[:]
    return group['categories'].asstr()[:][group['codes'][:]]

def dense_rows(group, indices, n_bins):
    arr = np.zeros((len(indices), n_bins), dtype=np.float64)
    ptr = group['indptr'][:]
    for i, row in enumerate(indices):
        a, b = int(ptr[row]), int(ptr[row + 1])
        arr[i, group['indices'][a:b]] = group['data'][a:b]
    return arr

cells = pd.read_csv(cohort)
assert len(cells) == 300 and cells.cellID.is_unique
ids = cells.cellID.astype(str).tolist()
with h5py.File(depth, 'r') as f, h5py.File(mu, 'r') as m:
    all_ids = f['obs/cellID'].asstr()[:]
    assert np.array_equal(all_ids, m['obs/_index'].asstr()[:])
    index = {v: i for i, v in enumerate(all_ids)}
    rows = np.array([index[c] for c in ids])
    seq, start, end = cat(f['var/seq']), f['var/start'][:], f['var/end'][:]
    n = len(start)
    raw = dense_rows(f['layers/raw'], rows, n)
    processed = dense_rows(f['X'], rows, n)
    labels = cat(m['obs/wnn_leiden_0.3'])
    cell_table = pd.DataFrame({'cell': ids, 'h5ad_row': rows, 'cluster': labels[rows],
                              'coverage_tercile': cells.tercile.astype(str),
                              'raw_total_retained_bins': raw.sum(axis=1),
                              'metadata_total_raw_depth': m['obs/total_raw_depth'][:][rows],
                              'metadata_snp_coverage': m['obs/snp_coverage'][:][rows]})
    assert np.all(raw >= 0) and np.all(raw == np.floor(raw))
mat = {}
for name in ['HA', 'HB', 'AT']:
    coords = pq.read_table(allelic / (name + '.parquet'), columns=['chr', 'start', 'end'])
    assert np.array_equal(coords['chr'].to_numpy(), seq)
    assert np.array_equal(coords['start'].to_numpy(), start)
    assert np.array_equal(coords['end'].to_numpy(), end)
    mat[name] = pq.read_table(allelic / (name + '.parquet'), columns=ids).to_pandas()[ids].to_numpy(dtype=np.float64).T
    assert mat[name].shape == raw.shape
    assert np.all(mat[name] >= 0) and np.all(mat[name] == np.floor(mat[name]))
assert np.array_equal(mat['AT'], mat['HA'] + mat['HB'])
cell_table['at_total_retained_bins'] = mat['AT'].sum(axis=1)
cell_table.to_csv(out / 'cells.tsv', sep='\t', index=False)

windows = [('chr7_8.3_10.7','chr7',8.3,10.7,'candidate'),
           ('chr7_51.4_53.5','chr7',51.4,53.5,'candidate'),
           ('chr7_67.4_70.0','chr7',67.4,70.0,'candidate'),
           ('chr7_70.0_73.3','chr7',70.0,73.3,'candidate'),
           ('chr9_20_27','chr9',20,27,'candidate_broad_window'),
           ('chr9_22.5_26.6','chr9',22.5,26.6,'candidate'),
           ('chr9_70_100','chr9',70,100,'historical_comparison_not_truth')]
region_rows, count_rows, bin_rows = [], [], []
for name, chrom, lo, hi, role in windows:
    wstart, wend = int(round(lo * 1e6)) + 1, int(round(hi * 1e6))
    idx = np.flatnonzero((seq == chrom) & (start <= wend) & (end >= wstart))
    assert len(idx)
    left = np.flatnonzero((seq == chrom) & (end < wstart))[-50:]
    right = np.flatnonzero((seq == chrom) & (start > wend))[:50]
    for suffix, ix in [('target',idx),('left_flank',left),('right_flank',right)]:
        rid = name + '__' + suffix
        if not len(ix): continue
        # Flanks are 50 retained bins, not necessarily a contiguous 5 Mb window.
        total_span = int(end[ix[-1]] - start[ix[0]] + 1)
        retained_bp = int(np.sum(end[ix] - start[ix] + 1))
        region_rows.append({'region_id':rid,'parent':name,'part':suffix,'role':role,
                            'chrom':chrom,'requested_start_1based':wstart if suffix=='target' else int(start[ix[0]]),
                            'requested_end_inclusive':wend if suffix=='target' else int(end[ix[-1]]),
                            'first_retained_start':int(start[ix[0]]),'last_retained_end':int(end[ix[-1]]),
                            'n_bins':len(ix),'retained_bp':retained_bp,'span_bp':total_span,
                            'gap_bp_within_retained_span':total_span-retained_bp})
        for j in ix:
            bin_rows.append({'region_id':rid,'bin_index':int(j),'chrom':chrom,'start':int(start[j]),'end':int(end[j])})
        rs, xs = raw[:,ix], processed[:,ix]
        aa, bb, tt = mat['HA'][:,ix], mat['HB'][:,ix], mat['AT'][:,ix]
        for i,c in enumerate(ids):
            count_rows.append({'cell':c,'cluster':str(labels[rows[i]]),'region_id':rid,'n_bins':len(ix),
                              'raw_sum':float(rs[i].sum()),'raw_mean':float(rs[i].mean()),
                              'raw_zero_bins':int((rs[i]==0).sum()),
                              'processed_sum':float(xs[i].sum()),'processed_mean':float(xs[i].mean()),
                              'processed_zero_bins':int((xs[i]==0).sum()),
                              'ha_sum':float(aa[i].sum()),'hb_sum':float(bb[i].sum()),'at_sum':float(tt[i].sum()),
                              'at_zero_bins':int((tt[i]==0).sum()),'at_ge2_bins':int((tt[i]>=2).sum()),
                              'at_ge10_bins':int((tt[i]>=10).sum())})
pd.DataFrame(region_rows).to_csv(out/'regions.tsv',sep='\t',index=False)
pd.DataFrame(bin_rows).to_csv(out/'region_bins.tsv',sep='\t',index=False)
pd.DataFrame(count_rows).to_csv(out/'region_counts.tsv',sep='\t',index=False)
selected = pd.read_csv(segments,sep='\t')
selected[selected.chromosome.isin(['chr7','chr9'])].to_csv(out/'existing_chr7_chr9_segments.tsv',sep='\t',index=False)
inputs = [depth,mu,cohort,segments]+[allelic/(x+'.parquet') for x in ['HA','HB','AT']]
manifest = {'date':'2026-10-29','sample':'MCF10-A; user supplied identity; not region-level CN truth',
            'repo_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            'repo_status':subprocess.check_output(['git','status','--short'],text=True).splitlines(),
            'panel':'Existing coverage-stratified 300 cells; legacy sampling seed unverified; exact membership frozen',
            'new_randomness':'None; these are deterministic descriptive exports',
            'coordinate_convention':'1-based inclusive; window starts converted from Mb boundaries by adding 1',
            'flanks':'50 retained bins on each side on same chromosome; may overlap other candidate windows',
            'allelic_caution':'Sums are descriptive; phase continuity and counting units require upstream audit; do not treat sums as LOH tests',
            'checks':{'unique_cells':True,'barcode_alignment':True,'coordinates_match_all_matrices':True,
                      'AT_equals_HA_plus_HB_all_300':True,'raw_nonnegative_integer':True,'no_refiltering':True},
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
                           'h5py':h5py.__version__,'pyarrow':pyarrow.__version__},
            'inputs':[identity(p) for p in inputs],
            'artifacts':[identity(p) for p in sorted(out.glob('*.tsv'))]}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
src=Path('/tmp/prepare_low_coverage_research_panel.py').read_text()
(out/'PREPARATION.md').write_text('# Deterministic panel preparation\n\nCommand from repository root: `python3 /tmp/prepare_low_coverage_research_panel.py`\n\nNo random sampling or statistical inference was performed. All 300 selected cells were checked for exact AT = HA + HB and matching bin coordinates. Exported regional sums do not establish phase coherence or common depth/AT counting units.\n\n```python\n'+src+'\n```\n')
print(json.dumps({'cells':len(ids),'regions_with_flanks':len(region_rows),'cell_region_rows':len(count_rows),
                  'cluster_sizes':cell_table.cluster.value_counts().sort_index().to_dict(),
                  'checks':manifest['checks'],'output_dir':str(out)},indent=2))

```
