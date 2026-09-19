"""Descriptive flank comparisons from fixed exported boundaries, not a caller."""
import hashlib
from pathlib import Path
import numpy as np


def summarize(values):
    x=np.asarray(values,dtype=float)
    if not len(x) or np.any(~np.isfinite(x)) or np.any(x<0):
        raise ValueError('Expected nonempty finite nonnegative pseudobulk values')
    q1,q3=np.percentile(x,[25,75]); iqr=q3-q1
    keep=x[(x>=q1-1.5*iqr)&(x<=q3+1.5*iqr)]
    mean=float(x.mean()); std=float(x.std(ddof=1)) if len(x)>1 else None
    return dict(n=len(x),mean=mean,median=float(np.median(x)),iqr_mean=float(keep.mean()),
                retained=len(keep),std=std,cv=std/mean if mean>0 and std is not None else None,
                zero_fraction=float(np.mean(x==0)))


def compare_profile(values,boundaries,max_bins):
    x=np.asarray(values,dtype=float)
    if any(not np.isfinite(b) or int(b)!=b or not 0<b<len(x) for b in boundaries):
        raise ValueError('Boundaries must be integer offsets strictly inside the chromosome')
    cuts=[0,*sorted(set(int(b) for b in boundaries)),len(x)]
    intervals=[dict(start=a,end=b,summary=summarize(x[a:b])) for a,b in zip(cuts[:-1],cuts[1:])]
    rows=[]
    for i in range(1,len(intervals)-1):
        seg=intervals[i]
        if seg['end']-seg['start']>max_bins:continue
        row=dict(seg)
        for name,j in [('left',i-1),('right',i+1)]:
            flank=intervals[j]; effects={}
            for metric in ['mean','median','iqr_mean']:
                a,b=seg['summary'][metric],flank['summary'][metric]
                effects[metric]=float(np.log2(a/b)) if a>0 and b>0 else None
            row[name]={**flank,'effects':effects}
        rows.append(row)
    return rows


def generate_view(store,chrom,sources,max_bins):
    if store.arrays is None:raise ValueError('Pseudobulk index is not ready')
    idx=np.flatnonzero(store.var_seq==chrom)
    if not len(idx):raise ValueError('Unknown chromosome')
    if not sources or any(s not in store.members for s in sources):raise ValueError('Unknown or empty source selection')
    if not np.array_equal(idx,np.arange(idx[0],idx[-1]+1)):raise ValueError('Chromosome bins are not contiguous')
    rows=[]; profiles={}
    for source in dict.fromkeys(sources):
        si=list(store.members).index(source)
        profile=store.arrays['pb_x'][si,idx]
        profiles[source]=profile.tolist()
        bp=store.tables['breakpoints']
        cuts=bp.loc[(bp.source==source)&(bp.chromosome==chrom),'absolute_bin'].tolist()
        for row in compare_profile(profile,cuts,max_bins):
            rows.append(dict(row,key=f"{source}|{chrom}|{row['start']}",source=source,
                             start_bp=int(store.var_start[idx[row['start']]]),end_bp=int(store.var_end[idx[row['end']-1]])))
    return dict(chrom=chrom,sources=list(dict.fromkeys(sources)),max_bins=max_bins,rows=rows,profiles=profiles,
                start_bp=store.var_start[idx].tolist(),end_bp=store.var_end[idx].tolist(),
                provenance=dict(dataset=store.identity,implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                aggregation='Mean of normalized X over original WNN res-0.3 members; float32 cached profile',
                                boundaries='Fixed exported accepted breakpoints; chromosome-local [start,end)',
                                iqr='Two-sided Q1-1.5*IQR through Q3+1.5*IQR, inclusive',
                                uncertainty='Per-bin sample standard deviation (ddof=1), not uncertainty of a cluster estimate',
                                zero_policy='No pseudocount; log2 ratio undefined if either level is zero',
                                scope='Descriptive recomputation; no segmentation, CN, SRD, bootstrap or merge decisions'))


def paired_bootstrap(cell_levels,replicates=200,seed=42):
    """Paired resampling of member cells; columns left, segment, right means."""
    x=np.asarray(cell_levels,dtype=float)
    if x.ndim!=2 or x.shape[1]!=3 or not len(x) or np.any(~np.isfinite(x)) or np.any(x<0):
        raise ValueError('Expected member-cell mean X for three regions')
    rng=np.random.default_rng(seed)
    samples=x[rng.integers(0,len(x),size=(replicates,len(x)))].mean(axis=1)
    original=x.mean(axis=0);out={}
    for side,j in [('left',0),('right',2)]:
        ok=(samples[:,1]>0)&(samples[:,j]>0)
        effects=np.log2(samples[ok,1]/samples[ok,j])
        observed=float(np.log2(original[1]/original[j])) if original[1]>0 and original[j]>0 else None
        out[side]=dict(observed=observed,valid=len(effects),invalid=replicates-len(effects),
                       ci95=np.percentile(effects,[2.5,97.5]).tolist() if len(effects) else None,
                       se=float(np.std(effects,ddof=1)) if len(effects)>1 else None,
                       sign_agreement=float(np.mean(np.sign(effects)==np.sign(observed))) if len(effects) and observed is not None else None,
                       effects=effects.tolist())
    return dict(out,n_cells=len(x),replicates=replicates,seed=seed,estimator='Arithmetic mean only',
                interpretation='Paired whole-cell percentile bootstrap with fixed selected boundaries and membership; exploratory, conditional on this selection, not a calibrated significance test or a merge rule')


def bootstrap_segment(store,source,chrom,start,end):
    import h5py
    if source not in store.members:raise ValueError('Unknown source')
    idx=np.flatnonzero(store.var_seq==chrom)
    if not len(idx):raise ValueError('Unknown chromosome')
    bp=store.tables['breakpoints']
    cuts=[0,*sorted(set(bp.loc[(bp.source==source)&(bp.chromosome==chrom),'absolute_bin'].astype(int))),len(idx)]
    pairs=list(zip(cuts[:-1],cuts[1:]))
    if (start,end) not in pairs:raise ValueError('Selection must match an accepted segment')
    i=pairs.index((start,end))
    if i==0 or i==len(pairs)-1:raise ValueError('Both adjacent flanks are required')
    bounds=pairs[i-1:i+2];members=store.members[source]
    levels=np.zeros((len(members),3))
    with h5py.File(store.h5ad_path,'r') as f:
        g=f['X']; ptr=g['indptr'][...]
        for r,cell in enumerate(members):
            lo,hi=int(ptr[cell]),int(ptr[cell+1])
            columns=g['indices'][lo:hi];values=g['data'][lo:hi]
            for j,(a,b) in enumerate(bounds):
                levels[r,j]=values[(columns>=idx[a])&(columns<=idx[b-1])].sum(dtype=np.float64)/(b-a)
    return dict(paired_bootstrap(levels),source=source,chrom=chrom,start=start,end=end,
                region_bounds=bounds,dataset=store.identity,
                implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
