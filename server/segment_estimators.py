"""Experimental segment summaries; production scientific modules stay unchanged.

Only the per-segment location estimate is exchanged. Genome normalization,
bin-count weighting, scale search, clipping and state mapping retain production
semantics. The arithmetic path delegates directly to assign_gainloss_new.
"""
import numpy as np

from . import pyepi_adapter as sci
from . import science

METHODS = ('arithmetic', 'iqr_upper', 'iqr_two_sided')


def segment_summary(values, method):
    x = np.asarray(values, dtype=float)
    if method not in METHODS:
        raise ValueError('Unknown segment estimator')
    if not len(x) or np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError('Segment requires finite nonnegative X values')
    if method == 'arithmetic':
        return {'value': float(np.mean(x)), 'retained': len(x), 'total': len(x)}
    q1, q3 = np.percentile(x, [25,75])
    lower = q1 - 1.5*(q3-q1) if method == 'iqr_two_sided' else 0
    upper = q3 + 1.5*(q3-q1)
    retained = int(np.count_nonzero((x >= lower) & (x <= upper)))
    return {'value': float(sci.trimmed_mean_iqr(x, lb=method == 'iqr_two_sided')),
            'retained': retained, 'total': len(x), 'lower': float(lower), 'upper': float(upper)}


def assign_with_estimator(xfull, labels, method='arithmetic'):
    if method not in METHODS:
        raise ValueError('Unknown segment estimator')
    xfull, labels = np.asarray(xfull,dtype=float), np.asarray(labels)
    if len(xfull) != len(labels) or np.any(~np.isfinite(xfull)) or np.any(xfull < 0):
        raise ValueError('Finite nonnegative X and matching segment labels are required')
    baseline = float(sci.trimmed_mean_iqr(xfull, lb=False))
    if not np.isfinite(baseline) or baseline <= 0:
        raise ValueError('CN fitting requires a positive full-genome baseline')
    ids, inverse, lengths = np.unique(labels,return_inverse=True,return_counts=True)
    normalized = xfull/baseline
    info = {int(i): segment_summary(normalized[labels==i],method) for i in ids}
    means = np.array([info[int(i)]['value'] for i in ids])
    if method == 'arithmetic':
        result = science.assign_cn(xfull, labels.tolist())
    else:
        median = float(np.median(means[inverse]))
        if not np.isfinite(median) or median <= 0:
            raise ValueError('IQR segment estimates do not yield a positive scale grid; CN is undefined for this option')
        # Use original lengths, not the number of bins retained by trimming.
        best_s, continuous, integer = sci.weighted_scale_search(means.copy(), lengths, median)
        ci, cc = integer[inverse], continuous[inverse]
        holmes = np.clip(ci,1,3)-1
        watson = np.where(cc<=1,0,np.where(cc>=3,2,1))
        combined = np.full(len(xfull),-1.0)
        for h,w,value in [(0,0,0),(1,1,1),(2,2,2),(0,1,0.5),(2,1,1.5)]:
            combined[(holmes==h)&(watson==w)] = value
        result = {'cn_int':ci, 'cn_cont':cc, 'cn5':combined, 'best_s':float(best_s)}
    result.update(baseline=baseline,segment_estimator=method,segment_info=info)
    return result
