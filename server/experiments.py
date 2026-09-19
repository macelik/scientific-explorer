"""Immutable interventions on normalized X, using the supplied scientific functions."""
from __future__ import annotations

import time
import numpy as np

from . import science
from . import pyepi_adapter as sci

PARAMETERS = {'k': 2, 'n_permutations': 1000, 'seed': 42, 'alpha': 0.001}


def apply_edit(source, start, end, operation, value, preserve_total):
    source = np.asarray(source, dtype=np.float64)
    if source.ndim != 1 or not np.all(np.isfinite(source)) or np.any(source < 0):
        raise ValueError('Source X must be a finite nonnegative vector')
    if not (isinstance(start, (int, np.integer)) and isinstance(end, (int, np.integer))
            and 0 <= start < end <= len(source)):
        raise ValueError('Select a nonempty half-open interval inside the chromosome')
    if not np.isfinite(value):
        raise ValueError('The edit value must be finite')
    if operation not in ('multiply', 'add', 'set'):
        raise ValueError('Unknown operation')
    if operation in ('multiply', 'set') and value < 0:
        raise ValueError('Multiplication factors and replacement values must be nonnegative')
    edited = source.copy()
    with np.errstate(over='ignore', invalid='ignore'):
        region = source[start:end] * value if operation == 'multiply' else (
            source[start:end] + value if operation == 'add' else np.full(end-start, value))
    if not np.all(np.isfinite(region)):
        raise ValueError('This edit produces nonfinite X values')
    clipped = int(np.sum(region < 0))
    edited[start:end] = np.maximum(region, 0)
    total_before, total_after = float(source.sum()), float(edited.sum())
    if not np.isfinite(total_after):
        raise ValueError('The edited genome total is not finite')
    factor = 1.0
    if preserve_total:
        if total_after <= 0 or total_before <= 0:
            raise ValueError('Cannot preserve genome total when the original or edited total is zero')
        factor = total_before / total_after
        edited *= factor
    if not np.all(np.isfinite(edited)):
        raise ValueError('Normalization produces nonfinite X values')
    return edited, {
        'operation': operation, 'value': float(value), 'preserve_total': bool(preserve_total),
        'clipped_bins': clipped, 'normalization_factor': factor,
        'original_total': total_before, 'edited_total': float(edited.sum()),
        'total_before_normalization': total_after,
        'changed_bins_genome': int(np.count_nonzero(edited != source)),
    }


def apply_region_edit(source, offset, n, request, production_cn=None):
    """Extend/duplicate in-place on the retained grid; never insert coordinates."""
    start, end = request['start'], request['end']
    if not (isinstance(start, int) and isinstance(end, int) and 0 <= start < end <= n):
        raise ValueError('Select a nonempty region inside the chromosome')
    operation = request['operation']
    if operation in ('multiply', 'add', 'set'):
        edited, info = apply_edit(source, offset+start, offset+end, operation,
                                 request['value'], request.get('preserve_total', False))
        return edited, {**info, 'event_start': start, 'event_end': end,
                        'target_start': start, 'target_end': end}
    source = np.asarray(source, dtype=float)
    if source.ndim != 1 or np.any(~np.isfinite(source)) or np.any(source < 0):
        raise ValueError('Source X must be finite and nonnegative')
    rng = np.random.default_rng(request.get('seed', 42))
    edited = source.copy()
    info = {'operation': operation, 'clipped_bins': 0, 'event_start': start, 'event_end': end}
    donor = source[offset+start:offset+end]
    if operation == 'simulate':
        state = {'loss': 0, 'base': 1, 'gain': 2}.get(request.get('cn_state'))
        if state is None or production_cn is None:
            raise ValueError('Select an available production loss/base/gain donor class')
        donor = source[np.asarray(production_cn) == state]
        if len(donor) < 8:
            raise ValueError('This cell has fewer than eight donor bins in the selected exact production class')
        edited[offset+start:offset+end] = rng.choice(donor, end-start, replace=True)
        info.update(target_start=start, target_end=end, donor_class=request['cn_state'])
    elif operation in ('extend', 'duplicate'):
        amount, direction = request.get('length', end-start), request.get('direction', 'right')
        if not isinstance(amount, int) or amount < 1 or direction not in ('left', 'right'):
            raise ValueError('Extension requires a positive integer length and left/right direction')
        lo, hi = (max(0, start-amount), start) if direction == 'left' else (end, min(n, end+amount))
        actual = hi-lo
        if not actual:
            raise ValueError('No retained bins are available in that direction')
        if operation == 'extend':
            values = rng.choice(donor, actual, replace=True)
        else:
            indices = np.arange(-actual, 0) if direction == 'left' else np.arange(actual)
            values = donor[indices % len(donor)]
        edited[offset+lo:offset+hi] = values
        info.update(target_start=lo, target_end=hi, event_start=min(start,lo), event_end=max(end,hi),
                    requested_extension=amount, actual_extension=actual, direction=direction)
    else:
        raise ValueError('Unknown operation')
    info.update(donor_count=len(donor), donor_mean=float(np.mean(donor)),
                donor_zero_fraction=float(np.mean(donor == 0)), seed=request.get('seed', 42),
                sampling=('Exact ordered cyclic copy' if operation == 'duplicate' else
                          'Independent empirical resampling; spatial correlation is not preserved'))
    before, after = float(source.sum()), float(edited.sum())
    factor = 1.0
    if request.get('preserve_total'):
        if before <= 0 or after <= 0:
            raise ValueError('Cannot preserve total when the original or edited total is zero')
        factor = before/after
        edited *= factor
    if not np.all(np.isfinite(edited)):
        raise ValueError('The edit produced nonfinite values')
    info.update(normalization_factor=factor, original_total=before, edited_total=float(edited.sum()),
                changed_bins_genome=int(np.count_nonzero(edited != source)),
                preserve_total=request.get('preserve_total', False))
    return edited, info


def inspect_tree(xfull, offset, n):
    """Root and both children; only gated k=2 splits enter the final segmentation."""
    nodes, boundaries = [], []

    def inspect(start, end, side, parent_accepted):
        x = xfull[offset+start:offset+end]
        node = science.compute_node(x)
        node.update(start=start, end=end, side=side, accepted=False,
                    hypothetical=not parent_accepted, tests=None)
        b = node['argmax_b']
        if b is not None:
            tests = science.run_split_tests(x, b, xfull, 1000, 42)
            node['tests'] = tests
            node['passes'] = bool(tests['p_local'] < 0.001 and tests['p_global'] < 0.001)
            node['accepted'] = bool(parent_accepted and node['passes'])
            if node['accepted']:
                boundaries.append(start+b)
        nodes.append(node)
        return node

    root = inspect(0, n, 'root', True)
    if root['argmax_b'] is not None:
        inspect(0, root['argmax_b'], 'left', root['accepted'])
        inspect(root['argmax_b'], n, 'right', root['accepted'])
    return {'nodes': nodes, 'boundaries': sorted(boundaries)}


def splice_labels(production, offset, n, boundaries):
    labels = np.asarray(production).copy()
    next_id = int(labels.max()) + 1
    edges = [0, *sorted(set(boundaries)), n]
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        labels[offset+a:offset+b] = next_id+i
    return labels


def fit_cn(xfull, labels, offset, n, chroms):
    baseline = float(sci.trimmed_mean_iqr(xfull, lb=False))
    if not np.isfinite(baseline) or baseline <= 0:
        return {'status': 'undefined', 'reason': 'Genome trimmed-mean baseline is zero or undefined', 'cn5': None}
    # The original fitter has no valid positive scale grid if the bin-weighted
    # segment median vanishes. Report it rather than returning NaN/inf in JSON.
    values = np.empty(len(xfull))
    for label in np.unique(labels):
        mask = labels == label
        values[mask] = np.mean(xfull[mask]) / baseline
    if np.median(values) <= 0:
        return {'status': 'undefined', 'reason': 'Bin-weighted segment median is zero; scale fitting is undefined', 'cn5': None}
    result = science.assign_cn(xfull, labels.tolist())
    if not np.isfinite(result['best_s']) or result['best_s'] <= 0:
        return {'status': 'undefined', 'reason': 'No finite positive fitted scale', 'cn5': None}
    return {
        'status': 'done', 'best_s': result['best_s'], 'baseline': baseline,
        'cn5': result['cn5'][offset:offset+n].tolist(),
        'cn_int': result['cn_int'][offset:offset+n].tolist(),
        'cn_cont': result['cn_cont'][offset:offset+n].tolist(),
        'holmes': (np.clip(result['cn_int'][offset:offset+n], 1, 3)-1).tolist(),
        'watson': np.where(result['cn_cont'][offset:offset+n] <= 1, 0,
                           np.where(result['cn_cont'][offset:offset+n] >= 3, 2, 1)).tolist(),
        '_genome_cn5': result['cn5'],
    }


def run_experiment(xfull, production, chroms, starts, ends, request, pkg_dir):
    science._ensure(pkg_dir)
    t0 = time.perf_counter()
    chrom = next(c for c in chroms if c['name'] == request['chrom'])
    offset, n = chrom['offset'], chrom['n']
    start, end = request['start'], request['end']
    if not 0 <= start < end <= n:
        raise ValueError('Region bounds are outside the selected chromosome')
    xfull = np.asarray(xfull, dtype=np.float64)
    donor_cn = science.assign_cn(xfull, production.tolist())['cn5'] if request['operation'] == 'simulate' else None
    edited, intervention = apply_region_edit(xfull, offset, n, request, donor_cn)
    original_tree = inspect_tree(xfull, offset, n)
    edited_tree = inspect_tree(edited, offset, n)
    root_b = original_tree['nodes'][0]['argmax_b']
    same_boundary = (science.run_split_tests(edited[offset:offset+n], root_b, edited, 1000, 42)
                     if root_b is not None else None)
    edited_labels = splice_labels(production, offset, n, edited_tree['boundaries'])
    event_start, event_end = intervention['event_start'], intervention['event_end']
    prod_bounds = (np.where(np.diff(production[offset:offset+n]) != 0)[0]+1).tolist()
    event_bounds = [b for b in [event_start, event_end] if 0 < b < n]
    isolated = splice_labels(production, offset, n, prod_bounds + event_bounds)
    cn = {
        'original': fit_cn(xfull, production, offset, n, chroms),
        'fixed': fit_cn(edited, production, offset, n, chroms),
        'resegmented': fit_cn(edited, edited_labels, offset, n, chroms),
        'event_isolated': fit_cn(edited, isolated, offset, n, chroms),
    }
    ref = cn['original'].get('_genome_cn5')
    for key, result in cn.items():
        full = result.get('_genome_cn5')
        if full is not None and ref is not None:
            result['changed_bins'] = int(np.count_nonzero(full[offset:offset+n] != ref[offset:offset+n]))
            result['changed_elsewhere'] = [
                {'chrom': c['name'], 'changed_bins': int(np.count_nonzero(
                    full[c['offset']:c['offset']+c['n']] != ref[c['offset']:c['offset']+c['n']]))}
                for c in chroms if c['name'] != request['chrom'] and np.any(
                    full[c['offset']:c['offset']+c['n']] != ref[c['offset']:c['offset']+c['n']])
            ]
    for result in cn.values():
        result.pop('_genome_cn5', None)
        if result['status'] == 'done':
            result['event_calls'] = {kind: {str(v): int(np.count_nonzero(np.asarray(result[kind][event_start:event_end]) == v))
                                          for v in np.unique(result[kind][event_start:event_end])}
                                     for kind in ('cn5', 'holmes', 'watson')}
    edges = []
    for label, b in [('left', event_start), ('right', event_end)]:
        edge = {'edge': label, 'boundary': b, 'chromosome_edge': b in (0,n),
                'accepted': b in edited_tree['boundaries'], 'nodes': []}
        for node in edited_tree['nodes']:
            if node['hypothetical'] or not node['start'] < b < node['end']:
                continue
            local = b-node['start']
            test = science.run_split_tests(edited[offset+node['start']:offset+node['end']], local, edited, 1000, 42)
            rank = 1+int(np.count_nonzero(np.asarray(node['ad']) > test['observed_ad']))
            edge['nodes'].append({'side': node['side'], 'rank': rank, 'is_winner': local == node['argmax_b'], **test})
        edges.append(edge)
    sl = slice(offset, offset+n)
    return {
        'request': request, 'parameters': PARAMETERS, 'intervention': intervention,
        'original': original_tree, 'edited': edited_tree,
        'edited_at_original_boundary': same_boundary,
        'start_bp': np.asarray(starts[sl], dtype=int).tolist(),
        'end_bp': np.asarray(ends[sl], dtype=int).tolist(),
        'x_original': xfull[sl].tolist(), 'x_edited': edited[sl].tolist(),
        'cn': cn, 'compute_ms': round(1000*(time.perf_counter()-t0), 1),
        'event_edges': edges,
        'provenance': {
            'signal': 'GC-corrected, normalized X; not raw fragment counts',
            'null_pool': 'Each scenario uses its own full-genome X for the global test',
            'segmentation_scope': 'Selected chromosome only, k=2; other chromosome boundaries remain production',
            'coordinates': 'Chromosome-local half-open bin interval; genomic bp are 1-based inclusive',
            'source_hash': request.get('dataset_hash'), 'science_hash': request.get('science_hash'),
        },
    }
