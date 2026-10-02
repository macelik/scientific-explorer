"""User-defined merge rules for the exploratory adjacent-segment merge.

A rule is a list of conditions on one internal segment's left/right flank
scores that must ALL hold; a segment is flagged when ANY enabled rule matches
(OR of ANDs). The topmost matching rule claims the segment, and that rule's
mergeInto metric picks the neighbour to merge into (smaller |score|, ties
left). Rows are the per-segment dicts produced by hatch_merge.score_rows, the
same rows /api/integration/hatch-scores returns, so web/src/mergeRules.ts can
evaluate identical rules client-side (see tests/fixtures/merge_rules_parity.json).

These are exploratory diagnostics, independent of the original
bootstrap/chain-guard merge prototype; nothing here recomputes CN.
"""
from __future__ import annotations

import math

from . import hatch_merge as hm

CONDITION_KINDS = ('weak_flank', 'sign_mismatch', 'ambiguous', 'max_size')
SRD_METRICS = ('srd', 'srd_phi')
HATCH_PATTERNS = ('/', '\\', 'x', '-', '|', '.')


def _number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _check_metric(metric, estimator, where):
    if metric not in hm.METRICS:
        raise ValueError(f'{where}: unknown metric {metric!r}')
    if estimator not in hm.ESTIMATORS:
        raise ValueError(f'{where}: unknown estimator {estimator!r}')


def _validate_condition(c, where):
    if not isinstance(c, dict):
        raise ValueError(f'{where}: each condition must be an object')
    kind = c.get('kind')
    if kind == 'weak_flank':
        _check_metric(c.get('metric'), c.get('estimator', 'mean'), where)
        t = c.get('threshold')
        if not (_number(t) and math.isfinite(t) and t > 0):
            raise ValueError(f'{where}: threshold must be a positive finite number')
    elif kind == 'sign_mismatch':
        _check_metric(c.get('metric'), c.get('estimator', 'mean'), where)
    elif kind == 'ambiguous':
        if c.get('srdMetric') not in SRD_METRICS:
            raise ValueError(f'{where}: srdMetric must be one of {SRD_METRICS}')
        if c.get('estimator', 'mean') not in hm.ESTIMATORS:
            raise ValueError(f'{where}: unknown estimator {c.get("estimator")!r}')
    elif kind == 'max_size':
        n = c.get('maxBins')
        if not (isinstance(n, int) and not isinstance(n, bool) and n >= 1):
            raise ValueError(f'{where}: maxBins must be an integer >= 1')
    else:
        raise ValueError(f'{where}: unknown condition kind {kind!r}')


def validate_rules(rules):
    """Validate and return the enabled rules, in list order.

    Disabled rules are ignored apart from the id checks, so a half-edited
    disabled rule never blocks a run. Raises ValueError on any problem."""
    if not isinstance(rules, list) or not all(isinstance(r, dict) for r in rules):
        raise ValueError('rules must be a list of objects')
    ids = [r.get('id') for r in rules]
    if any(not isinstance(i, str) or not i for i in ids):
        raise ValueError('every rule needs a non-empty string id')
    if len(set(ids)) != len(ids):
        raise ValueError('rule ids must be unique')
    enabled = [r for r in rules if r.get('enabled')]
    if not enabled:
        raise ValueError('At least one merge rule must be enabled')
    for r in enabled:
        where = f'rule {r["id"]!r}'
        conditions = r.get('conditions')
        if not isinstance(conditions, list) or not conditions:
            raise ValueError(f'{where}: needs at least one condition')
        for c in conditions:
            _validate_condition(c, where)
        merge_into = r.get('mergeInto') if isinstance(r.get('mergeInto'), dict) else {}
        _check_metric(merge_into.get('metric'), merge_into.get('estimator', 'mean'), f'{where} mergeInto')
        hatch = r.get('hatch') if isinstance(r.get('hatch'), dict) else {}
        if hatch.get('pattern') not in HATCH_PATTERNS:
            raise ValueError(f'{where}: unknown hatch pattern {hatch.get("pattern")!r}')
    return enabled


def row_score(row, side, metric, estimator='mean'):
    """The finite signed flank score for one side of a score row, else None."""
    key = f'delta_{estimator}_{side}' if metric == 'delta_log2fc' else f'{metric}_{side}'
    v = row.get(key)
    return float(v) if v is not None and math.isfinite(v) else None


def condition_holds(c, row):
    kind = c['kind']
    if kind == 'max_size':
        return (row['end'] - row['start']) <= c['maxBins']
    estimator = c.get('estimator', 'mean')
    if kind == 'ambiguous':
        return hm.classify_ambiguous(
            row_score(row, 'left', c['srdMetric']), row_score(row, 'right', c['srdMetric']),
            row_score(row, 'left', 'delta_log2fc', estimator), row_score(row, 'right', 'delta_log2fc', estimator),
        ) == 'ambiguous'
    left = row_score(row, 'left', c['metric'], estimator)
    right = row_score(row, 'right', c['metric'], estimator)
    if kind == 'sign_mismatch':
        return hm.classify_transition(left, right)
    if kind == 'weak_flank':
        return left is not None and right is not None and min(abs(left), abs(right)) < c['threshold']
    raise ValueError(f'unknown condition kind {kind!r}')


def evaluate_rules(rows, rules, allow_gap_crossing=False):
    """Per score row: which rules match, which claims it, and the merge
    direction/score or why it cannot merge. `rules` must be validated."""
    out = []
    for row in rows:
        matched = [r['id'] for r in rules if all(condition_holds(c, row) for c in r['conditions'])]
        ev = {'start': row['start'], 'end': row['end'], 'matched': matched,
              'claimed_by': None, 'side': None, 'score': None, 'blocked': None}
        if matched:
            rule = next(r for r in rules if r['id'] == matched[0])
            ev['claimed_by'] = rule['id']
            metric = rule['mergeInto']['metric']
            estimator = rule['mergeInto'].get('estimator', 'mean')
            left = row_score(row, 'left', metric, estimator)
            right = row_score(row, 'right', metric, estimator)
            if left is None or right is None:
                ev['blocked'] = 'nonfinite'
            else:
                side = 'left' if abs(left) <= abs(right) else 'right'
                ev['side'] = side
                ev['score'] = abs(left) if side == 'left' else abs(right)
                if not allow_gap_crossing and row[f'gap_{side}']:
                    ev['blocked'] = 'gap'
        out.append(ev)
    return out


def run_rule_merge(profile, boundaries, n_bins, var_start, var_end, chrom_offset,
                   rules, allow_gap_crossing=False, max_steps=10000):
    """Greedy exploratory merge driven by user-defined rules, one merge per
    step. Each step recomputes every score (pooled phi included) for the
    current partition, evaluates the enabled rules, and merges the candidate
    with the lowest (claiming-rule index, weaker-flank |score|, start bin).
    Gap-blocked and non-finite candidates are skipped, never redirected.
    Stops when no candidate remains."""
    enabled = validate_rules(rules)
    order = {r['id']: k for k, r in enumerate(enabled)}
    names = {r['id']: (r.get('name') or r['id']) for r in enabled}
    segs = hm.flat_partition(boundaries, n_bins)
    original = [{'start': s, 'end': e} for s, e in segs]
    steps = [{'step': 0, 'segments': list(original), 'removed_boundary': None, 'merged_interval': None,
              'segment': None, 'rule_id': None, 'rule_name': None, 'side': None, 'score': None, 'phi': None}]
    for _ in range(max_steps + 1):
        phi, rows = hm.score_rows(profile, segs, var_start, var_end, chrom_offset)
        evals = evaluate_rules(rows, enabled, allow_gap_crossing)
        candidates = [(k + 1, ev) for k, ev in enumerate(evals)
                      if ev['claimed_by'] is not None and ev['blocked'] is None]
        if not candidates:
            break
        i, ev = min(candidates, key=lambda t: (order[t[1]['claimed_by']], t[1]['score'], t[1]['start']))
        a = i - 1 if ev['side'] == 'left' else i
        removed = segs[a][1]
        merged = (segs[a][0], segs[a + 1][1])
        segs = segs[:a] + [merged] + segs[a + 2:]
        steps.append({
            'step': len(steps), 'segments': [{'start': s, 'end': e} for s, e in segs],
            'removed_boundary': removed, 'merged_interval': [merged[0], merged[1]],
            'segment': [ev['start'], ev['end']], 'rule_id': ev['claimed_by'],
            'rule_name': names[ev['claimed_by']], 'side': ev['side'], 'score': ev['score'], 'phi': phi,
        })
    else:
        raise RuntimeError('Merge did not terminate within max_steps')
    return {'original': original, 'final': steps[-1]['segments'], 'steps': steps,
            'settings': {'rules': rules, 'allow_gap_crossing': allow_gap_crossing}}
