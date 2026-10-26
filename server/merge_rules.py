"""User-defined merge rules for the exploratory adjacent-segment merge.

A rule is a list of conditions on one internal segment's left/right flank
scores that must ALL hold; a segment is flagged when ANY enabled rule matches
(OR of ANDs). The topmost matching rule claims the segment, and that rule's
mergeInto metric picks the neighbour to merge into (smaller |score|, ties
left). Rows are the per-segment dicts produced by hatch_merge.score_rows, the
same rows /api/integration/hatch-scores returns, so web/src/mergeRules.ts can
evaluate identical rules client-side (see tests/fixtures/merge_rules_parity.json).

rules[0] is Rule 0, the always-active rule built from the tracks' section-1
settings; it is evaluated first and claims before any user rule.

These are exploratory diagnostics, independent of the original
bootstrap/chain-guard merge prototype; nothing here recomputes CN.
"""
from __future__ import annotations

import math

from . import hatch_merge as hm

CONDITION_KINDS = ('weak_flank', 'sign_mismatch', 'ambiguous', 'max_size')
SRD_METRICS = ('srd', 'srd_phi')
HATCH_PATTERNS = ('/', '\\', 'x', '-', '|', '.')
RULE0_ID = 'rule0'
RULE0_TYPE = 'flank_threshold'
RULE0_CRITERIA = ('delta_median', 'delta_mean', 'delta_iqr_mean', 'srd_phi')
RULE0_COMBINE = ('all', 'any')


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


def _positive(v):
    return _number(v) and math.isfinite(v) and v > 0


def _check_hatch(r, where):
    hatch = r.get('hatch') if isinstance(r.get('hatch'), dict) else {}
    if hatch.get('pattern') not in HATCH_PATTERNS:
        raise ValueError(f'{where}: unknown hatch pattern {hatch.get("pattern")!r}')


def _validate_rule0(r):
    where = 'Rule 0'
    criteria = r.get('criteria')
    if not isinstance(criteria, list) or not criteria:
        raise ValueError(f'{where}: check at least one criterion')
    if not all(isinstance(c, str) for c in criteria) or len(set(criteria)) != len(criteria) or any(c not in RULE0_CRITERIA for c in criteria):
        raise ValueError(f'{where}: criteria must be distinct values from {RULE0_CRITERIA}')
    if r.get('combine') not in RULE0_COMBINE:
        raise ValueError(f'{where}: combine must be one of {RULE0_COMBINE}')
    for key in ('tau', 'srdRef'):
        if not _positive(r.get(key)):
            raise ValueError(f'{where}: {key} must be a positive finite number')
    n = r.get('maxBins')
    if not (isinstance(n, int) and not isinstance(n, bool) and n >= 1):
        raise ValueError(f'{where}: maxBins must be an integer >= 1')
    for key in ('mergeSignMismatch', 'mergeFlankConflict'):
        if not isinstance(r.get(key), bool):
            raise ValueError(f'{where}: {key} must be true or false')
    _check_hatch(r, where)


def validate_rules(rules):
    """Validate and return [Rule 0, *enabled user rules], in list order.

    rules[0] must be Rule 0 (always active). Disabled user rules are ignored
    apart from the id checks, so a half-edited disabled rule never blocks a
    run. Raises ValueError on any problem."""
    if not isinstance(rules, list) or not all(isinstance(r, dict) for r in rules):
        raise ValueError('rules must be a list of objects')
    ids = [r.get('id') for r in rules]
    if any(not isinstance(i, str) or not i for i in ids):
        raise ValueError('every rule needs a non-empty string id')
    if len(set(ids)) != len(ids):
        raise ValueError('rule ids must be unique')
    if not rules or rules[0].get('id') != RULE0_ID or rules[0].get('type') != RULE0_TYPE:
        raise ValueError('rules[0] must be Rule 0 (id "rule0", type "flank_threshold")')
    for r in rules[1:]:
        if r.get('type', 'conditions') != 'conditions':
            raise ValueError(f'rule {r["id"]!r}: only rules[0] may be Rule 0')
    _validate_rule0(rules[0])
    enabled = [r for r in rules[1:] if r.get('enabled') is True]
    for r in enabled:
        where = f'rule {r["id"]!r}'
        conditions = r.get('conditions')
        if not isinstance(conditions, list) or not conditions:
            raise ValueError(f'{where}: needs at least one condition')
        for c in conditions:
            _validate_condition(c, where)
        merge_into = r.get('mergeInto') if isinstance(r.get('mergeInto'), dict) else {}
        _check_metric(merge_into.get('metric'), merge_into.get('estimator', 'mean'), f'{where} mergeInto')
        _check_hatch(r, where)
    return [rules[0], *enabled]


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


def _criterion_value(row, side, criterion):
    v = row.get(f'{criterion}_{side}')
    return float(v) if v is not None and math.isfinite(v) else None


def _smaller_side(left, right):
    if left is None and right is None:
        return None
    if right is None:
        return 'left'
    if left is None:
        return 'right'
    return 'left' if abs(left) <= abs(right) else 'right'


def evaluate_rule0(row, r0, allow_gap_crossing=False):
    """Rule 0 on one score row: per-criterion flank/threshold checks, the
    all/any test, sign and flank-conflict checks, direction, ordering score
    and outcome. See docs/superpowers/specs/2026-10-26-rule0-visual-merge-design.md §3."""
    out = {'start': row['start'], 'end': row['end'], 'evaluated': False, 'criteria': [], 'passed': False,
           'sign_mismatch': False, 'flank_conflict': False, 'held_back': [], 'outcome': 'none',
           'side': None, 'score': None, 'blocked': None}
    if row['end'] - row['start'] > r0['maxBins']:
        return out
    out['evaluated'] = True
    for c in r0['criteria']:
        left, right = _criterion_value(row, 'left', c), _criterion_value(row, 'right', c)
        threshold = r0['srdRef'] if c == 'srd_phi' else r0['tau']
        side = _smaller_side(left, right)
        value = None if side is None else abs(left if side == 'left' else right)
        out['criteria'].append({'criterion': c, 'side': side, 'value': value, 'threshold': threshold,
                                'fraction': None if value is None else value / threshold,
                                'passed': value is not None and value < threshold})
    passing = [c for c in out['criteria'] if c['passed']]
    out['passed'] = len(passing) == len(out['criteria']) if r0['combine'] == 'all' else bool(passing)
    sl, sr = _criterion_value(row, 'left', 'srd_phi'), _criterion_value(row, 'right', 'srd_phi')
    out['sign_mismatch'] = sl is not None and sr is not None and sl * sr < 0
    if not out['passed']:
        if out['sign_mismatch']:
            out['outcome'] = 'signs_only'
        return out
    out['flank_conflict'] = len({c['side'] for c in passing}) > 1
    fractions = [c['fraction'] for c in passing]
    out['score'] = max(fractions) if r0['combine'] == 'all' else min(fractions)
    best = min(passing, key=lambda c: (c['fraction'], 0 if c['side'] == 'left' else 1))
    out['side'] = best['side'] if out['flank_conflict'] else passing[0]['side']
    if out['sign_mismatch'] and not r0['mergeSignMismatch']:
        out['held_back'].append('signs')
    if out['flank_conflict'] and not r0['mergeFlankConflict']:
        out['held_back'].append('flanks')
    if out['held_back']:
        out['outcome'] = 'held_back'
        return out
    out['outcome'] = 'would_merge'
    if not allow_gap_crossing and row[f'gap_{out["side"]}']:
        out['blocked'] = 'gap'
    return out


def evaluate_rules(rows, rules, allow_gap_crossing=False):
    """Per score row: which rules match, which claims it, and the merge
    direction/score or why it cannot merge. `rules` must be validated; a
    leading Rule 0 is evaluated first and claims before user rules."""
    r0 = rules[0] if rules and rules[0].get('type') == RULE0_TYPE else None
    user = rules[1:] if r0 is not None else rules
    out = []
    for row in rows:
        r0ev = evaluate_rule0(row, r0, allow_gap_crossing) if r0 is not None else None
        matched = ([RULE0_ID] if r0ev is not None and r0ev['outcome'] == 'would_merge' else []) + \
            [r['id'] for r in user if all(condition_holds(c, row) for c in r['conditions'])]
        ev = {'start': row['start'], 'end': row['end'], 'matched': matched,
              'claimed_by': None, 'side': None, 'score': None, 'blocked': None, 'rule0': r0ev}
        if matched and matched[0] == RULE0_ID:
            ev.update(claimed_by=RULE0_ID, side=r0ev['side'], score=r0ev['score'], blocked=r0ev['blocked'])
        elif matched:
            rule = next(r for r in user if r['id'] == matched[0])
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
