"""Compare the archived small-segment flank tables with live hatch-scores on
the pilot chromosomes (spec 2026-10-26 §2). Needs a running server:
    BASE_URL=http://127.0.0.1:8767 python3 tests/compare_archived_live.py
Writes validation/rule0/archived_vs_live.json and prints a summary."""
import json
import math
import os
import sys
import urllib.request
from pathlib import Path

BASE = os.environ.get('BASE_URL', 'http://127.0.0.1:8766')
OUT = Path(__file__).resolve().parents[1] / 'validation' / 'rule0' / 'archived_vs_live.json'

# (archived table, archived column stem, live column stem)
PAIRS = [
    ('flank_scores', 'delta', 'delta_median'),
    ('flank_scores', 'delta_mean', 'delta_mean'),
    ('flank_scores', 'delta_trim', 'delta_iqr_mean'),
    ('flank_scores', 'srd_phi', 'srd_phi'),
    ('screening', 'log2fc', 'delta_mean'),
    ('screening', 'srd', 'srd'),
]


def get(path):
    with urllib.request.urlopen(BASE + path, timeout=300) as r:
        return json.loads(r.read())


def post(path, body):
    req = urllib.request.Request(BASE + path, data=json.dumps(body).encode(), headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.loads(r.read())


def num(v):
    return None if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)


def main():
    meta = get('/api/integration/meta')
    pilots, sources = meta['pilot_chromosomes'], meta['sources']
    live = {}
    for chrom in pilots:
        for src in sources:
            res = post('/api/integration/hatch-scores', {'chrom': chrom, 'source': src})
            for row in res['rows']:
                live[(src, chrom, row['start'], row['end'])] = row
    report = {'base_url': BASE, 'pilot_chromosomes': pilots, 'pairs': []}
    for table, a_stem, l_stem in PAIRS:
        rows = meta.get(table) or []
        matched = unmatched = compared = 0
        max_diff, worst = 0.0, None
        for r in rows:
            key = (r['source'], r['chromosome'], int(r['bin_start']), int(r['bin_end']))
            lr = live.get(key)
            if lr is None:
                unmatched += 1
                continue
            matched += 1
            for side in ('left', 'right'):
                a, b = num(r.get(f'{a_stem}_{side}')), num(lr.get(f'{l_stem}_{side}'))
                if a is None or b is None:
                    continue
                compared += 1
                d = abs(a - b)
                if d > max_diff:
                    max_diff, worst = d, {'key': list(key), 'side': side, 'archived': a, 'live': b}
        report['pairs'].append({'archived': f'{table}.{a_stem}', 'live': l_stem, 'archived_rows': len(rows),
                                'matched_segments': matched, 'unmatched_segments': unmatched,
                                'compared_values': compared, 'max_abs_diff': max_diff, 'worst': worst})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + '\n')
    for p in report['pairs']:
        print(f"{p['archived']:>26} vs {p['live']:<15} matched {p['matched_segments']:>4}/{p['archived_rows']:<4} "
              f"values {p['compared_values']:>4} max|diff| {p['max_abs_diff']:.3g}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
