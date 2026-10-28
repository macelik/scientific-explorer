import numpy as np

from server import allelic_scan as sc

FAST = {'k': 1, 'n_permutations': 200, 'alpha': 0.01, 'min_child_bins': 20, 'min_terminal_bins': 20}


def _step(n_left=80, n_right=80, lo=0.05, hi=0.30, sd=0.01, seed=1):
    rng = np.random.default_rng(seed)
    return np.r_[rng.normal(lo, sd, n_left), rng.normal(hi, sd, n_right)]


def test_geometry_rules_match_the_depth_run():
    ok, why = sc.validate_candidate_geometry(10, 0, 100, 0, 100, 20, 50)
    assert not ok and why == ('left_child_too_short', 'left_terminal_too_short')
    ok, why = sc.validate_candidate_geometry(30, 20, 100, 0, 100, 20, 50)
    assert ok and why == ()          # internal left edge: no terminal rule
    ok, why = sc.validate_candidate_geometry(60, 0, 100, 0, 100, 20, 50)
    assert not ok and why == ('right_terminal_too_short',)


def test_step_gives_one_breakpoint_and_two_segments():
    x = _step()
    out = sc.scan_chromosome(x, np.ones(len(x), bool), FAST)
    assert [b['bin'] for b in out['breakpoints']] == [80]
    b = out['breakpoints'][0]
    assert b['gap_start'] is None and b['gap_end'] is None and b['recursion_depth'] == 0
    assert b['p_local'] < 0.01 and b['p_global'] < 0.01
    segs = out['segments']
    assert [(s['bin_start'], s['bin_end'], s['kept_bins']) for s in segs] == [(0, 79, 80), (80, 159, 80)]
    assert abs(segs[0]['mean'] - 0.05) < 0.01 and abs(segs[1]['mean'] - 0.30) < 0.01


def test_skipped_bins_map_back_with_a_gap():
    x = _step()
    kept = np.ones(len(x), bool)
    kept[75:85] = False                       # the change at 80 falls inside skipped bins
    out = sc.scan_chromosome(x, kept, FAST)
    assert len(out['breakpoints']) == 1
    b = out['breakpoints'][0]
    assert (b['bin'], b['gap_start'], b['gap_end']) == (85, 75, 85)
    segs = out['segments']
    assert [(s['bin_start'], s['bin_end'], s['n_bins'], s['kept_bins']) for s in segs] == [(0, 84, 85, 75), (85, 159, 75, 75)]


def test_terminal_rule_blocks_a_step_near_the_end():
    x = _step(n_left=10, n_right=150)
    out = sc.scan_chromosome(x, np.ones(len(x), bool), {**FAST, 'min_child_bins': 5, 'min_terminal_bins': 20})
    assert 10 not in [b['bin'] for b in out['breakpoints']]
    top = [r for r in out['audit'] if r['candidate_rank'] == 1][0]
    assert top['candidate_absolute_bin'] == 10 and 'left_terminal_too_short' in top['rejection_reasons']


def test_noise_stops_at_the_first_tested_candidate():
    x = np.random.default_rng(3).normal(0.1, 0.01, 160)
    out = sc.scan_chromosome(x, np.ones(len(x), bool), {**FAST, 'k': 3})
    assert out['breakpoints'] == []
    tested = [r for r in out['audit'] if r['selected_for_significance']]
    assert len(tested) == 1 and tested[0]['recursion_depth'] == 0 and tested[0]['p_local'] >= 0.01


def test_no_or_one_kept_bin():
    x = _step()
    none = sc.scan_chromosome(x, np.zeros(len(x), bool), FAST)
    assert none['breakpoints'] == [] and none['audit'] == []
    assert none['segments'] == [{'bin_start': 0, 'bin_end': 159, 'n_bins': 160, 'kept_bins': 0, 'mean': None, 'median': None}]
    one = np.zeros(len(x), bool); one[5] = True
    s = sc.scan_chromosome(x, one, FAST)['segments']
    assert len(s) == 1 and s[0]['kept_bins'] == 1 and abs(s[0]['mean'] - x[5]) < 1e-12


def test_sums_cover_kept_bins_only():
    x = _step()
    kept = np.ones(len(x), bool); kept[:5] = False
    ha = np.ones(len(x)); hb = np.full(len(x), 2.0)
    s = sc.scan_chromosome(x, kept, FAST, ha=ha, hb=hb)['segments'][0]
    assert (s['ha_sum'], s['hb_sum'], s['at_sum']) == (75.0, 150.0, 225.0)


def test_segments_partition_the_chromosome_and_run_is_deterministic():
    x = np.r_[_step(seed=4), _step(seed=5) + 0.1]
    settings = {**FAST, 'k': 3}
    a = sc.scan_chromosome(x, np.ones(len(x), bool), settings)
    b = sc.scan_chromosome(x, np.ones(len(x), bool), settings)
    assert a == b
    segs = a['segments']
    assert segs[0]['bin_start'] == 0 and segs[-1]['bin_end'] == len(x) - 1
    assert all(segs[i]['bin_end'] + 1 == segs[i + 1]['bin_start'] for i in range(len(segs) - 1))


def test_scan_task_round_trip():
    x = _step()
    src, chrom, out = sc.scan_task(('cluster1', 'chr2', x, np.ones(len(x), bool), None, None, FAST, None, None))
    assert (src, chrom) == ('cluster1', 'chr2') and [b['bin'] for b in out['breakpoints']] == [80]


def test_global_test_uses_the_given_pool(monkeypatch):
    sc._ensure()
    seen = []
    real = sc.sci.global_permutation_test_ad

    def spy(seq_data, *a, **k):
        seen.append(len(seq_data))
        return real(seq_data, *a, **k)
    monkeypatch.setattr(sc.sci, 'global_permutation_test_ad', spy)
    x = _step()
    kept = np.ones(len(x), bool); kept[:4] = False
    sc.scan_chromosome(x, kept, FAST)
    pool = np.random.default_rng(9).normal(0.1, 0.05, 500)
    sc.scan_chromosome(x, kept, FAST, pool=pool)
    assert seen == [156, 500]                 # chromosome kept bins, then the given genome pool
