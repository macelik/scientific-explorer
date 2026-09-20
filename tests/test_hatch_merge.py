from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from server import hatch_merge as hm


def test_flat_partition_half_open_and_drops_empty():
    assert hm.flat_partition([3, 3, 7], 10) == [(0, 3), (3, 7), (7, 10)]
    assert hm.flat_partition([], 5) == [(0, 5)]


def test_segment_summary_adds_start_end_sum_to_summarize():
    profile = np.array([2., 2., 4., 4.])
    s = hm.segment_summary(profile, 1, 3)
    assert s['start'] == 1 and s['end'] == 3
    assert s['sum'] == pytest.approx(6.0)
    assert s['mean'] == pytest.approx(3.0)
