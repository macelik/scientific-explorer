# server/merge_rounds.py
"""Merge rounds: the working segmentation is the original accepted breakpoints
with each round's genome-wide rule merge applied in turn. Every round merges on
the original pb_x; only the boundaries carry over.
Spec: docs/superpowers/specs/2026-10-27-merge-rounds-design.md §2."""
from __future__ import annotations

import hashlib
import json
import threading
from collections import OrderedDict
from dataclasses import dataclass, field, replace
from typing import List

import numpy as np
import pandas as pd

from . import cn_recall as cr
from . import merge_rules as mr


def normalize_rounds(rounds) -> List[dict]:
    if rounds is None:
        return []
    if not isinstance(rounds, list):
        raise ValueError('rounds must be a list')
    out = []
    for i, r in enumerate(rounds, 1):
        if not isinstance(r, dict) or not isinstance(r.get('rules'), list):
            raise ValueError(f'round {i}: rules must be a list')
        try:
            mr.validate_rules(r['rules'])
        except ValueError as e:
            raise ValueError(f'round {i}: {e}') from e
        out.append({'rules': r['rules'], 'allow_gap_crossing': bool(r.get('allow_gap_crossing', False))})
    return out


def rounds_hash(rounds) -> str:
    return hashlib.sha256(json.dumps(rounds, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


@dataclass
class BaseInputs:
    """What merge_genome needs; RecallInputs also satisfies this shape."""
    var_seq: np.ndarray
    var_start: np.ndarray
    var_end: np.ndarray
    pb_x: np.ndarray
    breakpoints: pd.DataFrame
    sources: List[str] = field(default_factory=list)

    def chromosomes(self) -> List[str]:
        return list(dict.fromkeys(self.var_seq.tolist()))


def fold(inp, rounds):
    """Apply each (normalized) round in turn; returns (breakpoints, per-round summaries)."""
    bp, per = inp.breakpoints, []
    for i, r in enumerate(rounds, 1):
        try:
            m = cr.merge_genome(replace(inp, breakpoints=bp), r['rules'], r['allow_gap_crossing'])
        except cr.RecallError as e:
            raise cr.RecallError(f'round {i}: {e}') from e
        bp = m['merged_breakpoints'].reset_index(drop=True)
        per.append({'pairs': m['pairs'], 'merges': sum(p['merges'] for p in m['pairs'])})
    return bp, per


class BaseCache:
    def __init__(self, size: int = 8):
        self.size = size
        self._d: OrderedDict = OrderedDict()
        self.lock = threading.Lock()

    def get(self, inp, rounds) -> dict:
        h = rounds_hash(rounds)
        key = (id(inp.breakpoints), h)
        with self.lock:
            e = self._d.get(key)
            if e is not None and e['source'] is inp.breakpoints:
                self._d.move_to_end(key)
                return e
        bp, per = fold(inp, rounds)
        e = {'source': inp.breakpoints, 'hash': h, 'breakpoints': bp, 'rounds': per, 'extra': {}}
        with self.lock:
            self._d[key] = e
            self._d.move_to_end(key)
            while len(self._d) > self.size:
                self._d.popitem(last=False)
        return e


CACHE = BaseCache(8)
