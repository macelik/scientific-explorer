"""Data store: cohort metadata, categorical CN matrix, bin coordinates,
production breakpoints/clusters, and on-demand `.X` rows read straight from
the h5ad's CSR arrays (no full-matrix load).

Coordinate systems used throughout the service:
  * genome-wide bin index  g   in [0, n_bins)              (h5ad var order)
  * chromosome-local index i   in [0, n_chrom)             (g = offset + i)
  * node-local index       j   in [0, n_node)              (i = node.start + j)
Genomic intervals from the h5ad are 1-based inclusive: [start, end].
"""
from __future__ import annotations

import hashlib
import json
import os
import threading
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import h5py
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from . import pyepi_adapter
from .config import Config


CN_LEVELS = [0.0, 0.5, 1.0, 1.5, 2.0]
CN_LABELS = {0.0: "loss", 0.5: "putative loss", 1.0: "base", 1.5: "putative gain", 2.0: "gain"}

# Seed exploration regions (Mb, chromosome). Starting intervals, not rules.
SEED_REGIONS = [
    {"id": "9p21", "label": "9p21 (CDKN2A) event", "chrom": "chr9", "start_mb": 19.0, "end_mb": 26.0,
     "note": "real focal deletion the pipeline under-calls"},
    {"id": "chr9-transition", "label": "chr9 accessibility transition", "chrom": "chr9", "start_mb": 27.0, "end_mb": 33.0,
     "note": "dominant technical artifact (~32.4 Mb edge)"},
    {"id": "chr9-shoulder", "label": "chr9 accessibility shoulder", "chrom": "chr9", "start_mb": 6.8, "end_mb": 7.4,
     "note": "secondary artifact, common depth-1 competitor"},
    {"id": "chr13-gain", "label": "chr13 gain", "chrom": "chr13", "start_mb": 94.0, "end_mb": 98.0, "note": ""},
    {"id": "chr1-distal-gain", "label": "chr1 distal gain", "chrom": "chr1", "start_mb": 195.0, "end_mb": 249.0, "note": ""},
]

REFERENCE_CELLS = [
    {"cell": "TTAGGCTAGGCCGGAA-1", "chrom": "chr9", "label": "typical chr9: depth-0 at 32.4 Mb artifact"},
    {"cell": "GAAGGATGTACGCGCA-1", "chrom": "chr9", "label": "chr9 shoulder case, small right-child segment"},
    {"cell": "TGTGGAGCAGGCTAGA-1", "chrom": "chr1", "label": "chr1 two-child relationship"},
]

DEFAULTS = {
    "n_permutations": 1000,
    "seed": 42,
    "alpha": 0.001,
    "rolling_window": 25,
    "rolling_min_periods": 8,
    "half_window_bins": 20,
    "n_competing_peaks": 5,
    "min_peak_separation_mb": 6.0,
    "log2_clip": 4.0,
    "production_k": 2,
}


def _file_identity(path: str) -> dict:
    st = os.stat(path)
    return {"path": os.path.abspath(path), "size": st.st_size, "mtime": int(st.st_mtime)}


@dataclass
class Chrom:
    name: str
    offset: int
    n: int


class DataStore:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        cfg.validate()
        self.sci = pyepi_adapter.load(cfg.package_dir)
        t0 = time.time()
        self._lock = threading.Lock()
        self._load_cohort()
        self._load_h5ad_meta()
        self._validate_alignment()
        self._load_breakpoints()
        self._cluster_order = None
        self._x_cache: Dict[str, np.ndarray] = {}
        self._cluster_cache: Dict[str, np.ndarray] = {}
        self.identity = self._compute_identity()
        self.load_seconds = time.time() - t0

    # ------------------------------------------------------------------ load
    def _load_cohort(self) -> None:
        cd = self.cfg.cohort_dir
        self.cells = pd.read_csv(os.path.join(cd, "cells.csv"))
        self.bins = pd.read_csv(os.path.join(cd, "bins.csv"))
        self.cn = np.load(os.path.join(cd, "cn_matrix.npy"))
        if self.cn.shape != (len(self.cells), len(self.bins)):
            raise ValueError(f"cn_matrix shape {self.cn.shape} != (cells {len(self.cells)}, bins {len(self.bins)})")
        vals = set(np.unique(self.cn).tolist())
        if not vals.issubset(set(CN_LEVELS)):
            raise ValueError(f"unexpected CN values in cn_matrix: {sorted(vals)}")
        self.cn_u8 = np.rint(self.cn * 2).astype(np.uint8)  # 0..4, categorical, never averaged
        self.cell_row = {c: k for k, c in enumerate(self.cells.cellID.values)}

    def _load_h5ad_meta(self) -> None:
        self.h5 = h5py.File(self.cfg.h5ad, "r")
        obs_ids = self.h5["obs/cellID"][...].astype(str)
        self.h5_row = {c: k for k, c in enumerate(obs_ids)}
        self.n_cells_total = len(obs_ids)
        cats = self.h5["var/seq/categories"][...].astype(str)
        codes = self.h5["var/seq/codes"][...]
        self.var_seq = cats[codes]
        self.var_start = self.h5["var/start"][...].astype(np.int64)
        self.var_end = self.h5["var/end"][...].astype(np.int64)
        self.n_bins = len(self.var_seq)
        self.indptr = self.h5["X/indptr"][...]
        self.X_data = self.h5["X/data"]
        self.X_indices = self.h5["X/indices"]
        # chromosome layout in h5ad order (first occurrence)
        order = list(dict.fromkeys(self.var_seq.tolist()))
        self.chroms: List[Chrom] = []
        for name in order:
            idx = np.where(self.var_seq == name)[0]
            if not np.array_equal(idx, np.arange(idx[0], idx[0] + len(idx))):
                raise ValueError(f"chromosome {name} bins are not contiguous in the h5ad")
            self.chroms.append(Chrom(name=name, offset=int(idx[0]), n=int(len(idx))))
        self.chrom_by_name = {c.name: c for c in self.chroms}

    def _validate_alignment(self) -> None:
        b = self.bins
        if len(b) != self.n_bins:
            raise ValueError("bins.csv row count differs from h5ad var")
        if not (np.array_equal(b.start.values, self.var_start) and np.array_equal(b.end.values, self.var_end)):
            raise ValueError("bins.csv coordinates do not match h5ad var")
        if not np.array_equal(b.seq.values.astype(str), self.var_seq):
            raise ValueError("bins.csv chromosomes do not match h5ad var")
        exp_local = np.concatenate([np.arange(c.n) for c in self.chroms])
        if not np.array_equal(b.bin_in_chrom.values, exp_local):
            raise ValueError("bins.csv bin_in_chrom is not consistent with chromosome layout")
        missing = [c for c in self.cells.cellID if c not in self.h5_row]
        if missing:
            raise ValueError(f"{len(missing)} cohort cells missing from h5ad, e.g. {missing[:3]}")
        widths = np.unique(self.var_end - self.var_start + 1)
        self.bin_width_bp = int(widths[0]) if len(widths) == 1 else None

    def _load_breakpoints(self) -> None:
        bp = pd.read_csv(os.path.join(self.cfg.cohort_dir, "breakpoints.csv"))
        self.breakpoints = bp
        self._bp_index: Dict[tuple, list] = {}
        for (cell, seq), grp in bp.groupby(["cell", "seq"]):
            self._bp_index[(cell, seq)] = sorted(
                [dict(bp=int(r.breakpoint), ad=float(r.ad_dist), p_local=float(r.p_value_local),
                      p_global=float(r.p_value_global)) for r in grp.itertuples()],
                key=lambda r: r["bp"])

    def _compute_identity(self) -> dict:
        parts = [
            _file_identity(self.cfg.h5ad),
            _file_identity(os.path.join(self.cfg.cohort_dir, "cn_matrix.npy")),
            _file_identity(os.path.join(self.cfg.cohort_dir, "cells.csv")),
            _file_identity(os.path.join(self.cfg.cohort_dir, "breakpoints.csv")),
        ]
        h = hashlib.sha1(json.dumps(parts, sort_keys=True).encode()).hexdigest()[:16]
        return {"hash": h, "files": parts, "n_cells": int(len(self.cells)), "n_bins": int(self.n_bins),
                "n_cells_total_h5ad": int(self.n_cells_total)}

    # ------------------------------------------------------------- accessors
    def production_breakpoints(self, cell: str, chrom: str) -> list:
        return list(self._bp_index.get((cell, chrom), []))

    def x_row(self, cell: str) -> np.ndarray:
        """Full-genome `.X` row (GC-corrected, normalised) as float64."""
        with self._lock:
            if cell in self._x_cache:
                return self._x_cache[cell]
        r = self.h5_row[cell]
        a, b = int(self.indptr[r]), int(self.indptr[r + 1])
        with self._lock:
            idx = self.X_indices[a:b]
            dat = self.X_data[a:b]
        row = np.zeros(self.n_bins, dtype=np.float64)
        row[idx] = dat
        with self._lock:
            self._x_cache[cell] = row
        return row

    def chrom_slice(self, chrom: str) -> slice:
        c = self.chrom_by_name[chrom]
        return slice(c.offset, c.offset + c.n)

    def production_clusters(self, cell: str) -> np.ndarray:
        """Segment id per genome-wide bin, rebuilt from production breakpoints
        with the exact loop pyEpiAneufinder uses (pyEpiAneufinder.py cluster
        construction). Verified equal to clusters.json in tests/."""
        if cell in self._cluster_cache:
            return self._cluster_cache[cell]
        counter = 1
        out: List[int] = []
        for c in self.chroms:
            bps = [r["bp"] for r in self.production_breakpoints(cell, c.name)]
            if not bps:
                out += [counter] * c.n
            else:
                edges = sorted([0, c.n] + bps)
                sizes = np.diff(edges)
                out += np.repeat(range(counter, len(sizes) + counter), sizes).tolist()
            counter = max(out) + 1
        arr = np.asarray(out, dtype=np.int32)
        self._cluster_cache[cell] = arr
        return arr

    def spliced_clusters(self, cell: str, chrom: str, boundaries: List[int]) -> np.ndarray:
        """Replace one chromosome's segmentation with the given chromosome-local
        boundaries, keeping every other chromosome's production segments and
        re-numbering so ids stay unique and monotonic (same convention)."""
        counter = 1
        out: List[int] = []
        for c in self.chroms:
            if c.name == chrom:
                bps = sorted(set(int(b) for b in boundaries if 0 < int(b) < c.n))
            else:
                bps = [r["bp"] for r in self.production_breakpoints(cell, c.name)]
            if not bps:
                out += [counter] * c.n
            else:
                edges = sorted([0, c.n] + bps)
                sizes = np.diff(edges)
                out += np.repeat(range(counter, len(sizes) + counter), sizes).tolist()
            counter = max(out) + 1
        return np.asarray(out, dtype=np.int32)

    def cluster_order(self) -> List[int]:
        """Ward/Euclidean hierarchical clustering on genome-wide categorical CN
        profiles; leaf order matches karyo_gainloss (dendrogram leaves, reversed)."""
        if self._cluster_order is None:
            cache = os.path.join(self.cfg.cache_dir, f"cluster_order_{self.identity['hash']}.json")
            if os.path.exists(cache):
                with open(cache) as f:
                    self._cluster_order = json.load(f)
            else:
                Z = linkage(pdist(self.cn, metric="euclidean"), method="ward")
                leaves = leaves_list(Z)[::-1]
                self._cluster_order = [int(v) for v in leaves]
                with open(cache, "w") as f:
                    json.dump(self._cluster_order, f)
        return self._cluster_order

    # ------------------------------------------------------------- metadata
    def meta(self) -> dict:
        cells = []
        for r in self.cells.itertuples():
            cells.append(dict(cellID=r.cellID, row=int(r.Index), raw_coverage=int(r.raw_coverage), tercile=r.tercile,
                              n_segments=int(r.n_segments), best_s=float(r.best_s), frac_loss=float(r.frac_loss),
                              frac_gain=float(r.frac_gain), frac_base=float(r.frac_base)))
        chroms = []
        for c in self.chroms:
            s = self.var_start[c.offset]
            e = self.var_end[c.offset + c.n - 1]
            chroms.append(dict(name=c.name, offset=c.offset, n=c.n, start_bp=int(s), end_bp=int(e)))
        return {
            "dataset": {
                "identity": self.identity,
                "h5ad": self.cfg.h5ad,
                "cohort_dir": self.cfg.cohort_dir,
                "package_dir": self.cfg.package_dir,
                "bin_width_bp": self.bin_width_bp,
                "signal": ".X (GC-corrected, integer-rounded, total-normalised to 1e5 per cell)",
            },
            "chromosomes": chroms,
            "cells": cells,
            "cluster_order": self.cluster_order(),
            "cn_levels": [{"value": v, "label": CN_LABELS[v]} for v in CN_LEVELS],
            "regions": SEED_REGIONS,
            "reference_cells": REFERENCE_CELLS,
            "defaults": DEFAULTS,
        }

    def bins_binary(self) -> bytes:
        arr = np.empty((self.n_bins, 2), dtype=np.int32)
        arr[:, 0] = self.var_start
        arr[:, 1] = self.var_end
        return arr.tobytes()

    def karyogram_binary(self) -> bytes:
        return np.ascontiguousarray(self.cn_u8).tobytes()
