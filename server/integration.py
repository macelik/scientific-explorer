"""Multimodal integration story: read-only adapters over the exported artifacts.

Inputs (all optional; the tabs report what is missing):
  * integrated.h5mu       - cells, cluster labels (several resolutions), WNN weights,
                            coverage fields, three UMAPs, PCA, joint/modality graphs,
                            allelic abs-dBAF matrix (mod/haplo/X)
  * count_matrix.h5ad     - the explorer's own X (GC-corrected, normalised, NOT log1p) used
                            for the per-cluster pseudobulk depth profiles
  * raw_accepted_breakpoints.tsv (per cluster + all)  - pseudobulk segmentation
  * data/derived/*        - cluster-level segments/CN and per-cell CN under cluster segmentation
  * data/prototype/*      - small-segment flank metrics, SRD screen, merge prototype outputs

Nothing here recomputes clustering, integration or the modified caller. Heavy arrays are
built once in a background thread and cached under the app cache directory keyed by the
input fingerprints.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import json
import os
import threading
import time
from typing import Dict, List, Optional

import h5py
import numpy as np
import pandas as pd

from .config import Config


SOURCES = ["cluster0", "cluster1", "cluster2", "cluster3", "cluster4", "all"]
CLUSTER_KEYS = {
    "wnn": ["wnn_leiden_0.3", "wnn_leiden_0.5", "wnn_leiden_1.0", "wnn_leiden_2.0"],
    "depth": ["depth_leiden_0.3", "depth_leiden_0.5", "depth_leiden_1.0", "depth_leiden_2.0"],
    "haplo": ["haplo_leiden_0.3", "haplo_leiden_0.5", "haplo_leiden_1.0", "haplo_leiden_2.0"],
}
CHR8_SCORE = {"chrom": "chr8", "num_mb": (46.0, 120.0), "den_mb": (10.0, 35.0)}


def _fid(path: str) -> dict:
    st = os.stat(path)
    return {"path": os.path.abspath(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def _read_cat(g: h5py.Group, key: str) -> np.ndarray:
    obj = g[key]
    if isinstance(obj, h5py.Group):
        cats = obj["categories"][...].astype(str)
        codes = obj["codes"][...]
        out = cats[codes]
        out[codes < 0] = ""
        return out
    v = obj[...]
    return v.astype(str) if v.dtype.kind in "OS" else v


def _csr_rows(g: h5py.Group, lo: int, hi: int, n_cols: int) -> np.ndarray:
    indptr = g["indptr"][lo:hi + 1]
    a, b = int(indptr[0]), int(indptr[-1])
    data = g["data"][a:b]
    idx = g["indices"][a:b]
    out = np.zeros((hi - lo, n_cols), dtype=np.float64)
    for r in range(hi - lo):
        s, e = int(indptr[r] - a), int(indptr[r + 1] - a)
        out[r, idx[s:e]] = data[s:e]
    return out


class IntegrationStore:
    def __init__(self, cfg: Config, var_seq: np.ndarray, var_start: np.ndarray, var_end: np.ndarray,
                 h5ad_path: str, cohort_cells: List[str]):
        self.cfg = cfg
        self.var_seq, self.var_start, self.var_end = var_seq, var_start, var_end
        self.n_bins = len(var_seq)
        self.h5ad_path = h5ad_path
        self.cohort_cells = set(cohort_cells)
        self.problems: List[str] = []
        self.notes: List[str] = []
        self.manifest: dict = {}
        self.h5mu_path = self._resolve_h5mu()
        self.available = self.h5mu_path is not None
        self.index_state = {"status": "idle", "step": "", "done": 0, "total": 6, "started": None, "finished": None, "error": None}
        self._lock = threading.RLock()
        self.arrays: Optional[dict] = None
        self._karyo_gz: Optional[bytes] = None
        self.tables: Dict[str, pd.DataFrame] = {}
        self.cells: Optional[pd.DataFrame] = None
        self.emb: Dict[str, np.ndarray] = {}
        self.graphs: Dict[str, dict] = {}
        if self.available:
            try:
                self._load_light()
            except Exception as e:  # noqa: BLE001
                self.available = False
                self.problems.append(f"failed to load integration inputs: {type(e).__name__}: {e}")
        self.identity = self._identity()

    # ------------------------------------------------------------ resolve
    def _resolve_h5mu(self) -> Optional[str]:
        cands = []
        if self.cfg.h5mu:
            cands.append(self.cfg.h5mu)
        cands.append(os.path.join(self.cfg.integration_dir, "data", "integrated.h5mu"))
        man = os.path.join(self.cfg.prototype_dir, "run_manifest.json")
        if os.path.exists(man):
            try:
                with open(man) as f:
                    self.manifest = json.load(f)
                p = self.manifest.get("inputs", {}).get("integrated")
                if p:
                    cands.append(p)
            except Exception as e:  # noqa: BLE001
                self.problems.append(f"run_manifest.json unreadable: {e}")
        for c in cands:
            if c and os.path.isfile(os.path.realpath(c)):
                if c != cands[0]:
                    self.notes.append(f"h5mu resolved via fallback: {os.path.realpath(c)}")
                return os.path.realpath(c)
        self.problems.append("integrated.h5mu not found (checked: " + "; ".join(x for x in cands if x) + "). Set PYEPI_SCIENTIFIC_H5MU.")
        return None

    def _identity(self) -> dict:
        files = []
        paths = [self.h5mu_path, self.h5ad_path, self.cfg.cluster_breakpoints]
        for folder in [self.cfg.derived_dir, self.cfg.prototype_dir]:
            paths.extend(str(p) for p in sorted(Path(folder).glob('*')) if p.suffix in ('.tsv', '.json', '.npy'))
        for p in dict.fromkeys(paths):
            if p and os.path.isfile(p):
                files.append(_fid(p))
        payload = {"files": files, "adapter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        h = hashlib.sha1(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]
        return {"hash": h, "files": files}

    # ------------------------------------------------------------- light
    def _load_light(self) -> None:
        f = h5py.File(self.h5mu_path, "r")
        self.h5 = f
        ids = f["obs/_index"][...].astype(str)
        self.cell_ids = ids
        self.cell_index = {c: i for i, c in enumerate(ids)}
        n = len(ids)
        if len(set(ids)) != n:
            raise ValueError("Duplicate MuData cell IDs")
        with h5py.File(self.h5ad_path, "r") as depth:
            if not np.array_equal(depth["obs/cellID"][...].astype(str), ids):
                raise ValueError("count_matrix.h5ad cell order differs from h5mu obs")
        # bin alignment against the explorer's h5ad
        dv = f["mod/depth/var"]
        dseq = _read_cat(dv, "seq")
        if len(dseq) != self.n_bins or not np.array_equal(dseq, self.var_seq) or not np.array_equal(dv["start"][...], self.var_start) or not np.array_equal(dv["end"][...], self.var_end):
            raise ValueError("h5mu depth bins do not match the explorer's count_matrix bins")
        nh = f["mod/haplo/var/_index"].shape[0]
        if nh != self.n_bins:
            raise ValueError(f"haplo modality has {nh} bins, depth has {self.n_bins}")
        self.notes.append("Allelic bins are assumed positionally aligned to depth. The supplied haplo var has no genomic coordinates, so equal bin counts do not independently verify this alignment.")
        obs = {"cell": ids}
        for k in ["gc_norm_depth", "snp_coverage", "total_raw_depth", "n_bins_snp_covered", "wnn_depth_weight", "wnn_haplo_weight"]:
            obs[k] = f["obs"][k][...]
        for k in CLUSTER_KEYS["wnn"]:
            obs[k] = _read_cat(f["obs"], k)
        for m in ["depth", "haplo"]:
            mo = f[f"mod/{m}/obs"]
            mids = _read_cat(mo, "_index") if "_index" in mo else _read_cat(mo, "cellID")
            if not np.array_equal(mids, ids):
                raise ValueError(f"mod/{m} obs order differs from mudata obs")
            for k in CLUSTER_KEYS[m]:
                if k in mo:
                    obs[k] = _read_cat(mo, k)
        self.cells = pd.DataFrame(obs)
        self.cells["in_explorer_cohort"] = self.cells["cell"].isin(self.cohort_cells)
        # embeddings
        self.emb = {
            "wnn": f["obsm/X_umap"][...].astype(np.float32),
            "depth": f["mod/depth/obsm/X_umap"][...].astype(np.float32),
            "haplo": f["mod/haplo/obsm/X_umap"][...].astype(np.float32),
            "depth_pc1": f["mod/depth/obsm/X_pca"][:, 0].astype(np.float32),
            "haplo_pc1": f["mod/haplo/obsm/X_pca"][:, 0].astype(np.float32),
            "depth_pc2": f["mod/depth/obsm/X_pca"][:, 1].astype(np.float32),
            "haplo_pc2": f["mod/haplo/obsm/X_pca"][:, 1].astype(np.float32),
        }
        # graphs (sparse, kept as CSR arrays)
        for name, grp in [("wnn", "obsp/connectivities"), ("depth", "mod/depth/obsp/connectivities"), ("haplo", "mod/haplo/obsp/connectivities")]:
            g = f[grp]
            self.graphs[name] = {"indptr": g["indptr"][...], "indices": g["indices"][...], "data": g["data"][...]}
        self.uns = {}
        try:
            self.uns["neighbors"] = {k: (f[f"uns/neighbors/params/{k}"][()].decode() if isinstance(f[f"uns/neighbors/params/{k}"][()], bytes) else (f[f"uns/neighbors/params/{k}"][()].item() if isinstance(f[f"uns/neighbors/params/{k}"], h5py.Dataset) else "group")) for k in f["uns/neighbors/params"]}
        except Exception:  # noqa: BLE001
            pass
        # tables
        bp = pd.read_csv(self.cfg.cluster_breakpoints, sep="\t")
        self.tables["breakpoints"] = bp
        for name in ["karyogram_cluster_segments", "percell_best_s"]:
            p = os.path.join(self.cfg.derived_dir, name + ".tsv")
            if os.path.exists(p):
                self.tables[name] = pd.read_csv(p, sep="\t")
            else:
                self.problems.append(f"derived table missing: {p}")
        for name in ["small_segment_flank_scores", "small_segment_scores", "small_segment_bootstrap", "srd_descriptive", "retained_boundaries",
                     "simplified_segments", "initial_boundary_scores", "unresolved_segments", "threshold_sensitivity", "review_summary"]:
            p = os.path.join(self.cfg.prototype_dir, name + ".tsv")
            if os.path.exists(p):
                self.tables[name] = pd.read_csv(p, sep="\t")
            else:
                self.notes.append(f"prototype table not available: {name}")
        audit = os.path.join(self.cfg.prototype_dir, "segmentation_candidate_audit.tsv")
        self._audit_path = audit if os.path.exists(audit) else None
        self._audit: Optional[pd.DataFrame] = None
        if "percell_best_s" in self.tables:
            t = self.tables["percell_best_s"].set_index("cell")
            self.cells["percell_best_s"] = self.cells["cell"].map(t["best_s"]).astype(float)
        # cluster membership per source, validated against manifest sizes
        self.members = {}
        lab = self.cells["wnn_leiden_0.3"].to_numpy()
        for s in SOURCES:
            self.members[s] = np.arange(n) if s == "all" else np.where(lab == s[-1])[0]
        ms = self.manifest.get("cluster_sizes", {})
        for s, v in ms.items():
            if s in self.members and len(self.members[s]) != v:
                self.problems.append(f"{s}: manifest size {v} != h5mu membership {len(self.members[s])}")
        nb = bp.groupby("source")["n_cells"].first().to_dict()
        for s, v in nb.items():
            if s in self.members and len(self.members[s]) != v:
                self.problems.append(f"{s}: breakpoints n_cells {v} != h5mu membership {len(self.members[s])}")

    # ------------------------------------------------------------ heavy
    def cache_path(self) -> str:
        return os.path.join(self.cfg.cache_dir, f"integration_{self.identity['hash']}.npz")

    def start_index(self) -> None:
        with self._lock:
            if not self.available or self.index_state["status"] in ("running", "done"):
                return
            self.index_state.update({"status": "running", "started": time.time(), "error": None, "done": 0})
        threading.Thread(target=self._run_index, daemon=True).start()

    def _step(self, label: str, done: int) -> None:
        with self._lock:
            self.index_state.update({"step": label, "done": done})

    def _run_index(self) -> None:
        try:
            cp = self.cache_path()
            if os.path.exists(cp):
                self._step("loading cache", 5)
                z = np.load(cp)
                self.arrays = {k: z[k] for k in z.files}
            else:
                self.arrays = self._build_arrays()
                np.savez(cp, **self.arrays)
            self._step("per-cell karyogram", 6)
            self._prepare_karyo()
            with self._lock:
                self.index_state.update({"status": "done", "finished": time.time(), "done": 6, "step": "done"})
        except Exception as e:  # noqa: BLE001
            with self._lock:
                self.index_state.update({"status": "error", "error": f"{type(e).__name__}: {e}"})

    def _build_arrays(self) -> dict:
        n = len(self.cell_ids)
        out: dict = {}
        # ---- pass over X (explorer h5ad; NOT the log1p copy inside the h5mu)
        self._step("pseudobulk depth: pass over X rows", 1)
        g = h5py.File(self.h5ad_path, "r")
        hid = g["obs/cellID"][...].astype(str)
        if not np.array_equal(hid, self.cell_ids):
            raise ValueError("count_matrix.h5ad cell order differs from h5mu obs")
        Xg = g["X"]
        sums = np.zeros((len(SOURCES), self.n_bins), dtype=np.float64)
        src_of = np.full(n, -1, dtype=int)
        for i, s in enumerate(SOURCES[:-1]):
            src_of[self.members[s]] = i
        c8 = np.where(self.var_seq == CHR8_SCORE["chrom"])[0]
        mb = self.var_start[c8] / 1e6
        num = c8[(mb >= CHR8_SCORE["num_mb"][0]) & (mb <= CHR8_SCORE["num_mb"][1])]
        den = c8[(mb >= CHR8_SCORE["den_mb"][0]) & (mb <= CHR8_SCORE["den_mb"][1])]
        chr8 = np.full(n, np.nan)
        zero_frac = np.zeros(n)
        step = 250
        for lo in range(0, n, step):
            hi = min(n, lo + step)
            block = _csr_rows(Xg, lo, hi, self.n_bins)
            for r in range(hi - lo):
                si = src_of[lo + r]
                if si >= 0:
                    sums[si] += block[r]
                sums[-1] += block[r]
            a = block[:, num].mean(axis=1); b = block[:, den].mean(axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                chr8[lo:hi] = np.where((a > 0) & (b > 0), np.log2(a / b), np.nan)
            zero_frac[lo:hi] = (block == 0).mean(axis=1)
            self._step(f"pseudobulk depth: {hi}/{n} cells", 1)
        g.close()
        counts = np.array([len(self.members[s]) for s in SOURCES], dtype=np.float64)
        out["pb_x"] = (sums / counts[:, None]).astype(np.float32)
        out["chr8_score"] = chr8.astype(np.float32)
        out["zero_frac"] = zero_frac.astype(np.float32)
        # ---- allelic pseudobulk (mean abs-dBAF over member cells; zero = no coverage OR balanced)
        self._step("pseudobulk allelic abs-dBAF", 2)
        Hg = self.h5["mod/haplo/X"]
        hs = np.zeros((len(SOURCES), self.n_bins)); hn = np.zeros((len(SOURCES), self.n_bins))
        for lo in range(0, n, step):
            hi = min(n, lo + step)
            block = _csr_rows(Hg, lo, hi, self.n_bins)
            nz = block != 0
            for r in range(hi - lo):
                si = src_of[lo + r]
                if si >= 0:
                    hs[si] += block[r]; hn[si] += nz[r]
                hs[-1] += block[r]; hn[-1] += nz[r]
        out["pb_haplo"] = (hs / counts[:, None]).astype(np.float32)
        out["pb_haplo_nonzero_frac"] = (hn / counts[:, None]).astype(np.float32)
        # ---- spectral component of the joint graph (scanpy diffmap, as in fig1a)
        self._step("WNN spectral component (diffmap)", 3)
        out["wnn_spectral1"] = self._diffmap1()
        # ---- per-cell CN (cluster segmentation, per-cell scaling) as uint8 = 2*state, 255 = missing
        self._step("per-cell CN matrix", 4)
        p = os.path.join(self.cfg.derived_dir, "percell_cn_states.npy")
        if os.path.exists(p):
            cn = np.load(p, mmap_mode="r")
            if cn.shape != (n, self.n_bins):
                raise ValueError(f"percell_cn_states shape {cn.shape} != ({n},{self.n_bins})")
            u8 = np.full(cn.shape, 255, dtype=np.uint8)
            for lo in range(0, n, 500):
                blk = np.asarray(cn[lo:lo + 500], dtype=np.float32)
                ok = np.isfinite(blk)
                v = np.rint(blk * 2)
                u8[lo:lo + 500][ok] = v[ok].astype(np.uint8)
            out["percell_cn_u8"] = u8
        else:
            self.problems.append("percell_cn_states.npy missing; per-cell karyogram unavailable")
        return out

    def _diffmap1(self) -> np.ndarray:
        try:
            import anndata as ad
            import scanpy as sc
            from scipy.sparse import csr_matrix
            g = self.graphs["wnn"]
            n = len(self.cell_ids)
            C = csr_matrix((g["data"], g["indices"], g["indptr"]), shape=(n, n))
            d = self.h5["obsp/distances"]
            D = csr_matrix((d["data"][...], d["indices"][...], d["indptr"][...]), shape=(n, n))
            a = ad.AnnData(obs=pd.DataFrame(index=self.cell_ids))
            a.obsp["connectivities"] = C; a.obsp["distances"] = D
            a.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances", "params": {"method": "umap", "n_neighbors": 20}}
            sc.tl.diffmap(a, n_comps=5, random_state=42)
            return np.asarray(a.obsm["X_diffmap"])[:, 1].astype(np.float32)
        except Exception as e:  # noqa: BLE001
            self.notes.append(f"WNN spectral component not computed ({type(e).__name__}: {e}); fig1a spectral panels unavailable")
            return np.full(len(self.cell_ids), np.nan, dtype=np.float32)

    def _prepare_karyo(self) -> None:
        import gzip
        if self.arrays and "percell_cn_u8" in self.arrays:
            raw = np.ascontiguousarray(self.arrays["percell_cn_u8"]).tobytes()
            self._karyo_gz = gzip.compress(raw, compresslevel=4)

    # ---------------------------------------------------------- outputs
    def status(self) -> dict:
        with self._lock:
            st = dict(self.index_state)
        return {"available": self.available, "h5mu": self.h5mu_path, "problems": self.problems, "notes": self.notes,
                "identity": self.identity, "index": st, "cache": self.cache_path() if self.available else None}

    def meta(self) -> dict:
        assert self.cells is not None
        cells = self.cells.copy()
        for k, v in self.emb.items():
            if v.ndim == 2:
                cells[f"umap_{k}_x"] = v[:, 0]; cells[f"umap_{k}_y"] = v[:, 1]
            else:
                cells[k] = v
        if self.arrays is not None:
            for k in ["chr8_score", "zero_frac", "wnn_spectral1"]:
                if k in self.arrays:
                    cells[k] = self.arrays[k]
        cells = cells.replace({np.nan: None})
        def tbl(name):
            t = self.tables.get(name)
            return None if t is None else t.replace({np.nan: None}).to_dict(orient="records")
        chroms = list(dict.fromkeys(self.var_seq.tolist()))
        sizes = {s: int(len(self.members[s])) for s in SOURCES}
        return {
            "identity": self.identity, "h5mu": self.h5mu_path, "notes": self.notes, "problems": self.problems,
            "sources": SOURCES, "source_sizes": sizes,
            "cluster_keys": CLUSTER_KEYS, "default_cluster_key": "wnn_leiden_0.3",
            "chromosomes": chroms, "pilot_chromosomes": self.manifest.get("params", {}).get("chromosomes", []),
            "manifest": self.manifest, "neighbors_params": self.uns.get("neighbors", {}),
            "chr8_score": CHR8_SCORE,
            "cells": cells.to_dict(orient="records"),
            "breakpoints": tbl("breakpoints"),
            "cluster_segments": tbl("karyogram_cluster_segments"),
            "flank_scores": tbl("small_segment_flank_scores"),
            "screening": tbl("small_segment_scores"),
            "bootstrap": tbl("small_segment_bootstrap"),
            "srd_descriptive": tbl("srd_descriptive"),
            "retained_boundaries": tbl("retained_boundaries"),
            "simplified_segments": tbl("simplified_segments"),
            "initial_boundary_scores": tbl("initial_boundary_scores"),
            "unresolved_segments": tbl("unresolved_segments"),
            "threshold_sensitivity": tbl("threshold_sensitivity"),
            "review_summary": tbl("review_summary"),
            "index_done": self.arrays is not None,
        }

    def pseudobulk_binary(self) -> bytes:
        assert self.arrays is not None
        a = np.concatenate([self.arrays["pb_x"], self.arrays["pb_haplo"], self.arrays["pb_haplo_nonzero_frac"]], axis=0).astype(np.float32)
        return np.ascontiguousarray(a).tobytes()

    def karyogram_gz(self) -> Optional[bytes]:
        return self._karyo_gz

    def cell_detail(self, cell: str, k: int = 15) -> dict:
        i = self.cell_index[cell]
        out = {"cell": cell, "index": int(i), "neighbors": {}}
        for name, g in self.graphs.items():
            a, b = int(g["indptr"][i]), int(g["indptr"][i + 1])
            idx = g["indices"][a:b]; w = g["data"][a:b]
            order = np.argsort(-w)[:k]
            out["neighbors"][name] = [{"cell": str(self.cell_ids[j]), "index": int(j), "w": float(w[o]),
                                       "cluster": str(self.cells["wnn_leiden_0.3"].iloc[j])} for o, j in zip(order, idx[order])]
        sets = {n: set(int(x) for x in self.graphs[n]["indices"][int(self.graphs[n]["indptr"][i]):int(self.graphs[n]["indptr"][i + 1])]) for n in self.graphs}
        out["overlap"] = {"depth_vs_haplo": len(sets["depth"] & sets["haplo"]), "wnn_vs_depth": len(sets["wnn"] & sets["depth"]),
                          "wnn_vs_haplo": len(sets["wnn"] & sets["haplo"]), "n_depth": len(sets["depth"]), "n_haplo": len(sets["haplo"]), "n_wnn": len(sets["wnn"])}
        return out

    def candidates(self, source: str, chrom: str) -> List[dict]:
        if self._audit_path is None:
            return []
        if self._audit is None:
            self._audit = pd.read_csv(self._audit_path, sep="\t")
        t = self._audit[(self._audit["source"] == source) & (self._audit["chromosome"] == chrom)]
        return t.replace({np.nan: None}).to_dict(orient="records")
