"""loader.py - load depth X, cluster labels, accepted breakpoints; build
observed and whole-cell-bootstrap per-bin pseudobulks for (source, chromosome).

Bootstrap draws are whole cells WITH replacement inside each cluster,
deterministically derived from a master seed per (cluster, replicate) pair.
The SAME replicate cell selections are reused for every chromosome of that
source and every sensitivity setting; the replicate pseudobulk is the per-bin
mean over the resampled cells (same scale as the observed mean pseudobulk).
"""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd


def replicate_draws(n_cells, n_reps, seed, cluster_pos):
    """Deterministic whole-cell draws for one cluster.

    One ``numpy.random.SeedSequence`` child per (cluster, replicate) via a
    hierarchical seed [master_seed, cluster_pos, replicate]; stable across
    calls and independent of how many clusters exist.
    """
    out = []
    for r in range(n_reps):
        ss = np.random.SeedSequence([int(seed), int(cluster_pos), int(r)])
        rng = np.random.default_rng(ss)
        out.append(rng.integers(0, n_cells, size=n_cells))
    return out


def draw_hash(reps):
    h = hashlib.sha256()
    for arr in reps:
        h.update(np.asarray(arr).astype(np.int64).tobytes())
    return h.hexdigest()


def load_annotations(adata, integrated, cluster_key):
    import muon as mu
    if "cellID" in adata.obs.columns:
        adata.obs.index = [str(c) for c in adata.obs["cellID"]]
    m = mu.read_h5mu(integrated)
    cl = m.obs[cluster_key].astype(str)
    cl.index = [str(c) for c in cl.index]
    cell_to_cluster = {str(c): f"cluster{l}" for c, l in cl.items()}
    missing = [c for c in adata.obs.index if c not in cell_to_cluster]
    if missing:
        raise ValueError(f"{len(missing)} cells absent from cluster assignment")
    labels = np.array([cell_to_cluster[c] for c in adata.obs.index])
    members = {f"cluster{i}": list(np.where(labels == f"cluster{i}")[0])
               for i in range(5)}
    return labels, members


def cluster_order():
    return [f"cluster{i}" for i in range(5)]


def chrom_bin_positions(adata):
    chrom_of = adata.var["seq"].astype(str).to_numpy()
    positions = adata.var["start"].to_numpy(dtype=float)
    return chrom_of, positions


def observed_pseudobulk(X, member_rows, chrom_cols):
    """Per-bin mean over a source's cells (mean-over-cells pseudobulk)."""
    sub = X[member_rows][:, chrom_cols]
    if hasattr(sub, "toarray"):
        return np.asarray(sub.toarray(), dtype=float).mean(axis=0)
    return np.asarray(sub, dtype=float).mean(axis=0)


def replicate_pseudobulks(X, member_rows, chrom_cols, reps, n_rep):
    """(n_rep, n_bins) per-bin means over the resampled whole cells.

    Per-bin mean over sampled rows equals a sparse->dense per-cell matrix
    times a replicate count matrix divided by cluster size; implicit zeros are
    preserved (no pseudocount).
    """
    sub = X[member_rows][:, chrom_cols]
    if hasattr(sub, "toarray"):
        block = np.asarray(sub.toarray(), dtype=float)
    else:
        block = np.asarray(sub, dtype=float)
    n_g = block.shape[0]
    C = np.zeros((n_rep, n_g), dtype=float)
    for r, idx in enumerate(reps):
        C[r, :] = np.bincount(np.asarray(idx, dtype=np.int64), minlength=n_g)
    return (C @ block) / n_g


def small_internal_segments(boundary_bins, n_bins, max_small):
    """Original flat partition and internal segments <= max_small bins."""
    cuts = [0] + sorted(int(b) for b in boundary_bins) + [int(n_bins)]
    segs = [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1)
            if cuts[i + 1] > cuts[i]]
    small = []
    for k in range(1, len(segs) - 1):
        s, e = segs[k]
        if e - s <= max_small:
            small.append((k, s, e))
    return segs, small


def read_breakpoints(path):
    return pd.read_csv(path, sep="\t")


def boundary_bins_for(raw, source, chromosome):
    sub = raw[(raw["source"] == source) & (raw["chromosome"] == chromosome)]
    return sorted(int(b) for b in sub["absolute_bin"])
