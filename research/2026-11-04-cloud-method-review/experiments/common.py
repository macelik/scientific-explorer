"""Shared, read-only helpers for the 2026-11-04 cloud method review experiments.

Nothing here writes to the data package, `reference/` or `validation/`. Every
loader reads by barcode name / genomic coordinate and checks alignment.
"""
from __future__ import annotations

import base64
import hashlib
import html
import importlib.util
import io
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REVIEW = HERE.parent
REPO = REVIEW.parent.parent
SHARED = REPO / "research" / "2026-10-29-low-coverage" / "shared"
SEED = 20261104


# ----------------------------------------------------------------------------- paths / identity
def default_paths():
    root = Path(os.environ.get("PYEPI_SCIENTIFIC_ROOT", REPO.parent / "cloud-review-data-2026-11-04"))
    pkg = Path(os.environ.get("PYEPI_SCIENTIFIC_PACKAGE", REPO.parent / "pyepi-upstream" / "pyEpiAneufinder"))
    return {
        "root": root,
        "package": pkg,
        "h5ad": root / "pyepi_results" / "count_matrix.h5ad",
        "h5mu": root / "integration-story" / "data" / "integrated.h5mu",
        "allelic": root / "aligned_matrices",
        "genome_pool": root / "integration-story" / "data" / "genome_pool",
        "scratch": Path(os.environ.get("REVIEW_SCRATCH", "/tmp/cnv-review-scratch")),
    }


def sha256(path, block=8 << 20):
    d = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            d.update(b)
    return d.hexdigest()


def environment():
    pk = {}
    for name in ["numpy", "pandas", "scipy", "h5py", "anndata", "pyarrow", "skmisc", "matplotlib", "statsmodels"]:
        try:
            mod = __import__(name)
            pk[name] = getattr(mod, "__version__", "?")
        except Exception:  # noqa: BLE001
            pk[name] = None
    git = None
    try:
        git = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:  # noqa: BLE001
        pass
    up = None
    try:
        up = subprocess.run(["git", "-C", str(default_paths()["package"].parent), "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    except Exception:  # noqa: BLE001
        pass
    return {"python": platform.python_version(), "platform": platform.platform(), "packages": pk,
            "workbench_commit_at_run": git, "pyepi_upstream_commit": up, "cpu_count": os.cpu_count()}


def load_pyepi_module(name, package_dir=None):
    """Import one pinned upstream module by file path without executing the package __init__."""
    package_dir = Path(package_dir or default_paths()["package"])
    pkg = "pyEpiAneufinder"
    if pkg not in sys.modules:
        import types
        m = types.ModuleType(pkg)
        m.__path__ = [str(package_dir)]
        sys.modules[pkg] = m
    full = f"{pkg}.{name}"
    if full in sys.modules:
        return sys.modules[full]
    spec = importlib.util.spec_from_file_location(full, package_dir / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[full] = mod
    spec.loader.exec_module(mod)
    return mod


# ----------------------------------------------------------------------------- AnnData (h5py, read-only)
class H5AD:
    def __init__(self, path):
        self.path = Path(path)
        self.f = h5py.File(self.path, "r")
        v = self.f["var"]
        cats = v["seq/categories"][...].astype(str)
        self.seq = cats[v["seq/codes"][...]]
        self.start = v["start"][...].astype(np.int64)
        self.end = v["end"][...].astype(np.int64)
        self.gc = v["GC"][...].astype(np.float64)
        self.cells = self.f["obs/cellID"][...].astype(str)
        self.n_cells, self.n_bins = len(self.cells), len(self.seq)
        self._ip = {k: self.f[k + "/indptr"][...] for k in ("X", "layers/raw")}

    def rows(self, key, lo, hi, dtype=np.float64):
        ip = self._ip[key]
        s, e = ip[lo], ip[hi]
        data = self.f[key + "/data"][s:e]
        ind = self.f[key + "/indices"][s:e]
        out = np.zeros((hi - lo, self.n_bins), dtype=dtype)
        for r in range(hi - lo):
            a, b = ip[lo + r] - s, ip[lo + r + 1] - s
            out[r, ind[a:b]] = data[a:b]
        return out

    def row_nnz(self, key):
        return np.diff(self._ip[key])

    def explicit_zeros(self, key, chunk=5_000_000):
        d = self.f[key + "/data"]
        n = 0
        for s in range(0, d.shape[0], chunk):
            blk = d[s:s + chunk]
            n += int((blk == 0).sum())
        return n


def read_cat(group, key):
    g = group[key]
    if isinstance(g, h5py.Group):
        cats = g["categories"][...].astype(str)
        return cats[g["codes"][...]]
    return g[...].astype(str)


def wnn_labels(h5mu_path, key="wnn_leiden_0.3"):
    with h5py.File(h5mu_path, "r") as f:
        obs = f["obs/_index"][...].astype(str)
        lab = read_cat(f["obs"], key)
    return pd.Series(lab, index=obs, name=key)


def group_means(d, key, masks, step=250):
    """Float64 mean of rows of `key` ('X' or 'layers/raw') for each boolean mask (dict name->mask)."""
    names = list(masks)
    sums = np.zeros((len(names), d.n_bins))
    for lo in range(0, d.n_cells, step):
        hi = min(d.n_cells, lo + step)
        blk = d.rows(key, lo, hi)
        for k, nm in enumerate(names):
            sel = masks[nm][lo:hi]
            if sel.any():
                sums[k] += blk[sel].sum(axis=0)
    return {nm: sums[k] / masks[nm].sum() for k, nm in enumerate(names)}


def cell_clusters(d, p=None):
    p = p or default_paths()
    lab = wnn_labels(p["h5mu"])
    missing = set(d.cells) - set(lab.index)
    if missing:
        raise ValueError(f"{len(missing)} AnnData cells missing from h5mu")
    return lab.loc[d.cells].to_numpy()


def stratified_halves(labels, seed_tuple):
    """Split cells into two halves within each label, using np.random.default_rng(seed_tuple)."""
    rng = np.random.default_rng(seed_tuple)
    a = np.zeros(len(labels), bool)
    for lab in sorted(set(labels)):
        idx = np.flatnonzero(labels == lab)
        perm = rng.permutation(idx)
        a[perm[: len(idx) // 2]] = True
    return a


# ----------------------------------------------------------------------------- shared frozen inputs
def shared_cells():
    return pd.read_csv(SHARED / "cells.tsv", sep="\t")


def shared_regions():
    return pd.read_csv(SHARED / "regions.tsv", sep="\t")


def shared_region_bins():
    return pd.read_csv(SHARED / "region_bins.tsv", sep="\t")


# ----------------------------------------------------------------------------- output helpers
def write_json(path, obj):
    def conv(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, Path):
            return str(o)
        raise TypeError(type(o))
    Path(path).write_text(json.dumps(obj, indent=1, default=conv) + "\n")


def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    return base64.b64encode(buf.getvalue()).decode()


CSS = """
:root{--bg:#fbfbfa;--fg:#1d1f21;--muted:#5d6166;--line:#d9dbde;--card:#ffffff;--accent:#2f5d8a;--warn:#8a5a00;--bad:#9b2c2c}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#17191c;--fg:#e6e7e9;--muted:#a2a7ad;--line:#33373c;--card:#1f2226;--accent:#8fb6dd;--warn:#e0b45c;--bad:#ef8f8f}}
:root[data-theme="dark"]{--bg:#17191c;--fg:#e6e7e9;--muted:#a2a7ad;--line:#33373c;--card:#1f2226;--accent:#8fb6dd;--warn:#e0b45c;--bad:#ef8f8f}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
main{max-width:1080px;margin:0 auto;padding:28px 16px 64px}h1{font-size:1.55rem;margin:.2rem 0 .3rem}h2{font-size:1.15rem;margin:2rem 0 .5rem;border-bottom:1px solid var(--line);padding-bottom:.25rem}
.sub{color:var(--muted);margin:0 0 1rem}.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px 16px;margin:12px 0}
.tag{display:inline-block;font-size:.78rem;padding:1px 8px;border-radius:10px;border:1px solid var(--line);color:var(--muted);margin-right:6px}
.warn{color:var(--warn)}.bad{color:var(--bad)}code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.86em}
pre{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:10px;overflow-x:auto;white-space:pre-wrap}
.tw{overflow-x:auto;margin:8px 0}table{border-collapse:collapse;font-size:.86rem;min-width:50%}th,td{border-bottom:1px solid var(--line);padding:4px 8px;text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left}th{color:var(--muted);font-weight:600}img{max-width:100%;height:auto;background:#fff;border-radius:6px;border:1px solid var(--line)}
ul,ol{padding-left:1.2rem}li{margin:.15rem 0}h3{font-size:1rem;margin:1.2rem 0 .4rem}
.toc ol{columns:2;margin:.3rem 0 0}.toc a{color:var(--accent);text-decoration:none}.verdict{border-left:4px solid var(--accent);font-size:1.02rem}
.cap{color:var(--muted);font-size:.86rem;margin:.2rem 0 1rem}.cols{display:grid;grid-template-columns:1fr 1fr;gap:12px}
@media (max-width:700px){.cols{grid-template-columns:1fr}.toc ol{columns:1}}
"""


def html_table(df, floatfmt="{:.4g}", max_rows=200):
    df = df.head(max_rows)
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    body = []
    for _, r in df.iterrows():
        cells = []
        for v in r.values:
            if isinstance(v, (float, np.floating)):
                s = "" if not np.isfinite(v) else floatfmt.format(v)
            else:
                s = str(v)
            cells.append(f"<td>{html.escape(s)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f'<div class="tw"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def write_html_report(path, title, subtitle, sections, provenance, key_findings=None, verdict=None, glossary=None,
                      supports=None, not_supports=None, decision_impact=None, next_steps=None):
    """Self-contained report. Backward compatible: (heading, html) sections plus optional structured blocks.

    key_findings: list of short HTML strings with numbers (shown first, in a summary box).
    verdict: one-sentence bottom line for the reader.
    glossary: dict column/term -> explanation (rendered as a definition table).
    supports / not_supports: lists of what the evidence does / does not support.
    decision_impact: HTML on how the result changes the method recommendation.
    next_steps: list of follow-up checks that would resolve what remains open.
    Sections may be 2-tuples (heading, html) or 3-tuples (heading, html, level) with level 2 or 3.
    """
    blocks = []
    if verdict:
        blocks.append(("Bottom line", f"<div class='card verdict'>{verdict}</div>", 2))
    if key_findings:
        blocks.append(("Key findings", "<div class='card'><ol>" + "".join(f"<li>{k}</li>" for k in key_findings) + "</ol></div>", 2))
    for sec in sections:
        h, frag = sec[0], sec[1]
        lvl = sec[2] if len(sec) > 2 else 2
        blocks.append((h, frag, lvl))
    if supports or not_supports:
        frag = "<div class='cols'>"
        frag += "<div class='card'><h3>Supports</h3><ul>" + "".join(f"<li>{x}</li>" for x in (supports or [])) + "</ul></div>"
        frag += "<div class='card'><h3>Does not support / cannot tell</h3><ul>" + "".join(f"<li>{x}</li>" for x in (not_supports or [])) + "</ul></div>"
        frag += "</div>"
        blocks.append(("What the evidence supports", frag, 2))
    if decision_impact:
        blocks.append(("Decision impact", f"<div class='card'>{decision_impact}</div>", 2))
    if next_steps:
        blocks.append(("Next checks", "<ul>" + "".join(f"<li>{x}</li>" for x in next_steps) + "</ul>", 2))
    if glossary:
        rows = "".join(f"<tr><td><code>{html.escape(k)}</code></td><td style='text-align:left;white-space:normal'>{v}</td></tr>" for k, v in glossary.items())
        blocks.append(("Glossary", f"<div class='tw'><table><thead><tr><th>term / column</th><th style='text-align:left'>meaning</th></tr></thead><tbody>{rows}</tbody></table></div>", 2))
    blocks.append(("Reproduce and provenance", "<pre>" + html.escape(json.dumps(provenance, indent=1, default=str)) + "</pre>", 2))
    ids, toc = [], []
    for i, (h, _, lvl) in enumerate(blocks):
        sid = f"s{i}"
        ids.append(sid)
        if lvl == 2:
            toc.append(f"<li><a href='#{sid}'>{html.escape(h)}</a></li>")
    parts = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
             f"<title>{html.escape(title)}</title><style>{CSS}</style></head><body><main>",
             f"<h1>{html.escape(title)}</h1><p class='sub'>{subtitle}</p>",
             "<p><span class='tag'>exploratory diagnostic</span><span class='tag'>MCF10-A development data</span>"
             "<span class='tag'>not CN truth</span></p>",
             "<nav class='toc card'><b>Contents</b><ol>" + "".join(toc) + "</ol></nav>"]
    for sid, (h, frag, lvl) in zip(ids, blocks):
        tag = "h2" if lvl == 2 else "h3"
        parts.append(f"<section id='{sid}'><{tag}>{html.escape(h)}</{tag}>{frag}</section>")
    parts.append("</main></body></html>")
    Path(path).write_text("".join(parts))


def caption(text):
    """Figure/table caption: what is plotted and how to read it."""
    return f"<p class='cap'>{text}</p>"


def img(b64, alt):
    return f"<img alt='{html.escape(alt)}' src='data:image/png;base64,{b64}'>"
