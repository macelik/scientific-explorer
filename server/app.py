"""FastAPI service for the interactive pyEpiAneufinder explorer."""
from __future__ import annotations

import asyncio
import json
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Literal

import numpy as np
from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .config import APP_DIR, Config
from .data import DEFAULTS, DataStore
from .jobs import JobManager
from .integration import IntegrationStore
from .integration_api import register_integration
from . import science


cfg = Config.from_env()
store: Optional[DataStore] = None
jobs: Optional[JobManager] = None
integration: Optional[IntegrationStore] = None
node_cache: Dict[str, dict] = {}
node_lock = threading.Lock()
node_pool = ThreadPoolExecutor(max_workers=4)

events_state = {"status": "idle", "done": 0, "total": 0, "rows": None, "started": None, "finished": None, "error": None}
events_lock = threading.RLock()

app = FastAPI(title="pyEpiAneufinder explorer")
app.add_middleware(GZipMiddleware, minimum_size=2048)


@app.on_event("startup")
def _startup() -> None:
    global store, jobs
    store = DataStore(cfg)
    ncpu = os.cpu_count() or 4
    jobs = JobManager(cfg.cache_dir, interactive_workers=min(4, ncpu), background_workers=max(1, min(6, ncpu - 4)))
    _load_events_cache()
    global integration
    try:
        integration = IntegrationStore(cfg, store.var_seq, store.var_start, store.var_end, cfg.h5ad, list(store.cells.cellID))
        if integration.available:
            integration.start_index()
        print(f"[explorer] integration: {'available' if integration.available else 'unavailable'} {integration.problems}")
    except Exception as e:  # noqa: BLE001
        print(f"[explorer] integration store failed: {e}")
    print(f"[explorer] data root: {cfg.data_root}")
    print(f"[explorer] cohort: {len(store.cells)} cells x {store.n_bins} bins; loaded in {store.load_seconds:.1f}s")


@app.on_event("shutdown")
def _shutdown() -> None:
    if jobs:
        jobs.shutdown()


def _store() -> DataStore:
    if store is None:
        raise HTTPException(503, "data not loaded")
    return store


def _cell_or_404(cell: str) -> None:
    if cell not in _store().cell_row:
        raise HTTPException(404, f"cell {cell} not in cohort")


def _chrom_or_404(chrom: str):
    s = _store()
    if chrom not in s.chrom_by_name:
        raise HTTPException(404, f"unknown chromosome {chrom}")
    return s.chrom_by_name[chrom]


# ------------------------------------------------------------------- meta
@app.get("/api/meta")
def get_meta():
    return _store().meta()


@app.get("/api/bins.bin")
def get_bins():
    return Response(content=_store().bins_binary(), media_type="application/octet-stream")


@app.get("/api/karyogram.bin")
def get_karyogram():
    s = _store()
    return Response(content=s.karyogram_binary(), media_type="application/octet-stream",
                    headers={"X-Shape": f"{s.cn_u8.shape[0]},{s.cn_u8.shape[1]}", "X-Encoding": "uint8 = 2*CN"})


# ------------------------------------------------------------------- cell
@app.get("/api/cell/{cell}")
def get_cell(cell: str):
    s = _store()
    _cell_or_404(cell)
    t0 = time.time()
    xfull = s.x_row(cell)
    row = s.cell_row[cell]
    meta = s.cells.iloc[row]
    bps = {c.name: s.production_breakpoints(cell, c.name) for c in s.chroms}
    clusters = s.production_clusters(cell)
    out = {
        "cellID": cell,
        "row": int(row),
        "raw_coverage": int(meta.raw_coverage),
        "tercile": meta.tercile,
        "best_s_export": float(meta.best_s),
        "n_segments": int(meta.n_segments),
        "genome": science.genome_baselines(xfull),
        "x": xfull.astype(np.float32).tolist(),
        "production_breakpoints": bps,
        "production_clusters": clusters.tolist(),
        "production_cn5": s.cn[row].astype(float).tolist(),
        "load_ms": round((time.time() - t0) * 1000, 1),
    }
    return out


# ------------------------------------------------------------------- node
class NodeReq(BaseModel):
    cell: str
    chrom: str
    start: int = Field(ge=0)
    end: int = Field(ge=1)


def _node_key(r: NodeReq) -> str:
    return f"node:{_store().identity['hash']}:{r.cell}:{r.chrom}:{r.start}:{r.end}"


def _compute_node_sync(r: NodeReq) -> dict:
    s = _store()
    c = s.chrom_by_name[r.chrom]
    if not (0 <= r.start < r.end <= c.n):
        raise HTTPException(400, f"node bounds [{r.start},{r.end}) outside chromosome [0,{c.n})")
    xfull = s.x_row(r.cell)
    x_chrom = xfull[c.offset:c.offset + c.n]
    x_node = x_chrom[r.start:r.end]
    t0 = time.time()
    res = science.compute_node(x_node)
    starts = s.var_start[c.offset + r.start:c.offset + r.end]
    ends = s.var_end[c.offset + r.start:c.offset + r.end]
    prod = s.production_breakpoints(r.cell, r.chrom)
    prod_set = {p["bp"] for p in prod}
    argmax_chrom = (r.start + res["argmax_b"]) if res["argmax_b"] is not None else None
    res.update({
        "cell": r.cell, "chrom": r.chrom, "start": r.start, "end": r.end, "offset_genome": int(c.offset),
        "x": x_node.astype(np.float32).tolist(),
        "start_bp": starts.astype(int).tolist(),
        "end_bp": ends.astype(int).tolist(),
        "argmax_chrom": argmax_chrom,
        "argmax_in_production": bool(argmax_chrom in prod_set) if argmax_chrom is not None else False,
        "production_breakpoints_in_node": [p for p in prod if r.start < p["bp"] < r.end],
        "compute_ms": round((time.time() - t0) * 1000, 1),
    })
    return res


@app.post("/api/node")
async def post_node(r: NodeReq):
    _cell_or_404(r.cell)
    _chrom_or_404(r.chrom)
    key = _node_key(r)
    with node_lock:
        if key in node_cache:
            out = dict(node_cache[key])
            out["cached"] = True
            return out
    loop = asyncio.get_running_loop()
    res = await loop.run_in_executor(node_pool, _compute_node_sync, r)
    with node_lock:
        node_cache[key] = res
    out = dict(res)
    out["cached"] = False
    return out


# ------------------------------------------------------------------ tests
class TestReq(BaseModel):
    cell: str
    chrom: str
    start: int = Field(ge=0)
    end: int = Field(ge=1)
    b: int = Field(ge=1)
    n_permutations: int = DEFAULTS["n_permutations"]
    seed: int = DEFAULTS["seed"]


@app.post("/api/test")
def post_test(r: TestReq):
    s = _store()
    _cell_or_404(r.cell)
    c = _chrom_or_404(r.chrom)
    if not (0 <= r.start < r.end <= c.n):
        raise HTTPException(400, "node bounds outside chromosome")
    if not (0 < r.b < r.end - r.start):
        raise HTTPException(400, "boundary must leave both children non-empty")
    key = f"test:{s.identity['hash']}:{r.cell}:{r.chrom}:{r.start}:{r.end}:{r.b}:{r.n_permutations}:{r.seed}"
    xfull = s.x_row(r.cell)
    x_node = xfull[c.offset + r.start:c.offset + r.end]
    assert jobs is not None
    return jobs.submit(key, science.run_split_tests, x_node, r.b, xfull, r.n_permutations, r.seed, cfg.package_dir,
                       kind="interactive", meta=r.model_dump())


@app.get("/api/jobs/{jid}")
def get_job(jid: str):
    assert jobs is not None
    j = jobs.get(jid)
    if j is None:
        raise HTTPException(404, "unknown job")
    return j


# ------------------------------------------------------------- CN assign
class AssignReq(BaseModel):
    cell: str
    chrom: str
    boundaries: List[int]
    segment_estimator: Literal['arithmetic', 'iqr_upper', 'iqr_two_sided'] = 'arithmetic'


@app.post("/api/assign_cn")
def post_assign(r: AssignReq):
    s = _store()
    _cell_or_404(r.cell)
    c = _chrom_or_404(r.chrom)
    bad = [b for b in r.boundaries if not (0 < b < c.n)]
    if bad:
        raise HTTPException(400, f"boundaries outside (0,{c.n}): {bad}")
    t0 = time.time()
    xfull = s.x_row(r.cell)
    prod_cl = s.production_clusters(r.cell)
    man_cl = s.spliced_clusters(r.cell, r.chrom, r.boundaries)
    prod = science.assign_cn(xfull, prod_cl.tolist())
    from .segment_estimators import assign_with_estimator
    arithmetic = science.assign_cn(xfull, man_cl.tolist())
    try:
        man = (arithmetic if r.segment_estimator == 'arithmetic' else
               assign_with_estimator(xfull, man_cl, r.segment_estimator))
    except ValueError as e:
        raise HTTPException(422, str(e)) from e
    row = s.cell_row[r.cell]
    export_cn = s.cn[row]
    sl = slice(c.offset, c.offset + c.n)
    changed_elsewhere = []
    for cc in s.chroms:
        if cc.name == r.chrom:
            continue
        a = prod["cn5"][cc.offset:cc.offset + cc.n]
        b = man["cn5"][cc.offset:cc.offset + cc.n]
        nchg = int(np.sum(a != b))
        if nchg:
            changed_elsewhere.append({"chrom": cc.name, "changed_bins": nchg, "n": cc.n})
    if len(set(man_cl.tolist())) != (man_cl.max()):
        raise HTTPException(500, "segment id collision")
    out = {
        "cell": r.cell, "chrom": r.chrom, "boundaries": sorted(set(r.boundaries)),
        "segment_estimator": r.segment_estimator,
        "estimator_scope": "All full-genome segments; original genome baseline and original bin-count weights",
        "arithmetic_manual": {
            "best_s": arithmetic['best_s'], "cn5": arithmetic['cn5'][sl].tolist(),
            "cn_int": arithmetic['cn_int'][sl].tolist(), "cn_cont": arithmetic['cn_cont'][sl].tolist(),
            "segments": science.segment_table(xfull, man_cl, arithmetic, c.offset, c.n, s.var_start, s.var_end),
        },
        "estimator_changed_bins": int(np.count_nonzero(arithmetic['cn5'][sl] != man['cn5'][sl])),
        "estimator_changed_elsewhere": [
            {'chrom':cc.name,'changed_bins':int(np.count_nonzero(arithmetic['cn5'][cc.offset:cc.offset+cc.n] != man['cn5'][cc.offset:cc.offset+cc.n]))}
            for cc in s.chroms if cc.name != r.chrom and np.any(arithmetic['cn5'][cc.offset:cc.offset+cc.n] != man['cn5'][cc.offset:cc.offset+cc.n])
        ],
        "production": {
            "best_s": prod["best_s"],
            "best_s_export": float(s.cells.iloc[row].best_s),
            "matches_export": bool(np.array_equal(prod["cn5"], export_cn)),
            "cn5": prod["cn5"][sl].tolist(),
            "cn_int": prod["cn_int"][sl].tolist(),
            "cn_cont": prod["cn_cont"][sl].tolist(),
            "segments": science.segment_table(xfull, prod_cl, prod, c.offset, c.n, s.var_start, s.var_end),
        },
        "manual": {
            "best_s": man["best_s"],
            "cn5": man["cn5"][sl].tolist(),
            "cn_int": man["cn_int"][sl].tolist(),
            "cn_cont": man["cn_cont"][sl].tolist(),
            "segments": science.segment_table(xfull, man_cl, man, c.offset, c.n, s.var_start, s.var_end),
            "n_segments_genome": int(len(np.unique(man_cl))),
        },
        "changed_bins_this_chrom": int(np.sum(prod["cn5"][sl] != man["cn5"][sl])),
        "changed_elsewhere": changed_elsewhere,
        "genome_baseline_trimmed": science.genome_baselines(xfull)["baseline_trimmed"],
        "compute_ms": round((time.time() - t0) * 1000, 1),
    }
    return out


# ------------------------------------------------------------ event index
def _events_cache_path() -> str:
    return os.path.join(cfg.cache_dir, f"events_{_store().identity['hash']}.json")


def _load_events_cache() -> None:
    p = _events_cache_path()
    if os.path.exists(p):
        with open(p) as f:
            rows = json.load(f)
        with events_lock:
            events_state.update({"status": "done", "rows": rows, "done": len(_store().cells), "total": len(_store().cells)})


def _run_event_index() -> None:
    s = _store()
    assert jobs is not None
    chroms = [{"name": c.name, "offset": c.offset, "n": c.n} for c in s.chroms]
    futures = []
    for cell in s.cells.cellID.values:
        xfull = s.x_row(cell)
        prod = {c.name: [p["bp"] for p in s.production_breakpoints(cell, c.name)] for c in s.chroms}
        futures.append(jobs.background.submit(science.index_cell, cell, xfull, chroms, prod, cfg.package_dir))
    rows: List[dict] = []
    try:
        for f in futures:
            rows.extend(f.result())
            with events_lock:
                events_state["done"] += 1
        with open(_events_cache_path(), "w") as fh:
            json.dump(rows, fh)
        with events_lock:
            events_state.update({"status": "done", "rows": rows, "finished": time.time()})
    except Exception as e:  # noqa: BLE001
        with events_lock:
            events_state.update({"status": "error", "error": f"{type(e).__name__}: {e}"})


@app.post("/api/events/index")
def start_event_index():
    with events_lock:
        if events_state["status"] in ("running", "done"):
            return _events_status()
        events_state.update({"status": "running", "done": 0, "total": len(_store().cells), "started": time.time(),
                             "error": None})
    threading.Thread(target=_run_event_index, daemon=True).start()
    return _events_status()


def _events_status() -> dict:
    with events_lock:
        return {k: v for k, v in events_state.items() if k != "rows"}


@app.get("/api/events/status")
def get_events_status():
    return _events_status()


@app.get("/api/events/candidates")
def get_event_candidates():
    with events_lock:
        if events_state["status"] != "done":
            return {"status": events_state["status"], "rows": None, "done": events_state["done"],
                    "total": events_state["total"]}
        return {"status": "done", "rows": events_state["rows"]}


# --------------------------------------------------------------- sessions
_SAFE = re.compile(r"^[A-Za-z0-9._-]{1,80}$")


class SessionReq(BaseModel):
    name: str
    payload: dict


@app.get("/api/sessions")
def list_sessions():
    out = []
    for fn in sorted(os.listdir(cfg.sessions_dir)):
        if fn.endswith(".json"):
            p = os.path.join(cfg.sessions_dir, fn)
            try:
                with open(p) as f:
                    d = json.load(f)
                out.append({"name": fn[:-5], "saved_at": d.get("saved_at"), "dataset_hash": d.get("dataset", {}).get("hash"),
                            "cell": d.get("state", {}).get("cell"), "chrom": d.get("state", {}).get("chrom"),
                            "size": os.path.getsize(p)})
            except Exception:
                out.append({"name": fn[:-5], "error": "unreadable"})
    return out


@app.post("/api/sessions")
def save_session(r: SessionReq):
    if not _SAFE.match(r.name):
        raise HTTPException(400, "session name must be 1-80 chars of letters, digits, . _ -")
    doc = {"format": "pyepi-explorer-session", "version": 1, "saved_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "dataset": _store().identity, "state": r.payload}
    with open(os.path.join(cfg.sessions_dir, r.name + ".json"), "w") as f:
        json.dump(doc, f)
    return {"ok": True, "name": r.name}


@app.get("/api/sessions/{name}")
def load_session(name: str):
    if not _SAFE.match(name):
        raise HTTPException(400, "bad name")
    p = os.path.join(cfg.sessions_dir, name + ".json")
    if not os.path.exists(p):
        raise HTTPException(404, "no such session")
    with open(p) as f:
        return json.load(f)


@app.get("/api/health")
def health():
    return {"ok": True, "loaded": store is not None, "events": _events_status(), "integration": (integration.status()["index"] if integration else None)}


# --------------------------------------------------------------- frontend
from .experiment_api import register_experiments
register_experiments(app, _store, lambda: jobs, cfg)
register_integration(app, lambda: integration)

DIST = os.path.join(APP_DIR, "web", "dist")
if os.path.isdir(DIST):
    app.mount("/assets", StaticFiles(directory=os.path.join(DIST, "assets")), name="assets")

    @app.get("/{path:path}")
    def spa(path: str):
        target = os.path.join(DIST, path)
        if path and os.path.isfile(target):
            return FileResponse(target)
        return FileResponse(os.path.join(DIST, "index.html"))
