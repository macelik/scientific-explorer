"""Read-only routes for the multimodal integration story tabs."""
from __future__ import annotations

import asyncio
import json
from concurrent.futures import ThreadPoolExecutor

from pydantic import BaseModel, Field
from typing import List
from .flank_view import generate_view, bootstrap_segment

from fastapi import HTTPException
from fastapi.responses import Response

from .integration import IntegrationStore

_pool = ThreadPoolExecutor(max_workers=2)
_meta_cache: dict = {}


class FlankViewRequest(BaseModel):
    chrom: str
    sources: List[str] = Field(min_length=1, max_length=6)
    max_bins: int = Field(default=100, ge=1, le=10000)


class BootstrapRequest(BaseModel):
    source: str
    chrom: str
    start: int = Field(ge=0)
    end: int = Field(gt=0)


def register_integration(app, get_integration):
    def _ig() -> IntegrationStore:
        ig = get_integration()
        if ig is None:
            raise HTTPException(503, "integration data not loaded")
        return ig

    @app.post("/api/integration/flank-view")
    def flank_view(req: FlankViewRequest):
        ig = _ig()
        if not ig.available or ig.arrays is None:
            raise HTTPException(503, "Pseudobulk index is not ready")
        try:
            return generate_view(ig, req.chrom, req.sources, req.max_bins)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e

    @app.post("/api/integration/flank-bootstrap")
    async def flank_bootstrap(req: BootstrapRequest):
        ig = _ig()
        if not ig.available:
            raise HTTPException(503, "Integration inputs unavailable")
        try:
            return await asyncio.get_running_loop().run_in_executor(
                _pool, bootstrap_segment, ig, req.source, req.chrom, req.start, req.end)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e

    @app.get("/api/integration/status")
    def status():
        ig = get_integration()
        if ig is None:
            return {"available": False, "problems": ["integration store not initialised"], "index": {"status": "idle"}}
        return ig.status()

    @app.post("/api/integration/index")
    def start_index():
        ig = _ig()
        ig.start_index()
        return ig.status()

    @app.get("/api/integration/meta")
    async def meta():
        ig = _ig()
        if not ig.available:
            raise HTTPException(503, "integration inputs unavailable: " + "; ".join(ig.problems))
        key = (ig.identity["hash"], ig.arrays is not None)
        if key not in _meta_cache:
            loop = asyncio.get_running_loop()
            body = await loop.run_in_executor(_pool, lambda: json.dumps(ig.meta()))
            _meta_cache.clear()
            _meta_cache[key] = body
        return Response(content=_meta_cache[key], media_type="application/json")

    @app.get("/api/integration/pseudobulk.bin")
    def pseudobulk():
        ig = _ig()
        if ig.arrays is None:
            raise HTTPException(503, "pseudobulk index not ready")
        return Response(content=ig.pseudobulk_binary(), media_type="application/octet-stream",
                        headers={"X-Layout": "float32 [3 blocks x 6 sources x n_bins]: pb_x, pb_haplo_abs_dbaf, pb_haplo_nonzero_frac; sources cluster0..4,all"})

    @app.get("/api/integration/percell_karyogram.bin")
    def percell_karyogram():
        ig = _ig()
        gz = ig.karyogram_gz()
        if gz is None:
            raise HTTPException(503, "per-cell karyogram not ready")
        return Response(content=gz, media_type="application/octet-stream",
                        headers={"Content-Encoding": "gzip", "X-Encoding": "uint8 = 2*CN state, 255 = missing; rows = h5mu obs order"})

    @app.get("/api/integration/cell/{cell}")
    def cell_detail(cell: str):
        ig = _ig()
        if cell not in ig.cell_index:
            raise HTTPException(404, f"{cell} not in the h5mu")
        return ig.cell_detail(cell)

    @app.get("/api/integration/candidates")
    def candidates(source: str, chrom: str):
        ig = _ig()
        return {"source": source, "chrom": chrom, "rows": ig.candidates(source, chrom)}
