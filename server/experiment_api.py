"""Experiment routes kept separate from production node/session APIs."""
import hashlib
import json
from pathlib import Path
from typing import Literal

from fastapi import HTTPException
from pydantic import BaseModel, Field

from .experiments import run_experiment, apply_region_edit


class ExperimentReq(BaseModel):
    cell: str
    chrom: str
    start: int = Field(ge=0, strict=True)
    end: int = Field(ge=1, strict=True)
    operation: Literal['multiply', 'add', 'set', 'simulate', 'extend', 'duplicate']
    value: float = Field(default=0.5, allow_inf_nan=False)
    preserve_total: bool = False
    cn_state: Literal['loss', 'base', 'gain'] = 'loss'
    direction: Literal['left', 'right'] = 'right'
    length: int = Field(default=40, ge=1, le=30000, strict=True)
    seed: int = Field(default=42, ge=0, le=4294967295, strict=True)


def register_experiments(app, get_store, get_jobs, cfg):
    @app.post('/api/experiment/preview')
    def preview_experiment(req: ExperimentReq):
        s = get_store()
        if req.cell not in s.cell_row or req.chrom not in s.chrom_by_name:
            raise HTTPException(404, 'Unknown cohort cell or chromosome')
        c = s.chrom_by_name[req.chrom]
        try:
            original = s.x_row(req.cell)
            edited, info = apply_region_edit(original, c.offset, c.n, req.model_dump(), s.cn[s.cell_row[req.cell]])
        except ValueError as e:
            raise HTTPException(400, str(e)) from e
        return {'original': original[c.offset:c.offset+c.n].tolist(),
                'edited': edited[c.offset:c.offset+c.n].tolist(), 'intervention': info}

    @app.post('/api/experiment')
    def submit_experiment(req: ExperimentReq):
        s = get_store()
        if req.cell not in s.cell_row or req.chrom not in s.chrom_by_name:
            raise HTTPException(404, 'Unknown cohort cell or chromosome')
        c = s.chrom_by_name[req.chrom]
        if not 0 <= req.start < req.end <= c.n:
            raise HTTPException(400, 'Region must contain retained bins inside the chromosome')
        if req.operation in ('multiply', 'set') and req.value < 0:
            raise HTTPException(400, 'Multiplier and replacement value must be nonnegative')
        scientific_files = [Path(cfg.package_dir)/name for name in
                            ['distance_statistics.py', 'get_breakpoints.py', 'assign_somy.py']]
        scientific_files.append(Path(__file__).with_name('experiments.py'))
        science_hash = hashlib.sha256(b''.join(p.read_bytes() for p in scientific_files)).hexdigest()[:20]
        payload = {**req.model_dump(), 'dataset_hash': s.identity['hash'], 'science_hash': science_hash}
        key = 'experiment:' + hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
        chroms = [{'name': v.name, 'offset': v.offset, 'n': v.n} for v in s.chroms]
        return get_jobs().submit(key, run_experiment, s.x_row(req.cell).copy(),
            s.production_clusters(req.cell).copy(), chroms, s.var_start, s.var_end,
            payload, cfg.package_dir, kind='interactive', meta=payload)
