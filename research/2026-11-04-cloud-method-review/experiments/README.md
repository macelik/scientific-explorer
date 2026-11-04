# Experiments

Every experiment has a `run.py`, small output tables, a `summary.json` with its parameters,
seeds and environment, and a self-contained `report.html`. Scripts only read the data package;
they write into their own folder (experiment 00 also writes `../input_audit.json`). Large
floating intermediates from experiment 01 go to `$REVIEW_SCRATCH` and are not committed.

| # | Folder | Question | Report |
|---|---|---|---|
| 00 | [00_input_audit](00_input_audit/) | Identities, units, alignment, stage availability | [report](00_input_audit/report.html) |
| 01 | [01_reconstruct_x](01_reconstruct_x/) | Does the pinned source reproduce X? Floating intermediates, amplification, environment sensitivity | [report](01_reconstruct_x/report.html) |
| 02 | [02_zero_clip_audit](02_zero_clip_audit/) | QC-removed bins (incl. CDKN2A/B), rounding/clipping zeros by coverage | [report](02_zero_clip_audit/report.html) |
| 03 | [03_normalization_aggregation](03_normalization_aggregation/) | Estimand, CN thresholds, scale-grid edge, denominator sensitivity (segmentation fixed) | [report](03_normalization_aggregation/report.html) |
| 04 | [04_fixed_region_contrasts](04_fixed_region_contrasts/) | Raw-count contrasts, homogeneity, dosage compatibility for saved windows, called gains/losses, literature loci, random windows | [report](04_fixed_region_contrasts/report.html) |
| 05 | [05_segmentation_calibration](05_segmentation_calibration/) | Exchangeable-null acceptance, sampling share, half-split boundary stability | [report](05_segmentation_calibration/report.html) |
| 06 | [06_allelic_information](06_allelic_information/) | Allele observability, effective imbalance, anchored power, evidence labels | [report](06_allelic_information/report.html) |
| 07 | [07_subwindow_profiles](07_subwindow_profiles/) | Per-bin shape, GC attenuation and gaps in the saved windows | [report](07_subwindow_profiles/report.html) |

## Environment

Python 3.10.20 (uv-managed venv) with numpy 1.24.4, pandas 2.2.0, scipy 1.8.0, h5py 3.6.0,
anndata 0.8.0, pyarrow 23.0.0, scikit-misc 0.5.2 (installed `--no-deps`: its metadata require
numpy ≥ 1.26.4, unlike the recorded handoff environment), natsort 8.0.2, matplotlib < 3.9.
Experiment 01's sensitivity check also used Python 3.13.18 / numpy 2.5.4 / scikit-misc 0.5.3
and Python 3.10 with scikit-misc 0.1.4, 0.2.0 and 0.3.1. Each `summary.json` records the
versions at run time.

```bash
uv venv -p 3.10 ../envs/cnv310
uv pip install -p ../envs/cnv310/bin/python numpy==1.24.4 pandas==2.2.0 scipy==1.8.0 h5py==3.6.0 anndata==0.8.0 mudata==0.2.1 pyarrow==23.0.0 natsort==8.0.2 "matplotlib<3.9" fastapi==0.119.1 uvicorn==0.38.0 pydantic==2.11.10 httpx==0.28.1 pytest==9.1.1 seaborn tqdm biopython
uv pip install -p ../envs/cnv310/bin/python --no-deps scikit-misc==0.5.2
```

## Commands (from the workbench root, after the handoff README's data/source setup)

```bash
export PYEPI_SCIENTIFIC_ROOT="$PWD/../cloud-review-data-2026-11-04"
export PYEPI_SCIENTIFIC_PACKAGE="$PWD/../pyepi-upstream/pyEpiAneufinder"
export REVIEW_SCRATCH=/path/to/scratch   # ~0.7 GB for experiment 01 intermediates
PY="../envs/cnv310/bin/python -I"
$PY research/2026-11-04-cloud-method-review/experiments/01_reconstruct_x/run.py
$PY research/2026-11-04-cloud-method-review/experiments/01_reconstruct_x/env_sensitivity.py --tag py310_np124_skmisc052
$PY research/2026-11-04-cloud-method-review/experiments/01_reconstruct_x/analyze.py
$PY research/2026-11-04-cloud-method-review/experiments/02_zero_clip_audit/run.py
$PY research/2026-11-04-cloud-method-review/experiments/03_normalization_aggregation/run.py
$PY research/2026-11-04-cloud-method-review/experiments/04_fixed_region_contrasts/run.py
$PY research/2026-11-04-cloud-method-review/experiments/05_segmentation_calibration/run.py
$PY research/2026-11-04-cloud-method-review/experiments/06_allelic_information/run.py
$PY research/2026-11-04-cloud-method-review/experiments/07_subwindow_profiles/run.py
$PY research/2026-11-04-cloud-method-review/experiments/00_input_audit/run.py   # reads summaries of 01, 02, 06
```

Runtime on 4 CPUs: under 2 minutes each (experiment 05 is about 80 s).

## Seeds

New-experiment seed 20261104 with numpy PCG64 `default_rng` and named child tuples:
`(20261104, 2)` stratified cell halves (exp 05); `(20261104, 1, L, trial)` null sequences (exp 05);
`(20261104, 4)` random windows (exp 04, reused by 06 and 07); `(20261104, 6, region[, 1])`
allelic simulations and bootstrap (exp 06); `20261104` subsample of 40 cells for a figure (exp 01).
The pinned pyEpi permutation tests use their own fixed `random_state=42`. The 300-cell panel is
the frozen `research/2026-10-29-low-coverage/shared/cells.tsv` membership; no new panel was drawn.
