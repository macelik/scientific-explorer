"""Experiment 01c: is the near-zero LOESS prediction environment-sensitive?

Refits the pinned get_loess_smoothed for the cells listed in mismatch_bins.tsv plus a fixed
set of 50 cells with any prediction <= 0.05 (first 50 by row order, deterministic) in
whatever interpreter runs this script, and writes fitted values at those cells' near-zero
bins to env_sensitivity_<tag>.tsv. Run once per environment, then compare the files.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common as C  # noqa: E402

OUT = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    p = C.default_paths()
    rf = C.load_pyepi_module("render_fragments", p["package"])
    d = C.H5AD(p["h5ad"])
    cells = pd.read_csv(OUT / "per_cell.tsv.gz", sep="\t")
    mism = pd.read_csv(OUT / "mismatch_bins.tsv", sep="\t")
    rows = sorted(set(mism.row) | set(cells.loc[cells.fit_min <= 0.05, "row"].head(50)))
    out = []
    for r in rows:
        raw = d.rows("layers/raw", r, r + 1)[0]
        fit = np.asarray(rf.get_loess_smoothed(raw, d.gc), float)
        for b in np.flatnonzero(fit <= 0.05):
            out.append({"row": r, "bin": int(b), "raw": raw[b], "fit": fit[b]})
        for b in mism.loc[mism.row == r, "bin"]:
            if fit[b] > 0.05:
                out.append({"row": r, "bin": int(b), "raw": raw[b], "fit": fit[b]})
    df = pd.DataFrame(out).drop_duplicates(["row", "bin"])
    df.to_csv(OUT / f"env_sensitivity_{a.tag}.tsv", sep="\t", index=False)
    env = C.environment()
    C.write_json(OUT / f"env_sensitivity_{a.tag}.json", {"tag": a.tag, "rows": rows, "n": len(df), "environment": env})
    print(a.tag, len(rows), len(df), env["python"], env["packages"]["numpy"], env["packages"]["skmisc"])


if __name__ == "__main__":
    main()
