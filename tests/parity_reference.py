"""Generate reference values with numpy/pandas + pyEpi functions for the
browser statistics module. Output: tests/parity_reference.json"""
import json, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from server.config import Config
from server.data import DataStore
from server import science

cfg = Config.from_env()
s = DataStore(cfg)
sci = s.sci
cases = []
for cell, chrom, start, end in [("TTAGGCTAGGCCGGAA-1", "chr9", 0, 1031), ("TTAGGCTAGGCCGGAA-1", "chr9", 0, 275),
                                ("GAAGGATGTACGCGCA-1", "chr9", 276, 1031), ("TGTGGAGCAGGCTAGA-1", "chr1", 0, 2232)]:
    c = s.chrom_by_name[chrom]
    xfull = s.x_row(cell)
    x = xfull[c.offset + start:c.offset + end]
    node = science.compute_node(x)
    d = np.asarray(node["ad"])
    n = len(x)
    b = node["argmax_b"]
    starts = s.var_start[c.offset + start:c.offset + end]
    ends = s.var_end[c.offset + start:c.offset + end]
    windows = []
    for (ws, we) in [(max(0, b - 20), min(n, b + 20)), (0, 40), (n - 40, n), (5, 8), (0, 1), (10, 12), (b - 3, b + 3)]:
        js = [j for j in range(max(ws, 1), min(we, n))]
        adv = d[[j - 1 for j in js]] if js else np.array([])
        xv = x[ws:we]
        def summ(v):
            v = np.asarray(v, dtype=float)
            out = {"n": int(len(v))}
            if len(v) == 0:
                return dict(out, mean=None, median=None, std=None, cv=None, lag1=None)
            mu = float(v.mean()); sd = float(v.std(ddof=0))
            out.update(mean=mu, median=float(np.median(v)), std=sd, cv=(sd / mu if mu != 0 else None))
            if len(v) - 1 >= 3 and np.std(v[:-1]) > 0 and np.std(v[1:]) > 0:
                out["lag1"] = float(np.corrcoef(v[:-1], v[1:])[0, 1])
            else:
                out["lag1"] = None
            return out
        acf = []
        for k in range(1, 11):
            if len(adv) - k >= 3 and np.std(adv[:-k]) > 0 and np.std(adv[k:]) > 0:
                acf.append(float(np.corrcoef(adv[:-k], adv[k:])[0, 1]))
            else:
                acf.append(None)
        windows.append({
            "s": ws, "e": we, "ad_j": js, "ad_stats": summ(adv), "x_stats": summ(xv), "ad_acf": acf,
            "span_bp": int(ends[we - 1] - starts[ws] + 1) if we > ws else 0,
            "retained_bp": int(np.sum(ends[ws:we] - starts[ws:we] + 1)),
        })
    roll = pd.Series(x).rolling(25, center=True, min_periods=8).mean().to_numpy()
    base_g = float(sci.trimmed_mean_iqr(xfull, lb=False))
    with np.errstate(divide="ignore", invalid="ignore"):
        lfc = np.log2(roll / base_g)
    lfc = np.clip(lfc, -4, 4)
    # peaks (static helper) with min_sep 6 Mb and half 20
    sys.path.insert(0, os.path.join(cfg.data_root, "2026-09-15-presentation-story", "scripts"))
    os.chdir(cfg.data_root)
    from _peak_panel import find_spaced_peaks
    mb = starts / 1e6
    peaks = find_spaced_peaks(d, mb, b, n_peaks=5, min_sep_mb=6.0, half_win=20)
    cases.append({
        "cell": cell, "chrom": chrom, "start": start, "end": end, "n": n, "argmax_b": b, "argmax_ad": node["argmax_ad"],
        "baseline_trimmed_node": node["baseline_trimmed"], "baseline_median_node": node["baseline_median"],
        "baseline_trimmed_genome": base_g,
        "windows": windows,
        "rolling": [None if np.isnan(v) else float(v) for v in roll],
        "log2fc": [None if np.isnan(v) else float(v) for v in lfc],
        "peaks": [{"b": int(p["idx"]) + 1, "ad": float(p["ad"]), "is_winner": bool(p["is_winner"])} for p in peaks],
    })
with open(os.path.join(HERE, "parity_reference.json"), "w") as f:
    json.dump(cases, f)
print("wrote", len(cases), "cases")
