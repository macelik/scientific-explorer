"""Numerical validation of the explorer service against direct pyEpiAneufinder calls.
Run with the server up on 127.0.0.1:8765:  python3 tests/validate.py
"""
import json, os, sys, time, urllib.request
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from server.config import Config
from server.data import DataStore
from server import science

API = os.environ.get("API", "http://127.0.0.1:8765")
def get(p):
    with urllib.request.urlopen(API + p) as r: return json.load(r)
def post(p, body):
    req = urllib.request.Request(API + p, data=json.dumps(body).encode(), headers={"content-type": "application/json"})
    with urllib.request.urlopen(req) as r: return json.load(r)

cfg = Config.from_env(); s = DataStore(cfg); sci = s.sci
out = []; fails = 0
def check(name, ok, detail=""):
    global fails
    fails += (not ok)
    out.append(f"[{'PASS' if ok else 'FAIL'}] {name}{(' -- ' + detail) if detail else ''}")
    print(out[-1], flush=True)

# 1. derived clusters == clusters.json for all cohort cells
t0 = time.time()
with open(cfg.clusters_json) as f: cj = json.load(f)
n_ok = sum(1 for c in s.cells.cellID if np.array_equal(np.asarray(cj[c]), s.production_clusters(c)))
check("production clusters rebuilt from breakpoints.csv equal clusters.json", n_ok == len(s.cells), f"{n_ok}/{len(s.cells)} cells, {time.time()-t0:.1f}s")
del cj

# 2. production CN + best_s reproduced for all 300 cells
t0 = time.time(); ok_cn = 0; ok_s = 0
for r in s.cells.itertuples():
    res = science.assign_cn(s.x_row(r.cellID), s.production_clusters(r.cellID).tolist())
    ok_cn += np.array_equal(res["cn5"], s.cn[r.Index]); ok_s += abs(res["best_s"] - r.best_s) < 1e-12
check("assign_gainloss_new reproduces exported five-state CN rows", ok_cn == len(s.cells), f"{ok_cn}/{len(s.cells)}")
check("assign_gainloss_new reproduces exported best_s", ok_s == len(s.cells), f"{ok_s}/{len(s.cells)}, {time.time()-t0:.1f}s")

# 3. node AD arrays and baselines via API vs direct calls
for cell, chrom, start, end in [("TTAGGCTAGGCCGGAA-1", "chr9", 0, 1031), ("TTAGGCTAGGCCGGAA-1", "chr9", 0, 275), ("GAAGGATGTACGCGCA-1", "chr9", 276, 1031), ("TGTGGAGCAGGCTAGA-1", "chr1", 0, 2232), ("TGTGGAGCAGGCTAGA-1", "chr1", 1245, 2232)]:
    c = s.chrom_by_name[chrom]; xfull = s.x_row(cell); x = xfull[c.offset + start:c.offset + end]
    d = np.asarray(sci.seq_dist_ad(x)); t0 = time.time(); node = post("/api/node", dict(cell=cell, chrom=chrom, start=start, end=end)); dt = time.time() - t0
    check(f"node AD {cell} {chrom}[{start},{end}) equals direct seq_dist_ad", np.allclose(node["ad"], d, rtol=0, atol=1e-9) and node["argmax_b"] == int(np.argmax(d)) + 1, f"argmax b={node['argmax_b']} AD={node['argmax_ad']:.4f}; API {dt*1000:.0f} ms (server compute {node['compute_ms']} ms)")
    check(f"node baselines {cell} {chrom}[{start},{end})", abs(node["baseline_trimmed"] - sci.trimmed_mean_iqr(x, lb=False)) < 1e-12 and abs(node["baseline_median"] - np.median(x)) < 1e-12, f"trimmed {node['baseline_trimmed']:.4f} median {node['baseline_median']:.4f}")

# 4. permutation tests via API vs direct calls, argmax and non-argmax splits
def wait_job(j):
    while j["status"] == "running": time.sleep(0.2); j = get(f"/api/jobs/{j['id']}")
    return j
for cell, chrom, start, end, b in [("TTAGGCTAGGCCGGAA-1", "chr9", 0, 1031, 275), ("TTAGGCTAGGCCGGAA-1", "chr9", 0, 1031, 200), ("TTAGGCTAGGCCGGAA-1", "chr9", 0, 275, 225), ("GAAGGATGTACGCGCA-1", "chr9", 276, 1031, 10), ("TGTGGAGCAGGCTAGA-1", "chr1", 0, 2232, 1245)]:
    c = s.chrom_by_name[chrom]; xfull = s.x_row(cell); x = xfull[c.offset + start:c.offset + end]
    obs = sci.dist_ad(x[:b], x[b:]); pl = sci.permutation_test_ad(x[:b], x[b:], obs, 1000, random_state=42); pg = sci.global_permutation_test_ad(xfull, b, len(x) - b, obs, 1000, random_state=42)
    t0 = time.time(); j = wait_job(post("/api/test", dict(cell=cell, chrom=chrom, start=start, end=end, b=b))); dt = time.time() - t0
    r = j["result"]
    check(f"tests {cell} {chrom}[{start},{end}) b={b}{' (non-argmax)' if b == 200 or b == 10 else ''}", r and abs(r["p_local"] - pl) < 1e-15 and abs(r["p_global"] - pg) < 1e-15 and abs(r["observed_ad"] - obs) < 1e-12, f"p_local={pl:.6g} p_global={pg:.6g} obs={obs:.4f}; {'cached' if j.get('cached') else f'{dt*1000:.0f} ms'}")
    prod = {p["bp"]: p for p in s.production_breakpoints(cell, chrom)}
    if start + b in prod:
        pp = prod[start + b]
        check(f"  matches production breakpoints.csv row for bp {start+b}", abs(pp["p_local"] - pl) < 1e-12 and abs(pp["p_global"] - pg) < 1e-12 and abs(pp["ad"] - obs) < 1e-9, f"production p_local={pp['p_local']:.6g}")

# 5. CN splicing integrity
cell, chrom = "TTAGGCTAGGCCGGAA-1", "chr9"
a = post("/api/assign_cn", dict(cell=cell, chrom=chrom, boundaries=[225, 275, 731]))
check("assign_cn with production boundaries reproduces production result", a["production"]["matches_export"] and a["changed_bins_this_chrom"] == 0 and abs(a["manual"]["best_s"] - a["production"]["best_s"]) < 1e-15 and not a["changed_elsewhere"])
b = post("/api/assign_cn", dict(cell=cell, chrom=chrom, boundaries=[190, 260]))
man = s.spliced_clusters(cell, chrom, [190, 260]); prod = s.production_clusters(cell)
same_structure = True
for c in s.chroms:
    if c.name == chrom: continue
    pa = prod[c.offset:c.offset + c.n]; ma = man[c.offset:c.offset + c.n]
    if not np.array_equal(np.diff(pa) != 0, np.diff(ma) != 0): same_structure = False
ids = np.unique(man)
check("spliced clusters keep every other chromosome's segment structure", same_structure)
check("spliced cluster ids unique, contiguous and monotonic", np.array_equal(ids, np.arange(1, len(ids) + 1)) and np.all(np.diff(man) >= 0))
check("manual boundaries [190,260] produce a full-genome refit", b["manual"]["segments"][0]["end_bin"] == 190 and len(b["manual"]["segments"]) == 3, f"best_s prod {b['production']['best_s']:.6f} manual {b['manual']['best_s']:.6f}; changed on chr9 {b['changed_bins_this_chrom']} bins; elsewhere {b['changed_elsewhere']}; {b['compute_ms']} ms")

# 6. literal median zero handling (this cell's genome-wide median is 0)
cd = get(f"/api/cell/{cell}")
check("literal median baseline reported as exactly 0 for a >50%-zero cell (UI shows log2FC undefined, no substitution)", cd["genome"]["baseline_median"] == 0.0 and cd["genome"]["zero_fraction"] > 0.5, f"zero fraction {cd['genome']['zero_fraction']:.3f}, trimmed baseline {cd['genome']['baseline_trimmed']:.4f}")
check("cell endpoint x row equals h5ad row (float32 round-trip)", np.allclose(cd["x"], s.x_row(cell), rtol=1e-6, atol=1e-6))

# 7. karyogram/clustered order sanity
m = get("/api/meta")
check("cluster order is a permutation of cohort rows", sorted(m["cluster_order"]) == list(range(len(s.cells))))
check("event index complete", get("/api/events/status")["status"] == "done", str(get("/api/events/status")))
ev = get("/api/events/candidates")
if ev["rows"]:
    rows = ev["rows"]; root9 = [r for r in rows if r["chrom"] == "chr9" and r["node"] == "root"]
    starts = s.var_start[s.chrom_by_name["chr9"].offset:][:1031] / 1e6
    n2034 = sum(1 for r in root9 if 20 <= starts[r["bp"]] <= 34)
    check("chr9 depth-0 candidates in 20-34 Mb reproduce the documented 203/300", n2034 == 203, f"{n2034}/300 ({100*n2034/300:.1f}%)")
    acc = [r for r in rows if r["accepted"]]
    prod_all = {(r.cell, r.seq, int(r.breakpoint)) for r in s.breakpoints.itertuples()}
    check("every 'accepted' candidate is a production breakpoint", all((r["cell"], r["chrom"], r["bp"]) in prod_all for r in acc), f"{len(acc)} accepted candidates of {len(rows)}")

out.append(f"\n{len(out)-fails}/{len(out)} checks passed")
print(out[-1])
os.makedirs(os.path.join(ROOT, "validation"), exist_ok=True)
with open(os.path.join(ROOT, "validation", "validate_output.txt"), "w") as f: f.write("\n".join(out) + "\n")
sys.exit(1 if fails else 0)
