# Allelic segmentation: validation record (2026-10-28)

Updated after the app's depth breakpoints switched to the genome-pool segmentation (see `../depth-genome-pool/README.md`). Sections 2 and 4 now show the genome-pool depth counts; the allelic counts and section 3 are unchanged. Section 1 still checks the copied recursion against the original chromosome-pool file.

What was measured on the real data for the "Integration · allelic segmentation" tab. Nothing here is interpreted beyond the numbers.

## 1. Depth reproduction gate

`tests/validate_allelic_depth_reproduction.py` runs the copied constrained `recursive_getbp` on the depth tracks with the **chromosome** global pool (how `raw_accepted_breakpoints.tsv` was made) and compares with that original file (`integration-story/data/raw_accepted_breakpoints.tsv`, pinned in the script; the app now loads the genome-pool file instead). Result in `depth_reproduction.json`:

| Reference breakpoints | Reproduced | Missing | Extra | Mismatched | Seconds |
|---|---|---|---|---|---|
| 814 | 814 | 0 | 0 | 0 | 63.4 (18 workers) |

Settings: k 6, 1000 permutations, alpha 0.001, min child bins 20, min terminal bins 50.

## 2. Walkthrough

`tests/browser_allelic_walkthrough.cjs` against a private server on port 8799: **PASS**.
- 10 plot panels (5 clusters, depth plus allelic).
- Agreement table has 5 rows (one per cluster) genome-wide.
- Export `allelic_vs_depth.tsv` starts with `# signal: abs_dbaf`.
- After switching the signal to signed dBAF the allelic y range is [-0.5, 0.5], and the "earlier run" note appears.
- No page errors or console errors.
- Recorded in `allelic_check_results.json`; screenshot in `allelic.png` (depth segments switched to "original").

`runSeconds` in the JSON is 0 because the |dBAF| job was already cached on disk from earlier attempts of the same script that session, so the Run click returned at once. It is not a compute time.

## 3. Per-signal runs (`signals.json`)

Each signal run at default settings (min AT 10, k 6, 1000 permutations, alpha 0.001, min child 20, min terminal 50, global pool genome) over all six sources and 22 chromosomes, through `POST /api/allelic/segment`. All four ended `done`.

| Signal | cluster0 | cluster1 | cluster2 | cluster3 | cluster4 | all | Total | Gap breakpoints | Seconds |
|---|---|---|---|---|---|---|---|---|---|
| `abs_dbaf` | 7 | 8 | 8 | 11 | 4 | 11 | 49 | 9 | 0 (cached) |
| `dbaf` | 12 | 12 | 8 | 9 | 10 | 11 | 62 | 19 | 16 |
| `shrunk_abs_dbaf` | 6 | 8 | 7 | 9 | 4 | 12 | 46 | 9 | 14 |
| `oriented_all` | 6 | 5 | 6 | 5 | 6 | 18 | 46 | 12 | 14 |

Every signal finds only a handful of breakpoints per source across the whole genome (4 to 18). For comparison, the depth tracks carry 117 to 156 breakpoints per cluster (genome-pool segmentation). The `abs_dbaf` job was served from the cache (see section 2), so its seconds are not a compute time.

## 4. Table 1 for |dBAF|, genome-wide

From `allelic_check_results.json` (`genome`), depth segments "current round" with 0 merge rounds, matching distance N = 3 bins. Median offset is over matched pairs.

| Source | Depth bps | Allelic bps | Matched | Depth only | Allelic only | Median offset |
|---|---|---|---|---|---|---|
| cluster0 | 137 | 7 | 1 | 136 | 6 | 0.0 |
| cluster1 | 143 | 8 | 1 | 142 | 7 | 1.0 |
| cluster2 | 128 | 8 | 2 | 126 | 6 | 1.0 |
| cluster3 | 156 | 11 | 1 | 155 | 10 | 1.0 |
| cluster4 | 117 | 4 | 2 | 115 | 2 | 0.5 |

The same table for chr8 only (`chr8` in the JSON) had 7 to 10 depth breakpoints and 1 allelic breakpoint per source, with a match in cluster4 only.

## 5. Settings used

- Signal `abs_dbaf` for the walkthrough; all four in section 3.
- min AT 10; k 6; permutations 1000; alpha 0.001; min child bins 20; min terminal bins 50.
- Global pool: **genome** (the tab default, as in original pyEpi) for sections 2 to 4. The depth reproduction gate (section 1) used the **chromosome** pool.
- Depth breakpoints compared against are the app's current ones, now made with the same genome-wide global pool as the allelic scan (original pyEpi's test). The matched counts still depend on the matching distance and on how few allelic breakpoints there are; they are not a measure of correctness.
