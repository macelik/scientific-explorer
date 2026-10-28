# Depth cluster segmentation with pyEpi's genome-wide global test (2026-10-28)

The app's depth cluster breakpoints were made by `run_bootstrap_diagnostics.py`, whose global permutation test sampled from each chromosome's own bins (`recursive_getbp(seq, seq, ...)`). Original pyEpi passes the whole cell (`counts.X[i, :]`). `tools/segment_depth_genome_pool.py` recomputes the six segmentations (cluster0..cluster4 and all) with each pseudobulk's whole genome (25,337 bins) as the pool, at the same settings (k 6, 1000 permutations, alpha 0.001, min child 20, min terminal 50), using the copied recursion that reproduces the original file 814/814 with the chromosome pool. The app defaults now point at the new files. The originals in `integration-story/data/` were not touched.

Run: 58 s on 18 workers; 823 breakpoints; 17,222 audit rows. Record in `integration-story/data/genome_pool/run.json`. The run is deterministic (a second run gave a byte-identical breakpoints file). Derived tables and figures were regenerated with `make_karyograms.py` (per-cell CN calling took about 200 s).

## What changed

| Source | Old (chromosome pool) | New (genome pool) | Shared | Old only | New only |
|---|---|---|---|---|---|
| cluster0 | 133 | 137 | 125 | 8 | 12 |
| cluster1 | 143 | 143 | 134 | 9 | 9 |
| cluster2 | 122 | 128 | 119 | 3 | 9 |
| cluster3 | 164 | 156 | 148 | 16 | 8 |
| cluster4 | 116 | 117 | 111 | 5 | 6 |
| all | 136 | 142 | 130 | 6 | 12 |
| Total | 814 | 823 | 767 | 47 | 56 |

Shared means the same source, chromosome and chromosome-local bin. 767 of 814 old breakpoints (94%) reappear; 47 disappear and 56 are new. Shared breakpoints keep their recursion depth in every case.

Chromosomes with the most changed breakpoints (old-only plus new-only; old to new total): chr17 13 (33 to 46), chr18 12 (34 to 22), chr12 12 (43 to 45), chr7 12 (82 to 86), chr5 11 (66 to 65), chr8 7 (41 to 48), chr10 6 (50 to 56), chr2 5 (54 to 55), chr4 5 (35 to 30), chr6 4 (56 to 56).

Median p_global of the shared breakpoints is 0.000999 (1/1001, the floor for 1000 permutations) both before and after; the mean is the same, so p_global does not separate the two segmentations at this alpha. Full numbers are in `comparison.json`.

This is a finding, not a pass/fail: the genome pool moves a small share of breakpoints, concentrated on a few chromosomes (chr17, chr18, chr12, chr7, chr5).

## Tests and walkthroughs

- `python3 -m pytest -q`: 245 passed. One test stub (`test_cache_identity_changes_with_depth_and_manifest`) needed the new `segmentation_audit` config field; nothing else changed.
- Depth reproduction gate (pinned to the original file): 814 reference, 814 reproduced, 0 missing, 0 extra, 0 mismatched.
- `browser_hatch_merge_walkthrough.cjs`, `browser_cn_recall_walkthrough.cjs` and `browser_allelic_walkthrough.cjs` all PASS on a private server (port 8799) with **no assertion edited**: their assertions compare against values computed in the page, not fixed depth counts. The recorded numbers in the output JSONs did change with the segmentation:

| Record | Value | Old | New |
|---|---|---|---|
| cn-recall | merges | 358 | 367 |
| cn-recall | binsChanged | 12949 | 11949 |
| cn-recall | changedPalette | 15834 | 14768 |
| cn-recall | cells with any change | 3,460 of 3,460 | 3,460 of 3,460 |
| cn-recall | median bins changed | 8.76% | 8.55% |
| cn-recall | best_s change median | -0.0022 [-0.0632, 0.0507] | -0.0029 [-0.0420, 0.0443] |
| cn-recall | bpBefore / bpUndo | 35 | 41 |
| cn-recall | bpAfter / bpExpectedAfter | 21 | 23 |
| cn-recall | zzOriginal | 33 | 36 |
| cn-recall | zzRound1 | 7 | 10 |
| hatch-merge | shapesBaseline | 65 | 90 |
| hatch-merge | shapesWithRule | 335 | 432 |
| hatch-merge | Rule 0 would merge now (r0High, wmOff) | 13 | 14 |
| hatch-merge | held back (hbOff / signs) | 2 | 5 |
| hatch-merge | wmOn | 15 | 19 |
| hatch-merge | signs disagree only / merges last (chr7) | 12 | 13 |
| hatch-merge | small (one chr7 row) | 4 | 3 |
| allelic | depth breakpoints, genome table 1 (cluster0, 2, 3, 4) | 133, 122, 164, 116 | 137, 128, 156, 117 |
| allelic | depth only (cluster0, 2, 3, 4) | 132, 120, 163, 114 | 136, 126, 155, 115 |
| allelic | chr8 depth breakpoints (cluster0, 1, 3) | 6, 6, 6 | 7, 10, 7 |
| allelic | "depth breakpoint" plot-element count | 35 | 41 |

Allelic breakpoint counts, matched counts and median offsets in the allelic tab did not change. `validation/allelic/signals.json` was re-recorded with the same four-signal command and the counts are identical; the jobs were served from `cache/allelic/`, so the seconds in the file are the ones from the first recording.
