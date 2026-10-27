# Genome-wide merge and CN re-call, with a before/after karyogram

**Status:** agreed in conversation 2026-10-27.
**Builds on:** `2026-10-26-rule0-revision-2-design.md` and `2026-10-02-merge-rules-design.md` (branch `merge-rules`, unmerged).

## 1. Goal

The merge in section 2 applies to one chromosome and one source. The user wants to see what merging does to the copy-number calls.

One action applies the current merge rules to every chromosome of every source and re-calls CN on the merged segmentation, per cluster and per cell. Both karyograms can then switch between the original calls, the re-called calls, and a view of only the changed bins.

The re-call uses the same pyEpi function and the same procedure as `integration-story/scripts/make_karyograms.py`, unmodified. The original files are never written.

**Agreed in conversation:**
- **Scope:** all chromosomes × all sources. The current chromosome is not special.
- **Comparison:** Before · After · Changed, plus an impact summary.
- **Results:** exportable as files. In memory only, so they do not survive a server restart.
- **Segmentation tracks:** unchanged. They keep showing the original breakpoints.

**Background facts:**
- **The current merge is correct.** At the default Rule 0 (τ 0.2, ref 3.3, `all`), most source/chromosome pairs have 0 or 1 merges. Genome-wide it is 10 merges of 946 segments. At τ 0.5, ref 8, `any` it is 383.
- **CN calling fits one scale factor `best_s` per row over all of that row's segments.** So a merge on one chromosome can change calls on other chromosomes. This is real behaviour, and the docs say so.
- **Cost:** a genome-wide merge takes under 1 s, and cluster re-calling is instant. Per-cell re-calling takes about 0.06 s per cell, about 3.5 min single-process for 3,460 cells.

## 2. Pipeline (server)

`server/cn_recall.py` holds three pure functions, plus a job driver.

### Stage 1: merge

For each source in `SOURCES` (`cluster0`…`cluster4`, `all`) and each chromosome in genome order:

- Run `merge_rules.run_rule_merge` on that pair, with the same inputs `/api/integration/hatch-merge` uses: `pb_x` profile, chromosome-local accepted boundaries, `var_start`/`var_end`, and `chrom_offset`. The result for any pair is therefore identical to a single-chromosome "Run merge".
- Record, per pair:
  - the original and merged segment counts;
  - the number of merges, and how many of them were phase 2;
  - the merged chromosome-local boundaries, which are a subset of the original ones.

The merged breakpoints table has the columns of `raw_accepted_breakpoints.tsv`. Each kept row is copied unchanged from the original row with the same `source`/`chromosome`/`absolute_bin`. Removed rows are dropped.

### Stage 2: cluster CN

For each source:

- **Segment labels.** Built from the merged boundaries with the logic of `make_karyograms.build_seg_labels`: chromosome-local cuts, with segment ids increasing genome-wide.
- **Pseudobulk.** `X[members].mean(axis=0)` in float64 from `count_matrix.h5ad`, as `make_karyograms.py` computes it. This is deliberately not the float32 `pb_x` used for merging.
- **CN call.** `assign_gainloss_new(pseudobulk, seg_labels)`, through `pyepi_adapter`.
- **Output rows** have exactly the columns of `karyogram_cluster_segments.tsv`:

      group, n_cells, best_s, segment, chromosome, bin_start, bin_end, n_bins, continuous_cn, integer_cn, watson, cn_call

  They use the same rounding as `make_karyograms.py`.

### Stage 3: per-cell CN

- **Calls.** Each cell in cluster k gets `assign_gainloss_new(cell_X, seg_labels of cluster k)`. Cells are read from `count_matrix.h5ad`, in the h5mu obs order that `IntegrationStore` already validates.
- **Execution.** The work runs in chunks on `JobManager.background` (a process pool), and progress is reported after each chunk.
- **Output:**
  - a uint8 matrix (n_cells × n_bins) encoded like `percell_cn_u8`: `rint(2·state)`, with 255 for not finite;
  - per-cell `best_s`.

### Diff (server, after stages 2 and 3)

- **Per source:**
  - the number of bins whose `cn_call` changed;
  - that number as a percentage of called bins;
  - `best_s` before and after;
  - segments before and after.
- **Per-cell totals:**
  - the number of cells with ≥ 1 changed bin;
  - the median % of bins changed;
  - `best_s` change (median, min and max).
- **Transitions:** a count of (before state → after state) over changed bins, kept separately for clusters (excluding `all`) and for cells.
- **Before** is the on-disk `karyogram_cluster_segments.tsv` and `percell_cn_u8`.

### Identity

The job hash is sha256 of the canonical JSON of `{rules, allow_gap_crossing}`, the dataset identity hash, and the `merge_rules.py`, `hatch_merge.py` and `cn_recall.py` sources.

- **Reuse:** a job with the same hash is reused whether it is running or finished.
- **Retention:** the server keeps the 3 most recent finished results in memory, and drops the oldest.

## 3. API

### `POST /api/integration/recall`

Body: `{rules, allow_gap_crossing}`, as for `/hatch-merge`. Validated with `merge_rules.validate_rules` (422 with the same messages).

Returns `{id, hash, status}`.

### `GET /api/integration/recall/{id}`

Returns:

```
{id, hash, status: running|done|error, stage: merge|cluster|percell|done, progress: {done, total}, error,
 settings, provenance,
 merge: {pairs: [{source, chrom, original, merged, merges, phase2}], merged_breakpoints: [...rows]} | null,
 cluster_segments: [...rows] | null,
 diff: {sources: [...], cells: {...} | null, transitions: {clusters: [...], cells: [...] | null}} | null,
 percell_ready: bool}
```

`merge` and `cluster_segments` are filled once stage 2 finishes. `diff.cells` and `percell_ready` are filled once stage 3 finishes.

**Provenance:** dataset identity, implementation sha256, and this scope sentence:

> "Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on the rule-merged segmentation; originals unchanged."

Unknown id → 404.

### `GET /api/integration/recall/{id}/percell.bin`

Gzipped uint8 matrix, in the same layout as `/api/integration/percell_karyogram.bin`. Returns 409 until per-cell results are ready.

### Errors

- **pyEpi failure.** Any exception in a stage, including in a worker, ends the job with `status: error` and `error: "<stage> <source>[/<chrom>]: <Type>: <message>"`.
- **Missing per-cell matrix.** If the original per-cell matrix is missing (`percell_cn_u8` absent), stages 1–2 still run, and stage 3 fails with a clear message.

## 4. UI

### Section 2 (Merge rules)

**The button.** "Merge genome & re-call CN" goes next to "Run merge". It sends the same rules as "Run merge", including Rule 0.

**Status line**, polled every 1 s while running:

    merging… → cluster CN ✓ → per-cell 1,240 / 3,460

- It shows done or the error text.
- The button is disabled while a job with the current hash runs.
- If the rules or the gap setting change after a finished run, a notice appears: "Rules changed since this re-call; showing the previous one."

### Karyograms (fig 3c and 3d)

**The switch.** A segmented control, **Before · After · Changed**, is placed above the cluster karyogram and applies to both karyograms.
- Its state lives in the integration store and is not saved in sessions.
- After and Changed are disabled until a re-call exists.
- When per-cell results are not ready, the per-cell panel shows "re-calling per-cell CN… N / M" in After and Changed.

**The three views:**
- **Before:** exactly as now.
- **After:** rows come from `cluster_segments` (cluster) and `percell.bin` (cells). The row sub-labels show the after `best_s` and segment count.
- **Changed:** a bin whose state differs is drawn in its after colour. Every other called bin is drawn `#f1f1f1`, and not-called bins `#e6e6e6`. This needs a per-row mask, or a reserved value in the existing `ClusterKaryogram` palette.

**Hover in After and Changed:**

    before <state> (<label>) → after <state> (<label>)

followed by the cluster hover's existing segment text from the after segmentation, and `best_s <before> → <after>`.

**Impact summary.** Under the switch, shown whenever a re-call exists:
- **Table, one row per source:** segments before → after, merges (of which phase 2), `best_s` before → after, and bins changed (n and %).
- **Cell line:** "cells with any change: N of M · median % bins changed · best_s change median [min, max]".
- **Transitions:** the top 5 transitions for clusters and for cells, as `1 → 1.5: N bins`.
- **Cross-chromosome note:** "A merge on one chromosome can change calls elsewhere: best_s is fitted over the whole genome."

**Export buttons**, beside the summary. Each file is TSV with `#`-prefixed header lines for settings (rules JSON, gap setting) and provenance:

| File | Content |
|---|---|
| `merged_breakpoints.tsv` | the raw-breakpoints columns |
| `recalled_cluster_segments.tsv` | the cluster-segments columns |
| `recalled_percell_best_s.tsv` | `cell, wnn_cluster, best_s_before, best_s_after, bins_changed` |
| `cn_change_summary.tsv` | the per-source table, then the transitions |

## 5. Testing

**Python (pytest):**
- **Identity: cluster.** With rules = Rule 0 at a setting that matches nothing (for example τ 1e-9 and srdRef 1e-9, combine `all`), the merged breakpoints equal the originals. The re-called cluster rows equal `karyogram_cluster_segments.tsv` field for field (best_s to 4 dp). This test needs the real data and is skipped if it is absent.
- **Identity: per-cell.** For the same no-op rules, a sample of 20 cells (fixed seed, across all clusters) re-calls to states equal to `percell_cn_states.npy`, NaN-aware. Skipped if the data is absent.
- **Stage 1 on synthetic data.** Merged boundaries ⊆ original boundaries, and the per-pair counts equal `len(steps) − 1`.
- **Segment labels.** Equal `make_karyograms.build_seg_labels` on a synthetic layout with gaps and empty chromosomes.
- **Diff.** Changed counts and transition tables on a hand-built before/after.
- **API:**
  - 422 for invalid rules;
  - 404 for an unknown id;
  - 409 for `percell.bin` before it is ready;
  - the same hash reuses the job;
  - an error in a stage gives `status: error` with the stage named.

**TypeScript helper suite:** the Changed-mode row builder (mask and colours) and the export TSV writers (header lines, columns).

**Browser walkthrough:** extend `tests/browser_hatch_merge_walkthrough.cjs`. At τ 0.5 / ref 8 / `any`:
- start a genome re-call and wait for done;
- check that After and Changed become enabled;
- check that Changed draws at least one changed cluster bin;
- check that the summary shows more than 0 merges;
- check that the per-cell panel loads after the re-call;
- switch back to Before and check that the rows equal the original;
- no console errors.

**Regression:** full pytest, the TS helper suites, and `npm run build`.

## 6. Out of scope

- Writing results to disk, or restoring them after a restart.
- Changing the segmentation tracks.
- Re-scoring the flank view from merged segments.
- τ calibration (see `docs/notes/2026-10-27-tau-calibration.md`).
