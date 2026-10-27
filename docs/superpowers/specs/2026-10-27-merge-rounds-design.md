# Merge rounds: adopt a merged segmentation as the new base, with a call-noise table

**Status:** agreed in conversation 2026-10-27.
**Builds on:** `2026-10-27-genome-merge-cn-recall-design.md` (branch `merge-rules`, unmerged).

## 1. Goal

A genome re-call shows the effect of one rule set. The user wants to apply a second, different rule set on top of the merged result, and a third after that.

After a re-call, the user can adopt its merged segmentation as the working base. Everything in section 2 then reads from that base:
- the chromosome tracks and the "Segments on chr" table;
- the hatch scores;
- "Run merge";
- the next re-call.

A new table under the tracks shows, round by round, how much call noise is left.

**Agreed in conversation:**
- **Representation:** a recipe of rounds (option 1). The base is "original + round 1 rules + round 2 rules + …", rebuilt on the server on request. It is saved with the session.
- **Karyogram switch:** Original · Before · After · Changed. Before means the current base.
- **Noise table:** under the tracks. It has a this-chromosome / genome-wide scope, one column per round plus After, and definitions shown as visible text above or below the table.
- **Not included:** the spike/dip rule. It is parked until the user has tried a few rounds.

**Background, measured on real data** at τ 0.5 / ref 8 / `any` with mergeSignMismatch off:

| | Original | After 1 round |
|---|---|---|
| Cluster rows: segments | 138–186 | 91–107 |
| Cluster rows: zig-zags | 33–40 | 5–10 |
| Cluster rows: zig-zag median length (bins) | ~38 | ~140 |
| Cluster rows: gain next to loss | 6–25 | 0–2 |
| Per cell (150 sampled): segments | 160.5 | 97.8 |
| Per cell: zig-zags | 31.0 | 8.9 |
| Per cell: gain next to loss | 13.4 | 1.7 |

The same rules applied a second time give 0 merges: the merge already runs until nothing is left to merge. A second round therefore only helps with different rules.

## 2. Data model and server

### Recipe

The integration store holds `rounds: Round[]`:

    Round = {rules, allow_gap_crossing, label, merges, noise}

- **`rules`** is a frozen copy of the full rules list sent with the re-call, including Rule 0's section-1 settings. Editing the rule editor afterwards never changes an adopted round.
- **`merges`** is the round's genome-wide merge count.
- **`noise`** is the snapshot described in §4.
- **Empty list:** `[]` means the original segmentation.
- **Saved with the session:** `rounds` goes in `serialize()`, so it survives reloads and server restarts. A session without the field loads as `[]`.
- **Removing rounds:** only from the end, via "Undo last round" and "Back to original". There is no editing or reordering.

The server receives rounds as `[{rules, allow_gap_crossing}]`. `label`, `merges` and `noise` are client-side display data and are never sent.

### Rebuilding a base

`server/merge_rounds.py` provides `base_breakpoints(inp, rounds)`:
- It starts from the original `raw_accepted_breakpoints.tsv` rows.
- For each round in turn, it runs `cn_recall.merge_genome` with that round's rules on the previous result's boundaries.
- It returns the final breakpoint rows plus, per round, the `pairs` summary.

Notes:
- `merge_genome` gets a breakpoints table to start from. It already reads `inp.breakpoints`, so each round is a shallow copy of the inputs with `breakpoints` replaced.
- Each round's rules are checked with `merge_rules.validate_rules`. An invalid round gives 422 with `round <i>: <message>`, where i counts from 1.
- Each round takes under 1 s. Results are cached by `rounds_hash(rounds)`, the sha256 of canonical JSON, with the 8 most recent bases kept.
- **The merge profile is always `pb_x`**, the original pseudobulk, as for `/hatch-merge`. Merging changes boundaries, not the signal.

### Routes

Each of `/api/integration/hatch-scores`, `/api/integration/hatch-merge` and `/api/integration/recall` gains an optional `rounds` field:
- **Left out or `[]`:** behaviour is byte-identical to today's.
- **Set:** the route reads the base's boundaries instead of the original ones.

#### `POST /api/integration/base`

Body: `{rounds}`. Returns:

```
{hash, n_rounds, breakpoints: [...rows], rounds: [{pairs, merges}],
 cluster_segments: [...rows], noise: {...}, percell: {status: original|ready|running|missing, job: id|null}}
```

- **`breakpoints`** has the columns of `raw_accepted_breakpoints.tsv`, as `merged_breakpoints` does today.
- **`cluster_segments`** are the base's cluster CN calls. They are computed as in re-call stage 2 (float64 member mean, `assign_gainloss_new`), which is instant. For `[]` they are the on-disk `karyogram_cluster_segments.tsv`.
- **`noise`:**
  - For `[]`, the full statistics computed from the on-disk files: cluster level and cells.
  - Otherwise, the cluster-level statistics for the base, plus the cell statistics when its per-cell calls are in memory.

### Per-cell calls of a base

- **For `[]`:** the on-disk `percell_cn_u8`, as now.
- **For a later base:** the per-cell result of the re-call job that produced it.
  - Adopting **pins** that job: `POST /api/integration/recall/{id}/pin`, body `{rounds_hash}`.
  - A pinned job is never evicted by the 3-job limit.
  - Only one base is pinned at a time; pinning a new one unpins the previous one.
- **After a restart, or when no pinned job matches:**
  - `/base` reports `percell.status: missing`.
  - It starts a **base recompute** job. This is a calls-only job: re-call stage 3 run directly on the base's boundaries, with no merge stage and no diff. It therefore never needs an earlier base's per-cell calls, and nothing recurses.
  - Its per-cell result equals that of the re-call that produced the base, because the segmentation and the calling are the same.
  - It pins that job and returns `status: running`.
  - The cell noise for the base is computed from the recompute's results.
  - The client polls the job as it polls re-calls now.

### Re-call changes

- **Hash:** sha256 over `rounds`, `rules`, `allow_gap_crossing`, dataset identity and implementation sources. `merge_rounds.py` joins the implementation sources.
- **Starting point:** stage 1 starts from the base's boundaries.
- **The diff's "before" is the base:**
  - Cluster calls are the base's `cluster_segments`.
  - Per-cell calls are the base's per-cell matrix. For `[]` this is the on-disk matrix; otherwise it is the pinned job's.
  - If the base's per-cell matrix is not available, stages 1–2 run, the job waits for the base recompute, and stage 3's diff then compares against it.
- **Noise:** the job records statistics for its own result (§4). Cluster level is filled after stage 2 and cells after stage 3.
- **Every export header** lists each base round's rules in order, then the job's own rules.

### Provenance

The re-call scope sentence becomes:

> "Exploratory: CN re-called with pyEpi assign_gainloss_new (unmodified) on a rule-merged segmentation (N earlier rounds); originals unchanged."

## 3. UI

### Round strip (section 2, above the rule editor)

    Segmentation base:  Original → Round 1 (τ 0.5 · ref 8 · any, 358 merges) → Round 2 (…, 41 merges)   [Undo last round] [Back to original]

- Each round's label is built from its Rule 0 settings, plus "+ N user rules" when enabled user rules exist. The title text lists the full rules.
- Undo and Back take effect immediately, with no dialog.
- After any change to `rounds`, the store:
  1. fetches `/base`;
  2. refetches hatch scores;
  3. clears a pending re-call that no longer matches the base.
- **Section-2 heading:** "Base: original" or "Base: round N (N rounds of merges on the original)".

### Adopting

The impact summary gets **"Adopt as round N"**, next to the export buttons.

- **When it is enabled:** the re-call is `done` and its `rounds` equal the current `rounds`. Otherwise the button is disabled with the text "computed on a different base".
- **What adopting does:**
  1. appends `{rules, allow_gap_crossing, label, merges, noise}` taken from the re-call's settings and diff;
  2. pins the job;
  3. clears the pending re-call;
  4. refreshes the base.
- **The rule editor is unchanged.** Re-running with the same rules shows 0 merges.

### Karyogram switch: Original · Before · After · Changed

- **Original:** the on-disk calls. Hidden when `rounds` is empty, because it is then the same as Before.
- **Before:** the current base's calls: cluster from `/base`, cells from the base's per-cell matrix.
- **After:** the pending re-call.
- **Changed:** bins that differ between Before and After, drawn as now.
- **When things are disabled or pending:**
  - After and Changed are disabled without a pending re-call.
  - While the base's per-cell calls are recomputing, the per-cell panel in Before shows "re-calling base per-cell CN… N / M".
- **Labels:** hover text and row sub-labels name the compared bases, for example `round 1 <state> → after <state>` and `round 1: best_s … · 94 segments`.

### Tracks

- The breakpoint ticks, the small-segment hatches, the "Segments on chr" table and the hatch scores all show the current base.
- The flank view and the bootstrap/archived views stay on the original archives, and their text says so where it mentions breakpoints.

## 4. Noise table: "Call noise by round"

### Placement

Directly below the "Segments on chr" table.

- **Scope toggle:** "this chromosome" / "genome-wide".
- **Sources:** the ones selected for the tracks.

### Layout

One block per source, with a cluster part and a cells part. There is one column per round (Original, Round 1, …) and a final **After** column while a re-call is pending.

| Level | Measure |
|---|---|
| cluster row | segments |
| cluster row | zig-zags |
| cluster row | zig-zag median length (bins) |
| cluster row | gain next to loss |
| cells (n = …) | segments per cell |
| cells | zig-zags per cell |
| cells | zig-zag median length (bins) |
| cells | gain next to loss per cell |

- Empty values show "—".
- Cell values that are not yet available show "re-calling…".

### Definitions

These are **visible text directly under the table**, not hover-only, in this wording:

- Calls are read per chromosome as the sequence of segment states. Segments with no call are skipped. Gain = state > 1 (weak gain included), loss = state < 1.
- **Zig-zag:** an internal segment whose state differs from both neighbours, while the two neighbours have the same state as each other.
- **Gain next to loss:** a boundary with a gain on one side and a loss on the other. This can be real (for example an unbalanced rearrangement), so it is counted but never merged automatically.
- **Cell values:** the mean per cell over all cells of the cluster. The zig-zag median length is taken over all of those cells' zig-zags.

### Computation

`server/call_noise.py` provides:

- **`seq_noise(states, lengths)`** for one chromosome's segment-state sequence. Non-finite states are dropped first. It returns `{segments, zigzags, zigzag_lengths, gain_loss_adjacent}`.
- **`cluster_noise(cluster_segments_df)`** returns `{source: {chrom | 'genome': {...}}}`.
- **`cells_noise(u8_rows, seg, var_seq)`** works from per-cell u8 rows:
  - Each cell's state per segment is read at the segment's first bin, since states are constant within a segment.
  - It returns, per chromosome and genome-wide, the sums needed for the mean per cell, plus the pooled zig-zag lengths.
  - Re-call workers return it next to `u8`, so the full matrix is not scanned again.

**Statistics format:**

```
noise = {source: {chrom | 'genome': {cluster: {segments, zigzags, zigzag_median_bins, gain_loss_adjacent},
                                      cells: {n_cells, segments_mean, zigzags_mean, zigzag_median_bins, gain_loss_adjacent_mean} | null}}}
```

- The median of an empty list is `null`.
- Genome-wide values sum the per-chromosome counts, and the median is taken over the pooled lengths.

**Original values.** They are computed once from the on-disk files and cached. The real-data check (§6) pins the measured cluster values.

**Snapshot.** Adopting copies the re-call's `noise` into the round entry. Earlier rounds' columns therefore stay filled after a restart, even when their per-cell calls are no longer in memory.

### Export

A fifth export button writes `noise_by_round.tsv`.

- **Columns:** `source, chromosome, level, measure`, then one column per round (`original`, `round1`, …) and `after` when a re-call is pending.
- **Rows:** `chromosome` is a chromosome name or `genome`, and `level` is `cluster` or `cells`.
- **Header:** the same `#` lines as the other exports, listing each round's rules.

## 5. Errors

- **Invalid round:** an invalid round in `rounds` gives 422 `round <i>: <message>` on every route that accepts `rounds`.
- **Failed recompute:** if a base recompute fails, the per-cell panel shows the job's error text. Section 2 keeps working, because the cluster calls, tracks and hatch scores do not depend on per-cell calls.
- **Session from a different dataset:** a session's rounds might not validate against the current server. In that case `/base` returns 422, and the store resets `rounds` to `[]` with the notice "Saved merge rounds could not be applied: <message>".

## 6. Testing

**Python (pytest):**
- **Fold:**
  - `base_breakpoints([])` equals the original breakpoints.
  - `[R]` equals the merged breakpoints of a single re-call with R.
  - `[R, R]` equals `[R]`.
- **Round 2 on synthetic data:** each round's boundaries are a subset of the previous round's, and the per-pair counts equal `len(steps) − 1`.
- **Routes:**
  - Without `rounds`, `/hatch-scores`, `/hatch-merge` and `/recall` give results identical to today's (existing tests unchanged).
  - With `rounds`, they read the base's boundaries.
  - An invalid round gives 422 naming the round.
- **Hash:** a different `rounds` gives a different job; the same `rounds` reuses it.
- **Pinning:**
  - A pinned job survives 3 later finished re-calls.
  - Pinning a second job unpins the first.
- **`seq_noise` on hand-built sequences:**
  - a zig-zag;
  - a staircase, which is not a zig-zag;
  - a weak-gain zig-zag (1 → 1.5 → 1);
  - gain next to loss;
  - skipped non-finite segments;
  - chromosome-edge segments, which are never zig-zags;
  - an empty chromosome.
- **`cells_noise`:** equals `seq_noise` applied per cell on a small hand-built matrix.
- **Real data:** the original's genome-wide cluster noise equals the measured values (cluster2: 144 segments, 33 zig-zags, 6 gain next to loss). Skipped if the data is absent.
- **Restart path:** with no pinned job, `/base` for `[R]` reports `missing`, then `running`. Once done, its per-cell result equals the earlier re-call's, on a fixture.

**TypeScript helper suite:**
- round labels;
- the switch's available modes per state;
- the noise-table row builder ("—" and "re-calling…" cases);
- the `noise_by_round.tsv` writer (header lines, columns);
- `serialize()` round-trip of `rounds`, and loading a session without `rounds`.

**Browser walkthrough**, extending `tests/browser_cn_recall_walkthrough.cjs`. At τ 0.5 / ref 8 / `any`:
1. Re-call and adopt as round 1.
   - The strip shows Round 1.
   - The tracks' breakpoint count drops.
   - The noise table gains a Round 1 column whose genome-wide cluster zig-zags are lower than Original.
   - The definitions text is visible.
2. Re-call again with the same rules: the summary shows 0 merges.
3. Undo the last round: the breakpoint count equals the original.
4. No console errors.

**Regression:** full pytest, all TS suites, `npm run build`, and the hatch-merge walkthrough.

## 7. Out of scope

- The spike/dip rule (parked).
- Editing or reordering earlier rounds.
- Writing merged segmentations to disk, other than exports.
- Moving the flank view or the bootstrap/archived views onto the base.
- τ calibration (see `docs/notes/2026-10-27-tau-calibration.md`).
