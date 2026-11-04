# Merge rules: separate visualization from user-defined merge rules

**Status:** design approved in conversation 2026-10-02; awaiting written-spec review.
**Supersedes:** the union-of-three-hatch-layers merge eligibility model (commit `3c86329`).

## 1. Goal

The Integration · segmentation & CN tab currently mixes two jobs in one
control row: display settings for the chromosome tracks, and three fixed
hatch layers (transition zone / merge proposal / ambiguous) whose union also
drives the exploratory merge. The merge rules are hard to read and cannot be
combined (e.g. "weak SRD/√φ flank AND small segment").

This redesign:

- Splits the tab into **1. Visualization** (display only, never affects
  merging) and **2. Merge rules** (an ordered, user-built list).
- Lets a user create any number of merge rules. Each rule is a set of
  conditions that must **all** hold (AND); a segment is flagged if **any**
  enabled rule matches (OR of ANDs).
- Gives each rule its own hatch pattern and colour, so what is hatched on the
  tracks is exactly what the merge would act on.

Out of scope: changing how scores (SRD, SRD/√φ, Δlog2FC, pooled φ) are
computed; recomputing CN on merged boundaries; the archived
flank-consistency / SRD-screen label layers (they stay as display layers).

## 2. Rule data model

Rules evaluate **segments**, not boundaries. A flagged segment is merged into
its more similar neighbour.

```ts
type Metric = 'srd' | 'srd_phi' | 'delta_log2fc'
type Estimator = 'mean' | 'median' | 'iqr_mean'   // used only when metric = delta_log2fc

type Condition =
  | { kind: 'weak_flank';    metric: Metric; estimator: Estimator; threshold: number } // min(|L|,|R|) < threshold
  | { kind: 'sign_mismatch'; metric: Metric; estimator: Estimator }                    // L * R < 0
  | { kind: 'ambiguous';     srdMetric: 'srd' | 'srd_phi'; estimator: Estimator }      // weaker SRD-variant flank != weaker Δlog2FC flank
  | { kind: 'max_size';      maxBins: number }                                         // segment length <= maxBins

type HatchPattern = '/' | '\\' | 'x' | '-' | '|' | '.'

interface MergeRule {
  id: string
  name: string
  enabled: boolean
  conditions: Condition[]                               // all must hold
  mergeInto: { metric: Metric; estimator: Estimator }   // defines the "weaker flank"
  hatch: { pattern: HatchPattern; color: string }
}
```

Semantics, exact:

- **Scores.** For an internal segment, L and R are its flank scores against
  the left and right neighbouring segments, as already returned per row by
  `POST /api/integration/hatch-scores` (`srd_left`, `srd_phi_left`,
  `delta_{estimator}_left`, and the `_right` equivalents).
- **Non-finite scores.** A condition whose required score is null or
  non-finite evaluates to false. `weak_flank` requires both L and R finite.
- **`weak_flank`**: `min(|L|, |R|) < threshold` (strict).
- **`sign_mismatch`**: `L * R < 0` (zero on either side is not a mismatch).
- **`ambiguous`**: weaker side by `srdMetric` (smaller absolute value) differs
  from weaker side by Δlog2FC (`estimator`); a tie on either side is not
  ambiguous. Identical to today's `classify_ambiguous` when `srdMetric = srd_phi`.
- **`max_size`**: `end - start <= maxBins` (segment length in bins).
- **Rule match**: an enabled rule with ≥ 1 condition matches a segment if all
  its conditions are true.
- **Claiming rule**: the topmost enabled rule (list order) that matches the
  segment.
- **Merge direction**: under the claiming rule's `mergeInto` metric, the flank
  with the smaller |score|. Ties (|L| = |R|) choose the left flank. If either
  `mergeInto` score is non-finite, the segment is flagged (hatched) but not
  merge-eligible; the tooltip says why.
- **`mergeInto` default**: when a rule is created, it is set to the metric and
  estimator of its first condition that has one; for a size-only rule,
  SRD/√φ. The user can change it.
- **Edge segments** (first/last on the chromosome) are never evaluated, as
  they have only one flank — consistent with `hatch-scores`' internal-segment
  rows.
- **Gaps.** If the chosen merge direction crosses a genomic gap
  (`gap_left`/`gap_right`) and gap crossing is off, the segment is not
  merge-eligible in that step. It is never redirected to the other
  neighbour. The tooltip shows "blocked by gap".

### Presets and session migration

`+ Add rule` offers: blank, and three presets equal to today's layers.
Saved sessions containing the old `hatchLayers` object are migrated on
restore into the corresponding preset rules, preserving metric, estimator,
threshold, and on/off state:

| Old layer | Preset rule | Hatch |
|---|---|---|
| transition zone | `sign_mismatch(metric, estimator)` | `/` |
| merge proposal | `weak_flank(metric, estimator, threshold)` | `x` |
| ambiguous flank preference | `ambiguous(srd_phi, estimator)` | `\` |

Each migrated rule's `mergeInto` = its condition's metric/estimator
(`srd_phi`/estimator for ambiguous). Rule order: transition, proposal,
ambiguous. Old `small_max_bins` was a merge-tool-local setting not stored in
sessions, so it has nothing to migrate.

## 3. UI layout

```
┌─ 1. Visualization ── (display only, never affects merging)
│  chromosome · sources/clusters
│  colour small segments by [metric] on flank [..] · scale [binned/continuous]
│  level basis · smoothing · small-segment size / τ
│  label layers (archived, unchanged): SRD-screen proposal · consistency · reason · values · CN strip · ...
└──
[ Chromosome tracks: colouring + one hatch overlay per enabled, complete rule ]
  legend: one entry per enabled rule — pattern swatch, colour, name, "flags N segments"
┌─ 2. Merge rules ── [+ Add rule ▾]
│  one card per rule: enable checkbox · priority number · name · hatch pattern + colour
│    ALL of: condition rows (kind dropdown + parameters + remove) · [+ condition]
│    merge into the flank with smaller |<metric>|
│    → flags N segments on <chrom> · <source>
│    controls: ↑ ↓ (reorder = priority) · duplicate · delete
│  ☐ allow merging across a genomic gap (off by default)
│  [Run merge] — disabled, with an explanation, until ≥ 1 enabled rule with ≥ 1 condition
│  results: step slider · original vs current tracks · history table (+ "rule" column) · JSON / CSV
└──
```

- Hatches and "flags N" counts update live (client-side) as rules are edited.
- A segment matched by several rules shows each matching rule's pattern
  overlaid. Patterns may repeat across rules.
- Clicking a segment lists, per enabled rule, each condition with its L/R
  values and pass/fail, the claiming rule, merge direction, and any
  "blocked by gap" / "non-finite mergeInto score" note.
- Reorder with ↑/↓ buttons only (no drag-and-drop, no new dependency).
- Any rule or gap-setting change after a run marks the result stale (existing
  frozen-settings signature pattern).
- Rules and the gap setting are saved in sessions; merge results are not
  (unchanged behaviour).

## 4. Merge algorithm and API

`POST /api/integration/hatch-merge`

```json
{ "chrom": "chr8", "source": "cluster4",
  "rules": [ MergeRule, ... ],
  "allow_gap_crossing": false }
```

Replaces the `transition` / `proposal` / `ambiguous` / `small_max_bins`
fields (clean break; the frontend is the only caller). Disabled rules may be
sent and are ignored.

Validation (HTTP 422): no enabled rule; an enabled rule with zero
conditions; unknown condition kind, metric, estimator, or hatch pattern;
`threshold` not finite or ≤ 0; `maxBins` < 1; duplicate rule ids.

Loop (one merge per step, `max_steps` guard retained):

1. Recompute segment summaries and pooled φ for the current partition
   (unconditionally, every step).
2. For each internal segment compute flank scores and evaluate rules top to
   bottom; record the claiming rule.
3. For each claimed segment, determine the merge direction under the claiming
   rule's `mergeInto`; drop it if gap-blocked or non-finite.
4. Pick one candidate by: claiming-rule index (ascending), then weaker-flank
   |score| under `mergeInto` (ascending), then segment start bin (ascending).
5. Remove the boundary on the chosen side; record
   `{step, segments, removed_boundary, merged_interval, rule_id, rule_name, side, score, phi}`.
6. Stop when no candidate remains.

Response shape unchanged apart from the step fields above and
`settings = {rules, allow_gap_crossing}` (echoed from the request, so
summaries and exports always describe the run being shown).

`POST /api/integration/hatch-scores` is unchanged.

## 5. Code layout

- **New `server/merge_rules.py`**: condition evaluation, rule claiming,
  direction choice, and `run_rule_merge(...)`. Reuses scoring primitives
  (`segment_summary`, `pooled_phi`, `flank_score`, `is_gap`, `flat_partition`)
  and the `classify_*`/`weaker_side` functions in `server/hatch_merge.py`.
- **`server/hatch_merge.py`**: remove the layer-based loop (`_segment_status`,
  `_layer_eligibility`, `_eligible_boundaries`, `_validate_layer`,
  `run_merge`); keep scoring, classification, and `hatch_scores`.
- **`server/hatch_merge_api.py`**: new Pydantic models for rules/conditions
  (discriminated union on `kind`); route calls `merge_rules.run_rule_merge`.
- **New `web/src/mergeRules.ts`**: types, `evaluateRules(rows, rules)`
  (per segment: matched rules, claiming rule, direction, block reason),
  presets, `migrateHatchLayers(old)`, default `mergeInto` derivation.
  `web/src/hatchMerge.ts` keeps types for scores/results and the API client;
  its request type changes to the new shape.
- **`web/src/integrationStore.ts`**: `mergeRules: MergeRule[]` and
  `allowGapCrossing: boolean` replace `hatchLayers`; `serialize`/`restore`
  updated, with restore migrating legacy `hatchLayers`.
- **New `web/src/components/MergeRulesEditor.tsx`**: section 2's rule cards.
- **`web/src/components/IntegrationSegmentationView.tsx`**: section 1
  regrouping; remove the old hatch-layer controls; draw per-rule hatch
  overlays, legend, and tooltip lines from `evaluateRules`.
- **`web/src/components/HatchMergeView.tsx`**: becomes the run/results panel
  driven by `mergeRules`/`allowGapCrossing`; history table gains a rule
  column; CSV includes rule id and name.

## 6. Testing

- **Python unit (`tests/test_merge_rules.py`)**: each condition kind including
  non-finite and tie cases; AND within a rule and OR across rules; claiming
  rule = topmost; `mergeInto` direction and left-tie rule; ordering (rule,
  score, start); gap-blocked segments skipped, not redirected; φ recomputed
  each step; each 422 validation case; each migration preset reproduces the
  old layer's classification on the existing `tests/test_hatch_merge.py`
  fixtures. Old layer-loop tests in `tests/test_hatch_merge.py` are removed
  with the code they test; scoring/classification tests stay.
- **Parity**: one shared JSON fixture (`tests/fixtures/merge_rules_parity.json`:
  score rows, rules, expected per-segment matched rules, claiming rule,
  direction, block reason) asserted by pytest against `merge_rules.py` and by
  a TypeScript test (esbuild + node) against `mergeRules.ts`.
- **Browser (`tests/browser_hatch_merge_walkthrough.cjs`, rewritten)**: add a
  rule → hatch appears and "flags N" shows; change a threshold → count
  changes; two rules OR-combined; reorder; run merge → history shows rule
  column; stale notice after edit; session export → restore round-trip;
  restore of a legacy-`hatchLayers` session migrates to three preset rules.
  This also closes the UI-level verification gap noted on commit `3c86329`.
- **Regression**: full backend suite and `npm run build` pass.
- **Docs**: update hatch-merge text in `docs/USAGE.md`, the API section of
  `README.md`, worked example 2 in `docs/SCIENTIFIC_WORKFLOW.md`, and
  `docs/ARCHITECTURE.md`.

## 7. Scientific framing

Unchanged in spirit: the merge remains "an independent exploratory rule, not
the original bootstrap/chain-guard merge prototype," and does not recompute
CN calls. Rules are user-defined diagnostics, not validated decision rules;
UI copy and docs must say so.
