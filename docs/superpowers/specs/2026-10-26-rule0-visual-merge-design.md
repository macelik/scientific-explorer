# Rule 0: merge what the tracks hatch

**Status:** design approved in conversation 2026-10-26; revised the same day (flank dropdown → criteria checklist); awaiting written-spec review.
**Builds on:** `2026-10-02-merge-rules-design.md` (branch `merge-rules`, unmerged). Where the two disagree, this spec wins.
**Changes from that spec:** section 1 is no longer display-only. Its criteria and thresholds define **Rule 0**, which is hatched on the tracks and is always part of the merge.

## 1. Goal

User testing of the Integration · segmentation & CN tab found:

1. Colouring by a log2 effect with flank "min |effect| flank" hatches nothing. Segments whose flanks are below the τ and/or SRD/√φ thresholds should be hatched, showing what would merge. "SRD ref" should read SRD/√φ.
2. The flank-consistency hatch should mark only segments whose flanks have opposite SRD/√φ signs, and overlaps with the would-merge hatch should look different.
3. The candidate audit shows no rejection reason for candidates rejected on p-value.
4. The "SRD-screen proposal: hatch + colour" layer is not needed.
5. The merge should act on what is hatched.

## 2. Score source

The tracks' small-segment colouring, hatching, tooltips and click card use **live** flank scores from `POST /api/integration/hatch-scores` (unchanged endpoint), for any chromosome and every internal segment with `end - start ≤ small ≤` bins. The archived small-segment tables are no longer the source for the tracks.

Colour-by options on the tracks:

| Option | Live column | Was (archived id) |
|---|---|---|
| Δlog2FC · median levels | `delta_median_*` | `delta_median` |
| Δlog2FC · mean levels | `delta_mean_*` | `delta_mean`, `log2fc` (SRD-screen) |
| Δlog2FC · two-sided IQR-mean levels | `delta_iqr_mean_*` | `delta_trim` |
| \|SRD\| | `srd_*` | `srd` |
| \|SRD/√φ\| | `srd_phi_*` | `srd_phi` |

A segment is coloured by the colour-by metric's value on that metric's own smaller-|value| flank (a missing side loses; both missing → grey "N/A"). The **flank** dropdown and the **effect metric for flank choice** dropdown are removed (replaced by the criteria checklist, §3).

Saved sessions with an archived colour id are mapped by the "Was" column. Before switching, the plan verifies on the pilot chromosomes that each archived column agrees with its live counterpart on matching segments (same source, chromosome, start/end), and records the agreement (max absolute difference, number of matched segments) in `validation/`. If a pair disagrees beyond floating-point noise, the record says so and the docs note that track values on pilot chromosomes changed.

Stays archived (unchanged, labelled archived): the small-segment diagnostics panel (figs 4, 6, 7, 8), the merge-prototype summary, the "retained reason" label layer (shown only on segments matching an archived row), the candidate audit, and the archived extras in the click card (bootstrap fraction, SRD-screen bootstrap SE), shown when an archived row matches.

## 3. Rule 0

### Inputs (all from section 1)

- `criteria`: a non-empty set, checked in a checklist, from
  - `delta_median` — "min |Δlog2FC median| flank < τ"
  - `delta_mean` — "min |Δlog2FC mean| flank < τ"
  - `delta_iqr_mean` — "min |Δlog2FC IQR-mean| flank < τ"
  - `srd_phi` — "min |SRD/√φ| flank < SRD/√φ ref"

  Default: `delta_mean` and `srd_phi`.
- `combine`: `all` | `any` (default `all`)
- `tau` (τ) > 0, the threshold for every Δlog2FC criterion
- `srdRef` > 0, the threshold for `srd_phi` (UI label "SRD/√φ ref"; default 3.3)
- `maxBins` = the "small ≤" value (integer ≥ 1)
- `mergeSignMismatch`: boolean, default `false` (UI: "merge would-merge segments whose flank signs disagree")
- `mergeFlankConflict`: boolean, default `false` (UI: "merge segments whose criteria pick different flanks, toward the one furthest below its threshold")

### Per internal segment with `end - start ≤ maxBins`

A missing or non-finite score is "missing".

1. **Each checked criterion c, metric X:** its flank `f_c` is the side with smaller |X| (a missing side loses; tie → left; both missing → c fails and has no flank). c **passes** if `|X_{f_c}| < threshold_c` (τ for Δlog2FC criteria, srdRef for `srd_phi`). Its **fraction** is `|X_{f_c}| / threshold_c`.
2. **Test.** `all`: every checked criterion passes. `any`: at least one passes. The criteria that count from here on are the **passing** criteria.
3. **Sign mismatch.** `srd_phi_left * srd_phi_right < 0` (zero or missing → no mismatch). Independent of the criteria.
4. **Flank conflict.** The passing criteria do not all pick the same flank.
5. **Direction.** No conflict → that shared flank. Conflict → the flank of the passing criterion with the smallest fraction (tie → left), used only if `mergeFlankConflict` is true.
6. **Ordering score.** `all`: the largest fraction across checked criteria. `any`: the smallest fraction among passing criteria.
7. **Outcome.**
   - Fails the test → not merged by Rule 0.
   - Passes, and is held back if either applies: sign mismatch with `mergeSignMismatch = false`, or flank conflict with `mergeFlankConflict = false`. A held-back segment is not merged, and it records which reason(s) apply.
   - Passes, not held back → **would merge**, flagged "signs disagree" and/or "flanks disagree" when those apply (their switches being on).
   - A would-merge segment whose direction side has a genomic gap while gap crossing is off → still hatched would-merge, marked **blocked by gap**, not merged in that step, never redirected.

Rule 0 has no on/off switch: everything hatched would-merge (not blocked) takes part in the merge.

### Merge

Rule 0 joins the existing rule engine as a rule with `type: "flank_threshold"` and fixed id `rule0`, always first in the list. User rule cards (`type: "conditions"`, the existing model) follow it and OR with it. Claiming, ordering and the loop are as before: one merge per step, candidate order (rule index, score, segment start), where Rule 0's score is the ordering score above; every score, φ and every criterion's flank choice are recomputed after each merge. A Rule 0 merge joins the segment to its direction flank.

## 4. Tracks: hatching and legend

Rule 0 drawing per segment:

| Case | Drawing | Merges |
|---|---|---|
| would merge | × in Rule 0 colour | yes |
| … signs disagree (switch on) | × plus amber vertical \| overlay | yes |
| … flanks disagree (switch on) | × plus violet horizontal ― overlay | yes |
| held back | amber \| for signs and/or violet ― for flanks, plus ⚠ on the marker | no |
| signs disagree only (fails the test) | amber \| | no |

- A segment that is "would merge" with both disagreements draws both overlays. Held back for both reasons draws both marks.
- Blocked by gap: dotted outline on the hatch and a tooltip note.
- User rule cards keep their own patterns and colours, as now.
- The "signs-disagree hatch" layer checkbox (replacing "flank-consistency hatch (×/‒/) + x/o/ label") toggles only the non-merging "signs disagree only" case. Merging and held-back cases are always drawn. The x/o labels are removed.
- Legend: one entry per case present, with counts, e.g. "would merge: 6 · signs disagree: 2 · flanks disagree: 1 · held back: 3 · signs disagree only: 4 · blocked by gap: 1".
- The marker tooltip lists, per checked criterion: its flank, value, threshold, fraction and pass/fail; then the combine result, sign check, flank-conflict check, direction, outcome and (if any) which user rules also match.

## 5. Section 1 and section 2 UI

**Section 1 (heading "1. Visualization and Rule 0")** keeps colour by, scale, track level, smoothing, small ≤ and τ. It renames SRD ref to "SRD/√φ ref", removes the flank and "effect metric for flank choice" dropdowns, and adds the criteria checklist with its all/any toggle and the two switches. It also removes the "SRD-screen proposal: hatch + colour by min-SRD flank" layer, and the text saying small-segment metrics exist only for the pilot chromosomes.

**Section 2 ("2. Merge rules")** shows a pinned, non-editable, non-deletable, non-movable **Rule 0 card** first: a one-line summary of the section-1 settings and "flags N segments" (would merge, across displayed rows). User rule cards follow. "Run merge" requires that Rule 0 is valid (≥ 1 criterion; τ, srdRef finite > 0; maxBins integer ≥ 1) and every enabled user rule is complete; the user rule list may be empty.

## 6. Candidate audit reasons (issue 3)

For a candidate with `selected_for_significance` true, `accepted` false and an empty `rejection_reasons`, the hover shows a derived reason:

- `p_local ≥ 0.001` → "local permutation p = <p_local> ≥ 0.001"
- else `p_global ≥ 0.001` → "global permutation p = <p_global> ≥ 0.001"
- else "not accepted (no recorded reason)"

Every derived reason ends with "(cutoff inferred: all accepted candidates have p < 0.001)". The archived table and server route are unchanged; this is client-side display.

## 7. API

`POST /api/integration/hatch-merge` body stays `{chrom, source, rules, allow_gap_crossing}`. `rules[0]` must be the Rule 0 object:

```json
{"id": "rule0", "type": "flank_threshold", "name": "Rule 0 — shown on the tracks",
 "criteria": ["delta_mean", "srd_phi"], "combine": "all",
 "tau": 0.2, "srdRef": 3.3, "maxBins": 100,
 "mergeSignMismatch": false, "mergeFlankConflict": false,
 "hatch": {"pattern": "x", "color": "#b91c1c"}}
```

Other rules keep the existing shape, with `type` defaulting to `"conditions"`. Validation (422): missing or misplaced Rule 0 (rules[0] must be `type: "flank_threshold"` with id `rule0`, and no other rule may use either); `criteria` empty, with duplicates, or containing an unknown value; invalid `combine`; `tau` / `srdRef` not finite or ≤ 0; `maxBins` not an integer ≥ 1; non-boolean `mergeSignMismatch` / `mergeFlankConflict`; plus all existing user-rule checks. Steps from Rule 0 report `rule_id: "rule0"`, `side` = the direction flank, `score` = the ordering score.

## 8. Code layout

- `server/merge_rules.py`: Rule 0 validation and evaluation (`evaluate_rule0(row, rule0, allow_gap_crossing)`), wired into `validate_rules`, `evaluate_rules` and `run_rule_merge`.
- `web/src/mergeRules.ts`: TS port of the same, plus the display outcome (case, held-back reasons, per-criterion detail), and `rule0FromSeg(seg)` building Rule 0 from section-1 settings.
- `web/src/integrationStore.ts`: section-1 settings gain `criteria`, `combine`, `mergeSignMismatch`, `mergeFlankConflict`, and drop `flank` and `effectMetric`. Colour-by ids become live. Restore maps archived colour ids (table in §2) and old `consistency` / `proposal` label flags.
- `web/src/components/IntegrationSegmentationView.tsx`: tracks drawn from live rows with Rule 0 cases; legend; tooltip; click card; candidate-audit derived reasons; removal of the SRD-screen proposal layer and the two flank dropdowns.
- `web/src/components/MergeRulesEditor.tsx` / `HatchMergeView.tsx`: pinned Rule 0 card; requests always lead with Rule 0.

## 9. Testing

- **Python unit tests:**
  - each criterion, including missing values and ties
  - `all` and `any`
  - ordering score in both modes
  - sign mismatch, including zero and missing values
  - flank conflict, held back vs. merged by its switch
  - direction toward the smallest fraction
  - both switches together
  - gap-blocked segments not redirected
  - Rule 0 first in priority
  - recompute between steps
  - every new 422 case
- **Parity:** extend `tests/fixtures/merge_rules_parity.json` with Rule 0 cases covering every criterion, both combine modes and every outcome. Assert them in pytest and in `tests/merge_rules_helpers.ts`.
- **TypeScript:**
  - candidate-audit derived reasons
  - session restore mapping archived colour ids, dropping the proposal flag and the old flank settings
  - colour value taken from the colour metric's own smaller-|value| flank
- **Validation record:** archived-vs-live agreement on the pilot chromosomes (§2).
- **Browser** (`tests/browser_hatch_merge_walkthrough.cjs`):
  - a non-pilot chromosome shows Rule 0 hatches
  - changing τ, the criteria, or all/any changes the Rule 0 count
  - each switch moves segments between held back and merging
  - a merge run's history names Rule 0
  - the candidate audit shows a derived p-value reason
  - no console errors
- **Regression:** full backend suite and `npm run build`.

## 10. Scientific framing

Rule 0 is an exploratory diagnostic built from the displayed settings, not a validated decision rule. The merge remains independent of the original bootstrap/chain-guard prototype and does not recompute CN. The candidate-audit cutoff is labelled as inferred, not recorded.
