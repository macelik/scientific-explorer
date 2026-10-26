# Rule 0: merge what the tracks hatch

**Status:** design approved in conversation 2026-10-26 (sections A and B); awaiting written-spec review.
**Builds on:** `2026-10-02-merge-rules-design.md` (branch `merge-rules`, unmerged). Where the two disagree, this spec wins.
**Changes from that spec:** section 1 ("Visualization") is no longer display-only. Its flank/threshold settings define **Rule 0**, which is hatched on the tracks and is always part of the merge.

## 1. Goal

User testing of the Integration · segmentation & CN tab found:

1. Colouring by a log2 effect with flank "min |effect| flank" hatches nothing. Segments whose chosen flank is below the τ and/or SRD/√φ thresholds should be hatched, showing what would merge. "SRD ref" should read SRD/√φ.
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

Saved sessions with an archived id are mapped by the "Was" column. Before switching, the plan verifies on the pilot chromosomes that each archived column agrees with its live counterpart on matching segments (same source, chromosome, start/end), and records the agreement (max absolute difference, number of matched segments) in `validation/`. If a pair disagrees beyond floating-point noise, the record says so and the docs note that track values on pilot chromosomes changed.

Stays archived (unchanged, labelled archived): the small-segment diagnostics panel (figs 4, 6, 7, 8), the merge-prototype summary, the "retained reason" label layer (shown only on segments matching an archived row), the candidate audit, and the archived extras in the click card (bootstrap fraction, SRD-screen bootstrap SE), shown when an archived row matches.

## 3. Rule 0

### Inputs (all from section 1)

- flank rule `F`: `min_effect` | `min_srd` | `max_effect` | `left` | `right` (UI labels unchanged: min |effect| flank, min |SRD/√φ| flank, max |effect| flank, left, right)
- effect estimator `E`: `median` | `mean` | `iqr_mean` — the colour-by option when it is a Δlog2FC option, otherwise the "effect metric for flank choice" selector. That selector now offers only the three Δlog2FC estimators (it previously also offered SRD and SRD/√φ): τ is a threshold on |Δlog2FC|, and choosing the flank by SRD/√φ is the `min_srd` flank rule.
- `tau` (τ) > 0, on |Δlog2FC|
- `srdRef` > 0, on |SRD/√φ| (UI label "SRD/√φ ref"; default 3.3)
- test mode `M`: `effect` | `srd` | `both` | `either` (default `both`)
- `maxBins` = the "small ≤" value (integer ≥ 1)
- `mergeSignMismatch`: boolean, default `false` (UI: "merge would-merge segments whose flank signs disagree")

### Per internal segment with `end - start ≤ maxBins`

Let `eL, eR` be Δlog2FC(E) and `sL, sR` be SRD/√φ for the left/right flank; a missing or non-finite value is "missing".

1. **Chosen flank f.**
   - `min_effect`: the side with smaller |e|; `min_srd`: the side with smaller |s|; `max_effect`: the side with larger |e|. A missing value loses to a present one for `min_*` (treated as +∞) and for `max_effect` (treated as −∞). Ties choose left. Both missing → no flank; segment not evaluated.
   - `left` / `right`: that side (not evaluated if both its e and s are missing).
2. **Test on f.** `effect`: |e_f| < τ. `srd`: |s_f| < srdRef. `both`: both. `either`: at least one. A missing value fails its part.
3. **Sign mismatch.** `sL * sR < 0` (zero or missing → no mismatch). Independent of f and of the test.
4. **Fraction (ordering score).** `fe = |e_f| / τ`, `fs = |s_f| / srdRef`. `effect` → `fe`; `srd` → `fs`; `both` → `max(fe, fs)`; `either` → the smaller of the parts that pass.
5. **Outcome.**
   - Fails the test → not merged by Rule 0.
   - Passes, no sign mismatch → **would merge**.
   - Passes, sign mismatch, `mergeSignMismatch = true` → **would merge (signs disagree)**.
   - Passes, sign mismatch, `mergeSignMismatch = false` → **held back**: not merged.
   - A would-merge segment whose side f has a genomic gap while gap crossing is off → still hatched would-merge, marked **blocked by gap**, not merged in that step, never redirected.

Rule 0 has no on/off switch: everything hatched would-merge (not blocked) takes part in the merge.

### Merge

Rule 0 joins the existing rule engine as a rule with `type: "flank_threshold"` and fixed id `rule0`, always first in the list. User rule cards (`type: "conditions"`, the existing model) follow it and OR with it. Claiming, ordering and the loop are as before: one merge per step, candidate order (rule index, score, segment start), where Rule 0's score is the fraction above; every score, φ and Rule 0's flank choice are recomputed after each merge. A Rule 0 merge joins the segment to flank f.

## 4. Tracks: hatching and legend

One drawing per segment, by case:

| Case | Drawing | Merges |
|---|---|---|
| would merge | × in Rule 0 colour | yes |
| would merge (signs disagree) | × in Rule 0 colour plus amber vertical \| overlay | yes |
| held back | amber \| plus ⚠ on the marker | no |
| signs disagree only (fails the test) | amber \| | no |

- Blocked by gap: dotted outline on the hatch and a tooltip note.
- User rule cards keep their own patterns and colours, as now.
- The "signs-disagree hatch" layer checkbox (replacing "flank-consistency hatch (×/‒/) + x/o/ label") toggles only the two non-merging amber cases. Merging cases are always drawn. The x/o labels are removed.
- Legend: one entry per case present, with counts, e.g. "would merge: 6 · would merge, signs disagree: 2 · held back (signs disagree): 1 · signs disagree: 4 · blocked by gap: 1".
- The marker tooltip lists: chosen flank and why, e/s on both flanks, each test part pass/fail, fraction, sign check, the outcome, and (if any) which user rules also match.

## 5. Section 1 and section 2 UI

**Section 1 (heading "1. Visualization and Rule 0")** keeps colour by, flank, effect metric for flank choice, scale, track level, smoothing, small ≤ and τ; renames SRD ref to "SRD/√φ ref"; adds "test mode" and the sign switch. The "SRD-screen proposal: hatch + colour by min-SRD flank" layer is removed. Text saying small-segment metrics exist only for the pilot chromosomes is removed.

**Section 2 ("2. Merge rules")** shows a pinned, non-editable, non-deletable, non-movable **Rule 0 card** first: a one-line summary of the section-1 settings and "flags N segments" (would merge + would merge (signs disagree), across displayed rows). User rule cards follow. "Run merge" requires that Rule 0 is valid (τ, srdRef finite > 0; maxBins integer ≥ 1) and every enabled user rule is complete; a user rule list may be empty.

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
 "flank": "min_effect", "effectEstimator": "mean", "tau": 0.2, "srdRef": 3.3,
 "mode": "both", "maxBins": 100, "mergeSignMismatch": false,
 "hatch": {"pattern": "x", "color": "#b91c1c"}}
```

Other rules keep the existing shape, with `type` defaulting to `"conditions"`. Validation (422): missing or misplaced Rule 0 (rules[0] must be `type: "flank_threshold"` with id `rule0`, and no other rule may use either); invalid `flank`, `effectEstimator`, or `mode`; `tau` / `srdRef` not finite or ≤ 0; `maxBins` not an integer ≥ 1; non-boolean `mergeSignMismatch`; plus all existing user-rule checks. Steps from Rule 0 report `rule_id: "rule0"`, `side` = f, `score` = the fraction.

## 8. Code layout

- `server/merge_rules.py`: Rule 0 validation and evaluation (`evaluate_rule0(row, rule0, allow_gap_crossing)`), wired into `validate_rules`, `evaluate_rules` and `run_rule_merge`.
- `web/src/mergeRules.ts`: TS port of the same, plus display outcome (`would_merge` | `would_merge_signs` | `held_back` | `signs_only` | `none`), and `rule0FromSeg(seg)` building Rule 0 from section-1 settings.
- `web/src/integrationStore.ts`: section-1 settings gain `mode` and `mergeSignMismatch`; colour-by ids become live; restore maps archived ids (table in §2) and old `consistency`/`proposal` label flags.
- `web/src/components/IntegrationSegmentationView.tsx`: tracks drawn from live rows with Rule 0 cases; legend; tooltip; click card; candidate-audit derived reasons; removal of the SRD-screen proposal layer.
- `web/src/components/MergeRulesEditor.tsx` / `HatchMergeView.tsx`: pinned Rule 0 card; requests always lead with Rule 0.

## 9. Testing

- Python unit: each flank rule incl. missing values and ties; each mode; fraction per mode; sign mismatch incl. zero/missing; held back vs merged by the switch; gap-blocked not redirected; Rule 0 first in priority; recompute between steps; every new 422 case.
- Parity: extend `tests/fixtures/merge_rules_parity.json` with Rule 0 cases covering every flank rule, mode and outcome; asserted by pytest and by `tests/merge_rules_helpers.ts`.
- TS: candidate-audit derived reasons; session restore mapping archived colour ids and dropping the proposal flag.
- Validation record: archived-vs-live agreement on pilot chromosomes (§2).
- Browser (`tests/browser_hatch_merge_walkthrough.cjs`): a non-pilot chromosome shows Rule 0 hatches; changing τ / mode changes the Rule 0 count; the sign switch moves segments between held back and merging; a merge run's history names Rule 0; candidate audit shows a derived p-value reason; no console errors.
- Regression: full backend suite and `npm run build`.

## 10. Scientific framing

Rule 0 is an exploratory diagnostic built from the displayed settings, not a validated decision rule; the merge remains independent of the original bootstrap/chain-guard prototype and does not recompute CN. The candidate-audit cutoff is labelled as inferred, not recorded.
