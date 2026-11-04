# Rule 0, revision 2: gap fix, fraction-based direction, final signs phase, Rule 0 copies

**Status:** agreed in conversation 2026-10-26 (after testing revision 1).
**Builds on:** `2026-10-26-rule0-visual-merge-design.md` and `2026-10-02-merge-rules-design.md` (branch `merge-rules`, unmerged). Where they disagree, this spec wins.

## 1. Feedback addressed

1. "Accepted breakpoints" and "colour by recursion depth" overlap. Accepted breakpoints should always be drawn as lines coloured by recursion depth.
2. The dotted red "genomic gap at that flank boundary" mark and the hover line "genomic gap: left right" appear on every segment.
3. The "merge segments whose criteria pick different flanks" switch is unnecessary. Always merge toward the criterion furthest below its threshold. Use the same rule to pick the direction for segments whose flank signs disagree.
4. Add an option to also merge signs-disagree segments that fail Rule 0, as the last category, drawn with the same × hatch in a different colour (black).
5. The hover box shows only the min-flank value. It should show left and right values with pass/fail.
6. The section-2 "+ Add rule…" menu should offer Rule 0 as its first, default entry: an editable copy of the section-1 settings.

## 2. Genomic-gap bug (issue 2)

Bin coordinates in `var/start` / `var/end` are 1-based and inclusive, so contiguous bins satisfy `start[i+1] = end[i] + 1`. `hatch_merge.is_gap` tested `start[right] != end[left]`, which is true at every boundary. On the real data, 549 of 25,336 boundaries are true gaps.

**Fix.** `is_gap` is true iff `start[right] > end[left] + 1`. This is also correct for half-open coordinates (`start = end`), which the existing unit test uses. It matches the frontend's `withGaps`, which tests `start > end + 1 bp`.

**Consequences.**
- Gap flags, gap-blocked merges, the dotted outline, the dotted red flank line and the hover/click-card gap text all become correct without further changes.
- Revision 1 reported "blocked by gap: 13 · would merge now: 0" on chr8 at τ 5 / ref 8. That count came from this bug.
- The validation record gains the before/after count of gap boundaries on the real data.

## 3. Breakpoints (issue 1)

- Accepted breakpoints are always drawn, coloured by recursion depth.
- The "accepted breakpoints" and "colour by recursion depth" layer checkboxes are removed.
- The depth legend is always shown.
- Saved sessions with the old flags restore without error, and the flags are dropped.

## 4. Rule 0 changes

### Removed: `mergeFlankConflict`

The switch and the violet ― marker are removed.

**Direction.** For a segment that passes, the direction is always the flank of the passing criterion with the smallest fraction (ties → left). When all passing criteria pick the same flank, this is that flank. A flank conflict is still computed and is reported in the tooltip only ("criteria pick different flanks → toward <side>"). It never holds a segment back.

`held_back` can now contain only `"signs"`.

The API ignores a `mergeFlankConflict` key if present, so old clients and saved settings do not 422.

### Kept: `mergeSignMismatch`

Behaviour is unchanged: a passing segment whose SRD/√φ flank signs disagree merges when this is on and is held back (⚠) when it is off. Its direction uses the fraction rule above.

### New: `mergeSignsOnly` (boolean, default `false`)

UI label: "also merge signs-disagree segments that fail Rule 0 (last, black ×)".

Applies to an evaluated segment (≤ `maxBins`) that fails the Rule 0 test and has a sign mismatch.
- **Off:** outcome is `signs_only`, display-only, the amber | as now.
- **On:** outcome is `signs_merge`.
  - **Direction:** the flank of the checked criterion with the smallest fraction among criteria that have a value, whether they pass or not (ties → left).
  - **Score:** that fraction.
  - **No criterion has a value:** `blocked: "nonfinite"`.
  - **Gap:** if the direction flank has a genomic gap and gap crossing is off, `blocked: "gap"`, never redirected.

Held-back segments (pass the test, signs disagree, `mergeSignMismatch` off) are not `signs_merge`. They stay held back.

### Merge order: phases

Every evaluated segment claim gets a **phase**:

| Phase | Claims |
|---|---|
| 1 | Rule 0 `would_merge`, and every user rule's match |
| 2 | Rule 0 `signs_merge`, only when no user rule matches the segment. A user rule that matches claims it in phase 1 instead. |

Each step merges the candidate with the lowest `(phase, rule index, score, start)`. Every score and φ are recomputed after each merge, as before. So phase 2 runs only once no phase-1 candidate is left. If a phase-2 merge creates a new phase-1 candidate, that candidate goes next.

A phase-2 step reports:
- `rule_id: "rule0"`
- `phase: 2`
- `rule_name: "<Rule 0 name> · signs disagree, fails Rule 0 (last phase)"`

Phase-1 steps report `phase: 1`.

## 5. Rule 0 copies in section 2 (issue 6)

- **Menu.** The "+ Add rule…" menu lists **"Rule 0 copy (current section-1 settings)"** first.
- **Picking it.** It adds a user rule of `type: "flank_threshold"`, with the same fields as Rule 0 except:
  - it has its own id (not `rule0`) and `name: "Rule 0 copy"`;
  - it has `enabled`;
  - it has no `mergeSignsOnly`, which is ignored on copies;
  - its hatch pattern is `x` and its colour comes from the rule-colour cycle.
- **Editing.** Its card edits criteria, all/any, τ, SRD/√φ ref, size limit (`maxBins`), the sign switch, name, hatch and enabled. It can be moved, duplicated and deleted like any card.
- **Evaluation.** A copy is evaluated with the same function as Rule 0 (`mergeSignsOnly` treated as false). It **matches** when its outcome is `would_merge`. When it claims a segment, its side, score and `blocked` come from that evaluation.
- **Validation (enabled copies, 422).** The same checks as Rule 0, except `mergeSignsOnly`. Its id may not be `rule0`. Only `rules[0]` may have id `rule0`.
- **Existing preset rules** (`type: "conditions"`) are unchanged.

## 6. Tracks

**Drawing per segment**

| Case | Drawing |
|---|---|
| Rule 0 would merge, not blocked | × `#b91c1c` |
| Rule 0 would merge, gap-blocked | dotted outline `#b91c1c` |
| … plus signs disagree | amber \| overlay (merging or held back), as now |
| held back | amber \| plus ⚠, as now |
| `signs_merge`, not blocked | × `#111827` (black) |
| `signs_merge`, blocked | dotted outline `#111827` |
| `signs_only` (option off) | amber \|, when the signs layer is on, as now |

- User rules, including copies, keep their own hatch pattern and colour.
- The dotted red flank-boundary line now appears only at real gaps.

**Legend.**
- Rule 0: "would merge now · blocked by gap · signs disagree · ⚠ held back".
- Then "× (black) signs disagree, fails Rule 0, merges last: N" when the option is on, or "| signs disagree only: N" when it is off.
- The "flanks disagree" entry is removed.

**Tooltip and click card, per checked criterion:**

`min |Δlog2FC mean| flank < 0.2: L −0.123 ✓ · R 0.401 ✗ → left (fraction 0.62) pass`

- Values are signed, and ✓/✗ is |value| < threshold for each side.
- An N/A side shows `N/A`.
- Then: the sign check (SRD/√φ L/R signs), "criteria pick different flanks → toward X" when there is a conflict, and the outcome line ("merges into left", "blocked by gap", "⚠ held back: signs disagree", "fails the test; signs disagree, merges last into right", …).
- Gap text appears only if a side really has a gap.

## 7. Section 1 controls

- Remove the flank-conflict checkbox.
- Add the `mergeSignsOnly` checkbox next to the sign switch.
- Keep the "signs-disagree hatch (amber |) on segments that do not pass Rule 0" layer, which affects only the display-only `signs_only` case.
- Session restore drops `mergeFlankConflict`, `labels.breakpoints` and `labels.depthColor`, and defaults `mergeSignsOnly` to false.
- The pinned Rule 0 card summary shows the sign switch and the signs-only option.
- Its count adds "merges last: N" when the option is on.

## 8. API

`rules[0]`:

```json
{"id":"rule0","type":"flank_threshold","name":"Rule 0 — shown on the tracks",
 "criteria":["delta_mean","srd_phi"],"combine":"all","tau":0.2,"srdRef":3.3,"maxBins":100,
 "mergeSignMismatch":false,"mergeSignsOnly":false,"hatch":{"pattern":"x","color":"#b91c1c"}}
```

- `mergeSignsOnly` must be a boolean on Rule 0 (422 otherwise).
- User rules may be `type: "conditions"` (the default) or `type: "flank_threshold"` (a copy, with `enabled`).
- Steps gain `phase`.

## 9. Testing

**Python:**
- `is_gap` with 1-based inclusive coordinates.
- Direction toward the smallest fraction on a conflict, with no held-back for flanks.
- `mergeFlankConflict` ignored.
- `signs_merge` direction and score from failing criteria, plus the nonfinite and gap cases.
- Phase ordering in a run: a phase-1 candidate goes before a lower-score phase-2 one, and phase-2 merges happen after phase 1 runs out.
- A user rule claiming a would-be `signs_merge` segment in phase 1.
- Copy validation and evaluation.
- 422 for non-boolean `mergeSignsOnly`, and for a copy with id `rule0`.

**Parity fixture:**
- Drop `mergeFlankConflict`, add `mergeSignsOnly`, and add cases for `signs_merge` and a copy.
- Also compare `phase`.

**TypeScript:** the same through parity, plus:
- restore migration (dropped keys);
- `makeRule0Copy`;
- the per-side tooltip line formatter.

**Browser walkthrough:**
- On chr8 at the defaults and at a looser setting, the number of gap-flagged segments is small, not all of them.
- No "flank conflict" control.
- Toggling `mergeSignsOnly` changes the black-× count, and a run's history shows a phase-2 step.
- Adding "Rule 0 copy" from the menu creates an editable card.
- Breakpoint layer checkboxes are absent.
- No console errors.

**Regression:** full pytest, the TS helper suites, and `npm run build`.

## 10. Scientific framing

Unchanged: these are exploratory diagnostics. The signs-only phase merges segments that fail every displayed threshold, solely because their flank SRD/√φ signs disagree. The docs state this plainly.
