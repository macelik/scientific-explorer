# Usage

Open **http://127.0.0.1:8766** (see the root README's [Quick start](../README.md#quick-start)
for how to launch it).

## X experiments

Choose a reference cell from the cohort, then open **X experiments**.

1. Select a named region, reuse an explorer window, draw on the X plot, or enter bin/Mb bounds.
2. Choose a transformation:
   - **Multiply / add / replace** changes the selected X values directly.
   - **Simulate loss / base / gain** samples values from the same cell's original genome-wide bins with the corresponding exact production call. The seed is editable.
   - **Extend with similar values** samples the selected source region into neighbouring bins.
   - **Duplicate ordered pattern** copies that region left or right, repeating its pattern when needed.
3. Inspect the X preview, then select **Run experiment**.
4. Use **Run width series** for 10, 25, 50, 100 and 150-bin scenarios. For extension/duplication, these are additional bins; for other edits, they are event widths centred on the selected region. Actual widths are shown when chromosome edges clip the interval.

Extension overwrites existing neighbouring bins. It does not insert genomic coordinates. Every experiment starts from original X; edits never accumulate silently. The optional total-preservation checkbox rescales the entire genome after editing.

### Reading the result

- The original/edited AD curves show whether the winning candidate moves. Both child curves and exact split tests are expandable.
- **Event-edge diagnostics** show whether each intended boundary was accepted, its rank among candidates in the tested node, and both permutation p-values.
- Four CN rows distinguish original calls, edited values with production boundaries, automatic resegmentation, and an explicitly diagnostic segmentation with event edges imposed.
- Expand **Holmes, Watson and combined calls** to see the number of event bins in each state under each boundary treatment.
- Run history retains the most recent 12 results per cell/chromosome. Inspect earlier runs or export the complete retained history. A changed draft marks the displayed result stale until rerun.

Save/download a normal session to preserve drafts, retained experiments, comparison windows, and manual trees. Reload and upload that JSON to restore them. Export a single experiment as JSON for provenance and arrays, or export regional values as CSV. Sessions are separate from Fable's sessions.

### Interpretation

These are interventions on **GC-corrected, normalized X**, not simulated raw fragments or a rerun of preprocessing. Empirical class sampling is conditioned on existing calls, retains the donor value distribution in expectation, and discards spatial correlation. It is not a validated biological CN simulator. Pattern duplication introduces repeated structure deliberately.

The copied production functions have no hard 50/100-bin size cutoff. Candidate competition, the two significance tests, recursion depth, and downstream fitting determine whether a particular event is resolved. Width-series results describe these specific realizations, not calibrated detection probabilities or a universal size threshold.

Automatic experiments offer **Recursion limit k** from 1 to 10 (default 2). A run tests depths 0 through k−1, with 1,000 permutations, seed 42 for tests, and strict `p < 0.001` for both tests. Increasing k continues only through accepted parents. The immediate children of a rejected root remain explicitly hypothetical and never expand further. The sampling seed is separate from the test seed. Edited global tests use the edited full-genome X row. Only the selected chromosome is automatically resegmented; other chromosomes retain their production boundaries. All CN rows are fitted genome-wide. Inside **Both children and exact split tests**, use **Inspection depth** to view deeper tested nodes. Orange shading marks the modified target clipped to the displayed child intervals; for extend/duplicate it highlights the overwritten bins, not the untouched donor region. Solid vertical lines show accepted breakpoints (including accepted descendants); dashed lines show unaccepted candidates. Blue represents original X and brown edited X. Changing k marks the last result stale; rerun to update it. Depth is retained in sessions and exports. Manual segmentation remains available for user-directed boundaries on original X.

## Manual segmentation / segment estimators

In **Manual segmentation**, place or accept breakpoints, then use **Segment summary**:

- **Arithmetic mean (production)** is the default.
- **IQR mean · upper tail only, keep zeros** keeps values from zero through `Q3 + 1.5 × IQR` in each segment, then averages them.
- **IQR mean · two-sided rule** keeps values between `Q1 − 1.5 × IQR` and `Q3 + 1.5 × IQR`. This can remove zeros when the lower fence is positive.

Press **Finish & assign CN** after choosing the method. The table shows each segment's arithmetic mean, selected estimate, and retained/total bins. For IQR, three CN strips compare production boundaries with arithmetic means, your manual boundaries with arithmetic means, and those same manual boundaries with IQR means. The latter pair isolates the effect of the estimator. Counts report changes to the combined five-state calls; continuous estimates can change even when that state stays the same.

The selected estimator applies to **all genome segments**. Full-genome normalization retains its original upper-tail IQR baseline; segment lengths still weight the scale search even when some values are trimmed. Consequently, both summaries and the shared scale can change calls on other chromosomes. AD curves and split tests remain based on original X. If the selected estimates cannot produce a positive scale grid, the app reports an undefined fit instead of inventing CN calls.

Changing the method marks the previous result outdated. The option and frozen result are saved in sessions, and CSV exports record the method and retained counts. This is an exploratory alternative to production CN assignment, not a validated replacement.

## Integration tabs

**Integration · clusters**: pick a cluster labelling (WNN / depth / allelic at several resolutions),
click a cluster chip or lasso cells in any UMAP; the same cells light up in all three embeddings
and the cluster is preselected in the segmentation tab. Click a cell for its memberships, coverage,
WNN weights, chr8 score and neighbour overlap; cohort cells can be opened in the cell explorer.
The collapsible panels reproduce the story diagnostics (coverage by cluster, PC1/spectral vs
coverage, weights, chr8 group and what a higher resolution adds).

**Integration · segmentation & CN**: the cluster-level and per-cell karyograms come first (click a
cluster row to jump its chromosome into the tracks; click a cell for its card). In the tracks,
choose the metric that colours small segments (read on its own smaller-|value| flank), switch between
binned bands and a continuous scale, and add the label layers you want to talk about (accepted
breakpoints, recursion-depth colour, small segments, signs-disagree hatch, retained reason (archived),
length + value labels, cluster CN strip, candidate audit, allelic). Click a shaded segment or a
point in the diagnostics to open the segment's full flank table. τ and SRD/√φ ref are Rule 0's
thresholds (see the section below).

### Visualization and merge rules

The Chromosome tracks panel starts with **1. Visualization and Rule 0**. It holds the colouring, scale, track level, smoothing, "small ≤" and τ settings, and the criteria that define Rule 0.

**Colouring.** Small internal segments (at most "small ≤" bins) are coloured by live flank scores computed on demand for any chromosome, not only the pilot chromosomes. The colour-by options are Δlog2FC (median, mean or two-sided IQR-mean levels), |SRD| and |SRD/√φ|. Each segment takes the chosen metric's value on that metric's own smaller-|value| flank; a missing side loses, and a segment with both sides missing is grey "N/A". The tracks used to read these values from the archived small-segment tables for the pilot chromosomes. The live and archived values agree to within about 5e-05 (the archived tables are rounded to four decimals) on all 162 flank and 196 screening rows of the pilot chromosomes; see [the agreement record](../validation/rule0/README.md). The archived diagnostics panel, merge-prototype summary, "retained reason" label layer, candidate audit and the archived extras in the click card (bootstrap fraction, SRD-screen bootstrap SE) stay archived and are labelled as such.

**Rule 0.** The same section defines Rule 0, the rule hatched on the tracks and always part of the merge. It has:

- A **criteria checklist**, at least one of: min |Δlog2FC median| flank < τ; min |Δlog2FC mean| flank < τ; min |Δlog2FC IQR-mean| flank < τ; min |SRD/√φ| flank < SRD/√φ ref. The defaults are the mean criterion and the SRD/√φ criterion, with SRD/√φ ref 3.3. Each criterion looks at its own smaller-|value| flank (ties go left).
- An **all / any** toggle: with *all* every checked criterion must pass; with *any* one is enough.
- Two switches, both off by default. **Merge would-merge segments whose flank signs disagree** (the SRD/√φ values on the two flanks have opposite signs). **Merge segments whose criteria pick different flanks, toward the one furthest below its threshold** (the passing criterion with the smallest value-to-threshold ratio wins; ties go left).

A segment that passes the test and is not held back is a *would merge* segment and is drawn with a × hatch in Rule 0's colour. The other cases:

- *Would merge, signs disagree* (switch on): × plus an amber vertical | overlay.
- *Would merge, flanks disagree* (switch on): × plus a violet horizontal ― overlay. A segment with both disagreements draws both overlays.
- *Held back*: the segment passes the test but a disagreement applies while its switch is off. It is not merged and is drawn with the amber | and/or violet ― marks and a ⚠ on its marker.
- *Signs disagree only*: the segment fails the test but its flank SRD/√φ signs differ; drawn with an amber | and never merged. The "signs-disagree hatch" layer checkbox toggles only this case; would-merge and held-back marks are always drawn.
- *Blocked by gap*: a would-merge segment whose chosen flank sits across a genomic gap is still hatched and counted as would merge, drawn with a dotted outline, and counted as "blocked by gap" in the legend. It merges only when **allow merging across a genomic gap** is ticked. On chr8 at default settings every would-merge segment is gap-adjacent, so enable gap crossing to see merges there.

The legend has one entry per case present, with counts. Hover over a segment's marker for the per-criterion detail: each checked criterion's flank, value, threshold, fraction and pass/fail, then the combined result, the sign check, the flank-conflict check, the direction, the outcome and which user rules also match. These hatches and settings only change what is drawn until you run a merge. The earlier "SRD-screen proposal" layer, the flank dropdown and the "effect metric for flank choice" dropdown were removed; the criteria checklist replaces them.

**2. Merge rules**, below the tracks, starts with a pinned **Rule 0 card**. It cannot be edited, deleted or moved; it summarises the section 1 settings on one line and shows how many segments it flags (would-merge segments across the displayed rows). Below it is an ordered list of rules you build yourself, which add to Rule 0. Add a blank rule or a preset (weak flank, flank signs disagree, ambiguous flank preference, small segment). Each rule has:

- **Conditions, all of which must hold** for an internal segment: *weaker flank below threshold* (min(|L|, |R|) < threshold, for SRD, SRD/√φ or Δlog2FC with an estimator), *flank signs disagree*, *SRD and Δlog2FC disagree on the weaker flank*, and *segment size at most N bins*.
- **Merge into**: the metric that decides which neighbour a flagged segment joins (the flank with the smaller |score|; ties go left).
- **A hatch pattern and colour.** Each enabled rule hatches exactly the segments it flags, live, and its card shows how many segments it flags on the displayed rows. Hover over a hatched segment's diamond marker for each condition's values and the merge outcome.

A segment is flagged when Rule 0 or any enabled rule matches. **Run merge** then merges one flagged segment per step: Rule 0 goes first, then the user rules in list order, then the smallest score, then the leftmost segment, recomputing every score, including the pooled dispersion φ and each criterion's flank choice, and re-evaluating every rule after each merge (so segments not hatched at the start can become eligible later). Run merge works with Rule 0 alone; the user rule list may be empty. It needs at least one criterion, positive τ and SRD/√φ ref, an integer "small ≤" of at least 1, and every enabled user rule complete. The φ shown in the merge history is the pooled dispersion of the partition the merge decision was made on, i.e. before that step's merge. Genomic gaps are never crossed unless you allow it; a segment whose chosen neighbour sits across a gap is skipped, not redirected. Use the step slider (or Reset/Undo) to inspect intermediate segmentations and download the run (JSON) or the merge history (CSV), which names the rule behind each merge (Rule 0 appears as `rule0`). Rules and the section 1 settings are saved in sessions. Sessions saved before this change are converted into preset rules (migrated merge-proposal rules do not flag a segment whose other flank score is missing, e.g. Δlog2FC with a zero level, unlike the old layer). Sessions that stored an archived colour choice are mapped to the live metric: `delta_median` to Δlog2FC median, `delta_mean` and `log2fc` (SRD-screen) to Δlog2FC mean, `delta_trim` to Δlog2FC IQR-mean, `srd` to |SRD|, `srd_phi` to |SRD/√φ|. The old flank and effect-metric settings and the old proposal layer flag are dropped, and the old flank-consistency label flag becomes the signs-disagree layer. Rule 0 is an exploratory diagnostic built from the displayed settings, not a validated decision rule. It is independent of the original bootstrap/chain-guard merge prototype, and CN calls are not recomputed on merged boundaries.

**Candidate audit.** For a candidate that was selected for significance but not accepted and has no recorded rejection reason, hovering shows a reason derived from the archived p-values: "local permutation p = … ≥ 0.001", otherwise "global permutation p = … ≥ 0.001", otherwise "not accepted (no recorded reason)". Each derived reason ends with "(cutoff inferred: all accepted candidates have p < 0.001)"; the cutoff is inferred, not recorded.

## Recompute integration flank comparisons

In **Integration · segmentation & CN**, use **Generate alternative flank view**. Choose chromosome/source and maximum internal-segment size, then compare arithmetic, median and two-sided IQR effects against both adjacent segments. Click a point or choose a segment to see its unsmoothed bin values and summaries.

**Compute paired cell bootstrap** provides a separate arithmetic-mean comparison using 200 paired whole-cell replicates (seed 42). It does not rerun clustering, segmentation, SRD or merging. Download the JSON/CSV to retain these generated results; they are not yet stored in normal app sessions.

The original generators and configuration were found upstream and successfully rerun: all 11 diagnostic tables match their archived versions byte-for-byte. The app now loads these verified reference tables automatically, restoring the original SRD/noise diagnostics, merge summary and candidate audit. They are stored separately under [reference/integration-prototype](../reference/integration-prototype/README.md); the deleted story prototype directory stays absent. The original median merge bootstrap uses seed 20260907 and remains distinct from the alternative arithmetic bootstrap above. See [the integration review](../validation/integration-review/REVIEW.md) for details and scientific caveats.

## Integration chart loading

The segmentation tab renders its chromosome and diagnostic charts as you scroll them into view. This avoids initializing all eight plots during the initial CN-panel load. Once displayed, plots retain their zoom while scrolling. See [responsiveness checks](../validation/integration-freeze/REPORT.md) for the reported freeze investigation and its current limits.
