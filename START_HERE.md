# Scientific explorer

Open **http://127.0.0.1:8766**. This is the separate version; Fable's `interactive-explorer` is untouched.

If the service is stopped, run from this directory:

```bash
python3 run.py
```

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

## Reading the result

- The original/edited AD curves show whether the winning candidate moves. Both child curves and exact split tests are expandable.
- **Event-edge diagnostics** show whether each intended boundary was accepted, its rank among candidates in the tested node, and both permutation p-values.
- Four CN rows distinguish original calls, edited values with production boundaries, automatic resegmentation, and an explicitly diagnostic segmentation with event edges imposed.
- Expand **Holmes, Watson and combined calls** to see the number of event bins in each state under each boundary treatment.
- Run history retains the most recent 12 results per cell/chromosome. Inspect earlier runs or export the complete retained history. A changed draft marks the displayed result stale until rerun.

Save/download a normal session to preserve drafts, retained experiments, comparison windows, and manual trees. Reload and upload that JSON to restore them. Export a single experiment as JSON for provenance and arrays, or export regional values as CSV. Sessions are separate from Fable's sessions.

## Interpretation

These are interventions on **GC-corrected, normalized X**, not simulated raw fragments or a rerun of preprocessing. Empirical class sampling is conditioned on existing calls, retains the donor value distribution in expectation, and discards spatial correlation. It is not a validated biological CN simulator. Pattern duplication introduces repeated structure deliberately.

The copied production functions have no hard 50/100-bin size cutoff. Candidate competition, the two significance tests, recursion depth, and downstream fitting determine whether a particular event is resolved. Width-series results describe these specific realizations, not calibrated detection probabilities or a universal size threshold.

Automatic experiments offer **Recursion limit k** from 1 to 10 (default 2). A run tests depths 0 through k−1, with 1,000 permutations, seed 42 for tests, and strict `p < 0.001` for both tests. Increasing k continues only through accepted parents. The immediate children of a rejected root remain explicitly hypothetical and never expand further. The sampling seed is separate from the test seed. Edited global tests use the edited full-genome X row. Only the selected chromosome is automatically resegmented; other chromosomes retain their production boundaries. All CN rows are fitted genome-wide. Inside **Both children and exact split tests**, use **Inspection depth** to view deeper tested nodes. Orange shading marks the modified target clipped to the displayed child intervals; for extend/duplicate it highlights the overwritten bins, not the untouched donor region. Solid vertical lines show accepted breakpoints (including accepted descendants); dashed lines show unaccepted candidates. Blue represents original X and brown edited X. Changing k marks the last result stale; rerun to update it. Depth is retained in sessions and exports. Manual segmentation remains available for user-directed boundaries on original X.

## Compare arithmetic and IQR segment means

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
choose the metric that colours small segments and which flank it is read from, switch between
binned bands and a continuous scale, and add the label layers you want to talk about (proposal,
x/o// consistency, reason, values, CN strip, candidate audit, allelic). Click a shaded segment or a
point in the diagnostics to open the segment's full flank table. τ and the SRD reference are
inspection settings, not thresholds.

## Recompute integration flank comparisons

In **Integration · segmentation & CN**, use **Generate alternative flank view**. Choose chromosome/source and maximum internal-segment size, then compare arithmetic, median and two-sided IQR effects against both adjacent segments. Click a point or choose a segment to see its unsmoothed bin values and summaries.

**Compute paired cell bootstrap** provides a separate arithmetic-mean comparison using 200 paired whole-cell replicates (seed 42). It does not rerun clustering, segmentation, SRD or merging. Download the JSON/CSV to retain these generated results; they are not yet stored in normal app sessions.

The original generators and configuration were found upstream and successfully rerun: all 11 diagnostic tables match their archived versions byte-for-byte. The app now loads these verified reference tables automatically, restoring the original SRD/noise diagnostics, merge summary and candidate audit. They are stored separately under [reference/integration-prototype](reference/integration-prototype/README.md); the deleted story prototype directory stays absent. The original median merge bootstrap uses seed 20260907 and remains distinct from the alternative arithmetic bootstrap above. See [the integration review](validation/integration-review/REVIEW.md) for details and scientific caveats.

## Integration chart loading

The segmentation tab renders its chromosome and diagnostic charts as you scroll them into view. This avoids initializing all eight plots during the initial CN-panel load. Once displayed, plots retain their zoom while scrolling. See [responsiveness checks](validation/integration-freeze/REPORT.md) for the reported freeze investigation and its current limits.
