# Integration review and alternative view

Scope: preserve story artifacts and original interactive-explorer. Audit the existing scientific-explorer integration; add an on-demand alternative without replacing the source figures or making method changes.

Findings driving the work: prototype metric/bootstrap/manifest tables are absent; original figure script paths are not portable; allelic nonzero fraction is mislabeled coverage; malformed Plotly allelic axes and raw-depth shape coordinates can distort tracks; integration fingerprints omit depth and several exported inputs; manifest loading is conditional and overwritten.

Visual direction: keep the existing quiet scientific workspace. Put one Generate alternative flank view action near the segmentation introduction. Show a paired left/right effect plot, selected segment with both neighboring regions, bin-value distributions and exact summary table. Use consistent left/segment/right colors and restrained selection highlighting. No ornamental animation.

Compute from available mean pseudobulk and fixed accepted boundaries, by selected source/chromosome and maximum internal-segment width. Compare arithmetic mean, median, and two-sided 1.5-IQR mean. Keep zero-level ratios undefined, expose original/retained lengths and per-bin spread, and make no SRD/bootstrap/merge claims. Export generated values and provenance. Frozen results become stale when scope changes.

Validation: synthetic numerical/coordinate tests; source identity/manifest tests; real artifact checks; TypeScript build; browser generation, selection, exports, stale handling and corrected tracks. Record missing evidence separately from confirmed code defects.
