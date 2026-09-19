# Experiment recursion and child annotations

Implemented 2026-09-19 in scientific-explorer only.

The experiment editor now exposes recursion k=1–10, default 2. Both scenarios
use the selected limit, as do width-series runs. The backend uses the same
argmax and strict local/global acceptance gates as the supplied recursive caller.
Depth 0 is the chromosome root; k allows testing depths 0 through k−1. Existing
immediate hypothetical children of a rejected root are retained for inspection,
but never enter the segmentation or expand further. Original data remain immutable.

The root and child plots show solid accepted-breakpoint lines by scenario color
and dashed unaccepted candidate lines. Child shading intersects the edit target
with the intervals of the displayed nodes. In extension/duplication experiments,
this target is the overwritten interval rather than the untouched source plus
extension. Displayed scenario intervals may differ; their clipped overlaps are
merged before shading so overlapping regions do not become artificially darker.
All positions use chromosome-local retained bins mapped to actual genomic starts
and ends. The inspection-depth selector displays deeper reached nodes.

Changing k makes the previous result stale. Effective k is shown in run history,
the result and provenance, and included in the cache key, exports and sessions.
Old saved drafts/results without a request k continue to mean k=2.

## Verification

- 44 Python checks passed, including comparison with the original
  `recursive_getbp` at k=1, 3 and 5, default k=2 parity, pruned-branch handling,
  and rejection of invalid/noninteger depths.
- TypeScript build passed. The overlay helper checks cover clipped targets,
  genomic boundary coordinates, exclusion of node-edge/outside breakpoints,
  accepted/candidate distinction, old-session depth inference and stale drafts.
- Browser check ran k=4 on TTAGGCTAGGCCGGAA-1 / chr9, checked every displayed
  child breakpoint against the returned accepted boundaries, confirmed deeper
  nodes, selected depth 2, changed k to mark results stale and restored the
  draft/result from a session. No page errors. See `browser-results.json` and
  the depth-1/depth-2 screenshots, both visually inspected.

## Reload investigation deferred by the user

The user reports Chrome becoming unresponsive on the initial cohort view, with
a second tab working while the first remains open. Twelve cold/warm loads and
reloads across four test-browser tabs did not reproduce it (0.89–1.62 seconds to
load and navigate). See `baseline-reload-probe.json`. The user asked to move on
to the experiment changes. No freeze fix is claimed, and no speculative change
to startup or Chrome configuration was made.
