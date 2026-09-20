# Segmentation-tab responsiveness

2026-09-20. The user reports that a fresh Integration · segmentation page shows
the cluster and per-cell CN panels, then becomes unresponsive before the
chromosome tracks load. No session restore is involved.

## Investigation

The persistent hang was not reproduced using either the test Chromium or the
installed Google Chrome, including an isolated visible Chrome window. Profiling
did confirm unnecessary synchronous rendering: all eight Plotly charts were
initialized on entry, including off-screen diagnostics. The track ready callback
updated parent React state, causing another parent render; its newly created
config object also retriggered the Plotly effect despite unchanged track data.

The baseline profile and viewport/heartbeat results are preserved as
`baseline-profile.json` and `baseline-browser.json`. These show active chart
rendering, not evidence of an infinite loop. No Chrome configuration, data,
scientific computation or original interactive-explorer source was changed.

## Changes

- The eight segmentation plots now reserve their layout space and initialize
  when they enter the viewport. A plot stays mounted after first display, so
  scrolling preserves zoom and selection.
- The chromosome plot's ready callback stores its element in a ref instead of
  rerendering the parent. Its configuration is memoized.
- Numerical tracks, boundaries, small-segment shading/hatches, diagnostics and
  export actions retain their existing definitions.

## Verification

`browser_integration_rendering.cjs` first failed against the previous build
because the off-screen chromosome chart was already mounted. It passes with
the change: zero off-screen track charts on initial load; chromosome tracks,
five allelic overlays, all seven diagnostic plots, controls and tab navigation
work with no page errors. See `rendering-results.json`.

The visible Chrome profile also passed loading, navigation, and resizing across
800–1920 pixels; only the one viewed chromosome chart was initialized during
that sequence, compared with eight in the baseline. The production TypeScript/
Vite build passed. All 66 recorded original-explorer source hashes still match.

The user's persistent hang still requires confirmation in their browser. This
change addresses verified rendering overhead; it does not establish that the
unreproduced hang had the same cause.
