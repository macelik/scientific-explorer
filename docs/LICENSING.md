# Licensing and provenance

This document records a licensing/provenance audit of this repository, performed
before adding a `LICENSE` file, so a license is never applied to code whose
origin hasn't been checked.

## What this repository contains

- **`server/`, `web/`, `tests/`, `docs/`** — application code, tests, and
  documentation authored for this project (with AI assistance; see the root
  README's "Research-led, AI-assisted engineering" section).
- **`reference/integration-prototype/`** — copied snapshots of prior research
  scripts, verified byte-for-byte reproducible against their upstream origin
  (see [`reference/integration-prototype/README.md`](../reference/integration-prototype/README.md)).
  No copyright/license header exists in the copied files; they are the
  repository owner's own prior research output, not a third-party import.
- **`web/`'s original scaffold** — adapted from a sibling internal prototype,
  `interactive-explorer`, itself unlicensed and undocumented as third-party
  open-source.

## What this repository does NOT contain

- **`pyEpiAneufinder`'s source.** `server/pyepi_adapter.py` loads it from an
  external, configurable path (`PYEPI_SCIENTIFIC_PACKAGE`) at runtime; its
  functions are called, never copied into this repository's source tree (see
  `server/science.py`'s module docstring: "Every statistic here is a direct
  call into the supplied pyEpiAneufinder functions... nothing is
  re-implemented"). This repository does not redistribute that package, so
  licensing this repository's own code does not purport to relicense it.

## Audit conclusion

No `LICENSE`, `COPYING`, or license metadata field was found anywhere in this
workspace — not in this repository, not in `interactive-explorer/`, not in
`pyEpiAneufinder/`. Nothing found indicates a third-party open-source license
applies to any code embedded in this repository.

## Decision

Confirmed by the repository owner (2026-09-30): they are the sole author/
rights-holder of the code under `server/`, `web/`, `tests/`, `docs/`, and
`reference/integration-prototype/`. An [MIT license](../LICENSE) now covers
this repository's own code. `pyEpiAneufinder` remains an external runtime
dependency with its own, separately-tracked license — this repository does
not relicense it. `reference/integration-prototype/` is covered by the same
MIT license as archival research material, since it is prior work by the
same research effort, not a third-party import.
