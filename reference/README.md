# reference/

Archival, read-only material used to verify reproducibility — never a
runtime dependency. See [`../docs/EXTENDING.md`](../docs/EXTENDING.md#when-something-belongs-in-reference-rather-than-runtime-code)
and [`../docs/LICENSING.md`](../docs/LICENSING.md).

## `integration-prototype/`

- `src/` — prototype source snapshots (e.g. the SRD/dispersion formulas that
  `server/hatch_merge.py` ports and cites), kept for exact reproducibility
  auditing. New code must not import from here.
- `archived/` — the historical, original tables from the upstream prototype
  run — the fixed reference point everything else is checked against.
- `regenerated/` / `screening-regenerated/` — freshly computed tables,
  checked byte-for-byte against `archived/`.
- `analysis/` — supporting analysis artifacts from the reproduction effort.

See [`integration-prototype/README.md`](integration-prototype/README.md) for
the detailed reproduction instructions and byte-for-byte verification
results.
