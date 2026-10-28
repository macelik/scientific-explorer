# History

Artifacts from earlier AI-assisted development sessions on this repository,
kept for provenance rather than as runtime dependencies:

- **`sonnet-changes.md`** — an explicit before/after changelog written during
  an early debugging session on the Integration segmentation tab, documenting
  exactly what was changed, why, and what was deliberately left alone. Written
  because this repository's contributor norms (see [AGENTS.md](../../AGENTS.md))
  require stating old behavior, new behavior, reason, and validation performed
  for any change to scientific/segmentation code — this file is the working
  example of that discipline in practice.
- **`2026-10-28-merge-rounds-allelic-segmentation.md`** — record of the
  session that added the genome merge / CN re-call / merge rounds, the
  allelic segmentation tab, and switched the depth cluster breakpoints to
  pyEpi's original genome-wide global permutation test (old behaviour, new
  behaviour, reason and validation for each scientific change).
- **`sonnet-backup-original/`** — byte-exact copies of source files as they
  stood immediately before that debugging session's edits, kept as a manual
  revert path from before this repository had its own git history to serve
  that purpose. Superseded by git history for any file touched since; kept
  here as a record of the session, not as an active revert mechanism.

None of these files is imported or referenced by any runtime code.
