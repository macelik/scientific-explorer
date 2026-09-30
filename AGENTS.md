# AGENTS.md

Guidance for AI agents and automated contributors working in this
repository — read this before making changes.

- Read [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and
  [docs/EXTENDING.md](docs/EXTENDING.md) before making a structural change;
  read [docs/SCIENTIFIC_WORKFLOW.md](docs/SCIENTIFIC_WORKFLOW.md) before
  making a scientific one.
- Do not change scientific/statistical behaviour just to make a UI or test
  pass. If a test fails because the underlying statistic changed, that is a
  scientific decision requiring the format below, not a quick fix.
- Do not silently substitute a formula, statistic, or significance test for
  a "simpler" or "more standard" one. Reuse existing functions (see
  docs/EXTENDING.md) rather than reimplementing.
- Preserve provenance and random seeds. Every stochastic result must remain
  reproducible from its recorded seed and parameters.
- Distinguish exploratory diagnostics from production/validated behaviour in
  both code comments and UI copy — see [docs/VALIDATION.md](docs/VALIDATION.md)
  for the existing convention (`consistent with`, `exploratory alternative`,
  never `proves` or `biologically validated`).
- Prefer small, inspectable changes over large rewrites. A prior
  whole-branch review in this repository (see
  [validation/hatch-merge/HANDOFF.md](validation/hatch-merge/HANDOFF.md))
  caught six real issues in one large feature; smaller changes are easier to
  verify.
- Run the relevant Python tests
  (`python3 -m pytest -q tests/test_hatch_merge.py tests/test_experiments.py tests/test_segment_estimators.py tests/test_integration_review.py`)
  and, for frontend changes, `npm run build` in `web/`, before considering a
  change done.
- Update the relevant `docs/*.md` file when adding a scientific feature — an
  undocumented capability is, for this repository's purposes, not finished.
- Never describe an exploratory transformation (simulation, empirical
  resampling, an alternative estimator) as biological ground truth.
- Do not rewrite or regenerate archived files under `reference/` or
  `validation/` — they are provenance records. If a check needs updating,
  add a new report; do not silently edit an old one.

If changing a scientific calculation, state explicitly in the PR/commit
message: the old behaviour, the new behaviour, the reason, and what
validation was performed. See [CONTRIBUTING.md](CONTRIBUTING.md) for the
full PR checklist.
