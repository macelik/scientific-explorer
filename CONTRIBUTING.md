# Contributing

This is scientific software. Correctness and provenance matter more here
than typical UI polish — a change that "looks right" but silently alters a
statistic, a significance test, or a segmentation boundary is a bigger
problem than a rough edge in the interface.

## Ways to contribute

- Bug reports (UI or operational)
- Algorithmic or scientific issues (a result you believe is wrong or
  mis-described)
- New diagnostics, experiments, or scoring methods (see
  [docs/EXTENDING.md](docs/EXTENDING.md))
- Code changes
- Tests and validation evidence
- Documentation

## Reporting an algorithmic or scientific issue

This is different from a UI bug report. Please include:

- What you expected and what you actually observed
- The exact settings used (metric, estimator, threshold, seed, chromosome/
  source, recursion depth — whatever applies)
- Whether it reproduces with a fresh session (rules out stale/cached state)

## Code changes

Read [AGENTS.md](AGENTS.md) first. Run the test suite before opening a PR:

```bash
python3 -m pytest -q tests/test_hatch_merge.py tests/test_experiments.py tests/test_segment_estimators.py tests/test_integration_review.py
cd web && npm run build
```

TypeScript helper logic that ports a Python classifier (like
`tests/hatch_merge_helpers.ts`) is bundled with esbuild and run under Node —
see [tests/README.md](tests/README.md) for the exact command.

## PR checklist

State which category your change affects:

- [ ] Visualization only
- [ ] Data transformation
- [ ] Statistical score
- [ ] Significance testing
- [ ] Segmentation
- [ ] CN assignment

If anything beyond "visualization only" is checked, state in the PR
description: the **old behaviour**, the **new behaviour**, the **reason**,
and **what validation was performed** (which tests, which script, and
whether results were compared against a prior run or archived reference).
