# Interactive CNV research guide

Open **[index.html](index.html)** in a browser. It is a single, self-contained
file: no server, package installation, fonts, or network connection are needed.
Opening a cited paper follows its external link. You can copy `index.html` alone
to another machine.

The guide explains all **36 entries** from the six research categories. These
are tools, model families, reference strategies and technical annotations, not
36 distinct CNV callers. Numbat, for example, appears from more than one
research perspective.

## Using the guide

- **Methods:** search or filter the library, click each inference step, and
  switch between the mechanism, a worked example, and adaptation to this
  repository. Each entry lists inputs, outputs, missing prerequisites, a first
  experiment, limitations and primary sources.
- **Interactive examples:** vary depth/accessibility, a reference level, allele
  observations, or count-model parameters. Equations and interpretation update
  with the plots. These are deterministic teaching calculations, not fitted
  models or executions of the named tools.
- **Our data:** inspect the saved 300-cell pilot across six candidate windows,
  three coverage terciles and five WNN clusters. Observed summaries remain
  separate from the teaching examples.
- **Workflow:** locate each entry in the proposed development sequence.
- **Compare:** add up to three entries to compare requirements and assumptions.
- **Save notes / Print view:** export the current method to Markdown or print
  the selected view. The URL fragment retains the selected method and panel,
  or the selected example; slider values and comparisons reset on reload.

## Scientific scope

MCF10-A is a development sample, not a normal/diploid truth label. The intended
algorithm must generalize across samples. Published behavior, proposed
adaptations, hypothetical examples, and observed pilot summaries are labelled
separately. The guide does not run a caller, estimate deletion prevalence,
change the application's scientific calculations, or read individual-cell
matrices.

The depth illustration shows a multiplicative ambiguity; its fixed background
curve is not sampled data. The allele illustration uses a uniform-prior
binomial/Beta posterior, **not** a beta-binomial emission or a calibrated LOH
test. Its numerical credible intervals use a 4,000-point midpoint grid. The
count illustration uses Poisson, NB2 (variance `mu + mu^2/r`), and zero-inflated
NB distributions, with a final bar collecting counts of 16 or more. The
extra-zero parameter is an assumption, not a deletion fraction. Distribution
parameterizations are documented by SciPy for
[negative binomial](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html),
[Poisson](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.poisson.html),
and [Beta](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.beta.html).

Pilot values are embedded directly from
[`zero_support_by_stratum.tsv`](../e-per-cell-evidence/artifacts/zero_support_by_stratum.tsv).
The existing panel membership is fixed; no new random sampling was performed.
AT counts selected SNP-site observations, not established independent
molecules. Historical seed and phase continuity remain unresolved. Whole-window
zeros, per-bin zeros and unavailable evidence have different meanings.

## Sources and rebuilding

`content/callers.json`, `content/models.json` and `content/context.json` preserve
the three GPT 6.1 Terra agents' structured explanations. The build merges them
with the original six `methods.json` files, retaining their detailed research
records. Existing research records and pilot tables are not rewritten.

From the repository root:

```bash
python3 research/2026-10-29-low-coverage/interactive-guide/build_guide.py
```

Edit `shell.html`, `guide.css`, `guide.js` or the structured explanation files,
then rebuild. `build_manifest.json` records the category counts and SHA-256
hashes of the source files, pilot table and generated HTML. The HTML has no
runtime dependencies.

## Verification

See [`verification/README.md`](verification/README.md) for the actual checks,
reproduction commands and limitations. Screenshots and browser reports live in
that new directory; archived repository validation records are unchanged.
