# Note (to implement later): calibrating the Δlog2FC merge threshold τ

**Status:** idea, not implemented. Discussed 2026-10-27.

## Why

Segments merge when |Δlog2FC| between a segment and its flank is < τ. τ should keep the smallest real copy-number step and merge anything below it. In scATAC data, steps are smaller than their clonal values. The causes are tumour fraction, background accessibility, smoothing and compression. A global normalisation factor largely cancels in an adjacent-segment ratio.

## Model

Lump all attenuation into one effective factor `a` (0 to 1):

    observed ratio = 1 + a · (CN/2 − 1)

At full clonality, the adjacent steps around a diploid baseline are:

| Adjacent pair | Clonal step (log2) |
|---|---|
| 2 → 1 (hemizygous loss) | −1.00 |
| 2 → 3 (single-copy gain) | +0.585 |
| 3 → 4 | +0.415 |
| 4 → 5 | +0.32 |

The smallest step to keep is the single-copy **gain**, not the deletion.

### How `a` shrinks the steps

| a | gain step | loss step |
|---|---|---|
| 1.0 | 0.585 | −1.00 |
| 0.8 | 0.485 | −0.74 |
| 0.6 | 0.379 | −0.51 |
| 0.5 | 0.322 | −0.42 |
| 0.4 | 0.263 | −0.32 |
| 0.3 | 0.202 | −0.23 |

- τ = 0.5 keeps a single-copy-gain boundary only if a ≳ 0.83, and a hemizygous-loss boundary only if a ≳ 0.59. For most scATAC data that is too strict.
- Gains merge before losses at any τ.
- A fixed fallback of τ ≈ 0.25–0.3 keeps gains down to a ≈ 0.4–0.5.

## Calibration from a known CNV (per cluster)

1. **Baseline.** Take a robust genome-wide centre of the cluster's signal, such as the interquartile mean, with the known CNV region left out.
2. **Observed step.** r_obs = log2(mean of the CNV region / baseline).
3. **Solve for `a`.** a = (2^r_obs − 1) / (c/2 − 1), where c is the region's known copy number:
   - loss (c = 1): a = 2·(1 − 2^r_obs)
   - gain (c = 3): a = 2·(2^r_obs − 1)
4. **Set τ.** τ ≈ f · log2(1 + a/2), with f ≈ 0.5. Check that τ is above the |Δ| noise in copy-neutral regions.

### Caveats

- **The baseline must be the neutral state.** The IQR mean is only neutral when most of the genome is diploid. For aneuploid clusters, prefer the most common segment-level value or known-neutral regions.
- **The known CNV must be clonal within the cluster, with the correct c.** Otherwise `a` also absorbs subclonality.
- **Measure on the same signal τ sees.** Use the same normalisation, smoothing and segment means of log2FC. Leave out the CNV's edge bins.
- **Cross-check at the boundaries.** Also measure the step directly at the CNV's boundaries, which is what τ actually sees. If the two disagree, look for local effects such as GC or edge smoothing.
- **`a` varies by cluster.** Calibrate τ per cluster, or use the smallest `a` across clusters.
- **Check gain and loss separately when both are known.** If they disagree, prefer the gain-based `a`.
- **Bootstrap.** Bootstrap the bins for a confidence interval on `a`, and therefore on τ.

**Alternative without a known event:** plot |Δ| at all accepted breakpoints and place τ at the trough between the noise mode and the real-step mode.
