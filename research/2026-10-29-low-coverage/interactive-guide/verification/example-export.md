# Numbat-multiome

A joint probabilistic model combines reference-normalized ATAC depth and phased SNP alleles along chromosomes.

## How it works
1. Prepare: Supply raw bin-by-cell ATAC counts, genomic bins, normalized ATAC reference and phased SNP-level allele counts.
2. Model: Use a Poisson-lognormal coverage model and a beta-binomial phased allele-count model to score candidate joint states.
3. Connect: A hidden Markov model (HMM) shares information across neighboring bins and models state transitions.
4. Pool: Iterative clone-aware pooling improves estimates and returns state probabilities and clones.

## Illustrative example — not a tool run
Joint evidence

Imagine a region with normalized depth near 0.6 and phased counts 18 versus 6.
- Compare depth against a defensible ATAC reference.
- Score allele imbalance using the phased counts and their uncertainty.
- Combine both signals with neighboring bins in the HMM.

A loss-like state may gain support, but this hypothetical input is not a tool result or proof of absolute copy number.

## Proposed adaptation
Stage: Choose a baseline

Available:
- Existing ATAC depth and site observations
- Published ATAC-bin mode is a relevant comparator

Needed:
- Suitable ATAC reference
- Compatible raw fragment counts, bin geometry, phased SNP-level counts and verified phase continuity

First experiment: Audit inputs and reference suitability before a small ATAC-bin feasibility run.

Success: Demonstrate input compatibility and report posterior uncertainty against independent evidence.

## Limitations
- No suitable normal reference is established here; MCF10-A is not assumed diploid.
- Current pyEpi pipeline uses AD segmentation, not this HMM.
- Missing alleles and uncertain phase restrict power; ATAC 0:0 sensitivity is unvalidated here.

## Sources
- [Numbat-multiome paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12531993/)
- [ATAC-bin input vignette](https://kharchenkolab.github.io/numbat/articles/numbat-multiome.html)

Scientific Explorer research guide · 2026-11-02. MCF10-A is development data, not CN ground truth.
