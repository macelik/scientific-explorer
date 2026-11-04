# Wave 1 findings for tasks B, C and D

Read the saved A, E and F findings for sources and limitations. These are
first-pass research results and exploratory diagnostics, not validated CN calls.

1. **Numbat-multiome is now directly relevant to scATAC.** Do not repeat the
   initial characterization of all Numbat ideas as RNA-only transfers. The
   [2025 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12531993/) and
   [official multiome vignette](https://kharchenkolab.github.io/numbat/articles/numbat-multiome.html)
   describe ATAC-bin analysis and comparison with WGS. Task A reports raw bin
   counts, phased SNP-level observations and a normalized ATAC reference as
   required inputs. The current no-reference MCF10-A setting is still a gap.
   Verify exact assumptions against these primary sources for your task.
2. **AT is not established as unique-fragment coverage.** The inspected builder
   sums selected heterozygous SNP-site DP. The same molecule may contribute at
   multiple SNPs; original counting settings remain unknown. Phase sets are
   not retained in the exported matrices. The builder's correspondence to the
   updated outputs is supported but its historical revision/run is unproven.
3. **The phasing setup supports hg38; alignment provenance is still open.**
   Do not perform assembly-dependent annotation joins as if exact alignment
   reference identity were established.
4. **The shared panel is ready.** `manifest.json`, `cells.tsv`, `regions.tsv`,
   `region_bins.tsv` and `region_counts.tsv` provide deterministic exports with
   input and output hashes. Every selected barcode matched, all bin coordinates
   matched, and AT=HA+HB was checked for all 300 cells. `PREPARATION.md` contains
   the complete preparation source. No random draw or refiltering was performed.
5. **Sparse bins differ from empty regions.** E's pilot finds many raw-zero
   bins but far fewer completely empty target windows. Whole-window allelic
   totals are very sparse. There was no null-model fitting, so no claim about
   excess structural zeros, loss probability or allelic power follows. Consult
   E's tables for exact values rather than treating these descriptive statements
   as calibrated tests.

Keep outputs in your own agent directory. Do not duplicate the count pilot or
re-run the full pipeline. For most wave-2 tasks, a source/input-compatibility
assessment is more useful than a speculative computation.
