# Pilot

No numerical pilot was run: the coordinator-owned frozen panel and manifest were not yet available during the initial audit; they are now available, but task E covers the count pilot, and the alignment build and count units remain unresolved. Running overlaps or ratios now would imply precision not supported by provenance. The read-only audit inspected the local VCF header, builder, pyEpi source and official method documentation. No cells were sampled, no seed was used, no source data modified, and no computation results are claimed. Reproducible inspection commands:

```bash
 gzip -cd /home/mami/Desktop/EGA/SNP_project/alpha/haplo_initial/updated_phasing/phasing/phased.vcf.gz | head -100
 sed -n '1,430p' /home/mami/Desktop/EGA/haplotype_aware_cnv/Build_Aligned_HA_HB.py
 rg -n -i 'blacklist|threshold_blacklist_bins|GC correction|process_fragments|process_count_matrix' /home/mami/Desktop/EGA/pyEpiAneufinder/pyEpiAneufinder/pyEpiAneufinder_v1.py
```

The gzip broken-pipe message from `head` is expected and does not imply a completed full-file integrity check. Paths are local provenance, not data exports.

## Provenance follow-up, 2026-10-29

Inspected script SHA-256 `b3892252f9198f14539d2a16be19b9a8fa8e346c2dc348731366f5cbb56475df`. Builder source line 60 performs sparse haplotype assignment; 170–205 reads GT and filters phased heterozygotes without preserving PS; 221–240 converts coordinates; 267–291 sets input and output defaults; 319–321 joins/selects genotypes; 370–392 aggregates and writes results. The planned output names match existing updated-phasing auxiliary files. No output-generation log or historical script hash was found, so the observed source-to-output relationship is conditional.

Relevant combined VCF header excerpt (first chromosome Eagle command; not proof of alignment reference):

```text
##contig=<ID=chr1,length=248956422>
##eagleVersion=2.4.1+htslib-1.9
##eagleCommand=eagle --numThreads 24 --vcfTarget pileup_outs/phasing/sctac_chr1.vcf.gz --vcfRef /work/project/ladcol_011/haplo/1000G_hg38//chr1.genotypes.bcf --geneticMapFile=/work/project/ladcol_011/haplo/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz --outPrefix pileup_outs/phasing/sctac_chr1.phased
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sctac
```

Header columns shown with spaces for readability; original is tab-separated. This command describes the chr1 phasing header, not verified per-chromosome commands. Follow-up used `sha256sum /home/mami/Desktop/EGA/haplotype_aware_cnv/Build_Aligned_HA_HB.py` and `nl -ba` on that source; no script executed on sample data.
