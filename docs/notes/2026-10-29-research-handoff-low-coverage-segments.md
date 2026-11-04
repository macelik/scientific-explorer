# Research handoff: low-coverage depth segments, allelic blind spots, and no normal reference

**Written:** 2026-10-29. **For:** a research agent. **From:** the scientific-explorer development session (user + Claude).

**Your job:**
1. Find methods we can adopt to tell real copy losses apart from low-coverage artefacts in scATAC CNV calling.
2. Find ways to infer normal (non-tumour) cells when none are labelled.

Look beyond CNV callers. Adjacent fields (image processing, segmentation, change-point detection, machine learning, single-cell foundation models) are explicitly in scope. Report methods with evidence and fit; don't write code in the repo.

---

## 1. Project goal

We are building a **haplotype-aware CNV analysis tool for single-cell ATAC data**. It combines two signals measured from the same fragments:

- **Read depth.** We use pyEpiAneufinder (pyEpi): 100 kb bins, GC correction and normalisation, with negative values clipped to 0. Breakpoints come from a recursive Anderson–Darling (AD) scan with local and global permutation tests. CN calls come from `assign_gainloss_new`.
- **Allelic signal.** HA/HB counts are reads on the two phased haplotypes at heterozygous SNPs, aggregated per 100 kb bin on the same grid. AT = HA + HB. The derived signals are dBAF = (HB − HA)/(2·AT), |dBAF|, shrunk |dBAF| = |HB − HA|/(2·(AT + 2)), and dBAF oriented to the all-cells direction.

**The long-term aim** is CN calls where allelic evidence confirms or rejects what depth suggests: one-copy loss, LOH, gain, allele-specific events. We currently segment **cluster pseudobulks** (5 WNN clusters plus all cells). Per-cell calls reuse the cluster segmentation.

## 2. The issue

The depth segmentation produces **short segments called as losses** in every cluster, for example on chr7. In those segments the depth really does drop, but the allelic track is almost empty there.

- **Why the allelic track is empty:** AT drops even more than depth. Most bins fall below min AT (10), so the allelic AD scan never sees them.
- **The user's reading:** allelic counts come from the same fragments, so low depth means low AT. Depth can flag a low-coverage region as a segment, while the allelic signal is blind to it.
- **The dilemma:** such a segment could be
  1. **a real one-copy loss.** Expected: about half depth, and the remaining reads all from one haplotype (LOH), so |dBAF| → 0.5;
  2. **a homozygous deletion**, possibly subclonal. Expected: depth → 0 in the cells that carry it, and no alleles left to measure. Pooled reads from cells without the deletion look balanced, so the allelic signal **cannot** confirm it. chr9 ~20–25 Mb (the CDKN2A region) is where the user expects a small focal event; depth catches it in most clusters, allelic misses it;
  3. **technical low coverage**: closed chromatin (ATAC is accessibility-driven), low SNP density, mappability. Expected: low depth, low AT, balanced alleles, the same in every cell type.

  **We have no normal reference cells.** So we cannot normalise depth against normal cells, and we cannot ask whether the dip is also present in normal cells.

### Evidence so far

The probe in Appendix B compares each segment with ±50 bins of flank, using the genome-pool depth segmentation without merge rounds:

| Region (all clusters) | Depth | AT | Observed \|dBAF\| | Expected if balanced |
|---|---|---|---|---|
| chr7 8.3–10.7 Mb | ×0.27–0.41 | ×0.09–0.15 | 0.09–0.23 | 0.08–0.20 |
| chr7 51.4–53.5 Mb | ×0.27–0.36 | ×0.08–0.11 | 0.10–0.25 | 0.09–0.19 |
| chr7 67.4–70.0 Mb | ×0.15–0.24 | ×0.03–0.09 | 0.13–0.25 | 0.14–0.25 |
| chr7 70.0–73.3 Mb | ×0.40–0.52 | ×0.22–0.32 | 0.08–0.12 | 0.05–0.10 |
| chr9 22.5–26.6 Mb (cluster0, cluster1, all) | ×0.55 | ×0.26–0.29 | 0.04–0.06 | 0.03–0.05 |

- **AT drops about 3× more than depth.** These regions are short of reads overlapping heterozygous SNPs, not just short of reads.
- **No LOH signature anywhere:** observed imbalance is close to what coverage noise gives with balanced alleles. That rules out a clonal one-copy loss in these segments.
- **The chr7 dips sit at identical coordinates in every cluster and in `all`.** That is consistent with a property of the region, but also with an event shared by all cells. We can't tell which without a reference.
- **Allelic AD at the defaults finds few breakpoints:** 4–18 per source genome-wide, against 117–156 depth breakpoints per cluster.

## 3. Questions for you

**A. How do similar tools handle low-coverage regions in loss calls, and how do they separate one-copy loss, homozygous deletion and technical dropout?**
- Start with scATAC and scRNA CNV tools: Numbat, CopyKAT, inferCNV, CaSpER, Alleloscope, epiAneufinder / pyEpiAneufinder, CHISEL / SIGNALS (scDNA, allele-specific), HoneyBADGER, XClone, scATAC-specific approaches (e.g. copyscAT).
- For each, report:
  - how it models coverage dropout;
  - how it uses allelic evidence when coverage is low (pooling, priors, HMM emission models);
  - whether it calls homozygous deletions;
  - what it uses as the reference.

**B. Adjacent fields.** Methods that separate a true signal change from missing or low-quality data in a 1-D sequence:
- change-point detection with missing or censored data and heteroscedastic noise (e.g. segmentation with count models, Poisson/negative-binomial change-points, Bayesian online change-point detection);
- image and signal processing: inpainting, edge detection with uncertainty masks, masked or partial convolution;
- joint segmentation of two coupled tracks: multi-track HMMs, joint change-point models such as FACETS/ASCAT-style allele-specific segmentation of logR + BAF from bulk tumour sequencing;
- machine-learning segmentation that outputs uncertainty or an explicit "insufficient evidence" state.

We want something that can say **"no allelic evidence possible here"** rather than **"balanced"**.

**C. Inferring normal cells without labels.**
- Look at scATAC foundation and representation models: the user mentioned atacformer, ChromFound and EpiFoundation. Also check scBasset, SCALE, PeakVI/scVI-family, Signac/ArchR label transfer, and cell-type annotation from accessibility atlases.
- Can any of them assign cell types or identify likely non-malignant cells (immune, stromal, endothelial) in a tumour scATAC sample?
- Check CNV-based approaches too: CopyKAT's diploid/aneuploid classification, Numbat's clone inference without a reference, inferCNV with an internal reference, using the most diploid-looking cluster.
- Check reference-free normalisation: a panel of normals from public scATAC data (e.g. healthy tissue atlases) used as an external reference for region-level accessibility and coverage expectations.

**D. Per-cell evidence for homozygous deletions.** Cells carrying a homozygous deletion should have structural zeros: exactly zero fragments in the region, which makes the cells bimodal. Cells without it should be Poisson-sparse. Find statistical tests or models for zero-inflation and bimodality in sparse single-cell counts that would work at 100 kb–5 Mb region sizes.

**E. External priors for technical low coverage:** blacklists, mappability tracks, constitutively closed or low-accessibility regions in public scATAC atlases, SNP-density maps. How do tools down-weight such regions?

### What a good answer looks like

- A ranked shortlist of methods: what each does, why it fits our setting (scATAC, 100 kb bins, pseudobulk plus per-cell, phased HA/HB, no normals), its inputs and outputs, licence, and code availability.
- A clear distinction between published, validated practice and ideas transferred from other fields.
- Concrete next experiments we could run on this dataset, for example a segment-level LOH test, a per-cell zero-bimodality test, a normal-cell inference route.
- Citations with links. Do not overclaim: the repo's rule is that exploratory results are described as "consistent with", never "proves" (see [`docs/VALIDATION.md`](../VALIDATION.md)).

## 4. Where everything is

All paths are relative to `/home/mami/Desktop/EGA/pyEpiAneufinder/analysis/2026-09-16-interactive-app-planning/` unless absolute.

### App and docs

- **App:** `scientific-explorer/` (FastAPI + React). Start a private server with `python3 run.py --port 8799` from `scientific-explorer/`. **Never stop or restart the server on port 8766; it is the user's.** Open `http://127.0.0.1:8799`.
- **Docs:**
  - `scientific-explorer/README.md`
  - `docs/USAGE.md` (Integration tabs, Merge rounds, Allelic segmentation)
  - `docs/SCIENTIFIC_WORKFLOW.md`
  - `docs/ARCHITECTURE.md`
- **Session record:** `docs/history/2026-10-28-merge-rounds-allelic-segmentation.md`
- **Specs:**
  - `docs/superpowers/specs/2026-10-28-allelic-segmentation-design.md`
  - `docs/superpowers/specs/2026-10-27-merge-rounds-design.md`
- **Validation records:**
  - `validation/allelic/` (allelic segmentation per signal, reproduction gate)
  - `validation/depth-genome-pool/` (depth breakpoints old vs new global test)
- **Deferred idea:** `docs/notes/2026-10-27-tau-calibration.md`, calibrating the merge threshold from a known CNV.

### Saved sessions

Load these in the app with the session bar: saved sessions → load.

- `scientific-explorer/sessions/session-2026-10-29-12-10.json` is the user's view of the issue:
  - Integration · allelic segmentation tab, signal **shrunk |dBAF|**, chromosome **chr7**, clusters 0–4;
  - AD defaults (min AT 10, k 6, 1000 permutations, alpha 0.001, min child 20, min terminal 50, genome global pool);
  - one merge round adopted.
  - Run segmentation; the result is cached and instant.
  - Then switch to **chr9** and look at 20–25 Mb.
- `scientific-explorer/sessions/demo-integration.json` is an older demo session.

### Data

- **Depth counts:** `integration-story/data/count_matrix.h5ad`, 3,460 cells × 25,337 bins of 100 kb. `X` is GC-corrected normalised depth; `layers['raw']` holds fragment counts; `var` has seq/start/end/GC.
- **Integrated object:** `integration-story/data/integrated.h5mu`. Cells, WNN clusters (`wnn_leiden_0.3` gives 5 clusters), depth and haplo modalities, and an abs-dBAF matrix in `mod/haplo/X`.
- **Allelic matrices:** `/home/mami/Desktop/EGA/haplotype_aware_cnv/hapcnv_outs_from_updated_phasing/aligned_matrices/{HA,HB,AT}.parquet`. Bins × barcodes on the same 25,337-bin grid. Columns are barcodes in a different order from the h5mu; match them by name.
- **Phasing inputs:** `/home/mami/Desktop/EGA/SNP_project/alpha/haplo_initial/updated_phasing/`.
- **Depth breakpoints used by the app:** `integration-story/data/genome_pool/raw_accepted_breakpoints.tsv`, from a genome-wide global permutation pool as in original pyEpi. Its cluster and per-cell CN tables are in `genome_pool/derived/`.
  - The earlier per-chromosome-pool breakpoints are in `integration-story/data/raw_accepted_breakpoints.tsv`, with tables in `data/derived/`.
- **Cached pseudobulks:**
  - depth (`pb_x`, 6 sources × bins): `scientific-explorer/cache/integration_*.npz`, the most recent file;
  - allelic sums (`ha`, `hb`, 6 × bins, sources in the order cluster0..cluster4, all): `scientific-explorer/cache/allelic/pseudobulk_*.npz`.

### Code

- **Depth pipeline:** `pyEpiAneufinder/` (package copy) and the original clone at `/home/mami/Desktop/EGA/pyEpiAneufinder/pyEpiAneufinder/`. See `get_breakpoints.py` (AD scan, permutation tests) and `assign_somy.py` (CN calling, `best_s` scale search).
- **Constrained depth scan used for the cluster breakpoints:** `/home/mami/Desktop/EGA/SNP_project/beta_revised/src/pyEpiAneufinder/get_breakpoints.py` and `.../src/segmentation/run_bootstrap_diagnostics.py`.
- **App modules:** `scientific-explorer/server/allelic_scan.py`, `allelic_data.py`, `allelic_api.py`, `cn_recall.py`, `merge_rounds.py`, `call_noise.py`.
- **Earlier allelic work:**
  - `/home/mami/Desktop/EGA/SNP_project/alpha/haplo_initial/`: PELT on pseudobulk |dBAF|; oriented dBAF for per-cell scores. See `AGENT_TASKS/implement_oriented_aware_ai.md` and `AGENT_TASKS/debate/oriented_ai_debate_conclusion.md`.
  - `/home/mami/Desktop/EGA/haplotype_aware_cnv/`: an older aligned-matrix pipeline on a different, 18,738-bin grid; not used now.

## 5. Constraints and open facts

- **No labelled normal cells.** Tumour type, genome build and expected known events beyond chr9 ~20–25 Mb are not recorded here; ask the user rather than assume. The coordinates are 1-based 100 kb bins; the build is probably hg38, but verify.
- **Keep exploratory results labelled as exploratory.** Do not modify the app, the data or the validation records. This is a read-and-research task.

## Appendix A: the expected-imbalance yardstick used above

For a bin with AT = n reads and balanced alleles, HB ~ Binomial(n, 0.5).

- **Expected noise:** E|HB/n − 0.5| is the |dBAF| expected from noise alone. We compare the AT-weighted mean of observed per-bin |dBAF| (bins with AT ≥ 2) with the AT-weighted mean of this expectation.
- **Under LOH** the observed value approaches 0.5 regardless of n ≥ 2.
- **Bins with AT = 1 are uninformative.** One read always gives |dBAF| = 0.5.

## Appendix B: probe script (read-only)

Run from `scientific-explorer/`. It reads the cached pseudobulks and the genome-pool segments.

```python
import numpy as np, pandas as pd, anndata as ad, glob, os
from scipy.stats import binom
S = ["cluster0","cluster1","cluster2","cluster3","cluster4","all"]
pb = np.load(max(glob.glob('cache/integration_*.npz'), key=os.path.getmtime))['pb_x'].astype(float)
al = np.load(max(glob.glob('cache/allelic/pseudobulk_*.npz'), key=os.path.getmtime)); HA, HB = al['ha'], al['hb']; AT = HA + HB
var = ad.read_h5ad('../integration-story/data/count_matrix.h5ad', backed='r').var.reset_index(drop=True)
seq, start = var.seq.astype(str).values, var.start.values
cs = pd.read_csv('../integration-story/data/genome_pool/derived/karyogram_cluster_segments.tsv', sep='\t')
expabs = lambda n: float((binom.pmf(np.arange(n+1), n, .5) * np.abs(np.arange(n+1)/n - .5)).sum())
def describe(si, b0, b1):
    fl = np.r_[max(0, b0-50):b0, b1+1:min(len(seq), b1+51)]; fl = fl[seq[fl] == seq[b0]]
    seg = np.arange(b0, b1+1); m = AT[si, seg] >= 2
    obs = np.average(np.abs(HB[si, seg][m]/AT[si, seg][m] - .5), weights=AT[si, seg][m])
    exp = np.average([expabs(int(n)) for n in AT[si, seg][m]], weights=AT[si, seg][m])
    return dict(depth_ratio=pb[si, seg].mean()/pb[si, fl].mean(), at_ratio=AT[si, seg].mean()/AT[si, fl].mean(),
                kept_frac=(AT[si, seg] >= 10).mean(), obs_abs_dbaf=obs, balanced_abs_dbaf=exp)
for r in cs[(cs.chromosome == 'chr7') & (cs.n_bins <= 60) & (cs.cn_call < 1)].itertuples():
    print(r.group, start[r.bin_start]/1e6, start[r.bin_end]/1e6, describe(S.index(r.group), r.bin_start, r.bin_end))
```
