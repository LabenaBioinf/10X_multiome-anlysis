# Single-nucleus multiome analysis of paired MIS-C flare and remission samples

Analysis code for the single-nucleus part of the article

> **Interferon-linked transcription and chromatin accessibility extend to naive T cells in multisystem
> inflammatory syndrome in children: a paired multi-omic cohort study.**
> Jenko Bizjan B, Kovač J, Tesovnik T, Šket R, Debeljak M, Zajc Avramovič M, Emeršič N, Vrhovšek B,
> Sabolić I, Smole A, Avčin T. (submitted)

Ten children with multisystem inflammatory syndrome in children (MIS-C) gave one PBMC sample at flare
and one at remission. Both samples were profiled with 10x Genomics Single Cell Multiome ATAC + Gene
Expression (RNA and chromatin accessibility in the same nuclei, 89,154 nuclei in total).
**The patient, not the nucleus, is the unit of analysis.** Counts are summed per patient × time point ×
population (pseudobulk) and compared with paired designs.

| | |
|---|---|
| Companion repository | [CISLD/MISC-bulk-cytokines-exome](https://github.com/CISLD/MISC-bulk-cytokines-exome): whole-blood RNA-seq, plasma cytokines, exome and HLA typing, public cohorts and **all figure scripts** |
| Data | Processed single-nucleus data: GEO [GSE301727](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301727). No data are included in this repository |
| Archive | Zenodo, all versions: [10.5281/zenodo.24120695](https://doi.org/10.5281/zenodo.24120695). The code used in the article is release **v2.1** |
| Licence | MIT |

## Workflow

```
Cell Ranger ARC outputs (per library)
        │
        ▼
SC_GEX-ATAC.R + scripts/            QC, normalisation, Harmony, reference mapping (Azimuth PBMC)
        │                           → harmonized_seurat.rds, per-population cell counts
        ▼
article_analyses/scripts/
  ├─ a6_composition_check.R          population proportions (exact paired Wilcoxon)
  ├─ pb/                             RNA pseudobulk: edgeR, camera, key-gene panel
  ├─ ambient/                        SoupX ambient-RNA sensitivity analysis
  └─ rev/                            TRBV11-2, chromVAR, FigR, sequencing depth, figure inputs
        │
        ▼
Figures: scripts/figures.sb in the companion repository
```

## Repository structure

```
SC_GEX-ATAC.R            main processing script (sources scripts/)
Singularity.def          container recipe for the processing pipeline (R 4.1.3)
scripts/                 functions used by SC_GEX-ATAC.R
article_analyses/
  scripts/               analyses reported in the article (R/Python + SLURM .sb job files)
  docs/                  written pre-declaration of the chromVAR extension
  logs/                  SLURM logs (empty)
CITATION.cff, LICENSE
```

## 1. Processing: from Cell Ranger ARC to the annotated object

Article Methods: *Single-nucleus data processing and composition*.

`SC_GEX-ATAC.R` expects the Cell Ranger ARC (v2.0) outputs of each library (`atac_peaks.bed`,
`per_barcode_metrics.csv`, `atac_fragments.tsv.gz`, `filtered_feature_bc_matrix.h5`) and the Azimuth
PBMC reference ([Hao et al. 2021](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat)).
It sources:

| Script | Step |
|---|---|
| `scripts/preprocess_and_QC.R` | Seurat/Signac objects per library; scDblFinder doublets; QC filters (1,000–100,000 fragments in peaks, TSS enrichment > 1, nucleosome signal < 4, > 500 UMIs, > 250 genes, < 20% mitochondrial reads) |
| `scripts/normalize.R` | SCTransform normalisation |
| `scripts/batch_effect.R` | Harmony on the RNA principal components |
| `scripts/annotate_and_integrate.R` | reference mapping to the Azimuth PBMC reference (30 populations), UMAP |
| `scripts/DE.R` | cell counts per population and condition (`cell_type_count.txt`). Its cell-level FindMarkers results are not used in the article |

The output, `harmonized_seurat.rds`, is the input of all analyses in `article_analyses/`.

## 2. Analyses reported in the article

Paths in the job files are written as `/path/to/...`. Set `WORKDIR` in each `.sb` file to the
`article_analyses` folder, which must also contain the `results/` tree and, in `tables_extracted/`, the
gene-panel tables that some scripts read (not included).

| Article item | Script (job file) | What it does |
|---|---|---|
| Figure 3C, Tables S4, S4b | `scripts/a6_composition_check.R` (`.sb`), `scripts/rev/composition_sensitivity_each20.py` | exact paired Wilcoxon tests of population proportions, BH over 30 populations; sensitivity analysis requiring ≥ 20 nuclei in each sample |
| Table S10 (TRBV11-2); ATAC pseudobulk | `scripts/rev/object_level_checks.R` (`.sb`) | *TRBV11-2* detection per sample; writes the ATAC pseudobulk used by chromVAR |
| Methods, sequencing depth | `scripts/rev/seq_depth.R` (`.sb`) | median raw GEX and ATAC read pairs per nucleus per library |
| Table S5, Figure 4B | `scripts/pb/pseudobulk_edgeR.R` (`.sb`) | RNA pseudobulk (≥ 10 nuclei per sample); edgeR quasi-likelihood `~ patient + condition`, `filterByExpr`, TMM, genome-wide BH; camera on the 50 Hallmark sets; populations with ≥ 4 complete pairs (20 of 30) |
| Tables S5b, S5c, S6; Figures 4C, 5A, 8 | `scripts/pb/pseudobulk_followup.R` (`.sb`) | 120-gene key-gene panel (selects the genes shown; estimates and genome-wide FDR from the primary model); targeted IL-12/STAT4 and IL-27/gp130 gene sets with camera |
| Tables S5e, S6b; genes marked ‡ | `scripts/pb/pseudobulk_designfilter.R`, `pseudobulk_designfilter_panel.R` (`.sb`) | sensitivity analysis with the alternative full-design expression filter; source of significance for panel genes below the primary filter |
| Tables S5d, S5f, S6c | `scripts/ambient/soupx_pseudobulk.R` (`.sb`), `soupx_concordance.py`, `soupx_text_numbers.py` | SoupX 1.6.2 ambient-RNA correction on the unfiltered Cell Ranger ARC matrices; reruns the pseudobulk models on the corrected counts; log2 fold-change concordance and the numbers cited in the text |
| Figure 6, Table S8e | `scripts/rev/gdt_chromvar.R` (`.sb`) | chromVAR deviations (35 cisBP motifs and their families) on the ATAC pseudobulk (≥ 20 nuclei per sample); limma `~ patient + condition`; BH within population |
| Figure 6 (naive and regulatory cells), Table S8f | `scripts/rev/chromvar_naive_extension.R` (`.sb`) | the same analysis in naive CD4⁺, naive CD8⁺, naive B and regulatory T cells; pre-declaration in `docs/` |
| Table S8g | `scripts/rev/chromvar_depth_sensitivity.R` (`.sb`) | sequencing-depth sensitivity: bias-corrected deviations and equal-depth binomial thinning (five seeds), same peaks, motifs, model and FDR families |
| Tables S8a–b, Supplementary Figures S4–S5 | `scripts/rev/figr_differential_edges.R` (`.sb`) | FigR in CD8⁺ TEM and γδ T cells, flare and remission nuclei analysed together; k-NN smoothing within library; motif-supported TF–DORC pairs; flare–remission difference tested by exact within-patient permutation |
| Table S8d | `scripts/rev/figr_module_test.R` (`.sb`), `figr_module_test_minlib.R` (`figr_module_test_min20.sb`) | transcription factor module tests with per-patient estimates |
| Table S8c | `scripts/rev/figr_nk.R` (`.sb`) | the same FigR tests in NK cells |
| Tables S8c–e (assembly) | `scripts/rev/consistency_rule.py` → `figr_supp_tables_v4.py` → `figr_modules_table_s8d.py` | collect module, chromVAR and NK results into table files; `figr_modules_table_s8d.py` writes Table S8d as reported (FDR over the modules reported in the article) |
| Figures 3, 5, 6 (inputs) | `scripts/rev/figure_data_sc.R` | per-sample values plotted in the figures |

### Run order

```
SC_GEX-ATAC.R
 ├─ a6_composition_check → composition_sensitivity_each20
 ├─ object_level_checks → gdt_chromvar → chromvar_naive_extension → chromvar_depth_sensitivity
 ├─ seq_depth
 ├─ pseudobulk_edgeR → pseudobulk_followup, pseudobulk_designfilter*
 │                  └─ ambient/soupx_pseudobulk
 ├─ figr_differential_edges → figr_module_test*, figr_nk → consistency_rule
 │                          → figr_supp_tables_v4 → figr_modules_table_s8d
 └─ figure_data_sc
then scripts/figures.sb in the companion repository
```

### Note on names in the FigR scripts

Variables and module names containing `published` or `S8` (for example `TBX21_published`) refer to the
TF–DORC list of an earlier version of the article, not to the final Supplementary Table S8.
These modules are not reported. `figr_modules_table_s8d.py` removes them before the FDR is computed.

## Software

The processing pipeline runs in the container built from `Singularity.def` (R 4.1.3, Seurat, Signac,
scDblFinder, FigR and the other packages listed in the recipe):

```bash
sudo singularity build SC_GEX-ATAC.sif Singularity.def
singularity exec --bind /path/to/data:/home/results --bind /path/to/resources:/home/resources \
  SC_GEX-ATAC.sif Rscript SC_GEX-ATAC.R
```

The article analyses are SLURM batch files (`sbatch <job>.sb`). Each step runs in an Apptainer container
(`apptainer exec`) with Seurat 5.1.0, Signac 1.13.0, edgeR, limma, chromVAR, FigR, msigdbr and ggplot2.
The Python scripts need Python ≥ 3.8 (standard library). The versions of all tools are given in the
article's Methods.

## Licence and citation

MIT licence (see `LICENSE`). If you use this code, please cite the article (see `CITATION.cff`).
The code used in the article is release v2.1, archived on Zenodo
(all versions: https://doi.org/10.5281/zenodo.24120695).
