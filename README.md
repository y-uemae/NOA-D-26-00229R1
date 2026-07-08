# TP53–LAG3 in WHO grade 4 glioma — analysis code

Analysis code and result tables for:

> **TP53 mutation is associated with selective upregulation of LAG3 independent of immune infiltration and antigen-presentation signatures in WHO grade 4 glioma**

Manuscript: *Neuro-Oncology Advances* (NOA-D-26-00229R1)

This repository contains the R scripts and the derived result tables needed to reproduce the analyses, figures, and supplementary tables in the manuscript. Raw sequencing/proteomic data are **not** redistributed here; they are obtained from the public sources described below.

---

## Repository layout

```
.
├── NOA-D-26-00229R1.Rproj      # open this in RStudio first (sets the project root)
├── README.md
├── LICENSE                     # MIT
├── scripts/                    # 59 analysis scripts (run in numeric step order)
├── results/
│   └── TP53/
│       ├── 20260221/           # main-analysis result tables (per step folder)
│       └── 20260630/           # revision analyses (CPTAC / ABSOLUTE / ComBat / log2FC)
└── data/
    └── raw/                    # (empty) place downloaded raw data here — see below
```

Paths inside every script are resolved with the [`here`](https://here.r-lib.org/) package, so the code is portable across machines **as long as the project is opened via the `.Rproj` file** (this fixes the project root). No absolute paths remain in the code.

---

## Software requirements

- R (≥ 4.3 recommended)
- RStudio (to open the `.Rproj`)
- Key R packages: `here`, `tidyverse`, `metafor`, `sva` (ComBat), `data.table`, `ggplot2`, `patchwork`, `ragg`. Bioconductor packages (e.g. `sva`) install via `BiocManager::install()`.

---

## Data acquisition

No raw data are included. To reproduce the pipeline from scratch, download the following into `data/raw/` as indicated. (The derived result tables in `results/` let you re-run most downstream/figure scripts without the raw data.)

**1. GDC (Genomic Data Commons)** — https://portal.gdc.cancer.gov/
WHO grade 4 glioma RNA-seq, WXS/MAF mutation, and clinical data (TCGA / CPTAC / HCMI). Cohort-selection queries are defined in `scripts/20260221_Step01 ...`. Place under `data/raw/GDC/glioma/`. (GDC open-access REST API; no authentication.)

**2. GLASS consortium** — Synapse `syn17038081` (https://www.synapse.org/#!Synapse:syn17038081)
Non-TCGA, WXS RNA-seq TPM matrix (`data_mrna_seq_tpm.txt`). Place at `data/raw/external_validation/difg_glass/data_mrna_seq_tpm.txt`.

**3. ABSOLUTE purity calls (PanCanAtlas)** — https://gdc.cancer.gov/about-data/publications/pancanatlas
File `TCGA_mastercalls.abs_tables_JSedit.fixed.txt`, used by the R1-3b ABSOLUTE sensitivity analysis. Place under the path referenced in `scripts/R1 3b step2 ...`.

---

## How to reproduce

1. Open `NOA-D-26-00229R1.Rproj` in RStudio (this sets the project root used by `here()`).
2. Install the required packages (above).
3. Run the scripts in `scripts/` in **numeric step order** (Step01 → Step02 → … → Step32, then the `_20260630` revision scripts). Each script reads from and writes to the corresponding `results/…` sub-folder.
4. To re-run only a downstream/figure step, ensure its input tables exist in `results/` (they are provided) and run that script directly.

### Important run-order notes

- **ComBat (R1-4, Table S7):** run `Revision combat_20260630.R` and then, **in the same R session**, `Revision combat patch_20260630.R`. The patch relies on in-memory objects (`Mcb`, `Msub`, `meta`) created by the first script and writes the final `ComBat_TP53_coefficients.csv`.
- **Checkpoint figure (Fig. 2):** `Step16c checkpoint figure.R` reads `step16a_v2_results.csv` (produced by `step16a_checkpoint_specificity_v2.R`). Run the v2 specificity script before the figure script.
- **ABSOLUTE (R1-3b, Table S6):** the analysis script is `R1 3b step2 absolute join refit_20260630.R`; it requires the PanCanAtlas ABSOLUTE table (see Data acquisition).

---

## Script → figure/table map

Main figures:

| Output | Script(s) |
|---|---|
| Fig. 1 (combined) | `step13_combined_figure1.R` |
| Fig. 1A (GDC WT/Mut, regression) | `step09_main.R`, `Step09b regression.R`, `step13c plot main wtmut v4.R` |
| Fig. 1B (GLASS) | `step10_glass.R` |
| Fig. 1C (meta-analysis forest) | `step11_meta.R`, `Step11c forest plot.R`, `Step13b forest plot v3.R` |
| Fig. 1D (TP53 subgroups) | `step12b_subgroup_classify.R`, `step12c_subgroup_analysis.R`, `step13a_plot_subgroup.R` |
| Fig. 2 (checkpoint specificity) | `step16a_checkpoint_specificity_v2.R`, `Step16b checkpoint regression.R`, `Step16c checkpoint figure.R` |
| Fig. 3A/3B (immune-score regression) | `Step17 immune score regression.R` |

Supplementary figures:

| Output | Script |
|---|---|
| Fig. S1 (grade 2/3) | `step_S1_grade23_violin.R` |
| Fig. S2 (quantile regression) | `Step25 quantile regression.R` |
| Fig. S3 (residual analysis) | `Step20 residual analysis.R` |
| Fig. S4 (immune score vs TP53) | `Step21 immune score tp53 diff.R`, `step31_visualization.R` |
| Fig. S5 (TP53 × immune interaction) | `Step18 interaction emm.R` |
| Fig. S6 (stratified robustness) | `Step19 robustness stratified_v3.R` |
| Fig. S8 (LAG3 distribution) | `Step22 lag3 distribution comparison.R` |
| Fig. S9 (TP53 × cohort interaction) | `Step24 interaction source.R` |
| Fig. S10 (permutation) | `Step26 permutation.R` |
| Fig. S11 (standardized effect) | `Step23 standardized effect.R` |
| Fig. S12 (GLASS visualization) | `step32c_glass_visualization.R` |

Supplementary tables (revision analyses):

| Output | Script |
|---|---|
| Table S5 (purity/hallmark) | `step32c_make_suppl_table_s5.R` |
| Table S6 (ABSOLUTE purity, R1-3b) | `R1 3b step2 absolute join refit_20260630.R` |
| Table S7 (ComBat batch correction, R1-4) | `Revision combat_20260630.R` + `Revision combat patch_20260630.R` |
| Table S8 (CPTAC protein detection, R1-1) | `Cptac lag3 tp53 analysis_20260630.R`, `Cptac checkpoint detection check_20260630.R` |
| Table S2 legend (log2FC, R2-1) | `Revision log2fc_20260630.R` |

The upstream data-construction scripts (Step01–Step08, Step27–Step32; cohort assembly, expression matrices, ESTIMATE purity, ssGSEA, screening) produce the intermediate tables that the figure/table scripts above consume.

---

## Notes on data handling

- Immune-score gene sets: T-cell (CD3D/E/G, CD8A/B, GZMA/B, PRF1), APM (B2M, TAP1/2, TAPBP, HLA-A/B/C, NLRC5), IFN-γ (STAT1, IRF1, IRF9, CXCL9/10/11, GBP1/2/4/5, IDO1).
- Tumor purity: ESTIMATE (expression-based) throughout; ABSOLUTE (DNA-based) used only for the TCGA-restricted R1-3b sensitivity analysis.
- A small number of diagnostic / structure-checking helper scripts used only during interactive development are **not** included, as they do not generate any manuscript figure or table.

---

## License

Code is released under the MIT License (see `LICENSE`). Third-party data (GDC/TCGA/CPTAC/HCMI, GLASS, PanCanAtlas) remain subject to their respective terms of use.

## Citation

If you use this code, please cite the manuscript (details to be added on acceptance) and this repository (DOI to be added after archiving).
