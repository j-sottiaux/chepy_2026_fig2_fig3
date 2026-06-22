# Reproducible analysis workflow for Figure 2 and Figure 3

This repository contains the R scripts and Shiny applications used to reproduce the proteomics, transcriptomics, and pathway-enrichment analyses associated with Figure 2 and Figure 3 of the manuscript.

The repository is intended for scientific reviewers who want to inspect, rerun, or audit the analysis workflow. The code is organized so that raw data are read from `data/00_raw/`, processed outputs are written to `data/`, and publication-style figures are written to `figures/`.

## Repository structure

```text
.
├── chepy_2026_fig2_fig3.Rproj
├── codex_shiny/
│   └── app.R
├── data/
│   ├── 00_raw/
│   ├── 01_proteo_toptables/
│   ├── 02_gsea_results_raw/
│   ├── 03_gsea_results_filtered/
│   ├── 04_ora_results_raw/
│   ├── 05_ora_results_filtered/
│   ├── 06_enrichment_merged/
│   ├── 07_enrichments_integration/
│   ├── 08_codex_data/
│   └── 09_proteo_toptables_extended/
├── diff_analysis_explorer_shiny/
│   └── app.R
├── figures/
│   ├── 01_exploratory_analysis/
│   ├── 02_differential_analysis/
│   └── 03_enrichments_integration/
├── renv/
├── renv.lock
└── scripts/
    ├── 01_functions.R
    ├── 02_proteomics_analysis.R
    └── 03_transcriptomics_analysis.R
```

## Software requirements

The project uses R and the package environment recorded in `renv.lock`. Reviewers should restore the environment before running the analysis.

```r
install.packages("renv")
renv::restore()
```

Main R packages used by the analysis include:

- `data.table`, `readxl`, `writexl`, `tidyverse`
- `ggrepel`, `ggthemes`, `pheatmap`, `grid`
- `limma`, `DESeq2`, `ashr`, `mixOmics`, `clusterProfiler`
- `AnnotationDbi`, `org.Hs.eg.db`

## Input data expected in `data/00_raw/`

The scripts assume that the following files are available in `data/00_raw/`:

| Input type | Expected filename pattern | Used by |
|---|---|---|
| Proteomics LFQ matrices | `proteo*.txt` | `scripts/02_proteomics_analysis.R` |
| Transcriptomics raw counts | `transcripto_raw_counts.txt` | `scripts/03_transcriptomics_analysis.R` |
| Transcriptomics background matrix | `transcripto_endopath_vst.csv` | `scripts/03_transcriptomics_analysis.R` |
| Gene-set databases | `*.gmt` | Proteomics and transcriptomics GSEA/ORA |

The `.gmt` gene-set databases (v2026_1_hs_symbols) can be downloaded from MSigDB: [GOBP](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2026.1.Hs/c5.go.bp.v2026.1.Hs.symbols.gmt) and [Reactome](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2026.1.Hs/c2.cp.reactome.v2026.1.Hs.symbols.gmt).

The proteomics workflow also reads and writes top tables in:

- `data/01_proteo_toptables/`
- `data/09_proteo_toptables_extended/`

## How to run the analysis

Open the RStudio project file:

```text
chepy_2026_fig2_fig3.Rproj
```

Then run the scripts in this order from the project root:

```r
source("scripts/02_proteomics_analysis.R")
source("scripts/03_transcriptomics_analysis.R")
```

`02_proteomics_analysis.R` should be run before `03_transcriptomics_analysis.R`, because the transcriptomics enrichment section reuses the GMT-derived objects created during the proteomics script, including `raw_gmt` and `ref_background_genes`.

## Analysis overview

### 1. Shared functions and metadata

`scripts/01_functions.R` defines reusable functions and global metadata used by both analysis scripts.

It includes:

- biological condition labels, plotting labels, and color palettes;
- proteomics sample-group definitions for cytoplasmic and nuclear fractions;
- matrix import and cleaning functions;
- PCA and PLS-DA helper functions;
- proteomics differential abundance functions;
- DESeq2 result extraction and gene-symbol mapping functions;
- GSEA and ORA helper functions;
- figure-generation functions for volcano plots, heatmaps, dimensionality reduction, and enrichment integration;
- Shiny helper functions used by the CODEX and differential-analysis explorer applications.

### 2. Proteomics workflow

`scripts/02_proteomics_analysis.R` performs the proteomics analysis.

Main steps:

1. Load proteomics LFQ matrices from `data/00_raw/`.
2. Clean protein identifiers and retain LFQ intensity columns.
3. Remove predefined outlier samples:
   - cytoplasm: column 45;
   - nucleus: column 2.
4. Run exploratory analyses with `mixOmics`:
   - PCA;
   - PLS-DA.
5. Run differential abundance analysis with `limma` for biological conditions versus healthy controls.
6. Generate volcano plots.
7. Run extended pairwise differential analyses across biological conditions.
8. Read GMT gene-set databases.
9. Run GSEA on ranked gene lists.
10. Run ORA on significant gene lists.
11. Merge and integrate enrichment results.
12. Generate CODEX-compatible enrichment data.

Default significance thresholds:

```r
padj_threshold <- 0.05
logFC_threshold <- 1
```

Default enrichment settings:

```r
padj_enrichment_cutoff <- 1
min_genecount_cutoff <- 15
max_genecount_cutoff <- 250
permutations_gsea <- 50000
multiple_testing_correction <- "BH"
eps_limit_gsea <- 1e-30
```

### 3. Transcriptomics workflow

`scripts/03_transcriptomics_analysis.R` performs the transcriptomics analysis.

Main steps:

1. Load raw count data from `data/00_raw/transcripto_raw_counts.txt`.
2. Define sample metadata for three conditions:
   - `ATA`;
   - `HC`;
   - `IFNa`.
3. Build a `DESeqDataSet` with design `~ condition`.
4. Filter low-count genes by retaining genes with counts greater than 10 in at least 3 samples.
5. Run DESeq2 differential expression analysis.
6. Apply variance-stabilizing transformation for exploratory visualization.
7. Run PCA and sample-distance heatmap analyses.
8. Extract contrasts:
   - `ATA` vs `HC`;
   - `ATA` vs `IFNa`;
   - `IFNa` vs `HC` for QC, although this contrast is not included in the final merged result lists.
9. Apply log-fold-change shrinkage with `ashr` for visualization-oriented results.
10. Map Ensembl identifiers to gene symbols using `org.Hs.eg.db`.
11. Generate volcano plots and heatmaps.
12. Run transcriptomics GSEA and ORA.
13. Merge and integrate transcriptomics enrichment results.

## Main output folders

| Output folder | Contents |
|---|---|
| `figures/01_exploratory_analysis/01_proteomics_datasets/` | Proteomics PCA and PLS-DA plots |
| `figures/01_exploratory_analysis/02_transcriptomics_dataset/` | Transcriptomics PCA plot |
| `figures/02_differential_analysis/01_proteomics_datasets/` | Proteomics volcano plots |
| `figures/02_differential_analysis/02_transcriptomics_dataset/` | Transcriptomics volcano plots and heatmaps |
| `figures/03_enrichments_integration/` | Integrated enrichment figures |
| `data/01_proteo_toptables/` | Proteomics differential abundance top tables |
| `data/02_gsea_results_raw/` | Raw GSEA outputs |
| `data/03_gsea_results_filtered/` | Filtered GSEA outputs |
| `data/04_ora_results_raw/` | Raw ORA outputs |
| `data/05_ora_results_filtered/` | Filtered ORA outputs |
| `data/06_enrichment_merged/` | Merged enrichment results |
| `data/07_enrichments_integration/` | Cross-enrichment integration tables |
| `data/08_codex_data/` | Data tables used by the CODEX Shiny application |
| `data/09_proteo_toptables_extended/` | Extended pairwise proteomics top tables used by the differential-analysis explorer |

## Shiny applications

Two Shiny applications are included.

### CODEX Shiny app

```r
shiny::runApp("codex_shiny")
```

This app uses processed enrichment and differential-analysis proteomic data, especially the outputs in `data/08_codex_data/`, to support interactive inspection of condition-specific enrichment results and volcano plots.

### Differential-analysis explorer Shiny app

```r
shiny::runApp("diff_analysis_explorer_shiny")
```

This app uses the extended proteomics differential abundance top tables in `data/09_proteo_toptables_extended/` to inspect pairwise comparisons across biological conditions.

## Reviewer notes

- The analysis is designed to be run from the project root. Relative paths assume the working directory is the directory containing the `.Rproj` file.
- The R environment should be restored with `renv::restore()` before rerunning the scripts.
- The scripts use deterministic seeding via `set.seed(12345)` in `scripts/01_functions.R`.
- The proteomics and transcriptomics analyses share enrichment helper objects. For that reason, run `scripts/02_proteomics_analysis.R` before `scripts/03_transcriptomics_analysis.R`.
- Some analysis products are generated as intermediate `.xlsx` files and then reused downstream; reviewers should not delete intermediate folders between steps unless they intend to regenerate the full workflow.
- Large Shiny-generated user plots, when exported interactively, are written to `figures/10_user_generated/` by the helper functions.

## Methods implemented in the scripts

The workflow uses established Bioconductor and R packages for omics analysis:

- `limma` is used for proteomics differential abundance modelling.
- `DESeq2` is used for transcriptomics count-based differential expression analysis.
- `ashr` is used for log-fold-change shrinkage in visualization-oriented DESeq2 outputs.
- `mixOmics` is used for PLS-DA-based exploratory multivariate analysis.
- `clusterProfiler` is used for GSEA and ORA enrichment analysis.

## Minimal reproducibility checklist

Before rerunning the analysis, verify that:

1. `renv::restore()` completes without package installation errors.
2. The working directory is the project root.
3. `data/00_raw/` contains the required proteomics, transcriptomics, and GMT files.
4. `scripts/01_functions.R` is sourced before either analysis script.
5. `scripts/02_proteomics_analysis.R` is run before `scripts/03_transcriptomics_analysis.R`.
6. Output folders under `data/` and `figures/` are writable.

## References

- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550.
- Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*. 2015;43(7):e47.
- Stephens M. False discovery rates: a new deal. *Biostatistics*. 2017;18(2):275-294.
- Rohart F, Gautier B, Singh A, Lê Cao KA. mixOmics: an R package for omics feature selection and multiple data integration. *PLOS Computational Biology*. 2017;13(11):e1005752.
- Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*. 2021;2(3):100141.
