# ----- 1. Script setup -----
## a. Load required libraries ----
### Core
library(data.table)
library(readxl)
library(writexl)

### Visualization / wrangling
library(tidyverse)
library(ggrepel)
library(ggthemes)
library(grid)

### Statistics
library(mixOmics)
library(limma)
library(ashr)
library(clusterProfiler)

## b. Reproducibility check ----
gc()
set.seed(12345)

## c. Significance thresholds ----
padj_threshold <- 0.05
logFC_threshold <- 1

## d. Enrichment specific settings ----
padj_enrichment_cutoff <- 1
min_genecount_cutoff <- 15
max_genecount_cutoff <- 250
permutations_gsea <- 50000
multiple_testing_correction <- "BH"
eps_limit_gsea <- 1e-30

## e. Biological conditions design ----
cytoATAp <- 1:10
nuclATAp <- c(1, 3:10)
ATAn <- 11:20
ACA <- 21:30
HC <- 31:40
cytoNS <- 41:44
nuclNS <- 41:45

cytoplasm_groups <- c(
  rep("dcSSc_ATAp", length(cytoATAp)),
  rep("dcSSc_ATAn", length(ATAn)),
  rep("lcSSc_ACA", length(ACA)),
  rep("HC", length(HC)),
  rep("NS", length(cytoNS))
)

nucleus_groups <- c(
  rep("dcSSc_ATAp", length(nuclATAp)),
  rep("dcSSc_ATAn", length(ATAn)),
  rep("lcSSc_ACA", length(ACA)),
  rep("HC", length(HC)),
  rep("NS", length(nuclNS))
)

## f. Source functions script ----
source(file = "scripts/01_functions.R")


# ---- 2. Process proteomic LFQ datasets -----
matrix_path <- list.files(
  path = "data/00_raw",
  pattern = "proteo.*\\.txt$",
  full.names = TRUE
)

proteo_raw <- lapply(matrix_path, import_matrix) %>%
  setNames(tools::file_path_sans_ext(basename(matrix_path)))

proteo_cytoplasm <- as.data.frame(proteo_raw$proteo_cytoplasm_lfq)
proteo_cytoplasm <- clean_proteo_df(proteo_cytoplasm)
proteo_cytoplasm <- proteo_cytoplasm[, -45] # Outlier removal
proteo_cyto_matrix <- as.matrix(proteo_cytoplasm)

proteo_nucleus <- as.data.frame(proteo_raw$proteo_nucleus_lfq)
proteo_nucleus <- clean_proteo_df(proteo_nucleus)
proteo_nucleus <- proteo_nucleus[, -2] # Outlier removal
proteo_nucl_matrix <- as.matrix(proteo_nucleus)


# ----- 3. Exploratory analysis -----
## a.PCA ----
proteo_cyto_pca <- compute_pca(df = proteo_cytoplasm, matrix_type = "cytoplasm")
plot_dim_reduction(proteo_cyto_pca)
proteo_nucl_pca <- compute_pca(df = proteo_nucleus, matrix_type = "nucleus")
plot_dim_reduction(proteo_nucl_pca)

## b.PLS-DA ----
proteo_cyto_plsda <- compute_plsda(df = proteo_cytoplasm, matrix_type = "cytoplasm")
plot_dim_reduction(proteo_cyto_plsda)
proteo_nucl_plsda <- compute_plsda(df = proteo_nucleus, matrix_type = "nucleus")
plot_dim_reduction(proteo_nucl_plsda)


# ----- 4. Differential analysis -----
## a. Biological conditions vs Healthy Controls (HC) ----
proteo_cyto_diff_analysis <- compute_diff_analysis(df = proteo_cyto_matrix, matrix_type = "cytoplasm")
proteo_nucl_diff_analysis <- compute_diff_analysis(df = proteo_nucl_matrix, matrix_type = "nucleus")

toptables_path <- list.files(
  path = "data/01_limma_toptables/",
  pattern = "*.xlsx",
  full.names = TRUE
)

proteo_toptables <- lapply(toptables_path, import_toptables) %>%
  setNames(tools::file_path_sans_ext(basename(toptables_path)))

create_volcano_plots(proteo_toptables)

## b. Pairwised comparisons across all biological conditions ----
proteo_cyto_extended_da <- compute_extended_diff_analysis(df = proteo_cyto_matrix, matrix_type = "cytoplasm")
proteo_nucl_extended_da <- compute_extended_diff_analysis(df = proteo_nucl_matrix, matrix_type = "nucleus")

toptables_path_extended <- list.files(
  path = "data/09_extended_toptables/",
  pattern = "*.xlsx",
  full.names = TRUE
)

### Toptables needed for the Differential Analysis Explorer (DAE) Shiny app ---
proteo_toptables_extended <- lapply(toptables_path_extended, import_toptables) %>%
  setNames(tools::file_path_sans_ext(basename(toptables_path_extended)))


# ----- 5. Gene Set Enrichment Analysis (GSEA) -----
## a. Process reference background genes ----
raw_gmt <- list.files(
  path = "data/00_raw",
  pattern = "*.gmt",
  full.names = TRUE
)

ref_background_genes <- lapply(raw_gmt, clean_gmt)
names(ref_background_genes) <- gsub(
  "_v2026_1_hs_symbols", "",
  tools::file_path_sans_ext(basename(raw_gmt))
)

## b. Process 'proteo_toptables' object into ranked lists (based on 't' value) ----
proteo_gsea_genelist <- clean_toptables_for_gsea(proteo_toptables)
names(proteo_gsea_genelist) <- gsub("_HC", "", tools::file_path_sans_ext(basename(toptables_path)))

## c. Run GSEA for each condition vs. 'ref_background_genes' ----
proteo_gsea_res <- compute_gsea(ref_background_genes, proteo_gsea_genelist)
proteo_gsea_res_filtered <- filter_gsea_results(proteo_gsea_res)
proteo_pathways_merged <- merge_enrichment_results(
  result_list = proteo_gsea_res_filtered,
  enrich_type = "gsea"
)


# ----- 6. Over-Representation Analysis (ORA) -----
## a. Process reference background genes ----
proteo_evapath <- as.data.frame(proteo_raw$proteo_evapath_lfq)
proteo_evapath <- clean_evapath(proteo_evapath)

proteo_filter <- unique(unlist(list(
  rownames(proteo_nucleus),
  rownames(proteo_cytoplasm),
  proteo_evapath$Genes
)))

xp_background_proteo <- filter_xp_proteome(ref_background_genes)
names(xp_background_proteo) <- gsub("_v2026_1_hs_symbols", "", tools::file_path_sans_ext(basename(raw_gmt)))

## b. Process 'proteo_toptables' object into unranked significant vectors (based on padj_threshold) ----
proteo_ora_genelist <- clean_toptables_for_ora(proteo_toptables)
names(proteo_ora_genelist) <- gsub("_HC", "", tools::file_path_sans_ext(basename(toptables_path)))

## c. Run ORA for each condition vs. 'xp_background_proteo' ----
proteo_ora_res <- compute_ora(xp_background_proteo, proteo_ora_genelist)
proteo_ora_res_filtered <- filter_ora_results(proteo_ora_res)
proteo_ora_pathways_merged <- merge_enrichment_results(result_list = proteo_ora_res_filtered, enrich_type = "ora")


# ----- 7. Cross-enrichment merging -----
proteo_enrich_integration <- integrate_proteo_enrich_res(proteo_gsea_res, proteo_ora_res)
plot_enrich_integration(proteo_enrich_integration)

### Data frame needed for the CODEX Shiny app ---
proteo_master_enrich_data <- process_codex_data(proteo_gsea_res, proteo_ora_res)
