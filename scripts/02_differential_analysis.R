### ** Project setup --------------------------------------------------------------

# Load required libraries
{
  library(data.table)
  library(readxl)
  library(writexl)
  
  library(tidyverse)
  library(ggrepel)
  library(ggthemes)
  
  library(mixOmics)
  library(limma)
  library(ashr)
  library(pheatmap)
  library(viridisLite)
}

# Reproducibility check
gc()
set.seed(12345)

# Source functions
source(file = "scripts/01_functions.R")

### 1/ Import & process proteomic LFQ data -----

# Import proteomic LFQ data
raw_matrix_path <- list.files(
  path = "data",
  pattern = "\\.txt$",
  full.names = TRUE
)

raw_matrix <- lapply(raw_matrix_path, import_raw_matrix) %>%
  setNames(tools::file_path_sans_ext(basename(raw_matrix_path)))

cytoplasm_df <- as.data.frame(raw_matrix$cytoplasm_matrix)
cytoplasm_df <- clean_proteo_df(cytoplasm_df)
cytoplasm_df <- cytoplasm_df[, -45]
cytoplasm_matrix <- as.matrix(cytoplasm_df)

nucleus_df <- as.data.frame(raw_matrix$nucleus_matrix)
nucleus_df <- clean_proteo_df(nucleus_df)
nucleus_df <- nucleus_df[, -2]
nucleus_matrix <- as.matrix(nucleus_df)

# Design groups for subsequent analysis
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

### 2/ Exploratory analysis / PCA, PLS-DA & heatmap clustering ----

compute_pca(df = cytoplasm_df, matrix_type = "cytoplasm")
compute_pca(df = nucleus_df, matrix_type = "nucleus")

compute_plsda(df = cytoplasm_df, matrix_type = "cytoplasm")
compute_plsda(df = nucleus_df, matrix_type = "nucleus")

create_heatmap(matrix = cytoplasm_matrix, matrix_type = "cytoplasm")
create_heatmap(matrix = nucleus_matrix, matrix_type = "nucleus")

### 3/ Differential analysis (biological conditions vs healthy controls) -----

cytoplasm_diff_analysis <- compute_diff_analysis(df = cytoplasm_matrix, matrix_type = "cytoplasm")
nucleus_diff_analysis <- compute_diff_analysis(df = nucleus_matrix, matrix_type = "nucleus")

### 4/ Data visualization -----

# a/ Volcano plots----
# Import limma top tables to plot each conditions
toptables_path <- list.files(
  path = "data/01_limma_toptables/",
  pattern = "*.xlsx",
  full.names = TRUE
)

toptables <- lapply(toptables_path, import_toptables) %>%
  setNames(tools::file_path_sans_ext(basename(toptables_path)))

# Set significance cutoff for the volcano plots
logFC_volcano_cutoff <- 1
padj_volcano_cutoff <- 0.05

create_volcano_plots(toptables)

# b/ MA plots ----
create_ma_plots(toptables)