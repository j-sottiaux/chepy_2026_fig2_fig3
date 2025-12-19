# Load required libraries
{
  library(data.table)
  library(readxl)
  library(writexl)

  library(tidyverse)
  library(ggrepel)
  library(ggthemes)

  library(mixOmics)
}

# Reproducibility check
gc()
set.seed(12345)

# Source functions
source(file = "scripts/00_functions.R")

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

# Overall roups design
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

### 2/ Exploratory analysis / PCA, PLS-DA ----
# PCA
cytoplasm_pca <- compute_pca(df = cytoplasm_df, matrix_type = "cytoplasm")
nucleus_pca <- compute_pca(df = nucleus_df, matrix_type = "nucleus")

# PLS-DA
cytoplasm_plsda <- compute_plsda(df = cytoplasm_df, matrix_type = "cytoplasm")
nucleus_plsda <- compute_plsda(df = nucleus_df, matrix_type = "nucleus")

### 3/ Dataviz ----
plot_dim_reduction(cytoplasm_pca)
plot_dim_reduction(nucleus_pca)
plot_dim_reduction(cytoplasm_plsda)
plot_dim_reduction(nucleus_plsda)
