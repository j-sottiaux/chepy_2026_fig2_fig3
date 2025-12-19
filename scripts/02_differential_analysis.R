# Load required libraries
{
  library(data.table)
  library(readxl)
  library(writexl)

  library(tidyverse)
  library(ggrepel)
  library(ggthemes)

  library(limma)
  library(ashr)
}

# Reproducibility check
gc()
set.seed(12345)

# Source functions
source(file = "scripts/00_functions.R")

### 1/ Differential analysis (biological conditions vs healthy controls) -----
cytoplasm_diff_analysis <- compute_diff_analysis(df = cytoplasm_matrix, matrix_type = "cytoplasm")
nucleus_diff_analysis <- compute_diff_analysis(df = nucleus_matrix, matrix_type = "nucleus")

### 2/ Data visualization -----
# a/ Volcano plots----
# Import limma toptables for plotting each conditions
toptables_path <- list.files(
  path = "data/01_limma_toptables/",
  pattern = "*.xlsx",
  full.names = TRUE
)

toptables <- lapply(toptables_path, import_toptables) %>%
  setNames(tools::file_path_sans_ext(basename(toptables_path)))

# Generate plots
logFC_volcano_cutoff <- 1
padj_volcano_cutoff <- 0.05

create_volcano_plots(toptables)
