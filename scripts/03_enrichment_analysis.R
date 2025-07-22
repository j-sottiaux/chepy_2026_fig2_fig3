# Project setup  -----

# Load required libraries
{
  # for the FCS analysis
  
  library(tidyverse)
  library(writexl)
  library(data.table)
  
  library(ggplot2)
  library(ggrepel)
  
  library(clusterProfiler)
}

# Reproducibility check
gc()
set.seed(12345)

# Source functions
source(file = "scripts/01_functions.R")

# significance & filtering thresholds
padj_threshold <- 0.05
logFC_threshold <- 1

# 1/ GSEA data pre-processing -----

# Process top tables into a gene list
gsea_genelist <- clean_toptables_for_gsea(toptables)
names(gsea_genelist) <- gsub("_HC", "", tools::file_path_sans_ext(basename(toptables_path)))

# Import & process reference bg genes

# List all raw databases
raw_gmt <- list.files(
  path = "data",
  pattern = "*.gmt",
  full.names = TRUE
)

print(raw_gmt)

# Import & process reference bg genes
ref_proteome <- lapply(raw_gmt, clean_gmt)
names(ref_proteome) <- gsub("_v2024_1_hs_symbols", "", tools::file_path_sans_ext(basename(raw_gmt)))


# 2/ ORA data pre-processing -----

ora_genelist <- clean_toptables_for_ora(toptables)
names(ora_genelist) <- gsub("_HC", "", tools::file_path_sans_ext(basename(toptables_path)))


# Import & process experiment bg genes

# Import & process whole cell proteome from evapass 2022
evapass_df <- as.data.frame(raw_matrix$evapass_matrix)
evapass_df <- clean_evapass(evapass_df)

# Create the xp_filter to subset ref_proteome
xp_filter <- unique(unlist(list(
  rownames(nucleus_df),
  rownames(cytoplasm_df),
  evapass_df$Genes
)))

# Create xp bg genes by "xp-filtering" gobp & reactome db
xp_proteome <- filter_xp_proteome(ref_proteome)
names(xp_proteome) <- gsub("_v2024_1_hs_symbols", "", tools::file_path_sans_ext(basename(raw_gmt)))

# 3/ Running enrichment analyses -----

## a/ GSEA ----

# Analysis settings
padj_enrichment_cutoff <- 1
min_genecount_cutoff <- 15
max_genecount_cutoff <- 300
permutations_gsea <- 50000
multiple_testing_correction <- "BH"
eps_limit_gsea <- 1e-30

# Run GSEA for each condition vs. ref_proteome
gsea_results <- compute_gsea(ref_proteome, gsea_genelist)
gsea_results_filtered <- filter_gsea_results(gsea_results)

gsea_pathways_merged <- merge_enrichment_results(result_list = gsea_results_filtered, enrich_type = "gsea")

## b/ ORA ----
# Run ORA for each condition vs. xp_proteome
ora_results <- compute_ora(xp_proteome, ora_genelist)
ora_results_filtered <- filter_ora_results(ora_results)

ora_pathways_merged <- merge_enrichment_results(result_list = ora_results_filtered, enrich_type = "ora")

## c/ Enrichments integration ----
enrich_integration <- integrate_gsea_ora_results(gsea_results, ora_results)

## d/ Data pre-processing for "Master volcano plot" ----
master_enrich_data <- create_enrichment_master_data(gsea_results, ora_results)

# 4/ Data visualization -----

## a/ Per condition lollipop charts ----
create_gsea_lpop_charts(gsea_results_filtered)
create_ora_lpop_charts(ora_results_filtered)

## b/ Cross-conditions scatter dotplots ----

create_xgsea_dotplots(pathways_merged = gsea_pathways_merged, cell_location = "cytoplasm")
create_xgsea_dotplots(pathways_merged = gsea_pathways_merged, cell_location = "nucleus")

create_xora_dotplots(pathways_merged = ora_pathways_merged, cell_location = "cytoplasm")
create_xora_dotplots(pathways_merged = ora_pathways_merged, cell_location = "nucleus")

## c/ Enrichments integration  plots ----
plot_enrich_integration(enrich_integration)