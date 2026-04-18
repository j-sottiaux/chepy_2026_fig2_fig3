# ----- 1. Script setup -----
## a. Load required libraries ----
### Annotation
library(AnnotationDbi)
library(org.Hs.eg.db)

### Visualization / wrangling
library(pheatmap)

### Statistics
library(DESeq2)
library(ashr)

## b. Reproducibility check ----
gc()
set.seed(12345)

## c. Source functions ----
source(file = "scripts/01_functions.R")


# ----- 2. Process transctiptomics dataset for DESeq2 -----
transcripto_counts <- as.matrix(
  read.table(
    file = "data/00_raw/transcripto_raw_counts.txt",
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    row.names = 1
  )
)

colnames(transcripto_counts) <- paste0("sample_", 1:9)

transcripto_metadata <- data.frame(
  row.names = paste0("sample_", 1:9),
  condition = factor(rep(c("ATA", "HC", "IFNa"), each = 3))
)

message(all(rownames(transcripto_metadata) == colnames(transcripto_counts))) # should be TRUE


# 2. Construct DESeqDataSet object ----
transcripto_dds <- DESeqDataSetFromMatrix(
  countData = transcripto_counts,
  colData = transcripto_metadata,
  design = ~condition
)
dim(transcripto_dds)

rna_small_groupsize <- 3
rna_keep <- rowSums(counts(transcripto_dds) > 10) >= rna_small_groupsize
transcripto_dds <- transcripto_dds[rna_keep, ] # Low counts pre-filtering to improve modeling speed
dim(transcripto_dds)


# ------ 3. Differential expression analysis -----
## a. Running analysis using DESeq 2 package ----
transcripto_dds <- DESeq(transcripto_dds)
transcripto_diff_analysis_res <- results(transcripto_dds)
transcripto_vsd <- vst(transcripto_dds, blind = TRUE)

transcripto_sample_dist <- dist(t(assay(transcripto_vsd)))
transcripto_sample_dist_matrix <- as.matrix(transcripto_sample_dist)
rownames(transcripto_sample_dist_matrix) <- paste(transcripto_vsd$condition, sep = "-")
colnames(transcripto_sample_dist_matrix) <- paste(transcripto_vsd$condition, sep = "-")

plotDispEsts(transcripto_dds) # Ok
plotPCA(transcripto_vsd, intgroup = "condition") # Ok, with caution due to low sample sizes
pheatmap(transcripto_sample_dist_matrix,
  clustering_distance_cols = transcripto_sample_dist,
  clustering_distance_rows = transcripto_sample_dist,
)

## b. Extract results + QC for visualizations and enrichment analysis ----
transcripto_ATA_vs_HC <- {
  res_raw <- results(transcripto_dds, contrast = c("condition", "ATA", "HC"))
  res_shrunk <- lfcShrink(
    transcripto_dds,
    contrast = c("condition", "ATA", "HC"),
    type = "ashr"
  )
  data.frame(
    as.data.frame(res_shrunk),
    stat = res_raw$stat
  )
}

transcripto_ATA_vs_IFNa <- {
  res_raw <- results(transcripto_dds, contrast = c("condition", "ATA", "IFNa"))
  res_shrunk <- lfcShrink(
    transcripto_dds,
    contrast = c("condition", "ATA", "IFNa"),
    type = "ashr"
  )
  data.frame(
    as.data.frame(res_shrunk),
    stat = res_raw$stat
  )
}

transcripto_IFNa_vs_HC <- {
  res_raw <- results(transcripto_dds, contrast = c("condition", "IFNa", "HC"))
  res_shrunk <- lfcShrink(
    transcripto_dds,
    contrast = c("condition", "IFNa", "HC"),
    type = "ashr"
  )
  data.frame(
    as.data.frame(res_shrunk),
    stat = res_raw$stat
  )
}

transcripto_ata_hc_qc <- hist(results(transcripto_dds, contrast = c("condition", "ATA", "HC"))$pvalue,
  main = "raw p-values histogram - ATA vs HC",
  xlab = "raw p-values"
) # Ok

transcripto_ata_ifna_qc <- hist(results(transcripto_dds, contrast = c("condition", "ATA", "IFNa"))$pvalue,
  main = "raw p-values histogram - ATA vs IFNa",
  xlab = "raw p-values"
) # Ok

transcripto_ifna_hc_qc <- hist(results(transcripto_dds, contrast = c("condition", "IFNa", "HC"))$pvalue,
  main = "raw p-values histogram - IFNa vs HC",
  xlab = "raw p-values"
) # Proceed with caution

## c. Merge results ----
transcripto_res <- list(
  transcripto_ATA_vs_HC = transcripto_ATA_vs_HC,
  transcripto_ATA_vs_IFNa = transcripto_ATA_vs_IFNa
  # transcripto_IFNa_vs_HC = transcripto_IFNa_vs_HC
)

## d. Gene symbol mapping ----
transcripto_toptables <- lapply(transcripto_res, map_gene_symbols)

## e. Visualization ----
create_volcano_transcripto(transcripto_toptables)
create_heatmaps_transcripto(
  toptables = transcripto_toptables,
  vsd = transcripto_vsd
)

# ----- 4. Gene Set Enrichment Analysis (GSEA) -----
## a. Process 'transcripto_toptables' object into ranked lists (based on 't' value) ----
transcripto_gsea_genelist <- clean_toptables_for_gsea(transcripto_toptables) %>%
  setNames(names(transcripto_toptables))

## b. Run GSEA for each condition vs. 'ref_background_genes' ----
transcripto_gsea_results <- compute_gsea(ref_background_genes, transcripto_gsea_genelist)
transcripto_gsea_results_filtered <- filter_gsea_results(transcripto_gsea_results)
transcripto_gsea_pathways_merged <- merge_enrichment_results(result_list = transcripto_gsea_results_filtered, enrich_type = "gsea")

# ----- 5. Over-Representation Analysis (ORA) -----
## a. Process reference background genes ----
transcripto_endopath <- fread("data/00_raw/transcripto_endopath_vst.csv")
rownames(transcripto_endopath) <- transcripto_endopath$V1
transcripto_endopath <- transcripto_endopath %>%
  dplyr::select(-V1)

transcripto_endopath <- clean_endopath(transcripto_endopath)

transcripto_filter <- unique(unlist(list(
  transcripto_toptables[["transcripto_ATA_vs_HC"]]$gene_id,
  transcripto_toptables[["transcripto_ATA_vs_IFNa"]]$gene_id,
  transcripto_endopath$gene_id
)))

xp_background_transcripto <- filter_xp_transcriptome(ref_background_genes)
names(xp_background_transcripto) <- gsub("_v2026_1_hs_symbols", "", tools::file_path_sans_ext(basename(raw_gmt)))

## b. Process 'transcripto_toptables' object into unranked significant vectors (based on padj_threshold) ----
transcripto_ora_genelist <- clean_toptables_for_ora(transcripto_toptables) %>%
  setNames(names(transcripto_toptables))


## c. Run ORA for each condition vs. 'xp_background_transcripto' ----
transcripto_ora_results <- compute_ora(xp_background_transcripto, transcripto_ora_genelist)
transcripto_ora_results_filtered <- filter_ora_results(transcripto_ora_results)
transcripto_ora_pathways_merged <- merge_enrichment_results(result_list = transcripto_ora_results_filtered, enrich_type = "ora")

# ----- 7. Cross-enrichment merging -----
transcripto_enrich_integration <- integrate_transcripto_enrich_res(transcripto_gsea_results, transcripto_ora_results)
plot_enrich_integration(transcripto_enrich_integration)
