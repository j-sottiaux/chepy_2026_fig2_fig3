# Load required libraries
{
  library(DESeq2)
  library(ashr)

  library(AnnotationDbi)
  library(org.Hs.eg.db)
}

# Reproducibility check
gc()
set.seed(12345)

# Source functions
source(file = "scripts/00_functions.R")


# 1.Import and prepare matrices for DESeq2 ----

## Raw counts matrix ---
rna_counts <- as.matrix(
  read.table(
    file = "data/00_aarhus_data/read_counts_rna.txt",
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    row.names = 1
  )
)

colnames(rna_counts) <- paste0("sample_", 1:9)

## Metadata data frame ---
rna_metadata <- data.frame(
  row.names = paste0("sample_", 1:9),
  condition = factor(rep(c("ATA", "HC", "IFNa"), each = 3))
)

## Sanity check ---
message(all(rownames(rna_metadata) == colnames(rna_counts))) # should be TRUE

# 2. Construct DESeqDataSet object ----
## Generate the object ---
aarhus_dds <- DESeqDataSetFromMatrix(
  countData = rna_counts,
  colData = rna_metadata,
  design = ~condition
)
dim(aarhus_dds)

## Low counts pre-filtering to improve modeling speed
rna_small_groupsize <- 3
rna_keep <- rowSums(counts(aarhus_dds) > 10) >= rna_small_groupsize
aarhus_dds <- aarhus_dds[rna_keep, ]

dim(aarhus_dds)


# 3. Differential expression analysis ----
aarhus_dds <- DESeq(aarhus_dds)
aarhus_diff_analysis_res <- results(aarhus_dds)

plotDispEsts(aarhus_dds) # Ok

aarhus_vsd <- vst(aarhus_dds, blind = TRUE)
plotPCA(aarhus_vsd, intgroup = "condition") # Ok, with caution due to low sample sizes

## Extract results + QC for visualizations and enrichment analysis ---
aarhus_ATA_vs_HC <- {
  res_raw <- results(aarhus_dds, contrast = c("condition", "ATA", "HC"))
  res_shrunk <- lfcShrink(
    aarhus_dds,
    contrast = c("condition", "ATA", "HC"),
    type = "ashr"
  )
  data.frame(
    as.data.frame(res_shrunk),
    stat = res_raw$stat
  )
}

ata_hc_qc <- hist(results(aarhus_dds, contrast = c("condition", "ATA", "HC"))$pvalue,
  main = "raw p-values histogram - ATA vs HC",
  xlab = "raw p-values"
) # Ok

aarhus_ATA_vs_IFNa <- {
  res_raw <- results(aarhus_dds, contrast = c("condition", "ATA", "IFNa"))
  res_shrunk <- lfcShrink(
    aarhus_dds,
    contrast = c("condition", "ATA", "IFNa"),
    type = "ashr"
  )
  data.frame(
    as.data.frame(res_shrunk),
    stat = res_raw$stat
  )
}

ata_ifna_qc <- hist(results(aarhus_dds, contrast = c("condition", "ATA", "IFNa"))$pvalue,
  main = "raw p-values histogram - ATA vs IFNa",
  xlab = "raw p-values"
) # Ok

aarhus_IFNa_vs_HC <- {
  res_raw <- results(aarhus_dds, contrast = c("condition", "IFNa", "HC"))
  res_shrunk <- lfcShrink(
    aarhus_dds,
    contrast = c("condition", "IFNa", "HC"),
    type = "ashr"
  )
  data.frame(
    as.data.frame(res_shrunk),
    stat = res_raw$stat
  )
}

ifna_hc_qc <- hist(results(aarhus_dds, contrast = c("condition", "IFNa", "HC"))$pvalue,
  main = "raw p-values histogram - IFNa vs HC",
  xlab = "raw p-values"
) # Proceed with caution

# 4. Map ensembl_id & gene_symbols for enrichment analysis ----
# Merge results toptables into a list object ---
aarhus_res <- list(
  aarhus_ATA_vs_HC = aarhus_ATA_vs_HC,
  aarhus_ATA_vs_IFNa = aarhus_ATA_vs_IFNa,
  aarhus_IFNa_vs_HC = aarhus_IFNa_vs_HC
)

# Gene Symbol mapping & visualization ---
aarhus_toptables <- lapply(aarhus_res, map_gene_symbols)
create_volcano_aarhus(aarhus_toptables)


# 5. Enrichment analysis ----
# Set genelist for GSEA & gene vector for ORA
aarhus_gsea_genelist <- clean_toptables_for_gsea(aarhus_toptables) %>%
  setNames(names(aarhus_toptables))

aarhus_ora_genelist <- clean_toptables_for_ora(aarhus_toptables) %>%
  setNames(names(aarhus_toptables))

# Run GSEA & filter results ---
aarhus_gsea_results <- compute_gsea(ref_proteome, aarhus_gsea_genelist)
aarhus_gsea_results_filtered <- filter_gsea_results(aarhus_gsea_results)
aarhus_gsea_pathways_merged <- merge_enrichment_results(result_list = aarhus_gsea_results_filtered, enrich_type = "gsea")

# Run ORA & filter results ---
aarhus_ora_results <- compute_ora(xp_proteome, aarhus_ora_genelist)
aarhus_ora_results_filtered <- filter_ora_results(aarhus_ora_results)
aarhus_ora_pathways_merged <- merge_enrichment_results(result_list = aarhus_ora_results_filtered, enrich_type = "ora")

# Merge enrichment results ---
aarhus_enrich_integration <- integrate_aarhus_gsea_ora_results(aarhus_gsea_results, aarhus_ora_results)
plot_enrich_integration(aarhus_enrich_integration)
