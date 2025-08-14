### 0/ Script setup --------------------------------------------------------------

# Load required libraries
{
  library(data.table)
  library(readxl)
  library(writexl)
  
  library(tidyverse)
  library(ggrepel)
  library(ggthemes)
  
}

# Reproducibility check
gc()
set.seed(12345)

#Gene signatures found in literature
EC_angiogenesis_gene_signature <-  c("PRND","LY6H","TNFRSF4","RGCC","CXCR4","PGF","IL32","COL4A1",
                                     "APLN","COL4A2","TP53I11","APOD","ACTB","LGALS1","SPARC","RGS3",
                                     "PXDN","LAMA4","GYPC","GPIHBP1","FSCN1","LXN","HSPG2","APLNR",
                                     "MCAM","TAGLN","VWA1","PFN1","TSPAN15","CCDC85B","CYTL1","SMTN",
                                     "IGFBP3","COL18A1","CALM1","C1orf54","ARHGAP18","RHOB","ACTG1","VSIR",
                                     "COL15A1","CTHRC1","IGFBP7","PRSS23","LDHB","MYL9","VWF","ODC1","CALD1", "ATP1B3")

EC_proliferative_gene_signature <- c("UBE2C","BIRC5","RRM2","PBK","AURKB","ASPM","NUSAP1","PCLAF",
                                     "TOP2A","TYMS","TK1","ZWINT","UBE2T","CENPF","CDK1","CENPW","PTTG1",
                                     "GGH","CDKN3","HMGB3","CKKS1B","DTYMK","PCNA","HMGN2","STMN1","SMC4",
                                     "TUBA1B","H2AFV","HMGB2","HMGB1","GAPDH","H2AFZ","TUBB","HIST1H4C","LGALS1",
                                     "NUCKS1","DUT","PFN1","COX8A","RANBP1","RHEB","MZT2B","PLP2","TUBA1C","SLC25A5",
                                     "SIVA1", "RAN", "TUBB4B","CAVIN3","UBE2S")

EC_endoMT_gene_signature <- read_xlsx("data/midstage_endoMT_gene_signature.xlsx") %>% 
  dplyr::pull(GeneID)

FB_pro_fibrotic_gene_signature <-  read_xlsx("data/zuh_fb_subpop.xlsx") %>% 
  dplyr::filter(cluster == "SFRP4/SFRP2+ Fib") %>% 
  dplyr::pull(gene)

cytoplasm_fibrosis_signature_genes <-  cytoplasm_df[rownames(cytoplasm_df) %in% FB_pro_fibrotic_gene_signature, ]
cytoplasm_fibrosis_heatmap <- pheatmap(cytoplasm_fibrosis_signature_genes,
                                           main = "Cytoplasm FB fibrotic gene signature expression levels",
                                           cluster_cols = TRUE,)

cytoplasm_proliferative_signature_genes <-  cytoplasm_df[rownames(cytoplasm_df) %in% EC_proliferative_gene_signature, ]
cytoplasm_proliferative_heatmap <- pheatmap(cytoplasm_proliferative_signature_genes,
                                       main = "Cytoplasm FB proliferative gene signature expression levels",
                                       cluster_cols = TRUE,)

cytoplasm_angiogenesis_signature_genes <-  cytoplasm_df[rownames(cytoplasm_df) %in% EC_angiogenesis_gene_signature, ]
cytoplasm_angiogenesis_heatmap <- pheatmap(cytoplasm_angiogenesis_signature_genes,
                                            main = "Cytoplasm FB angiogenesis gene signature expression levels",
                                            cluster_cols = TRUE,)

cytoplasm_endoMT_signature_genes <-  cytoplasm_df[rownames(cytoplasm_df) %in% EC_endoMT_gene_signature, ]
cytoplasm_endoMT_heatmap <- pheatmap(cytoplasm_endoMT_signature_genes,
                                           main = "Cytoplasm FB endoMT gene signature expression levels",
                                           cluster_cols = TRUE,)

  
nucleus_fibrosis_signature_genes <-  nucleus_df[rownames(nucleus_df) %in% FB_pro_fibrotic_gene_signature, ]
nucleus_fibrosis_heatmap <- pheatmap(nucleus_fibrosis_signature_genes,
                                       main = "Nucleus FB fibrotic gene signature expression levels",
                                       cluster_cols = TRUE,)

nucleus_proliferative_signature_genes <-  nucleus_df[rownames(nucleus_df) %in% EC_proliferative_gene_signature, ]
nucleus_proliferative_heatmap <- pheatmap(nucleus_proliferative_signature_genes,
                                            main = "Nucleus FB proliferative gene signature expression levels",
                                            cluster_cols = TRUE,)

nucleus_angiogenesis_signature_genes <-  nucleus_df[rownames(nucleus_df) %in% EC_angiogenesis_gene_signature, ]
nucleus_angiogenesis_heatmap <- pheatmap(nucleus_angiogenesis_signature_genes,
                                           main = "Nucleus FB angiogenesis gene signature expression levels",
                                           cluster_cols = TRUE,)

nucleus_endoMT_signature_genes <-  nucleus_df[rownames(nucleus_df) %in% EC_endoMT_gene_signature, ]
nucleus_endoMT_heatmap <- pheatmap(nucleus_endoMT_signature_genes,
                                     main = "nucleus FB endoMT gene signature expression levels",
                                     cluster_cols = TRUE,)