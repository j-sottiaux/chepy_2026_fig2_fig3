## Root reproducibility check
gc()
set.seed(12345)

# 1. Biological conditions management ----
## a. Proteomics biological condition setup ----
cytoATAp <- 1:10
nuclATAp <- c(1, 3:10)

ATAn <- 11:20
ACA  <- 21:30
HC   <- 31:40

cytoNS <- 41:44
nuclNS <- 41:45

proteo_sample_sets <- list(
  cytoplasm = list(
    ATAp = cytoATAp,
    ATAn = ATAn,
    ACA  = ACA,
    HC   = HC,
    NS   = cytoNS
  ),
  nucleus = list(
    ATAp = nuclATAp,
    ATAn = ATAn,
    ACA  = ACA,
    HC   = HC,
    NS   = nuclNS
  )
)

proteo_cond_levels <- c(
  "ATAp",
  "ATAn",
  "ACA",
  "HC",
  "NS"
)

proteo_cond_labels <- c(
  ATAp = "IgG ATA+",
  ATAn = "IgG ATA- / ACA-",
  ACA  = "IgG ACA+",
  HC   = "HC IgG",
  NS   = "NS"
)

proteo_cond_labels_expr <- parse(text = c(
  ATAp = "IgG^{ATA*'+'}",
  ATAn = "IgG^{ATA*'-'~'/'~ACA*'-'}",
  ACA  = "IgG^{ACA*'+'}",
  HC   = "HC~IgG",
  NS   = "NS"
))

proteo_cond_colors <- c(
  "ATAp" = "#DCA237",
  "ATAn" = "#469C76",
  "ACA"  = "#C17DA5",
  "HC"   = "#3070AD",
  "NS"   = "#787878"
)

transcripto_cond_levels <- c("ATA", "HC", "IFNa")

transcripto_cond_labels <- c(
  ATA  = "IgG ATA+",
  HC   = "HC IgG",
  IFNa = "IFNα"
)

transcripto_cond_labels_expr <- parse(text = c(
  ATA  = "ATA",
  HC   = "HC~IgG",
  IFNa = "IFNα"
))

transcripto_cond_colors <- c(
  ATA  = "#C66526",
  HC   = "#3070AD",
  IFNa = "#E2D418"
)

make_proteo_cond_vector <- function(matrix_type) {
  matrix_type <- tolower(matrix_type)
  
  if (!(matrix_type %in% names(proteo_sample_sets))) {
    stop("Error: 'matrix_type' must be either 'cytoplasm' or 'nucleus'.")
  }
  
  sample_set <- proteo_sample_sets[[matrix_type]]
  
  factor(
    rep(names(sample_set), lengths(sample_set)),
    levels = proteo_cond_levels
  )
}

make_proteo_cond_labels <- function(group) {
  factor(
    unname(proteo_cond_labels[as.character(group)]),
    levels = unname(proteo_cond_labels[proteo_cond_levels])
  )
}

check_group_matches_df <- function(df, group, matrix_type) {
  if (ncol(df) != length(group)) {
    stop(
      "Error: number of columns in df does not match number of samples in group.\n",
      "ncol(df) = ", ncol(df), "\n",
      "length(group) = ", length(group), "\n",
      "matrix_type = ", matrix_type
    )
  }
  invisible(TRUE)
}

# Condition display Shiny helpers ----
format_proteo_condition <- function(condition_key) {
  label <- unname(proteo_cond_labels[as.character(condition_key)])
  
  if (length(label) == 0 || is.na(label)) {
    as.character(condition_key)
  } else {
    label
  }
}

safe_proteo_condition <- function(condition_key) {
  format_proteo_condition(condition_key) |>
    stringr::str_replace_all("\\+", "pos") |>
    stringr::str_replace_all("-", "neg") |>
    stringr::str_replace_all("/", "_") |>
    stringr::str_replace_all("\\s+", "_") |>
    stringr::str_replace_all("[^A-Za-z0-9_]", "")
}

# 2. Data processing functions -----
import_matrix <- function(file_path) {
  raw_matrix <- fread(file_path) %>%
    return(list(raw_matrix))
}

import_toptables <- function(toptable_path) {
  toptable_list <- read_xlsx(toptable_path) %>%
    return(list(toptable_list))
}

clean_proteo_df <- function(df) {
  names(df)[names(df) == "T: T: Gene names"] <- "gene_id"

  missing_gene_ids <- which(df[, "gene_id"] == "")
  df[missing_gene_ids, "gene_id"] <- paste("missing", 1:length(missing_gene_ids), ";")

  geneID <- strsplit(df[, "gene_id"], ";")
  geneID <- sapply(geneID, "[[", 1)

  dups <- which(duplicated(geneID))
  geneID[dups] <- paste(geneID[dups], "_duplicate_", 1:length(dups), sep = "")

  row.names(df) <- geneID

  lfq_cols <- names(df)[grepl("LFQ intensity", names(df))]
  lfq_cols_ordered <- lfq_cols[order(as.numeric(gsub("\\D+", "", lfq_cols)))]

  clean_df <- df %>%
    dplyr::select(all_of(lfq_cols_ordered))

  return(clean_df)
}

extract_DESeq2_results <- function(dds_obj, condition_1, condition_2) {
  res <- results(dds_obj, contrast = c("condition", condition_1, condition_2))
  df <- as.data.frame(res@listData)
  rownames(df) <- rownames(res)
  df <- df %>%
    dplyr::select(baseMean, stat, pvalue, padj, log2FoldChange, lfcSE)

  return(df)
}

extract_DESeq2_dataviz <- function(dds_obj, condition_1, condition_2) {
  res <- results(dds_obj, contrast = c("condition", condition_1, condition_2))
  df <- as.data.frame(res@listData)
  res_shrunk <- lfcShrink(
    dds_obj,
    contrast = c("condition", condition_1, condition_2),
    type = "ashr"
  )


  res_shrunk <- data.frame(
    as.data.frame(res_shrunk),
    stat = res$stat
  )

  rownames(res_shrunk) <- rownames(res)
  res_shrunk <- res_shrunk %>%
    dplyr::select(baseMean, stat, pvalue, padj, log2FoldChange, lfcSE)

  return(res_shrunk)
}

map_gene_symbols <- function(df) {
  df$ensembl_gene_id <- rownames(df)
  annot <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = df$ensembl_gene_id,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )

  df <- df %>%
    dplyr::left_join(annot, by = c("ensembl_gene_id" = "ENSEMBL"))

  df$SYMBOL[trimws(df$SYMBOL) == ""] <- NA

  df <- df %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::rename(gene_id = SYMBOL) %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::rename(t = stat) %>%
    dplyr::rename(logFC = log2FoldChange) %>%
    dplyr::rename(ensembl_id = ensembl_gene_id) %>%
    dplyr::select(ensembl_id, gene_id, baseMean, logFC, lfcSE, t, pvalue, padj)

  return(df)
}

clean_toptables_for_gsea <- function(toptables) {
  lapply(names(toptables), function(name) {
    df <- toptables[[name]] # Extract each top table as a data frame

    df <- df %>%
      dplyr::distinct(gene_id, .keep_all = TRUE) %>%
      dplyr::select(gene_id, t)
    ranked_list <- df$t
    ranked_list <- setNames(df$t, as.character(df$gene_id))
    ranked_list <- ranked_list[order(-ranked_list)] # Ensure it remains sorted in decreasing order

    return(list("ranked_list" = ranked_list))
  })
}

clean_toptables_for_ora <- function(toptables) {
  lapply(toptables, function(df) {
    adjusted_pval_col <- if ("padj" %in% colnames(df)) {
      "padj"
    } else if ("adj.P.Val" %in% colnames(df)) {
      "adj.P.Val"
    } else {
      stop("No adjusted p-value column found: expected 'padj' or 'adj.P.Val'")
    }

    df$diffexpressed <- "no"
    df$diffexpressed[df$logFC > logFC_threshold & df[[adjusted_pval_col]] < padj_threshold] <- "up"
    df$diffexpressed[df$logFC < -logFC_threshold & df[[adjusted_pval_col]] < padj_threshold] <- "down"

    gene_vector <- df %>%
      dplyr::filter(diffexpressed != "no") %>%
      dplyr::distinct(gene_id, .keep_all = TRUE) %>%
      dplyr::pull(gene_id)

    list(gene_vector = gene_vector)
  })
}

clean_gmt <- function(file_path) {
  clean_gmt <- read.gmt(file_path)
  clean_gmt <- clean_gmt %>%
    mutate(term = str_replace_all(term, "_", " ") %>%
      str_replace_all("GOBP|KEGG MEDICUS REFERENCE|KEGG MEDICUS|REACTOME", "") %>%
      str_replace_all("MRNA", "mRNA") %>%
      str_replace_all("RRNA", "rRNA") %>%
      str_replace_all("SNNA", "snRNA") %>%
      str_replace_all("TRNA", "tRNA")) %>%
    return(list(clean_gmt))
}

clean_evapath <- function(df) {
  # Handle missing gene ids
  df <- df[df$Genes != "" & !is.na(df$Genes), ]

  # Create a new column to clean up labels
  geneID <- strsplit(df$Genes, ";")
  geneID <- sapply(geneID, "[[", 1)

  # Create the clean dataframe
  clean_df <- df %>%
    dplyr::select(Protein, Genes)

  return(clean_df)
}

clean_endopath <- function(df) {
  df$ensembl_id <- rownames(df)
  annot <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = df$ensembl_id,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )

  df <- df %>%
    dplyr::left_join(annot, by = c("ensembl_id" = "ENSEMBL"))

  df$SYMBOL[trimws(df$SYMBOL) == ""] <- NA

  df <- df %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::rename(gene_id = SYMBOL)

  df <- df[, c("ensembl_id", "gene_id")]

  return(df)
}

filter_xp_proteome <- function(ref_list) {
  lapply(names(ref_list), function(name) {
    df <- ref_list[[name]]
    df_filtered <- df[df$gene %in% proteo_filter, ]
    return(df_filtered)
  })
}

filter_xp_transcriptome <- function(ref_list) {
  lapply(names(ref_list), function(name) {
    df <- ref_list[[name]]
    df_filtered <- df[df$gene %in% transcripto_filter, ]
    return(df_filtered)
  })
}

# 3. Exploratory analysis -----
compute_pca <- function(df, matrix_type, save_path = "figures/01_exploratory_analysis") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  matrix_type <- tolower(matrix_type)
  if (!(matrix_type %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'matrix_type' must be either 'cytoplasm' or 'nucleus'.")
  }

  group <- make_proteo_cond_vector(matrix_type)
  group_label <- make_proteo_cond_labels(group)
  check_group_matches_df(df, group, matrix_type)
  
  matrix <- as.matrix(df)
  matrix_pca <- mixOmics::pca(t(matrix), 
                              ncomp = 3, 
                              center = TRUE,
                              scale = FALSE)

  explained_var <- 100 * matrix_pca$prop_expl_var$X[1:2]
  names(explained_var) <- paste0(matrix_type, "_PC", 1:2)

  coords <- as.data.frame(matrix_pca$variates$X[, 1:2, drop = FALSE])
  colnames(coords) <- c("x", "y")
  coords$group <- group
  coords$group_label <- group_label

  plot_list <- plotIndiv(
    matrix_pca,
    group = group_label,
    ind.names = FALSE,
    legend = TRUE,
    ellipse = FALSE,
    comp = 1:2,
    title = paste(matrix_type, "/", "PCA - comp. 1~2"),
    gg = TRUE
  )

  final_plot <- plot_list$graph

  return(list(
    plot = final_plot,
    data = coords,
    explained_var = explained_var,
    matrix_type = matrix_type,
    analysis_type = "pca"
  ))
}

compute_plsda <- function(df, matrix_type, save_path = "figures/01_exploratory_analysis") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  matrix_type <- tolower(matrix_type)
  if (!(matrix_type %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'matrix_type' must be either 'cytoplasm' or 'nucleus'.")
  }

  group <- make_proteo_cond_vector(matrix_type)
  group_label <- make_proteo_cond_labels(group)
  check_group_matches_df(df, group, matrix_type)

  matrix <- as.matrix(df)
  matrix_plsda <- mixOmics::plsda(t(matrix),
                                  group,
                                  ncomp = 3,
                                  scale = TRUE)

  explained_var <- NULL
  if (!is.null(matrix_plsda$prop_expl_var$X)) {
    explained_var <- 100 * matrix_plsda$prop_expl_var$X[1:2]
    names(explained_var) <- paste0(matrix_type, "_comp", 1:2)
  }

  coords <- as.data.frame(matrix_plsda$variates$X[, 1:2, drop = FALSE])
  colnames(coords) <- c("x", "y")
  coords$group <- group
  coords$group_label <- group_label

  plot_list <- plotIndiv(
    matrix_plsda,
    group = group_label,
    ind.names = FALSE,
    legend = TRUE,
    ellipse = FALSE,
    comp = 1:2,
    title = paste(matrix_type, "/", "PLS-DA - Conf. ellipses 95%"),
    gg = TRUE
  )

  final_plot <- plot_list$graph

  return(list(
    plot = final_plot,
    data = coords,
    explained_var = explained_var,
    matrix_type = matrix_type,
    analysis_type = "plsda"
  ))
}

# 4. Differential analysis functions -----
compute_diff_analysis <- function(df, matrix_type, save_path = "data/01_proteo_toptables/") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  matrix_type <- tolower(matrix_type)
  if (!(matrix_type %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'matrix_type' must be either 'cytoplasm' or 'nucleus'.")
  }

  group <- make_proteo_cond_vector(matrix_type)
  check_group_matches_df(df, group, matrix_type)

  group <- as.factor(group)
  design <- model.matrix(~ -1 + group)
  colnames(design) <- levels(group)

  contrast_matrix <- limma::makeContrasts(
    ATAp_HC = ATAp - HC,
    ATAn_HC = ATAn - HC,
    ACA_HC = ACA - HC,
    NS_HC = NS - HC,
    levels = design
  )

  fit <- limma::lmFit(df, design)
  fit_contrast <- limma::contrasts.fit(fit, contrast_matrix)
  fit_ebayes <- limma::eBayes(fit_contrast)

  comparisons <- colnames(contrast_matrix)
  top_tables <- list()

  for (comp in comparisons) {
    tt <- topTable(fit_ebayes, coef = comp, number = Inf, adjust.method = "BH", sort.by = "none")

    # Only apply ashr if t and logFC exist and are finite
    valid <- is.finite(tt$t) & is.finite(tt$logFC) & tt$t != 0

    if (any(valid)) {
      logFC <- tt$logFC[valid]
      se_logFC <- abs(tt$logFC[valid] / tt$t[valid])

      # Adaptive shrinkage
      ashr_fit <- ashr::ash(logFC, se_logFC)
      tt$logFC_shrunk <- tt$logFC # default to original
      tt$logFC_shrunk[valid] <- get_pm(ashr_fit) 
      tt$lfsr <- NA
      tt$lfsr[valid] <- get_lfsr(ashr_fit) # Add the LFSR values from ashr
    } else {
      tt$logFC_shrunk <- tt$logFC
      tt$lfsr <- NA
    }

    # Add additional columns for the shrunken logFC
    table_to_save <- tt %>% rownames_to_column(var = "gene_id")

    # Save the top table with the added shrunken logFC
    file_name <- paste0(save_path, paste0(matrix_type, "_", comp, ".xlsx"))
    write_xlsx(table_to_save, file_name)

    top_tables[[comp]] <- tt
  }

  print(paste("All top tables were saved in:", save_path))
  return(top_tables)
}

compute_extended_diff_analysis <- function(df, matrix_type, save_path = "data/09_proteo_toptables_extended/") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  matrix_type <- tolower(matrix_type)
  if (!(matrix_type %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'matrix_type' must be either 'cytoplasm' or 'nucleus'.")
  }

  group <- make_proteo_cond_vector(matrix_type)
  check_group_matches_df(df, group, matrix_type)
  
  group <- as.factor(group)
  design <- model.matrix(~ -1 + group)
  colnames(design) <- levels(group)

  contrast_matrix <- limma::makeContrasts(
    # ATAp
    ATAp_ATAn = ATAp - ATAn,
    ATAp_ACA = ATAp - ACA,
    ATAp_HC = ATAp - HC,
    ATAp_NS = ATAp - NS,

    # ATAn
    ATAn_ATAp = ATAn - ATAp,
    ATAn_ACA = ATAn - ACA,
    ATAn_HC = ATAn - HC,
    ATAn_NS = ATAn - NS,

    # ACA
    ACA_ATAp = ACA - ATAp,
    lcSSc_ACA_dcSSc_ATAn = ACA - ATAn,
    ACA_HC = ACA - HC,
    ACA_NS = ACA - NS,

    # HC
    HC_ATAp = HC - ATAp,
    HC_ATAn = HC - ATAn,
    HC_ACA = HC - ACA,
    HC_NS = HC - NS,

    # NS
    NS_ATAp = NS - ATAp,
    NS_ATAn = NS - ATAn,
    NS_ACA = NS - ACA,
    NS_HC = NS - HC,

    # levels
    levels = design
  )

  fit <- limma::lmFit(df, design)
  fit_contrast <- limma::contrasts.fit(fit, contrast_matrix)
  fit_ebayes <- limma::eBayes(fit_contrast)

  comparisons <- colnames(contrast_matrix)

  top_tables <- list()

  for (comp in comparisons) {
    tt <- limma::topTable(fit_ebayes, coef = comp, number = Inf, adjust.method = "BH", sort.by = "none")

    # Only apply ashr if t and logFC exist and are finite
    valid <- is.finite(tt$t) & is.finite(tt$logFC) & tt$t != 0

    if (any(valid)) {
      logFC <- tt$logFC[valid]
      se_logFC <- abs(tt$logFC[valid] / tt$t[valid])

      # Adaptive shrinkage
      ashr_fit <- ashr::ash(logFC, se_logFC)
      tt$logFC_shrunk <- tt$logFC # default to original
      tt$logFC_shrunk[valid] <- get_pm(ashr_fit) # Add the shrunken logFC values
      tt$lfsr <- NA
      tt$lfsr[valid] <- get_lfsr(ashr_fit) # Add the LFSR values from ashr
    } else {
      tt$logFC_shrunk <- tt$logFC
      tt$lfsr <- NA
    }

    # Add additional columns for the shrunken logFC
    table_to_save <- tt %>% rownames_to_column(var = "gene_id")

    # Save the top table with the added shrunken logFC
    file_name <- paste0(save_path, paste0(matrix_type, "_", comp, ".xlsx"))
    write_xlsx(table_to_save, file_name)

    top_tables[[comp]] <- tt
  }

  print(paste("All top tables were saved in:", save_path))
  return(top_tables)
}

# 5. Enrichment functions -----
## a. pathways categorization ----
RNA_pathways_labels <- paste(
  "mRNA", "rRNA", "tRNA", "RNA POLYMERASE", "SPLICING",
  "RNA PROCESSING", "RNA METABOLIC", "RNA STABILITY",
  "RNA LOCALIZATION", "TRANSCRIPTION",
  sep = "|"
)

DNA_pathways_labels <- paste(
  "DNA", "G2 M", "G2", "G1", "M PHASE", "MITOTIC", "CELL CYCLE",
  "TELOMERE MAINTENANCE", "TELOMERE ORGANIZATION", "TELOMERASE", "TELOMERIC REGION",
  "CHROMOSOME", "CHROMATIN", "CHROMATIDS", "DOUBLE STRAND",
  sep = "|"
)

translation_pathways_labels <- paste(
  "RIBOSOME", "AMINO ACIDS", "SELENOAMINO", "AMINO ACID",
  "TRANSLATION", "PEPTIDE BIOSYNTHESIS", "CYTOPLASMIC TRANSLATION",
  "RIBONUCLEOPROTEIN", "RIBOSOME BIOGENESIS",
  sep = "|"
)

immunity_pathways_labels <- paste(
  "ISG15 ANTIVIRAL MECHANISM", "IMMUNE", "ANTIVIRAL", "SLITS", "ROBO", "P53", "INTERFERON", "IFN", "INTERLEUKIN", "IL", "NF KB",
  "B CELL", "T CELL", "BCR", "TCR", "MHC", "ANTIGEN PROCESSING", "DEFENSE RESPONSE",
  "CELL KILLING", "RESPONSE TO VIRUS", "VIRUS", "VIRAL", "INFECTION", "COVID", "INFLUENZA", "STRESS RESPONSE", "STARVATION", "CELL DEATH",
  "FIBROSIS", "HEDGEHOG", "APOPTOSIS",
  sep = "|"
)


# Define the pathway categorization logic
categorize_pathway <- function(description) {
  description <- str_trim(description)
  description <- str_squish(description)

  case_when(
    str_detect(description, regex(RNA_pathways_labels, ignore_case = TRUE)) ~ "RNA metabolism",
    str_detect(description, regex(DNA_pathways_labels, ignore_case = TRUE)) ~ "DNA metabolism",
    str_detect(description, regex(translation_pathways_labels, ignore_case = TRUE)) ~ "Translation related",
    str_detect(description, regex(immunity_pathways_labels, ignore_case = TRUE)) ~ "Immune response",
    TRUE ~ "Other"
  )
}


## b. GSEA ----
compute_gsea <- function(ref_background, gene_list, base_path = "data/02_gsea_results_raw/") {
  obj_name <- deparse(substitute(gene_list))

  if (grepl("proteo", obj_name, ignore.case = TRUE)) {
    save_path <- paste0(base_path, "01_proteomics_datasets")
  } else {
    save_path <- paste0(base_path, "02_transcriptomics_dataset")
  }

  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  gsea_res <- lapply(names(ref_background), function(db_name) {
    term2gene <- ref_background[[db_name]]

    conditions_results <- lapply(names(gene_list), function(condition) {
      genelist <- gene_list[[condition]]$ranked_list

      # Ensure genelist is a named numeric vector
      if (!is.numeric(genelist) || is.null(genelist)) {
        stop(paste("Error: genelist for", condition, "is not a valid named numeric vector"))
      }

      # Run GSEA using clusterProfiler's GSEA() function
      # See settings in part 4/ of 02_gsea_analysis.R for parameters modifications
      message(paste("GSEA running for", condition, "with", db_name))

      gsea <- GSEA(
        geneList = genelist,
        TERM2GENE = term2gene,
        exponent = 1,
        eps = eps_limit_gsea,
        nPermSimple = permutations_gsea,
        minGSSize = min_genecount_cutoff,
        maxGSSize = max_genecount_cutoff,
        pvalueCutoff = padj_enrichment_cutoff,
        pAdjustMethod = multiple_testing_correction,
        verbose = TRUE,
        seed = TRUE,
        by = "fgsea"
      )

      df <- gsea@result


      if (nrow(df) == 0) {
        message(paste("No significant terms for", db_name, condition))
        return(NULL)
      }

      file_name <- file.path(save_path, paste0("gsea_", db_name, "_", condition, ".xlsx"))
      writexl::write_xlsx(df, file_name)
      message(paste(file_name, "saved"))

      return(gsea)
    })
    names(conditions_results) <- names(gene_list)
    return(conditions_results)
  })
  names(gsea_res) <- names(ref_background)
  return(gsea_res)
}

filter_gsea_results <- function(result_list, base_path = "data/03_gsea_results_filtered/") {
  obj_name <- deparse(substitute(result_list))

  if (grepl("proteo", obj_name, ignore.case = TRUE)) {
    save_path <- paste0(base_path, "01_proteomics_datasets")
  } else {
    save_path <- paste0(base_path, "02_transcriptomics_dataset")
  }

  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  res_filtered <- lapply(names(result_list), function(db_name) {
    db_data <- result_list[[db_name]]

    condition_results <- lapply(names(db_data), function(condition_name) {
      enrich_obj <- db_data[[condition_name]]

      # Check if enrich_obj is NULL
      if (is.null(enrich_obj)) {
        message(paste("No enrichment result for", db_name, condition_name))
        return(NULL)
      }

      # Check if enrich_obj has a @result slot
      if (!"result" %in% slotNames(enrich_obj)) {
        message(paste("Invalid enrichment object for", db_name, condition_name))
        return(NULL)
      }

      df <- enrich_obj@result

      # Filter
      table_to_save <- df %>%
        dplyr::filter(p.adjust < padj_threshold)

      if (nrow(table_to_save) == 0) {
        message(paste("No significant terms for", db_name, condition_name))
        return(NULL)
      }

      file_name <- file.path(save_path, paste0("gsea_", db_name, "_", condition_name, "_filtered.xlsx"))
      writexl::write_xlsx(table_to_save, file_name)
      message(paste(file_name, "saved"))

      return(table_to_save)
    })

    names(condition_results) <- names(db_data)
    return(condition_results)
  })

  names(res_filtered) <- names(result_list)
  return(res_filtered)
}

## c. ORA ----
compute_ora <- function(xp_background, gene_list, base_path = "data/04_ora_results_raw/") {
  obj_name <- deparse(substitute(gene_list))

  if (grepl("proteo", obj_name, ignore.case = TRUE)) {
    save_path <- paste0(base_path, "01_proteomics_datasets")
  } else {
    save_path <- paste0(base_path, "02_transcriptomics_dataset")
  }

  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  ora_res <- lapply(names(xp_background), function(db_name) {
    term2gene <- xp_background[[db_name]]

    conditions_results <- lapply(names(gene_list), function(condition) {
      gene_vector <- gene_list[[condition]]$gene_vector

      # if (is.null(gene_vector) || length(gene_vector) == 0) {
      #   message(paste("Skipping ORA for", condition, "with", db_name, ": gene_vector is NULL or empty"))
      #   return(NULL)
      # }

      # Run ORA using clusterProfiler's enricher() function
      # See settings  parameters modifications
      message(paste("ORA running for", condition, "with", db_name))

      ora <- enricher(
        gene = gene_vector,
        TERM2GENE = term2gene,
        pvalueCutoff = padj_enrichment_cutoff,
        qvalueCutoff = padj_enrichment_cutoff,
        pAdjustMethod = multiple_testing_correction,
        minGSSize = min_genecount_cutoff,
        maxGSSize = max_genecount_cutoff,
        TERM2NAME = NULL,
        gson = NULL
      )

      df <- ora@result


      if (nrow(df) == 0) {
        message(paste("No significant terms for", db_name, condition))
        return(NULL)
      }

      file_name <- file.path(save_path, paste0("ora_", db_name, "_", condition, ".xlsx"))
      writexl::write_xlsx(df, file_name)
      message(paste(file_name, "saved"))

      return(ora)
    })
    names(conditions_results) <- names(gene_list)
    return(conditions_results)
  })
  names(ora_res) <- names(xp_background)
  return(ora_res)
}

filter_ora_results <- function(result_list, base_path = "data/05_ora_results_filtered/") {
  obj_name <- deparse(substitute(result_list))

  if (grepl("proteo", obj_name, ignore.case = TRUE)) {
    save_path <- paste0(base_path, "01_proteomics_datasets")
  } else {
    save_path <- paste0(base_path, "02_transcriptomics_dataset")
  }

  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  res_filtered <- lapply(names(result_list), function(db_name) {
    db_data <- result_list[[db_name]]

    cond_results <- lapply(names(db_data), function(cond_name) {
      enrich_obj <- db_data[[cond_name]]

      # Check if enrich_obj is NULL
      if (is.null(enrich_obj)) {
        message(paste("No enrichment result for", db_name, cond_name))
        return(NULL)
      }

      # Check if enrich_obj has a @result slot
      if (!"result" %in% slotNames(enrich_obj)) {
        message(paste("Invalid enrichment object for", db_name, cond_name))
        return(NULL)
      }

      df <- enrich_obj@result

      # Filter
      table_to_save <- df %>%
        dplyr::filter(p.adjust < padj_threshold) %>%
        dplyr::filter(Count >= min_genecount_cutoff)

      if (nrow(table_to_save) == 0) {
        message(paste("No significant terms for", db_name, cond_name))
        return(NULL)
      }

      file_name <- file.path(save_path, paste0(db_name, "_", cond_name, "_filtered.xlsx"))
      writexl::write_xlsx(table_to_save, file_name)
      message(paste(file_name, "saved"))

      return(table_to_save)
    })

    names(cond_results) <- names(db_data)
    return(cond_results)
  })

  names(res_filtered) <- names(result_list)
  return(res_filtered)
}

## d. Combining analyses ----
merge_enrichment_results <- function(result_list, enrich_type, base_path = "data/06_enrichment_merged/") {
  obj_name <- deparse(substitute(result_list))

  if (grepl("proteo", obj_name, ignore.case = TRUE)) {
    save_path <- paste0(base_path, "01_proteomics_datasets")
  } else {
    save_path <- paste0(base_path, "02_transcriptomics_dataset")
  }

  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  enrich_type <- tolower(enrich_type)
  if (!(enrich_type %in% c("ora", "gsea"))) {
    stop("Error: 'enrich_type' must be either 'ora' or 'gsea'.")
  }

  merged_results <- list()

  for (db_name in names(result_list)) {
    message("Merging results for database: ", db_name)

    db_results <- list()

    for (condition in names(result_list[[db_name]])) {
      enrich_obj <- result_list[[db_name]][[condition]]

      # Check if valid enrich results
      if (is.data.frame(enrich_obj) && !is.null(enrich_obj) && nrow(enrich_obj) > 0) {
        enrich_obj$source <- condition
        db_results[[condition]] <- enrich_obj
      } else {
        warning(paste("Skipping empty or invalid result for:", db_name, condition))
      }
    }

    if (length(db_results) > 0) {
      merged_df <- do.call(rbind, db_results)
      merged_results[[db_name]] <- merged_df

      # Write to Excel
      file_name <- file.path(save_path, paste0(enrich_type, "_", db_name, "_results_merged.xlsx"))
      writexl::write_xlsx(merged_df, file_name)
      message("✅ Saved: ", file_name)
    } else {
      warning(paste("No data to merge for:", db_name))
    }
  }

  return(merged_results)
}

integrate_proteo_enrich_res <- function(gsea_obj, ora_obj,
                                        db_names = names(ref_background_genes),
                                        conditions = names(proteo_ora_genelist),
                                        save_path = "data/07_enrichments_integration/01_proteomics_dataset/",
                                        save_files = TRUE) {
  if (save_files && !dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  # Function to process one pair of db_name and condition
  process_pair <- function(db_name, condition) {
    gsea_df <- gsea_obj[[db_name]][[condition]]@result %>%
      dplyr::select(Description, NES, p.adjust) %>%
      dplyr::mutate(NES = abs(NES)) %>%
      dplyr::rename(gsea_padj = p.adjust)

    ora_df <- ora_obj[[db_name]][[condition]]@result %>%
      dplyr::select(Description, FoldEnrichment, p.adjust) %>%
      dplyr::rename(ora_padj = p.adjust)

    merged_df <- full_join(gsea_df, ora_df, by = "Description") %>%
      dplyr::mutate(
        gsea_padj = if_else(is.na(gsea_padj) & ora_padj < 0.05, 1, gsea_padj),
        ora_padj = if_else(is.na(ora_padj) & gsea_padj < 0.05, 1, ora_padj),
        pathway_relation = categorize_pathway(Description)
      ) %>%
      dplyr::select(Description, pathway_relation, NES, gsea_padj, FoldEnrichment, ora_padj)

    return(merged_df)
  }

  # Iterate over all combinations
  param_grid <- expand.grid(db_name = db_names, condition = conditions, stringsAsFactors = FALSE)

  results_list <- purrr::pmap(param_grid, function(db_name, condition) {
    tryCatch(
      {
        message("Processing: ", db_name, " / ", condition)
        result <- process_pair(db_name, condition)

        if (save_files) {
          file_path <- file.path(save_path, paste0(db_name, "_", condition, "_enrich_integration.xlsx"))
          writexl::write_xlsx(result, file_path)
          message("Saved: ", file_path)
        }

        return(result)
      },
      error = function(e) {
        warning("❌ Failed: ", db_name, " / ", condition, " → ", conditionMessage(e))
        return(NULL)
      }
    )
  })

  # Name the results
  names(results_list) <- paste(param_grid$db_name, param_grid$condition, sep = "_")
  return(results_list)
}

integrate_transcripto_enrich_res <- function(gsea_obj, ora_obj,
                                             db_names = names(ref_background_genes),
                                             conditions = names(transcripto_ora_genelist),
                                             save_path = "data/07_enrichments_integration/02_transcriptomics_dataset/",
                                             save_files = TRUE) {
  if (save_files && !dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  # Function to process one pair of db_name and condition
  process_pair <- function(db_name, condition) {
    gsea_df <- gsea_obj[[db_name]][[condition]]@result %>%
      dplyr::select(Description, NES, p.adjust) %>%
      dplyr::mutate(NES = abs(NES)) %>%
      dplyr::rename(gsea_padj = p.adjust)

    ora_df <- ora_obj[[db_name]][[condition]]@result %>%
      dplyr::select(Description, FoldEnrichment, p.adjust) %>%
      dplyr::rename(ora_padj = p.adjust)

    merged_df <- full_join(gsea_df, ora_df, by = "Description") %>%
      dplyr::mutate(
        gsea_padj = if_else(is.na(gsea_padj) & ora_padj < 0.05, 1, gsea_padj),
        ora_padj = if_else(is.na(ora_padj) & gsea_padj < 0.05, 1, ora_padj),
        pathway_relation = categorize_pathway(Description)
      ) %>%
      dplyr::select(Description, pathway_relation, NES, gsea_padj, FoldEnrichment, ora_padj)

    return(merged_df)
  }

  # Iterate over all combinations
  param_grid <- expand.grid(db_name = db_names, condition = conditions, stringsAsFactors = FALSE)

  results_list <- purrr::pmap(param_grid, function(db_name, condition) {
    tryCatch(
      {
        message("Processing: ", db_name, " / ", condition)
        result <- process_pair(db_name, condition)

        if (save_files) {
          file_path <- file.path(save_path, paste0(db_name, "_", condition, "_enrich_integration.xlsx"))
          writexl::write_xlsx(result, file_path)
          message("Saved: ", file_path)
        }

        return(result)
      },
      error = function(e) {
        warning("❌ Failed: ", db_name, " / ", condition, " → ", conditionMessage(e))
        return(NULL)
      }
    )
  })

  # Name the results
  names(results_list) <- paste(param_grid$db_name, param_grid$condition, sep = "_")
  return(results_list)
}

## e. Create "Master data" for the "Master volcano plot" ----
process_codex_data <- function(
  gsea_obj, ora_obj,
  db_names = names(ref_background_genes),
  conditions = names(proteo_ora_genelist),
  padj_threshold = 0.05,
  min_genecount_cutoff = 5, # used only for an ORA significance flag (no dropping)
  save_path = "data/08_codex_data",
  save_files = TRUE
) {
  if (isTRUE(save_files) && !dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  split_genes <- function(x) {
    x <- dplyr::coalesce(as.character(x), "")
    out <- strsplit(x, "/", fixed = TRUE)
    lapply(out, function(v) if (length(v) == 1 && v == "") character(0) else v)
  }

  process_pair <- function(db_name, condition) {
    gsea_df <- gsea_obj[[db_name]][[condition]]@result %>%
      dplyr::transmute(
        Description     = stringr::str_squish(Description),
        NES             = abs(NES),
        setSize         = setSize,
        gsea_padj       = p.adjust,
        gsea_core_genes = core_enrichment
      )

    ora_df <- ora_obj[[db_name]][[condition]]@result %>%
      dplyr::transmute(
        Description     = stringr::str_squish(Description),
        FoldEnrichment  = FoldEnrichment,
        Count           = Count,
        ora_padj        = p.adjust,
        ora_core_genes  = geneID
      )

    merged_df <- dplyr::full_join(gsea_df, ora_df, by = "Description") %>%
      dplyr::mutate(
        # replicate enrich_integration behavior:
        gsea_padj = dplyr::if_else(is.na(gsea_padj) & !is.na(ora_padj) & ora_padj < padj_threshold, 1, gsea_padj),
        ora_padj = dplyr::if_else(is.na(ora_padj) & !is.na(gsea_padj) & gsea_padj < padj_threshold, 1, ora_padj),

        # optional flags (do not filter rows)
        is_sig_gsea = !is.na(gsea_padj) & gsea_padj < padj_threshold,
        is_sig_ora = !is.na(ora_padj) & ora_padj < padj_threshold &
          (is.na(Count) | Count >= min_genecount_cutoff),

        # core gene overlap
        gsea_list = split_genes(gsea_core_genes),
        ora_list = split_genes(ora_core_genes),
        shared_list = purrr::map2(gsea_list, ora_list, intersect),
        n_shared_core_genes = purrr::map_int(shared_list, length),
        shared_core_genes = purrr::map_chr(
          shared_list,
          function(v) if (length(v) == 0) NA_character_ else paste(v, collapse = "/")
        ),
        denom = pmin(
          dplyr::coalesce(setSize, Count),
          dplyr::coalesce(Count, setSize)
        ),
        pct_core_overlap = dplyr::if_else(denom > 0, 100 * n_shared_core_genes / denom, 0),
        pathway_relation = categorize_pathway(Description)
      ) %>%
      dplyr::select(-gsea_list, -ora_list, -shared_list, -denom) %>%
      dplyr::arrange(dplyr::desc(n_shared_core_genes), gsea_padj)

    merged_df
  }

  param_grid <- expand.grid(db_name = db_names, condition = conditions, stringsAsFactors = FALSE)

  results_list <- purrr::pmap(param_grid, function(db_name, condition) {
    tryCatch(
      {
        message("Processing: ", db_name, " / ", condition)
        res <- process_pair(db_name, condition)

        if (isTRUE(save_files)) {
          file_path <- file.path(save_path, paste0(db_name, "_", condition, "_enrichment_master_data.xlsx"))
          writexl::write_xlsx(res, file_path)
          message("Saved: ", file_path)
        }

        res
      },
      error = function(e) {
        warning("Failed: ", db_name, " / ", condition, " -> ", conditionMessage(e))
        NULL
      }
    )
  })

  names(results_list) <- paste(param_grid$db_name, param_grid$condition, sep = "_")
  results_list
}

# 6. Data visualization functions -----
pathways_colors <- c(
  "DNA metabolism"       = "#BC272D",  
  "RNA metabolism"       = "#0D7D87",  
  "Translation related"  = "#50AD9F",
  "Immune response"      = "#C99B38",  
  "Other"                = "grey48"  
)

pathways_breaks <- c(
  "DNA metabolism",
  "RNA metabolism",
  "Immune response",
  "Translation related",
  "Other"
)

publication_theme <- {
  ggplot2::theme(
    # General aspect of the plot
    plot.background = element_rect(fill = "#FFFFFF"),
    panel.background = element_rect(fill = "#FFFFFF"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    plot.title = element_text(color = "black", family = "Helvetica Neue", face = "bold", hjust = 0, size = rel(1.5)),
    plot.subtitle = element_text(color = "black", family = "Helvetica Neue", hjust = 0, size = rel(1)),

    # Axis titles and texts
    axis.title.x = element_text(color = "black", family = "Helvetica Neue", size = rel(0.9)),
    axis.text.x = element_text(color = "black", family = "Helvetica Neue", size = rel(0.8), angle = 0, hjust = 0.5),
    axis.title.y = element_text(color = "black", family = "Helvetica Neue", size = rel(0.9)), ,
    axis.text.y = element_text(color = "black", family = "Helvetica Neue", size = rel(0.8), angle = 0, hjust = 0.5),
    axis.line = element_line(color = "#5f5f5f", linetype = "solid", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linetype = "solid", linewidth = 0.25),
    panel.grid.major = element_line(color = "#EAEAEA", linetype = "dotted", linewidth = 0.25),

    # Legend
    legend.box.margin = margin(0.1, 0.1, 0.1, 0.1),
    legend.position = c(0.8, 0.85),
    legend.background = element_rect(fill = "#FFFFFF", color = "#5f5f5f", linewidth = 0.1, linetype = "solid"),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(color = "black", family = "Helvetica Neue", face = "bold", size = rel(1)),
    legend.text = element_text(color = "black", family = "Helvetica Neue", face = "italic", size = rel(0.9))
  )
}

plot_dim_reduction <- function(df_obj, base_path = "figures/01_exploratory_analysis/") {
  obj_name <- deparse(substitute(df_obj))
  
  if (grepl("proteo", obj_name, ignore.case = TRUE)) {
    save_path <- paste0(base_path, "01_proteomics_datasets")
    data <- df_obj$data
    matrix_type <- df_obj$matrix_type
    analysis_type <- tolower(df_obj$analysis_type)
    
    group <- make_proteo_cond_vector(matrix_type)
    data$group <- group
    
    color_values <- proteo_cond_colors
    color_breaks <- proteo_cond_levels
    color_labels <- proteo_cond_labels_expr
    
  } else {
    save_path <- paste0(base_path, "02_transcriptomics_dataset")
    data <- df_obj$data
    
    matrix_type <- if (!is.null(df_obj$matrix_type)) {
      df_obj$matrix_type
    } else {
      "transcriptomics"
    }
    
    analysis_type <- if (!is.null(df_obj$analysis_type)) {
      tolower(df_obj$analysis_type)
    } else {
      "pca"
    }
    
    color_values <- transcripto_cond_colors
    color_breaks <- transcripto_cond_levels
    color_labels <- transcripto_cond_labels_expr
  }
  
  required_cols <- c("x", "y", "group")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop(
      "df_obj$data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  if (analysis_type == "pca" && !is.null(df_obj$explained_var)) {
    x_lab <- paste0("Principal Component 1 (", round(df_obj$explained_var[1], 1), "% var. explained)")
    y_lab <- paste0("Principal Component 2 (", round(df_obj$explained_var[2], 1), "% var. explained)")
  } else if (analysis_type == "plsda" && !is.null(df_obj$explained_var)) {
    x_lab <- paste0("Component 1 (", round(df_obj$explained_var[1], 1), "%)")
    y_lab <- paste0("Component 2 (", round(df_obj$explained_var[2], 1), "%)")
  } else if (analysis_type == "pca") {
    x_lab <- "Principal Component 1"
    y_lab <- "Principal Component 2"
  } else {
    x_lab <- "Component 1"
    y_lab <- "Component 2"
  }
  
  max_abs <- max(abs(c(data$x, data$y)), na.rm = TRUE)
  
  if (!is.finite(max_abs) || max_abs == 0) {
    max_abs <- 1
  }
  
  max_abs <- ceiling(max_abs / 10) * 10
  
  if (max_abs == 0) {
    max_abs <- 10
  }
  
  lims <- c(-max_abs, max_abs)
  axis_breaks <- seq(lims[1], lims[2], by = 10)
  
  final_plot <- ggplot(data, aes(x = x, y = y, color = group)) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = color_values,
      breaks = color_breaks,
      labels = color_labels,
      name = NULL
    ) +
    labs(
      title = paste(toupper(analysis_type), "-", matrix_type),
      x = x_lab,
      y = y_lab
    ) +
    scale_x_continuous(limits = lims, breaks = axis_breaks, expand = c(0, 0)) +
    scale_y_continuous(limits = lims, breaks = axis_breaks, expand = c(0, 0)) +
    coord_equal() +
    ggplot2::theme(
      plot.background = element_rect(fill = "#FFFFFF"),
      panel.background = element_rect(fill = "#F7F7F7"),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
      plot.title = element_text(
        color = "black", family = "Helvetica Neue", face = "bold",
        hjust = 0, size = rel(2), margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        color = "black", family = "Helvetica Neue",
        hjust = 0, size = rel(1.2), margin = margin(b = 2)
      ),
      axis.title.x = element_text(color = "black", family = "Helvetica Neue", size = rel(1)),
      axis.text.x = element_text(color = "black", family = "Helvetica Neue", size = rel(0.9)),
      axis.title.y = element_text(color = "black", family = "Helvetica Neue", size = rel(1)),
      axis.text.y = element_text(color = "black", family = "Helvetica Neue", size = rel(0.9)),
      axis.line = element_line(color = "#5f5f5f", linewidth = 0.25),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      panel.grid.major = element_line(color = "#EAEAEA", linetype = "solid", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.margin = margin(1, 1, 1, 1),
      legend.box.margin = margin(1, 1, 1, 1),
      legend.background = element_rect(fill = "#FFFFFF", linewidth = 0.1),
      legend.key.size = unit(0.35, "cm"),
      legend.title = element_text(color = "black", family = "Helvetica Neue", face = "bold", size = rel(1)),
      legend.text = element_text(color = "black", family = "Helvetica Neue", face = "italic", size = rel(0.9))
    )
  
  file_name <- paste0(matrix_type, "_", analysis_type, ".tiff")
  
  ggsave(
    filename = file_name,
    plot = final_plot,
    path = save_path,
    width = 7,
    height = 7,
    units = "in",
    dpi = 600,
    device = "tiff"
  )
  
  message(paste(file_name, "saved in:", save_path))
  
  return(final_plot)
}

create_volcano_plots <- function(toptables, save_path = "figures/02_differential_analysis/01_proteomics_datasets") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  lapply(names(toptables), function(name) {
    df <- toptables[[name]] # Extract each top table as a data frame

    # Create differential expression cutoffs
    df$diffexpressed <- "no"
    df$diffexpressed[df$logFC_shrunk > logFC_threshold & df$adj.P.Val < padj_threshold] <- "up"
    df$diffexpressed[df$logFC_shrunk < -logFC_threshold & df$adj.P.Val < padj_threshold] <- "down"

    # Select the top 20 genes
    df_sorted <- df[order(-abs(df$logFC_shrunk)), ]
    df_top20 <- df_sorted$gene_id[1:20]

    # Label top genes
    df$gene_label <- ifelse(df$gene_id %in% df_top20 & df$diffexpressed != "no",
      df$gene_id,
      NA
    )

    # Generate the volcano plot
    volcano_df <- ggplot(df, aes(x = logFC_shrunk, y = -log10(adj.P.Val), color = diffexpressed, label = gene_label)) +
      geom_point(size = 0.9, alpha = 0.5) +
      geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), col = "#dd9d6b", linetype = "dashed") +
      geom_hline(yintercept = -log10(padj_threshold), col = "#dd9d6b", linetype = "dashed") +
      scale_color_manual(
        name = "Differential abundance",
        values = c("down" = "#189392", "no" = "#dcdbc8", "up" = "#c43a50"),
        labels = c("Decreased", "No significant", "Increased")
      ) +
      geom_label_repel(
        aes(label = gene_label),
        size = 2.2,
        box.padding = 0.25,
        segment.color = "#d7d7d7",
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      coord_cartesian(ylim = c(0, 20), xlim = c(-8, 8)) +
      scale_x_continuous(breaks = seq(-8, 8, 2)) +
      scale_y_continuous(breaks = seq(0, 20, 2)) +
      labs(
        title = paste(name),
        subtitle = "Differential abundance of proteins (labeled as gene names)",
        x = "Shrunk logFC (ashr)",
        y = "-log10 adjusted p-value"
      ) +
      publication_theme

    # Save the plot
    file_name <- paste0("volcano_plot_", name, ".tiff")
    ggsave(
      filename = file_name,
      plot = volcano_df,
      path = save_path,
      width = 6, height = 6,
      units = "in",
      dpi = 600,
      device = "tiff"
    )

    message(paste(file_name, "saved in:", save_path))
    return(volcano_df)
  })
}

plot_enrich_integration <- function(enrich_list, base_path = "figures/03_enrichments_integration/") {
  obj_name <- deparse(substitute(enrich_list))

  if (grepl("proteo", obj_name, ignore.case = TRUE)) {
    save_path <- paste0(base_path, "01_proteomics_datasets")
  } else {
    save_path <- paste0(base_path, "02_transcriptomics_dataset")
  }

  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  lapply(names(enrich_list), function(integration_df) {
    df <- enrich_list[[integration_df]]
    thr <- -log10(padj_threshold)

    # Clean up labels + compute metrics
    df <- df %>%
      mutate(
        logGSEA = -log10(gsea_padj + 1e-10),
        logORA = -log10(ora_padj + 1e-10),
        co_signif = (logGSEA > thr & logORA > thr),
        pathway_relation = factor(pathway_relation, levels = pathways_breaks)
      )

    label_df <- df %>%
      filter(co_signif) %>%
      arrange(gsea_padj + ora_padj) %>%
      mutate(pathway_relation = factor(pathway_relation, levels = pathways_breaks)) %>% 
      slice_head(n = 10)

    xmax_val <- max(df$logGSEA, na.rm = TRUE)
    ymax_val <- max(df$logORA, na.rm = TRUE)

    x_lim <- round((xmax_val + 3))
    y_lim <- round((ymax_val + 3))

    # Extract clean title components
    name_parts <- strsplit(integration_df, "_", fixed = TRUE)[[1]]
    clean_subtitle <- integration_df
    if (length(name_parts) == 4) {
      db <- stringr::str_to_title(name_parts[1])
      comp <- name_parts[2]
      cond <- paste(name_parts[3:4], collapse = "_")
      clean_subtitle <- paste0(toupper(db), " against ", cond, "diff. analysis results", "(", comp, ")")
    } else if (length(name_parts) == 3) {
      db <- stringr::str_to_title(name_parts[1])
      comp <- name_parts[2]
      cond <- name_parts[3]
      clean_subtitle <- paste0(toupper(db), " against ", cond, " (", comp, ")")
    }

    enrichplot <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = subset(df, !co_signif),
        ggplot2::aes(logGSEA, logORA),
        color = "grey80", alpha = 0.5, size = 0.8, show.legend = FALSE
      ) +
      ggplot2::annotate(
        "rect",
        xmin = thr, xmax = x_lim, ymin = thr, ymax = y_lim,
        fill = "#FFFFFF", color = "grey60", linetype = "dashed", linewidth = 0.2, alpha = 0.3
      ) +
      ggplot2::annotate(
        "text",
        x = x_lim, y = thr - 0.2, label = "co-significance area",
        hjust = 1, size = 2, color = "grey60", fontface = "bold"
      ) +
      ggplot2::geom_point(
        data = subset(df, co_signif),
        ggplot2::aes(logGSEA, logORA, color = pathway_relation),
        alpha = 0.70, size = 1
      ) +
      ggrepel::geom_label_repel(
        data = label_df,
        ggplot2::aes(
          x = logGSEA,
          y = logORA,
          label = Description,
          color = pathway_relation
        ),

        # placement controls
        min.segment.length = 0.1,
        force = 1,
        force_pull = 0.01,
        seed = 123,
        nudge_x = (xmax_val),
        hjust = 0.5,
        max.overlaps = Inf,
        direction = "y",

        # visual harmony
        point.padding = 0.15,
        box.padding = 0.25,
        label.padding = grid::unit(0.15, "lines"),
        label.r = grid::unit(0.1, "lines"),
        label.size = 0.25,
        label.hjust = 0,
        size = 2.25,
        fontface = "bold",
        family = "Helvetica Neue",
        segment.color = "grey60",
        segment.size = 0.20,
        segment.curvature = 0,
        segment.ncp = 1,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(values = pathways_colors, 
                                  breaks = pathways_breaks, 
                                  drop = TRUE) +
      ggplot2::coord_cartesian(xlim = c(0, x_lim), ylim = c(0, y_lim), clip = "off") +
      ggplot2::labs(
        title = "Enrichment analyses integration",
        subtitle = bquote(.(clean_subtitle) ~ "/" ~ italic("top 10 co-significant pathways labeled")),
        x = expression(-log[10]("GSEA padj")),
        y = expression(-log[10]("ORA padj")),
        color = "Pathway Category"
      ) +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "#FFFFFF"),
        panel.background = ggplot2::element_rect(fill = "#F7F7F7"),
        plot.margin = ggplot2::margin(5, 25, 5, 5, "pt"),
        plot.title = ggplot2::element_text(family = "Helvetica Neue", face = "bold", size = ggplot2::rel(1.5)),
        plot.subtitle = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(1)),
        axis.title.x = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(0.8)),
        axis.title.y = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(0.8)),
        axis.text.x = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(0.6)),
        axis.text.y = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(0.6)),
        axis.line = ggplot2::element_line(color = "grey35", linetype = "solid", linewidth = 0.25),
        panel.grid.major = ggplot2::element_line(color = "#F0F0F0", linetype = "solid", linewidth = 0.3),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top",
        legend.background = ggplot2::element_rect(fill = "white"),
        legend.key.size = grid::unit(0.5, "cm"),
        legend.title = ggplot2::element_text(family = "Helvetica Neue", face = "bold", size = ggplot2::rel(0.6)),
        legend.text = ggplot2::element_text(family = "Helvetica Neue", face = "italic", size = ggplot2::rel(0.6))
      )

    file_name <- paste0(integration_df, "_enrich_integration_plot.tiff")
    ggplot2::ggsave(
      filename = file_name,
      plot = enrichplot,
      path = save_path,
      width = 8, height = 6, units = "in",
      dpi = 600, device = "tiff"
    )
    message(sprintf("saved: %s/%s", save_path, file_name))
    enrichplot
  })
}

create_codex_volcano <- function(enrich_category,
                                 db_name,
                                 cell_compartment,
                                 condition,
                                 pathway_description,
                                 save_path = "figures/10_user_generated/") {
  toptable_name <- paste0(cell_compartment, "_", condition, "_HC")
  enrich_name <- paste0(db_name, "_", cell_compartment, "_", condition)
  condition_label <- format_proteo_condition(condition)

  # Select enrichment data
  enrich_df <- proteo_master_enrich_data[[enrich_name]]

  if (enrich_category == "gsea") {
    enrich_df <- enrich_df %>%
      filter(!is.na(gsea_padj)) %>%
      select(Description, NES, setSize, gsea_padj, gsea_core_genes, pathway_relation) %>%
      rename(core_genes = gsea_core_genes)
  } else if (enrich_category == "ora") {
    enrich_df <- enrich_df %>%
      filter(!is.na(ora_padj)) %>%
      select(Description, FoldEnrichment, Count, ora_padj, ora_core_genes, pathway_relation) %>%
      rename(core_genes = ora_core_genes)
  } else {
    enrich_df <- enrich_df %>%
      filter(n_shared_core_genes != 0) %>%
      rename(core_genes = shared_core_genes)
  }

  enrich_df <- enrich_df %>%
    filter(str_detect(Description, fixed(pathway_description, ignore_case = TRUE)))

  validate(need(nrow(enrich_df) > 0, "No enrichment results found for this pathway."))

  # Extract genes
  core_gene_set <- enrich_df$core_genes %>%
    strsplit(split = "/") %>%
    unlist() %>%
    unique()

  toptable <- proteo_toptables[[toptable_name]] %>%
    mutate(
      diffexpressed = case_when(
        logFC_shrunk > logFC_threshold & adj.P.Val < padj_threshold ~ "up-regulated",
        logFC_shrunk < -logFC_threshold & adj.P.Val < padj_threshold ~ "down-regulated",
        TRUE ~ "no significant difference"
      ),
      gene_label = ifelse(gene_id %in% core_gene_set, gene_id, NA),
      highlight = ifelse(!is.na(gene_label), "Involved", "Not involved"),
      hover_text = paste0("<b>Gene name:</b> ", gene_id, "<br><b>Status:</b> ", diffexpressed)
    )

  volcano_plot <- ggplot(toptable, aes(
    x = logFC_shrunk,
    y = -log10(adj.P.Val),
    label = gene_label,
    color = highlight,
    text = hover_text,
    customdata = gene_id
  )) +
    geom_point(data = toptable %>% filter(is.na(gene_label)), size = 1, alpha = 0.20) +
    geom_point(data = toptable %>% filter(!is.na(gene_label)), size = 2, alpha = 1) +
    geom_vline(
      xintercept = c(-logFC_threshold, logFC_threshold),
      col = "#dd9d6b", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(padj_threshold),
      col = "#dd9d6b", linetype = "dashed"
    ) +
    geom_label_repel(
      data = subset(toptable, highlight == "Involved"),
      aes(label = gene_label),
      size = 3, box.padding = 0.25, segment.color = "#d7d7d7",
      max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_color_manual(
      name = "Genes that are:",
      values = c("Involved" = "#C43A4F", "Not involved" = "#CFCFCF")
    ) +
    coord_cartesian(ylim = c(0, 20), xlim = c(-8, 8)) +
    scale_x_continuous(breaks = seq(-8, 8, 2)) +
    scale_y_continuous(breaks = seq(0, 20, 2)) +
    labs(
      title = paste0("Volcano plot / ", cell_compartment, " / ", condition, " vs. HC"),
      subtitle = paste0(
        "Highlighted pathway: ", pathway_description, "\n",
        db_name, " database (", toupper(enrich_category), ")"
      ),
      x = "Shrunk logFC (ashr)",
      y = "-log10 adjusted p-value"
    ) +
    ggplot2::theme(
      plot.background = element_rect(fill = "#FFFFFF"),
      panel.background = element_rect(fill = "#FFFFFF"),
      plot.margin = margin(5, 25, 5, 5, "pt"),
      plot.title = element_text(family = "Helvetica Neue", face = "bold", size = rel(1.5)),
      plot.subtitle = element_text(family = "Helvetica Neue", size = rel(1)),
      axis.title.x = element_text(family = "Helvetica Neue", size = rel(0.8)),
      axis.title.y = element_text(family = "Helvetica Neue", size = rel(0.8)),
      axis.text.x = element_text(family = "Helvetica Neue", size = rel(0.6)),
      axis.text.y = element_text(family = "Helvetica Neue", size = rel(0.6)),
      axis.line = element_line(color = "grey35", linetype = "solid", linewidth = 0.25),
      panel.grid.major = element_line(color = "#F0F0F0", linetype = "solid", linewidth = 0.3),
      panel.grid.minor = element_line(color = "#F0F0F0", linetype = "dotted", linewidth = 0.2),
      legend.position = "top",
      legend.background = element_rect(fill = "white"),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = element_text(family = "Helvetica Neue", face = "bold", size = rel(0.6)),
      legend.text = element_text(family = "Helvetica Neue", face = "italic", size = rel(0.6))
    )

  # Save if path is provided
  if (!is.null(save_path)) {
    file_name <- paste0(
      toupper(db_name), "_", enrich_category, "_volcano_",
      condition, "_", str_replace_all(tolower(pathway_description), "\\s+", "_"), ".png"
    )
    ggsave(file.path(save_path, file_name),
      plot = volcano_plot, width = 8, height = 7.8, units = "in", dpi = 600
    )
    message(paste("Plot saved at:", file.path(save_path, file_name)))
  }

  return(volcano_plot)
}

create_codex_enrichment <- function(enrich_category,
                                    db_name,
                                    cell_compartment,
                                    condition,
                                    pathway_description,
                                    save_path = "figures/10_user_generated/",
                                    padj_threshold = 0.05) {
  enrich_name <- paste0(db_name, "_", cell_compartment, "_", condition)
  condition_label <- format_proteo_condition(condition)
  
  df <- proteo_master_enrich_data[[enrich_name]]
  validate(need(!is.null(df), "No enrichment table found for this selection."))

  thr <- -log10(padj_threshold)

  # Clamp padj to avoid -Inf/Inf issues; keep NA as NA
  clamp_p <- function(p) dplyr::if_else(is.na(p), NA_real_, pmax(pmin(p, 1), 1e-300))

  df <- df %>%
    dplyr::mutate(
      gsea_padj_c = clamp_p(gsea_padj),
      ora_padj_c = clamp_p(ora_padj),
      logGSEA = -log10(gsea_padj_c),
      logORA = -log10(ora_padj_c),
      co_signif = (logGSEA > thr & logORA > thr)
    )

  label_df <- df %>%
    dplyr::filter(co_signif) %>%
    dplyr::arrange(gsea_padj + ora_padj) %>%
    dplyr::slice_head(n = 10)

  xmax_val <- max(df$logGSEA, na.rm = TRUE)
  ymax_val <- max(df$logORA, na.rm = TRUE)
  x_lim <- round(xmax_val + 3)
  y_lim <- round(ymax_val + 3)

  # Subtitle should be based on enrich_name (not a pathway description)
  clean_subtitle <- paste0(
    toupper(db_name),
    " against ",
    condition_label,
    " (",
    cell_compartment,
    ")"
  )

  df_nonsig <- dplyr::filter(df, !co_signif)
  df_sig <- dplyr::filter(df, co_signif)

  enrich_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = df_nonsig,
      ggplot2::aes(logGSEA, logORA),
      color = "grey50", alpha = 0.4, size = 0.9, show.legend = FALSE
    ) +
    ggplot2::annotate(
      "rect",
      xmin = thr, xmax = x_lim, ymin = thr, ymax = y_lim,
      fill = "grey90", color = "grey30", linetype = "dashed", linewidth = 0.4, alpha = 0.4
    ) +
    ggplot2::annotate(
      "text",
      x = x_lim, y = thr - 0.2, label = "co-significance area",
      hjust = 1, size = 2.5, color = "grey30", fontface = "bold"
    ) +
    # Significant: color mapped to pathway_relation (legend comes only from this layer)
    ggplot2::geom_point(
      data = df_sig,
      ggplot2::aes(logGSEA, logORA, color = pathway_relation),
      alpha = 0.90, size = 2
    ) +
    ggrepel::geom_label_repel(
      data = label_df,
      ggplot2::aes(
        x = logGSEA,
        y = logORA,
        label = Description,
        color = pathway_relation
      ),
      min.segment.length = 0.1,
      force = 1,
      force_pull = 0.01,
      seed = 123,
      nudge_x = xmax_val,
      hjust = 0.5,
      max.overlaps = Inf,
      direction = "y",
      point.padding = 0.15,
      box.padding = 0.25,
      label.padding = grid::unit(0.15, "lines"),
      label.r = grid::unit(0.1, "lines"),
      label.size = 0.25,
      label.hjust = 0,
      size = 2.5,
      fontface = "bold",
      family = "Helvetica Neue",
      segment.color = "grey60",
      segment.size = 0.20,
      segment.curvature = 0,
      segment.ncp = 1,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = pathways_colors, drop = TRUE) +
    ggplot2::coord_cartesian(xlim = c(0, x_lim), ylim = c(0, y_lim), clip = "off") +
    ggplot2::labs(
      title = "Enrichment analyses integration",
      subtitle = bquote(.(clean_subtitle) ~ "/" ~ italic("top 10 co-significant pathways labeled")),
      x = expression(-log[10]("GSEA padj")),
      y = expression(-log[10]("ORA padj")),
      color = "Pathway Category"
    ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#FFFFFF"),
      panel.background = ggplot2::element_rect(fill = "#FFFFFF"),
      plot.margin = ggplot2::margin(5, 25, 5, 5, "pt"),
      plot.title = ggplot2::element_text(family = "Helvetica Neue", face = "bold", size = ggplot2::rel(1.5)),
      plot.subtitle = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(1)),
      axis.title.x = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(0.8)),
      axis.title.y = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(0.8)),
      axis.text.x = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(0.6)),
      axis.text.y = ggplot2::element_text(family = "Helvetica Neue", size = ggplot2::rel(0.6)),
      axis.line = ggplot2::element_line(color = "grey35", linetype = "solid", linewidth = 0.25),
      panel.grid.major = ggplot2::element_line(color = "#F0F0F0", linetype = "solid", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_line(color = "#F0F0F0", linetype = "dotted", linewidth = 0.2),
      legend.position = "top",
      legend.background = ggplot2::element_rect(fill = "white"),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_text(family = "Helvetica Neue", face = "bold", size = ggplot2::rel(0.6)),
      legend.text = ggplot2::element_text(family = "Helvetica Neue", face = "italic", size = ggplot2::rel(0.6))
    )

  if (!is.null(save_path)) {
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

    file_name <- paste0(
      toupper(db_name), "_enrich_plot_", condition, "_",
      stringr::str_replace_all(tolower(pathway_description), "\\s+", "_"),
      ".png"
    )

    ggplot2::ggsave(
      filename = file.path(save_path, file_name),
      plot = enrich_plot, width = 8, height = 7.8, units = "in", dpi = 600
    )
  }

  enrich_plot
}

create_dae_volcano <- function(toptables,
                               cell_compartment,
                               condition_a,
                               condition_b,
                               save_path = NULL,
                               use_shrunk = TRUE,
                               top_n = 20) {
  # Ensure expected global cutoffs exist (clear error if not)
  if (!exists("logFC_threshold", inherits = TRUE)) {
    stop("logFC_threshold is not defined in the environment.")
  }
  if (!exists("padj_threshold", inherits = TRUE)) {
    stop("padj_threshold is not defined in the environment.")
  }

  # Define axis limits used below
  xlim <- c(-8, 8)
  ylim <- c(0, 20)

  # Build name and fetch table
  toptable_name <- paste0(cell_compartment, "_", condition_a, "_", condition_b)
  condition_a_label <- format_proteo_condition(condition_a)
  condition_b_label <- format_proteo_condition(condition_b)
  
  tt <- toptables[[toptable_name]]
  if (is.null(tt)) stop("No toptable found for: ", toptable_name)

  # Ensure gene_id exists
  if (!("gene_id" %in% names(tt))) {
    tt <- tt %>% tibble::rownames_to_column(var = "gene_id")
  }

  # Choose LFC column
  lfc_col <- if (isTRUE(use_shrunk) && "logFC_shrunk" %in% names(tt)) "logFC_shrunk" else "logFC"
  if (!(lfc_col %in% names(tt))) stop("No logFC column found (expected ", lfc_col, ").")
  if (!("adj.P.Val" %in% names(tt))) stop("adj.P.Val column not found.")

  df <- tt %>%
    dplyr::mutate(
      lfc = .data[[lfc_col]],
      padj = adj.P.Val,
      padj_clamped = pmax(pmin(padj, 1), 1e-300),
      mlog10 = -log10(padj_clamped),
      diffexpressed = dplyr::case_when(
        is.finite(lfc) & is.finite(padj) & lfc > logFC_threshold & padj < padj_threshold ~ "up",
        is.finite(lfc) & is.finite(padj) & lfc < -logFC_threshold & padj < padj_threshold ~ "down",
        TRUE ~ "no"
      ),
      hover_text = paste0(
        "<b>Gene:</b> ", gene_id,
        "<br><b>", lfc_col, ":</b> ", signif(lfc, 4),
        "<br><b>adj.P.Val:</b> ", signif(padj, 4)
      )
    ) %>%
    dplyr::filter(is.finite(lfc), is.finite(mlog10))

  # Top N by absolute effect
  df_sorted <- df %>% dplyr::arrange(dplyr::desc(abs(lfc)))
  top_n <- min(top_n, nrow(df_sorted))
  top_ids <- df_sorted$gene_id[seq_len(top_n)]

  df <- df %>%
    dplyr::mutate(
      gene_label = dplyr::if_else(gene_id %in% top_ids & diffexpressed != "no", gene_id, NA_character_)
    )

  volcano_plot <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = lfc, y = mlog10, color = diffexpressed, label = gene_label)
  ) +
    ggplot2::geom_point(size = 0.9, alpha = 0.5) +
    ggplot2::geom_vline(
      xintercept = c(-logFC_threshold, logFC_threshold),
      col = "#dd9d6b", linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_threshold),
      col = "#dd9d6b", linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(
      name = "Differential abundance",
      values = c("down" = "#189392", "no" = "#dcdbc8", "up" = "#c43a50"),
      labels = c("down" = "Decreased", "no" = "No significant", "up" = "Increased")
    ) +
    ggrepel::geom_label_repel(
      ggplot2::aes(label = gene_label),
      size = 2.2,
      box.padding = 0.25,
      segment.color = "#d7d7d7",
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(ylim = ylim, xlim = xlim) +
    ggplot2::scale_x_continuous(breaks = seq(xlim[1], xlim[2], 2)) +
    ggplot2::scale_y_continuous(breaks = seq(ylim[1], ylim[2], 2)) +
    ggplot2::labs(
      title = paste(cell_compartment, "-", condition_a_label, "vs.", condition_b_label),
      subtitle = "Differential abundance (labeled = top effects among significant)",
      x = lfc_col,
      y = "-log10 adjusted p-value"
    ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#FFFFFF"),
      panel.background = ggplot2::element_rect(fill = "#FFFFFF"),
      plot.margin = ggplot2::margin(5, 25, 5, 5, "pt"),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.5)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.title.x = ggplot2::element_text(size = ggplot2::rel(0.8)),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.8)),
      axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.6)),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.6)),
      axis.line = ggplot2::element_line(color = "grey35", linewidth = 0.25),
      panel.grid.major = ggplot2::element_line(color = "#F0F0F0", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_line(color = "#F0F0F0", linetype = "dotted", linewidth = 0.2),
      legend.position = "top",
      legend.background = ggplot2::element_rect(fill = "white"),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.6)),
      legend.text = ggplot2::element_text(face = "italic", size = ggplot2::rel(0.6))
    )

  if (!is.null(save_path)) {
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    file_name <- paste0(
      cell_compartment, "_",
      safe_proteo_condition(condition_a), "_vs_",
      safe_proteo_condition(condition_b),
      "_volcano.png"
    )
    ggplot2::ggsave(
      filename = file.path(save_path, file_name),
      plot = volcano_plot,
      width = 8, height = 7.8, units = "in", dpi = 600
    )
  }

  list(plot = volcano_plot, data = df, lfc_col = lfc_col)
}

create_heatmaps_transcripto <- function(
  toptables,
  vsd,
  id_col = "ensembl_id",
  symbol_col = "gene_id",
  top_n = 25,
  save_path = "figures/02_differential_analysis/02_transcriptomics_dataset/"
) {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  vsd_mat <- assay(vsd)
  annotation_col <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])
  annotation_col <- annotation_col[colnames(vsd_mat), , drop = FALSE]

  group_sizes <- table(annotation_col$condition)
  gaps_col <- cumsum(group_sizes)[-length(group_sizes)]

  res_list <- lapply(names(toptables), function(name) {
    df <- toptables[[name]]

    fc_col <- if ("logFC" %in% colnames(df)) {
      "logFC"
    } else if ("log2FoldChange" %in% colnames(df)) {
      "log2FoldChange"
    } else {
      stop("No fold-change column found in ", name)
    }

    if (!id_col %in% colnames(df)) {
      stop("Column '", id_col, "' not found in ", name)
    }

    if (!"padj" %in% colnames(df)) {
      stop("Column 'padj' not found in ", name)
    }

    df_sig <- df %>%
      dplyr::filter(
        !is.na(.data[[fc_col]]),
        !is.na(.data$padj),
        !is.na(.data[[id_col]]),
        abs(.data[[fc_col]]) > logFC_threshold,
        .data$padj < padj_threshold
      ) %>%
      dplyr::arrange(.data$padj, dplyr::desc(abs(.data[[fc_col]]))) %>%
      dplyr::slice_head(n = top_n)

    if (nrow(df_sig) == 0) {
      message("No significant genes for ", name)
      return(NULL)
    }

    keep_ids <- df_sig[[id_col]][df_sig[[id_col]] %in% rownames(vsd_mat)]

    if (length(keep_ids) == 0) {
      message("No matching IDs found in vsd for ", name)
      return(NULL)
    }

    df_sig <- df_sig[match(keep_ids, df_sig[[id_col]]), , drop = FALSE]
    mat <- vsd_mat[keep_ids, , drop = FALSE]
    mat <- mat[, rownames(annotation_col), drop = FALSE]

    mat <- mat[apply(mat, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

    if (nrow(mat) == 0) {
      message("No variable genes left for ", name)
      return(NULL)
    }

    if (symbol_col %in% colnames(df_sig)) {
      rownames(mat) <- make.unique(as.character(df_sig[[symbol_col]][match(rownames(mat), df_sig[[id_col]])]))
    }

    mat_scaled <- t(scale(t(mat)))

    pheatmap::pheatmap(
      mat_scaled,
      annotation_col = annotation_col,
      annotation_names_col = FALSE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      gaps_col = gaps_col,
      border_color = "grey85",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 7,
      main = paste0(name, " - top ", nrow(mat_scaled), " DE genes"),
      filename = file.path(save_path, paste0("heatmap_", name, ".png")),
      width = 6,
      height = 8
    )

    message("heatmap_", name, ".png saved in: ", save_path)
    invisible(mat_scaled)
  })

  names(res_list) <- names(toptables)
  invisible(res_list)
}

create_volcano_transcripto <- function(toptables, save_path = "figures/02_differential_analysis/02_transcriptomics_dataset") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  lapply(names(toptables), function(name) {
    df <- toptables[[name]] # Extract each top table as a data frame
    nickname <- name

    # Create differential expression cutoffs
    df$diffexpressed <- "no"
    df$diffexpressed[df$logFC > logFC_threshold & df$padj < padj_threshold] <- "up"
    df$diffexpressed[df$logFC < -logFC_threshold & df$padj < padj_threshold] <- "down"

    # Select the top 20 genes
    df_sig <- df %>%
      dplyr::filter(
        !is.na(logFC),
        !is.na(padj),
        abs(logFC) > logFC_threshold,
        padj < padj_threshold
      ) %>%
      dplyr::arrange(desc(abs(logFC)))

    df_top20 <- head(df_sig$gene_id, 20)

    # Label top genes
    df$gene_label <- ifelse(df$gene_id %in% df_top20 & df$diffexpressed != "no",
      df$gene_id,
      NA
    )

    # Generate the volcano plot
    volcano_df <- ggplot(df, aes(x = logFC, y = -log10(padj), color = diffexpressed, label = gene_label)) +
      geom_point(size = 0.9, alpha = 0.5) +
      geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), col = "#dd9d6b", linetype = "dashed") +
      geom_hline(yintercept = -log10(padj_threshold), col = "#dd9d6b", linetype = "dashed") +
      scale_color_manual(
        name = "Differential abundance",
        values = c("down" = "#189392", "no" = "#dcdbc8", "up" = "#c43a50"),
        labels = c("Decreased", "No significant", "Increased")
      ) +
      geom_label_repel(
        aes(label = gene_label),
        size = 2.2,
        box.padding = 0.25,
        segment.color = "#d7d7d7",
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      coord_cartesian(ylim = c(0, 175), xlim = c(-10, 10)) +
      scale_x_continuous(breaks = seq(-10, 10, 2)) +
      scale_y_continuous(breaks = seq(0, 175, 25)) +
      labs(
        title = paste(nickname),
        subtitle = "Differential abundance of mRNA transcripts (labeled as gene names)",
        x = "log2 Fold Change",
        y = "-log10 adjusted p-value"
      ) +
      theme(
        # General aspect of the plot
        plot.background = element_rect(fill = "#FFFFFF"),
        panel.background = element_rect(fill = "#FFFFFF"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        plot.title = element_text(color = "black", family = "Helvetica Neue", face = "bold", hjust = 0, size = rel(1.5)),
        plot.subtitle = element_text(color = "black", family = "Helvetica Neue", hjust = 0, size = rel(1)),

        # Axis titles and texts
        axis.title.x = element_text(color = "black", family = "Helvetica Neue", size = rel(0.9)),
        axis.text.x = element_text(color = "black", family = "Helvetica Neue", size = rel(0.8), angle = 0, hjust = 0.5),
        axis.title.y = element_text(color = "black", family = "Helvetica Neue", size = rel(0.9)),
        axis.text.y = element_text(color = "black", family = "Helvetica Neue", size = rel(0.8), angle = 0, hjust = 0.5),
        axis.line = element_line(color = "#5f5f5f", linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(color = "black", linetype = "solid", linewidth = 0.25),
        panel.grid.major = element_line(color = "#EAEAEA", linetype = "dotted", linewidth = 0.25),

        # Legend
        legend.box.margin = margin(0.1, 0.1, 0.1, 0.1),
        legend.position = c(0.25, 0.85),
        legend.background = element_rect(fill = "#FFFFFF", color = "#5f5f5f", linewidth = 0.1, linetype = "solid"),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(color = "black", family = "Helvetica Neue", face = "bold", size = rel(1)),
        legend.text = element_text(color = "black", family = "Helvetica Neue", face = "italic", size = rel(0.9))
      )

    # Save the plot
    file_name <- paste0("volcano_", nickname, ".png")
    ggsave(
      filename = file_name,
      plot = volcano_df,
      path = save_path,
      width = 6, height = 6,
      units = "in",
      dpi = 600
      # device = "tiff"
    )

    message(paste(file_name, "saved in:", save_path))
    return(volcano_df)
  })
}
