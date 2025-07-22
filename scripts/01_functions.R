# 1/ Data processing functions -----
import_raw_matrix <- function(file_path) {
  raw_matrix <- fread(file_path) %>%
    return(list(raw_matrix))
}

import_toptables <- function(toptable_path) {
  toptable_list <- read_xlsx(toptable_path) %>%
    return(list(toptable_list))
}

clean_proteo_df <- function(df) {
  # Rename gene id column for ease of use
  names(df)[names(df) == "T: T: Gene names"] <- "gene_id"
  
  # Handle missing gene ids
  missing_gene_ids <- which(df[, "gene_id"] == "")
  df[missing_gene_ids, "gene_id"] <- paste("missing", 1:length(missing_gene_ids), ";")
  
  # Create a new column to clean up labels
  geneID <- strsplit(df[, "gene_id"], ";")
  geneID <- sapply(geneID, "[[", 1)
  
  # Label duplicated gene names to spot potential isoforms
  dups <- which(duplicated(geneID))
  geneID[dups] <- paste(geneID[dups], "_duplicate_", 1:length(dups), sep = "")
  
  # Clean labels and put gene ids as row names for clean df
  row.names(df) <- geneID
  
  # Select and reorder LFQ columns for subsequent analysis
  lfq_cols <- names(df)[grepl("LFQ intensity", names(df))]
  lfq_cols_ordered <- lfq_cols[order(as.numeric(gsub("\\D+", "", lfq_cols)))]
  
  # Create the clean dataframe
  clean_df <- df %>%
    dplyr::select(all_of(lfq_cols_ordered)) # Include the ens_gene_id column
  
  return(clean_df)
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
  lapply(names(toptables), function(name) {
    df <- toptables[[name]] # Extract each top table as a data frame
    
    df$diffexpressed <- "no"
    df$diffexpressed[df$logFC > logFC_threshold & df$adj.P.Val < padj_threshold] <- "up"
    df$diffexpressed[df$logFC < -logFC_threshold & df$adj.P.Val < padj_threshold] <- "down"
    
    df_filtered <- df %>%
      dplyr::filter(diffexpressed != "no") %>%
      dplyr::select(gene_id)
    
    gene_vector <- df_filtered$gene_id
    
    return(list("gene_vector" = gene_vector))
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

clean_evapass <- function(df) {
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

filter_xp_proteome <- function(ref_list) {
  lapply(names(ref_list), function(name) {
    df <- ref_list[[name]]
    df_filtered <- df[df$gene %in% xp_filter, ]
    return(df_filtered)
  })
}

# 2/ Exploratory analysis -----
compute_pca <- function(df, matrix_type, save_path = "figures/01_exploratory_analysis") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  matrix_type <- tolower(matrix_type)
  if (!(matrix_type %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'matrix_type' must be either 'cytoplasm' or 'nucleus'.")
  }
  
  if (matrix_type == "cytoplasm") {
    NS <- cytoNS
    ATAp <- cytoATAp
  } else if (matrix_type == "nucleus") {
    NS <- nuclNS
    ATAp <- nuclATAp
  }
  
  group <- factor(c(
    rep("dcSSc_ATAp", length(ATAp)),
    rep("dcSSc_ATAn", length(ATAn)),
    rep("lcSSc_ACA", length(ACA)),
    rep("HC", length(HC)),
    rep("NS", length(NS))
  ))
  
  matrix <- as.matrix(df)
  matrix_pca <- pca(t(matrix), ncomp = 3, center = TRUE, scale = TRUE)
  
  plot_list <- plotIndiv(matrix_pca,
                         group = group,
                         ind.names = FALSE,
                         legend = TRUE,
                         ellipse = TRUE,
                         comp = 1:2,
                         title = paste(matrix_type, "/", "PCA - comp. 1~2"),
                         gg = TRUE
  )
  
  pca_plot <- plot_list$graph
  
  final_plot <- pca_plot +
    theme(
      plot.background = element_rect(fill = "#F0F0F0", color = NA),
      panel.background = element_rect(fill = "#FFFFFF", color = NA),
      plot.margin = margin(15, 15, 15, 15),
      plot.title = element_text(color = "#304852", face = "bold", size = rel(1.6), hjust = 0),
      plot.subtitle = element_text(color = "#304852", size = rel(1.2), hjust = 0),
      axis.title = element_text(color = "#36595F", size = rel(1.2)),
      axis.text = element_text(color = "#36595F", size = rel(1)),
      axis.line = element_line(color = "#DADADA", linewidth = 0.25),
      axis.ticks = element_line(color = "#36595F", linewidth = 0.25),
      panel.grid.major = element_line(color = "#EAEAEA", linewidth = 0.25),
      legend.title = element_text(face = "bold", color = "#36595F"),
      legend.text = element_text(color = "#36595F"),
      strip.background = element_rect(fill = "#36595F"),
      strip.text = element_text(color = "white", face = "bold", size = rel(1))
    )
  
  # Save file
  file_name <- paste0("pca_", matrix_type, ".png")
  ggsave(
    filename = file_name,
    plot = final_plot,
    path = save_path,
    width = 7.7,
    height = 7.5,
    units = "in",
    dpi = 600
  )
  
  message(paste(file_name, "saved in:", save_path))
  return(final_plot)
}

compute_plsda <- function(df, matrix_type, save_path = "figures/01_exploratory_analysis") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  matrix_type <- tolower(matrix_type)
  if (!(matrix_type %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'matrix_type' must be either 'cytoplasm' or 'nucleus'.")
  }
  
  if (matrix_type == "cytoplasm") {
    NS <- cytoNS
    ATAp <- cytoATAp
  } else if (matrix_type == "nucleus") {
    NS <- nuclNS
    ATAp <- nuclATAp
  }
  
  group <- factor(c(
    rep("dcSSc_ATAp", length(ATAp)),
    rep("dcSSc_ATAn", length(ATAn)),
    rep("lcSSc_ACA", length(ACA)),
    rep("HC", length(HC)),
    rep("NS", length(NS))
  ))
  
  matrix <- as.matrix(df)
  matrix_plsda <- plsda(t(matrix), group, ncomp = 3, scale = TRUE)
  
  plot_list <- plotIndiv(matrix_plsda,
                         group = group,
                         ind.names = FALSE,
                         legend = TRUE,
                         ellipse = TRUE,
                         title = paste(matrix_type, "/", "PLS-DA - Conf. ellipses 95%"),
                         gg = TRUE
  )
  
  plsda_plot <- plot_list$graph
  
  final_plot <- plsda_plot +
    theme(
      plot.background = element_rect(fill = "#F0F0F0", color = NA),
      panel.background = element_rect(fill = "#FFFFFF", color = NA),
      plot.margin = margin(15, 15, 15, 15),
      plot.title = element_text(color = "#304852", face = "bold", size = rel(1.3), hjust = 0),
      plot.subtitle = element_text(color = "#304852", size = rel(1.2), hjust = 0),
      axis.title = element_text(color = "#36595F", size = rel(1.2)),
      axis.text = element_text(color = "#36595F", size = rel(1)),
      axis.line = element_line(color = "#DADADA", linewidth = 0.25),
      axis.ticks = element_line(color = "#36595F", linewidth = 0.25),
      panel.grid.major = element_line(color = "#EAEAEA", linewidth = 0.25),
      legend.title = element_text(face = "bold", color = "#36595F"),
      legend.text = element_text(color = "#36595F"),
      strip.background = element_rect(fill = "#36595F"),
      strip.text = element_text(color = "white", face = "bold", size = rel(1))
    )
  
  # Save file
  file_name <- paste0("plsda_", matrix_type, ".png")
  ggsave(
    filename = file_name,
    plot = final_plot,
    path = save_path,
    width = 7.7,
    height = 7.5,
    units = "in",
    dpi = 600
  )
  
  message(paste(file_name, "saved in:", save_path))
  return(plsda_plot)
}

# 3/ Differential analysis functions -----
compute_diff_analysis <- function(df, matrix_type, save_path = "data/01_limma_toptables/") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  matrix_type <- tolower(matrix_type)
  if (!(matrix_type %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'matrix_type' must be either 'cytoplasm' or 'nucleus'.")
  }
  
  if (matrix_type == "cytoplasm") {
    NS <- cytoNS
    ATAp <- cytoATAp
  } else if (matrix_type == "nucleus") {
    NS <- nuclNS
    ATAp <- nuclATAp
  }
  
  group <- factor(c(
    rep("dcSSc_ATAp", length(ATAp)),
    rep("dcSSc_ATAn", length(ATAn)),
    rep("lcSSc_ACA", length(ACA)),
    rep("HC", length(HC)),
    rep("NS", length(NS))
  ))
  
  design <- model.matrix(~ -1 + group)
  colnames(design) <- levels(group)
  
  contrast.matrix <- makeContrasts(
    dcSSc_ATAp_HC = dcSSc_ATAp - HC,
    dcSSc_ATAn_HC = dcSSc_ATAn - HC,
    lcSSc_ACA_HC = lcSSc_ACA - HC,
    NS_HC = NS - HC,
    levels = design
  )
  
  fit <- lmFit(df, design)
  fit_contrast <- contrasts.fit(fit, contrast.matrix)
  fit_ebayes <- eBayes(fit_contrast)
  
  comparisons <- c("dcSSc_ATAp_HC", "dcSSc_ATAn_HC", "lcSSc_ACA_HC", "NS_HC")
  top_tables <- list()
  
  for (comp in comparisons) {
    tt <- topTable(fit_ebayes, coef = comp, number = Inf, adjust.method = "BH", sort.by = "none")
    
    # Only apply ashr if t and logFC exist and are finite
    valid <- is.finite(tt$t) & is.finite(tt$logFC) & tt$t != 0
    
    if (any(valid)) {
      logFC <- tt$logFC[valid]
      se_logFC <- tt$logFC[valid] / tt$t[valid]
      
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
    file_name <- paste0(save_path, matrix_type, "_", comp, ".xlsx")
    write_xlsx(table_to_save, file_name)
    
    top_tables[[comp]] <- tt
  }
  
  print(paste("All top tables were saved in:", save_path))
  return(top_tables)
}

# 4/ Enrichment functions -----

# a/ pathways categorization ----
RNA_pathways_labels <- paste(
  "mRNA", "rRNA", "tRNA", "RNA POLYMERASE", "SPLICING",
  "RNA PROCESSING", "RNA METABOLIC", "RNA STABILITY",
  "RNA LOCALIZATION", "TRANSCRIPTION",
  sep = "|"
)

DNA_pathways_labels <- paste(
  "DNA", "G2 M", "G2", "G1", "M PHASE", "MITOTIC", "CELL CYCLE",
  "TELOMERE MAINTENANCE", "TELOMERE ORGANIZATION","TELOMERASE", "TELOMERIC REGION",
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
  "IMMUNE", "SLITS", "ROBO", "P53", "INTERFERON", "IFN", "INTERLEUKIN", "IL", "NF KB",
  "B CELL", "T CELL", "BCR", "TCR", "MHC", "ANTIGEN PROCESSING", "DEFENSE RESPONSE",
  "CELL KILLING", "RESPONSE TO VIRUS", "INFLUENZA", "STRESS RESPONSE", "STARVATION", "CELL DEATH",
  "FIBROSIS", "HEDGEHOG", "APOPTOSIS",
  sep = "|"
)


# Define the pathway categorization logic
categorize_pathway <- function(description) {
  description <- str_trim(description)
  description <- str_squish(description)
  
  case_when(
    str_detect(description, regex(RNA_pathways_labels, ignore_case = TRUE)) ~ "RNA metabolism & processing",
    str_detect(description, regex(DNA_pathways_labels, ignore_case = TRUE)) ~ "DNA metabolism & processing",
    str_detect(description, regex(translation_pathways_labels, ignore_case = TRUE)) ~ "Translation & ribosome biology",
    str_detect(description, regex(immunity_pathways_labels, ignore_case = TRUE)) ~ "Immune & stress response",
    TRUE ~ "Other"
  )
}


# b/ GSEA ----
compute_gsea <- function(ref_background, gene_list, save_path = "data/02_gsea_results_raw") {
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

filter_gsea_results <- function(result_list, save_path = "data/03_gsea_results_filtered") {
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

# c/ ORA ----
compute_ora <- function(xp_background, gene_list, save_path = "data/04_ora_results_raw") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  ora_res <- lapply(names(xp_background), function(db_name) {
    term2gene <- xp_background[[db_name]]
    
    conditions_results <- lapply(names(gene_list), function(condition) {
      gene_vector <- gene_list[[condition]]$gene_vector
      
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

filter_ora_results <- function(result_list, save_path = "data/05_ora_results_filtered") {
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


# d/ combining analyses ----
merge_enrichment_results <- function(result_list, enrich_type, save_path = "data/06_enrichment_results_merged") {
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
      file_name <- file.path(save_path, paste0(enrich_type, db_name, "_results_merged.xlsx"))
      writexl::write_xlsx(merged_df, file_name)
      message("✅ Saved: ", file_name)
    } else {
      warning(paste("No data to merge for:", db_name))
    }
  }
  
  return(merged_results)
}

integrate_gsea_ora_results <- function(gsea_obj, ora_obj,
                                       db_names = names(ref_proteome),
                                       conditions = names(ora_genelist),
                                       save_path = "data/07_enrichments_integration",
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

# e/ Create "Master data" for the "Master volcano plot" ----
create_enrichment_master_data <- function(gsea_obj, ora_obj,
                                          db_names = names(ref_proteome),
                                          conditions = names(ora_genelist),
                                          save_path = "data/08_master_enrichments",
                                          save_files = TRUE) {
  if (save_files && !dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  # Function to process one pair of db_name and condition
  process_pair <- function(db_name, condition) {
    gsea_df <- gsea_obj[[db_name]][[condition]]@result %>%
      dplyr::filter(p.adjust < padj_threshold) %>%
      dplyr::select(Description, NES, setSize, p.adjust, core_enrichment) %>%
      dplyr::mutate(NES = abs(NES)) %>%
      dplyr::rename(gsea_padj = p.adjust, gsea_core_genes = core_enrichment)
    
    ora_df <- ora_obj[[db_name]][[condition]]@result %>%
      dplyr::filter(p.adjust < padj_threshold) %>%
      dplyr::filter(Count >= min_genecount_cutoff) %>%
      dplyr::select(Description, FoldEnrichment, Count, p.adjust, geneID) %>%
      dplyr::rename(ora_padj = p.adjust, ora_core_genes = geneID)
    
    merged_df <- full_join(gsea_df, ora_df, by = "Description") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(shared_core_genes = paste(intersect(
        unlist(strsplit(gsea_core_genes, "/")),
        unlist(strsplit(ora_core_genes, "/"))
      ), collapse = "/")) %>%
      dplyr::mutate(n_shared_core_genes = length(intersect(
        unlist(strsplit(gsea_core_genes, "/")),
        unlist(strsplit(ora_core_genes, "/"))
      ))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(shared_core_genes = na_if(shared_core_genes, "")) %>%
      dplyr::mutate(pct_core_overlap = 100 * n_shared_core_genes / pmin(setSize, Count)) %>%
      dplyr::arrange(desc(n_shared_core_genes), gsea_padj) %>%
      dplyr::mutate(pathway_relation = categorize_pathway(Description))
  }
  
  # Iterate over all combinations
  param_grid <- expand.grid(db_name = db_names, condition = conditions, stringsAsFactors = FALSE)
  
  results_list <- purrr::pmap(param_grid, function(db_name, condition) {
    tryCatch(
      {
        message("Processing: ", db_name, " / ", condition)
        result <- process_pair(db_name, condition)
        
        if (save_files) {
          file_path <- file.path(save_path, paste0(db_name, "_", condition, "_enrichment_master_data.xlsx"))
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

# 5/ Data visualization functions -----
plot_theme <- {
  theme(
    # General aspect of the plot
    plot.background = element_rect(fill = "#F0F0F0"),
    panel.background = element_rect(fill = "#FFFFFF"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
    plot.title = element_text(color = "#1E5565", family = "Helvetica Neue", face = "bold", hjust = 0, size = rel(2)),
    plot.subtitle = element_text(color = "#304852", family = "Helvetica Neue", hjust = 0, size = rel(1.2)),
    
    # Axis titles and texts
    axis.title.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1)),
    axis.text.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(0.9), angle = 0, hjust = 0.5),
    axis.title.y = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1)), ,
    axis.text.y = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(0.9), angle = 0, hjust = 0.5),
    axis.line = element_line(color = "#DADADA", linetype = "solid", linewidth = 0.25),
    axis.ticks = element_line(color = "#36595F", linetype = "solid", linewidth = 0.25),
    panel.grid.major = element_line(color = "#F0F0F0", linetype = "solid", linewidth = 0.25),
    panel.grid.minor = element_line(color = "#F0F0F0", linetype = "dotted", linewidth = 0.15),
    
    # Legend
    legend.box.margin = margin(0.1, 0.1, 0.1, 0.1),
    legend.position = c(0.8, 0.85),
    legend.background = element_rect(fill = "#F0F0F0", color = "#b0b0b0", linewidth = 0.2, linetype = "solid"),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(color = "#36595F", family = "Helvetica Neue", face = "bold", size = rel(1)),
    legend.text = element_text(color = "#36595F", family = "Helvetica Neue", face = "italic", size = rel(0.9))
  )
}

create_heatmap <- function(matrix, matrix_type, save_path = "figures/01_exploratory_analysis") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  matrix_type <- tolower(matrix_type)
  if (!(matrix_type %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'matrix_type' must be either 'cytoplasm' or 'nucleus'.")
  }
  
  if (matrix_type == "cytoplasm") {
    NS <- cytoNS
    ATAp <- cytoATAp
  } else if (matrix_type == "nucleus") {
    NS <- nuclNS
    ATAp <- nuclATAp
  }
  
  # Build group labels
  group <- c(
    rep("dcSSc_ATAp", length(ATAp)),
    rep("dcSSc_ATAn", length(ATAn)),
    rep("lcSSc_ACA", length(ACA)),
    rep("HC", length(HC)),
    rep("NS", length(NS))
  )
  
  # Create unique names per sample
  group_labels <- unlist(lapply(split(seq_along(group), group), function(indices) {
    group_name <- group[indices[1]]
    paste0(group_name, "_", seq_along(indices))
  }))
  
  # Assign new colnames to matrix
  colnames(matrix) <- group_labels
  
  # Create heatmap
  p <- pheatmap(matrix,
                fontsize = 5,
                angle_col = 45,
                color = viridis::viridis(100),
                cutree_cols = 3,
                cutree_rows = 2,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                clustering_distance_cols = "euclidean",
                clustering_method = "ward.D2",
                show_rownames = FALSE,
                show_colnames = TRUE,
                main = paste("Samples clustering Heatmap for cellular compartment :", matrix_type)
  )
  
  # Save heatmap
  ggsave(
    filename = paste0("heatmap_", matrix_type, ".png"),
    plot = p$gtable,
    path = save_path,
    width = 9,
    height = 8,
    dpi = 300
  )
  
  message("Heatmap saved at: ", file.path(save_path, paste0("heatmap_", matrix_type, ".png")))
}

create_volcano_plots <- function(toptables, save_path = "figures/02_differential_analysis") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  lapply(names(toptables), function(name) {
    df <- toptables[[name]] # Extract each top table as a data frame
    
    # Create differential expression cutoffs
    df$diffexpressed <- "no"
    df$diffexpressed[df$logFC_shrunk > logFC_volcano_cutoff & df$adj.P.Val < padj_volcano_cutoff] <- "up"
    df$diffexpressed[df$logFC_shrunk < -logFC_volcano_cutoff & df$adj.P.Val < padj_volcano_cutoff] <- "down"
    
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
      geom_vline(xintercept = c(-logFC_volcano_cutoff, logFC_volcano_cutoff), col = "#dd9d6b", linetype = "dashed") +
      geom_hline(yintercept = -log10(padj_volcano_cutoff), col = "#dd9d6b", linetype = "dashed") +
      scale_color_manual(
        name = "Differential abundance",
        values = c("down" = "#189392", "no" = "#dcdbc8", "up" = "#c43a50"),
        labels = c("Decreased", "No significant", "Increased")
      ) +
      geom_label_repel(
        aes(label = gene_label),
        size = 3,
        box.padding = 0.25,
        segment.color = "#d7d7d7",
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      coord_cartesian(ylim = c(0, 20), xlim = c(-8, 8)) +
      scale_x_continuous(breaks = seq(-8, 8, 2)) +
      scale_y_continuous(breaks = seq(0, 20, 2)) +
      labs(
        title = paste("Volcano plot /", name),
        subtitle = "Differential abundance of proteins (labeled as gene names)",
        x = "Shrunk logFC (ashr)",
        y = "-log10 adjusted p-value"
      ) +
      plot_theme
    
    # Save the plot
    file_name <- paste0("volcano_plot_", name, ".png")
    ggsave(
      filename = file_name,
      plot = volcano_df,
      path = save_path,
      width = 7.7, height = 7.5,
      units = "in",
      dpi = 600
    )
    
    message(paste(file_name, "saved in:", save_path))
    return(volcano_df)
  })
}

create_ma_plots <- function(toptables, save_path = "figures/02_differential_analysis") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  lapply(names(toptables), function(name) {
    df <- toptables[[name]] # Extract each top table as a data frame
    
    # Create significance cutoff
    df$significant <- ifelse(df$adj.P.Val < padj_volcano_cutoff, "yes", "no")
    
    
    # Create differential expression cutoffs
    df$diffexpressed <- "no"
    df$diffexpressed[df$logFC_shrunk > logFC_volcano_cutoff & df$adj.P.Val < padj_volcano_cutoff] <- "up"
    df$diffexpressed[df$logFC_shrunk < -logFC_volcano_cutoff & df$adj.P.Val < padj_volcano_cutoff] <- "down"
    
    
    # Select the top 20 genes
    df_sorted <- df[order(-abs(df$logFC_shrunk)), ]
    df_top20 <- df_sorted$gene_id[1:20]
    
    # Label top genes
    df$gene_label <- ifelse(df$gene_id %in% df_top20 & df$diffexpressed != "no",
                            df$gene_id,
                            NA
    )
    
    # Generate the ma plot
    ma_df <- ggplot(df, aes(x = log(AveExpr), y = logFC_shrunk, color = significant, label = gene_label)) +
      geom_point(size = 0.9, alpha = 0.4) +
      geom_hline(yintercept = c(-logFC_volcano_cutoff, logFC_volcano_cutoff), col = "#dd9d6b", linetype = "dashed") +
      geom_label_repel(
        aes(label = gene_label),
        size = 3,
        box.padding = 0.25,
        segment.color = "#d7d7d7",
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      scale_color_manual(
        name = "Significant change",
        values = c("no" = "#dcdbc8", "yes" = "#dd9d6b"),
        labels = c("No", "Yes")
      ) +
      coord_cartesian(ylim = c(-8, 8), range(log(df$AveExpr), na.rm = TRUE)) +
      scale_x_continuous(breaks = seq(3, 4, 0.2)) +
      scale_y_continuous(breaks = seq(-8, 8, 2)) +
      labs(
        title = paste("MA plot /", name),
        subtitle = "Differential abundance of proteins (labeled as gene name)",
        x = "log Average Expression",
        y = "Shrunk log FC (ashr)"
      ) +
      plot_theme
    
    # Save the plot
    file_name <- paste0("ma_plot_", name, ".png")
    ggsave(
      filename = file_name,
      plot = ma_df,
      path = save_path,
      width = 7.7,
      height = 7.5,
      units = "in",
      dpi = 600
    )
    
    message(paste(file_name, "saved in:", save_path))
    return(ma_df)
  })
}

create_gsea_lpop_charts <- function(gsea_results, save_path = "figures/03_gsea_lollipop_charts") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  lapply(names(gsea_results), function(db_name) {
    db_data <- gsea_results[[db_name]]
    
    lapply(names(db_data), function(condition_name) {
      df <- db_data[[condition_name]]
      
      if (is.null(df) || nrow(df) == 0) {
        message(paste("Skipping", db_name, condition_name, "- empty dataframe"))
        return(NULL)
      }
      
      # Select top 30 pathways (balanced positive & negative)
      top_pos <- df %>%
        filter(NES > 0) %>%
        arrange(desc(NES)) %>%
        head(30)
      
      top_neg <- df %>%
        filter(NES < 0) %>%
        arrange(NES) %>%
        head(30)
      
      combined <- if (nrow(top_pos) + nrow(top_neg) > 30) {
        bind_rows(top_pos, top_neg) %>%
          arrange(desc(abs(NES))) %>%
          head(30)
      } else {
        bind_rows(top_pos, top_neg)
      }
      
      if (nrow(combined) == 0) {
        message(paste("Skipping", db_name, condition_name, "- no enriched pathways found"))
        return(NULL)
      }
      
      n_display <- nrow(combined)
      n_total <- nrow(df)
      
      # Plot
      lollipop_df <- ggplot(combined, aes(x = NES, y = fct_reorder(Description, NES))) +
        geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description),
                     color = "#36595F", size = 0.3, linetype = "dashed"
        ) +
        geom_point(aes(color = p.adjust, size = setSize)) +
        geom_vline(xintercept = 0, linetype = "solid", color = "#36595F", linewidth = 0.3) +
        scale_color_gradientn(
          colours = c("#604838", "#EB8A3E", "#EBB582"),
          guide = guide_colorbar(reverse = FALSE)
        ) +
        scale_size_continuous(range = c(2, 8)) +
        scale_x_continuous(expand = expansion(mult = 0.05)) +
        labs(
          title = paste("Top", n_display, "out of", n_total, "enriched pathways /", condition_name, "vs.", db_name),
          subtitle = "GSEA conducted with : padj cutoff = 0.05 (BH corr.) | 15 < Genes setSize < 300 | Permutations = 50k",
          x = "Normalized Enrichment Score • NES",
          y = NULL
        ) +
        theme(
          plot.background = element_rect(fill = "#F0F0F0"),
          panel.background = element_rect(fill = "#FFFFFF"),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          plot.title = element_text(color = "#304852", family = "Helvetica Neue", face = "bold", hjust = 0, size = rel(1.6)),
          plot.subtitle = ggtext::element_markdown(color = "#304852", family = "Helvetica Neue", hjust = 0, size = rel(1.2)),
          axis.title.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1)),
          axis.text.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(0.9), angle = 0, hjust = 0.5),
          axis.title.y = ggtext::element_markdown(family = "Helvetica Neue", size = rel(1)),
          axis.text.y = element_text(color = "#36595F", family = "Helvetica Neue", face = "italic", size = rel(0.9), angle = 0, hjust = 1),
          axis.line = element_line(color = "#DADADA", linetype = "solid", linewidth = 0.25),
          axis.ticks = element_line(color = "#36595F", linetype = "solid", linewidth = 0.25),
          panel.grid.major = element_line(color = "#F0F0F0", linetype = "solid", linewidth = 0.25),
        )
      
      file_name <- paste0("gsea_lpop_chart_", db_name, "_", condition_name, ".png")
      ggsave(
        filename = file_name,
        plot = lollipop_df,
        path = save_path,
        scale = 0.8,
        width = 18,
        height = 10,
        units = "in",
        dpi = 600
      )
      
      message(paste(file_name, "saved in:", save_path))
      return(lollipop_df)
    })
  })
}

create_ora_lpop_charts <- function(result_list, save_path = "figures/04_ora_lollipop_charts") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  lapply(names(result_list), function(db_name) {
    db_data <- result_list[[db_name]]
    
    lapply(names(db_data), function(condition_name) {
      df <- db_data[[condition_name]]
      
      if (is.null(df) || nrow(df) == 0) {
        message(paste("Skipping", db_name, condition_name, "- empty dataframe"))
        return(NULL)
      }
      
      # Select top 30 pathways (balanced positive & negative)
      top_pos <- df %>%
        filter(FoldEnrichment > 0) %>%
        arrange(desc(FoldEnrichment)) %>%
        head(30)
      
      top_neg <- df %>%
        filter(FoldEnrichment < 0) %>%
        arrange(FoldEnrichment) %>%
        head(30)
      
      combined <- if (nrow(top_pos) + nrow(top_neg) > 30) {
        bind_rows(top_pos, top_neg) %>%
          arrange(desc(abs(FoldEnrichment))) %>%
          head(30)
      } else {
        bind_rows(top_pos, top_neg)
      }
      
      if (nrow(combined) == 0) {
        message(paste("Skipping", db_name, condition_name, "- no enriched pathways found"))
        return(NULL)
      }
      
      n_display <- nrow(combined)
      n_total <- nrow(df)
      
      # Plot
      lollipop_df <- ggplot(combined, aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment))) +
        geom_segment(aes(x = 0, xend = FoldEnrichment, y = Description, yend = Description),
                     color = "#36595F", size = 0.3, linetype = "dashed"
        ) +
        geom_point(aes(color = p.adjust, size = Count)) +
        geom_vline(xintercept = 0, linetype = "solid", color = "#36595F", linewidth = 0.3) +
        scale_color_gradientn(
          colours = c("#604838", "#EB8A3E", "#EBB582"),
          guide = guide_colorbar(reverse = FALSE)
        ) +
        scale_size_continuous(range = c(2, 8)) +
        scale_x_continuous(expand = expansion(mult = 0.05)) +
        labs(
          title = paste("Top", n_display, "out of", n_total, "enriched pathways /", condition_name, "vs.", db_name),
          subtitle = "ORA conducted with : padj cutoff = 0.05 (BH corr.) | 15 < Genes setSize < 300",
          x = "Fold Enrichment",
          y = NULL
        ) +
        theme(
          plot.background = element_rect(fill = "#F0F0F0"),
          panel.background = element_rect(fill = "#FFFFFF"),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          plot.title = element_text(color = "#304852", family = "Helvetica Neue", face = "bold", hjust = 0, size = rel(1.8)),
          plot.subtitle = element_text(color = "#304852", family = "Helvetica Neue", hjust = 0, size = rel(1.2)),
          axis.title.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1)),
          axis.text.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(0.9), angle = 0, hjust = 0.5),
          axis.title.y = element_text(family = "Helvetica Neue", size = rel(1)),
          axis.text.y = element_text(color = "#36595F", family = "Helvetica Neue", face = "italic", size = rel(0.9), angle = 0, hjust = 1),
          axis.line = element_line(color = "#DADADA", linetype = "solid", linewidth = 0.25),
          axis.ticks = element_line(color = "#36595F", linetype = "solid", linewidth = 0.25),
          panel.grid.major = element_line(color = "#F0F0F0", linetype = "solid", linewidth = 0.25),
        )
      
      file_name <- paste0("ora_lpop_chart_", db_name, "_", condition_name, ".png")
      ggsave(
        filename = file_name,
        plot = lollipop_df,
        path = save_path,
        scale = 0.8,
        width = 18,
        height = 10,
        units = "in",
        dpi = 600
      )
      
      message(paste(file_name, "saved in:", save_path))
      return(lollipop_df)
    })
  })
}

create_xgsea_dotplots <- function(pathways_merged, cell_location, save_path = "figures/05_xgsea_dotplots") {
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  if (!(cell_location %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'location' must be either 'cytoplasm' or 'nucleus'.")
  }
  
  if (cell_location == "cytoplasm") {
    compartment <- "cytoplasm"
  } else if (cell_location == "nucleus") {
    compartment <- "nucleus"
  }
  
  # Loop over databases
  lapply(names(pathways_merged), function(db_name) {
    df <- as.data.frame(pathways_merged[[db_name]])
    
    # Extract location and filter
    df <- df %>%
      mutate(
        location = case_when(
          grepl("nucleus", source, ignore.case = TRUE) ~ "nucleus",
          grepl("cytoplasm", source, ignore.case = TRUE) ~ "cytoplasm",
          TRUE ~ NA_character_
        )
      )
    
    df$location <- factor(df$location, levels = c("nucleus", "cytoplasm"))
    df <- df %>% filter(location == compartment)
    
    df$source <- stringr::str_extract(df$source, "lcSSc_ACA|dcSSc_ATAn|dcSSc_ATAp|NS")
    
    # Select top 100 pathways (balanced positive & negative)
    top_pos <- df %>%
      filter(NES > 0) %>%
      arrange(desc(NES)) %>%
      head(100)
    
    top_neg <- df %>%
      filter(NES < 0) %>%
      arrange(NES) %>%
      head(100)
    
    combined_df <- if (nrow(top_pos) + nrow(top_neg) > 100) {
      bind_rows(top_pos, top_neg) %>%
        arrange(desc(abs(NES))) %>%
        head(100)
    } else {
      bind_rows(top_pos, top_neg)
    }
    
    if (nrow(combined_df) == 0) {
      message(paste("Skipping", db_name, "- no enriched pathways found"))
      return(NULL)
    }
    
    # Determine color scale dynamically
    min_nes <- min(combined_df$NES, na.rm = TRUE)
    max_nes <- max(combined_df$NES, na.rm = TRUE)
    
    if (min_nes < 0 && max_nes > 0) {
      # Mixed NES values – diverging viridis-friendly (red-blue approximation)
      color_scale <- scale_color_gradient2(
        low = "#A43820",
        mid = "#FFEB94",
        high = "#7CAA2D",
        midpoint = 0,
        guide = guide_colorbar(title = "NES", reverse = FALSE)
      )
    } else if (max_nes <= 0) {
      # Only negative NES – red sequential
      color_scale <- scale_color_gradient(
        low = "#A43820",
        high = "#FFEB94",
        guide = guide_colorbar(title = "NES", reverse = FALSE)
      )
    } else {
      # Only positive NES – green sequential
      color_scale <- scale_color_gradient(
        low = "#FFEB94",
        high = "#7CAA2D",
        guide = guide_colorbar(title = "NES", reverse = FALSE)
      )
    }
    
    n_display <- nrow(combined_df)
    n_total <- nrow(df)
    
    # Build dotplot
    dotplot <- ggplot(combined_df, aes(x = source, y = Description, color = NES)) +
      geom_point(aes(size = setSize), alpha = 0.8) +
      facet_grid(~location) +
      color_scale +
      labs(
        title = paste("Top", n_display, "out of", n_total, "enriched pathways /", db_name),
        subtitle = "GSEA conducted with: padj < 0.05 (BH corr.) | 15 < Genes setSize < 300 | Permutations = 50k",
        x = NULL,
        y = NULL
      ) +
      theme(
        plot.background = element_rect(fill = "#F0F0F0"),
        panel.background = element_rect(fill = "#FFFFFF"),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        plot.title = element_text(color = "#304852", family = "Helvetica Neue", face = "bold", hjust = 0, size = rel(2)),
        plot.subtitle = element_text(color = "#304852", family = "Helvetica Neue", hjust = 0, size = rel(1.3)),
        axis.title.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1), hjust = 0.5),
        axis.text.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1), hjust = 0.5),
        axis.title.y = element_text(family = "Helvetica Neue", size = rel(1)),
        axis.text.y = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(0.9), face = "italic", hjust = 1, margin = margin(t = 30, b = 30, r = 5, l = 5)),
        axis.line = element_line(color = "#DADADA", linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(color = "#36595F", linetype = "solid", linewidth = 0.25),
        panel.grid.major = element_line(color = "#F0F0F0", linetype = "solid", linewidth = 0.25),
        strip.background = element_rect(fill = "#36595F"),
        strip.text = element_text(size = 12, face = "bold", color = "#FFFFFF")
      )
    
    # Save the file
    file_name <- paste0("gsea_xconditions_", cell_location, "_", db_name, ".png")
    ggsave(
      filename = file_name,
      plot = dotplot,
      path = save_path,
      scale = 0.8,
      width = 20,
      height = 13,
      units = "in",
      dpi = 600
    )
    
    message(paste(file_name, "saved in:", save_path))
    return(dotplot)
  })
}

create_xora_dotplots <- function(pathways_merged, cell_location, save_path = "figures/06_xora_dotplots") {
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  if (!(cell_location %in% c("cytoplasm", "nucleus"))) {
    stop("Error: 'location' must be either 'cytoplasm' or 'nucleus'.")
  }
  
  if (cell_location == "cytoplasm") {
    compartment <- "cytoplasm"
  } else if (cell_location == "nucleus") {
    compartment <- "nucleus"
  }
  
  # Loop over databases
  lapply(names(pathways_merged), function(db_name) {
    df <- as.data.frame(pathways_merged[[db_name]])
    
    # Extract location and filter
    df <- df %>%
      mutate(
        location = case_when(
          grepl("nucleus", source, ignore.case = TRUE) ~ "nucleus",
          grepl("cytoplasm", source, ignore.case = TRUE) ~ "cytoplasm",
          TRUE ~ NA_character_
        )
      )
    
    df$location <- factor(df$location, levels = c("nucleus", "cytoplasm"))
    df <- df %>% filter(location == compartment)
    
    df$source <- stringr::str_extract(df$source, "lcSSc_ACA|dcSSc_ATAn|dcSSc_ATAp|NS")
    
    # Select top 100 pathways (balanced positive & negative)
    top_pos <- df %>%
      filter(FoldEnrichment > 0) %>%
      arrange(desc(FoldEnrichment)) %>%
      head(100)
    
    top_neg <- df %>%
      filter(FoldEnrichment < 0) %>%
      arrange(FoldEnrichment) %>%
      head(100)
    
    combined_df <- if (nrow(top_pos) + nrow(top_neg) > 100) {
      bind_rows(top_pos, top_neg) %>%
        arrange(desc(abs(FoldEnrichment))) %>%
        head(100)
    } else {
      bind_rows(top_pos, top_neg)
    }
    
    if (nrow(combined_df) == 0) {
      message(paste("Skipping", db_name, "- no enriched pathways found"))
      return(NULL)
    }
    
    # Determine color scale dynamically
    min_FoldEnrichment <- min(combined_df$FoldEnrichment, na.rm = TRUE)
    max_FoldEnrichment <- max(combined_df$FoldEnrichment, na.rm = TRUE)
    
    if (min_FoldEnrichment < 0 && max_FoldEnrichment > 0) {
      # Mixed FoldEnrichment values
      color_scale <- scale_color_gradient2(
        low = "#A43820",
        mid = "#FFEB94",
        high = "#7CAA2D",
        midpoint = 0,
        guide = guide_colorbar(title = "FoldEnrichment", reverse = FALSE)
      )
    } else if (max_FoldEnrichment <= 0) {
      # Only negative FoldEnrichment – red sequential
      color_scale <- scale_color_gradient(
        low = "#A43820",
        high = "#FFEB94",
        guide = guide_colorbar(title = "FoldEnrichment", reverse = FALSE)
      )
    } else {
      # Only positive FoldEnrichment – green sequential
      color_scale <- scale_color_gradient(
        low = "#FFEB94",
        high = "#7CAA2D",
        guide = guide_colorbar(title = "FoldEnrichment", reverse = FALSE)
      )
    }
    
    n_display <- nrow(combined_df)
    n_total <- nrow(df)
    
    # Build dotplot
    dotplot <- ggplot(combined_df, aes(x = source, y = Description, color = FoldEnrichment)) +
      geom_point(aes(size = Count), alpha = 0.8) +
      facet_grid(~location) +
      color_scale +
      labs(
        title = paste("Top", n_display, "out of", n_total, "enriched pathways /", db_name),
        subtitle = "ORA conducted with: padj < 0.05 (BH corr.) | 15 < Genes setSize < 300 | Permutations = 50k",
        x = NULL,
        y = NULL
      ) +
      theme(
        plot.background = element_rect(fill = "#F0F0F0"),
        panel.background = element_rect(fill = "#FFFFFF"),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        plot.title = element_text(color = "#304852", family = "Helvetica Neue", face = "bold", hjust = 0, size = rel(2)),
        plot.subtitle = element_text(color = "#304852", family = "Helvetica Neue", hjust = 0, size = rel(1.3)),
        axis.title.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1), hjust = 0.5),
        axis.text.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1), hjust = 0.5),
        axis.title.y = element_text(family = "Helvetica Neue", size = rel(1)),
        axis.text.y = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(0.9), face = "italic", hjust = 1, margin = margin(t = 30, b = 30, r = 5, l = 5)),
        axis.line = element_line(color = "#DADADA", linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(color = "#36595F", linetype = "solid", linewidth = 0.25),
        panel.grid.major = element_line(color = "#F0F0F0", linetype = "solid", linewidth = 0.25),
        strip.background = element_rect(fill = "#36595F"),
        strip.text = element_text(size = 12, face = "bold", color = "#FFFFFF")
      )
    
    # Save the file
    file_name <- paste0("gsea_xconditions_", cell_location, "_", db_name, ".png")
    ggsave(
      filename = file_name,
      plot = dotplot,
      path = save_path,
      scale = 0.8,
      width = 20,
      height = 13,
      units = "in",
      dpi = 600
    )
    
    message(paste(file_name, "saved in:", save_path))
    return(dotplot)
  })
}

plot_enrich_integration <- function(enrich_list, save_path = "figures/07_enrichments_integration") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  custom_colors <- c(
    "DNA metabolism & processing" = "#007191",
    "RNA metabolism & processing" = "#62c8d3",
    "Translation & ribosome biology" = "#f47a00",
    "Immune & stress response" = "#d31f11",
    "Other" = "#BBBBBB"
  )
  
  lapply(names(enrich_list), function(integration_df) {
    df <- enrich_list[[integration_df]]
    
    df <- df %>%
      dplyr::mutate(
        logGSEA = -log10(gsea_padj + 1e-10),
        logORA  = -log10(ora_padj + 1e-10)
      )
    
    # Label top pathways with highest significance
    label_df <- df %>%
      filter(logGSEA > -log10(0.05) & logORA > -log10(0.05)) %>%
      dplyr::arrange(gsea_padj + ora_padj) %>%
      dplyr::slice_head(n = 20)
    
    
    enrichplot <- ggplot(df, aes(x = logGSEA, y = logORA, color = pathway_relation)) +
      geom_point(size = 1.5, alpha = 0.4) +
      geom_hline(yintercept = -log10(0.05), col = "#dd9d6b", linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
      geom_vline(xintercept = -log10(0.05), col = "#dd9d6b", linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
      geom_label_repel(
        data = label_df,
        aes(label = Description),
        size = 2.5,
        box.padding = 0.4,
        segment.color = "#d7d7d7",
        force = 50,
        segment.size = 0.3,
        min.segment.length = 1,
        max.overlaps = 25,
        direction = "both",
        hjust = 0.5,
        show.legend = FALSE
      ) +
      scale_color_manual(values = custom_colors) +
      coord_cartesian(
        xlim = c(0, 11),
        ylim = c(0, 11)
      ) +
      labs(
        title = paste("Enrichments integration (GSEA & ORA) /", integration_df),
        subtitle = paste("Labeled pathways: top 20 of lowest GSEA & ORA adjusted p-values (sign. cutoff = 0.05 & BH correction)"),
        x = "-log10 GSEA adjusted pvalue",
        y = "-log10 ORA adjusted pvalue",
        color = "Pathway Category"
      ) +
      theme(
        plot.background = element_rect(fill = "#F0F0F0"),
        panel.background = element_rect(fill = "#FFFFFF"),
        plot.margin = margin(20, 20, 20, 20, unit = "pt"),
        plot.title = element_text(color = "#304852", family = "Helvetica Neue", face = "bold", hjust = 0, size = rel(2)),
        plot.subtitle = element_text(color = "#304852", family = "Helvetica Neue", hjust = 0, size = rel(1.3)),
        axis.title.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1), hjust = 0.5),
        axis.text.x = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1), hjust = 0.5),
        axis.title.y = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(1), hjust = 0.5),
        axis.text.y = element_text(color = "#36595F", family = "Helvetica Neue", size = rel(0.9)),
        axis.line = element_line(color = "#DADADA", linewidth = 0.25),
        axis.ticks = element_line(color = "#36595F", linewidth = 0.25),
        panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.25),
        strip.background = element_rect(fill = "#36595F"),
        strip.text = element_text(size = 12, face = "bold", color = "#FFFFFF"),
        legend.position = "top",
        legend.box.margin = margin(0.1, 0.1, 0.1, 0.1),
        legend.background = element_rect(fill = "#FFFFFF", color = "#b0b0b0", linewidth = 0.2),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(color = "#36595F", family = "Helvetica Neue", face = "bold", size = rel(1)),
        legend.text = element_text(color = "#36595F", family = "Helvetica Neue", face = "italic", size = rel(0.9))
      )
    
    file_name <- paste0(integration_df, "_enrich_integration_plot.png")
    ggsave(
      filename = file_name,
      plot = enrichplot,
      path = save_path,
      width = 12, height = 9,
      units = "in",
      dpi = 600
    )
    
    message(paste("saved:", file_name))
    return(enrichplot)
  })
}

create_master_volcano <- function(enrich_category,
                                  db_name,
                                  cell_compartment,
                                  condition,
                                  pathway_description,
                                  save_path = NULL) {
  
  toptable_name <- paste0(cell_compartment, "_", condition, "_HC")
  enrich_name <- paste0(db_name, "_", cell_compartment, "_", condition)
  
  # Select enrichment data
  enrich_df <- master_enrich_data[[enrich_name]]
  
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
  
  toptable <- toptables[[toptable_name]] %>%
    mutate(
      diffexpressed = case_when(
        logFC_shrunk > logFC_volcano_cutoff & adj.P.Val < padj_volcano_cutoff ~ "up-regulated",
        logFC_shrunk < -logFC_volcano_cutoff & adj.P.Val < padj_volcano_cutoff ~ "down-regulated",
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
    geom_point(data = toptable %>% filter(is.na(gene_label)), size = 0.2, alpha = 0.25) +
    geom_point(data = toptable %>% filter(!is.na(gene_label)), size = 0.8, alpha = 0.8) +
    geom_vline(xintercept = c(-logFC_volcano_cutoff, logFC_volcano_cutoff),
               col = "#dd9d6b", linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_volcano_cutoff),
               col = "#dd9d6b", linetype = "dashed") +
    geom_label_repel(data = subset(toptable, highlight == "Involved"),
                     aes(label = gene_label),
                     size = 3, box.padding = 0.25, segment.color = "#d7d7d7",
                     max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(name = "Genes that are:",
                       values = c("Involved" = "#C43A4F", "Not involved" = "#CFCFCF")) +
    coord_cartesian(ylim = c(0, 20), xlim = c(-8, 8)) +
    scale_x_continuous(breaks = seq(-8, 8, 2)) +
    scale_y_continuous(breaks = seq(0, 20, 2)) +
    labs(
      title = paste0("Volcano plot / ", cell_compartment, "_", condition, "_HC"),
      subtitle = paste0("Highlighted pathway: ", pathway_description, "\n",
                        db_name, " database (", toupper(enrich_category), ")"),
      x = "Shrunk logFC (ashr)",
      y = "-log10 adjusted p-value"
    ) +
    plot_theme
  
  # Save if path is provided
  if (!is.null(save_path)) {
    file_name <- paste0(
      toupper(db_name), "_", enrich_category, "_volcano_",
      condition, "_", str_replace_all(tolower(pathway_description), "\\s+", "_"), ".png"
    )
    ggsave(file.path(save_path, file_name),
           plot = volcano_plot, width = 8, height = 7.8, units = "in", dpi = 600)
    message(paste("Plot saved at:", file.path(save_path, file_name)))
  }
  
  return(volcano_plot)
}