create_volcano_aarhus <- function(df, save_path = "figures/aarhus_volcanos") {
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  nickname <- deparse(substitute(df))

  # Create differential expression cutoffs
  df$diffexpressed <- "no"
  df$diffexpressed[df$log2FoldChange > logFC_volcano_cutoff & df$padj < padj_volcano_cutoff] <- "up"
  df$diffexpressed[df$log2FoldChange < -logFC_volcano_cutoff & df$padj < padj_volcano_cutoff] <- "down"

  # Select the top 20 genes
  df_sig <- df |>
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(padj),
      padj > 0,
      abs(log2FoldChange) > logFC_volcano_cutoff,
      padj < padj_volcano_cutoff
    ) |>
    dplyr::arrange(desc(abs(log2FoldChange)))

  df_top20 <- head(df_sig$gene_symbol, 20)


  # Label top genes
  df$gene_label <- ifelse(df$gene_symbol %in% df_top20 & df$diffexpressed != "no",
    df$gene_symbol,
    NA
  )

  # Generate the volcano plot
  volcano_df <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed, label = gene_label)) +
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
      axis.title.y = element_text(color = "black", family = "Helvetica Neue", size = rel(0.9)), ,
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
    dpi = 600,
    # device = "tiff"
  )

  message(paste(file_name, "saved in:", save_path))
  return(volcano_df)
}

aarhus_ata_hc <- fread("data/aarhus_data/results_tables/ata_hc_igg_rna.csv")
aarhus_ata_ifna <- fread("data/aarhus_data/results_tables/ata_ifna_rna.csv")
aarhus_ifna_hc <- fread("data/aarhus_data/results_tables/ifna_hc_rna.csv")

create_volcano_aarhus(aarhus_ata_hc)
create_volcano_aarhus(aarhus_ata_ifna)
create_volcano_aarhus(aarhus_ifna_hc)

aarhus_rna_matrix <- read.table(
  file = "data/aarhus_data/results_tables/read_counts_rna.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

aarhus_rna_matrix_filt <- aarhus_rna_matrix[
  (
    (rowSums(aarhus_rna_matrix[, c("ATA_1", "ATA_2", "ATA_3")] == 0) > 0) +
      (rowSums(aarhus_rna_matrix[, c("HC_IgG_1", "HC_IgG_2", "HC_IgG_3")] == 0) > 0) +
      (rowSums(aarhus_rna_matrix[, c("IFNa_1", "IFNa_2", "IFNa_3")] == 0) > 0)
  ) < 2,
]
