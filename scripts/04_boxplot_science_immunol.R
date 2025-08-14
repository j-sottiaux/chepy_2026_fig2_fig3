### 0/ Script setup --------------------------------------------------------------

# Load required libraries
{
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)

}

# Reproducibility check
gc()
set.seed(12345)

# Define the save path
save_path <- "figures/bonus_lfq_intensities_boxplots"

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}


# Select the protein to be plotted here
gene_of_interest <- "HLA-C" 

# Function to extract and reshape for one dataset
extract_gene_data <- function(df, gene, source_name) {
  
  gene_row <- df[rownames(df) == gene, ]
  
  # Define LFQ column mapping (
  gene_long <- tibble(
    Condition = c(rep("dcSSc_ATAp", 10),   # LFQ 1–10
                  rep("dcSSc_ATAn", 10),   # LFQ 11–20
                  rep("lcSSc_ACA", 10),    # LFQ 21–30
                  rep("HC", 10)),          # LFQ 31–40
    LFQ = as.numeric(gene_row[1, 1:40])    # Only first 40 columns (drop NS but we can put it back if we need to)
  )
  
  gene_long$Compartment <- source_name
  return(gene_long)
}

# Apply to the two cellular compartments
cytoplasm_data <- NULL
nucleus_data   <- NULL

if (gene_of_interest %in% rownames(cytoplasm_df)) {
  cytoplasm_data <- extract_gene_data(cytoplasm_df, gene_of_interest, "Cytoplasm")
}

if (gene_of_interest %in% rownames(nucleus_df)) {
  nucleus_data <- extract_gene_data(nucleus_df, gene_of_interest, "Nucleus")
}

# Combine only available compartments
plot_data <- bind_rows(
  if (!is.null(cytoplasm_data)) cytoplasm_data,
  if (!is.null(nucleus_data)) nucleus_data
)

# If no compartments found, stop
if (nrow(plot_data) == 0) {
  stop(paste("Gene", gene_of_interest, "not found in any compartment."))
}


# Order and colors based on Marie-Elise IF outputs / can be changed here
condition_order <- c("dcSSc_ATAp", "dcSSc_ATAn", "lcSSc_ACA", "HC")
condition_colors <- c(
  "dcSSc_ATAp" = "#73A9E5",
  "dcSSc_ATAn" = "#88D0E5",
  "lcSSc_ACA"  = "#F71735",
  "HC"         = "#C8C6BD"
)
condition_comparisons <- list(
  c("dcSSc_ATAp", "dcSSc_ATAn"),
  c("dcSSc_ATAn", "lcSSc_ACA"),
  c("dcSSc_ATAp", "lcSSc_ACA"),
  c("lcSSc_ACA", "HC"),
  c("dcSSc_ATAn", "HC"),
  c("dcSSc_ATAp", "HC")
)

# Condition into ordered factor to plot in the order we like
plot_data$Condition <- factor(plot_data$Condition, levels = condition_order)

# Create the boxplot
lfq_boxplot <- ggplot(plot_data, aes(x = Condition, y = LFQ, fill = Condition, color = Condition)) +
  #geom_boxplot(outlier.shape = NA, alpha = 0.4, linewidth = 0.3) +             # shape set to boxplot by default here
  geom_violin(trim = FALSE, alpha = 0.4) +                                      # change shape to violin here
  geom_jitter(width = 0.1, size = 1, alpha = 0.5) +
  scale_fill_manual(values = condition_colors) +
  scale_color_manual(values = condition_colors) +
  stat_summary(fun = mean, geom = "point", shape  = 18, size = 3, alpha = 0.9, color = "#5f5f5f") +
  stat_compare_means(
    comparisons = condition_comparisons,
    method = "wilcox.test",      
    label = "p.signif",
  ) +
  facet_wrap(~Compartment, ncol = 1, scales = "free_y") +
  labs(
    title = paste("LFQ intensities for", gene_of_interest),
    y = "LFQ intensity",
    x = "Conditions"
  ) +
  theme(
    # General aspect of the plot
    plot.background = element_rect(fill = "#FFFFFF"),
    panel.background = element_rect(fill = "#FFFFFF"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
    plot.title = element_text(color = "black", family = "Helvetica Neue", face = "bold", hjust = 0, size = rel(2)),
    plot.subtitle = element_text(color = "black", family = "Helvetica Neue", hjust = 0, size = rel(1.2)),
    
    # Axis titles and texts
    axis.title.x = element_text(color = "black", family = "Helvetica Neue", size = rel(1)),
    axis.text.x = element_text(color = "black", family = "Helvetica Neue", size = rel(0.9), angle = 0, hjust = 0.5),
    axis.title.y = element_text(color = "black", family = "Helvetica Neue", size = rel(1)), ,
    axis.text.y = element_text(color = "black", family = "Helvetica Neue", size = rel(0.9), angle = 0, hjust = 0.5),
    axis.line = element_line(color = "#5f5f5f", linetype = "solid", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linetype = "solid", linewidth = 0.25),
    panel.grid.major = element_line(color = "#EAEAEA", linetype = "dotted", linewidth = 0.25),
    
    # Legend
    legend.box.margin = margin(0.1, 0.1, 0.1, 0.1),
    legend.position = "bottom",
    legend.background = element_rect(fill = "#FFFFFF", color = "#5f5f5f", linewidth = 0.1, linetype = "solid"),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(color = "black", family = "Helvetica Neue", face = "bold", size = rel(1)),
    legend.text = element_text(color = "black", family = "Helvetica Neue", face = "italic", size = rel(0.9))
  )


# Save the plot and ensure the name matches the shape 
file_name <- paste0("violin_plot_lfq_intensities_", gene_of_interest, ".png")
ggsave(
  filename = file_name,
  plot = lfq_boxplot,
  path = save_path,
  width = 7.7, height = 7.5,
  units = "in",
  dpi = 600
)

