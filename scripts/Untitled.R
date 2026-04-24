# Extract VST matrix
transcripto_vsd_mat <- assay(transcripto_vsd)

# PCA on samples, so transpose the gene x sample matrix
transcripto_pca <- prcomp(
  t(transcripto_vsd_mat),
  center = TRUE,
  scale. = FALSE
)

# Variance explained
transcripto_pca_var <- transcripto_pca$sdev^2 / sum(transcripto_pca$sdev^2)

# PCA dataframe
transcripto_pca_df <- as.data.frame(transcripto_pca$x) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::left_join(
    transcripto_metadata %>%
      tibble::rownames_to_column("sample_id"),
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("ATA", "HC", "IFNa")),
    condition_label = factor(
      dplyr::recode(
        as.character(condition),
        ATA  = "ATA",
        HC   = "HC IgG",
        IFNa = "IFNα"
      ),
      levels = c("ATA", "IFNα", "HC IgG")
    )
  )