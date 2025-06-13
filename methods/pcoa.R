# ==== Load Required Libraries ====
library(ggplot2)
library(patchwork)
library(readr)
library(dplyr)
library(cluster)
library(vegan)       # For vegdist
library(ggfortify)   # Optional, used in PCA

# ==== Read Metadata ====
metadata <- read_csv("metadata.csv")
metadata$sample_id <- as.character(metadata$sample_id)

# ==== File Paths ====
file_list <- list(
  Raw = "raw.csv",
  ALRA = "normalized_alra.csv",
  ComBat = "normalized_combat.csv",
  ConQuR = "normalized_conqur.csv",
  FSQN = "normalized_fsqn.csv",
  PLSDA = "normalized_plsda.csv"
)

# ==== PCoA Plot Function Using Metadata Grouping ====
pcoa_with_metadata_plot <- function(df, method_name, metadata, group_col = "batchid") {
  df$sample_id <- metadata$sample_id  # Add sample_id from metadata
  
  # Merge with metadata
  df_merged <- inner_join(df, metadata, by = "sample_id")
  
  # Extract numeric columns (features only)
  df_num <- df_merged %>%
    select(where(is.numeric)) %>%
    select(-any_of(group_col))
  
  # Remove constant columns
  df_num <- df_num[, apply(df_num, 2, function(x) sd(x) > 0), drop = FALSE]
  
  # Set rownames to sample IDs for distance matrix
  rownames(df_num) <- df_merged$sample_id
  
  # Compute Bray-Curtis distance
  dist_mat <- vegdist(df_num, method = "bray")
  
  # Run PCoA
  pcoa_result <- cmdscale(dist_mat, k = 2, eig = TRUE)
  pcoa_df <- as.data.frame(pcoa_result$points)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  pcoa_df$sample_id <- rownames(pcoa_df)
  
  # Add group info
  pcoa_df <- inner_join(pcoa_df, metadata, by = "sample_id")
  pcoa_df$group <- as.factor(pcoa_df[[group_col]])
  
  # Plot
  ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = group)) +
    geom_point(size = 2) +
    stat_ellipse(aes(group = group), type = "norm", linetype = 2, alpha = 0.5) +
    labs(title = paste("PCoA -", method_name), color = group_col) +
    theme_minimal()
}

# ==== Apply PCoA Plotting with Metadata Coloring ====
plots <- lapply(names(file_list), function(name) {
  df <- read_csv(file_list[[name]])
  pcoa_with_metadata_plot(df, name, metadata, group_col = "batchid")
})

# ==== Combine All Plots ====
combined_plot <- wrap_plots(plots, ncol = 2)
print(combined_plot)
