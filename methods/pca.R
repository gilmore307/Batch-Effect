# ==== Load Required Libraries ====
library(ggplot2)
library(patchwork)
library(readr)
library(dplyr)
library(ggfortify)
library(cluster)

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

# ==== PCA Plot Function Using Metadata Grouping ====
pca_with_metadata_plot <- function(df, method_name, metadata, group_col = "batchid") {
  df$sample_id <- metadata$sample_id  # Add sample_id from metadata
  
  # Merge with metadata
  df_merged <- inner_join(df, metadata, by = "sample_id")
  
  # Extract numeric columns excluding the grouping column
  df_num <- df_merged %>%
    select(where(is.numeric)) %>%
    select(-any_of(group_col))
  
  # Convert group column to factor
  df_group <- df_merged %>%
    mutate(group = as.factor(!!sym(group_col))) %>%
    select(group)
  
  # Remove constant columns
  df_num <- df_num[, apply(df_num, 2, function(x) sd(x) > 0), drop = FALSE]
  
  # PCA
  pca <- prcomp(df_num, scale. = TRUE)
  
  # Plot with group as color
  autoplot(
    pca,
    data = df_group,
    colour = "group",
    frame = TRUE,
    frame.type = "norm",
    main = paste("PCA -", method_name)
  ) +
    theme_minimal()
}


# ==== Apply PCA Plotting with Metadata Coloring ====
plots <- lapply(names(file_list), function(name) {
  df <- read_csv(file_list[[name]])
  pca_with_metadata_plot(df, name, metadata, group_col = "batchid")
})

# ==== Combine All Plots ====
combined_plot <- wrap_plots(plots, ncol = 2)
print(combined_plot)
