# ==== Load Required Libraries ====
library(ggplot2)
library(patchwork)
library(readr)
library(dplyr)
library(ggfortify)
library(cluster)

# ==== Read UID argument ====
#args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")
if (length(args) < 1) stop("Usage: Rscript pca.R <output_folder>")
output_folder <- args[1]

# ==== Read Metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"))
metadata$sample_id <- as.character(metadata$sample_id)

# ==== Find All Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))

# Add Raw file explicitly
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== PCA Plot Function Using Metadata Grouping ====
pca_with_metadata_plot <- function(df, method_name, metadata, group_col = "batchid") {
  df$sample_id <- metadata$sample_id
  df_merged <- inner_join(df, metadata, by = "sample_id")
  
  df_num <- df_merged %>%
    select(where(is.numeric)) %>%
    select(-any_of(group_col)) %>%
    select_if(~ sd(.) > 0)  # remove constant cols
  
  df_group <- df_merged %>%
    mutate(group = as.factor(!!sym(group_col))) %>%
    select(group)
  
  pca <- prcomp(df_num, scale. = TRUE)
  
  autoplot(
    pca,
    data = df_group,
    colour = "group",
    frame = TRUE,
    frame.type = "norm",
    main = paste("PCA -", method_name)
  ) + theme_minimal()
}

# ==== Generate PCA Plots ====
plots <- lapply(names(file_list), function(name) {
  df <- read_csv(file_list[[name]])
  pca_with_metadata_plot(df, name, metadata, group_col = "batchid")
})

# ==== Combine and Save ====
combined_plot <- wrap_plots(plots, ncol = 2)
ggsave(file.path(output_folder, "pca.tif"), combined_plot, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_folder, "pca.png"), combined_plot, width = 12, height = 10, dpi = 300)

