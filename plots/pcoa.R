# ==== Load Required Libraries ====
library(ggplot2)
library(patchwork)
library(readr)
library(dplyr)
library(vegan)
library(ggfortify)

# ==== Read UID argument ====
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")  # for testing; comment out in production
if (length(args) < 1) stop("Usage: Rscript pcoa.R <output_folder>")
output_folder <- args[1]

# ==== Read Metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"))
metadata$sample_id <- as.character(metadata$sample_id)

# ==== Find All Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))

# Add Raw file explicitly
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== PCoA Plot Function ====
pcoa_with_metadata_plot <- function(df, method_name, metadata, group_col = "batchid") {
  df$sample_id <- metadata$sample_id
  df_merged <- inner_join(df, metadata, by = "sample_id")
  
  df_num <- df_merged %>%
    select(where(is.numeric)) %>%
    select(-any_of(group_col))
  
  # Remove NA-only and constant columns
  df_num <- df_num[, colSums(!is.na(df_num)) > 0, drop = FALSE]
  df_num <- df_num[, apply(df_num, 2, function(x) sd(x, na.rm = TRUE)) > 0, drop = FALSE]
  
  # Remove rows with NA values
  complete_rows <- complete.cases(df_num)
  df_num <- df_num[complete_rows, , drop = FALSE]
  df_merged <- df_merged[complete_rows, , drop = FALSE]
  
  rownames(df_num) <- df_merged$sample_id
  
  # Compute Bray-Curtis distance and run PCoA
  dist_mat <- vegdist(df_num, method = "bray")
  pcoa_res <- cmdscale(dist_mat, k = 2, eig = TRUE)
  points_df <- as.data.frame(pcoa_res$points)
  colnames(points_df) <- c("PC1", "PC2")
  
  fake_prcomp <- list(
    x = as.matrix(points_df),
    sdev = sqrt(abs(pcoa_res$eig[1:2])),
    center = FALSE,
    scale = FALSE,
    rotation = matrix(0, nrow = ncol(df_num), ncol = 2)
  )
  rownames(fake_prcomp$rotation) <- colnames(df_num)
  colnames(fake_prcomp$rotation) <- c("PC1", "PC2")
  class(fake_prcomp) <- "prcomp"
  
  df_group <- df_merged %>%
    mutate(group = as.factor(!!sym(group_col))) %>%
    select(group)
  
  autoplot(fake_prcomp, data = df_group, colour = "group", frame = TRUE, frame.type = "norm",
           main = paste("PCoA -", method_name)) +
    theme_minimal()
}

# ==== Generate All PCoA Plots ====
pcoa_plots <- lapply(names(file_list), function(name) {
  df <- read_csv(file_list[[name]])
  pcoa_with_metadata_plot(df, name, metadata, group_col = "batchid")
})

# ==== Combine and Save ====
combined_pcoa <- wrap_plots(pcoa_plots, ncol = 2)
ggsave(file.path(output_folder, "pcoa.tif"), combined_pcoa, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_folder, "pcoa.png"), combined_pcoa, width = 12, height = 10, dpi = 300)
