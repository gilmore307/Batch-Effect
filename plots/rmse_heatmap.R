# ==== Load Required Libraries ====
library(pheatmap)
library(readr)
library(dplyr)
library(gridExtra)
library(grid)

# ==== Handle Argument ====
# args <- commandArgs(trailingOnly = TRUE)
args <- c("plots/output/example")  # For testing
if (length(args) < 1) stop("Usage: Rscript rmse_heatmap.R <output_folder>")
output_folder <- args[1]

# ==== Read Metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"))
metadata$sample_id <- as.character(metadata$sample_id)
metadata$batchid <- as.factor(metadata$batchid)

# ==== Find All Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== RMSE Matrix Computation ====
compute_rmse_matrix <- function(df, metadata, method_name, group_col = "batchid") {
  df$sample_id <- metadata$sample_id
  df_merged <- inner_join(df, metadata, by = "sample_id") %>%
    filter(complete.cases(.))
  
  batches <- df_merged[[group_col]]
  features <- df_merged %>%
    select(where(is.numeric)) %>%
    select_if(~ sd(., na.rm = TRUE) > 0) %>%
    as.matrix()
  
  # Log-transform Raw data only
  if (method_name == "Raw") {
    features <- log1p(features)
  }
  
  batch_levels <- sort(unique(batches))
  rmse_mat <- matrix(NA, nrow = length(batch_levels), ncol = length(batch_levels),
                     dimnames = list(batch_levels, batch_levels))
  
  for (i in batch_levels) {
    for (j in batch_levels) {
      X_i <- features[batches == i, , drop = FALSE]
      X_j <- features[batches == j, , drop = FALSE]
      
      if (nrow(X_i) > 0 && nrow(X_j) > 0) {
        X_i_mean <- colMeans(X_i, na.rm = TRUE)
        X_j_mean <- colMeans(X_j, na.rm = TRUE)
        rmse <- sqrt(mean((X_i_mean - X_j_mean)^2, na.rm = TRUE))
        rmse_mat[as.character(i), as.character(j)] <- round(rmse, 2)
      }
    }
  }
  
  return(rmse_mat)
}

# ==== RMSE Heatmap Plot Function ====
generate_rmse_heatmap_gtable <- function(file_path, method_name, metadata) {
  df <- read_csv(file_path)
  rmse_mat <- compute_rmse_matrix(df, metadata, method_name)
  
  if (all(is.na(rmse_mat))) {
    return(grid.text(paste("Insufficient data:", method_name)))
  }
  
  cat("\n--- RMSE matrix for method:", method_name, "---\n")
  print(rmse_mat)
  
  max_val <- max(rmse_mat, na.rm = TRUE)
  breaks <- seq(0, max_val, length.out = 101)
  
  p <- pheatmap(
    rmse_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    fontsize_number = 14,
    fontsize = 14,
    main = paste("Inter-Batch RMSE -", method_name),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = breaks,
    na_col = "gray",
    silent = TRUE
  )
  
  return(p$gtable)
}

# ==== Generate Heatmaps ====
rmse_heatmap_list <- lapply(names(file_list), function(name) {
  generate_rmse_heatmap_gtable(file_list[[name]], name, metadata)
})

# ==== Set Layout Parameters for Square Subplots ====
ncol_grid <- 3
nrow_grid <- ceiling(length(file_list) / ncol_grid)
plot_size <- 800  # pixels per square panel

# ==== Save Combined RMSE Heatmaps ====
tiff(file.path(output_folder, "rmse.tiff"),
     width = plot_size * ncol_grid,
     height = plot_size * nrow_grid,
     res = 150)
grid.arrange(grobs = rmse_heatmap_list, ncol = ncol_grid)
dev.off()

png(file.path(output_folder, "rmse.png"),
    width = plot_size * ncol_grid,
    height = plot_size * nrow_grid,
    res = 150)
grid.arrange(grobs = rmse_heatmap_list, ncol = ncol_grid)
dev.off()
