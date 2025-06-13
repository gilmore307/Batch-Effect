# ==== Load Required Libraries ====
library(vegan)
library(pheatmap)
library(readr)
library(dplyr)
library(gridExtra)
library(grid)

# ==== Handle Argument ====
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")  # For testing, replace with commandArgs in production
if (length(args) < 1) stop("Usage: Rscript bray-curtis_heatmap.R <output_folder>")
output_folder <- args[1]

# ==== Read Metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"))
metadata$sample_id <- as.character(metadata$sample_id)

# ==== Find All Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))

# Add Raw file explicitly
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== Bray-Curtis Batch-Level Heatmap Generator ====
generate_batch_heatmap_gtable <- function(file_path, method_name, metadata, group_col = "batchid") {
  df <- read_csv(file_path)
  df$sample_id <- metadata$sample_id
  df_merged <- inner_join(df, metadata, by = "sample_id")
  
  # Drop rows with NAs
  df_merged <- df_merged %>% filter(complete.cases(.))
  
  # Extract valid features
  taxa <- df_merged %>%
    select(where(is.numeric)) %>%
    select_if(~ sd(., na.rm = TRUE) > 0) %>%
    as.matrix()
  
  sample_ids <- df_merged$sample_id
  batches <- df_merged[[group_col]]
  
  if (nrow(taxa) < 2 || ncol(taxa) < 1) {
    warning(paste("Skipping method", method_name, "due to insufficient data."))
    return(grid.text(paste("Insufficient data:", method_name)))
  }
  
  dist_mat <- vegdist(taxa, method = "bray")
  dist_df <- as.matrix(dist_mat)
  rownames(dist_df) <- colnames(dist_df) <- sample_ids
  
  batch_levels <- sort(unique(batches))
  batch_dist_mat <- matrix(NA, length(batch_levels), length(batch_levels),
                           dimnames = list(batch_levels, batch_levels))
  
  # Compute inter-batch means
  for (i in batch_levels) {
    for (j in batch_levels) {
      samples_i <- sample_ids[batches == i]
      samples_j <- sample_ids[batches == j]
      
      if (length(samples_i) == 0 || length(samples_j) == 0) {
        batch_dist_mat[as.character(i), as.character(j)] <- NA
      } else {
        sub_dists <- dist_df[samples_i, samples_j, drop = FALSE]
        batch_dist_mat[as.character(i), as.character(j)] <- mean(sub_dists, na.rm = TRUE)
      }
    }
  }
  
  # Handle NA and constant matrix
  batch_dist_mat[is.na(batch_dist_mat)] <- NA
  if (length(unique(na.omit(c(batch_dist_mat)))) == 1) {
    batch_dist_mat <- batch_dist_mat + matrix(runif(length(batch_dist_mat), -1e-6, 1e-6),
                                              nrow = nrow(batch_dist_mat))
  }
  
  max_val <- max(batch_dist_mat, na.rm = TRUE)
  breaks <- seq(0, max_val, length.out = 101)
  
  p <- pheatmap(
    batch_dist_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    fontsize_number = 14,
    fontsize = 14,
    main = paste("Inter-Batch Bray-Curtis -", method_name),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = breaks,
    na_col = "gray",
    silent = TRUE
  )
  
  return(p$gtable)
}

# ==== Generate Heatmaps ====
heatmap_list <- lapply(names(file_list), function(name) {
  generate_batch_heatmap_gtable(file_list[[name]], name, metadata)
})

# ==== Save Combined Heatmap as TIFF ====
tiff(file.path(output_folder, "braycurtis.tif"), width = 1600, height = 1200, res = 150)
grid.arrange(grobs = heatmap_list, ncol = 2)
dev.off()

# ==== Save Combined Heatmap as PNG ====
png(file.path(output_folder, "braycurtis.png"), width = 1600, height = 1200, res = 150)
grid.arrange(grobs = heatmap_list, ncol = 2)
dev.off()
