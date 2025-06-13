# ==== Load Required Libraries ====
library(vegan)
library(pheatmap)
library(readr)
library(dplyr)
library(gridExtra)

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

# ==== Bray-Curtis Batch-Level Heatmap Generator (returns gtable) ====
generate_batch_heatmap_gtable <- function(file_path, method_name, metadata, group_col = "batchid") {
  df <- read_csv(file_path)
  
  # Align sample_id with metadata
  df$sample_id <- metadata$sample_id
  df_merged <- inner_join(df, metadata, by = "sample_id")
  
  # Extract taxa matrix
  taxa <- df_merged %>%
    select(where(is.numeric)) %>%
    select_if(~ sd(.) > 0) %>%
    as.matrix()
  
  sample_ids <- df_merged$sample_id
  batches <- df_merged[[group_col]]
  
  # Compute Bray-Curtis distance
  dist_mat <- vegdist(taxa, method = "bray")
  dist_df <- as.matrix(dist_mat)
  rownames(dist_df) <- colnames(dist_df) <- sample_ids
  
  # Get unique batches
  batch_levels <- sort(unique(na.omit(batches)))
  batch_dist_mat <- matrix(NA, length(batch_levels), length(batch_levels),
                           dimnames = list(batch_levels, batch_levels))
  
  for (i in batch_levels) {
    for (j in batch_levels) {
      samples_i <- sample_ids[batches == i]
      samples_j <- sample_ids[batches == j]
      sub_dists <- dist_df[samples_i, samples_j, drop = FALSE]
      batch_dist_mat[i, j] <- mean(sub_dists, na.rm = TRUE)
    }
  }
  
  # Replace NAs
  batch_dist_mat[is.na(batch_dist_mat)] <- 0
  
  # Avoid flat matrix issues
  if (length(unique(c(batch_dist_mat))) == 1) {
    batch_dist_mat <- batch_dist_mat + matrix(runif(length(batch_dist_mat), -1e-6, 1e-6),
                                              nrow = nrow(batch_dist_mat))
  }
  
  # Breaks
  max_val <- max(batch_dist_mat)
  breaks <- seq(0, max_val, length.out = 101)
  
  # Return gtable (not plot)
  p <- pheatmap(
    batch_dist_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    fontsize_number = 14,
    fontsize = 14,
    main = paste("Inter-Batch Bray-Curtis -", method_name),
    breaks = breaks,
    silent = TRUE  # <- return gtable object
  )
  return(p$gtable)
}

# ==== Generate All Heatmaps as gtable ====
heatmap_list <- lapply(names(file_list), function(name) {
  generate_batch_heatmap_gtable(file_list[[name]], name, metadata)
})

# ==== Combine and Display ====
grid.arrange(grobs = heatmap_list, ncol = 2)
