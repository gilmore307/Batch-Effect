# ==== Load Required Libraries ====
library(vegan)
library(pheatmap)
library(readr)
library(dplyr)
library(gridExtra)
library(grid)

# ==== Handle Argument ====
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")  # For testing
if (length(args) < 1) stop("Usage: Rscript bray-curtis_heatmap.R <output_folder>")
output_folder <- args[1]

# ==== Read Metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"))
metadata$sample_id <- as.character(metadata$sample_id)

# ==== Find All Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== Bray-Curtis Batch-Level Heatmap Generator ====
generate_batch_heatmap_gtable <- function(file_path,
                                          method_name,
                                          metadata,
                                          group_col = "batchid") {
  
  ## ---------- Load & join -------------------------------------------------
  df <- read_csv(file_path)
  df$sample_id <- metadata$sample_id
  df_merged <- inner_join(df, metadata, by = "sample_id") %>%
    filter(complete.cases(.))
  
  ## ---------- Extract numeric (taxa) matrix -------------------------------
  taxa <- df_merged %>%
    select(where(is.numeric)) %>%
    select_if(~ sd(., na.rm = TRUE) > 0) %>%
    as.matrix()
  
  ## ---------- NEW: normalise to relative abundance ------------------------
  if (any(taxa < 0, na.rm = TRUE))
    stop("Negative values detected – Bray-Curtis needs non-negative data.")
  
  row_totals <- rowSums(taxa, na.rm = TRUE)
  keep <- row_totals > 0        # drop empty samples
  taxa   <- taxa[keep, , drop = FALSE]
  batches <- df_merged[[group_col]][keep]
  sample_ids <- df_merged$sample_id[keep]
  
  taxa_rel <- sweep(taxa, 1, row_totals[keep], FUN = "/")
  
  ## ---------- Distance & aggregation -------------------------------------
  dist_mat <- vegdist(taxa_rel, method = "bray")
  dist_df  <- as.matrix(dist_mat)
  rownames(dist_df) <- colnames(dist_df) <- sample_ids
  
  batch_levels  <- sort(unique(batches))
  batch_dist_mat <- matrix(NA, length(batch_levels), length(batch_levels),
                           dimnames = list(batch_levels, batch_levels))
  
  for (i in batch_levels)
    for (j in batch_levels) {
      s_i <- sample_ids[batches == i]
      s_j <- sample_ids[batches == j]
      if (length(s_i) && length(s_j)) {
        batch_dist_mat[i, j] <- mean(dist_df[s_i, s_j], na.rm = TRUE)
      }
    }
  
  ## ---------- Cosmetic: avoid single-value colour scales ------------------
  if (length(unique(na.omit(c(batch_dist_mat)))) == 1)
    batch_dist_mat <- batch_dist_mat +
    matrix(runif(length(batch_dist_mat), -1e-6, 1e-6),
           nrow = nrow(batch_dist_mat))
  
  ## ---------- Plot --------------------------------------------------------
  max_val <- max(batch_dist_mat, na.rm = TRUE)
  breaks  <- seq(0, max_val, length.out = 101)
  
  p <- pheatmap(batch_dist_mat,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                display_numbers = TRUE,
                fontsize_number = 14,
                fontsize = 14,
                main = paste("Inter-Batch Bray-Curtis –", method_name),
                color = colorRampPalette(c("blue", "white", "red"))(100),
                breaks = breaks,
                na_col = "gray",
                silent = TRUE)
  
  p$gtable
}

# ==== Generate Heatmaps ====
heatmap_list <- lapply(names(file_list), function(name) {
  generate_batch_heatmap_gtable(file_list[[name]], name, metadata)
})

# ==== Layout for Square Subplots ====
ncol_grid <- 3
nrow_grid <- ceiling(length(heatmap_list) / ncol_grid)
plot_size <- 800  # pixels per subplot (square)

# ==== Save Combined Heatmaps ====
tiff(file.path(output_folder, "braycurtis.tiff"),
     width = plot_size * ncol_grid,
     height = plot_size * nrow_grid,
     res = 150)
grid.arrange(grobs = heatmap_list, ncol = ncol_grid)
dev.off()

png(file.path(output_folder, "braycurtis.png"),
    width = plot_size * ncol_grid,
    height = plot_size * nrow_grid,
    res = 150)
grid.arrange(grobs = heatmap_list, ncol = ncol_grid)
dev.off()
