# ==== Load Libraries ====
library(pROC)
library(pheatmap)
library(readr)
library(dplyr)
library(nnet)  # for multinom

# ==== Handle Argument ====
#args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")
if (length(args) < 1) stop("Usage: Rscript auc_heatmap.R <output_folder>")
output_folder <- args[1]

# ==== Read Metadata ====
metadata_path <- file.path(output_folder, "metadata.csv")
metadata <- read_csv(metadata_path)
metadata$sample_id <- as.character(metadata$sample_id)
metadata$batchid <- as.factor(metadata$batchid)

# ==== Find All Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))

# Add Raw file explicitly
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== Function to Compute AUC per Batch ====
compute_auc_matrix <- function(file_list, metadata, group_col = "batchid") {
  method_names <- names(file_list)
  batch_levels <- sort(unique(metadata[[group_col]]))
  auc_matrix <- matrix(NA, nrow = length(method_names), ncol = length(batch_levels),
                       dimnames = list(method_names, paste0("Batch_", batch_levels)))
  
  for (method in method_names) {
    df <- read_csv(file_list[[method]])
    df$sample_id <- metadata$sample_id
    df_merged <- inner_join(df, metadata, by = "sample_id")
    
    X <- df_merged %>% select(where(is.numeric)) %>% select_if(~ sd(.) > 0)
    y <- df_merged[[group_col]]
    
    for (b in batch_levels) {
      y_bin <- as.numeric(y == b)
      model <- try(glm(y_bin ~ ., data = cbind(y_bin, X), family = "binomial"), silent = TRUE)
      if (!inherits(model, "try-error")) {
        pred <- predict(model, newdata = X, type = "response")
        roc_obj <- roc(y_bin, pred, quiet = TRUE)
        auc_val <- as.numeric(auc(roc_obj))
        auc_matrix[method, paste0("Batch_", b)] <- auc_val
      }
    }
  }
  
  return(auc_matrix)
}

# ==== Compute AUC Matrix ====
auc_mat <- compute_auc_matrix(file_list, metadata, group_col = "batchid")

# ==== Save Heatmap as TIFF ====
tiff(file.path(output_folder, "auc_heatmap.tif"), width = 1000, height = 800, res = 150)
pheatmap(
  auc_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Population Effect AUC",
  fontsize_number = 12,
  fontsize = 14,
  color = colorRampPalette(c("white", "orange", "red"))(100),
  breaks = seq(0.5, 1, length.out = 101)
)
dev.off()

# ==== Save Heatmap as PNG ====
png(file.path(output_folder, "auc_heatmap.png"), width = 1000, height = 800, res = 150)
pheatmap(
  auc_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Population Effect AUC",
  fontsize_number = 12,
  fontsize = 14,
  color = colorRampPalette(c("white", "orange", "red"))(100),
  breaks = seq(0.5, 1, length.out = 101)
)
dev.off()
