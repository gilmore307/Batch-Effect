# ==== Load Libraries ====
library(pROC)
library(pheatmap)
library(readr)
library(dplyr)
library(nnet)  # for multinom

# ==== Read Metadata ====
metadata <- read_csv("metadata.csv")
metadata$sample_id <- as.character(metadata$sample_id)
metadata$batchid <- as.factor(metadata$batchid)

# ==== File Paths ====
file_list <- list(
  Raw = "raw.csv",
  ALRA = "normalized_alra.csv",
  ComBat = "normalized_combat.csv",
  ConQuR = "normalized_conqur.csv",
  FSQN = "normalized_fsqn.csv",
  PLSDA = "normalized_plsda.csv"
)

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
    
    # Get numeric features and batch
    X <- df_merged %>% select(where(is.numeric)) %>% select_if(~ sd(.) > 0)
    y <- df_merged[[group_col]]
    
    # One-vs-rest AUC per batch
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

# ==== Plot Heatmap ====
pheatmap(
  auc_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Population Effect AUC",
  fontsize_number = 12,
  fontsize = 14,
  color = colorRampPalette(c("white", "orange", "red"))(100),
  breaks = seq(0.5, 1, length.out = 101)  # AUC range [0.5, 1]
)
