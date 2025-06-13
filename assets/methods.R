# normalize_all_methods.R
# Normalize ConQuR sample data using various methods: ConQuR, ALRA, BDMMA, PLSDAbatch, ComBat, FSQN

# ---------------------------
# Load Libraries
# ---------------------------
library(ConQuR)
library(foreach)
library(doParallel)
library(dplyr)
library(Matrix)
library(SummarizedExperiment)
library(S4Vectors)
library(compositions)

# Optional: Load if installed
suppressWarnings({
  require(ALRA)
  require(BDMMAcorrect)
  require(PLSDAbatch)
  require(sva)
  require(FSQN)
})

# ---------------------------
# Load Sample Data
# ---------------------------
setwd('C:/Users/sunch/Desktop/Project/Batch-Effect/methods')
data("Sample_Data", package = "ConQuR")
taxa_mat <- Sample_Data[, 1:100]
metadata <- Sample_Data[, 101:ncol(Sample_Data)]
sample_ids <- rownames(Sample_Data)
metadata$sample_id <- rownames(Sample_Data)
covar <- Sample_Data[, c("sbp", "sex", "race", "age")]  # Metadata
batchid <- factor(Sample_Data$batchid)
write.csv(taxa_mat, file = "raw.csv", row.names = FALSE)
write.csv(metadata, file = "metadata.csv", row.names = FALSE)

# ---------------------------
# Normalize using ConQuR
# ---------------------------

conqur_result <- ConQuR(
  tax_tab = taxa_mat,
  batchid = batchid,
  covariates = covar,
  batch_ref = "0"  # Specify the reference batch
)

write.csv(conqur_result, "normalized_conqur.csv", row.names = FALSE)

# ---------------------------
# Normalize using ALRA
# ---------------------------
taxa_filtered <- taxa_mat[rowSums(taxa_mat) > 0, ]
taxa_norm <- sweep(taxa_filtered, 1, rowSums(taxa_filtered), "/")
taxa_norm[is.na(taxa_norm)] <- 0
taxa_norm_t <- t(taxa_norm)
alra_result <- alra(taxa_norm_t)
taxa_alra <- t(alra_result$A_norm_rank_k_cor_sc)
colnames(taxa_alra) <- colnames(taxa_mat)

write.csv(taxa_alra, "normalized_alra.csv", row.names = FALSE)

# ---------------------------
# Normalize using PLSDAbatch
# ---------------------------
# requires mixOmics to be installed
rownames(taxa_mat) <- rownames(Sample_Data)

taxa_clr <- clr(taxa_mat + 1)
Y.trt <- ifelse(metadata$sbp > 130, "high", "normal")    # binary treatment/class
Y.bat <- as.factor(metadata$batchid)   

result <- PLSDA_batch(X = taxa_clr, Y.trt = Y.trt, Y.bat = Y.bat, 
                      ncomp.trt = 1, ncomp.bat = 5)

corrected_data <- result$X.nobatch
write.csv(corrected_data, "normalized_plsda.csv", row.names = FALSE)

# ---------------------------
# Normalize using ComBat
# ---------------------------
taxa_rel_abund <- sweep(taxa_mat, 1, rowSums(taxa_mat), "/")
taxa_rel_abund[is.na(taxa_rel_abund)] <- 0
taxa_log <- log2(taxa_rel_abund + 1e-6)

taxa_log_t <- t(taxa_log)
combat_corrected <- ComBat(dat = taxa_log_t, batch = batchid, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
combat_corrected_final <- t(combat_corrected)
scaled_data <- scale(combat_corrected_final)
taxa_rel_corrected <- 2^scaled_data
write.csv(taxa_rel_corrected, "normalized_combat.csv", row.names = FALSE)

# ---------------------------
# Normalize using FSQN
# ---------------------------
taxa_cols <- which(sapply(Sample_Data, is.numeric))
batch_vector <- Sample_Data$batchid

reference_batch <- names(which.max(table(batch_vector)))
ref_matrix <- taxa_mat[batch_vector == reference_batch, , drop = FALSE]

reference_distribution <- apply(ref_matrix, 2, function(x) sort(x))
reference_target <- apply(reference_distribution, 1, median)  # vector of target values per rank

normalize_to_reference <- function(mat, target_vector) {
  apply(mat, 2, function(feature_column) {
    ranked <- rank(feature_column, ties.method = "min")
    sorted_original <- sort(feature_column)
    normalized <- target_vector[ranked]
    return(normalized)
  })
}

normalized_list <- lapply(unique(batch_vector), function(batch) {
  message("Normalizing batch: ", batch)
  mat <- taxa_mat[batch_vector == batch, , drop = FALSE]
  normed <- normalize_to_reference(mat, reference_target)
  rownames(normed) <- rownames(mat)
  return(normed)
})

normalized_matrix <- do.call(rbind, normalized_list)
normalized_matrix <- normalized_matrix[rownames(taxa_mat), ]

write.csv(normalized_matrix, "normalized_fsqn.csv", row.names = FALSE)
