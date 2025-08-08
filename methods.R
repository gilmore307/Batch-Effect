# ---------------------------
# Handle Arguments
# ---------------------------
#args <- commandArgs(trailingOnly = TRUE)
args <- c("MMUPHin,CQR,RUV,MetaDICT,SVD,PN,FAbatch,ComBatSeq", "output/example")

if (length(args) < 2) {
  stop("Usage: Rscript normalize_all_methods.R <method_list> <output_folder>")
}

method_list <- unlist(strsplit(args[1], ","))
output_folder <- args[2]

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# ---------------------------
# Load Libraries
# ---------------------------
library(dplyr)
library(Matrix)
library(SummarizedExperiment)
library(S4Vectors)

suppressWarnings({
  require(ALRA)
  require(BDMMAcorrect)
  require(PLSDAbatch)
  require(sva)
  require(FSQN)
  require(edgeR)
  require(DESeq2)
  require(metagenomeSeq)
  require(GUniFrac)
  require(preprocessCore)
  require(pamr)
  require(limma)
  require(huge)
  require(MMUPHin)
})

# ---------------------------
# Load Data
# ---------------------------
custom_matrix_path <- file.path(output_folder, "raw.csv")
custom_metadata_path <- file.path(output_folder, "metadata.csv")
default_matrix_path <- file.path("assets", "raw.csv")
default_metadata_path <- file.path("assets", "metadata.csv")

if (file.exists(custom_matrix_path) && file.exists(custom_metadata_path)) {
  cat("✅ Using uploaded user files\n")
  taxa_mat <- read.csv(custom_matrix_path, check.names = FALSE)
  metadata <- read.csv(custom_metadata_path)
} else {
  cat("⚠️ No uploaded files found — using default assets\n")
  taxa_mat <- read.csv(default_matrix_path, check.names = FALSE)
  metadata <- read.csv(default_metadata_path)
  write.csv(taxa_mat, file = file.path(output_folder, "raw.csv"), row.names = FALSE)
  write.csv(metadata, file = file.path(output_folder, "metadata.csv"), row.names = FALSE)
}

rownames(taxa_mat) <- metadata$sample_id
rownames(metadata) <- metadata$sample_id

# Ensure required metadata
required_columns <- c("sample_id", "batchid")
missing <- setdiff(required_columns, colnames(metadata))
if (length(missing) > 0) {
  stop(paste("Missing required metadata columns:", paste(missing, collapse = ", ")))
}
covar <- metadata[, !(colnames(metadata) %in% c("sample_id", "batchid")), drop = FALSE]
covar <- as.data.frame(lapply(covar, function(col) {
  if (is.numeric(col)) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
  } else if (is.factor(col)) {
    col[is.na(col)] <- levels(col)[1]
  }
  return(col)
}))

batchid <- factor(metadata$batchid)

# ---------------------------
# Define Reference vs Non-Reference
# ---------------------------
reference_batch <- 0
ref_idx <- which(batchid == reference_batch)
nonref_idx <- which(batchid != reference_batch)

ref_matrix <- taxa_mat[ref_idx, , drop = FALSE]
nonref_matrix <- taxa_mat[nonref_idx, , drop = FALSE]

# ---------------------------
# Helper
# ---------------------------
normalize_tss <- function(mat) {
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1  # avoid division by zero
  mat <- sweep(mat, 1, row_sums, "/")
  mat[is.na(mat)] <- 0
  return(mat)
}

if ("QN" %in% method_list) {
  cat("Running QN...\n")
  library(preprocessCore)
  
  # Ensure matrices are numeric with features in rows, samples in columns
  ref_abs  <- as.matrix(ref_matrix)
  all_abs  <- as.matrix(taxa_mat)
  
  # Determine target distribution from the reference matrix
  ref_target <- normalize.quantiles.determine.target(ref_abs)
  
  # Quantile-normalize all samples to the reference target
  qn_corrected <- normalize.quantiles.use.target(all_abs, target = ref_target)
  
  # Preserve dimnames
  dimnames(qn_corrected) <- dimnames(all_abs)
  
  write.csv(qn_corrected, file.path(output_folder, "normalized_qn.csv"), row.names = FALSE)
}

if ("BMC" %in% method_list) {
  cat("Running BMC...\n")
  require(pamr)
  
  norm <- normalize_tss(taxa_mat)
  norm[norm == 0] <- min(norm[norm != 0]) * 0.65
  log_mat <- log(norm)
  batch_factor <- factor(batchid)
  pam_input <- list(x = as.matrix(t(log_mat)), batchlabels = batch_factor)
  bmc_corrected <- pamr.batchadjust(pam_input)$x
  write.csv(t(bmc_corrected), file.path(output_folder, "normalized_bmc.csv"), row.names = FALSE)
}

if ("limma" %in% method_list) {
  cat("Running limma...\n")
  require(limma)
  
  mat <- as.matrix(taxa_mat)
  mat[mat == 0] <- min(mat[mat != 0]) * 0.65
  log_mat <- log(mat)
  batch_factor <- factor(batchid)
  limma_corrected <- removeBatchEffect(t(log_mat), batch = batch_factor, covariates = as.matrix(covar))
  write.csv(t(limma_corrected), file.path(output_folder, "normalized_limma.csv"), row.names = FALSE)
}

if ("ConQuR" %in% method_list) {
  cat("Running ConQuR with parallel backend...\n")
  
  require(ConQuR)
  require(doParallel)
  
  # Set up parallel backend
  num_cores <- parallel::detectCores(logical = TRUE) - 1
  
  # Use existing batchid and reference_batch
  batchid <- as.factor(metadata$batchid)
  batch_ref <- as.character(reference_batch)  # convert to character as required
  
  # Extract covariates
  covariates <- metadata[, colnames(covar), drop = FALSE]
  rownames(covariates) <- NULL  # ensure no rownames
  
  # Run ConQuR
  conqur_result <- ConQuR(
    tax_tab = taxa_mat,
    batchid = batchid,
    covariates = covariates,
    batch_ref = batch_ref,
    logistic_lasso = FALSE,
    quantile_type = "standard",
    simple_match = FALSE,
    lambda_quantile = "2p/n",
    interplt = FALSE,
    delta = 0.4999,
    taus = seq(0.005, 0.995, by = 0.005),
    num_core = num_cores
  )
  
  # Save result
  write.csv(conqur_result, file.path(output_folder, "normalized_conqur.csv"), row.names = FALSE)
}

if ("PLSDA" %in% method_list) {
  cat("Running PLSDAbatch...\n")
  require(PLSDAbatch)
  require(compositions)
  
  # Ensure rownames are set
  rownames(taxa_mat) <- rownames(metadata)
  
  # Apply CLR transformation using compositions::clr
  taxa_clr <- clr(taxa_mat + 1)  # small offset to avoid log(0)
  
  # Check if phenotype column is valid
  if (!("phenotype" %in% colnames(metadata))) {
    stop("❌ 'phenotype' column not found in metadata.")
  }
  if (length(unique(metadata$phenotype)) != 2) {
    stop("❌ 'phenotype' must have exactly 2 unique values (e.g., case/control).")
  }
  
  # Use phenotype for treatment group
  Y.trt <- as.factor(metadata$phenotype)
  
  # Use batchid for batch group
  Y.bat <- as.factor(metadata$batchid)
  
  # Run PLSDA_batch
  result <- PLSDA_batch(X = taxa_clr, Y.trt = Y.trt, Y.bat = Y.bat,
                        ncomp.trt = 1, ncomp.bat = 5)
  
  write.csv(result$X.nobatch, file.path(output_folder, "normalized_plsda.csv"), row.names = FALSE)
}


if ("ComBat" %in% method_list) {
  cat("Running ComBat...\n")
  require(sva)
  
  mat <- as.matrix(taxa_mat)
  
  pseudo <- min(mat[mat > 0]) * 0.65
  mat[mat == 0] <- pseudo
  
  log_mat <- log2(mat)
  log_mat_t <- t(log_mat)
  
  mod <- model.matrix(~ ., data = covar)
  
  combat_corrected <- ComBat(dat = log_mat_t, batch = batchid, mod = mod,
                             par.prior = FALSE, prior.plots = FALSE)
  combat_corrected_t <- t(combat_corrected)
  
  write.csv(combat_corrected_t, file.path(output_folder, "normalized_combat.csv"), row.names = FALSE)
}

if ("FSQN" %in% method_list) {
  cat("Running FSQN...\n")
  require(FSQN)
  
  normalized <- quantileNormalizeByFeature(taxa_mat, ref_matrix)
  write.csv(normalized, file.path(output_folder, "normalized_fsqn.csv"), row.names = FALSE)
}

if ("MMUPHin" %in% method_list) {
  cat("Running MMUPHin...\n")
  library(MMUPHin)
  
  # convert to integer counts
  count_mat <- round(as.matrix(taxa_mat))     # samples × features, whole numbers
  count_mat[count_mat < 0] <- 0               # just in case
  
  # transpose to features × samples, as required
  taxa_t <- t(count_mat)
  
  metadata$batchid <- factor(metadata$batchid)
  
  fit <- adjust_batch(
    feature_abd = taxa_t,
    batch       = "batchid",
    covariates  = colnames(covar),
    data        = metadata,
    control     = list(verbose = FALSE)
  )
  
  mmuphin_corrected <- t(fit$feature_abd_adj)   # back to samples × features
  write.csv(mmuphin_corrected,file.path(output_folder, "normalized_mmuphin.csv"),row.names = FALSE)
}


if ("CQR" %in% method_list) {
  cat("Running Composite Quantile Regression (CQR)...\n")
  
  # Prepare input objects
  samplebyotu_countmatrix <- as.data.frame(taxa_mat)
  samplebyotu_countmatrix$ID <- rownames(taxa_mat)
  samplebyotu_countmatrix <- samplebyotu_countmatrix[, c("ID", setdiff(colnames(samplebyotu_countmatrix), "ID"))]
  
  indep_vars <- intersect(colnames(metadata), colnames(covar))
  count_matrix_covar <- metadata[, indep_vars, drop = FALSE]
  count_matrix_covar$ID <- rownames(metadata)
  count_matrix_covar <- count_matrix_covar[, c("ID", indep_vars)]
  
  clinicalinfo_total <- data.frame(ID = rownames(metadata), Batch = metadata$batchid)
  
  # Source minimal code directly (no working directory change, assume already loaded or placed in environment)
  source("C:/Users/sunch/Desktop/Project/Batch-Effect/assets/CQR/Background_code.R", local = TRUE)
  source("C:/Users/sunch/Desktop/Project/Batch-Effect/assets/CQR/Data_Preprocessing.R", local = TRUE)
  source("C:/Users/sunch/Desktop/Project/Batch-Effect/assets/CQR/Processing.R", local = TRUE)
  
  # Save output
  write.csv(otu_final, file.path(output_folder, "normalized_cqr.csv"), row.names = FALSE)
}

if ("RUV" %in% method_list) {
  cat("Running fastRUV-III-NB...\n")
  library(ruvIIInb)
  library(DelayedArray)
  
  Y <- t(as.matrix(taxa_mat))  # genes x samples
  
  Y <- round(Y * 100)
  Y <- DelayedArray(Y)
  
  # Use ALL genes as controls (row names)
  control_genes <- rownames(Y)
  
  # Pseudo-replicate matrix
  M <- model.matrix(~ 0 + factor(batchid))
  
  # Numeric batch vector 1..B
  batch_numeric <- as.numeric(factor(batchid))
  
  ruv_result <- fastruvIII.nb(
    Y = Y,
    M = M,
    ctl = control_genes,
    batch = batch_numeric,
    k = 2,
    ncores = 6
  )
  
  corrected_log <- assay(ruv_result, "logPAC")
  
  # convert to a base matrix *after* the transpose
  out_mat <- as.matrix(t(corrected_log))
  
  write.csv(out_mat,file.path(output_folder, "normalized_ruv.csv"), row.names = FALSE)
}

if ("MetaDICT" %in% method_list) {
  cat("Running MetaDICT...\n")
  require(MetaDICT)
  require(vegan)
  
  O <- t(as.matrix(taxa_mat)) 
  meta <- metadata
  meta$batch <- meta$batchid 
  
  O <- O[rowSums(O) > 0, ]
  
  dist_mat <- as.matrix(vegdist(O, method = "bray"))
  
  metadict_res <- MetaDICT(O, meta, distance_matrix = dist_mat)
  
  corrected_mat <- t(metadict_res$count)
  write.csv(corrected_mat, file.path(output_folder, "normalized_metadict.csv"), row.names = FALSE)
}

if ("SVD" %in% method_list) {
  cat("Running SVD-based batch correction...\n")
  
  # CLR-like transform
  clr_matrix <- log2(normalize_tss(taxa_mat) + 1e-6)
  
  # Replace non-finite values with row means
  clr_matrix <- apply(clr_matrix, 1, function(x) {
    x[!is.finite(x)] <- NA
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  })
  clr_matrix <- t(clr_matrix)
  
  # Identify zero-SD features
  feature_sd <- apply(clr_matrix, 2, sd)
  zero_var_cols <- which(feature_sd == 0)
  variable_cols <- setdiff(seq_len(ncol(clr_matrix)), zero_var_cols)
  
  clr_var <- clr_matrix[, variable_cols, drop = FALSE]
  feature_mean <- apply(clr_var, 2, mean)
  feature_sd <- apply(clr_var, 2, sd)
  X_scaled <- scale(clr_var, center = TRUE, scale = TRUE)
  
  # Run SVD
  svd_decomp <- svd(crossprod(X_scaled))
  a1 <- svd_decomp$u[, 1]
  t1 <- X_scaled %*% a1 / sqrt(drop(crossprod(a1)))
  c1 <- crossprod(X_scaled, t1) / drop(crossprod(t1))
  X_deflated <- X_scaled - t1 %*% t(c1)
  
  # Rescale
  X_restored <- matrix(NA, nrow = nrow(X_deflated), ncol = ncol(X_deflated))
  for (i in 1:ncol(X_restored)) {
    X_restored[, i] <- X_deflated[, i] * feature_sd[i] + feature_mean[i]
  }
  
  # Build final output matrix with restored constant features
  full_matrix <- matrix(NA, nrow = nrow(clr_matrix), ncol = ncol(clr_matrix))
  colnames(full_matrix) <- colnames(clr_matrix)
  rownames(full_matrix) <- rownames(clr_matrix)
  
  # Fill variable features
  full_matrix[, variable_cols] <- X_restored
  
  # Fill constant features with original values
  if (length(zero_var_cols) > 0) {
    full_matrix[, zero_var_cols] <- clr_matrix[, zero_var_cols]
  }
  
  write.csv(full_matrix, file.path(output_folder, "normalized_svd.csv"), row.names = FALSE)
}

if ("PN" %in% method_list) {
  cat("Running Percentile Normalisation (PN)...\n")
  
  if (!("phenotype" %in% colnames(metadata))) {
    stop("❌ 'phenotype' column is required in metadata for PN.")
  }
  
  # Check if phenotype is binary
  pheno_vals <- unique(metadata$phenotype)
  if (length(pheno_vals) != 2) {
    stop("❌ 'phenotype' column must have exactly 2 unique values (e.g., case/control).")
  }
  
  # Convert to binary: 0 = control, 1 = case (alphabetically)
  trt <- as.numeric(factor(metadata$phenotype, levels = sort(pheno_vals))) - 1
  
  # TSS normalize
  tss <- normalize_tss(taxa_mat)
  
  # Apply PN
  pn_corrected <- percentile_norm(data = tss, batch = metadata$batchid, trt = trt, ctrl.grp = 0)
  
  write.csv(pn_corrected, file.path(output_folder, "normalized_pn.csv"), row.names = FALSE)
}

if ("FAbatch" %in% method_list) {
  cat("Running FAbatch...\n")
  library(bapred)
  
  if (!("phenotype" %in% colnames(metadata))) {
    stop("❌ 'phenotype' column is required in metadata for PN.")
  }
  
  # Check if phenotype is binary
  pheno_vals <- unique(metadata$phenotype)
  if (length(pheno_vals) != 2) {
    stop("❌ 'phenotype' column must have exactly 2 unique values (e.g., case/control).")
  }
  
  # Prepare inputs
  x <- as.matrix(log2(normalize_tss(taxa_mat) + 1e-6))
  y <- factor(metadata$phenotype, levels = sort(pheno_vals))
  batch <- factor(metadata$batchid)
  
  # Run fabatch (factors estimated automatically)
  fa_out <- fabatch(x = x, y = y, batch = batch)
  x_adj <- fa_out$xadj
  
  colnames(x_adj) <- colnames(taxa_mat)

  write.csv(x_adj, file.path(output_folder, "normalized_fabatch.csv"), row.names = FALSE)
}

if ("ComBatSeq" %in% method_list) {
  cat("Running ComBat-Seq...\n")
  require(sva)
  if (!("phenotype" %in% colnames(metadata))) {
    stop("❌ 'phenotype' column is required in metadata for ComBat-Seq.")
  }
  
  # Prepare raw count matrix
  count_matrix <- round(as.matrix(taxa_mat) * 100)  # ComBat-Seq requires integers
  
  # Check and remove samples with zero library size
  library_sizes <- rowSums(count_matrix)
  nonzero_samples <- library_sizes > 0
  if (any(!nonzero_samples)) {
    cat("⚠️ Removing", sum(!nonzero_samples), "samples with zero total counts.\n")
    cat("Removed sample(s):", paste(rownames(count_matrix)[!nonzero_samples], collapse = ", "), "\n")
    count_matrix <- count_matrix[nonzero_samples, , drop = FALSE]
    metadata <- metadata[nonzero_samples, , drop = FALSE]
  }
  
  # Ensure matching rownames
  if (!all(rownames(count_matrix) == rownames(metadata))) {
    stop("❌ Sample IDs in count matrix and metadata do not match after filtering.")
  }
  
  # Run ComBat-Seq
  combat_seq <- ComBat_seq(counts = t(count_matrix),  # genes x samples
                           batch = metadata$batchid,
                           group = metadata$phenotype)
  
  # Transpose back to samples x features
  combat_seq <- t(combat_seq)
  rownames(combat_seq) <- rownames(count_matrix)
  colnames(combat_seq) <- colnames(count_matrix)
  
  write.csv(combat_seq, file.path(output_folder, "normalized_combatseq.csv"), row.names = FALSE)
}
