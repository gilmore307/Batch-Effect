# ---------------------------
# Handle Arguments
# ---------------------------
#args <- commandArgs(trailingOnly = TRUE)
args <- c("BMC,ComBat,ConQuR,PLSDA,ALRA", "output/example")

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
library(ConQuR)
library(foreach)
library(doParallel)
library(dplyr)
library(Matrix)
library(SummarizedExperiment)
library(S4Vectors)
library(compositions)

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

# Ensure required metadata
required_columns <- c("sample_id", "batchid")
missing <- setdiff(required_columns, colnames(metadata))
if (length(missing) > 0) {
  stop(paste("Missing required metadata columns:", paste(missing, collapse = ", ")))
}
covar <- metadata[, !(colnames(metadata) %in% c("sample_id", "batchid")), drop = FALSE]
batchid <- factor(metadata$batchid)

normalize_tss <- function(mat) {
  mat <- sweep(mat, 1, rowSums(mat), "/")
  mat[is.na(mat)] <- 0
  return(mat)
}

# ---------------------------
# Normalize Methods
# ---------------------------

# TSS
if ("TSS" %in% method_list) {
  cat("Running TSS...\n")
  norm <- normalize_tss(taxa_mat)
  write.csv(norm, file.path(output_folder, "normalized_tss.csv"), row.names = FALSE)
}

# UQ
if ("UQ" %in% method_list) {
  cat("Running UQ...\n")
  uq_norm <- apply(taxa_mat, 1, function(x) x / quantile(x[x > 0], 0.75))
  write.csv(t(uq_norm), file.path(output_folder, "normalized_uq.csv"), row.names = FALSE)
}

# MED
if ("MED" %in% method_list) {
  cat("Running MED...\n")
  med_norm <- apply(taxa_mat, 1, function(x) x / median(x[x > 0]))
  write.csv(t(med_norm), file.path(output_folder, "normalized_med.csv"), row.names = FALSE)
}

# CSS
if ("CSS" %in% method_list) {
  cat("Running CSS...\n")
  css_obj <- cumNorm(newMRexperiment(t(taxa_mat)))
  css_norm <- t(MRcounts(css_obj, norm = TRUE))
  write.csv(css_norm, file.path(output_folder, "normalized_css.csv"), row.names = FALSE)
}

# TMM
if ("TMM" %in% method_list) {
  cat("Running TMM...\n")
  dge <- DGEList(counts = t(taxa_mat))
  dge <- calcNormFactors(dge, method = "TMM")
  tmm_norm <- t(cpm(dge))
  write.csv(tmm_norm, file.path(output_folder, "normalized_tmm.csv"), row.names = FALSE)
}

# GMPR
if ("GMPR" %in% method_list) {
  cat("Running GMPR...\n")
  gmpr_factor <- GMPR(taxa_mat)
  
  if (any(is.na(gmpr_factor))) {
    cat("⚠️ Warning: NA values detected in GMPR size factors — replacing with median\n")
    gmpr_factor[is.na(gmpr_factor)] <- median(gmpr_factor, na.rm = TRUE)
  }
  
  gmpr_norm <- sweep(taxa_mat, 1, gmpr_factor, "/")
  write.csv(gmpr_norm, file.path(output_folder, "normalized_gmpr.csv"), row.names = FALSE)
}


# CLR
if ("CLR" %in% method_list) {
  cat("Running CLR...\n")
  norm <- normalize_tss(taxa_mat)
  norm[norm == 0] <- min(norm[norm != 0]) * 0.65
  clr_norm <- apply(norm, 1, clr)
  write.csv(t(clr_norm), file.path(output_folder, "normalized_clr.csv"), row.names = FALSE)
}

# LOG
if ("LOG" %in% method_list) {
  cat("Running LOG...\n")
  norm <- normalize_tss(taxa_mat)
  norm[norm == 0] <- min(norm[norm != 0]) * 0.65
  log_norm <- log(norm)
  write.csv(log_norm, file.path(output_folder, "normalized_log.csv"), row.names = FALSE)
}

# AST
if ("AST" %in% method_list) {
  cat("Running AST...\n")
  norm <- normalize_tss(taxa_mat)
  ast_norm <- asin(sqrt(norm))
  write.csv(ast_norm, file.path(output_folder, "normalized_ast.csv"), row.names = FALSE)
}

# STD
if ("STD" %in% method_list) {
  cat("Running STD...\n")
  norm <- normalize_tss(taxa_mat)
  std_norm <- t(apply(norm, 1, function(x) (x - mean(x)) / sd(x)))
  write.csv(std_norm, file.path(output_folder, "normalized_std.csv"), row.names = FALSE)
}

# RANK
if ("Rank" %in% method_list) {
  cat("Running Rank...\n")
  norm <- normalize_tss(taxa_mat)
  rank_norm <- t(apply(norm, 1, rank))
  write.csv(rank_norm, file.path(output_folder, "normalized_rank.csv"), row.names = FALSE)
}

# BLOM
if ("Blom" %in% method_list) {
  cat("Running Blom...\n")
  norm <- normalize_tss(taxa_mat)
  c <- 3/8
  blom_norm <- t(apply(norm, 1, function(x) qnorm((rank(x) - c)/(length(x) - 2*c + 1))))
  write.csv(blom_norm, file.path(output_folder, "normalized_blom.csv"), row.names = FALSE)
}

# NPN
if ("NPN" %in% method_list) {
  cat("Running NPN...\n")
  norm <- normalize_tss(taxa_mat)
  npn_norm <- t(huge.npn(t(norm), npn.func = "truncation"))
  write.csv(npn_norm, file.path(output_folder, "normalized_npn.csv"), row.names = FALSE)
}

# logCPM
if ("logCPM" %in% method_list) {
  cat("Running logCPM...\n")
  taxa_mat[taxa_mat == 0] <- 1
  logcpm_norm <- t(cpm(t(taxa_mat), log = TRUE))
  write.csv(logcpm_norm, file.path(output_folder, "normalized_logcpm.csv"), row.names = FALSE)
}

# VST
if ("VST" %in% method_list) {
  cat("Running VST...\n")
  taxa_mat[taxa_mat == 0] <- 1
  
  # Transpose for DESeq2 (samples in columns, features in rows)
  count_matrix <- round(t(taxa_mat))  # must be integers
  
  # Make sure sample names match
  sample_ids <- colnames(count_matrix)
  col_data <- data.frame(condition = batchid)
  rownames(col_data) <- sample_ids
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~condition)
  
  # Apply VST
  vst_mat <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
  write.csv(t(vst_mat), file.path(output_folder, "normalized_vst.csv"), row.names = FALSE)
}


# QN
if ("QN" %in% method_list) {
  cat("Running QN...\n")
  norm <- normalize_tss(taxa_mat)
  norm[norm == 0] <- min(norm[norm != 0]) * 0.65
  log_mat <- log(norm)
  ref_quantiles <- normalize.quantiles.determine.target(as.matrix(t(log_mat)))
  qn_corrected <- normalize.quantiles.use.target(as.matrix(t(log_mat)), target = ref_quantiles)
  dimnames(qn_corrected) <- dimnames(t(log_mat))
  write.csv(t(qn_corrected), file.path(output_folder, "normalized_qn.csv"), row.names = FALSE)
}

# BMC
if ("BMC" %in% method_list) {
  cat("Running BMC...\n")
  norm <- normalize_tss(taxa_mat)
  norm[norm == 0] <- min(norm[norm != 0]) * 0.65
  log_mat <- log(norm)
  batch_factor <- factor(batchid)
  pam_input <- list(x = as.matrix(t(log_mat)), batchlabels = batch_factor)
  bmc_corrected <- pamr.batchadjust(pam_input)$x
  write.csv(t(bmc_corrected), file.path(output_folder, "normalized_bmc.csv"), row.names = FALSE)
}

# limma
if ("limma" %in% method_list) {
  cat("Running limma...\n")
  norm <- normalize_tss(taxa_mat)
  norm[norm == 0] <- min(norm[norm != 0]) * 0.65
  log_mat <- log(norm)
  batch_factor <- factor(batchid)
  limma_corrected <- removeBatchEffect(t(log_mat), batch = batch_factor)
  write.csv(t(limma_corrected), file.path(output_folder, "normalized_limma.csv"), row.names = FALSE)
}

if ("ConQuR" %in% method_list) {
  cat("Running ConQuR with parallel backend...\n")
  
  # Register parallel backend
  num_cores <- parallel::detectCores(logical = TRUE) - 1
  registerDoParallel(cores = num_cores)
  cat(paste0("Using ", num_cores, " cores\n"))
  
  # Run ConQuR
  conqur_result <- ConQuR(
    tax_tab = taxa_mat,
    batchid = batchid,
    covariates = covar,
    batch_ref = "0",
    simple_match = TRUE  # this reduces complexity
  )
  
  write.csv(conqur_result, file.path(output_folder, "normalized_conqur.csv"), row.names = FALSE)
  
  # Optional: shut down cluster afterward
  stopImplicitCluster()
}


if ("ALRA" %in% method_list) {
  cat("Running ALRA...\n")
  taxa_filtered <- taxa_mat[rowSums(taxa_mat) > 0, ]
  taxa_norm <- sweep(taxa_filtered, 1, rowSums(taxa_filtered), "/")
  taxa_norm[is.na(taxa_norm)] <- 0
  taxa_norm_t <- t(taxa_norm)
  alra_result <- alra(taxa_norm_t)
  taxa_alra <- t(alra_result$A_norm_rank_k_cor_sc)
  colnames(taxa_alra) <- colnames(taxa_mat)
  write.csv(taxa_alra, file.path(output_folder, "normalized_alra.csv"), row.names = FALSE)
}

if ("PLSDA" %in% method_list) {
  print("Running PLSDAbatch...")
  rownames(taxa_mat) <- rownames(metadata)
  taxa_clr <- clr(taxa_mat + 1)
  
  Y.trt <- ifelse(metadata$batchid == 0, "control", "treatment")
  Y.bat <- as.factor(metadata$batchid)
  
  result <- PLSDA_batch(X = taxa_clr, Y.trt = Y.trt, Y.bat = Y.bat, ncomp.trt = 1, ncomp.bat = 5)
  write.csv(result$X.nobatch, file.path(output_folder, "normalized_plsda.csv"), row.names = FALSE)
}


if ("ComBat" %in% method_list) {
  cat("Running ComBat...\n")
  taxa_rel_abund <- sweep(taxa_mat, 1, rowSums(taxa_mat), "/")
  taxa_rel_abund[is.na(taxa_rel_abund)] <- 0
  taxa_log <- log2(taxa_rel_abund + 1e-6)
  taxa_log_t <- t(taxa_log)
  combat_corrected <- ComBat(dat = taxa_log_t, batch = batchid, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
  scaled_data <- scale(t(combat_corrected))
  taxa_rel_corrected <- 2^scaled_data
  write.csv(taxa_rel_corrected, file.path(output_folder, "normalized_combat.csv"), row.names = FALSE)
}

if ("FSQN" %in% method_list) {
  cat("Running FSQN...\n")
  
  reference_batch <- 0  # explicitly use batchid == 0
  ref_matrix <- taxa_mat[batchid == reference_batch, , drop = FALSE]
  reference_distribution <- apply(ref_matrix, 2, function(x) sort(x))
  reference_target <- apply(reference_distribution, 1, median)
  
  normalize_to_reference <- function(mat, target_vector) {
    apply(mat, 2, function(feature_column) {
      ranked <- rank(feature_column, ties.method = "min")
      normalized <- target_vector[ranked]
      return(normalized)
    })
  }
  
  normalized_list <- lapply(unique(batchid), function(batch) {
    mat <- taxa_mat[batchid == batch, , drop = FALSE]
    normed <- normalize_to_reference(mat, reference_target)
    rownames(normed) <- rownames(mat)
    return(normed)
  })
  
  normalized_matrix <- do.call(rbind, normalized_list)
  normalized_matrix <- normalized_matrix[rownames(taxa_mat), ]
  write.csv(normalized_matrix, file.path(output_folder, "normalized_fsqn.csv"), row.names = FALSE)
}

