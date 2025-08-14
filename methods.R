# ===========================
# normalize_all_methods.R
# ===========================

# ---------------------------
# Handle Arguments
# ---------------------------
# args <- commandArgs(trailingOnly = TRUE)
args <- c("FSQN", "output/example")

if (length(args) < 2) stop("Usage: Rscript normalize_all_methods.R <method_list> <output_folder>")

method_list   <- unlist(strsplit(args[1], ","))
output_folder <- args[2]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

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
# Console Logger & Guards
# ---------------------------
.ts <- function() format(Sys.time(), "%H:%M:%S")
say <- function(...){ cat(sprintf("[%s] ", .ts()), paste0(..., collapse=""), "\n") }

start_step <- function(name){ say("▶️  START: ", name); proc.time() }
ok_step    <- function(name, t0){ dt <- proc.time()-t0; say("✅ DONE:  ", name, sprintf(" (%.2fs)", dt["elapsed"])) }
warn_step  <- function(name, msg){ say("⚠️ WARN:  ", name, " — ", msg) }
fail_step  <- function(name, msg){ say("❌ FAIL:  ", name, " — ", msg); stop(paste0(name, ": ", msg)) }

check_table <- function(x, name = "table", allow_negative = TRUE) {
  if (!is.data.frame(x) && !is.matrix(x)) fail_step(name, "Not a data.frame/matrix.")
  DF <- as.data.frame(x, stringsAsFactors = FALSE)
  rn <- rownames(DF); cn <- colnames(DF)
  
  num_cols <- vapply(DF, is.numeric, TRUE)
  if (any(num_cols)) {
    Mnum <- as.matrix(DF[num_cols])
    bad_fin <- which(!is.finite(Mnum), arr.ind = TRUE)
    if (nrow(bad_fin) > 0) {
      cells <- apply(bad_fin, 1, function(rc)
        sprintf("(%s, %s)", if(!is.null(rn)) rn[rc[1]] else rc[1],
                if(!is.null(cn)) cn[which(num_cols)][rc[2]] else rc[2]))
      fail_step(name, paste0("NA/NaN/Inf at ", paste(head(cells,10), collapse=", "),
                             if(nrow(bad_fin)>10) sprintf(" ... and %d more", nrow(bad_fin)-10) else ""))
    }
    if (!allow_negative) {
      negs <- which(Mnum < 0, arr.ind = TRUE)
      if (nrow(negs) > 0) {
        cells <- apply(negs, 1, function(rc)
          sprintf("(%s, %s)", if(!is.null(rn)) rn[rc[1]] else rc[1],
                  if(!is.null(cn)) cn[which(num_cols)][rc[2]] else rc[2]))
        fail_step(name, paste0("Negative values at ", paste(head(cells,10), collapse=", "),
                               if(nrow(negs)>10) sprintf(" ... and %d more", nrow(negs)-10) else ""))
      }
    }
  }
  if (any(!num_cols)) {
    Mchr <- DF[!num_cols]
    na_idx <- which(is.na(Mchr), arr.ind = TRUE)
    if (nrow(na_idx) > 0) {
      cells <- apply(na_idx, 1, function(rc)
        sprintf("(%s, %s)", if(!is.null(rn)) rn[rc[1]] else rc[1],
                if(!is.null(cn)) cn[which(!num_cols)][rc[2]] else rc[2]))
      fail_step(name, paste0("Missing values at ", paste(head(cells,10), collapse=", "),
                             if(nrow(na_idx)>10) sprintf(" ... and %d more", nrow(na_idx)-10) else ""))
    }
  }
  if (nrow(DF) < 2 || ncol(DF) < 2) fail_step(name, "Too few rows/cols.")
  invisible(TRUE)
}

write_out <- function(obj, path){
  nm <- basename(path)
  t0 <- start_step(paste0("Write ", nm))
  write.csv(obj, path, row.names = FALSE)
  ok_step(paste0("Write ", nm), t0)
}

# Simpler run_method: we LOG warnings but DO NOT muffle globally
run_method <- function(name, expr){
  t0 <- start_step(name)
  out <- withCallingHandlers(
    tryCatch({
      val <- force(expr)
      ok_step(name, t0)
      val
    }, error = function(e){
      fail_step(name, conditionMessage(e))
    }),
    warning = function(w){
      warn_step(name, conditionMessage(w))
    }
  )
  invisible(out)
}

post_summary <- function(M, name){
  M <- as.matrix(M)
  say(sprintf("%s summary: min=%.5f max=%.5f neg=%.2f%% NA=%d",
              name, min(M, na.rm=TRUE), max(M, na.rm=TRUE),
              100*mean(M<0, na.rm=TRUE), sum(!is.finite(M))))
}

normalize_tss <- function(mat){
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  out <- sweep(mat, 1, row_sums, "/")
  out[is.na(out)] <- 0
  out
}

# ---------------------------
# Load Data
# ---------------------------
custom_matrix_path    <- file.path(output_folder, "raw.csv")
custom_metadata_path  <- file.path(output_folder, "metadata.csv")
default_matrix_path   <- file.path("assets", "raw.csv")
default_metadata_path <- file.path("assets", "metadata.csv")

if (file.exists(custom_matrix_path) && file.exists(custom_metadata_path)) {
  say("✅ Using uploaded user files")
  uploaded_mat <- read.csv(custom_matrix_path, check.names = FALSE)
  metadata     <- read.csv(custom_metadata_path)
} else {
  say("⚠️ No uploaded files found — using default assets")
  uploaded_mat <- read.csv(default_matrix_path, check.names = FALSE)
  metadata     <- read.csv(default_metadata_path)
}

# ---------------------------
# Ensure INPUT is integer counts, and write origin.csv + raw.csv
# ---------------------------
source("Relative_to_Raw.R")
taxa_mat <- check_and_convert_counts(
  uploaded_mat,
  fixed_library_size   = 10000,
  force_row_normalize  = TRUE,
  output_folder        = output_folder,
  quiet = FALSE
)

# Attach rownames (expects metadata$sample_id)
rownames(taxa_mat) <- metadata$sample_id
rownames(metadata) <- metadata$sample_id

# Keep batchid as factor (0 allowed)
metadata$batchid <- factor(metadata$batchid, levels = unique(metadata$batchid))
batchid <- metadata$batchid

# Covariates cleanup
covar <- metadata[, !(colnames(metadata) %in% c("sample_id", "batchid","phenotype")), drop = FALSE]
covar <- as.data.frame(lapply(covar, function(col) {
  if (is.numeric(col))      col[is.na(col)] <- mean(col, na.rm = TRUE)
  else if (is.factor(col))  if (anyNA(col)) col[is.na(col)] <- levels(col)[1]
  else                      col[is.na(col)] <- as.character(stats::na.omit(col))[1]
  col
}))

# Reference batch for QN/FSQN
reference_batch <- 0
ref_idx <- which(batchid == reference_batch)
ref_matrix <- taxa_mat[ref_idx, , drop = FALSE]

# ---------------------------
# Global input checks
# ---------------------------
run_method("Input checks", {
  check_table(taxa_mat, "taxa_mat", allow_negative = FALSE)  # counts enforced
  check_table(metadata, "metadata", allow_negative = TRUE)
  if (!identical(rownames(taxa_mat), rownames(metadata))) {
    fail_step("Alignment", "rownames(taxa_mat) != rownames(metadata).")
  }
  say("Method list: ", paste(method_list, collapse=", "))
  say("Samples: ", nrow(taxa_mat), " | Features: ", ncol(taxa_mat))
  say("Reference batch level: ", as.character(reference_batch),
      " | #ref samples: ", length(ref_idx))
})

# ---------------------------
# Methods (no _pos.csv, CQR removed)
# ---------------------------

# QN
if ("QN" %in% method_list) {
  run_method("QN", {
    require(preprocessCore)
    all_abs <- as.matrix(taxa_mat)
    ref_abs <- as.matrix(ref_matrix)
    if (nrow(ref_abs) < 2) warn_step("QN", "Reference batch has <2 samples; results may be unstable.")
    ref_target <- normalize.quantiles.determine.target(ref_abs)
    qn_corrected <- normalize.quantiles.use.target(all_abs, target = ref_target)
    dimnames(qn_corrected) <- dimnames(all_abs)
    post_summary(qn_corrected, "QN")
    write_out(qn_corrected, file.path(output_folder, "normalized_qn.csv"))
  })
}

# BMC
if ("BMC" %in% method_list) {
  run_method("BMC", {
    suppressPackageStartupMessages(suppressWarnings(require(pamr)))
    
    # TSS first (rows sum to 1)
    norm <- normalize_tss(taxa_mat)
    if (all(norm == 0)) fail_step("BMC", "All-zero after TSS.")
    
    # tiny pseudo for zeros before log
    zmin <- min(norm[norm > 0])
    if (!is.finite(zmin)) fail_step("BMC", "No non-zero entries for log().")
    norm[norm == 0] <- zmin * 0.65
    
    # log-scale for pamr.batchadjust
    log_mat <- log(norm)
    pam_input <- list(x = as.matrix(t(log_mat)), batchlabels = factor(batchid))
    
    adj_log_t <- pamr.batchadjust(pam_input)$x      # features x samples (log-scale)
    adj_log   <- t(adj_log_t)                        # back to samples x features
    
    # ---- make it non-negative for distances ----
    # inverse of log() is exp(); clip numerical fuzz; re-TSS
    adj_pos <- exp(adj_log)
    adj_pos[!is.finite(adj_pos)] <- 0
    adj_pos[adj_pos < 0] <- 0
    adj_pos <- normalize_tss(adj_pos)
    
    post_summary(adj_pos, "BMC (positive)")
    write_out(adj_pos, file.path(output_folder, "normalized_bmc.csv"))
  })
}

# limma
if ("limma" %in% method_list) {
  run_method("limma", {
    suppressPackageStartupMessages(suppressWarnings(require(limma)))
    
    # 1) Pseudocount and log
    mat <- as.matrix(taxa_mat)
    nz <- mat[mat > 0]
    if (!length(nz)) fail_step("limma", "No positive entries to define pseudo.")
    pseudo <- min(nz) * 0.65
    mat[mat == 0] <- pseudo
    
    log_mat <- log(mat)  # natural log is fine; base doesn't matter if you invert with exp()
    
    # 2) Batch removal on log scale
    limma_corrected_log_t <- removeBatchEffect(
      t(log_mat),
      batch = factor(batchid),
      covariates = if (ncol(covar) > 0) as.matrix(covar) else NULL
    )
    limma_corrected_log <- t(limma_corrected_log_t)
    
    # 3) Back to positive space for distances / downstream tools
    limma_pos <- exp(limma_corrected_log)
    limma_pos[!is.finite(limma_pos)] <- 0
    limma_pos[limma_pos < 0] <- 0
    limma_pos <- normalize_tss(limma_pos)  # rows sum to 1
    
    post_summary(limma_pos, "limma (positive)")
    write_out(limma_pos, file.path(output_folder, "normalized_limma.csv"))
  })
}

# ConQuR
if ("ConQuR" %in% method_list) {
  run_method("ConQuR", {
    suppressPackageStartupMessages(suppressWarnings({
      library(ConQuR); library(doParallel)
    }))
    num_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
    say("ConQuR cores: ", num_cores)
    batchidF <- as.factor(metadata$batchid)
    batch_ref <- as.character(reference_batch)
    covariates <- metadata[, colnames(covar), drop = FALSE]
    rownames(covariates) <- NULL
    
    res <- withCallingHandlers(
      ConQuR(
        tax_tab = taxa_mat,
        batchid = batchidF,
        covariates = covariates,
        batch_ref = batch_ref,
        logistic_lasso = FALSE,
        quantile_type = "standard",
        simple_match = FALSE,
        lambda_quantile = "2p/n",
        interplt = FALSE,
        delta = 0.4999,
        taus = seq(0.05, 0.95, by = 0.05),
        num_core = num_cores
      ),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("Solution may be nonunique", msg) ||
            grepl("^glm\\.fit:", msg)) {
          # silence these recurring warnings
          invokeRestart("muffleWarning")
        } else {
          warn_step("ConQuR", msg)
        }
      }
    )
    
    post_summary(res, "ConQuR")
    write_out(res, file.path(output_folder, "normalized_conqur.csv"))
  })
}

# PLSDA
if ("PLSDA" %in% method_list) {
  run_method("PLSDAbatch", {
    suppressPackageStartupMessages(suppressWarnings({
      library(PLSDAbatch)
      library(compositions)
    }))
    
    if (!("phenotype" %in% colnames(metadata)))
      fail_step("PLSDAbatch", "'phenotype' not found.")
    if (length(unique(metadata$phenotype)) != 2)
      fail_step("PLSDAbatch", "'phenotype' must be binary.")
    
    # CLR transform
    taxa_clr <- compositions::clr(taxa_mat + 1)
    
    # Run PLSDA batch correction
    Y.trt <- as.factor(metadata$phenotype)
    Y.bat <- as.factor(metadata$batchid)
    result <- PLSDA_batch(
      X = taxa_clr,
      Y.trt = Y.trt,
      Y.bat = Y.bat,
      ncomp.trt = 1,
      ncomp.bat = 5
    )
    
    # Corrected matrix (still in CLR space)
    out <- result$X.nobatch
    
    # ---- Back-transform to positive ----
    out_pos <- compositions::clrInv(out)  # inverse CLR → positive ratios
    out_pos[!is.finite(out_pos)] <- 0
    out_pos[out_pos < 0] <- 0
    out_pos <- normalize_tss(out_pos)     # row sums = 1
    
    post_summary(out_pos, "PLSDAbatch (positive)")
    write_out(out_pos, file.path(output_folder, "normalized_plsda.csv"))
  })
}

# ComBat-Seq
if ("ComBatSeq" %in% method_list) {
  run_method("ComBat-Seq", {
    require(sva)
    if (!("phenotype" %in% colnames(metadata))) fail_step("ComBat-Seq", "'phenotype' is required.")
    count_matrix <- round(as.matrix(taxa_mat) * 100)
    library_sizes <- rowSums(count_matrix)
    nonzero_samples <- library_sizes > 0
    if (any(!nonzero_samples)) {
      say("ComBat-Seq: removing ", sum(!nonzero_samples), " samples with zero library size: ",
          paste(rownames(count_matrix)[!nonzero_samples], collapse = ", "))
      count_matrix <- count_matrix[nonzero_samples, , drop = FALSE]
      metadata <- metadata[nonzero_samples, , drop = FALSE]
    }
    if (!all(rownames(count_matrix) == rownames(metadata))) {
      fail_step("ComBat-Seq", "Sample IDs mismatch after filtering.")
    }
    combat_seq <- ComBat_seq(counts = t(count_matrix),  # genes x samples
                             batch = metadata$batchid,
                             group = metadata$phenotype)
    out <- t(combat_seq)
    rownames(out) <- rownames(count_matrix)
    colnames(out) <- colnames(count_matrix)
    post_summary(out, "ComBat-Seq")
    write_out(out, file.path(output_folder, "normalized_combatseq.csv"))
  })
}

# ComBat
if ("ComBat" %in% method_list) {
  run_method("ComBat", {
    suppressPackageStartupMessages(suppressWarnings(require(sva)))
    
    # 1) Add pseudocount
    mat <- as.matrix(taxa_mat)
    nz <- mat[mat > 0]
    if (!length(nz))
      fail_step("ComBat", "No positive entries to define pseudocount.")
    pseudo <- min(nz) * 0.65
    mat[mat == 0] <- pseudo
    
    # 2) Log-transform (base 2)
    log_mat_t <- t(log2(mat))
    
    # 3) Build model
    mod <- if (ncol(covar) > 0) model.matrix(~ ., data = covar) else NULL
    
    # 4) Run ComBat
    combat_log_t <- ComBat(
      dat = log_mat_t,
      batch = batchid,
      mod = mod,
      par.prior = FALSE,
      prior.plots = FALSE
    )
    combat_log <- t(combat_log_t)
    
    # 5) Back-transform to positive space
    combat_pos <- 2^combat_log
    combat_pos[!is.finite(combat_pos)] <- 0
    combat_pos[combat_pos < 0] <- 0
    combat_pos <- normalize_tss(combat_pos)  # row sums = 1
    
    post_summary(combat_pos, "ComBat (positive)")
    write_out(combat_pos, file.path(output_folder, "normalized_combat.csv"))
  })
}

# FSQN
if ("FSQN" %in% method_list) {
  run_method("FSQN", {
    require(FSQN)
    normalized <- quantileNormalizeByFeature(taxa_mat, ref_matrix)
    post_summary(normalized, "FSQN")
    write_out(normalized, file.path(output_folder, "normalized_fsqn.csv"))
  })
}

# MMUPHin
if ("MMUPHin" %in% method_list) {
  run_method("MMUPHin", {
    require(MMUPHin)
    count_mat <- round(as.matrix(taxa_mat))  # counts guard
    count_mat[count_mat < 0] <- 0
    taxa_t <- t(count_mat)
    metadata$batchid <- factor(metadata$batchid)
    fit <- adjust_batch(
      feature_abd = taxa_t,
      batch       = "batchid",
      covariates  = colnames(covar),
      data        = metadata,
      control     = list(verbose = FALSE,
                         diagnostic_plot = NULL)
    )
    out <- t(fit$feature_abd_adj)
    post_summary(out, "MMUPHin")
    write_out(out, file.path(output_folder, "normalized_mmuphin.csv"))
  })
}

# RUV (fastRUV-III-NB) — robust ctl handling
if ("RUV" %in% method_list) {
  run_method("fastRUV-III-NB", {
    suppressPackageStartupMessages(suppressWarnings({
      library(ruvIIInb)
    }))
    
    # start clean
    try(closeAllConnections(), silent = TRUE)
    try({ while (dev.cur() > 1) dev.off() }, silent = TRUE)
    invisible(gc())
    
    # genes x samples (counts already prepared upstream)
    Y <- t(as.matrix(taxa_mat))
    
    # enforce sample order and drop obvious all-zero genes
    stopifnot(!is.null(colnames(Y)))
    samp_ids <- colnames(Y)
    if (is.null(rownames(metadata)) || !all(samp_ids %in% rownames(metadata))) {
      fail_step("RUV", "metadata rownames must be sample IDs matching colnames(taxa_mat).")
    }
    keep_genes <- rowSums(Y) > 0
    if (!any(keep_genes)) fail_step("RUV", "All genes have zero counts.")
    Y <- Y[keep_genes, samp_ids, drop = FALSE]
    
    # ---- ctl as ALL feature names (so ruvIIInb can remap after its own filtering)
    ctl_names <- rownames(Y)
    
    # design matrices in the exact Y sample order
    batch_factor  <- factor(metadata[samp_ids, "batchid"])
    M             <- model.matrix(~ 0 + batch_factor)
    rownames(M)   <- samp_ids
    batch_numeric <- as.numeric(batch_factor)
    
    # sanity checks
    stopifnot(identical(rownames(M), samp_ids))
    stopifnot(length(batch_numeric) == ncol(Y))
    
    # modest parallelism
    workers <- max(1L, min(4L, parallel::detectCores(logical = TRUE) - 1L))
    
    ruv_result <- fastruvIII.nb(
      Y                = as.matrix(Y),  # in-memory
      M                = M,             # samples x groups
      ctl              = ctl_names,     # character vector (safe across internal filtering)
      batch            = batch_numeric, # length == ncol(Y)
      k                = 2,
      pCells.touse     = 0.05,          # 5% subsample for alpha (speed)
      use.pseudosample = FALSE,
      ncores           = workers
    )
    
    corrected_log <- assay(ruv_result, "logPAC")
    out <- as.matrix(t(corrected_log))  # back to samples x features
    
    post_summary(out, "RUV-III-NB")
    write_out(out, file.path(output_folder, "normalized_ruv.csv"))
    
    # clean after
    try(closeAllConnections(), silent = TRUE)
    try({ while (dev.cur() > 1) dev.off() }, silent = TRUE)
    invisible(gc())
  })
}

# PN
if ("PN" %in% method_list) {
  run_method("PN", {
    if (!("phenotype" %in% colnames(metadata))) fail_step("PN", "'phenotype' is required.")
    pheno_vals <- unique(metadata$phenotype)
    if (length(pheno_vals) != 2) fail_step("PN", "'phenotype' must be binary.")
    trt <- as.numeric(factor(metadata$phenotype, levels = sort(pheno_vals))) - 1
    tss <- normalize_tss(taxa_mat)
    if (all(tss == 0)) fail_step("PN", "All zero after TSS.")
    pn_corrected <- percentile_norm(data = tss, batch = metadata$batchid, trt = trt, ctrl.grp = 0)
    post_summary(pn_corrected, "PN")
    write_out(pn_corrected, file.path(output_folder, "normalized_pn.csv"))
  })
}

# MetaDICT
if ("MetaDICT" %in% method_list) {
  run_method("MetaDICT", {
    suppressPackageStartupMessages(suppressWarnings({
      library(MetaDICT); library(vegan)
    }))
    O <- t(as.matrix(taxa_mat))
    meta <- metadata
    meta$batch <- meta$batchid
    O <- O[rowSums(O) > 0, , drop=FALSE]
    if (nrow(O) < 2) fail_step("MetaDICT", "Too few non-zero samples after filtering.")
    dist_mat <- as.matrix(vegdist(O, method = "bray"))
    metadict_res <- MetaDICT(O, meta, distance_matrix = dist_mat)
    corrected_mat <- t(metadict_res$count)
    post_summary(corrected_mat, "MetaDICT")
    write_out(corrected_mat, file.path(output_folder, "normalized_metadict.csv"))
  })
}

# SVD
if ("SVD" %in% method_list) {
  run_method("SVD", {
    cat("Running SVD-based batch correction...\n")
    
    # CLR-like transform
    clr_matrix <- log2(normalize_tss(taxa_mat) + 1e-6)
    
    # Replace non-finite values with row means
    clr_matrix <- t(apply(clr_matrix, 1, function(x) {
      x[!is.finite(x)] <- NA
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    }))
    
    # Identify zero-SD features
    feature_sd <- apply(clr_matrix, 2, sd)
    zero_var_cols <- which(feature_sd == 0)
    variable_cols <- setdiff(seq_len(ncol(clr_matrix)), zero_var_cols)
    
    if (!length(variable_cols)) fail_step("SVD", "All features zero variance; nothing to correct.")
    
    # Center and scale variable features
    clr_var <- clr_matrix[, variable_cols, drop = FALSE]
    feature_mean <- apply(clr_var, 2, mean)
    feature_sd <- apply(clr_var, 2, sd)
    X_scaled <- scale(clr_var, center = TRUE, scale = TRUE)
    
    # Run SVD on covariance
    svd_decomp <- svd(crossprod(X_scaled))
    a1 <- svd_decomp$u[, 1]
    t1 <- X_scaled %*% a1 / sqrt(drop(crossprod(a1)))
    c1 <- crossprod(X_scaled, t1) / drop(crossprod(t1))
    X_deflated <- X_scaled - t1 %*% t(c1)
    
    # Rescale back to original mean/sd
    X_restored <- sweep(X_deflated, 2, feature_sd, `*`)
    X_restored <- sweep(X_restored, 2, feature_mean, `+`)
    
    # Build final matrix with restored constant features
    full_matrix <- matrix(NA, nrow = nrow(clr_matrix), ncol = ncol(clr_matrix))
    colnames(full_matrix) <- colnames(clr_matrix)
    rownames(full_matrix) <- rownames(clr_matrix)
    full_matrix[, variable_cols] <- X_restored
    if (length(zero_var_cols) > 0) {
      full_matrix[, zero_var_cols] <- clr_matrix[, zero_var_cols]
    }
    
    # ---- Back-transform to positive space ----
    positive_matrix <- 2^full_matrix - 1e-6          # still a matrix
    positive_matrix[positive_matrix < 0] <- 0        # clamp negatives to 0
    
    # TSS-normalize rows
    positive_matrix <- normalize_tss(positive_matrix)
    
    post_summary(positive_matrix, "SVD")
    write_out(positive_matrix, file.path(output_folder, "normalized_svd.csv"))
  })
}

# FAbatch (simple & robust)
if ("FAbatch" %in% method_list) {
  run_method("FAbatch", {
    suppressPackageStartupMessages(suppressWarnings(library(bapred)))
    
    if (!("phenotype" %in% colnames(metadata))) fail_step("FAbatch", "'phenotype' is required.")
    pheno_vals <- unique(metadata$phenotype)
    if (length(pheno_vals) != 2) fail_step("FAbatch", "'phenotype' must be binary.")
    
    # Build log-TSS features (continuous); replace non-finite per row
    X <- as.matrix(log2(normalize_tss(taxa_mat) + 1e-6))
    X <- t(apply(X, 1, function(r){
      r[!is.finite(r)] <- NA
      r[is.na(r)] <- mean(r, na.rm = TRUE)
      r
    }))
    
    y     <- factor(metadata$phenotype, levels = sort(pheno_vals))
    batch <- factor(metadata$batchid)
    
    # Variance filter (very lenient)
    v  <- apply(X, 2, var)
    keep_var <- is.finite(v) & v > 1e-12
    if (!any(keep_var)) fail_step("FAbatch", "All features ~zero variance.")
    Xv <- X[, keep_var, drop = FALSE]
    
    # Ensure p > n_b for *every* batch
    max_nb <- max(table(batch))
    # pick K = min(p, max_nb + margin)
    margin <- 5L
    K      <- min(ncol(Xv), max_nb + margin)
    if (K <= max_nb) {
      fail_step("FAbatch", sprintf("Not enough informative features to satisfy p > max batch size (have %d, need > %d).", ncol(Xv), max_nb))
    }
    
    # Select top-K variance features (no PCA, keeps life simple)
    ord <- order(apply(Xv, 2, var), decreasing = TRUE)
    sel <- ord[seq_len(K)]
    Xk  <- Xv[, sel, drop = FALSE]
    
    # Standardize for fabatch
    Xz <- scale(Xk, center = TRUE, scale = TRUE)
    
    # Try fabatch with conservative settings
    fa_out <- tryCatch(
      fabatch(
        x = Xz, y = y, batch = batch,
        nbf = NULL,          # let it choose up to maxnbf
        minerr = 1e-6,
        probcrossbatch = FALSE,
        maxiter = 100,
        maxnbf  = 8
      ),
      error = function(e) e
    )
    
    if (inherits(fa_out, "error")) {
      # Final fallback: tiny jitter to break collinearity
      warn_step("FAbatch", paste("Retrying with tiny jitter due to:", conditionMessage(fa_out)))
      Xz_jit <- Xz + matrix(rnorm(length(Xz), 0, 1e-8), nrow(Xz))
      fa_out <- fabatch(
        x = Xz_jit, y = y, batch = batch,
        nbf = NULL, minerr = 1e-6, probcrossbatch = FALSE, maxiter = 100, maxnbf = 8
      )
    }
    
    # Adjusted subset
    Xadj_sub <- fa_out$xadj
    # Un-standardize back to log space
    m <- attr(Xz, "scaled:center"); s <- attr(Xz, "scaled:scale")
    Xadj_sub <- sweep(Xadj_sub, 2, s, `*`)
    Xadj_sub <- sweep(Xadj_sub, 2, m, `+`)
    
    # Rebuild full matrix: adjusted for selected features, original for the rest
    out <- X
    colnames(out) <- colnames(X)  # already set
    out[, colnames(Xk)] <- Xadj_sub
    
    # Optional: map to positive row-normalized abundances for downstream distances
    pos <- 2^out - 1e-6
    pos[pos < 0] <- 0
    pos <- normalize_tss(pos)
    
    post_summary(pos, "FAbatch")
    write_out(pos, file.path(output_folder, "normalized_fabatch.csv"))
  })
}

# ---------------------------
# Output audit
# ---------------------------
expected_outputs <- list(
  QN="normalized_qn.csv", BMC="normalized_bmc.csv", limma="normalized_limma.csv",
  ConQuR="normalized_conqur.csv", PLSDA="normalized_plsda.csv", ComBat="normalized_combat.csv",
  FSQN="normalized_fsqn.csv", MMUPHin="normalized_mmuphin.csv",
  RUV="normalized_ruv.csv", MetaDICT="normalized_metadict.csv", SVD="normalized_svd.csv",
  PN="normalized_pn.csv", FAbatch="normalized_fabatch.csv", ComBatSeq="normalized_combatseq.csv"
)

audit_outputs <- function(selected, out_dir) {
  say("— Output audit —")
  for (m in selected) {
    fname <- expected_outputs[[m]]
    if (is.null(fname)) { say("  ", m, ": (no declared output)"); next }
    p <- file.path(out_dir, fname)
    if (file.exists(p)) {
      sz <- file.info(p)$size
      say("  ✅ ", m, " -> ", basename(p), " (", format(sz, big.mark=","), " bytes)")
    } else {
      say("  ❌ ", m, " -> ", basename(p), " (missing)")
    }
  }
}

say("🎉 All requested methods completed.")
audit_outputs(method_list, output_folder)
