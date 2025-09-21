# ---------------------------
# Handle Arguments
# ---------------------------
# All Methods: FSQN, QN, BMC, limma, ConQuR, PLSDA, ComBat, MMUPHin, RUV, MetaDICT, SVD, PN, FAbatch, ComBatSeq, DEBIAS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- c("FSQN, QN, BMC, limma, ConQuR, PLSDA, ComBat, MMUPHin, RUV, MetaDICT, SVD, PN, FAbatch, ComBatSeq, DEBIAS", "output/example")
}

# Split and trim whitespace so names match exactly
method_list   <- trimws(unlist(strsplit(args[1], ",")))
# Drop empty tokens
method_list   <- method_list[nzchar(method_list)]
output_folder <- args[2]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# ---------------------------
# Load Libraries
# ---------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(SummarizedExperiment)
  library(S4Vectors)
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
# Source preprocess.r (helpers: say(), write_tss_clr(), detect_input_form(),
# to_tss(), to_log(), to_clr(), to_counts(), check_table(), post_summary(), ...)
# Ensure preprocess.R uses our intended output/matrix paths (avoid commandArgs hijack)
# ---------------------------
assign("PREPROC_OUTPUT_DIR", output_folder, envir = .GlobalEnv)
assign("PREPROC_MATRIX_PATH", file.path(output_folder, "raw.csv"), envir = .GlobalEnv)
source("preprocess.R")

# If preprocess.r did not define the logger for some reason, provide a minimal one
if (!exists("say", mode = "function")) {
  .ts <- function() format(Sys.time(), "%H:%M:%S")
  say <- function(...){ cat(sprintf("[%s] ", .ts()), paste0(..., collapse=""), "\n") }
  start_step <- function(name){ say("▶️ START: ", name); proc.time() }
  ok_step    <- function(name, t0){ dt <- proc.time()-t0; say("✅ DONE:  ", name, sprintf(" (%.2fs)", dt["elapsed"])) }
  warn_step  <- function(name, msg){ say("⚠️ WARN:  ", name, " — ", msg) }
  fail_step  <- function(name, msg){ say("❌ FAIL:  ", name, " — ", msg); stop(paste0(name, ": ", msg)) }
}

 

# ---------------------------
# Helper to run a step with logging
# ---------------------------
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
    warning = function(w){ warn_step(name, conditionMessage(w)) }
  )
  invisible(out)
}

# ---------------------------
# Desired inputs per method (uses converters from preprocess.r)
# ---------------------------
expected_input <- list(
  QN        = "tss",
  BMC       = "log",
  limma     = "log",
  ConQuR    = "counts",
  PLSDA     = "clr",
  ComBat    = "log",
  FSQN      = "tss",
  MMUPHin   = "tss",
  RUV       = "counts",
  MetaDICT  = "tss",
  SVD       = "log",
  PN        = "tss",
  FAbatch   = "log",
  ComBatSeq = "counts",
  DEBIAS    = "counts"
)

# get input matrix for a method given a base matrix and its form
get_input_for <- function(method, base_M, base_form) {
  target <- expected_input[[method]]
  if (is.null(target)) stop("Unknown method in expected_input: ", method)
  if (target == base_form) return(base_M)
  say(sprintf("🔄 Convert for %s: %s → %s", method, base_form, target))
  out <- switch(
    target,
    "tss"    = to_tss(base_M, base_form),
    "log"    = to_log(base_M, base_form),
    "clr"    = to_clr(base_M, base_form),
    "counts" = to_counts(base_M, base_form),
    stop("Unhandled target: ", target)
  )
  say(sprintf("   ↳ done (%s → %s) | n=%d×%d", base_form, target, nrow(out), ncol(out)))
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
# Prepare INPUT
# ---------------------------
check_table(uploaded_mat, "uploaded_mat", allow_negative = TRUE)
input_form <- detect_input_form(uploaded_mat)
say("ℹ️ Detected input form: ", input_form)

if (!("sample_id" %in% colnames(metadata)))
  fail_step("Alignment", "'sample_id' column not found in metadata.")

# Make duplicate sample IDs unique (append _1, _2, ...)
metadata <- make_unique_sample_ids(metadata, id_col = "sample_id", sep = "_")
if (nrow(uploaded_mat) != nrow(metadata))
  fail_step("Alignment", "Row count mismatch between matrix and metadata.")
rownames(uploaded_mat) <- metadata$sample_id
rownames(metadata)     <- metadata$sample_id

base_M     <- as.matrix(uploaded_mat)
base_form  <- input_form

# Factor batch, covariates
metadata$batch_id <- factor(metadata$batch_id, levels = unique(metadata$batch_id))
batch_id <- metadata$batch_id

covar <- metadata[, !(colnames(metadata) %in% c("sample_id","batch_id","phenotype")), drop = FALSE]
covar <- as.data.frame(lapply(covar, function(col) {
  if (is.numeric(col))      col[is.na(col)] <- mean(col, na.rm = TRUE)
  else if (is.factor(col))  { if (anyNA(col)) col[is.na(col)] <- levels(col)[1] }
  else                      { if (anyNA(col)) { kept <- suppressWarnings(as.character(stats::na.omit(col))); if (length(kept)) col[is.na(col)] <- kept[1] } }
  col
}))

# Reference batch for methods needing it
reference_batch <- levels(batch_id)[1]
ref_idx <- which(batch_id == reference_batch)

# ---------------------------
# Global checks
# ---------------------------
run_method("Input checks", {
  check_table(base_M, "base_M (detected input)", allow_negative = TRUE)
  check_table(metadata, "metadata", allow_negative = TRUE)
  if (!identical(rownames(base_M), rownames(metadata))) {
    fail_step("Alignment", "rownames(base_M) != rownames(metadata).")
  }
  say("Method list: ", paste(method_list, collapse=", "))
  say("Samples: ", nrow(base_M), " | Features: ", ncol(base_M))
  say("Reference batch level: ", as.character(reference_batch),
      " | #ref samples: ", length(ref_idx))
  if (length(ref_idx) < 1) warn_step("Reference", "No samples in reference batch level.")
})

# ---------------------------
# Methods (uses get_input_for + write_tss_clr from preprocess.r)
# ---------------------------

# QN — expects TSS
if ("QN" %in% method_list) {
  run_method("QN", {
    require(preprocessCore)
    X_tss  <- get_input_for("QN", base_M, base_form)         # rows = samples
    ref_tss <- X_tss[ref_idx, , drop=FALSE]
    if (nrow(ref_tss) < 2) warn_step("QN", "Reference batch has <2 samples; results may be unstable.")
    target  <- normalize.quantiles.determine.target(ref_tss)
    qn_tss  <- normalize.quantiles.use.target(X_tss, target = target)
    dimnames(qn_tss) <- dimnames(X_tss)
    write_tss_clr("QN", qn_tss, "positive", "normalized_qn.csv")
  })
}

# BMC — expects log(TSS)
if ("BMC" %in% method_list) {
  run_method("BMC", {
    require(pamr)
    X_log  <- get_input_for("BMC", base_M, base_form)
    pam_in <- list(x = as.matrix(t(X_log)), batchlabels = factor(batch_id))
    adj_log <- t(pamr.batchadjust(pam_in)$x)
    write_tss_clr("BMC", adj_log, "log", "normalized_bmc.csv")
  })
}

# limma — expects log
if ("limma" %in% method_list) {
  run_method("limma", {
    require(limma)
    X_log <- get_input_for("limma", base_M, base_form)
    adj_t <- removeBatchEffect(
      t(X_log),
      batch = factor(batch_id),
      covariates = if (ncol(covar) > 0) as.matrix(covar) else NULL
    )
    adj <- t(adj_t)
    write_tss_clr("limma", adj, "log", "normalized_limma.csv")
  })
}

# ConQuR — expects counts
if ("ConQuR" %in% method_list) {
  run_method("ConQuR", {
    suppressPackageStartupMessages({ library(ConQuR); library(doParallel) })
    X_cnt <- get_input_for("ConQuR", base_M, base_form)
    covariates <- metadata[, colnames(covar), drop = FALSE]; rownames(covariates) <- NULL
    num_cores <- max(1, parallel::detectCores(TRUE) - 1)
    res_pos <- suppressWarnings(
      ConQuR(
        tax_tab = X_cnt,
        batchid = as.factor(metadata$batch_id),
        covariates = covariates,
        batch_ref = as.character(reference_batch),
        logistic_lasso = FALSE, quantile_type = "standard", simple_match = FALSE,
        lambda_quantile = "2p/n", interplt = FALSE, delta = 0.4999,
        taus = seq(0.05, 0.95, by = 0.05), num_core = num_cores
      )
    )
    write_tss_clr("ConQuR", res_pos, "positive", "normalized_conqur.csv")
  })
}

# PLSDA — expects CLR
if ("PLSDA" %in% method_list) {
  run_method("PLSDAbatch", {
    require(PLSDAbatch)
    if (!("phenotype" %in% colnames(metadata))) fail_step("PLSDAbatch", "'phenotype' not found.")
    if (length(unique(metadata$phenotype)) != 2) fail_step("PLSDAbatch", "'phenotype' must be binary.")
    X_clr <- get_input_for("PLSDA", base_M, base_form)
    res <- PLSDA_batch(
      X = X_clr,
      Y.trt = as.factor(metadata$phenotype),
      Y.bat = as.factor(metadata$batch_id),
      ncomp.trt = 1, ncomp.bat = 5
    )
    write_tss_clr("PLSDAbatch", res$X.nobatch, "clr", "normalized_plsda.csv")
  })
}

# ComBat — expects log
if ("ComBat" %in% method_list) {
  run_method("ComBat", {
    require(sva)
    X_log <- get_input_for("ComBat", base_M, base_form)
    adj_t <- ComBat(
      dat = t(X_log),
      batch = batch_id,
      mod = if (ncol(covar) > 0) model.matrix(~ ., data = covar) else NULL,
      par.prior = FALSE, prior.plots = FALSE
    )
    adj <- t(adj_t)
    write_tss_clr("ComBat", adj, "log", "normalized_combat.csv")
  })
}

# FSQN — expects TSS
if ("FSQN" %in% method_list) {
  run_method("FSQN", {
    require(FSQN)
    X_tss <- get_input_for("FSQN", base_M, base_form)
    ref_tss <- X_tss[ref_idx, , drop=FALSE]
    out_tss <- quantileNormalizeByFeature(X_tss, ref_tss)
    write_tss_clr("FSQN", out_tss, "positive", "normalized_fsqn.csv")
  })
}

# MMUPHin — expects TSS
if ("MMUPHin" %in% method_list) {
  run_method("MMUPHin", {
    require(MMUPHin)
    X_tss <- get_input_for("MMUPHin", base_M, base_form)
    feat_counts <- t(round(X_tss * 1e6))  # features x samples
    fit <- adjust_batch(
      feature_abd = feat_counts,
      batch       = "batch_id",
      covariates  = colnames(covar),
      data        = transform(metadata, batch_id=factor(batch_id)),
      control     = list(verbose = FALSE, diagnostic_plot = NULL)
    )
    out_pos <- t(fit$feature_abd_adj)
    write_tss_clr("MMUPHin", out_pos, "positive", "normalized_mmuphin.csv")
  })
}

# RUV (fastRUV-III-NB) — expects counts; use the *counts* assay for outputs
if ("RUV" %in% method_list) {
  run_method("fastRUV-III-NB", {
    suppressPackageStartupMessages(library(ruvIIInb))
    # ruviinb expects genes x samples; we build Y accordingly
    Y <- t(get_input_for("RUV", base_M, base_form))  # genes x samples
    if (is.null(colnames(Y))) colnames(Y) <- rownames(metadata)
    samp_ids <- colnames(Y)
    if (is.null(rownames(metadata)) || !all(samp_ids %in% rownames(metadata)))
      fail_step("RUV", "metadata rownames must match colnames.")
    keep <- rowSums(Y) > 0
    if (!any(keep)) fail_step("RUV", "All genes have zero counts.")
    Y <- Y[keep, samp_ids, drop=FALSE]
    ctl_names <- rownames(Y)  # using all genes as controls (fast variant)
    batch_factor <- factor(metadata[samp_ids, "batch_id"])
    M <- model.matrix(~ 0 + batch_factor); rownames(M) <- samp_ids
    
    fit <- fastruvIII.nb(
      Y = as.matrix(Y), M = M, ctl = ctl_names,
      batch = as.numeric(batch_factor), k = 2,
      pCells.touse = 0.05, use.pseudosample = FALSE,
      ncores = max(1L, min(4L, parallel::detectCores(TRUE) - 1L))
    )
    
    # Use the *counts* assay as native output for downstream TSS/CLR files
    out_counts <- t(assay(fit, "counts"))   # samples x genes/features
    write_tss_clr("RUV-III-NB", out_counts, "counts", "normalized_ruv.csv")
  })
}

# MetaDICT — expects TSS
if ("MetaDICT" %in% method_list) {
  run_method("MetaDICT", {
    suppressPackageStartupMessages({ library(MetaDICT); library(vegan) })
    O <- t(get_input_for("MetaDICT", base_M, base_form))  # samples x features for vegdist
    meta <- transform(metadata, batch = batch_id)
    O <- O[rowSums(O) > 0, , drop=FALSE]
    if (nrow(O) < 2) fail_step("MetaDICT", "Too few non-zero samples.")
    D <- as.matrix(vegdist(O, method = "bray"))
    res <- MetaDICT(O, meta, distance_matrix = D)
    out_pos <- t(res$count)
    write_tss_clr("MetaDICT", out_pos, "positive", "normalized_metadict.csv")
  })
}

# SVD — expects log
if ("SVD" %in% method_list) {
  run_method("SVD", {
    cat("Running SVD-based batch correction in log/CLR space...\n")
    X_log <- get_input_for("SVD", base_M, base_form)
    # Replace non-finite with row means
    X_log <- t(apply(X_log, 1, function(r){ r[!is.finite(r)] <- NA; r[is.na(r)] <- mean(r, na.rm = TRUE); r }))
    # Identify variable features
    feature_sd <- apply(X_log, 2, sd)
    zero_var_cols <- which(feature_sd == 0)
    variable_cols <- setdiff(seq_len(ncol(X_log)), zero_var_cols)
    if (!length(variable_cols)) fail_step("SVD", "All features zero variance.")
    Xv <- X_log[, variable_cols, drop = FALSE]
    mu <- colMeans(Xv); sdv <- apply(Xv, 2, sd)
    Z  <- scale(Xv, center = TRUE, scale = TRUE)
    # Remove first principal component (surrogate/batch component)
    s  <- svd(crossprod(Z))
    a1 <- s$u[,1]
    t1 <- Z %*% a1 / sqrt(drop(crossprod(a1)))
    c1 <- crossprod(Z, t1) / drop(crossprod(t1))
    Zdef <- Z - t1 %*% t(c1)
    # Rescale back to log space
    Xrest <- sweep(Zdef, 2, sdv, `*`)
    Xrest <- sweep(Xrest, 2, mu, `+`)
    full <- X_log
    full[, variable_cols] <- Xrest
    if (length(zero_var_cols) > 0) full[, zero_var_cols] <- X_log[, zero_var_cols]
    write_tss_clr("SVD", full, "log", "normalized_svd.csv")
  })
}

# PN — expects TSS
if ("PN" %in% method_list) {
  run_method("PN", {
    if (!("phenotype" %in% colnames(metadata))) fail_step("PN", "'phenotype' is required.")
    pheno_vals <- unique(metadata$phenotype)
    if (length(pheno_vals) != 2) fail_step("PN", "'phenotype' must be binary.")
    trt <- as.numeric(factor(metadata$phenotype, levels = sort(pheno_vals))) - 1
    X_tss <- get_input_for("PN", base_M, base_form)
    if (all(X_tss == 0)) fail_step("PN", "All zero after TSS.")
    pn_pos <- percentile_norm(data = X_tss, batch = metadata$batch_id, trt = trt, ctrl.grp = 0)
    write_tss_clr("PN", pn_pos, "positive", "normalized_pn.csv")
  })
}

# FAbatch — expects log
if ("FAbatch" %in% method_list) {
  run_method("FAbatch", {
    suppressPackageStartupMessages(library(bapred))
    if (!("phenotype" %in% colnames(metadata))) fail_step("FAbatch", "'phenotype' is required.")
    pheno_vals <- unique(metadata$phenotype)
    if (length(pheno_vals) != 2) fail_step("FAbatch", "'phenotype' must be binary.")
    X_log <- get_input_for("FAbatch", base_M, base_form)
    # replace non-finite per row
    X_log <- t(apply(X_log, 1, function(r){
      r[!is.finite(r)] <- NA
      r[is.na(r)] <- mean(r, na.rm = TRUE)
      r
    }))
    y     <- factor(metadata$phenotype, levels = sort(pheno_vals))
    batch <- factor(metadata$batch_id)
    v  <- apply(X_log, 2, var)
    keep_var <- is.finite(v) & v > 1e-12
    if (!any(keep_var)) fail_step("FAbatch", "All features ~zero variance.")
    Xv <- X_log[, keep_var, drop = FALSE]
    max_nb <- max(table(batch))
    K      <- min(ncol(Xv), max_nb + 5L)
    if (K <= max_nb) {
      fail_step("FAbatch", sprintf("Need p > max batch size (have %d, need > %d).", ncol(Xv), max_nb))
    }
    ord <- order(apply(Xv, 2, var), decreasing = TRUE)
    sel <- ord[seq_len(K)]
    Xk  <- Xv[, sel, drop = FALSE]
    Xz <- scale(Xk, center = TRUE, scale = TRUE)
    fa_out <- tryCatch(
      fabatch(
        x = Xz, y = y, batch = batch,
        nbf = NULL, minerr = 1e-6, probcrossbatch = FALSE, maxiter = 100, maxnbf = 8
      ),
      error = function(e) e
    )
    if (inherits(fa_out, "error")) {
      warn_step("FAbatch", paste("Retry with tiny jitter:", conditionMessage(fa_out)))
      Xz <- Xz + matrix(rnorm(length(Xz), 0, 1e-8), nrow(Xz))
      fa_out <- fabatch(
        x = Xz, y = y, batch = batch,
        nbf = NULL, minerr = 1e-6, probcrossbatch = FALSE, maxiter = 100, maxnbf = 8
      )
    }
    Xadj <- X_log
    m <- attr(Xz, "scaled:center"); s <- attr(Xz, "scaled:scale")
    Xadj_sub <- sweep(fa_out$xadj, 2, s, `*`)
    Xadj_sub <- sweep(Xadj_sub, 2, m, `+`)
    Xadj[, colnames(Xk)] <- Xadj_sub
    write_tss_clr("FAbatch", Xadj, "log", "normalized_fabatch.csv")
  })
}

# ComBat-Seq — expects counts
if ("ComBatSeq" %in% method_list) {
  run_method("ComBat-Seq", {
    require(sva)
    if (!("phenotype" %in% colnames(metadata))) fail_step("ComBat-Seq", "'phenotype' is required.")
    counts <- get_input_for("ComBatSeq", base_M, base_form)
    libsz <- rowSums(counts); keep <- libsz > 0
    if (any(!keep)) {
      say("ComBat-Seq: removing ", sum(!keep), " samples with zero library size.")
      counts  <- counts[keep, , drop = FALSE]
      metadata <- metadata[keep, , drop = FALSE]
    }
    if (!all(rownames(counts) == rownames(metadata))) fail_step("ComBat-Seq", "Sample IDs mismatch after filtering.")
    adj <- ComBat_seq(counts = t(counts), batch = metadata$batch_id, group = metadata$phenotype)
    out_counts <- t(adj)
    write_tss_clr("ComBat-Seq", out_counts, "counts", "normalized_combatseq.csv")
  })
}

# DEBIAS — Python package via reticulate; expects counts + phenotype; writes TSS+CLR
if ("DEBIAS" %in% method_list) {
  run_method("DEBIAS", {
    suppressPackageStartupMessages(library(reticulate))
    py <- Sys.getenv("RETICULATE_PYTHON")
    if (!nzchar(py)) py <- "C:/Users/sunch/ANACON~1/python.exe"  # your Anaconda Python
    reticulate::use_python(py, required = TRUE)
    
    have_np <- py_module_available("numpy")
    have_dm <- py_module_available("debiasm")
    if (!have_np || !have_dm) {
      system2(py, c("-m","pip","install","--upgrade","pip"), stdout = TRUE, stderr = TRUE)
      if (!have_np) system2(py, c("-m","pip","install","numpy"), stdout = TRUE, stderr = TRUE)
      if (!have_dm) system2(py, c("-m","pip","install","DEBIAS-M"), stdout = TRUE, stderr = TRUE)
    }
    
    np <- import("numpy", delay_load = TRUE)
    debiasm <- import("debiasm", delay_load = TRUE)
    
    X_cnt <- tryCatch(
      get_input_for("DEBIAS", base_M, base_form),
      error = function(e) { say("DEBIAS: expected_input missing; converting to counts."); to_counts(base_M, base_form) }
    )
    X_cnt[!is.finite(X_cnt)] <- 0
    X_cnt[X_cnt < 0] <- 0
    if (!("phenotype" %in% colnames(metadata))) {
      fail_step("DEBIAS", "'phenotype' column required in metadata.")
    }
    y_all <- metadata$phenotype
    
    # classifier vs regressor
    is_num <- is.numeric(y_all)
    uniq_vals <- unique(y_all[!is.na(y_all)])
    is_integerish <- is_num && all(abs(uniq_vals - round(uniq_vals)) < 1e-8)
    use_classifier <- (!is_num) || (is_integerish && length(uniq_vals) <= 10)
    
    # batch in first col (0-based), then counts
    b0 <- as.integer(factor(batch_id)) - 1L
    if (any(is.na(b0))) fail_step("DEBIAS", "Invalid batch IDs.")
    X_with_batch <- cbind(b0, round(X_cnt))
    
    # split (prefer a held-out batch); ensure >=2 train rows
    uniq_b <- sort(unique(b0))
    if (length(uniq_b) >= 2) {
      val_batch <- tail(uniq_b, 1)
      val_inds  <- (b0 == val_batch)
    } else {
      set.seed(123)
      val_inds  <- rep(FALSE, nrow(X_with_batch))
      val_inds[sample.int(nrow(X_with_batch), max(1L, floor(0.2 * nrow(X_with_batch))))] <- TRUE
    }
    if (sum(!val_inds) < 2) {
      set.seed(123)
      val_inds  <- rep(FALSE, nrow(X_with_batch))
      val_inds[sample.int(nrow(X_with_batch), max(1L, floor(0.2 * nrow(X_with_batch))))] <- TRUE
    }
    
    X_train_R <- X_with_batch[!val_inds, , drop = FALSE]
    X_val_R   <- X_with_batch[val_inds, , drop = FALSE]
    y_train_R <- y_all[!val_inds]
    
    if (use_classifier && length(unique(y_train_R[!is.na(y_train_R)])) < 2 && length(unique(y_all)) >= 2) {
      set.seed(42)
      val_inds  <- rep(FALSE, nrow(X_with_batch))
      val_inds[sample.int(nrow(X_with_batch), max(1L, floor(0.2 * nrow(X_with_batch))))] <- TRUE
      X_train_R <- X_with_batch[!val_inds, , drop = FALSE]
      X_val_R   <- X_with_batch[val_inds, , drop = FALSE]
      y_train_R <- y_all[!val_inds]
    }
    
    # numpy arrays
    X_train <- np$array(X_train_R, dtype = "int64")
    X_val   <- np$array(X_val_R,   dtype = "int64")
    X_full  <- np$array(X_with_batch, dtype = "int64")
    
    if (use_classifier) {
      yf <- as.factor(y_train_R)
      if (anyNA(yf)) { tab <- sort(table(yf), decreasing = TRUE); yf[is.na(yf)] <- names(tab)[1] }
      y_train <- np$array(as.integer(yf) - 1L, dtype = "int64")
      model <- debiasm$DebiasMClassifier(x_val = X_val)
    } else {
      yr <- as.numeric(y_train_R); if (anyNA(yr)) yr[is.na(yr)] <- mean(yr, na.rm = TRUE)
      y_train <- np$array(yr, dtype = "float32")
      model <- debiasm$DebiasMRegressor(x_val = X_val)
    }
    
    model$fit(X_train, y_train)
    X_debiased <- model$transform(X_full)
    Xd <- as.matrix(py_to_r(X_debiased))
    
    # robust shape handling
    p_full <- ncol(X_with_batch)   # 1 + features
    p_feat <- ncol(X_cnt)
    if (nrow(Xd) == p_full && ncol(Xd) == nrow(X_with_batch)) Xd <- t(Xd)
    
    if (ncol(Xd) == p_full) {
      out_counts <- Xd[, -1, drop = FALSE]          # includes batch -> drop it
    } else if (ncol(Xd) == p_feat) {
      out_counts <- Xd                               # features-only
    } else {
      fail_step("DEBIAS", sprintf("Unexpected transform shape: got %d cols (p_feat=%d, p_full=%d).",
                                  ncol(Xd), p_feat, p_full))
    }
    
    dimnames(out_counts) <- dimnames(X_cnt)
    write_tss_clr("DEBIAS", out_counts, "counts", "normalized_debias.csv")
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
  PN="normalized_pn.csv", FAbatch="normalized_fabatch.csv", ComBatSeq="normalized_combatseq.csv",
  DEBIAS = "normalized_debias.csv"
)

audit_outputs <- function(selected, out_dir) {
  say("— Output audit (TSS & CLR) —")
  for (m in selected) {
    fname <- expected_outputs[[m]]
    if (is.null(fname)) { say("  ", m, ": (no declared output)"); next }
    base <- sub("\\.csv$", "", fname, ignore.case = TRUE)
    p_tss <- file.path(out_dir, paste0(base, "_tss.csv"))
    p_clr <- file.path(out_dir, paste0(base, "_clr.csv"))
    
    if (file.exists(p_tss)) {
      sz <- file.info(p_tss)$size
      say("  ✅ ", m, " -> ", basename(p_tss), " (", format(sz, big.mark=","), " bytes)")
    } else {
      say("  ❌ ", m, " -> ", basename(p_tss), " (missing)")
    }
    
    if (file.exists(p_clr)) {
      sz <- file.info(p_clr)$size
      say("  ✅ ", m, " -> ", basename(p_clr), " (", format(sz, big.mark=","), " bytes)")
    } else {
      say("  ❌ ", m, " -> ", basename(p_clr), " (missing)")
    }
  }
}

say("🎉 All requested methods completed. (Outputs are TSS + CLR)")
audit_outputs(method_list, output_folder)
