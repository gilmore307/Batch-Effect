# ===========================
# normalize_all_methods.R  (Form-aware inputs, CLR outputs)
# ===========================

# ---------------------------
# Handle Arguments
# ---------------------------
# args <- commandArgs(trailingOnly = TRUE)
args <- c("FSQN,QN,BMC,limma,ConQuR,PLSDA,ComBat,MMUPHin,RUV,MetaDICT,SVD,PN,FAbatch,ComBatSeq", "output/example")

if (length(args) < 2) stop("Usage: Rscript normalize_all_methods.R <method_list_csv> <output_folder>")

method_list   <- unlist(strsplit(args[1], ","))
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
  if (nrow(DF) < 2 || ncol(DF) < 2) fail_step(name, "Too few rows/cols.")
  num_cols <- vapply(DF, is.numeric, TRUE)
  if (any(num_cols)) {
    Mnum <- as.matrix(DF[num_cols])
    bad_fin <- which(!is.finite(Mnum), arr.ind = TRUE)
    if (nrow(bad_fin) > 0) fail_step(name, "NA/NaN/Inf present.")
  }
  if (any(!num_cols)) {
    Mchr <- DF[!num_cols]
    na_idx <- which(is.na(Mchr), arr.ind = TRUE)
    if (nrow(na_idx) > 0) fail_step(name, "Missing values in non-numeric columns.")
  }
  invisible(TRUE)
}

write_out <- function(obj, path){
  nm <- basename(path)
  t0 <- start_step(paste0("Write ", nm))
  write.csv(obj, path, row.names = FALSE)
  ok_step(paste0("Write ", nm), t0)
}

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
  say(sprintf("%s summary: min=%.5f max=%.5f rowmean|avg=%.2e NA=%d",
              name, min(M, na.rm=TRUE), max(M, na.rm=TRUE),
              mean(abs(rowMeans(M))), sum(!is.finite(M))))
}

# ---------------------------
# Input form detection + converters
# ---------------------------
.approx_eq <- function(x, y, tol=1e-6) all(abs(x - y) <= tol, na.rm = TRUE)

is_counts_matrix <- function(M, int_tol = 1e-6) {
  M <- as.matrix(M)
  if (!all(is.finite(M))) return(FALSE)
  if (any(M < 0, na.rm = TRUE)) return(FALSE)
  v <- M[is.finite(M)]
  if (!length(v)) return(FALSE)
  frac <- mean(abs(v - round(v)) <= int_tol)
  frac >= 0.97
}

is_proportions_matrix <- function(M, row_tol = 1e-3, frac_rows = 0.9) {
  M <- as.matrix(M)
  if (!all(is.finite(M))) return(FALSE)
  if (any(M < 0, na.rm = TRUE)) return(FALSE)
  rs <- rowSums(M)
  mean(abs(rs - 1) <= row_tol, na.rm = TRUE) >= frac_rows
}

is_clr_matrix <- function(M, tol_mean = 1e-5, min_frac = 0.9, require_signmix = TRUE) {
  M <- as.matrix(M)
  if (!all(is.finite(M))) return(FALSE)
  if (nrow(M) < 2 || ncol(M) < 2) return(FALSE)
  rm <- rowMeans(M)
  mean_ok <- mean(abs(rm) <= tol_mean)
  signmix_ok <- 1
  if (require_signmix) {
    pos <- rowSums(M > 0, na.rm = TRUE)
    neg <- rowSums(M < 0, na.rm = TRUE)
    signmix_ok <- mean(pos > 0 & neg > 0)
  }
  (mean_ok >= min_frac) && (signmix_ok >= min_frac)
}

is_log_like_matrix <- function(M) {
  M <- as.matrix(M)
  if (!all(is.finite(M))) return(FALSE)
  mix_frac <- mean(rowSums(M > 0) > 0 & rowSums(M < 0) > 0)
  clrish   <- is_clr_matrix(M)
  (mix_frac >= 0.8) && !clrish
}

detect_input_form <- function(M) {
  if (is_counts_matrix(M))           return("counts")
  if (is_proportions_matrix(M))      return("proportions")
  if (is_clr_matrix(M))              return("clr")
  if (is_log_like_matrix(M))         return("log")
  if (all(as.matrix(M) >= 0, na.rm = TRUE)) return("positive")
  "log"
}

# Positive helpers
.normalize_tss <- function(mat){
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  out <- sweep(mat, 1, row_sums, "/")
  out[is.na(out)] <- 0
  out
}

to_positive <- function(M, from, pseudo_min = 1e-6) {
  M <- as.matrix(M)
  if (from == "positive" || from == "proportions") return(M)
  if (from == "counts")    return(M)
  if (from == "log" || from == "clr") {
    P <- exp(M)
    P[!is.finite(P)] <- 0
    P[P < 0] <- 0
    P[P == 0] <- pseudo_min
    return(P)
  }
  stop("Unknown 'from' in to_positive: ", from)
}

to_counts <- function(M, from, scale = 1e5) {
  M <- as.matrix(M)
  if (from == "counts") return(M)
  if (from %in% c("proportions","positive")) {
    C <- round(M * scale)
    return(C)
  }
  if (from %in% c("log","clr")) {
    P <- exp(M); P[P < 0] <- 0
    C <- round(P * scale)
    return(C)
  }
  stop("Unknown 'from' in to_counts: ", from)
}

to_log <- function(M, from) {
  M <- as.matrix(M)
  if (from == "log") return(M)
  if (from == "clr") return(M)
  if (from %in% c("proportions","positive","counts")) {
    P <- M
    nz <- P[P > 0]
    if (!length(nz)) stop("to_log: no positive entries")
    pc <- max(min(nz) * 0.65, 1e-6)
    P[P == 0] <- pc
    return(log(P))
  }
  stop("Unknown 'from' in to_log: ", from)
}

to_clr <- function(M, from) {
  M <- as.matrix(M)
  if (from == "clr") return(sweep(M, 1, rowMeans(M), "-"))
  if (from == "log") return(sweep(M, 1, rowMeans(M), "-"))
  if (from %in% c("proportions","positive","counts")) {
    P <- M
    nz <- P[P > 0]
    if (!length(nz)) stop("to_clr: no positive entries")
    pc <- max(min(nz) * 0.65, 1e-6)
    P[P == 0] <- pc
    L <- log(P)
    return(sweep(L, 1, rowMeans(L), "-"))
  }
  stop("Unknown 'from' in to_clr: ", from)
}

# Method → expected input form (native)
expected_input <- list(
  QN="positive",
  BMC="log",
  limma="log",
  ConQuR="positive",
  PLSDA="clr",
  ComBat="log",
  FSQN="positive",
  MMUPHin="positive",      # accepts counts or proportions; we give positive abundances
  RUV="counts",
  MetaDICT="positive",
  SVD="log",
  PN="positive",
  FAbatch="log",
  ComBatSeq="counts"
)

# Convert native outputs to CLR for writing
# native_type ∈ {"clr","log","positive","counts"}
finalize_to_clr <- function(native, native_type, method_name="") {
  if (native_type == "clr") {
    out <- sweep(as.matrix(native), 1, rowMeans(native), "-")
  } else if (native_type == "log") {
    L <- as.matrix(native); L[!is.finite(L)] <- NA
    if (any(is.na(L))) {
      L <- t(apply(L, 1, function(r){ r[is.na(r)] <- mean(r, na.rm = TRUE); r }))
    }
    out <- sweep(L, 1, rowMeans(L), "-")
  } else if (native_type %in% c("positive","counts")) {
    P <- as.matrix(native)
    P[!is.finite(P)] <- NA
    P[P < 0] <- 0
    if (all(P == 0, na.rm = TRUE)) fail_step(method_name, "All zeros before CLR conversion.")
    nz <- P[P > 0]
    if (!length(nz)) fail_step(method_name, "No positive entries before CLR conversion.")
    pc <- max(min(nz) * 0.65, 1e-6)
    P[P == 0] <- pc
    L <- log(P)
    out <- sweep(L, 1, rowMeans(L), "-")
  } else {
    fail_step(method_name, paste0("Unknown native_type: ", native_type))
  }
  out
}

write_clr <- function(method, native, native_type, filename) {
  clr <- finalize_to_clr(native, native_type, method)
  post_summary(clr, paste0(method, " (CLR)"))
  write_out(clr, file.path(output_folder, filename))
}

# Unified accessor: get input matrix for a method given a base matrix and its form
get_input_for <- function(method, base_M, base_form) {
  target <- expected_input[[method]]
  if (is.null(target)) stop("Unknown method in expected_input: ", method)
  if (target == base_form) return(base_M)
  if (target == "positive")   return(to_positive(base_M, base_form))
  if (target == "log")        return(to_log(base_M, base_form))
  if (target == "clr")        return(to_clr(base_M, base_form))
  if (target == "counts")     return(to_counts(base_M, base_form))
  stop("Unhandled target: ", target)
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
# Prepare INPUT — detect form and keep as base
# ---------------------------
check_table(uploaded_mat, "uploaded_mat", allow_negative = TRUE)
input_form <- detect_input_form(uploaded_mat)
say("ℹ️ Detected input form: ", input_form)

# Attach rownames (expects metadata$sample_id)
if (!("sample_id" %in% colnames(metadata)))
  fail_step("Alignment", "'sample_id' column not found in metadata.")
if (nrow(uploaded_mat) != nrow(metadata))
  fail_step("Alignment", "Row count mismatch between matrix and metadata.")
rownames(uploaded_mat) <- metadata$sample_id
rownames(metadata)     <- metadata$sample_id

# Keep this as the base matrix + base form
base_M     <- as.matrix(uploaded_mat)
base_form  <- input_form

# Factor batch, covariates
metadata$batchid <- factor(metadata$batchid, levels = unique(metadata$batchid))
batchid <- metadata$batchid

covar <- metadata[, !(colnames(metadata) %in% c("sample_id","batchid","phenotype")), drop = FALSE]
covar <- as.data.frame(lapply(covar, function(col) {
  if (is.numeric(col))      col[is.na(col)] <- mean(col, na.rm = TRUE)
  else if (is.factor(col))  if (anyNA(col)) col[is.na(col)] <- levels(col)[1]
  else                      col[is.na(col)] <- as.character(stats::na.omit(col))[1]
  col
}))

# Reference batch for methods needing it (QN/FSQN)
reference_batch <- 0
ref_idx <- which(batchid == reference_batch)

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
})

# ---------------------------
# Methods (each uses get_input_for to receive expected input form; outputs saved as CLR)
# ---------------------------

# QN — expects positive
if ("QN" %in% method_list) {
  run_method("QN", {
    require(preprocessCore)
    X_pos   <- get_input_for("QN", base_M, base_form)
    ref_pos <- X_pos[ref_idx, , drop=FALSE]
    if (nrow(ref_pos) < 2) warn_step("QN", "Reference batch has <2 samples; results may be unstable.")
    target  <- normalize.quantiles.determine.target(ref_pos)
    qn_pos  <- normalize.quantiles.use.target(X_pos, target = target)
    dimnames(qn_pos) <- dimnames(X_pos)
    write_clr("QN", qn_pos, "positive", "normalized_qn.csv")
  })
}

# BMC — expects log(proportions)
if ("BMC" %in% method_list) {
  run_method("BMC", {
    require(pamr)
    # Use PN-positive to ensure proportions then log
    X_pos  <- get_input_for("PN", base_M, base_form)
    X_tss  <- .normalize_tss(X_pos)
    zmin   <- min(X_tss[X_tss > 0]); if (!is.finite(zmin)) fail_step("BMC", "No non-zero entries.")
    X_tss[X_tss == 0] <- zmin * 0.65
    X_log  <- log(X_tss)
    pam_in <- list(x = as.matrix(t(X_log)), batchlabels = factor(batchid))
    adj_log <- t(pamr.batchadjust(pam_in)$x)
    write_clr("BMC", adj_log, "log", "normalized_bmc.csv")
  })
}

# limma — expects log
if ("limma" %in% method_list) {
  run_method("limma", {
    require(limma)
    X_log <- get_input_for("limma", base_M, base_form)
    adj_t <- removeBatchEffect(
      t(X_log),
      batch = factor(batchid),
      covariates = if (ncol(covar) > 0) as.matrix(covar) else NULL
    )
    adj <- t(adj_t)
    write_clr("limma", adj, "log", "normalized_limma.csv")
  })
}

# ConQuR — expects positive
if ("ConQuR" %in% method_list) {
  run_method("ConQuR", {
    suppressPackageStartupMessages({ library(ConQuR); library(doParallel) })
    X_pos <- get_input_for("ConQuR", base_M, base_form)
    covariates <- metadata[, colnames(covar), drop = FALSE]; rownames(covariates) <- NULL
    num_cores <- max(1, parallel::detectCores(TRUE) - 1)
    res_pos <- suppressWarnings(
      ConQuR(
        tax_tab = X_pos,
        batchid = as.factor(metadata$batchid),
        covariates = covariates,
        batch_ref = as.character(reference_batch),
        logistic_lasso = FALSE, quantile_type = "standard", simple_match = FALSE,
        lambda_quantile = "2p/n", interplt = FALSE, delta = 0.4999,
        taus = seq(0.05, 0.95, by = 0.05), num_core = num_cores
      )
    )
    write_clr("ConQuR", res_pos, "positive", "normalized_conqur.csv")
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
      Y.bat = as.factor(metadata$batchid),
      ncomp.trt = 1, ncomp.bat = 5
    )
    write_clr("PLSDAbatch", res$X.nobatch, "clr", "normalized_plsda.csv")
  })
}

# ComBat — expects log
if ("ComBat" %in% method_list) {
  run_method("ComBat", {
    require(sva)
    X_log <- get_input_for("ComBat", base_M, base_form)
    adj_t <- ComBat(
      dat = t(X_log),
      batch = batchid,
      mod = if (ncol(covar) > 0) model.matrix(~ ., data = covar) else NULL,
      par.prior = FALSE, prior.plots = FALSE
    )
    adj <- t(adj_t)
    write_clr("ComBat", adj, "log", "normalized_combat.csv")
  })
}

# FSQN — expects positive
if ("FSQN" %in% method_list) {
  run_method("FSQN", {
    require(FSQN)
    X_pos <- get_input_for("FSQN", base_M, base_form)
    ref_pos <- X_pos[ref_idx, , drop=FALSE]
    out_pos <- quantileNormalizeByFeature(X_pos, ref_pos)
    write_clr("FSQN", out_pos, "positive", "normalized_fsqn.csv")
  })
}

# MMUPHin — expects positive (or counts)
if ("MMUPHin" %in% method_list) {
  run_method("MMUPHin", {
    require(MMUPHin)
    X_pos <- get_input_for("MMUPHin", base_M, base_form)
    fit <- adjust_batch(
      feature_abd = t(round(X_pos)),        # integer-ish positive abundances
      batch       = "batchid",
      covariates  = colnames(covar),
      data        = transform(metadata, batchid=factor(batchid)),
      control     = list(verbose = FALSE, diagnostic_plot = NULL)
    )
    out_pos <- t(fit$feature_abd_adj)
    write_clr("MMUPHin", out_pos, "positive", "normalized_mmuphin.csv")
  })
}

# RUV (fastRUV-III-NB) — expects counts
if ("RUV" %in% method_list) {
  run_method("fastRUV-III-NB", {
    suppressPackageStartupMessages(library(ruvIIInb))
    Y <- t(get_input_for("RUV", base_M, base_form))  # genes x samples counts
    if (is.null(colnames(Y))) colnames(Y) <- rownames(metadata)
    samp_ids <- colnames(Y)
    if (is.null(rownames(metadata)) || !all(samp_ids %in% rownames(metadata)))
      fail_step("RUV", "metadata rownames must match colnames.")
    keep <- rowSums(Y) > 0
    if (!any(keep)) fail_step("RUV", "All genes have zero counts.")
    Y <- Y[keep, samp_ids, drop=FALSE]
    ctl_names <- rownames(Y)
    batch_factor <- factor(metadata[samp_ids, "batchid"])
    M <- model.matrix(~ 0 + batch_factor); rownames(M) <- samp_ids
    fit <- fastruvIII.nb(
      Y = as.matrix(Y), M = M, ctl = ctl_names,
      batch = as.numeric(batch_factor), k = 2,
      pCells.touse = 0.05, use.pseudosample = FALSE,
      ncores = max(1L, min(4L, parallel::detectCores(TRUE) - 1L))
    )
    out_log <- t(assay(fit, "logPAC"))
    write_clr("RUV-III-NB", out_log, "log", "normalized_ruv.csv")
  })
}

# MetaDICT — expects positive
if ("MetaDICT" %in% method_list) {
  run_method("MetaDICT", {
    suppressPackageStartupMessages({ library(MetaDICT); library(vegan) })
    O <- t(get_input_for("MetaDICT", base_M, base_form))
    meta <- transform(metadata, batch = batchid)
    O <- O[rowSums(O) > 0, , drop=FALSE]
    if (nrow(O) < 2) fail_step("MetaDICT", "Too few non-zero samples.")
    D <- as.matrix(vegdist(O, method = "bray"))
    res <- MetaDICT(O, meta, distance_matrix = D)
    out_pos <- t(res$count)
    write_clr("MetaDICT", out_pos, "positive", "normalized_metadict.csv")
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
    
    write_clr("SVD", full, "log", "normalized_svd.csv")
  })
}

# PN — expects positive
if ("PN" %in% method_list) {
  run_method("PN", {
    if (!("phenotype" %in% colnames(metadata))) fail_step("PN", "'phenotype' is required.")
    pheno_vals <- unique(metadata$phenotype)
    if (length(pheno_vals) != 2) fail_step("PN", "'phenotype' must be binary.")
    trt <- as.numeric(factor(metadata$phenotype, levels = sort(pheno_vals))) - 1
    X_pos <- get_input_for("PN", base_M, base_form)
    X_tss <- .normalize_tss(X_pos)
    if (all(X_tss == 0)) fail_step("PN", "All zero after TSS.")
    pn_pos <- percentile_norm(data = X_tss, batch = metadata$batchid, trt = trt, ctrl.grp = 0)
    write_clr("PN", pn_pos, "positive", "normalized_pn.csv")
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
    batch <- factor(metadata$batchid)
    
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
    
    write_clr("FAbatch", Xadj, "log", "normalized_fabatch.csv")
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
    adj <- ComBat_seq(counts = t(counts), batch = metadata$batchid, group = metadata$phenotype)
    out_counts <- t(adj)
    write_clr("ComBat-Seq", out_counts, "counts", "normalized_combatseq.csv")
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
  say("— Output audit (all outputs are CLR) —")
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

say("🎉 All requested methods completed. (Outputs are CLR)")
audit_outputs(method_list, output_folder)
