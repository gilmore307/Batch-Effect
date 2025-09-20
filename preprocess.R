# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
output_folder <- if (length(args) >= 1) args[1] else "output/example"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

matrix_path <- if (length(args) >= 2) args[2] else file.path(output_folder, "raw.csv")
if (!file.exists(matrix_path)) {
  alt <- file.path("assets", "raw.csv")
  if (file.exists(alt)) {
    message("No uploaded file found at ", matrix_path, "; using ", alt)
    matrix_path <- alt
  } else {
    stop("Neither ", matrix_path, " nor ", alt, " exists.")
  }
}

# ---------------------------
# Console Logger & Guards
# ---------------------------
.ts <- function() format(Sys.time(), "%H:%M:%S")
say <- function(...){ cat(sprintf("[%s] ", .ts()), paste0(..., collapse=""), "\n") }
start_step <- function(name){ say("▶️ START: ", name); proc.time() }
ok_step    <- function(name, t0){ dt <- proc.time()-t0; say("✅ DONE:  ", name, sprintf(" (%.2fs)", dt["elapsed"])) }
warn_step  <- function(name, msg){ say("⚠️ WARN:  ", name, " — ", msg) }
fail_step  <- function(name, msg){ say("❌ FAIL:  ", name, " — ", msg); stop(paste0(name, ": ", msg)) }

post_summary <- function(M, name){
  M <- as.matrix(M)
  say(sprintf("%s summary: min=%.5f max=%.5f rowmean|avg=%.2e NA=%d",
              name, min(M, na.rm=TRUE), max(M, na.rm=TRUE),
              mean(abs(rowMeans(M))), sum(!is.finite(M))))
}

check_table <- function(x, name = "table", allow_negative = TRUE) {
  if (!is.data.frame(x) && !is.matrix(x)) fail_step(name, "Not a data.frame/matrix.")
  DF <- as.data.frame(x, stringsAsFactors = FALSE)
  if (nrow(DF) < 2 || ncol(DF) < 2) fail_step(name, "Too few rows/cols.")
  num_cols <- vapply(DF, is.numeric, TRUE)
  if (any(num_cols)) {
    Mnum <- as.matrix(DF[num_cols])
    bad_fin <- which(!is.finite(Mnum), arr.ind = TRUE)
    if (nrow(bad_fin) > 0) fail_step(name, "NA/NaN/Inf present.")
    if (!allow_negative && any(Mnum < 0, na.rm = TRUE)) fail_step(name, "Negative values not allowed.")
  }
  if (any(!num_cols)) {
    Mchr <- DF[!num_cols]
    na_idx <- which(is.na(Mchr), arr.ind = TRUE)
    if (nrow(na_idx) > 0) fail_step(name, "Missing values in non-numeric columns.")
  }
  invisible(TRUE)
}

# ---------------------------
# Input form detection + converters (minimal set)
# ---------------------------

.normalize_tss <- function(mat){
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  out <- sweep(mat, 1, row_sums, "/")
  out[is.na(out)] <- 0
  out
}

is_counts_matrix <- function(M, tol = 1e-6, frac = 0.97, min_row_sum = 1) {
  M <- as.matrix(M)
  if (any(!is.finite(M))) return(FALSE)
  if (any(M < 0, na.rm = TRUE)) return(FALSE)
  intish <- mean(abs(M - round(M)) <= tol, na.rm = TRUE)
  if (is.na(intish) || intish < frac) return(FALSE)
  rs <- rowSums(M)
  if (!any(rs > min_row_sum, na.rm = TRUE)) return(FALSE)
  TRUE
}

is_tss_matrix <- function(M, rel_tol = 0.05) {
  M <- as.matrix(M)
  nr <- nrow(M); nc <- ncol(M)
  if (!nr || !nc) return(FALSE)
  if (any(!is.finite(M))) return(FALSE)
  if (min(M) < 0)         return(FALSE)
  storage.mode(M) <- "double"
  rs <- .rowSums(M, nr, nc)
  mu  <- mean(rs)
  tol <- max(1e-5, rel_tol * abs(mu))
  mean(abs(rs - mu) <= tol) >= 0.97
}

is_clr_matrix <- function(M, tol_mean = 1e-5, min_frac = 0.9) {
  M <- as.matrix(M)
  if (any(!is.finite(M))) return(FALSE)
  nr <- nrow(M); nc <- ncol(M)
  if (nr < 2 || nc < 2) return(FALSE)
  rm <- rowMeans(M)
  mean_ok <- mean(abs(rm) <= tol_mean)
  pos <- rowSums(M > 0)
  neg <- rowSums(M < 0)
  signmix_ok <- mean((pos > 0 & neg > 0) | (pos == 0 & neg == 0))
  (mean_ok >= min_frac) && (signmix_ok >= min_frac)
}

is_nonneg_log_matrix <- function(M, zero_pos_thresh = 0.8,
                                 backcheck_max = 50000L, int_tol = 1e-6) {
  M <- as.matrix(M)
  nr <- nrow(M); nc <- ncol(M)
  if (nr < 2 || nc < 2) return(FALSE)
  if (any(!is.finite(M))) return(FALSE)
  if (min(M) < 0) return(FALSE)
  zero_pos_rows <- (rowSums(M == 0) > 0) & (rowSums(M > 0) > 0)
  if (mean(zero_pos_rows) < zero_pos_thresh) return(FALSE)
  v <- as.vector(M)
  if (length(v) > backcheck_max) v <- sample(v, backcheck_max)
  back <- expm1(v); back[back < 0] <- 0
  mean(abs(back - round(back)) <= int_tol) >= 0.5
}

detect_input_form <- function(M) {
  if (is_tss_matrix(M))        return("tss")
  if (is_clr_matrix(M))        return("clr")
  if (is_counts_matrix(M))     return("counts")
  if (is_nonneg_log_matrix(M)) return("log")
  if (all(as.matrix(M) >= 0, na.rm = TRUE)) return("positive")
  "log"
}

to_tss <- function(M, from = NULL) {
  M <- as.matrix(M)
  if (is.null(from)) from <- detect_input_form(M)
  if (from %in% c("counts", "positive", "tss")) {
    M[!is.finite(M)] <- 0
    M[M < 0] <- 0
    P <- .normalize_tss(M)
    dimnames(P) <- dimnames(M)
    return(P)
  }
  if (from == "log") {
    P <- if (all(M >= 0, na.rm = TRUE)) expm1(M) else exp(M)
    P[!is.finite(P)] <- 0
    P[P < 0] <- 0
    P <- .normalize_tss(P)
    dimnames(P) <- dimnames(M)
    return(P)
  }
  if (from == "clr") {
    M[!is.finite(M)] <- -Inf
    row_max <- apply(M, 1, function(x) {
      xm <- max(x, na.rm = TRUE)
      if (!is.finite(xm)) 0 else xm
    })
    M <- sweep(M, 1, row_max, "-")
    P <- exp(M)
    P[!is.finite(P)] <- 0
    P <- .normalize_tss(P)
    dimnames(P) <- dimnames(M)
    return(P)
  }
  stop("to_tss: unknown 'from'=", from)
}

to_clr <- function(M, from = NULL, pseudo_min = 1e-6) {
  M <- as.matrix(M)
  if (is.null(from)) from <- detect_input_form(M)
  if (from == "clr") return(M)
  C <- to_tss(M, from = from)
  C[!is.finite(C)] <- 0
  C[C < 0] <- 0
  if (any(C == 0)) {
    C <- t(apply(C, 1, function(r) {
      if (all(r == 0)) return(rep(pseudo_min, length(r)))
      nz <- r[r > 0]
      eps <- max(pseudo_min, 0.65 * min(nz))
      r[r == 0] <- eps
      r
    }))
  }
  L <- log(C)
  out <- sweep(L, 1, rowMeans(L), "-")
  dimnames(out) <- dimnames(M)
  out
}

# ---------------------------
# Writer: emits BOTH TSS and CLR with suffixes into output_folder
# ---------------------------
write_tss_clr <- function(method, native, native_type, filename) {
  base <- sub("\\.csv$", "", filename, ignore.case = TRUE)
  # --- TSS ---
  tss <- to_tss(native, from = native_type)
  post_summary(tss, paste0("🔄 ", method, " (TSS)"))
  nm_tss <- basename(file.path(output_folder, paste0(base, "_tss.csv")))
  t0 <- start_step(paste0("Write ", nm_tss))
  write.csv(tss, file.path(output_folder, paste0(base, "_tss.csv")), row.names = FALSE)
  ok_step(paste0("Write ", nm_tss), t0)
  
  # --- CLR ---
  clr <- to_clr(native, from = native_type)
  post_summary(clr, paste0("🔄 ", method, " (CLR)"))
  nm_clr <- basename(file.path(output_folder, paste0(base, "_clr.csv")))
  t0 <- start_step(paste0("Write ", nm_clr))
  write.csv(clr, file.path(output_folder, paste0(base, "_clr.csv")), row.names = FALSE)
  ok_step(paste0("Write ", nm_clr), t0)
}

# ---------------------------
# Run
# ---------------------------
say("▶ Converting to TSS & CLR")
say("Output folder: ", normalizePath(output_folder, winslash = "/"))
say("Input matrix : ", normalizePath(matrix_path,  winslash = "/"))

uploaded_mat <- read.csv(matrix_path, check.names = FALSE)
check_table(uploaded_mat, "uploaded_mat", allow_negative = TRUE)
input_form <- detect_input_form(uploaded_mat)
say("ℹ️ Detected input form: ", input_form)

base_M <- as.matrix(uploaded_mat)
write_tss_clr("INPUT", base_M, input_form, basename(matrix_path))

say("🎉 Done. Wrote _tss and _clr files into: ", normalizePath(output_folder, winslash = "/"))
