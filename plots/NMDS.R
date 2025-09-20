# =========================
# NMDS figures (Aitchison on CLR; Bray-Curtis on TSS)
# =========================

# ==== Libraries ====
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(patchwork)   # layouts + legend collection

# Map method codes to short labels for figures
method_short_label <- function(x) {
  map <- c(
    qn = "QN", bmc = "BMC", limma = "Limma", conqur = "ConQuR",
    plsda = "PLSDA-batch", combat = "ComBat", fsqn = "FSQN", mmuphin = "MMUPHin",
    ruv = "RUV-III-NB", metadict = "MetaDICT", svd = "SVD", pn = "PN",
    fabatch = "FAbatch", combatseq = "ComBat-seq", debias = "DEBIAS-M"
  )
  sapply(x, function(v){ lv <- tolower(v); if (lv %in% names(map)) map[[lv]] else v })
}

  library(rlang)
  library(vegan)       # distances + NMDS
})

# ==== Helpers (robust ranges/padding) ====
safe_range <- function(v) { v <- v[is.finite(v)]; if (!length(v)) c(0,0) else range(v, na.rm = TRUE) }
finite_range <- function(...) { v <- unlist(list(...)); v <- v[is.finite(v)]; if (!length(v)) c(0,0) else range(v, na.rm = TRUE) }
safe_pad <- function(r, frac = 0.06) { dx <- diff(r); if (!is.finite(dx) || dx <= 0) 1e-6 else dx * frac }

# ==== Ellipse union bounds (falls back to point ranges if needed) ====
ellipse_union_bounds <- function(df_scores, group_var, level = 0.95, n = 240) {
  if (!nrow(df_scores)) return(list(x = c(0,0), y = c(0,0)))
  chi <- sqrt(qchisq(level, df = 2))
  xmins <- c(); xmaxs <- c(); ymins <- c(); ymaxs <- c()
  levs <- levels(df_scores[[group_var]]); if (is.null(levs)) levs <- unique(df_scores[[group_var]])
  for (lev in levs) {
    sub <- df_scores[df_scores[[group_var]] == lev, c("AX1","AX2"), drop = FALSE]
    if (nrow(sub) < 3 || any(!is.finite(as.matrix(sub)))) next
    S  <- tryCatch(stats::cov(sub, use = "complete.obs"), error = function(e) NULL)
    mu <- colMeans(sub, na.rm = TRUE)
    if (is.null(S) || any(!is.finite(S))) next
    eg <- eigen(S, symmetric = TRUE); if (any(!is.finite(eg$values))) next
    t  <- seq(0, 2*pi, length.out = n)
    R  <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)))
    pts <- t(chi * R %*% rbind(cos(t), sin(t))); pts <- sweep(pts, 2, mu, "+")
    xmins <- c(xmins, min(pts[,1], na.rm = TRUE)); xmaxs <- c(xmaxs, max(pts[,1], na.rm = TRUE))
    ymins <- c(ymins, min(pts[,2], na.rm = TRUE)); ymaxs <- c(ymaxs, max(pts[,2], na.rm = TRUE))
  }
  if (!length(xmins) || !length(ymins)) {
    return(list(x = safe_range(df_scores$AX1), y = safe_range(df_scores$AX2)))
  }
  list(x = c(min(xmins, na.rm = TRUE), max(xmaxs, na.rm = TRUE)),
       y = c(min(ymins, na.rm = TRUE), max(ymaxs, na.rm = TRUE)))
}

# ==== Compositional helpers ====
safe_closure <- function(X) {
  X[!is.finite(X)] <- 0; X[X < 0] <- 0
  rs <- rowSums(X); rs[!is.finite(rs) | rs == 0] <- 1
  sweep(X, 1, rs, "/")
}
clr_transform <- function(X) {
  Xt <- safe_closure(X)
  # replace zeros with a small replacement then renormalize
  for (i in seq_len(nrow(Xt))) {
    xi <- Xt[i, ]; pos <- xi > 0
    if (!any(pos)) { xi[] <- 1/length(xi); pos <- xi > 0 }
    if (any(!pos)) {
      eps <- min(xi[pos]) * 0.5
      xi[!pos] <- eps; xi <- xi / sum(xi)
    }
    Xt[i, ] <- xi
  }
  L <- log(Xt); sweep(L, 1, rowMeans(L), "-")
}

# ==== Args / IO ====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))
}
if (!("batch_id" %in% names(metadata)) && ("batch_id" %in% names(metadata))) {
  metadata$batch_id <- metadata$batch_id
}

# ---- Find normalized files ----
clr_paths <- list.files(output_folder, pattern = "^normalized_.*_clr\\.csv$", full.names = TRUE)
tss_paths <- list.files(output_folder, pattern = "^normalized_.*_tss\\.csv$", full.names = TRUE)

# Fallback: if no suffix-specific outputs, use any normalized_*.csv for both
if (!length(clr_paths) && !length(tss_paths)) {
  any_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
  clr_paths <- any_paths
  tss_paths <- any_paths
}

name_from <- function(paths, suffix) gsub(paste0("^normalized_|_", suffix, "\\.csv$"), "", basename(paths))
file_list_clr <- setNames(clr_paths, method_short_label(name_from(clr_paths, "clr")))
file_list_tss <- setNames(tss_paths, method_short_label(name_from(tss_paths, "tss")))

# Include raw_clr.csv / raw_tss.csv as "Before correction" if present
raw_clr_fp <- file.path(output_folder, "raw_clr.csv")
raw_tss_fp <- file.path(output_folder, "raw_tss.csv")
if (file.exists(raw_clr_fp)) file_list_clr <- c("Before correction" = raw_clr_fp, file_list_clr)
if (file.exists(raw_tss_fp)) file_list_tss <- c("Before correction" = raw_tss_fp, file_list_tss)

if (!length(file_list_clr) && !length(file_list_tss)) {
  stop("No normalized files found (expected raw_clr.csv/raw_tss.csv and/or normalized_*_clr.csv / normalized_*_tss.csv) in ", output_folder)
}

# ==== NMDS frames: Aitchison on CLR ====
compute_nmds_frames_aitch <- function(df, metadata, model.vars = c("batch_id","phenotype"), k = 2) {
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfm <- inner_join(df, metadata, by = "sample_id")
  
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  if (!ncol(X)) stop("No numeric features for NMDS (CLR).")
  # If negatives exist, assume CLR; else build CLR from non-negative data
  Xclr <- if (any(X < 0, na.rm = TRUE)) {
    sweep(X, 1, rowMeans(X, na.rm = TRUE), "-")
  } else {
    clr_transform(X)
  }
  Da <- dist(Xclr, method = "euclidean")
  
  set.seed(123)
  fit <- vegan::monoMDS(Da, k = k, maxit = 500, smin = 1e-12, trace = FALSE)
  sc  <- fit$points; colnames(sc) <- paste0("NMDS", seq_len(ncol(sc)))
  plot.df <- data.frame(sample_id = dfm$sample_id, as.data.frame(sc), check.names = FALSE)
  present <- intersect(model.vars, names(dfm))
  for (v in present) plot.df[[v]] <- factor(as.character(dfm[[v]]), levels = unique(as.character(metadata[[v]])))
  list(plot.df = plot.df, stress = fit$stress, used.vars = present)
}

# ==== NMDS frames: Bray-Curtis on TSS ====
compute_nmds_frames_bray <- function(df, metadata, model.vars = c("batch_id","phenotype"), k = 2) {
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfm <- inner_join(df, metadata, by = "sample_id")
  
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  Xtss <- safe_closure(X)
  Db <- vegan::vegdist(Xtss, method = "bray")
  
  set.seed(123)
  fit <- vegan::monoMDS(Db, k = k, maxit = 500, smin = 1e-12, trace = FALSE)
  sc  <- fit$points; colnames(sc) <- paste0("NMDS", seq_len(ncol(sc)))
  plot.df <- data.frame(sample_id = dfm$sample_id, as.data.frame(sc), check.names = FALSE)
  present <- intersect(model.vars, names(dfm))
  for (v in present) plot.df[[v]] <- factor(as.character(dfm[[v]]), levels = unique(as.character(metadata[[v]])))
  list(plot.df = plot.df, stress = fit$stress, used.vars = present)
}

# ==== NMDS panel (scatter with ellipses) with GLOBAL limit overrides ====
nmds_panel <- function(plot.df, stress, model.vars, axes = c(1,2),
                       label = NULL, xlim_override = NULL, ylim_override = NULL,
                       palette_name = "Batch") {
  mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                "#AEC7E8","#FFBB78","#98DF8A","#D62728","#FF9896","#C5B0D5",
                "#8C564B","#C49C94","#F7B6D2","#7F7F7F","#C7C7C7","#DBDB8D",
                "#17BECF","#9EDAE5")
  var.color <- model.vars[1]
  var.shape <- ifelse(length(model.vars) >= 2, model.vars[2], NA_character_)
  xcol <- paste0("NMDS", axes[1]); ycol <- paste0("NMDS", axes[2])
  
  # per-panel: ranges via 95% ellipses
  scores <- data.frame(
    AX1 = plot.df[[xcol]],
    AX2 = plot.df[[ycol]],
    batch = if (var.color %in% names(plot.df)) droplevels(plot.df[[var.color]]) else factor(1)
  )
  if (!is.factor(scores$batch)) scores$batch <- factor(scores$batch)
  eb <- ellipse_union_bounds(scores, "batch", level = 0.95, n = 240)
  xr <- safe_range(scores$AX1); yr <- safe_range(scores$AX2)
  xlim <- finite_range(xr, eb$x); ylim <- finite_range(yr, eb$y)
  pad_x <- safe_pad(xlim, 0.06); pad_y <- safe_pad(ylim, 0.06)
  xlim <- c(xlim[1] - pad_x, xlim[2] + pad_x)
  ylim <- c(ylim[1] - pad_y, ylim[2] + pad_y)
  if (!is.null(xlim_override)) xlim <- xlim_override
  if (!is.null(ylim_override)) ylim <- ylim_override
  
  title_text <- if (is.null(label)) sprintf("NMDS (stress=%.3f)", stress) else sprintf("%s - NMDS (stress=%.3f)", label, stress)
  
  p <- ggplot(plot.df, aes(x = !!sym(xcol), y = !!sym(ycol), colour = !!sym(var.color))) +
    {
      if (!is.na(var.shape) && var.shape %in% names(plot.df)) {
        geom_point(aes(shape = !!sym(var.shape)), size = 2, alpha = 0.9)
      } else {
        geom_point(size = 2, alpha = 0.9)
      }
    } +
    stat_ellipse(aes(group = !!sym(var.color)),
                 type = "norm", level = 0.95,
                 linewidth = 0.7, linetype = 1, show.legend = FALSE, na.rm = TRUE) +
    scale_color_manual(values = mbecCols, name = palette_name) +
    {
      if (!is.na(var.shape) && var.shape %in% names(plot.df)) {
        shape_vals <- c(0,1,2,3,6,8,15,16,17,23,25,4,5,9)
        nshape <- nlevels(plot.df[[var.shape]])
        if (nshape <= length(shape_vals)) {
          scale_shape_manual(values = shape_vals[seq_len(nshape)], name = "Phenotype")
        } else scale_shape_discrete(name = "Phenotype")
      }
    } +
    guides(colour = guide_legend(order = 1, nrow = 1, byrow = TRUE),
           shape  = guide_legend(order = 2, nrow = 1)) +
    labs(title = title_text, x = paste0("NMDS", axes[1]), y = paste0("NMDS", axes[2])) +
    scale_x_continuous(limits = xlim, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = ylim, expand = expansion(mult = c(0.02, 0.02))) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'bottom',
      legend.direction = 'horizontal',
      legend.box = 'vertical',
      plot.title = element_text(size = 10)
    )
  p
}

# ==== Params ====
batch_var  <- "batch_id"
shape_var  <- "phenotype"
model_vars <- if (is.na(shape_var)) c(batch_var) else c(batch_var, shape_var)
axes_to_plot <- c(1, 2)
ncol_grid <- 2
SYMMETRIC_AXES <- FALSE  # set TRUE to force symmetry about 0 (optional)

# =========================
# Set 1: NMDS - Aitchison (CLR)
# =========================
message("NMDS (Aitchison on CLR)")

frames_cache_clr <- list(); all_x <- NULL; all_y <- NULL
for (nm in names(file_list_clr)) {
  cat("Computing CLR NMDS:", nm, "\n")
  df <- read_csv(file_list_clr[[nm]], show_col_types = FALSE)
  fr <- compute_nmds_frames_aitch(df, metadata, model.vars = model_vars, k = max(axes_to_plot))
  frames_cache_clr[[nm]] <- fr
  
  xcol <- paste0("NMDS", axes_to_plot[1]); ycol <- paste0("NMDS", axes_to_plot[2])
  scores <- data.frame(
    AX1 = fr$plot.df[[xcol]], AX2 = fr$plot.df[[ycol]],
    batch = if (batch_var %in% names(fr$plot.df)) droplevels(fr$plot.df[[batch_var]]) else factor(1)
  )
  if (!is.factor(scores$batch)) scores$batch <- factor(scores$batch)
  eb <- ellipse_union_bounds(scores, "batch", level = 0.95, n = 240)
  all_x <- if (is.null(all_x)) finite_range(scores$AX1, eb$x) else finite_range(all_x, scores$AX1, eb$x)
  all_y <- if (is.null(all_y)) finite_range(scores$AX2, eb$y) else finite_range(all_y, scores$AX2, eb$y)
}
pad_x <- safe_pad(all_x, 0.06); pad_y <- safe_pad(all_y, 0.06)
xlim_clr <- c(all_x[1] - pad_x, all_x[2] + pad_x)
ylim_clr <- c(all_y[1] - pad_y, all_y[2] + pad_y)
if (isTRUE(SYMMETRIC_AXES)) { xh <- max(abs(xlim_clr)); yh <- max(abs(ylim_clr)); xlim_clr <- c(-xh, xh); ylim_clr <- c(-yh, yh) }

plots_clr <- lapply(names(file_list_clr), function(nm) {
  fr <- frames_cache_clr[[nm]]
  nmds_panel(fr$plot.df, fr$stress, model_vars,
             axes = axes_to_plot, label = nm,
             xlim_override = xlim_clr, ylim_override = ylim_clr, palette_name = "Batch")
})
names(plots_clr) <- names(file_list_clr)

# ---- Combine & save (CLR) ----
n_panels_clr <- length(plots_clr)
if (n_panels_clr == 1L) {
  combined_clr <- plots_clr[[1]] +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical",
          plot.margin = margin(8, 14, 8, 14))
  w_clr <- 9.5; h_clr <- 6
} else {
  combined_clr <- wrap_plots(plots_clr, ncol = ncol_grid) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical",
          plot.margin = margin(8, 14, 8, 14))
  w_clr <- 9.5 * min(ncol_grid, n_panels_clr)
  h_clr <- 6   * ceiling(n_panels_clr / ncol_grid)
}

ggsave(file.path(output_folder, "nmds_aitchison.png"),
       plot = combined_clr, width = w_clr, height = h_clr, dpi = 300)
ggsave(file.path(output_folder, "nmds_aitchison.tif"),
       plot = combined_clr, width = w_clr, height = h_clr, dpi = 300, compression = "lzw")

# =========================
# Set 2: NMDS - Bray-Curtis (TSS)
# =========================
message("NMDS (Bray-Curtis on TSS)")

frames_cache_tss <- list(); all_x <- NULL; all_y <- NULL
for (nm in names(file_list_tss)) {
  cat("Computing TSS NMDS:", nm, "\n")
  df <- read_csv(file_list_tss[[nm]], show_col_types = FALSE)
  fr <- compute_nmds_frames_bray(df, metadata, model.vars = model_vars, k = max(axes_to_plot))
  frames_cache_tss[[nm]] <- fr
  
  xcol <- paste0("NMDS", axes_to_plot[1]); ycol <- paste0("NMDS", axes_to_plot[2])
  scores <- data.frame(
    AX1 = fr$plot.df[[xcol]], AX2 = fr$plot.df[[ycol]],
    batch = if (batch_var %in% names(fr$plot.df)) droplevels(fr$plot.df[[batch_var]]) else factor(1)
  )
  if (!is.factor(scores$batch)) scores$batch <- factor(scores$batch)
  eb <- ellipse_union_bounds(scores, "batch", level = 0.95, n = 240)
  all_x <- if (is.null(all_x)) finite_range(scores$AX1, eb$x) else finite_range(all_x, scores$AX1, eb$x)
  all_y <- if (is.null(all_y)) finite_range(scores$AX2, eb$y) else finite_range(all_y, scores$AX2, eb$y)
}
pad_x <- safe_pad(all_x, 0.06); pad_y <- safe_pad(all_y, 0.06)
xlim_tss <- c(all_x[1] - pad_x, all_x[2] + pad_x)
ylim_tss <- c(all_y[1] - pad_y, all_y[2] + pad_y)
if (isTRUE(SYMMETRIC_AXES)) { xh <- max(abs(xlim_tss)); yh <- max(abs(ylim_tss)); xlim_tss <- c(-xh, xh); ylim_tss <- c(-yh, yh) }

plots_tss <- lapply(names(file_list_tss), function(nm) {
  fr <- frames_cache_tss[[nm]]
  nmds_panel(fr$plot.df, fr$stress, model_vars,
             axes = axes_to_plot, label = nm,
             xlim_override = xlim_tss, ylim_override = ylim_tss, palette_name = "Batch")
})
names(plots_tss) <- names(file_list_tss)

# ---- Combine & save (TSS) ----
n_panels_tss <- length(plots_tss)
if (n_panels_tss == 1L) {
  combined_tss <- plots_tss[[1]] +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical",
          plot.margin = margin(8, 14, 8, 14))
  w_tss <- 9.5; h_tss <- 6
} else {
  combined_tss <- wrap_plots(plots_tss, ncol = ncol_grid) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical",
          plot.margin = margin(8, 14, 8, 14))
  w_tss <- 9.5 * min(ncol_grid, n_panels_tss)
  h_tss <- 6   * ceiling(n_panels_tss / ncol_grid)
}

ggsave(file.path(output_folder, "nmds_braycurtis.png"),
       plot = combined_tss, width = w_tss, height = h_tss, dpi = 300)
ggsave(file.path(output_folder, "nmds_braycurtis.tif"),
       plot = combined_tss, width = w_tss, height = h_tss, dpi = 300, compression = "lzw")

# =========================
# NMDS ranking (with baseline-only assessment)
# =========================

# --- helpers ---
compute_centroids_nmds <- function(scores, batch_var = "batch_id") {
  scores %>%
    dplyr::group_by(!!rlang::sym(batch_var)) %>%
    dplyr::summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")
}
compute_centroid_distances <- function(centroids) {
  if (nrow(centroids) < 2) return(NA_real_)
  as.numeric(mean(dist(centroids[, c("NMDS1", "NMDS2")], method = "euclidean")))
}
compute_within_dispersion_nmds <- function(scores, batch_var = "batch_id") {
  if (!all(c("NMDS1","NMDS2", batch_var) %in% names(scores))) return(NA_real_)
  levs <- levels(scores[[batch_var]]); if (is.null(levs)) levs <- unique(scores[[batch_var]])
  ws <- c(); ns <- c()
  for (lev in levs) {
    sub <- scores[scores[[batch_var]] == lev, c("NMDS1","NMDS2"), drop = FALSE]
    n <- nrow(sub)
    if (n < 2) next
    d <- stats::dist(sub, method = "euclidean")
    ws <- c(ws, mean(d)); ns <- c(ns, n)
  }
  if (!length(ws)) return(NA_real_)
  stats::weighted.mean(ws, w = ns)
}

# then combine via geometric mean (higher = better).
nmds_metric_score <- function(batch_distance, stress, stress_cap = 0.30) {
  if (is.na(batch_distance) || is.na(stress)) return(NA_real_)
  S_batch  <- 1 / (1 + batch_distance)                           
  S_stress <- pmax(0, 1 - pmin(stress, stress_cap) / stress_cap) 
  sqrt(S_batch * S_stress)
}

methods_clr <- names(frames_cache_clr)
methods_tss <- names(frames_cache_tss)
all_methods <- union(methods_clr, methods_tss)

only_baseline <- (length(all_methods) == 1L) && identical(all_methods, "Before correction")

if (only_baseline) {
  # ===== Baseline-only assessment (no ranking) =====
  assess_rows <- list()
  
  if ("Before correction" %in% methods_clr) {
    fr <- frames_cache_clr[["Before correction"]]
    md <- metadata[match(fr$plot.df$sample_id, metadata$sample_id), , drop = FALSE]
    scores <- fr$plot.df %>%
      dplyr::select(sample_id, NMDS1, NMDS2) %>%
      dplyr::mutate(batch_id = factor(md$batch_id))
    cents <- compute_centroids_nmds(scores, "batch_id")
    D_between <- compute_centroid_distances(cents)
    W_within  <- compute_within_dispersion_nmds(scores, "batch_id")
    stress    <- fr$stress
    score     <- nmds_metric_score(D_between, stress)
    needs_correction <- is.finite(D_between) && is.finite(W_within) && (D_between > W_within)
    
    assess_rows[["CLR"]] <- tibble::tibble(
      Method            = "Before correction",
      Geometry          = "Aitchison (CLR)",
      Batch_Distance    = D_between,
      Within_Dispersion = W_within,
      NMDS_Stress       = stress,
      Score             = score,
      Needs_Correction  = needs_correction
    )
  }
  
  if ("Before correction" %in% methods_tss) {
    fr <- frames_cache_tss[["Before correction"]]
    md <- metadata[match(fr$plot.df$sample_id, metadata$sample_id), , drop = FALSE]
    scores <- fr$plot.df %>%
      dplyr::select(sample_id, NMDS1, NMDS2) %>%
      dplyr::mutate(batch_id = factor(md$batch_id))
    cents <- compute_centroids_nmds(scores, "batch_id")
    D_between <- compute_centroid_distances(cents)
    W_within  <- compute_within_dispersion_nmds(scores, "batch_id")
    stress    <- fr$stress
    score     <- nmds_metric_score(D_between, stress)
    needs_correction <- is.finite(D_between) && is.finite(W_within) && (D_between > W_within)
    
    assess_rows[["TSS"]] <- tibble::tibble(
      Method            = "Before correction",
      Geometry          = "Bray-Curtis (TSS)",
      Batch_Distance    = D_between,
      Within_Dispersion = W_within,
      NMDS_Stress       = stress,
      Score             = score,
      Needs_Correction  = needs_correction
    )
  }
  
  assess_df <- dplyr::bind_rows(assess_rows)
  
  # Combined view (geometric mean of available scores; OR on correction flags)
  comb_score <- dplyr::case_when(
    nrow(assess_df) >= 2 && all(is.finite(assess_df$Score)) ~ sqrt(prod(assess_df$Score)),
    TRUE                                                    ~ max(assess_df$Score, na.rm = TRUE)
  )
  needs_corr_any <- any(assess_df$Needs_Correction, na.rm = TRUE)
  
  assess_df <- dplyr::bind_rows(
    assess_df,
    tibble::tibble(
      Method            = "Before correction",
      Geometry          = "Combined",
      Batch_Distance    = NA_real_,
      Within_Dispersion = NA_real_,
      NMDS_Stress       = NA_real_,
      Score             = comb_score,
      Needs_Correction  = needs_corr_any
    )
  )
  
  print(assess_df, n = nrow(assess_df))
  readr::write_csv(assess_df, file.path(output_folder, "nmds_raw_assessment.csv"))
  
  message(if (isTRUE(needs_corr_any)) {
    "Assessment: Batch separation exceeds within-batch spread in at least one geometry - correction recommended."
  } else {
    "Assessment: Batch separation does not exceed within-batch spread - correction may not be necessary."
  })
  
} else {
  # ===== Multi-method ranking =====
  rank_tbl <- dplyr::tibble(
    Method = character(),
    Batch_Distance_Aitchison = numeric(),
    Batch_Distance_Bray      = numeric(),
    NMDS_Stress_Aitchison    = numeric(),
    NMDS_Stress_Bray         = numeric(),
    Score_Aitchison          = numeric(),
    Score_Bray               = numeric(),
    Combined_Score           = numeric()
  )
  
  for (m in all_methods) {
    # Aitchison (CLR NMDS)
    D_a <- NA_real_; S_a <- NA_real_; stress_a <- NA_real_
    if (m %in% methods_clr) {
      fr_a <- frames_cache_clr[[m]]
      md_a <- metadata[match(fr_a$plot.df$sample_id, metadata$sample_id), , drop = FALSE]
      scores_a <- fr_a$plot.df %>%
        dplyr::select(sample_id, NMDS1, NMDS2) %>%
        dplyr::mutate(batch_id = factor(md_a$batch_id))
      cents_a <- compute_centroids_nmds(scores_a, "batch_id")
      D_a <- compute_centroid_distances(cents_a)
      stress_a <- fr_a$stress
      S_a <- nmds_metric_score(D_a, stress_a)
    }
    
    # Bray-Curtis (TSS NMDS)
    D_b <- NA_real_; S_b <- NA_real_; stress_b <- NA_real_
    if (m %in% methods_tss) {
      fr_b <- frames_cache_tss[[m]]
      md_b <- metadata[match(fr_b$plot.df$sample_id, metadata$sample_id), , drop = FALSE]
      scores_b <- fr_b$plot.df %>%
        dplyr::select(sample_id, NMDS1, NMDS2) %>%
        dplyr::mutate(batch_id = factor(md_b$batch_id))
      cents_b <- compute_centroids_nmds(scores_b, "batch_id")
      D_b <- compute_centroid_distances(cents_b)
      stress_b <- fr_b$stress
      S_b <- nmds_metric_score(D_b, stress_b)
    }
    
    # Combine available NMDS metric scores (geometric mean; higher = better)
    Combined <- dplyr::case_when(
      !is.na(S_a) && !is.na(S_b) ~ sqrt(S_a * S_b),
      !is.na(S_a)                ~ S_a,
      !is.na(S_b)                ~ S_b,
      TRUE                       ~ NA_real_
    )
    
    rank_tbl <- dplyr::bind_rows(rank_tbl, dplyr::tibble(
      Method = m,
      Batch_Distance_Aitchison = D_a,
      Batch_Distance_Bray      = D_b,
      NMDS_Stress_Aitchison    = stress_a,
      NMDS_Stress_Bray         = stress_b,
      Score_Aitchison          = S_a,
      Score_Bray               = S_b,
      Combined_Score           = Combined
    ))
  }
  
  ranked_nmds_only <- rank_tbl %>%
    dplyr::arrange(dplyr::desc(Combined_Score)) %>%
    dplyr::mutate(Rank = dplyr::row_number())
  
  print(ranked_nmds_only, n = nrow(ranked_nmds_only))
  readr::write_csv(ranked_nmds_only, file.path(output_folder, "nmds_ranking.csv"))
}