# ==== Load Required Libraries ====
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(patchwork)  # legend collecting & layout
  library(rlang)
})

# ---- helpers ----
mbecUpperCase <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

guess_shape_var <- function(meta, batch_col = "batchid") {
  cand <- c("group","phenotype","condition","status","case_control","class","disease","label")
  hit <- cand[cand %in% names(meta)]
  if (length(hit)) return(hit[1])
  facs <- names(Filter(function(v) {
    !is.numeric(meta[[v]]) && !is.integer(meta[[v]]) &&
      length(unique(meta[[v]])) <= 12
  }, meta))
  facs <- setdiff(facs, c("sample_id", batch_col))
  if (length(facs)) return(facs[1])
  return(NA_character_)
}

# union bounds of 95% normal ellipses across groups
ellipse_union_bounds <- function(df_scores, group_var, level = 0.95, n = 240) {
  if (!nrow(df_scores)) return(list(x = c(0,0), y = c(0,0)))
  chi <- sqrt(qchisq(level, df = 2))
  xmins <- c(); xmaxs <- c(); ymins <- c(); ymaxs <- c()
  for (lev in levels(df_scores[[group_var]])) {
    sub <- df_scores[df_scores[[group_var]] == lev, c("PCX","PCY"), drop = FALSE]
    if (nrow(sub) < 3 || any(!is.finite(as.matrix(sub)))) next
    S <- tryCatch(stats::cov(sub, use = "complete.obs"), error = function(e) NULL)
    mu <- colMeans(sub, na.rm = TRUE)
    if (is.null(S) || any(!is.finite(S))) next
    eg <- eigen(S, symmetric = TRUE)
    if (any(!is.finite(eg$values))) next
    t <- seq(0, 2*pi, length.out = n)
    R <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)))
    pts <- t(chi * R %*% rbind(cos(t), sin(t)))
    pts <- sweep(pts, 2, mu, FUN = "+")
    xmins <- c(xmins, min(pts[,1], na.rm = TRUE))
    xmaxs <- c(xmaxs, max(pts[,1], na.rm = TRUE))
    ymins <- c(ymins, min(pts[,2], na.rm = TRUE))
    ymaxs <- c(ymaxs, max(pts[,2], na.rm = TRUE))
  }
  list(x = c(min(xmins, na.rm = TRUE), max(xmaxs, na.rm = TRUE)),
       y = c(min(ymins, na.rm = TRUE), max(ymaxs, na.rm = TRUE)))
}

# ==== Read UID argument ====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# ==== Read Metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))

# --------- Collect CLR files ---------
clr_paths <- list.files(output_folder, pattern = "^normalized_.*_clr\\.csv$", full.names = TRUE)

# include raw_clr.csv (as baseline) if present
raw_clr_fp <- file.path(output_folder, "raw_clr.csv")
if (file.exists(raw_clr_fp)) clr_paths <- c(raw_clr_fp, clr_paths)

if (!length(clr_paths)) stop("No CLR matrices found (expected 'raw_clr.csv' or 'normalized_*_clr.csv').")

method_names <- ifelse(basename(clr_paths) == "raw_clr.csv",
                       "Before correction",
                       gsub("^normalized_|_clr\\.csv$", "", basename(clr_paths)))
file_list <- setNames(clr_paths, method_names)

# ==== build plot.df + metric.df like mbecPCA would ====
# ensures consistent factor levels across panels so legend can be shared
compute_pca_frames <- function(df, metadata, model.vars = c("batchid","group"), n_pcs = 5) {
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfm <- inner_join(df, metadata, by = "sample_id")
  
  feature_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feature_cols)) %>% select(where(is.numeric))
  keep <- vapply(X, function(x) sd(x, na.rm = TRUE) > 0, logical(1))
  X <- as.matrix(X[, keep, drop = FALSE])
  if (!ncol(X)) stop("No variable numeric features remain for PCA.")
  
  # NOTE: If your inputs are CLR already, consider scale. = FALSE for strict Aitchison PCA.
  pc <- prcomp(X, center = TRUE, scale. = TRUE)
  k  <- min(n_pcs, ncol(pc$x))
  
  ve       <- (pc$sdev^2) / sum(pc$sdev^2)
  axis.min <- apply(pc$x[, seq_len(k), drop = FALSE], 2, min)
  axis.max <- apply(pc$x[, seq_len(k), drop = FALSE], 2, max)
  
  metric.df <- data.frame(
    axis.min      = as.numeric(axis.min),
    axis.max      = as.numeric(axis.max),
    var.explained = round(100 * ve[seq_len(k)], 2),
    stringsAsFactors = FALSE
  )
  
  plot.df <- data.frame(
    sample_id = dfm$sample_id,
    as.data.frame(pc$x[, seq_len(k), drop = FALSE]),
    stringsAsFactors = FALSE
  )
  colnames(plot.df)[2:(k + 1)] <- paste0("PC", seq_len(k))
  
  # attach *consistent* factor levels from metadata so legends match across panels
  present <- intersect(model.vars, names(dfm))
  for (v in present) {
    levs <- unique(as.character(metadata[[v]]))
    plot.df[[v]] <- factor(as.character(dfm[[v]]), levels = levs)
  }
  
  list(plot.df = plot.df, metric.df = metric.df, used.vars = present)
}

# ==== panel: scatter + marginal densities; legend kept (not collected here) ====
mbecPCAPlot <- function(plot.df, metric.df, model.vars, pca.axes, label=NULL) {
  
  mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                "#AEC7E8","#FFBB78","#98DF8A","#D62728","#FF9896","#C5B0D5",
                "#8C564B","#C49C94","#F7B6D2","#7F7F7F","#C7C7C7","#DBDB8D",
                "#17BECF","#9EDAE5")
  
  var.color <- model.vars[1]
  var.shape <- ifelse(length(model.vars) >= 2, model.vars[2], NA_character_)
  
  xcol <- colnames(plot.df)[pca.axes[1] + 1]
  ycol <- colnames(plot.df)[pca.axes[2] + 1]
  
  x.label <- paste0(xcol, ": ", metric.df$var.explained[pca.axes[1]], "% expl.var")
  y.label <- paste0(ycol, ": ", metric.df$var.explained[pca.axes[2]], "% expl.var")
  if (!is.null(label)) x.label <- paste(label, "-", x.label)
  
  # compute zoomed-out limits using ellipse bounds + points
  scores <- data.frame(
    PCX = plot.df[[xcol]],
    PCY = plot.df[[ycol]],
    batch = if (var.color %in% names(plot.df)) droplevels(plot.df[[var.color]]) else factor(1)
  )
  ell_bounds <- ellipse_union_bounds(scores, "batch", level = 0.95, n = 240)
  xr <- range(scores$PCX, na.rm = TRUE); yr <- range(scores$PCY, na.rm = TRUE)
  xlim <- range(c(xr, ell_bounds$x)); ylim <- range(c(yr, ell_bounds$y))
  pad_x <- diff(xlim) * 0.06; pad_y <- diff(ylim) * 0.06
  xlim <- c(xlim[1] - pad_x, xlim[2] + pad_x)
  ylim <- c(ylim[1] - pad_y, ylim[2] + pad_y)
  
  pmar <- margin(10, 16, 10, 16)
  
  # main scatter (legend source)
  pMain <- ggplot(plot.df, aes(x = !!sym(xcol), y = !!sym(ycol), colour = !!sym(var.color))) +
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
    scale_color_manual(values = mbecCols, name = "Batch") +
    {
      if (!is.na(var.shape) && var.shape %in% names(plot.df)) {
        shape_vals <- c(0,1,2,3,6,8,15,16,17,23,25,4,5,9)
        nshape <- nlevels(plot.df[[var.shape]])
        if (nshape <= length(shape_vals)) {
          scale_shape_manual(values = shape_vals[seq_len(nshape)], name = "Phenotype")
        } else scale_shape_discrete(name = "Phenotype")
      }
    } +
    guides(
      colour = guide_legend(order = 1, nrow = 1, byrow = TRUE),  # many batches on one row
      shape  = guide_legend(order = 2, nrow = 1)
    ) +
    labs(title = NULL) +
    scale_x_continuous(limits = xlim, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = ylim, expand = expansion(mult = c(0.02, 0.02))) +
    xlab(x.label) + ylab(y.label) + theme_bw() +
    theme(
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'bottom',
      legend.direction = 'horizontal',
      legend.box = 'vertical',
      plot.margin = pmar
    )
  
  # top density (PC1) — no legend
  pTop <- ggplot(plot.df, aes(x = !!sym(xcol))) +
    {
      if (!is.na(var.shape) && var.shape %in% names(plot.df)) {
        geom_density(aes(fill = !!sym(var.color), linetype = !!sym(var.shape)),
                     linewidth = 0.3, alpha = 0.5, show.legend = FALSE)
      } else {
        geom_density(aes(fill = !!sym(var.color)),
                     linewidth = 0.3, alpha = 0.5, show.legend = FALSE)
      }
    } +
    ylab("Density") +
    scale_fill_manual(values = mbecCols, guide = "none") +
    scale_x_continuous(limits = xlim, expand = expansion(mult = c(0.02, 0.02))) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'none',
      axis.title.x = element_blank(),
      plot.margin = pmar
    )
  
  # right density (PC2)
  pRight <- ggplot(plot.df, aes(y = !!sym(ycol))) +
    {
      if (!is.na(var.shape) && var.shape %in% names(plot.df)) {
        geom_density(
          aes(x = after_stat(density), fill = !!sym(var.color), linetype = !!sym(var.shape)),
          linewidth = 0.3, alpha = 0.5, orientation = "y", show.legend = FALSE
        )
      } else {
        geom_density(
          aes(x = after_stat(density), fill = !!sym(var.color)),
          linewidth = 0.3, alpha = 0.5, orientation = "y", show.legend = FALSE
        )
      }
    } +
    xlab(NULL) + ylab("Density") +
    scale_fill_manual(values = mbecCols, guide = "none") +
    scale_y_continuous(limits = ylim, expand = expansion(mult = c(0.02, 0.02))) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      axis.title.y = element_blank(),
      plot.margin = margin(10, 16, 10, 16)
    )
  
  # assemble (DON'T collect here; we'll collect once globally)
  design <- "
A#
CB
"
  (pTop + pRight + pMain) + plot_layout(design = design, widths = c(3, 1), heights = c(1.6, 3.2))
}

# ==== choose covariates (auto-detect shape var) ====
batch_var  <- "batchid"
shape_var  <- guess_shape_var(metadata, batch_var)
model_vars <- if (is.na(shape_var)) c(batch_var) else c(batch_var, shape_var)
message(sprintf("Using color=%s%s",
                batch_var,
                if (is.na(shape_var)) " (shape: none)" else paste0(", shape=", shape_var)))

# ==== build + save PCA panels for each matrix (single legend across ALL) ====
pcs_to_plot <- c(1, 2)
frames_cache <- list()   # <--- build this for ranking
plots <- list()

for (nm in names(file_list)) {
  cat("Processing:", nm, "\n")
  df <- read_csv(file_list[[nm]], show_col_types = FALSE)
  
  frames <- compute_pca_frames(df, metadata, model.vars = model_vars, n_pcs = max(pcs_to_plot))
  frames_cache[[nm]] <- frames
  
  plt <- mbecPCAPlot(
    plot.df   = frames$plot.df,
    metric.df = frames$metric.df,
    model.vars = frames$used.vars,
    pca.axes  = pcs_to_plot,
    label     = nm
  )
  
  plots[[nm]] <- plt
}

# ---- Combine ALL panels and keep ONLY ONE legend at the bottom (horizontal) ----
ncol_grid <- 2
n_panels  <- length(plots)
if (n_panels == 0L) stop("No PCA panels to plot.")

if (n_panels == 1L) {
  combined <- plots[[1]] +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "vertical",
      plot.margin      = margin(8, 14, 8, 14)
    )
  w <- 9.5; h <- 6
} else {
  combined <- wrap_plots(plots, ncol = ncol_grid) +
    plot_layout(guides = "collect") &
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "vertical",
      plot.margin      = margin(8, 14, 8, 14)
    )
  w <- 9.5 * min(ncol_grid, n_panels)
  h <- 6   * ceiling(n_panels / ncol_grid)
}

ggsave(file.path(output_folder, "pca.png"),
       plot = combined, width = w, height = h, dpi = 300)
ggsave(file.path(output_folder, "pca.tif"),
       plot = combined, width = w, height = h, dpi = 300, compression = "lzw")

# =========================
# PCA ranking / assessment
# =========================

compute_centroids_pca <- function(scores, batch_var = "batchid") {
  scores %>%
    dplyr::group_by(!!rlang::sym(batch_var)) %>%
    dplyr::summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
}
compute_centroid_distances <- function(centroids) {
  if (nrow(centroids) < 2) return(NA_real_)
  as.numeric(mean(dist(centroids[, c("PC1","PC2")], method = "euclidean")))
}
# Weighted average within-batch dispersion on PC1–PC2 (higher = batches internally more spread)
compute_within_dispersion <- function(scores, batch_var = "batchid") {
  if (!all(c("PC1","PC2", batch_var) %in% names(scores))) return(NA_real_)
  levs <- levels(scores[[batch_var]])
  if (is.null(levs)) levs <- unique(scores[[batch_var]])
  ws <- c(); ns <- c()
  for (lev in levs) {
    sub <- scores[scores[[batch_var]] == lev, c("PC1","PC2"), drop = FALSE]
    n <- nrow(sub)
    if (n < 2) next
    d <- stats::dist(sub, method = "euclidean")
    ws <- c(ws, mean(d))
    ns <- c(ns, n)
  }
  if (!length(ws)) return(NA_real_)
  stats::weighted.mean(ws, w = ns)
}
pca_metric_score <- function(batch_distance, coverage) {
  if (is.na(batch_distance) || is.na(coverage)) return(NA_real_)
  coverage <- max(0, min(1, coverage))
  (1 / (1 + batch_distance)) * coverage  # higher = better
}

only_baseline <- (length(frames_cache) == 1L) && identical(names(frames_cache), "Before correction")

if (only_baseline) {
  # ---- Assess RAW only (no ranking) ----
  m <- "Before correction"
  fr <- frames_cache[[m]]
  stopifnot(!is.null(fr))
  
  # Build scores with batch
  md <- metadata[match(fr$plot.df$sample_id, metadata$sample_id), , drop = FALSE]
  scores <- fr$plot.df %>%
    dplyr::select(sample_id, PC1, PC2) %>%
    dplyr::mutate(batchid = factor(md$batchid))
  
  cents <- compute_centroids_pca(scores, "batchid")
  D_between <- compute_centroid_distances(cents)
  W_within  <- compute_within_dispersion(scores, "batchid")
  
  ve <- fr$metric.df$var.explained
  coverage2 <- sum(ve[1:min(2, length(ve))], na.rm = TRUE) / 100
  score <- pca_metric_score(D_between, coverage2)
  
  # Simple, data-driven decision rule:
  # If batches are separated more than their average within-batch spread -> likely batch effect.
  needs_correction <- is.finite(D_between) && is.finite(W_within) && (D_between > W_within)
  
  assess_df <- tibble::tibble(
    Method            = m,
    Batch_Distance    = D_between,
    Within_Dispersion = W_within,
    Coverage_PC1_PC2  = coverage2,
    Score             = score,
    Needs_Correction  = needs_correction
  )
  
  print(assess_df)
  readr::write_csv(assess_df, file.path(output_folder, "pca_raw_assessment.csv"))
  
  msg <- if (isTRUE(needs_correction)) {
    "Assessment: Batch separation exceeds within-batch spread — correction recommended."
  } else {
    "Assessment: Batch separation does not exceed within-batch spread — correction may not be necessary."
  }
  message(msg)
  
} else {
  # ---- Rank multiple matrices (as before) ----
  rank_tbl <- dplyr::tibble(Method = character(),
                            Batch_Distance = numeric(),
                            Coverage = numeric(),
                            Score = numeric())
  
  for (m in names(frames_cache)) {
    fr <- frames_cache[[m]]
    if (!all(c("PC1","PC2","sample_id") %in% names(fr$plot.df))) next
    md <- metadata[match(fr$plot.df$sample_id, metadata$sample_id), , drop = FALSE]
    scores <- fr$plot.df %>%
      dplyr::select(sample_id, PC1, PC2) %>%
      dplyr::mutate(batchid = factor(md$batchid))
    cents <- compute_centroids_pca(scores, "batchid")
    D <- compute_centroid_distances(cents)
    ve <- fr$metric.df$var.explained
    cov2 <- sum(ve[1:min(2, length(ve))], na.rm = TRUE) / 100
    
    rank_tbl <- dplyr::bind_rows(rank_tbl, dplyr::tibble(
      Method = m,
      Batch_Distance = D,
      Coverage = cov2,
      Score = pca_metric_score(D, cov2)
    ))
  }
  
  ranked_pca <- rank_tbl %>%
    dplyr::arrange(dplyr::desc(Score)) %>%
    dplyr::mutate(Rank = dplyr::row_number())
  
  print(ranked_pca, n = nrow(ranked_pca))
  readr::write_csv(ranked_pca, file.path(output_folder, "pca_ranking.csv"))
}
