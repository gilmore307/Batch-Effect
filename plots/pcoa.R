# ==== Libraries ====
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(patchwork)  # layouts + legend collection
  library(rlang)
  library(vegan)      # Bray–Curtis
})

# ==== Helpers (robust ranges/padding) ====
safe_range <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) c(0, 0) else range(v, na.rm = TRUE)
}
finite_range <- function(...) {
  v <- unlist(list(...))
  v <- v[is.finite(v)]
  if (!length(v)) c(0, 0) else range(v, na.rm = TRUE)
}
safe_pad <- function(r, frac = 0.06) {
  dx <- diff(r)
  if (!is.finite(dx) || dx <= 0) return(1e-6) else dx * frac
}

# ==== Ellipse union bounds (falls back to point ranges if needed) ====
ellipse_union_bounds <- function(df_scores, group_var, level = 0.95, n = 240) {
  if (!nrow(df_scores)) return(list(x = c(0,0), y = c(0,0)))
  chi <- sqrt(qchisq(level, df = 2))
  xmins <- c(); xmaxs <- c(); ymins <- c(); ymaxs <- c()
  levs <- levels(df_scores[[group_var]])
  if (is.null(levs)) levs <- unique(df_scores[[group_var]])
  for (lev in levs) {
    sub <- df_scores[df_scores[[group_var]] == lev, c("AX1","AX2"), drop = FALSE]
    if (nrow(sub) < 3 || any(!is.finite(as.matrix(sub)))) next
    S  <- tryCatch(stats::cov(sub, use = "complete.obs"), error = function(e) NULL)
    mu <- colMeans(sub, na.rm = TRUE)
    if (is.null(S) || any(!is.finite(S))) next
    eg <- eigen(S, symmetric = TRUE)
    if (any(!is.finite(eg$values))) next
    t  <- seq(0, 2*pi, length.out = n)
    R  <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)))
    pts <- t(chi * R %*% rbind(cos(t), sin(t)))
    pts <- sweep(pts, 2, mu, FUN = "+")
    xmins <- c(xmins, min(pts[,1], na.rm = TRUE))
    xmaxs <- c(xmaxs, max(pts[,1], na.rm = TRUE))
    ymins <- c(ymins, min(pts[,2], na.rm = TRUE))
    ymaxs <- c(ymaxs, max(pts[,2], na.rm = TRUE))
  }
  if (!length(xmins) || !length(ymins)) {
    return(list(x = safe_range(df_scores$AX1), y = safe_range(df_scores$AX2)))
  }
  list(x = c(min(xmins, na.rm = TRUE), max(xmaxs, na.rm = TRUE)),
       y = c(min(ymins, na.rm = TRUE), max(ymaxs, na.rm = TRUE)))
}

# ==== Compositional helpers ====
safe_closure <- function(X) {
  rs <- rowSums(X, na.rm = TRUE)
  zero_rows <- which(rs == 0 | !is.finite(rs))
  if (length(zero_rows)) {
    X[zero_rows, ] <- 1 / ncol(X)
    rs <- rowSums(X, na.rm = TRUE)
  }
  sweep(X, 1, rs, "/")
}
clr_transform <- function(X) {
  Xc <- safe_closure(X)
  for (i in seq_len(nrow(Xc))) {
    xi <- Xc[i, ]
    pos <- xi > 0 & is.finite(xi)
    if (!any(pos)) {
      xi[] <- 1 / length(xi); pos <- xi > 0
    }
    if (any(!pos)) {
      minpos <- min(xi[pos], na.rm = TRUE)
      repl <- min(minpos * 0.5, 1e-8)
      xi[!pos] <- repl
      xi <- xi / sum(xi)
    }
    Xc[i, ] <- xi
  }
  L <- log(Xc)
  sweep(L, 1, rowMeans(L), "-")
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
file_list_clr <- setNames(clr_paths, if (length(clr_paths)) name_from(clr_paths, "clr") else character())
file_list_tss <- setNames(tss_paths, if (length(tss_paths)) name_from(tss_paths, "tss") else character())

# Include raw_clr.csv / raw_tss.csv as "Before correction" if present
raw_clr_fp <- file.path(output_folder, "raw_clr.csv")
raw_tss_fp <- file.path(output_folder, "raw_tss.csv")
if (file.exists(raw_clr_fp)) file_list_clr <- c("Before correction" = raw_clr_fp, file_list_clr)
if (file.exists(raw_tss_fp)) file_list_tss <- c("Before correction" = raw_tss_fp, file_list_tss)

if (!length(file_list_clr) && !length(file_list_tss)) {
  stop("No normalized files found (expected raw_clr.csv/raw_tss.csv and/or normalized_*_clr.csv / normalized_*_tss.csv) in ", output_folder)
}

# ==== PCoA frames (Aitchison on CLR) ====
compute_pcoa_frames_aitch <- function(df, metadata, model.vars = c("batchid","phenotype"),
                                      n_axes = 5) {
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfm <- inner_join(df, metadata, by = "sample_id")
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feat_cols)) %>% select(where(is.numeric))
  keep <- vapply(X, function(x) sd(x, na.rm = TRUE) > 0, logical(1))
  X <- as.matrix(X[, keep, drop = FALSE])
  if (!ncol(X)) stop("No variable numeric features remain for PCoA (CLR).")
  
  # If negatives exist, assume already CLR; else do CLR from positive data
  has_neg <- any(X < 0, na.rm = TRUE)
  Xclr <- if (has_neg) sweep(X, 1, rowMeans(X, na.rm = TRUE), "-") else clr_transform(X)
  d <- dist(Xclr, method = "euclidean")
  
  # PCoA
  if (requireNamespace("ape", quietly = TRUE)) {
    pc <- ape::pcoa(d)
    k  <- min(n_axes, ncol(pc$vectors))
    sc <- pc$vectors[, 1:k, drop = FALSE]
    vals <- pc$values
    rel  <- if ("Rel_corr_eig" %in% names(vals)) vals$Rel_corr_eig else vals$Relative_eig
    var_expl <- round(100 * rel[1:k], 2)
  } else {
    cm <- cmdscale(d, k = n_axes, eig = TRUE, add = TRUE)
    sc <- cm$points
    k  <- ncol(sc)
    rel_eig <- cm$eig / sum(cm$eig[cm$eig > 0], na.rm = TRUE)
    var_expl <- round(100 * rel_eig[1:k], 2)
  }
  
  axis.min <- apply(sc, 2, min, na.rm = TRUE)
  axis.max <- apply(sc, 2, max, na.rm = TRUE)
  metric.df <- data.frame(axis.min = as.numeric(axis.min),
                          axis.max = as.numeric(axis.max),
                          var.explained = as.numeric(var_expl))
  plot.df <- data.frame(sample_id = dfm$sample_id, as.data.frame(sc), check.names = FALSE)
  colnames(plot.df)[2:(k + 1)] <- paste0("PCo", seq_len(k))
  
  present <- intersect(model.vars, names(dfm))
  for (v in present) {
    levs <- unique(as.character(metadata[[v]]))
    plot.df[[v]] <- factor(as.character(dfm[[v]]), levels = levs)
  }
  list(plot.df = plot.df, metric.df = metric.df, used.vars = present)
}

# ==== PCoA frames (Bray–Curtis on TSS) ====
compute_pcoa_frames_bray <- function(df, metadata, model.vars = c("batchid","phenotype"),
                                     n_axes = 5) {
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfm <- inner_join(df, metadata, by = "sample_id")
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  
  # Ensure non-negative & close rows to TSS (proportions)
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0
  Xtss <- safe_closure(X)
  
  # Bray–Curtis (vegan::vegdist)
  d <- vegan::vegdist(Xtss, method = "bray")
  
  # PCoA with Cailliez correction (robust for non-Euclidean)
  if (requireNamespace("ape", quietly = TRUE)) {
    pc <- ape::pcoa(as.matrix(d), correction = "cailliez")
    k  <- min(n_axes, ncol(pc$vectors))
    sc <- pc$vectors[, 1:k, drop = FALSE]
    vals <- pc$values
    rel  <- if ("Rel_corr_eig" %in% names(vals)) vals$Rel_corr_eig else vals$Relative_eig
    var_expl <- round(100 * rel[1:k], 2)
  } else {
    cm <- cmdscale(as.matrix(d), k = n_axes, eig = TRUE, add = TRUE)
    sc <- cm$points
    k  <- ncol(sc)
    rel_eig <- cm$eig / sum(cm$eig[cm$eig > 0], na.rm = TRUE)
    var_expl <- round(100 * rel_eig[1:k], 2)
  }
  
  axis.min <- apply(sc, 2, min, na.rm = TRUE)
  axis.max <- apply(sc, 2, max, na.rm = TRUE)
  metric.df <- data.frame(axis.min = as.numeric(axis.min),
                          axis.max = as.numeric(axis.max),
                          var.explained = as.numeric(var_expl))
  plot.df <- data.frame(sample_id = dfm$sample_id, as.data.frame(sc), check.names = FALSE)
  colnames(plot.df)[2:(k + 1)] <- paste0("PCo", seq_len(k))
  
  present <- intersect(model.vars, names(dfm))
  for (v in present) {
    levs <- unique(as.character(metadata[[v]]))
    plot.df[[v]] <- factor(as.character(dfm[[v]]), levels = levs)
  }
  list(plot.df = plot.df, metric.df = metric.df, used.vars = present)
}

# ==== PCoA panel (scatter + marginals) with GLOBAL limit overrides ====
pcoa_panel <- function(plot.df, metric.df, model.vars, axes = c(1,2), label = NULL,
                       xlim_override = NULL, ylim_override = NULL, palette_name = "Batch") {
  mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                "#AEC7E8","#FFBB78","#98DF8A","#D62728","#FF9896","#C5B0D5",
                "#8C564B","#C49C94","#F7B6D2","#7F7F7F","#C7C7C7","#DBDB8D",
                "#17BECF","#9EDAE5")
  var.color <- model.vars[1]
  var.shape <- ifelse(length(model.vars) >= 2, model.vars[2], NA_character_)
  xcol <- paste0("PCo", axes[1]); ycol <- paste0("PCo", axes[2])
  x.label <- paste0(xcol, ": ", metric.df$var.explained[axes[1]], "% expl.var")
  y.label <- paste0(ycol, ": ", metric.df$var.explained[axes[2]], "% expl.var")
  if (!is.null(label)) x.label <- paste(label, "-", x.label)
  
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
  
  pmar <- margin(10, 16, 10, 16)
  
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
  
  design <- "
A#
CB
"
  (pTop + pRight + pMain) + plot_layout(design = design, widths = c(3, 1), heights = c(1.6, 3.2))
}

# ==== Params ====
batch_var  <- "batchid"
shape_var  <- "phenotype"
model_vars <- if (is.na(shape_var)) c(batch_var) else c(batch_var, shape_var)
axes_to_plot <- c(1, 2)
ncol_grid <- 2
SYMMETRIC_AXES <- FALSE  # set TRUE to force symmetry about 0 (optional)

# =========================
# Set 1: Aitchison (CLR)
# =========================
message(sprintf("PCoA (Aitchison on CLR): color=%s%s",
                batch_var,
                if (is.na(shape_var)) " (shape: none)" else paste0(", shape=", shape_var)))

frames_cache_clr <- list()
global_x_clr <- NULL; global_y_clr <- NULL

for (nm in names(file_list_clr)) {
  cat("Computing CLR frames:", nm, "\n")
  df <- read_csv(file_list_clr[[nm]], show_col_types = FALSE)
  fr <- compute_pcoa_frames_aitch(df, metadata, model.vars = model_vars, n_axes = max(axes_to_plot))
  frames_cache_clr[[nm]] <- fr
  
  xcol <- paste0("PCo", axes_to_plot[1]); ycol <- paste0("PCo", axes_to_plot[2])
  scores <- data.frame(
    AX1 = fr$plot.df[[xcol]],
    AX2 = fr$plot.df[[ycol]],
    batch = if (batch_var %in% names(fr$plot.df)) droplevels(fr$plot.df[[batch_var]]) else factor(1)
  )
  if (!is.factor(scores$batch)) scores$batch <- factor(scores$batch)
  xr <- safe_range(scores$AX1); yr <- safe_range(scores$AX2)
  eb <- ellipse_union_bounds(scores, "batch", level = 0.95, n = 240)
  
  global_x_clr <- if (is.null(global_x_clr)) finite_range(xr, eb$x) else finite_range(global_x_clr, xr, eb$x)
  global_y_clr <- if (is.null(global_y_clr)) finite_range(yr, eb$y) else finite_range(global_y_clr, yr, eb$y)
}

pad_x <- safe_pad(global_x_clr, 0.06); pad_y <- safe_pad(global_y_clr, 0.06)
xlim_global_clr <- c(global_x_clr[1] - pad_x, global_x_clr[2] + pad_x)
ylim_global_clr <- c(global_y_clr[1] - pad_y, global_y_clr[2] + pad_y)
if (isTRUE(SYMMETRIC_AXES)) {
  x_half <- max(abs(xlim_global_clr)); y_half <- max(abs(ylim_global_clr))
  xlim_global_clr <- c(-x_half, x_half); ylim_global_clr <- c(-y_half, y_half)
}

plots_clr <- lapply(names(file_list_clr), function(nm) {
  fr <- frames_cache_clr[[nm]]
  pcoa_panel(fr$plot.df, fr$metric.df, model_vars,
             axes = axes_to_plot, label = nm,
             xlim_override = xlim_global_clr, ylim_override = ylim_global_clr,
             palette_name = "Batch")
})
names(plots_clr) <- names(file_list_clr)

# ---- Combine & save (CLR) ----
n_panels_clr <- length(plots_clr)
if (n_panels_clr == 1L) {
  combined_clr <- plots_clr[[1]] +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "vertical",
      plot.margin = margin(8, 14, 8, 14)
    )
  w_clr <- 9.5; h_clr <- 6
} else {
  combined_clr <- wrap_plots(plots_clr, ncol = ncol_grid) +
    plot_layout(guides = "collect") &
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "vertical",
      plot.margin = margin(8, 14, 8, 14)
    )
  w_clr <- 9.5 * min(ncol_grid, n_panels_clr)
  h_clr <- 6   * ceiling(n_panels_clr / ncol_grid)
}
ggsave(file.path(output_folder, "pcoa_aitchison.png"),
       plot = combined_clr, width = w_clr, height = h_clr, dpi = 300)
ggsave(file.path(output_folder, "pcoa_aitchison.tif"),
       plot = combined_clr, width = w_clr, height = h_clr, dpi = 300, compression = "lzw")

# =========================
# Set 2: Bray–Curtis (TSS)
# =========================
message(sprintf("PCoA (Bray–Curtis on TSS): color=%s%s",
                batch_var,
                if (is.na(shape_var)) " (shape: none)" else paste0(", shape=", shape_var)))

frames_cache_tss <- list()
global_x_tss <- NULL; global_y_tss <- NULL

for (nm in names(file_list_tss)) {
  cat("Computing TSS/Bray frames:", nm, "\n")
  df <- read_csv(file_list_tss[[nm]], show_col_types = FALSE)
  fr <- compute_pcoa_frames_bray(df, metadata, model.vars = model_vars, n_axes = max(axes_to_plot))
  frames_cache_tss[[nm]] <- fr
  
  xcol <- paste0("PCo", axes_to_plot[1]); ycol <- paste0("PCo", axes_to_plot[2])
  scores <- data.frame(
    AX1 = fr$plot.df[[xcol]],
    AX2 = fr$plot.df[[ycol]],
    batch = if (batch_var %in% names(fr$plot.df)) droplevels(fr$plot.df[[batch_var]]) else factor(1)
  )
  if (!is.factor(scores$batch)) scores$batch <- factor(scores$batch)
  xr <- safe_range(scores$AX1); yr <- safe_range(scores$AX2)
  eb <- ellipse_union_bounds(scores, "batch", level = 0.95, n = 240)
  
  global_x_tss <- if (is.null(global_x_tss)) finite_range(xr, eb$x) else finite_range(global_x_tss, xr, eb$x)
  global_y_tss <- if (is.null(global_y_tss)) finite_range(yr, eb$y) else finite_range(global_y_tss, yr, eb$y)
}

pad_x <- safe_pad(global_x_tss, 0.06); pad_y <- safe_pad(global_y_tss, 0.06)
xlim_global_tss <- c(global_x_tss[1] - pad_x, global_x_tss[2] + pad_x)
ylim_global_tss <- c(global_y_tss[1] - pad_y, global_y_tss[2] + pad_y)
if (isTRUE(SYMMETRIC_AXES)) {
  x_half <- max(abs(xlim_global_tss)); y_half <- max(abs(ylim_global_tss))
  xlim_global_tss <- c(-x_half, x_half); ylim_global_tss <- c(-y_half, y_half)
}

plots_tss <- lapply(names(file_list_tss), function(nm) {
  fr <- frames_cache_tss[[nm]]
  pcoa_panel(fr$plot.df, fr$metric.df, model_vars,
             axes = axes_to_plot, label = nm,
             xlim_override = xlim_global_tss, ylim_override = ylim_global_tss,
             palette_name = "Batch")
})
names(plots_tss) <- names(file_list_tss)

# ---- Combine & save (TSS) ----
n_panels_tss <- length(plots_tss)
if (n_panels_tss == 1L) {
  combined_tss <- plots_tss[[1]] +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "vertical",
      plot.margin = margin(8, 14, 8, 14)
    )
  w_tss <- 9.5; h_tss <- 6
} else {
  combined_tss <- wrap_plots(plots_tss, ncol = ncol_grid) +
    plot_layout(guides = "collect") &
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "vertical",
      plot.margin = margin(8, 14, 8, 14)
    )
  w_tss <- 9.5 * min(ncol_grid, n_panels_tss)
  h_tss <- 6   * ceiling(n_panels_tss / ncol_grid)
}
ggsave(file.path(output_folder, "pcoa_braycurtis.png"),
       plot = combined_tss, width = w_tss, height = h_tss, dpi = 300)
ggsave(file.path(output_folder, "pcoa_braycurtis.tif"),
       plot = combined_tss, width = w_tss, height = h_tss, dpi = 300, compression = "lzw")

# =========================
# Unified PCoA ranking (Aitchison + Bray combined)
# =========================

compute_centroids_pcoa <- function(scores, batch_var = "batchid") {
  scores %>%
    dplyr::group_by(!!rlang::sym(batch_var)) %>%
    dplyr::summarise(PCo1 = mean(PCo1), PCo2 = mean(PCo2), .groups = "drop")
}
compute_centroid_distances <- function(centroids) {
  if (nrow(centroids) < 2) return(NA_real_)
  as.numeric(mean(dist(centroids[, c("PCo1", "PCo2")], method = "euclidean")))
}
# Weighted mean within-batch pairwise distance on PCo1–PCo2
compute_within_dispersion_pcoa <- function(scores, batch_var = "batchid") {
  if (!all(c("PCo1","PCo2", batch_var) %in% names(scores))) return(NA_real_)
  levs <- levels(scores[[batch_var]])
  if (is.null(levs)) levs <- unique(scores[[batch_var]])
  ws <- c(); ns <- c()
  for (lev in levs) {
    sub <- scores[scores[[batch_var]] == lev, c("PCo1","PCo2"), drop = FALSE]
    n <- nrow(sub)
    if (n < 2) next
    d <- stats::dist(sub, method = "euclidean")
    ws <- c(ws, mean(d))
    ns <- c(ns, n)
  }
  if (!length(ws)) return(NA_real_)
  stats::weighted.mean(ws, w = ns)
}
pcoa_metric_score <- function(batch_distance, coverage) {
  if (is.na(batch_distance) || is.na(coverage)) return(NA_real_)
  coverage <- max(0, min(1, coverage))  # clamp [0,1]
  (1 / (1 + batch_distance)) * coverage  # higher = better
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
      dplyr::select(sample_id, PCo1, PCo2) %>%
      dplyr::mutate(batchid = factor(md$batchid))
    cents <- compute_centroids_pcoa(scores, "batchid")
    D_between <- compute_centroid_distances(cents)
    W_within  <- compute_within_dispersion_pcoa(scores, "batchid")
    ve <- fr$metric.df$var.explained
    coverage2 <- sum(ve[1:min(2, length(ve))], na.rm = TRUE) / 100
    score <- pcoa_metric_score(D_between, coverage2)
    needs_correction <- is.finite(D_between) && is.finite(W_within) && (D_between > W_within)
    
    assess_rows[["CLR"]] <- tibble::tibble(
      Method            = "Before correction",
      Geometry          = "Aitchison (CLR)",
      Batch_Distance    = D_between,
      Within_Dispersion = W_within,
      Coverage_PC1_PC2  = coverage2,
      Score             = score,
      Needs_Correction  = needs_correction
    )
  }
  
  if ("Before correction" %in% methods_tss) {
    fr <- frames_cache_tss[["Before correction"]]
    md <- metadata[match(fr$plot.df$sample_id, metadata$sample_id), , drop = FALSE]
    scores <- fr$plot.df %>%
      dplyr::select(sample_id, PCo1, PCo2) %>%
      dplyr::mutate(batchid = factor(md$batchid))
    cents <- compute_centroids_pcoa(scores, "batchid")
    D_between <- compute_centroid_distances(cents)
    W_within  <- compute_within_dispersion_pcoa(scores, "batchid")
    ve <- fr$metric.df$var.explained
    coverage2 <- sum(ve[1:min(2, length(ve))], na.rm = TRUE) / 100
    score <- pcoa_metric_score(D_between, coverage2)
    needs_correction <- is.finite(D_between) && is.finite(W_within) && (D_between > W_within)
    
    assess_rows[["TSS"]] <- tibble::tibble(
      Method            = "Before correction",
      Geometry          = "Bray–Curtis (TSS)",
      Batch_Distance    = D_between,
      Within_Dispersion = W_within,
      Coverage_PC1_PC2  = coverage2,
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
      Coverage_PC1_PC2  = NA_real_,
      Score             = comb_score,
      Needs_Correction  = needs_corr_any
    )
  )
  
  print(assess_df, n = nrow(assess_df))
  readr::write_csv(assess_df, file.path(output_folder, "pcoa_raw_assessment.csv"))
  
  message(if (isTRUE(needs_corr_any)) {
    "Assessment: Batch separation exceeds within-batch spread in at least one geometry — correction recommended."
  } else {
    "Assessment: Batch separation does not exceed within-batch spread — correction may not be necessary."
  })
  
} else {
  # ===== Multi-method ranking (as before) =====
  rank_tbl <- dplyr::tibble(
    Method = character(),
    Batch_Distance_Aitchison = numeric(),
    Coverage_Aitchison       = numeric(),
    Score_Aitchison          = numeric(),
    Batch_Distance_Bray      = numeric(),
    Coverage_Bray            = numeric(),
    Score_Bray               = numeric(),
    Combined_Score           = numeric()
  )
  
  for (m in union(names(frames_cache_clr), names(frames_cache_tss))) {
    # ----- Aitchison/CLR -----
    D_a <- NA_real_; Cov_a <- NA_real_; S_a <- NA_real_
    if (m %in% names(frames_cache_clr)) {
      fr_a <- frames_cache_clr[[m]]
      md_a <- metadata[match(fr_a$plot.df$sample_id, metadata$sample_id), , drop = FALSE]
      scores_a <- fr_a$plot.df %>%
        dplyr::select(sample_id, PCo1, PCo2) %>%
        dplyr::mutate(batchid = factor(md_a$batchid))
      cents_a <- compute_centroids_pcoa(scores_a, "batchid")
      D_a <- compute_centroid_distances(cents_a)
      ve_a <- fr_a$metric.df$var.explained
      Cov_a <- sum(ve_a[1:min(2, length(ve_a))], na.rm = TRUE) / 100
      S_a <- pcoa_metric_score(D_a, Cov_a)
    }
    
    # ----- Bray–Curtis/TSS -----
    D_b <- NA_real_; Cov_b <- NA_real_; S_b <- NA_real_
    if (m %in% names(frames_cache_tss)) {
      fr_b <- frames_cache_tss[[m]]
      md_b <- metadata[match(fr_b$plot.df$sample_id, metadata$sample_id), , drop = FALSE]
      scores_b <- fr_b$plot.df %>%
        dplyr::select(sample_id, PCo1, PCo2) %>%
        dplyr::mutate(batchid = factor(md_b$batchid))
      cents_b <- compute_centroids_pcoa(scores_b, "batchid")
      D_b <- compute_centroid_distances(cents_b)
      ve_b <- fr_b$metric.df$var.explained
      Cov_b <- sum(ve_b[1:min(2, length(ve_b))], na.rm = TRUE) / 100
      S_b <- pcoa_metric_score(D_b, Cov_b)
    }
    
    Combined <- dplyr::case_when(
      !is.na(S_a) && !is.na(S_b) ~ sqrt(S_a * S_b),
      !is.na(S_a)                ~ S_a,
      !is.na(S_b)                ~ S_b,
      TRUE                       ~ NA_real_
    )
    
    rank_tbl <- dplyr::bind_rows(rank_tbl, dplyr::tibble(
      Method = m,
      Batch_Distance_Aitchison = D_a,
      Coverage_Aitchison       = Cov_a,
      Score_Aitchison          = S_a,
      Batch_Distance_Bray      = D_b,
      Coverage_Bray            = Cov_b,
      Score_Bray               = S_b,
      Combined_Score           = Combined
    ))
  }
  
  ranked_pcoa_unified <- rank_tbl %>%
    dplyr::arrange(dplyr::desc(Combined_Score)) %>%
    dplyr::mutate(Rank = dplyr::row_number())
  
  print(ranked_pcoa_unified, n = nrow(ranked_pcoa_unified))
  readr::write_csv(ranked_pcoa_unified, file.path(output_folder, "pcoa_ranking.csv"))
}
