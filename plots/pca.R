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
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")
if (length(args) < 1) stop("Usage: Rscript pca.R <output_folder>")
output_folder <- args[1]

# ==== Read Metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))

# ==== Find All Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list  <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list  <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

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
    # Put Batch first, long horizontal row; Phenotype under it (legend.box='vertical')
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
  
  # right density (PC2) — use coord_flip(), set limits with scale_y_continuous
  # ---- right density (PC2) — rotate so PC2 is vertical, density extends right ----
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
    # match the main plot's y-range exactly
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
plots <- list()

for (nm in names(file_list)) {
  cat("Processing:", nm, "\n")
  df <- read_csv(file_list[[nm]], show_col_types = FALSE)
  
  frames    <- compute_pca_frames(df, metadata, model.vars = model_vars, n_pcs = max(pcs_to_plot))
  used_vars <- frames$used.vars
  
  plt <- mbecPCAPlot(
    plot.df   = frames$plot.df,
    metric.df = frames$metric.df,
    model.vars = used_vars,
    pca.axes  = pcs_to_plot,
    label     = nm
  )
  
  plots[[nm]] <- plt

}

# ---- Combine ALL panels and keep ONLY ONE legend at the bottom (horizontal) ----
ncol_grid <- 2
combined <- wrap_plots(plots, ncol = ncol_grid) +
  plot_layout(guides = "collect") &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.box       = "vertical",
    plot.margin = margin(8, 14, 8, 14)
  )

ggsave(file.path(output_folder, "pca.png"),
       plot = combined, width = 9.5 * ncol_grid, height = 6 * ceiling(length(plots) / ncol_grid),
       dpi = 300)
ggsave(file.path(output_folder, "pca.tif"),
       plot = combined, width = 9.5 * ncol_grid, height = 6 * ceiling(length(plots) / ncol_grid),
       dpi = 300, compression = "lzw")


# Function to compute the centroid of each batch in PCA space
compute_centroids <- function(pca_scores, batch_var) {
  # Compute centroids (mean of each batch group in PCA space)
  centroids <- pca_scores %>%
    group_by(!!sym(batch_var)) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2))
  
  return(centroids)
}

# Function to compute the distance between centroids of different batches
compute_centroid_distances <- function(centroids) {
  # Calculate pairwise Euclidean distances between centroids
  dist_matrix <- dist(centroids[, c("PC1", "PC2")], method = "euclidean")
  return(mean(dist_matrix))  # Average distance between batch centroids
}

# Initialize a data frame to store results
rank_results <- tibble(Method = character(), Batch_Distance = numeric())

# Loop over all methods to compute the distance between batch centroids
for (method_name in names(file_list)) {
  cat("Processing:", method_name, "\n")
  
  # Read in the data for the current method
  df <- read_csv(file_list[[method_name]], show_col_types = FALSE)
  
  # Compute PCA on the data
  frames <- compute_pca_frames(df, metadata, model.vars = c("batchid"), n_pcs = 2)
  
  # Get PCA scores and batch information
  pca_scores <- frames$plot.df %>% select(sample_id, PC1, PC2, batchid)
  
  # Compute centroids for each batch
  centroids <- compute_centroids(pca_scores, batch_var = "batchid")
  
  # Compute the distance between batch centroids
  batch_distance <- compute_centroid_distances(centroids)
  
  # Store the result
  rank_results <- bind_rows(rank_results, tibble(Method = method_name, Batch_Distance = batch_distance))
}

# Rank the methods based on batch separation
ranked_methods <- rank_results %>%
  arrange(desc(Batch_Distance)) %>%
  mutate(Rank = row_number())

# Display the ranked methods based on batch separation
ranked_methods
