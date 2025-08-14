# ==== Libraries ====
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(patchwork)  # for layouts + legend collection
  library(rlang)
})

# ==== Helpers ====
guess_shape_var <- function(meta, batch_col = "batchid") {
  cand <- c("group","phenotype","condition","status","case_control","class","disease","label")
  hit <- cand[cand %in% names(meta)]
  if (length(hit)) return(hit[1])
  facs <- names(Filter(function(v) {
    !is.numeric(meta[[v]]) && !is.integer(meta[[v]]) && length(unique(meta[[v]])) <= 12
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
  list(x = c(min(xmins, na.rm = TRUE), max(xmaxs, na.rm = TRUE)),
       y = c(min(ymins, na.rm = TRUE), max(ymaxs, na.rm = TRUE)))
}

# --- Aitchison (CLR) helpers ---
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

# ==== Read args / IO ====
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")
if (length(args) < 1) stop("Usage: Rscript pcoa.R <output_folder>")
output_folder <- args[1]

metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))

file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list  <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list  <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== PCoA frames (scores + variance) ====
compute_pcoa_frames <- function(df, metadata, model.vars = c("batchid","group"),
                                n_axes = 5, dist_method = c("auto","bray","aitchison","euclidean")) {
  dist_method <- match.arg(dist_method)
  
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfm <- inner_join(df, metadata, by = "sample_id")
  
  # numeric feature matrix
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feat_cols)) %>% select(where(is.numeric))
  keep <- vapply(X, function(x) sd(x, na.rm = TRUE) > 0, logical(1))
  X <- as.matrix(X[, keep, drop = FALSE])
  if (!ncol(X)) stop("No variable numeric features remain for PCoA.")
  
  has_neg <- any(X < 0, na.rm = TRUE)
  bray_ok <- (!has_neg) && requireNamespace("vegan", quietly = TRUE)
  
  # ---- choose distance: Bray if possible, else Aitchison (CLR + Euclidean) ----
  method <- dist_method
  if (dist_method == "auto") {
    method <- if (bray_ok) "bray" else "aitchison"
  } else if (dist_method == "bray" && !bray_ok) {
    method <- "aitchison"
  }
  
  if (method == "bray") {
    d <- vegan::vegdist(X, method = "bray")
  } else if (method == "aitchison") {
    if (has_neg) {
      Xclr <- sweep(X, 1, rowMeans(X, na.rm = TRUE), "-") # treat as already log-ratio
    } else {
      Xclr <- clr_transform(X)
    }
    d <- dist(Xclr, method = "euclidean")
  } else {
    d <- dist(X, method = "euclidean")
  }
  
  # ---- PCoA ----
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
  metric.df <- data.frame(
    axis.min      = as.numeric(axis.min),
    axis.max      = as.numeric(axis.max),
    var.explained = as.numeric(var_expl),
    stringsAsFactors = FALSE
  )
  
  plot.df <- data.frame(sample_id = dfm$sample_id, as.data.frame(sc), check.names = FALSE)
  colnames(plot.df)[2:(k + 1)] <- paste0("PCo", seq_len(k))
  
  # attach covariates with consistent factor levels from metadata
  present <- intersect(model.vars, names(dfm))
  for (v in present) {
    levs <- unique(as.character(metadata[[v]]))
    plot.df[[v]] <- factor(as.character(dfm[[v]]), levels = levs)
  }
  
  list(plot.df = plot.df, metric.df = metric.df, used.vars = present)
}

# ==== PCoA panel (scatter + marginal densities; legend bottom, horizontal) ====
pcoa_panel <- function(plot.df, metric.df, model.vars, axes = c(1,2), label = NULL) {
  mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                "#AEC7E8","#FFBB78","#98DF8A","#D62728","#FF9896","#C5B0D5",
                "#8C564B","#C49C94","#F7B6D2","#7F7F7F","#C7C7C7","#DBDB8D",
                "#17BECF","#9EDAE5")
  
  var.color <- model.vars[1]
  var.shape <- ifelse(length(model.vars) >= 2, model.vars[2], NA_character_)
  
  xcol <- paste0("PCo", axes[1])
  ycol <- paste0("PCo", axes[2])
  x.label <- paste0(xcol, ": ", metric.df$var.explained[axes[1]], "% expl.var")
  y.label <- paste0(ycol, ": ", metric.df$var.explained[axes[2]], "% expl.var")
  if (!is.null(label)) x.label <- paste(label, "-", x.label)
  
  # zoomed-out limits based on union of ellipses + points
  scores <- data.frame(
    AX1 = plot.df[[xcol]],
    AX2 = plot.df[[ycol]],
    batch = if (var.color %in% names(plot.df)) droplevels(plot.df[[var.color]]) else factor(1)
  )
  ell_bounds <- ellipse_union_bounds(scores, "batch", level = 0.95, n = 240)
  xr <- range(scores$AX1, na.rm = TRUE); yr <- range(scores$AX2, na.rm = TRUE)
  xlim <- range(c(xr, ell_bounds$x)); ylim <- range(c(yr, ell_bounds$y))
  pad_x <- diff(xlim) * 0.06; pad_y <- diff(ylim) * 0.06
  xlim <- c(xlim[1] - pad_x, xlim[2] + pad_x)
  ylim <- c(ylim[1] - pad_y, ylim[2] + pad_y)
  
  pmar <- margin(10, 16, 10, 16)
  
  # main scatter
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
      colour = guide_legend(order = 1, nrow = 1, byrow = TRUE),  # batches on one long row
      shape  = guide_legend(order = 2, nrow = 1)                 # phenotype under it (after collect)
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
  
  # top density (Axis 1)
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
  
  # right density (Axis 2), rotated so density extends RIGHT; same vertical scale as main
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
      plot.margin = pmar
    )
  
  # assemble panel (legend will be collected globally)
  design <- "
A#
CB
"
  (pTop + pRight + pMain) + plot_layout(design = design, widths = c(3, 1), heights = c(1.6, 3.2))
}

# ==== Params ====
dist_method <- "auto"  # "auto" | "bray" | "aitchison" | "euclidean"
batch_var   <- "batchid"
shape_var   <- guess_shape_var(metadata, batch_var)
model_vars  <- if (is.na(shape_var)) c(batch_var) else c(batch_var, shape_var)
message(sprintf("PCoA settings: distance=%s, color=%s%s",
                dist_method, batch_var,
                if (is.na(shape_var)) " (shape: none)" else paste0(", shape=", shape_var)))

# ==== Build + save ====
axes_to_plot <- c(1, 2)
plots <- list()

for (nm in names(file_list)) {
  cat("Processing (PCoA):", nm, "\n")
  df <- read_csv(file_list[[nm]], show_col_types = FALSE)
  
  frames <- compute_pcoa_frames(df, metadata, model.vars = model_vars,
                                n_axes = max(axes_to_plot), dist_method = dist_method)
  plt <- pcoa_panel(frames$plot.df, frames$metric.df, frames$used.vars,
                    axes = axes_to_plot, label = nm)
  
  plots[[nm]] <- plt
}

# Combine ALL panels and keep ONE legend at bottom (horizontal)
ncol_grid <- 2
combined <- wrap_plots(plots, ncol = ncol_grid) +
  plot_layout(guides = "collect") &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.box       = "vertical",
    plot.margin = margin(8, 14, 8, 14)
  )

ggsave(file.path(output_folder, "pcoa.png"),
       plot = combined, width = 9.5 * ncol_grid, height = 6 * ceiling(length(plots) / ncol_grid),
       dpi = 300)
ggsave(file.path(output_folder, "pcoa.tif"),
       plot = combined, width = 9.5 * ncol_grid, height = 6 * ceiling(length(plots) / ncol_grid),
       dpi = 300, compression = "lzw")


# Function to compute the centroid of each batch in the PCoA space
compute_centroids <- function(pcoa_scores, batch_var) {
  # Compute centroids (mean of each batch group in PCoA space)
  centroids <- pcoa_scores %>%
    group_by(!!sym(batch_var)) %>%
    summarise(PCo1 = mean(PCo1), PCo2 = mean(PCo2))
  
  return(centroids)
}

# Function to compute the distance between centroids of different batches
compute_centroid_distances <- function(centroids) {
  # Calculate pairwise Euclidean distances between centroids
  dist_matrix <- dist(centroids[, c("PCo1", "PCo2")], method = "euclidean")
  return(mean(dist_matrix))  # Average distance between batch centroids
}

# Initialize a data frame to store results
rank_results <- tibble(Method = character(), Batch_Distance = numeric())

# Loop over all methods to compute the distance between batch centroids
for (method_name in names(file_list)) {
  cat("Processing:", method_name, "\n")
  
  # Read in the data for the current method
  df <- read_csv(file_list[[method_name]], show_col_types = FALSE)
  
  # Compute PCoA on the data
  frames <- compute_pcoa_frames(df, metadata, model.vars = c("batchid"),
                                n_axes = 2, dist_method = "auto")
  
  # Get PCoA scores and batch information
  pcoa_scores <- frames$plot.df %>% select(sample_id, PCo1, PCo2, batchid)
  
  # Compute centroids for each batch
  centroids <- compute_centroids(pcoa_scores, batch_var = "batchid")
  
  # Compute the distance between batch centroids
  batch_distance <- compute_centroid_distances(centroids)
  
  # Store the result
  rank_results <- bind_rows(rank_results, tibble(Method = method_name, Batch_Distance = batch_distance))
}

# Rank the methods based on batch separation (distance between centroids)
ranked_methods <- rank_results %>%
  arrange(desc(Batch_Distance)) %>%
  mutate(Rank = row_number())

# Display the ranked methods based on batch separation
ranked_methods
