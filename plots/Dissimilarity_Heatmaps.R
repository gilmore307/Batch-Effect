# ==== Dissimilarity Heatmaps: Aitchison RMSE (CLR) & Bray–Curtis (TSS) ====
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(vegan)  # Bray–Curtis
})

# ==== IO ====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# ==== Metadata ====
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

# Include raw.csv as "Before correction" in both lists if present
raw_fp <- file.path(output_folder, "raw.csv")
if (file.exists(raw_fp)) {
  file_list_clr <- c("Before correction" = raw_fp, file_list_clr)
  file_list_tss <- c("Before correction" = raw_fp, file_list_tss)
}
if (!length(file_list_clr) && !length(file_list_tss)) {
  stop("No normalized files found (expected normalized_*_clr.csv and/or normalized_*_tss.csv) in ", output_folder)
}

batch_var <- "batchid"

# ==== Helpers ====
sort_levels_numeric <- function(x) {
  x <- as.character(x)
  xn <- suppressWarnings(as.numeric(x))
  if (all(!is.na(xn))) as.character(sort(xn)) else sort(x, method = "radix")
}

# TSS closure (rows sum to 1); robust to zeros/NA/negatives
safe_closure <- function(X) {
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0
  rs <- rowSums(X, na.rm = TRUE)
  bad <- which(rs == 0 | !is.finite(rs))
  if (length(bad)) {
    X[bad, ] <- 1 / ncol(X)
    rs <- rowSums(X, na.rm = TRUE)
  }
  sweep(X, 1, rs, "/")
}

# CLR transform with multiplicative replacement for zeros
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

# If negatives present, treat as already CLR and row-center; else do CLR from TSS
to_clr_for_rmse <- function(X) {
  if (any(X < 0, na.rm = TRUE)) {
    sweep(X, 1, rowMeans(X, na.rm = TRUE), "-")
  } else {
    clr_transform(X)
  }
}

# Euclidean distances -> RMSE by dividing by sqrt(p)
euclidean_to_rmse <- function(D_eucl, p) as.matrix(D_eucl) / sqrt(p)

# Average pairwise sample distances into batch×batch matrix
batch_distance_matrix <- function(D_sample, batch_factor, diag_mode = c("zero","mean","NA")) {
  diag_mode <- match.arg(diag_mode)
  M <- as.matrix(D_sample)
  b_levels <- levels(batch_factor)
  B <- length(b_levels)
  Db <- matrix(NA_real_, B, B, dimnames = list(b_levels, b_levels))
  for (i in seq_len(B)) {
    idx_i <- which(batch_factor == b_levels[i])
    for (j in seq_len(B)) {
      idx_j <- which(batch_factor == b_levels[j])
      if (i == j) {
        if (diag_mode == "zero") {
          Db[i, j] <- 0
        } else if (diag_mode == "mean") {
          if (length(idx_i) >= 2) {
            subM <- M[idx_i, idx_i, drop = FALSE]
            Db[i, j] <- mean(subM[upper.tri(subM)], na.rm = TRUE)
          } else Db[i, j] <- 0
        } else {
          Db[i, j] <- NA_real_
        }
      } else {
        Db[i, j] <- mean(M[idx_i, idx_j, drop = FALSE], na.rm = TRUE)
      }
    }
  }
  Db
}

# ---- Build Aitchison RMSE batch×batch matrix (CLR + Euclidean) ----
rmse_batch_matrix_aitchison <- function(df, metadata, batch_var = "batchid", diag_mode = "zero") {
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfm <- inner_join(df, metadata, by = "sample_id")
  
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  keep <- apply(X, 2, function(col) all(is.finite(col)) && sd(col, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) stop("No variable numeric features remain.")
  
  Xclr <- to_clr_for_rmse(X)
  D_eucl <- dist(Xclr, method = "euclidean")
  D_rmse <- euclidean_to_rmse(D_eucl, p = ncol(Xclr))
  
  b_levels <- sort_levels_numeric(unique(dfm[[batch_var]]))
  bfac <- factor(as.character(dfm[[batch_var]]), levels = b_levels)
  Db <- batch_distance_matrix(D_rmse, bfac, diag_mode = diag_mode)
  list(Db = Db, order = b_levels)
}

# ---- Build Bray–Curtis batch×batch matrix (TSS + vegdist 'bray') ----
dissim_batch_matrix_bray <- function(df, metadata, batch_var = "batchid", diag_mode = "zero") {
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfm <- inner_join(df, metadata, by = "sample_id")
  
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  X <- safe_closure(X)  # proportions
  D_bray <- vegan::vegdist(X, method = "bray")
  
  b_levels <- sort_levels_numeric(unique(dfm[[batch_var]]))
  bfac <- factor(as.character(dfm[[batch_var]]), levels = b_levels)
  Db <- batch_distance_matrix(D_bray, bfac, diag_mode = diag_mode)
  list(Db = Db, order = b_levels)
}

# ---- Generic heatmap panel (upper triangle, shared global scale) ----
upper_heatmap_panel <- function(Db, ord, title_label, fill_label,
                                global_min, global_max,
                                label_digits = 3,
                                text_threshold_frac = 0.6) {
  stopifnot(length(ord) == nrow(Db), length(ord) == ncol(Db))
  idx_map <- setNames(seq_along(ord), ord)
  
  long <- as.data.frame(Db) |>
    mutate(batch1 = rownames(Db)) |>
    pivot_longer(-batch1, names_to = "batch2", values_to = "val") |>
    mutate(
      batch1 = factor(batch1, levels = ord),
      batch2 = factor(batch2, levels = ord),
      i = idx_map[as.character(batch1)],
      j = idx_map[as.character(batch2)]
    ) |>
    filter(i < j) |>
    mutate(
      label   = ifelse(is.na(val), "", formatC(val, format = "f", digits = label_digits)),
      txt_col = ifelse(!is.na(val) & val >= (global_min + text_threshold_frac * (global_max - global_min)),
                       "white", "black")
    )
  
  ggplot(long, aes(x = batch2, y = batch1, fill = val)) +
    geom_tile(width = 0.92, height = 0.92) +
    geom_text(aes(label = label, colour = txt_col), size = 3) +
    scale_colour_identity(guide = "none") +
    scale_fill_viridis_c(
      name = fill_label,
      option = "D",
      direction = -1,         # darker = larger dissimilarity
      limits = c(global_min, global_max),
      oob = scales::squish
    ) +
    coord_fixed() +
    labs(title = title_label, x = NULL, y = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.ticks  = element_blank(),
      panel.grid  = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.margin = margin(8, 10, 8, 10)
    ) +
    guides(fill = guide_colourbar(title.position = "top"))
}

# ==== A) Aitchison RMSE heatmaps ====
diag_mode <- "zero"
label_digits <- 2
text_threshold_frac <- 0.6

mat_list_ait <- list()
ord_list_ait <- list()
for (nm in names(file_list_clr)) {
  cat("Computing Aitchison RMSE batch matrix:", nm, "\n")
  df <- read_csv(file_list_clr[[nm]], show_col_types = FALSE)
  comp <- rmse_batch_matrix_aitchison(df, metadata, batch_var = batch_var, diag_mode = diag_mode)
  mat_list_ait[[nm]] <- comp$Db
  ord_list_ait[[nm]] <- comp$order
}
vals_ait <- unlist(lapply(mat_list_ait, function(M) M[upper.tri(M, diag = FALSE)]))
gmin_ait <- ifelse(is.finite(min(vals_ait, na.rm = TRUE)), min(vals_ait, na.rm = TRUE), 0)
gmax_ait <- ifelse(is.finite(max(vals_ait, na.rm = TRUE)), max(vals_ait, na.rm = TRUE), gmin_ait + 1e-8)

plots_ait <- list()
for (nm in names(mat_list_ait)) {
  plots_ait[[nm]] <- upper_heatmap_panel(
    Db = mat_list_ait[[nm]],
    ord = ord_list_ait[[nm]],
    title_label = paste0("Dissimilarity Heatmap — RMSE (Aitchison/CLR) — ", nm),
    fill_label = "RMSE (Aitchison/CLR)",
    global_min = gmin_ait,
    global_max = gmax_ait,
    label_digits = label_digits,
    text_threshold_frac = text_threshold_frac
  )
}

ncol_grid <- 2
combined_ait <- wrap_plots(plots_ait, ncol = ncol_grid) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave(file.path(output_folder, "dissimilarity_heatmaps_aitchison.png"),
       plot = combined_ait, width = 8.5 * ncol_grid, height = 6 * ceiling(length(plots_ait) / ncol_grid),
       dpi = 300)
ggsave(file.path(output_folder, "dissimilarity_heatmaps_aitchison.tif"),
       plot = combined_ait, width = 8.5 * ncol_grid, height = 6 * ceiling(length(plots_ait) / ncol_grid),
       dpi = 300, compression = "lzw")

# Helper: mean of upper-triangle (off-diagonal) entries
upper_mean <- function(M) {
  M <- as.matrix(M)
  ut <- upper.tri(M, diag = FALSE)
  vals <- M[ut]
  if (!length(vals)) return(NA_real_)
  mean(vals, na.rm = TRUE)
}

# ==== B) Bray–Curtis heatmaps ====
mat_list_bc <- list()
ord_list_bc <- list()
for (nm in names(file_list_tss)) {
  cat("Computing Bray–Curtis batch matrix:", nm, "\n")
  df <- read_csv(file_list_tss[[nm]], show_col_types = FALSE)
  comp <- dissim_batch_matrix_bray(df, metadata, batch_var = batch_var, diag_mode = diag_mode)
  mat_list_bc[[nm]] <- comp$Db
  ord_list_bc[[nm]] <- comp$order
}
vals_bc <- unlist(lapply(mat_list_bc, function(M) M[upper.tri(M, diag = FALSE)]))
gmin_bc <- ifelse(is.finite(min(vals_bc, na.rm = TRUE)), min(vals_bc, na.rm = TRUE), 0)
gmax_bc <- ifelse(is.finite(max(vals_bc, na.rm = TRUE)), max(vals_bc, na.rm = TRUE), gmin_bc + 1e-8)

plots_bc <- list()
for (nm in names(mat_list_bc)) {
  plots_bc[[nm]] <- upper_heatmap_panel(
    Db = mat_list_bc[[nm]],
    ord = ord_list_bc[[nm]],
    title_label = paste0("Dissimilarity Heatmap — Bray–Curtis (TSS) — ", nm),
    fill_label = "Bray–Curtis dissimilarity",
    global_min = gmin_bc,
    global_max = gmax_bc,
    label_digits = label_digits,
    text_threshold_frac = text_threshold_frac
  )
}

combined_bc <- wrap_plots(plots_bc, ncol = ncol_grid) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave(file.path(output_folder, "dissimilarity_heatmaps_braycurtis.png"),
       plot = combined_bc, width = 8.5 * ncol_grid, height = 6 * ceiling(length(plots_bc) / ncol_grid),
       dpi = 300)
ggsave(file.path(output_folder, "dissimilarity_heatmaps_braycurtis.tif"),
       plot = combined_bc, width = 8.5 * ncol_grid, height = 6 * ceiling(length(plots_bc) / ncol_grid),
       dpi = 300, compression = "lzw")

# ==== Unified ranking (Aitchison RMSE + Bray–Curtis) ====
mean_ait <- if (length(mat_list_ait)) sapply(mat_list_ait, upper_mean) else numeric()
mean_bc  <- if (length(mat_list_bc))  sapply(mat_list_bc,  upper_mean) else numeric()

all_methods <- sort(unique(c(names(mean_ait), names(mean_bc))))

# weights for combining (geometric mean)
weights <- c(aitchison = 0.5, bray = 0.5)
wa <- weights["aitchison"]; wb <- weights["bray"]
if (is.na(wa)) wa <- 0.5
if (is.na(wb)) wb <- 0.5
wsum <- wa + wb; wa <- wa / wsum; wb <- wb / wsum

unified_rows <- lapply(all_methods, function(m) {
  m_ait <- if (m %in% names(mean_ait)) mean_ait[[m]] else NA_real_
  m_bc  <- if (m %in% names(mean_bc))  mean_bc[[m]]  else NA_real_
  # Convert “lower is better” into scores (higher is better)
  S_ait <- if (is.na(m_ait)) NA_real_ else 1 / (1 + m_ait)
  S_bc  <- if (is.na(m_bc))  NA_real_ else 1 / (1 + m_bc)
  # Weighted geometric mean
  S_comb <- if (!is.na(S_ait) && !is.na(S_bc)) {
    (S_ait ^ wa) * (S_bc ^ wb)
  } else if (!is.na(S_ait)) {
    S_ait
  } else if (!is.na(S_bc)) {
    S_bc
  } else {
    NA_real_
  }
  data.frame(
    Method                   = m,
    MeanUpper_Aitchison_RMSE = m_ait,
    MeanUpper_Bray           = m_bc,
    Score_Aitchison          = S_ait,
    Score_Bray               = S_bc,
    Combined_Score           = S_comb,
    stringsAsFactors = FALSE
  )
})

ranking_unified <- dplyr::bind_rows(unified_rows) %>%
  dplyr::filter(!is.na(Combined_Score)) %>%
  dplyr::arrange(dplyr::desc(Combined_Score)) %>%
  dplyr::mutate(Rank = dplyr::row_number())

print(as.data.frame(ranking_unified), row.names = FALSE)
readr::write_csv(ranking_unified, file.path(output_folder, "dissimilarity_ranking.csv"))
