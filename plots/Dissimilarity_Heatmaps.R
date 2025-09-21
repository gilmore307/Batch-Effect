# =========================
# Dissimilarity Heatmaps: Aitchison RMSE (CLR) & Bray-Curtis (TSS)
# =========================
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(vegan)  # Bray-Curtis
  
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

# Include raw_clr.csv / raw_tss.csv as "Before correction" if present
raw_clr_fp <- file.path(output_folder, "raw_clr.csv")
raw_tss_fp <- file.path(output_folder, "raw_tss.csv")
if (file.exists(raw_clr_fp)) file_list_clr <- c("Before correction" = raw_clr_fp, file_list_clr)
if (file.exists(raw_tss_fp)) file_list_tss <- c("Before correction" = raw_tss_fp, file_list_tss)

if (!length(file_list_clr) && !length(file_list_tss)) {
  stop("No normalized files found (expected raw_clr.csv/raw_tss.csv and/or normalized_*_clr.csv / normalized_*_tss.csv) in ", output_folder)
}

batch_var <- "batch_id"

# ==== Helpers ====
sort_levels_numeric <- function(x) {
  x <- as.character(x)
  xn <- suppressWarnings(as.numeric(x))
  if (all(!is.na(xn))) as.character(sort(xn)) else sort(x, method = "radix")
}

# Ensure a 1:1 alignment between data frame rows and metadata by
# pairing duplicate sample_ids in their order of appearance.
align_samples_1to1 <- function(df, metadata) {
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df <- df %>% mutate(sample_id = as.character(sample_id))
  md <- metadata %>% mutate(sample_id = as.character(sample_id))

  df <- df %>% group_by(sample_id) %>% mutate(.dup_idx = dplyr::row_number()) %>% ungroup()
  md <- md %>% group_by(sample_id) %>% mutate(.dup_idx = dplyr::row_number()) %>% ungroup()

  # Join by sample_id + within-id duplicate index to avoid many-to-many expansion
  out <- suppressWarnings(dplyr::inner_join(df, md, by = c("sample_id", ".dup_idx")))
  out <- dplyr::select(out, -".dup_idx")
  out
}

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

clr_transform <- function(X) {
  Xc <- safe_closure(X)
  for (i in seq_len(nrow(Xc))) {
    xi <- Xc[i, ]
    pos <- xi > 0 & is.finite(xi)
    if (!any(pos)) { xi[] <- 1 / length(xi); pos <- xi > 0 }
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

to_clr_for_rmse <- function(X) {
  if (any(X < 0, na.rm = TRUE)) {
    sweep(X, 1, rowMeans(X, na.rm = TRUE), "-")
  } else {
    clr_transform(X)
  }
}

euclidean_to_rmse <- function(D_eucl, p) as.matrix(D_eucl) / sqrt(p)

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

rmse_batch_matrix_aitchison <- function(df, metadata, batch_var = "batch_id", diag_mode = "zero") {
  dfm <- align_samples_1to1(df, metadata)
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

dissim_batch_matrix_bray <- function(df, metadata, batch_var = "batch_id", diag_mode = "zero") {
  dfm <- align_samples_1to1(df, metadata)
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  X <- safe_closure(X)
  D_bray <- vegan::vegdist(X, method = "bray")
  b_levels <- sort_levels_numeric(unique(dfm[[batch_var]]))
  bfac <- factor(as.character(dfm[[batch_var]]), levels = b_levels)
  Db <- batch_distance_matrix(D_bray, bfac, diag_mode = diag_mode)
  list(Db = Db, order = b_levels)
}

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
      direction = -1,
      limits = c(global_min, global_max),
      oob = scales::squish
    ) +
    # Keep full set of batch levels on both axes
    scale_x_discrete(limits = ord, drop = FALSE) +
    scale_y_discrete(limits = ord, drop = FALSE) +
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

upper_mean <- function(M) {
  M <- as.matrix(M)
  ut <- upper.tri(M, diag = FALSE)
  vals <- M[ut]
  if (!length(vals)) return(NA_real_)
  mean(vals, na.rm = TRUE)
}

# ==== A) Aitchison RMSE heatmaps ====
diag_mode <- "zero"
label_digits <- 2
text_threshold_frac <- 0.6
ncol_grid <- 2

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
    title_label = paste0("Dissimilarity Heatmap - RMSE (Aitchison/CLR) - ", nm),
    fill_label = "RMSE (Aitchison/CLR)",
    global_min = gmin_ait,
    global_max = gmax_ait,
    label_digits = label_digits,
    text_threshold_frac = text_threshold_frac
  )
}

# ---- Combine & save (Aitchison) ----
n_panels_ait <- length(plots_ait)
if (n_panels_ait == 1L) {
  combined_ait <- plots_ait[[1]] +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  w_ait <- 8.5; h_ait <- 6
} else {
  combined_ait <- wrap_plots(plots_ait, ncol = ncol_grid) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal")
  w_ait <- 8.5 * min(ncol_grid, n_panels_ait)
  h_ait <- 6   * ceiling(n_panels_ait / ncol_grid)
}
ggsave(file.path(output_folder, "dissimilarity_heatmaps_aitchison.png"),
       plot = combined_ait, width = w_ait, height = h_ait, dpi = 300)
ggsave(file.path(output_folder, "dissimilarity_heatmaps_aitchison.tif"),
       plot = combined_ait, width = w_ait, height = h_ait, dpi = 300, compression = "lzw")

# ==== B) Bray-Curtis heatmaps ====
mat_list_bc <- list()
ord_list_bc <- list()
for (nm in names(file_list_tss)) {
  cat("Computing Bray-Curtis batch matrix:", nm, "\n")
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
    title_label = paste0("Dissimilarity Heatmap - RMSR (Bray-Curtis/TSS) - ", nm),
    fill_label = "RMSR (Bray-Curtis/TSS)",
    global_min = gmin_bc,
    global_max = gmax_bc,
    label_digits = label_digits,
    text_threshold_frac = text_threshold_frac
  )
}

# ---- Combine & save (Bray-Curtis) ----
n_panels_bc <- length(plots_bc)
if (n_panels_bc == 1L) {
  combined_bc <- plots_bc[[1]] +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  w_bc <- 8.5; h_bc <- 6
} else {
  combined_bc <- wrap_plots(plots_bc, ncol = ncol_grid) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal")
  w_bc <- 8.5 * min(ncol_grid, n_panels_bc)
  h_bc <- 6   * ceiling(n_panels_bc / ncol_grid)
}
ggsave(file.path(output_folder, "dissimilarity_heatmaps_braycurtis.png"),
       plot = combined_bc, width = w_bc, height = h_bc, dpi = 300)
ggsave(file.path(output_folder, "dissimilarity_heatmaps_braycurtis.tif"),
       plot = combined_bc, width = w_bc, height = h_bc, dpi = 300, compression = "lzw")

# ==== Unified ranking (Aitchison RMSE + Bray-Curtis) OR baseline-only assessment ====
mean_ait <- if (length(mat_list_ait)) sapply(mat_list_ait, upper_mean) else numeric()
mean_bc  <- if (length(mat_list_bc))  sapply(mat_list_bc,  upper_mean) else numeric()

all_methods <- sort(unique(c(names(mean_ait), names(mean_bc))))
only_baseline <- (length(all_methods) == 1L) && identical(all_methods, "Before correction")

if (only_baseline) {
  assess_rows <- list()
  
  # Aitchison baseline assessment
  if ("Before correction" %in% names(file_list_clr)) {
    df_raw_clr <- read_csv(file_list_clr[["Before correction"]], show_col_types = FALSE)
    comp_zero  <- rmse_batch_matrix_aitchison(df_raw_clr, metadata, batch_var = batch_var, diag_mode = "zero")
    comp_mean  <- rmse_batch_matrix_aitchison(df_raw_clr, metadata, batch_var = batch_var, diag_mode = "mean")
    between    <- upper_mean(comp_zero$Db)
    within     <- mean(diag(comp_mean$Db), na.rm = TRUE)
    score      <- 1 / (1 + between)
    needs_corr <- is.finite(between) && is.finite(within) && (between > within)
    assess_rows[["CLR"]] <- tibble::tibble(
      Method            = "Before correction",
      Geometry          = "Aitchison (CLR)",
      Mean_Between      = between,
      Mean_Within       = within,
      Score             = score,
      Needs_Correction  = needs_corr
    )
  }
  
  # Bray-Curtis baseline assessment
  if ("Before correction" %in% names(file_list_tss)) {
    df_raw_tss <- read_csv(file_list_tss[["Before correction"]], show_col_types = FALSE)
    comp_zero  <- dissim_batch_matrix_bray(df_raw_tss, metadata, batch_var = batch_var, diag_mode = "zero")
    comp_mean  <- dissim_batch_matrix_bray(df_raw_tss, metadata, batch_var = batch_var, diag_mode = "mean")
    between    <- upper_mean(comp_zero$Db)
    within     <- mean(diag(comp_mean$Db), na.rm = TRUE)
    score      <- 1 / (1 + between)
    needs_corr <- is.finite(between) && is.finite(within) && (between > within)
    assess_rows[["TSS"]] <- tibble::tibble(
      Method            = "Before correction",
      Geometry          = "Bray-Curtis (TSS)",
      Mean_Between      = between,
      Mean_Within       = within,
      Score             = score,
      Needs_Correction  = needs_corr
    )
  }
  
  assess_df <- dplyr::bind_rows(assess_rows)
  
  # Combined summary (geometric mean of available scores; OR on correction flags)
  comb_score <- dplyr::case_when(
    nrow(assess_df) >= 2 && all(is.finite(assess_df$Score)) ~ sqrt(prod(assess_df$Score)),
    TRUE                                                    ~ max(assess_df$Score, na.rm = TRUE)
  )
  needs_corr_any <- any(assess_df$Needs_Correction, na.rm = TRUE)
  
  assess_df <- dplyr::bind_rows(
    assess_df,
    tibble::tibble(
      Method           = "Before correction",
      Geometry         = "Combined",
      Mean_Between     = NA_real_,
      Mean_Within      = NA_real_,
      Score            = comb_score,
      Needs_Correction = needs_corr_any
    )
  )
  
  print(assess_df, n = nrow(assess_df))
  readr::write_csv(assess_df, file.path(output_folder, "dissimilarity_raw_assessment.csv"))
  
  message(if (isTRUE(needs_corr_any)) {
    "Assessment: Between-batch dissimilarity exceeds within-batch spread in at least one geometry - correction recommended."
  } else {
    "Assessment: Between-batch dissimilarity does not exceed within-batch spread - correction may not be necessary."
  })
  
} else {
  # ---- Ranking (multiple methods) ----
  weights <- c(aitchison = 0.5, bray = 0.5)
  wa <- weights["aitchison"]; wb <- weights["bray"]
  if (is.na(wa)) wa <- 0.5
  if (is.na(wb)) wb <- 0.5
  wsum <- wa + wb; wa <- wa / wsum; wb <- wb / wsum
  
  unified_rows <- lapply(all_methods, function(m) {
    m_ait <- if (m %in% names(mean_ait)) mean_ait[[m]] else NA_real_
    m_bc  <- if (m %in% names(mean_bc))  mean_bc[[m]]  else NA_real_
    S_ait <- if (is.na(m_ait)) NA_real_ else 1 / (1 + m_ait)
    S_bc  <- if (is.na(m_bc))  NA_real_ else 1 / (1 + m_bc)
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
}
