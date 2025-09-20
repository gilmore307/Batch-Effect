#!/usr/bin/env Rscript

# ======================= LISI (PCA-only) end-to-end script =======================
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(scales)

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

  library(FNN)        # k-NN
  library(purrr)
  library(patchwork)  # combine plots
})

has_repel <- requireNamespace("ggrepel", quietly = TRUE)

# ------------------------------- CLI options ------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# defaults (can be overridden via flags):
opt_k       <- 30          # local neighborhood size
opt_npcs    <- 50          # number of PCs (set to NA to use CLR directly)
opt_coords  <- "pca"       # "pca" or "clr" (alias for n_pcs = NA)
opt_kauto   <- FALSE       # auto-select k by trade-off plateau

# parse flags: --k=30 --npcs=50 --coords=pca|clr --kauto
for (a in args) {
  if (grepl("^--k=", a))       opt_k      <- as.integer(sub("^--k=", "", a))
  if (grepl("^--npcs=", a))    opt_npcs   <- as.integer(sub("^--npcs=", "", a))
  if (grepl("^--coords=", a))  opt_coords <- tolower(sub("^--coords=", "", a))
  if (a == "--kauto")          opt_kauto  <- TRUE
}
if (!is.na(opt_npcs) && opt_coords == "clr") opt_npcs <- NA           # honor --coords=clr
if (is.na(opt_npcs) && opt_coords == "pca")  opt_npcs <- 50            # honor --coords=pca

message("Output folder: ", output_folder)
message("k (neighbors): ", opt_k, if (opt_kauto) "  [auto if enabled]" else "")
message("Coordinates: ", if (is.na(opt_npcs)) "CLR (no PCA)" else paste0("PCA (", opt_npcs, " PCs)"))

# ------------------------------- Helpers ----------------------------------------
safe_closure <- function(X) {
  # Row-wise closure to simplex with NA/Inf/zero handling
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
  # Per-row positive floor to avoid log(0) and handle non-finite values
  for (i in seq_len(nrow(Xc))) {
    xi <- Xc[i, ]
    pos <- xi > 0 & is.finite(xi)
    if (!any(pos)) { xi[] <- 1/length(xi); pos <- xi > 0 }
    if (any(!pos)) {
      m <- min(xi[pos], na.rm = TRUE)
      xi[!pos] <- min(m * 0.5, 1e-8)
      xi <- xi / sum(xi)
    }
    Xc[i, ] <- xi
  }
  L <- log(Xc)
  sweep(L, 1, rowMeans(L), "-")
}

# Per-sample LISI and normalized LISI in a k-NN graph
# norm = (LISI - 1) / (L - 1), where L is #unique labels
compute_lisi_normalized <- function(coords, labels, k) {
  n <- nrow(coords)
  if (n < 2) return(tibble(lisi = rep(NA_real_, n), norm = rep(NA_real_, n)))
  k <- max(1, min(k, n - 1))
  kn <- FNN::get.knn(coords, k = k)$nn.index
  L <- length(unique(labels))
  if (L < 2) return(tibble(lisi = rep(1, n), norm = rep(0, n)))
  
  lisi <- numeric(n)
  for (i in seq_len(n)) {
    neigh <- labels[kn[i, ]]
    tab <- table(neigh)
    p <- as.numeric(tab) / k
    lisi[i] <- 1 / sum(p^2)
  }
  norm <- (lisi - 1) / (L - 1)
  tibble(lisi = lisi, norm = pmin(pmax(norm, 0), 1))
}

get_coords <- function(Xclr, n_pcs = 50) {
  if (is.na(n_pcs)) return(Xclr)  # CLR space
  pr <- prcomp(Xclr, center = FALSE, scale. = FALSE)
  k <- min(n_pcs, ncol(pr$x))
  pr$x[, seq_len(k), drop = FALSE]
}

# Compute iLISI (batch) + cLISI (phenotype) in chosen geometry for one matrix
compute_lisi_for_method <- function(df, meta, method_name,
                                    batch_col = "batch_id", treat_col = "phenotype",
                                    n_pcs = 50, k_nn = 30) {
  # Align on sample_id
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(meta) && "sample_id" %in% names(meta)) df$sample_id <- meta$sample_id
    else stop("Input lacks 'sample_id' and does not align with metadata.")
  }
  df   <- df   %>% mutate(sample_id = as.character(sample_id))
  meta <- meta %>% mutate(sample_id = as.character(sample_id))
  dfx  <- inner_join(df, meta, by = "sample_id")
  if (!nrow(dfx)) stop("No overlapping samples between data and metadata.")
  if (!(batch_col %in% names(dfx)) || !(treat_col %in% names(dfx)))
    stop("Metadata must contain batch and treatment columns.")
  
  # Numeric features
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfx %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  keep <- apply(X, 2, function(z) all(is.finite(z)) && sd(z, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) stop("No valid numeric features after filtering.")
  
  # CLR -> coords (PCA or CLR direct)
  Xclr       <- clr_transform(X)
  coords     <- get_coords(Xclr, n_pcs = n_pcs)
  n          <- nrow(coords)
  k          <- max(1, min(k_nn, n - 1))
  
  if (k >= n - 1 || k > 0.5 * n) {
    warning(sprintf("k=%d is large relative to n=%d; LISI may approximate global proportions.", k, n))
  }
  
  # Labels
  batch <- factor(dfx[[batch_col]])
  pheno <- factor(dfx[[treat_col]])
  
  # LISI in chosen geometry
  ilisi <- compute_lisi_normalized(coords, batch, k)
  clisi <- compute_lisi_normalized(coords, pheno, k)
  
  bind_rows(
    tibble(Method = method_name, Metric = "iLISI", Value = ilisi$norm, Raw = ilisi$lisi),
    tibble(Method = method_name, Metric = "cLISI", Value = clisi$norm, Raw = clisi$lisi)
  )
}

# Global-proportion baseline (for the "too-global k" check)
global_lisi_norm <- function(lbl) {
  tab <- table(lbl); p <- as.numeric(tab) / sum(tab)
  LISI <- 1 / sum(p^2)
  L <- length(tab)
  (LISI - 1) / (L - 1)
}

# -------- Your ranking function (unchanged): higher iLISI, lower cLISI ----------
rank_lisi_methods <- function(summary_df) {
  stopifnot(all(c("Method","median_iLISI","median_cLISI") %in% names(summary_df)))
  summary_df %>%
    transmute(
      Method,
      iLISI = median_iLISI,
      cLISI = median_cLISI,
      r_i = rank(-iLISI, ties.method = "average"),  # higher iLISI is better
      r_c = rank( cLISI, ties.method = "average"),  # lower  cLISI is better
      Score = (r_i + r_c) / 2
    ) %>%
    arrange(Score, Method) %>%
    mutate(Rank = dplyr::row_number()) %>%
    select(Method, iLISI, cLISI, Score, Rank)
}

# ---- Optional: Auto-select k by trade-off + plateau (from earlier message) -----
select_k_auto <- function(file_list, metadata,
                          ks = c(10, 15, 20, 30, 40, 50, 60),
                          n_pcs = 50,
                          batch_col = "batch_id",
                          treat_col = "phenotype",
                          plateau_frac = 0.98,
                          slope_tol = 0.005) {
  stopifnot(length(file_list) > 0)
  grid <- purrr::map_dfr(ks, function(kv) {
    lisi_k <- purrr::map_dfr(names(file_list), function(nm) {
      df <- readr::read_csv(file_list[[nm]], show_col_types = FALSE)
      compute_lisi_for_method(
        df, metadata, method_name = nm,
        batch_col = batch_col, treat_col = treat_col,
        n_pcs = n_pcs, k_nn = kv
      )
    })
    summ <- lisi_k |>
      dplyr::group_by(Method, Metric) |>
      dplyr::summarise(median = median(Value, na.rm = TRUE), .groups = "drop") |>
      tidyr::pivot_wider(names_from = Metric, values_from = median) |>
      dplyr::mutate(separation = 1 - cLISI)
    
    tibble::tibble(
      k = kv,
      mean_iLISI = mean(summ$iLISI, na.rm = TRUE),
      mean_sep   = mean(summ$separation, na.rm = TRUE)
    )
  }) |>
    dplyr::mutate(
      score = ifelse(mean_iLISI + mean_sep > 0,
                     2 * mean_iLISI * mean_sep / (mean_iLISI + mean_sep), 0)
    ) |>
    dplyr::arrange(k) |>
    dplyr::mutate(
      dscore = c(NA, diff(score) / diff(k))  # discrete slope
    )
  
  smax <- max(grid$score, na.rm = TRUE)
  cand <- grid |>
    dplyr::filter(score >= plateau_frac * smax) |>
    dplyr::mutate(ok = is.na(dscore) | abs(dscore) <= slope_tol) |>
    dplyr::filter(ok)
  
  k_star <- if (nrow(cand)) min(cand$k) else grid$k[which.max(grid$score)]
  list(k = k_star, grid = grid)
}

# ------------------------------------ IO ----------------------------------------
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

metadata <- readr::read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) %>%
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

# Optionally auto-select k
if (opt_kauto) {
  ks_scan <- c(10, 15, 20, 30, 40, 50, 60)
  auto <- select_k_auto(file_list, metadata,
                        ks = ks_scan,
                        n_pcs = if (is.na(opt_npcs)) NA_integer_ else opt_npcs,
                        batch_col = "batch_id",
                        treat_col = "phenotype")
  opt_k <- auto$k
  readr::write_csv(auto$grid, file.path(output_folder, "LISI_k_grid.csv"))
  message(sprintf("Auto-selected k=%d (harmonic-mean trade-off, plateau rule).", opt_k))
}

# -------------------------- Compute LISI per method -----------------------------
lisi_long <- lapply(names(file_list), function(nm) {
  message("Computing LISI for: ", nm)
  df <- readr::read_csv(file_list[[nm]], show_col_types = FALSE)
  compute_lisi_for_method(
    df, metadata, method_name = nm,
    batch_col = "batch_id", treat_col = "phenotype",
    n_pcs = if (is.na(opt_npcs)) NA_integer_ else opt_npcs,
    k_nn  = opt_k
  )
}) %>% bind_rows() %>%
  mutate(Method = factor(Method, levels = names(file_list)),
         Metric = factor(Metric, levels = c("iLISI","cLISI")))

# -------------------------- Summaries + Ranking ---------------------------------
summary_df <- lisi_long %>%
  group_by(Method, Metric) %>%
  summarise(
    median = median(Value, na.rm = TRUE),
    se     = sd(Value, na.rm = TRUE) / sqrt(sum(is.finite(Value))),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Metric, values_from = c(median, se))

ranked <- rank_lisi_methods(summary_df)   # <-- your ranking function
print(ranked)
write_csv(ranked, file.path(output_folder, "LISI_ranking.csv"))

# ------------------------------------ Plots -------------------------------------
# 1) Scatter: iLISI (y) vs 1-cLISI (x)
plot_df <- summary_df %>%
  transmute(
    Method,
    x = 1 - median_cLISI,               # 1 - cLISI on x (higher = better separation)
    y =      median_iLISI,              # iLISI on y (higher = better mixing)
    ylab = sprintf("%.4f", y)
  )

p_scatter <- ggplot(plot_df, aes(x = x, y = y, color = Method)) +
  geom_point(size = 3) +
  {
    if (has_repel) {
      ggrepel::geom_text_repel(aes(label = ylab),
                               size = 3, seed = 1,
                               box.padding = 0.2, point.padding = 0.2,
                               show.legend = FALSE, max.overlaps = Inf)
    } else {
      geom_text(aes(label = ylab), nudge_y = 0.03, size = 3, show.legend = FALSE)
    }
  } +
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.01)) +
  labs(
    x = expression("1 - cLISI (phenotype separation)"),
    y = "iLISI (batch mixing)",
    title = paste0("LISI (", if (is.na(opt_npcs)) "CLR" else paste0(opt_npcs, " PCs"),
                   ", k=", opt_k, ")")
  ) +
  theme_bw() +
  theme(
    legend.position    = "right",
    legend.title       = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

# 2) Boxplots (distributions)
p_i <- ggplot(lisi_long %>% filter(Metric=="iLISI"),
              aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.9) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "iLISI distribution", y = "iLISI", x = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

p_c <- ggplot(lisi_long %>% filter(Metric=="cLISI"),
              aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.9) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "cLISI distribution", y = "cLISI", x = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

# 3) Combined panel: scatter on top, two boxplots below
combined <- p_scatter /
  (p_i | p_c) +
  plot_layout(heights = c(1.6, 1), guides = "collect") &
  theme(legend.position = "right")  # collect legend to the right

ggsave(file.path(output_folder, "LISI.png"),
       combined, width = 9.0, height = 8.0, dpi = 300)
ggsave(file.path(output_folder, "LISI.tif"),
       combined, width = 9.0, height = 8.0, dpi = 300, compression = "lzw")
