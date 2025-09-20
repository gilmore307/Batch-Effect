# ===================== EBM (CLR-only, UMAP+knn) =====================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(uwot)     # UMAP
  library(FNN)      # kNN
})

# --------- Args / config ---------
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")  # change/remove for CLI
if (length(args) < 1) stop("Usage: Rscript ebm_clr_only.R <output_folder>")
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# UMAP (CLR/Aitchison geometry)
UMAP_NEIGHB   <- 15
UMAP_MIN_DIST <- 0.3
UMAP_METRIC   <- "euclidean"   # CLR -> Euclidean

# EBM (kNN) params — "his method"
KNN_K         <- 50
KNN_POOLS     <- 50
KNN_PER_LABEL <- 100

# EBM labels column (batch mixing uses batches)
LABEL_COL <- "batchid"

set.seed(42)  # uwot uses global RNG

# --------- Load metadata ---------
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))
if (!(LABEL_COL %in% names(metadata))) {
  stop(sprintf("metadata.csv must contain a '%s' column.", LABEL_COL))
}

# --------- Collect CLR files (matches your PCA loader) ---------
clr_paths <- list.files(output_folder, pattern = "^normalized_.*_clr\\.csv$", full.names = TRUE)
raw_fp <- file.path(output_folder, "raw.csv")
if (file.exists(raw_fp)) clr_paths <- c(raw_fp, clr_paths)
if (!length(clr_paths)) stop("No CLR files found (expected 'normalized_*_clr.csv' or raw.csv).")

name_fun <- function(x) gsub("^normalized_|_clr\\.csv$", "", basename(x))
file_list <- setNames(
  clr_paths,
  ifelse(basename(clr_paths) == "raw.csv", "Before Correction", name_fun(clr_paths))
)

# --------- Helpers ---------
`%||%` <- function(a,b) if (is.null(a)) b else a

safe_numeric_matrix <- function(df) {
  num <- dplyr::select(df, where(is.numeric))
  if (!ncol(num)) return(matrix(numeric(0), nrow = nrow(df)))
  keep <- vapply(num, function(z) all(is.finite(z)) && sd(z) > 0, logical(1))
  as.matrix(num[, keep, drop = FALSE])
}

# CLR transform (used for raw.csv or any non-CLR input)
safe_closure <- function(X) {
  X <- as.matrix(X)
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0
  rs <- rowSums(X)
  rs[rs == 0] <- 1
  sweep(X, 1, rs, "/")
}
clr_transform <- function(X, pseudocount = 1) {
  Xp <- safe_closure(X) + (pseudocount / max(1, ncol(X)))
  L  <- log(Xp)
  sweep(L, 1, rowMeans(L), "-")
}

get_umap2d <- function(X, n_neighbors = 15, min_dist = 0.3, metric = "euclidean") {
  if (!is.matrix(X) || nrow(X) < 3 || ncol(X) < 2) return(NULL)
  Xs <- scale(X)
  umap(Xs, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric,
       n_components = 2, verbose = FALSE)
}

# ---- kNN-based EBM (higher = better) ----
normalized_entropy <- function(p, B) {
  p <- p[p > 0]
  if (!length(p)) return(0)
  -sum(p * log(p)) / log(B)
}
batch_entropy_mixing_knn <- function(
    X2d, batch, k = 50, n_pools = 50, n_per_label = 100, seed = NULL, return_per_anchor = FALSE
) {
  stopifnot(is.matrix(X2d) || is.data.frame(X2d))
  X2d <- as.matrix(X2d); n <- nrow(X2d)
  if (length(batch) != n) stop("length(batch) must equal nrow(X2d)")
  batch <- droplevels(factor(batch)); B <- nlevels(batch)
  if (B <= 1) {
    return(list(mean_entropy = 0, EBM_Score = 0,
                pool_means = rep(0, max(1, n_pools)),
                per_anchor = if (return_per_anchor) numeric(0) else NULL))
  }
  k_eff <- min(max(1, k), n - 1)
  nn <- FNN::get.knn(X2d, k = k_eff)$nn.index
  idx_by_batch <- split(seq_len(n), batch)
  per_pool_means <- numeric(n_pools)
  all_anchor_ent <- if (return_per_anchor) numeric(0) else NULL
  if (!is.null(seed)) set.seed(seed)
  
  for (p in seq_len(n_pools)) {
    anchors <- unlist(lapply(idx_by_batch, function(ix) {
      if (length(ix) >= n_per_label) sample(ix, n_per_label, replace = FALSE)
      else                           sample(ix, n_per_label, replace = TRUE)
    }), use.names = FALSE)
    
    ent <- vapply(anchors, function(i) {
      nbrs <- nn[i, ]
      cnt  <- table(batch[nbrs])
      pvec <- rep(0, B); names(pvec) <- levels(batch)
      pvec[names(cnt)] <- as.numeric(cnt) / sum(cnt)
      normalized_entropy(pvec, B)
    }, numeric(1))
    
    per_pool_means[p] <- mean(ent, na.rm = TRUE)
    if (return_per_anchor) all_anchor_ent <- c(all_anchor_ent, ent)
  }
  
  mean_ebm <- mean(per_pool_means, na.rm = TRUE)
  res <- list(mean_entropy = mean_ebm, EBM_Score = mean_ebm, pool_means = per_pool_means)
  if (return_per_anchor) res$per_anchor <- all_anchor_ent
  res
}

rank_methods_ebm <- function(ebm_table) {
  ebm_table %>% filter(is.finite(EBM)) %>% arrange(desc(EBM), Method) %>% mutate(Rank = row_number())
}

pretty_metric <- function(m) {
  m <- tolower(m)
  if (m == "euclidean") return("Euclidean (Aitchison)")
  if (m == "cosine")    return("Cosine")
  mbecUpper <- paste0(toupper(substr(m,1,1)), substr(m,2,nchar(m)))
  return(mbecUpper)
}

# --------- Compute EBM per CLR matrix ---------
ebm_tbl <- tibble(Method = character(), EBM = numeric(), CI95_lo = numeric(), CI95_hi = numeric())

for (nm in names(file_list)) {
  cat("Processing:", nm, "\n")
  df <- read_csv(file_list[[nm]], show_col_types = FALSE)
  
  # Ensure sample_id present
  if (!("sample_id" %in% names(df))) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else { warning(sprintf("Skipping %s: no sample_id and row mismatch.", nm)); next }
  }
  df <- df |> mutate(sample_id = as.character(sample_id))
  merged <- inner_join(df, metadata, by = "sample_id")
  if (!nrow(merged)) { warning(sprintf("Skipping %s: no overlap with metadata.", nm)); next }
  
  # Features (CLR matrix or raw)
  X <- safe_numeric_matrix(merged[, setdiff(names(df), "sample_id"), drop = FALSE])
  if (nrow(X) < 3 || ncol(X) < 2) { warning(sprintf("Skipping %s: insufficient features.", nm)); next }
  b <- as.factor(merged[[LABEL_COL]])
  
  # If negatives exist, treat as already-CLR; else CLR-transform
  has_neg <- any(X < 0, na.rm = TRUE)
  Xclr <- if (has_neg) {
    # re-center rows (CLR-like input)
    sweep(X, 1, rowMeans(X, na.rm = TRUE), "-")
  } else {
    clr_transform(X, pseudocount = 1)
  }
  
  um <- tryCatch(
    get_umap2d(Xclr, n_neighbors = UMAP_NEIGHB, min_dist = UMAP_MIN_DIST, metric = UMAP_METRIC),
    error = function(e) NULL
  )
  if (is.null(um)) { warning(sprintf("Skipping %s: UMAP failed.", nm)); next }
  
  out <- tryCatch(
    batch_entropy_mixing_knn(
      um, b, k = KNN_K, n_pools = KNN_POOLS, n_per_label = KNN_PER_LABEL, seed = 42
    ),
    error = function(e) NULL
  )
  if (is.null(out)) { warning(sprintf("Skipping %s: EBM failed.", nm)); next }
  
  m  <- out$mean_entropy
  pm <- out$pool_means
  n  <- length(pm)
  se <- if (n > 1) stats::sd(pm) / sqrt(n) else 0
  ci_lo <- max(0, m - 1.96 * se)
  ci_hi <- min(1, m + 1.96 * se)
  
  ebm_tbl <- bind_rows(ebm_tbl, tibble(Method = nm, EBM = m, CI95_lo = ci_lo, CI95_hi = ci_hi))
}

# --------- Rank, save ---------
ebm_ranked <- rank_methods_ebm(ebm_tbl)
readr::write_csv(ebm_ranked, file.path(output_folder, "ebm_ranking.csv"))
print(ebm_ranked, n = nrow(ebm_ranked))

# --------- Plot (AS-style) with subtitle IN the ggplot block ---------
p_ebm <- ggplot(ebm_ranked, aes(x = reorder(Method, EBM), y = EBM, fill = Method)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.4, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", EBM)), vjust = -0.4, size = 3.2) +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Batch Entropy Mixing (kNN on UMAP)",
    subtitle = sprintf(
      "Labels=%s • UMAP(metric=%s, n_neighbors=%d, min_dist=%.2f)",
      LABEL_COL, pretty_metric(UMAP_METRIC), UMAP_NEIGHB, UMAP_MIN_DIST
    ),
    x = "Method",
    y = "EBM (0–1, higher = better mixing)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave(file.path(output_folder, "ebm.png"), p_ebm, width = 8.5, height = 5.2, dpi = 300)
ggsave(file.path(output_folder, "ebm.tif"), p_ebm, width = 8.5, height = 5.2, dpi = 300, compression = "lzw")
