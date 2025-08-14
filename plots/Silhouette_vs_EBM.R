# ==== Libraries ====
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(cluster)   # silhouette
  library(uwot)      # UMAP
  library(hexbin)    # hex-binning for entropy
})

# ==== Args / config ====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

PHENO_COL     <- "phenotype"  # numeric 0/1 or 2-class factor in metadata.csv
UMAP_NEIGHB   <- 15           # UMAP n_neighbors
UMAP_MIN_DIST <- 0.3          # UMAP min_dist
UMAP_METRIC   <- "euclidean"  # UMAP metric
UMAP_BINS     <- 30           # hexbin xbins for entropy

set.seed(42)                  # uwot::umap uses global RNG

# ==== Read metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

if (!PHENO_COL %in% names(metadata))
  stop(sprintf("Phenotype column '%s' not found in metadata.csv", PHENO_COL))

# Prepare .outcome as a two-class factor; positive class first (for consistency)
if (is.numeric(metadata[[PHENO_COL]]) && length(unique(metadata[[PHENO_COL]])) == 2) {
  metadata <- metadata %>%
    mutate(.outcome = factor(ifelse(.data[[PHENO_COL]] == 1, "pos", "neg"),
                             levels = c("pos","neg")))
} else {
  levs <- levels(factor(metadata[[PHENO_COL]]))
  if (length(levs) != 2) stop(sprintf("Phenotype column '%s' must have exactly 2 classes.", PHENO_COL))
  pos <- levs[2]
  metadata <- metadata %>%
    mutate(.outcome = relevel(factor(.data[[PHENO_COL]]), ref = pos))
}

# ==== Collect files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== Utils ====
`%||%` <- function(a,b) if (is.null(a)) b else a

safe_numeric_matrix <- function(df) {
  num <- dplyr::select(df, where(is.numeric))
  if (ncol(num) == 0) return(matrix(numeric(0), nrow = nrow(df)))
  num <- num[, colSums(is.na(num)) == 0, drop = FALSE]  # drop any-NA cols
  if (ncol(num) > 0) {
    sds <- apply(num, 2, sd)
    num <- num[, sds > 0, drop = FALSE]                 # drop zero-variance
  }
  as.matrix(num)
}

# --- Shared UMAP helper (2D embedding on scaled features) ---
get_umap2d <- function(X, n_neighbors = 15, min_dist = 0.3, metric = "euclidean") {
  if (!is.matrix(X) || nrow(X) < 3 || ncol(X) < 2) return(NULL)
  Xs <- scale(X)
  emb <- uwot::umap(
    Xs,
    n_neighbors = n_neighbors,
    min_dist    = min_dist,
    metric      = metric,
    n_components = 2,
    verbose     = FALSE
  )
  as.matrix(emb)
}

# --- Silhouette of phenotype on UMAP (mean over samples), rescaled to [0,1] ---
silhouette_on_umap <- function(X, labels, n_neighbors = 15, min_dist = 0.3, metric = "euclidean") {
  um <- get_umap2d(X, n_neighbors, min_dist, metric)
  if (is.null(um)) return(NA_real_)
  labs <- droplevels(as.factor(labels))
  if (nlevels(labs) < 2) return(NA_real_)
  d <- dist(um, method = "euclidean")
  sil <- tryCatch(cluster::silhouette(as.integer(labs), d), error = function(e) NULL)
  if (is.null(sil)) return(NA_real_)
  S <- mean(sil[, "sil_width"], na.rm = TRUE)  # in [-1, 1]
  (1 + S) / 2                                  # rescale to [0,1]
}

# --- Entropy of batch mixing on UMAP hex-bins (normalized to [0,1]) ---
# Fix: drop zero-probability components to avoid 0*log(0) -> NaN
batch_entropy_umap <- function(X, batch, nbins = 30,
                               n_neighbors = 15, min_dist = 0.3, metric = "euclidean") {
  um <- get_umap2d(X, n_neighbors, min_dist, metric)
  if (is.null(um)) return(NA_real_)
  
  batch <- droplevels(as.factor(batch))  # keep only existing batch levels overall
  B <- nlevels(batch)
  if (B <= 1) return(0)
  
  hb <- hexbin::hexbin(um[,1], um[,2], xbins = nbins, IDs = TRUE)
  cell_id <- hb@cID  # vector length = nrow(um), each entry is a bin id
  
  # per-cell normalized Shannon entropy (drop p==0)
  H <- tapply(seq_along(cell_id), cell_id, function(idx) {
    cnt <- table(droplevels(batch[idx]))        # only levels present in this bin
    p <- as.numeric(cnt) / sum(cnt)
    p <- p[p > 0]                                # critical: remove zeros
    Hraw <- if (length(p)) -sum(p * log(p)) else 0
    Hraw / log(B)                                # normalize to [0,1] using global B
  })
  
  sizes <- tapply(seq_along(cell_id), cell_id, length)
  if (length(H) == 0 || length(sizes) == 0) return(NA_real_)
  sum(unlist(H) * unlist(sizes), na.rm = TRUE) / sum(unlist(sizes), na.rm = TRUE)
}

# ==== Compute metrics per method ====
results <- tibble(
  Method = character(),
  Silhouette = numeric(),   # phenotype silhouette on UMAP, [0,1]
  Entropy = numeric()       # batch-mixing entropy on UMAP, [0,1]
)

for (name in names(file_list)) {
  fpath <- file_list[[name]]
  if (!file.exists(fpath)) next
  df_raw <- read_csv(fpath, show_col_types = FALSE)
  
  # Ensure sample_id present
  if (!"sample_id" %in% names(df_raw)) {
    if (nrow(df_raw) == nrow(metadata)) {
      df_raw <- df_raw %>% mutate(sample_id = metadata$sample_id)
    } else {
      warning(sprintf("File %s has no sample_id and row count doesn't match metadata; skipping.", fpath))
      next
    }
  }
  
  df_merged <- df_raw %>%
    mutate(sample_id = as.character(sample_id)) %>%
    inner_join(metadata, by = "sample_id")
  
  feature_cols <- setdiff(names(df_raw), "sample_id")
  X <- safe_numeric_matrix(df_merged[, feature_cols, drop = FALSE])
  
  y     <- df_merged$.outcome
  batch <- df_merged$batchid %||% rep(1, nrow(df_merged))
  batch <- as.factor(batch)
  
  sil <- tryCatch(
    silhouette_on_umap(X, y, n_neighbors = UMAP_NEIGHB, min_dist = UMAP_MIN_DIST, metric = UMAP_METRIC),
    error = function(e) NA_real_
  )
  ent <- tryCatch(
    batch_entropy_umap(X, batch, nbins = UMAP_BINS, n_neighbors = UMAP_NEIGHB, min_dist = UMAP_MIN_DIST, metric = UMAP_METRIC),
    error = function(e) NA_real_
  )
  
  results <- bind_rows(results, tibble(
    Method = name, Silhouette = sil, Entropy = ent
  ))
}

# ==== Plot: Silhouette (biology) vs Entropy (batch mixing) ====
results_plot <- results %>%
  filter(is.finite(Silhouette), is.finite(Entropy)) %>%
  mutate(
    Silhouette = pmin(pmax(Silhouette, 0), 1),
    Entropy    = pmin(pmax(Entropy,    0), 1)
  )

p_se <- ggplot(results_plot, aes(x = Entropy, y = Silhouette, label = Method, color = Method)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = Inf) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Phenotype Silhouette (UMAP) vs. Entropy of Batch Mixing (UMAP hex-bins)",
       x = "Entropy of Batch Mixing (normalized, 0–1)",
       y = "Phenotype Silhouette (rescaled, 0–1)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(output_folder, "silhouette_vs_entropy.png"),
       p_se, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_folder, "silhouette_vs_entropy.tif"),
       p_se, width = 8, height = 6, dpi = 300)

# Function to rank the methods based on Silhouette and Entropy
rank_methods_silhouette_entropy <- function(results) {
  # Rank methods based on Silhouette (higher is better) and Entropy (lower is better)
  ranked_results <- results %>%
    arrange(desc(Silhouette), Entropy) %>%  # Prioritize Silhouette (higher) and then Entropy (lower)
    mutate(Rank = row_number())  # Assign ranks
  
  return(ranked_results)
}

# Rank the methods based on Silhouette and Entropy values
ranked_methods <- rank_methods_silhouette_entropy(results)

# Display the ranked methods
ranked_methods
