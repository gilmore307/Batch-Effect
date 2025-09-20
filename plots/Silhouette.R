# ===================== Silhouette (CLR-only, UMAP-based) =====================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(uwot)
  library(cluster)   # silhouette
})

# ---- Config ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

UMAP_NEIGHB   <- 15
UMAP_MIN_DIST <- 0.3
UMAP_METRIC   <- "euclidean"  # CLR -> Euclidean
LABEL_COL     <- "phenotype"  # label to evaluate silhouette on (fallbacks below)

set.seed(42)

# ---- Metadata ----
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))

# Try a sensible fallback label if 'phenotype' is missing
if (!(LABEL_COL %in% names(metadata))) {
  fallback <- c("celltype","label","group","condition","status","class")
  LABEL_COL <- fallback[fallback %in% names(metadata)][1] %||% stop(
    "No label column found (looked for: phenotype, celltype, label, group, condition, status, class)."
  )
}

# ---- Files (CLR only, plus raw.csv as 'Before Correction') ----
clr_paths <- list.files(output_folder, pattern = "^normalized_.*_clr\\.csv$", full.names = TRUE)
raw_fp <- file.path(output_folder, "raw.csv")
if (file.exists(raw_fp)) clr_paths <- c(raw_fp, clr_paths)
if (!length(clr_paths)) stop("No CLR files found (expected 'normalized_*_clr.csv' or raw.csv).")

name_fun <- function(x) gsub("^normalized_|_clr\\.csv$", "", basename(x))
file_list <- setNames(
  clr_paths,
  ifelse(basename(clr_paths) == "raw.csv", "Before Correction", name_fun(clr_paths))
)

# ---- Helpers (reuse from your EBM script if present) ----
`%||%` <- function(a,b) if (is.null(a)) b else a

safe_numeric_matrix <- function(df) {
  num <- dplyr::select(df, where(is.numeric))
  if (!ncol(num)) return(matrix(numeric(0), nrow = nrow(df)))
  keep <- vapply(num, function(z) all(is.finite(z)) && sd(z) > 0, logical(1))
  as.matrix(num[, keep, drop = FALSE])
}
safe_closure <- function(X) {
  X <- as.matrix(X); X[!is.finite(X)] <- 0; X[X < 0] <- 0
  rs <- rowSums(X); rs[rs == 0] <- 1
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

# Silhouette on UMAP; returns mean in [0,1] via rescaling (1 + S)/2
silhouette_on_umap <- function(umap_coords, labels_factor) {
  labs <- droplevels(as.factor(labels_factor))
  if (nlevels(labs) < 2) return(NA_real_)
  d <- dist(as.matrix(umap_coords), method = "euclidean")
  sil <- tryCatch(cluster::silhouette(as.integer(labs), d), error = function(e) NULL)
  if (is.null(sil)) return(NA_real_)
  S <- mean(sil[, "sil_width"], na.rm = TRUE)  # [-1, 1]
  (1 + S) / 2                                  # -> [0, 1]
}

# ---- Compute Silhouette per method ----
sil_tbl <- tibble(Method = character(), Silhouette = numeric())

for (nm in names(file_list)) {
  cat("Processing (Silhouette):", nm, "\n")
  df <- read_csv(file_list[[nm]], show_col_types = FALSE)
  
  if (!("sample_id" %in% names(df))) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else { warning(sprintf("Skipping %s: no sample_id and row mismatch.", nm)); next }
  }
  df <- df |> mutate(sample_id = as.character(sample_id))
  merged <- inner_join(df, metadata, by = "sample_id")
  if (!nrow(merged)) { warning(sprintf("Skipping %s: no overlap with metadata.", nm)); next }
  
  X <- safe_numeric_matrix(merged[, setdiff(names(df), "sample_id"), drop = FALSE])
  if (nrow(X) < 3 || ncol(X) < 2) { warning(sprintf("Skipping %s: insufficient features.", nm)); next }
  
  # Ensure CLR; if values already look CLR (negatives), just row-recenter
  has_neg <- any(X < 0, na.rm = TRUE)
  Xclr <- if (has_neg) sweep(X, 1, rowMeans(X, na.rm = TRUE), "-") else clr_transform(X, 1)
  
  um <- tryCatch(
    get_umap2d(Xclr, n_neighbors = UMAP_NEIGHB, min_dist = UMAP_MIN_DIST, metric = UMAP_METRIC),
    error = function(e) NULL
  )
  if (is.null(um)) { warning(sprintf("Skipping %s: UMAP failed.", nm)); next }
  
  labs <- merged[[LABEL_COL]]
  sil  <- tryCatch(silhouette_on_umap(um, labs), error = function(e) NA_real_)
  sil_tbl <- bind_rows(sil_tbl, tibble(Method = nm, Silhouette = sil))
}

# ---- Rank, save, plot (AS-style) ----
pretty_metric <- function(m) {
  m <- tolower(m)
  if (m == "euclidean") return("Euclidean (Aitchison)")
  if (m == "cosine")    return("Cosine")
  mbecUpper <- paste0(toupper(substr(m,1,1)), substr(m,2,nchar(m)))
  return(mbecUpper)
}

sil_ranked <- sil_tbl %>%
  filter(is.finite(Silhouette)) %>%
  arrange(desc(Silhouette), Method) %>%
  mutate(Rank = row_number())

readr::write_csv(sil_ranked, file.path(output_folder, "silhouette_ranking.csv"))
print(sil_ranked, n = nrow(sil_ranked))

p_sil <- ggplot(sil_ranked, aes(x = reorder(Method, Silhouette), y = Silhouette, fill = Method)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.4, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", Silhouette)), vjust = -0.4, size = 3.2) +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Silhouette Score (UMAP)",
    subtitle = sprintf("Labels=%s • UMAP(metric=%s, n_neighbors=%d, min_dist=%.2f)",
                       LABEL_COL, pretty_metric(UMAP_METRIC), UMAP_NEIGHB, UMAP_MIN_DIST),
    x = "Method", y = "Silhouette (0–1, higher = tighter class separation)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave(file.path(output_folder, "silhouette.png"), p_sil, width = 8.5, height = 5.2, dpi = 300)
ggsave(file.path(output_folder, "silhouette.tif"), p_sil, width = 8.5, height = 5.2, dpi = 300, compression = "lzw")
