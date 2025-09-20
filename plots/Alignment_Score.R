# ================= Alignment Score (AS) — CLR only, with ranking =================
suppressPackageStartupMessages({
  library(FNN)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# --------- Args / config ---------
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")  # change/remove for CLI
if (length(args) < 1) stop("Usage: Rscript alignment_score_clr.R <output_folder>")
output_folder <- args[1]

K_NEIGHBORS   <- 10     # k in kNN
VAR_PROP_MIN  <- 0.95   # keep PCs explaining at least 95% variance
MAX_PCS       <- 10     # safety cap
set.seed(42)

# --------- Load metadata ---------
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))
if (!("batchid" %in% names(metadata)))
  stop("metadata.csv must contain a 'batchid' column.")

# --------- Collect CLR files ---------
clr_paths <- list.files(output_folder, pattern = "^normalized_.*_clr\\.csv$", full.names = TRUE)

# include raw.csv (as baseline) if present
raw_fp <- file.path(output_folder, "raw.csv")
if (file.exists(raw_fp)) clr_paths <- c(raw_fp, clr_paths)

if (!length(clr_paths)) stop("No CLR matrices found (expected 'normalized_*_clr.csv').")

method_names <- ifelse(basename(clr_paths) == "raw.csv",
                       "Before correction",
                       gsub("^normalized_|_clr\\.csv$", "", basename(clr_paths)))
file_list <- setNames(clr_paths, method_names)

# --------- Helpers ---------
safe_numeric_matrix <- function(df) {
  num <- dplyr::select(df, where(is.numeric))
  if (!ncol(num)) return(matrix(numeric(0), nrow = nrow(df)))
  keep <- vapply(num, function(z) all(is.finite(z)) && sd(z) > 0, logical(1))
  as.matrix(num[, keep, drop = FALSE])
}

compute_alignment_score <- function(X, batch, k = 10, var_prop = 0.95, max_pcs = 10) {
  # X: samples x features (numeric), batch: factor/character vector
  if (!is.matrix(X) || nrow(X) < 3 || ncol(X) < 2) return(NA_real_)
  pca <- prcomp(X, scale. = TRUE)
  var_cum <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  npc <- min(max(which(var_cum < var_prop)) + 1, ncol(pca$x))
  npc <- max(2, min(npc, max_pcs))
  coords <- pca$x[, seq_len(npc), drop = FALSE]
  k_eff <- max(1, min(k, nrow(coords) - 1))
  nn <- FNN::get.knn(coords, k = k_eff)$nn.index
  batch <- as.character(batch)
  si <- vapply(seq_len(nrow(coords)), function(i) {
    nb <- batch[nn[i, ]]
    1 - mean(nb == batch[i], na.rm = TRUE)  # fraction of neighbors from other batches
  }, numeric(1))
  mean(si, na.rm = TRUE)
}

# NEW: simple ranking function (higher AS = better)
rank_alignment_methods <- function(as_table) {
  as_table %>%
    filter(is.finite(AS)) %>%
    arrange(desc(AS), Method) %>%
    mutate(Rank = row_number())
}

# --------- Compute AS per method ---------
as_tbl <- tibble(Method = character(), AS = numeric())

for (nm in names(file_list)) {
  fp <- file_list[[nm]]
  df <- read_csv(fp, show_col_types = FALSE)
  if (!("sample_id" %in% names(df))) {
    if (nrow(df) == nrow(metadata)) {
      df$sample_id <- metadata$sample_id
    } else {
      warning(sprintf("Skipping %s: no sample_id and row count mismatch.", nm)); next
    }
  }
  df <- df |> mutate(sample_id = as.character(sample_id))
  merged <- inner_join(df, metadata, by = "sample_id")
  if (!nrow(merged)) { warning(sprintf("Skipping %s: no overlap with metadata.", nm)); next }
  X <- safe_numeric_matrix(merged[, setdiff(names(df), "sample_id"), drop = FALSE])
  b <- merged$batchid
  ascore <- tryCatch(
    compute_alignment_score(X, b, k = K_NEIGHBORS, var_prop = VAR_PROP_MIN, max_pcs = MAX_PCS),
    error = function(e) NA_real_
  )
  as_tbl <- bind_rows(as_tbl, tibble(Method = nm, AS = ascore))
}

# --------- Rank, save, and print ---------
as_ranked <- rank_alignment_methods(as_tbl)
readr::write_csv(as_ranked, file.path(output_folder, "alignment_score_ranking.csv"))
print(as_ranked)

# --------- Plot ---------
p_as <- ggplot(as_ranked, aes(x = reorder(Method, AS), y = AS, fill = Method)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.4, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", AS)), vjust = -0.4, size = 3.2) +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Alignment Score",
       x = "Method", y = "AS (0–1, higher = better mixing)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave(file.path(output_folder, "alignment_score.png"), p_as, width = 8.5, height = 5.2, dpi = 300)
ggsave(file.path(output_folder, "alignment_score.tif"), p_as, width = 8.5, height = 5.2, dpi = 300, compression = "lzw")
