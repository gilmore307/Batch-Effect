# ================= Alignment Score (AS) - CLR only, with ranking =================
suppressPackageStartupMessages({
  library(FNN)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# --------- Args / config ---------
args <- commandArgs(trailingOnly = TRUE)

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

if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

K_NEIGHBORS   <- 10     # k in kNN
VAR_PROP_MIN  <- 0.95   # keep PCs explaining at least 95% variance
MAX_PCS       <- 10     # safety cap
AS_THRESHOLD  <- 0.60   # if baseline AS < 0.60, recommend correction
set.seed(42)

# --------- Load metadata ---------
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
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
file_list <- setNames(clr_paths, method_short_label(method_names))

# lock plotting order to the discovered order
method_levels <- names(file_list)

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

# Simple ranking (for CSV only; plotting keeps original order)
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
  # only use original data columns (exclude metadata)
  X <- safe_numeric_matrix(merged[, setdiff(names(df), "sample_id"), drop = FALSE])
  b <- merged$batch_id
  ascore <- tryCatch(
    compute_alignment_score(X, b, k = K_NEIGHBORS, var_prop = VAR_PROP_MIN, max_pcs = MAX_PCS),
    error = function(e) NA_real_
  )
  as_tbl <- bind_rows(as_tbl, tibble(Method = nm, AS = ascore))
}

# --------- Rank (if applicable) + Save + Plot in original order ---------
only_baseline <- length(method_levels) == 1L && identical(method_levels, "Before correction")

if (only_baseline) {
  # ---- Baseline-only assessment (no ranking) ----
  base_row <- as_tbl %>% filter(Method == "Before correction")
  if (nrow(base_row) == 0) {
    stop("No Alignment Score computed for 'Before correction'.")
  }
  base_row <- base_row %>%
    mutate(Needs_Correction = AS < AS_THRESHOLD)
  
  # save + print assessment
  readr::write_csv(base_row, file.path(output_folder, "alignment_score_raw_assessment.csv"))
  print(base_row)
  
  # single-bar plot in original order
  plot_df <- base_row %>% mutate(Method = factor(Method, levels = method_levels))
  
  p_as <- ggplot(plot_df, aes(x = Method, y = AS, fill = Method)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.4, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.3f", AS)), vjust = -0.4, size = 3.2) +
    scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "Alignment Score (baseline)",
         x = "Method", y = "AS (0- , higher = better mixing)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  ggsave(file.path(output_folder, "alignment_score.png"), p_as, width = 6.5, height = 4.6, dpi = 300)
  ggsave(file.path(output_folder, "alignment_score.tif"), p_as, width = 6.5, height = 4.6, dpi = 300, compression = "lzw")
  
  if (isTRUE(base_row$Needs_Correction[1])) {
    message(sprintf("Baseline AS = %.3f < %.2f - correction recommended.", base_row$AS[1], AS_THRESHOLD))
  } else {
    message(sprintf("Baseline AS = %.3f >= %.2f - correction may not be necessary.", base_row$AS[1], AS_THRESHOLD))
  }
  
} else {
  # ---- Multiple methods: keep bar order as in file_list; still produce ranking CSV ----
  as_ranked <- rank_alignment_methods(as_tbl)
  readr::write_csv(as_ranked, file.path(output_folder, "alignment_score_ranking.csv"))
  print(as_ranked)
  
  # Plot bars in original method order (DO NOT reorder by AS)
  plot_df <- as_tbl %>%
    mutate(Method = factor(Method, levels = method_levels))
  
  p_as <- ggplot(plot_df, aes(x = Method, y = AS, fill = Method)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.4, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.3f", AS)), vjust = -0.4, size = 3.2) +
    scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "Alignment Score",
         x = "Method", y = "AS (0-1, higher = better mixing)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  ggsave(file.path(output_folder, "alignment_score.png"), p_as, width = 8.5, height = 5.2, dpi = 300)
  ggsave(file.path(output_folder, "alignment_score.tif"), p_as, width = 8.5, height = 5.2, dpi = 300, compression = "lzw")
}