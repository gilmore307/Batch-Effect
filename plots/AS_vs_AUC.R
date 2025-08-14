# ==== Libraries ====
suppressPackageStartupMessages({
  library(FNN)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(caret)
  library(pROC)
  library(glmnet)
})

# ==== Args / config ====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# Your metadata outcome column:
PHENO_COL <- "phenotype"   # <-- uses numeric 0/1 in your metadata.csv
CV_FOLDS  <- 5
CV_REPS   <- 5

set.seed(42)

# ==== Read metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

if (!PHENO_COL %in% names(metadata))
  stop(sprintf("Phenotype column '%s' not found in metadata.csv", PHENO_COL))

# Prepare .outcome as a two-class factor with the POSITIVE class as the FIRST level
# If numeric 0/1, treat 1 as "pos", 0 as "neg"
if (is.numeric(metadata[[PHENO_COL]]) && length(unique(metadata[[PHENO_COL]])) == 2) {
  metadata <- metadata %>%
    mutate(.outcome = factor(ifelse(.data[[PHENO_COL]] == 1, "pos", "neg"),
                             levels = c("pos","neg")))
} else {
  # If it's character/factor with exactly 2 levels, take the second level as positive and relevel it first
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

compute_alignment_score <- function(data, batch_vector, k = 10, pcs = NULL, min_var_prop = 0.8) {
  if (!is.matrix(data) || nrow(data) < 3 || ncol(data) < 2) return(NA_real_)
  pca <- prcomp(data, scale. = TRUE)
  if (is.null(pcs)) {
    varprop <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
    pcs <- max(1, min(which(varprop >= min_var_prop)))
    pcs <- min(pcs, 10, ncol(pca$x))
  } else {
    pcs <- min(pcs, ncol(pca$x))
  }
  coords <- pca$x[, seq_len(pcs), drop = FALSE]
  k_eff <- max(1, min(k, nrow(coords) - 1))
  nn <- get.knn(coords, k = k_eff)$nn.index
  diffs <- vapply(seq_len(nrow(coords)), function(i) {
    neighbor_batches <- batch_vector[nn[i, ]]
    mean(neighbor_batches != batch_vector[i], na.rm = TRUE)
  }, numeric(1))
  mean(diffs, na.rm = TRUE)
}

compute_auc_cv <- function(X, y, folds = 5, reps = 5) {
  # X: matrix (samples x features), y: factor with positive level FIRST
  if (!is.matrix(X) || nrow(X) < 10 || ncol(X) < 2) return(NA_real_)
  df <- as.data.frame(X)
  df$.outcome <- droplevels(y)
  if (nlevels(df$.outcome) != 2) return(NA_real_)
  
  ctrl <- trainControl(
    method = "repeatedcv",
    number = folds,
    repeats = reps,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  set.seed(42)
  fit <- tryCatch(
    train(.outcome ~ ., data = df, method = "glmnet",
          metric = "ROC", trControl = ctrl, tuneLength = 5),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  as.numeric(getTrainPerf(fit)$TrainROC)  # mean CV AUC
}

# ==== Compute AS and AUC per method ====
results <- tibble(Method = character(), AS = numeric(), AUC = numeric())

for (name in names(file_list)) {
  fpath <- file_list[[name]]
  if (!file.exists(fpath)) next
  df_raw <- read_csv(fpath, show_col_types = FALSE)
  
  # Ensure sample_id present (match your earlier behavior)
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
  
  # ---- CRITICAL: build X ONLY from the feature file columns (no metadata) ----
  feature_cols <- setdiff(names(df_raw), "sample_id")
  X <- safe_numeric_matrix(df_merged[, feature_cols, drop = FALSE])
  
  # Outcome and batch
  y <- df_merged$.outcome
  batch <- df_merged$batchid %||% rep(1, nrow(df_merged))  # fallback if batchid missing
  
  ascore <- tryCatch(compute_alignment_score(X, batch, k = 10), error = function(e) NA_real_)
  auc    <- tryCatch(compute_auc_cv(X, y, folds = CV_FOLDS, reps = CV_REPS), error = function(e) NA_real_)
  
  results <- bind_rows(results, tibble(Method = name, AS = ascore, AUC = auc))
}

# ==== Plot: ONLY AS vs AUC scatter (save figure, no tables) ====
p_scatter <- ggplot(results, aes(x = AS, y = AUC, label = Method, color = Method)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.5, 1)) +
  labs(title = "Alignment Score vs. AUC (5x5 CV, glmnet)",
       x = "Alignment Score",
       y = "AUC (mean CV)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(output_folder, "alignment_vs_auc.png"),
       p_scatter, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_folder, "alignment_vs_auc.tif"),
       p_scatter, width = 8, height = 6, dpi = 300)

# Normalize AUC and AS (inverting AS so that higher is better)
results <- results %>%
  mutate(
    Norm_AUC = (AUC - min(AUC)) / (max(AUC) - min(AUC)),        # Normalize AUC to [0, 1]
    Norm_AS = (1 - AS - min(1 - AS)) / (max(1 - AS) - min(1 - AS))  # Normalize AS to [0, 1], lower AS is better
  )

# Choose weights for AUC and AS (you can adjust these based on preference)
weight_auc <- 0.6
weight_as <- 0.4

# Compute combined score based on normalized values and weights
results <- results %>%
  mutate(Combined_Score = weight_auc * Norm_AUC + weight_as * Norm_AS)

# Rank methods based on the combined score
ranked_results <- results %>%
  arrange(desc(Combined_Score)) %>%
  mutate(Rank = row_number())

# Display the ranked methods with their combined scores
ranked_results

