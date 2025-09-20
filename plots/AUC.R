# ================= AUROC curves (CLR only) - per method, with AUC ranking =================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(caret)
  library(randomForest)
  library(pROC)
})


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

# --------- Args / config ---------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

PHENO_COL     <- "phenotype"  # binary outcome in metadata.csv (0/1 or 2-level factor)
CV_FOLDS      <- 5
CV_REPS       <- 5
AUC_THRESHOLD <- 0.70         # baseline-only: if AUC < threshold, recommend correction
set.seed(42)

# --------- Load metadata & make .outcome (positive class first) ---------
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))

if (!PHENO_COL %in% names(metadata))
  stop(sprintf("Phenotype column '%s' not found in metadata.csv", PHENO_COL))

if (is.numeric(metadata[[PHENO_COL]]) && dplyr::n_distinct(metadata[[PHENO_COL]]) == 2) {
  # numeric 0/1 -> 1 = positive
  metadata <- metadata |>
    mutate(.outcome = factor(ifelse(.data[[PHENO_COL]] == 1, "pos", "neg"),
                             levels = c("pos","neg")))
} else {
  levs <- levels(factor(metadata[[PHENO_COL]]))
  if (length(levs) != 2) stop(sprintf("'%s' must have exactly 2 classes.", PHENO_COL))
  # take second level as positive, put it first
  pos <- levs[2]
  metadata <- metadata |>
    mutate(.outcome = relevel(factor(.data[[PHENO_COL]]), ref = pos))
}

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
method_levels <- names(file_list)  # preserve original order for legend

# --------- Helpers ---------
safe_numeric_matrix <- function(df) {
  num <- dplyr::select(df, where(is.numeric))
  if (!ncol(num)) return(matrix(numeric(0), nrow = nrow(df)))
  ok  <- vapply(num, function(z) all(is.finite(z)) && sd(z) > 0, logical(1))
  as.matrix(num[, ok, drop = FALSE])
}

fit_rf_cv <- function(df_features, y_factor, folds = 5, reps = 5) {
  # df_features: data.frame of numeric features; y_factor: factor with positive level FIRST
  train_df <- as.data.frame(df_features)
  train_df$.outcome <- droplevels(y_factor)
  
  ctrl <- trainControl(
    method = "repeatedcv",
    number = folds,
    repeats = reps,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"   # keep out-of-fold predictions for ROC
  )
  
  set.seed(42)
  train(
    .outcome ~ ., data = train_df,
    method = "rf",
    metric = "ROC",
    trControl = ctrl,
    tuneLength = 5,
    importance = TRUE
  )
}

# --------- Train per method, compute AUROC & ROC coords ---------
roc_rows <- list()
auc_rows <- list()

for (nm in names(file_list)) {
  fp <- file_list[[nm]]
  dat <- read_csv(fp, show_col_types = FALSE)
  
  # ensure merge key
  if (!"sample_id" %in% names(dat)) {
    if (nrow(dat) == nrow(metadata)) {
      dat$sample_id <- metadata$sample_id
    } else {
      warning(sprintf("Skipping %s: no sample_id and row count mismatch.", nm))
      next
    }
  }
  
  merged <- dat |>
    mutate(sample_id = as.character(sample_id)) |>
    inner_join(metadata, by = "sample_id")
  
  if (!nrow(merged)) {
    warning(sprintf("Skipping %s: no overlap with metadata.", nm))
    next
  }
  
  feature_cols <- setdiff(names(dat), "sample_id")
  X <- safe_numeric_matrix(merged[, feature_cols, drop = FALSE])
  y <- merged$.outcome
  if (!is.matrix(X) || ncol(X) < 2 || nrow(X) < 10) {
    warning(sprintf("Skipping %s: not enough usable features/samples.", nm))
    next
  }
  
  fit <- tryCatch(fit_rf_cv(X, y, folds = CV_FOLDS, reps = CV_REPS),
                  error = function(e) NULL)
  if (is.null(fit)) {
    warning(sprintf("Skipping %s: model fitting failed.", nm))
    next
  }
  
  # out-of-fold predictions at best tuning
  bt <- fit$bestTune
  pred <- fit$pred
  for (col in names(bt)) pred <- pred[pred[[col]] == bt[[col]], , drop = FALSE]
  # probability column for positive class is named by its level ("pos")
  if (!all(c("obs","pos") %in% names(pred))) {
    pos_level <- levels(y)[1]  # should be "pos"
    if (!pos_level %in% names(pred)) {
      warning(sprintf("Skipping %s: could not find positive class prob column.", nm))
      next
    }
    names(pred)[names(pred) == pos_level] <- "pos"
  }
  
  roc_obj <- pROC::roc(response = pred$obs,
                       predictor = pred$pos,
                       levels = c("neg","pos"),
                       direction = "<",
                       quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  rcoords <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Method = nm
  )
  roc_rows[[nm]] <- rcoords
  auc_rows[[nm]] <- data.frame(Method = nm, AUC = auc_val, stringsAsFactors = FALSE)
}

roc_df <- dplyr::bind_rows(roc_rows)
auc_tbl <- dplyr::bind_rows(auc_rows)

# if nothing computed:
if (!nrow(auc_tbl)) stop("No AUROC results could be computed.")

# --------- Baseline-only path (no ranking) ---------
only_baseline <- (length(method_levels) == 1L && identical(method_levels, "Before correction"))
if (only_baseline) {
  base_auc <- auc_tbl %>% filter(Method == "Before correction")
  if (!nrow(base_auc)) stop("No AUROC computed for 'Before correction'.")
  base_auc <- base_auc %>%
    mutate(Needs_Correction = AUC < AUC_THRESHOLD)
  readr::write_csv(base_auc, file.path(output_folder, "auroc_raw_assessment.csv"))
  print(base_auc)
  
  # Plot single ROC, keep original labeling/order
  label_map <- setNames(sprintf("%s (AUC=%.3f)", base_auc$Method, base_auc$AUC),
                        base_auc$Method)
  roc_df <- roc_df %>% filter(Method == "Before correction")
  roc_df$Legend <- factor(label_map[as.character(roc_df$Method)], levels = unname(label_map))
  
  p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Legend)) +
    geom_path(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    labs(title = "ROC curve (baseline)",
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)",
         color = NULL) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_folder, "auroc.png"), p_roc, width = 6.5, height = 5, dpi = 300)
  ggsave(file.path(output_folder, "auroc.tif"), p_roc, width = 6.5, height = 5, dpi = 300, compression = "lzw")
  
  if (isTRUE(base_auc$Needs_Correction[1])) {
    message(sprintf("Baseline AUC = %.3f < %.2f - correction recommended.", base_auc$AUC[1], AUC_THRESHOLD))
  } else {
    message(sprintf("Baseline AUC = %.3f >= %.2f - correction may not be necessary.", base_auc$AUC[1], AUC_THRESHOLD))
  }
  
} else {
  # --------- Multi-method: keep legend order as file_list (no re-ranking in plot) ---------
  # Save AUC ranking CSV, but DO NOT reorder legend/lines in the plot
  auc_ranked <- auc_tbl %>% arrange(desc(AUC)) %>% mutate(Rank = row_number())
  readr::write_csv(auc_ranked, file.path(output_folder, "auroc_ranking.csv"))
  print(auc_ranked)
  
  # Keep only methods that actually produced ROC coords, in original discovery order
  methods_plotted <- intersect(method_levels, unique(roc_df$Method))
  
  # label text includes AUC values but order follows methods_plotted
  auc_vec <- setNames(auc_tbl$AUC, auc_tbl$Method)
  labels_in_order <- sapply(methods_plotted, function(m) sprintf("%s (AUC=%.3f)", m, auc_vec[[m]]))
  label_map <- setNames(labels_in_order, methods_plotted)
  
  roc_df <- roc_df %>% filter(Method %in% methods_plotted)
  roc_df$Legend <- factor(label_map[as.character(roc_df$Method)], levels = unname(label_map))
  
  p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Legend)) +
    geom_path(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    labs(title = "ROC curves (5 Folds * 5 Repeats)",
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)",
         color = NULL) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_folder, "auroc.png"), p_roc, width = 8.8, height = 6.2, dpi = 300)
  ggsave(file.path(output_folder, "auroc.tif"), p_roc, width = 8.8, height = 6.2, dpi = 300, compression = "lzw")
  
  # quick console summary
  auc_annot <- auc_ranked %>% mutate(Label = sprintf("%s (AUC=%.3f)", Method, AUC)) %>% pull(Label)
  cat("AUC ranking:\n", paste(sprintf("%2d. %s", seq_along(auc_annot), auc_annot), collapse = "\n"), "\n")
}