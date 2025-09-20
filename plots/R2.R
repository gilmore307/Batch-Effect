# ===================== One-way ANOVA R^2 boxplots (Batch vs Treatment) — CLR *and* TSS =====================
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(forcats)
})

# ----------------- Args / IO -----------------
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")  # <-- change/remove for CLI use
if (length(args) < 1) stop("Usage: Rscript r2_boxplot_anova_dual.R <output_folder>")
output_folder <- args[1]

metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))

# Collect CLR and TSS sets (+ raw.csv for both)
clr_paths <- list.files(output_folder, pattern = "^normalized_.*_clr\\.csv$", full.names = TRUE)
tss_paths <- list.files(output_folder, pattern = "^normalized_.*_tss\\.csv$", full.names = TRUE)
raw_fp    <- file.path(output_folder, "raw.csv")

name_fun <- function(x, suffix) gsub(paste0("^normalized_|_", suffix, "\\.csv$"), "", basename(x))
file_list_clr <- setNames(clr_paths, name_fun(clr_paths, "clr"))
file_list_tss <- setNames(tss_paths, name_fun(tss_paths, "tss"))

if (file.exists(raw_fp)) {
  file_list_clr <- c("Before correction" = raw_fp, file_list_clr)
  file_list_tss <- c("Before correction" = raw_fp, file_list_tss)
}
if (!length(file_list_clr)) stop("No CLR inputs found (expected 'normalized_*_clr.csv' and/or 'raw.csv').")
if (!length(file_list_tss)) stop("No TSS inputs found (expected 'normalized_*_tss.csv' and/or 'raw.csv').")

# ----------------- Helpers -----------------
safe_closure <- function(X) {
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0
  rs <- rowSums(X, na.rm = TRUE)
  bad <- which(rs == 0 | !is.finite(rs))
  if (length(bad)) { X[bad, ] <- 1 / ncol(X); rs <- rowSums(X, na.rm = TRUE) }
  sweep(X, 1, rs, "/")
}
clr_transform <- function(X) {
  Xc <- safe_closure(X)
  for (i in seq_len(nrow(Xc))) {
    xi <- Xc[i, ]; pos <- xi > 0 & is.finite(xi)
    if (!any(pos)) { xi[] <- 1/length(xi); pos <- xi > 0 }
    if (any(!pos)) {
      m <- min(xi[pos], na.rm = TRUE)
      xi[!pos] <- min(m*0.5, 1e-8)
      xi <- xi / sum(xi)
    }
    Xc[i, ] <- xi
  }
  L <- log(Xc)
  sweep(L, 1, rowMeans(L), "-")
}
# Unadjusted ANOVA R^2 = 1 - SSE/SST (clamped to [0,1])
anova_r2 <- function(y, g) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]; g <- droplevels(factor(g[ok]))
  if (length(y) < 3 || nlevels(g) < 2) return(NA_real_)
  fit <- tryCatch(lm(y ~ g), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  r2 <- summary(fit)$r.squared
  if (!is.finite(r2)) return(NA_real_)
  max(0, min(1, as.numeric(r2)))
}

# ----------------- Core: per-feature one-way ANOVA R^2 (Batch vs Treatment) -----------------
# geometry: "CLR" (Aitchison) or "TSS" (proportions)
compute_anova_r2_BT <- function(df, meta, batch_col = "batchid", treat_col = "phenotype", geometry = c("CLR","TSS")) {
  geometry <- match.arg(geometry)
  # align to metadata
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(meta)) df$sample_id <- meta$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfx <- inner_join(df, meta, by = "sample_id")
  
  if (!(batch_col %in% names(dfx))) stop(sprintf("Batch column '%s' not in metadata.", batch_col))
  if (!(treat_col %in% names(dfx))) stop(sprintf("Treatment column '%s' not in metadata.", treat_col))
  
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfx %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  
  # drop features that are non-finite or constant
  keep <- apply(X, 2, function(z) all(is.finite(z)) && sd(z, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) return(tibble(Feature = character(), Effect = character(), R2 = numeric()))
  
  if (geometry == "CLR") {
    Y <- if (any(X < 0, na.rm = TRUE)) sweep(X, 1, rowMeans(X, na.rm = TRUE), "-") else clr_transform(X)
  } else {
    Y <- safe_closure(X)  # TSS proportions
  }
  colnames(Y) <- colnames(X)
  
  # keep rows with complete covariates
  keep_row <- !is.na(dfx[[batch_col]]) & !is.na(dfx[[treat_col]])
  dfx <- dfx[keep_row, , drop = FALSE]
  Y   <- Y[keep_row, , drop = FALSE]
  
  dfx[[batch_col]] <- factor(dfx[[batch_col]])
  dfx[[treat_col]] <- factor(dfx[[treat_col]])
  
  # compute per-feature R^2
  res <- lapply(seq_len(ncol(Y)), function(j) {
    y <- Y[, j]
    r2_b <- anova_r2(y, dfx[[batch_col]])
    r2_t <- anova_r2(y, dfx[[treat_col]])
    tibble(Feature = colnames(Y)[j],
           Effect  = c("Batch","Treatment"),
           R2      = as.numeric(c(r2_b, r2_t)))
  }) %>% bind_rows()
  
  res %>% filter(is.finite(R2), R2 >= 0, R2 <= 1)
}

# ----------------- Build per-feature R^2 across all methods -----------------
batch_col <- "batchid"
treat_col <- "phenotype"
if (!("phenotype" %in% names(metadata))) stop("metadata.csv lacks 'phenotype'.")
if (dplyr::n_distinct(metadata$phenotype) < 2) stop("'phenotype' needs at least 2 levels.")
if (dplyr::n_distinct(metadata$batchid)   < 2) stop("'batchid' needs at least 2 levels.")

# CLR set
r2_long_clr <- lapply(names(file_list_clr), function(nm) {
  message("Per-feature ANOVA R^2 (CLR): ", nm)
  df <- read_csv(file_list_clr[[nm]], show_col_types = FALSE)
  out <- compute_anova_r2_BT(df, metadata, batch_col, treat_col, geometry = "CLR")
  out$Method <- nm
  out
}) %>% bind_rows()

# TSS set
r2_long_tss <- lapply(names(file_list_tss), function(nm) {
  message("Per-feature ANOVA R^2 (TSS): ", nm)
  df <- read_csv(file_list_tss[[nm]], show_col_types = FALSE)
  out <- compute_anova_r2_BT(df, metadata, batch_col, treat_col, geometry = "TSS")
  out$Method <- nm
  out
}) %>% bind_rows()

# tidy labels/order for plotting
method_levels_clr <- names(file_list_clr)
method_levels_tss <- names(file_list_tss)

tidy_long <- function(df, method_levels) {
  df %>%
    mutate(
      Effect = case_when(
        tolower(Effect) %in% c("batch","batchid","batch_id") ~ "Batch",
        tolower(Effect) %in% c("treatment","phenotype","group","trt") ~ "Treatment",
        TRUE ~ Effect
      ),
      Effect = factor(Effect, levels = c("Batch","Treatment")),
      Method = factor(Method, levels = method_levels)
    ) %>%
    filter(!is.na(R2), is.finite(R2), R2 >= 0, R2 <= 1)
}
r2_long_clr <- tidy_long(r2_long_clr, method_levels_clr)
r2_long_tss <- tidy_long(r2_long_tss, method_levels_tss)

# ----------------- Figure (Figure-7 style boxes) -----------------
box_cols <- c(Batch = "#FF7F0E", Treatment = "#BDBDBD")

make_boxplot <- function(r2_long_df, method_levels, title) {
  med_df <- r2_long_df %>%
    dplyr::group_by(Method, Effect) %>%
    dplyr::summarize(med = median(R2), .groups = "drop")
  
  ggplot(r2_long_df, aes(x = Effect, y = R2, fill = Effect)) +
    geom_boxplot(width = 0.7, outlier.size = 0.7) +   # <-- dot, not underscore
    stat_boxplot(geom = "errorbar", width = 0.35) +
    scale_fill_manual(values = c(Batch = "#FF7F0E", Treatment = "#BDBDBD"),
                      name = "Effect", drop = FALSE) +
    facet_grid(. ~ Method, scales = "free_x", space = "free_x") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(y = expression("One-way ANOVA "*R^2), x = NULL, title = title) +  # title can be expression()
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid   = element_blank(),
      axis.text    = element_text(size = 10),
      axis.title   = element_text(size = 12),
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text   = element_text(size = 10),
      plot.title   = element_text(hjust = 0.5, size = rel(1.2))
    ) +
    geom_text(
      data = med_df,
      aes(x = Effect, y = pmin(med + 0.03, 0.98), label = sprintf("%.3f", med)),
      inherit.aes = FALSE, size = 3
    )
}


p_clr <- make_boxplot(
  r2_long_clr, method_levels_clr,
  expression("Per-feature " * R^2 * " (one-way ANOVA) — CLR (Aitchison)")
)

p_tss <- make_boxplot(
  r2_long_tss, method_levels_tss,
  expression("Per-feature " * R^2 * " (one-way ANOVA) — TSS (proportions)")
)

ggsave(file.path(output_folder, "R2__aitchison.png"), p_clr, width = 10, height = 5.2, dpi = 300)
ggsave(file.path(output_folder, "R2_aitchison.tif"), p_clr, width = 10, height = 5.2, dpi = 300, compression = "lzw")
ggsave(file.path(output_folder, "R2_braycurtis.png"), p_tss, width = 10, height = 5.2, dpi = 300)
ggsave(file.path(output_folder, "R2_braycurtis.tif"), p_tss, width = 10, height = 5.2, dpi = 300, compression = "lzw")

# ----------------- Unified ranking (optional)
# Per geometry, summarize by median R2: want higher Treatment, lower Batch.
score_from_long <- function(df) {
  df %>%
    group_by(Method, Effect) %>%
    summarise(median_R2 = median(R2, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Effect, values_from = median_R2) %>%
    mutate(Score = pmax(0, Treatment) * pmax(0, 1 - Batch)) %>%
    select(Method, Score)
}
score_clr <- score_from_long(r2_long_clr) %>% rename(Score_CLR = Score)
score_tss <- score_from_long(r2_long_tss) %>% rename(Score_TSS = Score)

ranking_unified <- full_join(score_clr, score_tss, by = "Method") %>%
  mutate(
    Combined_Score = dplyr::case_when(
      !is.na(Score_CLR) & !is.na(Score_TSS) ~ sqrt(Score_CLR * Score_TSS),
      is.na(Score_TSS)                       ~ Score_CLR,
      is.na(Score_CLR)                       ~ Score_TSS,
      TRUE                                   ~ NA_real_
    )
  ) %>%
  arrange(desc(Combined_Score)) %>%
  mutate(Rank = row_number())

readr::write_csv(ranking_unified, file.path(output_folder, "r2_ranking.csv"))

print(ranking_unified, n = nrow(ranking_unified))
