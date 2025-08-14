# ================= pRDA per-feature R2 boxplots (Batch vs Treatment) =================
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
})
if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package 'vegan' is required. install.packages('vegan')")
}

# --------- Args / IO ---------
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")  # <-- change/remove for CLI use
if (length(args) < 1) stop("Usage: Rscript r2_boxplot_per_feature.R <output_folder>")
output_folder <- args[1]

metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))

# matrices: Raw + normalized_*
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list  <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list  <- c("Before correction" = file.path(output_folder, "raw.csv"), file_list)

# --------- Helpers (CLR + safe R^2) ---------
safe_closure <- function(X) {
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
safe_adjR2 <- function(fit) {
  out <- tryCatch(vegan::RsquareAdj(fit), error = function(e) list(r.squared=NA_real_, adj.r.squared=NA_real_))
  adj <- suppressWarnings(out$adj.r.squared)
  if (is.finite(adj)) return(max(0, adj))
  r <- suppressWarnings(out$r.squared)
  if (is.finite(r)) return(max(0, r))
  0
}

# --------- Core: per-feature pure R2 for Batch and Treatment ---------
# returns long tibble: Feature, Effect (Batch/Treatment), R2
compute_feature_r2_BT <- function(df, meta, batch_col = "batchid", treat_col = "phenotype") {
  # ensure sample_id present & align to metadata
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(meta)) df$sample_id <- meta$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfx <- inner_join(df, meta, by = "sample_id")
  
  # required covariates
  if (!(batch_col %in% names(dfx))) stop(sprintf("Batch column '%s' not in metadata.", batch_col))
  if (!(treat_col %in% names(dfx))) stop(sprintf("Treatment column '%s' not in metadata.", treat_col))
  
  # numeric features
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfx %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  
  # drop non-finite / constant features
  keep <- apply(X, 2, function(z) all(is.finite(z)) && sd(z, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) return(tibble(Feature = character(), Effect = character(), R2 = numeric()))
  
  # response in Aitchison space (or row-centered if already log-ratio)
  Y <- if (any(X < 0, na.rm = TRUE)) sweep(X, 1, rowMeans(X, na.rm = TRUE), "-") else clr_transform(X)
  colnames(Y) <- colnames(X)
  
  # filter rows with NA covariates and make factors
  keep_row <- !is.na(dfx[[batch_col]]) & !is.na(dfx[[treat_col]])
  dfx <- dfx[keep_row, , drop = FALSE]
  Y   <- Y[keep_row, , drop = FALSE]
  
  dfx[[batch_col]] <- factor(dfx[[batch_col]])
  dfx[[treat_col]] <- factor(dfx[[treat_col]])
  if (nlevels(dfx[[treat_col]]) < 2 || nlevels(dfx[[batch_col]]) < 2) {
    return(tibble(Feature = character(), Effect = character(), R2 = numeric()))
  }
  
  # per-feature R2 via pRDA with 1-column response
  get_r2_one <- function(yvec) {
    Y1 <- as.matrix(yvec)
    f_t <- as.formula(paste("Y1 ~", treat_col, "+ Condition(", batch_col, ")"))  # Treatment | Batch
    f_b <- as.formula(paste("Y1 ~", batch_col, "+ Condition(", treat_col, ")"))  # Batch | Treatment
    rt  <- tryCatch(safe_adjR2(vegan::rda(f_t, data = dfx)), error = function(e) NA_real_)
    rb  <- tryCatch(safe_adjR2(vegan::rda(f_b, data = dfx)), error = function(e) NA_real_)
    c(Batch = rb, Treatment = rt)
  }
  
  res <- lapply(seq_len(ncol(Y)), function(j) {
    r2 <- get_r2_one(Y[, j])
    tibble(Feature = colnames(Y)[j],
           Effect  = c("Batch","Treatment"),
           R2      = as.numeric(r2))
  }) %>% bind_rows()
  
  res %>% filter(is.finite(R2), R2 >= 0, R2 <= 1)
}

# --------- Build per-feature R2 across all methods ---------
batch_col <- "batchid"
treat_col <- "phenotype"
if (!("phenotype" %in% names(metadata))) stop("metadata.csv lacks 'phenotype'.")
if (dplyr::n_distinct(metadata$phenotype) < 2) stop("'phenotype' needs at least 2 levels.")

r2_long_df <- lapply(names(file_list), function(nm) {
  message("Per-feature R2 (Batch/Treatment): ", nm)
  df <- read_csv(file_list[[nm]], show_col_types = FALSE)
  out <- compute_feature_r2_BT(df, metadata, batch_col, treat_col)
  out$Method <- nm
  out
}) %>% bind_rows()

# tidy labels/order for plotting
r2_long_df <- r2_long_df %>%
  mutate(
    Effect = case_when(
      tolower(Effect) %in% c("batch","batchid","batch_id") ~ "Batch",
      tolower(Effect) %in% c("treatment","phenotype","group","trt") ~ "Treatment",
      TRUE ~ Effect
    ),
    Effect = factor(Effect, levels = c("Batch","Treatment")),
    Method = factor(Method, levels = names(file_list))
  ) %>%
  filter(!is.na(R2), is.finite(R2), R2 >= 0, R2 <= 1)

# --------- Pretty boxplot function (styled like your box_plot) ---------
r2_box_plot <- function(r2_long_df,
                        legend.title = "Effects",
                        ylab = expression("Explained variance ("*R^2*")"),
                        color.set = c(Batch = "#FF7F0E", Treatment = "#BDBDBD"),
                        x.angle = 45, x.hjust = 1, x.vjust = 1,
                        title = NULL) {
  
  # medians for labels
  med_df <- r2_long_df %>%
    group_by(Method, Effect) %>%
    summarize(med = median(R2), .groups = "drop")
  
  ggplot(r2_long_df, aes(x = Effect, y = R2, fill = Effect)) +
    geom_boxplot(width = 0.7, outlier_size = 0.7) +
    stat_boxplot(geom = "errorbar", width = 0.35) +
    scale_fill_manual(values = color.set, name = legend.title, drop = FALSE) +
    facet_grid(. ~ Method, scales = "free_x", space = "free_x") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(y = ylab, x = NULL, title = title) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = x.angle, hjust = x.hjust, vjust = x.vjust),
      panel.grid   = element_blank(),
      axis.text    = element_text(size = 10),
      axis.title   = element_text(size = 12),
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text   = element_text(size = 10),
      plot.title   = element_text(hjust = 0.5, size = rel(1.2))
    ) +
    # median labels just above each box
    geom_text(
      data = med_df,
      aes(x = Effect, y = pmin(med + 0.03, 0.98), label = sprintf("%.3f", med)),
      inherit.aes = FALSE, size = 3
    )
}

# --------- Plot & save ---------
p_r2 <- r2_box_plot(
  r2_long_df,
  legend.title = "Effects",
  title = "Explained variance by Batch vs Treatment (per feature)"
)

ggsave(file.path(output_folder, "R2_boxplot_per_feature.png"), p_r2, width = 9, height = 5.2, dpi = 300)
ggsave(file.path(output_folder, "R2_boxplot_per_feature.tif"), p_r2, width = 9, height = 5.2, dpi = 300, compression = "lzw")

# Function to rank the methods based on Treatment and Batch R²
rank_r2_methods <- function(r2_long_df) {
  # Compute the median R² for Batch and Treatment
  median_r2 <- r2_long_df %>%
    group_by(Method, Effect) %>%
    summarise(median_R2 = median(R2, na.rm = TRUE), .groups = "drop") %>%
    spread(Effect, median_R2)  # Wide format: one column for Treatment, one for Batch
  
  # Rank methods based on median R² for Treatment (higher is better) and Batch (lower is better)
  ranked_results <- median_r2 %>%
    arrange(desc(Treatment), Batch) %>%  # Prioritize Treatment (higher) and then Batch (lower)
    mutate(Rank = row_number())  # Assign ranks
  
  return(ranked_results)
}

# Rank the methods based on per-feature R² values for Batch and Treatment
ranked_r2_methods <- rank_r2_methods(r2_long_df)

# Display the ranked methods
ranked_r2_methods
