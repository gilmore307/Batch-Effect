# ==== Libraries ====
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ==== IO ====
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")
if (length(args) < 1) stop("Usage: Rscript bc_batch_heatmap.R <output_folder>")
output_folder <- args[1]

# ==== Read metadata & file list ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))

file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list  <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list  <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

batch_var <- "batchid"   # change here if your batch column is named differently

# ==== Helpers ====
# Sort character/factor levels numerically if possible, else alphabetically
sort_levels_numeric <- function(x) {
  x <- as.character(x)
  xn <- suppressWarnings(as.numeric(x))
  if (all(!is.na(xn))) as.character(sort(xn)) else sort(x, method = "radix")
}

# Row-wise closure to 1; if a row sums to 0, make it uniform
safe_closure <- function(X) {
  rs <- rowSums(X, na.rm = TRUE)
  bad <- which(rs == 0 | !is.finite(rs))
  if (length(bad)) {
    X[bad, ] <- 1 / ncol(X)
    rs <- rowSums(X, na.rm = TRUE)
  }
  sweep(X, 1, rs, "/")
}

# If matrix has negatives (CLR/log), back-transform to composition for Bray
to_composition_for_bray <- function(X) {
  if (any(X < 0, na.rm = TRUE)) {
    Xp <- exp(X)
    safe_closure(Xp)
  } else {
    safe_closure(X)
  }
}

# Bray–Curtis distance (rows = samples)
bray_dist <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  rs <- rowSums(X)
  d <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      num <- sum(abs(X[i, ] - X[j, ]))
      den <- rs[i] + rs[j]
      d[i, j] <- d[j, i] <- if (den > 0) num / den else 0
    }
  }
  stats::as.dist(d)
}

# Average pairwise sample distances to batch×batch matrix
batch_distance_matrix <- function(D_sample, batch_factor, diag_mode = c("zero","mean","NA")) {
  diag_mode <- match.arg(diag_mode)
  M <- as.matrix(D_sample)
  b_levels <- levels(batch_factor)
  B <- length(b_levels)
  Db <- matrix(NA_real_, B, B, dimnames = list(b_levels, b_levels))
  
  for (i in seq_len(B)) {
    idx_i <- which(batch_factor == b_levels[i])
    for (j in seq_len(B)) {
      idx_j <- which(batch_factor == b_levels[j])
      if (i == j) {
        if (diag_mode == "zero") {
          Db[i, j] <- 0
        } else if (diag_mode == "mean") {
          if (length(idx_i) >= 2) {
            subM <- M[idx_i, idx_i, drop = FALSE]
            Db[i, j] <- mean(subM[upper.tri(subM)], na.rm = TRUE)
          } else Db[i, j] <- 0
        } else {
          Db[i, j] <- NA_real_
        }
      } else {
        Db[i, j] <- mean(M[idx_i, idx_j, drop = FALSE], na.rm = TRUE)
      }
    }
  }
  Db
}

# ==== One panel: batch × batch BC heatmap with cell labels ====
bc_batch_heatmap_panel <- function(df, metadata, title_label = "Matrix",
                                   batch_var = "batchid", diag_mode = "zero",
                                   label_digits = 2, text_threshold = 0.6) {
  
  # attach sample_id if missing
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(metadata)) df$sample_id <- metadata$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfm <- inner_join(df, metadata, by = "sample_id")
  
  # numeric feature matrix (rows = samples)
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfm %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  # drop constant columns
  keep <- apply(X, 2, function(x) sd(x, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) stop("No variable numeric features remain.")
  
  # composition (Bray requirement)
  Xcomp <- to_composition_for_bray(X)
  
  # sample-level Bray
  D_samp <- bray_dist(Xcomp)
  
  # present batch levels (numeric order 0,1,2,...)
  present_batches <- sort_levels_numeric(unique(dfm[[batch_var]]))
  bfac <- factor(as.character(dfm[[batch_var]]), levels = present_batches)
  
  # batch×batch average distances
  Db <- batch_distance_matrix(D_samp, bfac, diag_mode = diag_mode)
  
  # enforce numeric order on axes
  ord <- present_batches
  Db  <- Db[ord, ord, drop = FALSE]
  
  # long form + labels and text color
  long <- as.data.frame(Db) |>
    mutate(batch1 = rownames(Db)) |>
    pivot_longer(-batch1, names_to = "batch2", values_to = "bc") |>
    mutate(
      batch1 = factor(batch1, levels = ord),
      batch2 = factor(batch2, levels = ord),
      label  = ifelse(is.na(bc), "", formatC(bc, format = "f", digits = label_digits)),
      txt_col = ifelse(!is.na(bc) & bc >= text_threshold, "white", "black")
    )
  
  ggplot(long, aes(x = batch2, y = batch1, fill = bc)) +
    geom_tile() +
    geom_text(aes(label = label, colour = txt_col), size = 3) +
    scale_colour_identity(guide = "none") +
    scale_fill_viridis_c(name = "Bray–Curtis", option = "C",
                         limits = c(0, 1), oob = scales::squish) +
    coord_fixed() +
    labs(title = paste0("BC (batch × batch) – ", title_label), x = NULL, y = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.ticks  = element_blank(),
      panel.grid  = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(8, 10, 8, 10)
    ) +
    guides(fill = guide_colourbar(title.position = "top"))
}

# ==== Build panels & save ====
diag_mode <- "zero"  # "zero" | "mean" | "NA" (for diagonal)
label_digits <- 2    # decimals shown in each cell
text_threshold <- 0.6  # switch text to white if bc >= 0.6

plots <- list()
for (nm in names(file_list)) {
  cat("Processing (BC batch heatmap):", nm, "\n")
  df <- read_csv(file_list[[nm]], show_col_types = FALSE)
  
  plt <- bc_batch_heatmap_panel(
    df, metadata,
    title_label   = nm,
    batch_var     = batch_var,
    diag_mode     = diag_mode,
    label_digits  = label_digits,
    text_threshold = text_threshold
  )
  plots[[nm]] <- plt
}

# ==== Combine all methods on one page with ONE colorbar at bottom ====
ncol_grid <- 2
combined <- wrap_plots(plots, ncol = ncol_grid) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave(file.path(output_folder, "bc_heatmap.png"),
       plot = combined, width = 8.5 * ncol_grid, height = 6 * ceiling(length(plots) / ncol_grid),
       dpi = 300)

ggsave(file.path(output_folder, "bc_heatmap.tif"),
       plot = combined, width = 8.5 * ncol_grid, height = 6 * ceiling(length(plots) / ncol_grid),
       dpi = 300, compression = "lzw")

# Function to calculate average Bray-Curtis dissimilarity for each method
average_bray_curtis <- function(batch_distance_matrix) {
  # Ensure the matrix is numeric
  batch_distance_matrix <- as.numeric(batch_distance_matrix)
  
  # Check if the matrix is not empty
  if (length(batch_distance_matrix) == 0) {
    return(NA)  # Return NA if the matrix is empty
  }
  
  # Exclude diagonal elements (self-comparisons) and calculate the mean
  return(mean(batch_distance_matrix[lower.tri(batch_distance_matrix)], na.rm = TRUE))
}

# Create a vector to store average BC dissimilarities for each method
avg_bray_curtis_values <- sapply(plots, function(plt) {
  # Extract the Bray-Curtis dissimilarity matrix from the plot object
  batch_dist_matrix <- plt$data %>%
    pivot_wider(names_from = batch2, values_from = bc) %>%
    select(-batch1) %>%
    as.matrix()
  
  # Ensure the matrix is numeric before processing
  batch_dist_matrix <- as.numeric(batch_dist_matrix)
  
  # Calculate the average dissimilarity for this method
  average_bray_curtis(batch_dist_matrix)
})

# Rank the methods based on average Bray-Curtis dissimilarity
method_ranking <- sort(avg_bray_curtis_values, na.last = TRUE)

# Display the ranking
method_ranking
