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
if (length(args) < 1) stop("Usage: Rscript rmse_batch_heatmap.R <output_folder>")
output_folder <- args[1]

# ==== Read metadata & file list ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))

file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list  <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list  <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

batch_var <- "batchid"   # change this if needed

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

# CLR transform with multiplicative replacement for zeros
clr_transform <- function(X) {
  Xc <- safe_closure(X)
  for (i in seq_len(nrow(Xc))) {
    xi <- Xc[i, ]
    pos <- xi > 0 & is.finite(xi)
    if (!any(pos)) {
      xi[] <- 1 / length(xi); pos <- xi > 0
    }
    if (any(!pos)) {
      minpos <- min(xi[pos], na.rm = TRUE)
      repl <- min(minpos * 0.5, 1e-8)
      xi[!pos] <- repl
      xi <- xi / sum(xi)
    }
    Xc[i, ] <- xi
  }
  L <- log(Xc)
  sweep(L, 1, rowMeans(L), "-")
}

# Prepare data for RMSE (CLR space). If negatives present, treat as already CLR and just row-center.
to_clr_for_rmse <- function(X) {
  if (any(X < 0, na.rm = TRUE)) {
    sweep(X, 1, rowMeans(X, na.rm = TRUE), "-")
  } else {
    clr_transform(X)
  }
}

# Euclidean distance matrix -> RMSE matrix (divide by sqrt(p))
euclidean_to_rmse <- function(D_eucl, p) {
  as.matrix(D_eucl) / sqrt(p)
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

# Compute RMSE batch×batch matrix and return both matrix & long data
rmse_batch_heatmap_data <- function(df, metadata, batch_var = "batchid",
                                    diag_mode = "zero") {
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
  
  # drop columns with NA or zero variance
  keep <- apply(X, 2, function(col) all(is.finite(col)) && sd(col, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) stop("No variable numeric features remain.")
  
  # CLR space for RMSE
  Xclr <- to_clr_for_rmse(X)
  
  # Euclidean -> RMSE
  D_eucl <- dist(Xclr, method = "euclidean")
  D_rmse <- euclidean_to_rmse(D_eucl, p = ncol(Xclr))
  
  # batch factor with numeric-sorted levels
  present_batches <- sort_levels_numeric(unique(dfm[[batch_var]]))
  bfac <- factor(as.character(dfm[[batch_var]]), levels = present_batches)
  
  # batch×batch RMSE (mean pairwise)
  Db <- batch_distance_matrix(D_rmse, bfac, diag_mode = diag_mode)
  
  # enforce numeric order on axes
  ord <- present_batches
  Db  <- Db[ord, ord, drop = FALSE]
  
  # long form for plotting
  long <- as.data.frame(Db) |>
    mutate(batch1 = rownames(Db)) |>
    pivot_longer(-batch1, names_to = "batch2", values_to = "rmse") |>
    mutate(batch1 = factor(batch1, levels = ord),
           batch2 = factor(batch2, levels = ord))
  
  list(Db = Db, long = long, order = ord)
}

# Build a plot from long data with a PANEL-SPECIFIC color scale + legend BELOW the heatmap
rmse_batch_heatmap_plot <- function(long, ord, title_label,
                                    label_digits = 3, cap_quantile = 1.0) {
  # per-panel limits (cap only the upper end if you like)
  panel_min <- min(long$rmse, na.rm = TRUE)
  panel_max <- if (cap_quantile >= 1) {
    max(long$rmse, na.rm = TRUE)
  } else {
    quantile(long$rmse, probs = cap_quantile, na.rm = TRUE, names = FALSE, type = 7)
  }
  if (!is.finite(panel_max) || panel_max <= panel_min) panel_max <- panel_min + 1e-8
  text_threshold <- (panel_min + panel_max) / 2
  
  long <- long %>%
    mutate(
      label   = ifelse(is.na(rmse), "", formatC(rmse, format = "f", digits = label_digits)),
      txt_col = ifelse(!is.na(rmse) & rmse >= text_threshold, "white", "black")
    )
  
  ggplot(long, aes(x = batch2, y = batch1, fill = rmse)) +
    geom_tile() +
    geom_text(aes(label = label, colour = txt_col), size = 3, na.rm = TRUE) +
    scale_colour_identity(guide = "none") +
    scale_fill_viridis_c(
      name = "RMSE (CLR space)", option = "C",
      limits = c(panel_min, panel_max),
      oob = scales::squish,
      na.value = "grey95"    # diagonal shows as light grey
    ) +
    coord_fixed() +
    labs(title = paste0("RMSE (batch × batch) – ", title_label), x = NULL, y = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.ticks  = element_blank(),
      panel.grid  = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.margin = margin(8, 10, 8, 10)
    ) +
    guides(fill = guide_colourbar(title.position = "top"))
}

# ==== Build panels & save (each with its own legend/scale) ====
diag_mode     <- "NA"   # "zero" | "mean" | "NA"
label_digits  <- 2        # decimals in tiles
cap_quantile  <- 1.0      # set <1 (e.g., 0.99) to cap extreme max within a panel

plots <- list()
for (nm in names(file_list)) {
  cat("Processing (RMSE batch heatmap):", nm, "\n")
  df <- read_csv(file_list[[nm]], show_col_types = FALSE)
  
  dat <- rmse_batch_heatmap_data(df, metadata, batch_var = batch_var, diag_mode = diag_mode)
  plt <- rmse_batch_heatmap_plot(dat$long, dat$order, title_label = nm,
                                 label_digits = label_digits, cap_quantile = cap_quantile)
  plots[[nm]] <- plt
}

# ==== Combine all panels (each retains its own legend at the bottom) ====
ncol_grid <- 2
combined <- wrap_plots(plots, ncol = ncol_grid)  # NOTE: no guides='collect' -> per-panel legends kept

ggsave(file.path(output_folder, "rmse_heatmap.png"),
       plot = combined, width = 8.5 * ncol_grid, height = 6 * ceiling(length(plots) / ncol_grid),
       dpi = 300)

ggsave(file.path(output_folder, "rmse_heatmap.tif"),
       plot = combined, width = 8.5 * ncol_grid, height = 6 * ceiling(length(plots) / ncol_grid),
       dpi = 300, compression = "lzw")

# Function to compute the mean RMSE for each method
compute_mean_rmse <- function(long_data) {
  # Calculate the mean RMSE across all batch × batch comparisons
  mean_rmse <- mean(long_data$rmse, na.rm = TRUE)
  return(mean_rmse)
}

# Initialize a data frame to store results
rank_results <- tibble(Method = character(), Mean_RMSE = numeric())

# Loop over all methods to compute the mean RMSE
for (method_name in names(file_list)) {
  cat("Processing:", method_name, "\n")
  
  # Read in the data for the current method
  df <- read_csv(file_list[[method_name]], show_col_types = FALSE)
  
  # Compute RMSE batch × batch matrix and long data for heatmap
  rmse_data <- rmse_batch_heatmap_data(df, metadata, batch_var = batch_var, diag_mode = "NA")
  
  # Compute the mean RMSE for the current method
  mean_rmse <- compute_mean_rmse(rmse_data$long)
  
  # Store the result
  rank_results <- bind_rows(rank_results, tibble(Method = method_name, Mean_RMSE = mean_rmse))
}

# Rank the methods based on the mean RMSE
ranked_methods <- rank_results %>%
  arrange(Mean_RMSE) %>%  # Sort in ascending order of Mean RMSE
  mutate(Rank = row_number())

# Display the ranked methods based on RMSE
ranked_methods

