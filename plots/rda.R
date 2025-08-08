# ==== Load Required Libraries ====
library(ggplot2)
library(patchwork)
library(readr)
library(dplyr)
library(ggfortify)   # still used for your PCA, not needed for RDA
library(cluster)
library(vegan)       # <— RDA

# ==== Read UID argument ====
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")
if (length(args) < 1) stop("Usage: Rscript rda.R <output_folder>")
output_folder <- args[1]

# ==== Read Metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"))
metadata$sample_id <- as.character(metadata$sample_id)

# ==== Find All Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== RDA Plot Function ====
# predictors: metadata columns used as constraints in the RDA formula
# group_col: which column to color/ellipse by in the plot (often one of predictors)
rda_with_metadata_plot <- function(df,
                                   method_name,
                                   metadata,
                                   predictors = c("batchid"),
                                   group_col = predictors[1],
                                   transform = c("hellinger", "none", "log1p")) {
  transform <- match.arg(transform)
  
  # merge by sample_id (assumes df has same row order as metadata used before)
  df$sample_id <- metadata$sample_id
  df_merged <- inner_join(df, metadata, by = "sample_id")
  
  # response matrix (features)
  Y <- df_merged %>%
    select(where(is.numeric)) %>%
    select_if(~ sd(.) > 0)
  
  if (ncol(Y) < 2 || nrow(Y) < 3) {
    return(ggplot() + ggtitle(paste("RDA -", method_name, "(not enough variation)")))
  }
  
  # optional transforms (recommended for compositional data)
  Yt <- switch(transform,
               hellinger = decostand(Y, method = "hellinger"),
               log1p     = log1p(Y),
               none      = as.matrix(Y))
  
  # build formula dynamically: Yt ~ pred1 + pred2 + ...
  # (Note: predictors can be factors or numeric; rda() handles both)
  if (!all(predictors %in% names(df_merged))) {
    missing <- paste(setdiff(predictors, names(df_merged)), collapse = ", ")
    stop(paste("Predictors missing from metadata:", missing))
  }
  form <- as.formula(paste("Yt ~", paste(predictors, collapse = " + ")))
  
  # fit RDA
  fit <- rda(form, data = df_merged)
  
  # site scores for constrained axes
  site_scores <- scores(fit, display = "sites", choices = 1:2)
  plot_df <- cbind(as.data.frame(site_scores),
                   group = as.factor(df_merged[[group_col]]))
  
  # variance explained by constrained axes
  eig <- fit$CCA$eig
  prop <- if (length(eig) >= 2) eig / sum(eig) else c(1, 0)
  xlab <- paste0("RDA1 (", sprintf("%.1f", 100 * prop[1]), "%)")
  ylab <- paste0("RDA2 (", sprintf("%.1f", 100 * prop[2]), "%)")
  
  ggplot(plot_df, aes(x = RDA1, y = RDA2, color = group)) +
    geom_point(size = 1.8, alpha = 0.9) +
    stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.15, color = NA, type = "norm") +
    coord_equal() +
    labs(title = paste("RDA -", method_name),
         x = xlab, y = ylab, color = group_col, fill = group_col) +
    theme_minimal() +
    theme(legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
}

# ==== Generate RDA Plots ====
# Choose your RDA constraints here (example: just batchid)
constraints <- c("batchid")
group_for_color <- "batchid"  # color/ellipses by this column

plots <- lapply(names(file_list), function(name) {
  df <- read_csv(file_list[[name]])
  rda_with_metadata_plot(
    df = df,
    method_name = name,
    metadata = metadata,
    predictors = constraints,
    group_col = group_for_color,
    transform = "hellinger"  # change to "none" or "log1p" if you prefer
  )
})

# ==== Combine and Save ====
ncol_grid <- 3
nrow_grid <- ceiling(length(plots) / ncol_grid)
subplot_size <- 6  # inches per subplot (square)

combined_plot <- wrap_plots(plots, ncol = ncol_grid)

ggsave(file.path(output_folder, "rda.tif"), combined_plot,
       width = subplot_size * ncol_grid,
       height = subplot_size * nrow_grid,
       dpi = 300)

ggsave(file.path(output_folder, "rda.png"), combined_plot,
       width = subplot_size * ncol_grid,
       height = subplot_size * nrow_grid,
       dpi = 300)
