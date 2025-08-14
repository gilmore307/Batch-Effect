# ================= pRDA variance partition — stacked bars with fixed-position on-bar labels =================
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
})
if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package 'vegan' is required. install.packages('vegan')")
}

# --------- Args / IO ---------
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")  # change/remove for CLI usage
if (length(args) < 1) stop("Usage: Rscript prda_varpart_plot.R <output_folder>")
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

# --------- Core: compute fractions (named numeric vector of length 4) ---------
compute_prda_parts <- function(df, meta, batch_col = "batchid", treat_col = "phenotype") {
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(meta)) df$sample_id <- meta$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfx <- inner_join(df, meta, by = "sample_id")
  
  if (!(batch_col %in% names(dfx))) stop(sprintf("Batch column '%s' not in metadata.", batch_col))
  if (!(treat_col %in% names(dfx))) stop(sprintf("Treatment column '%s' not in metadata.", treat_col))
  dfx <- dfx %>% filter(!is.na(.data[[batch_col]]), !is.na(.data[[treat_col]]))
  
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfx %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  keep <- apply(X, 2, function(z) all(is.finite(z)) && sd(z, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) return(c(Treatment=0, Intersection=0, Batch=0, Residuals=1))
  
  Y <- if (any(X < 0, na.rm = TRUE)) sweep(X, 1, rowMeans(X, na.rm = TRUE), "-") else clr_transform(X)
  
  dfx[[batch_col]] <- factor(dfx[[batch_col]])
  dfx[[treat_col]] <- factor(dfx[[treat_col]])
  if (nlevels(dfx[[treat_col]]) < 2 || nlevels(dfx[[batch_col]]) < 2) {
    return(c(Treatment=0, Intersection=0, Batch=0, Residuals=1))
  }
  
  Ymat <- Y
  f_both   <- as.formula(paste("Ymat ~", treat_col, "+", batch_col))
  f_t_pure <- as.formula(paste("Ymat ~", treat_col, "+ Condition(", batch_col, ")"))
  f_b_pure <- as.formula(paste("Ymat ~", batch_col, "+ Condition(", treat_col, ")"))
  
  fit_both  <- vegan::rda(f_both,   data = dfx)
  fit_t     <- vegan::rda(f_t_pure, data = dfx)
  fit_b     <- vegan::rda(f_b_pure, data = dfx)
  
  r2_both   <- safe_adjR2(fit_both)
  r2_t_pure <- safe_adjR2(fit_t)
  r2_b_pure <- safe_adjR2(fit_b)
  
  r2_inter  <- max(0, r2_both - r2_t_pure - r2_b_pure)
  r2_resid  <- max(0, 1 - (r2_t_pure + r2_b_pure + r2_inter))
  
  parts <- c(Treatment = r2_t_pure, Intersection = r2_inter,
             Batch = r2_b_pure, Residuals = r2_resid)
  parts[!is.finite(parts)] <- 0
  parts <- pmax(0, parts)
  s <- sum(parts)
  if (s > 0) parts <- parts / s
  parts
}

# --------- Compute for each method -> parts_df ---------
batch_col <- "batchid"
treat_col <- "phenotype"
if (!("phenotype" %in% names(metadata))) stop("metadata.csv has no 'phenotype'.")
if (dplyr::n_distinct(metadata$phenotype) < 2) stop("'phenotype' needs ≥2 levels.")

parts_df <- lapply(names(file_list), function(nm) {
  message("Computing pRDA varpart: ", nm)
  df <- read_csv(file_list[[nm]], show_col_types = FALSE)
  v  <- compute_prda_parts(df, metadata, batch_col = batch_col, treat_col = treat_col)
  tibble::tibble(Method = nm, Fraction = as.numeric(v))  # 4 rows per method
}) |> bind_rows()

# --------- Add Component & fixed-position on-bar labels ----------
component_order <- c("Treatment","Intersection","Batch","Residuals")
abbr <- c(Treatment="T", Intersection="I", Batch="B", Residuals="R")
label_pos_pct <- c(Treatment=95, Intersection=65, Batch=35, Residuals=5)  # % heights

stopifnot(all(dplyr::count(parts_df, Method)$n == 4))

parts_df <- parts_df %>%
  arrange(Method) %>%
  group_by(Method) %>%
  mutate(Component = factor(component_order[row_number()], levels = component_order)) %>%
  ungroup() %>%
  mutate(
    Method   = factor(Method, levels = names(file_list)),
    Fraction = pmin(pmax(Fraction, 0), 1),
    label    = paste0(abbr[as.character(Component)], " = ", round(Fraction*100)),
    y_pos    = label_pos_pct[as.character(Component)] / 100
  )

# --------- Plot (labels ON the bars at the fixed y positions) ---------
cols <- c(
  "Residuals"    = "#1F77B4",
  "Batch"        = "#FF7F0E",
  "Intersection" = "#FFD54F",
  "Treatment"    = "#BDBDBD"
)

p <- ggplot(parts_df, aes(x = Method, y = Fraction, fill = Component)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.4) +
  scale_fill_manual(values = cols,
                    breaks = c("Residuals","Batch","Intersection","Treatment"),
                    name = "Variation sources") +
  # fixed-position labels (on-bar)
  geom_text(aes(y = y_pos, label = label),
            size = 3.2, vjust = 0.5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.05),  # 105% headroom as requested
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Methods", y = "Explained variance (%)") +
  theme_bw() +
  theme(
    legend.position    = "right",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

# --------- Save ---------
ggsave(file.path(output_folder, "pRDA.png"),
       p, width = 7.2, height = 5.6, dpi = 300)
ggsave(file.path(output_folder, "pRDA.tif"),
       p, width = 7.2, height = 5.6, dpi = 300, compression = "lzw")

# Function to rank the methods based on Treatment and Batch variance
rank_prda_methods <- function(parts_df) {
  # Compute the rank based on Treatment and Batch fractions
  # We rank based on Treatment fraction first (higher is better), then Batch fraction (lower is better)
  
  ranked_results <- parts_df %>%
    group_by(Method) %>%
    summarise(
      Treatment = sum(Fraction[Component == "Treatment"], na.rm = TRUE),
      Batch = sum(Fraction[Component == "Batch"], na.rm = TRUE)
    ) %>%
    arrange(desc(Treatment), Batch) %>%  # Prioritize Treatment, then Batch
    mutate(Rank = row_number())  # Assign ranks based on sorted order
  
  return(ranked_results)
}

# Rank the methods based on pRDA variance partitioning
ranked_prda_methods <- rank_prda_methods(parts_df)

# Display the ranked methods
ranked_prda_methods
