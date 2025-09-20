# ==== Libraries ====
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)      # for pivot_wider()
library(gridExtra)
library(grid)       # for unit()

# ==== Args / config (for input and output folder) ====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# Your metadata outcome column:
PHENO_COL <- "phenotype"   # <-- uses numeric 0/1 in your metadata.csv

# ==== Read metadata (metadata.csv) ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

# Check if the phenotype column exists
if (!PHENO_COL %in% names(metadata))
  stop(sprintf("Phenotype column '%s' not found in metadata.csv", PHENO_COL))

# Prepare .outcome as a two-class factor with the POSITIVE class as the FIRST level
# If numeric 0/1, treat 1 as "Positive", 0 as "Negative"
if (is.numeric(metadata[[PHENO_COL]]) && length(unique(metadata[[PHENO_COL]])) == 2) {
  metadata <- metadata %>%
    mutate(.outcome = factor(ifelse(.data[[PHENO_COL]] == 1, "Positive", "Negative"),
                             levels = c("Positive","Negative")))
} else {
  # If it's character/factor with exactly 2 levels, take the second level as positive and relevel it first
  levs <- levels(factor(metadata[[PHENO_COL]]))
  if (length(levs) != 2) stop(sprintf("Phenotype column '%s' must have exactly 2 classes.", PHENO_COL))
  pos <- levs[2]
  metadata <- metadata %>%
    mutate(.outcome = relevel(factor(.data[[PHENO_COL]]), ref = pos))
}

# ==== Collect files (e.g., normalized files) ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
raw_path <- file.path(output_folder, "raw.csv")

# Choose primary feature file (prefer raw.csv, else first normalized_*.csv)
primary_file <- if (file.exists(raw_path)) {
  raw_path
} else if (length(file_paths) > 0) {
  file_paths[1]
} else {
  stop("No feature files found: expected 'raw.csv' or at least one 'normalized_*.csv' in ", output_folder)
}

# Keep a named list (optional; first element is the one we read)
file_list <- c(Raw = raw_path, setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths))))

# ==== Merge feature data with metadata (df_merged) ====
df_raw <- read_csv(primary_file, show_col_types = FALSE)

# Ensure sample_id is included and matches metadata
if (!"sample_id" %in% names(df_raw)) {
  if (nrow(df_raw) == nrow(metadata)) {
    df_raw <- df_raw %>% mutate(sample_id = metadata$sample_id)
  } else {
    warning(sprintf("File %s has no sample_id and row count doesn't match metadata; skipping.", primary_file))
  }
}

# Merge metadata and feature data
df_merged <- df_raw %>%
  mutate(sample_id = as.character(sample_id)) %>%
  inner_join(metadata, by = "sample_id")

# Ensure that batchid and .outcome are factors for proper plotting
df_merged <- df_merged %>%
  mutate(batchid = factor(batchid),
         .outcome = factor(.outcome))

# ==== Prepare the mosaic data ====
mosaic_data <- df_merged %>%
  group_by(batchid, .outcome) %>%
  summarise(count = n(), .groups = "drop")

# Calculate the total count to get proportions
total_count <- sum(mosaic_data$count)

# Add proportion column
mosaic_data <- mosaic_data %>%
  mutate(proportion = count / total_count)

# ==== Define the Mosaic Plot Function (separate palettes for batch vs outcome) ====
mbecMosaicPlot <- function(study.summary, model.vars) {
  
  # Color pools
  mbecCols <- c("#9467bd", "#BCBD22", "#2CA02C", "#E377C2", "#1F77B4", "#FF7F0E",
                "#AEC7E8", "#FFBB78", "#98DF8A", "#D62728", "#FF9896", "#C5B0D5",
                "#8C564B", "#C49C94", "#F7B6D2", "#7F7F7F", "#C7C7C7", "#DBDB8D",
                "#17BECF", "#9EDAE5")
  
  batchCols   <- mbecCols
  outcomeCols <- c(Positive = "#222222", Negative = "#BBBBBB")  # distinct from batch pool
  
  # Safety check: no overlap allowed
  if (length(intersect(toupper(batchCols), toupper(outcomeCols))) > 0) {
    stop("Outcome colors overlap with batch color pool; choose different outcome colors.")
  }
  
  main_color <- "#004B5A"
  legend.cex <- 0.7
  legend.title.cex <- 0.75
  
  # ----- Facet 1: bars colored by batch (legend = Batch) -----
  plot.v2 <- ggplot(
    study.summary, aes(x = batchid, y = proportion, group = .outcome, fill = batchid)
  ) +
    facet_grid(cols = vars(.outcome), scales = "free", space = "free_x", drop = TRUE) +
    geom_bar(stat = "identity", width = 0.9) +
    guides(fill = guide_legend(title = "Batch", reverse = TRUE, keywidth = 1, keyheight = 1)) +
    scale_fill_manual(values = batchCols) +
    ylab("Proportion of all observations") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(colour = main_color, size = 12),
          axis.ticks = element_blank(),
          axis.line = element_line(color = "#7F7F7F"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 90),
          legend.position = 'bottom', legend.box = 'horizontal',
          legend.direction = 'horizontal',
          legend.key.height = unit(0.2, 'cm'),
          legend.key.width  = unit(0.1, 'cm'),
          legend.title = element_text(size = rel(legend.title.cex)),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.text  = element_text(size = rel(legend.cex)),
          plot.margin = unit(c(0.2, 0.2, 0.05, 0.2), "cm"))
  
  # ----- Facet 2: bars colored by outcome (legend = Outcome) -----
  plot.v1 <- ggplot(
    study.summary, aes(x = .outcome, y = proportion, fill = .outcome)
  ) +
    facet_grid(cols = vars(batchid), scales = "free", space = "free_x", drop = TRUE) +
    geom_bar(stat = "identity", width = 0.9) +
    guides(fill = guide_legend(title = "Outcome", reverse = TRUE, keywidth = 1, keyheight = 1)) +
    scale_fill_manual(values = outcomeCols, breaks = names(outcomeCols)) +
    ylab("Proportion of all observations") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(colour = main_color, size = 12),
          axis.ticks = element_blank(),
          axis.line = element_line(color = "#7F7F7F"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 90),
          legend.position = 'bottom', legend.box = 'horizontal',
          legend.direction = 'horizontal',
          legend.key.height = unit(0.2, 'cm'),
          legend.key.width  = unit(0.1, 'cm'),
          legend.title = element_text(size = rel(legend.title.cex)),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.text  = element_text(size = rel(legend.cex)),
          plot.margin = unit(c(0.05, 0.2, 0.2, 0.2), "cm"))
  
  # Legend extractor
  g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    tmp$grobs[[which(vapply(tmp$grobs, function(x) x$name, FUN.VALUE = character(1)) == "guide-box")]]
  }
  legend.v2 <- g_legend(plot.v2)
  legend.v1 <- g_legend(plot.v1)
  
  gridExtra::grid.arrange(
    plot.v2 + theme(legend.position = "none"),
    plot.v1 + theme(legend.position = "none"),
    gridExtra::grid.arrange(legend.v1, legend.v2, ncol = 2, nrow = 1),
    ncol = 1, nrow = 3, widths = c(1), heights = c(1, 1, 0.2), padding = -10
  )
}

# ==== Plot the Mosaic Plot ====
plot.mosaic <- mbecMosaicPlot(study.summary = mosaic_data, model.vars = c('batchid', '.outcome'))
ggsave(file.path(output_folder, "mosaic_plot.png"), plot = plot.mosaic, width = 12, height = 8, dpi = 300)
ggsave(file.path(output_folder, "mosaic_plot.tif"), plot = plot.mosaic, width = 12, height = 8, dpi = 300, compression = "lzw")

# ==== Statistical Test to detect imbalance ====
# Conduct a Chi-Square test to assess if phenotype distribution differs across batches
contingency_table <- mosaic_data %>%
  select(batchid, .outcome, count) %>%
  tidyr::pivot_wider(names_from = .outcome, values_from = count, values_fill = 0) %>%
  as.data.frame()

rownames(contingency_table) <- contingency_table$batchid
contingency_table$batchid <- NULL

# Perform Chi-Square Test
chi_test_result <- chisq.test(as.matrix(contingency_table))

# Report
message(sprintf("Chi-square p-value: %.4g", chi_test_result$p.value))
if (chi_test_result$p.value < 0.05) {
  message("Batch effects are significant. Consider applying batch correction.")
} else {
  message("Batch effects are not significant. No batch correction needed.")
}
