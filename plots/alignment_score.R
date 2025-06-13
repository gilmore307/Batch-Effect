# ==== Load Libraries ====
library(FNN)
library(readr)
library(dplyr)
library(ggplot2)

# ==== Handle Argument ====
#args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")
if (length(args) < 1) stop("Usage: Rscript alignment_score.R <output_folder>")
output_folder <- args[1]

# ==== Read Metadata ====
metadata_path <- file.path(output_folder, "metadata.csv")
metadata <- read_csv(metadata_path)
metadata$sample_id <- as.character(metadata$sample_id)

# ==== Find All Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))

# Add Raw file explicitly
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== Alignment Score Function ====
compute_alignment_score <- function(data, batch_vector, k = 10) {
  pca <- prcomp(data, scale. = TRUE)
  coords <- pca$x[, 1:10]
  neighbors <- get.knn(coords, k = k)$nn.index
  mean(sapply(1:nrow(data), function(i) {
    neighbor_batches <- batch_vector[neighbors[i, ]]
    mean(neighbor_batches != batch_vector[i])
  }))
}

# ==== Compute Scores ====
alignment_scores <- data.frame(Method = character(), Score = numeric())

for (name in names(file_list)) {
  df <- read_csv(file_list[[name]])
  df$sample_id <- metadata$sample_id
  df_merged <- inner_join(df, metadata, by = "sample_id")
  
  X <- df_merged %>%
    select(where(is.numeric)) %>%
    select_if(~ sd(.) > 0) %>%
    as.matrix()
  
  batch <- df_merged$batchid
  score <- compute_alignment_score(X, batch, k = 10)
  alignment_scores <- rbind(alignment_scores, data.frame(Method = name, Score = score))
}

# ==== Plot and Save ====
plot <- ggplot(alignment_scores, aes(x = reorder(Method, Score), y = Score, fill = Method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "Comparison of Alignment Scores",
    x = "Method",
    y = "Alignment Score"
  ) +
  ylim(0, 1) +
  theme_minimal(base_size = 14)

# Save as TIFF (high-quality)
tiff(file.path(output_folder, "alignment_score.tif"), width = 1000, height = 800, res = 150)
print(plot)
dev.off()

# Save as PNG (lightweight alternative)
png(file.path(output_folder, "alignment_score.png"), width = 1000, height = 800, res = 150)
print(plot)
dev.off()

