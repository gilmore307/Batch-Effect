# ==== Load Libraries ====
library(FNN)
library(dplyr)
library(readr)
library(ggplot2)

# ==== Handle Argument ====
# args <- commandArgs(trailingOnly = TRUE)
args <- c("output/example")  # Replace for actual usage
output_folder <- args[1]

# ==== Read Metadata ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"))
metadata$sample_id <- as.character(metadata$sample_id)

# ==== Find Normalized Files ====
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== Entropy Function ====
compute_entropy_score <- function(data, batch_labels, k = 10) {
  # PCA first to reduce dimensions
  pca <- prcomp(data, scale. = TRUE)
  coords <- pca$x[, 1:min(nrow(data)-1, 10)]
  
  # Get k-nearest neighbors
  neighbors <- get.knn(coords, k = k)$nn.index
  
  # Compute entropy for each sample
  entropy_vec <- sapply(1:nrow(coords), function(i) {
    neighbor_batches <- batch_labels[neighbors[i, ]]
    freq <- table(neighbor_batches) / k
    -sum(freq * log2(freq + 1e-10))  # add small value to avoid log(0)
  })
  
  mean(entropy_vec)
}

# ==== Compute Entropy Scores for Each Method ====
entropy_scores <- data.frame(Method = character(), Entropy = numeric())

for (name in names(file_list)) {
  df <- read_csv(file_list[[name]])
  df$sample_id <- metadata$sample_id
  df_merged <- inner_join(df, metadata, by = "sample_id")
  
  X <- df_merged %>%
    select(where(is.numeric)) %>%
    select_if(~ sd(.) > 0) %>%
    as.matrix()
  
  batch <- as.factor(df_merged$batchid)
  entropy <- compute_entropy_score(X, batch, k = 10)
  
  entropy_scores <- rbind(entropy_scores, data.frame(Method = name, Entropy = entropy))
}

# ==== Plot ====
plot <- ggplot(entropy_scores, aes(x = reorder(Method, Entropy), y = Entropy, fill = Method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "Entropy of Batch Mixing",
    x = "Method",
    y = "Entropy"
  ) +
  ylim(0, log2(length(unique(metadata$batchid)))) +
  theme_minimal(base_size = 14)

# ==== Save Plots ====
# Save as TIFF (high-quality)
tiff(file.path(output_folder, "entropy_batch_mixing.tif"), width = 1000, height = 800, res = 150)
print(plot)
dev.off()

# Save as PNG (lightweight alternative)
png(file.path(output_folder, "entropy_batch_mixing.png"), width = 1000, height = 800, res = 150)
print(plot)
dev.off()
