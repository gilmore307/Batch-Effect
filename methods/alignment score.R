# ==== Load Libraries ====
library(FNN)
library(readr)
library(dplyr)
library(ggplot2)

# ==== Read Metadata ====
metadata <- read_csv("metadata.csv")
metadata$sample_id <- as.character(metadata$sample_id)

# ==== File Paths ====
file_list <- list(
  Raw = "raw.csv",
  ALRA = "normalized_alra.csv",
  ComBat = "normalized_combat.csv",
  ConQuR = "normalized_conqur.csv",
  FSQN = "normalized_fsqn.csv",
  PLSDA = "normalized_plsda.csv"
)

# ==== Alignment Score Function ====
compute_alignment_score <- function(data, batch_vector, k = 10) {
  # Run PCA
  pca <- prcomp(data, scale. = TRUE)
  coords <- pca$x[, 1:10]  # use first 10 PCs
  
  # Find k-nearest neighbors
  knn_result <- get.knn(coords, k = k)
  neighbors <- knn_result$nn.index
  
  # Compute alignment score
  alignment_scores <- sapply(1:nrow(data), function(i) {
    neighbor_batches <- batch_vector[neighbors[i, ]]
    mean(neighbor_batches != batch_vector[i])  # fraction of neighbors from other batches
  })
  
  return(mean(alignment_scores))
}

# ==== Compute Scores for All Methods ====
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

# ==== Plot Alignment Scores ====
ggplot(alignment_scores, aes(x = reorder(Method, Score), y = Score, fill = Method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "Comparison of alignment scores",
    x = "Method",
    y = "Alignment Score"
  ) +
  ylim(0, 1) +
  theme_minimal(base_size = 14)
