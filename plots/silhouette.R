# ==== Load Required Libraries ====
library(ggplot2)
library(cluster)
library(vegan)   
library(dplyr)
library(readr)   
library(patchwork)

# ==== Read UID argument ====
args <- c("output/example")  
if (length(args) < 1) stop("Usage: Rscript silhouette_plot.R <output_folder>")
output_folder <- args[1]

# ==== Read Data ====
metadata <- read_csv(file.path(output_folder, "metadata.csv"))
metadata$sample_id <- as.character(metadata$sample_id)

# Find All Normalized Files
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))

# Add Raw file explicitly
file_list <- c(Raw = file.path(output_folder, "raw.csv"), file_list)

# ==== Clustering and Silhouette Function ====
silhouette_plot <- function(df, method_name) {
  # Perform KMeans clustering
  kmeans_result <- kmeans(df, centers = 3) 
  
  # Calculate silhouette width for each cluster
  sil <- silhouette(kmeans_result$cluster, dist(df))
  
  # Calculate the overall silhouette width
  overall_silhouette_width <- mean(sil[, 3])
  
  # Group-wise silhouette width (average per cluster)
  group_silhouette_width <- aggregate(sil[, 3], by = list(kmeans_result$cluster), FUN = mean)
  colnames(group_silhouette_width) <- c("Cluster", "Silhouette_Width")
  
  # Prepare data for plotting
  sil_df <- data.frame(
    silhouette_width = sil[, 3],
    cluster = as.factor(sil[, 1]),
    sample_id = rownames(df)
  )
  
  # Create the plot
  p <- ggplot(sil_df, aes(x = silhouette_width, fill = cluster)) +
    geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.6) +
    labs(
      title = paste("Silhouette Plot -", method_name, 
                    "\n(Overall Silhouette: ", round(overall_silhouette_width, 3), ")"),
      x = "Silhouette Coefficient",
      y = "Counts",
      fill = "Cluster"  # Set legend title
    ) +
    scale_fill_manual(values = c("red", "green", "blue")) +  # Adjust colors as needed
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust the axis labels
  
  return(p)
}

# ==== Generate All Silhouette Plots ====
plots <- lapply(names(file_list), function(name) {
  df <- read_csv(file_list[[name]])
  df_num <- df %>%
    select(where(is.numeric)) %>%
    select_if(~ sd(.) > 0)  # Remove constant columns
  
  silhouette_plot(df_num, name)
})

# ==== Combine and Save with Square Subplots ====
ncol_grid <- 3
nrow_grid <- ceiling(length(plots) / ncol_grid)
subplot_size <- 6  # inches per subplot (square)

combined_plot <- wrap_plots(plots, ncol = ncol_grid)

ggsave(file.path(output_folder, "silhouette_plot.tif"), combined_plot,
       width = subplot_size * ncol_grid,
       height = subplot_size * nrow_grid,
       dpi = 300)

ggsave(file.path(output_folder, "silhouette_plot.png"), combined_plot,
       width = subplot_size * ncol_grid,
       height = subplot_size * nrow_grid,
       dpi = 300)
