# ---------------------------
# Handle Arguments
# ---------------------------
# Uncomment below if running from command line
# args <- commandArgs(trailingOnly = TRUE)

# Hardcoded for testing
args <- c("MMUPHin,CQR,RUV,MetaDICT,SVD,PN,FAbatch,ComBatSeq", "output/example")

if (length(args) < 2) {
  stop("Usage: Rscript normalize_all_methods.R <method_list> <output_folder>")
}

method_list <- unlist(strsplit(args[1], ","))
output_folder <- args[2]

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# ---------------------------
# Load Libraries
# ---------------------------
library(MBECS)
library(ggplot2)
library(gridExtra)
library(phyloseq)

# ---------------------------
# Load Data
# ---------------------------
custom_matrix_path <- file.path(output_folder, "raw.csv")
custom_metadata_path <- file.path(output_folder, "metadata.csv")
default_matrix_path <- file.path("assets", "raw.csv")
default_metadata_path <- file.path("assets", "metadata.csv")

if (file.exists(custom_matrix_path) && file.exists(custom_metadata_path)) {
  cat("✅ Using uploaded user files\n")
  taxa_mat <- read.csv(custom_matrix_path, check.names = FALSE)
  metadata <- read.csv(custom_metadata_path)
} else {
  cat("⚠️ No uploaded files found — using default assets\n")
  taxa_mat <- read.csv(default_matrix_path, check.names = FALSE)
  metadata <- read.csv(default_metadata_path)
  write.csv(taxa_mat, file = file.path(output_folder, "raw.csv"), row.names = FALSE)
  write.csv(metadata, file = file.path(output_folder, "metadata.csv"), row.names = FALSE)
}

# ---------------------------
# Format Data
# ---------------------------
# Rename for MBECS compatibility
colnames(metadata)[colnames(metadata) == "phenotype"] <- "group"
colnames(metadata)[colnames(metadata) == "batchid"] <- "batch"

# Convert group and batch to factors
metadata$group <- as.factor(metadata$group)
metadata$batch <- as.factor(metadata$batch)

rownames(taxa_mat) <- metadata$sample_id
rownames(metadata) <- metadata$sample_id

# ---------------------------
# Create phyloseq object (required for MBECS v1.12.0)
# ---------------------------
otu <- otu_table(as.matrix(taxa_mat), taxa_are_rows = FALSE)
smp <- sample_data(metadata)
physeq <- phyloseq(otu, smp)

# ---------------------------
# Run MBECS Batch Effect Diagnostics
# ---------------------------
mbec.obj <- mbecProcessInput(physeq)

mbec.obj <- mbecTransform(mbec.obj, method = "clr")

mbec.obj <- mbecReportPrelim(mbec.obj)

# ---------------------------
# Visualize Pre-Correction Diagnostics
# ---------------------------
pdf(file.path(output_folder, "MBECS_Preliminary_Report.pdf"), width = 14, height = 6)
gridExtra::grid.arrange(
  mbec.obj@report.pre$pca,
  mbec.obj@report.pre$rle,
  mbec.obj@report.pre$boxplot,
  ncol = 3
)
dev.off()

# Also print to console
print(mbec.obj@report.pre$stats)

# Save individual plots (optional)
ggsave(file.path(output_folder, "pca_plot.png"), mbec.obj@report.pre$pca)
ggsave(file.path(output_folder, "rle_plot.png"), mbec.obj@report.pre$rle)
ggsave(file.path(output_folder, "boxplot.png"), mbec.obj@report.pre$boxplot)

# Save report stats
write.csv(mbec.obj@report.pre$stats, file = file.path(output_folder, "mbecs_stats.csv"))
