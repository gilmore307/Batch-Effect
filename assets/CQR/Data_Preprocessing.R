###############################################
# Data Preprocessing (minimal, batch-focused)
###############################################

library(dplyr)
library(readxl)
library(tibble)

# --- Settings ---
criteria <- 0.005

# (1) Remove OTUs below frequency threshold
OTU_table <- samplebyotu_countmatrix
OTU_sum <- colSums(OTU_table[ , -1])
OTU_freq <- OTU_sum / sum(OTU_sum)
filtered_OTU <- which(OTU_freq > criteria * 1e-6)
samplebyotu_countmatrix <- OTU_table[, c(1, filtered_OTU + 1)]  # +1 for ID column

# (2) Remove zero samples and OTUs
otu_data <- samplebyotu_countmatrix
rownames(otu_data) <- otu_data$ID
otu_data <- otu_data[, -1]
otu_filtered <- rm_zero(otu_data)

# (3) Remove OTUs that appear in only one batch (confounded)
batch_info <- clinicalinfo_total
error_taxa <- Filter(function(taxon) {
  length(unique(batch_info$Batch[otu_filtered[, taxon] > 0])) == 1
}, colnames(otu_filtered))

otu_filtered <- otu_filtered[, !(colnames(otu_filtered) %in% error_taxa)]

# (4) Remove constant OTUs (all 0 or all same value)
otu_filtered <- otu_filtered[, apply(otu_filtered, 2, function(x) !(all(x == 0) || all(x == max(x))))]

# Final metadata match
meta_datasetselected <- count_matrix_covar %>%
  filter(ID %in% rownames(otu_filtered))

# Final OTU table
otu_final <- otu_filtered
