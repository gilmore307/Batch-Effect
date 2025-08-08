###############################################
# Background functions (batch effect only)
###############################################

library(vegan)
library(compositions)

# Remove samples and OTUs with all zero counts
rm_zero <- function(otu_table) {
  otu_1 <- otu_table
  zero_col_idx <- as.integer(which(apply(otu_1, 1, sum) == 0))
  if(length(zero_col_idx) > 0) {
    otu_with_ID_ch1 <- otu_table[,-zero_col_idx]
  } else {
    otu_with_ID_ch1 <- otu_table
  }
  otu_2 <- otu_with_ID_ch1
  zero_row_idx <- as.integer(which(apply(otu_2, 1, sum) == 0))
  if(length(zero_row_idx) > 0) {
    otu_with_ID_ch2 <- otu_with_ID_ch1[-zero_row_idx,]
  } else {
    otu_with_ID_ch2 <- otu_with_ID_ch1
  }
  otu_rm_zero <- otu_with_ID_ch2
  return(otu_with_ID_ch2)
}

# Compute PERMANOVA R² for various distance metrics
PERMANOVA_R2 <- function(table, batch) {
  pseudocount = 1e-06 
  res <- list()
  dissim <- c("bray", "aitchison", "manhattan", "canberra")
  for (d in dissim) {
    if (d == "bray") {
      bray_dist <- vegdist(table, method = "bray")
      res[[d]] <- adonis2(bray_dist ~ batch)$R2[1]
    } else {
      if (d == "aitchison") {
        table = table + pseudocount
      }
      res[[d]] <- adonis2(table ~ batch, method = d, na.rm = TRUE)$R2[1]
    }
  }
  return(res)
}
