#!/usr/bin/env Rscript
# Save EXACT inputs used in the vignette analyses:
#  1) abundance.csv  = CLR(ad.filter)  (samples x OTUs)
#  2) metadata.csv   = sample_id + sequencing_run_date + initial_phenol_concentration.regroup

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(mixOmics)                 # logratio.transfo()
  library(TreeSummarizedExperiment)
  library(PLSDAbatch)               # AD_data, PreFL()
})

# -- Load dataset from PLSDAbatch
data("AD_data", package = "PLSDAbatch")
se <- AD_data$FullData

# -- Start from raw counts (as in the vignette)
ad.count <- assays(se)$Count   # expected: samples x OTUs

# -- Pre-filter (vignette: PreFL on raw counts)
ad.filter.res <- PreFL(data = ad.count)
ad.filter     <- ad.filter.res$data.filter   # still samples x OTUs

# -- CLR transform (vignette uses offset = 1 on raw-count-derived matrix)
ad.clr <- logratio.transfo(X = ad.filter, logratio = "CLR", offset = 1)
class(ad.clr) <- "matrix"      # ensure plain matrix for downstream tools

# -- Get the metadata FIELDS used in analyses (from the vignette)
#    They come from rowData() in this dataset (rows are samples).
rd <- as.data.frame(rowData(se))

required_cols <- c("sequencing_run_date", "initial_phenol_concentration.regroup")
missing <- setdiff(required_cols, colnames(rd))
if (length(missing)) {
  stop("Missing required metadata columns in rowData(): ",
       paste(missing, collapse = ", "))
}

# -- Align samples between ad.clr (rownames) and metadata (rownames)
common <- intersect(rownames(rd), rownames(ad.clr))
if (length(common) == 0) stop("No overlapping sample IDs between abundance and metadata.")
# keep same order in both
ad.clr <- ad.clr[common, , drop = FALSE]
rd     <- rd[common, , drop = FALSE]

# -- Build the minimal metadata used in analyses
metadata <- data.frame(
  sample_id = common,
  sequencing_run_date = rd$sequencing_run_date,
  initial_phenol_concentration.regroup = rd$initial_phenol_concentration.regroup,
  check.names = FALSE
)

# -- Write EXACTLY two files
write.csv(as.data.frame(ad.clr), "raw.csv", quote = FALSE)
write.csv(rd, "metadata.csv", row.names = FALSE, quote = FALSE)

cat("Wrote files:\n  - raw.csv (", nrow(ad.clr), " samples x ",
    ncol(ad.clr), " OTUs)\n  - metadata.csv (", nrow(metadata),
    " samples x ", ncol(metadata), " columns)\n", sep = "")
