library(coda.base)
library(quantreg)
library(doParallel)
library(foreach)

cat("Starting CQR processing...\n")

# Parameters
pseudocount <- 1e-6

# Prepare matrices
otu_data <- as.matrix(otu_final + pseudocount)
otu_data <- log(otu_data)
otu_clr <- clr(otu_data)
X <- as.matrix(meta_datasetselected[, indep_vars, drop = FALSE])

# Parallel setup
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat(paste0("Using ", n_cores, " cores.\n"))

# Track time
start_time <- Sys.time()

corrected <- foreach(j = 1:ncol(otu_clr), .combine = 'cbind', .packages = "quantreg") %dopar% {
  print(j)
  y <- otu_clr[, j]
  fit <- rq(y ~ X, tau = seq(0.05, 0.95, by = 0.05))
  y_hat <- rowMeans(predict(fit, newdata = data.frame(X)))
  cat(sprintf("Finished column %d\n", j))
  y - y_hat
}

end_time <- Sys.time()
stopCluster(cl)

# Inverse CLR
corrected_exp <- exp(corrected)
corrected_exp <- corrected_exp / rowSums(corrected_exp)
otu_final <- as.data.frame(corrected_exp)
colnames(otu_final) <- colnames(otu_clr)
rownames(otu_final) <- rownames(otu_clr)

cat("CQR batch effect correction completed.\n")
