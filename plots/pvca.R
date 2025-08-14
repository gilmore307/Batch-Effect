# ================= PVCA (Principal Variance Component Analysis) — mbec-style, pRDA label style =================
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Package 'lme4' is required. install.packages('lme4')")
}
library(lme4)

# --------- Helpers ---------
pvca_zero <- function() {
  tibble(
    Component = factor(c("Treatment","Interaction","Batch","Residuals"),
                       levels = c("Treatment","Interaction","Batch","Residuals")),
    Fraction  = c(0,0,0,1)
  )
}

# --------- Core PVCA (mbec approach, robust) ---------
# Steps:
# 1) Center features (no scaling), compute sample–sample correlation.
# 2) Eigen-decompose; choose PCs: min 3, max 10, reaching cumvar_threshold.
# 3) For each selected PC: LMM random effects (batch, treat, interaction).
# 4) Standardize per-PC variances, weight by PC eigenvalue share, renormalize.
compute_pvca <- function(df, meta, batch_col = "batchid", treat_col = "phenotype",
                         cumvar_threshold = 0.60, scale_features = TRUE,  # kept for API compatibility
                         na_action = stats::na.omit, quiet = TRUE) {
  
  # align by sample_id
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(meta) && "sample_id" %in% names(meta)) {
      df$sample_id <- meta$sample_id
    } else {
      return(pvca_zero())
    }
  }
  df   <- dplyr::mutate(df,   sample_id = as.character(sample_id))
  meta <- dplyr::mutate(meta, sample_id = as.character(sample_id))
  
  dfx <- dplyr::inner_join(df, meta, by = "sample_id")
  if (!nrow(dfx)) return(pvca_zero())
  if (!(batch_col %in% names(dfx)) || !(treat_col %in% names(dfx))) return(pvca_zero())
  
  # numeric features (samples x features)
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfx %>%
    dplyr::select(dplyr::all_of(feat_cols)) %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix()
  if (!ncol(X)) return(pvca_zero())
  
  # keep finite, non-constant features
  keep <- apply(X, 2, function(z) all(is.finite(z)) && stats::sd(z, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) return(pvca_zero())
  
  # names & order
  s.names <- dfx$sample_id
  rownames(X) <- s.names
  
  # center-by-feature (no scaling), sample–sample correlation
  Xc <- sweep(X, 2, colMeans(X, na.rm = TRUE), "-")
  Xt <- t(Xc)
  if (ncol(Xt) < 3) return(pvca_zero())
  
  tmp.cor <- stats::cor(Xt, use = "pairwise.complete.obs")
  
  # eigen + PC selection (3..10 to hit threshold)
  eCor <- eigen(tmp.cor, symmetric = TRUE)
  eVal <- eCor$values
  eVec <- eCor$vectors
  if (any(!is.finite(eVal)) || length(eVal) < 3) return(pvca_zero())
  
  sum.eVal <- sum(eVal)
  if (!is.finite(sum.eVal) || sum.eVal <= 0) return(pvca_zero())
  prop.PCs <- eVal / sum.eVal
  
  k.hit <- which(cumsum(prop.PCs) >= cumvar_threshold)
  k.hit <- if (length(k.hit)) k.hit[1] else length(prop.PCs)
  n.samples <- nrow(tmp.cor)
  n.PCs <- min(10, n.samples, length(eVal), max(3, k.hit))
  if (n.samples < 3 || n.PCs < 1) return(pvca_zero())
  
  # per-PC LMM with random effects (batch, treat, interaction)
  rownames(eVec) <- s.names
  meta_ord <- dfx[, c("sample_id", batch_col, treat_col)]
  meta_ord <- meta_ord[match(s.names, meta_ord$sample_id), , drop = FALSE]
  meta_ord[[batch_col]] <- factor(meta_ord[[batch_col]])
  meta_ord[[treat_col]] <- factor(meta_ord[[treat_col]])
  if (nlevels(meta_ord[[batch_col]]) < 2 || nlevels(meta_ord[[treat_col]]) < 2) return(pvca_zero())
  
  scores  <- eVec[, seq_len(n.PCs), drop = FALSE]
  colnames(scores) <- paste0(".PC", seq_len(n.PCs))
  
  ctrl <- lmerControl(
    optimizer = "bobyqa",
    calc.derivs = FALSE,
    check.conv.singular = "ignore",
    check.conv.hess = "ignore",
    check.nobs.vs.nlev = "ignore",
    check.nobs.vs.rankZ = "ignore"
  )
  
  w_treat <- 0; w_batch <- 0; w_inter <- 0; w_res <- 0
  scaled.eVal <- eVal / sum.eVal
  
  for (j in seq_len(n.PCs)) {
    dat_j <- data.frame(
      .score = scores[, j],
      .batch = meta_ord[[batch_col]],
      .treat = meta_ord[[treat_col]]
    )
    dat_j$.inter <- interaction(dat_j$.batch, dat_j$.treat, drop = TRUE)
    
    # fit with interaction; if singular, drop interaction
    fit <- tryCatch(
      {
        if (quiet) suppressWarnings(suppressMessages(
          lmer(.score ~ 1 + (1|.batch) + (1|.treat) + (1|.inter),
               data = dat_j, REML = TRUE, na.action = na_action, control = ctrl)
        )) else
          lmer(.score ~ 1 + (1|.batch) + (1|.treat) + (1|.inter),
               data = dat_j, REML = TRUE, na.action = na_action, control = ctrl)
      },
      error = function(e) NULL
    )
    
    if (!is.null(fit) && isSingular(fit, tol = 1e-6)) {
      fit <- tryCatch(
        {
          if (quiet) suppressWarnings(suppressMessages(
            lmer(.score ~ 1 + (1|.batch) + (1|.treat),
                 data = dat_j, REML = TRUE, na.action = na_action, control = ctrl)
          )) else
            lmer(.score ~ 1 + (1|.batch) + (1|.treat),
                 data = dat_j, REML = TRUE, na.action = na_action, control = ctrl)
        },
        error = function(e) NULL
      )
    }
    
    pj <- scaled.eVal[j]
    if (!is.null(fit)) {
      vc <- suppressWarnings(as.data.frame(VarCorr(fit)))
      v_batch <- sum(vc$vcov[vc$grp == ".batch"],    na.rm = TRUE)
      v_treat <- sum(vc$vcov[vc$grp == ".treat"],    na.rm = TRUE)
      v_inter <- sum(vc$vcov[vc$grp == ".inter"],    na.rm = TRUE)
      v_res   <- sum(vc$vcov[vc$grp == "Residual"],  na.rm = TRUE)
      if (!is.finite(v_res)) v_res <- tryCatch(stats::sigma(fit)^2, error = function(e) NA_real_)
      
      v_tot <- v_batch + v_treat + v_inter + v_res
      if (!is.finite(v_tot) || v_tot <= 0) {
        if (is.finite(pj)) w_res <- w_res + pj
      } else {
        w_batch <- w_batch + pj * (if (is.finite(v_batch)) v_batch / v_tot else 0)
        w_treat <- w_treat + pj * (if (is.finite(v_treat)) v_treat / v_tot else 0)
        w_inter <- w_inter + pj * (if (is.finite(v_inter)) v_inter / v_tot else 0)
        w_res   <- w_res   + pj * (if (is.finite(v_res))   v_res   / v_tot else 0)
      }
    } else {
      if (is.finite(pj)) w_res <- w_res + pj
    }
  }
  
  # aggregate & renormalize — PRESERVE NAMES!
  parts <- c(Treatment = w_treat, Interaction = w_inter, Batch = w_batch, Residuals = w_res)
  parts[!is.finite(parts)] <- 0
  nm <- names(parts)
  parts <- pmax(parts, 0)  # pmax can drop names
  names(parts) <- nm
  
  s <- sum(parts)
  if (is.finite(s) && s > 0) {
    parts <- parts / s
  } else {
    parts <- c(Treatment = 0, Interaction = 0, Batch = 0, Residuals = 1)
  }
  
  Fraction <- as.numeric(parts[c("Treatment","Interaction","Batch","Residuals")])
  Fraction[!is.finite(Fraction)] <- 0
  Fraction <- pmin(pmax(Fraction, 0), 1)
  
  tibble(
    Component = factor(c("Treatment","Interaction","Batch","Residuals"),
                       levels = c("Treatment","Interaction","Batch","Residuals")),
    Fraction  = Fraction
  )
}

# --------- Use it (IO) ---------
args <- commandArgs(trailingOnly = TRUE)
output_folder <- if (length(args)) args[1] else "output/example"

metadata <- readr::read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(sample_id = as.character(sample_id))

# matrices: Raw + normalized_*
file_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
file_list  <- setNames(file_paths, gsub("^normalized_|\\.csv$", "", basename(file_paths)))
file_list  <- c("Before correction" = file.path(output_folder, "raw.csv"), file_list)

pvca_list <- lapply(names(file_list), function(nm) {
  message("Computing PVCA: ", nm)
  df <- readr::read_csv(file_list[[nm]], show_col_types = FALSE)
  out <- compute_pvca(df, metadata,
                      batch_col = "batchid",
                      treat_col = "phenotype",
                      cumvar_threshold = 0.60,
                      scale_features = TRUE,
                      quiet = TRUE)
  out$Method <- nm
  out
})

pvca_df <- dplyr::bind_rows(pvca_list) %>%
  dplyr::mutate(Method = factor(Method, levels = names(file_list)))

# --------- Plot PVCA in pRDA style: fixed on-bar labels (T/I/B/R) ---------
component_order <- c("Treatment","Interaction","Batch","Residuals")
abbr <- c(Treatment="T", Interaction="I", Batch="B", Residuals="R")
# fixed vertical label positions as fractions of bar height (tweak if desired)
label_pos_frac <- c(Treatment=0.95, Interaction=0.65, Batch=0.35, Residuals=0.08)

pvca_plot_df <- pvca_df %>%
  dplyr::mutate(
    Method    = factor(Method, levels = names(file_list)),
    Component = factor(as.character(Component), levels = component_order),
    Fraction  = ifelse(is.finite(Fraction), Fraction, 0)
  ) %>%
  dplyr::mutate(Fraction = pmin(pmax(Fraction, 0), 1)) %>%
  dplyr::arrange(Method, Component) %>%
  dplyr::mutate(
    label = paste0(abbr[as.character(Component)], " = ", round(Fraction*100)),
    y_pos = label_pos_frac[as.character(Component)]
  )

# sanity check: 4 rows per method
stopifnot(all(dplyr::count(pvca_plot_df, Method)$n == 4))

cols <- c(
  "Residuals"    = "#1F77B4",
  "Batch"        = "#FF7F0E",
  "Interaction"  = "#FFD54F",
  "Treatment"    = "#BDBDBD"
)

p <- ggplot(pvca_plot_df, aes(x = Method, y = Fraction, fill = Component)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.4) +
  scale_fill_manual(values = cols,
                    breaks = c("Residuals","Batch","Interaction","Treatment"),
                    name = "Variation sources") +
  # fixed-position labels on bars
  geom_text(aes(y = y_pos, label = label), size = 3.2, vjust = 0.5) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05),             # 105% headroom to avoid clipping labels
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Methods", y = "Explained variance (%)", title = "PVCA (weighted by PC variance)") +
  theme_bw() +
  theme(
    legend.position    = "right",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

# --------- Save ---------
ggsave(file.path(output_folder, "PVCA.png"), p, width = 7.2, height = 5.6, dpi = 300)
ggsave(file.path(output_folder, "PVCA.tif"), p, width = 7.2, height = 5.6, dpi = 300, compression = "lzw")


# Function to rank methods based on Treatment and Batch variance
rank_pvca_methods <- function(pvca_df) {
  # Compute the median fraction for Treatment and Batch components
  median_pvca <- pvca_df %>%
    group_by(Method, Component) %>%
    summarise(median_Fraction = median(Fraction, na.rm = TRUE), .groups = "drop") %>%
    spread(Component, median_Fraction)  # Wide format: one column for Treatment, one for Batch
  
  # Rank methods based on median Treatment fraction (higher is better) and Batch fraction (lower is better)
  ranked_results <- median_pvca %>%
    arrange(desc(Treatment), Batch) %>%  # Prioritize Treatment (higher) and then Batch (lower)
    mutate(Rank = row_number())  # Assign ranks
  
  return(ranked_results)
}

# Rank the methods based on PVCA results for Treatment and Batch components
ranked_pvca_methods <- rank_pvca_methods(pvca_plot_df)

# Display the ranked methods
ranked_pvca_methods
