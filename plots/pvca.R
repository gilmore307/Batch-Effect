# ================= PVCA (Principal Variance Component Analysis) — mbec-style, pRDA table style =================
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(gridExtra)  # tableGrob + arrangeGrob
  library(grid)       # rectGrob, segmentsGrob
  library(gtable)     # add borders/lines to table grob
})

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Package 'lme4' is required. install.packages('lme4')")
}
library(lme4)

# --------- Helpers ---------
pvca_zero <- function() {
  tibble(
    Component = factor(c("Treatment","Intersection","Batch","Residuals"),
                       levels = c("Treatment","Intersection","Batch","Residuals")),
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
    dat_j$.inter <- interaction(dat_j$.batch, dat_j$.treat, drop = TRUE)  # internal name can stay ".inter"
    
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
  parts <- c(Treatment = w_treat, Intersection = w_inter, Batch = w_batch, Residuals = w_res)
  parts[!is.finite(parts)] <- 0
  nm <- names(parts)
  parts <- pmax(parts, 0)  # pmax can drop names
  names(parts) <- nm
  
  s <- sum(parts)
  if (is.finite(s) && s > 0) {
    parts <- parts / s
  } else {
    parts <- c(Treatment = 0, Intersection = 0, Batch = 0, Residuals = 1)
  }
  
  Fraction <- as.numeric(parts[c("Treatment","Intersection","Batch","Residuals")])
  Fraction[!is.finite(Fraction)] <- 0
  Fraction <- pmin(pmax(Fraction, 0), 1)
  
  tibble(
    Component = factor(c("Treatment","Intersection","Batch","Residuals"),
                       levels = c("Treatment","Intersection","Batch","Residuals")),
    Fraction  = Fraction
  )
}

# --------- Use it (IO) ---------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

metadata <- readr::read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(sample_id = as.character(sample_id))

# matrices: Raw + normalized_*
file_paths <- list.files(
  output_folder,
  pattern = "^normalized_.*_clr\\.csv$",
  full.names = TRUE
)

# Clean names (remove prefix + suffix)
file_names <- gsub("^normalized_|_clr\\.csv$", "", basename(file_paths))
file_list  <- setNames(file_paths, file_names)

# Prepend raw.csv as "Before Correction" if present
raw_fp <- file.path(output_folder, "raw.csv")
if (file.exists(raw_fp)) {
  file_list <- c("Before Correction" = raw_fp, file_list)
}

if (length(file_list) == 0) {
  stop("No CLR files found. Expected files matching 'normalized_*_clr.csv' in ", output_folder)
}

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

# --------- Plot PVCA in pRDA style: no on-bar labels, table below ----------
component_order <- c("Treatment","Intersection","Batch","Residuals")
pvca_plot_df <- pvca_df %>%
  dplyr::mutate(
    Method    = factor(Method, levels = names(file_list)),
    Component = factor(as.character(Component), levels = component_order),
    Fraction  = ifelse(is.finite(Fraction), Fraction, 0)
  ) %>%
  dplyr::mutate(Fraction = pmin(pmax(Fraction, 0), 1)) %>%
  dplyr::arrange(Method, Component)

# sanity check: 4 rows per method
stopifnot(all(dplyr::count(pvca_plot_df, Method)$n == 4))

cols <- c(
  "Residuals"    = "#1F77B4",
  "Batch"        = "#FF7F0E",
  "Intersection" = "#FFD54F",
  "Treatment"    = "#BDBDBD"
)

p <- ggplot(pvca_plot_df, aes(x = Method, y = Fraction, fill = Component)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.4) +
  scale_fill_manual(values = cols,
                    breaks = c("Residuals","Batch","Intersection","Treatment"),
                    name = "Variation sources") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05),             # 105% headroom
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

# --------- Build the styled values table (percentages) ----------
tbl <- pvca_plot_df %>%
  dplyr::select(Method, Component, Fraction) %>%
  dplyr::mutate(`%` = scales::percent(Fraction, accuracy = 1)) %>%
  dplyr::select(-Fraction) %>%
  tidyr::pivot_wider(names_from = Component, values_from = `%`) %>%
  dplyr::arrange(Method) %>%
  dplyr::select(Method, Treatment, Intersection, Batch, Residuals)

nr <- nrow(tbl); nc <- ncol(tbl)
stripe_rows <- rep(c("#FBFCFF", "#F7F8FC"), length.out = nr)
fill_core <- matrix(rep(stripe_rows, each = nc), nrow = nr, ncol = nc)

tbl_theme <- gridExtra::ttheme_minimal(
  core = list(
    fg_params = list(cex = 0.9),
    bg_params = list(fill = fill_core, col = "#D0D7DE"),  # light inner borders
    padding   = unit(c(6, 4), "mm")
  ),
  colhead = list(
    fg_params = list(cex = 1.0, fontface = 2),
    bg_params = list(fill = "#ECEFF4", col = "#D0D7DE"),
    padding   = unit(c(7, 4), "mm")
  )
)

tbl_grob <- gridExtra::tableGrob(tbl, rows = NULL, theme = tbl_theme)

# Add a clean outer border around the table
tbl_grob <- gtable::gtable_add_grob(
  tbl_grob,
  grobs = grid::rectGrob(gp = grid::gpar(fill = NA, col = "#9AA0A6", lwd = 1)),
  t = 1, l = 1, b = nrow(tbl_grob), r = ncol(tbl_grob), name = "outer-border"
)

# Add a horizontal rule under the header (fix y0/y1 to avoid diagonal)
header_rows <- which(tbl_grob$layout$name %in% c("colhead-fg", "colhead-bg"))
if (length(header_rows)) {
  header_bottom <- max(tbl_grob$layout$b[header_rows])
  tbl_grob <- gtable::gtable_add_grob(
    tbl_grob,
    grobs = grid::segmentsGrob(
      x0 = unit(0, "npc"), x1 = unit(1, "npc"),
      y0 = unit(0, "npc"), y1 = unit(0, "npc"),
      gp = grid::gpar(col = "#9AA0A6", lwd = 1.2, lineend = "butt")
    ),
    t = header_bottom, l = 1, r = ncol(tbl_grob), name = "header-sep"
  )
}

# --------- Stack plot over table and save ----------
combined <- gridExtra::arrangeGrob(p, tbl_grob, ncol = 1, heights = c(3, 1.35))

ggsave(file.path(output_folder, "PVCA.png"),
       plot = combined, width = 7.2, height = 6.8, dpi = 300)

ggsave(file.path(output_folder, "PVCA.tif"),
       plot = combined, width = 7.2, height = 6.8, dpi = 300, compression = "lzw")

# --------- Rank the methods based on Treatment (desc) then Batch (asc) ----------
rank_pvca_methods <- function(df_long) {
  df_long %>%
    group_by(Method) %>%
    summarise(
      Treatment = sum(Fraction[Component == "Treatment"], na.rm = TRUE),
      Batch     = sum(Fraction[Component == "Batch"],     na.rm = TRUE)
    ) %>%
    arrange(desc(Treatment), Batch) %>%
    mutate(Rank = row_number())
}

ranked_pvca_methods <- rank_pvca_methods(pvca_plot_df)
print(ranked_pvca_methods)
readr::write_csv(ranked_pvca_methods, file.path(output_folder, "pvca_ranking.csv"))
