# ================= pRDA variance partition in two geometries =================
# A) Aitchison (CLR + RDA)    -> uses files: normalized_*_clr.csv  (+ raw.csv if present)
# B) Bray-Curtis (TSS + dbRDA)-> uses files: normalized_*_tss.csv  (+ raw.csv if present)
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(gridExtra)   # tableGrob + arrangeGrob
  library(grid)        # grobs

# Map method codes to short labels for figures
method_short_label <- function(x) {
  map <- c(
    qn = "QN", bmc = "BMC", limma = "Limma", conqur = "ConQuR",
    plsda = "PLSDA-batch", combat = "ComBat", fsqn = "FSQN", mmuphin = "MMUPHin",
    ruv = "RUV-III-NB", metadict = "MetaDICT", svd = "SVD", pn = "PN",
    fabatch = "FAbatch", combatseq = "ComBat-seq", debias = "DEBIAS-M"
  )
  sapply(x, function(v){ lv <- tolower(v); if (lv %in% names(map)) map[[lv]] else v })
}

  library(gtable)      # table tweaks
  library(vegan)       # rda, capscale, RsquareAdj
})

# --------- Args / IO ---------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- "output/example"  # default folder for quick runs
}
output_folder <- args[1]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

metadata <- read_csv(file.path(output_folder, "metadata.csv"), show_col_types = FALSE) |>
  mutate(sample_id = as.character(sample_id))
}
if (!("batch_id" %in% names(metadata)) && ("batch_id" %in% names(metadata))) {
  metadata$batch_id <- metadata$batch_id
}

# ---- Find normalized files ----
clr_paths <- list.files(output_folder, pattern = "^normalized_.*_clr\\.csv$", full.names = TRUE)
tss_paths <- list.files(output_folder, pattern = "^normalized_.*_tss\\.csv$", full.names = TRUE)

# Fallback: if no suffix-specific outputs, use any normalized_*.csv for both
if (!length(clr_paths) && !length(tss_paths)) {
  any_paths <- list.files(output_folder, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
  clr_paths <- any_paths
  tss_paths <- any_paths
}

name_from <- function(paths, suffix) gsub(paste0("^normalized_|_", suffix, "\\.csv$"), "", basename(paths))
file_list_clr <- setNames(clr_paths, method_short_label(name_from(clr_paths, "clr")))
file_list_tss <- setNames(tss_paths, method_short_label(name_from(tss_paths, "tss")))

# Include raw_clr.csv / raw_tss.csv as "Before correction" if present
raw_clr_fp <- file.path(output_folder, "raw_clr.csv")
raw_tss_fp <- file.path(output_folder, "raw_tss.csv")
if (file.exists(raw_clr_fp)) file_list_clr <- c("Before correction" = raw_clr_fp, file_list_clr)
if (file.exists(raw_tss_fp)) file_list_tss <- c("Before correction" = raw_tss_fp, file_list_tss)

if (!length(file_list_clr) && !length(file_list_tss)) {
  stop("No normalized files found (expected raw_clr.csv/raw_tss.csv and/or normalized_*_clr.csv / normalized_*_tss.csv) in ", output_folder)
}

# --------- Helpers (CLR + safe R^2) ---------
safe_closure <- function(X) {
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0
  rs <- rowSums(X, na.rm = TRUE)
  bad <- which(rs == 0 | !is.finite(rs))
  if (length(bad)) { X[bad, ] <- 1 / ncol(X); rs <- rowSums(X, na.rm = TRUE) }
  sweep(X, 1, rs, "/")
}
clr_transform <- function(X) {
  Xc <- safe_closure(X)
  for (i in seq_len(nrow(Xc))) {
    xi <- Xc[i, ]; pos <- xi > 0 & is.finite(xi)
    if (!any(pos)) { xi[] <- 1/length(xi); pos <- xi > 0 }
    if (any(!pos)) {
      m <- min(xi[pos], na.rm = TRUE)
      xi[!pos] <- min(m*0.5, 1e-8)
      xi <- xi / sum(xi)
    }
    Xc[i, ] <- xi
  }
  L <- log(Xc)
  sweep(L, 1, rowMeans(L), "-")
}
safe_adjR2 <- function(fit) {
  out <- tryCatch(vegan::RsquareAdj(fit), error = function(e) list(r.squared=NA_real_, adj.r.squared=NA_real_))
  adj <- suppressWarnings(out$adj.r.squared)
  if (is.finite(adj)) return(max(0, adj))
  r <- suppressWarnings(out$r.squared)
  if (is.finite(r)) return(max(0, r))
  0
}

# --------- Core calculators (return named numeric: Treatment, Intersection, Batch, Residuals) ---------
compute_prda_parts_aitch <- function(df, meta, batch_col = "batch_id", treat_col = "phenotype") {
  # ensure sample_id present & align
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(meta)) df$sample_id <- meta$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfx <- inner_join(df, meta, by = "sample_id") %>%
    filter(!is.na(.data[[batch_col]]), !is.na(.data[[treat_col]]))
  
  # features
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfx %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  keep <- apply(X, 2, function(z) all(is.finite(z)) && sd(z, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) return(c(Treatment=0, Intersection=0, Batch=0, Residuals=1))
  
  # Aitchison response (CLR; if already log-ratio, row-center)
  Y <- if (any(X < 0, na.rm = TRUE)) sweep(X, 1, rowMeans(X, na.rm = TRUE), "-") else clr_transform(X)
  
  dfx[[batch_col]] <- factor(dfx[[batch_col]])
  dfx[[treat_col]] <- factor(dfx[[treat_col]])
  if (nlevels(dfx[[treat_col]]) < 2 || nlevels(dfx[[batch_col]]) < 2) {
    return(c(Treatment=0, Intersection=0, Batch=0, Residuals=1))
  }
  
  Ymat <- Y
  f_both   <- as.formula(paste("Ymat ~", treat_col, "+", batch_col))
  f_t_pure <- as.formula(paste("Ymat ~", treat_col, "+ Condition(", batch_col, ")"))
  f_b_pure <- as.formula(paste("Ymat ~", batch_col, "+ Condition(", treat_col, ")"))
  
  fit_both  <- vegan::rda(f_both,   data = dfx)
  fit_t     <- vegan::rda(f_t_pure, data = dfx)
  fit_b     <- vegan::rda(f_b_pure, data = dfx)
  
  r2_both   <- safe_adjR2(fit_both)
  r2_t_pure <- safe_adjR2(fit_t)
  r2_b_pure <- safe_adjR2(fit_b)
  
  r2_inter  <- max(0, r2_both - r2_t_pure - r2_b_pure)
  r2_resid  <- max(0, 1 - (r2_t_pure + r2_b_pure + r2_inter))
  
  parts <- c(Treatment = r2_t_pure, Intersection = r2_inter,
             Batch = r2_b_pure, Residuals = r2_resid)
  parts[!is.finite(parts)] <- 0
  parts <- pmax(0, parts)
  s <- sum(parts)
  if (s > 0) parts <- parts / s
  parts
}

compute_prda_parts_bray <- function(df, meta, batch_col = "batch_id", treat_col = "phenotype") {
  # ensure sample_id present & align
  if (!"sample_id" %in% names(df)) {
    if (nrow(df) == nrow(meta)) df$sample_id <- meta$sample_id
    else stop("Input lacks 'sample_id' and row count != metadata; can't align samples.")
  }
  df  <- df %>% mutate(sample_id = as.character(sample_id))
  dfx <- inner_join(df, meta, by = "sample_id") %>%
    filter(!is.na(.data[[batch_col]]), !is.na(.data[[treat_col]]))
  
  # features (TSS proportions for Bray-Curtis)
  feat_cols <- setdiff(names(df), "sample_id")
  X <- dfx %>% select(all_of(feat_cols)) %>% select(where(is.numeric)) %>% as.matrix()
  keep <- apply(X, 2, function(z) all(is.finite(z)) && sd(z, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) return(c(Treatment=0, Intersection=0, Batch=0, Residuals=1))
  
  Xtss <- safe_closure(X)
  
  dfx[[batch_col]] <- factor(dfx[[batch_col]])
  dfx[[treat_col]] <- factor(dfx[[treat_col]])
  if (nlevels(dfx[[treat_col]]) < 2 || nlevels(dfx[[batch_col]]) < 2) {
    return(c(Treatment=0, Intersection=0, Batch=0, Residuals=1))
  }
  
  # distance-based RDA (dbRDA) via capscale with Bray-Curtis
  f_both   <- as.formula(paste("Xtss ~", treat_col, "+", batch_col))
  f_t_pure <- as.formula(paste("Xtss ~", treat_col, "+ Condition(", batch_col, ")"))
  f_b_pure <- as.formula(paste("Xtss ~", batch_col, "+ Condition(", treat_col, ")"))
  
  fit_both  <- vegan::capscale(f_both,   data = dfx, distance = "bray")
  fit_t     <- vegan::capscale(f_t_pure, data = dfx, distance = "bray")
  fit_b     <- vegan::capscale(f_b_pure, data = dfx, distance = "bray")
  
  r2_both   <- safe_adjR2(fit_both)
  r2_t_pure <- safe_adjR2(fit_t)
  r2_b_pure <- safe_adjR2(fit_b)
  
  r2_inter  <- max(0, r2_both - r2_t_pure - r2_b_pure)
  r2_resid  <- max(0, 1 - (r2_t_pure + r2_b_pure + r2_inter))
  
  parts <- c(Treatment = r2_t_pure, Intersection = r2_inter,
             Batch = r2_b_pure, Residuals = r2_resid)
  parts[!is.finite(parts)] <- 0
  parts <- pmax(0, parts)
  s <- sum(parts)
  if (s > 0) parts <- parts / s
  parts
}

# --------- Plot + styled table (reusable) ---------
plot_prda_with_table <- function(parts_df, file_list, title_prefix, outfile_prefix) {
  component_order <- c("Treatment","Intersection","Batch","Residuals")
  stopifnot(all(dplyr::count(parts_df, Method)$n == 4))
  
  parts_df <- parts_df %>%
    arrange(Method) %>%
    group_by(Method) %>%
    mutate(Component = factor(component_order[row_number()], levels = component_order)) %>%
    ungroup() %>%
    mutate(
      Method   = factor(Method, levels = names(file_list)),
      Fraction = pmin(pmax(Fraction, 0), 1)
    )
  
  cols <- c(
    "Residuals"    = "#1F77B4",
    "Batch"        = "#FF7F0E",
    "Intersection" = "#FFD54F",
    "Treatment"    = "#BDBDBD"
  )
  
  p <- ggplot(parts_df, aes(x = Method, y = Fraction, fill = Component)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.4) +
    scale_fill_manual(
      values = cols,
      breaks = component_order,   # c("Treatment","Intersection","Batch","Residuals")
      limits = component_order,   # lock the scale order
      name   = "Variation sources"
    )+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1.05), expand = expansion(mult = c(0, 0))) +
    labs(x = "Methods", y = "Explained variance (%)",
         title = title_prefix) +
    theme_bw() +
    theme(
      legend.position    = "right",
      legend.title       = element_text(size = 10),
      legend.text        = element_text(size = 9),
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(hjust = 0.5, face = "bold")
    )
  
  # Values table
  tbl <- parts_df %>%
    dplyr::select(Method, Component, Fraction) %>%
    dplyr::mutate(`%` = scales::percent(Fraction, accuracy = 1)) %>%
    dplyr::select(-Fraction) %>%
    tidyr::pivot_wider(names_from = Component, values_from = `%`) %>%
    dplyr::arrange(Method)
  
  nr <- nrow(tbl); nc <- ncol(tbl)
  stripe_rows <- rep(c("#FBFCFF", "#F7F8FC"), length.out = nr)
  fill_core <- matrix(rep(stripe_rows, each = nc), nrow = nr, ncol = nc)
  
  tbl_theme <- gridExtra::ttheme_minimal(
    core = list(
      fg_params = list(cex = 0.9),
      bg_params = list(fill = fill_core, col = "#D0D7DE"),
      padding   = grid::unit(c(6, 4), "mm")
    ),
    colhead = list(
      fg_params = list(cex = 1.0, fontface = 2),
      bg_params = list(fill = "#ECEFF4", col = "#D0D7DE"),
      padding   = grid::unit(c(7, 4), "mm")
    )
  )
  
  # Just build the table — no header rule (prevents the diagonal slash)
  tbl_grob <- gridExtra::tableGrob(tbl, rows = NULL, theme = tbl_theme)
  
  header_rows <- grep("^colhead", tbl_grob$layout$name)
  if (length(header_rows)) {
    header_bottom <- max(tbl_grob$layout$b[header_rows])
    tbl_grob <- gtable::gtable_add_grob(
      tbl_grob,
      grobs = grid::segmentsGrob(
        x0 = grid::unit(0, "npc"), x1 = grid::unit(1, "npc"),
        y0 = grid::unit(1, "npc"), y1 = grid::unit(1, "npc"),
        gp = grid::gpar(col = "#9AA0A6", lwd = 1.2)
      ),
      t = header_bottom, l = 1, r = ncol(tbl_grob), name = "header-sep"
    )
  }
  
  
  # Add a rule under the header
  header_rows <- which(tbl_grob$layout$name %in% c("colhead-fg", "colhead-bg"))
  if (length(header_rows)) {
    header_bottom <- max(tbl_grob$layout$b[header_rows])
    tbl_grob <- gtable::gtable_add_grob(
      tbl_grob,
      grobs = grid::segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                                 gp = grid::gpar(col = "#9AA0A6", lwd = 1.2)),
      t = header_bottom, l = 1, r = ncol(tbl_grob), name = "header-sep"
    )
  }
  
  combined <- gridExtra::arrangeGrob(p, tbl_grob, ncol = 1, heights = c(3, 1.35))
  
  ggsave(file.path(output_folder, paste0(outfile_prefix, ".png")),
         plot = combined, width = 7.6, height = 6.9, dpi = 300)
  ggsave(file.path(output_folder, paste0(outfile_prefix, ".tif")),
         plot = combined, width = 7.6, height = 6.9, dpi = 300, compression = "lzw")
  
  invisible(list(plot = p, table = tbl))
}

# --------- Compute & plot: A) Aitchison (CLR + RDA) ---------
if (length(file_list_clr)) {
  batch_col <- "batch_id"; treat_col <- "phenotype"
  if (!("phenotype" %in% names(metadata))) stop("metadata.csv has no 'phenotype'.")
  if (dplyr::n_distinct(metadata$phenotype) < 2) stop("'phenotype' needs >= 2 levels.")
  
  parts_df_aitch <- lapply(names(file_list_clr), function(nm) {
    message("pRDA (Aitchison) for: ", nm)
    df <- read_csv(file_list_clr[[nm]], show_col_types = FALSE)
    v  <- compute_prda_parts_aitch(df, metadata, batch_col, treat_col)
    tibble::tibble(Method = nm, Fraction = as.numeric(v))
  }) |> bind_rows()
  
  plot_prda_with_table(
    parts_df_aitch, file_list_clr,
    title_prefix  = "pRDA variance partition — Aitchison (CLR + RDA)",
    outfile_prefix = "pRDA_aitchison"
  )
}

# --------- Compute & plot: B) Bray-Curtis (TSS + dbRDA) ---------
if (length(file_list_tss)) {
  batch_col <- "batch_id"; treat_col <- "phenotype"
  if (!("phenotype" %in% names(metadata))) stop("metadata.csv has no 'phenotype'.")
  if (dplyr::n_distinct(metadata$phenotype) < 2) stop("'phenotype' needs ≥ 2 levels.")
  
  parts_df_bray <- lapply(names(file_list_tss), function(nm) {
    message("pRDA (Bray-Curtis) for: ", nm)
    df <- read_csv(file_list_tss[[nm]], show_col_types = FALSE)
    v  <- compute_prda_parts_bray(df, metadata, batch_col, treat_col)
    tibble::tibble(Method = nm, Fraction = as.numeric(v))
  }) |> bind_rows()
  
  plot_prda_with_table(
    parts_df_bray, file_list_tss,
    title_prefix  = "pRDA variance partition — Bray-Curtis (TSS + dbRDA)",
    outfile_prefix = "pRDA_braycurtis"
  )
}

# ==== Unified pRDA scoring (Aitchison + Bray-Curtis) =========================
# Per-geometry score (higher = better):
#   Score = Treatment / (Treatment + Batch + 1e-12)
# Unified score = geometric mean of available per-geometry scores.

component_order <- c("Treatment","Intersection","Batch","Residuals")

ensure_components <- function(parts_df) {
  if (!"Component" %in% names(parts_df)) {
    parts_df <- parts_df %>%
      group_by(Method) %>%
      mutate(Component = factor(component_order[row_number()],
                                levels = component_order)) %>%
      ungroup()
  } else {
    parts_df <- parts_df %>%
      mutate(Component = factor(as.character(Component),
                                levels = component_order))
  }
  parts_df
}

extract_TB <- function(parts_df) {
  parts_df %>%
    ensure_components() %>%
    group_by(Method) %>%
    summarise(
      Treatment = sum(Fraction[Component == "Treatment"], na.rm = TRUE),
      Batch     = sum(Fraction[Component == "Batch"],     na.rm = TRUE),
      .groups   = "drop"
    )
}

score_tb <- function(tb) {
  eps <- 1e-12
  tb %>%
    mutate(Score = Treatment / (Treatment + Batch + eps))
}

tb_ait <- if (exists("parts_df_aitch")) score_tb(extract_TB(parts_df_aitch)) %>%
  rename(Treatment_Aitchison = Treatment,
         Batch_Aitchison     = Batch,
         Score_Aitchison     = Score) else NULL

tb_bra <- if (exists("parts_df_bray")) score_tb(extract_TB(parts_df_bray)) %>%
  rename(Treatment_Bray = Treatment,
         Batch_Bray     = Batch,
         Score_Bray     = Score) else NULL

# Determine if we only have the baseline ("Before correction")
avail_methods <- unique(c(
  if (exists("parts_df_aitch")) unique(parts_df_aitch$Method) else character(),
  if (exists("parts_df_bray"))  unique(parts_df_bray$Method)  else character()
))
only_baseline <- length(avail_methods) == 1L && identical(avail_methods, "Before correction")

if (only_baseline) {
  # ---- Baseline-only assessment (no ranking) ----
  assess_rows <- list()
  
  if (exists("parts_df_aitch")) {
    pf <- ensure_components(parts_df_aitch) %>% filter(Method == "Before correction")
    tr <- sum(pf$Fraction[pf$Component == "Treatment"], na.rm = TRUE)
    bt <- sum(pf$Fraction[pf$Component == "Batch"],     na.rm = TRUE)
    it <- sum(pf$Fraction[pf$Component == "Intersection"], na.rm = TRUE)
    rs <- sum(pf$Fraction[pf$Component == "Residuals"],    na.rm = TRUE)
    needs <- (is.finite(bt) && is.finite(tr) && (bt >= tr)) || (bt > 0.05)
    assess_rows[["Aitchison (CLR)"]] <- tibble::tibble(
      Geometry = "Aitchison (CLR)",
      Treatment = tr, Batch = bt, Intersection = it, Residuals = rs,
      Score = if (is.finite(tr) && is.finite(bt)) tr / (tr + bt + 1e-12) else NA_real_,
      Needs_Correction = needs
    )
  }
  
  if (exists("parts_df_bray")) {
    pf <- ensure_components(parts_df_bray) %>% filter(Method == "Before correction")
    tr <- sum(pf$Fraction[pf$Component == "Treatment"], na.rm = TRUE)
    bt <- sum(pf$Fraction[pf$Component == "Batch"],     na.rm = TRUE)
    it <- sum(pf$Fraction[pf$Component == "Intersection"], na.rm = TRUE)
    rs <- sum(pf$Fraction[pf$Component == "Residuals"],    na.rm = TRUE)
    needs <- (is.finite(bt) && is.finite(tr) && (bt >= tr)) || (bt > 0.05)
    assess_rows[["Bray-Curtis (TSS)"]] <- tibble::tibble(
      Geometry = "Bray-Curtis (TSS)",
      Treatment = tr, Batch = bt, Intersection = it, Residuals = rs,
      Score = if (is.finite(tr) && is.finite(bt)) tr / (tr + bt + 1e-12) else NA_real_,
      Needs_Correction = needs
    )
  }
  
  assess_df <- dplyr::bind_rows(assess_rows)
  print(assess_df, n = nrow(assess_df))
  readr::write_csv(assess_df, file.path(output_folder, "pRDA_raw_assessment.csv"))
  
  if (any(assess_df$Needs_Correction, na.rm = TRUE)) {
    message("pRDA baseline: Batch fraction ≥ Treatment (and/or Batch > 0.05) in at least one geometry — correction recommended.")
  } else {
    message("pRDA baseline: Batch fraction appears modest relative to Treatment — correction may not be necessary.")
  }
  
} else {
  # ---- Normal multi-method ranking ----
  unified <- dplyr::full_join(tb_ait, tb_bra, by = "Method") %>%
    rowwise() %>%
    mutate(
      Combined_Score = {
        s <- c(Score_Aitchison, Score_Bray)
        s <- s[is.finite(s)]
        if (!length(s)) NA_real_ else exp(mean(log(s)))
      }
    ) %>%
    ungroup() %>%
    arrange(desc(Combined_Score)) %>%
    mutate(Rank = row_number())
  
  print(unified, n = nrow(unified))
  readr::write_csv(unified, file.path(output_folder, "pRDA_ranking.csv"))
}