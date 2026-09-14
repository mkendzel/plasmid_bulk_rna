# Overview, DE, GSEA and per-gene figures from the checkpoints written by
# 3_DE_limma.R and 4_GSEA.R. Every loop is driven by the contrast registry.
# contrast added in 0_config.R appears here without further edits.
# Run from the repo root.

# ---- Libraries ----
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)
library(pheatmap)

# ---- Load helper functions + project config ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("template_scripts/0_config.R")

# ---- Data load ----
# tt is named by contrast; logcpm (voom's log2 CPM) and sample_meta by tissue
tt          <- load_checkpoint("tt_annotated", dir = dir_rds)
logcpm      <- load_checkpoint("logcpm", dir = dir_rds)
sample_meta <- load_checkpoint("sample_meta", dir = dir_rds)

contrast_registry <- load_checkpoint("contrast_registry", dir = dir_rds)

# ---- Coverage check ----
# 3_DE_limma.R drops contrasts whose groups lost all samples but checkpoints the
# full registry.
coverage <- purrr::map_dfr(TISSUE_LEVELS, function(tis) {
  expect <- contrasts_for(tis)
  tibble::tibble(
    tissue     = tis,
    n_expected = length(expect),
    n_toptable = sum(expect %in% names(tt)),
    missing    = paste(setdiff(expect, names(tt)), collapse = ", ")
  )
})
print(coverage)

if (any(coverage$n_toptable != coverage$n_expected)) {
  warning("Contrasts missing from tt - re-run 3_DE_limma.R:\n",
          paste(coverage$missing[nzchar(coverage$missing)], collapse = "\n"),
          call. = FALSE)
}

# Contrasts present in tt for one tissue, in registry order
contrasts_present <- function(tis) intersect(contrasts_for(tis), names(tt))

# All toptables of one tissue stacked, with registry metadata and the Up/Down/NS
# call applied. Feeds the volcano grid and p-value histograms.
long_tt <- function(tis) {
  purrr::map_dfr(contrasts_present(tis), function(cn) {
    dplyr::mutate(tt[[cn]], contrast = cn)
  }) |>
    dplyr::left_join(
      dplyr::select(contrast_registry, name, label, type, min_n),
      by = c("contrast" = "name")
    ) |>
    dplyr::mutate(
      label     = factor(label, levels = contrast_labels(tis)),
      direction = de_direction(logFC, adj.P.Val)
    )
}

# Facet grid geometry scaled to the contrast count
facet_dims <- function(n, ncol = 4, w_per = 4, h_per = 3.2) {
  ncol <- min(ncol, n)
  list(ncol = ncol, width = max(7, ncol * w_per),
       height = max(5, ceiling(n / ncol) * h_per))
}

# Annotation frame and colours for the sample heatmaps
sample_annotation <- function(sm) {
  ann <- data.frame(treatment = sm$treatment, row.names = sm$sample)
  cols <- list(treatment = TREATMENT_COLORS[levels(droplevels(sm$treatment))])
  list(ann = ann, colors = cols)
}

# =============================================================================
# Overview
# =============================================================================
dir_overview <- "overview"

# ---- PCA ----
# prcomp once per tissue, shared by the PC1/2, PC3/4 and scree plots. Non-finite
# and zero-variance genes are dropped before ranking.
pca_fit <- function(lc, n_top = QC_PARAMS$pca_n_feats) {
  v <- matrixStats::rowVars(lc)
  names(v) <- rownames(lc)
  v <- v[is.finite(v) & v > 0]
  keep <- utils::head(names(sort(v, decreasing = TRUE)), min(n_top, length(v)))
  p <- stats::prcomp(t(lc[keep, , drop = FALSE]), scale. = FALSE)
  list(x = p$x, pct = round(100 * p$sdev^2 / sum(p$sdev^2), 1), n_used = length(keep))
}

pca_plot <- function(fit, sm, pcs = c(1, 2), title = NULL) {

  if (ncol(fit$x) < max(pcs)) {
    return(empty_plot(title, note = paste0("Fewer than PC", max(pcs), " available")))
  }

  d <- tibble::tibble(
    sample = rownames(fit$x),
    xx     = fit$x[, pcs[1]],
    yy     = fit$x[, pcs[2]]
  ) |>
    dplyr::left_join(dplyr::select(sm, sample, treatment, condition), by = "sample")

  centroids <- d |>
    dplyr::group_by(condition, treatment) |>
    dplyr::summarise(cx = mean(xx), cy = mean(yy), .groups = "drop")

  d <- dplyr::left_join(d, centroids, by = c("condition", "treatment"))

  ggplot(d, aes(xx, yy)) +
    # Legs from each sample to its group centroid
    geom_segment(aes(xend = cx, yend = cy, colour = treatment),
                 alpha = 0.5, linewidth = 0.4, show.legend = FALSE) +
    geom_point(aes(fill = treatment), shape = 21,
               size = 3.2, colour = "black", stroke = 0.5) +
    geom_point(data = centroids, aes(cx, cy, fill = treatment),
               shape = 21, size = 5, colour = "black", stroke = 1.1,
               inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = TREATMENT_COLORS, drop = TRUE) +
    scale_colour_manual(values = TREATMENT_COLORS, drop = TRUE) +
    labs(x = paste0("PC", pcs[1], ": ", fit$pct[pcs[1]], "% variance"),
         y = paste0("PC", pcs[2], ": ", fit$pct[pcs[2]], "% variance"),
         title = title) +
    theme_classic()
}

scree_plot <- function(fit, n_pc = 10, title = NULL) {
  n <- min(n_pc, length(fit$pct))
  d <- tibble::tibble(
    pc  = factor(paste0("PC", seq_len(n)), levels = paste0("PC", seq_len(n))),
    pct = fit$pct[seq_len(n)],
    cum = cumsum(fit$pct)[seq_len(n)]
  )

  ggplot(d, aes(pc, pct)) +
    geom_col(fill = "grey70", colour = "black", width = 0.7) +
    geom_line(aes(y = cum, group = 1), colour = "#D62728", linewidth = 0.7) +
    geom_point(aes(y = cum), colour = "#D62728", size = 2) +
    geom_text(aes(label = sprintf("%.0f", pct)), vjust = -0.4, size = 3) +
    labs(x = NULL, y = "% variance (bars)   |   cumulative (line)", title = title) +
    theme_classic()
}

for (tis in TISSUE_LEVELS) {

  fit <- pca_fit(logcpm[[tis]])

  save_fig(pca_plot(fit, sample_meta[[tis]], c(1, 2),
                    sprintf("PCA - %s (top %d variable genes)", tis, fit$n_used)),
           paste0("pca_", tis), width = 8, height = 6, subdir = dir_overview)

  save_fig(pca_plot(fit, sample_meta[[tis]], c(3, 4),
                    sprintf("PCA PC3/PC4 - %s", tis)),
           paste0("pca_pc34_", tis), width = 8, height = 6, subdir = dir_overview)

  save_fig(scree_plot(fit, 10, paste0("Scree - ", tis)),
           paste0("pca_scree_", tis), width = 7, height = 5, subdir = dir_overview)
}

# ---- Sample-sample correlation ----
for (tis in TISSUE_LEVELS) {

  lc <- logcpm[[tis]]
  sm <- sample_meta[[tis]]
  lc <- lc[matrixStats::rowSds(lc) > 0, , drop = FALSE]

  cm <- stats::cor(lc, method = "pearson")

  if (!all(is.finite(cm))) {
    message("Non-finite correlations for ", tis, " - skipping sample correlation heatmap")
    next
  }

  a <- sample_annotation(sm)

  ph <- pheatmap::pheatmap(
    cm, silent = TRUE,
    clustering_distance_rows = stats::as.dist(1 - cm),
    clustering_distance_cols = stats::as.dist(1 - cm),
    annotation_col    = a$ann,
    annotation_colors = a$colors,
    main = sprintf("Sample correlation (Pearson, %d genes) - %s", nrow(lc), tis)
  )

  save_fig(ph, paste0("sample_cor_", tis), width = 9, height = 8, subdir = dir_overview)
}

# ---- Heatmap: most variable genes ----
n_top_variable <- 50

for (tis in TISSUE_LEVELS) {

  lc  <- logcpm[[tis]]
  sm  <- sample_meta[[tis]]
  map <- symbol_map(tt[contrasts_present(tis)])

  v <- matrixStats::rowVars(lc)
  names(v) <- rownames(lc)
  v <- v[is.finite(v) & v > 0]
  top_ens <- utils::head(names(sort(v, decreasing = TRUE)), min(n_top_variable, length(v)))

  mat <- lc[top_ens, , drop = FALSE]
  rownames(mat) <- label_genes(top_ens, map)

  mat_scaled <- zscore_rows(mat)

  if (!heatmap_ready(mat_scaled)) {
    message("Too few variable genes for ", tis, " - skipping top-variable heatmap")
    next
  }

  desired_order <- sm |> dplyr::arrange(treatment, rep) |> dplyr::pull(sample)
  mat_scaled    <- mat_scaled[, desired_order, drop = FALSE]

  a <- sample_annotation(sm)
  a$ann <- a$ann[desired_order, , drop = FALSE]

  ph <- pheatmap::pheatmap(
    mat_scaled, silent = TRUE,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    annotation_col    = a$ann,
    annotation_colors = a$colors,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    main = paste0("Top ", nrow(mat_scaled), " variable genes - ", tis)
  )

  save_fig(ph, paste0("heatmap_top", n_top_variable, "_variable_", tis),
           width = 10, height = 9, subdir = dir_overview)
}

# ---- log2FC heatmap across contrasts ----
top_n_genes <- 250

for (tis in TISSUE_LEVELS) {

  cns     <- contrasts_present(tis)
  res_use <- tt[cns]

  sig_tbl <- purrr::map_dfr(res_use, function(df) {
    df[is_de(df$logFC, df$adj.P.Val), c("ensembl_id", "logFC")]
  })

  if (nrow(sig_tbl) == 0) {
    save_fig(empty_plot(paste0("Top DE genes, log2 fold change - ", tis),
                        note = sprintf("No gene passes adj.P.Val < %s and |logFC| > %s",
                                       P_CUT, LFC_CUT)),
             paste0("heatmap_log2FC_", tis), width = 8, height = 5, subdir = dir_overview)
    next
  }

  top_genes <- sig_tbl |>
    dplyr::group_by(ensembl_id) |>
    dplyr::summarise(max_abs_lfc = max(abs(logFC)), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(max_abs_lfc)) |>
    dplyr::slice_head(n = top_n_genes) |>
    dplyr::pull(ensembl_id)

  # Non-significant cells become NA and are drawn black
  lfc_mat <- vapply(res_use, function(df) {
    idx <- match(top_genes, df$ensembl_id)
    lfc <- df$logFC[idx]
    lfc[!is_de(lfc, df$adj.P.Val[idx])] <- NA_real_
    lfc
  }, numeric(length(top_genes)))

  reg_t <- contrast_registry |>
    dplyr::filter(name %in% cns) |>
    dplyr::arrange(type, treatment, ref_treatment)

  lfc_mat <- lfc_mat[, reg_t$name, drop = FALSE]
  colnames(lfc_mat) <- reg_t$label
  rownames(lfc_mat) <- label_genes(top_genes, symbol_map(res_use))

  ann <- data.frame(type = reg_t$type, row.names = reg_t$label)

  # NA -> 0 for the clustering distance only
  lfc_for_cluster <- lfc_mat
  lfc_for_cluster[is.na(lfc_for_cluster)] <- 0

  ph <- pheatmap::pheatmap(
    lfc_mat, silent = TRUE,
    cluster_rows      = if (heatmap_ready(lfc_for_cluster))
                          stats::hclust(stats::dist(lfc_for_cluster)) else FALSE,
    cluster_cols      = FALSE,
    annotation_col    = ann,
    annotation_colors = list(type = type_cols[levels(droplevels(reg_t$type))]),
    na_col            = "black",
    show_rownames     = FALSE,
    angle_col         = 45,
    main = paste0("Top ", nrow(lfc_mat), " DE genes, log2 fold change - ", tis,
                  "  (black = not significant)")
  )

  save_fig(ph, paste0("heatmap_log2FC_", tis),
           width = max(8, 3 + nrow(reg_t) * 0.42), height = 10, subdir = dir_overview)
}

# ---- DE gene counts across contrasts ----
for (tis in TISSUE_LEVELS) {

  # Built over the full registry so a contrast missing from tt still gets a
  # labelled slot.
  counts_t <- purrr::map_dfr(contrasts_for(tis), function(cn) {
    if (!cn %in% names(tt)) {
      return(tibble::tibble(contrast = cn, up = NA_integer_, down = NA_integer_))
    }
    d <- de_direction(tt[[cn]]$logFC, tt[[cn]]$adj.P.Val)
    tibble::tibble(contrast = cn, up = sum(d == "Up"), down = sum(d == "Down"))
  }) |>
    dplyr::left_join(dplyr::select(contrast_registry, name, label, type, min_n),
                     by = c("contrast" = "name")) |>
    dplyr::mutate(label = factor(label, levels = contrast_labels(tis)))

  plot_df <- counts_t |>
    tidyr::pivot_longer(c(up, down), names_to = "direction", values_to = "n") |>
    dplyr::mutate(
      signed    = ifelse(direction == "up", n, -n),
      direction = factor(ifelse(direction == "up", "Up", "Down"), levels = DE_LEVELS)
    )

  p <- ggplot(plot_df, aes(label, signed, fill = direction)) +
    geom_col(colour = "black", width = 0.72, linewidth = 0.25) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(n), "", format(n, big.mark = ",")),
                  vjust = ifelse(signed >= 0, -0.4, 1.3)),
              size = 2.5, na.rm = TRUE) +
    # Unreplicated contrasts flagged once, above the panel
    geom_text(data = dplyr::distinct(counts_t, label, type, min_n),
              aes(x = label, y = Inf, label = ifelse(min_n == 1, "n = 1", "")),
              inherit.aes = FALSE, vjust = 1.4, size = 2.6, colour = "#B22222") +
    facet_grid(~ type, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = DE_COLS, drop = FALSE, name = NULL) +
    scale_y_continuous(labels = function(z) format(abs(z), big.mark = ",")) +
    labs(x = NULL, y = "DE genes (up above, down below)",
         title = paste0("Differentially expressed genes - ", tis),
         subtitle = sprintf("adj.P.Val < %s and |log2FC| > %s", P_CUT, LFC_CUT)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          strip.background = element_rect(fill = "grey92", colour = NA))

  save_fig(p, paste0("de_counts_", tis),
           width = max(7, 2.2 + nrow(counts_t) * 0.62), height = 6,
           subdir = dir_overview)
}

# =============================================================================
# Per-contrast DE figures
# =============================================================================

# ---- Volcano ----
for (tis in TISSUE_LEVELS) {

  sub_v <- file.path("de", "volcano", tis)

  for (cn in contrasts_present(tis)) {
    r <- reg_of(cn)
    save_fig(
      volcano_plot(tt[[cn]], title = r$label, subtitle = cn,
                   caption = if (isTRUE(r$min_n == 1)) "Unreplicated (n = 1) - exploratory" else NULL,
                   label_n = 15),
      paste0("volcano_", cn), width = 7, height = 6, subdir = sub_v)
  }

  d  <- long_tt(tis)
  fd <- facet_dims(nlevels(droplevels(d$label)))

  p <- ggplot(d, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(colour = direction), alpha = 0.5, size = 0.5) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = -log10(P_CUT), linetype = "dashed", colour = "grey55") +
    facet_wrap(~ label, ncol = fd$ncol, drop = FALSE) +
    scale_colour_manual(values = DE_COLS, drop = FALSE, name = NULL) +
    labs(title = paste0("Volcano - ", tis), x = "log2 fold change",
         y = "-log10 adjusted p-value") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 7.5))

  save_fig(p, paste0("volcano_grid_", tis), width = fd$width, height = fd$height,
           subdir = file.path("de", "volcano"))
}

# ---- MA ----
ma_plot <- function(df, title = NULL, subtitle = NULL, caption = NULL) {
  d <- dplyr::mutate(df, direction = de_direction(logFC, adj.P.Val))

  ggplot(d, aes(AveExpr, logFC)) +
    geom_point(aes(colour = direction), alpha = 0.6, size = 1) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = DE_COLS, drop = FALSE, name = NULL) +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = "Average expression (log2 CPM)", y = "log2 fold change") +
    theme_classic() +
    theme(plot.caption = element_text(colour = "#B22222", hjust = 0))
}

for (tis in TISSUE_LEVELS) {

  sub_m <- file.path("de", "ma", tis)

  for (cn in contrasts_present(tis)) {
    r <- reg_of(cn)
    save_fig(
      ma_plot(tt[[cn]], title = r$label, subtitle = cn,
              caption = if (isTRUE(r$min_n == 1)) "Unreplicated (n = 1) - exploratory" else NULL),
      paste0("ma_", cn), width = 7, height = 6, subdir = sub_m)
  }
}

# ---- p-value histograms ----
for (tis in TISSUE_LEVELS) {

  d  <- long_tt(tis)
  fd <- facet_dims(nlevels(droplevels(d$label)))

  p <- ggplot(d, aes(P.Value)) +
    geom_histogram(binwidth = 0.02, fill = "grey70", colour = "black",
                   linewidth = 0.2) +
    facet_wrap(~ label, ncol = fd$ncol, scales = "free_y", drop = FALSE) +
    labs(title = paste0("Raw p-value distribution - ", tis),
         x = "P.Value", y = "genes") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 7.5))

  save_fig(p, paste0("pvalue_hist_", tis), width = fd$width, height = fd$height,
           subdir = file.path("de"))
}

# =============================================================================
# GSEA
# =============================================================================
# ---- Hallmark bubble plot ----
gsea_hallmark <- load_checkpoint("gsea_hallmark", dir = dir_rds)

for (tis in TISSUE_LEVELS) {

  d <- gsea_hallmark |>
    dplyr::filter(tissue == tis, !is.na(padj), padj <= P_CUT, !is.na(NES)) |>
    dplyr::mutate(
      label   = factor(label, levels = contrast_labels(tis)),
      pathway = pathway_label(pathway)
    )

  if (nrow(d) == 0) {
    save_fig(empty_plot(paste0("Hallmark GSEA - ", tis),
                        note = paste0("No pathway at padj <= ", P_CUT)),
             paste0("gsea_bubble_", tis), width = 8, height = 5, subdir = "gsea")
    next
  }

  # Pathways ordered by mean NES across the contrasts they are significant in
  d$pathway <- stats::reorder(d$pathway, d$NES, mean)

  p <- ggplot(d, aes(x = label, y = pathway, size = -log10(padj), fill = NES)) +
    geom_point(shape = 21, colour = "black") +
    scale_size(name = "-log10(padj)", range = c(2, 8)) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1, 1) * safe_absmax(d$NES)) +
    scale_x_discrete(position = "top", drop = FALSE) +
    labs(x = NULL, y = "Hallmark gene set",
         title = paste0("Hallmark GSEA - ", tis),
         subtitle = paste0("padj <= ", P_CUT)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black", angle = 45, hjust = 0))

  save_fig(p, paste0("gsea_bubble_", tis),
           width = max(8, 3 + dplyr::n_distinct(d$label) * 0.9),
           height = max(5, 1 + dplyr::n_distinct(d$pathway) * 0.22),
           subdir = "gsea")
}

# =============================================================================
# Per-gene expression
# =============================================================================
# Genes plotted one panel each. Empty runs the top DE genes of every contrast.
GENES_OF_INTEREST <- character(0)
n_top_genes <- 10

# Long expression frame for a set of symbols, one row per gene-sample. Symbols
# with several surviving loci collapse to the highest-mean-logCPM locus.
expr_long <- function(lc, sm, map, symbols) {

  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) == 0L) return(tibble::tibble())

  ens  <- map$ensembl_id[map$SYMBOL %in% symbols]
  rows <- which(strip_ens_version(rownames(lc)) %in% ens)
  if (length(rows) == 0L) return(tibble::tibble())
  lc <- lc[rows, , drop = FALSE]

  d <- tibble::tibble(
    ensembl_id = rep(strip_ens_version(rownames(lc)), times = ncol(lc)),
    sample     = rep(colnames(lc), each = nrow(lc)),
    logcpm     = as.vector(lc)
  ) |>
    dplyr::left_join(map, by = "ensembl_id") |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL %in% symbols)

  best <- d |>
    dplyr::group_by(SYMBOL, ensembl_id) |>
    dplyr::summarise(m = mean(logcpm, na.rm = TRUE), .groups = "drop_last") |>
    dplyr::slice_max(m, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(SYMBOL, ensembl_id)

  d |>
    dplyr::semi_join(best, by = c("SYMBOL", "ensembl_id")) |>
    dplyr::left_join(dplyr::select(sm, sample, treatment, tissue, condition),
                     by = "sample")
}

# Group means with SEM. n = 1 groups get NA so geom_errorbar draws nothing.
expr_summary <- function(d) {
  d |>
    dplyr::group_by(SYMBOL, condition, treatment, tissue) |>
    dplyr::summarise(
      n         = sum(!is.na(logcpm)),
      mean_expr = mean(logcpm, na.rm = TRUE),
      sd_expr   = stats::sd(logcpm, na.rm = TRUE),
      sem_expr  = ifelse(n < 2L, NA_real_, sd_expr / sqrt(n)),
      .groups   = "drop"
    )
}

gene_barplot <- function(pts, summ, gene) {
  d_pts  <- dplyr::filter(pts,  SYMBOL == gene)
  d_summ <- dplyr::filter(summ, SYMBOL == gene)

  ggplot(d_summ, aes(treatment, mean_expr, fill = treatment)) +
    geom_col(colour = "black", width = 0.72, linewidth = 0.25) +
    geom_errorbar(aes(ymin = mean_expr - sem_expr, ymax = mean_expr + sem_expr),
                  width = 0.25, na.rm = TRUE) +
    geom_point(data = d_pts, aes(treatment, logcpm, fill = treatment),
               position = position_jitter(width = 0.15, height = 0),
               shape = 21, size = 2.2, colour = "black", stroke = 0.4) +
    scale_fill_manual(values = TREATMENT_COLORS, drop = TRUE) +
    labs(x = NULL, y = "Expression (log2 CPM)", title = gene) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          plot.title  = element_text(face = "bold.italic"))
}

for (tis in TISSUE_LEVELS) {

  res_t <- tt[contrasts_present(tis)]
  map   <- symbol_map(res_t)

  symbols <- if (length(GENES_OF_INTEREST) > 0) {
    GENES_OF_INTEREST
  } else {
    unique(unlist(lapply(res_t, function(df) {
      df |>
        dplyr::filter(is_de(logFC, adj.P.Val), !is.na(SYMBOL)) |>
        dplyr::slice_min(adj.P.Val, n = n_top_genes, with_ties = FALSE) |>
        dplyr::pull(SYMBOL)
    })))
  }

  pts <- expr_long(logcpm[[tis]], sample_meta[[tis]], map, symbols)
  if (nrow(pts) == 0) {
    message("No expression rows for the requested genes in ", tis)
    next
  }

  summ <- expr_summary(pts)

  for (g in sort(unique(pts$SYMBOL))) {
    save_fig(gene_barplot(pts, summ, g), paste0("expr_", g),
             width = 6, height = 4, subdir = file.path("genes", tis))
  }
}
