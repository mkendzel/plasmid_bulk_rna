# Sample QC on the counts checkpoint from 1_import.R. Writes figures/qc/ and
# results/qc_report.md. Nothing is filtered and no checkpoint is written; samples
# to remove go in DROP_SAMPLES in 0_config.R and are dropped by 3_DE_limma.R.
# Run from the repo root.

# ---- Libraries ----
library(DESeq2)
library(RNAseqQC)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# ---- Load helper functions + project config ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("template_scripts/0_config.R")

# ---- Load from 1_import.R ----
counts       <- load_checkpoint("counts_raw", dir = dir_rds)
vendor_genes <- load_checkpoint("vendor_genes", dir = dir_rds)

cnt <- round(as.matrix(counts[, sample_map$sample, drop = FALSE]))

bt_idx  <- match(rownames(cnt), vendor_genes$gene_id)
biotype <- vendor_genes$gene_biotype[bt_idx]
symbol  <- vendor_genes$gene_name[bt_idx]

# ---- Count-derived metrics ----
lib_size <- colSums(cnt)

# Percent of each library on the selected genes. NA, not 0, when the annotation
# carries none of the requested biotypes.
read_share <- function(keep) {
  keep <- !is.na(keep) & keep
  if (!any(keep)) return(rep(NA_real_, ncol(cnt)))
  colSums(cnt[keep, , drop = FALSE]) / lib_size * 100
}

qc_summary <- tibble::tibble(
  sample             = colnames(cnt),
  total_counts       = lib_size,
  detected_genes     = colSums(cnt > 0),
  genes_min10        = colSums(cnt >= QC_PARAMS$detection_min),
  top100_frac        = apply(cnt, 2, function(x) {
    sum(sort(x, decreasing = TRUE)[seq_len(QC_PARAMS$complexity_top_n)]) / sum(x) * 100
  }),
  rRNA_pct           = read_share(biotype %in% c("rRNA", "rRNA_pseudogene")),
  mito_pct           = read_share(grepl(MITO_PREFIX, symbol)),
  protein_coding_pct = read_share(biotype == "protein_coding")
) |>
  dplyr::left_join(
    dplyr::select(sample_map, sample, treatment, tissue, condition),
    by = "sample"
  )

# ---- Vendor alignment metrics ----
# Per-sample mapping TSVs, two columns of Metric / Value pairs; each metric becomes
# a snake_case column. Skipped when the delivery archive is absent.
vk <- vendor_key()

if (!is.null(vk)) {
  align_stats <- do.call(rbind, lapply(seq_len(nrow(vk)), function(i) {
    m <- read_delivery(zip_sample_stats(vk$idx[i], vk$code[i]),
                       check.names = FALSE, stringsAsFactors = FALSE)
    v <- stats::setNames(as.list(m$Value),
                         tolower(gsub("[^A-Za-z0-9]+", "_", m$Metric)))
    tibble::as_tibble(c(list(idx = vk$idx[i]), v))
  })) |>
    dplyr::rename(dedup_rate = dedup_mapped_reads_rate) |>
    dplyr::left_join(dplyr::select(sample_map, idx, sample), by = "idx") |>
    dplyr::select(-idx)

  qc_summary <- dplyr::left_join(qc_summary, align_stats, by = "sample")
}

# ---- vst-derived metrics ----
dds <- DESeq2::DESeqDataSetFromMatrix(
  cnt,
  colData = as.data.frame(sample_map),
  design  = ~ 1
)

SummarizedExperiment::rowData(dds)$gene_name    <- symbol
SummarizedExperiment::rowData(dds)$gene_biotype <- biotype

vsd <- DESeq2::vst(
  RNAseqQC::filter_genes(dds,
                         min_count = QC_PARAMS$filter_min_count,
                         min_rep   = QC_PARAMS$filter_min_rep),
  blind = TRUE
)

# Correlation and replicate agreement measured within condition, against a
# sample's own replicates.
vmat <- SummarizedExperiment::assay(vsd)[, qc_summary$sample, drop = FALSE]
grp  <- as.character(qc_summary$condition)

cor_mat <- stats::cor(vmat, method = "pearson")
diag(cor_mat) <- NA

peers <- lapply(seq_along(grp), function(i) setdiff(which(grp == grp[i]), i))

qc_summary$within_condition_cor <- vapply(seq_along(peers), function(i) {
  if (length(peers[[i]]) == 0) NA_real_
  else stats::median(cor_mat[i, peers[[i]]], na.rm = TRUE)
}, numeric(1))

qc_summary$median_abs_M <- vapply(seq_along(peers), function(i) {
  p <- peers[[i]]
  if (length(p) == 0) return(NA_real_)
  ref <- if (length(p) == 1) vmat[, p] else matrixStats::rowMedians(vmat[, p, drop = FALSE])
  stats::median(abs(vmat[, i] - ref), na.rm = TRUE)
}, numeric(1))

# One PCA feeds both pca_centroid_dist and the PCA figure
top_var <- order(matrixStats::rowVars(vmat), decreasing = TRUE)[seq_len(QC_PARAMS$pca_n_feats)]
pca     <- stats::prcomp(t(vmat[top_var, ]), scale. = FALSE)
pca_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc_x    <- pca$x[, 1:2]

centroid <- t(vapply(split(as.data.frame(pc_x), grp), colMeans, numeric(2)))[grp, , drop = FALSE]
qc_summary$pca_centroid_dist <- sqrt(rowSums((pc_x - centroid)^2))
qc_summary$PC1 <- pc_x[, 1]
qc_summary$PC2 <- pc_x[, 2]

# ---- Evaluate the checks ----
qc_checks <- dplyr::filter(QC_THRESHOLDS, metric %in% colnames(qc_summary))

qc_status <- do.call(rbind, lapply(seq_len(nrow(qc_checks)), function(i) {
  ck <- qc_checks[i, ]
  v  <- as.numeric(qc_summary[[ck$metric]])

  status <- if (ck$direction == "higher") {
    ifelse(v < ck$fail, "FAIL", ifelse(v < ck$warn, "WARN", "PASS"))
  } else {
    ifelse(v > ck$fail, "FAIL", ifelse(v > ck$warn, "WARN", "PASS"))
  }

  tibble::tibble(
    sample    = qc_summary$sample,
    metric    = ck$label,
    category  = ck$category,
    unit      = ck$unit,
    direction = ck$direction,
    value     = v,
    warn_at   = ck$warn,
    fail_at   = ck$fail,
    status    = ifelse(is.na(v), "NA", status)
  )
})) |>
  dplyr::mutate(
    sample = factor(sample, levels = sample_map$sample),
    metric = factor(metric, levels = qc_checks$label),
    status = factor(status, levels = c("PASS", "WARN", "FAIL", "NA"))
  )

qc_fail_samples <- sort(unique(as.character(
  qc_status$sample[qc_status$status == "FAIL"]
)))

n_status <- table(qc_status$status)
message(sprintf("QC: %d pass, %d warn, %d fail over %d samples x %d checks",
                n_status[["PASS"]], n_status[["WARN"]], n_status[["FAIL"]],
                dplyr::n_distinct(qc_status$sample), nrow(qc_checks)))
if (length(qc_fail_samples) > 0) {
  message("FAIL on at least one check: ", paste(qc_fail_samples, collapse = ", "))
}

# ---- Plots ----
fig_qc <- "qc"

# Sample-labelled x axis, shared by the bar and PCA figures
x_samples <- theme(
  axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                    colour = "black", size = 7),
  panel.grid.major.x = element_blank()
)

# Per-sample bars, one facet per check. Threshold lines drawn only where they fall
# inside the observed range.
guides_df <- qc_status |>
  dplyr::group_by(metric) |>
  dplyr::summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE),
                   warn_at = warn_at[1], fail_at = fail_at[1], .groups = "drop") |>
  tidyr::pivot_longer(c(warn_at, fail_at), names_to = "level", values_to = "y") |>
  dplyr::filter(y >= lo, y <= hi)

qc_bars <- ggplot(qc_status, aes(x = sample, y = value, fill = status)) +
  geom_col() +
  geom_hline(data = guides_df, aes(yintercept = y, linetype = level),
             colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = qc_status_cols, drop = FALSE) +
  scale_linetype_manual(values = c(warn_at = "dashed", fail_at = "dotted"),
                        labels = c(warn_at = "warn", fail_at = "fail")) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = NULL, y = NULL, fill = NULL, linetype = NULL) +
  theme_bw() +
  x_samples +
  theme(legend.position = "top")

save_fig(qc_bars, "qc_metrics_by_sample",
         width = 13, height = max(6, 2 + 2.2 * ceiling(nrow(qc_checks) / 2)),
         subdir = fig_qc)

# PCA from the same prcomp as pca_centroid_dist. Only flagged samples are labelled.
pca_df <- qc_summary |>
  dplyr::mutate(flagged = sample %in% qc_status$sample[qc_status$status != "PASS"])

qc_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = treatment, shape = tissue)) +
  geom_point(size = 3, colour = "grey20", alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(pca_df, flagged),
    aes(x = PC1, y = PC2, label = sample), size = 2.6, colour = "grey15",
    min.segment.length = 0, max.overlaps = Inf, inherit.aes = FALSE,
    box.padding = 0.4
  ) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  scale_shape_manual(values = tissue_shapes) +
  labs(x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)"),
       fill = NULL, shape = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw()

save_fig(qc_pca, "qc_pca", width = 8, height = 6, subdir = fig_qc)

# RNAseqQC distribution views, run-level and unlabelled
save_fig(RNAseqQC::plot_total_counts(dds),       "qc_total_counts",       width = 9, height = 6, subdir = fig_qc)
save_fig(RNAseqQC::plot_library_complexity(dds, show_progress = FALSE),
                                                 "qc_library_complexity", width = 9, height = 6, subdir = fig_qc)
save_fig(RNAseqQC::plot_gene_detection(dds),     "qc_gene_detection",     width = 9, height = 6, subdir = fig_qc)
save_fig(RNAseqQC::plot_biotypes(dds),           "qc_biotypes",           width = 9, height = 6, subdir = fig_qc)
save_fig(RNAseqQC::mean_sd_plot(vsd),            "qc_mean_sd",            width = 7, height = 5, subdir = fig_qc)
save_fig(RNAseqQC::plot_pca_scatters(vsd, n_PCs = QC_PARAMS$pca_n_pcs,
                                     n_feats = QC_PARAMS$pca_n_feats,
                                     color_by = "treatment", shape_by = "tissue"),
         "qc_pca_scatters", width = 11, height = 10, subdir = fig_qc)

# ComplexHeatmap, not ggplot: drawn through a function so save_fig() gets a device
save_fig(
  function() {
    ComplexHeatmap::draw(
      RNAseqQC::plot_sample_clustering(vsd,
                                       n_feats   = QC_PARAMS$pca_n_feats,
                                       anno_vars = c("treatment", "tissue"),
                                       distance  = "euclidean")
    )
  },
  "qc_sample_clustering", width = 10, height = 9, subdir = fig_qc
)

# Status grid: one tile per sample x check, value in the tile, bold when flagged.
# Metrics reversed to match QC_THRESHOLDS order top to bottom.
qc_tile <- ggplot(qc_status, aes(x = sample, y = metric, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(
    aes(label = fmt_qc(value, unit, compact = TRUE),
        fontface = ifelse(status %in% c("WARN", "FAIL"), "bold", "plain")),
    size = 2.3, colour = "grey15"
  ) +
  scale_fill_manual(values = qc_status_cols, drop = FALSE) +
  scale_y_discrete(limits = rev(levels(qc_status$metric))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
    axis.text.y     = element_text(colour = "black"),
    panel.grid      = element_blank(),
    legend.position = "top"
  )

save_fig(qc_tile, "qc_status_grid",
         width = max(8, 2 + 0.45 * nrow(sample_map)),
         height = 1 + 0.42 * nrow(qc_checks), subdir = fig_qc)

# ---- Report ----
write_qc_report(
  path              = file.path(ensure_dir(dir_res), "qc_report.md"),
  title             = paste("QC report —", basename(PROJ)),
  qc_summary        = qc_summary,
  qc_status         = qc_status,
  thresholds        = QC_THRESHOLDS,
  params            = QC_PARAMS,
  prov              = capture_provenance(),
  counts_checkpoint = basename(tail(list.files(dir_rds,
                                               pattern = "^counts_raw_v\\d+\\.rds$",
                                               full.names = TRUE), 1))
)
