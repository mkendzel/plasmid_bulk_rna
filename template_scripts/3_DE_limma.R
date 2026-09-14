# limma-voom differential expression over the contrast registry in 0_config.R.
# Reads the checkpoints from 1_import.R and the QC verdict recorded in
# DROP_SAMPLES. Run from the repo root.

# ---- Libraries ----
# Bioconductor first so dplyr masks S4Vectors::rename/select, not the reverse
library(edgeR)
library(limma)
library(AnnotationDbi)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)

# ---- Load helper functions + project config ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("template_scripts/0_config.R")

library(ORG_DB, character.only = TRUE)

# ---- Load from 1_import.R ----
# load_checkpoint() takes the highest version; it should match the
# "Counts checkpoint" row of results/qc_report.md.
counts       <- load_checkpoint("counts_raw", dir = dir_rds)
vendor_genes <- load_checkpoint("vendor_genes", dir = dir_rds)

# ---- Sample filtering ----
counts <- counts[, setdiff(colnames(counts), DROP_SAMPLES), drop = FALSE]

# Replication recomputed from the surviving samples; scripts 4-5 read min_n to
# flag unreplicated contrasts
cond_n_kept <- sample_map |>
  dplyr::filter(sample %in% colnames(counts)) |>
  dplyr::count(condition, .drop = FALSE, name = "n") |>
  tibble::deframe()

contrast_registry$min_n <- contrast_min_n(contrast_registry, cond_n_kept)

unreplicated <- names(cond_n_kept)[cond_n_kept < 2]
if (length(DROP_SAMPLES) > 0) {
  message("Dropped ", length(DROP_SAMPLES), " sample(s): ",
          paste(DROP_SAMPLES, collapse = ", "))
}
if (length(unreplicated) > 0) {
  message("Groups below 2 replicates: ", paste(unreplicated, collapse = ", "))
}

save_checkpoint(counts, "counts_filtered", dir = dir_rds,
                notes = paste("dropped:", paste(DROP_SAMPLES, collapse = ", ")))
save_checkpoint(contrast_registry, "contrast_registry", dir = dir_rds,
                notes = "min_n recomputed after DROP_SAMPLES")

# ---- limma-voom fit, one tissue at a time ----
# Each tissue is subset and gets its own gene filter, TMM factors, voom trend and
# eBayes fit. Cell-means design (~ 0 + treatment); makeContrasts builds that
# tissue's comparisons from the registry expressions.
fit_tissue <- function(tis, counts, registry) {

  sm <- sample_map |>
    dplyr::filter(tissue == tis, sample %in% colnames(counts)) |>
    dplyr::mutate(treatment = droplevels(treatment),
                  condition = droplevels(condition))

  mat <- round(as.matrix(counts[, sm$sample, drop = FALSE]))

  design <- model.matrix(~ 0 + treatment, data = sm)
  colnames(design) <- levels(sm$treatment)

  dge <- edgeR::DGEList(counts = mat, group = sm$treatment)
  dge$samples <- cbind(
    dge$samples,
    as.data.frame(sm[, c("sample", "treatment", "tissue", "rep")])
  )

  keep <- edgeR::filterByExpr(dge, design = design)
  message(tis, ": keeping ", sum(keep), " of ", length(keep), " genes")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)

  v <- NULL
  save_fig(function() v <<- limma::voom(dge, design, plot = TRUE),
           paste0("voom_mean_variance_", tis),
           width = 7, height = 5, subdir = "qc")

  fit <- limma::lmFit(v, design)

  reg <- registry[registry$tissue == tis, ]

  # Drop contrasts whose groups lost all samples during QC
  has_groups <- vapply(reg$expr, function(e) {
    all(.groups_in_expr(e) %in% colnames(design))
  }, logical(1))

  if (any(!has_groups)) {
    message(tis, ": skipping ", sum(!has_groups),
            " contrast(s) with no remaining samples: ",
            paste(reg$name[!has_groups], collapse = ", "))
    reg <- reg[has_groups, ]
  }

  cm <- limma::makeContrasts(contrasts = reg$expr, levels = design)
  colnames(cm) <- reg$name

  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))

  tt <- stats::setNames(
    lapply(reg$name, function(cn) {
      limma::topTable(fit2, coef = cn, number = Inf, sort.by = "none")
    }),
    reg$name
  )

  list(sample_meta = sm, design = design, dge = dge,
       voom = v, fit = fit2, contrast_matrix = cm, toptables = tt)
}

fits <- stats::setNames(
  lapply(TISSUE_LEVELS, fit_tissue, counts = counts, registry = contrast_registry),
  TISSUE_LEVELS
)

lapply(fits, function(f) dim(f$design))
lapply(fits, function(f) length(f$toptables))

# Every tissue's toptables in one list, named by contrast
toptables <- do.call(c, unname(lapply(fits, function(f) f$toptables)))

# ---- Map gene IDs ----
all_ensembl <- toptables |>
  purrr::map(rownames) |>
  unlist(use.names = FALSE) |>
  strip_ens_version() |>
  unique()

# SYMBOL from the offline map in genesets/; ENTREZID from the org.db package,
# which the KEGG block in script 4 needs.
gene_map <- readRDS(GENE_MAP_PATH) |>
  tibble::as_tibble() |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(ensembl_id = ENSEMBL) |>
  dplyr::filter(ensembl_id %in% all_ensembl)

entrez_map <- suppressMessages(AnnotationDbi::select(
  get(ORG_DB),
  keys    = all_ensembl,
  keytype = "ENSEMBL",
  columns = "ENTREZID"
)) |>
  tibble::as_tibble() |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(ensembl_id = ENSEMBL, entrez_id = ENTREZID)

gene_map <- dplyr::full_join(gene_map, entrez_map, by = "ensembl_id")

# Vendor gene names fill the rows the map misses
if ("gene_name" %in% colnames(vendor_genes)) {
  vendor_sym <- vendor_genes |>
    dplyr::transmute(
      ensembl_id = strip_ens_version(gene_id),
      gene_name  = dplyr::na_if(gene_name, "")
    ) |>
    dplyr::distinct(ensembl_id, .keep_all = TRUE)

  gene_map <- gene_map |>
    dplyr::left_join(vendor_sym, by = "ensembl_id") |>
    dplyr::mutate(SYMBOL = dplyr::coalesce(SYMBOL, gene_name)) |>
    dplyr::select(-gene_name)
}

# ---- Append identifiers to every contrast ----
annotate_toptables <- function(tt_list, gene_map) {
  purrr::map(tt_list, function(df) {
    df |>
      as.data.frame() |>
      tibble::rownames_to_column("ensembl_id") |>
      dplyr::mutate(ensembl_id = strip_ens_version(ensembl_id)) |>
      dplyr::left_join(gene_map, by = "ensembl_id") |>
      tibble::as_tibble()
  })
}

tt_annotated <- annotate_toptables(toptables, gene_map)

# ---- Save checkpoints ----
# Downstream columns are limma's: logFC, AveExpr, t, P.Value, adj.P.Val, B.
# Expression values are voom log2-CPM (v$E). Per-tissue objects are lists named
# by tissue; tt_annotated is named by contrast.
save_checkpoint(lapply(fits, `[[`, "dge"),         "dge",          dir = dir_rds)
save_checkpoint(lapply(fits, `[[`, "voom"),        "voom",         dir = dir_rds)
save_checkpoint(lapply(fits, `[[`, "fit"),         "fit",          dir = dir_rds)
save_checkpoint(lapply(fits, function(f) f$voom$E), "logcpm",      dir = dir_rds)
save_checkpoint(lapply(fits, `[[`, "sample_meta"), "sample_meta",  dir = dir_rds)
save_checkpoint(tt_annotated,                      "tt_annotated", dir = dir_rds)

# ---- DE counts per contrast ----
de_counts <- purrr::imap_dfr(tt_annotated, function(df, cn) {
  d <- de_direction(df$logFC, df$adj.P.Val)
  tibble::tibble(contrast = cn, up = sum(d == "Up"), down = sum(d == "Down"))
}) |>
  dplyr::left_join(contrast_registry, by = c("contrast" = "name")) |>
  dplyr::select(tissue, type, contrast, label, min_n, up, down) |>
  dplyr::arrange(tissue, type, contrast)

de_counts

write.table(de_counts, file.path(ensure_dir(dir_res), "de_gene_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Every contrast's toptable, one sheet per contrast
wb <- openxlsx::createWorkbook()
for (cn in names(tt_annotated)) {
  sheet <- substr(gsub("[^A-Za-z0-9]", "_", cn), 1, 31)
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, tt_annotated[[cn]])
  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
}
openxlsx::saveWorkbook(wb, file.path(dir_res, "de_toptables.xlsx"), overwrite = TRUE)
