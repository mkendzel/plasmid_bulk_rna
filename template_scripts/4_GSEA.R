# Pathway analysis on the limma contrasts from 3_DE_limma.R: Hallmark GSEA, an
# optional custom gene set, and KEGG over-representation. Tables and checkpoints
# only; figures are drawn in 5_Graphs.R. Run from the repo root.

# ---- Libraries ----
library(fgsea)
library(clusterProfiler)
library(AnnotationDbi)
library(dplyr)
library(tibble)
library(purrr)

# ---- Load helper functions + project config ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("template_scripts/0_config.R")

library(ORG_DB, character.only = TRUE)

# ---- Load from 3_DE_limma.R ----
tt <- load_checkpoint("tt_annotated", dir = dir_rds)

contrast_registry <- load_checkpoint("contrast_registry", dir = dir_rds)

ensure_dir(dir_res)

# Registry columns joined onto every results table
reg_cols <- c("name", "type", "treatment", "ref_treatment", "tissue",
              "label", "min_n")

# ---- Hallmark GSEA ----
# run_gsea() (R/run_gsea.R) ranks by moderated t, de-duplicating symbols on
# highest AveExpr.
hallmark_sets <- fgsea::gmtPathways(HALLMARK_GMT)

gsea_results <- purrr::imap(tt, function(df, contrast_name) {
  message("fgsea: ", contrast_name)
  run_gsea(as.data.frame(df), hallmark_sets, gene_col = "SYMBOL")
})

# leadingEdge is a list column; collapse it for writing to disk
flatten_le <- function(x) {
  if (!"leadingEdge" %in% names(x)) return(x)
  x$leadingEdge <- vapply(x$leadingEdge, paste, character(1), collapse = ";")
  x
}

gsea_tbl_all <- purrr::imap_dfr(gsea_results, function(res, contrast_name) {
  tibble::as_tibble(flatten_le(as.data.frame(res))) |>
    dplyr::mutate(contrast = contrast_name, .before = 1)
}) |>
  dplyr::left_join(dplyr::select(contrast_registry, dplyr::all_of(reg_cols)),
                   by = c("contrast" = "name"), suffix = c("", "_reg")) |>
  dplyr::arrange(tissue, type, contrast, padj)

write.table(gsea_tbl_all, file.path(dir_res, "gsea_hallmark_all.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

save_checkpoint(gsea_tbl_all, "gsea_hallmark", dir = dir_rds,
                notes = "fgsea on Hallmark, ranked by moderated t")

# ---- Custom geneset GSEA ----
# Runs only when CUSTOM_GMT names a file. Ranked on moderated t, matching the
# Hallmark ranking above.
if (!is.na(CUSTOM_GMT)) {

  gmt_sym <- clusterProfiler::read.gmt(CUSTOM_GMT)

  sym2ent <- suppressMessages(AnnotationDbi::select(
    get(ORG_DB),
    keys    = unique(gmt_sym$gene),
    keytype = "SYMBOL",
    columns = "ENTREZID"
  ))

  gmt_ent <- gmt_sym |>
    dplyr::left_join(
      tibble::as_tibble(sym2ent) |>
        dplyr::rename(gene = SYMBOL, entrez_id = ENTREZID) |>
        dplyr::filter(!is.na(entrez_id)) |>
        dplyr::distinct(gene, .keep_all = TRUE),
      by = "gene"
    ) |>
    dplyr::filter(!is.na(entrez_id)) |>
    dplyr::transmute(term = term, gene = entrez_id)

  ranked_entrez <- purrr::map(tt, function(df) {
    df2 <- df |>
      dplyr::filter(!is.na(entrez_id), !is.na(t)) |>
      dplyr::mutate(entrez_id = as.character(entrez_id)) |>
      dplyr::arrange(dplyr::desc(abs(t))) |>
      dplyr::distinct(entrez_id, .keep_all = TRUE)

    sort(stats::setNames(df2$t, df2$entrez_id), decreasing = TRUE)
  })

  # pvalueCutoff/minGSSize = 1 keeps every result for downstream filtering
  gsea_custom <- purrr::imap(ranked_entrez, function(geneList, contrast_name) {
    message("custom GSEA: ", contrast_name)
    suppressWarnings(clusterProfiler::GSEA(
      geneList     = geneList,
      TERM2GENE    = gmt_ent,
      pvalueCutoff = 1,
      minGSSize    = 1
    ))
  })

  gsea_tbl_custom <- purrr::imap_dfr(gsea_custom, function(res, contrast_name) {

    out <- tibble::as_tibble(res@result)

    if (nrow(out) == 0) {
      return(tibble::tibble(
        contrast = contrast_name, term = NA_character_,
        setSize = NA_integer_, NES = NA_real_, pvalue = NA_real_,
        p.adjust = NA_real_
      ))
    }

    out |>
      dplyr::transmute(
        contrast = contrast_name, term = ID,
        setSize, NES, pvalue, p.adjust
      )
  }) |>
    dplyr::left_join(dplyr::select(contrast_registry, dplyr::all_of(reg_cols)),
                     by = c("contrast" = "name"), suffix = c("", "_reg")) |>
    dplyr::mutate(
      NES      = round(NES, 3),
      pvalue   = signif(pvalue, 3),
      p.adjust = signif(p.adjust, 3)
    ) |>
    dplyr::arrange(tissue, p.adjust)

  write.table(gsea_tbl_custom, file.path(dir_res, "gsea_custom_all.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  save_checkpoint(gsea_tbl_custom, "gsea_custom", dir = dir_rds,
                  notes = paste("clusterProfiler GSEA on", basename(CUSTOM_GMT)))
}

# ---- KEGG over-representation ----
# Significant genes per contrast, called with the same cutoffs as everything else.
sig_entrez <- purrr::map(tt, function(df) {
  d <- df[is_de(df$logFC, df$adj.P.Val) & !is.na(df$entrez_id), ]
  unique(as.character(d$entrez_id))
})

kegg_results <- purrr::imap(sig_entrez, function(genes, contrast_name) {
  if (length(genes) == 0) return(NULL)
  message("KEGG: ", contrast_name, " (", length(genes), " genes)")
  clusterProfiler::enrichKEGG(gene = genes, organism = KEGG_ORG)
})

kegg_tbl_all <- purrr::imap_dfr(kegg_results, function(res, contrast_name) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(tibble::tibble())
  tibble::as_tibble(as.data.frame(res)) |>
    dplyr::mutate(contrast = contrast_name, .before = 1)
})

if (nrow(kegg_tbl_all) > 0) {
  kegg_tbl_all <- kegg_tbl_all |>
    dplyr::left_join(dplyr::select(contrast_registry, dplyr::all_of(reg_cols)),
                     by = c("contrast" = "name"), suffix = c("", "_reg")) |>
    dplyr::arrange(tissue, contrast, p.adjust)

  write.table(kegg_tbl_all, file.path(dir_res, "kegg_all.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  save_checkpoint(kegg_tbl_all, "kegg", dir = dir_rds,
                  notes = paste("enrichKEGG on", KEGG_ORG, "significant genes"))
}

# ---- KEGG pathway diagrams ----
# pathview writes PNGs into the working directory, so each call runs inside the
# contrast's output folder. PATHWAY_IDS = NULL draws the top pathways per
# contrast; a character vector draws exactly those.
PATHWAY_IDS <- NULL
TOP_PATHWAYS <- 5

for (cn in names(kegg_results)) {

  res <- kegg_results[[cn]]
  if (is.null(res)) next

  ids <- if (is.null(PATHWAY_IDS)) {
    utils::head(as.data.frame(res)$ID, TOP_PATHWAYS)
  } else {
    PATHWAY_IDS
  }
  if (length(ids) == 0) next

  df <- tt[[cn]] |>
    dplyr::filter(!is.na(entrez_id), !is.na(logFC)) |>
    dplyr::distinct(entrez_id, .keep_all = TRUE)

  logFC <- stats::setNames(df$logFC, as.character(df$entrez_id))

  outdir <- ensure_dir(dir_fig, "pathways", cn)

  for (pid in ids) {
    withr::with_dir(outdir, {
      pathview::pathview(
        gene.data   = logFC,
        pathway.id  = pid,
        species     = KEGG_ORG,
        limit       = list(gene = 3, cpd = 1),
        out.suffix  = cn,
        kegg.native = TRUE
      )
    })
  }
}
