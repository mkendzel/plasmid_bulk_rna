# Ranked statistic vector for GSEA.
#
# Duplicate gene identifiers must be collapsed before ranking; `dedupe_by` picks
# which row survives.
#   "AveExpr"  keep the highest-expressed locus  - what run_gsea() does for Hallmark
#   "abs_stat" keep the most extreme statistic   - what 4_GSEA.R does for the custom set
#
# The two collections are therefore ranked by different rules for the same
# contrast. That divergence is real; unifying it would invalidate the
# checkpointed gsea_als, so it is documented rather than fixed.
#
# Figure scripts call this to rebuild the exact vector fgsea saw, so a running-ES
# curve reaches the NES reported in the checkpoint.
rank_stats <- function(tt, gene_col = "SYMBOL", stat_col = "t",
                       dedupe_by = c("AveExpr", "abs_stat"), quiet = FALSE) {

  dedupe_by <- match.arg(dedupe_by)
  tt <- as.data.frame(tt)

  if (!all(c(gene_col, stat_col) %in% names(tt))) {
    stop("rank_stats: missing column(s) ",
         paste(setdiff(c(gene_col, stat_col), names(tt)), collapse = ", "))
  }

  g    <- as.character(tt[[gene_col]])
  keep <- !is.na(g) & nzchar(g) & !is.na(tt[[stat_col]])
  tt   <- tt[keep, , drop = FALSE]
  if (nrow(tt) == 0L) return(setNames(numeric(0), character(0)))

  ord <- if (dedupe_by == "AveExpr") order(-tt$AveExpr) else order(-abs(tt[[stat_col]]))
  tt  <- tt[ord, , drop = FALSE]

  n_before <- nrow(tt)
  tt <- tt[!duplicated(as.character(tt[[gene_col]])), , drop = FALSE]
  if (!quiet) {
    message("rank_stats: removed ", n_before - nrow(tt),
            " duplicate gene names (kept by ", dedupe_by, ")")
  }

  sort(setNames(tt[[stat_col]], as.character(tt[[gene_col]])), decreasing = TRUE)
}
