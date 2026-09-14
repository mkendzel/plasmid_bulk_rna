# decoupleR ULM activity from a limma-voom topTable.
# Ranks on the moderated t statistic; positive score = higher in the numerator
# group of the contrast.
# net: decoupleR network, e.g. get_collectri() or get_progeny().
# mor_col: signed weight column, "mor" for CollecTRI, "weight" for PROGENy.
# minsize: minimum number of detected targets per source.
run_activity <- function(tt, net, contrast_name = "contrast",
                         gene_col = "SYMBOL", mor_col = "mor", minsize = 5) {
  tt <- tt[!is.na(tt[[gene_col]]) & tt[[gene_col]] != "", ]

  n_before <- nrow(tt)
  tt <- tt[order(-tt$AveExpr), ]
  tt <- tt[!duplicated(tt[[gene_col]]), ]
  message("Removed ", n_before - nrow(tt), " duplicate gene names (kept highest AveExpr)")

  # Guard against a symbol nomenclature mismatch between the data and the network
  shared <- length(intersect(unique(net$target), tt[[gene_col]]))
  message("Network targets found in data: ", shared, " of ", length(unique(net$target)))
  if (shared < 100) {
    stop("run_activity: only ", shared, " network targets matched the data; ",
         "check gene symbol nomenclature")
  }

  # Genes as rows, contrast as the single column
  mat <- matrix(tt$t, ncol = 1, dimnames = list(tt[[gene_col]], contrast_name))

  res <- decoupleR::run_ulm(
    mat     = mat,
    network = net,
    .source = "source",
    .target = "target",
    .mor    = mor_col,
    minsize = minsize
  )

  res <- as.data.frame(res)
  res$padj <- p.adjust(res$p_value, method = "BH")
  res <- res[order(res$p_value), ]
  return(res)
}
