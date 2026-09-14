# Long expression frame for a set of gene symbols: one row per gene-sample.
#
# Symbols with several surviving loci are collapsed to the highest-mean-logCPM
# locus. Without this the caller silently averages two loci into one bar, because
# the symbol map is deduplicated on ensembl_id, not on SYMBOL.
#
#   lc      logcpm matrix, rownames versioned Ensembl IDs
#   sm      sample_meta for the same experiment
#   map     symbol_map() output
#   symbols character vector of gene symbols to keep
gene_expr_long <- function(lc, sm, map, symbols) {

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

  if (nrow(d) == 0L) return(tibble::tibble())

  # One locus per symbol: the highest mean expression, matching the dedupe rule
  # rank_stats(dedupe_by = "AveExpr") uses.
  best <- d |>
    dplyr::group_by(SYMBOL, ensembl_id) |>
    dplyr::summarise(m = mean(logcpm, na.rm = TRUE), .groups = "drop_last") |>
    dplyr::slice_max(m, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(SYMBOL, ensembl_id)

  d |>
    dplyr::semi_join(best, by = c("SYMBOL", "ensembl_id")) |>
    dplyr::left_join(dplyr::select(sm, sample, line, stim, condition), by = "sample")
}

# Group means with SEM for a gene_expr_long() frame. n = 1 groups get NA rather
# than 0 so geom_errorbar draws nothing instead of warning.
gene_expr_summary <- function(d) {
  d |>
    dplyr::group_by(SYMBOL, condition, line, stim) |>
    dplyr::summarise(
      n         = sum(!is.na(logcpm)),
      mean_expr = mean(logcpm, na.rm = TRUE),
      sd_expr   = stats::sd(logcpm, na.rm = TRUE),
      sem_expr  = ifelse(n < 2L, NA_real_, sd_expr / sqrt(n)),
      .groups   = "drop"
    )
}
