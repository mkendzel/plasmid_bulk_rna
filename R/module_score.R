# Per-sample marker-module scores from logCPM.
#
# A bulk library is a blend, so a "cell type score" here is not a cell count. It
# is the average standardised expression of a hand-picked marker set, which moves
# with the fraction of that cell type but also with per-cell expression. Read it
# as composition evidence, not as a proportion.
#
# Scoring is z-per-gene-then-mean rather than mean-of-logCPM: without the z step
# a single high-abundance marker (GFAP, VIM) dominates its module and the score
# stops describing the set. z is taken across all samples of one experiment, so
# scores are only comparable within an experiment - the same rule the logFC and
# NES values follow.

# Standardise a vector, returning 0 rather than NaN for a constant gene. A marker
# with no variance carries no composition information; contributing 0 leaves it
# out of the module mean instead of poisoning it.
safe_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# Long frame of module marker expression, one row per gene-sample, with the
# z-score alongside the raw logCPM. Feeds both module_score() and the marker
# heatmap, so the two figures cannot disagree about which locus a symbol used.
#
#   modules  named list of symbol vectors
module_expr_z <- function(lc, sm, map, modules) {

  out <- purrr::imap_dfr(modules, function(genes, mod) {
    d <- gene_expr_long(lc, sm, map, genes)
    if (nrow(d) == 0L) return(tibble::tibble())
    d |>
      dplyr::group_by(SYMBOL) |>
      dplyr::mutate(z = safe_z(logcpm)) |>
      dplyr::ungroup() |>
      dplyr::mutate(module = mod)
  })

  # Every module empty gives a 0x0 tibble, which has no `module` column to
  # factor. Callers guard on nrow(), so hand back the empty frame untouched.
  if (nrow(out) == 0L) return(out)

  dplyr::mutate(out, module = factor(module, levels = names(modules)))
}

# One score per sample per module: the mean z across that module's detected
# markers. n_genes travels with the score because a module that lost half its
# markers to the detection filter is a weaker measurement, and the caller should
# be able to say so on the figure.
module_score <- function(z_long) {
  if (nrow(z_long) == 0L) return(tibble::tibble())
  z_long |>
    dplyr::group_by(sample, line, stim, condition, module) |>
    dplyr::summarise(
      n_genes = dplyr::n_distinct(SYMBOL),
      score   = mean(z, na.rm = TRUE),
      .groups = "drop"
    )
}

# Wide sample x module frame plus a single maturity axis.
#
# maturity = mean(mature modules) - mean(immature modules). One number per sample
# collapses the four-panel story into something that can go on an axis or into a
# regression as a covariate. `mature` and `immature` name entries of `modules`.
module_score_wide <- function(scores, mature, immature) {
  if (nrow(scores) == 0L) return(tibble::tibble())

  w <- scores |>
    dplyr::select(sample, line, stim, condition, module, score) |>
    tidyr::pivot_wider(names_from = module, values_from = score)

  present <- function(x) intersect(x, names(w))
  mat <- present(mature)
  imm <- present(immature)

  w |>
    dplyr::mutate(
      mature_mean   = rowMeans(dplyr::across(dplyr::all_of(mat)),  na.rm = TRUE),
      immature_mean = rowMeans(dplyr::across(dplyr::all_of(imm)), na.rm = TRUE),
      maturity      = mature_mean - immature_mean
    )
}
