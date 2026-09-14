# Grouped marker-expression bar plots from a count matrix.

# logCPM per marker panel: one row per panel, one column per sample.
# Counts are taken from the unfiltered matrix and summed within a panel.
# Effective library sizes are lib.size * norm.factors from the filterByExpr +
# calcNormFactors DGEList.
# A panel counts as detected when at least one treatment group has a median
# count of min_count or more. Undetected panels are NA across every sample.
#
#   counts       full count matrix, versioned Ensembl rownames
#   info         sample metadata for one tissue, with `sample` and `treatment`
#   map          symbol_map() output
#   panels       named list, panel label -> gene symbols summed into that panel
#   min_count    detection threshold on the per-group median count
#   prior_count  offset added before the log, as in voom
marker_panel_logcpm <- function(counts, info, map, panels,
                                min_count = 5, prior_count = 0.5) {

  cnt <- round(as.matrix(counts[, info$sample, drop = FALSE]))

  dge  <- edgeR::DGEList(counts = cnt, group = info$treatment)
  keep <- edgeR::filterByExpr(dge)
  dge  <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  eff_lib <- dge$samples$lib.size * dge$samples$norm.factors

  rn <- strip_ens_version(rownames(cnt))
  grp <- droplevels(factor(info$treatment))

  summed <- t(vapply(panels, function(syms) {
    rows <- which(rn %in% map$ensembl_id[map$SYMBOL %in% syms])
    if (length(rows) == 0L) return(rep(NA_real_, ncol(cnt)))
    tot <- colSums(cnt[rows, , drop = FALSE])
    detected <- any(tapply(tot, grp, stats::median) >= min_count)
    if (detected) tot else rep(NA_real_, ncol(cnt))
  }, numeric(ncol(cnt))))

  dimnames(summed) <- list(names(panels), colnames(cnt))

  log2((summed + prior_count) / (eff_lib + 1) * 1e6)
}

# Long frame for marker_panel_logcpm() output: one row per panel-sample.
# NA panels are dropped; facet_wrap(drop = FALSE) leaves them blank.
marker_panel_long <- function(m, info, group_col = "treatment") {
  tibble::tibble(
    SYMBOL = rep(rownames(m), times = ncol(m)),
    sample = rep(colnames(m), each = nrow(m)),
    logcpm = as.vector(m)
  ) |>
    dplyr::filter(!is.na(logcpm)) |>
    dplyr::left_join(
      dplyr::select(info, sample, group = dplyr::all_of(group_col)),
      by = "sample"
    ) |>
    dplyr::mutate(SYMBOL = factor(SYMBOL, levels = rownames(m))) |>
    dplyr::arrange(SYMBOL, group)
}

# Panels with no data.
missing_markers <- function(d, symbols) setdiff(symbols, as.character(unique(d$SYMBOL)))

# Group means with SEM. n = 1 groups get NA.
marker_group_summary <- function(d) {
  d |>
    dplyr::group_by(SYMBOL, group, .drop = FALSE) |>
    dplyr::summarise(
      n         = sum(!is.na(logcpm)),
      mean_expr = mean(logcpm, na.rm = TRUE),
      sem_expr  = ifelse(n < 2L, NA_real_, stats::sd(logcpm, na.rm = TRUE) / sqrt(n)),
      .groups   = "drop"
    ) |>
    dplyr::filter(n > 0L)
}

# One facet_wrap block: bar at the group mean, whiskers at +/- SEM, one point
# per sample. drop = FALSE keeps a labelled but empty panel for any symbol with
# no rows in `d`.
marker_block <- function(d, summ, symbols, block_label, ncol,
                         group_labels, group_colors, value_label,
                         y_expand, point_size, show_x) {

  d$SYMBOL    <- factor(as.character(d$SYMBOL), levels = symbols)
  summ$SYMBOL <- factor(as.character(summ$SYMBOL), levels = symbols)

  ggplot(summ, aes(x = group, y = mean_expr, fill = group)) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_expr - sem_expr, ymax = mean_expr + sem_expr),
                  width = 0.2, linewidth = 0.4, color = "grey30", na.rm = TRUE) +
    geom_jitter(data = d, aes(x = group, y = logcpm), inherit.aes = FALSE,
                width = 0.15, height = 0, size = point_size, shape = 21,
                fill = "white", color = "grey20", stroke = 0.4) +
    facet_wrap(~ SYMBOL, scales = "free_y", ncol = ncol, drop = FALSE) +
    scale_fill_manual(values = group_colors, labels = group_labels, drop = FALSE) +
    scale_x_discrete(labels = group_labels, drop = FALSE) +
    scale_y_continuous(expand = y_expand) +
    labs(x = NULL, y = value_label, fill = NULL, subtitle = block_label) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.subtitle = element_text(face = "bold", size = 13, color = "black"),
          strip.text = element_text(face = "bold", size = 13, color = "black"),
          axis.title.y = element_text(size = 13, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.text.x = if (show_x) element_text(size = 11, color = "black")
                        else element_blank())
}

# Marker panels stacked one block per gene set, blocks separated by a rule.
#
#   d             marker_expr_long() output
#   symbols       every panel to draw, in order; undetected ones stay blank
#   sets          named list of symbol vectors, one block per element, in
#                 display order. NULL draws a single unlabelled block.
#   group_levels  x-axis order; defaults to the levels present on d$group
#   group_labels  display labels, named by group level
#   group_colors  bar fills, named by group level
#   ncol          panels per row; defaults to the largest block
make_marker_barplot <- function(d, plot_title, symbols, sets = NULL,
                                group_levels = NULL, group_labels = NULL,
                                group_colors = NULL,
                                value_label = "Expression (log2 CPM)",
                                ncol = NULL, point_size = 1.8,
                                divider_color = "grey30") {

  if (is.null(sets)) sets <- stats::setNames(list(symbols), "")
  if (is.null(ncol)) ncol <- max(lengths(sets))

  if (is.null(group_levels)) group_levels <- levels(factor(d$group))
  if (is.null(group_labels)) group_labels <- stats::setNames(group_levels, group_levels)
  if (is.null(group_colors)) {
    group_colors <- stats::setNames(grDevices::hcl.colors(length(group_levels), "Set 2"),
                                    group_levels)
  }

  if (nrow(d) == 0L) {
    return(empty_plot(plot_title, note = "No marker genes detected"))
  }

  d$group  <- factor(as.character(d$group), levels = group_levels)
  d$SYMBOL <- factor(as.character(d$SYMBOL), levels = symbols)
  d <- d[!is.na(d$group), ]

  summ <- marker_group_summary(d)

  # Bottom expansion is added only when a panel reaches below zero.
  goes_negative <- any(d$logcpm < 0, na.rm = TRUE) ||
    any(summ$mean_expr - summ$sem_expr < 0, na.rm = TRUE)
  y_expand <- expansion(mult = c(if (goes_negative) 0.08 else 0, 0.08))

  blocks <- lapply(seq_along(sets), function(i) {
    block_syms <- sets[[i]]
    marker_block(
      d[as.character(d$SYMBOL) %in% block_syms, ],
      summ[as.character(summ$SYMBOL) %in% block_syms, ],
      block_syms, names(sets)[i], ncol,
      group_labels, group_colors, value_label,
      y_expand, point_size, show_x = TRUE
    )
  })

  block_rows <- ceiling(lengths(sets) / ncol)

  divider <- ggplot() +
    annotate("segment", x = 0, xend = 1, y = 0.5, yend = 0.5,
             color = divider_color, linewidth = 0.6) +
    scale_x_continuous(expand = expansion(0)) +
    theme_void()

  parts   <- vector("list", 2L * length(blocks) - 1L)
  heights <- numeric(length(parts))
  for (i in seq_along(blocks)) {
    parts[[2L * i - 1L]]   <- blocks[[i]]
    heights[2L * i - 1L]   <- block_rows[i]
    if (i < length(blocks)) {
      parts[[2L * i]] <- divider
      heights[2L * i] <- 0.10
    }
  }

  patchwork::wrap_plots(parts, ncol = 1) +
    patchwork::plot_layout(heights = heights) +
    patchwork::plot_annotation(
      title = plot_title,
      theme = theme(plot.title = element_text(face = "bold", size = 16, color = "black"))
    )
}
