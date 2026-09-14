# Bar plot of decoupleR ULM activity scores (CollecTRI TFs or PROGENy pathways).
# score < 0 is higher in group_neg, score > 0 in group_pos.
# act_df: tidy table with source, score, padj columns.
# n_show: number of sources to plot, taken by largest |score|. Inf plots all.
# padj_cut: keep only sources at or below this padj. NA skips the filter.
# keep_sources: character vector of exact source names. When given, those are
#   plotted and both n_show and padj_cut are skipped.
# Returns NULL when no source survives the filters.
make_activity_barplot <- function(act_df, plot_title, group_neg, group_pos,
                                  n_show = 20, padj_cut = 0.05,
                                  keep_sources = NULL,
                                  value_label = "TF activity (ULM score)") {

  if (is.null(keep_sources)) {
    if (!is.na(padj_cut)) act_df <- act_df[!is.na(act_df$padj) & act_df$padj <= padj_cut, ]
    act_df <- act_df[order(-abs(act_df$score)), ]
    act_df <- act_df[seq_len(min(n_show, nrow(act_df))), ]
  } else {
    act_df <- act_df[act_df$source %in% keep_sources, ]
  }

  if (nrow(act_df) == 0) {
    warning("make_activity_barplot: nothing to plot for '", plot_title,
            "', returning NULL", call. = FALSE)
    return(NULL)
  }

  act_df <- act_df[order(act_df$score), ]
  act_df$source <- factor(act_df$source, levels = act_df$source)
  act_df$direction <- factor(ifelse(act_df$score > 0, group_pos, group_neg),
                             levels = c(group_neg, group_pos))

  fill_values <- contrast_colors(group_neg, group_pos)

  x_lab <- direction_axis_label(value_label, group_neg, group_pos)

  ggplot(act_df, aes(x = score, y = source, fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_x_continuous(limits = max(abs(act_df$score)) * c(-1.05, 1.05)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = x_lab, y = NULL, fill = "Higher in", title = plot_title) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.major.y = element_blank(), legend.position = "right",
          plot.title = element_text(face = "bold", size = 17),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12))
}
