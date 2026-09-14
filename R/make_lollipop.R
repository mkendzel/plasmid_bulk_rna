# Lollipop graphs of the fgsea pathways that clear nes_cut and p_cut
# ordered by NES.
# Point size is -log10(padj). NES < 0 is higher in group_neg, NES > 0 in group_pos.
# nes_cut / p_cut: selection cutoffs, |NES| > nes_cut and p < p_cut.
# p_col: which p-value column the cutoff applies to, "pval" (raw) or "padj".
#   Bubble size is always -log10(padj) regardless of this choice.
# keep_pathways: character vector of exact pathway names. When given, those
#    pathways are plotted regardless of significance and both cutoffs are skipped.
# size_range: min/max bubble diameter passed to scale_size_continuous().
make_lollipop <- function(gsea_df, plot_title, group_neg, group_pos,
                          keep_pathways = NULL, size_range = c(2.5, 8),
                          nes_cut = 1.5, p_cut = 0.05,
                          p_col = c("pval", "padj")) {
  p_col <- match.arg(p_col)

  sig_paths <- gsea_df %>%
    {
      if (is.null(keep_pathways)) filter(., abs(NES) > nes_cut, .data[[p_col]] < p_cut)
      else filter(., pathway %in% keep_pathways)
    } %>%
    mutate(
      pathway_clean = str_remove(pathway, "^HALLMARK_") %>%
        str_replace_all("_", " ") %>%
        str_to_title(),
      direction = factor(ifelse(NES > 0, group_pos, group_neg),
                         levels = c(group_neg, group_pos)),
      neg_log10_padj = -log10(padj)
    ) %>%
    arrange(NES)

  sig_paths$pathway_clean <- factor(sig_paths$pathway_clean, levels = sig_paths$pathway_clean)

  fill_values <- contrast_colors(group_neg, group_pos)

  x_lab <- direction_axis_label("Normalized Enrichment Score (NES)", group_neg, group_pos)

  ggplot(sig_paths, aes(x = NES, y = pathway_clean)) +
    geom_segment(aes(x = 0, xend = NES, y = pathway_clean, yend = pathway_clean),
                 color = "grey60", linewidth = 0.4) +
    geom_point(aes(size = neg_log10_padj, fill = direction), shape = 21, stroke = 0.3) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_size_continuous(range = size_range, name = "-log10(padj)") +
    scale_x_continuous(limits = max(abs(sig_paths$NES)) * c(-1.05, 1.05)) +
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
