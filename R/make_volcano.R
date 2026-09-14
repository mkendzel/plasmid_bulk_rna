# General volcano plots

# Volcano of one topTable. Grey = not in pathway_name, blue = in it,
# red = in it and past the dashed cutoffs. Labels are thinned.
# label_size: text size of the repelled gene labels.
make_pathway_volcano <- function(res_df, pathway_name, set_label, contrast_label, tissue,
                                 group_neg, group_pos, n_bins = 8,
                                 pathway_list = pathways,
                                 padj_cut = padj_cutoff, lfc_cut = lfc_cutoff,
                                 label_size = 3) {
  target_genes <- pathway_list[[pathway_name]]

  volcano_df <- res_df %>%
    filter(!is.na(SYMBOL), !is.na(logFC), !is.na(adj.P.Val)) %>%
    mutate(
      neg_log10_padj = pmin(-log10(adj.P.Val), 50),
      is_target = SYMBOL %in% target_genes,
      sig = is_target & adj.P.Val < padj_cut & abs(logFC) > lfc_cut
    )

  n_detected <- n_distinct(volcano_df$SYMBOL[volcano_df$is_target])
  n_sig      <- n_distinct(volcano_df$SYMBOL[volcano_df$sig])
  x_limit    <- max(abs(volcano_df$logFC)) * 1.05

  ggplot(volcano_df, aes(x = logFC, y = neg_log10_padj)) +
    geom_point(data = filter(volcano_df, !is_target),
               color = "grey75", size = 0.5, alpha = 0.4) +
    geom_point(data = filter(volcano_df, is_target & !sig),
               color = "steelblue", size = 1.5, alpha = 0.7) +
    geom_point(data = filter(volcano_df, sig),
               color = "firebrick", size = 2, alpha = 0.9) +
    geom_text_repel(
      data = label_set(volcano_df, n_bins),
      aes(label = SYMBOL),
      size = label_size, max.overlaps = Inf, seed = 42,
      color = "#7F0000", fontface = "bold",
      bg.color = "white", bg.r = 0.15,
      box.padding = 0.7, point.padding = 0.2, force = 5,
      min.segment.length = 0, segment.color = "grey50", segment.size = 0.3
    ) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey40") +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    labs(
      x = direction_axis_label("log2 Fold Change", group_neg, group_pos),
      y = "-log10(adj.P.Val)",
      title = paste0(contrast_label, " (", tissue, ")"),
      subtitle = paste0(set_label, " — ", n_sig, "/", n_detected, " set genes significant")
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"))
}

# Volcano of every gene in a topTable, coloured by direction. Same cutoffs,
# labels, and y cap as make_ifn_volcano.
make_volcano <- function(res_df, contrast_label, tissue, group_neg, group_pos,
                         n_bins = 8,
                         padj_cut = padj_cutoff, lfc_cut = lfc_cutoff) {
  volcano_df <- res_df %>%
    filter(!is.na(SYMBOL), !is.na(logFC), !is.na(adj.P.Val)) %>%
    mutate(
      neg_log10_padj = pmin(-log10(adj.P.Val), 50),
      sig = adj.P.Val < padj_cut & abs(logFC) > lfc_cut,
      direction = factor(
        case_when(sig & logFC > 0 ~ group_pos,
                  sig & logFC < 0 ~ group_neg,
                  TRUE ~ "Not significant"),
        levels = c(group_neg, group_pos, "Not significant")
      )
    )

  n_up    <- n_distinct(volcano_df$SYMBOL[volcano_df$sig & volcano_df$logFC > 0])
  n_down  <- n_distinct(volcano_df$SYMBOL[volcano_df$sig & volcano_df$logFC < 0])
  x_limit <- max(abs(volcano_df$logFC)) * 1.05

  color_values <- c(contrast_colors(group_neg, group_pos),
                    "Not significant" = "grey75")

  ggplot(volcano_df, aes(x = logFC, y = neg_log10_padj)) +
    geom_point(data = filter(volcano_df, !sig),
               aes(color = direction), size = 0.5, alpha = 0.4) +
    geom_point(data = filter(volcano_df, sig),
               aes(color = direction), size = 1.5, alpha = 0.8) +
    geom_text_repel(
      data = label_set(volcano_df, n_bins),
      aes(label = SYMBOL, color = direction),
      size = 3, max.overlaps = Inf, segment.color = "grey40",
      show.legend = FALSE, seed = 42
    ) +
    scale_color_manual(values = color_values, drop = FALSE) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey40") +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    labs(
      x = direction_axis_label("log2 Fold Change", group_neg, group_pos),
      y = "-log10(adj.P.Val)",
      color = "Higher in",
      title = paste0(contrast_label, " (", tissue, ")"),
      subtitle = paste0("All genes — ", n_up, " higher in ", group_pos, ", ",
                        n_down, " higher in ", group_neg,
                        " (adj.P < ", padj_cut, ", |logFC| > ", lfc_cut, ")")
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"))
}
