# Volcano panel for one annotated toptable.
#
# Expects limma's columns (logFC, adj.P.Val) plus SYMBOL. The Up / Down / NS call
# comes from de_call.R, so the counts drawn on the panel match de_gene_counts.tsv.

volcano_plot <- function(df, title = NULL, subtitle = NULL, caption = NULL,
                         label_n = 0) {

  d <- dplyr::mutate(df,
    neglog10_padj = -log10(adj.P.Val),
    direction     = de_direction(logFC, adj.P.Val)
  )

  n_up   <- sum(d$direction == "Up",   na.rm = TRUE)
  n_down <- sum(d$direction == "Down", na.rm = TRUE)

  p <- ggplot2::ggplot(d, ggplot2::aes(logFC, neglog10_padj)) +
    ggplot2::geom_point(ggplot2::aes(colour = direction), alpha = 0.6, size = 1) +
    ggplot2::geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(P_CUT),
                        linetype = "dashed", colour = "grey40") +
    # Inf placement keeps the counts on-panel when the y axis is near-degenerate
    ggplot2::annotate("text", x = Inf, y = Inf, label = paste0("Up: ", n_up),
                      hjust = 1.1, vjust = 1.5, colour = DE_COLS[["Up"]], size = 3.2) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = paste0("Down: ", n_down),
                      hjust = -0.1, vjust = 1.5, colour = DE_COLS[["Down"]], size = 3.2) +
    ggplot2::scale_colour_manual(values = DE_COLS, drop = FALSE, name = NULL) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = "log2 fold change", y = "-log10 adjusted p-value") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.caption = ggplot2::element_text(colour = "#B22222", hjust = 0))

  if (label_n > 0) {
    lab <- d |>
      dplyr::filter(direction != "NS", !is.na(SYMBOL), SYMBOL != "") |>
      dplyr::slice_min(adj.P.Val, n = label_n, with_ties = FALSE)
    if (nrow(lab) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = lab, ggplot2::aes(label = SYMBOL), size = 2.6, fontface = "italic",
        max.overlaps = 20, min.segment.length = 0, segment.size = 0.3,
        segment.colour = "grey60")
    }
  }
  p
}
