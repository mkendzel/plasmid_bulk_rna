# Label thinning functions for crowded volcano panels.

# Keep the top-ranked row per cell of an n_bins x n_bins grid over the plotted
# range. Crowded regions get thinned to one label, isolated points survive.
thin_labels <- function(df, x, y, rank_by, n_bins = 8) {
  if (nrow(df) < 2) return(df)

  df %>%
    mutate(
      .x_bin = cut({{ x }}, breaks = n_bins, labels = FALSE),
      .y_bin = cut({{ y }}, breaks = n_bins, labels = FALSE)
    ) %>%
    group_by(.x_bin, .y_bin) %>%
    slice_max({{ rank_by }}, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-.x_bin, -.y_bin)
}

# Candidate labels for a volcano: significant genes, one row per symbol
# Rank favours genes that are both large-effect and confident.
# Expects the columns make_volcano() / make_pathway_volcano() 
#    build: sig, neg_log10_padj, logFC, SYMBOL.
label_set <- function(volcano_df, n_bins) {
  volcano_df %>%
    filter(sig) %>%
    mutate(label_rank = neg_log10_padj * abs(logFC)) %>%
    arrange(desc(label_rank)) %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    thin_labels(logFC, neg_log10_padj, label_rank, n_bins = n_bins)
}
