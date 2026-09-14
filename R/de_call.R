# The Up / Down / NS call, in one place.
#
# This decision used to appear in five spots with three different spellings
# (`> LFC_CUT` in the volcano, `>= LFC_CUT` in the heatmaps and Venn), so the
# counts annotated on a volcano disagreed with the counts in de_gene_counts.tsv.
# `>` is authoritative because that is what 3_DE_limma.R writes to the TSV.
#
# P_CUT and LFC_CUT come from 0_config.R; LFC_CUT is on the log2 scale.

DE_LEVELS <- c("Up", "Down", "NS")
DE_COLS   <- c(Up = "#D62728", Down = "#1F77B4", NS = "grey75")

# TRUE where a gene clears both cutoffs. NA logFC or NA adjusted p is FALSE.
is_de <- function(logFC, adj_p, p_cut = P_CUT, lfc_cut = LFC_CUT) {
  !is.na(logFC) & !is.na(adj_p) & adj_p < p_cut & abs(logFC) > lfc_cut
}

# Factor with fixed levels, so scale_*_manual(drop = FALSE) keeps a legend entry
# even for a contrast where nothing is significant.
de_direction <- function(logFC, adj_p, p_cut = P_CUT, lfc_cut = LFC_CUT) {
  sig <- is_de(logFC, adj_p, p_cut, lfc_cut)
  factor(ifelse(sig & logFC > 0, "Up", ifelse(sig, "Down", "NS")), levels = DE_LEVELS)
}

# Significant gene IDs from an annotated toptable, one direction or both.
de_genes <- function(df, direction = c("both", "Up", "Down"), id_col = "ensembl_id",
                     p_cut = P_CUT, lfc_cut = LFC_CUT) {
  direction <- match.arg(direction)
  d <- de_direction(df$logFC, df$adj.P.Val, p_cut, lfc_cut)
  keep <- if (direction == "both") d != "NS" else d == direction
  unique(df[[id_col]][which(keep)])
}
