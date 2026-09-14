# Treatment identity colours shared by every figure that colours by group.
# Control is dark enough to stay distinct from the grey75 used for
# non-significant points in make_volcano().
TREATMENT_COLORS <- c(
  Control = "#595959",
  CL      = "#D73027",
  DOPC    = "#4575B4"
)

# Named colour vector for a two-group contrast. Groups named in `palette` keep
# their identity colour; anything else falls back to blue for group_neg and red
# for group_pos.
contrast_colors <- function(group_neg, group_pos,
                            palette = TREATMENT_COLORS,
                            fallback = c("#4575B4", "#D73027")) {
  grps <- c(group_neg, group_pos)
  out  <- stats::setNames(fallback, grps)
  hit  <- grps %in% names(palette)
  out[hit] <- unname(palette[grps[hit]])
  out
}
