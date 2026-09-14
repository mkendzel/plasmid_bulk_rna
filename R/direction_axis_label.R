# X-axis label with arrows marking which group each side is higher in.
direction_axis_label <- function(value_label, group_neg, group_pos) {
  paste0(value_label, "\n← Higher in ", group_neg, "   |   Higher in ", group_pos, " →")
}
