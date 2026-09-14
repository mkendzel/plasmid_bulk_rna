# Guards for plotting 25 contrasts, six of which have no significant genes.
#
# The contract these enforce: a missing figure file looks like a pipeline
# failure, an annotated empty panel is a result. Nothing below ever returns NULL
# in place of a plot.

# range() on an empty or all-non-finite vector returns c(Inf, -Inf), which kills
# coord_fixed / scale_*_continuous downstream.
safe_range <- function(x, default = c(-1, 1)) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(default)
  r <- range(x)
  if (diff(r) == 0) r + c(-1, 1) * max(abs(r[1]) * 0.05, 1e-6) else r
}

safe_max <- function(x, default = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) default else max(x)
}

# Symmetric limit for a diverging scale; never zero, so breaks stay valid.
safe_absmax <- function(x, default = 1) {
  m <- safe_max(abs(x), 0)
  if (m == 0) default else m
}

# cor() errors on fewer than 3 complete pairs and on a constant vector.
safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  stats::cor(x[ok], y[ok], method = method)
}

# t.test() errors below 2 observations and on two constant vectors. Ctrl2 is
# n = 1 at IFNb and IFNg, so any per-stimulus or Mock-only split can reach that
# case and must degrade to NA rather than abort a table.
safe_t_p <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 2L || length(b) < 2L) return(NA_real_)
  if (stats::sd(a) == 0 && stats::sd(b) == 0) return(NA_real_)
  tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_)
}

# One-sample version: is this vector's mean different from mu?
safe_t_p_one <- function(x, mu = 0) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || stats::sd(x) == 0) return(NA_real_)
  tryCatch(stats::t.test(x, mu = mu)$p.value, error = function(e) NA_real_)
}

# Placeholder drawn in place of a figure that has no data. Always written.
empty_plot <- function(title = NULL, subtitle = NULL, note = "No data to plot") {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = note,
                      size = 4, colour = "grey35") +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey35")
    )
}

# The checkpointed gsea_hallmark stores leadingEdge as one ";"-collapsed string
# per row. Gene symbols never contain ";", so the split is unambiguous. An empty
# leading edge is stored as "" (paste(character(0), collapse = ";")), hence the
# nzchar filters.
split_leading_edge <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(character(0))
  out <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  unique(out[nzchar(out)])
}

# TRUE when a matrix is safe to hand to pheatmap with clustering on.
heatmap_ready <- function(m, min_rows = 2L) {
  if (is.null(m) || !is.matrix(m) || nrow(m) < min_rows || ncol(m) < 2L) return(FALSE)
  ok <- apply(m, 1, function(r) {
    r <- r[is.finite(r)]
    length(r) > 1L && stats::sd(r) > 0
  })
  sum(ok) >= min_rows
}

# Drop rows that would become all-NaN under t(scale(t(.))).
drop_constant_rows <- function(m) {
  keep <- apply(m, 1, function(r) {
    r <- r[is.finite(r)]
    length(r) > 1L && stats::sd(r) > 0
  })
  m[keep, , drop = FALSE]
}

# Z-score by row, after removing rows that cannot be scaled.
zscore_rows <- function(m) {
  m <- drop_constant_rows(m)
  if (nrow(m) == 0L) return(m)
  t(scale(t(m)))
}
