# Single-file PDF writer.
#
# Separate from save_fig()
# always uses cairo_pdf: allows for non-ASCII glyphs and custom arrows
save_pdf <- function(plot, filename, width = 7, height = 5, dir = fig_dir) {
  ggsave(file.path(dir, filename), plot,
         width = width, height = height, device = cairo_pdf)
}
