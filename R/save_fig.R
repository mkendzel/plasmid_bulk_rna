# The one figure writer. Everything in scripts 3 and 5-7 goes through this.
#
# PNG is the default: the pipeline emits ~160 figures and a 300 dpi LZW TIFF runs
# 3-8 MB apiece. Pass formats = FIG_FORMATS_PUB at the handful of call sites that
# are actual manuscript panels.
#
# bg = "white" is not optional. ggsave() defaults to a transparent background for
# themes that do not paint plot.background, which is why 2_QC.R grew its own
# png_out() wrapper; that wrapper is retired in favour of this.

FIG_FORMATS     <- c("png")
FIG_FORMATS_PUB <- c("pdf", "tiff")

save_fig <- function(p, name, width = 7, height = 6, subdir = NULL,
                     formats = FIG_FORMATS, dpi = 300, bg = "white") {

  outdir <- if (is.null(subdir)) ensure_dir(dir_fig) else ensure_dir(dir_fig, subdir)

  files <- vapply(formats, function(fmt) {
    f <- file.path(outdir, paste0(name, ".", fmt))

    switch(fmt,
      png  = grDevices::png(f, width = width, height = height, units = "in",
                            res = dpi, bg = bg),
      tiff = grDevices::tiff(f, width = width, height = height, units = "in",
                             res = dpi, compression = "lzw", bg = bg),
      pdf  = grDevices::pdf(f, width = width, height = height, bg = bg),
      stop("save_fig: unsupported format '", fmt, "'")
    )
    on.exit(grDevices::dev.off(), add = TRUE)

    # pheatmap and UpSetR draw to the device rather than returning a printable
    # ggplot, so they need their own branches.
    if (inherits(p, "pheatmap")) {
      grid::grid.draw(p$gtable)
    } else if (inherits(p, "upset")) {
      print(p)
    } else if (is.function(p)) {
      p()
    } else {
      print(p)
    }
    f
  }, character(1))

  invisible(unname(files))
}
