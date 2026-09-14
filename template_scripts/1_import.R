# Import the Plasmidsaurus expression matrix and rename columns to sample names.
# Nothing is filtered or dropped. Run once per delivery, from the repo root.

# ---- Libraries ----
library(dplyr)
library(tibble)

# ---- Load helper functions + project config ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("template_scripts/0_config.R")

# ---- Import ----
expr <- read_delivery(
  zip_expr,
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

# Gene annotation shipped with the matrix; used for biotype QC and symbol fallback
vendor_genes <- expr[, intersect(c("gene_id", "gene_name", "gene_biotype"),
                                 colnames(expr))]

# ---- Count matrix ----
count_cols <- grep("_count$", colnames(expr), value = TRUE)

counts <- expr[, c("gene_id", count_cols)]
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

colnames(counts) <- sub("_count$", "", colnames(counts))

# ---- Rename columns to sample names ----
# Three keys, tried in order: the vendor idx -> code map parsed from the archive,
# the submitted sample code, then the trailing index of the column name. The
# vendor map is checked against sample_map first, so a renumbered delivery errors
# instead of relabelling.
vk <- vendor_key()

rename_by_vendor <- character(0)
if (!is.null(vk)) {
  chk <- dplyr::inner_join(vk, sample_map[, c("idx", "vendor", "sample")],
                           by = "idx")
  stopifnot(
    nrow(chk) == nrow(sample_map),
    identical(chk$code, chk$vendor)
  )
  rename_by_vendor <- stats::setNames(chk$sample, paste0(RUN_ID, "_", chk$idx))
}

rename_by_code <- stats::setNames(sample_map$sample, sample_map$vendor)
rename_by_idx  <- stats::setNames(sample_map$sample, as.character(sample_map$idx))

new_names <- vapply(colnames(counts), function(cn) {
  if (cn %in% names(rename_by_vendor)) return(rename_by_vendor[[cn]])
  if (cn %in% names(rename_by_code))   return(rename_by_code[[cn]])
  idx <- sub("^.*_", "", cn)
  if (grepl("^[0-9]+$", idx) && idx %in% names(rename_by_idx)) return(rename_by_idx[[idx]])
  NA_character_
}, character(1), USE.NAMES = FALSE)

# List every failure before stopping
if (anyNA(new_names)) {
  message("Unmapped count columns:\n  ",
          paste(colnames(counts)[is.na(new_names)], collapse = "\n  "))
}
colnames(counts) <- new_names

stopifnot(
  !anyNA(colnames(counts)),
  !any(duplicated(colnames(counts))),
  length(setdiff(sample_map$sample, colnames(counts))) == 0
)

counts <- counts[, sample_map$sample, drop = FALSE]

# ---- Checkpoints ----
save_checkpoint(counts, "counts_raw", dir = dir_rds,
                notes = "all samples, vendor columns renamed, nothing filtered")
save_checkpoint(vendor_genes, "vendor_genes", dir = dir_rds,
                notes = "gene_id / gene_name / gene_biotype from the vendor matrix")
