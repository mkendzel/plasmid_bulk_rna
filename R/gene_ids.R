# Ensembl <-> symbol plumbing shared by every figure that labels genes.
#
# logcpm rownames carry Ensembl versions ("ENSG00000000003.15"); the annotated
# toptables carry the stripped form in `ensembl_id`. Everything here works in the
# stripped space.

strip_ens_version <- function(x) sub("\\..*$", "", x)

# ensembl_id -> SYMBOL, one row per Ensembl ID. Accepts a single annotated
# toptable or a list of them (the first element is used - annotation is identical
# across contrasts of one experiment).
symbol_map <- function(tt) {
  if (is.list(tt) && !is.data.frame(tt)) tt <- tt[[1L]]
  tt |>
    dplyr::select(ensembl_id, SYMBOL) |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
    dplyr::distinct(ensembl_id, .keep_all = TRUE)
}

# Row labels for a matrix keyed by (possibly versioned) Ensembl IDs. Falls back
# to the Ensembl ID where no symbol maps, and is always unique so it can be
# assigned to rownames().
label_genes <- function(ids, map) {
  s   <- map$SYMBOL[match(strip_ens_version(ids), map$ensembl_id)]
  bad <- is.na(s) | s == ""
  s[bad] <- ids[bad]
  make.unique(s)
}
