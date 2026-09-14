# Resolving CollecTRI regulons named with a UniProt accession ("A0A087WSP5")
# to their gene symbol ("Stat1").

# UniProtKB accession syntax.
is_uniprot_acc <- function(x) {
  grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", x)
}

# accession -> gene symbol, NA where nothing resolves. org.Mm.eg.db is queried
# first, then the UniProt REST service for whatever is left. `online = FALSE`
# and a failed request both skip the REST pass.
uniprot_symbol_map <- function(acc, online = TRUE) {
  acc <- unique(acc[!is.na(acc) & nzchar(acc)])
  out <- stats::setNames(rep(NA_character_, length(acc)), acc)
  if (length(acc) == 0L) return(out)

  local <- suppressMessages(tryCatch(
    AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys = acc,
                          column = "SYMBOL", keytype = "UNIPROT",
                          multiVals = "first"),
    error = function(e) stats::setNames(rep(NA_character_, length(acc)), acc)
  ))
  out[names(local)] <- unname(local)

  todo <- names(out)[is.na(out)]
  if (!online || length(todo) == 0L) return(out)

  # Blocks of 100 accessions per request
  for (chunk in split(todo, ceiling(seq_along(todo) / 100))) {
    q <- paste0(
      "https://rest.uniprot.org/uniprotkb/search?format=tsv&size=500",
      "&fields=accession,gene_primary",
      "&query=", paste0("accession:", chunk, collapse = "+OR+")
    )
    tsv <- tryCatch(utils::read.delim(url(q), stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (is.null(tsv) || !nrow(tsv)) next
    hit <- tsv[nzchar(tsv$Gene.Names..primary.), ]
    out[hit$Entry] <- hit$Gene.Names..primary.
  }

  out
}

# Swap accessions for symbols. Anything absent from `lookup`, and any accession
# that resolved to NA, keeps the label it came in with.
relabel_sources <- function(x, lookup) {
  hit <- !is.na(lookup[x]) & x %in% names(lookup)
  x[hit] <- unname(lookup[x[hit]])
  x
}

# Apply relabel_sources() to an activity table's `source` column. Warns when
# two regulons in one contrast collapse onto the same label.
relabel_activity <- function(act, lookup) {
  act$source <- relabel_sources(act$source, lookup)

  per_contrast <- act[c("source", "contrast")]
  dup <- unique(act$source[duplicated(per_contrast)])
  if (length(dup)) {
    warning("relabel_activity: label collision on ", paste(dup, collapse = ", "),
            call. = FALSE)
  }

  act
}
