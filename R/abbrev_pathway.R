# Applied in order to title-cased hallmark names: acronyms first, then words
PATHWAY_ABBREV <- c(
  "Epithelial Mesenchymal Transition" = "EMT",
  "Reactive Oxygen Species"           = "ROS",
  "Unfolded Protein Response"         = "UPR",
  "Oxidative Phosphorylation"         = "OxPhos",
  "Tnfa"                              = "TNFA",
  "Nfkb"                              = "NF-kB",
  "Il6"                               = "IL6",
  "Il2"                               = "IL2",
  "Jak"                               = "JAK",
  "Stat3"                             = "STAT3",
  "Stat5"                             = "STAT5",
  "E2f"                               = "E2F",
  "G2m"                               = "G2M",
  "Myc"                               = "MYC",
  "Kras"                              = "KRAS",
  "Mtorc1"                            = "mTORC1",
  "Pi3k"                              = "PI3K",
  "Akt"                               = "AKT",
  "Tgf"                               = "TGF",
  "Wnt"                               = "WNT",
  "Dna"                               = "DNA",
  "Uv"                                = "UV",
  "P53"                               = "p53",
  "Interferon"                        = "IFN",
  "Signaling"                         = "Sig.",
  "Response"                          = "Resp.",
  "Metabolism"                        = "Metab.",
  "Via"                               = "via"
)

# Helper to shorten hallmark pathway names for axis labels
abbrev_pathway <- function(x, max_chars = 30) {
  x <- str_to_title(str_replace_all(str_remove(x, "^HALLMARK_"), "_", " "))
  x <- str_replace_all(x, PATHWAY_ABBREV)
  str_trunc(x, max_chars)
}

# Hand-set labels that override abbrev_pathway() in plot titles and subtitles
PATHWAY_LABEL <- c(
  HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFN-α Response",
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-γ Response"
)

# Display name for a hallmark set: PATHWAY_LABEL if listed, else abbrev_pathway()
pathway_label <- function(x, max_chars = 60) {
  out <- abbrev_pathway(x, max_chars = max_chars)
  hit <- x %in% names(PATHWAY_LABEL)
  out[hit] <- unname(PATHWAY_LABEL[x[hit]])
  out
}

# Filename-safe slug for a hallmark set: "HALLMARK_IL6_JAK_STAT3_SIGNALING" ->
# "il6_jak_stat3_signaling"
pathway_slug <- function(x) {
  str_remove(x, "^HALLMARK_") |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_remove_all("^_|_$")
}
