# Shared config for PROJECT_NAME. Sourced by scripts 1-5

library(dplyr)
library(tibble)
library(tidyr)

# ---- Paths ----
PROJ    <- "analysis/PROJECT_NAME"
dir_raw <- file.path(PROJ, "data", "plasmidsaurus")
dir_rds <- file.path(PROJ, "data", "r_objects")
dir_fig <- file.path(PROJ, "figures")
dir_res <- file.path(PROJ, "results")

# Plasmidsaurus delivery code, e.g. "LSS7T8" or "FCMSVB", given for each experiment
RUN_ID <- "RUNID"

# The results archive. Downloaded from Plasmidsaurus directly.
# Everything scripts 1 and 2 needs
# read via unz() and nothing is extracted into the repo.
results_zip <- file.path(dir_raw, paste0(RUN_ID, "_results.zip"))

# Members of results_zip.
#   expr          gene_id / gene_name / gene_biotype, then <vendor>_cpm and
#                 <vendor>_count per sample (raw values)
#   map_stats     one row per sample of run-level alignment numbers
#   biotype       reads per gene biotype, genes with 5+ reads
#   sample_stats  per-sample mapping TSV, Metric / Value pairs
zip_expr       <- paste0(RUN_ID, "-expression-matrix.tsv")
zip_map_stats  <- paste0(RUN_ID, "-mapping-stats-reads.csv")
zip_biotype    <- paste0(RUN_ID, "-gene-biotype-5plus_reads-summary.csv")
zip_sample_stats <- function(idx, code) {
  sprintf("%s_mapping-stats/%s_%d_%s.tsv", RUN_ID, RUN_ID, idx, code)
}

# Read one member of results_zip. `inner` is a path inside the archive, `reader`
# the parser for it (read.csv for the .csv members).
read_delivery <- function(inner, reader = utils::read.delim, ...) {
  if (!file.exists(results_zip)) {
    stop("Delivery archive not found: ", results_zip, call. = FALSE)
  }
  reader(unz(results_zip, inner), ...)
}

# idx -> vendor sample code, parsed from the per-sample filenames in results_zip.
# NULL when the archive is absent.
vendor_key <- function(zip = results_zip) {
  if (!file.exists(zip)) return(NULL)
  nm <- utils::unzip(zip, list = TRUE)$Name
  m  <- regmatches(nm, regexec(sprintf("mapping-stats/%s_(\\d+)_(.+)\\.tsv$", RUN_ID), nm))
  m  <- m[lengths(m) == 3L]
  if (length(m) == 0L) return(NULL)
  tibble::tibble(
    idx  = as.integer(vapply(m, `[`, character(1), 2L)),
    code = vapply(m, `[`, character(1), 3L)
  ) |>
    dplyr::arrange(idx)
}

# Create a directory and return its path. save_fig() calls this with dir_fig.
ensure_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# ---- Organism ----
# Mouse settings; the commented block is the human equivalent.
ORG_DB        <- "org.Mm.eg.db"
KEGG_ORG      <- "mmu"
GENE_MAP_PATH <- "genesets/mouse_gene_map_ensembl_symbol.rds"
HALLMARK_GMT  <- "genesets/mh.all.v2026.1.Mm.symbols.gmt"
MITO_PREFIX   <- "^mt-"   # matched against gene_name for the mito read share

# ORG_DB        <- "org.Hs.eg.db"
# KEGG_ORG      <- "hsa"
# GENE_MAP_PATH <- "genesets/human_gene_map_ensembl_symbol.rds"
# HALLMARK_GMT  <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
# MITO_PREFIX   <- "^MT-"

# Second gene set for script 4. NA runs Hallmark only.
CUSTOM_GMT <- NA_character_

# ---- Cutoffs ----
# Read by de_call.R, volcano_plot.R and scripts 3-5. LFC_CUT is on the log2 scale.
P_CUT   <- 0.05
LFC_CUT <- 1.5

# ---- Sample map example ----
# Vendor sample -> sample name and design factors.
#   vendor     sample code submitted to Plasmidsaurus
#   idx        position in the delivery; matrix columns are <RUN_ID>_<idx>
#   sample     name used downstream, "<treatment>_<tissue>_<rep>"
#   treatment  primary design factor
#   tissue     second design factor
#   rep        replicate number within treatment x tissue
sample_map <- tibble::tribble(
  ~vendor,      ~idx, ~sample,            ~treatment, ~tissue,  ~rep,
  "CTRL-L1",      1L, "Control_liver_1",  "Control",  "liver",    1L,
  "CTRL-L2",      2L, "Control_liver_2",  "Control",  "liver",    2L,
  "CTRL-L3",      3L, "Control_liver_3",  "Control",  "liver",    3L,
  "A-L1",         4L, "TreatA_liver_1",   "TreatA",   "liver",    1L,
  "A-L2",         5L, "TreatA_liver_2",   "TreatA",   "liver",    2L,
  "A-L3",         6L, "TreatA_liver_3",   "TreatA",   "liver",    3L,
  "B-L1",         7L, "TreatB_liver_1",   "TreatB",   "liver",    1L,
  "B-L2",         8L, "TreatB_liver_2",   "TreatB",   "liver",    2L,
  "B-L3",         9L, "TreatB_liver_3",   "TreatB",   "liver",    3L,
  "CTRL-S1",     10L, "Control_spleen_1", "Control",  "spleen",   1L,
  "CTRL-S2",     11L, "Control_spleen_2", "Control",  "spleen",   2L,
  "CTRL-S3",     12L, "Control_spleen_3", "Control",  "spleen",   3L,
  "A-S1",        13L, "TreatA_spleen_1",  "TreatA",   "spleen",   1L,
  "A-S2",        14L, "TreatA_spleen_2",  "TreatA",   "spleen",   2L,
  "A-S3",        15L, "TreatA_spleen_3",  "TreatA",   "spleen",   3L,
  "B-S1",        16L, "TreatB_spleen_1",  "TreatB",   "spleen",   1L,
  "B-S2",        17L, "TreatB_spleen_2",  "TreatB",   "spleen",   2L,
  "B-S3",        18L, "TreatB_spleen_3",  "TreatB",   "spleen",   3L
)

# ---- Factor levels ----
# REF_TREATMENT is the denominator of every treatment contrast.
# Script 3 fits each tissue separately; a single-tissue project lists one level.
TREATMENT_LEVELS <- c("Control", "TreatA", "TreatB")
TISSUE_LEVELS    <- c("liver", "spleen")
REF_TREATMENT    <- "Control"

sample_map$treatment <- factor(sample_map$treatment, levels = TREATMENT_LEVELS)
sample_map$tissue    <- factor(sample_map$tissue,    levels = TISSUE_LEVELS)

# Treatment x tissue group, used for replicate counts, QC and plot colours
sample_map$condition <- paste(sample_map$treatment, sample_map$tissue, sep = "_")

# Condition order follows the factor levels above
cond_levels <- sample_map |>
  dplyr::arrange(treatment, tissue) |>
  dplyr::pull(condition) |>
  unique()

sample_map$condition <- factor(sample_map$condition, levels = cond_levels)

# Replicates per group. Groups at n = 1 are flagged by min_n below, not dropped.
cond_n <- sample_map |>
  dplyr::count(condition, name = "n") |>
  tibble::deframe()

stopifnot(
  !any(duplicated(sample_map$vendor)),
  !any(duplicated(sample_map$idx)),
  !any(duplicated(sample_map$sample)),
  !anyNA(sample_map$treatment),
  !anyNA(sample_map$tissue),
  REF_TREATMENT %in% TREATMENT_LEVELS
)

# ---- Samples dropped after QC ----
# Names listed here leave the analysis at the top of script 3.
DROP_SAMPLES <- character(0)
# DROP_SAMPLES <- c("TreatB_liver_2")

# ---- Contrast registry ----
# One row per contrast, fitted within its tissue. `expr` goes to
# limma::makeContrasts against that tissue's cell-means design, whose columns
# are the treatment levels. Scripts 3-5 join on `name`.
#   name           unique id, used in file names
#   type           grouping for facets and annotation colours
#   tissue         the tissue fit the contrast is taken from
#   expr           the contrast in treatment-level arithmetic
#   label          title and axis text
#   min_n          smallest group size in the contrast
.non_ref <- setdiff(TREATMENT_LEVELS, REF_TREATMENT)

# Each treatment vs the reference treatment
.treat <- tidyr::expand_grid(tissue = TISSUE_LEVELS, treatment = .non_ref) |>
  dplyr::transmute(
    name          = sprintf("%s_vs_%s_%s", treatment, REF_TREATMENT, tissue),
    type          = "treatment",
    treatment     = treatment,
    ref_treatment = REF_TREATMENT,
    tissue        = tissue,
    expr          = sprintf("%s - %s", treatment, REF_TREATMENT),
    label         = sprintf("%s: %s vs %s", tissue, treatment, REF_TREATMENT)
  )

# Every pair of non-reference treatments
.pairs <- if (length(.non_ref) >= 2) {
  stats::setNames(
    as.data.frame(t(utils::combn(.non_ref, 2)), stringsAsFactors = FALSE),
    c("ref_treatment", "treatment")
  )
} else {
  data.frame(ref_treatment = character(0), treatment = character(0))
}

.between <- tidyr::expand_grid(tissue = TISSUE_LEVELS, .pairs) |>
  dplyr::transmute(
    name          = sprintf("%s_vs_%s_%s", treatment, ref_treatment, tissue),
    type          = "between_treatment",
    treatment     = treatment,
    ref_treatment = ref_treatment,
    tissue        = tissue,
    expr          = sprintf("%s - %s", treatment, ref_treatment),
    label         = sprintf("%s: %s vs %s", tissue, treatment, ref_treatment)
  )

contrast_registry <- dplyr::bind_rows(.treat, .between)

# Treatment levels named in a contrast expression
.groups_in_expr <- function(e) {
  toks <- unlist(strsplit(e, "[^A-Za-z0-9_]+"))
  intersect(toks, TREATMENT_LEVELS)
}

# Smallest group size per contrast, from replicate counts per condition
contrast_min_n <- function(registry, n = cond_n) {
  vapply(seq_len(nrow(registry)), function(i) {
    g <- paste(.groups_in_expr(registry$expr[i]), registry$tissue[i], sep = "_")
    min(n[g])
  }, integer(1))
}

contrast_registry$min_n <- contrast_min_n(contrast_registry)

contrast_registry <- contrast_registry |>
  dplyr::mutate(
    type      = factor(type, levels = c("treatment", "between_treatment")),
    treatment = factor(treatment, levels = TREATMENT_LEVELS),
    tissue    = factor(tissue, levels = TISSUE_LEVELS)
  ) |>
  dplyr::arrange(tissue, type, treatment, ref_treatment)

rm(.treat, .between, .pairs, .non_ref)

stopifnot(
  nrow(contrast_registry) > 0,
  !any(duplicated(contrast_registry$name)),
  !anyNA(contrast_registry$min_n)
)

# Contrast names for one tissue, in registry order
contrasts_for <- function(tis) {
  contrast_registry$name[contrast_registry$tissue == tis]
}

# Contrast labels for one tissue, in registry order. The `registry` default
# resolves against globalenv(), so a script that overwrites contrast_registry
# with the checkpoint gets the min_n recomputed in 3_DE_limma.R after DROP_SAMPLES.
contrast_labels <- function(tis, registry = contrast_registry) {
  registry |>
    dplyr::filter(tissue == tis) |>
    dplyr::arrange(type, treatment, ref_treatment) |>
    dplyr::pull(label) |>
    unique()
}

# Registry row for a contrast name; an all-NA row for an unknown name
reg_of <- function(cn, registry = contrast_registry) {
  registry[match(cn, registry$name), ]
}

# ---- Palettes ----
# Overrides TREATMENT_COLORS from R/treatment_colors.R, sourced before this file.
TREATMENT_COLORS <- c(
  Control = "#595959",
  TreatA  = "#D73027",
  TreatB  = "#4575B4"
)

tissue_shapes <- stats::setNames(c(21, 22, 23, 24, 25)[seq_along(TISSUE_LEVELS)],
                                 TISSUE_LEVELS)

type_cols <- c(treatment = "#4C72B0", between_treatment = "#DD8452")

# QC status tiles in script 2
qc_status_cols <- c(PASS = "#EFF3EF", WARN = "#FBE0A6", FAIL = "#F2B4AC", "NA" = "#E6E6E6")

# ---- QC thresholds ----
# Every check in script 2 reads warn/fail from this table. Deleting a row drops
# that check; a metric the script cannot compute is skipped.
#   direction  "higher" = larger is better, "lower" = smaller is better
#   warn/fail  absolute values in the metric's own unit
#   unit       "count" | "pct" | "ratio", formatting only
QC_THRESHOLDS <- tibble::tribble(
  ~metric,                ~label,                    ~category,    ~direction, ~warn,  ~fail,  ~unit,
  "input_reads",          "Read count",              "alignment",  "higher",   10e6,   5e6,    "count",
  "mapping_rate",         "Mapping rate",            "alignment",  "higher",   95,     90,     "pct",
  "dedup_rate",           "Dedup rate",              "alignment",  "higher",   35,     25,     "pct",
  "total_counts",         "Total counts",            "libsize",    "higher",   5e6,    2.5e6,  "count",
  "top100_frac",          "Counts in top 100 genes", "complexity", "lower",    40,     50,     "pct",
  "detected_genes",       "Detected genes",          "detection",  "higher",   14000,  10000,  "count",
  "genes_min10",          "Genes with >= 10 counts", "detection",  "higher",   10000,  7000,   "count",
  "protein_coding_pct",   "Protein-coding reads",    "biotype",    "higher",   85,     80,     "pct",
  "rRNA_pct",             "rRNA reads",              "biotype",    "lower",    5,      10,     "pct",
  "mito_pct",             "Mitochondrial reads",     "biotype",    "lower",    10,     15,     "pct",
  "within_condition_cor", "Within-condition r",      "clustering", "higher",   0.95,   0.92,   "ratio",
  "pca_centroid_dist",    "PCA centroid distance",   "pca",        "lower",    5,      8,      "ratio",
  "median_abs_M",         "Replicate deviation",     "replicate",  "lower",    0.30,   0.45,   "ratio"
)

# Non-threshold QC knobs
QC_PARAMS <- list(
  filter_min_count = 10,   # gene kept when it has this many counts ...
  filter_min_rep   = 3,    # ... in at least this many samples, before vst
  complexity_top_n = 100,  # genes summed for top100_frac
  detection_min    = 10,   # count threshold for genes_min10
  pca_n_feats      = 500,  # top-variable genes for PCA and clustering
  pca_n_pcs        = 5     # PCs in the scatter matrix
)
