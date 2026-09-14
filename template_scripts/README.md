# Template scripts

Starting point for a bulk RNA-seq project built from a Plasmidsaurus delivery.
limma-voom throughout, with every project-specific value in `0_config.R`.

## Setting up a project

```powershell
Copy-Item -Recurse analysis\_template "analysis\YOUR_PROJECT"
Copy-Item template_scripts\*.R "analysis\YOUR_PROJECT\scripts"
```

Put the delivery archive in `analysis/YOUR_PROJECT/data/plasmidsaurus/`, then in
each script change the `source()` line to the project's own `0_config.R`.

All scripts run from the repo root, in order.

| Script | Reads | Writes |
|---|---|---|
| `0_config.R` | - | nothing; sourced by 1-5 |
| `1_import.R` | `<RUN_ID>_results.zip` | `counts_raw`, `vendor_genes` |
| `2_QC.R` | `counts_raw`, delivery archive | `figures/qc/`, `results/qc_report.md` |
| `3_DE_limma.R` | `counts_raw` | `dge`, `voom`, `fit`, `logcpm`, `sample_meta` (lists by tissue), `tt_annotated` (list by contrast), `results/de_*` |
| `4_GSEA.R` | `tt_annotated` | `gsea_hallmark`, `gsea_custom`, `kegg`, `results/*.tsv`, `figures/pathways/` |
| `5_Graphs.R` | `tt_annotated`, `logcpm`, `sample_meta`, `gsea_hallmark` | `figures/overview/`, `figures/de/`, `figures/gsea/`, `figures/genes/` |

Checkpoints are versioned RDS files in `data/r_objects/`, written by
`save_checkpoint()` and read by `load_checkpoint()`; each one gets a `.md`
sidecar recording the code that produced it.

## What to fill in

`0_config.R`, in order:

1. `PROJ` and `RUN_ID`.
2. Organism block: mouse is active, human is commented below it.
3. `sample_map` — one row per sample, transcribed from the submission sheet.
   `idx` is the position in the delivery; matrix columns are `<RUN_ID>_<idx>`.
4. `TREATMENT_LEVELS`, `TISSUE_LEVELS`, `REF_TREATMENT`.
5. `TREATMENT_COLORS` — one colour per treatment level.
6. `P_CUT`, `LFC_CUT`, `QC_THRESHOLDS`.
7. `DROP_SAMPLES`, after reading `results/qc_report.md`.

## Design

Each tissue in `TISSUE_LEVELS` is fitted on its own: samples are subset, then
get their own `filterByExpr`, TMM factors, voom trend and `eBayes` fit. The
design is cell-means on treatment (`~ 0 + treatment`). A single-tissue project
lists one level in `TISSUE_LEVELS`.

Contrasts are rows of `contrast_registry`. Each row names the tissue it is
fitted in, and its `expr` column is treatment-level arithmetic handed to
`limma::makeContrasts`:

| type | pattern |
|---|---|
| `treatment` | `<treatment> - <ref>` |
| `between_treatment` | `<treatment> - <other>` |

The two blocks that build it are a default. Rows can be added, deleted or
hand-written; scripts 3-5 loop over whatever the registry holds.

`condition` (`<treatment>_<tissue>`) is the replicate group used for `min_n`,
the QC metrics in script 2 and plot colours.

## Shared code

Helpers in `R/` are sourced by every script and are project-agnostic:
`checkpoint.R`, `save_fig.R`, `de_call.R` (the Up/Down/NS call), `plot_guards.R`
(empty-data guards), `volcano_plot.R`, `gene_ids.R`, `run_gsea.R`,
`write_qc_report.R`, `abbrev_pathway.R`.

Gene sets and Ensembl-symbol maps live in `genesets/` and are referenced from
scripts as `genesets/...`.
