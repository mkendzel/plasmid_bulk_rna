# plasmid_bulk_rna

Template pipeline for a bulk RNA-seq project sequenced through Plasmidsaurus.
limma-voom throughout, with every project-specific value collected in one
config file. See [template_scripts/README.md](template_scripts/README.md) for
the full walkthrough (what each script reads/writes, and what to fill in).

## Quick start

1. Copy the project skeleton and scripts:

   ```powershell
   Copy-Item -Recurse analysis\_template "analysis\YOUR_PROJECT"
   Copy-Item template_scripts\*.R "analysis\YOUR_PROJECT\scripts"
   ```

2. Download the **results archive** from Plasmidsaurus (a `.zip` file named
   with your run ID, e.g. `LSS7T8_results.zip`) and drop it, unmodified, into
   `analysis/YOUR_PROJECT/data/plasmidsaurus/`. The scripts read directly out
   of this zip — nothing needs to be extracted.

3. Open that project's `0_config.R` and edit it for your experiment:
   - Set `PROJ` and `RUN_ID` (the run ID in the zip's filename).
   - **Fill in `sample_map`** — one row per sample, transcribed from your
     Plasmidsaurus submission sheet (vendor sample code, its position `idx`
     in the delivery, your sample name, and its `treatment`/`tissue`/`rep`).
     This is the one step the archive can't do for you: it tells the
     pipeline which vendor code is which sample.
   - Update `TREATMENT_LEVELS`, `TISSUE_LEVELS`, `REF_TREATMENT`, and the
     other settings below `sample_map` to match.

4. Run scripts 1-5 in order, from the repo root.

## Layout

- `template_scripts/` — the numbered pipeline scripts (`0_config.R` through
  `5_Graphs.R`) and a shell script for aligning/quantifying raw FASTQs
  yourself (FastQC, STAR, featureCounts) — not needed if you're starting from
  the Plasmidsaurus results archive, since they already ran this for you.
- `R/` — project-agnostic helper functions sourced by every script.
- `genesets/` — Ensembl-to-symbol gene ID maps and Hallmark gene sets (mouse +
  human) used by the default config.
- `analysis/_template/` — empty project skeleton (`data/`, `figures/`,
  `results/`, `scripts/`) copied per new project.
