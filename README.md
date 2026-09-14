# plasmid_bulk_rna

Template pipeline for a bulk RNA-seq project built from a Plasmidsaurus
delivery. limma-voom throughout, with every project-specific value collected
in one config file. See [template_scripts/README.md](template_scripts/README.md)
for the full walkthrough (what each script reads/writes, and what to fill in).

## Quick start

```powershell
Copy-Item -Recurse analysis\_template "analysis\YOUR_PROJECT"
Copy-Item template_scripts\*.R "analysis\YOUR_PROJECT\scripts"
```

Put the delivery archive in `analysis/YOUR_PROJECT/data/plasmidsaurus/`, edit
that project's `0_config.R`, then run scripts 1-5 in order from the repo root.

## Layout

- `template_scripts/` — the numbered pipeline scripts (`0_config.R` through
  `5_Graphs.R`) and the shell script for alignment/quantification upstream of
  the delivery (FastQC, STAR, featureCounts).
- `R/` — project-agnostic helper functions sourced by every script.
- `genesets/` — Ensembl-to-symbol gene ID maps and Hallmark gene sets (mouse +
  human) used by the default config.
- `analysis/_template/` — empty project skeleton (`data/`, `figures/`,
  `results/`, `scripts/`) copied per new project.
