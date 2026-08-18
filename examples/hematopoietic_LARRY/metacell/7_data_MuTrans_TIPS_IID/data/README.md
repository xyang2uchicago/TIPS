# Cached outputs (Felix `data/` slot)

**Expression input:** read directly from Larry Seurat each run — no SCE file required.

`F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/larry_BioTIP/BioTIP_attractor/seu_attractor_MuTrans_HVG.rds`

This folder only stores **11.1 DEG outputs** (filenames include cutoffs):

- `DEG_perState_min.prop{min.prop}_lfc{logFC}_FDFR{FDR}.rds` — filtered gene lists per attractor
- `markers.up_ttest_min.prop{min.prop}_lfc{logFC}_FDFR{FDR}.rds` — full `findMarkers` tables

Defaults: `min.prop=0.25`, `lfc=0.5`, `FDR=0.05` — set **`logFC.cut` and `fdr.cut` in `11.1`** (USER INPUT block). Downstream steps read `results_core_4/deg_cutoffs_latest.rds`.

SCE helpers (`load_larry_sce`, `build_mutrans_sce`, …) live in `GSE140802_lineage_tracking/TIPS_extend.R` and are loaded via `tips_configure()`.
