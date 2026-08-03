# GSE140802 Larry MuTrans — TIPS STRING

When MuTrans can't handle >5000 cells effentially, we applied seacell to get 1200 metacells.

Note that in this datset, only total-count-normalized coutns are downloadable.  

Therefore, MuTrans pipeline skiped the Seurat::normalization line.

To benchmark MuTrans outcomes, we applied BioTIP on the saem metacell matrix, 

working on log2-transforemed 'data'  here.



**Start here.** Felix-style layout (mirrors `felix/TIPS/examples/GSE87038/GSE87038_STRING/code_core_13/`). Arm folders (`code_core_<TAG>/`) hold scripts and `run_TIPS_core.R`; workflow docs live in **this file only**.

```
7_data_MuTrans_TIPS_STRING/
  README.md              ← you are here
  ../TIPS_extend_v2.R    ← 24.1 dual-pull helper + lineage-lean reclaim (shared with IID)
  code_core_<TAG>/       # e.g. code_core_4_9vs11, code_core_4_10vs11, code_core_1_5vs7
  results_core_<TAG>/    # PPI_weight/, cisTarget_predicted_<CP>/, ...
  data/                  # shared caches (see below)
```

**IID variant:** sibling `7_data_MuTrans_TIPS_IID/` (Felix `GSE87038_IID` analogue).

## Comparison arms


| Folder               | CP           | CM      | CF     | Run                                                |
| -------------------- | ------------ | ------- | ------ | -------------------------------------------------- |
| `code_core_4_9vs11`  | 4 progenitor | 9 Baso  | 11 Meg | `source(".../code_core_4_9vs11/run_TIPS_core.R")`  |
| `code_core_4_10vs11` | 4            | 10 Mast | 11 Meg | `source(".../code_core_4_10vs11/run_TIPS_core.R")` |
| `code_core_1_5vs7`   | 1            | 5 Mono  | 7 Neut | `source(".../code_core_1_5vs7/run_TIPS_core.R")`   |


Edit **USER PARAMETERS** at the top of each `run_TIPS_core.R` (TAG, CP/CM/CF clusters, `motif_target_strategy`, etc.).

## Step order (`run_TIPS_core.R`)

`11.1` → `11.1.1` → `11.2.0` → `11.3` → `12.0` → `24.1` → `24.3` → `25_venn`

(`24.0` runs inside `24.1`.)

**Overlap guard:** `11.1`, `24.0`, and `24.1` stop if CTS has no overlap with either descendant HiG (CM/CF DEG).

## Recent changes (Holly / Larry MuTrans)



### Script naming — removed `_cardiac`

Larry MuTrans is not the cardiac study. Step scripts were renamed:


| Old                                         | Current                             |
| ------------------------------------------- | ----------------------------------- |
| `11.1_STRINGweighted_CTS_cardiac_network.R` | `11.1_STRINGweighted_CTS_network.R` |
| `11.3_CTS_cardiac_network_ANND_pagerank.R`  | `11.3_CTS_network_ANND_pagerank.R`  |




### STRING duplicate check (`11.1.1`) — STRING only

Duplicate vertex cleanup is **STRING-only**. The sibling IID pipeline skips `11.1.1` (IID gene symbols are already unique).

`11.1.1_check_vertex_duplication.R` resolves duplicate STRING peptide → gene symbol mappings. On the **first** successful run it calls biomaRt and writes a **shared cache** for all arms:

`data/unique_STRING_mapping_correction.txt`

When that file exists, re-runs skip biomaRt. Per-arm output: `results_core_<TAG>/correct_n_edges_HiG_STRING2.14.0.rds`.

### `TIPS_extend_v2.R` — `24.1` dual-pull helper + lineage-lean reclaim (all five Larry arms)

**Location:** `GSE140802_lineage_tracking/TIPS_extend_v2.R` (parent of `7_data_MuTrans_TIPS_STRING/` and `7_data_MuTrans_TIPS_IID/`). Sourced automatically from each `24.1_acat_CTS.cisTarget_dualpull_clean.R` via `tips_source_extend()` — this is now the **default** (`tips_extend_r_path()` in `00_configuration.R` resolves to `TIPS_extend_v2.R` unless `Sys.setenv(TIPS_EXTEND_R = "/path/to/file.R")` overrides it, e.g. to pin the older `TIPS_extend.R`, Part 1 only).

**Versioning:** `TIPS_extend_v2.R` is a superset — Part 1 (dual-pull helpers, below) is byte-identical to the original `TIPS_extend.R`; Part 2 (lineage-lean reclaim, below) is new. A future `TIPS_extend_v3.R` would follow the same pattern.

#### Part 1 — dual-pull candidate fallback (runs inside `24.1`, before `dualpull_final_table.tsv` exists)

Felix dual-pull uses `CM_candidate` / `CF_candidate` = cisTarget motif target **and** `CM_hi` / `CF_hi` (CTS gene also a DEG in the CM/CF cluster). Some arms (e.g. `code_core_1_5vs7`, CP1 → Mono vs Neut) have empty `CM_hi`/`CF_hi` on CTS rows, so auto PageRank TFs (e.g. ETS1/FOS/IRF7) produce empty pull subgraphs and `empty final_table`.

`TIPS_extend_v2.R` keeps Felix logic first and only adds fallbacks when needed:


| Step             | Arms that already succeeded (`4_9vs11`, `4_10vs11`, IID) | Sparse arms (`1_5vs7`-like)                                         |
| ---------------- | -------------------------------------------------------- | ------------------------------------------------------------------- |
| CM/CF candidates | Felix `CM_hi` & `CF_hi` (non-zero) → **same as before**  | Often zero → HiG fallback                                           |
| HiG fallback     | **Skipped** (not needed)                                 | Motif targets in `HiG_<CM>` / `HiG_<CF>` when `CM_hi`/`CF_hi` empty |
| PageRank re-pick | **Skipped** if current TFs already have candidates       | Walks PageRank when auto TFs lack dual-pull support                 |


**Do the other four completed pipelines need a rerun?** **No.** Their `results_core_`* outputs were not modified. On a future re-run of `24.1`, behavior should match the original Felix path when `CM_hi`/`CF_hi` are already populated (extend fallbacks do not run).

`24.1` **USER INPUT** (all arms):


| Setting                    | Default          | Purpose                                                           |
| -------------------------- | ---------------- | ----------------------------------------------------------------- |
| `tips_extend_auto_key_tfs` | `TRUE`           | Walk PageRank for TFs with dual-pull support when auto picks fail |
| `key_TFs_override`         | `NULL`           | Felix-style manual TF list (overrides auto + extend)              |
| `rebuild_heatmap`          | `TRUE` / `FALSE` | Set `TRUE` when changing TFs or after extend picks new key_TFs    |


**Only** `code_core_1_5vs7` **(if** `24.1` **failed with** `empty final_table`**):** sync `TIPS_extend_v2.R` + updated `24.1`, then:

```r
rebuild_heatmap <- TRUE
tips_extend_auto_key_tfs <- TRUE
source(".../code_core_1_5vs7/24.1_acat_CTS.cisTarget_dualpull_clean.R")
```

Expect `[TIPS_extend] skip ETS1 ...`, then `[TIPS_extend] key_TFs: ...` with different TFs and non-zero `vcount(g_CT_sub)`.

`motif_target_strategy` (`solo` / `fam` / `merge`) is unchanged — all five Larry arms use `merge` in `run_TIPS_core.R`; `merge` did not cause the `1_5vs7` failure.

#### Part 2 — lineage-lean gene modules (new, Holly 2026-07-29; runs AFTER `dualpull_final_table.tsv` exists)

`24.1`'s dual-pull output classifies every biased gene into **1 shared set** (bridge genes — on both CM and CF arms) and **2 lineage-lean sets** (CM-exclusive, CF-exclusive). Holly added a reclaim step: a shared/bridge gene is reassigned back to a lineage-lean set when its edges show a direction specific to one arm, matching how she manually curated the human (Felix) results. **This changes mouse results relative to the un-reclaimed sets** — expected, not a bug.

Two reclaim rules (both applied when `module_set = "lean"`, the default):

| Rule | Function | Logic |
|------|----------|-------|
| Opposite-direction bridge | `tips_reclaim_bridge()` | Gene-level: CM shows only `increase` and CF only `decrease` (or vice versa) → reclaim to the `increase` arm. Edge-level: a shared edge is CM-`increase`/CF-`decrease` (or vice versa) → reclaim the bridge gene(s) on that edge. |
| CF-increase majority | `tips_reclaim_cf_increase_only()` | For bridge genes with **no CM increase**: reclaim to CF if CF `increase` edges outnumber CF `decrease` edges. Captures endothelial-biased regulators (e.g. **HAND2, FOXA2**) while correctly excluding genes like **TAL1** (only one CF-increase edge among mostly-decrease edges — stays shared). |

**Entry points:**
- `tips_lineage_lean_genes(path, module_set)` — full detail: `CM`/`CF` gene vectors, `bridge_genes`, and a per-gene `gene_table` audit (`branch` + `membership` reason: `arm_exclusive` / `opposite_direction_reclaim` / `cf_increase_only_bridge`). `module_set`: `"lean"` (default, with reclaim), `"lean_exclusive_only"` (arm-exclusive only, no reclaim), or `"full"` (all lineage-biased-edge genes, no shared/lean split).
- `tips_lineage_lean_modules(path, module_set)` — thin wrapper mapping `CM → Blood`, `CF → Endothelium` for cross-atlas UCell scripts.

**Used by:** `HEP_cross_atlas_transfer_Ucell.R` (cross-atlas UCell scoring, not yet in this repo) and cross-atlas heatmaps. Not called from `24.1`/`24.3` itself — Part 2 is a separate post-processing step over `PPI_graph_GRN_prediction_CTS_<CP>_dualpull_final_table.tsv`.

## Expression / grouping

- Seurat: `.../larry_BioTIP/BioTIP_attractor/seu_attractor_MuTrans_HVG.rds`
- `load_larry_sce()` → SCE assay `logcounts` from Seurat RNA `data` layer (`group_col` → `sce$label`)
- BioTIP CTS: `.../BioTIP_attractor/CTS_Lib_data.RData` + `CTS_summary_data.csv` filter
- MuTrans metacells use `attractor` labels (not `leiden_r0_6`; that column is for single-cell BioTIP only)



## Larry vs Felix `code_core_13`

1. If seurat object existing, read normalziaed count and meta from seurat without saving redundant sce.rdata object
2. `TAG` + `tips_configure()` → `code_core_<TAG>/`, `results_core_<TAG>/`
3. Larry MuTrans Seurat + BioTIP CTS (not Felix `sce_E8.25`)
4. `db_species` 10090 → mm10 cisTarget + `motifAnnotations_mgi` (Felix GSE87038 example: hg38 + hgnc + `toupper`)
5. `11.1.1` — STRING ID de-duplication (cached in `data/unique_STRING_mapping_correction.txt`; not used for IID)
6. `12.0` — `seed_TF` auto from PageRank TFs in HiGCTS with BC > 0 (Felix: manual); falls back to CTS if no HiGCTS
7. `24.1` — top 3 PageRank TFs among top 20 in HiGCTS with cisTarget motifs; `motif_target_strategy` solo/fam/merge
8. Shared STRING / cisTarget / Felix `Shared_Data` paths
9. `Real_names` for 24.3 figure titles



## Felix `code_core_13` vs Holly full TIPS


|                        | Felix code_core_13                                    | Holly full TIPS                |
| ---------------------- | ----------------------------------------------------- | ------------------------------ |
| ATAC/ChIP              | None (dummy access = 1)                               | Real ATAC / ISL1 ChIP          |
| PPIN                   | STRING modular 11.1–11.3                              | Same idea, dataset-specific    |
| 11.3 rewiring          | Skipped (`p.PageRank = NA`) — slow, unused downstream | Full null testing              |
| seed_TF / key_TFs      | Manual in Felix examples                              | Holly auto-extends from motifs |
| cisTarget              | Ranking DB + `heatmap_blocked_CTS_*_v3.tsv`           | Often `*_scATAC_*` variants    |
| 11.1.1 duplicate check | STRING examples                                       | STRING only; IID skipped       |


Holly refs: `GSE175634_iPSC_CM_weighted_v9/24.1.1_...`, `IbarraSoria2018_E8.25_v9/24.1_...`

## Core inputs & shared paths


| Resource            | Path                                                                                                                                   |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Larry Seurat        | `F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/larry_BioTIP/BioTIP_attractor/seu_attractor_MuTrans_HVG.rds` |
| BioTIP CTS          | `.../larry_BioTIP/BioTIP_attractor/CTS_Lib_data.RData`                                                                                 |
| STRING PPIN         | `F:/projects/Share/PPIN/STRING/STRINGdb/mouse/full/`                                                                                   |
| cisTarget           | `F:/projects/Share/TF/cistarget/`                                                                                                      |
| Felix Shared_Data   | `F:/projects/TIPS/results/felix/TIPS/examples/Shared_Data/`                                                                            |
| Weinreb Table S3    | `F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/aaw3381-Weinreb-Table-S3.xlsx` (`25_venn` only)                                  |
| Felix GSE87038 PPIN | `F:/projects/TIPS/results/felix/TIPS/examples/GSE87038/data/PPIN/`                                                                     |




## `data/` folder (shared caches)


| File                                                      | From                                | Notes                                              |
| --------------------------------------------------------- | ----------------------------------- | -------------------------------------------------- |
| `sce_LARRY_mutrans_metacells.RData`                       | `00_prepare_inputs.R`               | optional                                           |
| `DEG_perState_*.rds`, `markers.up_ttest_min.prop0.25.rds` | `11.1_STRINGweighted_CTS_network.R` | per run                                            |
| `unique_STRING_mapping_correction.txt`                    | `11.1.1_check_vertex_duplication.R` | STRING only; first biomaRt run; shared across arms |




## Environment variables (optional overrides)


| Variable                     | Default                        | Purpose                                  |
| ---------------------------- | ------------------------------ | ---------------------------------------- |
| `TIPS_GROUP_COL`             | `attractor`                    | Seurat `meta.data` column for DEG groups |
| `TIPS_CP_CLUSTER`            | arm-specific (see table above) | CP attractor id                          |
| `TIPS_CM_CLUSTER`            | arm-specific                   | CM attractor id                          |
| `TIPS_CF_CLUSTER`            | arm-specific                   | CF attractor id                          |
| `TIPS_DB_SPECIES`            | `10090`                        | NCBI taxon (mouse)                       |
| `TIPS_MOTIF_TARGET_STRATEGY` | from runner                    | `solo`, `fam`, or `merge`                |
| `SEURAT_RDS`                 | Larry path above               | Override Seurat object                   |
| `STRING_PPIN_DIR`            | Share STRING path              | Override PPIN folder                     |




## Key outputs (per arm)

- `GSE140802_STRING_graph_perState_notsimplified.rds`
- `PPI_weight/GSE140802_STRING_graph_perState_simplified_combinedweighted.rds`
- `binary_annot_CTS_<CP>_DEG.tsv`
- `cisTarget_predicted_<CP>/PPI_graph_GRN_prediction_CTS_<CP>_dualpull_final_table.tsv`
- `venn_TD_CTS_HiG.pdf` (`25_venn`)

