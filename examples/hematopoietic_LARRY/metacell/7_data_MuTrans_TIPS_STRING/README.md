# GSE140802 Larry MuTrans — TIPS STRING

When MuTrans can't handle >5000 cells effentially, we applied seacell to get 1200 metacells.

Note that in this datset, only total-count-normalized coutns are downloadable.  

Therefore, MuTrans pipeline skiped the Seurat::normalization line.

To benchmark MuTrans outcomes, we applied BioTIP on the saem metacell matrix, 

working on log2-transforemed 'data'  here.



**Start here.** Felix-style layout (mirrors `felix/TIPS/examples/GSE87038/GSE87038_STRING/code_core_13/`). Arm folders (`code_core_<TAG>/`) hold scripts and `run_TIPS_core.R`; workflow docs live in **this file only**.

```
7_data_MuTrans_TIPS_STRING/
  README.md                ← you are here
  ../../TIPS_extend_v2.R   ← 24.1 dual-pull helper (shared with IID and single_cell/)
  code_core_4_9vs11/
  results_core_4_9vs11/    # PPI_weight/, cisTarget_predicted_<CP>/, ...
  data/                    # shared caches (see below)
```

**IID variant:** sibling `7_data_MuTrans_TIPS_IID/` (Felix `GSE87038_IID` analogue).

## Comparison arms

Only one arm is kept in this repo; the other comparison arms (`code_core_4_10vs11`, `code_core_4_9vs10`, `code_core_6_9vs10`, `code_core_1_5vs7`) were pruned as unnecessary. The table below documents the arm that remains; the tuning notes further down that reference the removed arms (e.g. `code_core_1_5vs7`) are kept for background on why `TIPS_extend.R`/`TIPS_extend_v2.R` has its HiG-fallback and PageRank re-pick logic.

| Folder               | CP           | CM      | CF     | Run                                                |
| -------------------- | ------------ | ------- | ------ | -------------------------------------------------- |
| `code_core_4_9vs11`  | 4 progenitor | 9 Baso  | 11 Meg | `source(".../code_core_4_9vs11/run_TIPS_core.R")`  |

Edit **USER PARAMETERS** at the top of `run_TIPS_core.R` (TAG, CP/CM/CF clusters, `motif_target_strategy`, etc.).

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

### `TIPS_extend_v2.R` — `24.1` dual-pull helper (shared across metacell and single_cell)

**Location:** `hematopoietic_LARRY/TIPS_extend_v2.R` (top level, sibling to both `metacell/` and `single_cell/`). `TIPS_extend.R` (the leaner, Part-1-only predecessor) has been retired — every arm, including this one, now sources `TIPS_extend_v2.R`. Sourced automatically from each `24.1_acat_CTS.cisTarget_dualpull_clean.R` via `tips_extend_r_path()` / `tips_source_extend()` in `00_configuration.R`.

Felix dual-pull uses `CM_candidate` / `CF_candidate` = cisTarget motif target **and** `CM_hi` / `CF_hi` (CTS gene also a DEG in the CM/CF cluster). The now-removed `code_core_1_5vs7` arm (CP1 → Mono vs Neut) had empty `CM_hi`/`CF_hi` on CTS rows, so auto PageRank TFs (e.g. ETS1/FOS/IRF7) produced empty pull subgraphs and `empty final_table` — this is why `TIPS_extend_v2.R` carries HiG-fallback and PageRank re-pick logic below, even though the only arm left in this repo (`code_core_4_9vs11`) doesn't hit it:


| Step             | Arms that already succeed (e.g. `4_9vs11`)               | Sparse arms (`1_5vs7`-like)                                         |
| ---------------- | -------------------------------------------------------- | ------------------------------------------------------------------- |
| CM/CF candidates | Felix `CM_hi` & `CF_hi` (non-zero) → **same as before**  | Often zero → HiG fallback                                           |
| HiG fallback     | **Skipped** (not needed)                                 | Motif targets in `HiG_<CM>` / `HiG_<CF>` when `CM_hi`/`CF_hi` empty |
| PageRank re-pick | **Skipped** if current TFs already have candidates       | Walks PageRank when auto TFs lack dual-pull support                 |


`24.1` **USER INPUT**:


| Setting                    | Default          | Purpose                                                           |
| -------------------------- | ---------------- | ----------------------------------------------------------------- |
| `tips_extend_auto_key_tfs` | `TRUE`           | Walk PageRank for TFs with dual-pull support when auto picks fail |
| `key_TFs_override`         | `NULL`           | Felix-style manual TF list (overrides auto + extend)              |
| `rebuild_heatmap`          | `TRUE` / `FALSE` | Set `TRUE` when changing TFs or after extend picks new key_TFs    |

`motif_target_strategy` (`solo` / `fam` / `merge`) — `run_TIPS_core.R` uses `merge`.

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
7. `24.1` — top 3 PageRank TFs among top 20 in HiGCTS with cisTarget motifs; `motif_target_strategy` solo/fam/merge, which mean: 
    solo (Felix default) --
    Use only the exact column cisTarget_<TF>.motif_target (TF name alone in the header). Strictest: only genes assigned to that TF’s own motif, not a shared family column.
    fam “family” ---
    If the TF also appears in multi-TF family columns (e.g. GATA1;GATA2;…), pick the smallest such family column that contains the TF (fewest TF names in the header = “tightest” family). Ignores the solo column when a family column exists.
    merge (union) --
    Union of all columns where the TF appears (solo + every family column). Creates cisTarget_<TF>.merged_motif_target = 1 if the gene is a target in any of those columns. Most inclusive target set.

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

