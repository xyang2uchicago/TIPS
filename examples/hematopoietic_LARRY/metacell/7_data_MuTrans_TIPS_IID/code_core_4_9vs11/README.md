# TIPS IID — Larry MuTrans (`code_core_4_9vs11`)

Felix `GSE87038_IID/code_core_13/` analogue. **Target:** in-house IID PPIN (not STRING).

> **Current state:** runners still source `11.1_STRINGweighted_...` and STRING paths in `00_configuration.R` until `11.1_IIDweighted_...` is wired. Use sibling `7_data_MuTrans_TIPS_STRING/` for production STRING runs.

## This arm

| Parameter | Value |
|-----------|--------|
| TAG | `4_9vs11` |
| CP / CM / CF | `4` / `9` (Baso) / `11` (Meg) |
| `group_col` | `attractor` |
| `motif_target_strategy` | `merge` |
| Results | `../results_core_4_9vs11/` |

## Run

```r
source("F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_IID/code_core_4_9vs11/run_TIPS_core.R")
```

## Step order (intended)

`11.1_IID` → `11.2.0` → `11.3` → `12.0` → `24.1` → `24.3` → `25_venn` (no `11.1.1` for IID)

## Expression / grouping

Same as STRING arms: `load_larry_sce()` from `seu_attractor_MuTrans_HVG.rds`, `group_col` default `attractor`.

## Larry vs Felix / Holly

Same dual-pull + cisTarget logic as STRING README (`code_core_4_9vs11/README.md` in STRING folder). IID replaces STRING 11.1–11.3 graph build with IID-weighted graphs (`*_IID_graph_*`).

## Environment variables

Same as STRING; `TIPS_WD` should point to `7_data_MuTrans_TIPS_IID/` when IID config is finalized.

## Key outputs (IID naming, when wired)

- `GSE140802_IID_graph_perState_notsimplified.rds`
- `PPI_weight/GSE140802_IID_graph_perState_simplified_combinedweighted.rds`
- `cisTarget_predicted_4/...dualpull_final_table.tsv`
