# hematopoietic_LARRY/data

Reference resources shared across the `hematopoietic_LARRY` subtrees (`metacell/`,
`single_cell/`, and top-level `8.x`/`2.x` scripts) but not used anywhere outside
`hematopoietic_LARRY/` — genuinely cross-tree resources live in `../../Shared_Data/`
instead. See `../../MISSING_DATA.md` for outstanding data gaps.

## MuTrans attractor-level BioTIP run

| File | Used by |
|---|---|
| `seu_attractor_MuTrans_HVG.rds` | Seurat object (MuTrans attractor grouping), default `seurat_rds`/`LARRY_SEURAT_RDS` across all three `code_core` arms (`metacell/7_data_MuTrans_TIPS_IID`, `metacell/7_data_MuTrans_TIPS_STRING`, `single_cell/code_core_11_10vs17`). |
| `CTS_Lib_data.RData` | BioTIP CTS module list for the attractor-level clustering, `biotip_cts_rdata` in the same three `code_core` arms. |
| `CTS_summary_data.csv` | BioTIP CTS significance summary (MCI_P, IC_g/IC_s local p, localmax) for the same clustering, `biotip_cts_summary_csv`. |

These three are **not** the cell-level leiden_r0_8 BioTIP CTS run used by
`single_cell/9_Held_out_clone_fate_prediction_cursor` (that one lives in that pipeline's
own `data/` folder) — attractor-level and leiden-cell-level are different clusterings.
See the header comment in `9_Held_out_clone_fate_prediction_cursor/00_configuration.R`.

## Held-out clone-fate prediction inputs (from Holly)

| File | Used by |
|---|---|
| `Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds` | Cell-level (leiden_r0_8) re-analyzed Seurat object, `single_cell/9_Held_out_clone_fate_prediction_cursor` (`seurat_rds`) and its own top-level build script `2_invitro_NMtrajectory_creatSeurat_v2.R`. |
| `h5ad_export/obs_metadata.csv` | h5ad-exported cell metadata, same pipeline (`obs_metadata_csv`) and `2.1_normcount_2_scanpy.R`. |
| `cell_cycle_signature_cor_3khvg.rds` | Cell-cycle signature correlation, same pipeline (`hvg_rds`) and top-level `8.2.1`/`8.2.2` figure scripts. |
| `larry_BioTIP/jaccard_STRING_vs_IID_4_9vs11.tsv` | STRING-vs-IID Jaccard comparison, same pipeline (`string_iid_jaccard`) and `metacell/8.3_SEACELL_MuTrans_summary_fig.R`. |
| `larry_figures/td_genes_scores_4_to_{9,11}_seacell.csv` | MuTrans transition-driver gene scores, same pipeline (`mutrans_td_dir`) and the `8.x` comparison/figure scripts across `metacell/` and top-level. |

## Other

| File | Used by |
|---|---|
| `ensembl_peptide_to_gene_mmusculus.rds` | Ensembl protein-ID → mouse gene symbol mapping, `11.1.1_check_vertex_duplication.R` (STRING PPIN construction) in both `metacell/7_data_MuTrans_TIPS_STRING` and `single_cell/code_core_11_10vs17`. |
