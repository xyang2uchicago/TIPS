# cardiac/data

Reference resources shared across the `cardiac/` dataset trees (`GSE175634`, `GSE87038`,
`IbarraSoria2018`) but not used anywhere outside `cardiac/` — genuinely cross-tree
resources (used by `hematoendothelial/` or `hematopoietic_LARRY/` too) live in
`../../Shared_Data/` instead. See `../../MISSING_DATA.md` for outstanding data gaps.

## Fetal/in-vitro heart annotation (Ameen et al. 2022, *Cell*)

| File | Used by |
|---|---|
| `Ameen2022cell-supplement-10.xlsx` | In-vitro (iPSC-CM) marker gene annotation, `24.0_acat_load_input_clean*.R` (sheet 4). |
| `Ameen2022cell_Table1_sheet3_markergene_Laksshman2026update.xlsx` | Fetal heart marker gene annotation, `24.0_acat_load_input_clean*.R` (sheet 1). |
| `Ameen2022cell_Table1_Sheet1_cellmetadata_Laksshman2026update.xlsx` | Not currently referenced by any script. |
| `Ameen2022cell_media-1.xlsx` | Not currently referenced live — only appears in commented-out lines in `24.0_acat_load_input_clean*.R`. |

## ISL1 gene sets

| File | Used by |
|---|---|
| `GSE195476_ISL1/ISL1_set.rds` | ISL1 ChIP/target gene sets, cardiac progenitor annotation — `24.0_acat_load_input_clean*.R`, `11.8_compared_ISL1KO_DEGs*.R`, `24_acat_CTS.CP_pioneer_heatmap*.R`, `11.10_CTS_table_report.R`. |
| `GSE80383_Isl1/ISL1_Mm2Hg.GS_uniqueSymbol_FC1.2_padj0.05.RData` | Same scripts — mouse-to-human-mapped ISL1 gene sets. |
| `ISL1_CP_bound_gene.csv` | Not currently referenced by any script. |

## Other files

| File | Used by |
|---|---|
| `readme_filename_map_xy.txt` | Filename → cell-type identity mapping for GSE175634 pseudobulk/peak files, `24.0_acat_load_input_clean*.R`. |
