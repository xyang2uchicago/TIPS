# Shared_Data

Species- and publication-level reference resources shared across multiple experiment
trees (`cardiac/`, `hematoendothelial/`, `hematopoietic_LARRY/`). Anything
here is read by more than one pipeline — dataset-specific inputs live in each
experiment's own `data/` folder instead. See `../MISSING_DATA.md` for what's still
missing and where each gap stands (get from Holly / published / code-generated / undecided).

## Gene annotation

| File | Contents |
|---|---|
| `coding_genes.rds` | Human (HGNC) protein-coding gene symbols. Code-generated — built by `build_coding_genes.R`, not stored here; run the script to regenerate. |
| `coding_genes_mouse.rds` | Mouse (MGI) protein-coding gene symbols. Code-generated — built by `build_coding_genes_mouse.R`, not stored here; run the script to regenerate. |

## cisTarget motif rankings (`cistarget/`)

| File | Species | Source |
|---|---|---|
| `mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` | mouse | https://resources.aertslab.org/cistarget/databases/ |
| `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` | human | same |

## PPI networks (`PPIN/`)

Expected here, per-species STRING protein-protein interactions and IID annotated PPIs:

| File | Species | Source |
|---|---|---|
| `<taxon>.protein.{aliases,info,links}.v12.0.txt.gz` | per-taxon (e.g. `10090` mouse, `9606` human) | STRING v12.0, https://string-db.org |
| `mouse_annotated_PPIs.txt.gz` | mouse | IID, https://iid.ophid.utoronto.ca/ |
| `human_annotated_PPIs.txt.gz` | human | IID, https://iid.ophid.utoronto.ca/ — **not present yet**, being downloaded. |

## Code-generated

| File | Built by | Depends on |
|---|---|---|
| `iid2025_human_mouse_conserved_global.rds` | `0_IIDweighted_PPIN_build.R` (this folder) | `PPIN/human_annotated_PPIs.txt.gz` — filters IID's human PPI set to mouse-conserved edges.

## Published reference datasets

| File | Source |
|---|---|
| `aaw3381-Weinreb-Table-S3.xlsx` | Weinreb, Rodriguez-Fraticelli, Camargo & Klein, *Science* 2020 — LARRY lineage-tracing paper (GEO GSE140802), Table S3. |
| `stateFate_inVitro_metadata.txt.gz`, `stateFate_inVitro_clone_matrix.mtx.gz` | Public in-vitro LARRY dataset from the same paper. |

## Other files

| File | Used by |
|---|---|
| `CHD_Cilia_Genelist.rds` | Congenital heart disease / cilia gene list — TF-CHD enrichment across both `cardiac/` and `hematoendothelial/` (`IbarraSoria2018_STRING`/`_IID` `code_core_endothelial.b`) pipelines. |
