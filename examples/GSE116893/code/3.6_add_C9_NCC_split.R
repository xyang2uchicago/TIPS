########## BEGINNING OF USER INPUT ##########
wd <- "C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/"
########## END OF USER INPUT ##########

library(Seurat)
library(SingleCellExperiment)

data_dir <- paste0(wd, "data/")

# ---- load Seurat object and compute module scores ----
# malignant_obs_scoreAdded.rds is not available locally; recompute scores here.
# Gene sets from PI's Umap_Exploration_xy_patient_level_NCCscore.R
seurat_obj <- readRDS(paste0(data_dir, "malignant_leiden.rds"))

migratoryNCC_genes  <- list(c("SOX9", "FOXD3", "SNAI2"))
notch_core_genes    <- list(c("HES1", "HES5", "HEY1", "HEY2", "HEYL", "NRARP", "DTX1", "DTX4"))
notch_feedfwd_genes <- list(c("NOTCH1", "NOTCH2", "NOTCH3", "JAG1", "MAML2", "HES1"))

seurat_obj <- AddModuleScore(seurat_obj, features = migratoryNCC_genes,  name = "NCCmigratory")
seurat_obj <- AddModuleScore(seurat_obj, features = notch_core_genes,    name = "NOTCH_core")
seurat_obj <- AddModuleScore(seurat_obj, features = notch_feedfwd_genes, name = "NOTCH_feedforward")
# Seurat appends "1" to the name → columns: NCCmigratory1, NOTCH_core1, NOTCH_feedforward1

range(seurat_obj@meta.data[["NCCmigratory1"]])

# ---- load existing SCE and add scores to colData ----
sce <- readRDS(paste0(data_dir, "sce.rds"))
dim(sce)
table(colData(sce)[["leiden_0.6"]])

stopifnot(all(colnames(sce) %in% rownames(seurat_obj@meta.data)))

score_map <- c(
    NCCmigratory_1      = "NCCmigratory1",
    NOTCH_core_1        = "NOTCH_core1",
    NOTCH_feedforward_1 = "NOTCH_feedforward1"
)
for (sce_col in names(score_map)) {
    seurat_col <- score_map[[sce_col]]
    colData(sce)[[sce_col]] <- seurat_obj@meta.data[colnames(sce), seurat_col]
    # NCCmigratory_1 range: -0.1899229 1.462688
}

cat("NCCmigratory_1 range:", range(colData(sce)[["NCCmigratory_1"]]), "\n")

# ---- create C9 split cluster label ----
# Cells in cluster 9:
#   NCCmigratory_1 >  0 → "9_hi"  (NCC-high = early-migratory progenitor / CP analog)
#   NCCmigratory_1 <= 0 → "9_lo"  (NCC-low  = committed / CF analog)
# All other clusters: keep leiden_0.6 value unchanged.
cluster <- as.character(colData(sce)[["leiden_0.6"]])
ncc     <- colData(sce)[["NCCmigratory_1"]]
is_c9   <- cluster == "9"

leiden_C9split                                          <- cluster
leiden_C9split[is_c9 & !is.na(ncc) & ncc >  0]        <- "9_hi"
leiden_C9split[is_c9 & !is.na(ncc) & ncc <= 0]        <- "9_lo"
leiden_C9split[is_c9 &  is.na(ncc)]                    <- "9_lo"  # NA → low

colData(sce)[["leiden_C9split"]] <- leiden_C9split

cat("\nCluster 9 split by NCCmigratory_1:\n")
print(table(leiden_C9split[is_c9]))
# Cluster 9 split by NCCmigratory_1:

#  9_hi  9_lo 
#  1001 10867 


cat("\nFull cluster table (leiden_C9split):\n")
print(table(leiden_C9split))
# leiden_C9split
#     1    10    11    12    13    14    15    16    17     2     3     4     5 
# 28394 10039  2911  1862  1834   642   348   189    73 28380 27401 25381 24098 
#     6     7     8  9_hi  9_lo 
# 19850 19743 18215  1001 10867

# ---- save ----
saveRDS(sce, file = paste0(data_dir, "sce_C9split.rds"))
cat("\nSaved:", paste0(data_dir, "sce_C9split.rds"), "\n")
