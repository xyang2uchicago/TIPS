### sce -> h5ad on midway3 ##############################################
# ssh xyang2@midway3.rcc.uchicago.edu
# module load python/anaconda-2022.05
# source activate /project/xyang2/software-packages/env/velocity_2025Feb_xy/

# rcchelp balance  
# rcchelp usage
# squeue -u xyang2
# squeue -p bigmem --state=PD | wc -l
# squeue -p caslake --state=PD | wc -l

# sinteractive -p bigmem  --mem=500G  --account=pi-xyang2 --time=1:00:00  -c 3    # requires n_cores  
# sinteractive -p caslake  --cpus-per-task=1 --mem=180G --time=6:00:00  --account=pi-xyang2 -c 3 
# sinteractive -p gpu --account=pi-xyang2 --gres=gpu:1 --mem=180GB  --time=6:00:00 -c 3  
# cd /project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/
# R
## To quit your interactive job:
# exit or Ctrl-D


setwd('/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/')


# library(SingleCellExperiment)
# library(slingshot)
library(dplyr) 
 

# larry_nm <- readRDS( "Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds")
# larry_nm
# # An object of class Seurat 
# # 25289 features across 92527 samples within 1 assay 
# # Active assay: RNA (25289 features, 3000 variable features)
# #  3 layers present: counts, data, scale.data
# #  3 dimensional reductions calculated: spring, pca, umap
# colnames(larry_nm@meta.data)
# #  [1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"
# #  [4] "Library"              "Cell.barcode"         "Time.point"
# #  [7] "Starting.population"  "Cell.type.annotation" "Well"
# # [10] "SPRING.x"             "SPRING.y"             "cell_index0"
# # [13] "cell_id"              "NM_trajectory"        "Fig6B_label"
# # [16] "CellCycleSig1"        "Cell.type.clean"      "leiden_r0_2"
# # [19] "leiden_r0_4"          "leiden_r0_6"          "leiden_r0_8"
# # [22] "leiden_r1"            "leiden_r1_2"          "leiden_r0_5"
# # [25] "leiden_r0_6_anno"


# ###### 1) extract raw counts matrix ##############################################
# ## vewcasue the only count matrix ever released is stateFate_inVitro_normed_counts.mtx.gz,
# # explicitly documented as reporting UMI transcript counts per gene per cell after total-counts normalization (L1 normalization on cells). 
# range(larry_nm$nCount_RNA)
# sd(larry_nm$nCount_RNA)  # Every one of your cells has the exact same total
# ## it is recomputed nCount_RNA after total-count normalization, not raw depth
# # [1] 2878.749 2878.749
# hist(larry_nm$nCount_RNA, breaks = 50)

library(zellkonverter); library(SingleCellExperiment); library(Seurat)
package.version('zellkonverter') #"1.12.5"
package.version('SingleCellExperiment') # "1.24.1"
package.version('Seurat') #"5.3.0"

seu = readRDS('Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds')
seu
# An object of class Seurat
# 25289 features across 92527 samples within 1 assay
# Active assay: RNA (25289 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  3 dimensional reductions calculated: spring, pca, umap

range(seu[['RNA']]$data) #[1] 0.00000 6.73168
range(seu[['RNA']]$scale.data) #[1] -2.159904 10.000000
dim(seu[['RNA']]$data) # 25289 92527
dim(seu[['RNA']]$scale.data) #  3000 92527


library(Matrix)
dir.create("h5ad_export", showWarnings = FALSE)

range(seu[['RNA']]$data) #[1] 0.00000 6.73168

normed <- seu[["RNA"]]$counts        # sparse, total-count normalized (not raw)
writeMM(normed, file = "h5ad_export/normed_counts.mtx")
write.csv(rownames(normed), "h5ad_export/genes.csv", row.names = FALSE)
write.csv(colnames(normed), "h5ad_export/barcodes.csv", row.names = FALSE)
write.csv(seu@meta.data, "h5ad_export/obs_metadata.csv", row.names = TRUE)
 

# mat <- as(seu[["RNA"]]$data , "dgCMatrix")
# writeMM(mat, file = "h5ad_export/normed_log1p.mtx")
mat <- as(seu[["RNA"]]$scale.data , "dgCMatrix")
writeMM(mat, file = "h5ad_export/scale.data.mtx")
write.csv(rownames(mat), "h5ad_export/scale.data.genes.csv", row.names = FALSE)


for (red in Seurat::Reductions(seu)) {
  write.csv(Seurat::Embeddings(seu, red),
            paste0("h5ad_export/", red, "_embedding.csv"), row.names = TRUE)
}

writeLines(paste(
  "'normed_counts'/X is total-count (L1) normalized UMI data from GSE140802's",
  "stateFate_inVitro_normed_counts.mtx.gz. NOT raw counts - none were ever",
  "publicly released for this dataset.",
  "spring embedding is downloaded, pca and umap is computed byHolly using Seurat 5.",
  "normed_log1p.mtx is the log1p-transformed normed counts matrix,",
  "scale.data.mtx is the scaled data matrix."
), "h5ad_export/data_provenance.txt")


celltype_pal_df <- data.frame(
  Cell.type.clean = names(celltype_pal),
  color = unname(celltype_pal)
)

write.csv(celltype_pal_df, "celltype_pal.csv", row.names = FALSE)


################# check the h5ad file on midway3 ##############################################

import anndata as ad, pandas as pd, scipy.io as sio, glob, os


# pthyeon set wd to /project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/
os.chdir('/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/')


norm_counts = sio.mmread("h5ad_export/normed_counts.mtx").T.tocsr()   # -> cells x genes

genes = pd.read_csv("h5ad_export/genes.csv")["x"].values
barcodes = pd.read_csv("h5ad_export/barcodes.csv")["x"].values
obs = pd.read_csv("h5ad_export/obs_metadata.csv", index_col=0)

adata = ad.AnnData(X=norm_counts, obs=obs, var=pd.DataFrame(index=genes))  #!!!!!!!!!!
adata.obs_names = barcodes

## add Seurat-exported embeddings
for f in glob.glob("h5ad_export/*_embedding.csv"):
    name = os.path.basename(f).replace("_embedding.csv", "")
    adata.obsm[f"X_{name}"] = pd.read_csv(f, index_col=0).values

adata.uns['data_provenance'] = open("h5ad_export/data_provenance.txt").read()

# read scaled matrix exported from R
scaled = sio.mmread("h5ad_export/scale.data.mtx").T.tocsr()  
genes = pd.read_csv("h5ad_export/scale.data.genes.csv")["x"].values
# make sure order matches adata
scaled_df = pd.DataFrame.sparse.from_spmatrix(
    scaled,
    index=adata.obs_names,
    columns=genes
)

# AnnData store it as obsm, not layers to avoid confusion of gene dimention (the current X is all genes alough L1 normalized)
adata.obsm["Seurat_scaled_data"] = scaled_df.sparse.to_coo().tocsr()
adata.uns["Seurat_scaled_variable_genes"] = genes


## save color palette for cell types
pal_df = pd.read_csv("celltype_pal.csv")
celltype_pal = dict(zip(pal_df["Cell.type.clean"], pal_df["color"]))

adata.obs["Cell.type.clean"] = adata.obs["Cell.type.clean"].astype("category")
cats = adata.obs["Cell.type.clean"].cat.categories

# check mismatch
print(set(cats) - set(celltype_pal.keys()))  # set()

adata.uns["Cell.type.clean_colors"] = [celltype_pal[x] for x in cats]


adata.write_h5ad("h5ad_export/larry_NM_trajectory_normedCounts.h5ad")


### verify

import anndata as ad, pandas as pd, scipy.io as sio, glob, os
import scanpy as sc

# pthyeon set wd to /project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/
data_path = '/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/'

adata = sc.read_h5ad(data_path + "h5ad_export/larry_NM_trajectory_normedCounts.h5ad")
adata
# AnnData object with n_obs × n_vars = 92527 × 25289
    # obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Library', 'Cell.barcode', 'Time.point', 'Starting.population', 'Cell.type.annotation', 'Well', 'SPRING.x', 'SPRING.y', 'cell_index0', 'cell_id', 'NM_trajectory', 'CellCycleSig1', 'S.Score', 'G2M.Score', 'Phase', 'leiden_r0_2', 'leiden_r0_4', 'leiden_r0_6', 'leiden_r0_8', 'leiden_r1', 'leiden_r1_2', 'leiden_r0_5', 'Cell.type.clean', 'leiden_r0_6_anno'
    # uns: 'data_provenance', 'Seurat_scaled_variable_genes', 'Cell.type.clean_colors'
    # obsm: 'X_pca', 'X_spring', 'X_umap', 'Seurat_scaled_data'

### this is the total-count normalized (L1) counts matrix, I can't get the real raw counts matrix for this dataset,
### However, all genes are present in the matrix. Seurat analysis based on it makes sense. 
adata.X.min() # 0.00000
adata.X.max() # 837.5551043206162


adata.obs.columns
# Index(['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Library', 'Cell.barcode',
#        'Time.point', 'Starting.population', 'Cell.type.annotation', 'Well',
#        'SPRING.x', 'SPRING.y', 'cell_index0', 'cell_id', 'NM_trajectory',
#        'CellCycleSig1', 'S.Score', 'G2M.Score', 'Phase', 'leiden_r0_2',
#        'leiden_r0_4', 'leiden_r0_6', 'leiden_r0_8', 'leiden_r1', 'leiden_r1_2',
#        'leiden_r0_5', 'Cell.type.clean', 'leiden_r0_6_anno'],
#       dtype='object')
adata.var.index
# Index(['0610006L08Rik', '0610007P14Rik', '0610009B22Rik', '0610009E02Rik',
#        '0610009L18Rik', '0610009O20Rik', '0610010F05Rik', '0610010K14Rik',
#        '0610011F06Rik', '0610012D04Rik',
#        ...
#        'mt-Co2', 'mt-Co3', 'mt-Cytb', 'mt-Nd1', 'mt-Nd2', 'mt-Nd3', 'mt-Nd4',
#        'mt-Nd4l', 'mt-Nd5', 'mt-Nd6'],
#       dtype='object', length=25289)







# ### then
# rm  normed_counts.mtx
# rm  genes.csv
# rm  barcodes.csv
# rm umap_embedding.csv
# rm pca_embedding.csv
# rm spring_embedding.csv
