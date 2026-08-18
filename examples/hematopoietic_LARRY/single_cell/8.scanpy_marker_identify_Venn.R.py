# scp -p -r E:\Git_Holly\BioTIP\R\BioTIP_update_06232025.R   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/seu_dev/GSE175634_iPSC_CM/BioTIP/.
# scp -p -r F:\projects\TIPS\source\GSE140802_lineage_tracking\5_BioTIP_NMtrajectoryCell_leiden0_6.R   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/code/.
# scp -p -r F:\projects\TIPS\source\GSE140802_lineage_tracking\5_BioTIP.sbatch   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/code/.

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

wd_path = '/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/h5ad_export/'

import scanpy as sc
import pandas as pd
import numpy as np
import os
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42  # saved figures maintain text as text
plt.rcParams['ps.fonttype'] = 42 

os.chdir(wd_path)
# -----------------------------
# 1. Load h5ad
# -----------------------------
data_path = '/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/'

adata = sc.read_h5ad(data_path + "h5ad_export/larry_NM_trajectory_normedCounts.h5ad")
adata
# AnnData object with n_obs × n_vars = 92527 × 25289
    # obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Library', 'Cell.barcode', 'Time.point', 'Starting.population', 'Cell.type.annotation', 'Well', 'SPRING.x', 'SPRING.y', 'cell_index0', 'cell_id', 'NM_trajectory', 'CellCycleSig1', 'S.Score', 'G2M.Score', 'Phase', 'leiden_r0_2', 'leiden_r0_4', 'leiden_r0_6', 'leiden_r0_8', 'leiden_r1', 'leiden_r1_2', 'leiden_r0_5', 'Cell.type.clean', 'leiden_r0_6_anno'
    # uns: 'data_provenance', 'Seurat_scaled_variable_genes', 'Cell.type.clean_colors'
    # obsm: 'X_pca', 'X_spring', 'X_umap', 'Seurat_scaled_data'

ct_col = "Cell.type.clean"

celltypes_use = list(adata.obs[ct_col].cat.categories)
celltypes_use
# ['Baso', 'Ccr7_DC', 'Eos', 'Erythroid', 'Lymphoid', 'Mast', 'Meg', 'Mono', 'Neut', 'Unddiff', 'pDC']

published_markers = {
    # from Table S3, DGE of progenitors in vitro
    "Baso": ["Ligp1", "Ly6a", "Cpa3",], #Fig 2J
    "Eos": ["Cpa3", "Ikzf2", "Cyp11a1", "Akr1c13", "Alox5"],
    "Meg": ["Podxl", "Itga2b", "Plscr2", "Pbx1", "Nckap1"],
    "Mast": ["Ikzf2", "Arnt", "1810022K09Rik"],
    "Neut": ["GP3", "Ly6a",#Fig 2J
                "S100a8", "Cd74", "Csf1r", 'Elane', 'H2-Aa', 'Klf4'],  # Fig 4A
    "Mono": ["Ly6a", "Tie1", "S100a4"], #Fig 2J
    "Lymphoid": ["Ms4a3", "Ctsg", "Gjb3", "Foxp1", "Mdm2"],  # Fig 2J

    # Table S3 has only one "Dendritic cell" category, not Ccr7_DC vs pDC
    "DC": ["Flt3", "H2afy", "Ighm", "Dnlt"],  # Fig 2J

    # from Table S3, DGE of progenitors in vivo
    "Erythroid": ["Gata2", 'Itga2b','Car2'],  # Fig 2J

    # no "Undifferentiated" category in Table S3;
    # use in vivo MPP/GMP progenitor-associated genes as fallback
    "Unddiff": ["Ly6a", "Meis1", "Flt3", "Grb10", "Lgals3bp",   # Fig S11
                'Procr',  # page 3
                'Rbp1', 'Ly6a', 'Esam', 'Hlf'], # Fig 2J

    "cellCycle": ["Ube2c", "Hmgb2", "Hmgn2", "Tuba1b", "Ccnb1", "Tubb5", "Top2a", "Tubb4b"],

}


# -----------------------------
# 2. Make sure cell type is categorical
# -----------------------------
adata.obs[ct_col] = adata.obs[ct_col].astype("category")

present_celltypes = list(adata.obs[ct_col].cat.categories)
missing = sorted(set(celltypes_use) - set(present_celltypes))

print("Missing cell types:", missing) # []

celltypes_use = [x for x in celltypes_use if x in present_celltypes]

# subset to requested cell types only
ad = adata[adata.obs[ct_col].isin(celltypes_use)].copy()
ad.obs[ct_col] = ad.obs[ct_col].cat.remove_unused_categories()
ad.obs[ct_col] = ad.obs[ct_col].cat.reorder_categories(celltypes_use, ordered=True)

# -----------------------------
# 3. Make sure X is log-normalized expression
# -----------------------------
# Since your workflow used:
# counts = total-count-normalized counts
# data   = log1p(total-count-normalized counts)
#
# For marker detection, use log expression, NOT scaled data.
# If your h5ad has log expression in X, do nothing.
#
# If instead you saved log expression in a layer, uncomment one of these:
ad.layers["norm_counts"] = ad.X.copy()
ad.X = np.log1p(ad.X)
# ad.X = ad.layers["data"].copy()

#####################################################
## plot the known markes based on the Seurat umap and the downloaded SPRING embedding
#####################################################

# handle missing genes
all_genes = sorted(set(sum(published_markers.values(), [])))
missing_genes = [g for g in all_genes if g not in adata.var_names]
print("Missing genes:", missing_genes) # []

published_markers_present = {
    ct: [g for g in genes if g in adata.var_names]
    for ct, genes in published_markers.items()
    if ct in present_celltypes
}

# optional: print final marker panel used
for ct, genes in published_markers_present.items():
    print(ct, genes)
# Baso ['Gata2', 'Akr1c13', 'Cpa3', 'Alox5', 'Slc18a2']
# Eos ['Cpa3', 'Ikzf2', 'Cyp11a1', 'Akr1c13', 'Alox5']
# Meg ['Podxl', 'Itga2b', 'Plscr2', 'Pbx1', 'Nckap1']
# Mast ['Ikzf2', 'Arnt', '1810022K09Rik']
# Neut ['Muc13', 'Srgn', 'Ccl9', 'Elane', 'Igfbp4']
# Mono ['Rbms1', 'Sirpa', 'Set', 'Tuba1b', 'Tk1']
# Lymphoid ['Ighm', 'Il12a', 'Ctr9', 'Gsn', 'Satb1']
# Ccr7_DC ['Ly86', 'Ighm', 'Irf8', 'Mpeg1', 'Olfm1']
# pDC ['Ly86', 'Ighm', 'Irf8', 'Mpeg1', 'Olfm1']
# Erythroid ['Podxl', 'F2r', 'Itga2b', 'Pf4', 'Car2']
# Unddiff ['Rbp1', 'Esam', 'Ly6a', 'Grb10', 'Lgals3bp']

sc.pl.dotplot(
    adata,
    var_names=published_markers_present,
    groupby=ct_col,
    categories_order=celltypes_use,
    standard_scale="var",
    dendrogram=False,
    # save="_published_TableS3_markers_Cell.type.clean.pdf"
    save="_published_markers_Cell.type.clean.pdf"
)



#####################################################
## data-driven marker identification
#####################################################
# -----------------------------
# 4. Differential marker test
# -----------------------------
sc.tl.rank_genes_groups(
    ad,
    groupby=ct_col,
    groups=celltypes_use,
    reference="rest",
    method="t-test_overestim_var",  # wilcox is not good
    use_raw=False,
    pts=True,
    corr_method="benjamini-hochberg"
)

sc.pl.rank_genes_groups(ad,
                        save="_rankgene_tov_Cell.type.clean.pdf")  # top 20 by default


# robust extraction across Scanpy versions
dfs = []
for ct in celltypes_use:
    tmp = sc.get.rank_genes_groups_df(ad, group=ct)
    tmp["group"] = ct
    dfs.append(tmp)

markers = pd.concat(dfs, axis=0)
markers.shape()
markers.to_csv("markers_per_" + ct_col +"tov.csv", index=False)  #!!!!!!!!!!

## save for future follow ups
ad.write_h5ad("larry_NM_trajectory_log1pNormedCounts_rankgene_tov.h5ad") #!!!!!!!!!!!!!


# -----------------------------
# 5. Filter and get top 5 per cell type
# -----------------------------
markers = markers.replace([np.inf, -np.inf], np.nan)

# keep positive markers only
if "logfoldchanges" in markers.columns:
    markers_pos = markers[markers["logfoldchanges"] > 0].copy()
else:
    markers_pos = markers.copy()


# require adjusted p-value when available
if "pvals_adj" in markers_pos.columns:
    markers_sig = markers_pos[markers_pos["pvals_adj"] < 0.05].copy()
else:
    markers_sig = markers_pos.copy()


# optional: require expressed in at least 10% of cells in the group
if "pct_nz_group" in markers_sig.columns:
    markers_sig = markers_sig[markers_sig["pct_nz_group"] >= 0.10].copy()

print("Number of markers passing filters:", markers_sig.shape[0]) # 1100
top5_list = []

for ct in celltypes_use:
    d = markers_sig[markers_sig["group"] == ct].copy()
    # fallback if too few pass filters
    if d.shape[0] < 5:
        d = markers_pos[markers_pos["group"] == ct].copy()
    # sort the markers by scores in descending order and get the top 5
    d = d.sort_values("scores", ascending=False).head(5)
    top5_list.append(d)

top5_markers = pd.concat(top5_list, axis=0)

# -----------------------------
# 6. Save results
# -----------------------------
cols_keep = [
    c for c in [
        "group",
        "names",
        "scores",
        "logfoldchanges",
        "pvals",
        "pvals_adj",
        "pct_nz_group",
        "pct_nz_reference"
    ]
    if c in top5_markers.columns
]

top5_markers = top5_markers[cols_keep]
top5_markers.to_csv("top5_markers_per_" + ct_col +".csv", index=False)

print(top5_markers)

# convenient dictionary
top5_dict = (
    top5_markers
    .groupby("group")["names"]
    .apply(list)
    .to_dict()
)

print(top5_dict)

# Dotplot of top 5 markers per cell type
sc.pl.dotplot(
    ad,
    var_names=top5_dict,
    groupby=ct_col,
    standard_scale="var",
    dendrogram=False,
    save="_top5_markers_Cell.type.clean.pdf"
)


###############################################
### in local R , compare the identified top 5 markers wit hpublished markers
###############################################
# scp -p -r   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/* F:\projects\TIPS\results\GSE140802_lineage_tracking\inVitro_NMtrajectory\.

setwd('F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/')

top5 = read.csv('top5_markers_per_Cell.type.clean.csv')
dim(top5) # 55  8
head(top5)
#     group   names    scores logfoldchanges        pvals    pvals_adj
# 1    Baso    Ccl6 168.49463       6.010245 0.000000e+00 0.000000e+00
# 2    Baso   Gata2 155.59550       5.120154 0.000000e+00 0.000000e+00
# 3    Baso Csf2rb2 152.87750       5.454580 0.000000e+00 0.000000e+00
# 4    Baso    Cpa3 147.98814       6.863646 0.000000e+00 0.000000e+00
# 5    Baso  Csf2rb 124.41495       3.916251 0.000000e+00 0.000000e+00
# 6 Ccr7_DC   H2-Aa  37.75036       8.377263 1.217643e-89 3.079296e-85
#   pct_nz_group pct_nz_reference
# 1    0.9809636       0.26550888
# 2    0.9715629       0.17231965
# 3    0.9760282       0.23527279
# 4    0.9819036       0.08787611
# 5    0.9706228       0.44916848
# 6    0.9747899       0.07064323
unique(top5$group)
#  [1] "Baso"      "Ccr7_DC"   "Eos"       "Erythroid" "Lymphoid"  "Mast"     
#  [7] "Meg"       "Mono"      "Neut"      "Unddiff"   "pDC"    
top5$group_clean = top5$group
top5$group_clean[top5$group_clean == "Ccr7_DC"] = "DC"
top5$group_clean[top5$group_clean == "pDC"] = "DC"
unique(top5$group_clean)
#  [1] "Baso"      "DC"        "Eos"       "Erythroid" "Lymphoid"  "Mast"     
#  [7] "Meg"       "Mono"      "Neut"      "Unddiff"      


merge(top5, pub_markers, by.x = "group_clean", by.y = "lineage_clean", all.x = TRUE)
#    group_clean   names    scores logfoldchanges        pvals    pvals_adj
# 1         Baso    Ccl6 168.49463       6.010245 0.000000e+00 0.000000e+00
# 2         Baso   Gata2 155.59550       5.120154 0.000000e+00 0.000000e+00
# 3         Baso Csf2rb2 152.87750       5.454580 0.000000e+00 0.000000e+00
# 4         Baso    Cpa3 147.98814       6.863646 0.000000e+00 0.000000e+00
# 5         Baso  Csf2rb 124.41495       3.916251 0.000000e+00 0.000000e+00

data_path = 'F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/'
pub_markers = readxl::read_excel(paste0(data_path, 'aaw3381-Weinreb-Table-S3.xlsx'), sheet = 1)
dim(pub_markers) # [1] 447   4
head(pub_markers)
# # A tibble: 6 × 4
#   Lineage       `Gene symbol` `"-log10(adj. p-val)"` "log2(fold-enrichment)…¹
#   <chr>         <chr>                          <dbl>                    <dbl>
# 1 Megakaryocyte Podxl                           9                        0.32
# 2 Megakaryocyte Itga2b                          9                        0.33
# 3 Megakaryocyte Plscr2                          9                        0.08
# 4 Megakaryocyte Pbx1                            9                        0.38
# 5 Megakaryocyte Nckap1                          9                        0.08
# 6 Megakaryocyte Mpl                             8.95                     0.14
# # ℹ abbreviated name: ¹​`"log2(fold-enrichment)"`

unique(pub_markers$Lineage)
# [1] "Megakaryocyte"  "Mast cell"      "Basophil"       "Eosinophil"    
# [5] "Neutrophil"     "Monocyte"       "Dendritic cell" "Lymphoid"      
pub_markers$lineage_clean = pub_markers$Lineage
pub_markers$lineage_clean[pub_markers$lineage_clean == "Megakaryocyte"] = "Meg"
pub_markers$lineage_clean[pub_markers$lineage_clean == "Mast cell"] = "Mast"
pub_markers$lineage_clean[pub_markers$lineage_clean == "Basophil"] = "Baso"
pub_markers$lineage_clean[pub_markers$lineage_clean == "Eosinophil"] = "Eos"
pub_markers$lineage_clean[pub_markers$lineage_clean == "Neutrophil"] = "Neut"
pub_markers$lineage_clean[pub_markers$lineage_clean == "Monocyte"] = "Mono"
pub_markers$lineage_clean[pub_markers$lineage_clean == "Dendritic cell"] = "DC"

unique(pub_markers$lineage_clean)
# [1] "Meg"      "Mast"     "Baso"     "Eos"      "Neut"     "Mono"    
# [7] "DC"       "Lymphoid"

library(dplyr)
setdiff(top5$group_clean, pub_markers$lineage_clean) %>% unique()  # [1] "Erythroid" "Unddiff"  


pub_markers2  = readxl::read_excel(paste0(data_path, 'aaw3381-Weinreb-Table-S3.xlsx'), sheet = 2)
dim(pub_markers2) # [1]  190   4
head(pub_markers2)
unique(pub_markers2$Lineage)
# [1] "Neutrophil"     "B cell"         "Monocyte"       "Erythroid"     
# [5] "T cell"         "MPP/GMP"        "NK cell"        "Dendritic cell"
# [9] "Basophil"
pub_markers2$lineage_clean = pub_markers2$Lineage
pub_markers2$lineage_clean[pub_markers2$lineage_clean == "Neutrophil"] = "Neut"
pub_markers2$lineage_clean[pub_markers2$lineage_clean == "B cell"] = "Lymphoid"
pub_markers2$lineage_clean[pub_markers2$lineage_clean == "Monocyte"] = "Mono"
pub_markers2$lineage_clean[pub_markers2$lineage_clean == "Erythroid"] = "Erythroid"
pub_markers2$lineage_clean[pub_markers2$lineage_clean == "T cell"] = "Lymphoid"
pub_markers2$lineage_clean[pub_markers2$lineage_clean == "MPP/GMP"] = "Unddiff"
pub_markers2$lineage_clean[pub_markers2$lineage_clean == "NK cell"] = "Lymphoid"
pub_markers2$lineage_clean[pub_markers2$lineage_clean == "Dendritic cell"] = "DC"
pub_markers2$lineage_clean[pub_markers2$lineage_clean == "Basophil"] = "Baso"
unique(pub_markers2$lineage_clean)
# [1] "Neut"      "Lymphoid"  "Mono"      "Erythroid" "Unddiff"   "DC"       
#[7] "Baso" 


### venn diagran of pub_markers2 vs pub_markers per lineage_clean
library(gplots)
x = intersect(pub_markers2$lineage_clean, pub_markers$lineage_clean)
y = intersect(top5$group_clean, pub_markers$lineage_clean)
z = intersect(top5$group_clean, pub_markers2$lineage_clean)
x = unique(c(x, y , z))
par(mfrow=c(2, round(length(x) / 2)))
for (i in x) {
    venn(list(Pub2 = pub_markers2$`Gene symbol`[which(pub_markers2$lineage_clean == i)], 
              Pub1 = pub_markers$`Gene symbol`[which(pub_markers$lineage_clean == i)],
              Top5 = top5$names[which(top5$group_clean == i)]))

    mtext(i, side = 3, line = -1, font = 2)
}


dev.off()
