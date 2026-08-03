
setwd("F:/projects/TIPS/results/GSE140802_lineage_tracking")

library(Seurat)
packageVersion("Seurat") # ‘5.4.0.9022’
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(patchwork)
library(curl)

library(dplyr)
library(ggrastr)
library(ggrepel)


## Use Seurat v5 assay style
options(Seurat.object.assay.version = "v5")
options(timeout = 3600)

## 1) dwonload the expriment 1 in vitro data and save locallya
# urls <- c(
#   counts = "https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVitro_normed_counts.mtx.gz",
#   genes  = "https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVitro_gene_names.txt.gz",
#   meta   = "https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVitro_metadata.txt.gz",
#   nm_traj = "https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVitro_neutrophil_monocyte_trajectory.txt.gz",
#   clone  = "https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/stateFate_inVitro_clone_matrix.mtx.gz"
# )

data_path <- "f:/projects/TIPS/data/GSE140802_lineage_tracking/Github/"
indir = 'inVitro/'
outdir = 'inVitro_NMtrajectory/'

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## 4. load the full Seurat object
## refer to 0_Experiment1_invitro_creatSeurat.R
## -----------------------------
# the authors already did QC: Kept abundant cell barcodes based on total transcript histograms; 
# removed cells with >20% mitochondrial transcripts; 
# used SCRUBLET with manual thresholding; doublet filtering removed ~10% of cells.
# Cells within each experiment were then normalized to have the same total number of transcripts

### Use these author-matching choices:
npcs = 30  # was 50, reduced for speeding up UMAP only (not PCA) according to elbow plot
visualization_original = "SPRING"  # not "UMAP"
SPRING_filter_genes = c(85, 3, 3)
# cell_cycle_cor_threshold = 0.1
# author: remove the paper's cell-cycle-correlated signature genes
cell_cycle_signature <- c(
  "Ube2c", "Hmgb2", "Hmgn2", "Tuba1b",
  "Ccnb1", "Tubb5", "Top2a", "Tubb4b"
)
kNN_for_SPRING = 4
ForceAtlas2_steps = 500

# # author: NeuMo trajectory smoothing parameters
# ## we don't need because the stateFate_inVitro_neutrophil_pseudotime.txt.gz is already provided 
# NeuMo_trajectory_smoothing_N1 = 250
# NeuMo_trajectory_threshold1 = 0.40
# NeuMo_trajectory_smoothing_N2 = 50
# NeuMo_trajectory_threshold2 = 0.40
# NeuMo_coefficients = c(Neu = 0.1, Mo = 0.1, MPP = 1, Meg = -2, Other = -1)

# author: FateID parameters
FateID_cells = 20000
FateID_HVGs = 1000
heldout_clone_fraction = 0.5
heldout_fate_KNN = 400

larry <- readRDS(file.path(indir, "reanalyzed_Seurat5_noCellCycle_PCA_UMAP.rds"))
DefaultAssay(larry) # "RNA"

############################################################
## 2.0. understand the data distribution
############################################################
meta = larry@meta.data
table(meta$NM_trajectory )
# FALSE  TRUE 
# 38360 92527 
table(meta$NM_trajectory, meta$Cell.type.annotation)
       
  #        Baso Ccr7_DC   Eos Erythroid Lymphoid  Mast   Meg Monocyte
  # FALSE  5900     122   200       196      668  1033   628    12269
  # TRUE   4255     119   355       490      189  1551  1468     6913
       
  #       Neutrophil   pDC Undifferentiated
  # FALSE          4    64            17276
  # TRUE       22233    37            54917
## Barplot: Fraction of NM trajectory cells per cell type
table(meta$NM_trajectory, meta$Time.point)
#            2     4     6
  # FALSE   774 19447 18139
  # TRUE  27475 29051 36001

table( meta$Time.point, meta$Cell.type.annotation,meta$NM_trajectory)
#       
#        Baso Ccr7_DC   Eos Erythroid Lymphoid  Mast   Meg Monocyte
# FALSE  5900     122   200       196      668  1033   628    12269
# TRUE   4255     119   355       490      189  1551  1468     6913
       
#       Neutrophil   pDC Undifferentiated
# FALSE          4    64            17276
# TRUE       22233    37            54917
#       
#        Baso Ccr7_DC   Eos Erythroid Lymphoid  Mast   Meg Monocyte
# FALSE  5900     122   200       196      668  1033   628    12269
# TRUE   4255     119   355       490      189  1551  1468     6913




############################################################
## 3. Subset to neutrophil-monocyte trajectory
############################################################

larry_nm <- subset(larry, subset = NM_trajectory)
rm(larry)
gc()

## Re-run normalization/PCA/UMAP inside the trajectory subset.
## This usually gives a cleaner bifurcation panel than using full-dataset UMAP.
# larry_nm <- NormalizeData(
#   larry_nm,
#   normalization.method = "LogNormalize",
#   scale.factor = 10000,
#   verbose = TRUE
# )

range(larry_nm[['RNA']]$counts)
#[1]   0.000 837.5551
larry_nm[['RNA']]$data <- log1p(larry_nm[['RNA']]$counts)
range(larry_nm[['RNA']]$data)
#[1] 0.00000 6.73168

larry_nm <- FindVariableFeatures(
  larry_nm,
  selection.method = "vst", # default
  nfeatures = 3000,
  verbose = TRUE
)

 
# Removed genes with correlation R > 0.1 to signature: Ube2c, Hmgb2, Hmgn2, Tuba1b, Ccnb1, Tubb5, Top2a, Tubb4b.
# calculate the correlation of the cell cycle signature genes with the eight cell cycle signature scores
# Get Seurat's default cell cycle gene sets for mouse
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
cell_cycle_genes_seurat <- unique(c(s_genes, g2m_genes))
length(cell_cycle_genes_seurat) # 97

cell_cycle_genes_mouse <- stringr::str_to_title(cell_cycle_genes_seurat)
cell_cycle_genes_mouse <- intersect(
  cell_cycle_genes_mouse,
  rownames(larry_nm)
)
length(cell_cycle_genes_mouse) # 93
head(cell_cycle_genes_mouse)

larry_nm <- CellCycleScoring(
  larry_nm,
  s.features = stringr::str_to_title(cc.genes$s.genes),
  g2m.features = stringr::str_to_title(cc.genes$g2m.genes),
  set.ident = FALSE
)
colnames(larry_nm@meta.data)
#  [1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"        
#  [4] "Library"              "Cell.barcode"         "Time.point"          
#  [7] "Starting.population"  "Cell.type.annotation" "Well"                
# [10] "SPRING.x"             "SPRING.y"             "cell_index0"         
# [13] "cell_id"              "NM_trajectory"        "CellCycleSig1"       
# [16] "S.Score"              "G2M.Score"            "Phase"        

hvg <- VariableFeatures(larry_nm)

expr_hvg <- GetAssayData(
  larry_nm,
  assay = "RNA",
  layer = "data"
)[hvg, , drop = FALSE]

cell_cycle_signature_cor <- vapply(hvg, function(gene) {  
  x <- as.numeric(expr_hvg[gene, ])  
  if (sd(x, na.rm = TRUE) == 0) {
    return(NA_real_)
  }  
  suppressWarnings(
    cor(
      x,
      larry_nm$S.Score,
      use = "pairwise.complete.obs",
      method = "spearman"
    )
  )
  
}, numeric(1))

cell_cycle_signature_cor <- sort(cell_cycle_signature_cor, decreasing = TRUE)

head(cell_cycle_signature_cor, 30)  # they all biologicall important, not to removed!!
    #       H2afy         Prtn3         Plac8          Ttf1           Mpo 
    # 0.3368575     0.2386996     0.2340606     0.2321832     0.2304042 
    #    Ccdc62          Ctsg        Unc13d           Myc         Elane 
    # 0.2260632     0.2094567     0.2062644     0.1989810     0.1970690 
    #      Cd34       Gm37214       Ptprcap 2810403D21Rik      Tmem176b 
    # 0.1825454     0.1772943     0.1678188     0.1555660     0.1522584 
    #       Fry         Muc13          Srgn          Car2        Bcl11a 
    # 0.1439423     0.1365795     0.1333781     0.1325042     0.1270497 
    #      Sox4         H2-Q7         Ccnd2         Mef2c           Ak4 
    # 0.1187358     0.1182859     0.1175382     0.1156264     0.1150816 
    #      Tcf4         Egfl7          Scin          Bcl2          Lmo4 
    # 0.1075007     0.1074148     0.1062428     0.1059850     0.1048503 

tail(cell_cycle_signature_cor, 30)
#     Tecpr1       Psap       Ctsd      Timp2      Anxa1     Wfdc17      Itgam 
# -0.2530740 -0.2554771 -0.2571720 -0.2591316 -0.2662451 -0.2662870 -0.2671655 
#      Gpnmb      Itgb2       Sat1     Clec4n       Ccr1       Mmp8      Mpeg1 
# -0.2673966 -0.2682269 -0.2712195 -0.2733006 -0.2755371 -0.2798000 -0.2819408 
#      Ahnak      Itm2b       Cd52       Klf6       Pirb    Lilrb4a        Vim 
# -0.2820763 -0.2835312 -0.2866918 -0.2867993 -0.3012919 -0.3062674 -0.3084360 
#    Alox5ap       Lyz2       Ccl6     Tyrobp        Grn        Gsn     Lgals3 
# -0.3095076 -0.3152042 -0.3206376 -0.3218464 -0.3220165 -0.3423992 -0.3642211 
#     Fcer1g     S100a6 
# -0.3789861 -0.3988530 
saveRDS(cell_cycle_signature_cor, file.path(outdir, "cell_cycle_signature_cor_3khvg.rds"))

hist(cell_cycle_signature_cor, 100, main = "Spearman correlation: S.Score vs 3k HVGs")
abline(v = 0.3, col = "red")
abline(v = -0.3, col = "red")
dev.copy2pdf(file = file.path(outdir, "hist_Sscore_cor_3khvg.pdf"))


cell_cycle_cor_threshold <- 0.3
hvg_correlated_with_cc <- names(
  cell_cycle_signature_cor[
    abs(cell_cycle_signature_cor) > cell_cycle_cor_threshold
  ]
)
length(hvg_correlated_with_cc) # [1] 13
# [1] "H2afy"   "Pirb"    "Lilrb4a" "Vim"     "Alox5ap" "Lyz2"    "Ccl6"   
# [8] "Tyrobp"  "Grn"     "Gsn"     "Lgals3"  "Fcer1g"  "S100a6" 
# Within NM-trajectory cells, none of the published 8 cell-cycle signature genes were retained among the 3,000 HVGs. 
# Correlation of HVGs with Seurat S.Score was centered near zero; only 13 genes exceeded |Spearman rho| > 0.3, 
# and most were myeloid/lineage-associated rather than canonical cell-cycle genes. 
# Therefore, we retained the NM-subset HVGs for the primary PCA/trajectory analysis and used removal of these 13 genes as a sensitivity analysis.

larry_nm <- ScaleData(
  larry_nm,
  features = setdiff(hvg, cell_cycle_signature),
  verbose = TRUE
)


larry_nm <- RunPCA(
  larry_nm,
  features = VariableFeatures(larry_nm),
  npcs = 50,
  verbose = TRUE
)

ElbowPlot(larry_nm, ndims = 50) +
  geom_hline(yintercept = 1,   linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
  annotate("text", x = 45, y = 1,   label = "1",   vjust = -0.5, color = "red") +
  annotate("text", x = 45, y = 1.5, label = "1.5", vjust = -0.5, color = "red")
dev.copy2pdf(file = file.path(outdir, "Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_ElbowPlot.pdf"))

larry_nm <- RunUMAP(
  larry_nm,
  dims = 1:npcs,
  reduction = "pca",
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = TRUE
)

saveRDS(
  larry_nm,
  file = file.path(outdir, "Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds")
)


# scale.data is centered/z-scored, can be negative, and often contains only HVGs. 
# In my case it was generated after selecting/scaling HVGs, so it is not a full expression matrix for all genes. 
# It is useful for PCA/heatmap-like standardized display, not for biological marker intensity.
# check the dimention of three layers:
# In Seurat v5, use the Assays() and GetAssayData() accessors to explore 'counts', 'data', and 'scale.data' layers:
dim(GetAssayData(larry_nm, layer = "data")) # 25289 92527
dim(GetAssayData(larry_nm, layer = "scale.data")) # 3000 92527
dim(GetAssayData(larry_nm, layer = "counts")) # 25289 92527

range(GetAssayData(larry_nm, layer = "data")) #  0.00000 6.73168
range(GetAssayData(larry_nm, layer = "scale.data")) # -2.159904 10.000000
range(GetAssayData(larry_nm, layer = "counts")) #  0.0000 837.5551


## double check that PCA is not dominated by cell cycle:
pc_embed <- Embeddings(larry_nm, "pca")[, 1:30, drop = FALSE]

pc_cc_cor <- apply(pc_embed, 2, function(x) {
  cor(x, larry_nm$CellCycleSig1, use = "pairwise.complete.obs", method = "spearman")
})

pc_cc_cor

barplot(
  abs(pc_cc_cor),
  las = 2,
  main = "Absolute Spearman correlation: PCs vs CellCycleSig1",
  ylab = "|rho|",
  ylim = c(0, 0.3)
)
abline(h = 0.3, lty = 2, col = "red")
dev.copy2pdf(file = paste0(outdir, "barplot_noCellCycle.PCA_CellCycleSig1_cor.pdf"))

