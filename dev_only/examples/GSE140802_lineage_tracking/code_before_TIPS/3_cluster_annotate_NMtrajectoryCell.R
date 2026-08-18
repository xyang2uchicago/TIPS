# -- work on the logcounts matrix for all downstream tasks

# 1. Start from downloaded total-count-normalized matrix, SPRING coordinates, and annotations.
# 2. Subset to NM_trajectory cells.
# 3. Re-normalize with Seurat LogNormalize.
# 4. Select 3,000 HVGs within the NM subset.
# 5. Explicitly exclude the published 8 cell-cycle genes.
# 6. Confirm that the broader all-cell cell-cycle-correlated genes are not retained among NM HVGs.
# 7. Scale final HVGs, run PCA, then UMAP.
# (done locally here)
# 8. Leiden clustering on PCA/SNN graph
# 9. compare Leiden clusters with downloaded cell-type annotation
# (moving to midway3 to run the following steps)
#         ↓
# 10. Slingshot / PAGA / MST on Leiden clusters
# 11. algorithmically identify branch node + Neu/Mono terminal paths


# ssh xyang2@midway3.rcc.uchicago.edu
# module load python/anaconda-2022.05
# source activate /project/xyang2/software-packages/env/velocity_2022.05_xy

# rcchelp balance  
# rcchelp usage
# squeue -u xyang2
# squeue -p bigmem --state=PD | wc -l
# squeue -p caslake --state=PD | wc -l

# sinteractive -p amd --mem=180G  --account=pi-xyang2 --time=6:00:00  -c 1  
# sinteractive -p caslake  --cpus-per-task=1 --mem=180G --time=6:00:00  --account=pi-xyang2 
# sinteractive -p gpu --account=pi-xyang2 --gres=gpu:1 --mem=180GB  --time=6:00:00 -c 1  (used)
# cd /project/imoskowitz/xyang2/heart_dev/Atlas2_results
# R
## To quit your interactive job:
# exit or Ctrl-D



################ local run ########################
setwd('F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/')
# setwd('/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/')  # if midway3 running

library('SingleCellExperiment')
library(Seurat)
library(dplyr) 
# devtools::install_github("xyang2uchicago/BioTIP")
# library(BioTIP)
# packageVersion('BioTIP')  # 1.19.0
library(scuttle)
library(igraph)
library(pheatmap)
library(mclust)

library(SingleCellExperiment)


larry_nm <- readRDS( "Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds")
larry_nm
# An object of class Seurat 
# 25289 features across 92527 samples within 1 assay 
# Active assay: RNA (25289 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  3 dimensional reductions calculated: spring, pca, umap
 
if(!'Cell.type.clean' %in% colnames(larry_nm@meta.data)) {
  larry_nm$Cell.type.clean <- as.character(larry_nm$Cell.type.annotation)
  larry_nm$Cell.type.clean[larry_nm$Cell.type.clean == "Neutrophil"] <- "Neut"
  larry_nm$Cell.type.clean[larry_nm$Cell.type.clean == "Monocyte"] <- "Mono"
  larry_nm$Cell.type.clean[larry_nm$Cell.type.clean == "Undifferentiated"] <- "Unddiff"
  larry_nm$Cell.type.clean[larry_nm$Cell.type.clean == "outside_NM_trajectory"] <- "Other"

  cell_lvls <- sort(unique(larry_nm$Cell.type.clean[larry_nm$NM_trajectory]))
  larry_nm$Cell.type.clean <- factor(larry_nm$Cell.type.clean, levels = cell_lvls)
}

if('Fig6B_label' %in% colnames(larry_nm@meta.data))  larry_nm$Fig6B_label <- NULL

### Leiden clustering on PCA/SNN graph
########################################
npcs <- 30
cluster_res <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 0.5)
cluster_col <- paste0("leiden_r", gsub("\\.", "_", cluster_res))

larry_nm <- FindNeighbors(
  larry_nm,
  reduction = "pca",
  dims = 1:npcs,
  k.param = 20, # by default, k.param = 20
  verbose = TRUE
)

for(i in 1:length(cluster_res)){
  larry_nm <- FindClusters(
    larry_nm,
    resolution = cluster_res[i],
    algorithm = 4,       # Leiden
    random.seed = 1000,
    cluster.name = cluster_col[i],
    verbose = TRUE
  )
}
colnames(larry_nm@meta.data)
#  [1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"        
#  [4] "Library"              "Cell.barcode"         "Time.point"          
#  [7] "Starting.population"  "Cell.type.annotation" "Well"                
# [10] "SPRING.x"             "SPRING.y"             "cell_index0"         
# [13] "cell_id"              "NM_trajectory"        "Fig6B_label"         
# [16] "CellCycleSig1"        "Cell.type.clean"      "leiden_r0_2"         
# [19] "leiden_r0_4"          "leiden_r0_6"          "leiden_r0_8"         
# [22] "leiden_r1"            "leiden_r1_2"          "leiden_r0_5"
all(cluster_col %in% colnames(larry_nm@meta.data)) #T

for (i in cluster_col){
   cat(i, ": ", length(unique(larry_nm@meta.data[[i]])), "Clusters: ", table(larry_nm@meta.data[[i]]), "\n")
}

# leiden_r0_2 :  8 Clusters:  32787 22175 14234 10301 6010 4674 1824 522 
# leiden_r0_4 :  12 Clusters:  18850 14134 10428 10093 9573 9532 6621 6416 4353 1868 521 138 
# leiden_r0_6 :  18 Clusters:  11623 10714 10385 8422 8267 8163 6707 6344 4687 4119 3947 3530 1924 1857 930 526 245 137 
# leiden_r0_8 :  19 Clusters:  10332 10234 8030 6954 6798 6583 6387 5020 4823 4637 4517 4077 3638 3630 2623 1865 1717 525 137 
# leiden_r1 :  22 Clusters:  8624 8188 7336 6707 6362 5591 5138 4868 4822 4678 3942 3821 3597 3573 3325 3132 2936 1986 1858 1376 530 137 
# leiden_r1_2 :  22 Clusters:  8432 8318 7271 7030 6687 5665 4813 4786 4656 4272 3940 3855 3806 3588 3573 3335 2602 1995 1858 1378 530 137 
# leiden_r0_5 :  15 Clusters:  11627 10638 10336 10179 9088 7819 6881 6563 5112 4350 3741 2429 1869 1371 524 



saveRDS(larry_nm, file = "Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds")  # !!!!!!!!!!

### compare Leiden clusters with downloaded cell-type annotation
########################################
## plot umaps for each leiden cluster into one figure

p <- DimPlot(larry_nm, reduction = "umap", group.by = c( cluster_col,"Cell.type.clean"), 
      label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.2, raster = TRUE) & NoLegend()

pdf(file = "Leiden_clusters_UMAP.pdf")
print(p)
dev.off()


## plot heatmap for each leiden cluster vs cell type
p = list()

for(r in paste0("leiden_r", c('0_2', '0_4', '0_5', '0_6', '0_8', '1'))){
  tab_ct <- table(
    Leiden = larry_nm@meta.data[[r]],
    CellType = larry_nm$Cell.type.clean
  )

  prop_ct <- prop.table(tab_ct, margin = 1) # row-wise proportions, per Leiden cluster
  round(prop_ct, 2)

  ph <- pheatmap::pheatmap(
      prop_ct,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = paste0(r, " vs cell types"),
      silent = TRUE
    )
    
    p[[r]] <- ph$gtable
}

g <- gridExtra::arrangeGrob(
  grobs = unname(p),
  ncol = 2
)


pdf("Leiden_clusters_vs_celltypes_heatmaps.pdf", width = 14, height = 16)
grid::grid.draw(g)
dev.off()


## choose the best leiden cluster r=0.6 for trajectory analysis
## becasue until 0.6, Meg can be distinguished until leiden_0.6;  

r <- "leiden_r0_6"

md <- larry_nm@meta.data %>%
  as.data.frame() %>%
  mutate(
    leiden = as.character(.data[[r]]),
    Cell.type.clean = as.character(Cell.type.clean),
    Time.point = as.character(Time.point)
  )

cluster_summary <- md %>%
  group_by(leiden) %>%
  summarise(
    n = n(),
    
    frac_Und = mean(Cell.type.clean == "Und", na.rm = TRUE),
    frac_Neutrophil = mean(Cell.type.clean == "Neutrophil", na.rm = TRUE),
    frac_Monocyte = mean(Cell.type.clean == "Monocyte", na.rm = TRUE),
    
    frac_D2 = mean(Time.point == "2", na.rm = TRUE),
    frac_D4 = mean(Time.point == "4", na.rm = TRUE),
    frac_D6 = mean(Time.point == "6", na.rm = TRUE),
    
    major_celltype = names(sort(table(Cell.type.clean), decreasing = TRUE))[1],
    major_time = names(sort(table(Time.point), decreasing = TRUE))[1],
    purity_celltype = max(prop.table(table(Cell.type.clean))),
    purity_time = max(prop.table(table(Time.point))),
    
    .groups = "drop"
  ) %>%
  arrange(leiden)

cluster_summary$cluster_annotation <- ifelse(
  cluster_summary$purity_celltype >= 0.6,
  cluster_summary$major_celltype,
  ifelse(cluster_summary$frac_Und >= 0.4 & cluster_summary$frac_Monocyte >= 0.3, "Und/Mono",
  ifelse(cluster_summary$frac_Und >= 0.4 & cluster_summary$frac_Neutrophil >= 0.3, "Und/Neu",
  paste0("Mixed_", cluster_summary$major_celltype)))
)   
cluster_summary$cluster_annotation = paste0(cluster_summary$cluster_annotation, "_", cluster_summary$major_time)

cluster_summary %>% as.data.frame()
#    leiden     n frac_Und frac_Neutrophil frac_Monocyte     frac_D2    frac_D4    frac_D6 major_celltype major_time purity_celltype purity_time cluster_annotation
# 1       1 11623        0               0             0 0.540996300 0.25406522 0.20493848        Unddiff          2       0.9989676   0.5409963          Unddiff_2
# 2      10  4119        0               0             0 0.107793154 0.24496237 0.64724448           Mono          6       0.9065307   0.6472445             Mono_6
# 3      11  3947        0               0             0 0.013934634 0.27033190 0.71573347           Baso          6       0.9569293   0.7157335             Baso_6
# 4      12  3530        0               0             0 0.010764873 0.75552408 0.23371105           Neut          4       0.5056657   0.7555241       Mixed_Neut_4
# 5      13  1924        0               0             0 0.500000000 0.21673597 0.28326403        Unddiff          2       0.7811850   0.5000000          Unddiff_2
# 6      14  1857        0               0             0 0.051157781 0.17662897 0.77221325           Mast          6       0.8341411   0.7722132             Mast_6
# 7      15   930        0               0             0 0.000000000 0.11397849 0.88602151           Neut          6       0.8107527   0.8860215             Neut_6
# 8      16   526        0               0             0 0.102661597 0.31368821 0.58365019            Eos          6       0.5988593   0.5836502        Mixed_Eos_6
# 9      17   245        0               0             0 0.004081633 0.38367347 0.61224490        Unddiff          6       0.9510204   0.6122449          Unddiff_6
# 10     18   137        0               0             0 0.489051095 0.18248175 0.32846715        Ccr7_DC          2       0.8321168   0.4890511          Ccr7_DC_2
# 11      2 10714        0               0             0 0.003920105 0.22092589 0.77515400           Neut          6       0.9621990   0.7751540             Neut_6
# 12      3 10385        0               0             0 0.607896004 0.28483390 0.10727010        Unddiff          2       0.9996148   0.6078960          Unddiff_2
# 13      4  8422        0               0             0 0.107100451 0.31346474 0.57943481           Neut          6       0.9576110   0.5794348             Neut_6
# 14      5  8267        0               0             0 0.083464376 0.70509254 0.21144309        Unddiff          4       0.9655256   0.7050925          Unddiff_4
# 15      6  8163        0               0             0 0.249908122 0.35097391 0.39911797        Unddiff          6       0.8809261   0.3991180          Unddiff_6
# 16      7  6707        0               0             0 0.947368421 0.03951096 0.01312062        Unddiff          2       0.9664530   0.9473684          Unddiff_2
# 17      8  6344        0               0             0 0.230611602 0.38634931 0.38303909           Mono          4       0.4985813   0.3863493       Mixed_Mono_4
# 18      9  4687        0               0             0 0.355664604 0.18028590 0.46404950        Unddiff          6       0.5824621   0.4640495    Mixed_Unddiff_6


table(larry_nm@meta.data[["leiden_r0_6"]],larry_nm$Cell.type.clean)
  #     Baso Ccr7_DC   Eos Erythroid Lymphoid  Mast   Meg  Mono  Neut   pDC Unddiff
  # 1      0       0     0         0        1     0     0     1    10     0   11611
  # 2      0       0     0         0        0     0     0     0 10309     0     405
  # 3      0       0     0         0        4     0     0     0     0     0   10381
  # 4      0       0     0         0        0     0     0     0  8065     0     357
  # 5      0       0     0         0        0     0     0     0   285     0    7982
  # 6      0       0     0         0        0     0     0     0   972     0    7191
  # 7      1       3     0         0      184     0     0     0     1    36    6482
  # 8      0       0     0         0        0     0     0  3163    22     0    3159
  # 9      2       0     0       486        0     0  1467     0     2     0    2730
  # 10     0       0     0         1        0     0     0  3734     1     1     382
  # 11  3777       2    38         0        0     0     1     4     6     0     119
  # 12     0       0     0         0        0     0     0     0  1785     0    1745
  # 13   415       0     0         0        0     0     0     0     6     0    1503
  # 14    12       0     0         2        0  1549     0     1     1     0     292
  # 15     0       0     0         0        0     0     0    10   754     0     166
  # 16    47       0   315         1        0     2     0     0     5     0     156
  # 17     1       0     2         0        0     0     0     0     9     0     233
  # 18     0     114     0         0        0     0     0     0     0     0      23


## manually update according to the % of cell.type.clean >20%
cluster_summary[which(cluster_summary$leiden == 7), 'cluster_annotation'] <- "Unddiff/Lym_2"
cluster_summary[which(cluster_summary$leiden == 13), 'cluster_annotation'] <- "Unddiff/Baso_2"

larry_nm$leiden_r0_6_anno <- cluster_summary$cluster_annotation[
  match(
    larry_nm$leiden_r0_6,
    cluster_summary$leiden
  )
] 

table( larry_nm$leiden_r0_6_anno,  larry_nm$Cell.type.clean)
               
  #                  Baso Ccr7_DC   Eos Erythroid Lymphoid  Mast   Meg  Mono  Neut   pDC Unddiff
  # Baso_6           3777       2    38         0        0     0     1     4     6     0     119
  # Ccr7_DC_2           0     114     0         0        0     0     0     0     0     0      23
  # Mast_6             12       0     0         2        0  1549     0     1     1     0     292
  # Mixed_Eos_6        47       0   315         1        0     2     0     0     5     0     156
  # Mixed_Mono_4        0       0     0         0        0     0     0  3163    22     0    3159
  # Mixed_Neut_4        0       0     0         0        0     0     0     0  1785     0    1745
  # Mixed_Unddiff_6     2       0     0       486        0     0  1467     0     2     0    2730
  # Mono_6              0       0     0         1        0     0     0  3734     1     1     382
  # Neut_6              0       0     0         0        0     0     0    10 19128     0     928
  # Unddiff/Baso_2    415       0     0         0        0     0     0     0     6     0    1503
  # Unddiff/Lym_2       1       3     0         0      184     0     0     0     1    36    6482
  # Unddiff_2           0       0     0         0        5     0     0     1    10     0   21992
  # Unddiff_4           0       0     0         0        0     0     0     0   285     0    7982
  # Unddiff_6           1       0     2         0        0     0     0     0   981     0    7424



saveRDS(larry_nm, file = "Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds")   # !!!!!!!!!!



############## spring map , umap, river plot of cell-type composition (USED) , 
## and stacked area proportion plot (NOT used) #####################
########################################
library(dplyr)
library(ggplot2)
library(scales)
library(ggalluvial)

for(r in c("leiden_r0_6", "Cell.type.clean","leiden_r0_6_anno", "leiden_r1")){

  ## 1. define ONE shared cluster order and palette
  cluster_levels <- sort(unique(as.character(larry_nm@meta.data[[r]])))
  cluster_levels <- as.character(cluster_levels)

  cluster_pal <- setNames(
    scales::hue_pal()(length(cluster_levels)),
    cluster_levels
  )

 plot_df <- larry_nm@meta.data %>%
  as.data.frame() %>%
  dplyr::mutate(
    Time.point = as.numeric(as.character(Time.point)),
    cluster = factor(as.character(.data[[r]]), levels = names(cluster_pal))
  ) %>%
  dplyr::count(Time.point, cluster, name = "n") %>%
  dplyr::group_by(Time.point) %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::ungroup()

  p_area <- ggplot(
    plot_df,
    aes(x = Time.point, y = freq, fill = cluster, group = cluster)
  ) +
    geom_area(position = "stack", color = NA, alpha = 0.95) +
    scale_fill_manual(values = cluster_pal, drop = FALSE) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    scale_x_continuous(
      breaks = sort(unique(plot_df$Time.point)),
      labels = paste0("D", sort(unique(plot_df$Time.point)))
    ) +
    labs(
      x = "Time point",
      y = "Fraction of cells",
      fill = r,
      title = paste0("Cluster composition over time: ", r)
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  pdf(paste0("stacked_area_proportion_plot_", r, ".pdf"))
  print(p_area)
  dev.off()
  

  #if(r=="leiden_r0_6_anno"){
    p_umap <- DimPlot(
      larry_nm,
      reduction = "umap",
      group.by = r,
      cols = cluster_pal,
      label = TRUE,
      repel = TRUE,
      pt.size = 0.2,
      raster = TRUE
    ) + ggtitle("UMAP") +
      NoLegend()

    p_string <- DimPlot(
      larry_nm,
      reduction = "spring",
      group.by = r,
      cols = cluster_pal,
      label = TRUE,
      repel = TRUE,
      pt.size = 0.2,
      raster = TRUE
    ) + ggtitle("Spring map") +
      NoLegend()


  river_time_df <- larry_nm@meta.data %>%
    as.data.frame() %>%
    dplyr::mutate(
      Time.point = factor(
        paste0("D", Time.point),
        levels = c("D2", "D4", "D6")
      ),
      cluster = factor(
      as.character(.data[[r]]),
      levels = cluster_levels
    )
    ) %>%
    dplyr::count(Time.point, cluster, name = "n") %>%
    dplyr::group_by(Time.point) %>%
    dplyr::mutate(freq = n / sum(n)) %>%
    dplyr::ungroup()

  p_river_time <- ggplot(
    river_time_df,
    aes(
      x = Time.point,
      stratum = cluster,
      alluvium = cluster,
      y = freq,
      fill = cluster,
      label = cluster
    )
  ) +
    geom_flow(alpha = 0.85, color = "grey80", linewidth = 0.15) +
    geom_stratum(width = 0.22, color = "grey30") +
    geom_text(
      stat = "stratum",
      size = 3,
      color = "black"
    ) +
    scale_fill_manual(values = cluster_pal, drop = FALSE) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    labs(
      title = paste0("Cluster composition over time: ", r),
      x = "Time point",
      y = "Fraction of cells",
      fill = r
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

#  p_river_time

#   ggsave(
#     paste0("riverplot_", r, "_over_timepoint.pdf"),
#     p_river_time,
#     width = 7,
#     height = 5,
#     device = cairo_pdf
#   )
# }


    pdf(paste0("spring_umap_riverplot_", r, ".pdf"), width=15, height=5)
    print(p_string | p_umap |  p_river_time)
    dev.off()

}


######### river plot of cell-type composition over leiden_r0_6 clusters #########
library(dplyr)
library(ggplot2)
library(ggalluvial)

r <- "leiden_r0_6"

## same cluster color palette as UMAP/Slingshot
cluster_levels <- as.character(sort(as.numeric(unique(as.character(larry_nm@meta.data[[r]])))))

cluster_pal <- setNames(
  scales::hue_pal()(length(cluster_levels)),
  cluster_levels
)

river_df <- larry_nm@meta.data %>%
  as.data.frame() %>%
  dplyr::mutate(
    CellType = as.character(Cell.type.clean),
    Cluster = factor(as.character(.data[[r]]), levels = cluster_levels)
  ) %>%
  dplyr::count(CellType, Cluster, name = "n")

p_river_ct_cluster <- ggplot(
  river_df,
  aes(axis1 = CellType, axis2 = Cluster, y = n)
) +
  geom_alluvium(
    aes(fill = Cluster),
    width = 0.18,
    alpha = 0.85
  ) +
  geom_stratum(
    width = 0.18,
    fill = "grey95",
    color = "grey40"
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3
  ) +
  scale_fill_manual(values = cluster_pal, drop = FALSE) +
  scale_x_discrete(
    limits = c("Cell type", "Leiden cluster"),
    expand = c(0.08, 0.08)
  ) +
  labs(
    title = paste0("Cell-type composition of ", r, " clusters"),
    y = "Number of cells",
    fill = r
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_river_ct_cluster

ggsave(
  paste0("riverplot_celltype_over_", r, ".pdf"),
  p_river_ct_cluster,
  width = 8,
  height = 5,
  device = cairo_pdf
)
dev.off()