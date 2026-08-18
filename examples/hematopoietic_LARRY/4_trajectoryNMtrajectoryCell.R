# -- work on the logcounts matrix for all downstream tasks

# 1. Start from downloaded total-count-normalized matrix, SPRING coordinates, and annotations.
# 2. Subset to NM_trajectory cells.
# 3. Re-normalize with Seurat LogNormalize.
# 4. Select 3,000 HVGs within the NM subset.
# 5. Explicitly exclude the published 8 cell-cycle genes.
# 6. Confirm that the broader all-cell cell-cycle-correlated genes are not retained among NM HVGs.
# 7. Scale final HVGs, run PCA, then UMAP.
# (done locally )
# 8. Leiden clustering on PCA/SNN graph
# 9. compare Leiden clusters with downloaded cell-type annotation -> select leiden_0.6
#  # (done locally , now , moveing to midway3 to run the following steps)
# 10. convert Seurat object to SingleCellExperiment object -> root cluster = '4' -> Slingshot
# 11 / PAGA / MST on Leiden clusters
#         ↓
# 11. algorithmically identify branch node + Neu/Mono terminal paths


########################################
### Slingshot / PAGA / MST on Leiden clusters
## running on midway3 due to "matrixStats" version issue
########################################

# scp -p -r F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/*.rds xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/.
# scp -p -r F:\projects\TIPS\source\GSE140802_lineage_tracking\4*.R   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/code/

# ssh xyang2@midway3.rcc.uchicago.edu
# module load python/anaconda-2022.05
# source activate /project/xyang2/software-packages/env/velocity_2025Feb_xy/

# rcchelp balance  
# rcchelp usage
# squeue -u xyang2
# squeue -p bigmem --state=PD | wc -l
# squeue -p caslake --state=PD | wc -l

# sinteractive -p bigmem  --mem=500G  --account=pi-xyang2 --time=1:00:00  -c 1  
# sinteractive -p caslake  --cpus-per-task=1 --mem=180G --time=6:00:00  --account=pi-xyang2 
# sinteractive -p gpu --account=pi-xyang2 --gres=gpu:1 --mem=180GB  --time=6:00:00 -c 1  (used)
# cd /project/imoskowitz/xyang2/heart_dev/Atlas2_results
# R
## To quit your interactive job:
# exit or Ctrl-D

#setwd('F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/')
setwd('/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/')

pkgs <- c(
   "slingshot", "TrajectoryUtils", "MatrixGenerics",
   "DelayedMatrixStats", "sparseMatrixStats", "matrixStats"
)

sapply(pkgs, packageVersion)  # local / midway3
# $slingshot
# [1] 2 4 0 /2 10  0

# $TrajectoryUtils
# [1] 1 4 0 /  1 10  1

# $MatrixGenerics
# [1]  1 10  0/ 1 14  0

# $DelayedMatrixStats
# [1]  1 20  0 /1 24  0

# $sparseMatrixStats
# [1]  1 10  0 / 1 14 0

# $matrixStats
# [1] 1 4 1  / 1 5 0 



library(SingleCellExperiment)
library(slingshot)
library(dplyr) 
library(Seurat)

step1 = FALSE
if(step1){
    larry_nm <- readRDS( "Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds")
    larry_nm
    # An object of class Seurat 
    # 25289 features across 92527 samples within 1 assay 
    # Active assay: RNA (25289 features, 3000 variable features)
    #  3 layers present: counts, data, scale.data
    #  3 dimensional reductions calculated: spring, pca, umap
    colnames(larry_nm@meta.data)
#      [1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"
#  [4] "Library"              "Cell.barcode"         "Time.point"
#  [7] "Starting.population"  "Cell.type.annotation" "Well"
# [10] "SPRING.x"             "SPRING.y"             "cell_index0"
# [13] "cell_id"              "NM_trajectory"        "CellCycleSig1"
# [16] "S.Score"              "G2M.Score"            "Phase"
# [19] "leiden_r0_2"          "leiden_r0_4"          "leiden_r0_6"
# [22] "leiden_r0_8"          "leiden_r1"            "leiden_r1_2"
# [25] "leiden_r0_5"          "Cell.type.clean"      "leiden_r0_6_anno"


  ################################################################################
    ### first of all, convert Seurat object to SingleCellExperiment object ##########
    ################################################################################
    # scale.data is centered/z-scored, can be negative, and often contains only HVGs. 
    # In my case it was generated after selecting/scaling HVGs, so it is not a full expression matrix for all genes. 
    # It is useful for PCA/heatmap-like standardized display, not for biological marker intensity.
    # check the dimention of three layers:
    # In Seurat v5, use the Assays() and GetAssayData() accessors to explore 'counts', 'data', and 'scale.data' layers:
    dim(GetAssayData(larry_nm, layer = "data")) # 25289 92527
    dim(GetAssayData(larry_nm, layer = "scale.data")) # 3000 92527
    dim(GetAssayData(larry_nm, layer = "counts")) # 25289 92527

    range(GetAssayData(larry_nm, layer = "data")) #  0.00000  6.73168
    range(GetAssayData(larry_nm, layer = "scale.data")) # -2.159904 10.000000
    range(GetAssayData(larry_nm, layer = "counts")) #  0.0000 837.5551


    npcs <- 30
    DefaultAssay(larry_nm) <- "RNA"

    counts_mat <- GetAssayData(
    larry_nm,
    assay = "RNA",
    layer = "counts" # note this is the total-count-normalized matrix downladed from GEO!
    )

    logcounts_mat <- GetAssayData(
    larry_nm,
    assay = "RNA",
    layer = "data"
    )

    sce<- SingleCellExperiment(
    assays = list(
        counts = counts_mat,
        logcounts = logcounts_mat
    ),
    colData = larry_nm@meta.data
    )

    reducedDim(sce, "PCA") <- Embeddings(larry_nm, "pca")[, 1:npcs, drop = FALSE]

    if ("umap" %in% Reductions(larry_nm)) {
    reducedDim(sce, "UMAP") <- Embeddings(larry_nm, "umap")
    }

    if ("spring" %in% Reductions(larry_nm)) {
    reducedDim(sce, "SPRING") <- Embeddings(larry_nm, "spring")
    }

    colnames(colData(sce))
    #  [1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"
    #  [4] "Library"              "Cell.barcode"         "Time.point"
    #  [7] "Starting.population"  "Cell.type.annotation" "Well"
    # [10] "SPRING.x"             "SPRING.y"             "cell_index0"
    # [13] "cell_id"              "NM_trajectory"        "Fig6B_label"
    # [16] "CellCycleSig1"        "Cell.type.clean"      "leiden_r0_2"
    # [19] "leiden_r0_4"          "leiden_r0_6"          "leiden_r0_8"
    # [22] "leiden_r1"            "leiden_r1_2"          "leiden_r0_5"
    # [25] "leiden_r0_6_anno"      
    saveRDS(sce, file = "sce.rds")   #!!!!!!!!!!!!!



    sce
    # class: SingleCellExperiment 
    # dim: 25289 92527 
    # metadata(0):
    # assays(2): counts logcounts
    # rownames(25289): 0610006L08Rik 0610007P14Rik ... mt-Nd5 mt-Nd6
    # rowData names(0):
    # colnames(92527): cell_0 cell_1 ... cell_130883 cell_130886
    # colData names(27): orig.ident nCount_RNA ... leiden_r0_6_anno leiden_alg
    # reducedDimNames(3): PCA UMAP SPRING
    # mainExpName: NULL
    # altExpNames(0):

    cat("step1 done\n")
} # end of step1

################################################################################
### Option 1: slingshot algorithmically identify branch node + Neu/Mono terminal paths
################################################################################
step2 = TRUE
rerunslingshot = TRUE

# # find the earlest undiff root
# x = table(colData(sce)$leiden_r0_8, colData(sce)$Time.point)
# y = x[,1]/apply(x,1,sum)
# y
           # # 1            2            3            4            5            6
# # 0.5595238095 0.6052374438 0.0860523039 0.9469370147 0.2230067667 0.0059243506
           # # 7            8            9           10           11           12
# # 0.2071395021 0.0005976096 0.0491395397 0.0715980160 0.5160504760 0.1758646063
          # # 13           14           15           16           17           18
# # 0.1069268829 0.0126721763 0.3804803660 0.0552278820 0.0407687828 0.1009523810
          # # 19
# # 0.4890510949
# x = table(colData(sce)$leiden_r0_8, colData(sce)$Cell.type.clean)
# z = x[,'Unddiff']/apply(x,1,sum)
# z
         # # 1          2          3          4          5          6          7
# # 0.99970964 0.99960915 0.98468244 0.96764452 0.46734334 0.03691326 0.86268984
         # # 8          9         10         11         12         13         14
# # 0.06553785 0.02716152 0.08561570 0.87381005 0.10694138 0.09400770 0.50303030
        # # 15         16         17         18         19
# # 0.99008769 0.16139410 0.17880023 0.29714286 0.16788321
# y+z
         # # 1          2          3          4          5          6          7
# # 1.55923345 1.60484659 1.07073474 1.91458154 0.69035010 0.04283761 1.06982934
         # # 8          9         10         11         12         13         14
# # 0.06613546 0.07630106 0.15721372 1.38986053 0.28280598 0.20093458 0.51570248
        # # 15         16         17         18         19
# # 1.37056805 0.21662198 0.21956902 0.39809524 0.65693431

roots = c("4", "Unddiff", "3", "Unddiff_2")
names(roots) = c("leiden_r0_8","Cell.type.clean", "leiden_r0_6", "leiden_r0_6_anno")

if(step2){
  sce   = readRDS("sce.rds")
  range(counts(sce)) #  0.0000 837.5551
  range(logcounts(sce)) #  0.00000 6.73168
  for( r in c("leiden_r0_8"))  #,"Cell.type.clean", "leiden_r0_6", "leiden_r0_6_anno") )
  { 
    md <- colData(sce) %>%
      as.data.frame() %>%
      mutate(
        leiden = as.character(.data[[r]]),
        Cell.type.clean = as.character(.data[['Cell.type.clean']]),
        Time.point = as.character(.data[["Time.point"]])
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
    write.table(cluster_summary, file = paste0("cluster_summary_", r, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Choose the root algorithmically as the Leiden cluster most enriched for D2 undifferentiated cells:
    # root_cluster <- cluster_summary %>%
    #   filter(n >= 200) %>%
    #   slice_max(frac_D2, n = 1, with_ties = FALSE) %>%
    #   pull(leiden)
    root_cluster = roots[r]
    cat("root_cluster:", root_cluster, "\n") # "4"   Und    / Und/Lym/pDC_2

    rd <- as.matrix(reducedDim(sce, "PCA"))
    cl <- as.character(sce[[r]])
    ok <- !is.na(sce[[r]]) &
    rowSums(!is.finite(rd)) == 0
    rd2 <- rd[ok, , drop = FALSE]
    cl2 <- cl[ok]
    
    if(rerunslingshot){  #!!!!!!!!
        sds <- slingshot(
            data = rd2,
            clusterLabels = cl2,
            start.clus = as.character(root_cluster),
            approx_points = 150
        )
    } else sds <- readRDS(paste0("sds_slingshot_", r, ".rds"))

    lineages <- lapply(slingLineages(sds), as.character)
    lineages

    if(rerunslingshot){
        saveRDS(lineages, file = paste0("slingshot_lineages_", r, ".rds"))
        saveRDS(sds, file = paste0("sds_slingshot_", r, ".rds"))
        } else lineages <- readRDS(paste0("slingshot_lineages_", r, ".rds"))
    #    set.seed(1)
    #    keep <- sample(seq_len(nrow(rd2)), min(20000, nrow(rd2)))
    
    # define color palette for the clusters
    cluster_levels <- sort(unique(as.character(colData(sce)[[r]])))
    cluster_levels <- as.character(cluster_levels)
    cluster_pal <- setNames(scales::hue_pal()(length(cluster_levels)),cluster_levels)

    sds_plot <- as.SlingshotDataSet(sds)
    cl_fac <- factor(cl2, levels = names(cluster_pal))
    cols <- setNames(seq_along(levels(cl_fac)), levels(cl_fac))
   
    pdf(paste0("Slingshot_all_lineages_", r, "_PCA.pdf"))
    plot(
        rd2[, 1:2],
        col = cluster_pal[as.character(cl_fac)], ## use the color palette for the clusters
        pch = 16,
        cex = 0.25,
        asp = 1,
        main = paste0("Slingshot all lineages: ", r),
        xlab = "PC1",
        ylab = "PC2"
    )
    lines(sds_plot, type = "lineages", lwd = 3, col = "black")
    lines(sds_plot, type = "curves", lwd = 2, col = "red")
    legend(
        "bottomright",
        legend = names(cols),
        col = cols,
        pch = 16,
        pt.cex = 1.2,
        cex = 0.75,
        bty = "n",
        title = r
    )
    dev.off()

    cat("step2 done\n")
  } # end of for( r in c("Cell.type.clean", "leiden_r0_6", "leiden_r0_6_anno") )
} # end of step2



for(r in c("leiden_r0_8")){ #, "Cell.type.clean", "leiden_r0_6", "leiden_r0_6_anno") ) {
  print( slingLineages(sds))
}
# $Lineage1
# [1] "4"  "2"  "1"  "3"  "7"  "12" "9"  "6"  "8"

# $Lineage2
# [1] "4"  "2"  "1"  "3"  "7"  "14" "5"  "13"

# $Lineage3
# [1] "4"  "2"  "11" "16" "10" "18"

# $Lineage4
# [1] "4"  "2"  "1"  "15"

# $Lineage5
# [1] "4"  "2"  "11" "17"

# $Lineage6
# [1] "4"  "19"

# $Lineage1
# [1] "Unddiff_2"    "Unddiff_4"    "Unddiff_6"    "Mixed_Neut_4" "Mixed_Mono_4"
# [6] "Mono_6"

# $Lineage2
# [1] "Unddiff_2" "Unddiff_4" "Unddiff_6" "Neut_6"

# $Lineage3
# [1] "Unddiff_2"       "Unddiff/Baso_2"  "Mixed_Unddiff_6"

# $Lineage4
# [1] "Unddiff_2"      "Unddiff/Baso_2" "Mast_6"

# $Lineage5
# [1] "Unddiff_2"      "Unddiff/Baso_2" "Baso_6"

# $Lineage6
# [1] "Unddiff_2"     "Unddiff/Lym_2" "Ccr7_DC_2"

# $Lineage7
# [1] "Unddiff_2"      "Unddiff/Baso_2" "Mixed_Eos_6"


# $Lineage1
# [1] "Unddiff_2"    "Unddiff_4"    "Unddiff_6"    "Mixed_Neut_4" "Mixed_Mono_4"
# [6] "Mono_6"

# $Lineage2
# [1] "Unddiff_2" "Unddiff_4" "Unddiff_6" "Neut_6"

# $Lineage3
# [1] "Unddiff_2"       "Unddiff/Baso_2"  "Mixed_Unddiff_6"

# $Lineage4
# [1] "Unddiff_2"      "Unddiff/Baso_2" "Mast_6"

# $Lineage5
# [1] "Unddiff_2"      "Unddiff/Baso_2" "Baso_6"

# $Lineage6
# [1] "Unddiff_2"     "Unddiff/Lym_2" "Ccr7_DC_2"

# $Lineage7
# [1] "Unddiff_2"      "Unddiff/Baso_2" "Mixed_Eos_6"


# $Lineage1
# [1] "Unddiff_2"    "Unddiff_4"    "Unddiff_6"    "Mixed_Neut_4" "Mixed_Mono_4"
# [6] "Mono_6"

# $Lineage2
# [1] "Unddiff_2" "Unddiff_4" "Unddiff_6" "Neut_6"

# $Lineage3
# [1] "Unddiff_2"       "Unddiff/Baso_2"  "Mixed_Unddiff_6"

# $Lineage4
# [1] "Unddiff_2"      "Unddiff/Baso_2" "Mast_6"

# $Lineage5
# [1] "Unddiff_2"      "Unddiff/Baso_2" "Baso_6"

# $Lineage6
# [1] "Unddiff_2"     "Unddiff/Lym_2" "Ccr7_DC_2"

# $Lineage7
# [1] "Unddiff_2"      "Unddiff/Baso_2" "Mixed_Eos_6"

step3 = TRUE  # this step can be rerun locally
color_strategy=c('Neu_mon_focused', 'all_cleuster')
if (step3) {
  
  library(ggplot2)
  library(dplyr)
  library(ggrastr)
  library(grid)
  library(ggrepel)
  library(slingshot)

  sce <- readRDS("sce.rds")
  
  common_prefix2 <- function(a, b) {
    n <- min(length(a), length(b))
    same <- a[seq_len(n)] == b[seq_len(n)]
    if (all(same)) return(a[seq_len(n)])
    k <- which(!same)[1] - 1
    if (k == 0) return(character(0))
    a[seq_len(k)]
  }

  for (r in c("leiden_r0_8")){ #c("leiden_r0_6", "Cell.type.clean", "leiden_r0_6_anno")) {
    
    message("Plotting: ", r)
    
    cd <- colData(sce) %>%
      as.data.frame() %>%
      mutate(
        leiden = as.character(.data[[r]]),
        Cell.type.clean = as.character(Cell.type.clean),
        Time.point = as.character(Time.point)
      )
    
    ## -------------------------------------------------------
    ## 1. cluster summary for current r
    ## -------------------------------------------------------
    cluster_summary <- cd %>%
      group_by(leiden) %>%
      summarise(
        n = n(),
        frac_Und = mean(Cell.type.clean == "Unddiff", na.rm = TRUE),
        frac_Neutrophil = mean(Cell.type.clean == "Neut", na.rm = TRUE),
        frac_Monocyte = mean(Cell.type.clean == "Mono", na.rm = TRUE),
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
    
    ## -------------------------------------------------------
    ## 2. load Slingshot output
    ## -------------------------------------------------------
    sds <- readRDS(paste0("sds_slingshot_", r, ".rds"))
    lineages <- readRDS(paste0("slingshot_lineages_", r, ".rds"))
    
    terminal_clusters <- vapply(
      lineages,
      function(x) tail(x, 1),
      character(1)
    )
    
    lineage_tbl <- data.frame(
      lineage = seq_along(lineages),
      terminal = terminal_clusters,
      path = vapply(lineages, paste, collapse = " -> ", character(1))
    ) %>%
      left_join(
        cluster_summary,
        by = c("terminal" = "leiden")
      ) %>%
      mutate(
        terminal_class = case_when(
          frac_Neutrophil >= 0.50 & frac_D6 > 0.50 ~ "Neutrophil_terminal",
          frac_Monocyte >= 0.50 & frac_D6 > 0.50 ~ "Monocyte_terminal",
          TRUE ~ "Other_terminal"
        ),
        lineage_id = as.integer(lineage),
        terminal = as.character(terminal),
        lineage_label = ifelse(
          terminal_class == "Other_terminal",
          paste0("Other_", major_celltype),
          terminal_class
        )
      )
    
    print(lineage_tbl)
    write.csv(lineage_tbl, file = paste0("lineage_tbl_", r, ".csv"), row.names = FALSE)
    ## -------------------------------------------------------
    ## 3. identify Neu/Mono branch node
    ## -------------------------------------------------------
    candidate_ids <- lineage_tbl$lineage_id[
      lineage_tbl$terminal_class %in% c(
        "Neutrophil_terminal",
        "Monocyte_terminal"
      )
    ]
    
    has_neu <- any(lineage_tbl$terminal_class == "Neutrophil_terminal")
    has_mono <- any(lineage_tbl$terminal_class == "Monocyte_terminal")
    
    if (length(candidate_ids) >= 2 && has_neu && has_mono) {
      common_path <- Reduce(common_prefix2, lineages[candidate_ids])
      branch_node <- if (length(common_path) > 0) tail(common_path, 1) else NA_character_
    } else {
      branch_node <- NA_character_
    }
    
    branch_node <- as.character(branch_node)
    message("Branch node: ", branch_node)
    
    bifurcation_detected <- has_neu & has_mono & !is.na(branch_node)
    message("Bifurcation detected: ", bifurcation_detected)
    
    ## -------------------------------------------------------
    ## 4. prepare UMAP dataframe
    ## -------------------------------------------------------
    um <- as.matrix(reducedDim(sce, "UMAP"))[, 1:2]
    
    plot_df <- data.frame(
      UMAP_1 = um[, 1],
      UMAP_2 = um[, 2],
      leiden = as.character(cd[[r]]),
      celltype = as.character(cd$Cell.type.clean)
    )
     
    neu_terms <- lineage_tbl$terminal[
      lineage_tbl$terminal_class == "Neutrophil_terminal"
    ]
    
    mono_terms <- lineage_tbl$terminal[
      lineage_tbl$terminal_class == "Monocyte_terminal"
    ]
    
    other_terms <- lineage_tbl$terminal[
      lineage_tbl$terminal_class == "Other_terminal"
    ]
    
    # set the consistent cluster color palette
    cluster_levels <- if (r == "leiden_r0_6") {
      as.character(sort(as.numeric(unique(as.character(cd[[r]])))))
    } else  sort(unique(as.character(cd[[r]])))

    if(color_strategy != 'Neu_mon_focused') {
        pdf_name <- paste0("p_sling_umap_", r, "_alllabels.pdf")
        color_var <- "leiden"
        cluster_pal <- setNames(
          scales::hue_pal()(length(cluster_levels)),
          cluster_levels
        )
    } else {
      pdf_name <- paste0("p_sling_umap_", r, "_Neu_mon_focused.pdf")
      color_var <- "state"
      cluster_pal <- c(
        "Other" = "grey85",
        "Other terminal" = "grey40",
        "Unddiff" = "#E77325",
        "Unddiff_2" = "grey70",
        "Unddiff_4" = "grey55",
        "Unddiff_6" = "grey40",   
        "Branch cluster" = "#D55E00",
        "Neu terminal" = "#8E77E6",
        "Mono terminal" = "#1F78B4"
      )
    }
    plot_df$leiden <- factor(plot_df$leiden, levels = cluster_levels)

    plot_df$state <- case_when(
      plot_df$leiden == branch_node ~ "Branch cluster",
      plot_df$leiden %in% neu_terms ~ "Neu terminal",
      plot_df$leiden %in% mono_terms ~ "Mono terminal",
      plot_df$leiden %in% other_terms ~ "Other terminal",
      plot_df$celltype == "Unddiff" ~ "Unddiff",
      TRUE ~ "Other"
    )
    
    ## -------------------------------------------------------
    ## 5. cluster centroids and lineage paths on UMAP
    ## -------------------------------------------------------
    cent <- plot_df %>%
      group_by(leiden) %>%
      summarise(
        UMAP_1 = median(UMAP_1, na.rm = TRUE),
        UMAP_2 = median(UMAP_2, na.rm = TRUE),
        .groups = "drop"
      )
    
    path_df <- bind_rows(lapply(seq_len(nrow(lineage_tbl)), function(i) {
      lin_i <- lineage_tbl$lineage_id[i]
      data.frame(
        lineage_id = lin_i,
        lineage_label = lineage_tbl$lineage_label[i],
        order = seq_along(lineages[[lin_i]]),
        leiden = as.character(lineages[[lin_i]])
      )
    })) %>%
      left_join(cent, by = "leiden") %>%
      arrange(lineage_id, order)
    
    label_df <- cent %>%
      filter(leiden %in% unique(path_df$leiden)) %>%
      mutate(label = leiden)
    
    ## -------------------------------------------------------
    ## 6. plot
    ## -------------------------------------------------------
    p_sling_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2)) +
      
      ggrastr::rasterise(
        geom_point(aes(color = .data[[color_var]]), size = 0.22, alpha = 0.75),
        dpi = 300
      ) +
      
      geom_path(
        data = path_df,
        aes(
          x = UMAP_1,
          y = UMAP_2,
          group = lineage_id
        ),
        inherit.aes = FALSE,
        linewidth = 1.2,
        color = "black",
        arrow = arrow(length = unit(0.18, "inches"), type = "closed")
      ) +
      
      geom_point(
        data = cent %>% filter(leiden %in% unique(path_df$leiden)),
        aes(x = UMAP_1, y = UMAP_2),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        color = "black",
        size = 2.5,
        stroke = 0.7
      ) +
      
      geom_label_repel(
        data = label_df,
        aes(x = UMAP_1, y = UMAP_2, label = label),
        inherit.aes = FALSE,
        size = 4.5,
        fontface = "bold",
        fill = "white",
        label.size = 0,
        box.padding = 0.6,
        point.padding = 0.5,
        max.overlaps = Inf
      ) +
      
      scale_color_manual(
        values = cluster_pal,
        drop = FALSE,
        name = r
      ) +
      
      coord_equal() +
      theme_classic(base_size = 16) +
      theme(
        legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      ggtitle(
        paste0(
          "Slingshot trajectories: ",
          ncol(sce), " cells; ",
          r,
          "; branch = ", branch_node
        )
      )
    
    pdf(
      pdf_name,
      width = 10,
      height = 7
    )
    print(p_sling_umap)
    dev.off()
  }
  
} # end of for( r in c("Cell.type.clean", "leiden_r0_6", "leiden_r0_6_anno") )
  cat("step3 done\n")
} # end of step3
# scp -p -r  xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/*.rds F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/.

# Plotting: leiden_r0_8
  # lineage terminal                                       path    n   frac_Und
# 1       1        8 4 -> 2 -> 1 -> 3 -> 7 -> 12 -> 9 -> 6 -> 8 5020 0.06553785
# 2       2       13     4 -> 2 -> 1 -> 3 -> 7 -> 14 -> 5 -> 13 3638 0.09400770
# 3       3       18             4 -> 2 -> 11 -> 16 -> 10 -> 18  525 0.29714286
# 4       4       15                          4 -> 2 -> 1 -> 15 2623 0.99008769
# 5       5       17                         4 -> 2 -> 11 -> 17 1717 0.17880023
# 6       6       19                                    4 -> 19  137 0.16788321
  # frac_Neutrophil frac_Monocyte      frac_D2   frac_D4   frac_D6 major_celltype
# 1    0.9322709163  0.0021912351 0.0005976096 0.1400398 0.8593625           Neut
# 2    0.0002748763  0.9051676745 0.1069268829 0.2399670 0.6531061           Mono
# 3    0.0095238095  0.0000000000 0.1009523810 0.3123810 0.5866667            Eos
# 4    0.0091498284  0.0003812429 0.3804803660 0.2981319 0.3213877        Unddiff
# 5    0.0058241118  0.0000000000 0.0407687828 0.2382062 0.7210250            Meg
# 6    0.0000000000  0.0000000000 0.4890510949 0.1824818 0.3284672        Ccr7_DC
  # major_time purity_celltype purity_time      terminal_class lineage_id
# 1          6       0.9322709   0.8593625 Neutrophil_terminal          1
# 2          6       0.9051677   0.6531061   Monocyte_terminal          2
# 3          6       0.6000000   0.5866667      Other_terminal          3
# 4          2       0.9900877   0.3804804      Other_terminal          4
# 5          6       0.8089691   0.7210250      Other_terminal          5
# 6          2       0.8321168   0.4890511      Other_terminal          6
        # lineage_label
# 1 Neutrophil_terminal
# 2   Monocyte_terminal
# 3           Other_Eos
# 4       Other_Unddiff
# 5           Other_Meg
# 6       Other_Ccr7_DC

# Plotting: leiden_r0_6
#   lineage terminal                              path    n frac_Und
# 1       1       10 3 -> 1 -> 5 -> 6 -> 12 -> 8 -> 10 4119        0
# 2       2       15  3 -> 1 -> 5 -> 6 -> 4 -> 2 -> 15  930        0
# 3       3       17 3 -> 1 -> 5 -> 6 -> 12 -> 8 -> 17  245        0
# 4       4        9                 3 -> 1 -> 13 -> 9 4687        0
# 5       5       14                3 -> 1 -> 13 -> 14 1857        0
# 6       6       11                3 -> 1 -> 13 -> 11 3947        0
# 7       7       18                 3 -> 1 -> 7 -> 18  137        0
# 8       8       16                3 -> 1 -> 13 -> 16  526        0
#   frac_Neutrophil frac_Monocyte     frac_D2   frac_D4   frac_D6 major_celltype
# 1               0             0 0.107793154 0.2449624 0.6472445           Mono
# 2               0             0 0.000000000 0.1139785 0.8860215           Neut
# 3               0             0 0.004081633 0.3836735 0.6122449        Unddiff
# 4               0             0 0.355664604 0.1802859 0.4640495        Unddiff
# 5               0             0 0.051157781 0.1766290 0.7722132           Mast
# 6               0             0 0.013934634 0.2703319 0.7157335           Baso
# 7               0             0 0.489051095 0.1824818 0.3284672        Ccr7_DC
# 8               0             0 0.102661597 0.3136882 0.5836502            Eos
#   major_time purity_celltype purity_time terminal_class lineage_id
# 1          6       0.9065307   0.6472445 Other_terminal          1
# 2          6       0.8107527   0.8860215 Other_terminal          2
# 3          6       0.9510204   0.6122449 Other_terminal          3
# 4          6       0.5824621   0.4640495 Other_terminal          4
# 5          6       0.8341411   0.7722132 Other_terminal          5
# 6          6       0.9569293   0.7157335 Other_terminal          6
# 7          2       0.8321168   0.4890511 Other_terminal          7
# 8          6       0.5988593   0.5836502 Other_terminal          8
#   lineage_label
# 1    Other_Mono
# 2    Other_Neut
# 3 Other_Unddiff
# 4 Other_Unddiff
# 5    Other_Mast
# 6    Other_Baso
# 7 Other_Ccr7_DC
# 8     Other_Eos
# Branch node: 6
# Bifurcation detected: TRUE
# Plotting: Cell.type.clean
#   lineage  terminal                           path     n frac_Und
# 1       1      Mast Unddiff -> Eos -> Baso -> Mast  1551        0
# 2       2 Erythroid    Unddiff -> Meg -> Erythroid   490        0
# 3       3   Ccr7_DC Unddiff -> Lymphoid -> Ccr7_DC   119        0
# 4       4       pDC     Unddiff -> Lymphoid -> pDC    37        0
# 5       5      Mono                Unddiff -> Mono  6913        0
# 6       6      Neut                Unddiff -> Neut 22233        0
#   frac_Neutrophil frac_Monocyte    frac_D2   frac_D4   frac_D6 major_celltype
# 1               0             0 0.02450032 0.1753707 0.8001289           Mast
# 2               0             0 0.17959184 0.1612245 0.6591837      Erythroid
# 3               0             0 0.45378151 0.1848739 0.3613445        Ccr7_DC
# 4               0             0 0.43243243 0.2432432 0.3243243            pDC
# 5               0             0 0.12541588 0.3409518 0.5336323           Mono
# 6               0             0 0.04803670 0.3022984 0.6496649           Neut
#   major_time purity_celltype purity_time terminal_class lineage_id
# 1          6               1   0.8001289 Other_terminal          1
# 2          6               1   0.6591837 Other_terminal          2
# 3          2               1   0.4537815 Other_terminal          3
# 4          2               1   0.4324324 Other_terminal          4
# 5          6               1   0.5336323 Other_terminal          5
# 6          6               1   0.6496649 Other_terminal          6
#     lineage_label
# 1      Other_Mast
# 2 Other_Erythroid
# 3   Other_Ccr7_DC
# 4       Other_pDC
# 5      Other_Mono
# 6      Other_Neut
# Branch node: Unddiff
# Bifurcation detected: TRUE
# Plotting: leiden_r0_6_anno
#   lineage        terminal
# 1       1          Mono_6
# 2       2          Neut_6
# 3       3 Mixed_Unddiff_6
# 4       4          Mast_6
# 5       5          Baso_6
# 6       6       Ccr7_DC_2
# 7       7     Mixed_Eos_6
#                                                                            path
# 1 Unddiff_2 -> Unddiff_4 -> Unddiff_6 -> Mixed_Neut_4 -> Mixed_Mono_4 -> Mono_6
# 2                                 Unddiff_2 -> Unddiff_4 -> Unddiff_6 -> Neut_6
# 3                                Unddiff_2 -> Unddiff/Baso_2 -> Mixed_Unddiff_6
# 4                                         Unddiff_2 -> Unddiff/Baso_2 -> Mast_6
# 5                                         Unddiff_2 -> Unddiff/Baso_2 -> Baso_6
# 6                                       Unddiff_2 -> Unddiff/Lym_2 -> Ccr7_DC_2
# 7                                    Unddiff_2 -> Unddiff/Baso_2 -> Mixed_Eos_6
#       n frac_Und frac_Neutrophil frac_Monocyte    frac_D2   frac_D4   frac_D6
# 1  4119        0               0             0 0.10779315 0.2449624 0.6472445
# 2 20066        0               0             0 0.04704475 0.2548091 0.6981461
# 3  4687        0               0             0 0.35566460 0.1802859 0.4640495
# 4  1857        0               0             0 0.05115778 0.1766290 0.7722132
# 5  3947        0               0             0 0.01393463 0.2703319 0.7157335
# 6   137        0               0             0 0.48905109 0.1824818 0.3284672
# 7   526        0               0             0 0.10266160 0.3136882 0.5836502
#   major_celltype major_time purity_celltype purity_time terminal_class
# 1           Mono          6       0.9065307   0.6472445 Other_terminal
# 2           Neut          6       0.9532543   0.6981461 Other_terminal
# 3        Unddiff          6       0.5824621   0.4640495 Other_terminal
# 4           Mast          6       0.8341411   0.7722132 Other_terminal
# 5           Baso          6       0.9569293   0.7157335 Other_terminal
# 6        Ccr7_DC          2       0.8321168   0.4890511 Other_terminal
# 7            Eos          6       0.5988593   0.5836502 Other_terminal
#   lineage_id lineage_label
# 1          1    Other_Mono
# 2          2    Other_Neut
# 3          3 Other_Unddiff
# 4          4    Other_Mast
# 5          5    Other_Baso
# 6          6 Other_Ccr7_DC
# 7          7     Other_Eos
# Branch node: Unddiff_6
# Bifurcation detected: TRUE

#############################################################
# Option 2: MST trajectory using PCA space
## running locally for convenience
#############################################################
step4 = TRUE
if(step4){
    setwd('F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/')
    # # setwd('/project/imoskowitz/xyang2/TIPS/GSE140802_lineage_tracking/results/')

    pkgs <- c(
    "slingshot", "TrajectoryUtils", "MatrixGenerics",
    "DelayedMatrixStats", "sparseMatrixStats", "matrixStats"
    )
    sapply(pkgs, packageVersion) 

    

    library(Seurat)
    library(igraph)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(ggrastr)
    library(grid)
    library(igraph)

    larry_nm <- readRDS( "Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds")

    npcs <- 30

    X <- SeuratObject::Embeddings(larry_nm, "pca")[, 1:npcs]

    for( r in c("leiden_r0_8")){ #c("Cell.type.clean", "leiden_r0_6_anno", "leiden_r0_6") ){
        cl <- as.character(larry_nm@meta.data[[r]])

        cent <- aggregate(as.data.frame(X), by = list(cluster = cl), FUN = median)
        rownames(cent) <- cent$cluster

        g <- graph_from_adjacency_matrix(as.matrix(dist(cent[, -1])), mode = "undirected", weighted = TRUE, diag = FALSE)
        mst_g <- mst(g, weights = E(g)$weight)

        if(r != "leiden_r0_6_anno")  set.seed(50)  else set.seed(100)
        lay <- layout_with_fr(mst_g) 
        p_mst_pca <- plot(mst_g, layout = lay,
            vertex.label = V(mst_g)$name, vertex.size = 25, edge.width = 2,
            main = paste0("Trajectory: ", nrow(larry_nm), ' cells; ', r, '; PCA; MST'))

        pdf(file = paste0('MST_trajectory_', nrow(larry_nm), '_', r, '_PCA.pdf'))
        plot(mst_g, layout = lay,
            vertex.label = V(mst_g)$name, vertex.size = 25, edge.width = 2,
            main = paste0("Trajectory: ", nrow(larry_nm), ' cells; ', r, '; PCA; MST'))
        dev.off()
    }

  cat("step4 done\n")
}

    ############################################    
    ## finally, slingshot running on three major clsuters: undiff. Mono, Neut (NOT used!!) ######
    ############################################
    # https://htmlpreview.github.io/?https://github.com/xyang2uchicago/BioTIP/blob/master/tutorial/Gastrulation.html

    # NOTE: MST trajectory using PCA space is more make sense !!!!!!!!
step5 = FALSE
if(step5){
  library(scuttle)
  outdir = 'F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/three_cluster_only/'
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


    input =  'logcounts'
    if (input == 'counts') {
        by.cluster <- scuttle::aggregateAcrossCells(sce, ids=colData(sce)$Cell.type.clean, use.assay.type='counts') # use.assay.type='counts' by default, here after SCT-normalization, scale.data is used !!!
    }
    if (input == 'logcounts') {
        by.cluster <- scuttle::aggregateAcrossCells(sce, ids=colData(sce)$Cell.type.clean, use.assay.type='logcounts') # use.assay.type='counts' by default, here after SCT-normalization, scale.data is used !!!
    }
    centroids.cluster <- reducedDim(by.cluster, "PCA")
    dmat.cluster <- as.matrix(dist(centroids.cluster))
    set.seed(1000)
    trajectory.cluster <- graph.adjacency(dmat.cluster, mode = "undirected", weighted = TRUE)
    pdf(file=paste0(outdir, 'MST_trajectory_aggregatedBy',input,'_PCA.pdf'))
      plot(minimum.spanning.tree(trajectory.cluster))
      mtext(paste0('MST_trajectory_aggregatedBy ',input,'+PCA'), side = 3, line = 2, font = 2, cex = 1.5) 
    dev.off()


    input =  'counts'

    if (input == 'counts') {
        by.cluster <- aggregateAcrossCells(sce, ids=colData(sce)$Cell.type.clean, use.assay.type='counts') # use.assay.type='counts' by default, here after SCT-normalization, scale.data is used !!!
    }
    if (input == 'logcounts') {
        by.cluster <- aggregateAcrossCells(sce, ids=colData(sce)$Cell.type.clean, use.assay.type='logcounts') # use.assay.type='counts' by default, here after SCT-normalization, scale.data is used !!!
    }
    centroids.cluster <- reducedDim(by.cluster, "SPRING")
    dmat.cluster <- as.matrix(dist(centroids.cluster))
    set.seed(1000)
    trajectory.cluster <- graph.adjacency(dmat.cluster, mode = "undirected", weighted = TRUE)
    pdf(file=paste0(outdir, 'MST_trajectory_aggregatedBy',input,'_spring.pdf'))
      plot(minimum.spanning.tree(trajectory.cluster))
      mtext(paste0('MST_trajectory_aggregatedBy ',input,'+SPRING'), side = 3, line = 2, font = 2, cex = 1.5) 
    dev.off()

    
    library(SingleCellExperiment)
    library(slingshot)

    ## Keep only the relevant NM branch cells
    sce_nm <- sce[, sce$NM_trajectory &
                    sce$Cell.type.clean %in% c("Unddiff", "Neut", "Mono")]

    tp <- as.character(sce_nm$Time.point)
    ct <- as.character(sce_nm$Cell.type.clean)

    sce_nm$nm_branch_state <- NA_character_

    sce_nm$nm_branch_state[
    tp == "2" & ct == "Unddiff"
    ] <- "D2_Und_parent"

    sce_nm$nm_branch_state[
    tp %in% c("4", "6") & ct == "Neut"
    ] <- "D4D6_Neutrophil"

    sce_nm$nm_branch_state[
    tp %in% c("4", "6") & ct == "Mono"
    ] <- "D4D6_Monocyte"

    sce_nm <- sce_nm[, !is.na(sce_nm$nm_branch_state)]

    sce_nm$nm_branch_state <- factor(
    sce_nm$nm_branch_state,
    levels = c(
        "D2_Und_parent",
        "D4D6_Neutrophil",
        "D4D6_Monocyte"
    )
    )

    sce_nm
    # class: SingleCellExperiment 
    # dim: 25289 52129 
    # metadata(0):
    # assays(2): counts logcounts
    # rownames(25289): 0610006L08Rik 0610007P14Rik ... mt-Nd5 mt-Nd6
    # rowData names(0):
    # colnames(52129): cell_2 cell_3 ... cell_130877 cell_130881
    # colData names(18): orig.ident nCount_RNA ... Cell.type.clean
    #   nm_branch_state
    # reducedDimNames(3): PCA UMAP SPRING
    # mainExpName: NULL
    # altExpNames(0):
    colnames(colData(sce_nm))
    #  [1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"        
    #  [4] "Library"              "Cell.barcode"         "Time.point"          
    #  [7] "Starting.population"  "Cell.type.annotation" "Well"                
    # [10] "SPRING.x"             "SPRING.y"             "cell_index0"         
    # [13] "cell_id"              "NM_trajectory"        "Fig6B_label"         
    # [16] "CellCycleSig1"        "Cell.type.clean"      "nm_branch_state"     
    table(colData(sce_nm)$nm_branch_state)
    # D2_Und_parent D4D6_Neutrophil   D4D6_Monocyte 
    #         24918           21165            6046 

    ## Use PCA for inference
    sce_nm <- slingshot(
    sce_nm,
    clusterLabels = "nm_branch_state",
    reducedDim = "PCA",
    start.clus = "D2_Und_parent",
    end.clus = c("D4D6_Neutrophil", "D4D6_Monocyte")
    )

    colnames(label_df)
    #[1] "nm_branch_state" "DR1"             "DR2"             "n"              
    
    label_from <- label_df %>%
    dplyr::rename(
        from_DR1 = DR1,
        from_DR2 = DR2
    )

    label_to <- label_df %>%
    dplyr::rename(
        to_DR1 = DR1,
        to_DR2 = DR2
    )

    edge_df <- data.frame(
    from = "D2_Und_parent",
    to = c("D4D6_Neutrophil", "D4D6_Monocyte")
    ) %>%
    dplyr::left_join(
        label_from,
        by = c("from" = "nm_branch_state")
    ) %>%
    dplyr::left_join(
        label_to,
        by = c("to" = "nm_branch_state")
    )

    edge_df
    #            from              to  from_DR1  from_DR2   n.x   to_DR1   to_DR2
    # 1 D2_Und_parent D4D6_Neutrophil -4.416771 -3.481628 24918 7.649006 1.435865
    # 2 D2_Und_parent   D4D6_Monocyte -4.416771 -3.481628 24918 2.202951 5.038877
    #     n.y
    # 1 21165
    # 2  6046

    p_branch <- ggplot(plot_df, aes(x = DR1, y = DR2)) +

    ggrastr::rasterise(
        geom_point(
        aes(color = nm_branch_state),
        size = 0.25,
        alpha = 0.75
        ),
        dpi = 300
    ) +

    geom_curve(
        data = edge_df,
        aes(
        x = from_DR1,
        y = from_DR2,
        xend = to_DR1,
        yend = to_DR2
        ),
        inherit.aes = FALSE,
        curvature = 0.20,
        linewidth = 0.9,
        color = "grey25",
        arrow = arrow(length = unit(0.20, "inches"), type = "closed")
    ) +

    geom_label_repel(
        data = label_df,
        aes(
        x = DR1,
        y = DR2,
        label = nm_branch_state,
        color = nm_branch_state
        ),
        inherit.aes = FALSE,
        size = 5,
        fontface = "bold",
        fill = "white",
        alpha = 0.92,
        label.size = 0,
        box.padding = 0.6,
        point.padding = 0.4,
        max.overlaps = Inf
    ) +

    scale_color_manual(values = state_pal, drop = FALSE) +
    coord_equal() +
    theme_classic(base_size = 16) +
    theme(
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    ) +
    ggtitle(paste0(dimred_use, ": time-guided NM bifurcation, ", nrow(sce_nm), " cells"))
    
    pdf(file=paste0(outdir, 'timeGuided_NMT_NMbranch_', nrow(sce_nm),'cells.pdf'))
    print(p_branch)
    dev.off()
}
