setwd("f:/projects/TIPS/results/GSE140802_lineage_tracking")

library(Seurat)
packageVersion("Seurat") # ‘5.4.0.9022’
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(patchwork)
library(curl)

library(dplyr)
library(ggrepel)
library(ggrastr)
library(scales)

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

outdir = 'inVitro_NMtrajectory/'

# refer to 2_invitro_NMtrajectory_creatSeurat_v2.R
larry_nm <- readRDS(file.path(outdir, "Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds"))
larry_nm 
# An object of class Seurat 
# 25289 features across 92527 samples within 1 assay 
# Active assay: RNA (25289 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  3 dimensional reductions calculated: spring, pca, umap


# Check the column names of the metadata in the Seurat object
colnames(larry_nm@meta.data)
#  [1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"        
#  [4] "Library"              "Cell.barcode"         "Time.point"          
#  [7] "Starting.population"  "Cell.type.annotation" "Well"                
# [10] "SPRING.x"             "SPRING.y"             "cell_index0"         
# [13] "cell_id"              "NM_trajectory"        "Fig6B_label"         
# [16] "CellCycleSig1"        "Cell.type.clean"     
 
table(larry_nm@meta.data$Cell.type.clean)
#        Baso    Ccr7_DC        Eos  Erythroid   Lymphoid       Mast        Meg 
#       4255        119        355        490        189       1551       1468 
#   Monocyte Neutrophil        pDC        Und 
#       6913      22233         37      54917 

library(Seurat)
library(ggplot2)
library(patchwork)
library(ggrastr)

DefaultAssay(larry_nm) <- "RNA"

## Choose embedding
reduction_use <- "spring"
# reduction_use <- "umap"

stopifnot(reduction_use %in% Reductions(larry_nm))

plot_titles <- c(
  "Cell.type.clean" = "cell type",
  "Cebpa" = "Cebpa (HSP)",
  "Runx1" = "Runx1 (HSP)",
  "Gfi1"  = "Gfi1 (Neu)",
  "Cebpe" = "Cebpe (Neu)",
  "Irf8"  = "Irf8 (Mono)",
  "Klf4"  = "Klf4 (Mono)"
)

celltype_pal <- c(
  "Baso"       = "#FB8072",
  "Ccr7_DC"    = "#E69F00",
  "Eos"        = "#B3A400",
  "Erythroid"  = "#66A61E",
  "Lymphoid"   = "#00BA38",
  "Mast"       = "#00BFC4",
  "Meg"        = "#00AEEF",
  "Monocyte"   = "#0099E6",
  "Neutrophil" = "#B58CFF",
  "pDC"        = "#FF61CC",
  "Und"        = "#F05AA6"
)

## Make sure cell-type order matches the palette
larry_nm$Cell.type.clean <- factor(
  larry_nm$Cell.type.clean,
  levels = names(celltype_pal)
)

theme_panel <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

make_seurat_plot <- function(obj, var, reduction_use, plot_title) {
  
  is_gene <- var %in% rownames(obj)
  is_meta <- var %in% colnames(obj@meta.data)
  
  if (!is_gene && !is_meta) {
    stop("Variable not found as gene or metadata column: ", var)
  }
  
  ## =====================================================
  ## MODIFICATION 1:
  ## Special handling for Cell.type.clean.
  ## This puts cell-type annotation directly on the map.
  ## =====================================================
  if (var == "Cell.type.clean") {
    
    plot_obj <- DimPlot(
      object = obj,
      reduction = reduction_use,
      group.by = var,
      pt.size = 0.3,
      raster = FALSE
    ) +
      scale_color_manual(values = celltype_pal, drop = FALSE) +
      ggtitle(plot_title) +
      theme_panel
    
    ## =====================================================
    ## MODIFICATION 2:
    ## Add on-site labels instead of relying on legend.
    ## =====================================================
    plot_obj <- LabelClusters(
      plot = plot_obj,
      id = var,
      repel = TRUE,
      size = 5,
      fontface = "bold"
    ) +
      NoLegend()
    
    ## =====================================================
    ## MODIFICATION 3:
    ## Rasterize the plot to reduce file size.
    ## =====================================================
    plot_obj <- ggrastr::rasterise(
      plot_obj,
      dpi = 300
    )  # IMPORTANT to reduce the figure size for large dataset
    
    return(plot_obj)
  }
  
  ## =====================================================
  ## Gene-expression panels
  ## =====================================================
  if (is_gene) {
    
    plot_obj <- FeaturePlot(
      object = obj,
      features = var,
      reduction = reduction_use,
      order = TRUE,
      pt.size = 0.3,
      raster = FALSE,
      cols = c("lightgrey", "darkred"),
      min.cutoff = "q05",
      max.cutoff = "q95"
    ) +
      ggtitle(plot_title) +
      theme_panel
    
    plot_obj <- ggrastr::rasterise(
      plot_obj,
      dpi = 300
    )  # IMPORTANT to reduce the figure size for large dataset
    
    return(plot_obj)
  }
  
  ## =====================================================
  ## Optional branch: numeric metadata, e.g. entropy/module score
  ## =====================================================
  meta_vec <- obj@meta.data[[var]]
  
  if (is.numeric(meta_vec)) {
    
    plot_obj <- FeaturePlot(
      object = obj,
      features = var,
      reduction = reduction_use,
      order = TRUE,
      pt.size = 0.3,
      raster = FALSE,
      cols = c("lightgrey", "darkred"),
      min.cutoff = "q05",
      max.cutoff = "q95"
    ) +
      ggtitle(plot_title) +
      theme_panel
    
    plot_obj <- ggrastr::rasterise(
      plot_obj,
      dpi = 300
    )  # IMPORTANT to reduce the figure size for large dataset
    
    return(plot_obj)
  }
  
  ## =====================================================
  ## Optional branch: other discrete metadata
  ## =====================================================
  plot_obj <- DimPlot(
    object = obj,
    reduction = reduction_use,
    group.by = var,
    pt.size = 0.3,
    raster = FALSE
  ) +
    ggtitle(plot_title) +
    theme_panel +
    NoLegend()
  
  plot_obj <- ggrastr::rasterise(
    plot_obj,
    dpi = 300
  )  # IMPORTANT to reduce the figure size for large dataset
  
  return(plot_obj)
}


for(reduction_use in c("spring", "umap")) {
    ## ---------------------------------------------------------
    ## Build plots
    ## ---------------------------------------------------------

    celltype_var <- "Cell.type.clean"
    marker_genes <- setdiff(names(plot_titles), celltype_var)

    marker_genes_present <- marker_genes[marker_genes %in% rownames(larry_nm)]
    marker_genes_missing <- setdiff(marker_genes, marker_genes_present)

    if (length(marker_genes_missing) > 0) {
    warning(
        "These genes were not found in the Seurat object and will be skipped: ",
        paste(marker_genes_missing, collapse = ", ")
    )
    }

    p_celltype <- make_seurat_plot(
    obj = larry_nm,
    var = celltype_var,
    reduction_use = reduction_use,
    plot_title = plot_titles[[celltype_var]]
    )

    marker_plots <- lapply(marker_genes_present, function(gene) {
    make_seurat_plot(
        obj = larry_nm,
        var = gene,
        reduction_use = reduction_use,
        plot_title = plot_titles[[gene]]
    )
    })

    names(marker_plots) <- marker_genes_present

    ## ---------------------------------------------------------
    ## 1 + 3x2 arrangement
    ## ---------------------------------------------------------

    p_markers_3x2 <- wrap_plots(marker_plots, ncol = 2)

    p_all <- p_celltype + p_markers_3x2 # +
    #   plot_layout(heights = c(1.05, 2.0)) +
    #   plot_annotation(tag_levels = "A") &
    #   theme(
    #     plot.tag = element_text(size = 18, face = "bold")
    #   )

    #print(p_all)

    ggsave(
    filename = paste0(outdir, toupper(reduction_use), "_NMtrajectory_celltype_plus_markers.pdf"),
    plot = p_all,
    width = 12,
    height = 6,
    units = "in",
    device = cairo_pdf,
    dpi = 300,
    limitsize = FALSE
    )
}


# # for a Seurat object, use:
# FeaturePlot(larry_nm, features = gene, reduction = reduction_use)
# # for genes
# DimPlot(larry_nm, group.by = "Cell.type.clean", reduction = reduction_use)
## marker genes
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggrastr)
library(ggrepel)
library(dplyr)
library(patchwork)
 
genes_show <- c("Cebpa", "Runx1", "Gfi1", "Cebpe", "Irf8", "Klf4")

for(reduction_use in c("spring", "umap")) {
  ## large cell-type panel
  p_celltype <- DimPlot(
    larry_nm,
    reduction = reduction_use,
    group.by = "Cell.type.clean",
    pt.size = 0.25,
    label = FALSE
  ) +
    scale_color_manual(values = celltype_pal) +
    ggtitle("cell type") +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  ## rasterize cell points, before adding labels
  p_celltype <- ggrastr::rasterise(p_celltype, layers = "Point", dpi = 300)

  ## add cell-type labels manually at median position
  emb <- Embeddings(larry_nm, reduction_use) %>% as.data.frame()
  colnames(emb)[1:2] <- c("Dim1", "Dim2")
  emb$Cell.type.clean <- larry_nm$Cell.type.clean

  label_df <- emb %>%
    group_by(Cell.type.clean) %>%
    summarise(
      Dim1 = median(Dim1, na.rm = TRUE),
      Dim2 = median(Dim2, na.rm = TRUE),
      .groups = "drop"
    )

  p_celltype <- p_celltype +
    geom_text_repel(
      data = label_df,
      aes(x = Dim1, y = Dim2, label = Cell.type.clean),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold",
      color = "black",
      max.overlaps = Inf
    )

  ## gene panels
  gene_plots <- lapply(genes_show, function(g) {
    p <- FeaturePlot(
      larry_nm,
      features = g,
      reduction = reduction_use,
      pt.size = 0.25,
      order = TRUE
    ) +
      scale_colour_gradient(low = "lightgrey", high = "darkred") +
      ggtitle(g) +
      theme_classic(base_size = 13) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )

    ggrastr::rasterise(p, layers = "Point", dpi = 300)
  })

  p_genes <- patchwork::wrap_plots(gene_plots, ncol = 2)

  p_final <- p_celltype + p_genes +
    patchwork::plot_layout(widths = c(1.6, 1))

  pdf(
    paste0(reduction_use, "_NMtrajectoryCells_lineageMarkers.pdf"),
    width = 14,
    height = 7
  )
  print(p_final)
  dev.off()
}

## if sce object, use:

# Reductions(larry_nm) # [1] "spring" "pca" "umap" 
# for(reduction_use in c("spring", "umap")) { 
#   stopifnot(reduction_use %in% Reductions(larry_nm)) 
#   # If "Cell.type.clean" is included as one of the genes, it is discrete and cannot use scale_colour_gradient. 
#   # We need to handle louvain separately, with a suitable discrete color scale. 
#   plot_list <- lapply(genes_present, function(gene) 
#   { 
#     plotUMAP(larry_nm, colour_by = gene, point_size = 0.3, by_exprs_values = "X") + 
#     scale_colour_gradient(low = "lightgrey", high = "darkred") + 
#     ggtitle(plot_titles[[gene]]) + theme_classic() 
#     }
#   ) 
    
#   plots <- mapply(function(gene, plot_obj) { 
#     if (gene == "Cell.type.clean") { 
#         p <- plotUMAP(larry_nm, colour_by = "Cell.type.clean", point_size = 0.3, by_exprs_values = "X") + 
#         scale_colour_manual(values = celltype_pal) + 
#         ggtitle(plot_titles[["Cell.type"]]) + 
#         theme_classic() ggrastr::rasterise(p, dpi = 300) 
#     } else { 
#           ggrastr::rasterise(plot_obj, dpi = 300) #   IMPORTSNT to reduce the figure size for large dataset 
#           } 
#    }, genes, plot_list, SIMPLIFY = FALSE) 
          
#   g1 <- do.call(gridExtra::arrangeGrob, c(plots, ncol = 3)) 
          
#   pdf(file= paste0(reduction_use,'_NMtrajectoryCells_lineageMarkers.pdf')) 
#   grid.draw(g1) 
#   dev.off() 
#   }