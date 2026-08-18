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

data_path <- "F:/projects/TIPS/results/GSE140802_lineage_tracking/"
outdir = 'inVitro/'

larry <- readRDS(file.path(data_path, outdir, "reanalyzed_Seurat5_noCellCycle_PCA_UMAP.rds"))
larry
# An object of class Seurat 
# 25289 features across 130887 samples within 1 assay 
# Active assay: RNA (25289 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  3 dimensional reductions calculated: spring, pca, umap


genes_to_plot <- c(
  "Cebpa",
  "Runx1",
  "Gfi1",
  "Cebpe", #Elane,
  "Irf8",
  "Klf4"   #, "Csf1r"
)

# Check the column names of the metadata in the Seurat object
colnames(larry@meta.data)
#  [1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"        
#  [4] "Library"              "Cell.barcode"         "Time.point"          
#  [7] "Starting.population"  "Cell.type.annotation" "Well"                
# [10] "SPRING.x"             "SPRING.y"             "cell_index0"         
# [13] "cell_id"              "NM_trajectory"        "Fig6B_label"         
# [16] "CellCycleSig1"  

## Optional: rename for cleaner plot labels
table(larry@meta.data$Cell.type.annotation)
    #         Baso          Ccr7_DC              Eos        Erythroid 
    #        10155              241              555              686 
    #     Lymphoid             Mast              Meg         Monocyte 
    #          857             2584             2096            19182 
    #   Neutrophil              pDC Undifferentiated 
    #        22237              101            72193 
larry@meta.data$Cell.type.clean <- as.character(larry@meta.data$Cell.type.annotation)
larry@meta.data$Cell.type.clean[
  larry@meta.data$Cell.type.clean == "Undifferentiated"
] <- "Und"
 
celltype_pal <- c(
  "Baso"      = "#FB8072",
  "Ccr7_DC"   = "#E69F00",
  "Eos"       = "#B3A400",
  "Erythroid" = "#66A61E",
  "Lymphoid"  = "#00BA38",
  "Mast"      = "#00BFC4",
  "Meg"       = "#00AEEF",
  "Monocyte"      = "#0099E6",
  "Neutrophil"       = "#B58CFF",
  "pDC"       = "#FF61CC",
  "Und"       = "#F05AA6"
)
larry@meta.data$Cell.type.clean <- factor(larry@meta.data$Cell.type.clean, levels = names(celltype_pal))

time_pal <- c(
  "2" = "yellow",
  "4" = "orange",
  "6" = "brown"
)

## Visualize the celltype_pal palette
barplot(
  rep(1, length(celltype_pal)),
  col = celltype_pal,
  names.arg = names(celltype_pal),
  main = "celltype_pal color palette",
  ylab = "",
  yaxt = "n"    
)

state_pal <- c(
  "Other" = "grey85",
  "Parent_D2_Undiff" = "#E77325",
  "Desc_Neu_D4D6" = "#8E77E6",
  "Desc_Mo_D4D6" = "#1F78B4"
)


## ############################################################
## A) Plot Fig. based on SPRING coordinates
## ############################################################

## -----------------------------
## 5. Plot Fig. 6B draft by cell type
## -----------------------------

emb <- as.data.frame(Embeddings(larry, reduction = "spring"))
emb$cell_id <- rownames(emb)
emb$NM_trajectory <- larry$NM_trajectory
emb$Cell.type.annotation <- larry$Cell.type.annotation
emb$Time.point <- larry$Time.point
emb$Starting.population <- larry$Starting.population
emb <- emb %>%
  mutate(
    Cell.type.annotation = as.character(Cell.type.annotation),
    Time.point = as.character(Time.point)
  )
## simplify labels if they are too long
emb$Cell.type.clean <- emb$Cell.type.annotation
emb$Cell.type.clean[emb$Cell.type.clean == "Neutrophil"] <- "Neu"
emb$Cell.type.clean[emb$Cell.type.clean == "Monocyte"] <- "Mono"
emb$Cell.type.clean[emb$Cell.type.clean == "Undifferentiated"] <- "Und"
emb$Cell.type.clean[emb$Cell.type.clean == "outside_NM_trajectory"] <- "Other"

cell_lvls <- sort(unique(emb$Cell.type.clean[emb$NM_trajectory]))
emb$Cell.type.clean <- factor(emb$Cell.type.clean, levels = cell_lvls)

cell_pal <- setNames(scales::hue_pal()(length(cell_lvls)), cell_lvls)


label_df <- emb %>%
  filter(NM_trajectory) %>%
  filter(!is.na(Cell.type.clean), Cell.type.clean != "") %>%
  group_by(Cell.type.clean) %>%
  summarise(
    SPRING_1 = median(SPRING_1, na.rm = TRUE),
    SPRING_2 = median(SPRING_2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 20)   ## increase to 50 if too many labels

label_df$Cell.type.clean <- factor(label_df$Cell.type.clean, levels = cell_lvls)

## =========================================================
## 0. Prepare data
## =========================================================

## make sure emb has cell_id
if (!"cell_id" %in% colnames(emb)) {
  emb$cell_id <- rownames(emb)
}

## use Seurat normalized data if available, otherwise counts
rna_layer <- if ("data" %in% Layers(larry[["RNA"]])) "data" else "counts"
expr_mat <- GetAssayData(larry, assay = "RNA", layer = rna_layer)

## keep common cells only
common_cells <- intersect(emb$cell_id, colnames(expr_mat))
emb_plot <- emb %>%
  filter(cell_id %in% common_cells) %>%
  mutate(
    NMtraj = factor(ifelse(NM_trajectory, "Y", "N"), levels = c("N", "Y"))
  )

expr_mat <- expr_mat[, common_cells, drop = FALSE]

## fixed axis ranges so all panels align visually
xlims <- range(emb_plot$SPRING_1, na.rm = TRUE)
ylims <- range(emb_plot$SPRING_2, na.rm = TRUE)

## marker genes
genes_to_plot <- c(
  "Cebpa",
  "Runx1",
  "Gfi1",
  "Cebpe",
  "Irf8",
  "Klf4"
)

genes_present <- genes_to_plot[genes_to_plot %in% rownames(expr_mat)]
genes_missing <- setdiff(genes_to_plot, genes_present)

if (length(genes_missing) > 0) {
  message("Missing genes: ", paste(genes_missing, collapse = ", "))
}

## =========================================================
## 1. Common themes
## =========================================================

theme_big <- theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

theme_small <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

## =========================================================
## 2. Panel 1: original cell-type-labeled plot
## =========================================================

p_celltype_onplot <- ggplot() +

  ## background
  ggrastr::rasterise(
    geom_point(
      data = emb_plot %>% filter(!NM_trajectory),
      aes(x = SPRING_1, y = SPRING_2),
      color = "grey85",
      size = 0.12,
      alpha = 0.35
    ),
    dpi = 300
  ) +

  ## highlighted NM trajectory cells
  ggrastr::rasterise(
    geom_point(
      data = emb_plot %>% filter(NM_trajectory), #!!!!!!!!!!!!!!
      aes(x = SPRING_1, y = SPRING_2, color = Cell.type.clean),
      size = 0.45,
      alpha = 0.90
    ),
    dpi = 300
  ) +

  ## on-site labels
  geom_label_repel(
    data = label_df,
    aes(
      x = SPRING_1,
      y = SPRING_2,
      label = Cell.type.clean,
      color = Cell.type.clean
    ),
    size = 6,
    fontface = "bold",
    fill = "white",
    alpha = 0.90,
    label.size = 0,
    box.padding = 0.65,
    point.padding = 0.45,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +

  scale_color_manual(values = celltype_pal, drop = FALSE) +
  coord_equal(xlim = xlims, ylim = ylims, expand = FALSE) +
  labs(
    title = "NM trajectory in vitro landscape",
    subtitle = "NM trajectory cells colored by annotated cell type",
    x = "SPRING 1",
    y = "SPRING 2"
  ) +
  theme_big

## =========================================================
## 3. Panel 2: NM trajectory yes/no
## =========================================================

p_nmtraj <- ggplot() +

  ## non-NM
  ggrastr::rasterise(
    geom_point(
      data = emb_plot %>% filter(NMtraj == "N"),
      aes(x = SPRING_1, y = SPRING_2),
      color = "grey85",
      size = 0.12,
      alpha = 0.40
    ),
    dpi = 300
  ) +

  ## NM trajectory
  ggrastr::rasterise(
    geom_point(
      data = emb_plot %>% filter(NMtraj == "Y"),
      aes(x = SPRING_1, y = SPRING_2),
      color = "#A6CEE3",   ## light blue
      size = 0.40,
      alpha = 0.95
    ),
    dpi = 300
  ) +

  coord_equal(xlim = xlims, ylim = ylims, expand = FALSE) +
  labs(
    title = "NM trajectory assignment",
    subtitle = "Y = light blue; N = light grey",
    x = "SPRING 1",
    y = "SPRING 2"
  ) +
  theme_big

## =========================================================
## 4. Marker-plot helper
## =========================================================

make_marker_plot <- function(gene, emb_df, expr_mat) {

  plot_df <- emb_df %>%
    select(cell_id, SPRING_1, SPRING_2, NM_trajectory)

  plot_df$expr <- as.numeric(expr_mat[gene, plot_df$cell_id])

  ## clip to q05-q95 for a FeaturePlot-like display
  q_low  <- quantile(plot_df$expr, 0.05, na.rm = TRUE)
  q_high <- quantile(plot_df$expr, 0.95, na.rm = TRUE)

  if (is.na(q_low))  q_low  <- min(plot_df$expr, na.rm = TRUE)
  if (is.na(q_high)) q_high <- max(plot_df$expr, na.rm = TRUE)

  if (q_high <= q_low) {
    q_low  <- min(plot_df$expr, na.rm = TRUE)
    q_high <- max(plot_df$expr, na.rm = TRUE)
  }

  plot_df$expr_clip <- pmin(pmax(plot_df$expr, q_low), q_high)

  ggplot() +

    ## background
    ggrastr::rasterise(
      geom_point(
        data = plot_df,
        aes(x = SPRING_1, y = SPRING_2),
        color = "grey90",
        size = 0.12,
        alpha = 0.35
      ),
      dpi = 300
    ) +

    ## expression overlay
    ggrastr::rasterise(
      geom_point(
        data = plot_df %>% filter(expr_clip > q_low),
        aes(x = SPRING_1, y = SPRING_2, color = expr_clip),
        size = 0.28,
        alpha = 0.95
      ),
      dpi = 300
    ) +

    scale_color_gradient(
      low = "lightgrey",
      high = "darkred"
    ) +

    coord_equal(xlim = xlims, ylim = ylims, expand = FALSE) +
    labs(title = gene) +
    theme_small
}

marker_plots <- lapply(genes_present, function(g) {
  make_marker_plot(gene = g, emb_df = emb_plot, expr_mat = expr_mat)
})

names(marker_plots) <- genes_present

## =========================================================
## 5. Combine into 1 + 1 + 6 layout
## =========================================================

top_row <- p_celltype_onplot | p_nmtraj
bottom_row <- wrap_plots(marker_plots, ncol = 3)

p_all <- (top_row / bottom_row) +
  plot_layout(heights = c(1.0, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold")
  )

## =========================================================
## 6. Save
## =========================================================

ggsave(
  filename = file.path(outdir, "SPRING_LARRY_NM_states_plus_markers.pdf"),
  plot = p_all,
  width = 16,
  height = 11,
  units = "in",
  device = cairo_pdf,
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(outdir, "SPRING_LARRY_NM_states_plus_markers.pdf"),
  plot = p_all,
  width = 16,
  height = 11,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

#print(p_all)

####################################
## B) UMAP 
####################################
larry_nm <- subset(larry, subset = NM_trajectory)
dim(larry_nm) #25289 92527


#library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggrastr)
library(patchwork)
library(scales)

## ---------------------------------------------------------
## 1. choose UMAP reduction
## ---------------------------------------------------------

umap_red <- if ("umap_nocc" %in% names(larry_nm@reductions)) {
  "umap_nocc"
} else if ("umap" %in% names(larry_nm@reductions)) {
  "umap"
} else {
  stop("No UMAP reduction found in larry_nm.")
}

## ---------------------------------------------------------
## 2. build plotting dataframe from UMAP + metadata
## ---------------------------------------------------------

emb_umap <- as.data.frame(Embeddings(larry_nm, reduction = umap_red))
colnames(emb_umap)[1:2] <- c("UMAP_1", "UMAP_2")
emb_umap$cell_id <- rownames(emb_umap)

meta_umap <- larry_nm@meta.data %>%
  as.data.frame() %>%
  mutate(cell_id = rownames(.))

plot_df <- left_join(emb_umap, meta_umap, by = "cell_id")

## clean labels
if (!"Cell.type.clean" %in% colnames(plot_df)) {
  plot_df$Cell.type.clean <- as.character(plot_df$Cell.type.annotation)
}
plot_df$Cell.type.clean <- as.character(plot_df$Cell.type.clean)
plot_df$Cell.type.clean[plot_df$Cell.type.clean == "Undifferentiated"] <- "Undiff"

plot_df$Time.point <- factor(as.character(plot_df$Time.point), levels = c("2", "4", "6"))

## if TIPS_state does not exist yet, define it
if (!"TIPS_state" %in% colnames(plot_df)) {
  plot_df$TIPS_state <- "Other"

  plot_df$TIPS_state[
    plot_df$Time.point == "2" &
      plot_df$Cell.type.annotation == "Undifferentiated"
  ] <- "Parent_D2_Undiff"

  plot_df$TIPS_state[
    plot_df$Time.point %in% c("4", "6") &
      plot_df$Cell.type.annotation == "Neutrophil"
  ] <- "Desc_Neu_D4D6"

  plot_df$TIPS_state[
    plot_df$Time.point %in% c("4", "6") &
      plot_df$Cell.type.annotation == "Monocyte"
  ] <- "Desc_Mo_D4D6"
}

plot_df$TIPS_state <- factor(
  plot_df$TIPS_state,
  levels = c("Other", "Parent_D2_Undiff", "Desc_Neu_D4D6", "Desc_Mo_D4D6")
)



## ---------------------------------------------------------
## 4. label positions
## ---------------------------------------------------------
dim(plot_df) 
label_cell_df <- plot_df %>%
  group_by(Cell.type.clean) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 20)

label_time_df <- plot_df %>%
  group_by(Time.point) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

label_state_df <- plot_df %>%
  group_by(TIPS_state) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

## ---------------------------------------------------------
## 5. common theme
## ---------------------------------------------------------

theme_top <- theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

## ---------------------------------------------------------
## 6. panel A: cell type with colored label boxes
## ---------------------------------------------------------

pA <- ggplot() +

  ggrastr::rasterise(
    geom_point(
      data = plot_df,
      aes(x = UMAP_1, y = UMAP_2, color = Cell.type.clean),
      size = 0.26,
      alpha = 0.92
    ),
    dpi = 300
  ) +

  geom_label_repel(
    data = label_cell_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = Cell.type.clean,
      color = Cell.type.clean
    ),
    size = 6,
    fontface = "bold",
    fill = "white",
    alpha = 0.92,
    label.size = 0,
    box.padding = 0.60,
    point.padding = 0.40,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +

  scale_color_manual(values = celltype_pal, drop = FALSE) +
  coord_equal() +
  labs(
    title = "LARRY UMAP",
    subtitle = "NM trajectory cells colored by cell type",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_top

## ---------------------------------------------------------
## 7. panel B: time point with colored label boxes
## ---------------------------------------------------------

pB <- ggplot() +

  ggrastr::rasterise(
    geom_point(
      data = plot_df,
      aes(x = UMAP_1, y = UMAP_2, color = Time.point),
      size = 0.26,
      alpha = 0.88
    ),
    dpi = 300
  ) +

  geom_label_repel(
    data = label_time_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = Time.point,
      color = Time.point
    ),
    size = 7,
    fontface = "bold",
    fill = "white",
    alpha = 0.92,
    label.size = 0,
    box.padding = 0.70,
    point.padding = 0.50,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +

  scale_color_manual(values = time_pal, drop = FALSE) +
  coord_equal() +
  labs(
    title = "LARRY UMAP",
    subtitle = "Colored by time point",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_top

## ---------------------------------------------------------
## 8. panel C: three-state definition
##     draw Other first in grey, then colored states on top
## ---------------------------------------------------------

plot_df_other <- plot_df %>% filter(TIPS_state == "Other")
plot_df_focus <- plot_df %>% filter(TIPS_state != "Other")

pC <- ggplot() +

  ## Other first
  ggrastr::rasterise(
    geom_point(
      data = plot_df_other,
      aes(x = UMAP_1, y = UMAP_2),
      color = "grey85",
      size = 0.24,
      alpha = 0.55
    ),
    dpi = 300
  ) +

  ## highlighted 3 states
  ggrastr::rasterise(
    geom_point(
      data = plot_df_focus,
      aes(x = UMAP_1, y = UMAP_2, color = TIPS_state),
      size = 0.28,
      alpha = 0.95
    ),
    dpi = 300
  ) +

  geom_label_repel(
    data = label_state_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = TIPS_state,
      color = TIPS_state
    ),
    size = 5.6,
    fontface = "bold",
    fill = "white",
    alpha = 0.92,
    label.size = 0,
    box.padding = 0.65,
    point.padding = 0.45,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +

  scale_color_manual(values = state_pal, drop = FALSE) +
  coord_equal() +
  labs(
    title = "Three-state definition",
    subtitle = "Parent and two descendants shown explicitly",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_top

## ---------------------------------------------------------
## 9. keep your current marker FeaturePlots as before
## ---------------------------------------------------------

genes_to_plot <- c(
  "Cebpa",
  "Runx1",
  "Gfi1",
  "Cebpe",
  "Irf8",
  "Klf4"
)

genes_present <- genes_to_plot[genes_to_plot %in% rownames(larry_nm)]

p_features <- FeaturePlot(
  object = larry_nm,
  features = genes_present,
  reduction = umap_red,
  ncol = 3,
  order = TRUE,
  raster = TRUE,
  cols = c("lightgrey", "darkred"),
  min.cutoff = "q05",
  max.cutoff = "q95"
) &
  theme_classic(base_size = 14) &
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

## ---------------------------------------------------------
## 10. combine and save
## ---------------------------------------------------------

p_top <- pA | pB | pC

p_all <- (p_top / p_features) +
  plot_layout(heights = c(1.0, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold")
  )

ggsave(
  filename = file.path(outdir, "UMAP_LARRY_states_plus_markers_labelbox.pdf"),
  plot = p_all,
  width = 16,
  height = 11,
  units = "in",
  device = cairo_pdf,
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(outdir, "UMAP_LARRY_states_plus_markers_labelbox.jpg"),
  plot = p_all,
  width = 16,
  height = 11,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

print(p_all)