
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
outdir = 'inVitro/'
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


## 2) read the metadata
meta <- read.delim(gzfile(paste0(data_path, "stateFate_inVitro_metadata.txt.gz")), 
    header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
dim(meta)  # 130887      8
head(meta, 3)
#   Library      Cell barcode Time point Starting population
# 1  d6_2_2 GCGTGCAA-AGAAGTTA          6       Lin-Kit+Sca1-
# 2  d6_2_2 AAGGGACC-CTCGATGC          6       Lin-Kit+Sca1-
# 3  d6_2_2 CGTACCGA-AGCGCCTT          6       Lin-Kit+Sca1-
#   Cell type annotation Well SPRING-x SPRING-y
# 1     Undifferentiated    2  411.496  -96.190
# 2     Undifferentiated    2 -587.462 -306.925
# 3             Monocyte    2 1429.805 -429.300

## Make metadata column names R-friendly:
## "Time point" -> "Time.point"
## "Cell type annotation" -> "Cell.type.annotation"
## "SPRING-x" -> "SPRING.x"
names(meta) <- make.names(names(meta))
colnames(meta)
# [1] "Library"              "Cell.barcode"         "Time.point"          
# [4] "Starting.population"  "Cell.type.annotation" "Well"                
# [7] "SPRING.x"             "SPRING.y"     

## Add stable cell IDs based on matrix/metadata row order.
## The LARRY trajectory indices are base-0, so cell_0 corresponds to row 1.
meta$cell_index0 <- seq_len(nrow(meta)) - 1
meta$cell_id <- paste0("cell_", meta$cell_index0)
rownames(meta) <- meta$cell_id

## Check required columns
required_cols <- c("Time.point", "Starting.population", "Cell.type.annotation",
                   "SPRING.x", "SPRING.y")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

## Clean plotting variables
meta$Time.point <- factor(meta$Time.point, levels = sort(unique(meta$Time.point)))
meta$Cell.type.annotation <- factor(meta$Cell.type.annotation)

## -----------------------------
## 3. Read neutrophil/monocyte trajectory indices
## -----------------------------

read_base0_index_file <- function(path) {
  ## Robust reader:
  ## accepts one-column files, files with headers, or whitespace-separated values.
  raw_lines <- readLines(gzfile(paste0(data_path, "stateFate_inVitro_neutrophil_pseudotime.txt.gz")))
  tokens <- unlist(strsplit(raw_lines, "\\s+"))
  tokens <- tokens[nzchar(tokens)]
  idx <- suppressWarnings(as.integer(tokens))
  idx <- idx[!is.na(idx)]
  sort(unique(idx))
}

nm_idx0 <- read_base0_index_file(files["nm_traj"])
nm_cells <- paste0("cell_", nm_idx0)

message("Total cells in metadata: ", nrow(meta)) # Total cells in metadata: 130,887
message("Neutrophil/monocyte trajectory cells: ", length(nm_cells)) # Neutrophil/monocyte trajectory cells: 92,527

meta$NM_trajectory <- rownames(meta) %in% nm_cells

## Useful display label:
meta$Fig6B_label <- ifelse(
  meta$NM_trajectory,
  as.character(meta$Cell.type.annotation),
  "outside_NM_trajectory"
)
meta$Fig6B_label <- factor(meta$Fig6B_label)

## -----------------------------
## 4. Create a Seurat object
##    using normalized counts  ('after total-count normalization')('rows are cells, columns are genes')
## -----------------------------
# the authors already did QC: Kept abundant cell barcodes based on total transcript histograms; 
# removed cells with >20% mitochondrial transcripts; 
# used SCRUBLET with manual thresholding; doublet filtering removed ~10% of cells.
# Cells within each experiment were then normalized to have the same total number of transcripts


## a Dummy Seurat object is for plotting metadata/Spring coordinates only.
## It avoids loading the full expression matrix unless needed.
# dummy_counts <- Matrix::sparseMatrix(
#   i = rep(1, nrow(meta)),
#   j = seq_len(nrow(meta)),
#   x = rep(1, nrow(meta)),
#   dims = c(1, nrow(meta)),
#   dimnames = list("dummy", rownames(meta))
# )
# larry <- CreateSeuratObject(
#   counts = dummy_counts,
#   meta.data = meta,
#   assay = "RNA",
#   min.cells = 0,
#   min.features = 0
# )

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


message("Reading expression matrix...")
mat_cells_by_genes <- Matrix::readMM(gzfile(paste0(data_path, "stateFate_inVitro_normed_counts.mtx.gz")))
mat_cells_by_genes <- as(mat_cells_by_genes, "dgCMatrix")
dim(mat_cells_by_genes) #[1] 130887  25289

gene_names <- readLines(gzfile(paste0(data_path, "stateFate_inVitro_gene_names.txt.gz")))
gene_names <- make.unique(gene_names)
head(gene_names)
# [1] "0610006L08Rik" "0610007P14Rik" "0610009B22Rik" "0610009E02Rik"
# [5] "0610009L18Rik" "0610009O20Rik"
length(gene_names) #[1] 25289   

stopifnot(nrow(mat_cells_by_genes) == nrow(meta))
stopifnot(ncol(mat_cells_by_genes) == length(gene_names))

## README: rows = cells, columns = genes.
## Seurat expects rows = genes, columns = cells.
expr_genes_by_cells <- Matrix::t(mat_cells_by_genes)

rownames(expr_genes_by_cells) <- gene_names
colnames(expr_genes_by_cells) <- rownames(meta)
dim(expr_genes_by_cells) #[1] 25289 130887
rm(mat_cells_by_genes)
gc()
range(expr_genes_by_cells)
# [1]   0.000 870.983   not yet log-normalized


larry <- CreateSeuratObject(
  counts = expr_genes_by_cells,
  meta.data = meta,
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  project = "LARRY_inVitro"
)

## Add provided SPRING coordinates too, useful for comparison with paper/browser
## They used SPRING force-directed layouts built from total-count-normalized expression, 
## filtered HVGs, 50 PCs, a kNN graph, and ForceAtlas2.

spring_mat <- as.matrix(larry@meta.data[, c("SPRING.x", "SPRING.y")])
mode(spring_mat) <- "numeric"
rownames(spring_mat) <- colnames(larry)
colnames(spring_mat) <- c("SPRING_1", "SPRING_2")

larry[["spring"]] <- CreateDimReducObject(
  embeddings = spring_mat,
  key = "SPRING_",
  assay = DefaultAssay(larry)
)

## 4. logtransform + scale (NOT re Normalize), variable genes, PCA, UMAP
# data (log-normalized) → expression plotting and DE: FeaturePlot, VlnPlot, DotPlot, FindMarkers
# scale.data (z-scored) → PCA → which then feeds UMAP and clustering
############################################################

DefaultAssay(larry) <- "RNA"

## The downloaded matrix is already total-count normalized.
range(larry[['RNA']]$counts)
#[1]   0.000 870.983

## This still creates a Seurat-normalized data layer for PCA/UMAP plotting.
# larry <- NormalizeData(
#   object = larry,
#   normalization.method = "LogNormalize",
#   scale.factor = 10000,
#   verbose = TRUE
# )

larry[['RNA']]$data <- log1p(larry[['RNA']]$counts)
range(larry[['RNA']]$data)
#[1] 0.00000 6.77077

larry <- FindVariableFeatures(
  object = larry,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)

## 4. Add cell-cycle signature score
larry <- AddModuleScore(
  object = larry,
  features = list(cell_cycle_signature),
  name = "CellCycleSig"
)

cc_score <- larry$CellCycleSig1

## 5. Correlate HVGs with cell-cycle score
expr_hvg <- GetAssayData(
  larry,
  assay = "RNA",
  layer = "data"
)[VariableFeatures(larry), , drop = FALSE]

cc_cor <- apply(
  expr_hvg,
  1,
  function(x) cor(
    as.numeric(x),
    cc_score,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
)

## Weinreb used R > 0.1.
# ## however, we can't repeat this setting
# cellcycle_correlated_genes <- names(cc_cor)[cc_cor > 0.2]
# setdiff(cellcycle_correlated_genes, cell_cycle_signature)
# #[1] "Elane"  "Prtn3"  "Ctsg"   "Plac8"  "H2afy"  "Ccdc62" "Ttf1"   "Igfbp4"
# intersect(cellcycle_correlated_genes, cell_cycle_signature)
# #[1] "Top2a"
# cellcycle_correlated_genes <- names(cc_cor)[cc_cor > 0.1]
# setdiff(cellcycle_correlated_genes, cell_cycle_signature)
# #  [1] "Elane"         "Mpo"           "Prtn3"         "Ctsg"         
# #  [5] "Muc13"         "Cd34"          "Plac8"         "H2afy"        
# #  [9] "Gm37214"       "Fry"           "Lipg"          "Myc"          
# # [13] "Lmo4"          "H1f0"          "Ccdc62"        "Ttf1"         
# # [17] "2810403D21Rik" "Igfbp4"        "Gfi1"     
# intersect(cellcycle_correlated_genes, cell_cycle_signature)
# #[1] "Top2a"
# message("Original HVGs: ", length(VariableFeatures(larry))) # 3000
# message("Removed cell-cycle-correlated HVGs: ", length(cellcycle_correlated_genes))  # 24
# write.table(cellcycle_correlated_genes, 
#   file = file.path(outdir, paste0(length(cellcycle_correlated_genes), "_cellcycle_correlated_genes_fullSeurat5.txt")), 
#   quote = FALSE, row.names = FALSE, col.names = FALSE)

hvg_no_cellcycle <- setdiff(
  VariableFeatures(larry),
  cell_cycle_signature
)
message("Remaining HVGs for PCA/UMAP: ", length(hvg_no_cellcycle)) #  2999

## 6. Scale only no-cell-cycle HVGs
larry <- ScaleData(
  larry,
  features = hvg_no_cellcycle,
  verbose = TRUE
)

larry <- RunPCA(
  object = larry,
  features = VariableFeatures(larry),
  npcs = 50,
  verbose = TRUE
)

## inspect this and adjust dims if needed
ElbowPlot(larry, ndims = 50) +
  geom_hline(yintercept = 1,   linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
  annotate("text", x = 45, y = 1,   label = "1",   vjust = -0.5, color = "red") +
  annotate("text", x = 45, y = 1.5, label = "1.5", vjust = -0.5, color = "red")
dev.copy2pdf(file = file.path(outdir, "Seurat5_noCellCycle_PCA_ElbowPlot.pdf"))

larry <- RunUMAP(
  object = larry,
  dims = 1:npcs,   ## was 50 in the original code without UMAP running 
  reduction = "pca",
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = TRUE
)

saveRDS(
  larry,
  file = file.path(outdir, "reanalyzed_Seurat5_noCellCycle_PCA_UMAP.rds")
)

############################################################
## 5.0. understand the data distribution
############################################################
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


