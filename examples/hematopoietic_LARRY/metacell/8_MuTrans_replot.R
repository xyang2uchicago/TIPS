
MuTrans_path <- "F:/projects/TIPS/results/GSE140802_lineage_tracking/larry/data/processed/larry/mutrans/"

setwd("F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/larry_BioTIP/")

pkgs <- c(
   "slingshot", "TrajectoryUtils", "MatrixGenerics",
   "DelayedMatrixStats", "sparseMatrixStats", "matrixStats"
)

sapply(pkgs, packageVersion)  # local / midway3
# $slingshot
# [1] 2 4 0  

# $TrajectoryUtils
# [1] 1 4 0  

# $MatrixGenerics
# [1]  1 10  0 

# $DelayedMatrixStats
# [1]  1 20  0 

# $sparseMatrixStats
# [1]  1 10  0  

# $matrixStats
# [1] 1 4 1   



library(SingleCellExperiment)
library(dplyr) 
library(Seurat)

 
seu <- readRDS( "../BioTIP_attractor/seu_attractor_MuTrans_HVG.rds")
seu
# An object of class Seurat 
# 3000 features across 1200 samples within 1 assay 
# Active assay: RNA (3000 features, 0 variable features)
#  2 layers present: counts, data
#  1 dimensional reduction calculated: umap

colnames(seu@meta.data)
#  [1] "orig.ident"                  "nCount_RNA"                 
#  [3] "nFeature_RNA"                "Cell.type.clean"            
#  [5] "Cell.type.clean_purity"      "n_genes"                    
#  [7] "n_genes_by_counts"           "log1p_n_genes_by_counts"    
#  [9] "total_counts"                "log1p_total_counts"         
# [11] "pct_counts_in_top_50_genes"  "pct_counts_in_top_100_genes"
# [13] "pct_counts_in_top_200_genes" "pct_counts_in_top_500_genes"
# [15] "total_counts_mt"             "log1p_total_counts_mt"      
# [17] "pct_counts_mt"               "total_counts_ribo"          
# [19] "log1p_total_counts_ribo"     "pct_counts_ribo"            
# [21] "total_counts_hb"             "log1p_total_counts_hb"      
# [23] "pct_counts_hb"               "leiden_0.5"                 
# [25] "land"                        "entropy"                    
# [27] "attractor"       

## Violin plot of land score across attractor

library(ggplot2)
quantile(seu$land, seq(0.1,1,0.1))
# assuming 'seu' object and that 'land' and 'attractor' are columns in @meta.data
#       10%       20%       30%       40%       50%       60%       70% 
#  3.336326  3.450830  3.789678  3.952399  4.311432  4.396681  5.153804 
#       80%       90%      100% 
#  5.399334  6.108439 14.438444 

summary(seu$land[seu$attractor == '4'])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   4.311   5.409   5.409   5.973   6.481  11.584 

hist(seu@meta.data$land, 100)
abline(v=4.5, lty=2, col='red')
quantile(seu$land, seq(0.1,1,0.1))

library(ggpubr)
land_cut <- 5.4
land_ylim <- c(2, 14)

# Histogram plot of 'land'
p1 <- ggplot(seu@meta.data, aes(x = land)) +
  geom_histogram(bins = 100,  color = "black") +  #fill = "#69b3a2",
  geom_vline(xintercept = land_cut, linetype = "dashed", color = "red") +
  labs(x = "Land Score", y = "Count", title = "Histogram of Land Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate(
    "text",
    x = land_cut + 1,
    y = 0.5 * max(table(base::cut(seu$land, 100))),
    label = paste0(land_cut, " (80%) "),
    color = "red",
    vjust = -0.5,
    size = 5
  )

# Violin plot of 'land' by 'attractor' (publication-style)
# attractor is a factor in Seurat: NEVER as.integer(seu$attractor) — that uses codes 1..14, not labels 0..13
attractor_levels <- as.character(sort(unique(as.numeric(as.character(seu$attractor)))))
seu$attractor_f <- factor(as.character(seu$attractor), levels = attractor_levels)
if (any(is.na(seu$attractor_f))) {
  stop("attractor_f has NA — check seu$attractor: ", paste(utils::head(which(is.na(seu$attractor_f))), collapse = ", "))
}

p2 <- ggplot(seu@meta.data, aes(x = attractor_f, y = land, fill = attractor_f, group = attractor_f)) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    color = "black",
    linewidth = 0.35,
    alpha = 1,
    bw = 0.4
  ) +
  geom_boxplot(
    width = 0.08,
    outlier.shape = NA,
    fill = NA,
    color = "black",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = stats::setNames(
      ifelse(attractor_levels %in% c("4", "12"), "#f4b6b6", "white"),
      attractor_levels
    ),
    guide = "none"
  ) +
  geom_hline(yintercept = land_cut, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_y_continuous(limits = land_ylim, breaks = seq(2, 14, by = 2)) +
  labs(x = "Attractor", y = "Land score", title = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = ifelse(attractor_levels %in% c("4", "12"), "red", "black")
    )
  )

# Combine the two plots in two rows
ggarrange(
  p1, p2,
  nrow = 2,
  heights = c(1, 1)
)

dev.copy2pdf(file='MuTrans_land_across_attractor_hist_violin.pdf')