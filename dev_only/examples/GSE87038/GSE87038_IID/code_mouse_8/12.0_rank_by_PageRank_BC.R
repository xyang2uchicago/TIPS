## dorothea_mm$tf is native case, not uppercase, so toupper(TF_mouse) keeps the
## mouse-native TF list case-matched to the uppercased gene column.
## IID has no HiGCTS_8 in graph_list (didn't pass the edge-count threshold) -- uses
## CTS_8 directly for seed_TF.

source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "GSE87038", "GSE87038_IID/")
dir.create(paste0(wd, "results_mouse_8/"), showWarnings = FALSE, recursive = TRUE)
setwd(paste0(wd, "results_mouse_8/"))

library(clusterProfiler) # clusterProfiler v4.6.0
library(msigdbr)
packageVersion("msigdbr") # 25.1.1
library(dplyr)
library(ggplot2)
library(patchwork)

## get the TF database
library(dorothea)
# human or mouse
data(dorothea_mm, package = "dorothea")
TF_mouse <- unique(dorothea_mm$tf)
length(TF_mouse) # [1] 1113
TF_human <- unique(dorothea_hs$tf) # gene symbols
length(TF_human) # 1333

shared_path <- paste0(shared_data_path, "/")

coding_genes <- readRDS(file = paste0(shared_path, "coding_genes.rds")) %>% unique()
length(coding_genes) # 19930
names(coding_genes) <- NULL
CHD <- readRDS(paste0(shared_path, "CHD_Cilia_Genelist.rds"))
CHD <- unlist(CHD[c("Griffin2023_PCGC_AllCurated")])

source(here::here("R", "celltype_specific_weight_v10.R"))

#############################
db_specifc_output_path <- paste0(wd, "results/PPI_weight/")

df_PageRank <- readRDS(file = paste0(db_specifc_output_path, "IID_df_PAGERANK_strength_ANND.rewiring.P.rds"))

dim(df_PageRank)
colnames(df_PageRank)
table(df_PageRank$signature)

## Original ATAC/ChIP pipeline keeps the uppercase-mouse-as-human-proxy convention
## (matches Ameen2022/Maven2023/Gao2019 human-symbol sets in 24.0); toupper(TF_mouse)
## here so the mouse-native TF list still matches the uppercased gene column.
df_PageRank$gene <- toupper(df_PageRank$gene)

res_pr <- rank_TF_CHD_in_PPIN(df_PageRank, CHD, toupper(TF_mouse),
  signatures = c("HiG_8", "CTS_8"),
  key = "PageRank",
  top_TF_rank = 3, gene_top_n = 20, saveFigure = TRUE
)
#  => CP_rank_gene_by_pageRank.pdf

################################
df_betweenness <- read.table(file = paste0(db_specifc_output_path, "df_betweeness.tsv"), sep = "\t", header = T)
dim(df_betweenness)
colnames(df_betweenness)
table(df_betweenness$signature)

df_betweenness$gene <- toupper(df_betweenness$gene)

res_bw <- rank_TF_CHD_in_PPIN(df_betweenness, CHD, toupper(TF_mouse),
  signatures = c("HiG_8", "CTS_8"),
  key = "BetweennessCentrality",
  top_TF_rank = 3, gene_top_n = 20, saveFigure = TRUE
)
#  => CP_rank_gene_by_BetweennessCentrality.pdf

##################################################
## identify among the top_TF_rank (=3) TFs BetweennessCentrality > 0

(keyTF_cardiac.a <- subset(res_pr[["CTS_8"]])$gene)

seed_TF <- keyTF_cardiac.a
