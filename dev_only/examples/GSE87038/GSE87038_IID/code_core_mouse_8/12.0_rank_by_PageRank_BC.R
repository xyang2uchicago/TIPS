## 12.0_rank_by_PageRank_BC.R — code_core_mouse_8 (IID, C8 arm, CORE)
## IID has no HiGCTS_8 in graph_list -- uses CTS_8 directly for seed_TF (see original
## GSE87038_IID/code/12.0...R note).
##
## IMPORTANT: unlike STRING's raw data, IID's saved df_PageRank/
## df_betweenness/graph_list are ALREADY uppercase in the source rds files (confirmed via
## R: gene column is "ALX1" not "Alx1"). Native-case TF_mouse therefore matches nothing
## here -- use toupper(TF_mouse) instead, same as the code_mouse_8 pattern. This does NOT
## affect 24.0/24.1's mm10 cisTarget matching, which uses CTS/DEG from the shared
## GSE87038/data/ files -- those ARE native case regardless of IID vs STRING (they don't
## depend on which PPI database was used).

source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "GSE87038", "GSE87038_IID/")
dir.create(paste0(wd, "results_core_mouse_8/"), showWarnings = FALSE, recursive = TRUE)
setwd(paste0(wd, "results_core_mouse_8/"))

library(clusterProfiler) # clusterProfiler v4.6.0
library(msigdbr)
packageVersion("msigdbr") # 25.1.1
library(dplyr)
library(ggplot2)
library(patchwork)

## get the TF database
library(dorothea)
data(dorothea_mm, package = "dorothea")
TF_mouse <- unique(dorothea_mm$tf) # gene symbols, native case
length(TF_mouse) # [1] 1113

shared_path <- paste0(shared_data_path, "/")

coding_genes <- readRDS(file = paste0(shared_path, "coding_genes.rds")) %>% unique()
length(coding_genes) # 19930
names(coding_genes) <- NULL
CHD <- NULL

source(here::here("R", "celltype_specific_weight_v10.R"))

#############################
db_specifc_output_path <- paste0(wd, "results/PPI_weight/")

df_PageRank <- readRDS(file = paste0(db_specifc_output_path, "IID_df_PAGERANK_strength_ANND.rewiring.P.rds"))

dim(df_PageRank)
colnames(df_PageRank)
table(df_PageRank$signature)

## df_PageRank$gene is already uppercase in IID's raw data (see header note) --
## toupper(TF_mouse) to match.
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

## df_betweenness$gene is already uppercase too — see note above res_pr.
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
