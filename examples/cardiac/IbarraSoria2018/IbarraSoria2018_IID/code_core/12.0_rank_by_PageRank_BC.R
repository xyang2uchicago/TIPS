########## BEGINNING OF USER INPUT ##########
source(here::here("examples", "config.R"))
wd          <- tips_path("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID/")
shared_path <- paste0(shared_data_path, "/")
celltype_specific_weight_version <- '10'
CP_cluster  <- "cardiac.a"   # focal cluster; drives HiG_X, CTS_X, HiGCTS_X signatures
top_TF_rank <- 3
gene_top_n  <- 20
########## END OF USER INPUT ##########
setwd(paste0(wd, "results_core/PPI_weight/"))

library(clusterProfiler)
library(msigdbr)
packageVersion("msigdbr") # 25.1.1
library(dplyr)
library(ggplot2)
library(patchwork)

library(dorothea)
data(dorothea_mm, package = "dorothea")
TF_mouse <- unique(dorothea_mm$tf)
length(TF_mouse) # 1113
TF_human <- unique(dorothea_hs$tf)
length(TF_human) # 1333

coding_genes <- readRDS(file = paste0(shared_path, "coding_genes.rds")) %>% unique()
length(coding_genes) # 19930
names(coding_genes) <- NULL
CHD <- NULL

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

#############################
db_specifc_output_path <- paste0(wd, "results_core/PPI_weight/")

# IbarraSoria2018_IID uses no prefix on the PageRank file (vs GSE87038_IID which has "IID_" prefix)
df_PageRank <- readRDS(file = paste0(db_specifc_output_path, "df_PAGERANK_strength_ANND.rewiring.P.rds"))

dim(df_PageRank) # 5626   16
colnames(df_PageRank)
# [1] "signature"                   "gene"
# [3] "PageRank"                    "PPI_cat"
# [5] "EigenCentrality"             "p.PageRank"
# [7] "rank_by_p.PR"                "rank_by_PR"
# [9] "annd"                        "p.annd"
# [11] "strength"                    "rank_by_strength"
# [13] "normalized.strength"         "rank_by_normalized.strength"
# [15] "rank_by_ANND"                "rank_by_p.ANND"
table(df_PageRank$signature)
#              CTS_cardiac.a          CTS_endothelial.b
#                         19                          7
#                  HiG_blood              HiG_cardiac.a
#                        360                        360
# ...
#           HiGCTS_cardiac.a       HiGCTS_endothelial.b
#                          8                          2

signatures_found <- intersect(
    c(paste0("HiG_", CP_cluster), paste0("CTS_", CP_cluster), paste0("HiGCTS_", CP_cluster)),
    unique(df_PageRank$signature)
)
# HiGCTS_cardiac.a EXISTS in IID (unlike GSE87038_IID which lacks HiGCTS_8)
HiGCTS_sig <- paste0("HiGCTS_", CP_cluster)

df_PageRank$gene <- toupper(df_PageRank$gene)

res_pr <- rank_TF_CHD_in_PPIN(df_PageRank, CHD, TF_human,
    signatures = signatures_found,
    key = "PageRank",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
)
#  => CP_rank_gene_by_pageRank.pdf

################################
df_betweenness <- read.table(file = paste0(db_specifc_output_path, "df_betweeness.tsv"), sep = "\t", header = T)
dim(df_betweenness) # 5982    5
colnames(df_betweenness)
# [1] "signature"             "BetweennessCentrality" "gene"
# [4] "rank_by_BC"            "PPI_cat"

df_betweenness$gene <- toupper(df_betweenness$gene)

res_bw <- rank_TF_CHD_in_PPIN(df_betweenness, CHD, TF_human,
    signatures = signatures_found,
    key = "BetweennessCentrality",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
)
#  => CP_rank_gene_by_BetweennessCentrality.pdf

##################################################
## HiGCTS_cardiac.a exists but original 12.0 used CTS_cardiac.a directly for seed_TF
if (HiGCTS_sig %in% signatures_found) {
    message("HiGCTS_cardiac.a found in IID — matching original 12.0 behavior: using CTS_cardiac.a for seed_TF")
}
seed_TF <- subset(res_pr[[paste0("CTS_", CP_cluster)]])$gene
# "GATA4" "TBX5"  "MSX2"
