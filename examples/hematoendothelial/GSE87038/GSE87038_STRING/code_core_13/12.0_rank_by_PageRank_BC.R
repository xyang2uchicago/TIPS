########## BEGINNING OF USER INPUT ##########
code_dir <- here::here("examples", "hematoendothelial", "GSE87038", "GSE87038_STRING", "code_core_13")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
########## END OF USER INPUT ##########
setwd(results_dir)

library(clusterProfiler) # clusterProfiler v4.6.0
library(msigdbr)
packageVersion("msigdbr") # 25.1.1
library(dplyr)
library(ggplot2)
library(patchwork)

## get the TF database (mouse)
library(dorothea)
data(dorothea_mm, package = "dorothea")
TF_mouse <- unique(dorothea_mm$tf)
length(TF_mouse) # [1] 1113

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

#############################
df_PageRank <- readRDS(file = paste0(db_specifc_output_path, "df_PAGERANK_strength_ANND.rewiring.P.rds"))

dim(df_PageRank) # 10549   16
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
#      CTS_11      CTS_13      CTS_15      CTS_16    CTS_16.1       CTS_7
#          51          60          66          39          79          31
#       CTS_8       HiG_1      HiG_10      HiG_11      HiG_12      HiG_13
#          54         418         574         570         450         529
#      HiG_14      HiG_15      HiG_16      HiG_17      HiG_18      HiG_19
#         472         533         654         646         574         671
#       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6       HiG_7
#         562         545         419         424         634         408
#       HiG_8       HiG_9   HiGCTS_11   HiGCTS_13   HiGCTS_15   HiGCTS_16
#         445         477          23          15          44          17
#  HiGCTS_16.1    HiGCTS_7    HiGCTS_8
#          36          17          12

signatures_found <- intersect(
    c(paste0("HiG_", CP_cluster), paste0("CTS_", CP_cluster), paste0("HiGCTS_", CP_cluster)),
    unique(df_PageRank$signature)
)
HiGCTS_sig <- paste0("HiGCTS_", CP_cluster)

df_PageRank$gene <- toupper(df_PageRank$gene)

res_pr <- rank_TF_CHD_in_PPIN(df_PageRank, CHD, toupper(TF_mouse),
    signatures = signatures_found,
    key = "PageRank",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
)
#  => CP_rank_gene_by_pageRank.pdf

################################
df_betweenness <- read.table(file = paste0(db_specifc_output_path, "df_betweeness.tsv"), sep = "\t", header = T)
dim(df_betweenness) # 10549    5
colnames(df_betweenness)
# [1] "signature"             "BetweennessCentrality" "gene"
# [4] "rank_by_BC"            "PPI_cat"
table(df_betweenness$signature)
#      CTS_11      CTS_13      CTS_15      CTS_16    CTS_16.1       CTS_7
#          51          60          66          39          79          31
#       CTS_8       HiG_1      HiG_10      HiG_11      HiG_12      HiG_13
#          54         418         574         570         450         529
#      HiG_14      HiG_15      HiG_16      HiG_17      HiG_18      HiG_19
#         472         533         654         646         574         671
#       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6       HiG_7
#         562         545         419         424         634         408
#       HiG_8       HiG_9   HiGCTS_11   HiGCTS_13   HiGCTS_15   HiGCTS_16
#         445         477          23          15          44          17
#  HiGCTS_16.1    HiGCTS_7    HiGCTS_8
#          36          17          12

df_betweenness$gene <- toupper(df_betweenness$gene)

res_bw <- rank_TF_CHD_in_PPIN(df_betweenness, CHD, toupper(TF_mouse),
    signatures = signatures_found,
    key = "BetweennessCentrality",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
)
#  => CP_rank_gene_by_BetweennessCentrality.pdf

##################################################
## identify among the top_TF_rank (=3) TFs BetweennessCentrality > 0

if (HiGCTS_sig %in% signatures_found) {
    (seed_TF <- subset(res_pr[[HiGCTS_sig]])$gene) # character(0)  (HiGCTS_13 TFs have BC=0 with fully corrected HiG graphs)
    (seed_TF <- intersect(
        seed_TF,
        subset(res_bw[[HiGCTS_sig]], BetweennessCentrality > 0)$gene
    )) # character(0)
} else {
    message("No ", HiGCTS_sig, " in data — seed_TF not determined from HiGCTS")
}
