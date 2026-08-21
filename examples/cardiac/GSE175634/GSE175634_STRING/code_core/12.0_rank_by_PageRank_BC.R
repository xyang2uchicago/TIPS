########## BEGINNING OF USER INPUT ##########
code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_STRING", "code_core")
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

## get the TF database
library(dorothea)
data(dorothea_hs, package = "dorothea")
TF_human <- unique(dorothea_hs$tf) # gene symbols
length(TF_human) # [1] 1333

source(here::here("R", paste0("celltype_specific_weight_v", celltype_specific_weight_version, ".R")))

#############################
df_PageRank <- readRDS(file = paste0(db_specifc_output_path, "df_PAGERANK_strength_ANND.rewring.P.rds"))

dim(df_PageRank)
# [1] 4495   16
colnames(df_PageRank)
# [1] "signature"                    "gene"
# [3] "PageRank"                     "PPI_cat"
# [5] "EigenCentrality"              "p.PageRank"
# [7] "rank_by_p.PR"                 "rank_by_PR"
# [9] "annd"                         "p.annd"
# [11] "strength"                    "rank_by_strength"
# [13] "normalized.strength"         "rank_by_normalized.strength"
# [15] "rank_by_ANND"                "rank_by_p.ANND"
table(df_PageRank$signature)
# CTS_CP:51  CTS_CP.1:61  CTS_endoderm:60  CTS_muscle:57
# HiG_0:325  HiG_1:1154   HiG_10:70        HiG_12:20
# HiG_2:240  HiG_4:102    HiG_5:337        HiG_6:102
# HiG_7:104  HiG_9:1466   HiG_CP:129       HiG_endoderm:37
# HiG_muscle:118  HiGCTS_CP:18  HiGCTS_CP.1:20
# HiGCTS_endoderm:6  HiGCTS_muscle:18

signatures_found <- intersect(
    c(paste0("HiG_", CP_cluster), paste0("CTS_", CP_cluster), paste0("HiGCTS_", CP_cluster)),
    unique(df_PageRank$signature)
)
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
dim(df_betweenness)
# [1] 4583    5
colnames(df_betweenness)
# [1] "signature"             "BetweennessCentrality" "gene"
# [4] "rank_by_BC"            "PPI_cat"
table(df_betweenness$signature)
# CTS_CP:72  CTS_CP.1:77  CTS_endoderm:64  CTS_muscle:61
# HiG_0:325  HiG_1:1155   HiG_10:70        HiG_12:22
# HiG_2:242  HiG_4:104    HiG_5:339        HiG_6:106
# HiG_7:104  HiG_9:1467   HiG_CP:132       HiG_endoderm:42
# HiG_muscle:118  HiGCTS_CP:29  HiGCTS_CP.1:27
# HiGCTS_endoderm:7  HiGCTS_muscle:20

df_betweenness$gene <- toupper(df_betweenness$gene)

res_bw <- rank_TF_CHD_in_PPIN(df_betweenness, CHD, TF_human,
    signatures = signatures_found,
    key = "BetweennessCentrality",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
)
#  => CP_rank_gene_by_BetweennessCentrality.pdf

##################################################
## identify among the top_TF_rank (=3) TFs BetweennessCentrality > 0

if (HiGCTS_sig %in% signatures_found) {
    (seed_TF <- subset(res_pr[[HiGCTS_sig]])$gene)
    (seed_TF <- intersect(
        seed_TF,
        subset(res_bw[[HiGCTS_sig]], BetweennessCentrality > 0)$gene
    ))
} else {
    message("No ", HiGCTS_sig, " in data — seed_TF not determined from HiGCTS")
}
# seed_TF = "HOXB2"
# (top PageRank in HiGCTS_CP, only TF with BetweennessCentrality > 0)
