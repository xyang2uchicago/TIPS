########## BEGINNING OF USER INPUT ##########
code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_IID", "code_core")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
########## END OF USER INPUT ##########
setwd(results_dir)

library(clusterProfiler)
library(msigdbr)
packageVersion("msigdbr")
library(dplyr)
library(ggplot2)
library(patchwork)

library(dorothea)
data(dorothea_hs, package = "dorothea")
TF_human <- unique(dorothea_hs$tf)
length(TF_human) # 1333

source(here::here("R", paste0("celltype_specific_weight_v", celltype_specific_weight_version, ".R")))

#############################
df_PageRank <- readRDS(file = paste0(db_specifc_output_path, "df_PAGERANK_strength_ANND.rewring.P.rds"))

dim(df_PageRank)
# [1] 3941   16
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
#          CTS_CP        CTS_CP.1    CTS_endoderm      CTS_muscle
#              35              38              50              45
#           HiG_0           HiG_1          HiG_10          HiG_12
#             298            1068              55              13
#           HiG_2           HiG_4           HiG_5           HiG_6
#             207              72             293              77
#           HiG_7           HiG_9          HiG_CP    HiG_endoderm
#              85            1328              94              27
#      HiG_muscle       HiGCTS_CP     HiGCTS_CP.1 HiGCTS_endoderm
#             110              14              10               6
#   HiGCTS_muscle
#              16

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
colnames(df_betweenness)

df_betweenness$gene <- toupper(df_betweenness$gene)

res_bw <- rank_TF_CHD_in_PPIN(df_betweenness, CHD, TF_human,
    signatures = signatures_found,
    key = "BetweennessCentrality",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
)
#  => CP_rank_gene_by_BetweennessCentrality.pdf

##################################################
if (HiGCTS_sig %in% signatures_found) {
    (seed_TF <- subset(res_pr[[HiGCTS_sig]])$gene)
    # [1] "PRRX1" "HOXB2" "HOXB3"
    (seed_TF <- intersect(
        seed_TF,
        subset(res_bw[[HiGCTS_sig]], BetweennessCentrality > 0)$gene
    ))
    # seed_TF = c("PRRX1", "HOXB2")
    # (top PageRank in HiGCTS_CP with BetweennessCentrality > 0; HOXB3 excluded: BC = 0)
} else if (paste0("CTS_", CP_cluster) %in% signatures_found) {
    message("No ", HiGCTS_sig, " in IID data — using CTS_", CP_cluster, " for seed_TF")
    seed_TF <- subset(res_pr[[paste0("CTS_", CP_cluster)]])$gene
} else {
    message("No suitable signature found for seed_TF")
    seed_TF <- character(0)
}
