########## BEGINNING OF USER INPUT ##########
wd          <- "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/GSE87038_IID/"
shared_path <- paste0(wd, "../../Shared_Data/")
celltype_specific_weight_version <- '10'
CP_cluster  <- "13"  # focal cluster; drives HiG_X, CTS_X, HiGCTS_X signatures
top_TF_rank <- 3     # top N TFs to report per signature
gene_top_n  <- 20    # top N genes to label in plot
########## END OF USER INPUT ##########
setwd(paste0(wd, "results_core_13/"))

library(clusterProfiler)
library(msigdbr)
packageVersion("msigdbr") # 25.1.1
library(dplyr)
library(ggplot2)
library(patchwork)

## get the TF database
library(dorothea)
data(dorothea_mm, package = "dorothea")
TF_mouse <- unique(dorothea_mm$tf)
length(TF_mouse) # [1] 1113
TF_human <- unique(dorothea_hs$tf)
length(TF_human) # [1] 1333
setdiff(TF_human, toupper(TF_mouse))
setdiff(toupper(TF_mouse), TF_human)

coding_genes <- readRDS(file = paste0(shared_path, "coding_genes.rds")) %>% unique()
length(coding_genes) # 19930
names(coding_genes) <- NULL
CHD <- NULL

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

#############################
db_specifc_output_path <- paste0(wd, "results_core_13/PPI_weight/")

# IID uses IID-prefixed PageRank file (vs STRING's df_PAGERANK_strength_ANND.rewiring.P.rds)
df_PageRank <- readRDS(file = paste0(db_specifc_output_path, "IID_df_PAGERANK_strength_ANND.rewiring.P.rds"))

dim(df_PageRank) # 9574   16
colnames(df_PageRank)
# [1] "signature"                   "gene"
# [3] "PageRank"                    "EigenCentrality"
# [5] "PPI_cat"                     "strength"
# [7] "normalized.strength"         "annd"
# [9] "rank_by_PR"                  "rank_by_strength"
# [11] "rank_by_normalized.strength" "rank_by_ANND"
# [13] "p.PageRank"                  "rank_by_p.PR"
# [15] "p.annd"                      "rank_by_p.ANND"
table(df_PageRank$signature)
#    CTS_13    CTS_15  CTS_16.1     CTS_8     HiG_1    HiG_10    HiG_11    HiG_12
#        55        61        74        52       387       535       537       418
#    HiG_13    HiG_14    HiG_15    HiG_16    HiG_17    HiG_18    HiG_19     HiG_2
#       495       443       496       610       595       537       621       524
#     HiG_3     HiG_4     HiG_5     HiG_6     HiG_7     HiG_8     HiG_9 HiGCTS_15
#       509       385       388       584       374       411       442        41

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
dim(df_betweenness) # 9574    5
colnames(df_betweenness)
# [1] "signature"             "BetweennessCentrality" "gene"
# [4] "rank_by_BC"            "PPI_cat"

df_betweenness$gene <- toupper(df_betweenness$gene)

res_bw <- rank_TF_CHD_in_PPIN(df_betweenness, CHD, toupper(TF_mouse),
    signatures = signatures_found,
    key = "BetweennessCentrality",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
)
#  => CP_rank_gene_by_BetweennessCentrality.pdf

##################################################
## IID lacks HiGCTS_8; fall back to CTS_8 for seed_TF

if (HiGCTS_sig %in% signatures_found) {
    (seed_TF <- subset(res_pr[[HiGCTS_sig]])$gene)
    (seed_TF <- intersect(seed_TF, subset(res_bw[[HiGCTS_sig]], BetweennessCentrality > 0)$gene))
} else if (paste0("CTS_", CP_cluster) %in% signatures_found) {
    message("No ", HiGCTS_sig, " in IID data — using CTS_", CP_cluster, " for seed_TF")
    seed_TF <- subset(res_pr[[paste0("CTS_", CP_cluster)]])$gene
    # "BCL6B" "ELK3" "ETS1"
} else {
    message("No suitable signature found for seed_TF")
    seed_TF <- character(0)
}
