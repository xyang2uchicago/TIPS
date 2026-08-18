########## BEGINNING OF USER INPUT ##########
source(here::here("examples", "config.R"))
wd          <- tips_path("examples", "hematoendothelial", "IbarraSoria2018", "IbarraSoria2018_STRING/")
shared_path <- paste0(shared_data_path, "/")
celltype_specific_weight_version <- '10'
CP_cluster  <- "endothelial.b"   # focal cluster; drives HiG_X, CTS_X, HiGCTS_X signatures
top_TF_rank <- 3             # top N TFs to report per signature
gene_top_n  <- 20            # top N genes to label in plot
########## END OF USER INPUT ##########
setwd(paste0(wd, "results_core_endothelial.b/"))

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
length(TF_human) # [1] 1333
setdiff(TF_human, toupper(TF_mouse))
# [1] "AHRR"    "CASZ1"   "DUX4"    "DUXA"    "ETV7"    "FOXP1"   "GLI4"    "HES4"    "HHEX"    "HMGA2"
#  [11] "NKX2-8"  "PITX2"   "PITX3"   "POU5F1B" "RAX2"    "RHOXF1"  "SETBP1"  "SHOX"    "SMAD9"   "SP140"
# ....
setdiff(toupper(TF_mouse), TF_human)
# [1] "4932411N23RIK" "AW822073"      "B020011L13RIK" "BC025920"      "DUXBL1"        "DUXBL2"
# ....

coding_genes <- readRDS(file = paste0(shared_path, "coding_genes.rds")) %>% unique()
length(coding_genes) # 19930
names(coding_genes) <- NULL
CHD <- NULL

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

#############################
db_specifc_output_path <- paste0(wd, "results_core_endothelial.b/PPI_weight/")

df_PageRank <- readRDS(file = paste0(db_specifc_output_path, "df_PAGERANK_strength_ANND.rewiring.P.rds"))

dim(df_PageRank) # 6366   16
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
#                         30                         28
#                  HiG_blood              HiG_cardiac.a
#                        417                        403
#              HiG_cardiac.b              HiG_cardiac.c
#                        407                        443
#          HiG_endothelial.a          HiG_endothelial.b
#                        455                        416
#          HiG_endothelial.c          HiG_endothelial.d
#                        478                        425
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                        361                        351
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                        328                        365
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                        369                        335
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                        409                        323
#           HiGCTS_cardiac.a       HiGCTS_endothelial.b
#                         12                         11

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
dim(df_betweenness) # 6427    5
colnames(df_betweenness)
# [1] "signature"             "BetweennessCentrality" "gene"
# [4] "rank_by_BC"            "PPI_cat"
table(df_betweenness$signature)
#              CTS_cardiac.a          CTS_endothelial.b
#                         37                         33
#                  HiG_blood              HiG_cardiac.a
#                        420                        407
#              HiG_cardiac.b              HiG_cardiac.c
#                        409                        444
# ...
#           HiGCTS_cardiac.a       HiGCTS_endothelial.b
#                         13                         16

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
    (seed_TF <- subset(res_pr[[HiGCTS_sig]])$gene) # "TAL1" "GATA2" "ETV2"
    (seed_TF <- intersect(
        seed_TF,
        subset(res_bw[[HiGCTS_sig]], BetweennessCentrality > 0)$gene
    )) # "ETV2"
} else {
    message("No ", HiGCTS_sig, " in data — seed_TF not determined from HiGCTS")
}
