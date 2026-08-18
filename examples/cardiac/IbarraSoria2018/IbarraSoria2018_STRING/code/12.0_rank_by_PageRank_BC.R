########## BEGINNING OF USER INPUT ##########
source(here::here("examples", "config.R"))
wd          <- tips_path("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_STRING/")
shared_path <- paste0(shared_data_path, "/")
celltype_specific_weight_version <- '10'
CP_cluster  <- "cardiac.a"   # focal cluster; drives HiG_X, CTS_X, HiGCTS_X signatures
top_TF_rank <- 3              # top N TFs to report per signature
gene_top_n  <- 20             # top N genes to label in plot
########## END OF USER INPUT ##########
setwd(paste0(wd, "results/"))

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
#  [21] "SP140L"  "SPZ1"    "THAP10"  "THAP5"   "THAP6"   "THAP8"   "THAP9"   "TP53"    "TP63"    "TP73"
#  [31] "TRERF1"  "VENTX"   "YBX1"    "YY2"     "ZBED1"   "ZBED2"   "ZBED5"   "ZBED6"   "ZBTB47"  "ZFP69B"
#  [41] "ZFY"     "ZIM3"    "ZNF10"   "ZNF100"  "ZNF101" .....
setdiff(toupper(TF_mouse), TF_human)
#  [1] "4932411N23RIK" "AW822073"      "B020011L13RIK" "BC025920"      "DUXBL1"        "DUXBL2"
#   [7] "DUXBL3"        "DUXF3"         "E430018J23RIK" "HMGA1B"        "NKX2-9"        "PIN1RT1"
#  [13] "REX2"          "RHOX13"        "RSL1"          "RSLCAN18"      "SRCAP"         "TRP53"
#  [19] "TRP63"         "TRP73"         "ZFP105"        "ZFP108"        "ZFP109"        "ZFP11"
#  [25] "ZFP110"        "ZFP112"        "ZFP114"        "ZFP128"        "ZFP13"         "ZFP131"
#  [31] "ZFP141"        "ZFP142"        "ZFP143"        "ZFP146"        "ZFP148"        "ZFP160"
#  [37] "ZFP169"        "ZFP174"        "ZFP180"        "ZFP182"        "ZFP189"        "ZFP202"
#  [43] "ZFP207"        "ZFP212"        "ZFP213"        "ZFP217"        "ZFP219"        "ZFP235"
#  [49] "ZFP236"        "ZFP239"        "ZFP24"         "ZFP248"        "ZFP251"        "ZFP260"
#  [55] "ZFP263"        "ZFP266"        "ZFP273"        "ZFP276"        "ZFP277"        "ZFP280B"
#  [61] "ZFP280C"       "ZFP280D"       "ZFP281"        "ZFP282"        "ZFP286"        "ZFP287"
#  [67] "ZFP292"        "ZFP296"        "ZFP316"        "ZFP317"        "ZFP318"        "ZFP319"
#  [73] "ZFP322A"       "ZFP324"        "ZFP329"        "ZFP334"        "ZFP335"        "ZFP341"
#  [79] "ZFP354A"       "ZFP354B"       "ZFP354C"       "ZFP358"        "ZFP362"        "ZFP366"
#  [85] "ZFP369"        "ZFP382"        "ZFP383"        "ZFP384"        "ZFP386"        "ZFP39"
#  [91] "ZFP397"        "ZFP398"        "ZFP407"        "ZFP408"        "ZFP410"        "ZFP414"
#  [97] "ZFP420"        "ZFP422"        "ZFP423"        "ZFP426"        "ZFP429"        "ZFP438"
# ....

coding_genes <- readRDS(file = paste0(shared_path, "coding_genes.rds")) %>% unique()
length(coding_genes) # 19930
names(coding_genes) <- NULL
CHD <- readRDS(paste0(shared_path, "CHD_Cilia_Genelist.rds"))
CHD <- unlist(CHD[c("Griffin2023_PCGC_AllCurated")])

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

#############################
db_specifc_output_path <- paste0(wd, "results/PPI_weight/")

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


df_PageRank$gene <- toupper(df_PageRank$gene)

res_pr <- rank_TF_CHD_in_PPIN(df_PageRank, CHD, TF_human,
    signatures = c(paste0("HiG_", CP_cluster), paste0("CTS_", CP_cluster), paste0("HiGCTS_", CP_cluster)),
    key = "PageRank",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
)
#  => CP_rank_gene_by_pageRank.pdf

#############################
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
#          HiG_endothelial.a          HiG_endothelial.b
#                        455                        418
#          HiG_endothelial.c          HiG_endothelial.d
#                        480                        426
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                        367                        356
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                        329                        369
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                        371                        338
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                        412                        327
#           HiGCTS_cardiac.a       HiGCTS_endothelial.b
#                         13                         16

df_betweenness$gene <- toupper(df_betweenness$gene)

res_bw <- rank_TF_CHD_in_PPIN(df_betweenness, CHD, TF_human,
    signatures = c(paste0("HiG_", CP_cluster), paste0("CTS_", CP_cluster), paste0("HiGCTS_", CP_cluster)),
    key = "BetweennessCentrality",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
)
#  => CP_rank_gene_by_BetweennessCentrality.pdf

##################################################
## identify among the top_TF_rank (=3) TFs BetweennessCentrality > 0
## the output will pass to code 24.0xxxx

(seed_TF <- subset(res_pr[[paste0("HiGCTS_", CP_cluster)]])$gene) #  "MEF2C" "GATA4" "MSX2"

(seed_TF <- intersect(seed_TF, subset(res_bw[[paste0("HiGCTS_", CP_cluster)]], BetweennessCentrality > 0)$gene)) #  "MEF2C" "GATA4" "MSX2"
