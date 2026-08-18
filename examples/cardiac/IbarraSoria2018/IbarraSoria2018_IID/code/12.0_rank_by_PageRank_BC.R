source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID/")
setwd(paste0(wd, "results/PPI_weight/"))


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

shared_path <- paste0(shared_data_path, "/")

coding_genes <- readRDS(file = paste0(shared_path, "coding_genes.rds")) %>% unique()
length(coding_genes) # 19930
names(coding_genes) <- NULL
CHD <- readRDS(paste0(shared_path, "CHD_Cilia_Genelist.rds"))
CHD <- unlist(CHD[c("Griffin2023_PCGC_AllCurated")])

source(here::here("R", "celltype_specific_weight_v10.R"))

#############################
db_specifc_output_path <- paste0(wd, "results/PPI_weight/")

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
#              HiG_cardiac.b              HiG_cardiac.c
#                        374                        405
#          HiG_endothelial.a          HiG_endothelial.b
#                        402                        373
#          HiG_endothelial.c          HiG_endothelial.d
#                        424                        368
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                        322                        309
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                        291                        325
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                        331                        296
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                        361                        289
#           HiGCTS_cardiac.a       HiGCTS_endothelial.b
#                          8                          2


df_PageRank$gene <- toupper(df_PageRank$gene)

res_pr <- rank_TF_CHD_in_PPIN(df_PageRank, CHD, TF_human,
  signatures = c("HiG_cardiac.a", "CTS_cardiac.a", "HiGCTS_cardiac.a"),
  key = "PageRank",
  top_TF_rank = 3, gene_top_n = 20, saveFigure = TRUE
)
#  => CP_rank_gene_by_pageRank.pdf

#############################
df_betweenness <- read.table(file = paste0(db_specifc_output_path, "df_betweeness.tsv"), sep = "\t", header = T)
dim(df_betweenness) # 5982    5
colnames(df_betweenness)
# [1] "signature"             "BetweennessCentrality" "gene"
# [4] "rank_by_BC"            "PPI_cat"
table(df_betweenness$signature)
#              CTS_cardiac.a          CTS_endothelial.b
#                         37                         27
#                  HiG_blood              HiG_cardiac.a
#                        374                        384
#              HiG_cardiac.b              HiG_cardiac.c
#                        391                        421
#          HiG_endothelial.a          HiG_endothelial.b
#                        424                        389
#          HiG_endothelial.c          HiG_endothelial.d
#                        447                        392
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                        345                        330
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                        309                        343
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                        345                        315
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                        379                        305
#           HiGCTS_cardiac.a       HiGCTS_endothelial.b
#                         13                         12

df_betweenness$gene <- toupper(df_betweenness$gene)

res_bw <- rank_TF_CHD_in_PPIN(df_betweenness, CHD, TF_human,
  signatures = c("HiG_cardiac.a", "CTS_cardiac.a", "HiGCTS_cardiac.a"),
  key = "BetweennessCentrality",
  top_TF_rank = 3, gene_top_n = 20, saveFigure = TRUE
)
#  => CP_rank_gene_by_BetweennessCentrality.pdf

##################################################
## identify among the top_TF_rank (=5) TFs BetweennessCentrality > 0
## the output will pass to code 24.0xxxx

(keyTF_cardiac.a <- subset(res_pr[["CTS_cardiac.a"]])$gene) #  [1] "GATA4" "TBX5"  "MSX2"

# (keyTF_cardiac.a <- intersect(keyTF_cardiac.a, subset(res_bw[["HiGCTS_cardiac.a"]], BetweennessCentrality > 0)$gene)) #  "MEF2C" "GATA4" "MSX2"

seed_TF_2018 <- keyTF_cardiac.a

# rank_TF_CHD_in_PPIN(df_betweenness, CHD, TF_mouse,
# signatures=c('HiG_8','CTS_8', 'HiGCTS_8'),
# key = 'BetweennessCentrality',
# gene_top_n = 5, int_top_n = 20, saveFigure=TRUE)

# combined_table <- df_PageRank %>%
# filter(grepl("_8", signature)) %>%
# select(signature, gene, PageRank) %>%
# left_join(
# df_betweenness %>%
# filter(grepl("_8", signature)) %>%
# select(signature, gene, BetweennessCentrality),
# by = c("signature", "gene")
# ) %>%
# mutate(
# is_CHD = tolower(gene) %in% tolower(CHD),
# is_TF  = gene %in% TF_mouse
# )

# write.table(combined_table, "pagerank_betweenness_combined_8.txt", sep = "\t", row.names = FALSE, quote = FALSE)
