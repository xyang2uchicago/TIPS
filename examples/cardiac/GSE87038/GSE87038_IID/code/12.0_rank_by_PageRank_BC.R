source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "GSE87038", "GSE87038_IID/")
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

df_PageRank <- readRDS(file = paste0(db_specifc_output_path, "IID_df_PAGERANK_strength_ANND.rewiring.P.rds"))

dim(df_PageRank) # 6786   16
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
#    CTS_11    CTS_15    CTS_16  CTS_16.1     CTS_7     CTS_8     HiG_1    HiG_10
#        18        37        19        43         9        38       261       363
#    HiG_11    HiG_12    HiG_13    HiG_14    HiG_15    HiG_16    HiG_17    HiG_18
#       379       292       348       315       342       458       453       399
#    HiG_19     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_7     HiG_8
#       446       373       357       253       264       449       278       278
#     HiG_9 HiGCTS_15
#       303        11


df_PageRank$gene <- toupper(df_PageRank$gene)

res_pr <- rank_TF_CHD_in_PPIN(df_PageRank, CHD, TF_human,
  signatures = c("HiG_8", "CTS_8"),
  key = "PageRank",
  top_TF_rank = 3, gene_top_n = 20, saveFigure = TRUE
)
#  => CP_rank_gene_by_pageRank.pdf

################################
df_betweenness <- read.table(file = paste0(db_specifc_output_path, "df_betweeness.tsv"), sep = "\t", header = T)
dim(df_betweenness) # 7265    5
colnames(df_betweenness)
# [1] "signature"             "BetweennessCentrality" "gene"
# [4] "rank_by_BC"            "PPI_cat"
table(df_betweenness$signature)
#    CTS_11    CTS_15    CTS_16  CTS_16.1     CTS_7     CTS_8     HiG_1    HiG_10
#        46        60        34        73        26        51       273       380
#    HiG_11    HiG_12    HiG_13    HiG_14    HiG_15    HiG_16    HiG_17    HiG_18
#       400       309       367       333       359       477       475       420
#    HiG_19     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_7     HiG_8
#       470       393       376       270       281       465       291       292
#     HiG_9 HiGCTS_15
#       316        28

df_betweenness$gene <- toupper(df_betweenness$gene)

res_bw <- rank_TF_CHD_in_PPIN(df_betweenness, CHD, TF_human,
  signatures = c("HiG_8", "CTS_8"),
  key = "BetweennessCentrality",
  top_TF_rank = 3, gene_top_n = 20, saveFigure = TRUE
)
#  => CP_rank_gene_by_BetweennessCentrality.pdf

# combined_table <- df_PageRank %>%
#   filter(grepl("_cardiac.a", signature)) %>%
#   select(signature, gene, PageRank) %>%
#   left_join(
#     df_betweenness %>%
#       filter(grepl("_cardiac.a", signature)) %>%
#       select(signature, gene, BetweennessCentrality),
#     by = c("signature", "gene")
#   ) %>%
#   mutate(
#     is_CHD = tolower(gene) %in% tolower(CHD),
#     is_TF  = gene %in% TF_mouse
#   )

# write.table(combined_table, "pagerank_betweenness_combined_cardiac.a.txt", sep = "\t", row.names = FALSE, quote = FALSE)


##################################################
## identify among the top_TF_rank (=3) TFs BetweennessCentrality > 0

(keyTF_cardiac.a <- subset(res_pr[["CTS_8"]])$gene) #  "NR2F1" "IRX5"  "ALX1"

# use the above 3

# (keyTF_cardiac.a <- intersect(keyTF_cardiac.a, subset(res_bw[["CTS_8"]], BetweennessCentrality > 0)$gene)) #  "IRX5" "ALX1"

seed_TF <- keyTF_cardiac.a

####
# STRING originally used HiGCTS_8, IID doesnt have that so we are using CTS_8 instead. Need to double check with Xinan
# if this is a safe replacement, or if I should try using 2018 instead.
# If HiG would be used instead, we would get "NKX2-5" for 136 and 138

library(dplyr)
library(tidyr)

all_sigs <- intersect(unique(df_PageRank$signature), unique(df_betweenness$signature))

res_pr <- rank_TF_CHD_in_PPIN(
  df_PageRank,
  CHD,
  TF_human,
  signatures = all_sigs,
  key = "PageRank",
  top_TF_rank = 3,
  gene_top_n = 20,
  saveFigure = FALSE
)

res_bw <- rank_TF_CHD_in_PPIN(
  df_betweenness,
  CHD,
  TF_human,
  signatures = all_sigs,
  key = "BetweennessCentrality",
  top_TF_rank = 3,
  gene_top_n = 20,
  saveFigure = FALSE
)

common_sigs <- intersect(names(res_pr), names(res_bw))

keyTF_pr_summary <- bind_rows(lapply(common_sigs, function(sig) {
  pr_genes <- res_pr[[sig]]$gene

  data.frame(
    signature = sig,
    genes = if (length(pr_genes) == 0) "" else paste(pr_genes, collapse = ", "),
    stringsAsFactors = FALSE
  )
}))

keyTF_overlap_summary <- bind_rows(lapply(common_sigs, function(sig) {
  pr_genes <- res_pr[[sig]]$gene
  bw_genes <- res_bw[[sig]]$gene[res_bw[[sig]]$BetweennessCentrality > 0]
  overlap_genes <- intersect(pr_genes, bw_genes)

  data.frame(
    signature = sig,
    genes = if (length(overlap_genes) == 0) "" else paste(overlap_genes, collapse = ", "),
    stringsAsFactors = FALSE
  )
}))

print(as_tibble(keyTF_pr_summary), n = Inf)
# A tibble: 26 × 2
#    signature genes
#    <chr>     <chr>
#  1 CTS_11    ""
#  2 CTS_15    "GATA1, TAL1, CEBPD"
#  3 CTS_16    "HOXB7, CEBPD, ZEB1"
#  4 CTS_16.1  "RUNX1"
#  5 CTS_7     "SPI1, KLF1, GATA2"
#  6 CTS_8     "NR2F1, IRX5, ALX1"
#  7 HiG_1     "MYCN"
#  8 HiG_10    "JUN"
#  9 HiG_11    "JUN"
# 10 HiG_12    ""
# 11 HiG_13    "MYC, JUN"
# 12 HiG_14    "MYC"
# 13 HiG_15    "JUN"
# 14 HiG_16    ""
# 15 HiG_17    "JUN, FOS"
# 16 HiG_18    ""
# 17 HiG_19    ""
# 18 HiG_2     "JUN, FOS"
# 19 HiG_3     ""
# 20 HiG_4     ""
# 21 HiG_5     ""
# 22 HiG_6     "JUN"
# 23 HiG_7     "MYC"
# 24 HiG_8     "NKX2-5"
# 25 HiG_9     ""
# 26 HiGCTS_15 ""
print(as_tibble(keyTF_overlap_summary), n = Inf)
# A tibble: 26 × 2
#    signature genes
#    <chr>     <chr>
#  1 CTS_11    ""
#  2 CTS_15    "GATA1"
#  3 CTS_16    "HOXB7, ZEB1"
#  4 CTS_16.1  "RUNX1"
#  5 CTS_7     "SPI1"
#  6 CTS_8     "IRX5, ALX1"
#  7 HiG_1     ""
#  8 HiG_10    ""
#  9 HiG_11    "JUN"
# 10 HiG_12    ""
# 11 HiG_13    "MYC, JUN"
# 12 HiG_14    "MYC"
# 13 HiG_15    "JUN"
# 14 HiG_16    ""
# 15 HiG_17    "JUN"
# 16 HiG_18    ""
# 17 HiG_19    ""
# 18 HiG_2     ""
# 19 HiG_3     ""
# 20 HiG_4     ""
# 21 HiG_5     ""
# 22 HiG_6     ""
# 23 HiG_7     "MYC"
# 24 HiG_8     "NKX2-5"
# 25 HiG_9     ""
# 26 HiGCTS_15 ""
