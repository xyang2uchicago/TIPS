require(dplyr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(scales) # for color gradient
library(patchwork) # to arrange plots
library(gridExtra) # to arrange plots

########## BEGINNING OF USER INPUT ##########

wd = "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- paste0(wd, "data/")

celltype_specific_weight_version <- '10'
# source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))
source(paste0('../../code/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

db <- "GSE87038"

s <- "combined" # specificity method

CT_id <- c("7", "8", "11", "13", "15", "16", "16.1") # critical transition clusters
noncardiac_id <- c("1", "5", "6", "7", "9", "10", "13", "15") # Non-cardiac clusters

excluded <- c("HiGCTS_13", "CTS_13") # network too small, don't graph

########## END OF USER INPUT ##########

file <- paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"       "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"      "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"
# [15] "HiG_11"      "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"
# [29] "CTS_16"      "CTS_16.1"    "CTS_13"      "CTS_8"       "HiGCTS_13"

ids <- unique(sub(".*_", "", names(graph_list)))
other_cardiac_id <- setdiff(ids, c(CT_id, noncardiac_id))

###################
## assess CHD scores across PPI_cat
#################
CHD <- readRDS(file = paste0(inputdir, "CHD_Cilia_Genelist.rds"))
names(CHD)
CHD <- CHD$Griffin2023_PCGC_AllCurated

############################################
#### plot network for HiGCTS, each page is one HiG&CTS in 4 PPI thresholds
# refer to 11.2.1_weighted_CTS_cardiac_network.R & 11.2.0_update_network_weights.R to see howthe graph is built
############################################
p_listoflist = list()

for (int in c(paste0("HiGCTS_", CT_id))) {
    if(int %in% excluded){
        next
    }
    g <- graph_list[[int]]
    graph_attr(g, "name") <- int
	
	p_listoflist[[int]]  = plot_weighted_PPIN(g, layout = "fr", 
		CHD = CHD, node_size_title = "|Wilcox score|")
}

(names(p_listoflist))
# [1] "HiGCTS_7"    "HiGCTS_8"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1"

n <- length(p_listoflist)
pdf(file = paste0("network_view_PPI_HiGCTS_", db, ".pdf"), width = 12, height = 12)
for (i in seq(1, n, by = 2)) {
  grobs <- p_listoflist[i:min(i + 1, n)]
  print(grid.arrange(grobs = grobs, ncol = 2))
}
dev.off()

p_listoflist = list()
for (k in c(CT_id)) {

    int <- paste0("CTS_", k)

    if (int %in% excluded) {
        next
    }

    g <- graph_list[[int]]
    graph_attr(g, "name") <- int

    HiG_int <- paste0("HiGCTS_", k)
    HiGCTS  <- graph_list[[HiG_int]]

    HiG <- V(g)$name %in% V(HiGCTS)$name

    V(g)$is_HiG <- HiG

    p_listoflist[[int]] = plot_weighted_PPIN(
        g,
        layout = "fr",
        CHD = CHD,
        node_size_title = "|Wilcox score|"
    )
}

(names(p_listoflist))
# [1] "HiGCTS_7"    "HiGCTS_8"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1"

n <- length(p_listoflist)
pdf(file = paste0("network_view_PPI_CTS_", db, ".pdf"), width = 12, height = 12)
for (i in seq(1, n, by = 2)) {
  grobs <- p_listoflist[i:min(i + 1, n)]
  print(grid.arrange(grobs = grobs, ncol = 2))
}
dev.off()


#########################################################
## OPTION: plot for HiG PPIN, too dense to show
#########################################################

graph_list <- readRDS(file = paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds"))


p_list_HiG <- list()
for (int in grep("^HiG_", names(graph_list), value = TRUE)) {
		g = graph_list[[int]]
		graph_attr(g, "name") <- int
		p_list_HiG[[int]]  = plot_weighted_PPIN(g, layout = "fr", 
		CHD = CHD, node_size_title = "|Wilcox score|")
}

(names(p_list_HiG))
# [1] "HiG_1"  "HiG_2"  "HiG_3"  "HiG_4"  "HiG_5"  "HiG_6"  "HiG_9"  "HiG_10" "HiG_12" "HiG_14" "HiG_17" "HiG_18" "HiG_19" "HiG_7"  "HiG_11" "HiG_15" "HiG_16" "HiG_13" "HiG_8"

groups <- list(
    CT = paste0("HiG_", CT_id),
    Other = paste0("HiG_", other_cardiac_id),
    NonC = paste0("HiG_", noncardiac_id)
)

pdf(file = paste0("network_HiG_view_weight_shrink_", db, ".pdf"), width = 20, height = 20)

for (g in groups) {
    # filter out subclusters not in p_list_HiG
    valid <- intersect(g, names(p_list_HiG))
    if (length(valid) > 0)
        print(grid.arrange(grobs = p_list_HiG[valid], ncol = 4))
}

dev.off()
