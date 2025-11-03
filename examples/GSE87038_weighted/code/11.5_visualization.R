require(dplyr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(scales) # for color gradient
library(patchwork) # to arrange plots
library(gridExtra) # to arrange plots

wd = "/Users/felixyu/Documents/GSE87038_weighted/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- paste0(wd, "data/")

db <- "GSE87038"

# refer to 11.2.0_weighted_graph_attack_robustness.R
s <- "combined"

CT_id <- c(7, 8, 11, 13, 15, 16, 16.1)
noncardiac_id <- c(1, 5, 6, 7, 9, 10, 13, 15) # Non-cardiac clusters
other_cardiac_id <- setdiff(1:19, c(CT_id, noncardiac_id))

file <- paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

celltype_specific_weight_version <- '10'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

(names(graph_list))
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"       "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"      "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"
# [15] "HiG_11"      "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"
# [29] "CTS_16"      "CTS_16.1"    "CTS_13"      "CTS_8"       "HiGCTS_13"

###################
## assess CHD scores across PPI_cat
#################
CHD <- readRDS(file = paste0(inputdir, "CHD_Cilia_Genelist.rds"))
names(CHD)
CHD <- CHD$Griffin2023_PCGC_AllCurated

############################################
#### plot network for HiGCTS, each page is one HiG&CTS in 4 PPI thresholds
# refer to 11.2.1_weighted_CTS_cardiac_network.R & 11.2.0_update_network_weights.R to  see the the graph is built
############################################
g = graph_list[["HiGCTS_8"]]
range(E(g)$weight)
# [1] 0.0006951995 0.1414810923

p_listoflist = list()

# HiGCTS_13 too small, not included here
for (int in c(paste0("HiGCTS_", CT_id))) {
    if(int == "HiGCTS_13"){
        next
    }
    g <- graph_list[[int]]
	g = graph_list[[int]]
	
	p_listoflist[[int]]  = plot_weighted_PPIN(g, layout = "fr", 
		CHD = CHD, node_size_title = "|Wilcox score|")
}

(names(p_listoflist))
# [1] "HiGCTS_7"    "HiGCTS_8"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1"

pdf(file = "network_view_PPI_GSE87038.pdf", width = 12, height = 12)
print(grid.arrange(grobs = list(p_listoflist[[1]], p_listoflist[[2]]), ncol = 2))
print(grid.arrange(grobs = list(p_listoflist[[3]], p_listoflist[[4]]), ncol = 2))
print(grid.arrange(grobs = list(p_listoflist[[5]], p_listoflist[[6]]), ncol = 2))
dev.off()


#########################################################
## OPTION: plot for HiG PPIN, too dense to show
#########################################################

score_threshold <- "weight_shrink"
input_path <- paste0("/Users/felixyu/Documents/GSE87038_weighted/results")

# refer to 11.1_CTS_cardiac_network_degreeDistribution.R
if (grepl("score", score_threshold)) {
    graph_list <- readRDS(file = paste0(input_path, "/GSE87038_STRING_graph_perState_notsimplified.rds"))
    graph_list <- lapply(graph_list, simplify, edge.attr.comb ='max') # !!!!!!!!!!!!!!!!!!! # FIXED
} else {
    graph_list <- readRDS(file = paste0(input_path, "/PPI_weight/GSE87038_STRING_graph_perState_simplified_", s, "weighted.rds"))
}

p_list_HiG <- list()
for (int in grep("^HiG_", names(graph_list), value = TRUE)) {
		g = graph_list[[int]]
		
		p_list_HiG[[int]]  = plot_weighted_PPIN(g, layout = "fr", 
		CHD = CHD, node_size_title = "|Wilcox score|")
}


(names(p_list_HiG))
# [1] "HiG_1"  "HiG_2"  "HiG_3"  "HiG_4"  "HiG_5"  "HiG_6"  "HiG_9"  "HiG_10" "HiG_12" "HiG_14" "HiG_17" "HiG_18" "HiG_19" "HiG_7"  "HiG_11" "HiG_15" "HiG_16" "HiG_13" "HiG_8"

pdf(file = paste0("network_HiG_view_", score_threshold, "_", db, ".pdf"), width = 20, height = 20)
print(grid.arrange(grobs = p_list_HiG[CT_id], ncol = 4))
print( grid.arrange(grobs = p_list_HiG[other_cardiac_id], ncol=4)) 
print( grid.arrange(grobs = p_list_HiG[noncardiac_id], ncol=4)) 

dev.off()
