require(dplyr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(scales) # for color gradient
library(patchwork) # to arrange plots
library(gridExtra) # to arrange plots

########## BEGINNING OF USER INPUT ##########

wd = "/Users/felixyu/Documents/IbarraSoria2018/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- paste0(wd, "data/")

celltype_specific_weight_version <- '10'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

db <- "IbarraSoria2018"

s <- "combined" # specificity method

# Critical transition clusters
CT_id <- c(
  "cardiac.a", "endothelial.b"
)

# Non-cardiac HiG clusters
noncardiac_id <- c(
    "blood", "endothelial.a", "endothelial.b", "endothelial.c",
    "endothelial.d", "mesodermProgenitors"
)

########## END OF USER INPUT ##########

file <- paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [5] "HiG_blood"                  "HiG_mesodermProgenitors"    "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"  
#  [9] "HiG_somiticMesoderm"        "HiG_mixedMesoderm.a"        "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"              "HiG_endothelial.b"          "HiG_cardiac.a"             
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"           "CTS_endothelial.b"          "CTS_cardiac.a"    

# Filter cardiac, non-critical transition clusters
ids <- unique(sub(".*_", "", names(graph_list)))
other_cardiac_id <- setdiff(ids, c(CT_id, noncardiac_id))

###################
## assess CHD scores across three PPI_cat
#################
CHD <- readRDS(file = paste0(inputdir, "CHD_Cilia_Genelist.rds"))
names(CHD)
CHD <- CHD$Griffin2023_PCGC_AllCurated

############################################
#### plot network for HiGCTS, each page is one HiG&CTS in 4 PPI threshiolds
# refer to 11.1-weighted_CTS_cardiac_network.R & 11.2.0_update_network_weights.R to see how the graph is built
############################################
p_listoflist = list()

for (int in c(paste0("HiGCTS_", CT_id))) {
    g <- graph_list[[int]]
	g = graph_list[[int]]
	
	p_listoflist[[int]]  = plot_weighted_PPIN(g, layout = "fr", 
		CHD = CHD, node_size_title = "|Wilcox score|")
}

(names(p_listoflist))
# "HiGCTS_endothelial.b" "HiGCTS_cardiac.a" 


pdf(file = paste0("network_view_PPI_", db, ".pdf"), width = 12)
grid.arrange(p_listoflist[[1]], p_listoflist[[2]], ncol = 2)
dev.off()

#########################################################
## OPTION: plot for Hig PPIN, too dense to show
#########################################################

graph_list <- readRDS(file = paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds"))

p_list_HiG <- list()
for (int in grep("^HiG_", names(graph_list), value = TRUE)) {
		g = graph_list[[int]]
		
		p_list_HiG[[int]]  = plot_weighted_PPIN(g, layout = "fr", 
		CHD = CHD, node_size_title = "|Wilcox score|")
}


(names(p_list_HiG))
# [1] "HiGCTS_endothelial.b" "HiGCTS_cardiac.a"    
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [5] "HiG_blood"                  "HiG_mesodermProgenitors"    "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"  
#  [9] "HiG_somiticMesoderm"        "HiG_mixedMesoderm.a"        "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"              "HiG_endothelial.b"          "HiG_cardiac.a" 

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
