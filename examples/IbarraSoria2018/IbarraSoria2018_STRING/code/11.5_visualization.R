require(dplyr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(scales) # for color gradient
library(patchwork) # to arrange plots
library(gridExtra) # to arrange plots

########## BEGINNING OF USER INPUT ##########

wd = "/Users/felixyu/Documents/GitHub/TIPS/examples/IbarraSoria2018/IbarraSoria2018_STRING/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- paste0(wd, "../data/")
shared_path <- paste0(wd, "../../Shared_Data/")

celltype_specific_weight_version <- '10'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))
# source(paste0('../../code/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

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
#  [1] "HiG_blood"                  "HiG_cardiac.b"               "HiG_cardiac.c"              "HiG_endothelial.a"
#  [5] "HiG_endothelial.c"          "HiG_endothelial.d"           "HiG_extraembryonicMesoderm" "HiG_mesodermProgenitors"
#  [9] "HiG_mixedMesoderm.a"        "HiG_mixedMesoderm.b"         "HiG_pharyngealMesoderm"     "HiG_presomiticMesoderm.a"
# [13] "HiG_presomiticMesoderm.b"   "HiG_somiticMesoderm"         "HiG_endothelial.b"          "HiG_cardiac.a"
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"            "CTS_endothelial.b"          "CTS_cardiac.a"

ids <- unique(sub(".*_", "", names(graph_list)))
other_cardiac_id <- setdiff(ids, c(CT_id, noncardiac_id))

###################
## assess CHD scores across PPI_cat
#################
CHD <- readRDS(file = paste0(shared_path, "CHD_Cilia_Genelist.rds"))
names(CHD)
CHD <- CHD$Griffin2023_PCGC_AllCurated

############################################
#### plot network for HiGCTS, each page is one HiG&CTS in 4 PPI thresholds
# refer to 11.2.1_weighted_CTS_cardiac_network.R & 11.2.0_update_network_weights.R to see howthe graph is built
############################################
p_listoflist = list()

for (int in c(paste0("HiGCTS_", CT_id))) {
    g <- graph_list[[int]]
    graph_attr(g, "name") <- int
	
	p_listoflist[[int]]  = plot_weighted_PPIN(g, layout = "fr", 
		CHD = CHD, node_size_title = "|Wilcox score|")
}

(names(p_listoflist))
# [1] "HiGCTS_cardiac.a"     "HiGCTS_endothelial.b"

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

    g <- graph_list[[int]]
    graph_attr(g, "name") <- int

    HiG_int <- paste0("HiGCTS_", k)
    HiGCTS  <- graph_list[[HiG_int]]

    HiG <- V(g)$name %in% V(HiGCTS)$name

    V(g)$is_HiG <- HiG

    # CTS graphs have no DEG-based vertex weight; use uniform placeholder
    if (!"weight" %in% vertex_attr_names(g)) V(g)$weight <- 1

    p_listoflist[[int]] = plot_weighted_PPIN(
        g,
        layout = "fr",
        CHD = CHD,
        node_size_title = "|Wilcox score|"
    )
}

(names(p_listoflist))
# [1] "CTS_cardiac.a"     "CTS_endothelial.b"

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
#  [1] "HiG_blood"                  "HiG_cardiac.b"               "HiG_cardiac.c"              "HiG_endothelial.a"
#  [5] "HiG_endothelial.c"          "HiG_endothelial.d"           "HiG_extraembryonicMesoderm" "HiG_mesodermProgenitors"
#  [9] "HiG_mixedMesoderm.a"        "HiG_mixedMesoderm.b"         "HiG_pharyngealMesoderm"     "HiG_presomiticMesoderm.a"
# [13] "HiG_presomiticMesoderm.b"   "HiG_somiticMesoderm"         "HiG_endothelial.b"          "HiG_cardiac.a"

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