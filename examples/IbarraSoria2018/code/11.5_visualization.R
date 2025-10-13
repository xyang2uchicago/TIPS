require(dplyr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(scales) # for color gradient
library(patchwork) # to arrange plots
library(gridExtra) # to arrange plots

wd = "/Users/felixyu/Documents/IbarraSoria2018/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- paste0(wd, "data/")

db <- "IbarraSoria2018"

# refer to 11.2.0_weighted_graph_attack_robustness.R
s <- "combined"
file <- paste0("2018_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [5] "HiG_blood"                  "HiG_mesodermProgenitors"    "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"  
#  [9] "HiG_somiticMesoderm"        "HiG_mixedMesoderm.a"        "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"              "HiG_endothelial.b"          "HiG_cardiac.a"             
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"           "CTS_endothelial.b"          "CTS_cardiac.a"    

###################
## assess CHD scores across three PPI_cat
#################
CHD <- readRDS(file = paste0(inputdir, "CHD_Cilia_Genelist.rds"))
names(CHD)
CHD <- CHD$Griffin2023_PCGC_AllCurated

############################################
#### plot network for HiGCTS, each page is one HiG&CTS in 4 PPI threshiolds
# refer to 11.1-weighted_CTS_cardiac_network.R & 11.2.0_update_network_weights.R to  see the the graph is built
############################################
p_listoflist <- list()

for (int in c("HiGCTS_endothelial.b", "HiGCTS_cardiac.a")) {
    g <- graph_list[[int]]
    (vertex_attr_names(g)) # [1] "name"   "weight" "FDR"    "size"   "color"
    (edge_attr_names(g))
    # [1] "combined_score"  "weight"          "width"
    # [4] "original_weight" "corexp_sign"     "coexp_target"

    (V(g)$FDR)
    #  [1] 2.443857e-19 2.222633e-36 7.662580e-27 6.934620e-29 1.167551e-33 3.281019e-15 2.230434e-12 1.255132e-13
    #  [9] 1.831379e-28 3.471063e-22 4.148845e-33 6.991152e-25 1.564250e-16 1.415299e-17 2.583114e-37 3.994026e-14

    V(g)$is_CHD <- V(g)$name %in% CHD

    (table(E(g)$corexp_sign))
    # negative positive
    #        1       19

    # Remove nodes with no connections
    g_connected <- delete_vertices(g, which(degree(g) == 0))

    set.seed(1234) # Use to ensure layout repeatable
    layout_coords <- create_layout(g_connected, layout = "fr") # "fr" = Fruchterman-Reingold layout

    p_listoflist[[int]] <- ggraph(layout_coords) +
        geom_edge_link(aes(
            width = log10(weight), # Edge width = combined score
            color = corexp_sign # Edge color = pos/neg correlation
        ), alpha = 0.7) +

        geom_node_point(aes(
            size = weight, #  3 ,
            color = is_CHD # logFC_ave +1,
        )) +
        geom_node_text(aes(label = name), repel = TRUE, size = 3) + # <--- Add labels
        scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "gray70")) +
        scale_edge_color_manual(values = c("positive" = "orange", "negative" = "blue")) +
        scale_size_continuous(range = c(2, 5), name = "|Wilcox score|") +
        scale_edge_width_continuous(range = c(0.5, 3), name = "log10(E weights)") +

        theme_void() +
        labs(
            edge_color = "Specific co-exp",
            color = "Curated CHD genes"
        ) +
        theme(legend.position = "right") +
        ggtitle(paste0(db, ": ", int, " ", vcount(g_connected), "/", vcount(g), " PPI genes"))
}

(names(p_listoflist))
# "HiGCTS_endothelial.b" "HiGCTS_cardiac.a" 


pdf(file = "network_view_PPI_IbarraSoria2018.pdf", width = 12)
grid.arrange(p_listoflist[[1]], p_listoflist[[2]], ncol = 2)
dev.off()

############################################
#### Ideally, I want to plot network for each tipping point, each page is score 200, 400, 700, weighted in row
##  clusters in columns which are the HiG for parent, CT, and decenders
############################################

####################################
## however, the trajectory is more than lineage
## so, first plot HiG across all clusters, one score threshold per page
## non-cardiac and cardiac seperately
##############################################

score_threshold <- "weight_shrink"
input_path <- paste0(wd, "results/")

# refer to 11.1_CTS_cardiac_network_degreeDistribution.R
if (grepl("score", score_threshold)) {
    graph_list <- readRDS(file = paste0(input_path, "2018_STRING_graph_perState_notsimplified.rds"))
    graph_list <- lapply(graph_list, simplify) # !!!!!!!!!!!!!!!!!!!
} else {
    graph_list <- readRDS(file = paste0(input_path, "PPI_weight/2018_STRING_graph_perState_simplified_", s, "weighted.rds"))
}

p_list_HiG <- list()
for (int in grep("^HiG_", names(graph_list), value = TRUE)) {
    g <- graph_list[[int]]
    (V(g)$FDR)
    # [1] 1.267287e-121 7.654658e-103  7.693735e-51  3.771167e-51 9.032465e-132  6.385503e-99  3.121662e-88  1.649808e-33
    V(g)$is_CHD <- V(g)$name %in% CHD

    # Remove nodes with no connections
    g_connected <- delete_vertices(g, which(degree(g) == 0))

    set.seed(1234) # Use to ensure layout repeatable
    # layout_coords <- create_layout(g_connected, layout = "fr")  # "fr" = Fruchterman-Reingold layout
    layout_coords <- create_layout(g_connected, layout = "kk") # 'kk' = Kamada-Kawai (also stochastic, needs seed)

    p_list_HiG[[int]] <- ggraph(layout_coords) +
        geom_edge_link(aes(
            width = log10(weight), # Edge width = combined score
            color = corexp_sign # Edge color = pos/neg correlation
        ), alpha = 0.7) +

        geom_node_point(aes(
            size = weight, #  3 ,
            color = is_CHD # logFC_ave +1,
        )) +
        geom_node_text(aes(label = name), repel = TRUE, size = 3) + # <--- Add labels
        scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "gray70")) +
        scale_edge_color_manual(values = c("positive" = "orange", "negative" = "blue")) +
        scale_size_continuous(range = c(2, 5), name = "|Wilcox score|") +
        scale_edge_width_continuous(range = c(0.5, 3), name = "log10(E weights)") +

        theme_void() +
        labs(
            edge_color = "Specific co-exp",
            color = "Curated CHD genes"
        ) +
        theme(legend.position = "right") +
        ggtitle(paste0(db, ": ", int, " ", vcount(g_connected), "/", vcount(g), " PPI genes"))
    (p_list_HiG[[int]])
}


(names(p_list_HiG))
# [1] "HiGCTS_endothelial.b" "HiGCTS_cardiac.a"    
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [5] "HiG_blood"                  "HiG_mesodermProgenitors"    "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"  
#  [9] "HiG_somiticMesoderm"        "HiG_mixedMesoderm.a"        "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"              "HiG_endothelial.b"          "HiG_cardiac.a" 

# Cardiac-related HiG clusters
CT_id <- c(
  "cardiac.b", "cardiac.c", "cardiac.a"
)

# Non-cardiac HiG clusters
noncardiac_id <- c(
  "extraembryonicMesoderm", "mesodermProgenitors", "presomiticMesoderm.b", 
  "presomiticMesoderm.a", "somiticMesoderm", "mixedMesoderm.a", 
  "pharyngealMesoderm", "mixedMesoderm.b", "blood"
)

# Endothelial or other clusters
other_id <- c(
  "endothelial.a", "endothelial.c", "endothelial.d", "endothelial.b"
)

pdf(file = paste0("network_HiG_view_", score_threshold, "_IbarraSoria2018.pdf"), width = 20, height = 20)
print(grid.arrange(grobs = p_list_HiG[paste0("HiG_", CT_id)], ncol = 4))
print(grid.arrange(grobs = p_list_HiG[paste0("HiG_", other_id)], ncol = 4))
print(grid.arrange(grobs = p_list_HiG[paste0("HiG_", noncardiac_id)], ncol = 4))
dev.off()
