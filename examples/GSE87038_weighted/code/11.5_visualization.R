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
file <- paste0("GSE87038_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"       "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"      "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"
# [15] "HiG_11"      "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"
# [29] "CTS_16"      "CTS_16.1"    "CTS_13"      "CTS_8"       "HiGCTS_13"

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

# HiGCTS_13 too small, not included here
for (int in c("HiGCTS_7", "HiGCTS_8", "HiGCTS_8", "HiGCTS_11", "HiGCTS_15", "HiGCTS_16", "HiGCTS_16.1")) {
    g <- graph_list[[int]]
    (vertex_attr_names(g)) #  "name"   "weight" "FDR"
    (edge_attr_names(g)) # "weight"         "norm_PPI_score" "corexp_sign"    "coexp_target"

    (V(g)$FDR)
    # [1] 4.493309e-27 7.710597e-13 1.640461e-20 4.920997e-41 4.551240e-31
    # [6] 2.721066e-26 4.574606e-37 1.505239e-25 9.554066e-41 2.646641e-53
    # [11] 3.793807e-10 3.130554e-14 4.711156e-31 6.517444e-25
    V(g)$is_CHD <- V(g)$name %in% CHD

    (table(E(g)$corexp_sign))
    # negative positive
    #        3       7

    # Remove nodes with no connections
    g_connected <- delete_vertices(g, which(degree(g) == 0))

    set.seed(1234) # Use to ensure layout repeatable
    # layout_coords <- create_layout(g_connected, layout = "fr")  # "fr" = Fruchterman-Reingold layout
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
        scale_edge_width_continuous(range = c(0.1, 3), limits = range(log10(E(g)$weight), na.rm = TRUE), name = "log10(E weights)") +

        theme_void() +
        labs(
            edge_color = "Specific co-exp",
            color = "Curated CHD genes"
        ) +
        theme(legend.position = "right") +
        ggtitle(paste0(db, ": ", int, " ", vcount(g_connected), "/", vcount(g), " PPI genes"))
}

(names(p_listoflist))
# [1] "HiGCTS_7"    "HiGCTS_8"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1"

pdf(file = "network_view_PPI_GSE87038.pdf", width = 12, height = 12)
print(grid.arrange(grobs = list(p_listoflist[[1]], p_listoflist[[2]]), ncol = 2))
print(grid.arrange(grobs = list(p_listoflist[[3]], p_listoflist[[4]]), ncol = 2))
print(grid.arrange(grobs = list(p_listoflist[[5]], p_listoflist[[6]]), ncol = 2))
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
    g <- graph_list[[int]]
    (V(g)$FDR)
    # [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    V(g)$is_CHD <- V(g)$name %in% CHD

    # Remove nodes with no connections
    g_connected <- delete_vertices(g, which(degree(g) == 0))

    set.seed(1234) # Use to ensure layout repeatable
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
        scale_edge_width_continuous(range = c(0.1, 3), limits = range(log10(E(g)$weight), na.rm = TRUE), name = "log10(E weights)") +

        theme_void() +
        labs(
            edge_color = "Specific co-exp",
            color = "Curated CHD genes"
        ) +
        theme(legend.position = "right") +
        ggtitle(paste0(db, ": ", int, " ", vcount(g_connected), "/", vcount(g), " PPI genes"))
    print(p_list_HiG[[int]])
}


(names(p_list_HiG))
# [1] "HiG_1"  "HiG_2"  "HiG_3"  "HiG_4"  "HiG_5"  "HiG_6"  "HiG_9"  "HiG_10" "HiG_12" "HiG_14" "HiG_17" "HiG_18" "HiG_19" "HiG_7"  "HiG_11" "HiG_15" "HiG_16" "HiG_13" "HiG_8"

CT_id <- c(2, 4, 5, 8, 9, 17, 18, 19) # Cardiac-related clusters
noncardiac_id <- c(6, 7, 10, 11, 13, 14, 15) # Non-cardiac clusters
other_id <- setdiff(1:19, c(CT_id, noncardiac_id)) # Remaining clusters: 1, 3, 12, 16


pdf(file = paste0("network_HiG_view_", score_threshold, "_GSE87038.pdf"), width = 20, height = 20)
print(grid.arrange(grobs = p_list_HiG[CT_id], ncol = 4))
print(grid.arrange(grobs = p_list_HiG[other_id], ncol = 4))
print(grid.arrange(grobs = p_list_HiG[noncardiac_id], ncol = 4))

dev.off()
