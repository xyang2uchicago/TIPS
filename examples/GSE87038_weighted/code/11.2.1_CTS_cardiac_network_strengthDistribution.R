library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(igraph)

########## BEGINNING OF USER INPUT ##########

wd = "/Users/felixyu/Documents/GSE87038_weighted/"
setwd(paste0(wd, "results/"))
PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_palette <- c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)

db <- "GSE87038"

celltype_specific_weight_version <- '10'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

CT_id <- c("7", "8", "11", "13", "15", "16") # critical transition clusters
CT_id_formatted <- paste0("(_", CT_id, ")") %>% paste(collapse="|")

########## END OF USER INPUT ##########


graph_list <- readRDS(file = paste0(db, "_STRING_graph_perState_notsimplified.rds"))
(N0 <- sapply(graph_list, vcount))
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6
#         303         435         411         304         322         512
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#         358         422         341         364         524         457
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#         529         336         441         406         524         403
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#         332          14          19          30          13          30
#   HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#          13          10          31          51          66          39
#    CTS_16.1      CTS_13       CTS_8
#          79          60          54

########## remove name-duplicated Vertex due to inconsistence in STRING.db ###########
## refer to 11.1.0_correct_vertex_duplication.R
correct_n_edges <- readRDS("correct_n_edges_HiG_STRING2.14.0.rds")
for (g_name in unique(correct_n_edges$graph_id)) {
    # Subset the rows for this graph
    rows <- subset(correct_n_edges, graph_id == g_name)

    # Extract the vertex_index_to_remove as a character vector
    vertices_str <- rows$vertex_index_to_remove
    vertices_str <- vertices_str[!is.na(vertices_str)]

    # Split each string by "," and convert to numeric, then combine all
    vertices_to_remove <- unlist(
        lapply(vertices_str, function(s) as.numeric(strsplit(s, ",")[[1]]))
    )

    # Remove duplicates and sort decreasingly to avoid index shifting when deleting
    vertices_to_remove <- sort(unique(vertices_to_remove), decreasing = TRUE)
    if (length(vertices_to_remove) > 0) {
        graph_list[[g_name]] <- delete_vertices(graph_list[[g_name]], vertices_to_remove)
    }
}
N <- sapply(graph_list, vcount)
((N0 - N)[which(N0 - N > 0)])
#  HiG_5 HiG_17 HiG_18 HiG_19
#      2      1      2      2
graph_list <- lapply(graph_list, simplify, edge.attr.comb ='max') # !!!!!!!!!!!!!!!!!!! # FIXED
N2 <- sapply(graph_list, vcount)
(all(N == N2)) # [1] TRUE


(names(graph_list))
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_13"
# [26] "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"
# [31] "CTS_16.1"    "CTS_13"      "CTS_8"
edge_counts <- sapply(graph_list, ecount)
(edge_counts)
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6
#        4978        9544        7120        5362        5811       10993
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#        6937        8655        5828        7004       11724        8549
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#        7834        7572        9301        8632       10581        8508
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#        5880          17          15          54          10          31
#   HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#           9          10          52          54         207          70
#    CTS_16.1      CTS_13       CTS_8
#         210         147         165

graphs_with_duplicates <- sapply(graph_list, function(g) {
    vertex_names <- V(g)$name
    if (is.null(vertex_names)) {
        # If no names, use vertex indices
        vertex_names <- V(g)
    }
    any(duplicated(vertex_names))
})

# See which graphs have duplicates
(which(graphs_with_duplicates))
# named integer(0)

V_deg <- lapply(graph_list, function(x) {
    {
        igraph::strength(x, weights = E(x)$weight) / (vcount(x) - 1)
    } %>% sort(., decreasing = T)
}) %>%
    lapply(., function(x) {
        x %>%
            as.data.frame(strength = x) %>%
            mutate(gene = names(x), id = seq_along(x))
    }) %>%
    rbindlist(., idcol = names(.))
colnames(V_deg)[1:2] <- c("signature", "nor_strength")
V_deg$PPI_cat <- lapply(names(V_deg$signature), function(x) unlist(strsplit(x, split = "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))

# Aggregate strength Calculation:  this calculates the average strength for each signature (graph)
df <- aggregate(V_deg$nor_strength, by = list(V_deg$signature), FUN = mean) %>%
    mutate(k = aggregate(V_deg$id, by = list(V_deg$signature), FUN = max)[, 2]) %>%
    arrange(desc(x))
df$PPI_cat <- lapply(df$Group.1, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))

g_strength <- ggplot(data = df, aes(x = k, y = x, col = PPI_cat)) +
    scale_color_manual(values = PPI_color_palette) +
    geom_point(shape = 18, size = 5) +
    xlab("number of nodes per PPI_cat") +
    theme(legend.position = c(1, 1), legend.justification = c(0, 1)) +
    ylab("average GRN normalized strength") +
    ggtitle(db)



# cumulative (normalized= FALSE!!)strength distribution to a power law fit ########################
# PPI_cat: PPI network category: A: CTS; B: CTS&hiG; C; HiG  ############
#### normalized = FALSE by default

V_deg_dis <- lapply(graph_list, function(x) strength_distribution(x, normalized = FALSE, cumulative = TRUE)) %>%
    lapply(., function(x) {
        x %>%
            as.data.frame(strength_distribution = x) %>%
            mutate(k = seq_along(x))
    }) %>%
    rbindlist(., idcol = names(.))

colnames(V_deg_dis)[1:2] <- c("signature", "strength_distribution")
V_deg_dis$PPI_cat <- lapply(V_deg_dis$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))

(table(V_deg_dis$PPI_cat))
#    CTS HiGCTS    HiG 
#     69     30   1953 
V_deg_dis$cluster <- lapply(V_deg_dis$signature, function(x) unlist(strsplit(x, "_"))[2]) %>% unlist()

all(V_deg_dis$signature %in% names(graph_list))
V_deg_dis$n_nodes <- 0
for (i in seq_along(graph_list)) {
    j <- which(V_deg_dis$signature == names(graph_list)[i])
    V_deg_dis$n_nodes[j] <- vcount(graph_list[[i]])
}

## To provide insights into the distribution of node strengths in each signature and how the strength distribution varies across "transitory" and "steady" PPI_cats,
## we plot cumulative normalized strength distribution on log scale.
# showing how many (the cumulative fraction of) nodes in each signature  having a strength greater than or equal to k (the strength),
# with lines shaped based on the signature's PPI_cat ("transitory" or "steady")
g_strength_dis <- ggplot(
    data = V_deg_dis %>% filter(strength_distribution > 0),
    aes(x = k, y = strength_distribution, color = cluster, type = PPI_cat, size = PPI_cat)
) + #
    geom_line(aes(linetype = PPI_cat)) +
    xlab("cumulative  strength distribution") +
    scale_size_manual(values = PPI_size_palette) + # Set line width
    geom_text( # data=V_deg_dis_text,
        aes(label = n_nodes, color = cluster), # interaction(PPI_cat, cluster)),
        hjust = 1.1, vjust = 0.5, check_overlap = TRUE, size = 3
    ) + # Adding text for n_nodes
    theme(legend.position = c(0.2, 0.75), legend.justification = c(1, 1), legend.text = element_text(size = 5)) +
    coord_trans(x = "log10", y = "log10")

ggsave("strength_distribution_w_vsize.pdf", g_strength_dis, width = 11, height = 10)
(n_nodes <- lapply(graph_list, vcount) %>% unlist())
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6 
#         303         435         411         304         320         512 
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18 
#         358         422         341         364         523         455 
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13 
#         527         336         441         406         524         403 
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1 
#         332          14          19          30          13          30 
#   HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16 
#          13          10          31          51          66          39 
#    CTS_16.1      CTS_13       CTS_8 
#          79          60          54 

V_deg_nor_dis <- lapply(graph_list, function(g) {
    deg <- strength(g, weights = E(g)$weight)
    deg_table <- table(deg)
    df <- data.frame(
        k = as.integer(names(deg_table)), # the strength values
        freq = as.numeric(deg_table) # how many vertices have each strength value
    )
    df$nor_strength <- df$freq / sum(df$freq) # normalized frequency (probability); sum(df$freq)= vcount(g)
    df$nor_strength_cum <- rev(cumsum(rev(df$nor_strength))) # cumulative probability distribution
    df
}) %>%
    data.table::rbindlist(idcol = "signature")

V_deg_nor_dis$PPI_cat <- lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))

(table(V_deg_nor_dis$PPI_cat))
#    CTS HiGCTS    HiG
#    312     96   7618
V_deg_nor_dis$cluster <- lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x, "_"))[2]) %>% unlist()

all(V_deg_nor_dis$signature %in% names(graph_list))
V_deg_nor_dis$n_nodes <- 0
for (i in seq_along(graph_list)) {
    j <- which(V_deg_nor_dis$signature == names(graph_list)[i])
    V_deg_nor_dis$n_nodes[j] <- vcount(graph_list[[i]])
}

## To provide insights into the distribution of node strengths in each signature and how the strength distribution varies across "transitory" and "steady" PPI_cats,
## we plot cumulative normalized strength distribution on log scale.
# showing how many (the cumulative fraction of) nodes in each signature  having a strength greater than or equal to k (the strength),
# with lines shaped based on the signature's PPI_cat ("transitory" or "steady")
g_strength_dis <- V_deg_nor_dis %>%
    filter(k > 0, nor_strength_cum > 0) %>% # !!!!!!!!! NEW !!!!!!
    ggplot(aes(x = k, y = nor_strength_cum, color = cluster, linetype = PPI_cat, size = PPI_cat)) +
    geom_line() +
    scale_size_manual(values = PPI_size_palette) + # Set line width
    xlab("Normalized strength level") +
    ylab("cumulative normalized  strength distribution") +
    geom_text(aes(label = n_nodes), hjust = 1.1, vjust = 0.5, check_overlap = TRUE, size = 3) +
    theme(
        legend.position = c(0.2, 0.75),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 5)
    ) +
    coord_trans(x = "log10", y = "log10")
print(g_strength_dis)

g_strength_dis2 <- V_deg_nor_dis %>%
    filter(k > 0, nor_strength_cum > 0) %>% # !!!!!!!!! NEW !!!!!!
    ggplot(aes(x = k, y = nor_strength_cum, color = PPI_cat, linetype = PPI_cat, size = PPI_cat)) +
    geom_line(aes(group = signature, linetype = PPI_cat)) +
    scale_color_manual(values = PPI_color_palette) +
    scale_size_manual(values = PPI_size_palette) + # Set line width
    xlab("Normalized strength level") +
    ylab("Fraction of nodes having a normalized strength ≥ x") +
    theme(
        legend.position = c(0.2, 0.75),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 5)
    ) +
    coord_trans(x = "log10", y = "log10")
print(g_strength_dis2)

pdf(file = "normalized_strength_distribution.pdf")
print(g_strength)
print(g_strength_dis)
print(g_strength_dis2)
dev.off()




g_strength_dis <- ggplot(
    data = V_deg_dis %>% filter(strength_distribution > 0),
    aes(x = k, y = strength_distribution, color = cluster, type = PPI_cat, size = PPI_cat)
) + #
    geom_line(aes(linetype = PPI_cat)) +
    xlab(" strength level (x)") +
    scale_size_manual(values = PPI_size_palette) + # Set line width
    ylab("Fraction of nodes having a strength ≥ x") +
    theme(legend.position = c(0.2, 0.75), legend.justification = c(1, 1), legend.text = element_text(size = 5)) +
    coord_trans(x = "log10", y = "log10")
print(g_strength_dis)

g_strength_dis2 <- ggplot(
    data = V_deg_dis %>% filter(strength_distribution > 0),
    aes(x = k, y = strength_distribution, color = PPI_cat, size = PPI_cat)
) + #
    geom_line(aes(group = signature, linetype = PPI_cat)) +
    xlab(" strength level (x)") +
    scale_color_manual(values = PPI_color_palette) +
    scale_size_manual(values = PPI_size_palette) + # Set line width
    ylab("Fraction of nodes having a strength ≥ x") +
    theme(legend.position = c(0.2, 0.75), legend.justification = c(1, 1), legend.text = element_text(size = 5)) +
    coord_trans(x = "log10", y = "log10")
print(g_strength_dis2)

pdf(file = paste0("strength_", db, ".pdf"))
print(g_strength)
print(g_strength_dis2)
print(g_strength_dis)
dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

######## Compare node-strength distribution across three categories ################
print(dim(V_deg_dis)) # [1] 2052    6
print(head(V_deg_dis, 3))
#    signature strength_distribution     k PPI_cat cluster n_nodes
#       <char>                 <num> <int>  <fctr>  <char>   <num>
# 1:     HiG_1             1.0000000     1     HiG       1     303
# 2:     HiG_1             0.9669967     2     HiG       1     303
# 3:     HiG_1             0.9141914     3     HiG       1     303


###################################################
## 4.2) evaluate the strength distribution of each PPI_cat; NORMALIZED the cumulative strength is WRONG!!!!!!

V_deg_dis$normalized_strength_distribution <- V_deg_dis$strength_distribution / (V_deg_dis$n_nodes - 1)

# Density Plot / Kernel Density Estimate (KDE)
ggplot(V_deg_dis, aes(x = normalized_strength_distribution, fill = PPI_cat)) +
    geom_density(alpha = 0.5) + # Create density plot
    labs(x = "Normalized strength Distribution", y = "Density", title = "Density Plot: strength Distribution by Width Category") +
    theme_minimal() +
    scale_fill_manual(values = PPI_color_palette)
# Boxplot by Categories (NOT USED)
g1 <- ggplot(V_deg_dis, aes(x = factor(PPI_cat), y = normalized_strength_distribution, fill = PPI_cat)) +
    geom_boxplot() +
    labs(x = "PPIN Category", y = "Normalized strength Distribution", title = "PPINs for all clusters") +
    theme_minimal() +
    scale_fill_manual(values = PPI_color_palette) +
    stat_compare_means(
        method = "wilcox",
        comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")),
        p.adjust.method = "BH", # Adjust p-values using Benjamini-Hochberg (BH) method
        label = "p.signif"
    )
# Boxplot by Categories (NEW,  USED)
g2 <- ggplot(
    subset(V_deg_dis, grepl(CT_id_formatted, signature)),
    aes(x = factor(PPI_cat), y = normalized_strength_distribution, fill = PPI_cat)
) +
    geom_boxplot() +
    labs(x = "PPIN Category", y = "Normalized strength Distribution", title = "PPINs for transition clusters") +
    theme_minimal() +
    scale_fill_manual(values = PPI_color_palette) +
    stat_compare_means(
        method = "wilcox",
        comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")),
        p.adjust.method = "BH", # Adjust p-values using Benjamini-Hochberg (BH) method
        label = "p.signif"
    )
library(gridExtra)
pdf(file = paste0("boxplot_normalized_strength_", db, ".pdf"), height = 4)
print(grid.arrange(g1, g2, ncol = 2))
dev.off()

tmp <- subset(V_deg_dis, grepl(CT_id_formatted, signature))

print(table(V_deg_dis$PPI_cat))
#    CTS HiGCTS    HiG 
#     69     30   1953

print(table(tmp$PPI_cat))
#    CTS HiGCTS    HiG 
#     69     30    594
saveRDS(V_deg_dis, file = "V_deg_dis.rds") ## note that 'PPI_cast' was originally named as 'width'
