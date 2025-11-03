library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(patchwork)
library(igraph)

wd = "/Users/felixyu/Documents/GSE87038_weighted/"
setwd(paste0(wd, "results/"))
PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

db <- "GSE87038"

CT_id <- c(7, 8, 11, 13, 15, 16, 16.1)
CT_id_formatted <- paste0("_(", paste(CT_id, collapse = "|"), ")")

graph_list <- readRDS(file = paste0(db, "_STRING_graph_perState_notsimplified.rds"))
N0 <- sapply(graph_list, vcount)

# Check which graphs have duplicate vertex names
graphs_with_duplicates <- sapply(graph_list, function(g) {
    vertex_names <- V(g)$name
    if (is.null(vertex_names)) {
        # If no names, use vertex indices
        vertex_names <- V(g)
    }
    any(duplicated(vertex_names))
})

# See which graphs have duplicates
print(which(graphs_with_duplicates))
#  HiG_5 HiG_17 HiG_18 HiG_19
#      5     11     12     13

########## remove name-duplicated Vertex due to inconsistency in STRING.db ###########
## refer to 11.1.1_check_vertex_duplication.R
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

    cat("Removing from", g_name, ":", paste(vertices_to_remove, collapse = ", "), "\n")

    if (length(vertices_to_remove) > 0) {
        graph_list[[g_name]] <- delete_vertices(graph_list[[g_name]], vertices_to_remove)
    }
}
N <- sapply(graph_list, vcount)
print((N0 - N)[which(N0 - N > 0)])
#  HiG_5 HiG_17 HiG_18 HiG_19
#      2      1      2      2

graph_list <- lapply(graph_list, simplify, edge.attr.comb ='max') # !!!!!!!!!!!!!!!!!!! # FIXED
N2 <- sapply(graph_list, vcount)
all(N == N2) # [1] TRUE


names(graph_list)
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_13"
# [26] "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"
# [31] "CTS_16.1"    "CTS_13"      "CTS_8"
edge_counts <- sapply(graph_list, ecount)
edge_counts
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6
#        4978        9544        7120        5362        5808       10993
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#        6937        8655        5828        7004       11719        8543
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#        7827        7572        9301        8632       10581        8508
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#        5880          17          15          54          10          31
#   HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#           9          10          52          54         207          70
#    CTS_16.1      CTS_13       CTS_8
#         210         147         165

# Identify and print duplicate vertex indices for each graph, should return nothing if all duplicates are removed
for (g_name in names(graph_list)) {
    g <- graph_list[[g_name]]
    vertex_names <- V(g)$name

    if (is.null(vertex_names)) {
        message(g_name, ": No vertex names present.")
        next
    }

    duplicated_names <- vertex_names[duplicated(vertex_names)]

    if (length(duplicated_names) > 0) {
        cat("Graph:", g_name, "\n")
        dup_indices <- which(vertex_names %in% duplicated_names)
        print(dup_indices)
        cat("\n")
    }
}


V_deg <- lapply(graph_list, function(x) igraph::degree(x, normalized = T) %>% sort(., decreasing = T)) %>%
    lapply(., function(x) {
        x %>%
            as.data.frame(degree = x) %>%
            mutate(gene = names(x), id = seq_along(x))
    }) %>%
    rbindlist(., idcol = names(.))
colnames(V_deg)[1:2] <- c("signature", "nor_degree")
V_deg$PPI_cat <- lapply(names(V_deg$signature), function(x) unlist(strsplit(x, split = "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))

# Aggregate Degree Calculation:  this calculates the average degree for each signature (graph)
df <- aggregate(V_deg$nor_degree, by = list(V_deg$signature), FUN = mean) %>%
    mutate(k = aggregate(V_deg$id, by = list(V_deg$signature), FUN = max)[, 2]) %>%
    arrange(desc(x))
df$PPI_cat <- lapply(df$Group.1, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))

g_degree <- ggplot(data = df, aes(x = k, y = x, col = PPI_cat)) +
    scale_color_manual(values = PPI_color_palette) +
    geom_point(shape = 18, size = 5) +
    xlab("number of nodes per PPI_cat") +
    theme(legend.position = c(1, 1), legend.justification = c(0, 1)) +
    ylab("average GRN normalized degree") +
    ggtitle(db)



# cumulative (normalized= FALSE!!)degree distribution to a power law fit ########################
# PPI_cat: PPI network category: A: CTS; B: CTS&hiG; C; HiG  ############
#### normalized = FALSE by default
V_deg_dis <- lapply(graph_list, function(x) igraph::degree_distribution(x, cumulative = TRUE)) %>%
    lapply(., function(x) {
        x %>%
            as.data.frame(degree_distribution = x) %>%
            mutate(k = seq_along(x))
    }) %>%
    rbindlist(., idcol = names(.))

colnames(V_deg_dis)[1:2] <- c("signature", "degree_distribution")
V_deg_dis$PPI_cat <- lapply(V_deg_dis$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))

table(V_deg_dis$PPI_cat)
#    CTS HiGCTS    HiG
#    145     53   4370
V_deg_dis$cluster <- lapply(V_deg_dis$signature, function(x) unlist(strsplit(x, "_"))[2]) %>% unlist()

all(V_deg_dis$signature %in% names(graph_list))
V_deg_dis$n_nodes <- 0
for (i in seq_along(graph_list)) {
    j <- which(V_deg_dis$signature == names(graph_list)[i])
    V_deg_dis$n_nodes[j] <- vcount(graph_list[[i]])
}

## To provide insights into the distribution of node degrees in each signature and how the degree distribution varies across "transitory" and "steady" PPI_cats,
## we plot cumulative normalized degree distribution on log scale.
# showing how many (the cumulative fraction of) nodes in each signature  having a degree greater than or equal to k (the degree),
# with lines shaped based on the signature's PPI_cat ("transitory" or "steady")
g_degree_dis <- ggplot(
    data = V_deg_dis,
    aes(x = k, y = degree_distribution, color = cluster, type = PPI_cat)
) + #
    geom_line(aes(linetype = PPI_cat)) +
    xlab("cumulative degree distribution") +
    geom_text( # data=V_deg_dis_text,
        aes(label = n_nodes, color = cluster), # interaction(PPI_cat, cluster)),
        hjust = 1.1, vjust = 0.5, check_overlap = TRUE, size = 3
    ) + # Adding text for n_nodes
    theme(legend.position = c(0.2, 0.75), legend.justification = c(1, 1), legend.text = element_text(size = 5)) +
    coord_trans(x = "log10", y = "log10") #+ #xlab('degree k') +
# scale_x_reverse()  # Flip the x-axis from highest to smallest

print(g_degree_dis)
dev.copy2pdf(file = "degree_distribution_w_vsize.pdf")
(n_nodes <- lapply(graph_list, vcount) %>% unlist())
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6
#         303         435         411         304         320         512
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#         358         422         341         364         524         455
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#         527         336         441         406         524         403
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#         332          14          19          30          13          30
#   HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#          13          10          31          51          66          39
#    CTS_16.1      CTS_13       CTS_8
#          79          60          54


# cumulative NORMALIZED degree distribution to a power law fit ########################
# PPI_cat: PPI network category: A: CTS; B: CTS&hiG; C; HiG  ############
#### normalized = FALSE by default


V_deg_nor_dis <- lapply(graph_list, function(g) {
    deg <- igraph::degree(g)
    deg_table <- table(deg)
    df <- data.frame(
        k = as.integer(names(deg_table)),
        freq = as.numeric(deg_table)
    )
    df$nor_degree <- df$freq / sum(df$freq)
    df$nor_degree_cum <- rev(cumsum(rev(df$nor_degree)))
    df
}) %>%
    data.table::rbindlist(idcol = "signature")

# colnames(V_deg_nor_dis)[1:2]=c('signature','nor_degree_distribution')
V_deg_nor_dis$PPI_cat <- lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))

table(V_deg_nor_dis$PPI_cat)
# CTS HiGCTS    HiG
# 80     33   1033
V_deg_nor_dis$cluster <- lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x, "_"))[2]) %>% unlist()

all(V_deg_nor_dis$signature %in% names(graph_list))
V_deg_nor_dis$n_nodes <- 0
for (i in seq_along(graph_list)) {
    j <- which(V_deg_nor_dis$signature == names(graph_list)[i])
    V_deg_nor_dis$n_nodes[j] <- vcount(graph_list[[i]])
}

## To provide insights into the distribution of node degrees in each signature and how the degree distribution varies across "transitory" and "steady" PPI_cats,
## we plot cumulative normalized degree distribution on log scale.
# showing how many (the cumulative fraction of) nodes in each signature  having a degree greater than or equal to k (the degree),
# with lines shaped based on the signature's PPI_cat ("transitory" or "steady")
g_degree_dis <- V_deg_nor_dis %>%
    filter(k > 0, nor_degree_cum > 0) %>% # !!!!!!!!! NEW !!!!!!
    ggplot(aes(x = k, y = nor_degree_cum, color = cluster, linetype = PPI_cat)) +
    geom_line() +
    xlab("Normalized degree level") +
    ylab("cumulative normalized  degree distribution") +
    geom_text(aes(label = n_nodes), hjust = 1.1, vjust = 0.5, check_overlap = TRUE, size = 3) +
    theme(
        legend.position = c(0.2, 0.75),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 5)
    ) +
    coord_trans(x = "log10", y = "log10")

g_degree_dis2 <- V_deg_nor_dis %>%
    filter(k > 0, nor_degree_cum > 0) %>% # !!!!!!!!! NEW !!!!!!
    ggplot(aes(x = k, y = nor_degree_cum, color = PPI_cat, linetype = PPI_cat)) +
    geom_line(aes(group = signature, linetype = PPI_cat)) +
    scale_color_manual(values = PPI_color_palette) +
    xlab("Normalized degree level") +
    ylab("Fraction of nodes having a normalized degree ≥ x") +
    theme(
        legend.position = c(0.2, 0.75),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 5)
    ) +
    coord_trans(x = "log10", y = "log10")


pdf(file = "normalized_degree_distribution.pdf")
print(g_degree)
print(g_degree_dis)
print(g_degree_dis2)
dev.off()





g_degree_dis <- ggplot(
    data = V_deg_dis,
    aes(x = k, y = degree_distribution, color = cluster, type = PPI_cat)
) + #
    geom_line(aes(linetype = PPI_cat)) +
    xlab("degree level (x)") +
    ylab("Fraction of nodes having a degree ≥ x") +
    theme(legend.position = c(0.2, 0.75), legend.justification = c(1, 1), legend.text = element_text(size = 5)) +
    coord_trans(x = "log10", y = "log10")
g_degree_dis2 <- ggplot(
    data = V_deg_dis,
    aes(x = k, y = degree_distribution, color = PPI_cat)
) + #
    geom_line(aes(group = signature, linetype = PPI_cat)) +
    xlab("degree level (x)") +
    scale_color_manual(values = PPI_color_palette) +
    ylab("Fraction of nodes having a degree ≥ x") +
    theme(legend.position = c(0.2, 0.75), legend.justification = c(1, 1), legend.text = element_text(size = 5)) +
    coord_trans(x = "log10", y = "log10")

pdf(file = paste0("degree_", db, ".pdf"))
print(g_degree)
print(g_degree_dis2)
print(g_degree_dis)
dev.off()

######## Compare node-degree distribution across three categories ################
dim(V_deg_dis) # [1] 4568    6
head(V_deg_dis, 3)
#    signature degree_distribution     k PPI_cat cluster n_nodes
#       <char>               <num> <int>  <fctr>  <char>   <num>
# 1:     HiG_1           1.0000000     1     HiG       1     303
# 2:     HiG_1           0.9900990     2     HiG       1     303
# 3:     HiG_1           0.9735974     3     HiG       1     303



###################################################
## 4.2) evaluate the degree distribution of each PPI_cat; NORMALIZED the cumulative degree is WRONG!!!!!!

V_deg_dis$normalized_degree_distribution <- V_deg_dis$degree_distribution / (V_deg_dis$n_nodes - 1)

######## Compare node-degree distribution across three categories ################
dim(V_deg_dis) # [1] 4568   7
head(V_deg_dis, 3)
#    signature degree_distribution     k PPI_cat cluster n_nodes
#       <char>               <num> <int>  <fctr>  <char>   <num>
# 1:     HiG_1           1.0000000     1     HiG       1     303
# 2:     HiG_1           0.9900990     2     HiG       1     303
# 3:     HiG_1           0.9735974     3     HiG       1     303
#    normalized_degree_distribution
#                             <num>
# 1:                    0.003311258
# 2:                    0.003278474
# 3:                    0.003223832



# Density Plot / Kernel Density Estimate (KDE)
ggplot(V_deg_dis, aes(x = normalized_degree_distribution, fill = PPI_cat)) +
    geom_density(alpha = 0.5) + # Create density plot
    labs(x = "Normalized Degree Distribution", y = "Density", title = "Density Plot: Degree Distribution by Width Category") +
    theme_minimal() +
    scale_fill_manual(values = PPI_color_palette)
# Boxplot by Categories (NOT USED)
g1 <- ggplot(V_deg_dis, aes(x = factor(PPI_cat), y = normalized_degree_distribution, fill = PPI_cat)) +
    geom_boxplot() +
    labs(x = "PPIN Category", y = "Normalized Degree Distribution", title = "PPINs for all clusters") +
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
    aes(x = factor(PPI_cat), y = normalized_degree_distribution, fill = PPI_cat)
) +
    geom_boxplot() +
    labs(x = "PPIN Category", y = "Normalized Degree Distribution", title = "PPINs for transition clusters") +
    theme_minimal() +
    scale_fill_manual(values = PPI_color_palette) +
    stat_compare_means(
        method = "wilcox",
        comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")),
        p.adjust.method = "BH", # Adjust p-values using Benjamini-Hochberg (BH) method
        label = "p.signif"
    )


combined <- g1 + g2
ggsave(paste0("boxplot_normalized_degree_", db, ".pdf"), combined, width = 12, height = 4)

tmp <- subset(V_deg_dis, grepl(CT_id_formatted, signature))

table(V_deg_dis$PPI_cat)
#    CTS HiGCTS    HiG
#    145     53   4370

table(tmp$PPI_cat)
#    CTS HiGCTS    HiG
#    145     53   1321
saveRDS(V_deg_dis, file = "V_deg_dis.rds") ## note that 'PPI_cast' was originally named as 'width'
