library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(patchwork)
library(igraph)

########## BEGINNING OF USER INPUT ##########

wd = "/Users/felixyu/Documents/IbarraSoria2018/"
setwd(paste0(wd, "results/"))
PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

db <- "IbarraSoria2018"

CT_id <- c("endothelial.b", "cardiac.a")  # critical transition clusters
CT_id_formatted <- paste0("(_", CT_id, ")") %>% paste(collapse="|")

########## END OF USER INPUT ##########

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
(which(graphs_with_duplicates))
# named integer(0)


N <- sapply(graph_list, vcount)
((N0 - N)[which(N0 - N > 0)])
# named numeric(0)

graph_list <- lapply(graph_list, simplify, edge.attr.comb ='max') # !!!!!!!!!!!!!!!!!!! # FIXED
N2 <- sapply(graph_list, vcount)
(all(N == N2)) # [1] TRUE


(names(graph_list))
#  [1] "HiG_blood"                  "HiG_cardiac.b"             
#  [3] "HiG_cardiac.c"              "HiG_endothelial.a"         
#  [5] "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [7] "HiG_extraembryonicMesoderm" "HiG_mesodermProgenitors"   
#  [9] "HiG_mixedMesoderm.a"        "HiG_mixedMesoderm.b"       
# [11] "HiG_pharyngealMesoderm"     "HiG_presomiticMesoderm.a"  
# [13] "HiG_presomiticMesoderm.b"   "HiG_somiticMesoderm"       
# [15] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"          
# [17] "CTS_endothelial.b"          "CTS_cardiac.a"  
edge_counts <- sapply(graph_list, ecount)
(edge_counts)
#                  HiG_blood              HiG_cardiac.b 
#                      10986                       7164 
#              HiG_cardiac.c          HiG_endothelial.a 
#                       8853                       9643 
#          HiG_endothelial.c          HiG_endothelial.d 
#                       8915                       7825 
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors 
#                       5666                       4991 
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b 
#                       4527                       5930 
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a 
#                       6110                       4881 
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm 
#                       7797                       4687 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a 
#                         20                         28 
#          CTS_endothelial.b              CTS_cardiac.a 
#                         82                         79 

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

(table(V_deg_dis$PPI_cat))
#    CTS HiGCTS    HiG 
#     30     18   3034
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
(n_nodes)
#                  HiG_blood              HiG_cardiac.b 
#                        420                        409 
#              HiG_cardiac.c          HiG_endothelial.a 
#                        444                        455 
#          HiG_endothelial.c          HiG_endothelial.d 
#                        480                        426 
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors 
#                        367                        356 
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b 
#                        329                        369 
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a 
#                        371                        338 
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm 
#                        412                        327 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a 
#                         16                         13 
#          CTS_endothelial.b              CTS_cardiac.a 
#                         33                         37 


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
#    CTS HiGCTS    HiG 
#     26     15   1393 
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
(dim(V_deg_dis)) # [1] 3082    7
(head(V_deg_dis, 3))
#    signature degree_distribution     k PPI_cat cluster n_nodes
#       <char>               <num> <int>  <fctr>  <char>   <num>
# 1: HiG_blood           1.0000000     1     HiG   blood     420
# 2: HiG_blood           0.9928571     2     HiG   blood     420
# 3: HiG_blood           0.9857143     3     HiG   blood     420
#    normalized_degree_distribution
#                             <num>
# 1:                    0.002386635
# 2:                    0.002369587
# 3:                    0.002352540



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

(table(V_deg_dis$PPI_cat))
#    CTS HiGCTS    HiG 
#     30     18   3034

(table(tmp$PPI_cat))
#    CTS HiGCTS    HiG 
#     30     18      0 
saveRDS(V_deg_dis, file = "V_deg_dis.rds") ## note that 'PPI_cast' was originally named as 'width'
