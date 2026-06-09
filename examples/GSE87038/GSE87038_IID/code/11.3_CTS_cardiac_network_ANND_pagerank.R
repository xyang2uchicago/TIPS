library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(igraph)
library(rstatix)
library(brainGraph)

########## BEGINNING OF USER INPUT ##########

wd <- "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/GSE87038_IID/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- file.path(wd, "data/")

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

db <- "GSE87038"

s <- "combined" # specificity method

########## END OF USER INPUT ##########

file <- file.path(inputdir, paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds"))
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"
#  [7] "HiG_9"     "HiG_10"    "HiG_12"    "HiG_14"    "HiG_17"    "HiG_18"
# [13] "HiG_19"    "HiG_7"     "HiG_11"    "HiG_15"    "HiG_16"    "HiG_13"   
# [19] "HiG_8"     "HiGCTS_15" "CTS_7"     "CTS_11"    "CTS_15"    "CTS_16"
# [25] "CTS_16.1"  "CTS_8"

edge_counts <- sapply(graph_list, ecount)
(edge_counts)
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_9    HiG_10 
#      2254      3987      2959      2249      2056      4706      2758      3608
#    HiG_12    HiG_14    HiG_17    HiG_18    HiG_19     HiG_7    HiG_11    HiG_15
#      2430      2891      4503      3948      3399      2523      4158      3660
#    HiG_16    HiG_13     HiG_8 HiGCTS_15     CTS_7    CTS_11    CTS_15    CTS_16
#      5008      3490      2279        10        10        18        50        28 
#  CTS_16.1     CTS_8
#        57        65

(sapply(graph_list, vcount))
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_9    HiG_10
#       273       393       376       270       281       465       316       380
#    HiG_12    HiG_14    HiG_17    HiG_18    HiG_19     HiG_7    HiG_11    HiG_15
#       309       333       475       420       470       291       400       359 
#    HiG_16    HiG_13     HiG_8 HiGCTS_15     CTS_7    CTS_11    CTS_15    CTS_16
#       477       367       292        28        26        46        60        34
#  CTS_16.1     CTS_8
#        73        51

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
(which(graphs_with_duplicates)) # named integer(0)


#######################################################################
## fr. https://www.nature.com/articles/s41598-021-03625-w#Abs1
# 'Knn, page rank, and strength are the most relevant GRN features'
## First, evaluate the pageRanks,
## PageRank’s main difference from EigenCentrality is that it accounts for link direction.Like EigenCentrality,
## PageRank can help uncover influential or important nodes whose reach extends beyond just their direct connections.
## It’s especially useful in scenarios where link direction is important:
#
# https://igraph.org/r/doc/page_rank.html

page <- lapply(graph_list, function(x) page_rank(x, directed = FALSE, weights = E(x)$weight)$vector)

df <- lapply(page, function(x) data.frame(PageRank = x, gene = names(x)) %>% arrange(desc(PageRank))) %>%
    rbindlist(., idcol = names(.))
colnames(df)[1] <- "signature"
df$PPI_cat <- lapply(df$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))
(dim(df)[1]) # 7186

ic <- lapply(graph_list, function(x) eigen_centrality(x, weights = E(x)$weight)$vector)
IC <- lapply(ic, function(x) data.frame(EigenCentrality = x, gene = names(x)) %>% arrange(desc(EigenCentrality))) %>%
    rbindlist(., idcol = names(.))
colnames(IC)[1] <- "signature"
(dim(IC)) # [1] 7186    3
df <- merge(df, IC, by = c("signature", "gene"))
(dim(df)) # [1] 7186    5



## A) estimate the random PageRank by rewiring the edges while keeping the pro; this loop takes a while, Do NOT repeat !!
vn <- sapply(graph_list, vcount) # lengths(graph_list)
# Step 1: Calculate all pc values first
all_pc <- numeric(length(graph_list))
names(all_pc) <- names(graph_list)
for (i in names(graph_list)) {
    g <- graph_list[[i]]
    all_pc[i] <- mean(igraph::strength(g, weights = E(g)$weight)) / vn[i]
}

if (max(all_pc) > 1) { # ** v5 new
    # Step 2: Normalize to [0.01, 0.99] to preserve variability
    all_pc <- 0.01 + 0.98 * (all_pc - min(all_pc)) / (max(all_pc) - min(all_pc))
}

N <- 1000
set.seed(1234)
pr_P <- list()
for (i in names(graph_list)) {
    (cat(i))
    g <- graph_list[[i]]
    pc <- all_pc[i]
    g.random <- list()
    for (j in 1:N) {
        g.random[[j]] <- rewire(graph_list[[i]], each_edge(prob = pc)) ## rewiring the edges while keeping the pro
    }

    pr_random <- lapply(g.random, function(x) page_rank(x, directed = FALSE, weights = E(x)$weight)$vector)

    tmp <- lapply(pr_random, function(x) data.frame(PageRank = x, gene = names(x)) %>% arrange(desc(PageRank))) %>%
        rbindlist(., idcol = names(.))
    head(tmp)
#      PageRank   gene
#         <num> <char>
# 1: 0.08593170  EPHA2
# 2: 0.05744567   FGF8
# 3: 0.05480651  FGFR2
# 4: 0.05405320  NR2F1
# 5: 0.03860998   IRX5
# 6: 0.03855780   NPM3
    for (j in V(graph_list[[i]])$name) {
        pr_P[[i]][j] <- length(which(subset(tmp, gene == j)$PageRank >= page[[i]][j])) / N
    }
}
saveRDS(pr_P, file = "GSE87038_IID_PageRank_Pvalue_by_rewiring.rds")


pr_P <- readRDS(file = "GSE87038_IID_PageRank_Pvalue_by_rewiring.rds")
tmp <- lapply(pr_P, function(x) data.frame(p.PageRank = x, gene = names(x))) %>%
    rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"

df <- merge(df, tmp, by = c("signature", "gene")) %>%
    group_by(signature) %>%
    mutate(rank_by_p.PR = rank(p.PageRank)) %>%
    mutate(rank_by_PR = rank(-PageRank)) %>%
    ungroup()

(head(df))
#   signature gene     PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>     <chr>       <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 CTS_11    AARD      0.00676 CTS             0            1               31
# 2 CTS_11    ABHD4     0.00676 CTS             0            1               31
# 3 CTS_11    ASB4      0.00676 CTS             0            1               31
# 4 CTS_11    BAIAP2L1  0.0112  CTS             0.00682      0.999            8
# 5 CTS_11    BCAR3     0.0117  CTS             0.00750      0.999            8
# 6 CTS_11    CCDC80    0.00676 CTS             0            1               31
(dim(df)) # [1] 7265    8

(subset(df, tolower(gene) == "isl1"))
#   signature gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>     <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 CTS_8     ISL1  0.0219   CTS          0.00384         0.994          25
# 2 HiG_14    ISL1  0.000834 HiG          0.0000614       0.912          41 
# 3 HiG_2     ISL1  0.000400 HiG          0.00000174      0.828         298.
# 4 HiG_8     ISL1  0.000552 HiG          0.00000604      0.935         169

# number of significantly high pagerank per PPI_cats, too much control !
n.pr.high <- lapply(names(graph_list), function(x) nrow(subset(df, signature == x & p.PageRank < 0.05))) %>% unlist()
names(n.pr.high) <- names(graph_list)
(n.pr.high)
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_9    HiG_10
#         0         0         0         0         0         0         0         0
#    HiG_12    HiG_14    HiG_17    HiG_18    HiG_19     HiG_7    HiG_11    HiG_15
#         0         0         0         0         0         0         0         0
#    HiG_16    HiG_13     HiG_8 HiGCTS_15     CTS_7    CTS_11    CTS_15    CTS_16
#         0         0         0         0         0         0         0         0
#  CTS_16.1     CTS_8
#         0         0 

# write.table(df, file='df_PAGERANK.tsv',sep='\t', row.names=F)  #!!!!!!!!



## ANND (Average Nearest Neighbor strength) measures the strength-strength dependence adjacent to a vertex
##
## here, we evaluate the ann, which only works with simple graphs,
# is often used to characterize dependencies between strengths of a node and its neighbors in a network.
# a non-simple graph is to have multiple edges connecting two nodes or for there to be a self-edge.
# igraph::knn():
# res$knn:
# A numeric vector giving the average nearest neighbor strength for all vertices in the graph.
# res$knnk :
# Calculate the ANND (average nearest neighbor strength) of the given vertices and the same quantity in the function of vertex strength
# A numeric vector, its length is the maximum (total) vertex strength in the graph.
# The first element is the average nearest neighbor strength of vertices with strength one, etc.
# for zero strength vertices the answer in ‘knn’ is NaN
annd_observed <- list()
for (i in names(graph_list)) {
    G <- graph_list[[i]]
    # remove unconnected nodes
    V_Isolated <- which(igraph::degree(G) == 0)
    G <- delete_vertices(G, V_Isolated) # !!!!!!!!!
    annd_observed[[i]] <- knn(G, weights = E(G)$weight)$knn
}
# annd_observed = lapply(graph_list, function(x) knn(x, weights = E(x)$weight)$knn )   # ** update
(any(is.na(annd_observed[["CTS_8"]]))) # FALSE
rm(G)

## A) estimate the random annd by rewiring the edges while keeping the pro; this loop takes a while, Do NOT repeat !!
annd_P <- list()
N <- 1000
set.seed(1234)
pr_P <- list()
for (i in names(graph_list)) {
    (cat(i))
    g <- graph_list[[i]]
    pc <- all_pc[i]
    g.random <- list()
    for (j in 1:N) {
        g.random[[j]] <- rewire(graph_list[[i]], each_edge(prob = pc)) ## rewiring the edges while keeping the pro
        # cat(range(E(g.random[[1]])$weight)) # [1] 0.404 1.866
    }

    annd_random <- lapply(g.random, function(x) knn(x, weights = E(x)$weight)$knn) # ** update

    tmp <- lapply(annd_random, function(x) data.frame(annd = x, gene = names(x)) %>% arrange(desc(annd))) %>%
        rbindlist(., idcol = names(.))
    head(tmp, 3)
#        annd      gene
#       <num>    <char>
# 1: 13.00000     PTPRE
# 2: 10.65799 RAB11FIP1
# 3:  9.59835     EPHA7
    for (j in V(graph_list[[i]])$name) {
        annd_P[[i]][j] <- length(which(subset(tmp, gene == j)$knn >= annd_observed[[i]][j])) / N
    }
}
saveRDS(annd_P, file = "GSE87038_IID_annd_Pvalue_by_rewiring.rds")

annd_P <- readRDS(file = "GSE87038_IID_annd_Pvalue_by_rewiring.rds")

(unique(df$PPI_cat))
# [1] CTS    HiGCTS HiG
# Levels: CTS HiGCTS HiG

df <- rbind(
    subset(df, PPI_cat == "CTS"),
    subset(df, PPI_cat == "HiGCTS"),
    subset(df, PPI_cat == "HiG")
)
(unique(df$PPI_cat))
# [1] CTS    HiGCTS HiG
# Levels: CTS HiGCTS HiG

tmp <- lapply(annd_observed, function(x) data.frame(annd = x, gene = names(x)) %>% arrange(desc(annd))) %>%
    rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"
df <- merge(df, tmp, by = c("signature", "gene"))
(dim(df)) # [1] 6786    9

annd_P[["CTS_8"]]
tmp <- lapply(annd_P, function(x) data.frame(p.annd = x, gene = names(x))) %>%
    rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"
(dim(tmp)) # [1] 7265   3

df <- merge(df, tmp, by = c("signature", "gene"))
df[which(is.na(df$knn)), "p.annd"] <- NA ## due to nrow(df) 4878 > nrow(tmp)
(dim(df)) # [1] 6786   11

## merge back the normalized strength of vertex
# normalized_strength <- strength(g) / (vcount(g) - 1)
# refer to 11.1_CTS_cardiac_network_strengthDistribution.R
V_strength <- lapply(graph_list, function(g) {
    # Calculate strength and sort
    strength <- strength(g, weights = E(g)$weight)
    strength_sorted <- sort(strength, decreasing = TRUE)
    # Create data frame
    data.frame(
        strength = strength_sorted,
        gene = names(strength_sorted),
        id = seq_along(strength_sorted)
    )
}) %>%
    rbindlist(., idcol = "signature") %>%
    dplyr::rename("rank_by_strength" = "id")

V_strength_norm <- lapply(graph_list, function(g) {
    # Calculate normalized strength and sort
    norm_strength <- strength(g, weights = E(g)$weight) / (vcount(g) - 1)
    norm_strength_sorted <- sort(norm_strength, decreasing = TRUE)
    # Create data frame
    data.frame(
        normalized.strength = norm_strength_sorted,
        gene = names(norm_strength_sorted),
        id = seq_along(norm_strength_sorted)
    )
}) %>%
    rbindlist(., idcol = "signature") %>%
    dplyr::rename("rank_by_normalized.strength" = "id")

V_strength <- merge(V_strength, V_strength_norm, by = c("signature", "gene"))
(dim(V_strength)) # 7265    6

## add the V_strength & V_strength_norm infor
df <- merge(df, V_strength, by = c("signature", "gene"))

(dim(df)) # 6786   14
(head(df, 3)) 
#   signature     gene   PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
# 1    CTS_11 BAIAP2L1 0.01124110     CTS     0.006817923      0.999            8
# 2    CTS_11    BCAR3 0.01169290     CTS     0.007504827      0.999            8
# 3    CTS_11    DHRS7 0.01892042     CTS     0.073443392      0.999            8
#   rank_by_PR annd p.annd     strength rank_by_strength normalized.strength
# 1         18    4      0 0.0002066144               16        4.591432e-06
# 2         17    4      0 0.0002274308               15        5.054017e-06
# 3         16    3      0 0.0009656290               13        2.145842e-05
#   rank_by_normalized.strength
# 1                          16
# 2                          15
# 3                          13

(table(df$normalized.strength >= df$annd))
# FALSE
#  6786

## # Add rank_by_ANND column and rerank strength, normalized strength by considering ties !!!
df <- df %>%
    group_by(signature) %>% # Group by 'signature'
    mutate(rank_by_strength = rank(-strength, na.last = "keep")) %>% # highest to smallest
    # mutate(rank_by_normalized.strength = rank(-normalized.strength, na.last = "keep")) %>%
    mutate(rank_by_ANND = rank(-annd, na.last = "keep")) %>% # Rank the 'annd' values, ignoring NA values
    mutate(rank_by_PR = rank(-PageRank, na.last = "keep")) %>%
    mutate(rank_by_p.PR = rank(p.PageRank, na.last = "keep")) %>% # smallest to highest
    mutate(rank_by_p.ANND = rank(p.annd, na.last = "keep")) %>% # smallest to highest
    ungroup() # Ungroup after the operation
(colnames(df))
#IID PIPELINE RESULTS
#  [1] "signature"                   "gene"
#  [3] "PageRank"                    "PPI_cat"
#  [5] "EigenCentrality"             "p.PageRank"
#  [7] "rank_by_p.PR"                "rank_by_PR"
#  [9] "annd"                        "p.annd"
# [11] "strength"                    "rank_by_strength"
# [13] "normalized.strength"         "rank_by_normalized.strength"
# [15] "rank_by_ANND"                "rank_by_p.ANND"

saveRDS(df, file = "IID_df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
write.table(df, file = "IID_df_PAGERANK_strength_ANND.rewring.P.tsv", sep = "\t", row.names = F) # !!!!!!!!

##########################
## add the column of betweenness centrality
##########################
df <- readRDS(file = "IID_df_PAGERANK_strength_ANND.rewring.P.rds")

# igraph::betweenness() uses distance graph weights, but E(g) uses connection weights, thus we invert it.
betweenness_list <- lapply(graph_list, function(x) betweenness(x, weights = 1/E(x)$weight))
bc.median <- lapply(betweenness_list, median) %>% unlist()

for (i in seq_along(betweenness_list)) {
    betweenness_list[[i]] <- data.frame(BetweennessCentrality = betweenness_list[[i]], gene = names(betweenness_list[[i]])) %>%
        mutate(rank_by_BC = rank(-BetweennessCentrality, na.last = "keep")) # Rank the 'annd' values, ignoring NA values
}
df_BC <- betweenness_list %>% rbindlist(., idcol = names(.))
colnames(df_BC)[1] <- "signature"
df_BC$PPI_cat <- lapply(df_BC$signature, function(x) unlist(strsplit(x, split = "_"))[1]) %>% unlist()

(dim(df_BC)) # [1] 7265    5

write.table(df_BC, file = "df_betweeness.tsv", sep = "\t", row.names = F) # !!!!!!!!


########### betweenness centrality ############
df_BC <- read.table(file = "df_betweeness.tsv", sep = "\t", header = T)
df_BC$PPI_cat <- factor(df_BC$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

CHD <- readRDS(file.path(inputdir, "CHD_Cilia_Genelist.rds"))
df_BC$PCGC_AllCurated <- toupper(df_BC$gene) %in% toupper(unlist(CHD["Griffin2023_PCGC_AllCurated"]))

# Calculate top 5 significant genes within each box
df5 <- df_BC %>%
    filter(rank_by_BC <= 5 & BetweennessCentrality > 0) %>%
    ungroup()

write.table(df5[, c(
    "signature", "gene", "BetweennessCentrality", "PPI_cat", "rank_by_BC",
    "PCGC_AllCurated"
)], file = "table_top5_Betweenness_perPPI.tsv", sep = "\t", row.names = FALSE, quote = FALSE) # !!!!!!!!!!!!!!

(subset(df5, PPI_cat == "HiGCTS"))#
#     signature BetweennessCentrality   gene rank_by_BC PPI_cat PCGC_AllCurated
# 96  HiGCTS_15                    17   CDH5        1.0  HiGCTS           FALSE
# 97  HiGCTS_15                     6   RND2        3.5  HiGCTS           FALSE
# 98  HiGCTS_15                     1  CLDN5        5.0  HiGCTS           FALSE
# 99  HiGCTS_15                     6   GMFG        3.5  HiGCTS           FALSE
# 100 HiGCTS_15                    10 PTPN18        2.0  HiGCTS           FALSE

(df5_CHD <- subset(df5, PCGC_AllCurated == TRUE))

(df5_CHD)
#     signature BetweennessCentrality  gene rank_by_BC PPI_cat PCGC_AllCurated
# 9       HiG_2                 12502 ACTC1        4.0     HiG            TRUE
# 17      HiG_4                  5072 HDAC1        2.0     HiG            TRUE
# 21      HiG_5                  5198 HDAC1        3.0     HiG            TRUE
# 32      HiG_9                  6666 HDAC1        3.0     HiG            TRUE
# 53     HiG_17                  9734 ACTC1        5.0     HiG            TRUE
# 74     HiG_11                  6772 ACTC1        5.0     HiG            TRUE
# 106    CTS_11                    48 FGFR2        1.0     CTS            TRUE
# 125     CTS_8                   128 NR2F2        4.5     CTS            TRUE
# 126     CTS_8                   128 FGFR2        4.5     CTS            TRUE
# 127     CTS_8                   130  IRX5        3.0     CTS            TRUE



df <- df %>% filter(PPI_cat %in% c("CTS", "HiG")) %>% droplevels()

df_BC$signature <- factor(df_BC$signature, levels = unique(df_BC$signature))
pr <- ggplot(df_BC, aes(x = signature, y = log10(BetweennessCentrality + 1), colour = PPI_cat)) +
    geom_boxplot(position = "dodge2") +
    theme(
        legend.position = "none",
        legend.justification = c(1, 1), # Place legend at top-right corner
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    scale_color_manual(values = PPI_color_palette) +
    geom_text(data = df5_CHD, aes(label = gene), size = 2, hjust = -0.1, vjust = 0, check_overlap = TRUE) + # Adjust text labels
    labs(color = "PPI cat") # Optional: label for the color legend
pr_repel <- ggplot(df_BC, aes(x = signature, y = log10(BetweennessCentrality + 1), colour = PPI_cat)) +
    geom_boxplot(position = "dodge2") +
    theme(
        legend.position = "none", # c(1,1)
        legend.justification = c(1, 1), # Place legend at top-right corner
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    scale_color_manual(values = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")) +
    geom_text_repel(
        data = df5_CHD, # df5
        aes(label = gene),
        size = 2, # Adjust the size of the text labels
        box.padding = 0.5, # Add space between the text and the data points
        point.padding = 0.5, # Add space between the text and the points
        segment.color = "grey50", # Color for the line connecting the text to the points
        max.overlaps = 40, # Max number of overlaps before labels stop being placed
        show.legend = FALSE # Do not show text labels in the legend
    ) +
    # scale_x_discrete(limits = unique(df$signature)) +
    labs(color = "PPI cat") # Optional: label for the color legend


density_bc_plot <- ggplot(df_BC, aes(x = log10(BetweennessCentrality + 1), color = PPI_cat, fill = PPI_cat)) +
    geom_density(alpha = 0.3) + # Density lines with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
    coord_flip() + # Flip the axes to rotate the density plot
    labs(x = "Density of BetweennessCentrality", y = "") # Label the x-axis and remove the y-axis label
# Add statistical comparisons using stat_compare_means manyally

# Violin plot with statistical comparisons
violin_wilcox <- ggplot(df_BC, aes(x = PPI_cat, y = log10(BetweennessCentrality + 1), color = PPI_cat, fill = PPI_cat)) +
    geom_violin(alpha = 0.3) + # Violin plot with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "PPI category") + # Label the axes
    # Add statistical comparisons using stat_compare_means
    stat_compare_means(
        aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
        comparisons = list(c("HiG", "HiGCTS"), c("HiGCTS", "CTS"), c("HiG", "CTS")), # Specify comparisons
        method = "wilcox.test", # Non-parametric test (Wilcoxon)
        label = "p.signif", # Show significance labels (e.g., **, *, ns)
        label.x = 1.5, # Adjust x-position of the p-value text
        size = 4 # Adjust size of the p-value text
        , tip.length = 0
    ) +
    ggtitle("wilcox-test, all PPINs")

violin_t <- ggplot(df_BC, aes(x = PPI_cat, y = log10(BetweennessCentrality + 1), color = PPI_cat, fill = PPI_cat)) +
    geom_violin(alpha = 0.3) + # Violin plot with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "PPI category") + # Label the axes
    # Add statistical comparisons using stat_compare_means
    stat_compare_means(
        aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
        comparisons = list(c("HiG", "HiGCTS"), c("HiGCTS", "CTS"), c("HiG", "CTS")), # Specify comparisons
        method = "t.test", # Non-parametric test (Wilcoxon)
        label = "p.signif", # Show significance labels (e.g., **, *, ns)
        label.x = 1.5, # Adjust x-position of the p-value text
        size = 4 # Adjust size of the p-value text
        , tip.length = 0
    ) +
    ggtitle("t-test, all PPINs")



df_median <- df_BC %>%
    group_by(signature) %>%
    summarise(bc.median = median(BetweennessCentrality, na.rm = TRUE)) %>%
    as.data.frame()
df_median$PPI_cat <- lapply(df_median$signature %>% as.vector(), function(x) unlist(strsplit(x, split = "_"))[1]) %>% unlist()
df_median$PPI_cat <- factor(df_median$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

## it is more manke sense to access each signature by its median betweness rank !!!!
a <- grepl("^HiG_", df_median$signature)
b <- grepl("^HiGCTS_", df_median$signature)
c <- grepl("^CTS_", df_median$signature)
ks.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.1  HiG vs HiGCTS
ks.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 1.129e-05	HiG vs CTS
ks.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value =  1	HiGCTS vs CTS
wilcox.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.1183  HiG vs HiGCTS
wilcox.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 0.0002936	HiG vs CTS
wilcox.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = NA	HiGCTS vs CTS
# t.test(df_median$bc.median[a], df_median$bc.median[b]) # not enough HiGCTS
t.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value =  3.73e-05	HiG vs CTS
t.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = 3.73e-05	HiGCTS vs CTS


density_median_bc_plot <- ggplot(df_median, aes(x = log10(bc.median), color = PPI_cat, fill = PPI_cat)) +
    geom_density(alpha = 0.3) + # Density lines with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "Density of the median of BetweennessCentralitys per PPI", y = "")
# Add statistical comparisons using stat_compare_means manyally

x <- which(df_median$bc.median == 0)

violin_median_bc_wilcox <- ggplot(df_median, aes(x = PPI_cat, y = bc.median, color = PPI_cat, fill = PPI_cat)) +
    geom_violin(alpha = 0.3, drop = FALSE) + # Violin plot with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "PPI category", y = "median of BC per PPI") + # Label the axes
    # Add statistical comparisons using stat_compare_means
    stat_compare_means(
        aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
        comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
        method = "wilcox.test", # Non-parametric test (Wilcoxon)
        label = "p.signif", # Show significance labels (e.g., **, *, ns)
        label.x = 1.5, # Adjust x-position of the p-value text
        size = 4 # Adjust size of the p-value text
        , tip.length = 0
    ) +
    ylim(0, NA) + # Start from 0, let ggplot choose upper limit
    ggtitle("wilcox-test, median BC")

violin_median_bc_wilcox_ln <- ggplot(df_median, aes(x = PPI_cat, y = log10(bc.median + 1), color = PPI_cat, fill = PPI_cat)) +
    geom_violin(alpha = 0.3, drop = FALSE) + # Violin plot with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "PPI category", y = "log10(median of BC per PPI +1)") + # Label the axes
    # Add statistical comparisons using stat_compare_means
    stat_compare_means(
        aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
        comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
        method = "wilcox.test", # Non-parametric test (Wilcoxon)
        label = "p.signif", # Show significance labels (e.g., **, *, ns)
        label.x = 1.5, # Adjust x-position of the p-value text
        size = 4 # Adjust size of the p-value text
        , tip.length = 0
    ) +
    ylim(0, NA) + # Start from 0, let ggplot choose upper limit
    ggtitle("wilcox-test, median BC+1")

# Combine the boxplot and density plot
pdf(file = "BetweennessCentrality_GSE87038_v2.pdf", height = 10)
print(grid.arrange(pr_repel, density_median_bc_plot + coord_flip(), ncol = 2, widths = c(3, 1)))
print(grid.arrange(violin_median_bc_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_t, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_median_bc_wilcox, violin_median_bc_wilcox_ln, nrow = 2))
dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


########### plot PageRank ############
df <- readRDS(file = "IID_df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
(dim(df)) # [1] 6786   16

df <- rbind(
    subset(df, PPI_cat == "CTS"),
    subset(df, PPI_cat == "HiGCTS"),
    subset(df, PPI_cat == "HiG")
)

CHD <- readRDS(file.path(inputdir, "CHD_Cilia_Genelist.rds"))
df$PCGC_AllCurated <- toupper(df$gene) %in% toupper(unlist(CHD["Griffin2023_PCGC_AllCurated"]))

# Calculate top 5 significant genes within each box
df5 <- df %>%
    filter(rank_by_PR <= 5) %>%
    ungroup()
(subset(df5, PPI_cat == "HiGCTS"))
#   signature gene   PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>     <chr>     <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 HiGCTS_15 CDH5     0.168  HiGCTS                0          1          6.5
# 2 HiGCTS_15 CLDN5    0.108  HiGCTS                1          1          6.5
# 3 HiGCTS_15 GMFG     0.0907 HiGCTS                0          1          6.5
# 4 HiGCTS_15 PTPN18   0.102  HiGCTS                0          1          6.5
# 5 HiGCTS_15 RND2     0.0717 HiGCTS                0          1          6.5


write.table(df5[, c(
    "signature", "gene", "PageRank", "PPI_cat", "rank_by_PR",
    "normalized.strength", "rank_by_normalized.strength", "PCGC_AllCurated"
)], file = "table_top5_PageRank_perPPI.tsv", sep = "\t", row.names = FALSE, quote = FALSE) # !!!!!!!!!!!!!!

df5_CHD <- subset(df5, PCGC_AllCurated == TRUE)
(dim(df5)) # [1] 130  17
(dim(df5_CHD)) # [1] 8 17
(df5_CHD %>% as.data.frame())
#   signature  gene   PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
# 1    CTS_11 FGFR2 0.10066248     CTS      1.00000000      0.999          8.0
# 2     CTS_8 FGFR2 0.05480651     CTS      0.89790544      0.992          6.5
# 3     CTS_8  IRX5 0.03860998     CTS      0.01175613      0.996         37.0
# 4    HiG_11 HDAC1 0.01613842     HiG      1.00000000      0.901        248.0
# 5    HiG_12 HDAC1 0.01524682     HiG      0.01282003      0.934        283.5
# 6    HiG_19  FLNA 0.01112584     HiG      0.02054271      0.973        232.5
# 7     HiG_4 HDAC1 0.01701202     HiG      0.01740795      0.911         46.5
# 8     HiG_9 HDAC1 0.02178261     HiG      0.04833902      0.931        227.0
#   rank_by_PR      annd p.annd    strength rank_by_strength normalized.strength
# 1          1  2.291794      0 0.007387852                1        1.641745e-04
# 2          3  4.745451      0 0.014769619                3        2.953924e-04
# 3          5  3.763917      0 0.008777736                7        1.755547e-04
# 4          2 35.913179      0 0.081584708                2        2.044730e-04
# 5          3 37.695558      0 0.068916721                4        2.237556e-04
# 6          1 31.026887      0 0.024182650                7        5.156215e-05
# 7          4 37.705799      0 0.062556947                5        2.325537e-04
# 8          1 41.560844      0 0.086337173                1        2.740863e-04
#   rank_by_normalized.strength rank_by_ANND rank_by_p.ANND PCGC_AllCurated
# 1                           1           11            9.5            TRUE
# 2                           3           24           19.5            TRUE
# 3                           7           33           19.5            TRUE
# 4                           2          248          190.0            TRUE
# 5                           4           73          146.5            TRUE
# 6                           7          187          223.5            TRUE
# 7                           5          110          127.0            TRUE
# 8                           1          101          152.0            TRUE

df$signature <- factor(df$signature, levels = levels(df_BC$signature))
pr <- ggplot(df, aes(x = signature, y = PageRank, colour = PPI_cat)) +
    geom_boxplot(show.legend = TRUE) + # Enable legend for the boxplot
    scale_color_manual(values = PPI_color_palette) +
    geom_text(
        data = df5_CHD, aes(label = gene), # data=df5
        size = 2, # Adjust the size of the text labels
        hjust = -0.1, vjust = 0,
        check_overlap = TRUE
    ) + # Avoid text overlap
    theme(
        legend.position = "none", # c(1, 1),
        legend.justification = c(1, 1), # Place legend at top-right corner
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    scale_x_discrete(limits = unique(df$signature)) +
    labs(color = "PPI cat") # Optional: label for the color legend
pr_repel <- ggplot(df, aes(x = signature, y = PageRank, colour = PPI_cat)) +
    geom_boxplot(show.legend = TRUE) + # Enable legend for the boxplot
    scale_color_manual(values = PPI_color_palette) +
    geom_text_repel(
        data = df5_CHD, # df5
        aes(label = gene),
        size = 2, # Adjust the size of the text labels
        box.padding = 0.5, # Add space between the text and the data points
        point.padding = 0.5, # Add space between the text and the points
        segment.color = "grey50", # Color for the line connecting the text to the points
        max.overlaps = 20, # Max number of overlaps before labels stop being placed
        show.legend = FALSE # Do not show text labels in the legend
    ) +
    theme(
        legend.position = "none",
        # legend.justification = 'none', #c(1, 1),  # Place legend at top-right corner
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    # scale_x_discrete(limits = unique(df$signature)) +
    labs(color = "PPI cat") # Optional: label for the color legend


density_page_plot <- ggplot(df, aes(x = PageRank, color = PPI_cat, fill = PPI_cat)) +
    geom_density(alpha = 0.3) + # Density lines with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
    coord_flip() + # Flip the axes to rotate the density plot
    labs(x = "Density of PageRank", y = "") # Label the x-axis and remove the y-axis label
# Add statistical comparisons using stat_compare_means manyally
# Violin plot with statistical comparisons
violin_wilcox <- ggplot(df, aes(x = PPI_cat, y = PageRank, color = PPI_cat, fill = PPI_cat)) +
    geom_violin(alpha = 0.3) + # Violin plot with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "PPI category", y = "PageRank") + # Label the axes
    # Add statistical comparisons using stat_compare_means
    stat_compare_means(
        aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
        comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
        method = "wilcox.test", # Non-parametric test (Wilcoxon)
        label = "p.signif", # Show significance labels (e.g., **, *, ns)
        label.x = 1.5, # Adjust x-position of the p-value text
        size = 4 # Adjust size of the p-value text
        , tip.length = 0
    ) +
    ggtitle("wilcox test, all PPINs ")

violin_t <- ggplot(df, aes(x = PPI_cat, y = PageRank, color = PPI_cat, fill = PPI_cat)) +
    geom_violin(alpha = 0.3) + # Violin plot with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "PPI category", y = "PageRank") + # Label the axes
    # Add statistical comparisons using stat_compare_means
    stat_compare_means(
        aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
        comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
        method = "t.test", # Non-parametric test (Wilcoxon)
        label = "p.signif", # Show significance labels (e.g., **, *, ns)
        label.x = 1.5, # Adjust x-position of the p-value text
        size = 4 # Adjust size of the p-value text
        , tip.length = 0
    ) +
    ggtitle("t test, all PPINs ")

## it is more manke sense to access each signature by its median bc rank !!!!
# pg.median = lapply(page, median) %>% unlist
pg.median <- df %>%
    group_by(signature) %>%
    summarise(median_PageRank = median(PageRank, na.rm = TRUE))
a <- grepl("^HiG_", pg.median$signature)
b <- grepl("^HiGCTS_", pg.median$signature)
c <- grepl("^CTS_", pg.median$signature)
ks.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value =  0.1  HiG vs HiGCTS
ks.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  1.129e-05	HiG vs CTS
ks.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.2857	HiGCTS vs CTS
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value = 0.1  HiG vs HiGCTS
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  1.129e-05 	HiG vs CTS
wilcox.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value = 0.2857	HiGCTS vs CTS


df_median <- data.frame(
    pg.median = pg.median$median_PageRank,
    PPI_cat = lapply(pg.median$signature %>% as.vector(), function(x) unlist(strsplit(x, split = "_", fixed = T))[1]) %>% unlist()
)
df_median$PPI_cat <- factor(df_median$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))
density_median_page_plot <- ggplot(df_median, aes(x = pg.median, color = PPI_cat, fill = PPI_cat)) +
    geom_density(alpha = 0.3) + # Density lines with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "Density of the median of PageRanks per PPI", y = "")
# Add statistical comparisons using stat_compare_means manyally

violin_median_page_wilcox <- ggplot(df_median, aes(x = PPI_cat, y = pg.median, color = PPI_cat, fill = PPI_cat)) +
    geom_violin(alpha = 0.3) + # Violin plot with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "PPI category", y = "median of PageRanks per PPI") + # Label the axes
    # Add statistical comparisons using stat_compare_means
    stat_compare_means(
        aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
        comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
        method = "wilcox.test", # Non-parametric test (Wilcoxon)
        label = "p.signif", # Show significance labels (e.g., **, *, ns)
        label.x = 1.5, # Adjust x-position of the p-value text
        size = 4 # Adjust size of the p-value text
        , tip.length = 0
    ) +
    ggtitle("wilcox test, median PA")


# Combine the boxplot and density plot
pdf(file = "PageRank_GSE87038_v2.pdf", height = 10)
print(grid.arrange(pr, density_median_page_plot + coord_flip(), ncol = 2, widths = c(3, 1)))
print(grid.arrange(violin_median_page_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_t, pr, nrow = 2, heights = c(3, 3)))
dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


########### plot ANND (NOT USED   ) ############
{
    df <- readRDS(file = "IID_df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!

    df$label <- df$gene

    subset(df, signature == "HiGCTS_8")

    df <- rbind(
        subset(df, PPI_cat == "CTS"),
        subset(df, PPI_cat == "HiGCTS"),
        subset(df, PPI_cat == "HiG")
    )
    df$signature <- factor(df$signature, levels = unique(df$signature))

    # Step 1: Filter top 5 genes by annd within each signature
    top_genes <- df %>%
        group_by(signature) %>%
        arrange(desc(annd)) %>% # Sort in descending order of annd
        slice_head(n = 5) # Take the top 5 rows for each signature
    # subset the CHD genes within top 5
    CHD <- readRDS(file.path(inputdir, "CHD_Cilia_Genelist.rds"))

    top_genes_CHD <- subset(top_genes, toupper(gene) %in% toupper(unlist(CHD[c("Griffin2023_PCGC_AllCurated")])))
    (dim(top_genes)) # [1] 130  17
    (dim(top_genes_CHD)) # [1]  2 17
    print(top_genes_CHD, n = Inf, width = Inf)
#   signature gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <fct>     <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 CTS_8     NR2F2 0.0266   CTS          0.0333          0.994           25
# 2 HiG_8     ISL1  0.000552 HiG          0.00000604      0.935          169
#   rank_by_PR  annd p.annd  strength rank_by_strength normalized.strength
#        <dbl> <dbl>  <dbl>     <dbl>            <dbl>               <dbl>
# 1         16  7.53      0 0.00473                 18         0.0000947
# 2        276 74         0 0.0000599              276         0.000000206
#   rank_by_normalized.strength rank_by_ANND rank_by_p.ANND label
#                         <int>        <dbl>          <dbl> <chr>
# 1                          18            5           19.5 NR2F2
# 2                         276            2          140.  ISL1

    # Step 2: Create ggplot with boxplot and labels for top 5 genes
    pr <- ggplot(df[!is.na(df$annd), ], aes(x = signature, y = annd, colour = PPI_cat)) +
        geom_boxplot(show.legend = TRUE) + # Enable legend for the boxplot
        scale_color_manual(values = PPI_color_palette) +
        # Use ggrepel to avoid overlap and label top 5 genes based on annd
        geom_text_repel(
            data = top_genes_CHD, # Label only the top 5 genes
            aes(label = gene),
            size = 2, # Adjust the size of the text labels
            box.padding = 0.5, # Add space between the text and the data points
            point.padding = 0.5, # Add space between the text and the points
            segment.color = "grey50", # Color for the line connecting the text to the points
            max.overlaps = 20, # Max number of overlaps before labels stop being placed
            show.legend = FALSE # Do not show text labels in the legend
        ) +
        theme(
            legend.position = "none", # c(1, 1),
            legend.justification = c(0, 1), # Place legend at top-right corner
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        ) +
        # scale_x_discrete(limits = unique(df$signature)) + # Ensure x-axis respects the order of 'signature'
        labs(color = "PPI cat") # Optional: label for the color legend

    ## it is more manke sense to access each signature by its median page rank !!!!
    annd.median <- lapply(annd_observed, median, na.rm = TRUE) %>% unlist()

    A <- which(grepl("^CTS_", names(annd.median)))
    B <- which(grepl("^HiGCTS_", names(annd.median)))
    C <- which(grepl("^HiG_", names(annd.median)))
    ks.test(annd.median[C], annd.median[B]) # p-value = 0.1  HiG vs HiGCTS
    ks.test(annd.median[C], annd.median[A]) # p-value = 1.129e-05    HiG vs CTS
    ks.test(annd.median[B], annd.median[A]) # p-value = 0.2857	HiGCTS vs CTS
    wilcox.test(annd.median[C], annd.median[B]) # p-value = 0.1  HiG vs HiGCTS
    wilcox.test(annd.median[C], annd.median[A]) # p-value = 1.129e-05	HiG vs CTS
    wilcox.test(annd.median[B], annd.median[A]) # p-value = 0.2857	HiGCTS vs CTS


    df_median <- data.frame(
        annd.median = annd.median,
        PPI_cat = lapply(names(annd.median), function(x) unlist(strsplit(x, split = "_"))[1]) %>% unlist()
    )
    df_median$PPI_cat <- factor(df_median$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))
    density_median_annd_wilcox <- ggplot(df_median, aes(x = annd.median, color = PPI_cat, fill = PPI_cat)) +
        geom_density(alpha = 0.3) + # Density lines with transparency
        scale_color_manual(values = PPI_color_palette) +
        scale_fill_manual(values = PPI_color_palette) +
        theme_minimal() +
        theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
        labs(x = "Density of the median of ANND per PPI", y = "")
    # Add statistical comparisons using stat_compare_means manyally
    violin_median_annd_plot <- ggplot(df_median, aes(x = PPI_cat, y = annd.median, color = PPI_cat, fill = PPI_cat)) +
        geom_violin(alpha = 0.3) + # Violin plot with transparency
        scale_color_manual(values = PPI_color_palette) +
        scale_fill_manual(values = PPI_color_palette) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(x = "PPI category", y = "median of ANND per PPI") + # Label the axes
        # Add statistical comparisons using stat_compare_means
        stat_compare_means(
            aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
            comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
            method = "wilcox.test", # Non-parametric test (Wilcoxon)
            label = "p.signif", # Show significance labels (e.g., **, *, ns)
            label.x = 1.5, # Adjust x-position of the p-value text
            size = 4 # Adjust size of the p-value text
            , tip.length = 0
        ) +
        ggtitle("wilcox test, median ANND")

    # Combine the boxplot and density plot
    pdf(file = "annd_GSE87038_v2.pdf", height = 10)
    print(grid.arrange(pr, density_median_annd_wilcox + coord_flip(), ncol = 2, widths = c(3, 1)))
    print(grid.arrange(violin_median_annd_plot, pr, nrow = 2, heights = c(3, 3)))

    dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

}

############# plot normalized.strength (NO between-category difference! ) ###########
{
    df <- readRDS(file = "IID_df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
    df <- rbind(
        subset(df, PPI_cat == "CTS"),
        subset(df, PPI_cat == "HiGCTS"),
        subset(df, PPI_cat == "HiG")
    )
    df$signature <- factor(df$signature, levels = unique(df$signature))

    df$label <- df$gene
    subset(df, signature == "HiGCTS_8")

    CHD <- readRDS(file.path(inputdir, "CHD_Cilia_Genelist.rds"))
    df$PCGC_AllCurated <- toupper(df$gene) %in% toupper(unlist(CHD["Griffin2023_PCGC_AllCurated"]))

    # Step 1: Filter top 5 genes by normalized.strength within each signature
    top_genes <- df %>%
        group_by(signature) %>%
        arrange(desc(normalized.strength)) %>% # Sort in descending order of normalized.strength
        slice_head(n = 5) # Take the top 5 rows for each signature

    write.table(top_genes[, c(
        "signature", "gene", "normalized.strength", "PPI_cat", "rank_by_normalized.strength",
        "PCGC_AllCurated"
    )], file = "table_top5_strength_perPPI.tsv", sep = "\t", row.names = FALSE, quote = FALSE) # !!!!!!!!!!!!!!


    # subset the CHD genes within top 5
    top_genes_CHD <- subset(top_genes, PCGC_AllCurated == TRUE)
    (dim(top_genes)) # [1] 130   18
    (dim(top_genes_CHD)) # [1] 8  18

    # Violin plot with statistical comparisons
    violin_wilcox <- ggplot(df, aes(x = PPI_cat, y = log10(normalized.strength), color = PPI_cat, fill = PPI_cat)) +
        geom_violin(alpha = 0.3) + # Violin plot with transparency
        scale_color_manual(values = PPI_color_palette) +
        scale_fill_manual(values = PPI_color_palette) +
        theme_minimal() +
        theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
        labs(x = "PPI category", y = "log10 normalized.strength") + # Label the axes
        # Add statistical comparisons using stat_compare_means
        stat_compare_means(
            aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
            comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
            method = "wilcox.test", # Non-parametric test (Wilcoxon)
            label = "p.signif", # Show significance labels (e.g., **, *, ns)
            label.x = 1.5, # Adjust x-position of the p-value text
            size = 4 # Adjust size of the p-value text
            , tip.length = 0
        ) +
        ggtitle("wilcox test, all PPINs ")

    violin_t <- ggplot(df, aes(x = PPI_cat, y = log10(normalized.strength), color = PPI_cat, fill = PPI_cat)) +
        geom_violin(alpha = 0.3) + # Violin plot with transparency
        scale_color_manual(values = PPI_color_palette) +
        scale_fill_manual(values = PPI_color_palette) +
        theme_minimal() +
        theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
        labs(x = "PPI category", y = "normalized.strength") + # Label the axes
        # Add statistical comparisons using stat_compare_means
        stat_compare_means(
            aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
            comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
            method = "t.test", # Non-parametric test (Wilcoxon)
            label = "p.signif", # Show significance labels (e.g., **, *, ns)
            label.x = 1.5, # Adjust x-position of the p-value text
            size = 4 # Adjust size of the p-value text
            , tip.length = 0
        ) +
        ggtitle("t test, all PPINs ")

    # Step 2: Create ggplot with boxplot and labels for top 5 genes
    pr <- ggplot(df, aes(x = signature, y = log10(normalized.strength), colour = PPI_cat)) +
        geom_boxplot(show.legend = TRUE) + # Enable legend for the boxplot
        scale_color_manual(values = PPI_color_palette) +
        # Use ggrepel to avoid overlap and label top 5 genes based on normalized.strength
        geom_text_repel(
            data = top_genes_CHD, # Label only the top 5 genes
            aes(label = gene),
            size = 2, # Adjust the size of the text labels
            box.padding = 0.5, # Add space between the text and the data points
            point.padding = 0.5, # Add space between the text and the points
            segment.color = "grey50", # Color for the line connecting the text to the points
            max.overlaps = 20, # Max number of overlaps before labels stop being placed
            show.legend = FALSE # Do not show text labels in the legend
        ) +
        theme(
            legend.position = "none", # c(1, 1),
            legend.justification = c(0, 1), # Place legend at top-right corner
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        ) +
        scale_x_discrete(limits = unique(df$signature)) + # Ensure x-axis respects the order of 'signature'
        labs(color = "PPI cat") # Optional: label for the color legend

    ## it is more manke sense to access each signature by its median page rank !!!!
    df_median <- df %>%
        group_by(signature) %>%
        summarise(median_normalized_strength = median(normalized.strength, na.rm = TRUE))

    A <- which(grepl("^CTS_", df_median$signature))
    B <- which(grepl("^HiGCTS_", df_median$signature))
    C <- which(grepl("^HiG_", df_median$signature))
    ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value = 0.1  HiG vs HiGCTS
    ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value = 7.905e-05	HiG vs CTS
    ks.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value =  0.2857	HiGCTS vs CTS
    wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value = 0.1  HiG vs HiGCTS
    wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value = 2.259e-05	HiG vs CTS
    wilcox.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value =  0.2857	HiGCTS vs CTS


    df_median$PPI_cat <- lapply(df_median$signature %>% as.vector(), function(x) unlist(strsplit(x, split = "_"))[1]) %>% unlist()
    df_median$PPI_cat <- factor(df_median$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

    density_median_normalized.strength_plot <- ggplot(df_median, aes(
        x = log10(median_normalized_strength),
        color = PPI_cat, fill = PPI_cat
    )) +
        geom_density(alpha = 0.3) + # Density lines with transparency
        scale_color_manual(values = PPI_color_palette) +
        scale_fill_manual(values = PPI_color_palette) +
        theme_minimal() +
        theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
        labs(x = "Density of the median of normalzied node strength per PPI", y = "")
    # Add statistical comparisons using stat_compare_means manyally

    violin_median_normalized.strength_wilcox <- ggplot(df_median, aes(
        x = PPI_cat,
        y = log10(median_normalized_strength), color = PPI_cat, fill = PPI_cat
    )) +
        geom_violin(alpha = 0.3) + # Violin plot with transparency
        scale_color_manual(values = PPI_color_palette) +
        scale_fill_manual(values = PPI_color_palette) +
        theme_minimal() +
        theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
        labs(x = "PPI category", y = "log10. median of normalized node strength per PPI") + # Label the axes
        # Add statistical comparisons using stat_compare_means
        stat_compare_means(
            aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
            comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
            method = "wilcox.test", # Non-parametric test (Wilcoxon)
            label = "p.signif", # Show significance labels (e.g., **, *, ns)
            label.x = 1.5, # Adjust x-position of the p-value text
            size = 4 # Adjust size of the p-value text
            , tip.length = 0
        ) +
        ggtitle("wilcox, median nr_strength")

    # Combine the boxplot and density plot
    pdf(file = "normalized.node.strength_GSE87038_v2.pdf", height = 10)
    print(grid.arrange(pr, density_median_normalized.strength_plot + coord_flip(), ncol = 2, widths = c(3, 1)))
    print(grid.arrange(violin_median_normalized.strength_wilcox, pr, nrow = 2, heights = c(3, 3)))
    print(grid.arrange(violin_wilcox, pr, nrow = 2, heights = c(3, 3)))
    print(grid.arrange(violin_t, pr, nrow = 2, heights = c(3, 3)))

    dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}




# number of significantly high annd per PPI_cats
n.annd.high <- lapply(names(graph_list), function(x) nrow(subset(df, signature == x & round(p.annd, 2) <= 0.05))) %>% unlist()
names(n.annd.high) <- names(graph_list)
(n.annd.high)
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_9    HiG_10
#       261       373       357       253       264       449       303       363
#    HiG_12    HiG_14    HiG_17    HiG_18    HiG_19     HiG_7    HiG_11    HiG_15
#       292       315       453       399       446       278       379       342
#    HiG_16    HiG_13     HiG_8 HiGCTS_15     CTS_7    CTS_11    CTS_15    CTS_16
#       458       348       278        11         9        18        37        19
#  CTS_16.1     CTS_8
#        43        38

df_compare <- data.frame(
    signature = names(graph_list),
    n_sig.pagerank = n.pr.high,
    n_sig.annd = n.annd.high
)
df_compare$PPI_cat <- lapply(df_compare$signature, function(x) unlist(strsplit(x, split = "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))
ggplot(df_compare, aes(x = n_sig.pagerank, y = n_sig.annd)) +
    geom_point(aes(shape = PPI_cat, colour = PPI_cat), show.legend = FALSE) +
    scale_color_manual(values = PPI_color_palette) +
    geom_text_repel(aes(label = signature), hjust = -0.1, vjust = 0) +
    theme(legend.position = c(0, 0), legend.justification = c(0, 0))

#####################


##########################

###########################################################################################################################################
## Given a transitional state, CTS&HiG genes exhibit higher betweenness centrality in the CTS-derived network and the HiG-derived network
###########################################################################################################################################
bc <- read.table(file = "df_betweeness.tsv", header = TRUE)
(dim(bc)) # [1] 7265    5
(colnames(bc))
# [1] "signature"             "BetweennessCentrality" "gene"                  "rank_by_BC"            "PPI_cat"
## find out the cluster with CTSHiG
x <- grep("HiGCTS_", bc$signature, value = TRUE) %>% unique()
CTS <- lapply(x, function(x) unlist(strsplit(x, split = "_"))[2]) %>% unlist()
y <- grepl(".", CTS, fixed = TRUE)
if (any(y)) CT <- CTS[!y] else CT <- CTS

(CTS)
# [1] "15"
(CT)
# [1] "15"

plot_CTS_bc <- plot_HiG_bc <- list()
## compare the CTS&HiG genes to non-hiG genes per CTS-derived PPI
for (id in CTS) {
    if (grepl(".", id, fixed = TRUE)) id2 <- unlist(strsplit(id, split = ".", fixed = TRUE))[1] else id2 <- id
    CTS_PPI <- subset(bc, signature == paste("CTS", id, sep = "_"))
    HiG_PPI <- subset(bc, signature == paste("HiG", id2, sep = "_"))
    CTS_PPI$isHiG <- factor(CTS_PPI$gene %in% HiG_PPI$gene, levels = c("FALSE", "TRUE"))
    HiG_PPI$isCTS <- factor(HiG_PPI$gene %in% CTS_PPI$gene, levels = c("FALSE", "TRUE"))

    ############## ranked by betweenness
    # Compute exact p-value
    pval <- wilcox.test(BetweennessCentrality ~ isHiG, data = CTS_PPI)$p.value

    # Reorder gene factor levels by BetweennessCentrality (high to low)
    CTS_PPI <- CTS_PPI %>%
        mutate(gene = factor(gene, levels = gene[order(-BetweennessCentrality)]))

    plot_CTS_bc[[id]] <- ggplot(CTS_PPI, aes(x = gene, y = BetweennessCentrality, fill = isHiG)) +
        geom_boxplot(aes(group = isHiG), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) + # Add boxplot first
        geom_point(aes(color = isHiG), size = 3) +
        scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E7298A")) +
        geom_text_repel(aes(label = gene, color = isHiG), hjust = -0.1, vjust = 0) +
        theme(
            legend.position = c(0, 0),
            legend.justification = c(1, 1),
            axis.text.x = element_blank(), # Remove x-axis labels
            axis.ticks.x = element_blank()
        ) +
        annotate("text",
            x = (length(unique(CTS_PPI$gene)) + 1) / 2,
            y = max(CTS_PPI$BetweennessCentrality) * 0.8,
            label = paste0("wilcox p = ", signif(pval, 2), " F,T: ", table(CTS_PPI$isHiG) %>% toString()), size = 4
        )


    # Compute exact p-value
    pval <- wilcox.test(BetweennessCentrality ~ isCTS, data = HiG_PPI)$p.value

    # Reorder gene factor levels by BetweennessCentrality (high to low)
    HiG_PPI <- HiG_PPI %>%
        mutate(gene = factor(gene, levels = gene[order(-BetweennessCentrality)]))

    plot_HiG_bc[[id]] <- ggplot(HiG_PPI, aes(x = gene, y = BetweennessCentrality, fill = isCTS)) +
        geom_boxplot(aes(group = isCTS), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) + # Add boxplot first
        geom_point(aes(color = isCTS), size = 3) +
        scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E6AB02")) +
        geom_text_repel(aes(label = gene, color = isCTS), hjust = -0.1, vjust = 0) +
        theme(
            legend.position = c(0, 0),
            legend.justification = c(1, 1),
            axis.text.x = element_blank(), # Remove x-axis labels
            axis.ticks.x = element_blank()
        ) +
        annotate("text",
            x = (length(unique(HiG_PPI$gene)) + 1) / 2,
            y = max(HiG_PPI$BetweennessCentrality) * 0.8,
            label = paste0("wilcox p = ", signif(pval, 2), " F,T: ", table(HiG_PPI$isCTS) %>% toString()), size = 4
        )
}
#####

df <- readRDS(file = "IID_df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
(colnames(df))
#  [1] "signature"                   "gene"                        "PageRank"                    "PPI_cat"                     "EigenCentrality"             "p.PageRank"
#  [7] "rank_by_p.PR"                "rank_by_PR"                  "annd"                        "p.annd"                      "strength"                    "rank_by_strength"
# [13] "normalized.strength"         "rank_by_normalized.strength" "rank_by_ANND"                "rank_by_p.ANND"

plot_CTS_pr <- plot_HiG_pr <- list()

## compare the CTS&HiG genes to non-HiG genes per CTS-derived PPI
for (id in CTS) {
    if (grepl(".", id, fixed = TRUE)) id2 <- unlist(strsplit(id, split = ".", fixed = TRUE))[1] else id2 <- id
    CTS_PPI <- subset(df, signature == paste("CTS", id, sep = "_"))
    HiG_PPI <- subset(df, signature == paste("HiG", id2, sep = "_"))
    CTS_PPI$isHiG <- factor(CTS_PPI$gene %in% HiG_PPI$gene, levels = c("FALSE", "TRUE"))
    HiG_PPI$isCTS <- factor(HiG_PPI$gene %in% CTS_PPI$gene, levels = c("FALSE", "TRUE"))

    ############## ranked by pageRank
    # Compute exact p-value
    pval <- wilcox.test(PageRank ~ isHiG, data = CTS_PPI)$p.value

    # Reorder gene factor levels by PageRank (high to low)
    CTS_PPI <- CTS_PPI %>%
        mutate(gene = factor(gene, levels = gene[order(-PageRank)]))

    plot_CTS_pr[[id]] <- ggplot(CTS_PPI, aes(x = gene, y = PageRank, fill = isHiG)) +
        geom_boxplot(aes(group = isHiG), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) + # Add boxplot first
        geom_point(aes(color = isHiG), size = 3) +
        scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E7298A")) +
        geom_text_repel(aes(label = gene, color = isHiG), hjust = -0.1, vjust = 0) +
        theme(
            legend.position = c(0, 0),
            legend.justification = c(1, 1),
            axis.text.x = element_blank(), # Remove x-axis labels
            axis.ticks.x = element_blank()
        ) +
        annotate("text",
            x = (length(unique(CTS_PPI$gene)) + 1) / 2,
            y = max(CTS_PPI$PageRank) * 0.8,
            label = paste0("wilcox p = ", signif(pval, 2), " F,T: ", table(CTS_PPI$isHiG) %>% toString()), size = 4
        )


    # Compute exact p-value
    pval <- wilcox.test(PageRank ~ isCTS, data = HiG_PPI)$p.value

    # Reorder gene factor levels by PageRank (high to low)
    HiG_PPI <- HiG_PPI %>%
        mutate(gene = factor(gene, levels = gene[order(-PageRank)]))

    plot_HiG_pr[[id]] <- ggplot(HiG_PPI, aes(x = gene, y = PageRank, fill = isCTS)) +
        geom_boxplot(aes(group = isCTS), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) + # Add boxplot first
        geom_point(aes(color = isCTS), size = 3) +
        scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E6AB02")) +
        geom_text_repel(aes(label = gene, color = isCTS), hjust = -0.1, vjust = 0) +
        theme(
            legend.position = c(0, 0),
            legend.justification = c(1, 1),
            axis.text.x = element_blank(), # Remove x-axis labels
            axis.ticks.x = element_blank()
        ) +
        annotate("text",
            x = (length(unique(HiG_PPI$gene)) + 1) / 2,
            y = max(HiG_PPI$PageRank) * 0.8,
            label = paste0("wilcox p = ", signif(pval, 2), " F,T: ", table(HiG_PPI$isCTS) %>% toString()), size = 4
        )
}

plot_CTS_annd <- plot_HiG_annd <- list()

## compare the CTS&HiG genes to non-hiG genes per CTS-derived PPI
for (id in CTS) {
    if (grepl(".", id, fixed = TRUE)) id2 <- unlist(strsplit(id, split = ".", fixed = TRUE))[1] else id2 <- id
    CTS_PPI <- subset(df, signature == paste("CTS", id, sep = "_"))
    HiG_PPI <- subset(df, signature == paste("HiG", id2, sep = "_"))
    CTS_PPI$isHiG <- factor(CTS_PPI$gene %in% HiG_PPI$gene, levels = c("FALSE", "TRUE"))
    HiG_PPI$isCTS <- factor(HiG_PPI$gene %in% CTS_PPI$gene, levels = c("FALSE", "TRUE"))

    ############## ranked by annd
    # Compute exact p-value
    pval <- wilcox.test(annd ~ isHiG, data = CTS_PPI)$p.value

    # Reorder gene factor levels by annd (high to low)
    CTS_PPI <- CTS_PPI %>%
        mutate(gene = factor(gene, levels = gene[order(-annd)]))

    plot_CTS_annd[[id]] <- ggplot(CTS_PPI, aes(x = gene, y = annd, fill = isHiG)) +
        geom_boxplot(aes(group = isHiG), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) + # Add boxplot first
        geom_point(aes(color = isHiG), size = 3) +
        scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E7298A")) +
        geom_text_repel(aes(label = gene, color = isHiG), hjust = -0.1, vjust = 0) +
        theme(
            legend.position = c(0, 0),
            legend.justification = c(1, 1),
            axis.text.x = element_blank(), # Remove x-axis labels
            axis.ticks.x = element_blank()
        ) +
        annotate("text",
            x = (length(unique(CTS_PPI$gene)) + 1) / 2,
            y = max(CTS_PPI$annd, na.rm = T) * 0.8,
            label = paste0("wilcox p = ", signif(pval, 2), " F,T: ", table(CTS_PPI$isHiG) %>% toString()), size = 4
        )


    # Compute exact p-value
    pval <- wilcox.test(annd ~ isCTS, data = HiG_PPI)$p.value

    # Reorder gene factor levels by annd (high to low)
    HiG_PPI <- HiG_PPI %>%
        mutate(gene = factor(gene, levels = gene[order(-annd)]))

    plot_HiG_annd[[id]] <- ggplot(HiG_PPI, aes(x = gene, y = annd, fill = isCTS)) +
        geom_boxplot(aes(group = isCTS), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) + # Add boxplot first
        geom_point(aes(color = isCTS), size = 3) +
        scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E6AB02")) +
        geom_text_repel(aes(label = gene, color = isCTS), hjust = -0.1, vjust = 0) +
        theme(
            legend.position = c(0, 0),
            legend.justification = c(1, 1),
            axis.text.x = element_blank(), # Remove x-axis labels
            axis.ticks.x = element_blank()
        ) +
        annotate("text",
            x = (length(unique(HiG_PPI$gene)) + 1) / 2,
            y = max(HiG_PPI$annd, na.rm = T) * 0.8,
            label = paste0("wilcox p = ", signif(pval, 2), " F,T: ", table(HiG_PPI$isCTS) %>% toString()), size = 4
        )
}


pdf(file = "gene_ranked_by_importance_dotBoxPlot_6panel.pdf", height = 10.5)

for (id in CTS) {
  # handle possible suffixes in HiG names
  id2 <- if (grepl("\\.", id)) sub("\\..*", "", id) else id
  
  # CTS plots
  p_cts_annd <- plot_CTS_annd[[id]]
  p_cts_bc   <- plot_CTS_bc[[id]]
  p_cts_pr   <- plot_CTS_pr[[id]]
  
  # HiG plots: find matching names containing the CTS id
  hiG_matches <- grep(id2, names(plot_HiG_annd), value = TRUE)
  p_hig_annd <- if(length(hiG_matches) > 0) plot_HiG_annd[[hiG_matches[1]]] else NULL
  
  hiG_matches_bc <- grep(id2, names(plot_HiG_bc), value = TRUE)
  p_hig_bc   <- if(length(hiG_matches_bc) > 0) plot_HiG_bc[[hiG_matches_bc[1]]] else NULL
  
  hiG_matches_pr <- grep(id2, names(plot_HiG_pr), value = TRUE)
  p_hig_pr   <- if(length(hiG_matches_pr) > 0) plot_HiG_pr[[hiG_matches_pr[1]]] else NULL
  
  # Combine all 6 plots
  plot_list <- list(p_cts_annd, p_hig_annd,
                    p_cts_bc, p_hig_bc,
                    p_cts_pr, p_hig_pr)
  
  # Remove any NULLs in case some HiG plots don't exist
  plot_list <- Filter(Negate(is.null), plot_list)
  
  # Labels
  labels <- c(paste0("CTS_", id, " network"),
              paste0("HiG_", id2, " network"),
              paste0("CTS_", id, " network"),
              paste0("HiG_", id2, " network"),
              paste0("CTS_", id, " network"),
              paste0("HiG_", id2, " network"))
  labels <- labels[seq_along(plot_list)]
  
  # Arrange in 2 columns, 3 rows (adjust automatically if some plots are missing)
  x <- ggarrange(plotlist = plot_list,
                 labels = labels,
                 ncol = 2,
                 nrow = ceiling(length(plot_list)/2))
  
  print(x)
}

dev.off()