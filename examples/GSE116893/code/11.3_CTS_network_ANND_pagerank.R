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

########## BEGINNING OF USER INPUT ##########

wd = "C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- paste0(wd, "data/")
shared_path <- paste0(wd, "../Shared_Data/")

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

db <- "GSE116893"

s <- "combined" # specificity method

########## END OF USER INPUT ##########

file <- paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"
#  [7] "HiG_7"     "HiG_8"     "HiG_10"    "HiG_11"    "HiG_12"    "HiG_13"
# [13] "HiG_14"    "HiG_16"    "HiG_17"    "HiG_15"    "HiG_9"     "HiGCTS_15"
# [19] "HiGCTS_16" "CTS_15"    "CTS_16"    "CTS_9"
edge_counts <- sapply(graph_list, ecount)
(edge_counts)
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_7     HiG_8
#     13101      4272      3520      4372      4623      9872      4363      3946
#    HiG_10    HiG_11    HiG_12    HiG_13    HiG_14    HiG_16    HiG_17    HiG_15
#      5037     10588     10459     16636     24134     10239     37844     23162
#     HiG_9 HiGCTS_15 HiGCTS_16    CTS_15    CTS_16     CTS_9
#      2407       381       413       931       493       454
(sapply(graph_list, vcount))
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_7     HiG_8
#       372       354       284       173       353       579       335       311
#    HiG_10    HiG_11    HiG_12    HiG_13    HiG_14    HiG_16    HiG_17    HiG_15
#       278       520       548       378      1044       555      1385       728
#     HiG_9 HiGCTS_15 HiGCTS_16    CTS_15    CTS_16     CTS_9
#       221        50       124        75       138       110
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
(dim(df)[1]) # 8915

ic <- lapply(graph_list, function(x) eigen_centrality(x, weights = E(x)$weight)$vector)
IC <- lapply(ic, function(x) data.frame(EigenCentrality = x, gene = names(x)) %>% arrange(desc(EigenCentrality))) %>%
    rbindlist(., idcol = names(.))
colnames(IC)[1] <- "signature"
(dim(IC)) # [1] 8915    3
df <- merge(df, IC, by = c("signature", "gene"))
(dim(df)) # [1] 8915    5



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
    #      PageRank     gene
    #         <num>   <char>
    # 1: 0.01578567 Hsp90aa1
    # 2: 0.01442298    H3f3b
    # 3: 0.01400916     Npm1
    # 4: 0.01393744    Hspa8
    # 5: 0.01310504      Ncl
    # 6: 0.01188350    Hspd1
    for (j in V(graph_list[[i]])$name) {
        pr_P[[i]][j] <- length(which(subset(tmp, gene == j)$PageRank >= page[[i]][j])) / N
    }
}
saveRDS(pr_P, file = paste0(db, "_PageRank_Pvalue_by_rewiring.rds"))


pr_P <- readRDS(file = paste0(db, "_PageRank_Pvalue_by_rewiring.rds"))
tmp <- lapply(pr_P, function(x) data.frame(p.PageRank = x, gene = names(x))) %>%
    rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"

df <- merge(df, tmp, by = c("signature", "gene")) %>%
    group_by(signature) %>%
    mutate(rank_by_p.PR = rank(p.PageRank)) %>%
    mutate(rank_by_PR = rank(-PageRank)) %>%
    ungroup()

(head(df))
# # A tibble: 6 × 8
#   signature gene     PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>     <chr>       <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 CTS_15    ADD3      0.00339 CTS            1.19e- 2      0.598         44.5
# 2 CTS_15    ANP32B    0.00461 CTS            4.58e- 3      0.369         23
# 3 CTS_15    ANP32E    0.00960 CTS            1.40e- 1      0.322         18
# 4 CTS_15    APOLD1    0.00215 CTS            3.94e-19      0.635         60.5
# 5 CTS_15    ARHGAP35  0.00253 CTS            3.42e- 3      0.598         44.5
# 6 CTS_15    ASPM      0.0355  CTS            9.48e- 1      0.615         49.5
# # ℹ 1 more variable: rank_by_PR <dbl>
(dim(df)) # [1] 8915    8

(subset(df, tolower(gene) == "isl1"))
# # A tibble: 1 × 8
#   signature gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>     <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 HiG_12    ISL1   0.00166 HiG             0.00229      0.275          318
# # ℹ 1 more variable: rank_by_PR <dbl>

# number of significantly high pagerank per PPI_cats, too much control !
n.pr.high <- lapply(names(graph_list), function(x) nrow(subset(df, signature == x & p.PageRank < 0.05))) %>% unlist()
names(n.pr.high) <- names(graph_list)
(n.pr.high)
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_7     HiG_8
#        77         0         0        55         0         0         0         0
#    HiG_10    HiG_11    HiG_12    HiG_13    HiG_14    HiG_16    HiG_17    HiG_15
#         2         0         0       118         0         3         0         0
#     HiG_9 HiGCTS_15 HiGCTS_16    CTS_15    CTS_16     CTS_9
#         0         0         0         1         0         0

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
    V_Isolated <- which(degree(G) == 0)
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
    #        annd   gene
    #       <num> <char>
    # 1: 85.12626  Fkbp4
    # 2: 83.20335 Dnaja1
    # 3: 80.58264 Fam60a
    for (j in V(graph_list[[i]])$name) {
        annd_P[[i]][j] <- length(which(subset(tmp, gene == j)$knn >= annd_observed[[i]][j])) / N
    }
}
saveRDS(annd_P, file = paste0(db, "_annd_Pvalue_by_rewiring.rds"))

annd_P <- readRDS(file = paste0(db, "_annd_Pvalue_by_rewiring.rds"))

(unique(df$PPI_cat)) # CTS    HiGCTS HiG

df <- rbind(
    subset(df, PPI_cat == "CTS"),
    subset(df, PPI_cat == "HiGCTS"),
    subset(df, PPI_cat == "HiG")
)
(unique(df$PPI_cat)) # CTS    HiGCTS HiG

tmp <- lapply(annd_observed, function(x) data.frame(annd = x, gene = names(x)) %>% arrange(desc(annd))) %>%
    rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"
df <- merge(df, tmp, by = c("signature", "gene"))
(dim(df)) # [1] 8878    9

annd_P[["CTS_8"]]
# NULL  (CTS_8 does not exist in GSE116893 graph_list)
tmp <- lapply(annd_P, function(x) data.frame(p.annd = x, gene = names(x))) %>%
    rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"
(dim(tmp)) # [1] 8915    3

df <- merge(df, tmp, by = c("signature", "gene"))
df[which(is.na(df$knn)), "p.annd"] <- NA ## due to nrow(df) 4878 > nrow(tmp)
(dim(df)) # [1] 8878   10

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
(dim(V_strength)) # 8915    6

## add the V_strength & V_strength_norm infor
df <- merge(df, V_strength, by = c("signature", "gene"))
(dim(df)) # 8878   14
(head(df, 3))
#   signature   gene    PageRank PPI_cat EigenCentrality p.PageRank
# 1    CTS_15   ADD3 0.003393514     CTS     0.011929784      0.598
# 2    CTS_15 ANP32B 0.004606449     CTS     0.004580529      0.369
# 3    CTS_15 ANP32E 0.009597013     CTS     0.139993484      0.322
#   rank_by_p.PR rank_by_PR     annd p.annd   strength rank_by_strength
# 1         44.5         61 38.95112      0 0.1665578               57
# 2         23.0         54 19.16907      0 0.1851355               56
# 3         18.0         41 37.96229      0 1.0066784               39
#   normalized.strength rank_by_normalized.strength
# 1         0.002250781                          57
# 2         0.002501831                          56
# 3         0.013603762                          39

(table(df$normalized.strength >= df$annd))
# FALSE
#  8878

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
#  [1] "signature"                   "gene"
#  [3] "PageRank"                    "PPI_cat"
#  [5] "EigenCentrality"             "p.PageRank"
#  [7] "rank_by_p.PR"                "rank_by_PR"
#  [9] "annd"                        "p.annd"
# [11] "strength"                    "rank_by_strength"
# [13] "normalized.strength"         "rank_by_normalized.strength"
# [15] "rank_by_ANND"                "rank_by_p.ANND"
saveRDS(df, file = "df_PAGERANK_strength_ANND.rewiring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
write.table(df, file = "df_PAGERANK_strength_ANND.rewiring.P.tsv", sep = "\t", row.names = F) # !!!!!!!!

##########################
## add the column of betweenness centrality
##########################
df <- readRDS(file = "df_PAGERANK_strength_ANND.rewiring.P.rds")

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

(dim(df_BC)) # [1] 8915    5

write.table(df_BC, file = "df_betweeness.tsv", sep = "\t", row.names = F) # !!!!!!!!


########### betweenness centrality ############
df_BC <- read.table(file = "df_betweeness.tsv", sep = "\t", header = T)
df_BC$PPI_cat <- factor(df_BC$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

CHD <- readRDS(file = paste0(shared_path, "CHD_Cilia_Genelist.rds"))
df_BC$PCGC_AllCurated <- toupper(df_BC$gene) %in% toupper(unlist(CHD["Griffin2023_PCGC_AllCurated"]))

# Calculate top 5 significant genes within each box
df5 <- df_BC %>%
    filter(rank_by_BC <= 5 & BetweennessCentrality > 0) %>%
    ungroup()

write.table(df5[, c(
    "signature", "gene", "BetweennessCentrality", "PPI_cat", "rank_by_BC",
    "PCGC_AllCurated"
)], file = "table_top5_Betweenness_perPPI.tsv", sep = "\t", row.names = FALSE, quote = FALSE) # !!!!!!!!!!!!!!

(subset(df5, PPI_cat == "HiGCTS"))
#    signature BetweennessCentrality    gene rank_by_BC PPI_cat PCGC_AllCurated
# 86 HiGCTS_15                   122    PFN1          5  HiGCTS           FALSE
# 87 HiGCTS_15                   147   KIF11          3  HiGCTS           FALSE
# 88 HiGCTS_15                   146    ASPM          4  HiGCTS           FALSE
# 89 HiGCTS_15                   193  DIAPH3          2  HiGCTS           FALSE
# 90 HiGCTS_15                   260   TOP2A          1  HiGCTS           FALSE
# 91 HiGCTS_16                  3461  SCARB1          1  HiGCTS           FALSE
# 92 HiGCTS_16                  1130  ATP1B3          3  HiGCTS           FALSE
# 93 HiGCTS_16                   916   ITPR1          5  HiGCTS           FALSE
# 94 HiGCTS_16                   930 CYP17A1          4  HiGCTS           FALSE
# 95 HiGCTS_16                  1186    ERN1          2  HiGCTS           FALSE

(df5_CHD <- subset(df5, PCGC_AllCurated == TRUE))

(df5_CHD)
#   signature BetweennessCentrality    gene rank_by_BC PPI_cat PCGC_AllCurated
# 6     HiG_2                  4750 CACNA1C          5     HiG            TRUE
# 11    HiG_3                  3412    MYCN          4     HiG            TRUE
# 27    HiG_6                 17932 CACNA1C          1     HiG            TRUE
# 28    HiG_6                 12684    MYCN          2     HiG            TRUE
# 31    HiG_7                  7349    MYCN          1     HiG            TRUE
# 36    HiG_8                  4973 CACNA1C          5     HiG            TRUE
# 61   HiG_14                 51307  SEMA3A          3     HiG            TRUE
# 72   HiG_17                 55448     NF1          1     HiG            TRUE



df_BC <- rbind(
    subset(df_BC, PPI_cat == "CTS"),
    subset(df_BC, PPI_cat == "HiGCTS"),
    subset(df_BC, PPI_cat == "HiG")
)
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

# Caused by error in `kruskal.test.default()`:
# ! 'x' and 'g' must have the same length


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
ks.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.1754   HiG vs HiGCTS
ks.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 0.0614   HiG vs CTS
ks.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = 0.6      HiGCTS vs CTS
wilcox.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.1106   HiG vs HiGCTS
wilcox.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 0.08049  HiG vs CTS
wilcox.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = 0.5536   HiGCTS vs CTS
t.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.003025  HiG vs HiGCTS
t.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 0.008534  HiG vs CTS
t.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = 0.4034    HiGCTS vs CTS


density_median_bc_plot <- ggplot(df_median, aes(x = log10(bc.median), color = PPI_cat, fill = PPI_cat)) +
    geom_density(alpha = 0.3) + # Density lines with transparency
    scale_color_manual(values = PPI_color_palette) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "Density of the median of BetweennessCentralitys per PPI", y = "")
# Add statistical comparisons using stat_compare_means manyally

x <- which(df_median$bc.median == 0) # 5

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
pdf(file = paste0("BetweennessCentrality_", db, "_v2.pdf"), height = 10)
print(grid.arrange(pr_repel, density_median_bc_plot + coord_flip(), ncol = 2, widths = c(3, 1)))
print(grid.arrange(violin_median_bc_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_t, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_median_bc_wilcox, violin_median_bc_wilcox_ln, nrow = 2))
dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


########### plot PageRank ############
df <- readRDS(file = "df_PAGERANK_strength_ANND.rewiring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
(dim(df)) # [1] 8878   16

df <- rbind(
    subset(df, PPI_cat == "CTS"),
    subset(df, PPI_cat == "HiGCTS"),
    subset(df, PPI_cat == "HiG")
)

CHD <- readRDS(file = paste0(shared_path, "CHD_Cilia_Genelist.rds"))
df$PCGC_AllCurated <- toupper(df$gene) %in% toupper(unlist(CHD["Griffin2023_PCGC_AllCurated"]))

# Calculate top 5 significant genes within each box
df5 <- df %>%
    filter(rank_by_PR <= 5) %>%
    ungroup()
(subset(df5, PPI_cat == "HiGCTS"))
# # A tibble: 10 × 17
#    signature gene    PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#    <chr>     <chr>      <dbl> <fct>             <dbl>      <dbl>        <dbl>
#  1 HiGCTS_15 ASPM      0.0528 HiGCTS            0.922      0.659         37
#  2 HiGCTS_15 KIF11     0.0573 HiGCTS            1          0.62          33.5
#  3 HiGCTS_15 KIF15     0.0446 HiGCTS            0.820      0.665         39
#  4 HiGCTS_15 TOP2A     0.0579 HiGCTS            0.919      0.354         15
#  5 HiGCTS_15 TPX2      0.0536 HiGCTS            0.948      0.589         30
#  6 HiGCTS_16 CYP11B1   0.0285 HiGCTS            0.879      0.708        115
#  7 HiGCTS_16 CYP17A1   0.0328 HiGCTS            0.972      0.713        117
#  8 HiGCTS_16 FDX1      0.0310 HiGCTS            0.990      0.691        110
#  9 HiGCTS_16 SCARB1    0.0493 HiGCTS            0.754      0.62          98
# 10 HiGCTS_16 STAR      0.0325 HiGCTS            1          0.701        112.
# # ℹ 10 more variables: rank_by_PR <dbl>, annd <dbl>, p.annd <dbl>,
# #   strength <dbl>, rank_by_strength <dbl>, normalized.strength <dbl>,
# #   rank_by_normalized.strength <int>, rank_by_ANND <dbl>,
# #   rank_by_p.ANND <dbl>, PCGC_AllCurated <lgl>

write.table(df5[, c(
    "signature", "gene", "PageRank", "PPI_cat", "rank_by_PR",
    "normalized.strength", "rank_by_normalized.strength", "PCGC_AllCurated"
)], file = "table_top5_PageRank_perPPI.tsv", sep = "\t", row.names = FALSE, quote = FALSE) # !!!!!!!!!!!!!!

df5_CHD <- subset(df5, PCGC_AllCurated == TRUE)
(dim(df5)) # [1] 110  17
(dim(df5_CHD)) # [1]  5 17
(df5_CHD %>% as.data.frame())
#   signature    gene    PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
# 1    HiG_15    PTEN 0.005393929     HiG      0.19125682      0.169         55.5
# 2    HiG_17     NF1 0.003123284     HiG      0.68691156      0.426        779.0
# 3     HiG_2 CACNA1C 0.009962531     HiG      0.03000444      0.605        323.0
# 4     HiG_6 CACNA1C 0.007879795     HiG      1.00000000      0.763        535.0
# 5     HiG_8 CACNA1C 0.011616202     HiG      0.82370629      0.682        264.5
#   rank_by_PR      annd p.annd strength rank_by_strength normalized.strength
# 1          4 100.52816      0 3.197282               11         0.004397912
# 2          5  94.39445      0 2.433973                8         0.001758651
# 3          4  44.70995      0 3.078285                3         0.008720355
# 4          2  65.03190      0 5.753122                1         0.009953498
# 5          5  48.49607      0 2.657005                4         0.008570985
#   rank_by_normalized.strength rank_by_ANND rank_by_p.ANND PCGC_AllCurated
# 1                          11          385          364.5            TRUE
# 2                           8          446          691.5            TRUE
# 3                           3           84          177.5            TRUE
# 4                           1          115          290.0            TRUE
# 5                           4           54          155.5            TRUE

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
ks.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value = 0.0117    HiG vs HiGCTS
ks.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value = 0.001754  HiG vs CTS
ks.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value = 0.9       HiGCTS vs CTS
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value = 0.0117    HiG vs HiGCTS
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value = 0.001754  HiG vs CTS
wilcox.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value = 0.8       HiGCTS vs CTS


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
pdf(file = paste0("PageRank_", db, "_v2.pdf"), height = 10)
print(grid.arrange(pr, density_median_page_plot + coord_flip(), ncol = 2, widths = c(3, 1)))
print(grid.arrange(violin_median_page_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_t, pr, nrow = 2, heights = c(3, 3)))
dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

########### plot ANND (NOT USED   ) ############
{
    df <- readRDS(file = "df_PAGERANK_strength_ANND.rewiring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!

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
    CHD <- readRDS(file = paste0(shared_path, "CHD_Cilia_Genelist.rds"))

    top_genes_CHD <- subset(top_genes, toupper(gene) %in% toupper(unlist(CHD[c("Griffin2023_PCGC_AllCurated")])))
    (dim(top_genes)) # [1] 110  17  (22 graphs × 5; not auto-printed inside {})
    (dim(top_genes_CHD)) # not captured in Rscript output (inside {} block)
    (top_genes_CHD)
    # not captured in Rscript output (inside {} block)

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
    ks.test(annd.median[C], annd.median[B]) # p-value = 1.129e-05  HiG vs HiGCTS
    ks.test(annd.median[C], annd.median[A]) # p-value = 3.04e-06	HiG vs CTS
    ks.test(annd.median[B], annd.median[A]) # p-value = 0.01515	HiGCTS vs CTS
    wilcox.test(annd.median[C], annd.median[B]) # p-value = 1.129e-05  HiG vs HiGCTS
    wilcox.test(annd.median[C], annd.median[A]) # p-value = 3.04e-06	HiG vs CTS
    wilcox.test(annd.median[B], annd.median[A]) # p-value = 0.004662	HiGCTS vs CTS


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
    pdf(file = paste0("annd_", db, "_v2.pdf"), height = 10)
    print(grid.arrange(pr, density_median_annd_wilcox + coord_flip(), ncol = 2, widths = c(3, 1)))
    print(grid.arrange(violin_median_annd_plot, pr, nrow = 2, heights = c(3, 3)))

    dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

############# plot normalized.strength (NO between-category difference! ) ###########
{
    df <- readRDS(file = "df_PAGERANK_strength_ANND.rewiring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
    df <- rbind(
        subset(df, PPI_cat == "CTS"),
        subset(df, PPI_cat == "HiGCTS"),
        subset(df, PPI_cat == "HiG")
    )
    df$signature <- factor(df$signature, levels = unique(df$signature))

    df$label <- df$gene
    subset(df, signature == "HiGCTS_8")

    CHD <- readRDS(file = paste0(shared_path, "CHD_Cilia_Genelist.rds"))
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
    (dim(top_genes)) # [1] 110   18  (22 graphs × 5; not auto-printed inside {})
    (dim(top_genes_CHD)) # not captured in Rscript output (inside {} block)

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
    ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value = 1.129e-05  HiG vs HiGCTS
    ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value = 0.00228	HiG vs CTS
    ks.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value =  0.09091	HiGCTS vs CTS
    wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value = 1.129e-05  HiG vs HiGCTS
    wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value =  0.0002007	HiG vs CTS
    wilcox.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value =  0.03497	HiGCTS vs CTS


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
    pdf(file = paste0("normalized.node.strength_", db, "_v2.pdf"), height = 10)
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
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_7     HiG_8
#       372       354       282       173       351       579       334       310
#    HiG_10    HiG_11    HiG_12    HiG_13    HiG_14    HiG_16    HiG_17    HiG_15
#       276       520       547       378      1043       555      1382       728
#     HiG_9 HiGCTS_15 HiGCTS_16    CTS_15    CTS_16     CTS_9
#       220        46       119        69       132       108

df_compare <- data.frame(
    signature = names(graph_list),
    n_sig.pagerank = n.pr.high,
    n_sig.annd = n.annd.high
)
df_compare$PPI_cat <- lapply(df_compare$signature, function(x) unlist(strsplit(x, split = "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))



###########################################################################################################################################
## Given a transitional state, CTS&HiG genes exhibit higher betweenness centrality in the CTS-derived network and the HiG-derived network
###########################################################################################################################################
bc <- read.table(file = "df_betweeness.tsv", header = TRUE)
(dim(bc)) # [1] 8915    5
(colnames(bc))
# [1] "signature"             "BetweennessCentrality" "gene"                  "rank_by_BC"            "PPI_cat"
## find out the cluster with CTSHiG
x <- grep("HiGCTS_", bc$signature, value = TRUE) %>% unique()
CTS <- lapply(x, function(x) unlist(strsplit(x, split = "_"))[2]) %>% unlist()
y <- grepl(".", CTS, fixed = TRUE)
if (any(y)) CT <- CTS[!y] else CT <- CTS

(CTS)
# [1] "15" "16"
(CT)
# [1] "15" "16"

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

df <- readRDS(file = "df_PAGERANK_strength_ANND.rewiring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
(colnames(df))
#  [1] "signature"                   "gene"
#  [3] "PageRank"                    "PPI_cat"
#  [5] "EigenCentrality"             "p.PageRank"
#  [7] "rank_by_p.PR"                "rank_by_PR"
#  [9] "annd"                        "p.annd"
# [11] "strength"                    "rank_by_strength"
# [13] "normalized.strength"         "rank_by_normalized.strength"
# [15] "rank_by_ANND"                "rank_by_p.ANND"

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