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

wd = "/Users/felixyu/Documents/GSE87038_weighted/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- paste0(wd, "data/")

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

db <- "GSE87038"

s <- "combined" # specificity method

########## END OF USER INPUT ##########

file <- paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_8"
# [26] "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"      "CTS_16.1"
# [31] "CTS_13"      "CTS_8"       "HiGCTS_13"
edge_counts <- sapply(graph_list, ecount)
(edge_counts)
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6 
#        4978        9544        7120        5362        5808       10993 
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18 
#        6937        8655        5828        7004       11719        8546 
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13 
#        7828        7572        9301        8632       10581        8508 
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1 
#        5880          17          15          54          10          31 
#    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16    CTS_16.1 
#          10          52          54         207          70         210 
#      CTS_13       CTS_8 
#         147         165
(sapply(graph_list, vcount))
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6 
#         303         435         411         304         320         512 
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18 
#         358         422         341         364         523         455 
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13 
#         527         336         441         406         524         403 
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1 
#         332          14          19          30          13          30 
#    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16    CTS_16.1 
#          10          31          51          66          39          79 
#      CTS_13       CTS_8 
#          60          54
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
CTS.ID <- c("7", "8", "11", "13", "15", "16", "16.1")

page <- lapply(graph_list, function(x) page_rank(x, directed = FALSE, weights = E(x)$weight)$vector)

df <- lapply(page, function(x) data.frame(PageRank = x, gene = names(x)) %>% arrange(desc(PageRank))) %>%
    rbindlist(., idcol = names(.))
colnames(df)[1] <- "signature"
df$PPI_cat <- lapply(df$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))
(dim(df)[1]) # 8213

ic <- lapply(graph_list, function(x) eigen_centrality(x, weights = E(x)$weight)$vector)
IC <- lapply(ic, function(x) data.frame(EigenCentrality = x, gene = names(x)) %>% arrange(desc(EigenCentrality))) %>%
    rbindlist(., idcol = names(.))
colnames(IC)[1] <- "signature"
(dim(IC)) # [1] 8213    3
df <- merge(df, IC, by = c("signature", "gene"))
(dim(df)) # [1] 8213    5



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
saveRDS(pr_P, file = "GSE87038_PageRank_Pvalue_by_rewiring.rds")


pr_P <- readRDS(file = "GSE87038_PageRank_Pvalue_by_rewiring.rds")
tmp <- lapply(pr_P, function(x) data.frame(p.PageRank = x, gene = names(x))) %>%
    rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"

df <- merge(df, tmp, by = c("signature", "gene")) %>%
    group_by(signature) %>%
    mutate(rank_by_p.PR = rank(p.PageRank)) %>%
    mutate(rank_by_PR = rank(-PageRank)) %>%
    ungroup()

(head(df))
#   signature gene        PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>     <chr>          <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 CTS_11    1810011O10…  0.00430 CTS             0            0.952           14
# 2 CTS_11    Aard         0.00430 CTS             0            0.956           33
# 3 CTS_11    Abhd4        0.00430 CTS             0            0.952           14
# 4 CTS_11    Akr1b8       0.00430 CTS             0            0.952           14
# 5 CTS_11    Asb4         0.0525  CTS             0.00366      0.947            6
# 6 CTS_11    Baiap2l1     0.00913 CTS             0.133        0.961           39
(dim(df)) # [1] 8213    8

(subset(df, tolower(gene) == "isl1"))
#   signature gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>     <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 CTS_8     Isl1   0.0578  CTS              0.681       0.74           45 
# 2 HiGCTS_8  Isl1   0.187   HiGCTS           0.769       0.786           3 
# 3 HiG_14    Isl1   0.00378 HiG              0.0235      0.297         129 
# 4 HiG_2     Isl1   0.00349 HiG              0.133       0.3           116.
# 5 HiG_8     Isl1   0.00703 HiG              0.242       0.237          30 

# number of significantly high pagerank per PPI_cats, too much control !
n.pr.high <- lapply(names(graph_list), function(x) nrow(subset(df, signature == x & p.PageRank < 0.05))) %>% unlist()
names(n.pr.high) <- names(graph_list)
(n.pr.high)
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6 
#           0           0           0           0           0           0 
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18 
#           0           0           0           0           0           0 
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13 
#           0           0           0           0           0           0 
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1 
#           0           0           0           0           0           0 
#    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16    CTS_16.1 
#           0           0           0           0           0           0 
#      CTS_13       CTS_8 
#           0           0 

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
saveRDS(annd_P, file = "GSE87038_annd_Pvalue_by_rewiring.rds")

annd_P <- readRDS(file = "GSE87038_annd_Pvalue_by_rewiring.rds")

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
(dim(df)) # [1] 8060    9

annd_P[["CTS_8"]]
tmp <- lapply(annd_P, function(x) data.frame(p.annd = x, gene = names(x))) %>%
    rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"
(dim(tmp)) # [1] 8213    3

df <- merge(df, tmp, by = c("signature", "gene"))
df[which(is.na(df$knn)), "p.annd"] <- NA ## due to nrow(df) 4878 > nrow(tmp)
(dim(df)) # [1] 8060   10

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
(dim(V_strength)) # 8213    6

## add the V_strength & V_strength_norm infor
df <- merge(df, V_strength, by = c("signature", "gene"))
(dim(df)) # 8060   14
(head(df, 3))
#   signature     gene    PageRank PPI_cat EigenCentrality p.PageRank
# 1    CTS_11     Asb4 0.052483905     CTS     0.003662169      0.947
# 2    CTS_11 Baiap2l1 0.009130091     CTS     0.132614464      0.961
# 3    CTS_11   Ccdc80 0.022049838     CTS     0.102663878      0.962
#   rank_by_p.PR rank_by_PR     annd p.annd   strength rank_by_strength
# 1            6          4 1.892216      0 0.20650263                5
# 2           39         30 6.457612      0 0.02879440               25
# 3           42         21 6.692259      0 0.08354285               15
#   normalized.strength rank_by_normalized.strength
# 1         0.004130053                           5
# 2         0.000575888                          25
# 3         0.001670857                          15

(table(df$normalized.strength >= df$annd))
# FALSE 
#  8060 

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
saveRDS(df, file = "df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
write.table(df, file = "df_PAGERANK_strength_ANND.rewring.P.tsv", sep = "\t", row.names = F) # !!!!!!!!

##########################
## add the column of betweenness centrality
##########################
df <- readRDS(file = "df_PAGERANK_strength_ANND.rewring.P.rds")

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

(dim(df_BC)) # [1] 8213    5

write.table(df_BC, file = "df_betweeness.tsv", sep = "\t", row.names = F) # !!!!!!!!


########### betweenness centrality ############
df_BC <- read.table(file = "df_betweeness.tsv", sep = "\t", header = T)
df_BC$PPI_cat <- factor(df_BC$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

CHD <- readRDS(file = paste0(inputdir, "CHD_Cilia_Genelist.rds"))
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
#       signature BetweennessCentrality    gene rank_by_BC PPI_cat
# 96     HiGCTS_7                     7     F10        3.0  HiGCTS
# 97     HiGCTS_7                     5   Blvrb        4.0  HiGCTS
# 98     HiGCTS_7                     1    Klf1        5.0  HiGCTS
# 99     HiGCTS_7                    12 Hbb-bh1        1.0  HiGCTS
# 100    HiGCTS_7                     9  Hba-a1        2.0  HiGCTS
# 101   HiGCTS_11                    24  Maged2        3.0  HiGCTS
# 102   HiGCTS_11                    37  Pcolce        1.0  HiGCTS
# 103   HiGCTS_11                    19    Asb4        4.0  HiGCTS
# 104   HiGCTS_11                    28  Efemp2        2.0  HiGCTS
# 105   HiGCTS_11                    18  Twist2        5.0  HiGCTS
# 106   HiGCTS_15                   129    Cd34        1.0  HiGCTS
# 107   HiGCTS_15                   116   Gng11        2.0  HiGCTS
# 108   HiGCTS_15                    51    Tie1        5.0  HiGCTS
# 109   HiGCTS_15                   114   Cldn5        3.0  HiGCTS
# 110   HiGCTS_15                    62   Gpsm3        4.0  HiGCTS
# 111   HiGCTS_16                     4  Col1a1        1.5  HiGCTS
# 112   HiGCTS_16                     4  Col5a2        1.5  HiGCTS
# 113 HiGCTS_16.1                    71   Wisp1        2.0  HiGCTS
# 114 HiGCTS_16.1                    38    Bmp2        4.5  HiGCTS
# 115 HiGCTS_16.1                    68    Gas6        3.0  HiGCTS
# 116 HiGCTS_16.1                    38    Sdpr        4.5  HiGCTS
# 117 HiGCTS_16.1                   138   Postn        1.0  HiGCTS
# 118    HiGCTS_8                     5    Fgf8        2.5  HiGCTS
# 119    HiGCTS_8                     5 Igfbpl1        2.5  HiGCTS
# 120    HiGCTS_8                     8    Isl1        1.0  HiGCTS
#     PCGC_AllCurated
# 96            FALSE
# 97            FALSE
# 98            FALSE
# 99            FALSE
# 100           FALSE
# 101           FALSE
# 102           FALSE
# 103           FALSE
# 104           FALSE
# 105           FALSE
# 106           FALSE
# 107           FALSE
# 108           FALSE
# 109           FALSE
# 110           FALSE
# 111           FALSE
# 112           FALSE
# 113           FALSE
# 114            TRUE
# 115           FALSE
# 116           FALSE
# 117           FALSE
# 118           FALSE
# 119           FALSE
# 120            TRUE

(df5_CHD <- subset(df5, PCGC_AllCurated == TRUE))

(df5_CHD)
#       signature BetweennessCentrality   gene rank_by_BC PPI_cat PCGC_AllCurated
# 6         HiG_2                 21943  Acta2        1.0     HiG            TRUE
# 14        HiG_3                  9565  Acta2        2.0     HiG            TRUE
# 51       HiG_17                 11396  Acta2        2.0     HiG            TRUE
# 58       HiG_18                  8432  Acta2        5.0     HiG            TRUE
# 62       HiG_19                 18035  Acta2        2.0     HiG            TRUE
# 73       HiG_11                 20880  Acta2        1.0     HiG            TRUE
# 84       HiG_16                 14197  Acta2        3.0     HiG            TRUE
# 91        HiG_8                  5914 Nkx2-5        3.0     HiG            TRUE
# 94        HiG_8                  7894 Cdkn1c        1.0     HiG            TRUE
# 114 HiGCTS_16.1                    38   Bmp2        4.5  HiGCTS            TRUE
# 120    HiGCTS_8                     8   Isl1        1.0  HiGCTS            TRUE
# 129      CTS_11                   110   Fbn1        3.0     CTS            TRUE
# 130      CTS_11                    96  Fgfr2        4.0     CTS            TRUE
# 140      CTS_16                    84 Plagl1        5.0     CTS            TRUE
# 153       CTS_8                   171   Isl1        5.0     CTS            TRUE
# 154       CTS_8                   213   Osr1        4.0     CTS            TRUE
# 155       CTS_8                   254  Fgfr2        2.0     CTS            TRUE



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
ks.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value =  1.12930547713156e-05  HiG vs HiGCTS
ks.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value =   3.04043782304652e-06	HiG vs CTS
ks.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value =  1	HiGCTS vs CTS
wilcox.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.000295226244881333  HiG vs HiGCTS
wilcox.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value =  0.000123880353215404	HiG vs CTS
wilcox.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = 0.440400698139003	HiGCTS vs CTS
t.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 1.31031354436124e-05  HiG vs HiGCTS
t.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value =  1.48497712000408e-05	HiG vs CTS
t.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = 0.355917683749582	HiGCTS vs CTS


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
pdf(file = "BetweennessCentrality_GSE87038_v2.pdf", height = 10)
print(grid.arrange(pr_repel, density_median_bc_plot + coord_flip(), ncol = 2, widths = c(3, 1)))
print(grid.arrange(violin_median_bc_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_t, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_median_bc_wilcox, violin_median_bc_wilcox_ln, nrow = 2))
dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


########### plot PageRank ############
df <- readRDS(file = "df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
(dim(df)) # [1] 8060   16

df <- rbind(
    subset(df, PPI_cat == "CTS"),
    subset(df, PPI_cat == "HiGCTS"),
    subset(df, PPI_cat == "HiG")
)

CHD <- readRDS(file = paste0(inputdir, "CHD_Cilia_Genelist.rds"))
df$PCGC_AllCurated <- toupper(df$gene) %in% toupper(unlist(CHD["Griffin2023_PCGC_AllCurated"]))

# Calculate top 5 significant genes within each box
df5 <- df %>%
    filter(rank_by_PR <= 5) %>%
    ungroup()
(subset(df5, PPI_cat == "HiGCTS"))
# A tibble: 30 × 17
#    signature gene   PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#    <chr>     <chr>     <dbl> <fct>             <dbl>      <dbl>        <dbl>
#  1 HiGCTS_11 Asb4     0.119  HiGCTS           0.0941      0.958          4.5
#  2 HiGCTS_11 Ccdc80   0.0742 HiGCTS           0.608       0.967         10  
#  3 HiGCTS_11 Efemp2   0.0824 HiGCTS           0.612       0.969         12  
#  4 HiGCTS_11 Pcolce   0.188  HiGCTS           1           0.963          7  
#  5 HiGCTS_11 Sgce     0.0711 HiGCTS           0.0706      0.961          6  
#  6 HiGCTS_15 Cd34     0.115  HiGCTS           1           0.721         20  
#  7 HiGCTS_15 Cdh5     0.113  HiGCTS           0.996       0.711         19  
#  8 HiGCTS_15 Cldn5    0.0851 HiGCTS           0.826       0.734         22  
#  9 HiGCTS_15 Gng11    0.0715 HiGCTS           0.285       0.587          5  
# 10 HiGCTS_15 Tie1     0.0520 HiGCTS           0.507       0.763         24  
#    rank_by_PR  annd p.annd strength rank_by_strength normalized.strength
#         <dbl> <dbl>  <dbl>    <dbl>            <dbl>               <dbl>
#  1          2  1.24      0    0.171                2             0.00953
#  2          4  4.80      0    0.103                5             0.00572
#  3          3  4.28      0    0.119                3             0.00660
#  4          1  2.03      0    0.269                1             0.0150 
#  5          5  3         0    0.103                4             0.00574
#  6          1  8.22      0    1.26                 1             0.0433 
#  7          2  6.75      0    1.23                 2             0.0425 
#  8          3  8.13      0    0.938                3             0.0323 
#  9          4  5.61      0    0.625                4             0.0215 
# 10          5  7.64      0    0.554                5             0.0191 
#    rank_by_normalized.strength rank_by_ANND rank_by_p.ANND PCGC_AllCurated
#                          <int>        <dbl>          <dbl> <lgl>          
#  1                           2           12            7.5 FALSE          
#  2                           5            4            7.5 FALSE          
#  3                           3            5            7.5 FALSE          
#  4                           1           11            7.5 FALSE          
#  5                           4            8            7.5 FALSE          
#  6                           1            5           12.5 FALSE          
#  7                           2           11           12.5 FALSE          
#  8                           3            7           12.5 FALSE          
#  9                           4           15           12.5 FALSE          
# 10                           5            8           12.5 FALSE    

write.table(df5[, c(
    "signature", "gene", "PageRank", "PPI_cat", "rank_by_PR",
    "normalized.strength", "rank_by_normalized.strength", "PCGC_AllCurated"
)], file = "table_top5_PageRank_perPPI.tsv", sep = "\t", row.names = FALSE, quote = FALSE) # !!!!!!!!!!!!!!

df5_CHD <- subset(df5, PCGC_AllCurated == TRUE)
(dim(df5)) # [1] 160  17
(dim(df5_CHD)) # [1] 15 17
(df5_CHD %>% as.data.frame())
#    signature   gene    PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
# 1     CTS_11   Fbn1 0.055065672     CTS       0.4869366      0.973         30.5
# 2     CTS_11  Fgfr2 0.064327535     CTS       1.0000000      0.960         17.5
# 3      CTS_8  Fgfr2 0.041574332     CTS       0.3937628      0.458         13.0
# 4      CTS_8   Isl1 0.057804144     CTS       0.6809899      0.740         42.0
# 5  HiGCTS_16 Plagl1 0.114285714  HiGCTS       0.0000000      0.929          2.0
# 6   HiGCTS_8   Isl1 0.187381852  HiGCTS       0.7691712      0.786          3.0
# 7     HiG_11  Acta2 0.013871703     HiG       0.9973502      0.534        308.0
# 8     HiG_16  Acta2 0.011768339     HiG       1.0000000      0.770        517.0
# 9     HiG_17  Acta2 0.009240550     HiG       0.7519698      0.627        451.5
# 10    HiG_17  Actc1 0.008124728     HiG       0.8657376      0.723        498.5
# 11    HiG_18  Acta2 0.010387360     HiG       0.2678098      0.549        366.0
# 12    HiG_19  Acta2 0.011703470     HiG       1.0000000      0.704        519.0
# 13     HiG_2  Acta2 0.013565393     HiG       0.7456318      0.720        377.0
# 14     HiG_2  Actc1 0.010526542     HiG       0.8855112      0.810        402.0
# 15     HiG_3  Acta2 0.011665908     HiG       0.4713817      0.527        282.5
#    rank_by_PR      annd p.annd   strength rank_by_strength normalized.strength
# 1         3.0  5.873181      0  0.2466545              3.0         0.004933090
# 2         1.0  5.606124      0  0.3290442              1.0         0.006580885
# 3         4.0 12.525672      0  0.8116738              7.0         0.015314599
# 4         2.0 12.579043      0  1.4398707              2.0         0.027167371
# 5         4.5  1.000000      0  0.1212670              4.5         0.010105587
# 6         2.0  3.485632      0  0.4855130              2.0         0.053945886
# 7         1.0 74.216718      0  4.3194136              1.0         0.009816849
# 8         2.0 84.605783      0  4.5640296              1.0         0.008726634
# 9         2.0 96.246639      0  3.3943301              3.0         0.006502548
# 10        5.0 93.382448      0  3.1545475              5.0         0.006043194
# 11        4.0 77.000112      0  1.5957888              5.0         0.003514953
# 12        2.0 70.456011      0  2.1186153              1.0         0.004027786
# 13        2.0 77.610174      0 11.3800587              2.0         0.026221333
# 14        5.0 75.379152      0  9.8349919              4.0         0.022661272
# 15        4.0 65.374572      0  5.1496331              3.0         0.012560081
#    rank_by_normalized.strength rank_by_ANND rank_by_p.ANND PCGC_AllCurated
# 1                            3          8.0           16.5            TRUE
# 2                            1         11.0           16.5            TRUE
# 3                            7         21.0           24.5            TRUE
# 4                            2         20.0           24.5            TRUE
# 5                            5          7.5            4.5            TRUE
# 6                            2          3.0            4.0            TRUE
# 7                            1        198.0          219.0            TRUE
# 8                            1        151.0          260.0            TRUE
# 9                            3         96.0          261.5            TRUE
# 10                           5        122.0          261.5            TRUE
# 11                           5        134.0          226.0            TRUE
# 12                           1        118.0          262.0            TRUE
# 13                           2        222.0          218.0            TRUE
# 14                           4        245.0          218.0            TRUE
# 15                           3        170.0          203.5            TRUE

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
ks.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value =  1.12930547713156e-05  HiG vs HiGCTS
ks.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  3.04043782304652e-06	HiG vs CTS
ks.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.00116550116550117	HiGCTS vs CTS
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value = 1.12930547713156e-05  HiG vs HiGCTS
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  3.04043782304652e-06 	HiG vs CTS
wilcox.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.00116550116550117	HiGCTS vs CTS


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
    df <- readRDS(file = "df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!

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
    CHD <- readRDS(file = paste0(inputdir, "CHD_Cilia_Genelist.rds"))

    top_genes_CHD <- subset(top_genes, toupper(gene) %in% toupper(unlist(CHD[c("Griffin2023_PCGC_AllCurated")])))
    (dim(top_genes)) # [1] 160  17
    (dim(top_genes_CHD)) # [1]  7 17
    (top_genes_CHD)
    # # A tibble: 7 × 17
    # # Groups:   signature [7]
    #   signature   gene   PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
    #   <fct>       <chr>     <dbl> <fct>             <dbl>      <dbl>        <dbl>
    # 1 CTS_16.1    Tbx20  0.00446  CTS              0.0677      0.752         57  
    # 2 CTS_8       Tbx1   0.0134   CTS              0.222       0.753         43  
    # 3 HiGCTS_16.1 Bmp2   0.0343   HiGCTS           0.207       0.902          6.5
    # 4 HiGCTS_8    Isl1   0.187    HiGCTS           0.769       0.786          3  
    # 5 HiG_18      Ece1   0.000802 HiG              0.0118      0.474        217  
    # 6 HiG_5       Cited2 0.00131  HiG              0.0288      0.432        192. 
    # 7 HiG_7       Cited2 0.000970 HiG              0.0276      0.223         68  
    #   rank_by_PR   annd p.annd strength rank_by_strength normalized.strength
    #        <dbl>  <dbl>  <dbl>    <dbl>            <dbl>               <dbl>
    # 1         60  21.4       0   0.0374               55            0.000479
    # 2         25  18.1       0   0.296                21            0.00558 
    # 3         12   5.40      0   0.0920               12            0.00317 
    # 4          2   3.49      0   0.486                 2            0.0539  
    # 5        357 128.        0   0.0467              370            0.000103
    # 6        236  91.8       0   0.155               231            0.000487
    # 7        286 103.        0   0.0951              273            0.000284
    #   rank_by_normalized.strength rank_by_ANND rank_by_p.ANND label 
    #                         <int>        <dbl>          <dbl> <chr> 
    # 1                          55            3           34   Tbx20 
    # 2                          21            3           24.5 Tbx1  
    # 3                          12            3           11.5 Bmp2  
    # 4                           2            3            4   Isl1  
    # 5                         370            1          226   Ece1  
    # 6                         231            5          158.  Cited2
    # 7                         273            2          168   Cited2

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
    pdf(file = "annd_GSE87038_v2.pdf", height = 10)
    print(grid.arrange(pr, density_median_annd_wilcox + coord_flip(), ncol = 2, widths = c(3, 1)))
    print(grid.arrange(violin_median_annd_plot, pr, nrow = 2, heights = c(3, 3)))

    dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

############# plot normalized.strength (NO between-category difference! ) ###########
{
    df <- readRDS(file = "df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
    df <- rbind(
        subset(df, PPI_cat == "CTS"),
        subset(df, PPI_cat == "HiGCTS"),
        subset(df, PPI_cat == "HiG")
    )
    df$signature <- factor(df$signature, levels = unique(df$signature))

    df$label <- df$gene
    subset(df, signature == "HiGCTS_8")

    CHD <- readRDS(file = paste0(inputdir, "CHD_Cilia_Genelist.rds"))
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
    (dim(top_genes)) # [1] 160   18
    (dim(top_genes_CHD)) # [1] 15  18

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
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6 
#         300         435         406         301         316         510 
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18 
#         355         418         338         358         522         451 
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13 
#         523         335         437         406         519         403 
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1 
#         329           9          14          24           8          22 
#    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16    CTS_16.1 
#           7          23          32          60          34          67 
#      CTS_13       CTS_8 
#          50          48

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



###########################################################################################################################################
## Given a transitional state, CTS&HiG genes exhibit higher betweenness centrality in the CTS-derived network and the HiG-derived network
###########################################################################################################################################
bc <- read.table(file = "df_betweeness.tsv", header = TRUE)
(dim(bc)) # [1] 8213    5
(colnames(bc))
# [1] "signature"             "BetweennessCentrality" "gene"                  "rank_by_BC"            "PPI_cat"
## find out the cluster with CTSHiG
x <- grep("HiGCTS_", bc$signature, value = TRUE) %>% unique()
CTS <- lapply(x, function(x) unlist(strsplit(x, split = "_"))[2]) %>% unlist()
y <- grepl(".", CTS, fixed = TRUE)
if (any(y)) CT <- CTS[!y] else CT <- CTS

(CTS)
# "7"    "11"   "15"   "16"   "16.1" "8"    "13"
(CT)
# "7"  "11" "15" "16" "8"  "13"

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

df <- readRDS(file = "df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
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