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

wd = "/Users/felixyu/Documents/IbarraSoria2018/"

setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- paste0(wd, "data/")

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

db <- "IbarraSoria2018"

s <- "combined" # specificity method

########## END OF USER INPUT ##########

file <- paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [5] "HiG_blood"                  "HiG_mesodermProgenitors"    "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"  
#  [9] "HiG_somiticMesoderm"        "HiG_mixedMesoderm.a"        "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"              "HiG_endothelial.b"          "HiG_cardiac.a"             
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"           "CTS_endothelial.b"          "CTS_cardiac.a"
edge_counts <- sapply(graph_list, ecount)
(edge_counts)
# HiG_extraembryonicMesoderm          HiG_endothelial.a          HiG_endothelial.c          HiG_endothelial.d 
#                       6743                       6743                      15247                      15247 
#                  HiG_blood    HiG_mesodermProgenitors   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                      11724                      11724                      12662                      12662 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                      20602                      20602                       7557                       7557 
#              HiG_cardiac.b              HiG_cardiac.c          HiG_endothelial.b              HiG_cardiac.a 
#                      12912                      12912                       6608                       6608 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a          CTS_endothelial.b              CTS_cardiac.a 
#                         20                         28                         82                         79

(sapply(graph_list, vcount))
# HiG_extraembryonicMesoderm          HiG_endothelial.a          HiG_endothelial.c          HiG_endothelial.d 
#                        414                        414                        573                        573 
#                  HiG_blood    HiG_mesodermProgenitors   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                        559                        559                        537                        537 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                        552                        552                        432                        432 
#              HiG_cardiac.b              HiG_cardiac.c          HiG_endothelial.b              HiG_cardiac.a 
#                        549                        549                        392                        392 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a          CTS_endothelial.b              CTS_cardiac.a 
#                         16                         13                         33                         37 

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
## PageRank’s main difference from EigenCentrality is that it accounts for link direction. Like EigenCentrality,
## PageRank can help uncover influential or important nodes whose reach extends beyond just their direct connections.
## It’s especially useful in scenarios where link direction is important:
#
# https://igraph.org/r/doc/page_rank.html
CTS.ID <- c("endothelial.b", "cardiac.a")

page <- lapply(graph_list, function(x) page_rank(x, directed = FALSE, weights = E(x)$weight)$vector)

df <- lapply(page, function(x) data.frame(PageRank = x, gene = names(x)) %>% arrange(desc(PageRank))) %>%
    rbindlist(., idcol = names(.))
colnames(df)[1] <- "signature"
df$PPI_cat <- lapply(df$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))
(dim(df)[1]) # 8115

ic <- lapply(graph_list, function(x) eigen_centrality(x, weights = E(x)$weight)$vector)
IC <- lapply(ic, function(x) data.frame(EigenCentrality = x, gene = names(x)) %>% arrange(desc(EigenCentrality))) %>%
    rbindlist(., idcol = names(.))
colnames(IC)[1] <- "signature"
(dim(IC)) # [1] 8115    3
df <- merge(df, IC, by = c("signature", "gene"))
(dim(df)) # [1] 8115    5



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
    #      PageRank    gene
    #         <num>  <char>
    # 1: 0.01588449   Grwd1
    # 2: 0.01443920    Imp4
    # 3: 0.01325567 Mybbp1a
    # 4: 0.01302575    Atic
    # 5: 0.01295354   Hmgn2
    # 6: 0.01283500 Rsl24d1
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
#   signature     gene     PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR
#   <chr>         <chr>       <dbl> <fct>             <dbl>      <dbl>        <dbl>      <dbl>
# 1 CTS_cardiac.a Anxa2     0.0212  CTS             0.0207       0.626         19           18
# 2 CTS_cardiac.a Arhgap29  0.0252  CTS             0.00623      0.614         17.5         15
# 3 CTS_cardiac.a Arl4c     0.00483 CTS             0            0.766         37           34
# 4 CTS_cardiac.a Cgnl1     0.00483 CTS             0            0.753         32           34
# 5 CTS_cardiac.a Dact1     0.00836 CTS             0.00467      0.582          9.5         27
# 6 CTS_cardiac.a Dbt       0.00483 CTS             0            0.749         28           34
(dim(df)) # [1] 8115    8

(subset(df, tolower(gene) == "isl1"))
#   signature                  gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR
#   <chr>                      <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>      <dbl>
# 1 CTS_cardiac.a              Isl1   0.0440  CTS            0.597         0.756         33            9
# 2 HiGCTS_cardiac.a           Isl1   0.0914  HiGCTS         0.691         0.669          8.5          6
# 3 HiG_endothelial.a          Isl1   0.00166 HiG            0.000399      0.118         31          178
# 4 HiG_extraembryonicMesoderm Isl1   0.00194 HiG            0.00104       0.121         20          167

# number of significantly high pagerank per PPI_cats, too much control !
n.pr.high <- lapply(names(graph_list), function(x) nrow(subset(df, signature == x & p.PageRank < 0.05))) %>% unlist()
names(n.pr.high) <- names(graph_list)
(n.pr.high)
# HiG_extraembryonicMesoderm          HiG_endothelial.a          HiG_endothelial.c          HiG_endothelial.d 
#                          0                          0                         19                         16 
#                  HiG_blood    HiG_mesodermProgenitors   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                          0                          0                         11                         11 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                         28                         33                          4                          6 
#              HiG_cardiac.b              HiG_cardiac.c          HiG_endothelial.b              HiG_cardiac.a 
#                         12                         13                          3                          5 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a          CTS_endothelial.b              CTS_cardiac.a 
#                          0                          0                          0                          0 


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
(any(is.na(annd_observed[["CTS_cardiac.a"]]))) # FALSE
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
    #         annd     gene
    #        <num>   <char>
    # 1: 109.97744 Eif4ebp1
    # 2:  92.55861    Smyd2
    # 3:  90.53871     Plp2

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
(dim(df)) # [1] 8067    9

annd_P[["CTS_cardiac.a"]]
tmp <- lapply(annd_P, function(x) data.frame(p.annd = x, gene = names(x))) %>%
    rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"
(dim(tmp)) # [1] 8115    3

df <- merge(df, tmp, by = c("signature", "gene"))
df[which(is.na(df$knn)), "p.annd"] <- NA ## due to nrow(df) 4878 > nrow(tmp)
(dim(df)) # [1] 8067   10

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
(dim(V_strength)) # 8115    6

## add the V_strength & V_strength_norm info
df <- merge(df, V_strength, by = c("signature", "gene"))
(dim(df)) # 8067   14
(head(df, 3))
#       signature     gene    PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR     annd p.annd
# 1 CTS_cardiac.a    Anxa2 0.021185615     CTS     0.020682887      0.626         19.0         18 5.085841      0
# 2 CTS_cardiac.a Arhgap29 0.025220980     CTS     0.006230479      0.614         17.5         15 4.740786      0
# 3 CTS_cardiac.a    Dact1 0.008357241     CTS     0.004668790      0.582          9.5         27 4.125428      0
#    strength rank_by_strength normalized.strength rank_by_normalized.strength
# 1 0.1617582               17        0.0044932825                          17
# 2 0.1787239               16        0.0049645531                          16
# 3 0.0270023               27        0.0007500638                          27

(table(df$normalized.strength >= df$annd))
# FALSE 
#  8067 

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
#  [1] "signature"                   "gene"                        "PageRank"                   
#  [4] "PPI_cat"                     "EigenCentrality"             "p.PageRank"                 
#  [7] "rank_by_p.PR"                "rank_by_PR"                  "annd"                       
# [10] "p.annd"                      "strength"                    "rank_by_strength"           
# [13] "normalized.strength"         "rank_by_normalized.strength" "rank_by_ANND"               
# [16] "rank_by_p.ANND"  
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

(dim(df_BC)) # 8115    5

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
#               signature BetweennessCentrality  gene rank_by_BC PPI_cat PCGC_AllCurated
# 80 HiGCTS_endothelial.b                     9  Hhex        3.5  HiGCTS           FALSE
# 81 HiGCTS_endothelial.b                    17 Sox17        1.0  HiGCTS            TRUE
# 82 HiGCTS_endothelial.b                     8  Tal1        5.0  HiGCTS           FALSE
# 83 HiGCTS_endothelial.b                     9  Rhoj        3.5  HiGCTS           FALSE
# 84 HiGCTS_endothelial.b                    16  Etv2        2.0  HiGCTS           FALSE
# 85     HiGCTS_cardiac.a                    10  Klf6        4.0  HiGCTS           FALSE
# 86     HiGCTS_cardiac.a                    10 Sfrp5        4.0  HiGCTS           FALSE
# 87     HiGCTS_cardiac.a                    16  Msx2        2.0  HiGCTS            TRUE
# 88     HiGCTS_cardiac.a                    10 Gata4        4.0  HiGCTS            TRUE
# 89     HiGCTS_cardiac.a                    28 Mef2c        1.0  HiGCTS            TRUE

(df5_CHD <- subset(df5, PCGC_AllCurated == TRUE))
(df5_CHD)
#                 signature BetweennessCentrality  gene rank_by_BC PPI_cat PCGC_AllCurated
# 15      HiG_endothelial.d              9735.000 Sox17          5     HiG            TRUE
# 41    HiG_somiticMesoderm              9943.833 Rps19          4     HiG            TRUE
# 46    HiG_mixedMesoderm.a             11241.333 Rps19          4     HiG            TRUE
# 53 HiG_pharyngealMesoderm             10986.000   Eed          2     HiG            TRUE
# 58    HiG_mixedMesoderm.b             10749.000   Eed          2     HiG            TRUE
# 81   HiGCTS_endothelial.b                17.000 Sox17          1  HiGCTS            TRUE
# 87       HiGCTS_cardiac.a                16.000  Msx2          2  HiGCTS            TRUE
# 88       HiGCTS_cardiac.a                10.000 Gata4          4  HiGCTS            TRUE
# 89       HiGCTS_cardiac.a                28.000 Mef2c          1  HiGCTS            TRUE
# 91      CTS_endothelial.b                47.000 Sox17          5     CTS            TRUE
# 97          CTS_cardiac.a               109.000  Msx2          3     CTS            TRUE
# 99          CTS_cardiac.a               214.000 Mef2c          1     CTS            TRUE


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
ks.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.366013071895425
ks.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 0.366013071895425
ks.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value =  1
wilcox.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.162552609232464
wilcox.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 0.162552609232464
wilcox.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = NAN
t.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.000672926504386247
t.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 0.000672926504386247
t.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = NAN

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
pdf(file = paste0("BetweennessCentrality_", db,"_v2.pdf"), height = 10)
print(grid.arrange(pr_repel, density_median_bc_plot + coord_flip(), ncol = 2, widths = c(3, 1)))
print(grid.arrange(violin_median_bc_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_wilcox, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_t, pr, nrow = 2, heights = c(3, 3)))
print(grid.arrange(violin_median_bc_wilcox, violin_median_bc_wilcox_ln, nrow = 2))
dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


########### plot PageRank ############
df <- readRDS(file = "df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
(dim(df)) # [1] 8067   16

## reorder df$df$signature to be #E7298A #7570B3 and #E6AB02
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
# # A tibble: 10 × 17
#    signature            gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR  annd p.annd strength
#    <chr>                <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>      <dbl> <dbl>  <dbl>    <dbl>
#  1 HiGCTS_cardiac.a     Gata4   0.130  HiGCTS            0.715      0.53             3          2  6.31      0    0.747
#  2 HiGCTS_cardiac.a     Gata6   0.0939 HiGCTS            0.542      0.583            5          5  6.40      0    0.540
#  3 HiGCTS_cardiac.a     Mef2c   0.188  HiGCTS            1          0.576            4          1  5.69      0    1.14 
#  4 HiGCTS_cardiac.a     Msx2    0.119  HiGCTS            0.702      0.649            7          3  6.25      0    0.720
#  5 HiGCTS_cardiac.a     Sfrp5   0.0947 HiGCTS            0.393      0.512            1          4  5.15      0    0.511
#  6 HiGCTS_endothelial.b Etv2    0.155  HiGCTS            0.832      0.693            8          3  5.58      0    0.920
#  7 HiGCTS_endothelial.b Gata2   0.167  HiGCTS            0.864      0.745            9          2  4.62      0    1.01 
#  8 HiGCTS_endothelial.b Hhex    0.136  HiGCTS            0.767      0.76            11          4  5.60      0    0.829
#  9 HiGCTS_endothelial.b Sox17   0.108  HiGCTS            0.484      0.646            3          5  5.02      0    0.558
# 10 HiGCTS_endothelial.b Tal1    0.180  HiGCTS            1          0.75            10          1  5.67      0    1.15 
#    rank_by_strength normalized.strength rank_by_normalized.strength rank_by_ANND rank_by_p.ANND PCGC_AllCurated
#               <dbl>               <dbl>                       <int>        <dbl>          <dbl> <lgl>          
#  1                2              0.0622                           2            6            6.5 TRUE           
#  2                5              0.0450                           5            5            6.5 TRUE           
#  3                1              0.0946                           1            9            6.5 TRUE           
#  4                3              0.0600                           3            7            6.5 TRUE           
#  5                6              0.0426                           6           10            6.5 FALSE          
#  6                3              0.0613                           3            7            6   FALSE          
#  7                2              0.0670                           2            9            6   FALSE          
#  8                4              0.0552                           4            6            6   FALSE          
#  9                5              0.0372                           5            8            6   TRUE           
# 10                1              0.0765                           1            5            6   FALSE   

write.table(df5[, c(
    "signature", "gene", "PageRank", "PPI_cat", "rank_by_PR",
    "normalized.strength", "rank_by_normalized.strength", "PCGC_AllCurated"
)], file = "table_top5_PageRank_perPPI.tsv", sep = "\t", row.names = FALSE, quote = FALSE) # !!!!!!!!!!!!!!

df5_CHD <- subset(df5, PCGC_AllCurated == TRUE)
(dim(df5)) # 100  17
(dim(df5_CHD)) # 15 17
(df5_CHD %>% as.data.frame())
#                 signature  gene    PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR       annd p.annd   strength rank_by_strength normalized.strength rank_by_normalized.strength
# 1           CTS_cardiac.a Gata4 0.069606922     CTS      0.74310419      0.714         23.0          3   9.527368      0  0.9704233                3          0.02695620                           3
# 2           CTS_cardiac.a Gata6 0.053805790     CTS      0.54934872      0.722         24.0          5   9.076687      0  0.7298393                5          0.02027331                           5
# 3           CTS_cardiac.a Mef2c 0.118217061     CTS      1.00000000      0.713         22.0          1   7.779273      0  1.6446475                1          0.04568465                           1
# 4           CTS_cardiac.a  Msx2 0.061308985     CTS      0.61452550      0.761         30.0          4   8.789915      0  0.8951569                4          0.02486547                           4
# 5           CTS_cardiac.a  Tbx5 0.070150905     CTS      0.85442984      0.744         27.0          2  10.813689      0  1.0522409                2          0.02922891                           2
# 6        HiGCTS_cardiac.a Gata4 0.129870192  HiGCTS      0.71502275      0.530          3.0          2   6.312456      0  0.7467978                2          0.06223315                           2
# 7        HiGCTS_cardiac.a Gata6 0.093866630  HiGCTS      0.54225187      0.583          5.0          5   6.398557      0  0.5397180                5          0.04497650                           5
# 8        HiGCTS_cardiac.a Mef2c 0.187525656  HiGCTS      1.00000000      0.576          4.0          1   5.691685      0  1.1353058                1          0.09460882                           1
# 9        HiGCTS_cardiac.a  Msx2 0.119106075  HiGCTS      0.70190349      0.649          7.0          3   6.252759      0  0.7200779                3          0.06000649                           3
# 10   HiGCTS_endothelial.b Sox17 0.108237782  HiGCTS      0.48403374      0.646          3.0          5   5.023638      0  0.5584749                5          0.03723166                           5
# 11      HiG_endothelial.d Sox17 0.010926074     HiG      0.01487872      0.000          2.5          2  74.264148      0 25.7280000               50          0.04497902                          50
# 12    HiG_mixedMesoderm.a Rps19 0.008973646     HiG      0.82133962      0.003         10.0          5 130.184987      0 80.4920000                7          0.14608348                           7
# 13    HiG_mixedMesoderm.b   Eed 0.020582951     HiG      0.55868063      0.123         17.0          2  62.511422      0 41.5000000                3          0.09628770                           3
# 14 HiG_pharyngealMesoderm   Eed 0.020991274     HiG      0.55800736      0.100         14.0          2  62.511422      0 41.5000000                3          0.09628770                           3
# 15    HiG_somiticMesoderm Rps19 0.009070583     HiG      0.81991280      0.003          8.5          5 130.184987      0 80.4920000                7          0.14608348                           7
#    rank_by_ANND rank_by_p.ANND PCGC_AllCurated
# 1             8           15.5            TRUE
# 2             9           15.5            TRUE
# 3            15           15.5            TRUE
# 4            10           15.5            TRUE
# 5             4           15.5            TRUE
# 6             6            6.5            TRUE
# 7             5            6.5            TRUE
# 8             9            6.5            TRUE
# 9             7            6.5            TRUE
# 10            8            6.0            TRUE
# 11          337          287.0            TRUE
# 12           92          276.0            TRUE
# 13          250          214.5            TRUE
# 14          250          214.5            TRUE
# 15           95          276.0            TRUE

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
ks.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value =  0.0130718954248366
ks.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.0130718954248366
ks.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.333333333333333
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value = 0.0130718954248366
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.0130718954248366
wilcox.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.333333333333333

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
# Add statistical comparisons using stat_compare_means manually

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
    df <- readRDS(file = "df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!

    df$label <- df$gene
    subset(df, signature == "HiGCTS_cardiac.a")

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
    (dim(top_genes)) # 100  17
    (dim(top_genes_CHD)) # 8 17
    (top_genes_CHD)
    # # A tibble: 8 × 17
    # # Groups:   signature [5]
    #   signature             gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR   annd p.annd strength
    #   <fct>                 <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>      <dbl>  <dbl>  <dbl>    <dbl>
    # 1 CTS_cardiac.a         Isl1  0.0440   CTS             0.597        0.756         29            9  11.2       0    0.652
    # 2 CTS_cardiac.a         Tbx5  0.0702   CTS             0.854        0.744         27            2  10.8       0    1.05 
    # 3 HiGCTS_cardiac.a      Isl1  0.0914   HiGCTS          0.691        0.669          8.5          6   7.18      0    0.557
    # 4 HiGCTS_cardiac.a      Twis… 0.0735   HiGCTS          0.527        0.669          8.5          8   6.68      0    0.435
    # 5 HiGCTS_cardiac.a      Gata6 0.0939   HiGCTS          0.542        0.583          5            5   6.40      0    0.540
    # 6 HiG_endothelial.c     Cite… 0.000390 HiG             0.00199      0.318        237          519 151.        0    0.388
    # 7 HiG_endothelial.d     Cite… 0.000440 HiG             0.00182      0.223        106.         472 151.        0    0.334
    # 8 HiG_presomiticMesode… Cite… 0.000497 HiG             0.00530      0.282        144.         441 151.        0    0.601
    # # ℹ 6 more variables: rank_by_strength <dbl>, normalized.strength <dbl>, rank_by_normalized.strength <int>,
    # #   rank_by_ANND <dbl>, rank_by_p.ANND <dbl>, label <chr>

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
    ks.test(annd.median[C], annd.median[B]) # p-value = 0.0130718954248366
    ks.test(annd.median[C], annd.median[A]) # p-value = 0.0130718954248366
    ks.test(annd.median[B], annd.median[A]) # p-value = 0.333333333333333
    wilcox.test(annd.median[C], annd.median[B]) # p-value = 0.0130718954248366
    wilcox.test(annd.median[C], annd.median[A]) # p-value = 0.0130718954248366
    wilcox.test(annd.median[B], annd.median[A]) # p-value = 0.333333333333333

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
    df <- readRDS(file = "df_PAGERANK_strength_ANND.rewring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
    df <- rbind(
        subset(df, PPI_cat == "CTS"),
        subset(df, PPI_cat == "HiGCTS"),
        subset(df, PPI_cat == "HiG")
    )
    df$signature <- factor(df$signature, levels = unique(df$signature))

    df$label <- df$gene
    subset(df, signature == "HiGCTS_cardiac.a")

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
    (dim(top_genes)) # 100  18
    (dim(top_genes_CHD)) # 13 18

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
    ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value = 0.07843137
    ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value = 0.3660131
    ks.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value = 0.3333333 
    wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value = 0.05228758 
    wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value = 0.3921569
    wilcox.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value = 0.3333333 

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
        ggtitle("wilcox, medina nr_strength")

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
# HiG_extraembryonicMesoderm          HiG_endothelial.a          HiG_endothelial.c          HiG_endothelial.d 
#                        409                        409                        573                        573 
#                  HiG_blood    HiG_mesodermProgenitors   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                        557                        557                        536                        536 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                        551                        551                        428                        428 
#              HiG_cardiac.b              HiG_cardiac.c          HiG_endothelial.b              HiG_cardiac.a 
#                        548                        548                        391                        391 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a          CTS_endothelial.b              CTS_cardiac.a 
#                         11                         12                         28                         30  

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

# dev.copy2pdf(file='n.sig.pageRank_vs_n.sig.annd.pdf')



###########################################################################################################################################
## Given a transitional state, CTS&HiG genes exhibit higher betweenness centrality in the CTS-derived network and the HiG-derived network
###########################################################################################################################################
bc <- read.table(file = "df_betweeness.tsv", header = TRUE)
(dim(bc)) # 8115    5
(colnames(bc))
# [1] "signature"             "BetweennessCentrality" "gene"                  "rank_by_BC"            "PPI_cat"
## find out the cluster with CTSHiG
x <- grep("HiGCTS_", bc$signature, value = TRUE) %>% unique()
CTS <- lapply(x, function(x) unlist(strsplit(x, split = "_"))[2]) %>% unlist()
CT <- sub("\\..*", "", CTS)

(CTS)
# "endothelial.b" "cardiac.a"
(CT)
# "endothelial" "cardiac" 

plot_CTS_bc <- plot_HiG_bc <- list()
## compare the CTS&HiG genes to non-hiG genes per CTS-derived PPI
for (id in CTS) {

    if (grepl(".", id, fixed = TRUE)) {
        id2 <- sub("\\..*", "", id)  # strip suffix after dot
    } else {
        id2 <- id
    }

    # CTS cluster
    CTS_PPI <- subset(bc, signature == paste("CTS", id, sep = "_"))

    # HiG clusters containing this CTS id (handles suffixes like .a, .b, etc.)
    HiG_PPI <- subset(bc, grepl(id2, signature) & grepl("^HiG_", signature))

    # Label genes as in HiG or not
    CTS_PPI$isHiG <- factor(CTS_PPI$gene %in% HiG_PPI$gene, levels = c(FALSE, TRUE))
    HiG_PPI$isCTS <- factor(HiG_PPI$gene %in% CTS_PPI$gene, levels = c(FALSE, TRUE))



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
        mutate(gene = factor(gene, levels = unique(gene[order(-BetweennessCentrality)])))

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
    if (grepl(".", id, fixed = TRUE)) {
        id2 <- sub("\\..*", "", id)  # strip suffix after dot
    } else {
        id2 <- id
    }

    # CTS cluster
    CTS_PPI <- subset(df, signature == paste("CTS", id, sep = "_"))

    # HiG clusters containing this CTS id (handles suffixes like .a, .b, etc.)
    HiG_PPI <- subset(df, grepl(id2, signature) & grepl("^HiG_", signature))

    # Label genes as in HiG or not
    CTS_PPI$isHiG <- factor(CTS_PPI$gene %in% HiG_PPI$gene, levels = c(FALSE, TRUE))
    HiG_PPI$isCTS <- factor(HiG_PPI$gene %in% CTS_PPI$gene, levels = c(FALSE, TRUE))

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
        mutate(gene = factor(gene, levels = unique(gene[order(-PageRank)])))

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
    if (grepl(".", id, fixed = TRUE)) {
        id2 <- sub("\\..*", "", id)  # strip suffix after dot
    } else {
        id2 <- id
    }

    # CTS cluster
    CTS_PPI <- subset(df, signature == paste("CTS", id, sep = "_"))

    # HiG clusters containing this CTS id (handles suffixes like .a, .b, etc.)
    HiG_PPI <- subset(df, grepl(id2, signature) & grepl("^HiG_", signature))

    # Label genes as in HiG or not
    CTS_PPI$isHiG <- factor(CTS_PPI$gene %in% HiG_PPI$gene, levels = c(FALSE, TRUE))
    HiG_PPI$isCTS <- factor(HiG_PPI$gene %in% CTS_PPI$gene, levels = c(FALSE, TRUE))
    ############## ranked by annd
    # Compute exact p-value
    pval <- wilcox.test(annd ~ isHiG, data = CTS_PPI)$p.value

    # Reorder gene factor levels by annd (high to low)
    CTS_PPI <- CTS_PPI %>%
        mutate(gene = factor(gene, levels = unique(gene[order(-annd)])))

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
        mutate(gene = factor(gene, levels = unique(gene[order(-annd)])))

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

