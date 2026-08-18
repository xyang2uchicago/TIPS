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

source(here::here("examples", "config.R"))
wd <- tips_path("examples", "hematoendothelial", "IbarraSoria2018", "IbarraSoria2018_IID/")
outdir <- paste0(wd, "results_core_endothelial.b/PPI_weight/")
inputdir <- paste0(wd, "../data/")
shared_path <- paste0(shared_data_path, "/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

db <- "IbarraSoria2018"

s <- "combined" # specificity method

########## END OF USER INPUT ##########

file <- paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_blood"                  "HiG_cardiac.b"
#  [3] "HiG_cardiac.c"              "HiG_endothelial.a"
#  [5] "HiG_endothelial.c"          "HiG_endothelial.d"
#  [7] "HiG_extraembryonicMesoderm" "HiG_mesodermProgenitors"
#  [9] "HiG_mixedMesoderm.a"        "HiG_mixedMesoderm.b"
# [11] "HiG_pharyngealMesoderm"     "HiG_presomiticMesoderm.a"
# [13] "HiG_presomiticMesoderm.b"   "HiG_somiticMesoderm"
# [15] "HiG_endothelial.b"          "HiG_cardiac.a"
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"
# [19] "CTS_endothelial.b"          "CTS_cardiac.a"
edge_counts <- sapply(graph_list, ecount)
(edge_counts)
#                  HiG_blood              HiG_cardiac.b
#                       3118                       3498
#              HiG_cardiac.c          HiG_endothelial.a
#                       3764                       4382
#          HiG_endothelial.c          HiG_endothelial.d
#                       3966                       3247
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                       2654                       2187
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                       2134                       2750
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                       2717                       2129
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                       3462                       2162
#          HiG_endothelial.b              HiG_cardiac.a
#                       3156                       3072
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a
#                          1                          5
#          CTS_endothelial.b              CTS_cardiac.a
#                          5                         17

(sapply(graph_list, vcount))
#                  HiG_blood              HiG_cardiac.b
#                        374                        391
#              HiG_cardiac.c          HiG_endothelial.a
#                        421                        424
#          HiG_endothelial.c          HiG_endothelial.d
#                        447                        392
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                        345                        330
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                        309                        343
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                        345                        315
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                        379                        305
#          HiG_endothelial.b              HiG_cardiac.a
#                        389                        384
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a
#                         12                         13
#          CTS_endothelial.b              CTS_cardiac.a
#                         27                         37

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

page <- lapply(graph_list, function(x) page_rank(x, directed = FALSE, weights = E(x)$weight)$vector)

df <- lapply(page, function(x) data.frame(PageRank = x, gene = names(x)) %>% arrange(desc(PageRank))) %>%
  rbindlist(., idcol = names(.))
colnames(df)[1] <- "signature"
df$PPI_cat <- lapply(df$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
  unlist() %>%
  factor(., levels = c("CTS", "HiGCTS", "HiG"))
(dim(df)[1]) # 5982

ic <- lapply(graph_list, function(x) eigen_centrality(x, weights = E(x)$weight)$vector)
IC <- lapply(ic, function(x) data.frame(EigenCentrality = x, gene = names(x)) %>% arrange(desc(EigenCentrality))) %>%
  rbindlist(., idcol = names(.))
colnames(IC)[1] <- "signature"
(dim(IC)) # [1] 5982    3
df <- merge(df, IC, by = c("signature", "gene"))
(dim(df)) # [1] 5982    5


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
  cat(i)
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
  # 1: 0.10329182  GATA4
  # 2: 0.09851825   RND3
  # 3: 0.08672068   TBX5
  # 4: 0.07355501 PARD6B
  # 5: 0.06725620 DPYSL3
  # 6: 0.05204086  PTH1R
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
#   signature     gene    PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>         <chr>      <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 CTS_cardiac.a ANXA2    0.0461  CTS                   0      0.999          3.5
# 2 CTS_cardiac.a ARHGAP…  0.0278  CTS                   0      1             21.5
# 3 CTS_cardiac.a ARL4C    0.0300  CTS                   0      1             21.5
# 4 CTS_cardiac.a CGNL1    0.00691 CTS                   0      1             21.5
# 5 CTS_cardiac.a DACT1    0.00691 CTS                   0      1             21.5
# 6 CTS_cardiac.a DBT      0.00691 CTS                   0      1             21.5
(dim(df)) # [1] 5982    8

(subset(df, tolower(gene) == "isl1"))
#   signature       gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>           <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 CTS_cardiac.a   ISL1  0.00691  CTS         0                1             21.5
# 2 HiGCTS_cardiac… ISL1  0.0171   HiGCTS      0                1              7.5
# 3 HiG_cardiac.a   ISL1  0.000766 HiG         0.000143         0.922        194.
# 4 HiG_cardiac.b   ISL1  0.000515 HiG         0.0000159        0.871        198.
# 5 HiG_extraembry… ISL1  0.000511 HiG         0.0000217        0.888        191
# 6 HiG_mixedMesod… ISL1  0.000460 HiG         0.000000646      0.94         320.

# number of significantly high pagerank per PPI_cats, too much control !
n.pr.high <- lapply(names(graph_list), function(x) nrow(subset(df, signature == x & p.PageRank < 0.05))) %>% unlist()
names(n.pr.high) <- names(graph_list)
(n.pr.high)
#                  HiG_blood              HiG_cardiac.b
#                          0                          0
#              HiG_cardiac.c          HiG_endothelial.a
#                          0                          0
#          HiG_endothelial.c          HiG_endothelial.d
#                          0                          0
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                          0                          0
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                          0                          0
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                          0                          0
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                          0                          0
#          HiG_endothelial.b              HiG_cardiac.a
#                          0                          0
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a
#                          0                          0
#          CTS_endothelial.b              CTS_cardiac.a
#                          0                          0


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
  cat(i)
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
  #     annd   gene
  #    <num> <char>
  # 1:     5  MYOCD
  # 2:     5  MEF2C
  # 3:     5 TWIST1

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
(dim(df)) # [1] 5626    9

annd_P[["CTS_cardiac.a"]]
tmp <- lapply(annd_P, function(x) data.frame(p.annd = x, gene = names(x))) %>%
  rbindlist(., idcol = names(.))
colnames(tmp)[1] <- "signature"
(dim(tmp)) # [1] 5982    3

df <- merge(df, tmp, by = c("signature", "gene"))
df[which(is.na(df$knn)), "p.annd"] <- NA ## due to nrow(df) 4878 > nrow(tmp)
(dim(df)) # [1] 5626   10

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
(dim(V_strength)) # 5982    6

## add the V_strength & V_strength_norm info
df <- merge(df, V_strength, by = c("signature", "gene"))
(dim(df)) # 5626   14
(head(df, 3))
#       signature     gene   PageRank PPI_cat EigenCentrality p.PageRank
# 1 CTS_cardiac.a    ANXA2 0.04608295     CTS               0      0.999
# 2 CTS_cardiac.a ARHGAP29 0.02784757     CTS               0      1.000
# 3 CTS_cardiac.a    ARL4C 0.02999776     CTS               0      1.000
#   rank_by_p.PR rank_by_PR     annd p.annd    strength rank_by_strength
# 1          3.5        7.5 1.000000      0 0.002388866               16
# 2         21.5       15.5 4.000000      0 0.010000000               10
# 3         21.5       14.0 3.837341      0 0.010885294                9
#   normalized.strength rank_by_normalized.strength
# 1        6.635739e-05                          16
# 2        2.777778e-04                          10
# 3        3.023693e-04                           9

(table(df$normalized.strength >= df$annd))
# FALSE
#  5626

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
betweenness_list <- lapply(graph_list, function(x) betweenness(x, weights = 1 / E(x)$weight))
bc.median <- lapply(betweenness_list, median) %>% unlist()

for (i in seq_along(betweenness_list)) {
  betweenness_list[[i]] <- data.frame(BetweennessCentrality = betweenness_list[[i]], gene = names(betweenness_list[[i]])) %>%
    mutate(rank_by_BC = rank(-BetweennessCentrality, na.last = "keep")) # Rank the 'annd' values, ignoring NA values
}
df_BC <- betweenness_list %>% rbindlist(., idcol = names(.))
colnames(df_BC)[1] <- "signature"
df_BC$PPI_cat <- lapply(df_BC$signature, function(x) unlist(strsplit(x, split = "_"))[1]) %>% unlist()

(dim(df_BC)) # 5982    5

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
#           signature BetweennessCentrality  gene rank_by_BC PPI_cat
# 81 HiGCTS_cardiac.a                     3 GATA4          1  HiGCTS
#    PCGC_AllCurated
# 81            TRUE

(df5_CHD <- subset(df5, PCGC_AllCurated == TRUE))
(df5_CHD)
#              signature BetweennessCentrality   gene rank_by_BC PPI_cat
# 9        HiG_cardiac.b                  8118  ACTA2          4     HiG
# 13       HiG_cardiac.c                 16346  ACTC1          3     HiG
# 44 HiG_mixedMesoderm.a                  5415   BMP4          4     HiG
# 48 HiG_mixedMesoderm.b                  6181 NKX2-5          4     HiG
# 80       HiG_cardiac.a                 10190  ACTA2          3     HiG
# 81    HiGCTS_cardiac.a                     3  GATA4          1  HiGCTS
# 85       CTS_cardiac.a                     9  GATA4          3     CTS
# 89       CTS_cardiac.a                     4   TBX5          5     CTS
#    PCGC_AllCurated
# 9             TRUE
# 13            TRUE
# 44            TRUE
# 48            TRUE
# 80            TRUE
# 81            TRUE
# 85            TRUE
# 89            TRUE


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
ks.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.01307
ks.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 0.01307
ks.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value =  1
wilcox.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.02935
wilcox.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 0.02935
wilcox.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = NA
t.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 9.595e-06
t.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value = 9.595e-06
t.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = NA

density_median_bc_plot <- ggplot(df_median, aes(x = log10(bc.median), color = PPI_cat, fill = PPI_cat)) +
  geom_density(alpha = 0.3) + # Density lines with transparency
  scale_color_manual(values = PPI_color_palette) +
  scale_fill_manual(values = PPI_color_palette) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(x = "Density of the median of BetweennessCentralitys per PPI", y = "")
# Add statistical comparisons using stat_compare_means manyally

x <- which(df_median$bc.median == 0) # 1 2 3 4

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
(dim(df)) # [1] 5626   16

## reorder df$df$signature to be #E7298A #7570B3 and #E6AB02
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
# # A tibble: 7 × 17
#   signature       gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
#   <chr>           <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>
# 1 HiGCTS_cardiac… ANXA2    0.114 HiGCTS         1   e+ 0          1          5
# 2 HiGCTS_cardiac… ARL4C    0.114 HiGCTS         8.83e-17          1          5
# 3 HiGCTS_cardiac… GATA4    0.219 HiGCTS         9.42e-16          1          5
# 4 HiGCTS_cardiac… MSX2     0.114 HiGCTS         1   e+ 0          1          5
# 5 HiGCTS_cardiac… PODXL    0.114 HiGCTS         5.89e-17          1          5
# 6 HiGCTS_endothe… GATA2    0.286 HiGCTS         1   e+ 0          1          1.5
# 7 HiGCTS_endothe… TAL1     0.286 HiGCTS         1   e+ 0          1          1.5
# # ℹ 10 more variables: rank_by_PR <dbl>, annd <dbl>, p.annd <dbl>,
# #   strength <dbl>, rank_by_strength <dbl>, normalized.strength <dbl>,
# #   rank_by_normalized.strength <int>, rank_by_ANND <dbl>,
# #   rank_by_p.ANND <dbl>, PCGC_AllCurated <lgl>

write.table(df5[, c(
  "signature", "gene", "PageRank", "PPI_cat", "rank_by_PR",
  "normalized.strength", "rank_by_normalized.strength", "PCGC_AllCurated"
)], file = "table_top5_PageRank_perPPI.tsv", sep = "\t", row.names = FALSE, quote = FALSE) # !!!!!!!!!!!!!!

df5_CHD <- subset(df5, PCGC_AllCurated == TRUE)
(dim(df5)) # 97  17
(dim(df5_CHD)) # 4 17
(df5_CHD %>% as.data.frame())
#          signature  gene   PageRank PPI_cat EigenCentrality p.PageRank
# 1    CTS_cardiac.a GATA4 0.10329182     CTS    1.000000e+00      1.000
# 2    CTS_cardiac.a  TBX5 0.08672068     CTS    9.961792e-01      0.999
# 3 HiGCTS_cardiac.a GATA4 0.21930502  HiGCTS    9.420555e-16      1.000
# 4 HiGCTS_cardiac.a  MSX2 0.11428571  HiGCTS    1.000000e+00      1.000
#   rank_by_p.PR rank_by_PR     annd p.annd    strength rank_by_strength
# 1         12.5        1.0 1.651796      0 0.032801167              2.0
# 2          3.5        3.0 4.000000      0 0.030000000              3.5
# 3          5.0        1.0 1.000000      0 0.002801167              1.0
# 4          5.0        3.5 1.000000      0 0.002388866              2.5
#   normalized.strength rank_by_normalized.strength rank_by_ANND rank_by_p.ANND
# 1        0.0009111435                           2           16           10.0
# 2        0.0008333333                           4            4           10.0
# 3        0.0002334306                           1            6            4.5
# 4        0.0001990722                           3            6            4.5
#   PCGC_AllCurated
# 1            TRUE
# 2            TRUE
# 3            TRUE
# 4            TRUE

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
ks.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value =  0.01307
ks.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.01307
ks.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value = 0.3333
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[b, ]$median_PageRank) # p-value = 0.01307
wilcox.test(pg.median[a, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.01307
wilcox.test(pg.median[b, ]$median_PageRank, pg.median[c, ]$median_PageRank) # p-value =  0.3333

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
  df <- readRDS(file = "df_PAGERANK_strength_ANND.rewiring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!

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
  CHD <- readRDS(file = paste0(shared_path, "CHD_Cilia_Genelist.rds"))

  top_genes_CHD <- subset(top_genes, toupper(gene) %in% toupper(unlist(CHD[c("Griffin2023_PCGC_AllCurated")])))
  (dim(top_genes)) # 97  17
  (dim(top_genes_CHD)) # 7 17
  (top_genes_CHD)
  #     # A tibble: 7 × 17
  # # Groups:   signature [3]
  #   signature       gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
  #   <fct>           <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>
  # 1 CTS_cardiac.a   MEF2C 0.0107   CTS            5.70e- 2      1             12.5
  # 2 CTS_cardiac.a   TWIS… 0.00695  CTS            5.48e- 4      1             12.5
  # 3 CTS_cardiac.a   TBX5  0.0867   CTS            9.96e- 1      0.999          3.5
  # 4 HiGCTS_cardiac… MEF2C 0.111    HiGCTS         8.24e-16      0.999          1
  # 5 HiGCTS_cardiac… GATA6 0.109    HiGCTS         7.56e-16      1              5
  # 6 HiGCTS_cardiac… TWIS… 0.0180   HiGCTS         7.36e-18      1              5
  # 7 HiG_mixedMesod… ISL1  0.000460 HiG            6.46e- 7      0.94         302.
  # # ℹ 10 more variables: rank_by_PR <dbl>, annd <dbl>, p.annd <dbl>,
  # #   strength <dbl>, rank_by_strength <dbl>, normalized.strength <dbl>,
  # #   rank_by_normalized.strength <int>, rank_by_ANND <dbl>,
  # #   rank_by_p.ANND <dbl>, label <chr>

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
  ks.test(annd.median[C], annd.median[B]) # p-value = 0.01307
  ks.test(annd.median[C], annd.median[A]) # p-value = 0.01307
  ks.test(annd.median[B], annd.median[A]) # p-value = 0.3333
  wilcox.test(annd.median[C], annd.median[B]) # p-value = 0.02935
  wilcox.test(annd.median[C], annd.median[A]) # p-value = 0.01307
  wilcox.test(annd.median[B], annd.median[A]) # p-value = 0.2207

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
  subset(df, signature == "HiGCTS_cardiac.a")

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
  (dim(top_genes)) # 97  18
  (dim(top_genes_CHD)) # 7 18

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
  ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value = 0.01307
  ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value = 0.01307
  ks.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value = 0.3333
  wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value = 0.01307
  wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value = 0.01307
  wilcox.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value = 0.3333

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
#                  HiG_blood              HiG_cardiac.b
#                        360                        374
#              HiG_cardiac.c          HiG_endothelial.a
#                        405                        402
#          HiG_endothelial.c          HiG_endothelial.d
#                        424                        368
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                        322                        309
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                        291                        325
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                        331                        296
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                        361                        289
#          HiG_endothelial.b              HiG_cardiac.a
#                        373                        360
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a
#                          2                          8
#          CTS_endothelial.b              CTS_cardiac.a
#                          7                         19

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
(dim(bc)) # 5982    5
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
    id2 <- sub("\\..*", "", id) # strip suffix after dot
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

df <- readRDS(file = "df_PAGERANK_strength_ANND.rewiring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!
(colnames(df))
#  [1] "signature"                   "gene"                        "PageRank"                    "PPI_cat"                     "EigenCentrality"             "p.PageRank"
#  [7] "rank_by_p.PR"                "rank_by_PR"                  "annd"                        "p.annd"                      "strength"                    "rank_by_strength"
# [13] "normalized.strength"         "rank_by_normalized.strength" "rank_by_ANND"                "rank_by_p.ANND"

plot_CTS_pr <- plot_HiG_pr <- list()

## compare the CTS&HiG genes to non-HiG genes per CTS-derived PPI
for (id in CTS) {
  if (grepl(".", id, fixed = TRUE)) {
    id2 <- sub("\\..*", "", id) # strip suffix after dot
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
    id2 <- sub("\\..*", "", id) # strip suffix after dot
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
  p_cts_bc <- plot_CTS_bc[[id]]
  p_cts_pr <- plot_CTS_pr[[id]]

  # HiG plots: find matching names containing the CTS id
  hiG_matches <- grep(id2, names(plot_HiG_annd), value = TRUE)
  p_hig_annd <- if (length(hiG_matches) > 0) plot_HiG_annd[[hiG_matches[1]]] else NULL

  hiG_matches_bc <- grep(id2, names(plot_HiG_bc), value = TRUE)
  p_hig_bc <- if (length(hiG_matches_bc) > 0) plot_HiG_bc[[hiG_matches_bc[1]]] else NULL

  hiG_matches_pr <- grep(id2, names(plot_HiG_pr), value = TRUE)
  p_hig_pr <- if (length(hiG_matches_pr) > 0) plot_HiG_pr[[hiG_matches_pr[1]]] else NULL

  # Combine all 6 plots
  plot_list <- list(
    p_cts_annd, p_hig_annd,
    p_cts_bc, p_hig_bc,
    p_cts_pr, p_hig_pr
  )

  # Remove any NULLs in case some HiG plots don't exist
  plot_list <- Filter(Negate(is.null), plot_list)

  # Labels
  labels <- c(
    paste0("CTS_", id, " network"),
    paste0("HiG_", id2, " network"),
    paste0("CTS_", id, " network"),
    paste0("HiG_", id2, " network"),
    paste0("CTS_", id, " network"),
    paste0("HiG_", id2, " network")
  )
  labels <- labels[seq_along(plot_list)]

  # Arrange in 2 columns, 3 rows (adjust automatically if some plots are missing)
  x <- ggarrange(
    plotlist = plot_list,
    labels = labels,
    ncol = 2,
    nrow = ceiling(length(plot_list) / 2)
  )

  print(x)
}

dev.off()
# 1: ggrepel: 1559 unlabeled data points (too many overlaps). Consider increasing max.overlaps
# 2: ggrepel: 19 unlabeled data points (too many overlaps). Consider increasing max.overlaps
# 3: ggrepel: 1641 unlabeled data points (too many overlaps). Consider increasing max.overlaps
# 4: ggrepel: 1562 unlabeled data points (too many overlaps). Consider increasing max.overlaps
# 5: ggrepel: 1124 unlabeled data points (too many overlaps). Consider increasing max.overlaps
# 6: ggrepel: 30 unlabeled data points (too many overlaps). Consider increasing max.overlaps
# 7: ggrepel: 1188 unlabeled data points (too many overlaps). Consider increasing max.overlaps
# 8: ggrepel: 1132 unlabeled data points (too many overlaps). Consider increasing max.overlaps
