library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr)

library(igraph)
library(rstatix)

library(brainGraph)
library(MLmetrics)
library(sm)

########## BEGINNING OF USER INPUT ##########

wd = "/Users/felixyu/Documents/GSE87038_weighted/"

setwd(paste0(wd, "results/PPI_weight/"))

celltype_specific_weight_version <- '10'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_palette <- c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)

db <- "GSE87038"

s <- "combined"  # specificity method

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

###################################################
# See if weights have been updated
# combined_score: the
(names(edge_attr(graph_list[[1]])))
# "weight"         "norm_PPI_score" "corexp_sign"    "coexp_target"  
(all(E(graph_list[[1]])$weight == E(graph_list[[1]])$norm_PPI_score)) # FALSE


###################################################
## 4) evaluate the GRN robustness (stability) of each PPI_cat;
## two trials: node-normalized or not
# https://cwatson.github.io/files/brainGraph_UserGuide.pdf
################################################################
library(brainGraph)
library(data.table)
library(ggplot2)



# Targeted attack ###########################
set.seed(2020)
# Maximal component size as a function of vertices removed.
attack.vertex.strength <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "vertex", "degree"), idcol = names(graph_list))
attack.vertex.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "vertex", "btwn.cent"), idcol = names(graph_list)) # used !!
attack.vertex.strength %>% head()
attack.vertex.btwn %>% head()

attack.vertex <- rbind(attack.vertex.strength, attack.vertex.btwn)
colnames(attack.vertex)[1] <- "signature"
attack.vertex$PPI_cat <- lapply(attack.vertex$signature, function(x) unlist(strsplit(x, split = "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))
head(attack.vertex, 3)
attack.vertex$cluster <- lapply(attack.vertex$signature, function(x) unlist(strsplit(x, split = "_"))[2]) %>% unlist()

## plot attack measured by btwn.cent
g_attack <- ggplot(
    data = subset(attack.vertex, measure == "btwn.cent"),
    aes(
        x = removed.pct, y = comp.pct, color = PPI_cat, size = PPI_cat, # interaction(PPI_cat, signature),
        linetype = PPI_cat # comment this line if highlighting one subtype when color=signature
        # shape=PPI_cat  #
    )
) +
    # geom_line(size=ifelse(subset(attack.vertex ,measure == 'strength')$signature=='endothelial.a', 2, 1)) + # uisng this line if highlighting one subtype
    geom_line(aes(group = signature)) + # group ensures each line is drawn independently
    scale_color_manual(values = PPI_color_palette) +
    scale_size_manual(values = PPI_size_palette) +
    geom_abline(slope = -1, intercept = 1, col = "gray", lty = 2) +
    theme(legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1)) +
    ylab("the remaining maximal component size / the initial maximal component size") +
    ggtitle("vertex robustness by betweenness centrality")
print(g_attack)
gsignature_attack <- ggplot(
    data = subset(attack.vertex, measure == "btwn.cent"), # !!!!!!!!!!!!!
    aes(
        x = removed.pct, y = comp.pct, color = cluster, size = PPI_cat, # interaction(PPI_cat, signature),
        linetype = PPI_cat # comment this line if highlighting one subtype when color=signature
        # shape=PPI_cat  #
    )
) +
    # geom_line(size=ifelse(subset(attack.vertex ,measure == 'strength')$signature=='endothelial.a', 2, 1)) + # use this line if highlighting one subtype
    geom_line() +
    geom_abline(slope = -1, intercept = 1, col = "gray", lty = 2) +
    theme(legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1)) +
    scale_size_manual(values = PPI_size_palette) +
    ylab("the remaining maximal component size / the initial maximal component size") +
    ggtitle("vertex robustness by betweenness centrality")
print(gsignature_attack)

## plot attack measured by strength
g_attack2 <- ggplot(
    data = subset(attack.vertex, measure == "degree"),
    aes(
        x = removed.pct, y = comp.pct, col = PPI_cat, size = PPI_cat,
        linetype = PPI_cat # comment this line if highlighting one subtype when color=signature
        # shape=PPI_cat
    )
) +
    # geom_line(size=ifelse(subset(attack.vertex ,measure == 'strength')$signature=='endothelial.a', 2, 1)) + # use this line if highlighting one subtype
    geom_line(aes(group = signature)) + # comment this line if highlighting one subtype
    scale_color_manual(values = PPI_color_palette) +
    scale_size_manual(values = PPI_size_palette) +
    geom_abline(slope = -1, intercept = 1, col = "gray", lty = 2) +
    theme(legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1)) +
    ggtitle("vertex robustness measured by strength")
print(g_attack2)

gsignature_attack2 <- ggplot(
    data = subset(attack.vertex, measure == "degree"),
    aes(
        x = removed.pct, y = comp.pct, col = cluster, size = PPI_cat,
        linetype = PPI_cat
    )
) +
    geom_line(size = 1.2) + # <- thicker lines
    scale_size_manual(values = PPI_size_palette) +
    geom_abline(slope = -1, intercept = 1, col = "gray", lty = 2) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1),
        text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1.2, "cm")
    ) +
    ggtitle("Vertex robustness measured by strength")

   pdf(file=paste0('vertex_attack_', db,'.pdf')) 
   print(gridExtra::grid.arrange(g_attack ,  g_attack2 , 
					gsignature_attack ,  gsignature_attack2 , 
                      ncol = 2, nrow = 2))
   dev.off()


############## to restart here anytime #################################
library(brainGraph)
require(dplyr)
library(data.table)
library(ggplot2)

# refer to 11.1.1_CTS_cardiac_network_strengthDistribution.R
file <- paste0(db, "_STRING_graph_perState_simplified_combinedweighted.rds")
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_8"
# [26] "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"      "CTS_16.1"
# [31] "CTS_13"      "CTS_8"       "HiGCTS_13"


attack.vertex.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "vertex", "btwn.cent"), idcol = names(graph_list))
colnames(attack.vertex.btwn)[1] <- "signature"
(head(attack.vertex.btwn, 3))
#    signature                   type   measure comp.size  comp.pct removed.pct
#       <char>                 <char>    <char>     <num>     <num>       <num>
# 1:     HiG_1 Targeted vertex attack btwn.cent       300 1.0000000  0.00000000
# 2:     HiG_1 Targeted vertex attack btwn.cent       299 0.9966667  0.00330033
# 3:     HiG_1 Targeted vertex attack btwn.cent       298 0.9933333  0.00660066

(dim(attack.vertex.btwn)) #  8245    6
saveRDS(attack.vertex.btwn, file = "attack.vertex.btwn.rds")

# calculate the 'egde' btwn for each G first
# igraph::edge_betweenness() uses distance graph weights, but E(g) uses connection weights, thus we invert it.
for (j in seq_along(graph_list)) E(graph_list[[j]])$btwn <- edge_betweenness(graph_list[[j]], weights = 1/E(graph_list[[j]])$weight)
## then do the edge-attack analysis
attack.edge.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "edge", "btwn.cent"), idcol = names(graph_list)) # !!!!!!!!!!!!!
colnames(attack.edge.btwn)[1] <- "signature"
(dim(attack.edge.btwn)) #  151870      6
(head(attack.edge.btwn, 3))
#    signature                 type   measure comp.size comp.pct  removed.pct
#       <char>               <char>    <char>     <num>    <num>        <num>
# 1:     HiG_1 Targeted edge attack btwn.cent       300        1 0.0000000000
# 2:     HiG_1 Targeted edge attack btwn.cent       300        1 0.0002008839
# 3:     HiG_1 Targeted edge attack btwn.cent       300        1 0.0004017678

(table(attack.edge.btwn$signature))
#      CTS_11      CTS_13      CTS_15      CTS_16    CTS_16.1       CTS_7 
#          55         148         208          71         211          53 
#       CTS_8       HiG_1      HiG_10      HiG_11      HiG_12      HiG_13 
#         166        4979        8656        9302        5829        8509 
#      HiG_14      HiG_15      HiG_16      HiG_17      HiG_18      HiG_19 
#        7005        8633       10582       11720        8547        7829 
#       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6       HiG_7 
#        9545        7121        5363        5809       10994        7573 
#       HiG_8       HiG_9   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1 
#        5881        6938          16          55          11          32 
#    HiGCTS_7    HiGCTS_8 
#          18          11 

saveRDS(attack.edge.btwn, file = "attack.edge.btwn.rds")

# In a random failure analysis, you choose a vertex/edge at random, remove it, and calculate the maximum
#  component size until all elements have been removed.
####################################################
# refer to 11.2_CTS_cardiac_network_robustness_Miudway3.simulation.R
## run on Midway3, the following is long. DO NOT REPEAT !!!)


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

attack.edge.btwn <- readRDS(file = "attack.edge.btwn.rds")
attack.vertex.btwn <- readRDS(file = "attack.vertex.btwn.rds")
failure.vertex <- readRDS(paste0("failure.vertex_100_simplified_", s, "weighted.rds"))
(table(attack.vertex.btwn$signature, attack.vertex.btwn$type))
#               Targeted vertex attack
#   CTS_11                          52
#   CTS_13                          61
#   CTS_15                          67
#   CTS_16                          40
#   CTS_16.1                        80
#   CTS_7                           32
#   CTS_8                           55
#   HiG_1                          304
#   HiG_10                         423
#   HiG_11                         442
#   HiG_12                         342
#   HiG_13                         404
#   HiG_14                         365
#   HiG_15                         407
#   HiG_16                         525
#   HiG_17                         524
#   HiG_18                         456
#   HiG_19                         528
#   HiG_2                          436
#   HiG_3                          412
#   HiG_4                          305
#   HiG_5                          321
#   HiG_6                          513
#   HiG_7                          337
#   HiG_8                          333
#   HiG_9                          359
#   HiGCTS_11                       20
#   HiGCTS_15                       31
#   HiGCTS_16                       14
#   HiGCTS_16.1                     31
#   HiGCTS_7                        15
#   HiGCTS_8                        11

(subset(attack.vertex.btwn, signature == "HiGCTS_8"))
#     signature                   type   measure comp.size  comp.pct removed.pct
#        <char>                 <char>    <char>     <num>     <num>       <num>
#  1:  HiGCTS_8 Targeted vertex attack btwn.cent         7 1.0000000         0.0
#  2:  HiGCTS_8 Targeted vertex attack btwn.cent         4 0.5714286         0.1
#  3:  HiGCTS_8 Targeted vertex attack btwn.cent         3 0.4285714         0.2
#  4:  HiGCTS_8 Targeted vertex attack btwn.cent         3 0.4285714         0.3
#  5:  HiGCTS_8 Targeted vertex attack btwn.cent         2 0.2857143         0.4
#  6:  HiGCTS_8 Targeted vertex attack btwn.cent         2 0.2857143         0.5
#  7:  HiGCTS_8 Targeted vertex attack btwn.cent         1 0.1428571         0.6
#  8:  HiGCTS_8 Targeted vertex attack btwn.cent         1 0.1428571         0.7
#  9:  HiGCTS_8 Targeted vertex attack btwn.cent         1 0.1428571         0.8
# 10:  HiGCTS_8 Targeted vertex attack btwn.cent         1 0.1428571         0.9
# 11:  HiGCTS_8 Targeted vertex attack btwn.cent         0 0.0000000         1.0


failure.edge <- readRDS(paste0("failure.edge_100_simplified_", s, "weighted.rds"))
failure.dt <- rbind(failure.edge, failure.vertex)
(head(failure.dt, 3))
#     HiG_1                type measure comp.size comp.pct  removed.pct
#    <char>              <char>  <char>     <num>    <num>        <num>
# 1:  HiG_1 Random edge removal  random       300        1 0.0000000000
# 2:  HiG_1 Random edge removal  random       300        1 0.0002008839
# 3:  HiG_1 Random edge removal  random       300        1 0.0004017678

colnames(failure.dt)[1] <- "signature"
(table(failure.dt$signature, failure.dt$type))
#               Random edge removal Random vertex removal
#   CTS_11                       55                    52
#   CTS_13                      148                    61
#   CTS_15                      208                    67
#   CTS_16                       71                    40
#   CTS_16.1                    211                    80
#   CTS_7                        53                    32
#   CTS_8                       166                    55
#   HiG_1                      4979                   304
#   HiG_10                     8656                   423
#   HiG_11                     9302                   442
#   HiG_12                     5829                   342
#   HiG_13                     8509                   404
#   HiG_14                     7005                   365
#   HiG_15                     8633                   407
#   HiG_16                    10582                   525
#   HiG_17                    11720                   524
#   HiG_18                     8547                   456
#   HiG_19                     7829                   528
#   HiG_2                      9545                   436
#   HiG_3                      7121                   412
#   HiG_4                      5363                   305
#   HiG_5                      5809                   321
#   HiG_6                     10994                   513
#   HiG_7                      7573                   337
#   HiG_8                      5881                   333
#   HiG_9                      6938                   359
#   HiGCTS_11                    16                    20
#   HiGCTS_15                    55                    31
#   HiGCTS_16                    11                    14
#   HiGCTS_16.1                  32                    31
#   HiGCTS_7                     18                    15
#   HiGCTS_8                     11                    11


colnames(attack.vertex.btwn)[1] <- "signature"
(dim(failure.dt)) #  160115      6

robustness.dt <- rbind(failure.dt, attack.vertex.btwn, attack.edge.btwn)
(dim(robustness.dt)) #  320230      6
robustness.dt$PPI_cat <- lapply(robustness.dt$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))
(head(robustness.dt, 3))
#    signature                type measure comp.size comp.pct  removed.pct
#       <char>              <char>  <char>     <num>    <num>        <num>
# 1:     HiG_1 Random edge removal  random       300        1 0.0000000000
# 2:     HiG_1 Random edge removal  random       300        1 0.0002008839
# 3:     HiG_1 Random edge removal  random       300        1 0.0004017678

robustness.dt$experiment <- ifelse(grepl("edge", robustness.dt$type), "edge", "vertex")
robustness.dt$measure <- factor(robustness.dt$measure, levels = c("random", "btwn.cent"))

(table(robustness.dt$type, robustness.dt$measure))
#                          random btwn.cent
#   Random edge removal    151870         0
#   Random vertex removal    8245         0
#   Targeted edge attack        0    151870
#   Targeted vertex attack      0      8245


robustness.dt$type <- factor(robustness.dt$type,
    levels = c("Random edge removal", "Targeted edge attack", "Random vertex removal", "Targeted vertex attack")
)

robustness.dt$cluster <- lapply(robustness.dt$signature, function(x) unlist(strsplit(x, "_"))[2]) %>% unlist()

p_attack3 <- ggplot(
    robustness.dt,
    aes(x = removed.pct, y = comp.pct, col = cluster, linetype = PPI_cat, size = PPI_cat)
) +
    geom_line(show.legend = TRUE) +
    scale_size_manual(values = PPI_size_palette, name = "PPI Category") +
    facet_wrap(~type) +
    geom_abline(slope = -1, intercept = 1, col = "gray", lty = 2) +
    labs(
        x = "% edges/vertex removed",
        y = "% of max. component remaining",
        title = "Robustness by Cluster (Random vs. Targeted Attacks)"
    ) +
    guides(
        color = guide_legend(title = "Cluster"),
        linetype = guide_legend(title = "PPI Category"),
        size = guide_legend(title = "PPI Category")
    ) +
    theme(
        legend.position = "right", # vertical legend on the right
        legend.box = "vertical",
        legend.box.just = "center",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        panel.spacing.x = unit(0.5, "cm"),
        text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)
    )



print(p_attack3)
p_attack4 <- ggplot(
    robustness.dt,
    aes(x = removed.pct, y = comp.pct, col = PPI_cat, linetype = PPI_cat, size = PPI_cat)
) +
    geom_line(aes(group = signature), show.legend = TRUE) + # Ensure legend is shown
    scale_color_manual(values = PPI_color_palette, name = "PPI Category") +
    scale_size_manual(values = PPI_size_palette, name = "PPI Category") +
    facet_wrap(~type) +
    geom_abline(slope = -1, intercept = 1, col = "gray", lty = 2) +
    theme(
        legend.position = "right", # moves legend outside the plot
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.6, "cm"), # shrink legend box
        legend.spacing.x = unit(0.3, "cm"),
        text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.spacing.x = unit(0.5, "cm")
    ) +
    labs(
        x = "% edges/vertex removed", y = "% of max. component remaining",
        title = "Robustness under Random and Targeted Attacks"
    )

print(p_attack4)
pdf(file = paste0("attack_", db, ".pdf"), width = 8, height = 7)
plot(p_attack4)
plot(p_attack3)
dev.off() # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


### compare targeted attack vs random removal across categories  ################
#########################################################

# First, Run Shapiro-Wilk test, ensuring that we have enough non-NA observations per group
normality_tests <- robustness.dt %>%
    group_by(experiment, PPI_cat, measure) %>%
    filter(sum(!is.na(comp.pct)) >= 3) %>% # Ensure there are at least 3 non-NA values
    summarise(
        shapiro_test = list(
            tryCatch(
                expr = shapiro.test(comp.pct),
                error = function(e) NA # If the test fails, return NA instead
            )
        ),
        .groups = "drop"
    ) %>%
    mutate(
        p_value = sapply(shapiro_test, function(x) if (is.list(x)) x$p.value else NA) # Safely extract p-value
    ) %>%
    as.data.frame()

# Print the normality test results
(normality_tests)
#    experiment PPI_cat   measure shapiro_test      p_value
# 1        edge     CTS    random c(W = 0..... 6.162992e-22
# 2        edge     CTS btwn.cent c(W = 0..... 2.755804e-18
# 3        edge  HiGCTS    random c(W = 0..... 2.546931e-07
# 4        edge  HiGCTS btwn.cent c(W = 0..... 9.895911e-07
# 5        edge     HiG    random           NA           NA
# 6        edge     HiG btwn.cent           NA           NA
# 7      vertex     CTS    random c(W = 0..... 6.105198e-13
# 8      vertex     CTS btwn.cent c(W = 0..... 2.510821e-22
# 9      vertex  HiGCTS    random c(W = 0..... 7.115017e-06
# 10     vertex  HiGCTS btwn.cent c(W = 0..... 5.624286e-13
# 11     vertex     HiG    random           NA           NA
# 12     vertex     HiG btwn.cent           NA           NA

## then, Plot the Wilcox results for visualization and manually add fold changes !!!
# Wrap long x-axis labels to prevent overlap
# Combine experiment (edge vs vertex) into the x-axis labels
robustness.dt$type2 <- paste0(robustness.dt$experiment, "\n", robustness.dt$type)

# Wrap long labels for readability
robustness.dt$type2 <- factor(robustness.dt$type2,
    levels = c(
        "edge\nRandom edge removal", "edge\nTargeted edge attack",
        "vertex\nRandom vertex removal", "vertex\nTargeted vertex attack"
    ),
    labels = c(
        "Edge\nRandom", "Edge\nTargeted",
        "Vertex\nRandom", "Vertex\nTargeted"
    )
)

# Get y-position for text (adjust if needed)
y_max <- max(robustness.dt$comp.pct, na.rm = TRUE)

# Calculate fold change for each PPI_cat
fold_change <- robustness.dt %>%
    group_by(PPI_cat) %>%
    summarise(
        fold_change_edge = mean(comp.pct[type == "Random edge removal"], na.rm = TRUE) /
            mean(comp.pct[type == "Targeted edge attack"], na.rm = TRUE),
        fold_change_vertex = mean(comp.pct[type == "Random vertex removal"], na.rm = TRUE) /
            mean(comp.pct[type == "Targeted vertex attack"], na.rm = TRUE)
    )
(fold_change)
#   PPI_cat fold_change_edge fold_change_vertex
#   <fct>              <dbl>              <dbl>
# 1 CTS                1.04                1.66
# 2 HiGCTS             1.00                1.68
# 3 HiG                0.913               1.09

# Add annotation label and y-position to fold_change table
fold_change <- fold_change %>%
    mutate(
        label = paste0("Edge FC: ", round(fold_change_edge, 2), "\nVertex FC: ", round(fold_change_vertex, 2)),
        x = 2.5, # midpoint between the 4 attack types on x-axis
        y = y_max * 0.95 # place near the top, adjust multiplier if too close
    )

# Plot: Single combined boxplot with all types
g <- ggplot(robustness.dt, aes(x = type2, y = comp.pct, fill = measure, color = PPI_cat, size = PPI_cat)) +
    geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) +
    facet_wrap(~PPI_cat, ncol = 3) +
    scale_color_manual(values = PPI_color_palette) +
    scale_size_manual(values = PPI_size_palette) +
    scale_fill_manual(values = c("random" = "grey", "btwn.cent" = "white")) +
    theme_minimal(base_size = 14) +
    labs(
        title = "Combined Robustness Measures: Edge and Vertex Attacks",
        x = "Attack Type",
        y = "Component Percentage Remaining"
    ) +
    theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        strip.text = element_text(size = 14)
    ) +
    geom_label(
        data = fold_change,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = 4,
        hjust = 0.5,
        fill = "white",
        alpha = 0.7,
        label.size = 0
    )

# Filter the data for 'Targeted vertex attack'
targeted_vertex_data <- robustness.dt[robustness.dt$type == "Targeted vertex attack", ]
targeted_vertex_data <- targeted_vertex_data[!is.na(targeted_vertex_data$comp.pct), ]

g2 <- ggplot(targeted_vertex_data, aes(x = PPI_cat, y = comp.pct, fill = PPI_cat)) +
    geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) + # Dodge the boxes for each type per PPI_cat
    scale_fill_manual(values = PPI_color_palette) +
    geom_signif(
        comparisons = list(
            c("CTS", "HiGCTS"), c("HiGCTS", "HiG"), c("CTS", "HiG")
        ),
        map_signif_level = TRUE,
        step_increase = 0.1, # Adjusts spacing between the lines
        test = "wilcox.test", # Perform a t-test to calculate significance
        p.adjust.method = "holm" # Adjust p-values using the Holm method, default is "bonferroni"
    ) +
    theme_minimal() +
    labs(
        title = "Comparison of Robustness Measures by cent.btw",
        x = "PPI Category",
        y = "Component Percentage Remaining"
    ) +
    theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5) # center title
    )

print(g2)

pdf(file = paste0("box_wilcox-test_attack_", db, ".pdf"), width = 12)
print(g)
print(g2)
dev.off()

# confirm that the figure showing adjusted p-values
wilcox.test(
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "CTS")],
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiGCTS")]
)

# W = 21728, p-value = 0.1847
wilcox.test(
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "CTS")],
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiG")]
)

# W = 1031924, p-value < 2.2e-16
wilcox.test(
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiGCTS")],
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiG")]
)

# W = 320038, p-value = 1.008e-09

## finally, manually add the threshold of fold changes  ###############
robustness.dt <- robustness.dt %>%
    group_by(PPI_cat, type) %>%
    mutate(
        mean_comp_pct = mean(comp.pct, na.rm = TRUE) # Calculate mean comp.pct for each group
    ) %>%
    ungroup()

## for each signature, compare the comp.pct between targeted attach to its random removal to access the significance of 'hub'
#######################################
Hub_effect <- array(dim = length(graph_list))
names(Hub_effect) <- names(graph_list)
for (j in names(graph_list)) {
    tmp <- subset(robustness.dt, signature == j)
    dim(tmp) # 10566    11
    table(tmp$type)
    #    Random edge removal   Targeted edge attack  Random vertex removal Targeted vertex attack
    #                   4979                   4979                   304				304
    x <- subset(tmp, type == "Random edge removal")$comp.pct
    y <- subset(tmp, type == "Targeted edge attack")$comp.pct
    Hub_effect[j] <- wilcox.test(x, y)$p.value
}
df <- data.frame(
    edge_p = Hub_effect,
    edge_p.adj = p.adjust(Hub_effect, method = "bonferroni")
)

Hub_effect <- array(dim = length(graph_list))
names(Hub_effect) <- names(graph_list)
for (j in names(graph_list)) {
    tmp <- subset(robustness.dt, signature == j)
    x <- subset(tmp, type == "Random vertex removal")$comp.pct
    y <- subset(tmp, type == "Targeted vertex attack")$comp.pct
    Hub_effect[j] <- wilcox.test(x, y)$p.value
}
df$node_p <- Hub_effect
df$node_p.adj <- p.adjust(df$node_p, method = "bonferroni")

(df[which(df$edge_p.adj < 0.05), ])
#               edge_p    edge_p.adj      node_p node_p.adj
# HiG_1   3.637413e-51  1.163972e-49 0.087581455  1.0000000
# HiG_2   6.798073e-87  2.175383e-85 0.073828047  1.0000000
# HiG_3   1.487083e-65  4.758664e-64 0.003209893  0.1027166
# HiG_4   7.306022e-38  2.337927e-36 0.080204149  1.0000000
# HiG_5   5.755739e-56  1.841836e-54 0.214137044  1.0000000
# HiG_6  3.266294e-112 1.045214e-110 0.010388843  0.3324430
# HiG_9   3.612884e-77  1.156123e-75 0.097339992  1.0000000
# HiG_10 2.323184e-119 7.434188e-118 0.204113846  1.0000000
# HiG_12  1.453618e-62  4.651576e-61 0.073777875  1.0000000
# HiG_14  1.328686e-49  4.251795e-48 0.008231469  0.2634070
# HiG_17 5.077442e-136 1.624782e-134 0.034433999  1.0000000
# HiG_18 1.659211e-114 5.309475e-113 0.014798231  0.4735434
# HiG_19 1.377912e-141 4.409319e-140 0.008057704  0.2578465
# HiG_7   1.205950e-82  3.859038e-81 0.060765374  1.0000000
# HiG_11 2.869255e-100  9.181615e-99 0.044193528  1.0000000
# HiG_15 3.511180e-118 1.123578e-116 0.007319833  0.2342347
# HiG_16 1.026564e-138 3.285003e-137 0.014489861  0.4636756
# HiG_13  9.053573e-88  2.897143e-86 0.030944561  0.9902260
# HiG_8   2.234313e-53  7.149801e-52 0.215643811  1.0000000


(df[which(df$node_p.adj < 0.05), ])
#                  edge_p edge_p.adj       node_p   node_p.adj
# HiGCTS_16.1 0.892811555  1.0000000 5.747580e-04 0.0183922549
# CTS_11      0.222390265  1.0000000 1.198383e-04 0.0038348252
# CTS_16.1    0.013243249  0.4237840 4.469871e-06 0.0001430359
# CTS_13      0.007453426  0.2385096 7.639771e-06 0.0002444727


df$clust <- lapply(rownames(df), function(x) unlist(strsplit(x, split = "_"))[2]) %>% unlist()
df$PPI_cat <- lapply(rownames(df), function(x) unlist(strsplit(x, split = "_"))[1]) %>% unlist()
df$PPI_cat <- factor(df$PPI_cat, levels = c("HiG", "HiGCTS", "CTS"))

# Create a new column to differentiate edge and node
tmp <- df[, c(1, 2, 5, 6)] %>%
    mutate(type = "edge")
tmp2 <- df[, c(3, 4, 5, 6)] %>%
    mutate(type = "vertex")
colnames(tmp)[1:2] <- colnames(tmp2)[1:2] <- c("p", "p.adj")
tmp <- rbind(tmp, tmp2)
tmp$PPI_cat <- factor(tmp$PPI_cat, levels = c("HiG", "HiGCTS", "CTS"))

# Plot the results with subpanels
g <- ggplot(tmp, aes(x = as.factor(clust), y = -log10(p.adj), color = PPI_cat)) +
    geom_point(aes(shape = PPI_cat), size = 3) + # Points for each PPI_cat
    geom_line(aes(group = PPI_cat, linetype = PPI_cat)) + # Lines for each PPI_cat
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Horizontal line at adj.p = 0.05
    labs(
        x = "Cluster",
        y = "-log10(Adjusted p-values)",
        color = "State",
        shape = "State"
    ) +
    scale_color_manual(values = PPI_color_palette) +
    scale_shape_manual(values = c("HiG" = 16, "CTS" = 17, "HiGCTS" = 18)) +
    facet_wrap(~type, scales = "free", nrow = 2) + # Split into two rows for edge and node
    theme_minimal() +
    theme(legend.position = "top") +
    ggtitle("Adjusted p-values of Wilcox test comparing targeted attack vs random removal")
print(g)
## discuss the irrelevent of the p-values and node strength levels
tmp$count <- c(
    sapply(graph_list, ecount),
    sapply(graph_list, vcount)
)
g2 <- ggplot(tmp, aes(x = log10(count), y = -log10(p.adj), color = PPI_cat)) +
    geom_point(aes(shape = PPI_cat), size = 3) + # Points for each PPI_cat
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Horizontal line at adj.p = 0.05
    labs(
        x = "log10(PPI size)",
        y = "-log10(Adjusted p-values)",
        color = "State",
        shape = "State"
    ) +
    scale_color_manual(values = PPI_color_palette) +
    scale_shape_manual(values = c("HiG" = 16, "CTS" = 17, "HiGCTS" = 18)) +
    facet_wrap(~type, scales = "free", ncol = 2) + # Split into two rows for edge and node
    theme_minimal() +
    theme(legend.position = "top") +
    ggtitle("Adjusted p-values of Wilcox test comparing targeted attack vs random removal")
print(g2)
pdf(file = paste0("line.adjp_wilcox_attack_", db, ".pdf"))
print(g)
print(g2)
dev.off()


#############################################################
# 4.2) estimate the significance of robustness using AUC (????can't repeat mopuse observation!!)
# the observed lin always sit left to the simulated distributions,
# suggesting gene in the observed PPI is mroe important than randomly selected genes,
# becasue targeted attack at them significantly reduced the component size of the network
#
# for the
#############################################################


(names(graph_list))
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_8"
# [26] "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"      "CTS_16.1"
# [31] "CTS_13"      "CTS_8"       "HiGCTS_13"

observed_auc_list <- list()
for (j in names(graph_list)) {
    observed_auc_list[[j]] <- Area_Under_Curve(
        subset(robustness.dt, signature == j & type == "Targeted vertex attack")$removed.pct,
        subset(robustness.dt, signature == j & type == "Targeted vertex attack")$comp.pct
    )
}
(observed_auc_list %>% unlist())
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6 
#   0.4488999   0.4559308   0.4282119   0.4488711   0.4633703   0.4450483 
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18 
#   0.4548706   0.4681002   0.4488626   0.4328458   0.4536677   0.4446846 
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13 
#   0.4432373   0.4499778   0.4533461   0.4376471   0.4493834   0.4465732 
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1 
#   0.4632466   0.2817460   0.1973684   0.2388889   0.3076923   0.1833333 
#    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16    CTS_16.1 
#   0.3071429   0.2377279   0.1840959   0.3123737   0.2508013   0.2466465 
#      CTS_13       CTS_8 
#   0.2101852   0.2943673 

df_AUC <- data.frame(
    auc = observed_auc_list %>% unlist(),
    signature = names(observed_auc_list),
    PPI_cat = lapply(names(observed_auc_list), function(x) unlist(strsplit(x, split = "_"))[1]) %>% unlist()
)
df_AUC$PPI_cat <- factor(df_AUC$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

g3 <- ggplot(df_AUC, aes(x = PPI_cat, y = auc, fill = PPI_cat)) +
    geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = PPI_color_palette) +
    geom_signif(
        comparisons = list(
            c("CTS", "HiGCTS"), c("HiGCTS", "HiG"), c("CTS", "HiG")
        ),
        map_signif_level = TRUE,
        step_increase = 0.1,
        test = "wilcox.test",
        p.adjust.method = "holm"
    ) +
    theme_minimal() +
    labs(
        title = "Comparison of Robustness Measures by cent.btw",
        x = "PPI Category",
        y = "AUC of targeted vertex attack"
    ) +
    theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5)
    )
ggsave(
    filename = paste0("box_wilcox-test_attack_AUC_", db, ".pdf"),
    plot = g3,
    width = 10,
    height = 7
)


vn <- sapply(graph_list, vcount)
pn <- array(dim = length(graph_list))
names(pn) <- names(graph_list)
for (j in names(vn)) {
    pn[j] <- igraph::strength(graph_list[[j]], weights = E(graph_list[[j]])$weight) %>% mean() / vn[j]
}
(vn)
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


(sum(vn[grep("^HiG_", names(vn))]))     # 7717
(sum(vn[grep("^HiGCTS_", names(vn))]))  # 116
(sum(vn[grep("^CTS_", names(vn))]))     # 380

#### significance evidence 1, comp.pct (the ratio of maximal component size after each random removal to the observed graph's maximal component size)
### was calculated after	random vertex_attaching
## then calculate the AUC of CTS_muscle at muscle cluster vs HiG at nearby stady clusters (C1 & C5)
## this simulation runs longer than expected, therefore go to midway3 !!!
# N = 1000  for vertex removal N=100 for edge removal

# refer to 11.2_CTS_cardiac_network_robustness_Miudway3.simulation.R
## run on Midway3, the following is long. DO NOT REPEAT !!!)

attac_V_random <- readRDS(paste0("AUC_compt.pct_attac_V_random_1000runs_", s, "weighted.RDS"))
N <- nrow(attac_V_random)

col_exists <- function(df, col) col %in% colnames(df)

# Get all cluster IDs based on CTS_ columns present
cluster_ids <- unique(sub("^CTS_", "", grep("^CTS_", colnames(attac_V_random), value = TRUE)))

# Precompute base cluster IDs without suffixes like ".1"
cluster_base <- sapply(cluster_ids, function(x) {
    if (grepl(".", x, fixed = TRUE)) strsplit(x, ".", fixed = TRUE)[[1]][1] else x
})
names(cluster_base) <- cluster_ids

pdf("Vertex_All_AUC_density_plots.pdf", width = 20, height = 4)
par(
    mfrow = c(1, 7),
    mar = c(4.5, 4.5, 4.2, 2.5), # bottom, left, top, right (increased)
    oma = c(2, 1, 2, 1) # outer margin around full page
)


for (j in cluster_ids) {
    j0 <- cluster_base[j]

    # Construct column names for each category
    cts_col <- paste0("CTS_", j)
    higcts_col <- paste0("HiGCTS_", j)
    hig_col <- paste0("HiG_", j0)

    # Skip if none of these columns exist in the data
    if (!any(c(cts_col, higcts_col, hig_col) %in% colnames(attac_V_random))) {
        message("Skipping cluster ", j, ": no relevant columns found.")
        next
    }

    # Gather all existing columns for plotting
    obs_cols <- c()
    lf_groups <- c()
    x_lines <- numeric(0)

    if (col_exists(attac_V_random, cts_col)) {
        obs_cols <- c(obs_cols, cts_col)
        lf_groups <- c(lf_groups, rep("CTS", N))
        x_lines <- c(x_lines, attac_V_random[, cts_col])
    }
    if (col_exists(attac_V_random, higcts_col)) {
        obs_cols <- c(obs_cols, higcts_col)
        lf_groups <- c(lf_groups, rep("HiGCTS", N))
        x_lines <- c(x_lines, attac_V_random[, higcts_col])
    }
    if (col_exists(attac_V_random, hig_col)) {
        obs_cols <- c(obs_cols, hig_col)
        lf_groups <- c(lf_groups, rep("HiG", N))
        x_lines <- c(x_lines, attac_V_random[, hig_col])
    }

    lf <- factor(lf_groups, levels = c("CTS", "HiGCTS", "HiG"))

    # Calculate xlim range covering data and observed values
    observed_values <- sapply(obs_cols, function(col) {
        if (!is.null(observed_auc_list[[col]])) observed_auc_list[[col]] else NA
    })
    xlim_range <- range(c(x_lines, observed_values), na.rm = TRUE)

    # Plot density comparison
    sm.density.compare(
        x_lines,
        lf,
        xlab = "AUC of maximal component size",
        col = PPI_color_palette,
        xlim = xlim_range,
        lwd = 2,
        main = paste0("Cluster ", j)
    )

    # Add vertical lines for observed values
    for (col in obs_cols) {
        if (!is.null(observed_auc_list[[col]])) {
            label <- sub("^(CTS|HiGCTS|HiG)_.*", "\\1", col)
            abline(v = observed_auc_list[[col]], col = PPI_color_palette[label], lty = 2, lwd = 2)
        }
    }

    # KS test between CTS and HiG if both exist, add p-value as mtext
    if (col_exists(attac_V_random, cts_col) && col_exists(attac_V_random, hig_col)) {
        ks_p <- ks.test(attac_V_random[, cts_col], attac_V_random[, hig_col])$p.value
        mtext(sprintf("KS p = %.2g", ks_p), side = 1, line = 5)
    }

    # Add empirical p-values above plot for each group
    for (label in c("CTS", "HiGCTS", "HiG")) {
        colname <- paste0(label, "_", ifelse(label == "HiG", j0, j))
        if (col_exists(attac_V_random, colname) && !is.null(observed_auc_list[[colname]])) {
            above <- sum(attac_V_random[, colname] > observed_auc_list[[colname]])
            below <- sum(attac_V_random[, colname] < observed_auc_list[[colname]])
            txt_cluster <- ifelse(label == "HiG", j0, j)
            mtext(sprintf("%s_%s p= %.2f ; %.2f", label, txt_cluster, above / N, below / N),
                side = 3, line = switch(label,
                    "CTS" = 1,
                    "HiGCTS" = 2.5,
                    "HiG" = 4
                ),
                col = PPI_color_palette[label]
            )
        }
    }
}

dev.off()



#############################################################
# 4.2) estimate the significance of robustness   using AUC , edge removal
#############################################################

observed_auc_list <- list()
for (j in names(graph_list)) {
    observed_auc_list[[j]] <- Area_Under_Curve(
        subset(robustness.dt, signature == j & type == "Targeted edge attack")$removed.pct,
        subset(robustness.dt, signature == j & type == "Targeted edge attack")$comp.pct
    )
}
(observed_auc_list %>% unlist())
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6 
#   0.7059991   0.7104721   0.7235955   0.7093491   0.7138381   0.7281261 
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18 
#   0.7170029   0.7375547   0.7165139   0.7142507   0.7294497   0.7416438 
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13 
#   0.7513450   0.7326760   0.7323650   0.7349069   0.7455150   0.7239132 
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1 
#   0.7064360   0.6176471   0.4444444   0.6442901   0.6333333   0.5967742 
#    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16    CTS_16.1 
#   0.6071429   0.6300167   0.5586420   0.6465378   0.5522321   0.6590263 
#      CTS_13       CTS_8 
#   0.6808012   0.6021465

# ks.test(subset(robustness.dt, signature=='cardiac.a')$comp.pc, subset(robustness.dt, signature=='cardiac.c')$comp.pc)
# p-value < 2.2e-16

en <- mn <- numeric(length(graph_list)) # safer than array(dim=)
names(en) <- names(mn) <- names(graph_list)

for (j in names(graph_list)) {
    g <- graph_list[[j]]
    deg <- igraph::strength(g, weights = E(g)$weight)

    # Get edgelist as character matrix
    el <- as_edgelist(g, names = TRUE)

    # Sum node degrees for each edge, match names explicitly
    edge_strengths <- apply(el, 1, function(x) {
        if (x[1] %in% names(deg) && x[2] %in% names(deg)) {
            deg[x[1]] + deg[x[2]]
        } else {
            NA
        }
    })

    en[j] <- mean(edge_strengths, na.rm = TRUE) # avg edge strength
    mn[j] <- en[j]
}
options(scipen = 999)
(mn)
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6 
#   1.9879669   5.3753351   3.3102531   2.5006391   2.1490349   2.1797991 
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18 
#   1.9661476   4.5060989   2.6051531   2.4243945   2.1932094   1.0463892 
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13 
#   1.0355188   1.8499141   2.1243071   2.2028725   2.2447135   2.4200731 
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1 
#   2.5406947   1.8618209   0.2297711   1.0335435   0.2959712   0.3919628 
#    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16    CTS_16.1 
#   0.7694620   1.0199373   0.2860569   1.4769951   0.5877318   0.6082691 
#      CTS_13       CTS_8 
#   1.6047819   1.4260508 
options(scipen = 0)

attac_E_random <- readRDS(paste0("AUC_compt.pct_attac_E_random_100runs_", s, "weighted.RDS"))
N <- nrow(attac_E_random)

col_exists <- function(df, col) col %in% colnames(df)

# Get all cluster IDs from CTS_ columns
cluster_ids <- unique(sub("^CTS_", "", grep("^CTS_", colnames(attac_E_random), value = TRUE)))

# Base IDs without suffix (e.g., remove ".1")
cluster_base <- sapply(cluster_ids, function(x) {
    if (grepl(".", x, fixed = TRUE)) strsplit(x, ".", fixed = TRUE)[[1]][1] else x
})
names(cluster_base) <- cluster_ids

pdf("Edge_All_AUC_density_plots.pdf", width = 20, height = 4)
par(
    mfrow = c(1, 7),
    mar = c(4.5, 4.5, 4.2, 2.5), # bottom, left, top, right
    oma = c(2, 1, 2, 1) # increased top and bottom margins
)

for (j in cluster_ids) {
    j0 <- cluster_base[j]

    cts_col <- paste0("CTS_", j)
    higcts_col <- paste0("HiGCTS_", j)
    hig_col <- paste0("HiG_", j0)

    if (!any(c(cts_col, higcts_col, hig_col) %in% colnames(attac_E_random))) {
        message("Skipping cluster ", j, ": no relevant columns found.")
        next
    }

    obs_cols <- c()
    lf_groups <- c()
    x_lines <- numeric(0)

    if (col_exists(attac_E_random, cts_col)) {
        obs_cols <- c(obs_cols, cts_col)
        lf_groups <- c(lf_groups, rep("CTS", N))
        x_lines <- c(x_lines, attac_E_random[, cts_col])
    }
    if (col_exists(attac_E_random, higcts_col)) {
        obs_cols <- c(obs_cols, higcts_col)
        lf_groups <- c(lf_groups, rep("HiGCTS", N))
        x_lines <- c(x_lines, attac_E_random[, higcts_col])
    }
    if (col_exists(attac_E_random, hig_col)) {
        obs_cols <- c(obs_cols, hig_col)
        lf_groups <- c(lf_groups, rep("HiG", N))
        x_lines <- c(x_lines, attac_E_random[, hig_col])
    }

    lf <- factor(lf_groups, levels = c("CTS", "HiGCTS", "HiG"))

    observed_values <- sapply(obs_cols, function(col) {
        if (!is.null(observed_auc_list[[col]])) observed_auc_list[[col]] else NA
    })
    xlim_range <- range(c(x_lines, observed_values), na.rm = TRUE)

    sm.density.compare(
        x_lines,
        lf,
        xlab = "AUC of maximal component size",
        col = PPI_color_palette,
        xlim = xlim_range,
        lwd = 2,
        main = paste0("Cluster ", j)
    )

    for (col in obs_cols) {
        if (!is.null(observed_auc_list[[col]])) {
            label <- sub("^(CTS|HiGCTS|HiG)_.*", "\\1", col)
            abline(v = observed_auc_list[[col]], col = PPI_color_palette[label], lty = 2, lwd = 2)
        }
    }

    if (col_exists(attac_E_random, cts_col) && col_exists(attac_E_random, hig_col)) {
        ks_p <- ks.test(attac_E_random[, cts_col], attac_E_random[, hig_col])$p.value
        mtext(sprintf("KS p = %.2g", ks_p), side = 1, line = 5) # moved lower for spacing
    }

    for (label in c("CTS", "HiGCTS", "HiG")) {
        colname <- paste0(label, "_", ifelse(label == "HiG", j0, j))
        if (col_exists(attac_E_random, colname) && !is.null(observed_auc_list[[colname]])) {
            above <- sum(attac_E_random[, colname] > observed_auc_list[[colname]])
            below <- sum(attac_E_random[, colname] < observed_auc_list[[colname]])
            txt_cluster <- ifelse(label == "HiG", j0, j)
            mtext(sprintf("%s_%s p= %.2f ; %.2f", label, txt_cluster, above / N, below / N),
                side = 3,
                line = switch(label,
                    "CTS" = 1,
                    "HiGCTS" = 2.5,
                    "HiG" = 4
                ),
                col = PPI_color_palette[label]
            )
        }
    }
}

dev.off()
