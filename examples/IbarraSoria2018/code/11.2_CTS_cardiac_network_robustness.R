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

wd = "/Users/felixyu/Documents/IbarraSoria2018/"
source(paste0(wd, "code/celltype_specific_weight_v9.R"))

PPI_color_platte <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_platte <- c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)

setwd(paste0(wd, "/results/PPI_weight/"))

# refer to 11.1_CTS_cardiac_network_strengthDistribution.R
s <- "combined"
file <- paste0("2018_STRING_graph_perState_simplified_", s, "weighted.rds")
graph_list <- readRDS(file)
(names(graph_list))
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"         
#  [3] "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [5] "HiG_blood"                  "HiG_mesodermProgenitors"   
#  [7] "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"  
#  [9] "HiG_somiticMesoderm"        "HiG_mixedMesoderm.a"       
# [11] "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"             
# [15] "HiG_endothelial.b"          "HiG_cardiac.a"             
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"          
# [19] "CTS_endothelial.b"          "CTS_cardiac.a" 
edge_counts <- sapply(graph_list, ecount)
(edge_counts)
# HiG_extraembryonicMesoderm          HiG_endothelial.a 
#                       6743                       6743 
#          HiG_endothelial.c          HiG_endothelial.d 
#                      15247                      15247 
#                  HiG_blood    HiG_mesodermProgenitors 
#                      11724                      11724 
#   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                      12662                      12662 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a 
#                      20602                      20602 
#     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                       7557                       7557 
#              HiG_cardiac.b              HiG_cardiac.c 
#                      12912                      12912 
#          HiG_endothelial.b              HiG_cardiac.a 
#                       6608                       6608 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a 
#                         20                         28 
#          CTS_endothelial.b              CTS_cardiac.a 
#                         82                         79 

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
# [1] "combined_score"  "weight"          "width"           "original_weight"
# [5] "corexp_sign"     "coexp_target"
(all(E(graph_list[[1]])$combined_score / 1000 == E(graph_list[[1]])$original_weight)) # TRUE
(all(E(graph_list[[1]])$weight == E(graph_list[[1]])$original_weight)) # FALSE

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
    scale_color_manual(values = PPI_color_platte) +
    scale_size_manual(values = PPI_size_platte) +
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
    scale_size_manual(values = PPI_size_platte) +
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
    scale_color_manual(values = PPI_color_platte) +
    scale_size_manual(values = PPI_size_platte) +
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
    scale_size_manual(values = PPI_size_platte) +
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

   pdf(file='vertex_attack_IbarraSoria2018.pdf') 
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
file <- "2018_STRING_graph_perState_simplified_combinedweighted.rds"
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"         
#  [3] "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [5] "HiG_blood"                  "HiG_mesodermProgenitors"   
#  [7] "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"  
#  [9] "HiG_somiticMesoderm"        "HiG_mixedMesoderm.a"       
# [11] "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"             
# [15] "HiG_endothelial.b"          "HiG_cardiac.a"             
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"          
# [19] "CTS_endothelial.b"          "CTS_cardiac.a"   


attack.vertex.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "vertex", "btwn.cent"), idcol = names(graph_list))
colnames(attack.vertex.btwn)[1] <- "signature"
(head(attack.vertex.btwn, 3))
#                     signature                   type   measure comp.size comp.pct removed.pct
#                        <char>                 <char>    <char>     <num>    <num>       <num>
# 1: HiG_extraembryonicMesoderm Targeted vertex attack btwn.cent       409 1.000000 0.000000000
# 2: HiG_extraembryonicMesoderm Targeted vertex attack btwn.cent       408 0.997555 0.002415459
# 3: HiG_extraembryonicMesoderm Targeted vertex attack btwn.cent       407 0.995110 0.004830918

(dim(attack.vertex.btwn)) #  8135    6
saveRDS(attack.vertex.btwn, file = "attack.vertex.btwn.rds")

# calculate the 'edge' btwn for each G first
# igraph::edge_betweenness() uses distance graph weights, but E(g) uses connection weights, thus we invert it.
for (j in seq_along(graph_list)) E(graph_list[[j]])$btwn <- edge_betweenness(graph_list[[j]], weights = 1/E(graph_list[[j]])$weight)
## then do the edge-attack analysis
attack.edge.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "edge", "btwn.cent"), idcol = names(graph_list)) # !!!!!!!!!!!!!
colnames(attack.edge.btwn)[1] <- "signature"
(dim(attack.edge.btwn)) #  188339      6
(head(attack.edge.btwn, 3))
#                     signature                 type   measure comp.size comp.pct  removed.pct
#                        <char>               <char>    <char>     <num>    <num>        <num>
# 1: HiG_extraembryonicMesoderm Targeted edge attack btwn.cent       409        1 0.0000000000
# 2: HiG_extraembryonicMesoderm Targeted edge attack btwn.cent       409        1 0.0001483019
# 3: HiG_extraembryonicMesoderm Targeted edge attack btwn.cent       409        1 0.0002966039

(table(attack.edge.btwn$signature))
#                     signature                 type   measure comp.size comp.pct  removed.pct
#                        <char>               <char>    <char>     <num>    <num>        <num>
# 1: HiG_extraembryonicMesoderm Targeted edge attack btwn.cent       409        1 0.0000000000
# 2: HiG_extraembryonicMesoderm Targeted edge attack btwn.cent       409        1 0.0001483019
# 3: HiG_extraembryonicMesoderm Targeted edge attack btwn.cent       409        1 0.0002966039

saveRDS(attack.edge.btwn, file = "attack.edge.btwn.rds")

# In a random failure analysis, you choose a vertex/edge at random, remove it, and calculate the maximum
#  component size until all elements have been removed.
####################################################
# refer to 11.2_CTS_cardiac_network_robustness_Miudway3.simulation.R
## run on Midway3, the following is long. DO NOT REPEAT !!!)


(sapply(graph_list, vcount))
#              CTS_cardiac.a          CTS_endothelial.b                  HiG_blood              HiG_cardiac.a 
#                         80                         83                      11725                       6609 
#              HiG_cardiac.b              HiG_cardiac.c          HiG_endothelial.a          HiG_endothelial.b 
#                      12913                      12913                       6744                       6609 
#          HiG_endothelial.c          HiG_endothelial.d HiG_extraembryonicMesoderm    HiG_mesodermProgenitors 
#                      15248                      15248                       6744                      11725 
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a 
#                      20603                       7558                       7558                      12663 
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm           HiGCTS_cardiac.a       HiGCTS_endothelial.b 
#                      12663                      20603                         29                         21 
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

attack.edge.btwn <- readRDS(file = "attack.edge.btwn.rds")
attack.vertex.btwn <- readRDS(file = "attack.vertex.btwn.rds")
failure.vertex <- readRDS(paste0("failure.vertex_100_simplified_", s, "weighted.rds"))
(table(attack.vertex.btwn$signature, attack.vertex.btwn$type))
#                              Targeted vertex attack
#   CTS_cardiac.a                                  38
#   CTS_endothelial.b                              34
#   HiG_blood                                     560
#   HiG_cardiac.a                                 393
#   HiG_cardiac.b                                 550
#   HiG_cardiac.c                                 550
#   HiG_endothelial.a                             415
#   HiG_endothelial.b                             393
#   HiG_endothelial.c                             574
#   HiG_endothelial.d                             574
#   HiG_extraembryonicMesoderm                    415
#   HiG_mesodermProgenitors                       560
#   HiG_mixedMesoderm.a                           553
#   HiG_mixedMesoderm.b                           433
#   HiG_pharyngealMesoderm                        433
#   HiG_presomiticMesoderm.a                      538
#   HiG_presomiticMesoderm.b                      538
#   HiG_somiticMesoderm                           553
#   HiGCTS_cardiac.a                               14
#   HiGCTS_endothelial.b                           17

(subset(attack.vertex.btwn, signature == "HiGCTS_cardiac.a"))
#            signature                   type   measure comp.size   comp.pct removed.pct
#               <char>                 <char>    <char>     <num>      <num>       <num>
#  1: HiGCTS_cardiac.a Targeted vertex attack btwn.cent        12 1.00000000  0.00000000
#  2: HiGCTS_cardiac.a Targeted vertex attack btwn.cent        10 0.83333333  0.07692308
#  3: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         9 0.75000000  0.15384615
#  4: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         8 0.66666667  0.23076923
#  5: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         6 0.50000000  0.30769231
#  6: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         5 0.41666667  0.38461538
#  7: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         3 0.25000000  0.46153846
#  8: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         3 0.25000000  0.53846154
#  9: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         3 0.25000000  0.61538462
# 10: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         3 0.25000000  0.69230769
# 11: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         2 0.16666667  0.76923077
# 12: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1 0.08333333  0.84615385
# 13: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1 0.08333333  0.92307692
# 14: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         0 0.00000000  1.00000000


failure.edge <- readRDS(paste0("failure.edge_100_simplified_", s, "weighted.rds"))
failure.dt <- rbind(failure.edge, failure.vertex)
(head(failure.dt, 3))
#    HiG_extraembryonicMesoderm                type measure comp.size comp.pct  removed.pct
#                        <char>              <char>  <char>     <num>    <num>        <num>
# 1: HiG_extraembryonicMesoderm Random edge removal  random       409        1 0.0000000000
# 2: HiG_extraembryonicMesoderm Random edge removal  random       409        1 0.0001483019
# 3: HiG_extraembryonicMesoderm Random edge removal  random       409        1 0.0002966039

colnames(failure.dt)[1] <- "signature"
(table(failure.dt$signature, failure.dt$type))
#                              Random edge removal Random vertex removal
#   CTS_cardiac.a                               80                    38
#   CTS_endothelial.b                           83                    34
#   HiG_blood                                11725                   560
#   HiG_cardiac.a                             6609                   393
#   HiG_cardiac.b                            12913                   550
#   HiG_cardiac.c                            12913                   550
#   HiG_endothelial.a                         6744                   415
#   HiG_endothelial.b                         6609                   393
#   HiG_endothelial.c                        15248                   574
#   HiG_endothelial.d                        15248                   574
#   HiG_extraembryonicMesoderm                6744                   415
#   HiG_mesodermProgenitors                  11725                   560
#   HiG_mixedMesoderm.a                      20603                   553
#   HiG_mixedMesoderm.b                       7558                   433
#   HiG_pharyngealMesoderm                    7558                   433
#   HiG_presomiticMesoderm.a                 12663                   538
#   HiG_presomiticMesoderm.b                 12663                   538
#   HiG_somiticMesoderm                      20603                   553
#   HiGCTS_cardiac.a                            29                    14
#   HiGCTS_endothelial.b                        21                    17


colnames(attack.vertex.btwn)[1] <- "signature"
(dim(failure.dt)) #  196474      6

robustness.dt <- rbind(failure.dt, attack.vertex.btwn, attack.edge.btwn)
(dim(robustness.dt)) #  392948      6
robustness.dt$PPI_cat <- lapply(robustness.dt$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))
(head(robustness.dt, 3))
#                     signature                type measure comp.size comp.pct  removed.pct PPI_cat
#                        <char>              <char>  <char>     <num>    <num>        <num>  <fctr>
# 1: HiG_extraembryonicMesoderm Random edge removal  random       409        1 0.0000000000     HiG
# 2: HiG_extraembryonicMesoderm Random edge removal  random       409        1 0.0001483019     HiG
# 3: HiG_extraembryonicMesoderm Random edge removal  random       409        1 0.0002966039     HiG

robustness.dt$experiment <- ifelse(grepl("edge", robustness.dt$type), "edge", "vertex")
robustness.dt$measure <- factor(robustness.dt$measure, levels = c("random", "btwn.cent"))

(table(robustness.dt$type, robustness.dt$measure))
#                          random btwn.cent
#   Random edge removal    188339         0
#   Random vertex removal    8135         0
#   Targeted edge attack        0    188339
#   Targeted vertex attack      0      8135


robustness.dt$type <- factor(robustness.dt$type,
    levels = c("Random edge removal", "Targeted edge attack", "Random vertex removal", "Targeted vertex attack")
)

robustness.dt$cluster <- lapply(robustness.dt$signature, function(x) unlist(strsplit(x, "_"))[2]) %>% unlist()

p_attack3 <- ggplot(
    robustness.dt,
    aes(x = removed.pct, y = comp.pct, col = cluster, linetype = PPI_cat, size = PPI_cat)
) +
    geom_line(show.legend = TRUE) +
    scale_size_manual(values = PPI_size_platte, name = "PPI Category") +
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
    scale_color_manual(values = PPI_color_platte, name = "PPI Category") +
    scale_size_manual(values = PPI_size_platte, name = "PPI Category") +
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
pdf(file = "attack_IbarraSoria2018.pdf", width = 8, height = 7)
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
# 1        edge     CTS    random c(W = 0..... 1.505201e-08
# 2        edge     CTS btwn.cent c(W = 0..... 2.501683e-06
# 3        edge  HiGCTS    random c(W = 0..... 4.524175e-04
# 4        edge  HiGCTS btwn.cent c(W = 0..... 2.634337e-03
# 5        edge     HiG    random           NA           NA
# 6        edge     HiG btwn.cent           NA           NA
# 7      vertex     CTS    random c(W = 0..... 3.336774e-04
# 8      vertex     CTS btwn.cent c(W = 0..... 5.054868e-08
# 9      vertex  HiGCTS    random c(W = 0..... 6.787811e-02
# 10     vertex  HiGCTS btwn.cent c(W = 0..... 1.073263e-04
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
# 1 CTS                1.05                1.36
# 2 HiGCTS             1.01                1.34
# 3 HiG                0.956               1.04

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
    scale_color_manual(values = PPI_color_platte) +
    scale_size_manual(values = PPI_size_platte) +
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
    scale_fill_manual(values = PPI_color_platte) +
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

pdf(file = "box_wilcox-test_attack_IbarraSoria2018.pdf", width = 12)
print(g)
print(g2)
dev.off()

# confirm that the figure showing adjusted p-values
wilcox.test(
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "CTS")],
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiGCTS")]
)
# W = 990, p-value = 0.366

wilcox.test(
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "CTS")],
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiG")]
)
# W = 209550, p-value = 5.632e-05

wilcox.test(
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiGCTS")],
    targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiG")]
)
# W = 93293, p-value = 0.01586

## finally, manually add the threshold of fold changes  ###############
robustness.dt <- robustness.dt %>%
    group_by(PPI_cat, type) %>%
    mutate(
        mean_comp_pct = mean(comp.pct, na.rm = TRUE) # Calculate mean comp.pct for each group
    ) %>%
    ungroup()

## for each signature, compare the comp.pct between targeted attack to its random removal to access the significance of 'hub'
#######################################
Hub_effect <- array(dim = length(graph_list))
names(Hub_effect) <- names(graph_list)
for (j in names(graph_list)) {
    tmp <- subset(robustness.dt, signature == j)##NEW
    (dim(tmp)) # 14318    11
    (table(tmp$type))
    #    Random edge removal   Targeted edge attack  Random vertex removal Targeted vertex attack 
    #                   6744                   6744                    415                    415 
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
#                                  edge_p   edge_p.adj    node_p node_p.adj
# HiG_extraembryonicMesoderm 1.395868e-18 2.791735e-17 0.1409423          1
# HiG_endothelial.a          1.039399e-15 2.078797e-14 0.1538366          1
# HiG_endothelial.c          1.983482e-44 3.966964e-43 0.3560730          1
# HiG_endothelial.d          1.285423e-46 2.570847e-45 0.3404282          1
# HiG_blood                  6.022794e-31 1.204559e-29 0.3069206          1
# HiG_mesodermProgenitors    6.039746e-32 1.207949e-30 0.2672300          1
# HiG_presomiticMesoderm.b   3.298374e-36 6.596748e-35 0.3866887          1
# HiG_presomiticMesoderm.a   6.683827e-42 1.336765e-40 0.2701798          1
# HiG_somiticMesoderm        1.787474e-47 3.574947e-46 0.4347686          1
# HiG_mixedMesoderm.a        6.914545e-54 1.382909e-52 0.4030829          1
# HiG_pharyngealMesoderm     5.026050e-40 1.005210e-38 0.1429150          1
# HiG_mixedMesoderm.b        8.289912e-37 1.657982e-35 0.1841956          1
# HiG_cardiac.b              1.945315e-42 3.890630e-41 0.2279994          1
# HiG_cardiac.c              7.699959e-52 1.539992e-50 0.2917282          1
# HiG_endothelial.b          1.381454e-23 2.762908e-22 0.2490499          1
# HiG_cardiac.a              3.647988e-19 7.295975e-18 0.2035765          1

(df[which(df$node_p.adj < 0.05), ]) # No networks had significant hub effects at the node level
# [1] edge_p     edge_p.adj node_p     node_p.adj
# <0 rows> (or 0-length row.names)




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
    scale_color_manual(values = PPI_color_platte) +
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
    scale_color_manual(values = PPI_color_platte) +
    scale_shape_manual(values = c("HiG" = 16, "CTS" = 17, "HiGCTS" = 18)) +
    facet_wrap(~type, scales = "free", ncol = 2) + # Split into two rows for edge and node
    theme_minimal() +
    theme(legend.position = "top") +
    ggtitle("Adjusted p-values of Wilcox test comparing targeted attack vs random removal")
print(g2)
pdf(file = "line.adjp_wilcox_attack_IbarraSoria2018.pdf")
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

IDs_of_CTS <- c("endothelial.b", "cardiac.a")
(names(graph_list))
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [5] "HiG_blood"                  "HiG_mesodermProgenitors"    "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"  
#  [9] "HiG_somiticMesoderm"        "HiG_mixedMesoderm.a"        "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"              "HiG_endothelial.b"          "HiG_cardiac.a"             
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"           "CTS_endothelial.b"          "CTS_cardiac.a"

observed_auc_list <- list()
for (j in names(graph_list)) {
    observed_auc_list[[j]] <- Area_Under_Curve(
        subset(robustness.dt, signature == j & type == "Targeted vertex attack")$removed.pct,
        subset(robustness.dt, signature == j & type == "Targeted vertex attack")$comp.pct
    )
}
(observed_auc_list %>% unlist())
# HiG_extraembryonicMesoderm          HiG_endothelial.a          HiG_endothelial.c          HiG_endothelial.d 
#                  0.4608064                  0.4616981                  0.4790393                  0.4785429 
#                  HiG_blood    HiG_mesodermProgenitors   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                  0.4738400                  0.4720769                  0.4776363                  0.4736443 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                  0.4812117                  0.4801168                  0.4621463                  0.4656023 
#              HiG_cardiac.b              HiG_cardiac.c          HiG_endothelial.b              HiG_cardiac.a 
#                  0.4728770                  0.4752303                  0.4668530                  0.4644586 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a          CTS_endothelial.b              CTS_cardiac.a 
#                  0.3096591                  0.3461538                  0.2983683                  0.3396396 

df_AUC <- data.frame(
    auc = observed_auc_list %>% unlist(),
    signature = names(observed_auc_list),
    PPI_cat = lapply(names(observed_auc_list), function(x) unlist(strsplit(x, split = "_"))[1]) %>% unlist()
)
df_AUC$PPI_cat <- factor(df_AUC$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

g3 <- ggplot(df_AUC, aes(x = PPI_cat, y = auc, fill = PPI_cat)) +
    geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) + # Dodge the boxes for each type per PPI_cat
    scale_fill_manual(values = PPI_color_platte) +
    geom_signif(
        comparisons = list(
            c("CTS", "HiGCTS"), c("HiGCTS", "HiG"), c("CTS", "HiG")
        ),
        map_signif_level = TRUE,
        step_increase = 0.1, # Adjusts spacing between the lines
        test = "wilcox.test", # Perform a t-test to calculate significance
        p.adjust.method = "holm" # default is
    ) +
    theme_minimal() +
    labs(
        title = "Comparison of Robustness Measures by cent.btw",
        x = "PPI Category",
        y = "AUC of targeted vertex attack"
    ) +
    theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5) # center title
    )
print(g3)
dev.copy2pdf(file = "box_wilcox-test_attack_AUC_IbarraSoria2018.pdf", width = 10)


vn <- sapply(graph_list, vcount)
pn <- array(dim = length(graph_list))
names(pn) <- names(graph_list)
for (j in names(vn)) {
    pn[j] <- igraph::strength(graph_list[[j]], weights = E(graph_list[[j]])$weight) %>% mean() / vn[j]
}
(vn)
# HiG_extraembryonicMesoderm          HiG_endothelial.a          HiG_endothelial.c          HiG_endothelial.d                  HiG_blood    HiG_mesodermProgenitors   HiG_presomiticMesoderm.b 
#                  0.4608064                  0.4616981                  0.4790393                  0.4785429                  0.4738400                  0.4720769                  0.4776363 
#   HiG_presomiticMesoderm.a        HiG_somiticMesoderm        HiG_mixedMesoderm.a     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b              HiG_cardiac.b              HiG_cardiac.c 
#                  0.4736443                  0.4812117                  0.4801168                  0.4621463                  0.4656023                  0.4728770                  0.4752303 
#          HiG_endothelial.b              HiG_cardiac.a       HiGCTS_endothelial.b           HiGCTS_cardiac.a          CTS_endothelial.b              CTS_cardiac.a 
#                  0.4668530                  0.4644586                  0.3096591                  0.3461538                  0.2983683                  0.3396396 
# HiG_extraembryonicMesoderm          HiG_endothelial.a          HiG_endothelial.c          HiG_endothelial.d                  HiG_blood    HiG_mesodermProgenitors   HiG_presomiticMesoderm.b 
#                        414                        414                        573                        573                        559                        559                        537 
#   HiG_presomiticMesoderm.a        HiG_somiticMesoderm        HiG_mixedMesoderm.a     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b              HiG_cardiac.b              HiG_cardiac.c 
#                        537                        552                        552                        432                        432                        549                        549 
#          HiG_endothelial.b              HiG_cardiac.a       HiGCTS_endothelial.b           HiGCTS_cardiac.a          CTS_endothelial.b              CTS_cardiac.a 
#                        392                        392                         16                         13                         33                         37 


(sum(vn[grep("^HiG_", names(vn))]))     # 8016
(sum(vn[grep("^HiGCTS_", names(vn))]))  # 29
(sum(vn[grep("^CTS_", names(vn))]))     # 70

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
        col = PPI_color_platte,
        xlim = xlim_range,
        lwd = 2,
        main = paste0("Cluster ", j)
    )

    # Add vertical lines for observed values
    for (col in obs_cols) {
        if (!is.null(observed_auc_list[[col]])) {
            label <- sub("^(CTS|HiGCTS|HiG)_.*", "\\1", col)
            abline(v = observed_auc_list[[col]], col = PPI_color_platte[label], lty = 2, lwd = 2)
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
                col = PPI_color_platte[label]
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
# HiG_extraembryonicMesoderm          HiG_endothelial.a          HiG_endothelial.c          HiG_endothelial.d 
#                  0.6925460                  0.6909190                  0.7079443                  0.7094227 
#                  HiG_blood    HiG_mesodermProgenitors   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                  0.6988889                  0.6998049                  0.6997179                  0.7023177 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                  0.7053271                  0.7075465                  0.7206489                  0.7177318 
#              HiG_cardiac.b              HiG_cardiac.c          HiG_endothelial.b              HiG_cardiac.a 
#                  0.6967657                  0.6999990                  0.6972624                  0.6929412 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a          CTS_endothelial.b              CTS_cardiac.a 
#                  0.6113636                  0.7023810                  0.6308630                  0.6126582 

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
# HiG_extraembryonicMesoderm          HiG_endothelial.a          HiG_endothelial.c          HiG_endothelial.d 
#                   8.665859                   7.146630                  28.794115                  28.121675 
#                  HiG_blood    HiG_mesodermProgenitors   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                  17.625680                  16.383867                  30.028392                  29.370176 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                  48.241201                  48.712559                  17.298067                  17.401475 
#              HiG_cardiac.b              HiG_cardiac.c          HiG_endothelial.b              HiG_cardiac.a 
#                  26.305750                  27.133499                  12.910587                  12.530216 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a          CTS_endothelial.b              CTS_cardiac.a 
#                   1.359660                   1.187108                   1.566885                   1.129949  

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
        col = PPI_color_platte,
        xlim = xlim_range,
        lwd = 2,
        main = paste0("Cluster ", j)
    )

    for (col in obs_cols) {
        if (!is.null(observed_auc_list[[col]])) {
            label <- sub("^(CTS|HiGCTS|HiG)_.*", "\\1", col)
            abline(v = observed_auc_list[[col]], col = PPI_color_platte[label], lty = 2, lwd = 2)
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
                col = PPI_color_platte[label]
            )
        }
    }
}
dev.off()
