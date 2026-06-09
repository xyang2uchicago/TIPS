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

wd <- "/Users/felixyu/Documents/GitHub/TIPS/examples/IbarraSoria2018/IbarraSoria2018_IID/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- file.path(wd, "data/")

celltype_specific_weight_version <- "10"
source("/Users/felixyu/Documents/GitHub/TIPS/R/celltype_specific_weight_v10.R")

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_palette <- c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)

db <- "IbarraSoria2018"

s <- "combined" # specificity method

########## END OF USER INPUT ##########

file <- file.path(inputdir, paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds"))
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
# [1] "weight"         "norm_PPI_score" "corexp_sign"    "coexp_focal"
(all(E(graph_list[[1]])$combined_score / 1000 == E(graph_list[[1]])$original_weight)) # TRUE
(all(E(graph_list[[1]])$weight == E(graph_list[[1]])$original_weight)) # TRUE

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

pdf(file = "vertex_attack_IbarraSoria2018.pdf")
print(gridExtra::grid.arrange(g_attack, g_attack2,
  gsignature_attack, gsignature_attack2,
  ncol = 2, nrow = 2
))
dev.off()

############## to restart here anytime #################################
library(brainGraph)
require(dplyr)
library(data.table)
library(ggplot2)

# refer to 11.1.1_CTS_cardiac_network_strengthDistribution.R
file <- file.path(inputDir, paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds"))
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


attack.vertex.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "vertex", "btwn.cent"), idcol = names(graph_list))
colnames(attack.vertex.btwn)[1] <- "signature"
(head(attack.vertex.btwn, 5))
#    signature                   type   measure comp.size  comp.pct removed.pct
#       <char>                 <char>    <char>     <num>     <num>       <num>
# 1: HiG_blood Targeted vertex attack btwn.cent       360 1.0000000 0.000000000
# 2: HiG_blood Targeted vertex attack btwn.cent       356 0.9888889 0.002673797
# 3: HiG_blood Targeted vertex attack btwn.cent       355 0.9861111 0.005347594
# 4: HiG_blood Targeted vertex attack btwn.cent       354 0.9833333 0.008021390
# 5: HiG_blood Targeted vertex attack btwn.cent       353 0.9805556 0.010695187

(dim(attack.vertex.btwn)) #  6002    6
saveRDS(attack.vertex.btwn, file = "attack.vertex.btwn.rds")

# calculate the 'edge' btwn for each G first
# igraph::edge_betweenness() uses distance graph weights, but E(g) uses connection weights, thus we invert it.
for (j in seq_along(graph_list)) E(graph_list[[j]])$btwn <- edge_betweenness(graph_list[[j]], weights = 1 / E(graph_list[[j]])$weight)
## then do the edge-attack analysis
attack.edge.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "edge", "btwn.cent"), idcol = names(graph_list)) # !!!!!!!!!!!!!
colnames(attack.edge.btwn)[1] <- "signature"
(dim(attack.edge.btwn)) #  48446      6
(head(attack.edge.btwn, 3))
#    signature                 type   measure comp.size comp.pct  removed.pct
#       <char>               <char>    <char>     <num>    <num>        <num>
# 1: HiG_blood Targeted edge attack btwn.cent       360        1 0.0000000000
# 2: HiG_blood Targeted edge attack btwn.cent       360        1 0.0003207184
# 3: HiG_blood Targeted edge attack btwn.cent       360        1 0.0006414368

(table(attack.edge.btwn$signature))
#              CTS_cardiac.a          CTS_endothelial.b
#                         18                          6
#                  HiG_blood              HiG_cardiac.a
#                       3119                       3073
#              HiG_cardiac.b              HiG_cardiac.c
#                       3499                       3765
#          HiG_endothelial.a          HiG_endothelial.b
#                       4383                       3157
#          HiG_endothelial.c          HiG_endothelial.d
#                       3967                       3248
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                       2655                       2188
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                       2135                       2751
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                       2718                       2130
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                       3463                       2163
#           HiGCTS_cardiac.a       HiGCTS_endothelial.b
#                          6                          2

saveRDS(attack.edge.btwn, file = "attack.edge.btwn.rds")

# In a random failure analysis, you choose a vertex/edge at random, remove it, and calculate the maximum
#  component size until all elements have been removed.
####################################################
# refer to 11.2_CTS_cardiac_network_robustness_Miudway3.simulation.R
## run on Midway3, the following is long. DO NOT REPEAT !!!)


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

attack.edge.btwn <- readRDS(file = "attack.edge.btwn.rds")
attack.vertex.btwn <- readRDS(file = "attack.vertex.btwn.rds")
failure.vertex <- readRDS(paste0("failure.vertex_100_simplified_", s, "weighted.rds"))
(table(attack.vertex.btwn$signature, attack.vertex.btwn$type))
#                              Targeted vertex attack
#   CTS_cardiac.a                                  38
#   CTS_endothelial.b                              28
#   HiG_blood                                     375
#   HiG_cardiac.a                                 385
#   HiG_cardiac.b                                 392
#   HiG_cardiac.c                                 422
#   HiG_endothelial.a                             425
#   HiG_endothelial.b                             390
#   HiG_endothelial.c                             448
#   HiG_endothelial.d                             393
#   HiG_extraembryonicMesoderm                    346
#   HiG_mesodermProgenitors                       331
#   HiG_mixedMesoderm.a                           310
#   HiG_mixedMesoderm.b                           344
#   HiG_pharyngealMesoderm                        346
#   HiG_presomiticMesoderm.a                      316
#   HiG_presomiticMesoderm.b                      380
#   HiG_somiticMesoderm                           306
#   HiGCTS_cardiac.a                               14
#   HiGCTS_endothelial.b                           13

(subset(attack.vertex.btwn, signature == "HiGCTS_cardiac.a"))
#            signature                   type   measure comp.size comp.pct
#               <char>                 <char>    <char>     <num>    <num>
#  1: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         4     1.00
#  2: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         2     0.50
#  3: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         2     0.50
#  4: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         2     0.50
#  5: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         2     0.50
#  6: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1     0.25
#  7: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1     0.25
#  8: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1     0.25
#  9: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1     0.25
# 10: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1     0.25
# 11: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1     0.25
# 12: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1     0.25
# 13: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         1     0.25
# 14: HiGCTS_cardiac.a Targeted vertex attack btwn.cent         0     0.00
#     removed.pct
#           <num>
#  1:  0.00000000
#  2:  0.07692308
#  3:  0.15384615
#  4:  0.23076923
#  5:  0.30769231
#  6:  0.38461538
#  7:  0.46153846
#  8:  0.53846154
#  9:  0.61538462
# 10:  0.69230769
# 11:  0.76923077
# 12:  0.84615385
# 13:  0.92307692
# 14:  1.00000000


failure.edge <- readRDS(paste0("failure.edge_100_simplified_", s, "weighted.rds"))
failure.dt <- rbind(failure.edge, failure.vertex)
(head(failure.dt, 3))
#    HiG_blood                type measure comp.size comp.pct  removed.pct
#       <char>              <char>  <char>     <num>    <num>        <num>
# 1: HiG_blood Random edge removal  random       360        1 0.0000000000
# 2: HiG_blood Random edge removal  random       360        1 0.0003207184
# 3: HiG_blood Random edge removal  random       360        1 0.0006414368

colnames(failure.dt)[1] <- "signature"
(table(failure.dt$signature, failure.dt$type))
#                              Random edge removal Random vertex removal
#   CTS_cardiac.a                               18                    38
#   CTS_endothelial.b                            6                    28
#   HiG_blood                                 3119                   375
#   HiG_cardiac.a                             3073                   385
#   HiG_cardiac.b                             3499                   392
#   HiG_cardiac.c                             3765                   422
#   HiG_endothelial.a                         4383                   425
#   HiG_endothelial.b                         3157                   390
#   HiG_endothelial.c                         3967                   448
#   HiG_endothelial.d                         3248                   393
#   HiG_extraembryonicMesoderm                2655                   346
#   HiG_mesodermProgenitors                   2188                   331
#   HiG_mixedMesoderm.a                       2135                   310
#   HiG_mixedMesoderm.b                       2751                   344
#   HiG_pharyngealMesoderm                    2718                   346
#   HiG_presomiticMesoderm.a                  2130                   316
#   HiG_presomiticMesoderm.b                  3463                   380
#   HiG_somiticMesoderm                       2163                   306
#   HiGCTS_cardiac.a                             6                    14
#   HiGCTS_endothelial.b                         0                    13


colnames(attack.vertex.btwn)[1] <- "signature"
(dim(failure.dt)) #  54446      6

robustness.dt <- rbind(failure.dt, attack.vertex.btwn, attack.edge.btwn)
(dim(robustness.dt)) #  108894      6
robustness.dt$PPI_cat <- lapply(robustness.dt$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
  unlist() %>%
  factor(., levels = c("CTS", "HiGCTS", "HiG"))
(head(robustness.dt, 3))
#    signature                type measure comp.size comp.pct  removed.pct
#       <char>              <char>  <char>     <num>    <num>        <num>
# 1: HiG_blood Random edge removal  random       360        1 0.0000000000
# 2: HiG_blood Random edge removal  random       360        1 0.0003207184
# 3: HiG_blood Random edge removal  random       360        1 0.0006414368
#    PPI_cat
#     <fctr>
# 1:     HiG
# 2:     HiG
# 3:     HiG

robustness.dt$experiment <- ifelse(grepl("edge", robustness.dt$type), "edge", "vertex")
robustness.dt$measure <- factor(robustness.dt$measure, levels = c("random", "btwn.cent"))

(table(robustness.dt$type, robustness.dt$measure))
#                          random btwn.cent
#   Random edge removal     48444         0
#   Random vertex removal    6002         0
#   Targeted edge attack        0     48446
#   Targeted vertex attack      0      6002


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
# 1        edge     CTS    random c(W = 0..... 9.856623e-03
# 2        edge     CTS btwn.cent c(W = 0..... 1.190103e-01
# 3        edge  HiGCTS    random c(W = 0..... 3.887670e-01
# 4        edge  HiGCTS btwn.cent c(W = 0..... 2.020049e-01
# 5        edge     HiG    random           NA           NA
# 6        edge     HiG btwn.cent           NA           NA
# 7      vertex     CTS    random c(W = 0..... 3.397621e-02
# 8      vertex     CTS btwn.cent c(W = 0..... 1.480490e-11
# 9      vertex  HiGCTS    random c(W = 0..... 5.386657e-01
# 10     vertex  HiGCTS btwn.cent c(W = 0..... 1.845329e-04
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
# 1 CTS                1.16                2.08
# 2 HiGCTS             1.18                1.32
# 3 HiG                0.950               1.36

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

pdf(file = "box_wilcox-test_attack_IbarraSoria2018.pdf", width = 12)
print(g)
print(g2)
dev.off()

# confirm that the figure showing adjusted p-values
wilcox.test(
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "CTS")],
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiGCTS")]
)
# W = 304, p-value = 4.071e-07

wilcox.test(
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "CTS")],
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiG")]
)
# W = 185913, p-value = 0.5145

wilcox.test(
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiGCTS")],
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiG")]
)
# W = 92551, p-value = 0.1503

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
  tmp <- subset(robustness.dt, signature == j) ## NEW
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
#                                  edge_p   edge_p.adj       node_p   node_p.adj
# HiG_blood                  7.467968e-12 1.194875e-10 5.565045e-08 1.113009e-06
# HiG_cardiac.b              2.440438e-14 3.904701e-13 1.005866e-08 2.011732e-07
# HiG_endothelial.a          1.772570e-45 2.836111e-44 1.058956e-09 2.117911e-08
# HiG_endothelial.c          5.575150e-23 8.920240e-22 2.559150e-13 5.118300e-12
# HiG_endothelial.d          2.819900e-15 4.511841e-14 5.149142e-13 1.029828e-11
# HiG_extraembryonicMesoderm 1.857953e-31 2.972724e-30 4.516141e-09 9.032282e-08
# HiG_mesodermProgenitors    2.354842e-09 3.767748e-08 1.261726e-08 2.523452e-07
# HiG_mixedMesoderm.a        2.151169e-12 3.441870e-11 7.221433e-10 1.444287e-08
# HiG_mixedMesoderm.b        3.040934e-18 4.865494e-17 1.504632e-07 3.009264e-06
# HiG_pharyngealMesoderm     9.862547e-20 1.578007e-18 5.633255e-08 1.126651e-06
# HiG_presomiticMesoderm.a   5.235415e-14 8.376665e-13 1.361927e-10 2.723855e-09
# HiG_presomiticMesoderm.b   2.625025e-23 4.200039e-22 1.594908e-07 3.189817e-06
# HiG_somiticMesoderm        2.568387e-15 4.109419e-14 4.237431e-07 8.474863e-06
# HiG_endothelial.b          8.761090e-15 1.401774e-13 1.339490e-10 2.678980e-09
# HiG_cardiac.a              4.693108e-26 7.508972e-25 5.394127e-12 1.078825e-10

(df[which(df$node_p.adj < 0.05), ])
#                                  edge_p   edge_p.adj       node_p   node_p.adj
# HiG_blood                  7.467968e-12 1.194875e-10 5.565045e-08 1.113009e-06
# HiG_cardiac.b              2.440438e-14 3.904701e-13 1.005866e-08 2.011732e-07
# HiG_cardiac.c              6.185585e-02 9.896936e-01 2.293858e-11 4.587716e-10
# HiG_endothelial.a          1.772570e-45 2.836111e-44 1.058956e-09 2.117911e-08
# HiG_endothelial.c          5.575150e-23 8.920240e-22 2.559150e-13 5.118300e-12
# HiG_endothelial.d          2.819900e-15 4.511841e-14 5.149142e-13 1.029828e-11
# HiG_extraembryonicMesoderm 1.857953e-31 2.972724e-30 4.516141e-09 9.032282e-08
# HiG_mesodermProgenitors    2.354842e-09 3.767748e-08 1.261726e-08 2.523452e-07
# HiG_mixedMesoderm.a        2.151169e-12 3.441870e-11 7.221433e-10 1.444287e-08
# HiG_mixedMesoderm.b        3.040934e-18 4.865494e-17 1.504632e-07 3.009264e-06
# HiG_pharyngealMesoderm     9.862547e-20 1.578007e-18 5.633255e-08 1.126651e-06
# HiG_presomiticMesoderm.a   5.235415e-14 8.376665e-13 1.361927e-10 2.723855e-09
# HiG_presomiticMesoderm.b   2.625025e-23 4.200039e-22 1.594908e-07 3.189817e-06
# HiG_somiticMesoderm        2.568387e-15 4.109419e-14 4.237431e-07 8.474863e-06
# HiG_endothelial.b          8.761090e-15 1.401774e-13 1.339490e-10 2.678980e-09
# HiG_cardiac.a              4.693108e-26 7.508972e-25 5.394127e-12 1.078825e-10
# CTS_endothelial.b                    NA           NA 5.658740e-07 1.131748e-05
# CTS_cardiac.a                        NA           NA 1.005473e-06 2.010946e-05




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

observed_auc_list <- list()
for (j in names(graph_list)) {
  observed_auc_list[[j]] <- Area_Under_Curve(
    subset(robustness.dt, signature == j & type == "Targeted vertex attack")$removed.pct,
    subset(robustness.dt, signature == j & type == "Targeted vertex attack")$comp.pct
  )
}
(observed_auc_list %>% unlist())
#                  HiG_blood              HiG_cardiac.b
#                  0.3579397                  0.3573040
#              HiG_cardiac.c          HiG_endothelial.a
#                  0.3402745                  0.3627558
#          HiG_endothelial.c          HiG_endothelial.d
#                  0.3388626                  0.3305859
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                  0.3463768                  0.3468226
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                  0.3315818                  0.3614308
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                  0.3540829                  0.3291077
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                  0.3717722                  0.3577911
#          HiG_endothelial.b              HiG_cardiac.a
#                  0.3448417                  0.3364005
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a
#                  0.5000000                  0.3461538
#          CTS_endothelial.b              CTS_cardiac.a
#                  0.2333333                  0.2162162

df_AUC <- data.frame(
  auc = observed_auc_list %>% unlist(),
  signature = names(observed_auc_list),
  PPI_cat = lapply(names(observed_auc_list), function(x) unlist(strsplit(x, split = "_"))[1]) %>% unlist()
)
df_AUC$PPI_cat <- factor(df_AUC$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

g3 <- ggplot(df_AUC, aes(x = PPI_cat, y = auc, fill = PPI_cat)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) + # Dodge the boxes for each type per PPI_cat
  scale_fill_manual(values = PPI_color_palette) +
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


(sum(vn[grep("^HiG_", names(vn))])) # 5893
(sum(vn[grep("^HiGCTS_", names(vn))])) # 25
(sum(vn[grep("^CTS_", names(vn))])) # 64

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
#                  HiG_blood              HiG_cardiac.b
#                  0.7124768                  0.7456660
#              HiG_cardiac.c          HiG_endothelial.a
#                  0.7214780                  0.7638871
#          HiG_endothelial.c          HiG_endothelial.d
#                  0.7592952                  0.7414829
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                  0.7482729                  0.7157201
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                  0.7192493                  0.7490758
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                  0.7466733                  0.7358549
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                  0.7490958                  0.7423322
#          HiG_endothelial.b              HiG_cardiac.a
#                  0.7418726                  0.7497369
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a
#                  0.5000000                  0.5500000
#          CTS_endothelial.b              CTS_cardiac.a
#                  0.4600000                  0.5661765

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
#                  HiG_blood              HiG_cardiac.b
#                0.053814315                0.062127011
#              HiG_cardiac.c          HiG_endothelial.a
#                0.104277596                0.055827720
#          HiG_endothelial.c          HiG_endothelial.d
#                0.056035722                0.024415885
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                0.068484363                0.031505378
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                0.037001219                0.039702684
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                0.036529114                0.025196840
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                0.046664297                0.027767327
#          HiG_endothelial.b              HiG_cardiac.a
#                0.041299671                0.043511523
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a
#                0.003468301                0.003550597
#          CTS_endothelial.b              CTS_cardiac.a
#                0.028693660                0.040382680

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
