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

source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "GSE87038", "GSE87038_IID/")
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- file.path(wd, "data/")

celltype_specific_weight_version <- "10"
source(paste0("https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v", celltype_specific_weight_version, ".R"))

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_palette <- c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)

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
# "weight"         "norm_PPI_score" "corexp_sign"    "coexp_focal"
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

pdf(file = paste0("vertex_attack_", db, ".pdf"))
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
file <- file.path(inputdir, paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds"))
graph_list <- readRDS(file)

(names(graph_list))
#  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"
#  [7] "HiG_9"     "HiG_10"    "HiG_12"    "HiG_14"    "HiG_17"    "HiG_18"
# [13] "HiG_19"    "HiG_7"     "HiG_11"    "HiG_15"    "HiG_16"    "HiG_13"
# [19] "HiG_8"     "HiGCTS_15" "CTS_7"     "CTS_11"    "CTS_15"    "CTS_16"
# [25] "CTS_16.1"  "CTS_8"


attack.vertex.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "vertex", "btwn.cent"), idcol = names(graph_list))
colnames(attack.vertex.btwn)[1] <- "signature"
(head(attack.vertex.btwn, 3))
#    signature                   type   measure comp.size  comp.pct removed.pct
#       <char>                 <char>    <char>     <num>     <num>       <num>
# 1:     HiG_1 Targeted vertex attack btwn.cent       261 1.0000000 0.000000000
# 2:     HiG_1 Targeted vertex attack btwn.cent       259 0.9923372 0.003663004
# 3:     HiG_1 Targeted vertex attack btwn.cent       258 0.9885057 0.007326007

(dim(attack.vertex.btwn)) #  7291    6
saveRDS(attack.vertex.btwn, file = "attack.vertex.btwn.rds")

# calculate the 'egde' btwn for each G first
# igraph::edge_betweenness() uses distance graph weights, but E(g) uses connection weights, thus we invert it.
for (j in seq_along(graph_list)) E(graph_list[[j]])$btwn <- edge_betweenness(graph_list[[j]], weights = 1 / E(graph_list[[j]])$weight)
## then do the edge-attack analysis
attack.edge.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, "edge", "btwn.cent"), idcol = names(graph_list)) # !!!!!!!!!!!!!
colnames(attack.edge.btwn)[1] <- "signature"
(dim(attack.edge.btwn)) #  63130      6
(head(attack.edge.btwn, 3))
#    signature                 type   measure comp.size comp.pct  removed.pct
#       <char>               <char>    <char>     <num>    <num>        <num>
# 1:     HiG_1 Targeted edge attack btwn.cent       261        1 0.0000000000
# 2:     HiG_1 Targeted edge attack btwn.cent       261        1 0.0004436557
# 3:     HiG_1 Targeted edge attack btwn.cent       261        1 0.0008873114

(table(attack.edge.btwn$signature))
#    CTS_11    CTS_15    CTS_16  CTS_16.1     CTS_7     CTS_8     HiG_1    HiG_10
#        19        51        29        58        11        66      2255      3609
#    HiG_11    HiG_12    HiG_13    HiG_14    HiG_15    HiG_16    HiG_17    HiG_18
#      4159      2431      3491      2892      3661      5009      4504      3949
#    HiG_19     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_7     HiG_8
#      3400      3988      2960      2250      2057      4707      2524      2280
#     HiG_9 HiGCTS_15
#      2759        11

saveRDS(attack.edge.btwn, file = "attack.edge.btwn.rds")

# In a random failure analysis, you choose a vertex/edge at random, remove it, and calculate the maximum
#  component size until all elements have been removed.
####################################################
# refer to 11.2_CTS_cardiac_network_robustness_Miudway3.simulation.R
## run on Midway3, the following is long. DO NOT REPEAT !!!)


(sapply(graph_list, vcount))
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_9    HiG_10
#       273       393       376       270       281       465       316       380
#    HiG_12    HiG_14    HiG_17    HiG_18    HiG_19     HiG_7    HiG_11    HiG_15
#       309       333       475       420       470       291       400       359
#    HiG_16    HiG_13     HiG_8 HiGCTS_15     CTS_7    CTS_11    CTS_15    CTS_16
#       477       367       292        28        26        46        60        34
#  CTS_16.1     CTS_8
#        73        51

attack.edge.btwn <- readRDS(file = "attack.edge.btwn.rds")
attack.vertex.btwn <- readRDS(file = "attack.vertex.btwn.rds")
failure.vertex <- readRDS(paste0("failure.vertex_100_simplified_", s, "weighted.rds"))
(table(attack.vertex.btwn$signature, attack.vertex.btwn$type))
#           Targeted vertex attack
#   CTS_11                        47
#   CTS_15                        61
#   CTS_16                        35
#   CTS_16.1                      74
#   CTS_7                         27
#   CTS_8                         52
#   HiG_1                        274
#   HiG_10                       381
#   HiG_11                       401
#   HiG_12                       310
#   HiG_13                       368
#   HiG_14                       334
#   HiG_15                       360
#   HiG_16                       478
#   HiG_17                       476
#   HiG_18                       421
#   HiG_19                       471
#   HiG_2                        394
#   HiG_3                        377
#   HiG_4                        271
#   HiG_5                        282
#   HiG_6                        466
#   HiG_7                        292
#   HiG_8                        293
#   HiG_9                        317
#   HiGCTS_15                     29

(subset(attack.vertex.btwn, signature == "HiGCTS_8"))
# empty


failure.edge <- readRDS(paste0("failure.edge_100_simplified_", s, "weighted.rds"))
failure.dt <- rbind(failure.edge, failure.vertex)
(head(failure.dt, 3))
#     HiG_1                type measure comp.size comp.pct  removed.pct
#    <char>              <char>  <char>     <num>    <num>        <num>
# 1:  HiG_1 Random edge removal  random       261        1 0.0000000000
# 2:  HiG_1 Random edge removal  random       261        1 0.0004436557
# 3:  HiG_1 Random edge removal  random       261        1 0.0008873114

colnames(failure.dt)[1] <- "signature"
(table(failure.dt$signature, failure.dt$type))
#             Random edge removal Random vertex removal
#   CTS_11                     19                    47
#   CTS_15                     51                    61
#   CTS_16                     29                    35
#   CTS_16.1                   58                    74
#   CTS_7                      11                    27
#   CTS_8                      66                    52
#   HiG_1                    2255                   274
#   HiG_10                   3609                   381
#   HiG_11                   4159                   401
#   HiG_12                   2431                   310
#   HiG_13                   3491                   368
#   HiG_14                   2892                   334
#   HiG_15                   3661                   360
#   HiG_16                   5009                   478
#   HiG_17                   4504                   476
#   HiG_18                   3949                   421
#   HiG_19                   3400                   471
#   HiG_2                    3988                   394
#   HiG_3                    2960                   377
#   HiG_4                    2250                   271
#   HiG_5                    2057                   282
#   HiG_6                    4707                   466
#   HiG_7                    2524                   292
#   HiG_8                    2280                   293
#   HiG_9                    2759                   317
#   HiGCTS_15                  11                    29

colnames(attack.vertex.btwn)[1] <- "signature"
(dim(failure.dt)) #  70421      6

robustness.dt <- rbind(failure.dt, attack.vertex.btwn, attack.edge.btwn)
(dim(robustness.dt)) #  140842      6
robustness.dt$PPI_cat <- lapply(robustness.dt$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
  unlist() %>%
  factor(., levels = c("CTS", "HiGCTS", "HiG"))
(head(robustness.dt, 3))
#    signature                type measure comp.size comp.pct  removed.pct
#       <char>              <char>  <char>     <num>    <num>        <num>
# 1:     HiG_1 Random edge removal  random       261        1 0.0000000000
# 2:     HiG_1 Random edge removal  random       261        1 0.0004436557
# 3:     HiG_1 Random edge removal  random       261        1 0.0008873114

robustness.dt$experiment <- ifelse(grepl("edge", robustness.dt$type), "edge", "vertex")
robustness.dt$measure <- factor(robustness.dt$measure, levels = c("random", "btwn.cent"))

(table(robustness.dt$type, robustness.dt$measure))
#                          random btwn.cent
#   Random edge removal     63130         0
#   Random vertex removal    7291         0
#   Targeted edge attack        0     63130
#   Targeted vertex attack      0      7291


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
# 1        edge     CTS    random c(W = 0..... 1.500078e-09
# 2        edge     CTS btwn.cent c(W = 0..... 3.817775e-08
# 3        edge  HiGCTS    random c(W = 0..... 4.143536e-02
# 4        edge  HiGCTS btwn.cent c(W = 0..... 6.582988e-01
# 5        edge     HiG    random           NA           NA
# 6        edge     HiG btwn.cent           NA           NA
# 7      vertex     CTS    random c(W = 0..... 4.010969e-12
# 8      vertex     CTS btwn.cent c(W = 0..... 9.018697e-25
# 9      vertex  HiGCTS    random c(W = 0..... 1.137946e-01
# 10     vertex  HiGCTS btwn.cent c(W = 0..... 4.002824e-09
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
# 1 CTS                1.05                2.34
# 2 HiGCTS             1                   2.28
# 3 HiG                0.954               1.37

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
# W = 2480, p-value = 0.0001696


wilcox.test(
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "CTS")],
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiG")]
)
# W = 906080, p-value = 0.0004074


wilcox.test(
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiGCTS")],
  targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat == "HiG")]
)
# W = 92587, p-value = 0.4379

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
#              edge_p   edge_p.adj       node_p   node_p.adj
# HiG_1  1.506519e-16 3.916948e-15 3.627120e-07 9.430512e-06
# HiG_2  4.577468e-31 1.190142e-29 7.902987e-08 2.054777e-06
# HiG_3  1.967117e-10 5.114504e-09 3.055115e-10 7.943299e-09
# HiG_4  1.798962e-09 4.677302e-08 8.502522e-09 2.210656e-07
# HiG_5  1.325743e-12 3.446932e-11 1.673486e-09 4.351064e-08
# HiG_6  1.121175e-17 2.915056e-16 1.842558e-12 4.790652e-11
# HiG_9  5.870875e-17 1.526428e-15 5.951356e-09 1.547353e-07
# HiG_10 2.719330e-14 7.070257e-13 6.899793e-09 1.793946e-07
# HiG_12 7.090312e-17 1.843481e-15 6.850439e-09 1.781114e-07
# HiG_14 1.234290e-14 3.209154e-13 2.157660e-10 5.609915e-09
# HiG_17 4.124328e-10 1.072325e-08 1.418392e-14 3.687819e-13
# HiG_18 7.328695e-06 1.905461e-04 1.271245e-11 3.305236e-10
# HiG_19 2.120857e-10 5.514228e-09 3.173928e-17 8.252212e-16
# HiG_7  1.034385e-06 2.689400e-05 4.018828e-07 1.044895e-05
# HiG_11 4.624692e-32 1.202420e-30 8.385677e-09 2.180276e-07
# HiG_15 4.501129e-22 1.170294e-20 1.861238e-10 4.839218e-09
# HiG_16 3.224559e-14 8.383853e-13 7.560742e-12 1.965793e-10
# HiG_13 3.748390e-09 9.745813e-08 7.947369e-10 2.066316e-08
# HiG_8  1.513182e-08 3.934273e-07 1.944746e-08 5.056339e-07


(df[which(df$node_p.adj < 0.05), ])
#                 edge_p   edge_p.adj       node_p   node_p.adj
# HiG_1     1.506519e-16 3.916948e-15 3.627120e-07 9.430512e-06
# HiG_2     4.577468e-31 1.190142e-29 7.902987e-08 2.054777e-06
# HiG_3     1.967117e-10 5.114504e-09 3.055115e-10 7.943299e-09
# HiG_4     1.798962e-09 4.677302e-08 8.502522e-09 2.210656e-07
# HiG_5     1.325743e-12 3.446932e-11 1.673486e-09 4.351064e-08
# HiG_6     1.121175e-17 2.915056e-16 1.842558e-12 4.790652e-11
# HiG_9     5.870875e-17 1.526428e-15 5.951356e-09 1.547353e-07
# HiG_10    2.719330e-14 7.070257e-13 6.899793e-09 1.793946e-07
# HiG_12    7.090312e-17 1.843481e-15 6.850439e-09 1.781114e-07
# HiG_14    1.234290e-14 3.209154e-13 2.157660e-10 5.609915e-09
# HiG_17    4.124328e-10 1.072325e-08 1.418392e-14 3.687819e-13
# HiG_18    7.328695e-06 1.905461e-04 1.271245e-11 3.305236e-10
# HiG_19    2.120857e-10 5.514228e-09 3.173928e-17 8.252212e-16
# HiG_7     1.034385e-06 2.689400e-05 4.018828e-07 1.044895e-05
# HiG_11    4.624692e-32 1.202420e-30 8.385677e-09 2.180276e-07
# HiG_15    4.501129e-22 1.170294e-20 1.861238e-10 4.839218e-09
# HiG_16    3.224559e-14 8.383853e-13 7.560742e-12 1.965793e-10
# HiG_13    3.748390e-09 9.745813e-08 7.947369e-10 2.066316e-08
# HiG_8     1.513182e-08 3.934273e-07 1.944746e-08 5.056339e-07
# HiGCTS_15 7.620045e-01 1.000000e+00 2.778606e-06 7.224375e-05
# CTS_11    9.179205e-01 1.000000e+00 6.037989e-07 1.569877e-05
# CTS_15    2.103000e-01 1.000000e+00 4.170469e-12 1.084322e-10
# CTS_16    1.817796e-01 1.000000e+00 4.060433e-05 1.055713e-03
# CTS_16.1  8.922635e-01 1.000000e+00 8.234858e-11 2.141063e-09
# CTS_8     1.623003e-01 1.000000e+00 1.847193e-04 4.802703e-03


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
#  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"
#  [7] "HiG_9"     "HiG_10"    "HiG_12"    "HiG_14"    "HiG_17"    "HiG_18"
# [13] "HiG_19"    "HiG_7"     "HiG_11"    "HiG_15"    "HiG_16"    "HiG_13"
# [19] "HiG_8"     "HiGCTS_15" "CTS_7"     "CTS_11"    "CTS_15"    "CTS_16"
# [25] "CTS_16.1"  "CTS_8"

observed_auc_list <- list()
for (j in names(graph_list)) {
  observed_auc_list[[j]] <- Area_Under_Curve(
    subset(robustness.dt, signature == j & type == "Targeted vertex attack")$removed.pct,
    subset(robustness.dt, signature == j & type == "Targeted vertex attack")$comp.pct
  )
}
(observed_auc_list %>% unlist())
#      HiG_1      HiG_2      HiG_3      HiG_4      HiG_5      HiG_6      HiG_9
# 0.35266585 0.37107491 0.34105109 0.33473137 0.32942413 0.34581559 0.34463905
#     HiG_10     HiG_12     HiG_14     HiG_17     HiG_18     HiG_19      HiG_7
# 0.35905104 0.33777763 0.33423423 0.33575694 0.34418487 0.31693064 0.35689387
#     HiG_11     HiG_15     HiG_16     HiG_13      HiG_8  HiGCTS_15      CTS_7
# 0.36370383 0.34561566 0.35259491 0.34835886 0.34097517 0.17410714 0.41153846
#     CTS_11     CTS_15     CTS_16   CTS_16.1      CTS_8
# 0.14285714 0.09318182 0.14619377 0.12512529 0.16873065

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
#    HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_9    HiG_10
#       273       393       376       270       281       465       316       380
#    HiG_12    HiG_14    HiG_17    HiG_18    HiG_19     HiG_7    HiG_11    HiG_15
#       309       333       475       420       470       291       400       359
#    HiG_16    HiG_13     HiG_8 HiGCTS_15     CTS_7    CTS_11    CTS_15    CTS_16
#       477       367       292        28        26        46        60        34
#  CTS_16.1     CTS_8
#        73        51


(sum(vn[grep("^HiG_", names(vn))])) # 6947
(sum(vn[grep("^HiGCTS_", names(vn))])) # 28
(sum(vn[grep("^CTS_", names(vn))])) # 290

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
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_9    HiG_10
# 0.7320600 0.7523153 0.7277549 0.7382122 0.7353146 0.7494847 0.7418282 0.7316825
#    HiG_12    HiG_14    HiG_17    HiG_18    HiG_19     HiG_7    HiG_11    HiG_15
# 0.7337153 0.7275919 0.7296007 0.7504656 0.7368885 0.7265802 0.7425071 0.7619156
#    HiG_16    HiG_13     HiG_8 HiGCTS_15     CTS_7    CTS_11    CTS_15    CTS_16
# 0.7413575 0.7255237 0.7133935 0.4500000 0.6900000 0.4206349 0.4930303 0.4233193
#  CTS_16.1     CTS_8
# 0.4798887 0.5643725

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
# 0.037039995 0.088795202 0.056675185 0.046265892 0.034635309 0.039960566
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
# 0.036374545 0.076035932 0.046981202 0.042265343 0.033066542 0.019491901
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
# 0.017443164 0.030291971 0.040172395 0.041246032 0.042065084 0.046655542
#       HiG_8   HiGCTS_15       CTS_7      CTS_11      CTS_15      CTS_16
# 0.044048547 0.009352490 0.006229994 0.005943201 0.013914544 0.007211264
#    CTS_16.1       CTS_8
# 0.009086724 0.014347385
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
