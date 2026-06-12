library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library("gridExtra")
library(ggrepel)
library(ggpubr) # resuired by stat_compare_means()
library(igraph)

########## BEGINNING OF USER INPUT ##########

wd <- "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/GSE87038_IID/"
setwd(paste0(wd, "results/"))
inputdir <- file.path(wd, "results/")

db <- "GSE87038"

celltype_specific_weight_version <- "10"
source(paste0("https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v", celltype_specific_weight_version, ".R"))

cardiac_clusters <- c("2", "3", "4", "8", "11", "12", "14", "16", "17", "18", "19")
CT_id <- c("7", "8", "11", "13", "15", "16", "16.1") # critical transition clusters

redo <- FALSE # compute edge-betweenness communities. Only needed once!

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

########## END OF USER INPUT ##########


graph_list <- readRDS(file = file.path(inputdir, paste0(db, "_IID_graph_perState_notsimplified.rds")))
graph_list <- graph_list[names(graph_list) != "CTS_13"]
(N0 <- sapply(graph_list, vcount))
(N0)
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6
#         273         393         376         270         281         465
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#         316         380         309         333         475         420
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#         470         291         400         359         477         367
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#         292          10          17          28          11          28
#   HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#          12           9          26          46          60          34
#    CTS_16.1       CTS_8
#          73          51
########## remove name-duplicated Vertex due to inconsistence in STRING.db ###########

graph_list <- lapply(graph_list, simplify, edge.attr.comb = "max") # !!!!!!!!!!!!!!!!!!!
N2 <- sapply(graph_list, vcount)
(all(N0 == N2)) # [1] TRUE


(names(graph_list))
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_13"
# [26] "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"
# [31] "CTS_16.1"    "CTS_8"
edge_counts <- sapply(graph_list, ecount)
(edge_counts)
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6
#        2254        3987        2959        2249        2056        4706
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#        2758        3608        2430        2891        4503        3948
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#        3399        2523        4158        3660        5008        3490
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#        2279           0           2          10           6           6
#   HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#           3           1          10          18          50          28
#    CTS_16.1       CTS_8
#          57          65

###############################################################################
# 1) 'inspecting the interquartile range of betweenness centrality across cell states highlights a group of genes with high state-to-state variability'
# -> To quantify the topological importance of genes in GRN, the betweenness centrality of each gene is calculated, using the built-in NetworkX centrality method
# NetworkX → R igraph Equivalents:
# betweenness_centrality() -> betweenness()
###############################################################################
(is_weighted(graph_list[[1]])) # [1] TRUE
tmp_list <- graph_list # !!!!
for (i in seq_along(tmp_list)) {
  g <- tmp_list[[i]]
  E(g)$weight <- 1 / E(g)$weight
  tmp_list[[i]] <- g
}
BC <- sapply(tmp_list, betweenness, directed = FALSE, normalized = TRUE)

# 1) calculate the IQR per gene across all clusters
calculate_gene_iqr <- function(BC) {
  # Convert list to long-format data frame
  gene_data <- stack(BC)
  names(gene_data) <- c("value", "sample")

  # Add gene names by extracting from the original list
  gene_names <- unlist(lapply(BC, names))
  gene_data$gene <- gene_names

  # Calculate IQR for each gene
  result <- tapply(gene_data$value, gene_data$gene, IQR, na.rm = TRUE)

  return(result)
}

## all clusters
IQR <- calculate_gene_iqr(BC)
(range(IQR)) # [1] 0.00000000 0.06187858
(table(IQR > 0))
# FALSE  TRUE
#   746   894

# Sort by IQR value in descending order (like your plot)
plot_data <- data.frame(
  cell_type = names(IQR),
  iqr_value = as.numeric(IQR)
) %>%
  arrange(desc(iqr_value))

# Create the bar plot
n <- 15
g1 <- ggplot(plot_data[1:n, ], aes(x = reorder(cell_type, -iqr_value), y = iqr_value)) +
  geom_col(fill = "#7FCDCD", color = "white", width = 0.8) + # Light teal color
  labs(
    title = "Betweenness Centrality\ninterquartile range across cell clusters",
    x = paste("Top", n, "genes"),
    y = "IQR Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12, lineheight = 1.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "grey60"),
    axis.line.x = element_line(color = "grey60")
  ) +
  scale_y_continuous(
    breaks = seq(0, max(plot_data$iqr_value) * 1.1, by = 0.01),
    limits = c(0, max(plot_data$iqr_value) * 1.05),
    expand = c(0, 0)
  )

## cardiac-lineage clusters
cardiac <- setdiff(names(BC), paste0("HiG_", cardiac_clusters))
IQR <- calculate_gene_iqr(BC[cardiac])
(range(IQR)) # [1] 0.0000000 0.0646313
(table(IQR > 0))
# FALSE  TRUE
#   669   576

# Sort by IQR value in descending order (like your plot)
plot_data <- data.frame(
  cell_type = names(IQR),
  iqr_value = as.numeric(IQR)
) %>%
  arrange(desc(iqr_value))

# Create the bar plot
n <- 15
g <- ggplot(plot_data[1:n, ], aes(x = reorder(cell_type, -iqr_value), y = iqr_value)) +
  geom_col(fill = "#7FCDCD", color = "white", width = 0.8) + # Light teal color
  labs(
    title = "Betweenness Centrality\ninterquartile range across cardiac-lineage clusters",
    x = paste("Top", n, "genes"),
    y = "IQR Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12, lineheight = 1.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "grey60"),
    axis.line.x = element_line(color = "grey60")
  ) +
  scale_y_continuous(
    breaks = seq(0, max(plot_data$iqr_value) * 1.1, by = 0.01),
    limits = c(0, max(plot_data$iqr_value) * 1.05),
    expand = c(0, 0)
  )
pdf(file = "IQR_crossCluster_BetweennessCentrality.pdf")
grid.arrange(g1, g, ncol = 1)
dev.off()


###############################################################################
# 2) showing ‘the number of communities in the GRN graph is minimized in the intermediate, transitory state.’
# -> To estimate the communities, we employ NetworkX’s built-in functions girvan_newman() and greedy_modularity_communities() using default parameters
# -> the number of communities in the GRN graph defined as the highly connected sets of nodes
#
# NetworkX → R igraph Equivalents:
# girvan_newman() -> cluster_edge_betweenness()
# greedy_modularity_communities() -> cluster_fast_greedy()  Fast greedy modularity
# -> cluster_louvain()   Louvain method (also greedy modularity)
###############################################################################
# Community structure detection based on edge betweenness, this takes a while
## DO NOT repeat ########
if (redo) {
  EB <- sapply(tmp_list, cluster_edge_betweenness, directed = FALSE) # Run on midway3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  class(EB[[1]]) # "communities"
  saveRDS(EB, "cluster_edge_betweenness_list.rds")
}

EB <- readRDS(file = "cluster_edge_betweenness_list.rds")
# HiG line with other categories as dots at corresponding positions
nEB <- lengths(EB)


## fast greedy calculation gives different results, not USED
# EB = sapply(graph_list, cluster_fast_greedy )
# class(EB[[1]]) #"communities"
# # saveRDS(EB, 'cluster_fast_greedy_list.rds')  # NOT used !!
# plots = plot_nEB_ggplot(lengths(EB), PPI_color_palette, method='wilcox.test')
# grid.arrange(plots$line_plot, plots$boxplot, ncol = 2)

EB <- readRDS(file = "cluster_edge_betweenness_list.rds")
# HiG line with other categories as dots at corresponding positions
nEB <- lengths(EB)
plots <- plot_nEB_ggplot(nEB, PPI_color_palette, method = "wilcox.test")
grid.arrange(plots$line_plot, plots$boxplot, ncol = 2)
dev.copy2pdf(file = "community_number.pdf", width = 10)


############################################
# edge weigh distribution analysis across PPI categories (HiG, HiGCTS, CTS)
# showing no categry difference
############################################

# Required libraries
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(igraph)
# library(ggpubr)

# ==============================================================================
# MAIN ANALYSIS EXECUTION
# ==============================================================================
extract_edge_weights_by_category <- function(graph_list, PPI_color_palette, CT_id) {
  dfs <- lapply(names(graph_list), function(graph_name) {
    g <- graph_list[[graph_name]]

    # ---- skip graphs with no edges ----
    if (ecount(g) == 0) {
      return(NULL)
    }

    w <- E(g)$weight
    # ---- skip graphs with missing/empty weights ----
    if (is.null(w) || length(w) == 0) {
      return(NULL)
    }

    # infer PPI category from name prefix (HiG, HiGCTS, CTS)
    PPI_cat <- sub("_.*$", "", graph_name)

    # cluster id is whatever comes after the first underscore
    cluster_ID <- sub("^[^_]+_", "", graph_name)

    # label CT vs stable (customize if you have more labels)
    cluster_cat <- ifelse(cluster_ID %in% CT_id, "transitory", "stable")

    data.frame(
      sample = graph_name,
      PPI_cat = PPI_cat,
      edge_weight = as.numeric(w),
      num_edges = ecount(g),
      cluster_ID = cluster_ID,
      cluster_cat = cluster_cat,
      stringsAsFactors = FALSE
    )
  })

  edge_data <- dplyr::bind_rows(dfs)

  # keep factor order consistent with your palette
  edge_data$PPI_cat <- factor(edge_data$PPI_cat, levels = names(PPI_color_palette))

  edge_data
}

# Extract edge weight data by PPI category
edge_data <- extract_edge_weights_by_category(graph_list, PPI_color_palette, as.character(CT_id))
(head(edge_data, 3))
#   sample PPI_cat edge_weight num_edges cluster_ID cluster_cat
# 1  HiG_1     HiG        0.01      2254          1      stable
# 2  HiG_1     HiG        0.01      2254          1      stable
# 3  HiG_1     HiG        0.01      2254          1      stable

# Create plots for PPI category analysis
category_plots <- plot_edge_weight_distributions(edge_data, PPI_color_palette)

# Display all plots
("=== EDGE WEIGHT DISTRIBUTIONS BY PPI CATEGORY ===")
(category_plots$density_plot)
(category_plots$boxplot)

("=== SUMMARY STATISTICS ===")
edge_data %>%
  group_by(PPI_cat) %>%
  summarise(
    mean_weight = mean(edge_weight),
    median_weight = median(edge_weight),
    total_edges = sum(num_edges)
  )
# # A tibble: 3 × 4
#   PPI_cat mean_weight median_weight total_edges
#   <fct>         <dbl>         <dbl>       <int>
# 1 CTS           0.367         0.292      146303
# 2 HiGCTS        0.382         0.322        4672
# 3 HiG           0.366         0.295  1265162595
grid.arrange(category_plots$density_plot, category_plots$boxplot, ncol = 2)
dev.copy2pdf(file = "edge_weight.pdf")
