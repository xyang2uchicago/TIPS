# this code verify if it is true that 
# 1) 'inspecting the interquartile range of betweenness centrality across cell states highlights a group of genes with high state-to-state variability'
# -> To quantify the topological importance of genes in GRN, the betweenness centrality of each gene is calculated, using the built-in NetworkX centrality method
# 2) showing ‘the number of communities in the GRN graph is minimized in the intermediate, transitory state.’
# -> To quantify the topological importance of genes in GRN, the betweenness centrality of each gene is calculated, using the built-in NetworkX centrality method
# -> To estimate the communities, we employ NetworkX’s built-in functions girvan_newman() and greedy_modularity_communities() using default parameters 
# -> the number of communities in the GRN graph defined as the highly connected sets of nodes
#  suggesting ‘Interactions between genes are more compartmentalized in the stable states, with GRNs exhibiting a modular structure. Conversely, genes are more highly connected in the intermediate, unstable state.’
# 3) showing ‘edge weight distribution is narrower for pseudotime points corresponding to the stable states—indicting that fewer gene-gene connections are significant—while being broader for intermediate pseudotime points corresponding to the tipping point—indicating a larger number of connections between genes.’

# NetworkX → R igraph Equivalents:
# girvan_newman() -> cluster_edge_betweenness()
# greedy_modularity_communities() -> cluster_fast_greedy()  Fast greedy modularity
				# -> cluster_louvain()   Louvain method (also greedy modularity)

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

wd = "/Users/felixyu/Documents/IbarraSoria2018/"
setwd(paste0(wd, "results/"))

celltype_specific_weight_version <- '10'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

cardiac_clusters <- c(
    "cardiac.a", "cardiac.b", "cardiac.c",
    "extraembryonicMesoderm", "mixedMesoderm.a",
    "mixedMesoderm.b", "pharyngealMesoderm",
    "presomiticMesoderm.a", "presomiticMesoderm.b",
    "somiticMesoderm"
)

CT_id <- c("endothelial.b", "cardiac.a")  # critical transition clusters

redo = TRUE # compute edge-betweenness communities. Only needed once!

PPI_color_palette = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

########## END OF USER INPUT ##########
 
 
graph_list <- readRDS( file= paste0(db, '_STRING_graph_perState_notsimplified.rds'))  
(N0 = sapply(graph_list, vcount))
(N0)
#                  HiG_blood              HiG_cardiac.b 
#                        420                        409 
#              HiG_cardiac.c          HiG_endothelial.a 
#                        444                        455 
#          HiG_endothelial.c          HiG_endothelial.d 
#                        480                        426 
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors 
#                        367                        356 
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b 
#                        329                        369 
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a 
#                        371                        338 
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm 
#                        412                        327 
#          HiG_endothelial.b              HiG_cardiac.a 
#                        418                        407 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a 
#                         16                         13 
#          CTS_endothelial.b              CTS_cardiac.a 
#                         33                         37 

(names(graph_list))
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"          "HiG_endothelial.c"         
#  [4] "HiG_endothelial.d"          "HiG_blood"                  "HiG_mesodermProgenitors"   
#  [7] "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"   "HiG_somiticMesoderm"       
# [10] "HiG_mixedMesoderm.a"        "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"              "HiG_endothelial.b"         
# [16] "HiG_cardiac.a"              "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"          
# [19] "CTS_endothelial.b"          "CTS_cardiac.a"   
edge_counts <- sapply(graph_list, ecount)
(edge_counts)
#                  HiG_blood              HiG_cardiac.b 
#                      10986                       7164 
#              HiG_cardiac.c          HiG_endothelial.a 
#                       8853                       9643 
#          HiG_endothelial.c          HiG_endothelial.d 
#                       8915                       7825 
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors 
#                       5666                       4991 
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b 
#                       4527                       5930 
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a 
#                       6110                       4881 
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm 
#                       7797                       4687 
#          HiG_endothelial.b              HiG_cardiac.a 
#                       6855                       6513 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a 
#                         20                         28 
#          CTS_endothelial.b              CTS_cardiac.a 
#                         82                         79 

###############################################################################
# 1) 'inspecting the interquartile range of betweenness centrality across cell states highlights a group of genes with high state-to-state variability'
# -> To quantify the topological importance of genes in GRN, the betweenness centrality of each gene is calculated, using the built-in NetworkX centrality method
# NetworkX → R igraph Equivalents:
# betweenness_centrality() -> betweenness()
###############################################################################
(is_weighted(graph_list[[1]]))  # [1] TRUE
tmp_list = graph_list  #!!!!
for (i in seq_along(tmp_list)) {
  g = tmp_list[[i]]
  E(g)$weight = 1/E(g)$weight
  tmp_list[[i]] = g
}
BC = sapply(tmp_list, betweenness, directed=FALSE, normalized = TRUE)

#1) calculate the IQR per gene across all clusters
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
IQR = calculate_gene_iqr(BC)
(range(IQR))  # [1] 0.0000000 0.1818125
(table(IQR>0))
# FALSE  TRUE 
#   556   942
  
# Sort by IQR value in descending order (like your plot)
plot_data <- data.frame(
      cell_type = names(IQR),
      iqr_value = as.numeric(IQR)
    ) %>%
    arrange(desc(iqr_value))
  
  # Create the bar plot
n = 15
g1 = ggplot(plot_data[1:n,], aes(x = reorder(cell_type, -iqr_value), y = iqr_value)) +
    geom_col(fill = "#7FCDCD", color = "white", width = 0.8) +  # Light teal color
    labs(
      title = "Betweenness Centrality\ninterquartile range across cell clusters",
      x = paste("Top",n, "genes"),
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
cardiac = setdiff(names(BC), paste0("HiG_", cardiac_clusters))
IQR =calculate_gene_iqr(BC[cardiac])
(range(IQR))  # [1] 0.0000000 0.0950376
(table(IQR>0))
# FALSE  TRUE 
#   550   641 
  
# Sort by IQR value in descending order (like your plot)
plot_data <- data.frame(
      cell_type = names(IQR),
      iqr_value = as.numeric(IQR)
    ) %>%
    arrange(desc(iqr_value))
  
  # Create the bar plot
n = 15
g = ggplot(plot_data[1:n,], aes(x = reorder(cell_type, -iqr_value), y = iqr_value)) +
    geom_col(fill = "#7FCDCD", color = "white", width = 0.8) +  # Light teal color
    labs(
      title = "Betweenness Centrality\ninterquartile range across cardiac-lineage clusters",
      x = paste("Top",n, "genes"),
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
pdf(file= 'IQR_crossCluster_BetweennessCentrality.pdf' ) 
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
if(redo){
	EB = sapply(tmp_list, cluster_edge_betweenness, directed=FALSE)
	class(EB[[1]]) #"communities"
	saveRDS(EB, 'cluster_edge_betweenness_list.rds')
	
}

EB = readRDS(file='cluster_edge_betweenness_list.rds')
# HiG line with other categories as dots at corresponding positions
nEB = lengths(EB)


## fast greedy calculation gives different results, not USED
# EB = sapply(graph_list, cluster_fast_greedy )
# class(EB[[1]]) #"communities"
# # saveRDS(EB, 'cluster_fast_greedy_list.rds')  # NOT used !!
# plots = plot_nEB_ggplot(lengths(EB), PPI_color_palette, method='wilcox.test')
# grid.arrange(plots$line_plot, plots$boxplot, ncol = 2)

EB = readRDS(file='cluster_edge_betweenness_list.rds')
# HiG line with other categories as dots at corresponding positions
nEB = lengths(EB)
plots = plot_nEB_ggplot(nEB, PPI_color_palette, method='wilcox.test')
grid.arrange(plots$line_plot, plots$boxplot, ncol = 2)
dev.copy2pdf(file='community_number.pdf', width=10)
 
 
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
# Extract edge weight data by PPI category
edge_data <- extract_edge_weights_by_category(graph_list, PPI_color_palette, CT_id)
(head(edge_data, 3))
#      sample PPI_cat edge_weight num_edges cluster_ID cluster_cat
# 1 HiG_blood     HiG 0.016409879     10986      blood      stable
# 2 HiG_blood     HiG 0.024549195     10986      blood      stable
# 3 HiG_blood     HiG 0.006696524     10986      blood      stable

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
# 1 CTS          0.0746        0.0542       12965
# 2 HiGCTS       0.110         0.101         1184
# 3 HiG          0.0229        0.0135   830085899
grid.arrange(category_plots$density_plot, category_plots$boxplot, ncol = 2)
dev.copy2pdf(file='edge_weight.pdf')

