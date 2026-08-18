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

source('E:/Git_Holly/TIPS-main/celltype_specific_weight_v10.R')
 
setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9')
PPI_color_palette = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
 
 
graph_list <- readRDS( file= 'GSE175634_STRING_graph_perState_notsimplified.rds')  
(N0 = sapply(graph_list, vcount))
          # HiG_0           HiG_1           HiG_2          HiG_CP 
            # 325            1156             246             132 
          # HiG_4           HiG_5           HiG_6           HiG_7 
            # 105             339             108             108 
     # HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
            # 118            1468              70              42 
         # HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP 
             # 22              20               7              29 
    # HiGCTS_CP.1      CTS_muscle    CTS_endoderm          CTS_CP 
             # 27              61              64              72 
       # CTS_CP.1 
             # 77 
########## remove name-duplicated Vertex due to inconsistence in STRING.db ###########
## refer to 11.1.0_correct_vertex_duplication.R 
correct_n_edges = readRDS('correct_n_edges_HiG_STRING2.14.0.rds')
for(g_name in unique(correct_n_edges$graph_id)){
	vertices_to_remove = subset(correct_n_edges, graph_id == g_name)$vetex_index_to_remove
	if(any(is.na(vertices_to_remove))) vertices_to_remove = vertices_to_remove[!is.na(vertices_to_remove)]
	graph_list[[g_name]] = delete_vertices(graph_list[[g_name]], vertices_to_remove)
}
N = sapply(graph_list, vcount)
(N0-N)[which(N0-N>0)]
# HiG_1 HiG_2 HiG_4 HiG_6 HiG_7 HiG_9 
    # 1     4     1     2     4     1 

graph_list <- lapply(graph_list, simplify, edge.attr.comb = 'max') #!!!!!!!!!!!!!!!!!!!
N2 = sapply(graph_list, vcount)
all(N==N2)   # [1] TRUE 


names(graph_list)
 # [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"           "HiG_5"          
 # [7] "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"   
# [13] "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"     
# [19] "CTS_endoderm"    "CTS_CP"          "CTS_CP.1"    
edge_counts <- sapply(graph_list, ecount)
edge_counts
          # HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4           HiG_5 
          # 10066           31169            7129            1306            2269            4638 
          # HiG_6           HiG_7      HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
            # 894            1587            1989           42812            1213             141 
         # HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP     HiGCTS_CP.1      CTS_muscle 
             # 70              68              11              38              41             316 
   # CTS_endoderm          CTS_CP        CTS_CP.1 
            # 390             128             159 

###############################################################################
# 1) 'inspecting the interquartile range of betweenness centrality across cell states highlights a group of genes with high state-to-state variability'
# -> To quantify the topological importance of genes in GRN, the betweenness centrality of each gene is calculated, using the built-in NetworkX centrality method
# NetworkX → R igraph Equivalents:
# betweenness_centrality() -> betweenness()
###############################################################################
is_weighted(graph_list[[1]])  # [1] TRUE
tmp_list = graph_list  #!!!!
for (i in 1:length(tmp_list)) {
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
range(IQR)  # [1] 0.00000000 0.1946385
table(IQR>0)
# FALSE  TRUE 
  # 732  1400 
  
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
cardiac = setdiff(names(BC), c("HiG_0" , "HiG_2" , "HiG_7", "HiG_10" , "HiG_11" ))
IQR =calculate_gene_iqr(BC[cardiac])
range(IQR)  # [1] 0.00000000 0.1818138
table(IQR>0)
# FALSE  TRUE 
  # 631  1172 
  
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
redo= FALSE
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
dev.copy2pdf(file='community_number.pdf', width=10, height=5)
 
 
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
unstable_cluster_ID = c('muscle','endoderm','CP')

# ==============================================================================
# MAIN ANALYSIS EXECUTION
# ==============================================================================
# Extract edge weight data by PPI category
edge_data <- extract_edge_weights_by_category(graph_list, PPI_color_palette, unstable_cluster_ID)
head(edge_data, 3)
  # sample PPI_cat edge_weight num_edges cluster_ID cluster_cat
# 1  HiG_0     HiG       0.309     10066          0      stable
# 2  HiG_0     HiG       0.214     10066          0      stable
# 3  HiG_0     HiG       0.504     10066          0      stable

# Create plots for PPI category analysis
category_plots <- plot_edge_weight_distributions(edge_data, PPI_color_palette)

# Display all plots
print("=== EDGE WEIGHT DISTRIBUTIONS BY PPI CATEGORY ===")
print(category_plots$density_plot)
print(category_plots$boxplot)


category_plots
# $summary_stats
# # A tibble: 3 × 4
  # PPI_cat mean_weight median_weight total_edges
  # <fct>         <dbl>         <dbl>       <int>
# 1 CTS           0.437         0.348         993
# 2 HiGCTS        0.435         0.344         158
# 3 HiG           0.386         0.302      105283
grid.arrange(category_plots$density_plot, category_plots$boxplot, ncol = 2)
dev.copy2pdf(file='edge_weight.pdf', height=3.5, width=7)

