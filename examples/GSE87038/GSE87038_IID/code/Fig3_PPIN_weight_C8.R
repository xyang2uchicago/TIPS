library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(igraph)
library(rstatix)
library(scales)
library(pracma)
library(MLmetrics)
library(sm)

########## BEGINNING OF USER INPUT ##########

wd <- "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/GSE87038_IID/"
setwd(paste0(wd, "results/PPI_weight/"))
inputdir <- file.path(wd, "data/")
celltype_specific_weight_version <- "10"
source(paste0("https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v", celltype_specific_weight_version, ".R"))

db <- "GSE87038"

PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_palette <- c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)


s <- "combined" # specificity method

########## END OF USER INPUT ##########

###################################################
# Fig B) boxplot of normalized strength for all clusters per category
# original code: 11.3_CTS_cardiac_network_ANND_pagerank.R
# original pdf: normalized.node.strength_GSE87038_v2.pdf
################################################################
{
  CHD <- readRDS(file = file.path(inputdir, (paste0("CHD_Cilia_Genelist.rds"))))

  df <- readRDS(file = "IID_df_PAGERANK_strength_ANND.rewiring.P.rds")

  df <- df[grepl("_8$", df$signature), ]

  # Keep only desired PPI categories in the correct order
  df <- rbind(
    subset(df, PPI_cat == "CTS"),
    subset(df, PPI_cat == "HiGCTS"),
    subset(df, PPI_cat == "HiG")
  )

  df$signature <- as.character(df$signature) # <- important
  df$PPI_cat <- factor(df$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

  # lock x order explicitly (character vector)
  sig_order <- c("CTS_8", "HiGCTS_8", "HiG_8")
  df$signature <- factor(df$signature, levels = sig_order)


  # Enforce the factor order
  df$PPI_cat <- factor(df$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))
  df$signature <- factor(df$signature, levels = unique(df$signature))

  # Add CHD gene annotations
  df$PCGC_AllCurated <- toupper(df$gene) %in% toupper(unlist(CHD["Griffin2023_PCGC_AllCurated"]))

  # Identify top 5 genes per signature by normalized strength
  top_genes <- df %>%
    group_by(signature) %>%
    arrange(desc(normalized.strength)) %>%
    slice_head(n = 5) %>%
    ungroup()

  # Subset top CHD genes
  top_genes_CHD <- subset(top_genes, PCGC_AllCurated == TRUE)
  (dim(top_genes_CHD)) # Optional: check how many CHD genes were top 5

  # Optional: write out table of top 5 genes
  tb <- top_genes[, c("signature", "gene", "PPI_cat", "normalized.strength", "PCGC_AllCurated")]
  write.table(tb, file = "table_top5_strength_perPPI.tsv", sep = "\t", row.names = FALSE)

  # Plot
  boxplot_strength <- ggplot(df, aes(x = signature, y = normalized.strength, colour = PPI_cat)) +
    geom_boxplot(show.legend = TRUE, position = position_dodge2(preserve = "single")) +
    scale_color_manual(values = PPI_color_palette) +
    geom_text_repel(
      data = top_genes_CHD,
      aes(label = gene),
      size = 2,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = 20,
      show.legend = FALSE
    ) +
    theme(
      legend.position = "none",
      legend.justification = c(1, 1),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    stat_compare_means(
      comparisons = list(
        c("CTS_8", "HiGCTS_8"),
        c("CTS_8", "HiG_8"),
        c("HiGCTS_8", "HiG_8")
      ),
      method = "wilcox.test",
      label = "p.signif",
      p.adjust.method = "bonferroni",
      hide.ns = TRUE
    ) +
    scale_x_discrete(limits = unique(df$signature)) +
    labs(color = "PPI cat")

  vertex(boxplot_strength)
}

###################################################
# Fig F) ARC plot of vertex attack
# original code: 11.2_CTS_cardiac_network_robustness.R
# original pdf: attack_GSE87038.pdf
################################################################
{
  robustness.dt <- rbind(failure.dt, attack.vertex.btwn[, 1:6])
  failure.dt <- readRDS(file = paste0("failure.vertex_100_simplified_", s, "weighted.rds"))

  setnames(failure.dt, old = "HiG_1", new = "signature")

  # attack.edge.btwn = readRDS( file=paste0('attack.edge.btwn.rds'))
  attack.vertex.btwn <- readRDS(file = paste0("attack.vertex.btwn.rds"))

  failure.dt <- failure.dt[grepl("_8$", failure.dt$signature), ]
  attack.vertex.btwn <- attack.vertex.btwn[grepl("_8$", attack.vertex.btwn$signature), ]

  robustness.dt <- rbind(failure.dt, attack.vertex.btwn[, 1:6]) # , attack.edge.btwn[,1:6])
  (dim(robustness.dt)) #  345   9
  robustness.dt$PPI_cat <- lapply(robustness.dt$signature, function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>%
    factor(., levels = c("CTS", "HiGCTS", "HiG"))
  head(robustness.dt, 3)


  robustness.dt$measure %>% unique()
  # Levels: random btwn.cent
  robustness.dt <- subset(robustness.dt, measure != "degree")
  robustness.dt$experiment <- ifelse(grepl("edge", robustness.dt$type), "edge", "vertex")
  robustness.dt$measure <- factor(robustness.dt$measure, levels = c("random", "btwn.cent"))

  robustness.dt$type <- factor(robustness.dt$type,
    levels = c("Random edge removal", "Targeted edge attack", "Random vertex removal", "Targeted vertex attack")
  )
  robustness.dt$cluster <- lapply(robustness.dt$signature, function(x) unlist(strsplit(x, "_"))[2]) %>% unlist()

  p_attack4 <- ggplot(
    robustness.dt,
    aes(x = removed.pct, y = comp.pct, col = PPI_cat, linetype = PPI_cat)
  ) +
    geom_line(aes(group = signature, size = PPI_cat, shape = PPI_cat), show.legend = FALSE) + # colored by PPI_cat
    scale_color_manual(values = PPI_color_palette) +
    scale_size_manual(values = PPI_size_palette) + # Set line width
    # scale_shape_manual(values = c("HiG" = 16, "CTS" = 17, "HiGCTS" = 18)) +
    facet_wrap(~type) +
    geom_abline(slope = -1, intercept = 1, col = "gray", lty = 2) +
    theme(legend.position = c(0, 0), legend.justification = c(0, 0)) +
    labs(x = "% edges/vertex removed", y = "% of max. component remaining")
  vertex(p_attack4)
}

###################################################
# Fig J) boxplot of PageRank per PPIN
# original code: Code: 11.3_CTS_cardiac_network_ANND_pagerank.R
# original pdf: PageRank_GSE870383_v2.pdf
################################################################
{
  CHD <- readRDS(file = file.path(inputdir, (paste0("CHD_Cilia_Genelist.rds"))))

  df <- readRDS(file = "df_PAGERANK_strength_ANND.rewiring.P.rds") # !!!!!!!!!!!!!!!!!!!!!!!

  df <- df[grepl("_8$", df$signature), ]

  df <- rbind(
    subset(df, PPI_cat == "CTS"),
    subset(df, PPI_cat == "HiGCTS"),
    subset(df, PPI_cat == "HiG")
  )
  ## ensure the same order along x-axis
  # df$signature <- factor(df$signature, levels = signature_levels)
  df$label <- df$gene
  df$PCGC_AllCurated <- toupper(df$gene) %in% toupper(unlist(CHD["Griffin2023_PCGC_AllCurated"]))
  # Calculate top 5 significant genes within each box
  df5 <- df %>%
    filter(rank_by_PR <= 5) %>%
    ungroup()
  tb <- df5[, c("signature", "gene", "PPI_cat", "rank_by_PR", "PCGC_AllCurated")]
  write.table(tb, file = "table_top5_pagerank_perPPI.tsv", sep = "\t", row.names = F)

  df5_CHD <- subset(df5, PCGC_AllCurated == TRUE)
  (dim(df5_CHD)) # 0 18

  boxplot_pagerank <- ggplot(df, aes(x = signature, y = PageRank, colour = PPI_cat)) +
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
    stat_compare_means(
      comparisons = list(
        c("CTS_8", "HiGCTS_8"),
        c("CTS_8", "HiG_8"),
        c("HiGCTS_8", "HiG_8")
      ),
      method = "wilcox.test",
      label = "p.signif",
      p.adjust.method = "bonferroni",
      hide.ns = TRUE
    ) +
    scale_x_discrete(limits = unique(df$signature)) +
    labs(color = "PPI cat") # Optional: label for the color legend

  vertex(boxplot_pagerank)
}

###################################################
# Fig L) Edge weight density
# original code: 11.6_communities_edge_weight_distribution.R
# original pdf: edge_weight.pdf
###################################################
{
  graph_list_notsimplified <- readRDS(file = file.path(wd, "results", paste0(db, "_IID_graph_perState_notsimplified.rds")))

  # Remove duplicate vertices
  # correct_n_edges = readRDS(file.path(inputdir, '../correct_n_edges_HiG_STRING2.14.0.rds'))
  # for(g_name in unique(correct_n_edges$graph_id)){
  # 	vertices_to_remove = subset(correct_n_edges, graph_id == g_name)$vetex_index_to_remove
  # 	if(any(is.na(vertices_to_remove))) vertices_to_remove = vertices_to_remove[!is.na(vertices_to_remove)]
  # 	graph_list_notsimplified[[g_name]] = delete_vertices(graph_list_notsimplified[[g_name]], vertices_to_remove)
  # }
  N <- sapply(graph_list_notsimplified, vcount)

  graph_list <- lapply(graph_list_notsimplified, simplify, edge.attr.comb = "max") # !!!!!!!!!!!!!!!!!!!
  N2 <- sapply(graph_list, vcount)
  (all(N == N2)) #  TRUE

  graph_list <- graph_list[grepl("_8$", names(graph_list))]

  edge_data <- extract_edge_weights_by_category(graph_list, PPI_color_palette, "8")
  (head(edge_data, 3))
  #   sample PPI_cat edge_weight num_edges cluster_ID cluster_cat
  # 1  HiG_8     HiG        0.01      2279          8    unstable
  # 2  HiG_8     HiG        0.01      2279          8    unstable
  # 3  HiG_8     HiG        0.01      2279          8    unstable

  # Create plots for PPI category analysis
  category_plots <- plot_edge_weight_distributions(edge_data, PPI_color_palette)

  plot_density <- category_plots$density_plot

  edge_counts <- edge_data %>%
    group_by(PPI_cat) %>%
    summarise(
      edge_count = n(),
      .groups = "drop"
    )

  (edge_counts)
  #   PPI_cat edge_count
  #   <fct>        <int>
  # 1 CTS             65
  # 2 HiGCTS           1
  # 3 HiG           2279

  pairwise_pvals <- edge_data %>%
    pairwise_wilcox_test(
      formula = edge_weight ~ PPI_cat,
      p.adjust.method = "BH"
    )

  (pairwise_pvals)
  #   .y.         group1 group2    n1    n2 statistic     p p.adj p.adj.signif
  # * <chr>       <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>
  # 1 edge_weight CTS    HiGCTS    65     1       35  0.819 0.819 ns
  # 2 edge_weight CTS    HiG       65  2279    72188. 0.511 0.819 ns
  # 3 edge_weight HiGCTS HiG        1  2279     1019  0.733 0.819 ns
}

### Save Plots to Folder
plots <- list(
  boxplot_strength, # B
  p_attack4, # F
  boxplot_pagerank, # J
  plot_density # L
)

sizes <- list(
  c(5, 4), # B
  c(6, 3), # F
  c(5, 4), # J
  c(5, 5) # L
)

letters <- c("B", "F", "J", "L")

dir.create("plots_C8", showWarnings = FALSE)

file_names <- paste0(letters[seq_along(plots)], ".pdf")

for (i in seq_along(plots)) {
  filename <- file.path("plots_C8", file_names[i])
  plot_obj <- plots[[i]]

  if (!is.null(sizes[[i]])) {
    ggsave(
      filename = filename, plot = plot_obj,
      width = sizes[[i]][1], height = sizes[[i]][2]
    )
  } else {
    ggsave(filename = filename, plot = plot_obj)
  }
}
