require(dplyr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(scales)  # for color gradient
library(patchwork)  # to arrange plots
library(gridExtra)  # to arrange plots
library(pracma)
library(data.table)
library(ggpubr)
 
########## BEGINNING OF USER INPUT ##########

wd = "/Users/felixyu/Documents/IbarraSoria2018/"
setwd(paste0(wd, "results/PPI_weight/"))

celltype_specific_weight_version <- '10'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

db <- 'IbarraSoria2018'

CT_id <- c("endothelial.b", "cardiac.a") # critical transition clusters

CP_CTS <- "cardiac.a" # cardiac progenitor CTS

s = "combined" # specificity method

########## END OF USER INPUT ##########

file = paste0(db, '_STRING_graph_perState_simplified_',s,'weighted.rds')
graph_list <- readRDS( file)  
	
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


g_real = graph_list[[CP_CTS]]

pdf(file='simulation.pdf', width=15, height=5)
df = NULL
for(i in seq_along(graph_list)){
	g_real = graph_list[[i]]
	res = synthetic_simulation(g_real, main= names(graph_list)[i])
	grid.arrange(res$p_weights, res$p_line, res$p_AUC, ncol = 3)
	df = rbind(df,res$auc_summary) 
}
dev.off()

write.csv(df, "Simulation_AUC_summary.csv", row.names = FALSE)
network_colors = res$network_colors

# The real_PPİN is fragile — but still somewhat buffered by biological modularity.
# The degree-preserving null loses that structure and becomes pathologically brittle.

## wilcox test by category for all clusters ###
df$category = lapply(df$ID, function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist
df$cluster = lapply(df$ID, function(x) unlist(strsplit(x, split='_'))[2]) %>% unlist


df_summary <- df %>%
  group_by(category, network) %>%
  summarise(
    mean_AUC = mean(normalized_AUC, na.rm = TRUE),
    sd_AUC   = sd(normalized_AUC, na.rm = TRUE),
    n        = n(),
    se_AUC   = sd_AUC / sqrt(n)
  )
  
p1 <- ggplot(df, aes(x = network, y = normalized_AUC, fill = network)) +
  geom_boxplot(alpha = 0.6) +
  facet_wrap(~category, scales = "free_x") +
  #stat_compare_means(method = "anova", label.y = 1.1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(c("real_PPIN", "random"), 
						c("real_PPIN", "deg_preserving"),
						c("real_PPIN", "scale_free"),
						c("real_PPIN", "small_world")),
    label.y = c(1.05, 1.0)
  ) +
  labs(
    title = "Comparison of normalized AUC across networks within each category",
    x = "synthetic simulated network",
    y = "Normalized AUC"
  ) +
  scale_color_manual(values = network_colors) +
	  scale_fill_manual(values = network_colors)  +	
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9))
  
ggsave("bar_normalized_AUC_across_simulatedNetworks_Per_category.pdf", plot = p1, width = 10, height = 6)

### pairwise test for only transition clusters ##
df_summary <- df %>% 
  filter(cluster %in% as.character(CT_id)) %>%
  group_by(category, network) %>%
  summarise(
    mean_AUC = mean(normalized_AUC, na.rm = TRUE),
    sd_AUC   = sd(normalized_AUC, na.rm = TRUE),
    n        = n(),
    se_AUC   = sd_AUC / sqrt(n),
    .groups = "drop"   # <— this removes grouping and silences the message
  )  

df_plot <- df %>%
  filter(cluster %in% as.character(CT_id))

p2 <- ggplot(df_plot, aes(x = network, y = normalized_AUC, fill = network)) +
  geom_boxplot(alpha = 0.3, width = 0.6, outlier.shape = NA) +
  geom_line(aes(group = ID), color = "gray60", alpha = 0.6, linewidth = 0.6) +
  geom_point(size = 2.5, shape = 21, color = "black", 
             aes(fill = network), position = position_dodge(width = 0.4)) +
  facet_wrap(~category, scales = "free_x") +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,  alternative='greater', 
    label = "p.signif",
    comparisons = list(c("real_PPIN", "random"),
                       c("real_PPIN", "deg_preserving"),
                       c("real_PPIN", "scale_free"),
                       c("real_PPIN", "small_world")),
    label.y = c(1.05, 1.0)
  ) +
  labs(
    title = "Comparison of normalized AUC for CT clusters",
    x = "Synthetic simulated network",
    y = "Normalized AUC"
  ) +
  scale_fill_manual(values = network_colors) +
   theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)
  )
 ggsave("bar_normalized_AUC_across_simulatedNetworks_CTonly.pdf", plot = p2, width = 10, height = 6)

 
#   ## manually check ##
  
# HiG category
df_subset <- df_plot %>% filter(category == "HiG", network %in% c("real_PPIN", "deg_preserving")) %>% arrange(ID)
vals1 <- df_subset$normalized_AUC[df_subset$network == "real_PPIN"]
vals2 <- df_subset$normalized_AUC[df_subset$network == "deg_preserving"]
t.test(vals1, vals2, paired = TRUE, alternative = "greater") # p-value = 0.01927

# CTS category
df_subset <- df_plot %>% filter(category == "CTS", network %in% c("real_PPIN", "deg_preserving")) %>% arrange(ID)
vals1 <- df_subset$normalized_AUC[df_subset$network == "real_PPIN"]
vals2 <- df_subset$normalized_AUC[df_subset$network == "deg_preserving"]
t.test(vals1, vals2, paired = TRUE, alternative = "greater") # p-value = 0.09374

# HiGCTS category
df_subset <- df_plot %>% filter(category == "HiGCTS", network %in% c("real_PPIN", "deg_preserving")) %>% arrange(ID)
vals1 <- df_subset$normalized_AUC[df_subset$network == "real_PPIN"]
vals2 <- df_subset$normalized_AUC[df_subset$network == "deg_preserving"]
t.test(vals1, vals2, paired = TRUE, alternative = "greater") # p-value = 0.2963

# All categories together
df_subset <- df_plot %>% filter(network %in% c("real_PPIN", "deg_preserving")) %>% arrange(ID)
vals1 <- df_subset$normalized_AUC[df_subset$network == "real_PPIN"]
vals2 <- df_subset$normalized_AUC[df_subset$network == "deg_preserving"]
t.test(vals1, vals2, paired = TRUE, alternative = "greater") # p-value = 0.006846
t.test(vals1, vals2, paired = TRUE) # p-value = 0.01369
wilcox.test(vals1, vals2, paired = TRUE) # p-value = 0.0625
