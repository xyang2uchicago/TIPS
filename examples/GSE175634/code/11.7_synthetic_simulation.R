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
 

source('E:/Git_Holly/TIPS-main/celltype_specific_weight_v10.R')

setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9')
db = 'GSE175634_iPSC'

# refer to 11.2.0_weighted_graph_attack_robustness.R
s = "combined"
file = paste0('GSE175634_STRING_graph_perState_simplified_',s,'weighted.rds')
graph_list <- readRDS( file)  
	
names(graph_list)
 # [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"           "HiG_5"          
 # [7] "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"   
# [13] "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"     
# [19] "CTS_endoderm"    "CTS_CP"          "CTS_CP.1"    


g_real = graph_list[["CTS_CP.1"]]
# save(g_real, file = 'E:/Git_Holly/TIPS-main/data/graph_CTS_CP.1.RData')
pdf(file='simulation.pdf', width=15, height=5)
df = NULL
for(i in 1:length(graph_list)){
	g_real = graph_list[[i]]
	res = synthetic_simulation(g_real, main= names(graph_list)[i])
	grid.arrange(res$p_weights, res$p_line, res$p_AUC, ncol = 3)
	df = rbind(df,res$auc_summary) 
}
dev.off()

write.csv(df, "Simulation_AUC_summary.csv", row.names = FALSE)
network_colors = res$network_colors

# The real CP_PPİN is fragile — but still somewhat buffered by biological modularity.
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
  
ggplot(df, aes(x = network, y = normalized_AUC, fill = network)) +
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
  
dev.copy2pdf(file='bar_normalized_AUC_across_simulatedNetworks_Per_category.pdf')

### pairwise test for only 3 transitin clusters ##
df_summary <- df %>% 
  filter(cluster %in% c("endoderm", "CP", "muscle")) %>%
  group_by(category, network) %>%
  summarise(
    mean_AUC = mean(normalized_AUC, na.rm = TRUE),
    sd_AUC   = sd(normalized_AUC, na.rm = TRUE),
    n        = n(),
    se_AUC   = sd_AUC / sqrt(n),
    .groups = "drop"   # <— this removes grouping and silences the message
  )  

df_plot <- df %>%
  filter(cluster %in% c("endoderm", "CP", "muscle"))

ggplot(df_plot, aes(x = network, y = normalized_AUC, fill = network)) +
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
    title = "Comparison of normalized AUC for 3 CT clusters",
    x = "Synthetic simulated network",
    y = "Normalized AUC"
  ) +
  scale_fill_manual(values = network_colors) +
   theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)
  )
 dev.copy2pdf(file='bar_normalized_AUC_across_simulatedNetworks_CTonly.pdf')

 
  ## manually check ##
 df_subset <- df_plot %>% filter(category == "HiG", network %in% c("CP_PPIN", "deg_preserving"))
 t.test(normalized_AUC ~ network, data = df_subset, paired = TRUE, alternative='greater') # 0.078
 df_subset <- df_plot %>% filter(category == "CTS", network %in% c("CP_PPIN", "deg_preserving"))
 t.test(normalized_AUC ~ network, data = df_subset, paired = TRUE, alternative='greater') # 0.078
 df_subset <- df_plot %>% filter(category == "HiGCTS", network %in% c("CP_PPIN", "deg_preserving"))
 t.test(normalized_AUC ~ network, data = df_subset, paired = TRUE, alternative='greater') # 0.04966
 
  df_subset <- df_plot %>% filter(network %in% c("CP_PPIN", "deg_preserving"))
  t.test(normalized_AUC ~ network, data = df_subset, paired = TRUE, alternative='greater') # 0.001256
  t.test(normalized_AUC ~ network, data = df_subset, paired = TRUE ) #  p-value = 0.002512
  wilcox.test(normalized_AUC ~ network, data = df_subset, paired = TRUE )   # 0.003906