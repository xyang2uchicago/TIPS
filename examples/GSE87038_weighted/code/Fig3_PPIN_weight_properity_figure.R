library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)  # ggarrange()
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(igraph)
library(rstatix)
library(brainGraph)

wd = "/Users/felixyu/Documents/GSE87038_weighted/"
source(paste0(wd, "code/celltype_specific_weight_v10.R"))

PPI_color_palette = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_palette = c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)
signature_levels = c(
  "HiG_1", "HiG_2", "HiG_3", "HiG_4", "HiG_5", "HiG_6", "HiG_9", "HiG_10",
  "HiG_12", "HiG_14", "HiG_17", "HiG_18", "HiG_19", "HiG_7", "HiG_11",
  "HiG_15", "HiG_16", "HiG_13", "HiG_8", "HiGCTS_7", "HiGCTS_11", "HiGCTS_15",
  "HiGCTS_16", "HiGCTS_16.1", "HiGCTS_8", "CTS_7", "CTS_11", "CTS_15", "CTS_16",
  "CTS_16.1", "CTS_13", "CTS_8", "HiGCTS_13"
)

setwd(paste0(wd, 'results/PPI_weight/'))
inputdir = "../../data/"


# refer to 11.2.0_weighted_graph_attack_robustness.R
s = "combined"
file = paste0('GSE87038_STRING_graph_perState_simplified_',s,'weighted.rds')
graph_list <- readRDS(file)  
	
(names(graph_list))
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"       "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"      "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"      
# [15] "HiG_11"      "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"     
# [29] "CTS_16"      "CTS_16.1"    "CTS_13"      "CTS_8"       "HiGCTS_13"    

CHD = readRDS( file=paste0(inputdir, 'CHD_Cilia_Genelist.rds'))

###################################################
# Fig A ) normalized node strength analysis
# original code: 11.3_CTS_cardiac_network_ANND_pagerank.R
# original pdf: normalized.node.strength_GSE87038_v2.pdf
################################################################
{
df = readRDS(file='df_PAGERANK_strength_ANND.rewring.P.rds')  #!!!!!!!!!!!!!!!!!!!!!!!
df = rbind(subset(df, PPI_cat=='CTS'),
					subset(df, PPI_cat=='HiGCTS'),
					subset(df, PPI_cat=='HiG')
					)
df$label=df$gene

## Fig A top) boxplot of normalized strength ########
df_median = df %>% group_by(signature) %>%
					summarise(median_normalized_strength = median(normalized.strength, na.rm = TRUE))
df_median$PPI_cat = lapply(df_median$signature, function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist 
df_median$PPI_cat = factor(df_median$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 				
	
violin_median_normalized.strength_wilcox = ggplot(df_median, aes(x = PPI_cat, y = median_normalized_strength, color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_palette) +
	  scale_fill_manual(values = PPI_color_palette) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category", y = "median of normalized node strength per PPI") +  # Label the axes
	  # Add statistical comparisons using stat_compare_means
	  stat_compare_means(
		aes(group = PPI_cat),  # Grouping by the 'PPI_cat' column
		comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")),  # Specify comparisons
		method = "wilcox.test",  # Non-parametric test (Wilcoxon)
		label = "p.signif",  # Show significance labels (e.g., **, *, ns)
		label.x = 1.5,  # Adjust x-position of the p-value text
		size = 4  # Adjust size of the p-value text
		,tip.length =0
	  ) + 
ggtitle('wilcox, median nor_strength')
vertex(violin_median_normalized.strength_wilcox)

##  Fig A bottom) violin plot of median of normalized strength per category	
################################################################
df$PCGC_AllCurated = df$gene %in% unlist(CHD[c('Griffin2023_PCGC_AllCurated')])
 
top_genes <- df %>%
  group_by(signature) %>%
  arrange(desc(normalized.strength)) %>%   # Sort in descending order of normalized.strength
  slice_head(n = 5)  # Take the top 5 rows for each signature
  	
# subset the CHD genes within top 5 
top_genes_CHD = subset(top_genes, PCGC_AllCurated==TRUE)
(dim(top_genes))  # [1] 165  18
(dim(top_genes_CHD))  # [1] 12 18

df$PPI_cat = lapply(df$signature %>% as.vector, function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist 
df$PPI_cat = factor(df$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 				
# df$signature = factor(df$signature, levels = signature_levels)  

boxplot_strength = ggplot(df, aes(x = signature, y = log10(normalized.strength), colour = PPI_cat)) + 
	  geom_boxplot(show.legend = TRUE) +  # Enable legend for the boxplot
	  scale_color_manual(values = PPI_color_palette) +
	  # Use ggrepel to avoid overlap and label top 5 genes based on normalized.strength
	  geom_text_repel(
		data = top_genes_CHD,  # Label only the top 5 genes
		aes(label = gene),
		size = 2 ,  # Adjust the size of the text labels
		box.padding = 0.5,  # Add space between the text and the data points
		point.padding = 0.5,  # Add space between the text and the points
		segment.color = 'grey50',  # Color for the line connecting the text to the points
		max.overlaps = 20,  # Max number of overlaps before labels stop being placed
		show.legend = FALSE  # Do not show text labels in the legend
	  ) +
	  theme(
		legend.position = 'none', # c(1, 1), 
		#legend.justification = c(0, 1),  # Place legend at top-right corner
		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
	  ) +
	  scale_x_discrete(limits = signature_levels) + # Ensure x-axis respects the order of 'signature'
	  labs(color = "PPI cat")   # Optional: label for the color legend
    
vertex(boxplot_strength)

}

###################################################
# Fig B. degree) boxplot of normalized strength for only transition clusters per category
# original code: 11.2.1_CTS_cardiac_network_strengthDistribution.R
# original pdf: boxplot_normalized_strength_GSE87038.pdf
################################################################
{
V_deg_dis = lapply(graph_list, function(x) igraph::degree_distribution(x, cumulative=TRUE )) %>% 
		  lapply(., function(x) x %>% 
				  as.data.frame(degree_distribution=x)  %>% 
				  mutate(k=seq_along(x))) %>% 
rbindlist(.,idcol=names(.))

colnames(V_deg_dis)[1:2]=c('signature','degree_distribution')
V_deg_dis$PPI_cat = lapply(V_deg_dis$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
		factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 

(table(V_deg_dis$PPI_cat))
#    CTS HiGCTS    HiG 
#    145     53   4370
V_deg_dis$cluster = lapply(V_deg_dis$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 

all(V_deg_dis$signature %in% names(graph_list))
V_deg_dis$n_nodes <- 0
for(i in seq_along(graph_list)){
	j = which(V_deg_dis$signature == names(graph_list)[i])
	V_deg_dis$n_nodes[j] <- vcount(graph_list[[i]])
	}
V_deg_dis$normalized_degree_distribution =  V_deg_dis$degree_distribution / (V_deg_dis$n_nodes-1)
 
boxplot_transition_degree = ggplot(subset(V_deg_dis, grepl("_7|_8|_11|_13|_15|_16|_16.1",signature)), 
				aes(x = factor(PPI_cat), y = normalized_degree_distribution, fill = PPI_cat)) + 
		  geom_boxplot() +
		  labs(x = "PPIN Category", y = "Normalized Degree Distribution", title = "PPINs for all transition clusters") +
		  theme_minimal() +
		  scale_fill_manual(values = PPI_color_palette) +
		  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")), 
                     p.adjust.method = "BH",  # Adjust p-values using Benjamini-Hochberg (BH) method
                     label = "p.signif")
vertex(boxplot_transition_degree)
}

V_deg_dis$normalized_strength_distribution =  V_deg_dis$strength_distribution / (V_deg_dis$n_nodes-1)

boxplot_transition_strength = ggplot(subset(V_deg_dis, grepl("_7|_8|_11|_13|_15|_16|_16.1",signature)), 
				aes(x = factor(PPI_cat), y = normalized_strength_distribution, fill = PPI_cat)) + 
		  geom_boxplot() +
		  labs(x = "PPIN Category", y = "Normalized Strength Distribution", title = "PPINs for all transition clusters") +
		  theme_minimal() +
		  scale_fill_manual(values = PPI_color_palette) +
		  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")), 
                     p.adjust.method = "BH",  # Adjust p-values using Benjamini-Hochberg (BH) method
                     label = "p.signif")
vertex(boxplot_transition_strength)


###################################################
# Fig B. strength) boxplot of normalized strength for only transition clusters per category
# original code: 11.2.1_CTS_cardiac_network_strengthDistribution.R
# original pdf: boxplot_normalized_strength_GSE87038.pdf
################################################################
{
V_deg_dis = lapply(graph_list, function(x) strength_distribution(x, cumulative=TRUE )) %>% 
		  lapply(., function(x) x %>% 
				  as.data.frame(strength_distribution=x)  %>% 
				  mutate(k=seq_along(x))) %>% 
rbindlist(.,idcol=names(.))

colnames(V_deg_dis)[1:2]=c('signature','strength_distribution')
V_deg_dis$PPI_cat = lapply(V_deg_dis$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
		factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 

table(V_deg_dis$PPI_cat)
 # CTS HiGCTS    HiG 
# 400    132   1214
V_deg_dis$cluster = lapply(V_deg_dis$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 

all(V_deg_dis$signature %in% names(graph_list)) #T
V_deg_dis$n_nodes <- 0
for(i in seq_along(graph_list)){
	j = which(V_deg_dis$signature == names(graph_list)[i])
	V_deg_dis$n_nodes[j] <- vcount(graph_list[[i]])
	}
V_deg_dis$normalized_strength_distribution =  V_deg_dis$strength_distribution / (V_deg_dis$n_nodes-1)
 
boxplot_transition_strength = ggplot(subset(V_deg_dis, grepl("_7|_8|_11|_13|_15|_16|_16.1",signature)), 
				aes(x = factor(PPI_cat), y = normalized_strength_distribution, fill = PPI_cat)) + 
		  geom_boxplot() +
		  labs(x = "PPIN Category", y = "Normalized Strength Distribution", title = "PPINs for all transition clusters") +
        theme_minimal() +
		  scale_fill_manual(values = PPI_color_palette) +
		  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")), 
                     p.adjust.method = "BH",  # Adjust p-values using Benjamini-Hochberg (BH) method
                     label = "p.signif")
vertex(boxplot_transition_strength)
}

###################################################
# Fig C) boxplopt of median PageRAnk per category
# original code: 11.3_CTS_cardiac_network_ANND_pagerank.R
# original pdf: PageRank_GSE870383_v2.pdf
################################################################
{
df = readRDS(file='df_PAGERANK_strength_ANND.rewring.P.rds')  #!!!!!!!!!!!!!!!!!!!!!!!
df = rbind(subset(df, PPI_cat=='CTS'),
					subset(df, PPI_cat=='HiGCTS'),
					subset(df, PPI_cat=='HiG')
					)
df$label=df$gene

df_median = df %>% group_by(signature) %>%
					summarise(pg.median = median(PageRank, na.rm = TRUE))
df_median$PPI_cat = lapply(df_median$signature, function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist 
df_median$PPI_cat = factor(df_median$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 				
		
violin_median_pagerank = ggplot(df_median, aes(x = PPI_cat, y = pg.median, color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_palette) +
	  scale_fill_manual(values = PPI_color_palette) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category", y = "median of PageRanks per PPI") +  # Label the axes
	  # Add statistical comparisons using stat_compare_means
	  stat_compare_means(
		aes(group = PPI_cat),  # Grouping by the 'PPI_cat' column
		comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")),  # Specify comparisons
		method = "wilcox.test",  # Non-parametric test (Wilcoxon)
		label = "p.signif",  # Show significance labels (e.g., **, *, ns)
		label.x = 1.5,  # Adjust x-position of the p-value text
		size = 4  # Adjust size of the p-value text
		,tip.length =0
	  ) +
	  ggtitle('wilcox test, median PA')

vertex(violin_median_pagerank)
}

###################################################
# Fig D) boxplot of %remained_fraction of targeted attack vs random attack
# original code: 11.2_CTS_cardiac_network_robustness.R
# original pdf: box_wilcox-test_attack_GSE87038.pdf
################################################################
{
   attack.edge.btwn = readRDS(file='attack.edge.btwn.rds')
   attack.vertex.btwn = readRDS( file='attack.vertex.btwn.rds')
   failure.vertex = readRDS(paste0('failure.vertex_100_simplified_',s,'weighted.rds') )

   failure.edge = readRDS(paste0('failure.edge_100_simplified_',s,'weighted.rds')) 
   failure.dt <- rbind(failure.edge, failure.vertex)   
 
   colnames(failure.dt)[1] ='signature'
   colnames(attack.vertex.btwn)[1] = 'signature'
   (dim(failure.dt))  #[1] 160135      6
   
   robustness.dt <- rbind(failure.dt, attack.vertex.btwn, attack.edge.btwn)  
   (dim(robustness.dt))  #[1] 320270      6
   robustness.dt$PPI_cat = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
			factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 
 
    robustness.dt$experiment = ifelse(grepl('edge', robustness.dt$type), 'edge', 'vertex')
	robustness.dt$measure= factor(robustness.dt$measure, levels = c("random" ,  "btwn.cent"))

    (table(robustness.dt$type, robustness.dt$measure))
    #                          random btwn.cent
    #   Random edge removal    151870         0
    #   Random vertex removal    8245         0
    #   Targeted edge attack        0    151870
    #   Targeted vertex attack      0      8245


   robustness.dt$type = factor(robustness.dt$type,
			levels = c("Random edge removal" , "Targeted edge attack"  ,  "Random vertex removal" , "Targeted vertex attack") ) 
   robustness.dt$cluster = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 
                              

## then, Plot the Wilcox results for visualization and manually add fold chagnes !!!
boxplot_remained_fraction = ggplot(
  robustness.dt,
  aes(
    x = type, 
    y = comp.pct,
    fill = measure,
    color = PPI_cat,   # keep color
    size = experiment  # bold vs thin outline
  )
) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +
  facet_wrap(~ PPI_cat, ncol = 3) +
  scale_color_manual(values = PPI_color_palette, guide = "none") +  # keep color but hide legend
  scale_fill_manual(values = c("random" = "grey", "btwn.cent" = "white")) +
  scale_size_manual(values = c("edge" = 1.2, "vertex" = 0.4), name = "Experiment") +
  geom_signif(
    comparisons = list(
      c("Random edge removal", "Targeted edge attack"),
      c("Random vertex removal", "Targeted vertex attack")
    ),
    map_signif_level = TRUE,
    step_increase = 0.1,  # Adjusts spacing between the lines
    aes(group = type),
    test = "wilcox.test"  # Perform a t-test to calculate significance
  ) +
  theme_minimal() +
  labs(
    title = "Comparison of Robustness Measures by Type and State",
    x = "PPI Category", 
    y = "Component Percentage Remaining"
  ) +
  theme(legend.position = "top",
    axis.text.x = element_blank(),   # remove x-axis tick labels
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold")
  )

vertex(boxplot_remained_fraction)

## finally, manually add the threshold of fold changes  ###############
robustness.dt <- robustness.dt %>%
  group_by(PPI_cat, type) %>%
  mutate(
    mean_comp_pct = mean(comp.pct, na.rm = TRUE)  # Calculate mean comp.pct for each group
  ) %>%
  ungroup() 

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
# 1 CTS                1.06                1.54
# 2 HiGCTS             0.991               1.50
# 3 HiG                0.964               1.04

# TODO: significant if >1.5


### additionally, wilcox-test the between-group changes among observed PPINs
tmp = subset(robustness.dt,measure=='btwn.cent')
dim(tmp) #[1] 4604    9
wilcox.test(subset(tmp,PPI_cat=='HiG')$comp.pct, subset(tmp,PPI_cat=='CTS')$comp.pct)
# W = 138842816, p-value < 2.2e-16

wilcox.test(subset(tmp,PPI_cat=='HiG')$comp.pct, subset(tmp,PPI_cat=='HiGCTS')$comp.pct)
# W = 30162248, p-value < 2.2e-16

wilcox.test(subset(tmp,PPI_cat=='HiGCTS')$comp.pct, subset(tmp,PPI_cat=='CTS')$comp.pct)
# W = 151165, p-value = 0.001761
}

###################################################
# Fig E left) ARC plot of vertex attack
# original code: 11.2_CTS_cardiac_network_robustness.R
# original pdf: attack_GSE87038.pdf
################################################################
{  
# failure.vertex = readRDS(file=paste0('failure.vertex_100_simplified_',s,'weighted.rds'))  
# failure.edge = readRDS(paste0('failure.edge_100_simplified_',s,'weighted.rds'))
# failure.dt <- rbind(failure.edge, failure.vertex)   
failure.dt = readRDS(file=paste0('failure.vertex_100_simplified_',s,'weighted.rds'))  
colnames(failure.dt)[1] ='signature'

# attack.edge.btwn = readRDS( file=paste0('attack.edge.btwn.rds'))
attack.vertex.btwn = readRDS( file=paste0('attack.vertex.btwn.rds'))
colnames(attack.vertex.btwn)[1] = 'signature'	 

robustness.dt <- rbind(failure.dt, attack.vertex.btwn[,1:6])  #, attack.edge.btwn[,1:6])  
(dim(robustness.dt)) #  16490     6
robustness.dt$PPI_cat = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
		factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 
head(robustness.dt, 3)


robustness.dt$measure %>% unique
#[1] "random"     "btwn.cent"
robustness.dt = subset(robustness.dt ,measure != 'degree')
robustness.dt$experiment = ifelse(grepl('edge', robustness.dt$type), 'edge', 'vertex')
robustness.dt$measure= factor(robustness.dt$measure, levels = c("random" ,  "btwn.cent"))

robustness.dt$type = factor(robustness.dt$type,
		levels = c("Random edge removal" , "Targeted edge attack"  ,  "Random vertex removal" , "Targeted vertex attack") )
robustness.dt$cluster = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist                                
 
p_attack4 <- ggplot(robustness.dt,
                  aes(x=removed.pct, y=comp.pct, col=PPI_cat, linetype=PPI_cat)) +  
				geom_line(aes(group=signature,  size = PPI_cat, shape=PPI_cat), show.legend = FALSE) +  # colored by PPI_cat
				scale_color_manual(values = PPI_color_palette) +
				scale_size_manual(values = PPI_size_palette) +  # Set line width
				#scale_shape_manual(values = c("HiG" = 16, "CTS" = 17, "HiGCTS" = 18)) +  
				facet_wrap(~ type) + 
				geom_abline(slope=-1, intercept=1, col='gray', lty=2) +
				 theme(legend.position=c(0, 0), legend.justification=c(0, 0)) +
				 labs(x='% edges/vertex removed', y='% of max. component remaining')
vertex(p_attack4)
}

###################################################
# Fig E right) Boxplot of AUC for vertex attack
# original code: 11.2_CTS_cardiac_network_robustness.R
# original pdf: box_wilcox-test_attack_AUC_GSE87038.pdf
################################################################
{
library(MLmetrics)
library(sm)

IDs_of_CTS = c('7', '8', '11', '13', '15', '16', '16.1')

observed_auc_list = list()
for(j in names(graph_list)){
	observed_auc_list[[j]] = Area_Under_Curve(subset(robustness.dt, signature==j & type=='Targeted vertex attack')$removed.pct, 
				 subset(robustness.dt, signature==j & type=='Targeted vertex attack')$comp.pct ) 
} 
df_AUC = data.frame(auc=observed_auc_list %>% unlist, 
			signature = names(observed_auc_list), 
			PPI_cat = lapply(names(observed_auc_list), function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist)
df_AUC$PPI_cat  = factor(df_AUC$PPI_cat , levels = c('CTS',"HiGCTS" ,"HiG"))

boxplot_AUC_vertex_attack = ggplot(df_AUC, aes(x = PPI_cat, y = auc, fill = PPI_cat)) + 
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) +  # Dodge the boxes for each type per PPI_cat
  #facet_wrap(experiment ~ PPI_cat, ncol = 3) +  # Facet by PPI_cat, each PPI_cat gets a row
  scale_fill_manual(values = PPI_color_palette) +
  #scale_fill_manual(values = c("random" = "grey", "btwn.cent" = "orange")) +
  geom_signif(
    comparisons = list(
      c("CTS", "HiGCTS"), c("HiGCTS", "HiG"),c("CTS", "HiG")
    ),
    map_signif_level = TRUE,
    step_increase = 0.1,  # Adjusts spacing between the lines
    #aes(group = type),
    test = "wilcox.test"  # Perform a t-test to calculate significance
	#, p.adjust.method = "holm"   # default is 
  ) + 
  theme_minimal() +
  labs(
    title = "Comparison of Robustness Measures by cent.btw",
    x = "PPI Category", 
    y = "AUC of targeted vertex attack"
  ) +
  theme(legend.position = "top")
vertex(boxplot_AUC_vertex_attack)
}

###################################################
# Fig F) violin plot of median betweenness per category
# original code: Code: 11.3_CTS_cardiac_network_ANND_pagerank.R ??
# original pdf: BetweennessCentrality_GSE870383_v2.pdf
################################################################
{
df_BC = read.table(file='df_betweeness.tsv',sep='\t', header=T) 
df_BC$PPI_cat = factor(df_BC$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 			
## reorder df$df$signature to be #E7298A #7570B3 and #E6AB02
df5 <- df_BC %>%
	  filter(rank_by_BC <= 5 & BetweennessCentrality>0) %>%
	  ungroup()

df_median = df_BC %>% group_by(signature) %>%
					summarise(bc.median = median(BetweennessCentrality, na.rm = TRUE)) %>%
					as.data.frame()
df_median$PPI_cat = lapply(df_median$signature %>% as.vector, function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist 
df_median$PPI_cat = factor(df_median$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 				

violin_median_bc_wilcox = ggplot(df_median, aes(x = PPI_cat, y = bc.median , color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3, drop = FALSE) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_palette) +
	  scale_fill_manual(values = PPI_color_palette) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category", y = "median of BC per PPI") +  # Label the axes
	  # Add statistical comparisons using stat_compare_means
	  stat_compare_means(
		aes(group = PPI_cat),  # Grouping by the 'PPI_cat' column
		comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")),  # Specify comparisons
		method = "wilcox.test",  # Non-parametric test (Wilcoxon)
		label = "p.signif",  # Show significance labels (e.g., **, *, ns)
		label.x = 1.5,  # Adjust x-position of the p-value text
		size = 4  # Adjust size of the p-value text
		,tip.length =0
	  ) +
	  ylim(0, NA)  +  # Start from 0, let ggplot choose upper limit
	  ggtitle('wilcox-test, median BC')
 
plot(violin_median_bc_wilcox)
}

###################################################
# Fig G _new) cumulative distribution of normalized strength
# original code: 11.6_communities_edge_weight_distribution.R
# original pdf: edge_weight.pdf
################################################################
# EB = readRDS(file='cluster_edge_betweenness_list.rds')
# # HiG line with other categories as dots at corresponding positions
# nEB = lengths(EB)
# plots = plot_nEB_ggplot(nEB, PPI_color_palette, method='wilcox.test')
# grid.arrange(plots$line_plot, plots$boxplot, ncol = 2)
# dev.copy2pdf(file='community_number.pdf', width=10)

unstable_cluster_ID = c('7', '8', '11', '13', '15', '16', '16.1')

# Extract edge weight data by PPI category
edge_data <- extract_edge_weights_by_category(graph_list, PPI_color_palette, unstable_cluster_ID)
(head(edge_data, 3))
#   sample PPI_cat edge_weight num_edges cluster_ID cluster_cat
# 1  HiG_1     HiG 0.007203960      4978          1      stable
# 2  HiG_1     HiG 0.007968395      4978          1      stable
# 3  HiG_1     HiG 0.064959052      4978          1      stable

# Create plots for PPI category analysis
category_plots <- plot_edge_weight_distributions(edge_data, PPI_color_palette)

(category_plots$summary_stats)
# # A tibble: 3 × 4
#   PPI_cat mean_weight median_weight total_edges
#   <fct>         <dbl>         <dbl>       <int>
# 1 CTS           0.049         0.033         905
# 2 HiGCTS        0.082         0.049         137
# 3 HiG           0.017         0.01       150796

EB = readRDS(file='../cluster_edge_betweenness_list.rds')
# HiG line with other categories as dots at corresponding positions
nEB = lengths(EB)

community_plots = plot_nEB_ggplot(nEB, PPI_color_palette, method='wilcox.test')

edge_weight_community = grid.arrange(category_plots$density_plot,
							community_plots$line_plot, 
							ncol = 2)



###################################################
# Fig H) boxplot of betweenness per PPIN
# original code: Code: 11.3_CTS_cardiac_network_ANND_pagerank.R ??
# original pdf: BetweennessCentrality_GSE870383_v2.pdf
################################################################
{
CHD = readRDS( file=paste0(inputdir, 'CHD_Cilia_Genelist.rds'))

df_BC = read.table(file='df_betweeness.tsv',sep='\t', header=T) 
df_BC$PPI_cat = factor(df_BC$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG'))
df_BC$PCGC_AllCurated = toupper(df_BC$gene) %in% toupper(unlist(CHD['Griffin2023_PCGC_AllCurated']))
 			
df5 <- df_BC %>%
	  filter(rank_by_BC <= 5 & BetweennessCentrality>0) %>%
	  ungroup()
 	
df5_CHD = subset(df5, PCGC_AllCurated==TRUE)
(dim(df5_CHD))  # [1] 17  6

df_BC = rbind(subset(df_BC, PPI_cat=='CTS'),
					subset(df_BC, PPI_cat=='HiGCTS'),
					subset(df_BC, PPI_cat=='HiG')
					)
# df_BC$signature <- factor(df_BC$signature, levels = signature_levels)
boxplot_bc_log10 <- ggplot(df_BC, aes(
				x = signature, 
				y = log10(BetweennessCentrality + 1), 
			colour = PPI_cat)) +
  geom_boxplot(position = "dodge2") +
  scale_color_manual(values = PPI_color_palette) +
  geom_text(
    data = df5_CHD,
    aes(label = gene),
    size = 3.5,             # larger text
    vjust = -0.5,           # move above
    nudge_y = 0.2,          # space above box
    check_overlap = TRUE
  ) +
  theme(
    legend.position = 'none',
				#legend.justification = c(1, 1),  # Place legend at top-right corner
				axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
			)  + 
			geom_text(data = subset(df_BC, PCGC_AllCurated == TRUE) %>% filter(rank_by_BC <= 5),
					aes(label = gene), size = 2, 
					hjust = -0.1, vjust = 0, check_overlap = TRUE)  + # Adjust text labels
			scale_x_discrete(limits = signature_levels) +
			scale_color_manual(values = PPI_color_palette) 

vertex(boxplot_bc_log10)
}

###################################################
# Fig I) boxplot of pagerange per PPIN
# original code: Code: 11.3_CTS_cardiac_network_ANND_pagerank.R
# original pdf: PageRank_GSE870383_v2.pdf
################################################################
{
CHD = readRDS( file=paste0(inputdir, 'CHD_Cilia_Genelist.rds'))

df = readRDS(file='df_PAGERANK_strength_ANND.rewring.P.rds')  #!!!!!!!!!!!!!!!!!!!!!!!

df = rbind(subset(df, PPI_cat=='CTS'),
					subset(df, PPI_cat=='HiGCTS'),
					subset(df, PPI_cat=='HiG')
					)
## ensure the same order along x-axis					
# df$signature <- factor(df$signature, levels = signature_levels)
df$label=df$gene
df$PCGC_AllCurated = toupper(df$gene) %in% toupper(unlist(CHD['Griffin2023_PCGC_AllCurated']))	
# Calculate top 5 significant genes within each box
df5 <- df %>%
	  filter(rank_by_PR <= 5) %>%
	  ungroup()	
tb = df5[, c('signature','gene','PPI_cat','rank_by_PR','PCGC_AllCurated')]
write.table(tb, file= 'table_top5_pagerank_perPPI.tsv', sep='\t', row.names=F)
 	  
df5_CHD = subset(df5, PCGC_AllCurated==TRUE)
(dim(df5_CHD))  # [1] 13 18

boxplot_pagerank = ggplot(df, aes(x = signature,y = PageRank, colour = PPI_cat)) +
			  geom_boxplot(show.legend = TRUE) +  # Enable legend for the boxplot
  scale_color_manual(values = PPI_color_palette) +
			  geom_text(data=df5_CHD, aes(label = gene),  # data=df5
						size = 2,  # Adjust the size of the text labels
						hjust = -0.1, vjust = 0, 
						check_overlap = TRUE) +  # Avoid text overlap
			  theme(legend.position = 'none', #c(1, 1), 
					#legend.justification = c(1, 1),  # Place legend at top-right corner
					axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
			  scale_x_discrete(limits = signature_levels) +
			  labs(color = "PPI cat")  # Optional: label for the color legend

vertex(boxplot_pagerank)
}

### Save Plots to Folder
plots <- list(
	violin_median_normalized.strength_wilcox,  # A
	boxplot_strength,  # B   boxplot_degree
    boxplot_transition_strength + theme(legend.position = "none"),  # C   boxplot_transition_degree
    violin_median_pagerank,   # D 
    boxplot_remained_fraction, # E
    p_attack4, # F
    boxplot_AUC_vertex_attack, # G
    violin_median_bc_wilcox, # H
    edge_weight_community,   #  I
    boxplot_bc_log10, # J
    boxplot_pagerank # K
)

sizes <- list(
  c(4, 2), # A
  c(7, 4), # B
  c(3, 3), # C
  c(3, 3), # D
  NULL, # E
  c(6, 3), # F
  c(2, 3), # G
  c(4, 2), # H
  c(6,3), # I
  c(7, 4), # J
  c(7, 4) # K
)

dir.create("plots", showWarnings = FALSE)

file_names <- paste0(LETTERS[seq_along(plots)], ".pdf")

for (i in seq_along(plots)) {
  filename <- file.path("plots", file_names[i])
  plot_obj <- plots[[i]]
  
  if (!is.null(sizes[[i]])) {
    ggsave(filename = filename, plot = plot_obj,
           width = sizes[[i]][1], height = sizes[[i]][2])
  } else {
    ggsave(filename = filename, plot = plot_obj)
  }
}

ggsave(filename = 'plots/I.pdf', plot = edge_weight_community,
           width = 6, height = 3)