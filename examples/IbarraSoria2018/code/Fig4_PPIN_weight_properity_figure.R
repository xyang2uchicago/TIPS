library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(igraph)
library(rstatix)
 
library(brainGraph)
wd = "/Users/felixyu/Documents/IbarraSoria2018/"
source(paste0(wd, "code/celltype_specific_weight_v9.R"))

PPI_color_platte = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_platte = c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)

setwd(paste0(wd, "results/PPI_weight/"))
inputdir = "../"

# refer to 11.2.0_weighted_graph_attack_robustness.R
s = "combined"
file = paste0('2018_STRING_graph_perState_simplified_',s,'weighted.rds')
graph_list <- readRDS( file)  
	
(names(graph_list))
#  [1] "HiG_extraembryonicMesoderm" "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [5] "HiG_blood"                  "HiG_mesodermProgenitors"    "HiG_presomiticMesoderm.b"   "HiG_presomiticMesoderm.a"  
#  [9] "HiG_somiticMesoderm"        "HiG_mixedMesoderm.a"        "HiG_pharyngealMesoderm"     "HiG_mixedMesoderm.b"       
# [13] "HiG_cardiac.b"              "HiG_cardiac.c"              "HiG_endothelial.b"          "HiG_cardiac.a"             
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"           "CTS_endothelial.b"          "CTS_cardiac.a" 

CHD = readRDS( file=paste0(inputdir, 'CHD_Cilia_Genelist.rds'))

###################################################
# Fig A) Boxplot of normalized strength per PPI signature with top CHD genes labeled
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

df$PCGC_AllCurated = toupper(df$gene) %in% toupper(unlist(CHD['Griffin2023_PCGC_AllCurated']))
 
top_genes <- df %>%
  group_by(signature) %>%
  arrange(desc(normalized.strength)) %>%   # Sort in descending order of normalized.strength
  slice_head(n = 5)  # Take the top 5 rows for each signature
  	
# subset the CHD genes within top 5 
top_genes_CHD = subset(top_genes, PCGC_AllCurated==TRUE)
(dim(top_genes))  # [1] 100  18
(dim(top_genes_CHD))  # [1]  9 18

pr = ggplot(df, aes(x = signature, y = log10(normalized.strength), colour = PPI_cat)) + 
	  geom_boxplot(show.legend = TRUE) +  # Enable legend for the boxplot
	  scale_color_manual(values = PPI_color_platte) +
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
		legend.position = 'none',
		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
	  ) +
	  scale_x_discrete(limits = unique(df$signature)) + # Ensure x-axis respects the order of 'signature'
	  labs(color = "PPI cat")   # Optional: label for the color legend
    
(pr)
}

###################################################
# Fig B) Boxplot of median normalized degree for transition clusters per PPI category
# original code: 11.1.2_CTS_cardiac_network_degreeDistribution.R
# original pdf: boxplot_normalized_degree_GSE87038.pdf
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
#     30     18   4000 
V_deg_dis$cluster = lapply(V_deg_dis$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 

all(V_deg_dis$signature %in% names(graph_list))
V_deg_dis$n_nodes <- 0
for(i in seq_along(graph_list)){
	j = which(V_deg_dis$signature == names(graph_list)[i])
	V_deg_dis$n_nodes[j] <- vcount(graph_list[[i]])
	}
V_deg_dis$normalized_degree_distribution =  V_deg_dis$degree_distribution / (V_deg_dis$n_nodes-1)
 
g2 = ggplot(subset(V_deg_dis, grepl("_cardiac.a|_endothelial.b",signature)), 
				aes(x = factor(PPI_cat), y = normalized_degree_distribution, fill = PPI_cat)) + 
		  geom_boxplot() +
		  labs(x = "PPIN Category", y = "Normalized Degree Distribution", title = "PPINs for all transition clusters") +
		  theme_minimal() +
		  scale_fill_manual(values = PPI_color_platte) +
		  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")), 
                     p.adjust.method = "BH",  # Adjust p-values using Benjamini-Hochberg (BH) method
                     label = "p.signif")
(g2)
}

###################################################
# Fig C) Violin plot of median normalized strength per network
# original code: 11.3_CTS_cardiac_network_ANND_pagerank.R
# original pdf: normalized.node.strength_GSE87038_v2.pdf
################################################################
{
    df_median = df %>% group_by(signature) %>%
                        summarise(median_normalized_strength = median(normalized.strength, na.rm = TRUE))
    df_median$PPI_cat = lapply(df_median$signature, function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist 
    df_median$PPI_cat = factor(df_median$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 		


    g2_violin <- ggplot(df_median, aes(
        x = PPI_cat,
        y = median_normalized_strength, color = PPI_cat, fill = PPI_cat
    )) +
        geom_violin(alpha = 0.3) + # Violin plot with transparency
        scale_color_manual(values = PPI_color_platte) +
        scale_fill_manual(values = PPI_color_platte) +
        theme_minimal() +
        theme(legend.position = "none") + # , axis.text.y = element_blank(), axis.title.y = element_blank()) +
        labs(x = "PPI category", y = "log10. median of normalized node strength per PPI") + # Label the axes
        # Add statistical comparisons using stat_compare_means
        stat_compare_means(
            aes(group = PPI_cat), # Grouping by the 'PPI_cat' column
            comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")), # Specify comparisons
            method = "wilcox.test", # Non-parametric test (Wilcoxon)
            label = "p.signif", # Show significance labels (e.g., **, *, ns)
            label.x = 1.5, # Adjust x-position of the p-value text
            size = 4 # Adjust size of the p-value text
            , tip.length = 0
        ) +
        ggtitle("wilcox, median nr_strength")

    (g2_violin)
}

###################################################
# Fig D) Violin plot of median PageRank per PPI category
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
		
violin_median_page_wilcox = ggplot(df_median, aes(x = PPI_cat, y = pg.median, color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
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

(violin_median_page_wilcox)
}

###################################################
# Fig E) boxplot of %remained_fraction of targeted attack vs random attack
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
   (dim(failure.dt))  # [1] 196474      6
   
   robustness.dt <- rbind(failure.dt, attack.vertex.btwn, attack.edge.btwn)  
   (dim(robustness.dt))  # [1] 392948      6
   robustness.dt$PPI_cat = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
			factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 
 
    robustness.dt$experiment = ifelse(grepl('edge', robustness.dt$type), 'edge', 'vertex')
	robustness.dt$measure= factor(robustness.dt$measure, levels = c("random" ,  "btwn.cent"))

    (table(robustness.dt$type, robustness.dt$measure))
    #                          random btwn.cent
    #   Random edge removal    188339         0
    #   Random vertex removal    8135         0
    #   Targeted edge attack        0    188339
    #   Targeted vertex attack      0      8135


   robustness.dt$type = factor(robustness.dt$type,
			levels = c("Random edge removal" , "Targeted edge attack"  ,  "Random vertex removal" , "Targeted vertex attack") ) 
   robustness.dt$cluster = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 
                              

## Plot the Wilcox results for visualization and manually add fold chagnes !!!
g = ggplot(
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
  scale_color_manual(values = PPI_color_platte, guide = "none") +  # keep color but hide legend
  scale_fill_manual(values = c("random" = "grey", "btwn.cent" = "white")) +
  scale_size_manual(values = c("edge" = 1.2, "vertex" = 0.4), name = "Experiment") +
  geom_signif(
    comparisons = list(
      c("Random edge removal", "Targeted edge attack"),
      c("Random vertex removal", "Targeted vertex attack")
    ),
    map_signif_level = TRUE,
    step_increase = 0.1,
    aes(group = type),
    test = "wilcox.test"
  ) +
  theme_minimal() +
  labs(
    title = "Comparison of Robustness Measures",
    x = "PPI category", 
    y = "Component Percentage Remaining"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),   # remove x-axis tick labels
    axis.ticks.x = element_blank(),  # remove x-axis ticks
    strip.text = element_text(face = "bold")
  )

(g)

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
# 1 CTS                0.957               1.40
# 2 HiGCTS             0.956               1.35
# 3 HiG                0.921               1.04
}

###################################################
# Fig F) attack–robustness curve plot of vertex attack
# original code: 11.2_CTS_cardiac_network_robustness.R
# original pdf: attack_GSE87038.pdf
################################################################
{  
failure.dt = readRDS(file=paste0('failure.vertex_100_simplified_',s,'weighted.rds'))  
colnames(failure.dt)[1] ='signature'

attack.vertex.btwn = readRDS( file=paste0('attack.vertex.btwn.rds'))
colnames(attack.vertex.btwn)[1] = 'signature'	 

robustness.dt <- attack.vertex.btwn[,1:6]
(dim(robustness.dt))  # [1] 8135    6
robustness.dt$PPI_cat = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
		factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 
(head(robustness.dt, 3))

(robustness.dt$measure %>% unique)
#[1] "btwn.cent"

robustness.dt = subset(robustness.dt ,measure == 'btwn.cent')

robustness.dt$experiment = 'vertex'

robustness.dt$measure= factor(robustness.dt$measure, levels = c("random" ,  "btwn.cent"))

robustness.dt$type = factor(robustness.dt$type,
		levels = c("Random edge removal" , "Targeted edge attack"  ,  "Random vertex removal" , "Targeted vertex attack") )
robustness.dt$cluster = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist                                
 
p_attack4 <- ggplot(
  subset(robustness.dt, type == "Targeted vertex attack"),
  aes(
    x = removed.pct,
    y = comp.pct,
    color = PPI_cat,
    linetype = PPI_cat,
    size = PPI_cat,
    group = signature
  )
) +
  geom_line() +
  scale_color_manual(
    values = PPI_color_platte,
    name = "PPI_category"
  ) +
  scale_linetype_manual(
    values = c(
      "CTS" = "solid",         # Blue solid
      "HiGCTS" = "longdash",   # Pink long dash
      "HiG" = "33"             # Yellow short dash
    ),
    name = "PPI_category"
  ) +
  scale_size_manual(
    values = c(
      "CTS" = 1.3,      # thickest
      "HiGCTS" = 0.9,   # thinner pink
      "HiG" = 0.9       # thinner yellow
    ),
    guide = "none"  # hide size legend (optional)
  ) +
  geom_abline(
    slope = -1, intercept = 1,
    col = "gray", linetype = "dashed", size = 1.0
  ) +
  labs(
    x = "Fraction of vertices removed",
    y = "Fraction of largest component remaining"
  ) +
  theme_gray() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

(p_attack4)
}

###################################################
# Fig G) Boxplot of AUC for vertex attack
# original code: 11.2_CTS_cardiac_network_robustness.R
# original pdf: box_wilcox-test_attack_AUC_GSE87038.pdf
################################################################
{
library(MLmetrics)
library(sm)

IDs_of_CTS = c('cardiac.a', 'endothelial.b')

observed_auc_list = list()
for(j in names(graph_list)){
	observed_auc_list[[j]] = Area_Under_Curve(subset(robustness.dt, signature==j & type=='Targeted vertex attack')$removed.pct, 
				 subset(robustness.dt, signature==j & type=='Targeted vertex attack')$comp.pct ) 
} 
df_AUC = data.frame(auc=observed_auc_list %>% unlist, 
			signature = names(observed_auc_list), 
			PPI_cat = lapply(names(observed_auc_list), function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist)
df_AUC$PPI_cat  = factor(df_AUC$PPI_cat , levels = c('CTS',"HiGCTS" ,"HiG"))

g3 = ggplot(df_AUC, aes(x = PPI_cat, y = auc, fill = PPI_cat)) + 
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) +  # Dodge the boxes for each type per PPI_cat
  scale_fill_manual(values = PPI_color_platte) +
  geom_signif(
    comparisons = list(
      c("CTS", "HiGCTS"), c("HiGCTS", "HiG"),c("CTS", "HiG")
    ),
    map_signif_level = TRUE,
    step_increase = 0.1,  # Adjusts spacing between the lines
    test = "wilcox.test"  # Perform a t-test to calculate significance
  ) + 
  theme_minimal() +
  labs(
    title = "Comparison of Robustness Measures by cent.btw",
    x = "PPI Category", 
    y = "AUC"
  ) +
  theme(legend.position = "top")
(g3)
}

###################################################
# Fig H) violin plot of median betweenness per category
# original code: Code: 11.3_CTS_cardiac_network_ANND_pagerank.R ??
# original pdf: BetweennessCentrality_GSE870383_v2.pdf
################################################################
{
df_BC = read.table(file='df_betweeness.tsv',sep='\t', header=T) 
df_BC$PPI_cat = factor(df_BC$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 			

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
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
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
# Fig I) cumulative distribution of normalized strength
# original code: 11.2.1_CTS_cardiac_network_strengthDistribution.R
# original pdf: normalized_strength_distribution.pdf
################################################################
{
V_deg_nor_dis <- lapply(graph_list, function(g) {
	  deg <- strength(g)
	  deg_table <- table(deg)
	  df <- data.frame(
		k = as.integer(names(deg_table)),
		freq = as.numeric(deg_table)
	  )
	  df$nor_strength = df$freq / sum(df$freq)
	  df$nor_strength_cum = rev(cumsum(rev(df$nor_strength)))
	  df
	}) %>%
data.table::rbindlist(idcol = "signature")  

V_deg_nor_dis$PPI_cat = lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
		factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 

(table(V_deg_nor_dis$PPI_cat))
#    CTS HiGCTS    HiG 
#     59     25   7994 
V_deg_nor_dis$cluster = lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 

(all(V_deg_nor_dis$signature %in% names(graph_list))) #TRUE
V_deg_nor_dis$n_nodes <- 0
for(i in 1:length(graph_list)){
	j = which(V_deg_nor_dis$signature == names(graph_list)[i])
	V_deg_nor_dis$n_nodes[j] <- vcount(graph_list[[i]])
	}
g_strength_dis2 <- V_deg_nor_dis %>%
		  filter(k > 0, nor_strength_cum > 0) %>% #!!!!!!!!! NEW !!!!!!

ggplot(aes(x = k, y = nor_strength_cum, color = PPI_cat, linetype = PPI_cat, size = PPI_cat)) +
		  geom_line(aes(group=signature, linetype=PPI_cat)) +
		  scale_color_manual(values = PPI_color_platte) +
		  scale_size_manual(values = PPI_size_platte) +
		  xlab('Normalized strength level') + ylab('Fraction of nodes having a normalized strength ≥ x') +
			  theme(legend.position = c(0.2, 0.75),
				legend.justification = c(1, 1),
				legend.text = element_text(size = 5)) +
		  coord_trans(x = "log10", y = "log10")
(g_strength_dis2)
}

###################################################
# Fig J) boxplot of betweenness per PPIN
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
(dim(df5_CHD))  # [1] 13  6

df_BC = rbind(subset(df_BC, PPI_cat=='CTS'),
					subset(df_BC, PPI_cat=='HiGCTS'),
					subset(df_BC, PPI_cat=='HiG')
					)
df_BC$signature <- factor(df_BC$signature, levels = unique(df_BC$signature))
pr1 <- ggplot(df_BC, aes(x = signature, y = log10(BetweennessCentrality + 1), colour = PPI_cat)) +
  geom_boxplot(position = "dodge2") +
  scale_color_manual(values = PPI_color_platte) +
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
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

(pr1)
}

###################################################
# Fig K) boxplot of pagerank per PPIN
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
df$label=df$gene
df$PCGC_AllCurated = toupper(df$gene) %in% toupper(unlist(CHD['Griffin2023_PCGC_AllCurated']))	
# Calculate top 5 significant genes within each box
df5 <- df %>%
	  filter(rank_by_PR <= 5) %>%
	  ungroup()	
df5_CHD = subset(df5, PCGC_AllCurated==TRUE)
(dim(df5_CHD))  # [1] 8 18

pr2 <- ggplot(df, aes(x = signature, y = PageRank, colour = PPI_cat)) +
  geom_boxplot(show.legend = TRUE) +
  scale_color_manual(values = PPI_color_platte) +
  geom_text(
    data = df5_CHD,
    aes(label = gene),
    size = 3.5,
    vjust = -0.5,
    nudge_y = 0.005,  # adjust depending on PageRank range
    check_overlap = TRUE
  ) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  scale_x_discrete(limits = unique(df$signature)) +
  labs(color = "PPI cat")

(pr2)
}

### Save Plots to Folder
plots <- list(
    pr,
    g2,
    g2_violin,
    violin_median_page_wilcox,
    g,
    p_attack4,
    g3,
    violin_median_bc_wilcox,
    g_strength_dis2,
    pr1,
    pr2
)

sizes <- list(
  NULL,
  NULL,
  NULL,
  c(3, 3),
  NULL,
  c(5, 5),
  c(3, 5),
  c(3, 3),
  NULL,
  c(9, 5),
  c(9, 5)
)

dir.create("plots", showWarnings = FALSE)

file_names <- paste0(LETTERS[1:length(plots)], ".pdf")

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