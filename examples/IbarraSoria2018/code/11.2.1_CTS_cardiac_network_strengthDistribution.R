library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(igraph)

########## BEGINNING OF USER INPUT ##########

wd = '/Users/felixyu/Documents/IbarraSoria2018/'
setwd(paste0(wd, 'results/'))
PPI_color_palette <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_palette <- c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)

db <- "IbarraSoria2018"

celltype_specific_weight_version <- '10'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

CT_id <- c("cardiac.a", "endothelial.b")  # critical transition clusters
CT_id_formatted <- paste0("_(", paste(CT_id, collapse = "|"), ")")

########## END OF USER INPUT ##########
 
graph_list <- readRDS( file= paste0(db, '_STRING_graph_perState_notsimplified.rds'))
(N0 = sapply(graph_list, vcount))
(N0)
# HiG_extraembryonicMesoderm          HiG_endothelial.a 
#                        414                        414 
#          HiG_endothelial.c          HiG_endothelial.d 
#                        573                        573 
#                  HiG_blood    HiG_mesodermProgenitors 
#                        559                        559 
#   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                        537                        537 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a 
#                        552                        552 
#     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                        432                        432 
#              HiG_cardiac.b              HiG_cardiac.c 
#                        549                        549 
#          HiG_endothelial.b              HiG_cardiac.a 
#                        392                        392 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a 
#                         16                         13 
#          CTS_endothelial.b              CTS_cardiac.a 
#                         33                         37 


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
# HiG_extraembryonicMesoderm          HiG_endothelial.a 
#                      13486                      13486 
#          HiG_endothelial.c          HiG_endothelial.d 
#                      30494                      30494 
#                  HiG_blood    HiG_mesodermProgenitors 
#                      23448                      23448 
#   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                      25324                      25324 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a 
#                      41204                      41204 
#     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                      15114                      15114 
#              HiG_cardiac.b              HiG_cardiac.c 
#                      25824                      25824 
#          HiG_endothelial.b              HiG_cardiac.a 
#                      13216                      13216 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a 
#                         40                         56 
#          CTS_endothelial.b              CTS_cardiac.a 
#                        164                        158 

graphs_with_duplicates <- sapply(graph_list, function(g) {
  vertex_names <- V(g)$name
  if(is.null(vertex_names)) {
    # If no names, use vertex indices
    vertex_names <- V(g)
  }
  any(duplicated(vertex_names))
})

# See which graphs have duplicates
(which(graphs_with_duplicates))
# named integer(0)

V_deg = lapply(graph_list, function(x) {igraph::strength(x, weights = E(x)$weight )/ (vcount(x) - 1)} %>% sort(., decreasing=T)) %>%
           lapply(., function(x) x %>% 
                   as.data.frame(strength=x) %>% 
                   mutate(gene=names(x), id=seq_along(x))) %>% 
          rbindlist(.,idcol=names(.))
colnames(V_deg)[1:2]=c('signature','nor_strength')
V_deg$PPI_cat = lapply(names(V_deg$signature), function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist %>%
    factor(.,levels=c('CTS', 'HiGCTS', 'HiG'))
	
  # Aggregate strength Calculation:  this calculates the average strength for each signature (graph)	
df =  aggregate(V_deg$nor_strength, by =list(V_deg$signature), FUN=mean)  %>%  
			mutate(k=aggregate(V_deg$id, by =list(V_deg$signature), FUN=max)[,2] )  %>% 
			arrange(desc(x))
df$PPI_cat = lapply(df$Group.1, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
			factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 

g_strength <- ggplot(data=df, aes(x=k, y=x, col=PPI_cat)) +
			scale_color_manual(values = PPI_color_palette) +
           geom_point(shape=18, size=5) + xlab('number of nodes per PPI_cat') +
           theme(legend.position=c(1, 1), legend.justification=c(0, 1)) +
           ylab('average GRN normalized strength') + ggtitle(db)


  
  # cumulative (normalized= FALSE!!)strength distribution to a power law fit ########################
  # PPI_cat: PPI network category: A: CTS; B: CTS&hiG; C; HiG  ############
  #### normalized = FALSE by default
  V_deg_dis = lapply(graph_list, function(x) strength_distribution(x, normalized=FALSE, cumulative=TRUE )) %>% 
              lapply(., function(x) x %>% 
                      as.data.frame(strength_distribution=x)  %>% 
                      mutate(k=seq_along(x))) %>% 
    rbindlist(.,idcol=names(.))
  
  colnames(V_deg_dis)[1:2]=c('signature','strength_distribution')
  V_deg_dis$PPI_cat = lapply(V_deg_dis$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
			factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 

  (table(V_deg_dis$PPI_cat))
#    CTS HiGCTS    HiG 
#     30     20   3526 
  V_deg_dis$cluster = lapply(V_deg_dis$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 
  
    all(V_deg_dis$signature %in% names(graph_list))
	V_deg_dis$n_nodes <- 0
	for(i in seq_along(graph_list)){
		j = which(V_deg_dis$signature == names(graph_list)[i])
		V_deg_dis$n_nodes[j] <- vcount(graph_list[[i]])
	}
 
    ## To provide insights into the distribution of node strengths in each signature and how the strength distribution varies across "transitory" and "steady" PPI_cats,
	## we plot cumulative normalized strength distribution on log scale. 
	# showing how many (the cumulative fraction of) nodes in each signature  having a strength greater than or equal to k (the strength),
	# with lines shaped based on the signature's PPI_cat ("transitory" or "steady")
	g_strength_dis <- ggplot(data= V_deg_dis %>% filter(strength_distribution > 0),
                         aes(x=k, y= strength_distribution, color= cluster, type=PPI_cat, size=PPI_cat )) +  #
		geom_line(aes(linetype=PPI_cat)) + xlab('cumulative  strength distribution') +
		#scale_color_manual(values = PPI_color_palette) +  # Set line width
		scale_size_manual(values = PPI_size_palette) +  # Set line width
		geom_text(#data=V_deg_dis_text, 
				 aes(label=n_nodes, color=cluster),  # interaction(PPI_cat, cluster)), 
				 hjust=1.1, vjust=0.5, check_overlap = TRUE, size=3)   + # Adding text for n_nodes
		theme(legend.position=c(0.2, 0.75), legend.justification=c(1, 1), legend.text = element_text(size=5)) + 
		coord_trans(x="log10", y="log10")  
	   # scale_x_reverse()  # Flip the x-axis from highest to smallest
  
    ggsave("strength_distribution_w_vsize.pdf", g_strength_dis, width = 11, height = 10)
  (n_nodes =  lapply(graph_list, vcount) %>% unlist)
  (n_nodes)
# HiG_extraembryonicMesoderm          HiG_endothelial.a 
#                        414                        414 
#          HiG_endothelial.c          HiG_endothelial.d 
#                        573                        573 
#                  HiG_blood    HiG_mesodermProgenitors 
#                        559                        559 
#   HiG_presomiticMesoderm.b   HiG_presomiticMesoderm.a 
#                        537                        537 
#        HiG_somiticMesoderm        HiG_mixedMesoderm.a 
#                        552                        552 
#     HiG_pharyngealMesoderm        HiG_mixedMesoderm.b 
#                        432                        432 
#              HiG_cardiac.b              HiG_cardiac.c 
#                        549                        549 
#          HiG_endothelial.b              HiG_cardiac.a 
#                        392                        392 
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a 
#                         16                         13 
#          CTS_endothelial.b              CTS_cardiac.a 
#                         33                         37
  
  V_deg_nor_dis <- lapply(graph_list, function(g) {
		  deg <- strength(g, weights = E(g)$weight)
		  deg_table <- table(deg)
		  df <- data.frame(
			k = as.integer(names(deg_table)), #the strength values
			freq = as.numeric(deg_table)   # how many attices have each strength value
		    )
		  df$nor_strength = df$freq / sum(df$freq) # normalized frequency (probability); sum(df$freq)= vcount(g)
		  df$nor_strength_cum = rev(cumsum(rev(df$nor_strength))) # cumulative probability distribution
		  df
		}) %>%
  data.table::rbindlist(idcol = "signature")  
  
  V_deg_nor_dis$PPI_cat = lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
			factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 

  (table(V_deg_nor_dis$PPI_cat))
#    CTS HiGCTS    HiG 
#     59     25   7968 
  V_deg_nor_dis$cluster = lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 
  
    all(V_deg_nor_dis$signature %in% names(graph_list))
	V_deg_nor_dis$n_nodes <- 0
	for(i in seq_along(graph_list)){
		j = which(V_deg_nor_dis$signature == names(graph_list)[i])
		V_deg_nor_dis$n_nodes[j] <- vcount(graph_list[[i]])
	}
 
    ## To provide insights into the distribution of node strengths in each signature and how the strength distribution varies across "transitory" and "steady" PPI_cats,
	## we plot cumulative normalized strength distribution on log scale. 
	# showing how many (the cumulative fraction of) nodes in each signature  having a strength greater than or equal to k (the strength),
	# with lines shaped based on the signature's PPI_cat ("transitory" or "steady")
	g_strength_dis <- V_deg_nor_dis %>%
		  filter(k > 0, nor_strength_cum > 0) %>% #!!!!!!!!! NEW !!!!!!
		  ggplot(aes(x = k, y = nor_strength_cum, color = cluster, linetype = PPI_cat, size=PPI_cat)) +
		  geom_line() +
		  scale_size_manual(values = PPI_size_palette) +  # Set line width
		  xlab('Normalized strength level') + ylab('cumulative normalized  strength distribution') +
		  geom_text(aes(label = n_nodes), hjust = 1.1, vjust = 0.5, check_overlap = TRUE, size = 3) +
		theme(legend.position = c(0.2, 0.75),
				legend.justification = c(1, 1),
				legend.text = element_text(size = 5)) +
		  coord_trans(x = "log10", y = "log10")
	print(g_strength_dis)
	
	g_strength_dis2 <- V_deg_nor_dis %>%
		  filter(k > 0, nor_strength_cum > 0) %>% #!!!!!!!!! NEW !!!!!!
		  ggplot(aes(x = k, y = nor_strength_cum, color = PPI_cat, linetype = PPI_cat, size=PPI_cat)) +
		  geom_line(aes(group=signature, linetype=PPI_cat)) +
		  scale_color_manual(values = PPI_color_palette) +
		  scale_size_manual(values = PPI_size_palette) +  # Set line width
		  xlab('Normalized strength level') + ylab('Fraction of nodes having a normalized strength ≥ x') +
			  theme(legend.position = c(0.2, 0.75),
				legend.justification = c(1, 1),
				legend.text = element_text(size = 5)) +
		  coord_trans(x = "log10", y = "log10")
	print(g_strength_dis2)
	  
  pdf(file='normalized_strength_distribution.pdf')
  print(g_strength)
  print(g_strength_dis) 
  print(g_strength_dis2) 
  dev.off()



 
  g_strength_dis <- ggplot(data=V_deg_dis %>% filter(strength_distribution > 0),
                         aes(x= k, y= strength_distribution, color= cluster, type=PPI_cat, size=PPI_cat )) +  #
    geom_line(aes(linetype=PPI_cat)) + xlab(' strength level (x)') +
	scale_size_manual(values = PPI_size_palette) +  # Set line width
    ylab('Fraction of nodes having a strength ≥ x') + 
    theme(legend.position=c(0.2, 0.75), legend.justification=c(1, 1), legend.text = element_text(size=5)) + 
    coord_trans(x="log10", y="log10")
 print(g_strength_dis)
 
 g_strength_dis2 <- ggplot(data=V_deg_dis %>% filter(strength_distribution > 0),
                         aes(x= k, y= strength_distribution, color= PPI_cat , size=PPI_cat)) +  #
    geom_line(aes(group=signature, linetype=PPI_cat)) + xlab(' strength level (x)') +
    scale_color_manual(values = PPI_color_palette) +
	scale_size_manual(values = PPI_size_palette) +  # Set line width
	ylab('Fraction of nodes having a strength ≥ x') + 
    theme(legend.position=c(0.2, 0.75), legend.justification=c(1, 1), legend.text = element_text(size=5)) + 
    coord_trans(x="log10", y="log10") 
 print(g_strength_dis2)
 
 pdf(file=paste0('strength_', db, '.pdf'))
  print(g_strength)
  print( g_strength_dis2)
  print(g_strength_dis ) 
  dev.off()  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  ######## Compare node-strength distribution across three categories ################
  print(dim(V_deg_dis)) # [1] 3576    6
  print(head(V_deg_dis, 3))
#                     signature strength_distribution     k PPI_cat                cluster n_nodes
#                        <char>                 <num> <int>  <fctr>                 <char>   <num>
# 1: HiG_extraembryonicMesoderm             1.0000000     1     HiG extraembryonicMesoderm     414
# 2: HiG_extraembryonicMesoderm             0.9830918     2     HiG extraembryonicMesoderm     414
# 3: HiG_extraembryonicMesoderm             0.9541063     3     HiG extraembryonicMesoderm     414


###################################################
 ## 4.2) evaluate the strength distribution of each PPI_cat; NORMALIZED the cumulative strength is WRONG!!!!!! 
  
 V_deg_dis$normalized_strength_distribution =  V_deg_dis$strength_distribution / (V_deg_dis$n_nodes-1)


	 # Density Plot / Kernel Density Estimate (KDE)
	 ggplot(V_deg_dis, aes(x = normalized_strength_distribution, fill = PPI_cat)) +
		  geom_density(alpha = 0.5) +  # Create density plot
		  labs(x = "Normalized strength Distribution", y = "Density", title = "Density Plot: strength Distribution by Width Category") +
		  theme_minimal() +
		  scale_fill_manual(values = PPI_color_palette)
	# Boxplot by Categories (NOT USED)
	 g1 = ggplot(V_deg_dis, aes(x = factor(PPI_cat), y = normalized_strength_distribution, fill = PPI_cat)) + 
		  geom_boxplot() +
		  labs(x = "PPIN Category", y = "Normalized strength Distribution", title =  "PPINs for all clusters") +
		  theme_minimal() +
		  scale_fill_manual(values = PPI_color_palette) +
		  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")), 
                     p.adjust.method = "BH",  # Adjust p-values using Benjamini-Hochberg (BH) method
                     label = "p.signif")
	 # Boxplot by Categories (NEW,  USED)
	 g2 = ggplot(subset(V_deg_dis, grepl("_endothelial.b|_cardiac.a",signature)), 
				aes(x = factor(PPI_cat), y = normalized_strength_distribution, fill = PPI_cat)) + 
		  geom_boxplot() +
		  labs(x = "PPIN Category", y = "Normalized strength Distribution", title = "PPINs for transition clusters") +
		  theme_minimal() +
		  scale_fill_manual(values = PPI_color_palette) +
		  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")), 
                     p.adjust.method = "BH",  # Adjust p-values using Benjamini-Hochberg (BH) method
                     label = "p.signif")
library(gridExtra)
pdf(file=paste0('boxplot_normalized_strength_', db, '.pdf'), height=4)	
print(grid.arrange(g1, g2, ncol = 2))
dev.off()
 
    tmp = subset(V_deg_dis, grepl(CT_id_formatted, signature))
	
	print(table(V_deg_dis$PPI_cat))
#    CTS HiGCTS    HiG 
#     30     20   3526
	
	print(table(tmp$PPI_cat))
#    CTS HiGCTS    HiG 
#     30     20    350
    saveRDS(V_deg_dis, file='V_deg_dis.rds')  ## note that 'PPI_cast' was originally named as 'width'
	
	
	
  