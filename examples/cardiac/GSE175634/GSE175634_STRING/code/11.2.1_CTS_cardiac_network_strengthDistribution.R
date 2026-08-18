# this version simplify duplicated edges between any two vertices 
# Q1: Does the GRN of 37 CTS genes in E8.25_2018 dataset fragile?
# # to build PPI_cat-specific GRN for each PPI_cat 
# identify PPI_cat-specially highly expressed genes for steady PPI_cat but using the CTS for transitory PPI_cats

# Q2: Does the GRN of 19 CTS genes in hESC dataset fragile? 
# -- using string physicial links
# -- using HumanBase 'cardiac muscle' links
# to build PPI_cat-specific GRN for C5, C9 (CT), C7, C8, and C10 (CT) 
# 1-3) identify PPI_cat-specially highly expressed genes for C5, C7, and C8, and using the identified CTS for C9, C10
# 4) evaluate the GRN robustness (stability) of each PPI_cat
# 4.1) estimate the significance of strength distribution
# 4.2) estimate the significane of robustness
# 5) Does the most individual contributors matter 
#  ('specialized processes (e.g., cell differentiation) are mainly related to regulators with low Knn,
#    whereas essential processes are mainly related to regulators with high page rank or strength')
# 6) Does the GRN of  transitory PPI_cat exhibit bisignatures but not the steady PPI_cat?
# 7) Will the GRN of tranistory PPI_cat (smalle number of nodes) be compareable to the Subset the GRN of steady PPI_cat with the same number of nodes?


## compare the cardiac CTS from hESC:C10, E2018 2018 cardiac.a, and E8.25 2019 C8 (not the in vivo and EB C4)
## compare HEP CTS from E8.25 2019 C13, E8.25 2018 endo.b, EB C11,
# and E8.25 2019 C17 (or C8)
library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr) # resuired by stat_compare_means()
library(igraph)

setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weight_shrink')
PPI_color_platte = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_platte = c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)
 
 
graph_list <- readRDS( file= 'GSE175634_STRING_graph_perState_notsimplified.rds')  
(N0 = sapply(graph_list, vcount))
          # HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4           HiG_5 
            # 325            1156             246             132             105             339 
          # HiG_6           HiG_7      HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
            # 108             108             118            1468              70              42 
         # HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP     HiGCTS_CP.1      CTS_muscle 
             # 22              20               7              29              27              61 
   # CTS_endoderm          CTS_CP        CTS_CP.1 
             # 64              72              77 

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

graph_list <- lapply(graph_list, simplify) #!!!!!!!!!!!!!!!!!!!
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

# norm_strength <- strength(g, weights = E(g)$weight) / (vcount(g) - 1) 
V_deg = lapply(graph_list, function(x) {igraph::strength(x, weights = E(x)$weight )/ (vcount(x) - 1)} %>% sort(., decreasing=T)) %>%
           lapply(., function(x) x %>% 
                   as.data.frame(strength=x) %>% 
                   mutate(gene=names(x), id=1:length(x))) %>% 
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
			scale_color_manual(values = PPI_color_platte) +
           geom_point(shape=18, size=5) + xlab('number of nodes per PPI_cat') +
           theme(legend.position=c(1, 1), legend.justification=c(0, 1)) +
           ylab('average GRN normalized strength') + ggtitle('GSE175634')


  
  # cumulative (normalized= FALSE!!)strength distribution to a power law fit ########################
  # PPI_cat: PPI network category: A: CTS; B: CTS&hiG; C; HiG  ############
  #### normalized = FALSE by default
source('F:/projects/scRNA/source/cardiac_CTS_GRN/celltype_specific_weight_v3.R')
  V_deg_dis = lapply(graph_list, function(x) strength_distribution(x, normalized=FALSE, cumulative=TRUE )) %>% 
              lapply(., function(x) x %>% 
                      as.data.frame(strength_distribution=x)  %>% 
                      mutate(k=1:length(x))) %>% 
    rbindlist(.,idcol=names(.))
  
  colnames(V_deg_dis)[1:2]=c('signature','strength_distribution')
  V_deg_dis$PPI_cat = lapply(V_deg_dis$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
			factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 

  table(V_deg_dis$PPI_cat)
     # CTS HiGCTS    HiG 
    # 109     45   1951
  V_deg_dis$cluster = lapply(V_deg_dis$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 
  
    all(V_deg_dis$signature %in% names(graph_list))
	V_deg_dis$n_nodes <- 0
	for(i in 1:length(graph_list)){
		j = which(V_deg_dis$signature == names(graph_list)[i])
		V_deg_dis$n_nodes[j] <- vcount(graph_list[[i]])
		}
 
    ## To rovide insights into the distribution of node strengths in each signature and how the strength distribution varies across "transitory" and "steady" PPI_cats,
	## we plot cumulative normalized strength distribution on log scale. 
	# showing how many (the cumulative fraction of) nodes in each signature  having a strength greater than or equal to k (the strength),
	# with lines shaped based on the signature's PPI_cat ("transitory" or "steady")
	g_strength_dis <- ggplot(data= V_deg_dis %>% filter(strength_distribution > 0),
                         aes(x=k, y= strength_distribution, color= cluster, type=PPI_cat, size=PPI_cat )) +  #
		geom_line(aes(linetype=PPI_cat)) + xlab('cumulative  strength distribution') +
		#scale_color_manual(values = PPI_color_platte) +  # Set line width
		scale_size_manual(values = PPI_size_platte) +  # Set line width
		geom_text(#data=V_deg_dis_text, 
				 aes(label=n_nodes, color=cluster),  # interaction(PPI_cat, cluster)), 
				 hjust=1.1, vjust=0.5, check_overlap = TRUE, size=3)   + # Adding text for n_nodes
		theme(legend.position=c(0.2, 0.75), legend.justification=c(1, 1), legend.text = element_text(size=5)) + 
		coord_trans(x="log10", y="log10")  
	   # scale_x_reverse()  # Flip the x-axis from highest to smallest
  
  print(g_strength_dis) 
  dev.copy2pdf(file='strength_distribution_w_vsize.pdf')
  (n_nodes =  lapply(graph_list, vcount) %>% unlist)
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


 # cumulative NORMALIZED strength distribution to a power law fit ########################
  # PPI_cat: PPI network category: A: CTS; B: CTS&hiG; C; HiG  ############
  #### normalized = FALSE by default
    # V_deg_nor_dis = lapply(graph_list_nonempty, function(x) igraph::strength_distribution(x, cumulative=TRUE, normalized=TRUE )) %>% 
              # lapply(., function(x) x %>% 
                      # as.data.frame(strength_distribution=x)  %>% 
                      # mutate(k=1:length(x))) %>% 
    # rbindlist(.,idcol=names(.))
  #graph_list_nonempty <- Filter(function(g) ecount(g) > 0, graph_list)
  
  V_deg_nor_dis <- lapply(graph_list, function(g) {
		  deg <- strength(g, weights = E(g)$weight)
		  deg_table <- table(deg)
		  df <- data.frame(
			k = as.integer(names(deg_table)), #the strength values
			freq = as.numeric(deg_table)   # how many vertices have each strength value
		    )
		  df$nor_strength = df$freq / sum(df$freq) # normalized frequency (probability); sum(df$freq)= vcount(g)
		  df$nor_strength_cum = rev(cumsum(rev(df$nor_strength))) # cumulative probability distribution
		  df
		}) %>%
  data.table::rbindlist(idcol = "signature")  
  
  #colnames(V_deg_nor_dis)[1:2]=c('signature','nor_strength_distribution')
  V_deg_nor_dis$PPI_cat = lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
			factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 

  table(V_deg_nor_dis$PPI_cat)
     # CTS HiGCTS    HiG 
    #  227     66   4160
  V_deg_nor_dis$cluster = lapply(V_deg_nor_dis$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 
  
    all(V_deg_nor_dis$signature %in% names(graph_list))
	V_deg_nor_dis$n_nodes <- 0
	for(i in 1:length(graph_list)){
		j = which(V_deg_nor_dis$signature == names(graph_list)[i])
		V_deg_nor_dis$n_nodes[j] <- vcount(graph_list[[i]])
		}
 
    ## To rovide insights into the distribution of node strengths in each signature and how the strength distribution varies across "transitory" and "steady" PPI_cats,
	## we plot cumulative normalized strength distribution on log scale. 
	# showing how many (the cumulative fraction of) nodes in each signature  having a strength greater than or equal to k (the strength),
	# with lines shaped based on the signature's PPI_cat ("transitory" or "steady")
	g_strength_dis <- V_deg_nor_dis %>%
		  filter(k > 0, nor_strength_cum > 0) %>% #!!!!!!!!! NEW !!!!!!
		  ggplot(aes(x = k, y = nor_strength_cum, color = cluster, linetype = PPI_cat, size=PPI_cat)) +
		  geom_line() +
		  scale_size_manual(values = PPI_size_platte) +  # Set line width
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
		  scale_color_manual(values = PPI_color_platte) +
		  scale_size_manual(values = PPI_size_platte) +  # Set line width
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



 
   
  #dev.new(width=14)
  g_strength_dis <- ggplot(data=V_deg_dis %>% filter(strength_distribution > 0),
                         aes(x= k, y= strength_distribution, color= cluster, type=PPI_cat, size=PPI_cat )) +  #
    geom_line(aes(linetype=PPI_cat)) + xlab(' strength level (x)') +
	scale_size_manual(values = PPI_size_platte) +  # Set line width
    ylab('Fraction of nodes having a strength ≥ x') + 
    theme(legend.position=c(0.2, 0.75), legend.justification=c(1, 1), legend.text = element_text(size=5)) + 
    coord_trans(x="log10", y="log10") #+ #xlab('strength k') +
	#scale_x_reverse()  # Flip the x-axis from highest to smallest  NOT WORK!!! manually flip
 print(g_strength_dis)
 
 g_strength_dis2 <- ggplot(data=V_deg_dis %>% filter(strength_distribution > 0),
                         aes(x= k, y= strength_distribution, color= PPI_cat , size=PPI_cat)) +  #
    geom_line(aes(group=signature, linetype=PPI_cat)) + xlab(' strength level (x)') +
    scale_color_manual(values = PPI_color_platte) +
	scale_size_manual(values = PPI_size_platte) +  # Set line width
	ylab('Fraction of nodes having a strength ≥ x') + 
    theme(legend.position=c(0.2, 0.75), legend.justification=c(1, 1), legend.text = element_text(size=5)) + 
    coord_trans(x="log10", y="log10") 
 print(g_strength_dis2)
 
 pdf(file='strength_GSE175634.pdf')
  print(g_strength)
  print( g_strength_dis2)
  print(g_strength_dis ) 
  dev.off()  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  ######## Compare node-strength distribution across three categories ################
  dim(V_deg_dis) # [1] 2241    7
  head(V_deg_dis, 3)
   # signature strength_distribution     k PPI_cat cluster n_nodes
      # <char>                 <num> <int>  <fctr>  <char>   <num>
# 1:     HiG_0             1.0000000     1     HiG       0     325
# 2:     HiG_0             1.0000000     2     HiG       0     325
# 3:     HiG_0             0.9876923     3     HiG       0     325


###################################################
 ## 4.2) evaluate the strength distribution of each PPI_cat; NORMALIZED the cumulative strength is WRONG!!!!!! 
  
  # # cumulative NORMALIZED strength distribution to a power law fit ########################
  # # If TRUE then the result is divided by n−1, where n is the number of vertices in the graph.
 V_deg_dis$normalized_strength_distribution =  V_deg_dis$strength_distribution / (V_deg_dis$n_nodes-1)
  # #dev.new(width=14)
  # g_strength_dis <- ggplot(data=V_deg_dis,
                         # aes(x= k, y= normalized_strength_distribution, color= cluster, type=PPI_cat )) +  #
   # # #                  aes(x=strength_dis, color= PPI_cat, fill=PPI_cat)) +
   # # # geom_histogram(aes(y=..density..), position="identity" ,alpha=0.4) +
   # # # geom_density(alpha=0.3 ) + xlab('strength distribution') + 
    # geom_line(aes(linetype=PPI_cat)) + xlab('strength level (x)') +
    # ylab('Fraction of nodes having a normalized strength ≥ x') + 
    # theme(legend.position=c(0.2, 0.75), legend.justification=c(1, 1), legend.text = element_text(size=5)) + 
    # coord_trans(x="log10", y="log10") #+ #xlab('strength k') +
	# #scale_x_reverse()  # Flip the x-axis from highest to smallest  NOT WORK!!! manually flip
 
 # pdf(file='normalized_strength_GSE175634.pdf')
  # print(g_strength_dis ) 
  # dev.off()  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  ######## Compare node-strength distribution across three categories ################
  dim(V_deg_dis) # [1] 2105    7
  head(V_deg_dis, 3)
   # signature strength_distribution     k PPI_cat cluster n_nodes normalized_strength_distribution
      # <char>                 <num> <int>  <fctr>  <char>   <num>                            <num>
# 1:     HiG_0             1.0000000     1     HiG       0     325                      0.003086420
# 2:     HiG_0             1.0000000     2     HiG       0     325                      0.003086420
# 3:     HiG_0             0.9876923     3     HiG       0     325                      0.003048433

	 # Density Plot / Kernel Density Estimate (KDE)
	 ggplot(V_deg_dis, aes(x = normalized_strength_distribution, fill = PPI_cat)) +
		  geom_density(alpha = 0.5) +  # Create density plot
		  labs(x = "Normalized strength Distribution", y = "Density", title = "Density Plot: strength Distribution by Width Category") +
		  theme_minimal() +
		  scale_fill_manual(values = PPI_color_platte)
	# Boxplot by Categories (NOT USED)
	 g1 = ggplot(V_deg_dis, aes(x = factor(PPI_cat), y = normalized_strength_distribution, fill = PPI_cat)) + 
		  geom_boxplot() +
		  labs(x = "PPIN Category", y = "Normalized strength Distribution", title =  "PPINs for all clusters") +
		  theme_minimal() +
		  scale_fill_manual(values = PPI_color_platte) +
		  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")), 
                     p.adjust.method = "BH",  # Adjust p-values using Benjamini-Hochberg (BH) method
                     label = "p.signif")
	 # Boxplot by Categories (NEW,  USED)
	 g2 = ggplot(subset(V_deg_dis, grepl("_CP|_muscle|_endoderm",signature)), 
				aes(x = factor(PPI_cat), y = normalized_strength_distribution, fill = PPI_cat)) + 
		  geom_boxplot() +
		  labs(x = "PPIN Category", y = "Normalized strength Distribution", title = "PPINs for 3 transition clusters") +
		  theme_minimal() +
		  scale_fill_manual(values = PPI_color_platte) +
		  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CTS", "HiGCTS"), c("CTS", "HiG"), c("HiGCTS", "HiG")), 
                     p.adjust.method = "BH",  # Adjust p-values using Benjamini-Hochberg (BH) method
                     label = "p.signif")
library(gridExtra)
pdf(file='boxplot_normalized_strength_GSE175634.pdf', height=4)	
print(grid.arrange(g1, g2, ncol = 2))
dev.off()
 
    tmp = subset(V_deg_dis, grepl("_CP|_muscle|_endoderm",signature))
	
	table(V_deg_dis$PPI_cat)
	 # CTS HiGCTS    HiG 
   # 109     45   1951 
	
	table(tmp$PPI_cat)
	# CTS HiGCTS    HiG 
   # 109     45    176 
    saveRDS(V_deg_dis, file='V_deg_dis.rds')  ## note taht 'PPI_cast' was originally named as 'width'
	
	
	
  