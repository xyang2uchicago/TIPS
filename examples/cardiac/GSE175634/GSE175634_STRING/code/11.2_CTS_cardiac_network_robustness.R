
library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr) # resuired by stat_compare_means()

library(igraph)
library(rstatix)
 
library(brainGraph)

# source('F:/projects/scRNA/source/cardiac_CTS_GRN/celltype_specific_weight_v6.R')
source('E:/Git_Holly/TIPS-main/celltype_specific_weight_v10.R')

PPI_color_platte = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_platte = c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)

setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9')

 # refer to 11.1_CTS_cardiac_network_strengthDistribution.R 
s = "combined"
file = paste0('GSE175634_STRING_graph_perState_simplified_',s,'weighted.rds')
graph_list <- readRDS( file)  
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

# Check which graphs have duplicate vertex names
graphs_with_duplicates <- sapply(graph_list, function(g) {
  vertex_names <- V(g)$name
  if(is.null(vertex_names)) {
    # If no names, use vertex indices
    vertex_names <- V(g)
  }
  any(duplicated(vertex_names))
})

# See which graphs have duplicates
which(graphs_with_duplicates) #named integer(0)

###################################################
# See if weights have been updated
# combined_score: the 
names(edge_attr(graph_list[[1]]))
# [1] "combined_score"  "weight"          "width"           "original_weight" "corexp_sign"    
# [6] "coexp_target"  
all(E(graph_list[[1]])$combined_score/1000 == E(graph_list[[1]])$original_weight)  # TRUE
all(E(graph_list[[1]])$weight == E(graph_list[[1]])$original_weight)  # FALSE

###################################################
### ensure shrinked and non-shrined weights are different  #######
tmp = readRDS(file=paste0('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weight_PPImax/',file))
names(edge_attr(tmp[[1]]))
# [1] "combined_score"  "weight"          "width"           "original_weight" "corexp_sign"    
# [6] "coexp_target"  
all(E(graph_list[[1]])$original_weight == E(tmp[[1]])$original_weight)  # TRUE  
all(E(graph_list[[1]])$weight == E(tmp[[1]])$weight)  # FALSE  
all(E(graph_list[[1]])$coexp_target == E(tmp[[1]])$coexp_target)  # FALSE  
all(E(tmp[[1]])$weight == E(tmp[[1]])$original_weight)  # FALSE

range(E(graph_list[[1]])$weight - E(tmp[[1]])$weight)
#[1] -3.803734e+05 -2.676686e+00
range(E(graph_list[[1]])$coexp_target - E(tmp[[1]])$coexp_target)
# [1] -4.251198e-03 -1.149097e-07 
rm(tmp)


###################################################
 ## 4) evaluate the GRN robustness (stability) of each PPI_cat; 
 ## two trials: node-normalized or not  
  # https://cwatson.github.io/files/brainGraph_UserGuide.pdf
################################################################
 library(brainGraph)
 library(data.table)
 library(ggplot2)
  
 # refer to 11.1.1_CTS_cardiac_network_strengthDistribution.R 
 # V_deg_dis = readRDS('V_deg_dis.rds')
		 
  
  
 # Targeted attack ###########################
	set.seed(2020)
 # Maximal component size as a function of vertices removed.
	attack.vertex.strength <- rbindlist(lapply(graph_list, robustness_MonteCarlo, 'vertex', 'degree'), idcol=names(graph_list))
	attack.vertex.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, 'vertex', 'btwn.cent'), idcol=names(graph_list))  # used !!
	attack.vertex.strength %>% head
	attack.vertex.btwn %>% head

	attack.vertex <- rbind(attack.vertex.strength, attack.vertex.btwn)  
	colnames(attack.vertex)[1] = 'signature'
	attack.vertex$PPI_cat= lapply(attack.vertex$signature, function(x) unlist(strsplit(x, split='_'))[1] ) %>% unlist %>% 
				  factor(.,levels=c('CTS', 'HiGCTS', 'HiG'))
	head(attack.vertex, 3)
    attack.vertex$cluster= lapply(attack.vertex$signature, function(x) unlist(strsplit(x, split='_'))[2] ) %>% unlist 
   
   
   
   #mylabels <- sub('\\.', ', ', attack.vertex[, signature])
   ## plot attack measured by btwn.cent
    g_attack <- ggplot(data= subset(attack.vertex, measure == 'btwn.cent'),
                    aes(x=removed.pct, y=comp.pct, color=PPI_cat, size=PPI_cat, #interaction(PPI_cat, signature),
                    linetype=PPI_cat # masking this line if highlighting one subtype when color=signature
					#shape=PPI_cat  #  
                    )) +
     # geom_line(size=ifelse(subset(attack.vertex ,measure == 'strength')$signature=='endothelial.a', 2, 1)) + # uisng this line if highlighting one subtype
	 geom_line(aes(group=signature)) +  # group ensures each line is drawn independently
     scale_color_manual(values = PPI_color_platte) +
	 scale_size_manual(values = PPI_size_platte) +
     geom_abline(slope=-1, intercept=1, col='gray', lty=2) +
     theme(legend.position="inside", legend.position.inside = c(1, 1), legend.justification=c(1, 1)) + 
     ylab('the remaining maximal component size / the initial maximal component size') +
     ggtitle('vertex robustness by betweenness centrality')
   print(g_attack)
   
    gsignature_attack <- ggplot(data= subset(attack.vertex, measure == 'btwn.cent'),  #!!!!!!!!!!!!!
                    aes(x=removed.pct, y=comp.pct, color=cluster, size=PPI_cat, #interaction(PPI_cat, signature),
                    linetype=PPI_cat # masking this line if highlighting one subtype when color=signature
					#shape=PPI_cat  #  
                    )) +
     # geom_line(size=ifelse(subset(attack.vertex ,measure == 'strength')$signature=='endothelial.a', 2, 1)) + # uisng this line if highlighting one subtype
	 geom_line() +   
     #scale_color_manual(values = PPI_color_platte) +
     geom_abline(slope=-1, intercept=1, col='gray', lty=2) +
     theme(legend.position="inside", legend.position.inside = c(1, 1), legend.justification=c(1, 1)) + 
	 scale_size_manual(values = PPI_size_platte) +
     ylab('the remaining maximal component size / the initial maximal component size') +
     ggtitle('vertex robustness by betweenness centrality')
   print(gsignature_attack) 
  
   ## plot attack measured by strength
   g_attack2 <- ggplot(data= subset(attack.vertex ,measure == 'degree'),
                    aes(x=removed.pct, y=comp.pct, col=PPI_cat, size=PPI_cat,
                        linetype=PPI_cat # masking this line if highlighting one subtype when color=signature
                        #shape=PPI_cat 
						)) +
     # geom_line(size=ifelse(subset(attack.vertex ,measure == 'strength')$signature=='endothelial.a', 2, 1)) + # uisng this line if highlighting one subtype
     geom_line(aes(group=signature)) +   # masking this line if highlighting one subtype
     scale_color_manual(values = PPI_color_platte) +
	 scale_size_manual(values = PPI_size_platte) +
     geom_abline(slope=-1, intercept=1, col='gray', lty=2) +
     theme(legend.position="inside", legend.position.inside = c(1, 1), legend.justification=c(1, 1)) +
     ggtitle('vertex robustness measured by strength')
   print(g_attack2) 

   gsignature_attack2 <- ggplot(data= subset(attack.vertex ,measure == 'degree'),
                    aes(x=removed.pct, y=comp.pct, col=cluster, size=PPI_cat,
                        linetype=PPI_cat # masking this line if highlighting one subtype when color=signature
                        #shape=PPI_cat 
						)) +
     # geom_line(size=ifelse(subset(attack.vertex ,measure == 'strength')$signature=='endothelial.a', 2, 1)) + # uisng this line if highlighting one subtype
     geom_line() +    
     #scale_color_manual(values = PPI_color_platte) +
	 scale_size_manual(values = PPI_size_platte) +
     geom_abline(slope=-1, intercept=1, col='gray', lty=2) +
     theme(legend.position="inside", legend.position.inside = c(1, 1), legend.justification=c(1, 1)) +
     ggtitle('vertex robustness measured by strength')
   print(gsignature_attack2) 
   
   pdf(file='vertex_attack_GSE175634.pdf') 
   print(gridExtra::grid.arrange(g_attack ,  g_attack2 , 
					gsignature_attack ,  gsignature_attack2 , 
                      ncol = 2, nrow = 2))
   #plot(attack + attack2) 
   dev.off()
   
     
   ############## to restart here anytime ################################# 
    library(brainGraph)
    require(dplyr)
    library(data.table)
   # library(ggplot2)
   # 
   # # refer to 11.1.1_CTS_cardiac_network_strengthDistribution.R 
   # file = "GSE175634_STRING_graph_perState_simplified_combinedweighted.rds"
   # graph_list <- readRDS( file )  
  
   names(graph_list)
 # # [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"           "HiG_5"          
 # # [7] "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"   
# # [13] "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"     
# # [19] "CTS_endoderm"    "CTS_CP"          "CTS_CP.1"    

   
   	attack.vertex.btwn <- rbindlist(lapply(graph_list, robustness_MonteCarlo, 'vertex', 'btwn.cent'), idcol=names(graph_list))
	colnames(attack.vertex.btwn)[1] = 'signature'
	head(attack.vertex.btwn, 3)
   # signature                   type   measure comp.size  comp.pct removed.pct
      # <char>                 <char>    <char>     <num>     <num>       <num>
# 1:     HiG_0 Targeted vertex attack btwn.cent       325 1.0000000 0.000000000
# 2:     HiG_0 Targeted vertex attack btwn.cent       324 0.9969231 0.003076923
# 3:     HiG_0 Targeted vertex attack btwn.cent       323 0.9938462 0.006153846

	dim(attack.vertex.btwn) #  4604    6
   saveRDS(attack.vertex.btwn, file='attack.vertex.btwn.rds')  #!!!
	
   # calculate the distance-based weight and saved into a tmp_list
   tmp_list = graph_list  #!!!!
   for (i in 1:length(tmp_list)) {
    g = tmp_list[[i]]
    E(g)$weight = 1/E(g)$weight
    tmp_list[[i]] = g
   }
    ## then do the edge-attack analysis
   attack.edge.btwn <- rbindlist(lapply(tmp_list, robustness_MonteCarlo, 'edge', 'btwn.cent'), idcol=names(tmp_list))  #!!!!!!!!!!!!!
   colnames(attack.edge.btwn)[1] = 'signature'
   dim(attack.edge.btwn) #  106455     6
   head(attack.edge.btwn, 3)
   # signature                 type   measure comp.size comp.pct  removed.pct
      # <char>               <char>    <char>     <num>    <num>        <num>
# 1:     HiG_0 Targeted edge attack btwn.cent       325        1 0.000000e+00
# 2:     HiG_0 Targeted edge attack btwn.cent       325        1 9.934433e-05
# 3:     HiG_0 Targeted edge attack btwn.cent       325        1 1.986887e-04

  table(attack.edge.btwn$signature)
         # CTS_CP        CTS_CP.1    CTS_endoderm      CTS_muscle           HiG_0           HiG_1 
            # 129             160             391             317           10067           31170 
         # HiG_10          HiG_12           HiG_2           HiG_4           HiG_5           HiG_6 
           # 1214              71            7130            2270            4639             895 
          # HiG_7           HiG_9          HiG_CP    HiG_endoderm      HiG_muscle       HiGCTS_CP 
           # 1588           42813            1307             142            1990              39 
    # HiGCTS_CP.1 HiGCTS_endoderm   HiGCTS_muscle 
             # 42              12              69 


	saveRDS(attack.edge.btwn, file='attack.edge.btwn.rds')  #!!!
 
  # In a random failure analysis, you choose a vertex/edge at random, remove it, and calculate the maximum
  #  component size until all elements have been removed.  
  ####################################################
  # refer to 11.2_CTS_cardiac_network_robustness_Miudway3.simulation.R
  ## run on Midway3 following is teh log. DO NOT REPEAT !!!)
# scp -p -r  xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weight/failure*.rds  F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weight_shrink/.
	sapply(graph_list, vcount)
         # HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4           HiG_5 
            # 325            1155             242             132             104             339 
          # HiG_6           HiG_7      HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
            # 106             104             118            1467              70              42 
         # HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP     HiGCTS_CP.1      CTS_muscle 
             # 22              20               7              29              27              61 
   # CTS_endoderm          CTS_CP        CTS_CP.1 
             # 64              72              77 

	attack.edge.btwn = readRDS(file='attack.edge.btwn.rds')
	attack.vertex.btwn = readRDS( file='attack.vertex.btwn.rds')
    failure.vertex = readRDS(paste0('failure.vertex_100_simplified_',s,'weighted.rds') )
	table(attack.vertex.btwn$signature, attack.vertex.btwn$type)                
                  # Targeted vertex attack
  # CTS_CP                              73
  # CTS_CP.1                            78
  # CTS_endoderm                        65
  # CTS_muscle                          62
  # HiG_0                              326
  # HiG_1                             1156
  # HiG_10                              71
  # HiG_12                              23
  # HiG_2                              243
  # HiG_4                              105
  # HiG_5                              340
  # HiG_6                              107
  # HiG_7                              105
  # HiG_9                             1468
  # HiG_CP                             133
  # HiG_endoderm                        43
  # HiG_muscle                         119
  # HiGCTS_CP                           30
  # HiGCTS_CP.1                         28
  # HiGCTS_endoderm                      8
  # HiGCTS_muscle                       21

 subset(attack.vertex.btwn, signature=='HiGCTS_endoderm')
         # signature                   type   measure comp.size  comp.pct removed.pct
            # <char>                 <char>    <char>     <num>     <num>       <num>
# 1: HiGCTS_endoderm Targeted vertex attack btwn.cent         6 1.0000000   0.0000000
# 2: HiGCTS_endoderm Targeted vertex attack btwn.cent         4 0.6666667   0.1428571
# 3: HiGCTS_endoderm Targeted vertex attack btwn.cent         3 0.5000000   0.2857143
# 4: HiGCTS_endoderm Targeted vertex attack btwn.cent         2 0.3333333   0.4285714
# 5: HiGCTS_endoderm Targeted vertex attack btwn.cent         2 0.3333333   0.5714286
# 6: HiGCTS_endoderm Targeted vertex attack btwn.cent         1 0.1666667   0.7142857
# 7: HiGCTS_endoderm Targeted vertex attack btwn.cent         1 0.1666667   0.8571429
# 8: HiGCTS_endoderm Targeted vertex attack btwn.cent         0 0.0000000   1.0000000


   failure.edge = readRDS(paste0('failure.edge_100_simplified_',s,'weighted.rds')) 
   failure.dt <- rbind(failure.edge, failure.vertex)   
   #failure.dt <- failure.vertex
   head(failure.dt,3)
    # HiG_0                type measure comp.size comp.pct  removed.pct
   # <char>              <char>  <char>     <num>    <num>        <num>
# 1:  HiG_0 Random edge removal  random       325        1 0.000000e+00
# 2:  HiG_0 Random edge removal  random       325        1 9.934433e-05
# 3:  HiG_0 Random edge removal  random       325        1 1.986887e-04

   colnames(failure.dt)[1] ='signature'
	table(failure.dt$signature, failure.dt$type)                
                  # Random edge removal Random vertex removal
  # CTS_CP                          129                    73
  # CTS_CP.1                        160                    78
  # CTS_endoderm                    391                    65
  # CTS_muscle                      317                    62
  # HiG_0                         10067                   326
  # HiG_1                         31170                  1156
  # HiG_10                         1214                    71
  # HiG_12                           71                    23
  # HiG_2                          7130                   243
  # HiG_4                          2270                   105
  # HiG_5                          4639                   340
  # HiG_6                           895                   107
  # HiG_7                          1588                   105
  # HiG_9                         42813                  1468
  # HiG_CP                         1307                   133
  # HiG_endoderm                    142                    43
  # HiG_muscle                     1990                   119
  # HiGCTS_CP                        39                    30
  # HiGCTS_CP.1                      42                    28
  # HiGCTS_endoderm                  12                     8
  # HiGCTS_muscle                    69                    21


   colnames(attack.vertex.btwn)[1] = 'signature'
   dim(failure.dt)  #  111059     6
   
   robustness.dt <- rbind(failure.dt, attack.vertex.btwn, attack.edge.btwn)  
   dim(robustness.dt)  #[1] 222118      6
   robustness.dt$PPI_cat = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
			factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 
    head(robustness.dt, 3)
  # signature                type measure comp.size comp.pct  removed.pct PPI_cat
      # <char>              <char>  <char>     <num>    <num>        <num>  <fctr>
# 1:     HiG_0 Random edge removal  random       325        1 0.000000e+00     HiG
# 2:     HiG_0 Random edge removal  random       325        1 9.934433e-05     HiG
# 3:     HiG_0 Random edge removal  random       325        1 1.986887e-04     HiG

    robustness.dt$experiment = ifelse(grepl('edge', robustness.dt$type), 'edge', 'vertex')
	robustness.dt$measure= factor(robustness.dt$measure, levels = c("random" ,  "btwn.cent"))

    table(robustness.dt$type, robustness.dt$measure)
                      # random btwn.cent
  # Random edge removal    106455         0
  # Random vertex removal    4604         0
  # Targeted edge attack        0    106455
  # Targeted vertex attack      0      4604


   robustness.dt$type = factor(robustness.dt$type,
			levels = c("Random edge removal" , "Targeted edge attack"  ,  "Random vertex removal" , "Targeted vertex attack") )
  
   robustness.dt$cluster = lapply(robustness.dt$signature, function(x) unlist(strsplit(x , '_'))[2]) %>% unlist 
                               
    p_attack3 <- ggplot(robustness.dt,
                  aes(x=removed.pct, y=comp.pct, col=cluster, linetype=PPI_cat, size=PPI_cat)) +  # colored by cluster
				geom_line(show.legend = FALSE) +
				scale_size_manual(values = PPI_size_platte) +
				facet_wrap(~ type) + 
				geom_abline(slope=-1, intercept=1, col='gray', lty=2) +
				 theme(legend.position=c(0, 0), legend.justification=c(0, 0)) +
				 labs(x='% edges/vetex removed', y='% of max. component remaining')
	print(p_attack3)		   
	p_attack4 <- ggplot(robustness.dt,
                  aes(x=removed.pct, y=comp.pct, col=PPI_cat, linetype=PPI_cat, size=PPI_cat)) +  # colored by PPI_cat
				geom_line(aes(group=signature), show.legend = FALSE) +
				scale_color_manual(values = PPI_color_platte) +
				scale_size_manual(values = PPI_size_platte) +
				facet_wrap(~ type) + 
				geom_abline(slope=-1, intercept=1, col='gray', lty=2) +
				 theme(legend.position=c(0, 0), legend.justification=c(0, 0)) +
				 labs(x='% edges/vetex removed', y='% of max. component remaining')
	print(p_attack4)
    pdf(file='attack_GSE175634.pdf') 
	plot(p_attack4) 
	plot(p_attack3) 
    dev.off()  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

 ### comapre targeted attack vs random removal across categories  ################
#########################################################

# First, Run Shapiro-Wilk test, ensuring that we have enough non-NA observations per group
normality_tests <- robustness.dt %>%
  group_by(experiment, PPI_cat, measure) %>%
  filter(sum(!is.na(comp.pct)) >= 3) %>%  # Ensure there are at least 3 non-NA values
  summarise(
    shapiro_test = list(
      tryCatch(
        expr = shapiro.test(comp.pct),
        error = function(e) NA  # If the test fails, return NA instead
      )
    ), 
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(shapiro_test, function(x) if (is.list(x)) x$p.value else NA)  # Safely extract p-value
  ) %>%
  as.data.frame()

# Print the normality test results
print(normality_tests)
#    experiment PPI_cat   measure shapiro_test      p_value
# 1        edge     CTS    random c(W = 0..... 3.254787e-23
# 2        edge     CTS btwn.cent c(W = 0..... 3.112271e-17
# 3        edge  HiGCTS    random c(W = 0..... 3.670656e-07
# 4        edge  HiGCTS btwn.cent c(W = 0..... 4.018832e-08
# 5        edge     HiG    random           NA           NA
# 6        edge     HiG btwn.cent           NA           NA
# 7      vertex     CTS    random c(W = 0..... 4.595903e-10
# 8      vertex     CTS btwn.cent c(W = 0..... 9.181926e-18
# 9      vertex  HiGCTS    random c(W = 0..... 2.834539e-04
# 10     vertex  HiGCTS btwn.cent c(W = 0..... 4.835290e-09
# 11     vertex     HiG    random c(W = 0..... 9.933774e-35
# 12     vertex     HiG btwn.cent c(W = 0..... 1.094213e-36



## then, Plot the Wilcox results for visualization and manually add fold chagnes !!!
g = ggplot(robustness.dt, aes(x = type, y = comp.pct, fill = measure, color = PPI_cat)) + 
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) +  # Dodge the boxes for each type per PPI_cat
  facet_wrap(experiment ~ PPI_cat, ncol = 3) +  # Facet by PPI_cat, each PPI_cat gets a row
  scale_color_manual(values = PPI_color_platte) +
  # scale_size_manual(values = PPI_size_platte) +
  scale_fill_manual(values = c("random" = "grey", "btwn.cent" = "white")) +
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
  theme(legend.position = "top")
print(g)

# Filter the data for 'Targeted vertex attack'
targeted_vertex_data <- robustness.dt[robustness.dt$type == 'Targeted vertex attack', ]
targeted_vertex_data <- targeted_vertex_data[!is.na(targeted_vertex_data$comp.pct), ]

g2 = ggplot(targeted_vertex_data, aes(x = PPI_cat, y = comp.pct, fill = PPI_cat)) + 
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) +  # Dodge the boxes for each type per PPI_cat
  #facet_wrap(experiment ~ PPI_cat, ncol = 3) +  # Facet by PPI_cat, each PPI_cat gets a row
  scale_fill_manual(values = PPI_color_platte) +
  #scale_fill_manual(values = c("random" = "grey", "btwn.cent" = "white")) +
  geom_signif(
    comparisons = list(
      c("CTS", "HiGCTS"), c("HiGCTS", "HiG"),c("CTS", "HiG")
    ),
    map_signif_level = TRUE,
    step_increase = 0.1,  # Adjusts spacing between the lines
    #aes(group = type),
    test = "wilcox.test",  # Perform a t-test to calculate significance
	p.adjust.method = "holm"  # Adjust p-values using the Holm method, default is "bonferroni"
  ) + 
  theme_minimal() +
  labs(
    title = "Comparison of Robustness Measures by cent.btw",
    x = "PPI Category", 
    y = "Component Percentage Remaining"
  ) +
  theme(legend.position = "top")
print(g2)

pdf(file='box_wilcox-test_attack_GSE175634.pdf', width=3.5)
print(g)
print(g2)
dev.off()

# confirm taht the figure showing adjusted p-values 
wilcox.test(targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat=='CTS')],
			targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat=='HiGCTS')])
# W = 10518, p-value = 0.06667
wilcox.test(targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat=='CTS')],
			targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat=='HiG')])
# W = 380137, p-value < 2.2e-16
wilcox.test(targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat=='HiGCTS')],
			targeted_vertex_data$comp.pct[which(targeted_vertex_data$PPI_cat=='HiG')])
# 122960, p-value = 9.951e-08
 
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
fold_change
#    PPI_cat fold_change_edge fold_change_vertex
#   <fct>              <dbl>              <dbl>
# 1 CTS                1.06                1.54
# 2 HiGCTS             0.991               1.50
# 3 HiG                0.964               1.04


## for each signature, comapre the comp.pct between targeted attach to its random removal to access the significance of 'hub'
#######################################
	Hub_effect = array(dim=length(graph_list))
	names(Hub_effect) = names(graph_list)
	for(j in names(graph_list)){
		tmp = subset(robustness.dt, signature==j)
		dim(tmp)  # 230  9
		table(tmp$type)
		#    Random edge removal   Targeted edge attack  Random vertex removal Targeted vertex attack 
        #                   37                     37                     78 				78
		x = subset(tmp, type=='Random edge removal')$comp.pct
		y = subset(tmp, type=='Targeted edge attack')$comp.pct
		Hub_effect[j] = wilcox.test(x,y)$p.value
	}
	df=data.frame(edge_p = Hub_effect,
		edge_p.adj = p.adjust(Hub_effect, method='bonferroni'))
	
	Hub_effect = array(dim=length(graph_list))
	names(Hub_effect) = names(graph_list)
	for(j in names(graph_list)){
		tmp = subset(robustness.dt, signature==j)
		x = subset(tmp, type=='Random vertex removal')$comp.pct
		y = subset(tmp, type=='Targeted vertex attack')$comp.pct
		Hub_effect[j] = wilcox.test(x,y)$p.value
	}
	df$node_p = Hub_effect 
	df$node_p.adj = p.adjust(df$node_p, method='bonferroni')
 
  df[which(df$edge_p.adj<0.05),] 
#                   edge_p    edge_p.adj       node_p  node_p.adj
# HiG_0       3.233735e-11  6.790843e-10 0.4535332327 1.000000000
# HiG_1       5.609095e-99  1.177910e-97 0.4873374322 1.000000000
# HiG_2       1.775622e-14  3.728806e-13 0.7325582510 1.000000000
# HiG_CP      2.577084e-05  5.411876e-04 0.1173190552 1.000000000
# HiG_4       2.119556e-03  4.451068e-02 0.3270656926 1.000000000
# HiG_5       1.002736e-03  2.105745e-02 0.2013649833 1.000000000
# HiG_6       6.632838e-04  1.392896e-02 0.2476923942 1.000000000
# HiG_7       4.199996e-04  8.819991e-03 0.4446911224 1.000000000
# HiG_muscle  4.576384e-08  9.610405e-07 0.6311089349 1.000000000
# HiG_9      5.494467e-130 1.153838e-128 0.3688965383 1.000000000
# CTS_CP.1    2.889371e-05  6.067680e-04 0.0001465217 0.003076956

df[which(df$node_p.adj<0.05),] 
               # edge_p  edge_p.adj       node_p   node_p.adj
# CTS_CP   0.0328342394 0.689519027 3.307307e-08 6.945346e-07
# CTS_CP.1 0.0002442544 0.005129342 1.552792e-04 3.260863e-03


  df$clust = lapply(rownames(df), function(x) unlist(strsplit(x, split='_'))[2]) %>% unlist
	df$PPI_cat = lapply(rownames(df), function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist  
	df$PPI_cat = factor(df$PPI_cat, levels = c("HiG" ,   "HiGCTS", "CTS"))
    df$clust = factor(df$clust, levels = c("0" , "1",  "2"  ,  "CP" ,  "CP.1" ,   "4" ,  "5"  ,
					"6"  ,  "7" ,  "muscle" ,  "9"  , "10"   ,  "endoderm","12"        ))
   
 # Create a new column to differentiate edge and node
tmp = df[,c(1,2,5,6)] %>% 
		mutate(type = 'edge')
tmp2 = df[,c(3,4,5,6)] %>% 
		mutate(type = 'vertex')
colnames(tmp)[1:2] = colnames(tmp2)[1:2] = c('p','p.adj')		
tmp = rbind(tmp, tmp2)
tmp$PPI_cat = factor(tmp$PPI_cat, levels = c("HiG" ,   "HiGCTS", "CTS"))

# Plot the results with subpanels
g = ggplot(tmp, aes(x = as.factor(clust), y = -log10(p.adj), color = PPI_cat)) + 
  geom_point(aes(shape = PPI_cat), size = 3) +  # Points for each PPI_cat
  geom_line(aes(group = PPI_cat, linetype = PPI_cat)) +  # Lines for each PPI_cat
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Horizontal line at adj.p = 0.05
  labs(
    x = "Cluster", 
    y = "-log10(Adjusted p-values)", 
    color = "State", 
    shape = "State"
  ) +  
  scale_color_manual(values = PPI_color_platte) +  
  scale_shape_manual(values = c("HiG" = 16, "CTS" = 17, "HiGCTS" = 18)) +  
  facet_wrap(~ type, scales = "free", nrow = 2) +  # Split into two rows for edge and node
  theme_minimal() +
  theme(legend.position = "top") +
  ggtitle("Adjusted p-values of Wilcox test comparing targeted attack vs random removal") 
print(g)
## discuss the irrelevent of the p-values and node strength levels
  tmp$count = c(sapply(graph_list, ecount), 
				sapply(graph_list, vcount))
  g2 = ggplot(tmp, aes(x = log10(count), y = -log10(p.adj), color = PPI_cat)) + 
		  geom_point(aes(shape = PPI_cat), size = 3) +  # Points for each PPI_cat
		  #geom_line(aes(group = PPI_cat, linetype = PPI_cat)) +  # Lines for each PPI_cat
		  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Horizontal line at adj.p = 0.05
		  labs(
			x = "log10(PPI size)", 
			y = "-log10(Adjusted p-values)", 
			color = "State", 
			shape = "State"
		  ) +  
		  scale_color_manual(values = PPI_color_platte) +  
		  scale_shape_manual(values = c("HiG" = 16, "CTS" = 17, "HiGCTS" = 18)) +  
		  facet_wrap( ~ type , scales = "free", ncol = 2) +  # Split into two rows for edge and node
		  theme_minimal() +
		  theme(legend.position = "top") +
		  ggtitle("Adjusted p-values of Wilcox test comparing targeted attack vs random removal") 
print(g2)
    pdf(file='line.adjp_wilcox_attack_GSE175634.pdf')
	print(g)
	print(g2)
	dev.off()

	
    #############################################################
    # 4.2) estimate the significance of robustness using AUC (????can't repeat mopuse observation!!)
	# the observed lin always sit left to the simulated distributions,
	# suggesting gene in the observed PPI is mroe important than randomly selected genes, 
	# becasue targeted attack at them significantly reduced the component size of the network
	#
	# for the 
    #############################################################
    library(MLmetrics)
    library(sm)
	IDs_of_CTS = c("muscle" ,  "endoderm" ,"CP"   ,    "CP.1"  )
	names(graph_list)
 # [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"           "HiG_5"          
 # [7] "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"   
# [13] "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"     
# [19] "CTS_endoderm"    "CTS_CP"          "CTS_CP.1" 

     observed_auc_list = list()
	 for(j in names(graph_list)){
		observed_auc_list[[j]] = Area_Under_Curve(subset(robustness.dt, signature==j & type=='Targeted vertex attack')$removed.pct, 
                     subset(robustness.dt, signature==j & type=='Targeted vertex attack')$comp.pct ) 
    } 
	observed_auc_list %>% unlist
         # HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4           HiG_5 
      # 0.4780639       0.4887626       0.4884125       0.4093552       0.4501320       0.4548244 
          # HiG_6           HiG_7      HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
      # 0.4317425       0.4583025       0.4729244       0.4877886       0.4473469       0.3010204 
         # HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP     HiGCTS_CP.1      CTS_muscle 
      # 0.3459596       0.4111111       0.3809524       0.2318008       0.2314815       0.3209332 
   # CTS_endoderm          CTS_CP        CTS_CP.1 
      # 0.4159483       0.1801304       0.2327409 
	
	df_AUC = data.frame(auc=observed_auc_list %>% unlist, 
				signature = names(observed_auc_list), 
				PPI_cat = lapply(names(observed_auc_list), function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist)
	df_AUC$PPI_cat  = factor(df_AUC$PPI_cat , levels = c('CTS',"HiGCTS" ,"HiG"))
	
g3 = ggplot(df_AUC, aes(x = PPI_cat, y = auc, fill = PPI_cat)) + 
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.75)) +  # Dodge the boxes for each type per PPI_cat
  #facet_wrap(experiment ~ PPI_cat, ncol = 3) +  # Facet by PPI_cat, each PPI_cat gets a row
  scale_fill_manual(values = PPI_color_platte) +
  #scale_fill_manual(values = c("random" = "grey", "btwn.cent" = "orange")) +
  geom_signif(
    comparisons = list(
      c("CTS", "HiGCTS"), c("HiGCTS", "HiG"),c("CTS", "HiG")
    ),
    map_signif_level = TRUE,
    step_increase = 0.1,  # Adjusts spacing between the lines
    #aes(group = type),
    test = "wilcox.test" , # Perform a t-test to calculate significance
	p.adjust.method = "holm"   # default is 
  ) + 
  theme_minimal() +
  labs(
    title = "Comparison of Robustness Measures by cent.btw",
    x = "PPI Category", 
    y = "AUC of targeted vertext attack"
  ) +
  theme(legend.position = "top")
print(g3)
dev.copy2pdf(file='box_wilcox-test_attack_AUC_GSE175634.pdf', width=3.5)
	
	
	# ks.test(subset(robustness.dt, signature=='cardiac.a')$comp.pc, subset(robustness.dt, signature=='cardiac.c')$comp.pc)
    # p-value < 2.2e-16
    
	vn = sapply(graph_list, vcount)
	pn = array(dim=length(graph_list))
	names(pn) = names(graph_list)
	for(j in names(vn)){
		pn[j] = igraph::strength(graph_list[[j]], weights = E(graph_list[[j]])$weight)  %>% mean/vn[j]
    }
	vn
     # HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4           HiG_5 
            # 325            1155             242             132             104             339 
          # HiG_6           HiG_7      HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
            # 106             104             118            1467              70              42 
         # HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP     HiGCTS_CP.1      CTS_muscle 
             # 22              20               7              29              27              61 
   # CTS_endoderm          CTS_CP        CTS_CP.1 
             # 64              72              77 



	sum(vn[1:13])  #[1] 4226  #!!!
	sum(vn[14:17])  #[1] 83
	sum(vn[18:21])  # [1] 274
	
	#### significance evidence 1, comp.pct (the ratio of maximal component size after each random removal to the observed graph's maximal component size)
    ### was calculated after	random vertex_attaching 
	## then calculate the AUC of CTS_muscle at muscle cluster vs HiG at nearby stady clusters (C1 & C5) 
    ## this simulation runs longer than expected, therefore go to midway3 !!! 
	# N = 1000  for vertex removal N=100 for edge removal 
	
     # refer to 11.2_CTS_cardiac_network_robustness_Miudway3.simulation.R
  ## run on Midway3 following is teh log. DO NOT REPEAT !!!)
# scp -p -r  xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weight/AUC*.rds  F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weight_shrink/.

attac_V_random = readRDS(paste0('AUC_compt.pct_attac_V_random_1000runs_',s,'weighted.RDS'))
(N = nrow(attac_V_random))  # 1000
		    # plot density
    # https://www.statmethods.net/graphs/density.html

# Open PDF device for all plots
# # Create folder if it doesn't exist
# if (!dir.exists("AUC_vertex")) {
  # dir.create("AUC_vertex")
# }
pdf("Vetex_All_AUC_density_plots.pdf", width=20, height=7)
par(mfrow=c(2,7)) 	
	
for(j in c('muscle','endoderm', 'CP', 'CP.1')){	
    if(grepl('.', j, fixed=TRUE)) j0= unlist(strsplit(j, '.', fixed=TRUE))[1]  else j0=j# there is only one 'CP' cluster to calculate HiG
	
	lf = factor(c(rep('CTS', N), rep('HiGCTS', N), rep('HiG', N)),
            levels = c('CTS', 'HiGCTS', 'HiG'))  # Specify order explicitly  # ** corrected
    ## set xlim to plot all three line togewther
	x_lines = c(attac_V_random[, paste0('CTS_',j)], attac_V_random[, paste0('HiGCTS_',j)], attac_V_random[,paste0('HiG_', j0)])
	tmp = sm.density.compare(x_lines, lf, xlab='auc', 
			col= PPI_color_platte,  # ** corrected 
			xlim = range(c(x_lines, observed_auc_list[[paste0('CTS_',j)]],observed_auc_list[[paste0('HiGCTS_',j)]], observed_auc_list[[paste0('HiG_',j0)]])))
	abline(v=observed_auc_list[[paste0('CTS_',j)]], col=PPI_color_platte['CTS'], lty=2)
    abline(v=observed_auc_list[[paste0('HiGCTS_',j)]], col=PPI_color_platte['HiGCTS'], lty=2)
	abline(v=observed_auc_list[[paste0('HiG_',j0)]], col=PPI_color_platte['HiG'], lty=2)
	ks.test(attac_V_random[,paste0('CTS_',j)], attac_V_random[,paste0('HiG_',j0)])
	# D = 0.991, p-value < 2.2e-16
    ks.test(attac_V_random[,paste0('CTS_',j)], attac_V_random[,paste0('HiG_',j0)])
	# D = 0.998, p-value < 2.2e-16
    
    x = length(which(attac_V_random[,paste0('CTS_',j)]>observed_auc_list[[paste0('CTS_',j)]]))
	y = length(which(attac_V_random[,paste0('CTS_',j)]<observed_auc_list[[paste0('CTS_',j)]]))
	mtext(paste('CTS',j,'p= ', x/N,";", y/N),side=3, line=0,
          col=PPI_color_platte['CTS'])
	x = length(which(attac_V_random[,paste0('HiGCTS_',j)]>observed_auc_list[[paste0('HiGCTS_',j)]]))
	y = length(which(attac_V_random[,paste0('HiGCTS_',j)]<observed_auc_list[[paste0('HiGCTS_',j)]]))
    mtext(paste('HiGCTS',j,'p=', x/N, ";", y/N),side=3, line=1,
          col=PPI_color_platte['HiGCTS'])
	x = length(which(attac_V_random[,paste0('HiG_',j0)]>observed_auc_list[[paste0('HiG_',j0)]]))
	y = length(which(attac_V_random[,paste0('HiG_',j0)]<observed_auc_list[[paste0('HiG_',j0)]]))
	mtext(paste('HiG',j,'p=', x/N,";", y/N),side=3, line=2,
          col=PPI_color_platte['HiG'])
    #dev.copy2pdf(file=paste0('AUC_vertex/AUC_density_',j,'_random.vertex.removal_btwn.cent.pdf') )  
}    

clusterIDs = lapply(names(graph_list), function(x) unlist(strsplit(x, split='_'))[2]) %>% unlist
for(j in setdiff(clusterIDs,c('muscle','endoderm', 'CP', 'CP.1'))){	
	tmp= density( attac_V_random[,paste0('HiG_', j)])
	v = observed_auc_list[[paste0('HiG_',j)]]
    plot(tmp, col='grey',  xlim = range(c(v, attac_V_random[,paste0('HiG_', j)])),
		  xlab = '% maximal component size after removal to the observed maximal component size', 
		  ylab = 'Density', main = j)
	abline(v=v, col=PPI_color_platte['HiG'], lty=2)
    
 	x = length(which(attac_V_random[,paste0('HiG_',j)]>observed_auc_list[[paste0('HiG_',j)]]))
	y = length(which(attac_V_random[,paste0('HiG_',j)]<observed_auc_list[[paste0('HiG_',j)]]))
	mtext(paste('HiG',j,'p=', x/N,";", y/N),side=3, line=2,
          col=PPI_color_platte['HiG'])
    #dev.copy2pdf(file=paste0('AUC_vertex/AUC_density_',j,'_random.vertex.removal_btwn.cent.pdf') )  
}
dev.off()	
	
	
   #############################################################
    # 4.2) estimate the significance of robustness   using AUC , edge removal  
    #############################################################
   
     observed_auc_list = list()
	 for(j in names(graph_list)){
		observed_auc_list[[j]] = Area_Under_Curve(subset(robustness.dt, signature==j & type=='Targeted edge attack')$removed.pct, 
                     subset(robustness.dt, signature==j & type=='Targeted edge attack')$comp.pct ) 
    } 
	observed_auc_list %>% unlist
  #         HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4           HiG_5 
  #     0.6858468       0.7203489       0.6887157       0.5597807       0.6184507       0.6724290 
  #         HiG_6           HiG_7      HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
  #     0.6910997       0.6878484       0.7011913       0.7171298       0.6163821       0.7020263 
  #        HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP     HiGCTS_CP.1      CTS_muscle 
  #     0.6738095       0.7132353       0.6818182       0.5818713       0.5829268       0.6145935 
  #  CTS_endoderm          CTS_CP        CTS_CP.1 
  #     0.6812997       0.5692761       0.5367428 
	
	# ks.test(subset(robustness.dt, signature=='cardiac.a')$comp.pc, subset(robustness.dt, signature=='cardiac.c')$comp.pc)
    # p-value < 2.2e-16
    
	en = mn = array(dim=length(graph_list))
	names(en) = names(pn) = names(graph_list)
	for(j in names(vn)){
		# Edge strength is simply the sum of the strengths of the two nodes (vertices) connected by each edge
		deg <- igraph::strength(graph_list[[j]], weights = E(graph_list[[j]])$weight)
		en[j] =  apply(as_edgelist(graph_list[[j]]), 1, function(x) deg[x[1]] + deg[x[2]])
		mn[j] = (en[j] %>% mean) 
    }
	mn # averaged edge strength
    #                                                       HiG_0           HiG_1           HiG_2 
    #          NA              NA              NA      2.89182135      1.38987941      3.04511638 
    #      HiG_CP           HiG_4           HiG_5           HiG_6           HiG_7      HiG_muscle 
    #  0.50428237      7.26557084      4.89424044      2.41631360      0.96833129      1.59429487 
    #       HiG_9          HiG_10    HiG_endoderm          HiG_12   HiGCTS_muscle HiGCTS_endoderm 
    #  0.54808290     11.22242436      2.15378565      0.04041787      1.82650318      3.00357796 
    #   HiGCTS_CP     HiGCTS_CP.1      CTS_muscle    CTS_endoderm          CTS_CP        CTS_CP.1 
    #  0.41419400      0.33387924      0.14418649      0.25173236      0.61078904      0.34583759 

attac_E_random = readRDS(paste0('AUC_compt.pct_attac_E_random_100runs_',s,'weighted.RDS'))
(N = nrow(attac_E_random))  # 100
		
pdf("edge_All_AUC_density_plots.pdf", width=20, height=7)
par(mfrow=c(2,7)) 	
    # plot density
    # https://www.statmethods.net/graphs/density.html
for(j in c('muscle','endoderm', 'CP', 'CP.1')){	
    if(grepl('.', j, fixed=TRUE)) j0= unlist(strsplit(j, '.', fixed=TRUE))[1]  else j0=j# there is only one 'CP' cluster to calculate HiG
	if(!all(is.na(attac_E_random[, paste0('HiGCTS_',j)]))) {
		lf = factor(c(rep('CTS', N), rep('HiGCTS', N), rep('HiG', N)),
            levels = c('CTS', 'HiGCTS', 'HiG'))  # Specify order explicitly  # ** corrected
		## set xlim to plot all three line togewther
		x_lines = c(attac_E_random[, paste0('CTS_',j)], attac_E_random[, paste0('HiGCTS_',j)], attac_E_random[,paste0('HiG_', j0)])
		tmp = sm.density.compare(x_lines, lf, xlab='auc',
			col= PPI_color_platte,  # ** corrected 
			xlim = range(c(x_lines, observed_auc_list[[paste0('CTS_',j)]],observed_auc_list[[paste0('HiGCTS_',j)]], observed_auc_list[[paste0('HiG_',j0)]])))
		} else {
		lf = factor(c(rep('CTS', N), rep('HiG', N)),
            levels = c('CTS', 'HiG'))  # Specify order explicitly  # ** corrected
				## set xlim to plot all three line togewther
		x_lines = c(attac_E_random[, paste0('CTS_',j)], attac_E_random[,paste0('HiG_', j0)])
		tmp = sm.density.compare(x_lines, lf, xlab='auc', 	
			col= PPI_color_platte,  # ** corrected 		
			xlim = range(c(x_lines, observed_auc_list[[paste0('CTS_',j)]],observed_auc_list[[paste0('HiGCTS_',j)]], observed_auc_list[[paste0('HiG_',j0)]])))
		}
	abline(v=observed_auc_list[[paste0('CTS_',j)]], col=PPI_color_platte['CTS'], lty=2)
    abline(v=observed_auc_list[[paste0('HiGCTS_',j)]], col=PPI_color_platte['HiGCTS'], lty=2)
	abline(v=observed_auc_list[[paste0('HiG_',j0)]], col=PPI_color_platte['HiG'], lty=2)
	ks.test(attac_E_random[,paste0('CTS_',j)], attac_E_random[,paste0('HiG_',j0)])
	# D = 0.991, p-value < 2.2e-16
    ks.test(attac_E_random[,paste0('CTS_',j)], attac_E_random[,paste0('HiG_',j0)])
	# D = 0.998, p-value < 2.2e-16
    
    x = length(which(attac_E_random[,paste0('CTS_',j)]>observed_auc_list[[paste0('CTS_',j)]]))
	y = length(which(attac_E_random[,paste0('CTS_',j)]<observed_auc_list[[paste0('CTS_',j)]]))
	mtext(paste('CTS',j,'p= ', x/N,";", y/N),side=3, line=0,
          col=PPI_color_platte['CTS'])
	x = length(which(attac_E_random[,paste0('HiGCTS_',j)]>observed_auc_list[[paste0('HiGCTS_',j)]]))
	y = length(which(attac_E_random[,paste0('HiGCTS_',j)]<observed_auc_list[[paste0('HiGCTS_',j)]]))
    mtext(paste('HiGCTS',j,'p=', x/N, ";", y/N),side=3, line=1,
          col=PPI_color_platte['HiGCTS'])
	x = length(which(attac_E_random[,paste0('HiG_',j0)]>observed_auc_list[[paste0('HiG_',j0)]]))
	y = length(which(attac_E_random[,paste0('HiG_',j0)]<observed_auc_list[[paste0('HiG_',j0)]]))
	mtext(paste('HiG',j,'p=', x/N,";", y/N),side=3, line=2,
          col=PPI_color_platte['HiG'])
    #dev.copy2pdf(file=paste0('AUC_edge/AUC_density_',j,'_random.edge.removal_btwn.cent.pdf') )  
}    

clusterIDs = lapply(names(graph_list), function(x) unlist(strsplit(x, split='_'))[2]) %>% unlist
for(j in setdiff(clusterIDs,c('muscle','endoderm', 'CP', 'CP.1'))){	
	tmp= density( attac_E_random[,paste0('HiG_', j)] )
	v = observed_auc_list[[paste0('HiG_',j)]]
    plot(tmp, col='grey',  xlim = range(c(v, attac_E_random[,paste0('HiG_', j)])),
		  xlab = '% maximal component size after removal to the observed maximal component size', 
		  ylab = 'Density', main = 'Density plots for different muscle types')
	abline(v=v,  lty=2)
    
 	x = length(which(attac_E_random[,paste0('HiG_',j)]>observed_auc_list[[paste0('HiG_',j)]]))
	y = length(which(attac_E_random[,paste0('HiG_',j)]<observed_auc_list[[paste0('HiG_',j)]]))
	mtext(paste('HiG',j,'p=', x/N,";", y/N),side=3, line=2,
          col=PPI_color_platte['HiG'])
   # dev.copy2pdf(file=paste0('AUC_edge/AUC_density_',j,'_random.edge.removal_btwn.cent.pdf') )  
}
	
dev.off()		
	
	