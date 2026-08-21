# this version tests v9 code 
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
library(ggpubr)  # ggarrange()
library("gridExtra")
library(ggrepel)
library(ggpubr) # resui#E7298A by stat_compare_means()
library(igraph)
library(rstatix)
 
library(brainGraph)

code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_STRING", "code_core")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
PPI_size_platte = c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)
setwd(results_dir)


 # refer to 11.2.0_weighted_graph_attack_robustness.R
 #for(s in c("combined","ratio", "zscore", "diff")){
s = "combined"
file = paste0('GSE175634_STRING_graph_perState_simplified_',s,'weighted.rds')
graph_list <- readRDS( file)  
	
	#graph_list <- lapply(graph_list, simplify) #!!!!!!!!!!!!!!!!!!!
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
sapply(graph_list, vcount)
          # HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4           HiG_5 
            # 325            1155             242             132             104             339 
          # HiG_6           HiG_7      HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
            # 106             104             118            1467              70              42 
         # HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP     HiGCTS_CP.1      CTS_muscle 
             # 22              20               7              29              27              61 
   # CTS_endoderm          CTS_CP        CTS_CP.1 
             # 64              72              77 
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


 
    
    #######################################################################  
    ## fr. https://www.nature.com/articles/s41598-021-03625-w#Abs1 
    # 'Knn, page rank, and strength are the most relevant GRN features'
    ## First, evaluate the pageRanks,
    ## PageRank’s main difference from EigenCentrality is that it accounts for link direction.Like EigenCentrality,
    ## PageRank can help uncover influential or important nodes whose reach extends beyond just their direct connections. 
    ## It’s especially useful in scenarios where link direction is important:
    # 
    # https://igraph.org/r/doc/page_rank.html
	CTS.ID = c('muscle', 'essential','CP','CP.1')
	
    # tmp <- page_rank(graph_list[[1]], weights = E(graph_list[[1]])$weight)$vector 
    page = lapply(graph_list, function(x) page_rank(x, directed = FALSE, weights = E(x)$weight)$vector )   

    df = lapply(page, function(x) data.frame(PageRank=x, gene=names(x)) %>% arrange(desc(PageRank)))  %>%
      rbindlist(., idcol=names(.))   
    colnames(df)[1] ='signature'
    df$PPI_cat = lapply(df$signature, function(x) unlist(strsplit(x , '_'))[1]) %>% unlist %>%
			factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 
	dim(df)[1] # 4583    

    ic = lapply(graph_list, function(x) eigen_centrality(x, weights = E(x)$weight)$vector)
	IC = lapply(ic, function(x) data.frame(EigenCentrality=x, gene=names(x)) %>% arrange(desc(EigenCentrality)))  %>%
      rbindlist(., idcol=names(.))  
    colnames(IC)[1] ='signature'   
	dim(IC)  # [1] 4583    3
	df = merge(df, IC, by=c("signature","gene")) 
	dim(df)  # [1] 4583    5  # was nrow increased !!! now, equal !!


    ## A) estimate the random PageRank by rewiring the edges while keeping the pro; this loop takes a while, Do NOT repeat !!
    vn = sapply(graph_list, vcount)   #lengths(graph_list)
	# Step 1: Calculate all pc values first
	all_pc <- numeric(length(graph_list))
	names(all_pc) <- names(graph_list)
	for(i in names(graph_list)){
	  g <- graph_list[[i]]
	  all_pc[i] <- mean(igraph::strength(g, weights = E(g)$weight)) / vn[i]
	}
	
	if(max(all_pc)>1) {  # ** v5 new 
		# Step 2: Normalize to [0.01, 0.99] to preserve variability
		all_pc <- 0.01 + 0.98 * (all_pc - min(all_pc)) / (max(all_pc) - min(all_pc))
	}
	
    N = 1000
    set.seed(1234)
    pr_P = list() 
    for(i in names(graph_list)){
	  g = graph_list[[i]]
      pc = all_pc[i]
      g.random = list()
      for(j in 1:N){
        g.random[[j]] = rewire(graph_list[[i]], each_edge(prob = pc)) ## rewiring the edges while keeping the pro
		# cat(range(E(g.random[[1]])$weight)) # [1] 0.404 1.866 
      }
	  
      pr_random = lapply(g.random, function(x) page_rank(x, weights = E(x)$weight)$vector )   
   
      tmp= lapply(pr_random, function(x) data.frame(PageRank=x, gene=names(x)) %>% arrange(desc(PageRank)))  %>%
          rbindlist(., idcol=names(.)) 
      head(tmp)
      #     PageRank   gene
      # 1: 0.05167547   Rnd3
      # 2: 0.05151246 Twist1
      # 3: 0.04754937  Mef2c
      for(j in V(graph_list[[i]])$name){
        pr_P[[i]][j] <- length(which(subset(tmp, gene==j)$PageRank >= page[[i]][j]))/N
      }
    }	
    saveRDS(pr_P, file='GSE175634_PageRank_Pvalue_by_rewiring.rds')
    
   
	pr_P = readRDS(file='GSE175634_PageRank_Pvalue_by_rewiring.rds')
    tmp = lapply(pr_P, function(x) data.frame(p.PageRank=x, gene=names(x))) %>%
      rbindlist(., idcol=names(.))   
    colnames(tmp)[1] ='signature'  
	
    df = merge(df, tmp, by=c("signature","gene")) %>%
	  group_by(signature) %>%
	  mutate(rank_by_p.PR = rank(p.PageRank)) %>%
	  mutate(rank_by_PR = rank(-PageRank)) %>%
	  ungroup()
    
	head(df)
 #  A tibble: 6 × 8
 #   signature gene     PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR
#   <chr>     <chr>       <dbl> <fct>             <dbl>      <dbl>        <dbl>      <dbl>
# 1 CTS_CP    AFAP1L2   0.00277 CTS           0              0.855         48           62
# 2 CTS_CP    ANKRD26   0.0116  CTS           0.0000657      0.759          8           34
# 3 CTS_CP    ARHGAP28  0.00277 CTS           0              0.858         56.5         62
# 4 CTS_CP    ARL4C     0.0185  CTS           0              0.839         38           22
# 5 CTS_CP    ARRDC3    0.00277 CTS           0              0.861         63.5         62
# 6 CTS_CP    ASH2L     0.0212  CTS           0.0676         0.822         32.5         17
dim(df)  #[1] 4583     8

subset(df, gene=='ISL1')
#   signature   gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR
#   <chr>       <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>      <dbl>
# 1 CTS_CP.1    ISL1   0.0366  CTS            0.00729       0.872           67          4
# 2 HiGCTS_CP.1 ISL1   0.0794  HiGCTS         0.733         0.883           16          3
# 3 HiG_1       ISL1   0.00181 HiG            0.0539        0.782         1023        113
# 4 HiG_5       ISL1   0.00524 HiG            0.0843        0.281          180         51
# 5 HiG_9       ISL1   0.00197 HiG            0.116         0.767         1375         67
# 6 HiG_CP      ISL1   0.0136  HiG            0.000546      0.003           13         19

    
    # number of significantly high pagerank per PPI_cats, too much control !
    n.pr.high <- lapply(names(graph_list), function(x) nrow(subset(df, signature==x & p.PageRank<0.05))) %>% unlist
    names(n.pr.high) = names(graph_list)
    n.pr.high 
    #       HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4 
    #           1               0               2              44              25 
    #       HiG_5           HiG_6           HiG_7      HiG_muscle           HiG_9 
    #           0               0               0               0               0 
    #      HiG_10    HiG_endoderm          HiG_12   HiGCTS_muscle HiGCTS_endoderm 
    #          15               0               0               0               0 
    #   HiGCTS_CP     HiGCTS_CP.1      CTS_muscle    CTS_endoderm          CTS_CP 
    #           0               0               0               0               0 
    #    CTS_CP.1 
    #           0 



	# write.table(df, file='df_PAGERANK.tsv',sep='\t', row.names=F)  #!!!!!!!!



	## ANND (Average Nearest Neighbor strength) measures the strength-strength dependence adjacent to a vertex
    ##
	## here, we evaluate the ann, which only works with simple graphs,
    # is often used to characterize dependencies between strengths of a node and its neighbors in a network.
    # a non-simple graph is to have multiple edges connecting two nodes or for there to be a self-edge.
    # igraph::knn():
    # res$knn:
	# A numeric vector giving the average nearest neighbor strength for all vertices in the graph.
    # res$knnk :	
	# Calculate the ANND (average nearest neighbor strength) of the given vertices and the same quantity in the function of vertex strength
    # A numeric vector, its length is the maximum (total) vertex strength in the graph. 
    # The first element is the average nearest neighbor strength of vertices with strength one, etc.
    # for zero strength vertices the answer in ‘knn’ is NaN
    annd_observed = list()
    for(i in names(graph_list)){
         G = graph_list[[i]]
         # remove unconnected nodes
         V_Isolated = which(degree(G)==0)
         G = delete_vertices(G, V_Isolated)  #!!!!!!!!!  
         annd_observed[[i]] = knn(G, weights = E(G)$weight)$knn  
    }
    #annd_observed = lapply(graph_list, function(x) knn(x, weights = E(x)$weight)$knn )   # ** update
    any(is.na( annd_observed[['CTS_CP.1']])) #FALSE
    rm(G)
	
    ## A) estimate the random annd by rewiring the edges while keeping the pro; this loop takes a while, Do NOT repeat !!
    annd_P= list()
    N = 1000
    set.seed(1234)
    pr_P = list() 
    for(i in names(graph_list)){
	  g = graph_list[[i]]
      pc = all_pc[i]
      g.random = list()
      for(j in 1:N){
        g.random[[j]] = rewire(graph_list[[i]], each_edge(prob = pc)) ## rewiring the edges while keeping the pro
		# cat(range(E(g.random[[1]])$weight)) # [1] 0.404 1.866 
      }
	  
	  annd_random = lapply(g.random, function(x) knn(x, weights = E(x)$weight )$knn )   # ** update
      
      tmp= lapply(annd_random, function(x) data.frame(annd=x, gene=names(x)) %>% arrange(desc(annd)))  %>%
        rbindlist(., idcol=names(.)) 
      head(tmp, 3)
       # annd    gene
      # <num>  <char>
# 1: 19.00000     WLS
# 2: 15.18162   CHIC2
# 3: 14.00000   PLCE1
    
      for(j in V(graph_list[[i]])$name){
        annd_P[[i]][j] <- length(which(subset(tmp, gene==j)$knn >= annd_observed[[i]][j]))/N
      }
    }
    saveRDS(annd_P, file='GSE175634_annd_Pvalue_by_rewiring.rds')
    
	annd_P = readRDS( file='GSE175634_annd_Pvalue_by_rewiring.rds')
      ## reorder df$df$signature to be #E7298A #7570B3 and #E6AB02  
	unique(df$PPI_cat)  # CTS   HiGCTS   HiG  

	df = rbind(subset(df, PPI_cat=='CTS'),
						subset(df, PPI_cat=='HiGCTS'),
						subset(df, PPI_cat=='HiG')
						)
	unique(df$PPI_cat)  # CTS  HiGCTS   HiG   
	
	tmp = lapply(annd_observed, function(x) data.frame(annd=x, gene=names(x)) %>% arrange(desc(annd)))  %>%
      rbindlist(., idcol=names(.))  
    colnames(tmp)[1] ='signature'
    df = merge(df, tmp, by=c("signature","gene"))
     dim(df)  #[1] 4495    9

    annd_P[['CTS_CP.1']]
    tmp = lapply(annd_P, function(x) data.frame(p.annd=x, gene=names(x))) %>%
      rbindlist(., idcol=names(.))  
    colnames(tmp)[1] ='signature'
    dim(tmp)  #[1] 4583    3
	
	df = merge(df, tmp, by=c("signature","gene"))
    df[which(is.na(df$knn)),'p.annd'] = NA ## due to nrow(df) 4878 > nrow(tmp) 
    dim(df)  #[1] 4495   10

    ## merge back the normalized strength of vertex
	# normalized_strength <- strength(g) / (vcount(g) - 1)
    # refer to 11.1_CTS_cardiac_network_strengthDistribution.R 
	V_strength = lapply(graph_list, function(g) {
		# Calculate strength and sort
		strength <- strength(g, weights = E(g)$weight) 
		strength_sorted <- sort(strength, decreasing = TRUE)	
		# Create data frame
		data.frame(
			strength = strength_sorted,
			gene = names(strength_sorted), 
			id = 1:length(strength_sorted)
		)
	}) %>% 
	rbindlist(., idcol = "signature") %>%
	dplyr::rename('rank_by_strength' = 'id')
	
	V_strength_norm <- lapply(graph_list, function(g) {
		# Calculate normalized strength and sort
		norm_strength <- strength(g, weights = E(g)$weight) / (vcount(g) - 1)
		norm_strength_sorted <- sort(norm_strength, decreasing = TRUE)		
		# Create data frame
		data.frame(
			normalized.strength = norm_strength_sorted,
			gene = names(norm_strength_sorted), 
			id = 1:length(norm_strength_sorted)
		)
	}) %>% 
	rbindlist(., idcol = "signature") %>%
	dplyr::rename('rank_by_normalized.strength' = 'id')
	
	V_strength = merge(V_strength,V_strength_norm, by=c("signature","gene"))
	dim(V_strength ) # 4583    6

	## add the V_strength & V_strength_norm infor
	df <- merge(df, V_strength, by=c("signature","gene"))
	dim(df) # 4495   14
	head(df,3)
#   signature    gene   PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR
# 1    CTS_CP ANKRD26 0.01159497     CTS    6.574061e-05      0.759          8.0
# 2    CTS_CP   ARL4C 0.01846722     CTS    0.000000e+00      0.839         38.0
# 3    CTS_CP   ASH2L 0.02116792     CTS    6.761320e-02      0.822         32.5
#   rank_by_PR     annd p.annd   strength rank_by_strength normalized.strength
# 1         34 4.000000      0 0.01632206               39        0.0002298882
# 2         22 1.000000      0 0.01012606               44        0.0001426206
# 3         17 8.864765      0 0.18858257               17        0.0026560926
#   rank_by_normalized.strength
# 1                          39
# 2                          44
# 3                          17

	
 #   table(df$strength>=df$annd)   # all FALSE
 	 
  ## # Add rank_by_ANND column and rerank strength, normalized strength by considering ties !!!
	df <- df %>%
	  group_by(signature) %>%  # Group by 'signature'
	  mutate(rank_by_strength = rank(-strength, na.last = "keep")) %>%  # highest to smallest
	  # mutate(rank_by_normalized.strength = rank(-normalized.strength, na.last = "keep")) %>% 
	  mutate(rank_by_ANND = rank(-annd, na.last = "keep")) %>%  # Rank the 'annd' values, ignoring NA values
	  mutate(rank_by_PR = rank(-PageRank, na.last = "keep")) %>% 
	  mutate(rank_by_p.PR = rank(p.PageRank, na.last = "keep")) %>%   # smallest to highest 
	  mutate(rank_by_p.ANND = rank(p.annd, na.last = "keep")) %>%   # smallest to highest 
	  ungroup()  # Ungroup after the operation
   colnames(df)
 # [1] "signature"                   "gene"                        "PageRank"                    "PPI_cat"                    
 # [5] "EigenCentrality"             "p.PageRank"                  "rank_by_p.PR"                "rank_by_PR"                 
 # [9] "annd"                        "p.annd"                      "strength"                    "rank_by_strength"           
# [13] "normalized.strength"         "rank_by_normalized.strength" "rank_by_ANND"                "rank_by_p.ANND"   
    saveRDS(df, file='df_PAGERANK_strength_ANND.rewring.P.rds')  #!!!!!!!!!!!!!!!!!!!!!!!
    write.table(df, file='df_PAGERANK_strength_ANND.rewring.P.tsv',sep='\t', row.names=F)  #!!!!!!!!

##########################
## add the column of betweenness centrality
##########################
    df = readRDS( file='df_PAGERANK_strength_ANND.rewring.P.rds')
	
	
	betweenness_list = lapply(graph_list, function(x) betweenness(x, weights = 1/E(x)$weight))
	bc.median = lapply(betweenness_list, median) %>% unlist
	
	for(i in 1:length(betweenness_list)){
		 betweenness_list[[i]] = data.frame(BetweennessCentrality=betweenness_list[[i]], gene = names(betweenness_list[[i]])) %>%
				mutate(rank_by_BC = rank(-BetweennessCentrality, na.last = "keep"))   # Rank the 'annd' values, ignoring NA values
		 }
	df_BC = betweenness_list %>% rbindlist(., idcol=names(.))  
	colnames(df_BC)[1] <- "signature"
	df_BC$PPI_cat = lapply(df_BC$signature, function(x) unlist(strsplit(x, split="_"))[1]) %>% unlist
	 
	dim(df_BC)  # [1] 4583    5
 
 write.table(df_BC, file='df_betweeness.tsv',sep='\t', row.names=F)  #!!!!!!!!


########### betweenness centrality ############
 df_BC = read.table(file='df_betweeness.tsv',sep='\t', header=T) 
 df_BC$PPI_cat = factor(df_BC$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 				
 
	CHD = readRDS( file=here::here("examples", "Shared_Data", "CHD_Cilia_Genelist.rds"))
    df_BC$PCGC_AllCurated = df_BC$gene %in% unlist(CHD[c('Griffin2023_PCGC_AllCurated')])

	# Calculate top 5 significant genes within each box
	df5 <- df_BC %>%
	  filter(rank_by_BC <= 5 & BetweennessCentrality>0) %>%
	  ungroup()
 
 write.table(df5[, c("signature",     "gene"   , "BetweennessCentrality"  , "PPI_cat" , "rank_by_BC",
				  "PCGC_AllCurated"
				)], file='table_top5_Betweenness_perPPI.tsv', sep='\t',  row.names =FALSE, quote=FALSE) #!!!!!!!!!!!!!!

	
	subset(df5, PPI_cat=='HiGCTS')
#           signature BetweennessCentrality   gene rank_by_BC PPI_cat PCGC_AllCurated
# 66   HiGCTS_muscle                    16 FILIP1          4  HiGCTS           FALSE
# 67   HiGCTS_muscle                    51 PDLIM3          2  HiGCTS           FALSE
# 68   HiGCTS_muscle                    45  ACTA1          3  HiGCTS           FALSE
# 69   HiGCTS_muscle                     7    NES          5  HiGCTS           FALSE
# 70   HiGCTS_muscle                    61  TNNT2          1  HiGCTS           FALSE
# 71 HiGCTS_endoderm                     2  APOA1          2  HiGCTS           FALSE
# 72 HiGCTS_endoderm                     4    TTR          1  HiGCTS           FALSE
# 73       HiGCTS_CP                    55   BMP4          2  HiGCTS            TRUE
# 74       HiGCTS_CP                    73 PDGFRA          1  HiGCTS            TRUE
# 75       HiGCTS_CP                    45    VIM          3  HiGCTS           FALSE
# 76     HiGCTS_CP.1                    48   ISL1          2  HiGCTS            TRUE
# 77     HiGCTS_CP.1                   118  FGF10          1  HiGCTS           FALSE
# 78     HiGCTS_CP.1                    35 HAPLN1          4  HiGCTS           FALSE
# 79     HiGCTS_CP.1                    42   HAS2          3  HiGCTS            TRUE
# 80     HiGCTS_CP.1                    34 CITED2          5  HiGCTS            TRUE


	(df5_CHD = subset(df5, PCGC_AllCurated==TRUE))
      # signature BetweennessCentrality   gene rank_by_BC PPI_cat PCGC_AllCurated
# 16        HiG_CP                  1925 PDGFRA          3     HiG            TRUE
# 18        HiG_CP                   680   HAS2          5     HiG            TRUE
# 21         HiG_4                  1032   BMP4          3     HiG            TRUE
# 26         HiG_5                  6805  ACTC1          5     HiG            TRUE
# 28         HiG_5                 14209   MYH6          1     HiG            TRUE
# 41    HiG_muscle                   634  ACTA2          5     HiG            TRUE
# 57  HiG_endoderm                   122   CDH2          4     HiG            TRUE
# 73     HiGCTS_CP                    55   BMP4          2  HiGCTS            TRUE
# 74     HiGCTS_CP                    73 PDGFRA          1  HiGCTS            TRUE
# 76   HiGCTS_CP.1                    48   ISL1          2  HiGCTS            TRUE
# 79   HiGCTS_CP.1                    42   HAS2          3  HiGCTS            TRUE
# 80   HiGCTS_CP.1                    34 CITED2          5  HiGCTS            TRUE
# 90  CTS_endoderm                   264   GPC3          3     CTS            TRUE
# 91        CTS_CP                   227   BMP4          4     CTS            TRUE
# 92        CTS_CP                   343 PDGFRA          2     CTS            TRUE
# 100     CTS_CP.1                   477   HAS2          2     CTS            TRUE



	## reorder df$df$signature to be #E7298A #7570B3 and #E6AB02
	df_BC = rbind(subset(df_BC, PPI_cat=='CTS'),
						subset(df_BC, PPI_cat=='HiGCTS'),
						subset(df_BC, PPI_cat=='HiG')
						)
	df_BC$signature <- factor(df_BC$signature, levels = unique(df_BC$signature))
	pr <- ggplot(df_BC, aes(x = signature, y = log10(BetweennessCentrality + 1), colour = PPI_cat)) +
			  geom_boxplot(position = "dodge2") +  
			  theme(
				legend.position = 'none', 
				legend.justification = c(1, 1),  # Place legend at top-right corner
				axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
			  ) +
			  scale_color_manual(values = PPI_color_platte) + 
			  scale_size_manual(values = PPI_size_platte) +
			  geom_text(data = df5_CHD, aes(label = gene), size = 2, hjust = -0.1, vjust = 0, check_overlap = TRUE) +  # Adjust text labels
			  labs(color = "PPI cat")  # Optional: label for the color legend
	pr_repel <- ggplot(df_BC, aes(x = signature, y = log10(BetweennessCentrality + 1), colour = PPI_cat)) +
			  geom_boxplot(position = "dodge2") +  
			  theme(
				legend.position = 'none', # c(1,1)
				legend.justification = c(1, 1),  # Place legend at top-right corner
				axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
			  ) +
			  scale_color_manual(values = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")) + 
			  geom_text_repel(data = df5_CHD,  # df5
					aes(label = gene),
					size = 2,  # Adjust the size of the text labels
					box.padding = 0.5,  # Add space between the text and the data points
					point.padding = 0.5,  # Add space between the text and the points
					segment.color = 'grey50',  # Color for the line connecting the text to the points
					max.overlaps = 40,  # Max number of overlaps before labels stop being placed
					show.legend = FALSE  # Do not show text labels in the legend
				  ) +
			  #scale_x_discrete(limits = unique(df$signature)) +
			  labs(color = "PPI cat")  # Optional: label for the color legend

				
    density_bc_plot <- ggplot(df_BC, aes(x = log10(BetweennessCentrality+1), color = PPI_cat, fill = PPI_cat)) +
		  geom_density(alpha = 0.3) +  # Density lines with transparency
		  scale_color_manual(values = PPI_color_platte) +
		  scale_fill_manual(values = PPI_color_platte) +
		  scale_size_manual(values = PPI_size_platte) +
		  theme_minimal() +
		  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
		  coord_flip() +  # Flip the axes to rotate the density plot
		  labs(x = "Density of BetweennessCentrality", y = "")     # Label the x-axis and remove the y-axis label		  
		  # Add statistical comparisons using stat_compare_means manyally

	# Violin plot with statistical comparisons
	violin_wilcox = ggplot(df_BC, aes(x = PPI_cat, y = log10(BetweennessCentrality+1), color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category") +  # Label the axes
	  # Add statistical comparisons using stat_compare_means
	  stat_compare_means(
		aes(group = PPI_cat),  # Grouping by the 'PPI_cat' column
		comparisons = list( c("HiG", "HiGCTS"), c("HiGCTS", "CTS"),c("HiG", "CTS")),  # Specify comparisons
		method = "wilcox.test",  # Non-parametric test (Wilcoxon)
		label = "p.signif",  # Show significance labels (e.g., **, *, ns)
		label.x = 1.5,  # Adjust x-position of the p-value text
		size = 4  # Adjust size of the p-value text
		,tip.length =0
	  ) +
	  ggtitle('wilcox-test, all PPINs')
	
	violin_t = ggplot(df_BC, aes(x = PPI_cat, y = log10(BetweennessCentrality+1), color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
	  scale_size_manual(values = PPI_size_platte) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category") +  # Label the axes
	  # Add statistical comparisons using stat_compare_means
	  stat_compare_means(
		aes(group = PPI_cat),  # Grouping by the 'PPI_cat' column
		comparisons = list( c("HiG", "HiGCTS"), c("HiGCTS", "CTS"),c("HiG", "CTS")),  # Specify comparisons
		method = "t.test",  # Non-parametric test (Wilcoxon)
		label = "p.signif",  # Show significance labels (e.g., **, *, ns)
		label.x = 1.5,  # Adjust x-position of the p-value text
		size = 4  # Adjust size of the p-value text
		,tip.length =0
	  ) +
	  ggtitle('t-test, all PPINs')
	  
	# Caused by error in `kruskal.test.default()`:
    # ! 'x' and 'g' must have the same length 

	
	df_median = df_BC %>% group_by(signature) %>%
						summarise(bc.median = median(BetweennessCentrality, na.rm = TRUE)) %>%
						as.data.frame()
	df_median$PPI_cat = lapply(df_median$signature %>% as.vector, function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist 
	df_median$PPI_cat = factor(df_median$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 				

		## it is more manke sense to access each signature by its median betweness rank !!!!
		a = grepl('^HiG_', df_median$signature)
		b = grepl('^HiGCTS_', df_median$signature)
		c = grepl('^CTS_', df_median$signature)
    ks.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value =  0.01261  HiG vs HiGCTS
    ks.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value =   0.01261	HiG vs CTS
	ks.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value =  1	HiGCTS vs CTS
	wilcox.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.01278  HiG vs HiGCTS
    wilcox.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value =  0.01278	HiG vs CTS
	wilcox.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = NA	HiGCTS vs CTS
	t.test(df_median$bc.median[a], df_median$bc.median[b]) # p-value = 0.1945  HiG vs HiGCTS
    t.test(df_median$bc.median[a], df_median$bc.median[c]) # p-value =  0.1945	HiG vs CTS
	t.test(df_median$bc.median[b], df_median$bc.median[c]) # p-value = NA	HiGCTS vs CTS

	density_median_bc_plot <- ggplot(df_median, aes(x = log10(bc.median), color = PPI_cat, fill = PPI_cat)) +
		  geom_density(alpha = 0.3) +  # Density lines with transparency
		  scale_color_manual(values = PPI_color_platte) +  
		  scale_fill_manual(values = PPI_color_platte) +
		  theme_minimal() +
		  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
		  #coord_flip() +  # Flip the axes to rotate the density plot
		  labs(x = "Density of the median of BetweennessCentralitys per PPI", y = "")   
		# Add statistical comparisons using stat_compare_means manyally
	
	x = which(df_median$bc.median == 0)  # 5
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
	
	violin_median_bc_wilcox_ln = ggplot(df_median, aes(x = PPI_cat, y = log10(bc.median+1) , color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3, drop = FALSE) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category", y = "log10(median of BC per PPI +1)") +  # Label the axes
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
	  ggtitle('wilcox-test, median BC+1')
	
		# Combine the boxplot and density plot
	pdf(file='BetweennessCentrality_GSE1756343_v2.pdf')
    print(grid.arrange(pr_repel, density_median_bc_plot+ coord_flip() , ncol = 2, widths = c(3, 1)))
	print(grid.arrange(violin_median_bc_wilcox, pr,  nrow = 2, heights = c(1, 3)))
	print(grid.arrange(violin_wilcox, pr,  nrow = 2, heights = c(1, 3)))
	print(grid.arrange(violin_t, pr,  nrow = 2, heights = c(1, 3)))
	print(grid.arrange(violin_median_bc_wilcox, violin_median_bc_wilcox_ln, nrow=2))
    dev.off()  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


########### plot PageRank ############
    df = readRDS(file='df_PAGERANK_strength_ANND.rewring.P.rds')  #!!!!!!!!!!!!!!!!!!!!!!!
	dim(df) #[1] 4495   16

	## reorder df$df$signature to be #E7298A #7570B3 and #E6AB02
	df = rbind(subset(df, PPI_cat=='CTS'),
						subset(df, PPI_cat=='HiGCTS'),
						subset(df, PPI_cat=='HiG')
						)
						
	CHD = readRDS( file=here::here("examples", "Shared_Data", "CHD_Cilia_Genelist.rds"))
    df$PCGC_AllCurated = df$gene %in% unlist(CHD[c('Griffin2023_PCGC_AllCurated')])
	
    # Calculate top 5 significant genes within each box
	df5 <- df %>%
	  filter(rank_by_PR <= 5) %>%
	  ungroup()
    subset(df5, PPI_cat=='HiGCTS')
#  A tibble: 20 × 17
#    signature    gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR  annd p.annd strength rank_by_strength normalized.strength rank_by_normalized.s…¹
#    <chr>        <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>      <dbl> <dbl>  <dbl>    <dbl>            <dbl>               <dbl>                  <int>
#  1 HiGCTS_CP    BMP4    0.0898 HiGCTS            0.654      0.894         16            3  6.58      0    0.370                3             0.0132                       3
#  2 HiGCTS_CP    PDGF…   0.115  HiGCTS            0.878      0.902         17            2  5.59      0    0.510                2             0.0182                       2
#  3 HiGCTS_CP    RGS5    0.0546 HiGCTS            0.303      0.838         12            4  3.70      0    0.255                4             0.00909                      4
#  4 HiGCTS_CP    SFRP1   0.0520 HiGCTS            0.482      0.889         14.5          5  8.05      0    0.216                7             0.00770                      7
#  5 HiGCTS_CP    VIM     0.133  HiGCTS            1          0.884         13            1  6.40      0    0.583                1             0.0208                       1
#  6 HiGCTS_CP.1  FGF10   0.142  HiGCTS            1          0.866         12            1  5.33      0    0.532                1             0.0205                       1
#  7 HiGCTS_CP.1  HAND2   0.0718 HiGCTS            0.687      0.899         19            5  7.60      0    0.275                5             0.0106                       5
#  8 HiGCTS_CP.1  HAPL…   0.0782 HiGCTS            0.547      0.857          8            4  4.59      0    0.282                4             0.0108                       4
#  9 HiGCTS_CP.1  HAS2    0.0924 HiGCTS            0.853      0.867         13            2  6.29      0    0.368                2             0.0141                       2
# 10 HiGCTS_CP.1  ISL1    0.0794 HiGCTS            0.733      0.883         16            3  7.39      0    0.294                3             0.0113                       3
# 11 HiGCTS_endo… APOA1   0.224  HiGCTS            1          0.559          4            1  4.11      0    2.15                 1             0.358                        1
# 12 HiGCTS_endo… APOA2   0.207  HiGCTS            0.927      0.544          2            2  4.13      0    1.96                 3             0.326                        3
# 13 HiGCTS_endo… APOC1   0.205  HiGCTS            0.958      0.596          5            3  4.07      0    1.97                 2             0.329                        2
# 14 HiGCTS_endo… APOE    0.190  HiGCTS            0.873      0.55           3            4  4.11      0    1.79                 4             0.299                        4
# 15 HiGCTS_endo… TTR     0.122  HiGCTS            0.421      0.172          1            5  3.91      0    0.857                5             0.143                        5
# 16 HiGCTS_musc… ACTA1   0.130  HiGCTS            0.965      0.278          2            1  9.18      0    2.15                 1             0.113                        1
# 17 HiGCTS_musc… ACTC1   0.0937 HiGCTS            0.794      0.511         11            3  9.80      0    1.60                 3             0.0842                       3
# 18 HiGCTS_musc… ANKR…   0.0685 HiGCTS            0.574      0.544         12            5  9.80      0    1.14                 5             0.0600                       5
# 19 HiGCTS_musc… CNN1    0.0818 HiGCTS            0.641      0.507         10            4 10.0       0    1.36                 4             0.0715                       4
# 20 HiGCTS_musc… TNNT2   0.119  HiGCTS            1          0.56          15            2 10.3       0    2.07                 2             0.109                        2
# # ℹ abbreviated name: ¹​rank_by_normalized.strength
# # ℹ 3 more variables: rank_by_ANND <dbl>, rank_by_p.ANND <dbl>, PCGC_AllCurated <lgl>

# #   PCGC_AllCurated <lgl>

	write.table(df5[, c("signature",     "gene"   , "PageRank"  , "PPI_cat" , "rank_by_PR",
				 "normalized.strength"  , "rank_by_normalized.strength", "PCGC_AllCurated"
				)], file='table_top5_PageRank_perPPI.tsv', sep='\t',  row.names =FALSE, quote=FALSE) #!!!!!!!!!!!!!!

	df5_CHD = subset(df5, PCGC_AllCurated==TRUE)
	dim(df5)  #[1] 105   17
    dim(df5_CHD)  # [1] 21  17
	df5_CHD %>% as.data.frame
#         signature   gene    PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR       annd p.annd  strength rank_by_strength normalized.strength
# 1         CTS_CP   BMP4 0.041700069     CTS     0.572207811      0.891         51.0          3  11.692485      0 0.4803492                3         0.006765481
# 2         CTS_CP PDGFRA 0.056173592     CTS     0.768558299      0.876         50.0          2  10.088557      0 0.6704548                2         0.009443025
# 3       CTS_CP.1  HAND2 0.035485514     CTS     0.006987591      0.876         55.0          5  14.394795      0 0.4644797                4         0.006111575
# 4       CTS_CP.1   HAS2 0.047658130     CTS     0.005561831      0.867         49.0          2   9.526100      0 0.5591797                2         0.007357628
# 5       CTS_CP.1   ISL1 0.036617407     CTS     0.007291816      0.872         51.0          4  13.505401      0 0.4615943                5         0.006073609
# 6     CTS_muscle  ACTC1 0.046885181     CTS     0.767320162      0.789         56.0          3  19.892036      0 3.1027595                3         0.051712658
# 7      HiGCTS_CP   BMP4 0.089757747  HiGCTS     0.653888312      0.894         16.0          3   6.581750      0 0.3699930                3         0.013214034
# 8      HiGCTS_CP PDGFRA 0.114700169  HiGCTS     0.877634525      0.902         17.0          2   5.585933      0 0.5104018                2         0.018228637
# 9    HiGCTS_CP.1  HAND2 0.071814955  HiGCTS     0.686677772      0.899         19.0          5   7.603309      0 0.2748660                5         0.010571771
# 10   HiGCTS_CP.1   HAS2 0.092419408  HiGCTS     0.853227308      0.867         13.0          2   6.293260      0 0.3675556                2         0.014136753
# 11   HiGCTS_CP.1   ISL1 0.079388353  HiGCTS     0.733376893      0.883         16.0          3   7.386761      0 0.2940990                3         0.011311501
# 12 HiGCTS_muscle  ACTC1 0.093733482  HiGCTS     0.793547722      0.511         11.0          3   9.798229      0 1.5992411                3         0.084170584
# 13        HiG_10  RPS19 0.023161038     HiG     0.975281786      0.962         63.5          4  47.496526      0 6.5246634                2         0.094560339
# 14        HiG_12   MYH7 0.086967350     HiG     0.000000000      0.874         15.0          4   7.185046      0 0.1158520                4         0.005516760
# 15         HiG_4   BMP4 0.019128183     HiG     0.003072583      0.000          8.5          2  13.155925      0 0.8688782               65         0.008435711
# 16         HiG_5  ACTC1 0.018664896     HiG     0.999615338      0.880        324.0          2  65.900958      0 8.6298816                1         0.025532194
# 17         HiG_9  ACTC1 0.005650393     HiG     0.945915426      0.706       1291.5          3 125.775415      0 1.9932407                2         0.001359646
# 18         HiG_9   MYH6 0.005577847     HiG     1.000000000      0.749       1350.5          4 127.432256      0 2.0423901                1         0.001393172
# 19        HiG_CP   BMP4 0.023460806     HiG     0.005476214      0.000          3.0          3  22.850643      0 1.2698019               38         0.009693144
# 20        HiG_CP PDGFRA 0.024419816     HiG     0.002431687      0.000          3.0          1  21.982299      0 1.2525291               39         0.009561291
# 21    HiG_muscle  ACTA2 0.023564100     HiG     0.864854224      0.589         84.0          4  52.410521      0 5.1837620                4         0.044305658
#    rank_by_normalized.strength rank_by_ANND rank_by_p.ANND PCGC_AllCurated
# 1                            3           14           26.0            TRUE
# 2                            2           19           26.0            TRUE
# 3                            4            5           31.0            TRUE
# 4                            2           24           31.0            TRUE
# 5                            5            8           31.0            TRUE
# 6                            3           23           29.0            TRUE
# 7                            3            7            9.5            TRUE
# 8                            2            9            9.5            TRUE
# 9                            5            5           10.5            TRUE
# 10                           2           10           10.5            TRUE
# 11                           3            6           10.5            TRUE
# 12                           3           13            9.5            TRUE
# 13                           2            8           35.5            TRUE
# 14                           4           16           10.5            TRUE
# 15                          65           96           51.5            TRUE
# 16                           1           57          169.0            TRUE
# 17                           2          366          733.5            TRUE
# 18                           1          351          733.5            TRUE
# 19                          38           78           65.0            TRUE
# 20                          39           84           65.0            TRUE
# 21                           4           49           59.5            TRUE

	
	df$signature <- factor(df$signature, levels = levels(df_BC$signature))
    pr <- ggplot(df, aes(x = signature,y = PageRank, colour = PPI_cat)) +
			  geom_boxplot(show.legend = TRUE) +  # Enable legend for the boxplot
			  scale_color_manual(values = PPI_color_platte) +
			  scale_size_manual(values = PPI_size_platte) +
			  geom_text(data=df5_CHD, aes(label = gene),  # data=df5
						size = 2,  # Adjust the size of the text labels
						hjust = -0.1, vjust = 0, 
						check_overlap = TRUE) +  # Avoid text overlap
			  theme(legend.position = 'none', #c(1, 1), 
					legend.justification = c(1, 1),  # Place legend at top-right corner
					axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
			  scale_x_discrete(limits = unique(df$signature)) +
			  labs(color = "PPI cat")  # Optional: label for the color legend
	pr_repel <- ggplot(df, aes(x = signature, y = PageRank, colour = PPI_cat)) +
			  geom_boxplot(show.legend = TRUE) +  # Enable legend for the boxplot
			  scale_color_manual(values = PPI_color_platte) +
			  scale_size_manual(values = PPI_size_platte) +
			  geom_text_repel(data = df5_CHD,  # df5
					aes(label = gene),
					size = 2,  # Adjust the size of the text labels
					box.padding = 0.5,  # Add space between the text and the data points
					point.padding = 0.5,  # Add space between the text and the points
					segment.color = 'grey50',  # Color for the line connecting the text to the points
					max.overlaps = 20,  # Max number of overlaps before labels stop being placed
					show.legend = FALSE  # Do not show text labels in the legend
				  ) +
			theme(legend.position ='none', 
					#legend.justification = 'none', #c(1, 1),  # Place legend at top-right corner
					axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
			  #scale_x_discrete(limits = unique(df$signature)) +
			  labs(color = "PPI cat")  # Optional: label for the color legend

				
    density_page_plot <- ggplot(df, aes(x = PageRank, color = PPI_cat, fill = PPI_cat)) +
		  geom_density(alpha = 0.3) +  # Density lines with transparency
		  scale_color_manual(values = PPI_color_platte) +
		  scale_fill_manual(values = PPI_color_platte) +
		  theme_minimal() +
		  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
		  coord_flip() +  # Flip the axes to rotate the density plot
		  labs(x = "Density of PageRank", y = "")     # Label the x-axis and remove the y-axis label		  
		  # Add statistical comparisons using stat_compare_means manyally
	# Violin plot with statistical comparisons
	violin_wilcox = ggplot(df, aes(x = PPI_cat, y = PageRank, color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category", y = "PageRank") +  # Label the axes
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
	  ggtitle('wilcox test, all PPINs ')

    violin_t = ggplot(df, aes(x = PPI_cat, y = PageRank, color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category", y = "PageRank") +  # Label the axes
	  # Add statistical comparisons using stat_compare_means
	  stat_compare_means(
		aes(group = PPI_cat),  # Grouping by the 'PPI_cat' column
		comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")),  # Specify comparisons
		method = "t.test",  # Non-parametric test (Wilcoxon)
		label = "p.signif",  # Show significance labels (e.g., **, *, ns)
		label.x = 1.5,  # Adjust x-position of the p-value text
		size = 4  # Adjust size of the p-value text
		,tip.length =0
	  ) + 
	  ggtitle('t test, all PPINs ')
 	  
	## it is more manke sense to access each signature by its median bc rank !!!!
	# pg.median = lapply(page, median) %>% unlist
	pg.median = df %>% group_by(signature) %>%
					summarise(median_PageRank = median(PageRank, na.rm = TRUE))
    	a = grepl('^HiG_', pg.median$signature)
		b = grepl('^HiGCTS_', pg.median$signature)
		c = grepl('^CTS_', pg.median$signature)
	ks.test(pg.median[a,]$median_PageRank, pg.median[b,]$median_PageRank) # p-value =  0.0008403  HiG vs HiGCTS
    ks.test(pg.median[a,]$median_PageRank, pg.median[c,]$median_PageRank) # p-value =  0.01261	HiG vs CTS
	ks.test(pg.median[b,]$median_PageRank, pg.median[c,]$median_PageRank) # p-value =  0.02857	HiGCTS vs CTS
	wilcox.test(pg.median[a,]$median_PageRank, pg.median[b,]$median_PageRank) # p-value = 0.0008403  HiG vs HiGCTS
    wilcox.test(pg.median[a,]$median_PageRank, pg.median[c,]$median_PageRank) # p-value =  0.04454 	HiG vs CTS
	wilcox.test(pg.median[b,]$median_PageRank, pg.median[c,]$median_PageRank) # p-value =  0.02857	HiGCTS vs CTS
	
	df_median = data.frame(pg.median = pg.median$median_PageRank, 
					PPI_cat = lapply(pg.median$signature %>% as.vector, function(x) unlist(strsplit(x, split='_', fixed=T))[1]) %>% unlist)
	df_median$PPI_cat = factor(df_median$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 				
	density_median_page_plot <- ggplot(df_median, aes(x = pg.median, color = PPI_cat, fill = PPI_cat)) +
		  geom_density(alpha = 0.3) +  # Density lines with transparency
		  scale_color_manual(values = PPI_color_platte) +
		  scale_fill_manual(values = PPI_color_platte) +
		  theme_minimal() +
		  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
		  #coord_flip() +  # Flip the axes to rotate the density plot
		  labs(x = "Density of the median of PageRanks per PPI", y = "")   
		# Add statistical comparisons using stat_compare_means manyally
		
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
	  
 
	# Combine the boxplot and density plot
	pdf(file='PageRank_GSE1756343_v2.pdf')
    print(grid.arrange(pr, density_median_page_plot+ coord_flip() , ncol = 2, widths = c(3, 1)))
	print(grid.arrange(violin_median_page_wilcox, pr,  nrow = 2, heights = c(1, 3)))
	print(grid.arrange(violin_wilcox, pr,  nrow = 2, heights = c(1, 3)))
	print(grid.arrange(violin_t, pr,  nrow = 2, heights = c(1, 3)))
    dev.off()  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

########### plot ANND (NOT USED   ) ############
 {
	df = readRDS(file='df_PAGERANK_strength_ANND.rewring.P.rds')  #!!!!!!!!!!!!!!!!!!!!!!!
   
    df$label=df$gene
    # df$label[which(df$signature!='HiGCTS_CP.1')] = ''
    # df$label[which(round(df$p.annd,2)>0.05)] = ''
    # df$label[which(is.na(df$annd))] = ''
    subset(df, signature=='HiGCTS_CP.1')
    
	## reorder df$df$signature to be #E7298A #7570B3 and #E6AB02
	df = rbind(subset(df, PPI_cat=='CTS'),
						subset(df, PPI_cat=='HiGCTS'),
						subset(df, PPI_cat=='HiG')
						)
	df$signature <- factor(df$signature, levels = unique(df$signature))
	
  	# Step 1: Filter top 5 genes by annd within each signature
	top_genes <- df %>%
	  group_by(signature) %>%
	  arrange(desc(annd)) %>%   # Sort in descending order of annd
	  slice_head(n = 5)  # Take the top 5 rows for each signature
	# subset the CHD genes within top 5 
	CHD = readRDS( file=here::here("examples", "Shared_Data", "CHD_Cilia_Genelist.rds"))
 
	top_genes_CHD = subset(top_genes, gene %in% unlist(CHD[c('Griffin2023_PCGC_AllCurated')]))
	dim(top_genes)  #[1] 105   16
    dim(top_genes_CHD)  # [1] 5  16
	top_genes_CHD
# # A tibble: 5 × 17
# # Groups:   signature [5]
#   signature   gene  PageRank PPI_cat EigenCentrality p.PageRank rank_by_p.PR rank_by_PR  annd p.annd strength rank_by_strength
#   <fct>       <chr>    <dbl> <fct>             <dbl>      <dbl>        <dbl>      <dbl> <dbl>  <dbl>    <dbl>            <dbl>
# 1 CTS_CP      HEY1   0.00934 CTS             0.111        0.859           44         35 14.0       0   0.0819               27
# 2 CTS_CP.1    HAND2  0.0355  CTS             0.00699      0.876           55          5 14.4       0   0.464                 4
# 3 HiGCTS_CP   HEY1   0.0208  HiGCTS          0.132        0.816            3         16  8.33      0   0.0627               16
# 4 HiGCTS_CP.1 HAND2  0.0718  HiGCTS          0.687        0.899           19          5  7.60      0   0.275                 5
# 5 HiG_12      MGP    0.0357  HiG             0            0.821            2         13 11.9       0   0.0246               13
# # ℹ 5 more variables: normalized.strength <dbl>, rank_by_normalized.strength <int>, rank_by_ANND <dbl>, rank_by_p.ANND <dbl>, label <chr>

	# Step 2: Create ggplot with boxplot and labels for top 5 genes
	pr = ggplot(df[!is.na(df$annd),], aes(x = signature, y = annd, colour = PPI_cat)) + 
	  geom_boxplot(show.legend = TRUE) +  # Enable legend for the boxplot
	  scale_color_manual(values = PPI_color_platte) +
	  # Use ggrepel to avoid overlap and label top 5 genes based on annd
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
		legend.justification = c(0, 1),  # Place legend at top-right corner
		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
	  ) +
	  #scale_x_discrete(limits = unique(df$signature)) + # Ensure x-axis respects the order of 'signature'
	  labs(color = "PPI cat")   # Optional: label for the color legend
    
	## it is more manke sense to access each signature by its median page rank !!!!
	annd.median = lapply(annd_observed, median, na.rm=TRUE) %>% unlist
	
	A = which(grepl('^CTS_', names(annd.median)))
	B = which(grepl('^HiGCTS_', names(annd.median)))
	C = which(grepl('^HiG_', names(annd.median)))
    ks.test(annd.median[C], annd.median[B]) # p-value =  0.004202  HiG vs HiGCTS
    ks.test(annd.median[C], annd.median[A]) # p-value =  0.0126	HiG vs CTS
	ks.test(annd.median[B], annd.median[A]) # p-value =  0.2286	HiGCTS vs CTS
	wilcox.test(annd.median[C], annd.median[B]) # p-value = 0.001681 HiG vs HiGCTS
    wilcox.test(annd.median[C], annd.median[A]) # p-value =  0.01008	HiG vs CTS
	wilcox.test(annd.median[B], annd.median[A]) # p-value =  0.1143	HiGCTS vs CTS
	
	df_median = data.frame(annd.median = annd.median, 
					PPI_cat = lapply(names(annd.median), function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist)
	df_median$PPI_cat = factor(df_median$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 				
	density_median_annd_wilcox <- ggplot(df_median, aes(x = annd.median, color = PPI_cat, fill = PPI_cat)) +
		  geom_density(alpha = 0.3) +  # Density lines with transparency
		  scale_color_manual(values = PPI_color_platte) +
		  scale_fill_manual(values = PPI_color_platte) +
		  theme_minimal() +
		  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
		  #coord_flip() +  # Flip the axes to rotate the density plot
		  labs(x = "Density of the median of ANND per PPI", y = "")   
		# Add statistical comparisons using stat_compare_means manyally
	violin_median_annd_plot = ggplot(df_median, aes(x = PPI_cat, y = annd.median, color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
	  theme_minimal() +
	  theme(legend.position = "none") +  
	  labs(x = "PPI category", y = "median of ANND per PPI") +  # Label the axes
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
	  ggtitle('wilcox test, median ANND')

	# Combine the boxplot and density plot
	pdf(file='annd_GSE175634_v2.pdf')
    print(grid.arrange(pr, density_median_annd_wilcox+ coord_flip() , ncol = 2, widths = c(3, 1)))
	print(grid.arrange(violin_median_annd_plot, pr,  nrow = 2, heights = c(1, 3)))
	# print(violin_plot + ggtitle('wilcox test of all PPI'))

    dev.off()  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 }   
  	
 ############# plot normalized.strength (NO between-category difference! ) ###########
 {
 df = readRDS(file='df_PAGERANK_strength_ANND.rewring.P.rds')  #!!!!!!!!!!!!!!!!!!!!!!!
 	## reorder df$df$signature to be #E7298A #7570B3 and #E6AB02
	df = rbind(subset(df, PPI_cat=='CTS'),
						subset(df, PPI_cat=='HiGCTS'),
						subset(df, PPI_cat=='HiG')
						)
	df$signature <- factor(df$signature, levels = unique(df$signature))

    df$label=df$gene
    subset(df, signature=='HiGCTS_CP.1')
      
	CHD = readRDS( file=here::here("examples", "Shared_Data", "CHD_Cilia_Genelist.rds"))
	df$PCGC_AllCurated = df$gene %in% unlist(CHD[c('Griffin2023_PCGC_AllCurated')])

	# Step 1: Filter top 5 genes by normalized.strength within each signature
	top_genes <- df %>%
	  group_by(signature) %>%
	  arrange(desc(normalized.strength)) %>%   # Sort in descending order of normalized.strength
	  slice_head(n = 5)  # Take the top 5 rows for each signature
	  
	write.table(top_genes[, c("signature",     "gene"   , "normalized.strength"  , "PPI_cat" , "rank_by_normalized.strength",
				  "PCGC_AllCurated"
				)], file='table_top5_strength_perPPI.tsv', sep='\t',  row.names =FALSE, quote=FALSE) #!!!!!!!!!!!!!!
  
	    
	# subset the CHD genes within top 5 
	top_genes_CHD = subset(top_genes, PCGC_AllCurated==TRUE)
	dim(top_genes)  #[1] 105   18
    dim(top_genes_CHD)  # [1] 13  18
	
	# Violin plot with statistical comparisons
	violin_wilcox = ggplot(df, aes(x = PPI_cat, y = log10(normalized.strength), color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category", y = "log10 normalized.strength") +  # Label the axes
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
	  ggtitle('wilcox test, all PPINs ')

    violin_t = ggplot(df, aes(x = PPI_cat, y = log10(normalized.strength), color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category", y = "normalized.strength") +  # Label the axes
	  # Add statistical comparisons using stat_compare_means
	  stat_compare_means(
		aes(group = PPI_cat),  # Grouping by the 'PPI_cat' column
		comparisons = list(c("HiG", "CTS"), c("HiG", "HiGCTS"), c("HiGCTS", "CTS")),  # Specify comparisons
		method = "t.test",  # Non-parametric test (Wilcoxon)
		label = "p.signif",  # Show significance labels (e.g., **, *, ns)
		label.x = 1.5,  # Adjust x-position of the p-value text
		size = 4  # Adjust size of the p-value text
		,tip.length =0
	  ) + 
	  ggtitle('t test, all PPINs ')

	# Step 2: Create ggplot with boxplot and labels for top 5 genes
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
		legend.position = 'none', # c(1, 1), 
		legend.justification = c(0, 1),  # Place legend at top-right corner
		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
	  ) +
	  scale_x_discrete(limits = unique(df$signature)) + # Ensure x-axis respects the order of 'signature'
	  labs(color = "PPI cat")   # Optional: label for the color legend
    
	## it is more make sense to access each signature by its median strength !!!!
	df_median = df %>% group_by(signature) %>%
					summarise(median_normalized_strength = median(normalized.strength, na.rm = TRUE))
	
	A = which(grepl('^CTS_', df_median$signature))
	B = which(grepl('^HiGCTS_', df_median$signature))
	C = which(grepl('^HiG_', df_median$signature))
    ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value =  0.23  HiG vs HiGCTS
    ks.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value =  0.94	HiG vs CTS
	ks.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value =  0.77	HiGCTS vs CTS
	wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[B]) # p-value = 0.163  HiG vs HiGCTS
    wilcox.test(df_median$median_normalized_strength[C], df_median$median_normalized_strength[A]) # p-value =  1	HiG vs CTS
	wilcox.test(df_median$median_normalized_strength[B], df_median$median_normalized_strength[A]) # p-value =  0.3429	HiGCTS vs CTS

	df_median$PPI_cat = lapply(df_median$signature %>% as.vector(), function(x) unlist(strsplit(x, split='_'))[1]) %>% unlist 
	df_median$PPI_cat = factor(df_median$PPI_cat,levels=c('CTS', 'HiGCTS', 'HiG')) 				
	
	density_median_normalized.strength_plot <- ggplot(df_median, aes(x = log10(median_normalized_strength), 
						color = PPI_cat, fill = PPI_cat)) +
		  geom_density(alpha = 0.3) +  # Density lines with transparency
		  scale_color_manual(values = PPI_color_platte) +
		  scale_fill_manual(values = PPI_color_platte) +
		  theme_minimal() +
		  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank()) +
		  #coord_flip() +  # Flip the axes to rotate the density plot
		  labs(x = "Density of the median of normalzied node strength per PPI", y = "")   
		# Add statistical comparisons using stat_compare_means manyally
	
	violin_median_normalized.strength_wilcox = ggplot(df_median, aes(x = PPI_cat, 
						y = log10(median_normalized_strength), color = PPI_cat, fill = PPI_cat)) +
	  geom_violin(alpha = 0.3) +  # Violin plot with transparency
	  scale_color_manual(values = PPI_color_platte) +
	  scale_fill_manual(values = PPI_color_platte) +
	  theme_minimal() +
	  theme(legend.position = "none") + #, axis.text.y = element_blank(), axis.title.y = element_blank()) +
	  labs(x = "PPI category", y = "log10. median of normalized node strength per PPI") +  # Label the axes
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
	  ggtitle('wilcox, medina nr_strength')

	# Combine the boxplot and density plot
	pdf(file='normalized.node.strength_GSE175634_v2.pdf')
    print(grid.arrange(pr, density_median_normalized.strength_plot+ coord_flip() , ncol = 2, widths = c(3, 1)))
	print(grid.arrange(violin_median_normalized.strength_wilcox, pr,  nrow = 2, heights = c(1, 3)))
	print(grid.arrange(violin_wilcox, pr,  nrow = 2, heights = c(1, 3)))
	print(grid.arrange(violin_t, pr,  nrow = 2, heights = c(1, 3)))
	# grid.arrange(pr, density_normalized.strength_plot, ncol = 2, widths = c(3, 1))
 
    dev.off()  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
 } 	
     
    
    
    
    # number of significantly high annd per PPI_cats
    n.annd.high <- lapply(names(graph_list), function(x) nrow(subset(df, signature==x & round(p.annd,2)<=0.05))) %>% unlist
    names(n.annd.high) = names(graph_list)
    n.annd.high 
         # HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4           HiG_5           HiG_6           HiG_7 
            # 325            1154             240             129             102             337             102             104 
     # HiG_muscle           HiG_9          HiG_10    HiG_endoderm          HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP 
            # 118            1466              70              37              20              18               6              18 
    # HiGCTS_CP.1      CTS_muscle    CTS_endoderm          CTS_CP        CTS_CP.1 
             # 20              57              60              51              61 
    
    df_compare = data.frame(signature=names(graph_list), 
                  n_sig.pagerank = n.pr.high, 
                  n_sig.annd = n.annd.high) 
    df_compare$PPI_cat = lapply(df_compare$signature, function(x) unlist(strsplit(x , split='_'))[1]) %>% unlist %>%
						 factor(.,levels=c('CTS', 'HiGCTS', 'HiG')) 
    pdf(file='n.sig.pageRank_vs_n.sig.annd.pdf')
    print(ggplot(df_compare,aes(x= n_sig.pagerank, y=n_sig.annd)) +
      geom_point(aes(shape= PPI_cat, colour=PPI_cat), show.legend = FALSE) +
	  scale_color_manual(values = PPI_color_platte) +
      geom_text_repel(aes(label=signature), hjust= -0.1, vjust=0) +
      theme(legend.position=c(0, 0), legend.justification=c(0, 0)))
   dev.off()
    
    
    
###########################################################################################################################################
## Given a transitional state, CTS&HiG genes exhibit higher betweenness centrality in the CTS-derived network and the HiG-derived network  
###########################################################################################################################################
 bc = read.table(file='df_betweeness.tsv', header=TRUE) 
 dim(bc)  # [1] 4583    5
 colnames(bc)
 # [1] "signature"             "BetweennessCentrality" "gene"                 
# [4] "rank_by_BC"            "PPI_cat"     
 ## find out the cluster with CTSHiG 
 x = grep('HiGCTS_', bc$signature, value=TRUE) %>% unique
 CTS = lapply(x, function(x) unlist(strsplit(x, split='_'))[2]) %>% unlist
 y = grepl('.',CTS, fixed=TRUE)
 if(any(y)) CT = CTS[!y]  else CT = CTS

 CTS
 #["muscle"   "endoderm" "CP"       "CP.1" 
 CT
 #"muscle"   "endoderm" "CP" 

plot_CTS_bc = plot_HiG_bc = list()   
 ## compare the CTS&HiG genes to non-hiG genes per CTS-derived PPI
for(id in CTS){
   if(grepl('.', id, fixed=TRUE)) id2 = unlist(strsplit(id, split='.', fixed=TRUE))[1] else id2=id
   CTS_PPI = subset(bc, signature == paste('CTS',id, sep='_'))
   HiG_PPI = subset(bc, signature == paste('HiG',id2, sep='_'))
   CTS_PPI$isHiG = factor(CTS_PPI$gene %in%  HiG_PPI$gene , levels = c('FALSE', 'TRUE'))
   HiG_PPI$isCTS = factor(HiG_PPI$gene %in%  CTS_PPI$gene , levels = c('FALSE', 'TRUE'))
   
    ############## ranked by betweenness 
   # Compute exact p-value
	pval <- wilcox.test(BetweennessCentrality ~ isHiG, data = CTS_PPI)$p.value
	
	# Reorder gene factor levels by BetweennessCentrality (high to low)
	CTS_PPI <- CTS_PPI %>%
		mutate(gene = factor(gene, levels = gene[order(-BetweennessCentrality)]))

	plot_CTS_bc[[id]] = ggplot(CTS_PPI, aes(x = gene, y = BetweennessCentrality, fill = isHiG)) +
	  geom_boxplot(aes(group = isHiG), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) +  # Add boxplot first
	  geom_point(aes(color = isHiG), size = 3) +
	  scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E7298A")) +
	  geom_text_repel(aes(label = gene, color=isHiG), hjust = -0.1, vjust = 0) +
	  theme(
		legend.position = c(0, 0),
		legend.justification = c(1, 1),
		axis.text.x = element_blank(),  # Remove x-axis labels
		axis.ticks.x = element_blank()
	  ) +
	  # stat_compare_means(  doesn't work becasue the x-asix os genes rather than isHiG grouping
		# aes(group = isHiG),method = "wilcox.test",label = "p.signif",label.x = 1.5,size = 4) +
	 annotate("text", x =  (length(unique(CTS_PPI$gene)) + 1) / 2, 
			y = max(CTS_PPI$BetweennessCentrality) * 0.8, 
			label = paste0('wilcox p = ', signif(pval, 2), ' F,T: ', table(CTS_PPI$isHiG) %>% toString()), size = 4)  
	  # ggtitle(paste0('CTS genes of ' ,id2, ' Wilcox test ', pval_label))
	  
	  
	  # Compute exact p-value
	pval <- wilcox.test(BetweennessCentrality ~ isCTS, data = HiG_PPI)$p.value
	
	# Reorder gene factor levels by BetweennessCentrality (high to low)
	HiG_PPI <- HiG_PPI %>%
		mutate(gene = factor(gene, levels = gene[order(-BetweennessCentrality)]))

	plot_HiG_bc[[id]] = ggplot(HiG_PPI, aes(x = gene, y = BetweennessCentrality, fill = isCTS)) +
	  geom_boxplot(aes(group = isCTS), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) +  # Add boxplot first
	  geom_point(aes(color = isCTS), size = 3) +
	  scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E6AB02")) +
	  geom_text_repel(aes(label = gene, color=isCTS), hjust = -0.1, vjust = 0) +
	  theme(
		legend.position = c(0, 0),
		legend.justification = c(1, 1),
		axis.text.x = element_blank(),  # Remove x-axis labels
		axis.ticks.x = element_blank()
	  ) +
	   annotate("text", x = (length(unique(HiG_PPI$gene)) + 1) / 2, 
			y = max(HiG_PPI$BetweennessCentrality) * 0.8, 
			label = paste0('wilcox p = ', signif(pval, 2), ' F,T: ', table(HiG_PPI$isCTS) %>% toString()), size = 4)  
	   #ggtitle(paste0('CTS genes of ' ,id2, ' Wilcox test ', pval_label))

}   
 #####
 
 df = readRDS(file='df_PAGERANK_strength_ANND.rewring.P.rds')  #!!!!!!!!!!!!!!!!!!!!!!!
 colnames(df)
 # [1] "signature"                 "gene"                     
 # [3] "PageRank"                  "PPI_cat"                  
 # [5] "EigenCentrality"           "p.PageRank"               
 # [7] "rank_by_p.PR"              "rank_by_PR"               
 # [9] "annd"                      "p.annd"                   
# [11] "normalized.strength"         "rank_by_normalized.strength"
# [13] "strength"                    "rank_by_strength"           
# [15] "rank_by_ANND"              "rank_by_p.ANND"      

plot_CTS_pr = plot_HiG_pr = list()   

## compare the CTS&HiG genes to non-hiG genes per CTS-derived PPI
for(id in CTS){
   if(grepl('.', id, fixed=TRUE)) id2 = unlist(strsplit(id, split='.', fixed=TRUE))[1] else id2=id
   CTS_PPI = subset(df, signature == paste('CTS',id, sep='_'))
   HiG_PPI = subset(df, signature == paste('HiG',id2, sep='_'))
   CTS_PPI$isHiG = factor(CTS_PPI$gene %in%  HiG_PPI$gene , levels = c('FALSE', 'TRUE'))
   HiG_PPI$isCTS = factor(HiG_PPI$gene %in%  CTS_PPI$gene , levels = c('FALSE', 'TRUE')) 
   
   ############## ranked by pageRank 
   # Compute exact p-value
	pval <- wilcox.test(PageRank ~ isHiG, data = CTS_PPI)$p.value
	
	# Reorder gene factor levels by PageRank (high to low)
	CTS_PPI <- CTS_PPI %>%
		mutate(gene = factor(gene, levels = gene[order(-PageRank)]))

	plot_CTS_pr[[id]] = ggplot(CTS_PPI, aes(x = gene, y = PageRank, fill = isHiG)) +
	  geom_boxplot(aes(group = isHiG), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) +  # Add boxplot first
	  geom_point(aes(color = isHiG), size = 3) +
	  scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E7298A")) +
	  geom_text_repel(aes(label = gene, color=isHiG), hjust = -0.1, vjust = 0) +
	  theme(
		legend.position = c(0, 0),
		legend.justification = c(1, 1),
		axis.text.x = element_blank(),  # Remove x-axis labels
		axis.ticks.x = element_blank()
	  ) +
	  # stat_compare_means(  doesn't work becasue the x-asix os genes rather than isHiG grouping
		# aes(group = isHiG),method = "wilcox.test",label = "p.signif",label.x = 1.5,size = 4) +
	 annotate("text", x =  (length(unique(CTS_PPI$gene)) + 1) / 2, 
			y = max(CTS_PPI$PageRank) * 0.8, 
			label = paste0('wilcox p = ', signif(pval, 2), ' F,T: ', table(CTS_PPI$isHiG) %>% toString()), size = 4)  
	  # ggtitle(paste0('CTS genes of ' ,id2, ' Wilcox test ', pval_label))
	  
	  
	  # Compute exact p-value
	pval <- wilcox.test(PageRank ~ isCTS, data = HiG_PPI)$p.value
	
	# Reorder gene factor levels by PageRank (high to low)
	HiG_PPI <- HiG_PPI %>%
		mutate(gene = factor(gene, levels = gene[order(-PageRank)]))

	plot_HiG_pr[[id]] = ggplot(HiG_PPI, aes(x = gene, y = PageRank, fill = isCTS)) +
	  geom_boxplot(aes(group = isCTS), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) +  # Add boxplot first
	  geom_point(aes(color = isCTS), size = 3) +
	  scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E6AB02")) +
	  geom_text_repel(aes(label = gene, color=isCTS), hjust = -0.1, vjust = 0) +
	  theme(
		legend.position = c(0, 0),
		legend.justification = c(1, 1),
		axis.text.x = element_blank(),  # Remove x-axis labels
		axis.ticks.x = element_blank()
	  ) +
	   annotate("text", x = (length(unique(HiG_PPI$gene)) + 1) / 2, 
			y = max(HiG_PPI$PageRank) * 0.8, 
			label = paste0('wilcox p = ', signif(pval, 2), ' F,T: ', table(HiG_PPI$isCTS) %>% toString()), size = 4)  
	   #ggtitle(paste0('CTS genes of ' ,id2, ' Wilcox test ', pval_label))
  
}

plot_CTS_annd = plot_HiG_annd = list()   

## compare the CTS&HiG genes to non-hiG genes per CTS-derived PPI
for(id in CTS){
   if(grepl('.', id, fixed=TRUE)) id2 = unlist(strsplit(id, split='.', fixed=TRUE))[1] else id2=id
   CTS_PPI = subset(df, signature == paste('CTS',id, sep='_'))
   HiG_PPI = subset(df, signature == paste('HiG',id2, sep='_'))
   CTS_PPI$isHiG = factor(CTS_PPI$gene %in%  HiG_PPI$gene , levels = c('FALSE', 'TRUE'))
   HiG_PPI$isCTS = factor(HiG_PPI$gene %in%  CTS_PPI$gene , levels = c('FALSE', 'TRUE'))
  
   ############## ranked by annd 
   # Compute exact p-value
	pval <- wilcox.test(annd ~ isHiG, data = CTS_PPI)$p.value
	
	# Reorder gene factor levels by annd (high to low)
	CTS_PPI <- CTS_PPI %>%
		mutate(gene = factor(gene, levels = gene[order(-annd)]))

	plot_CTS_annd[[id]] = ggplot(CTS_PPI, aes(x = gene, y = annd, fill = isHiG)) +
	  geom_boxplot(aes(group = isHiG), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) +  # Add boxplot first
	  geom_point(aes(color = isHiG), size = 3) +
	  scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E7298A")) +
	  geom_text_repel(aes(label = gene, color=isHiG), hjust = -0.1, vjust = 0) +
	  theme(
		legend.position = c(0, 0),
		legend.justification = c(1, 1),
		axis.text.x = element_blank(),  # Remove x-axis labels
		axis.ticks.x = element_blank()
	  ) +
	  # stat_compare_means(  doesn't work becasue the x-asix os genes rather than isHiG grouping
		# aes(group = isHiG),method = "wilcox.test",label = "p.signif",label.x = 1.5,size = 4) +
	 annotate("text", x =  (length(unique(CTS_PPI$gene)) + 1) / 2, 
			y = max(CTS_PPI$annd, na.rm=T) * 0.8, 
			label = paste0('wilcox p = ', signif(pval, 2), ' F,T: ', table(CTS_PPI$isHiG) %>% toString()), size = 4)  
	  # ggtitle(paste0('CTS genes of ' ,id2, ' Wilcox test ', pval_label))
	  
	  
	  # Compute exact p-value
	pval <- wilcox.test(annd ~ isCTS, data = HiG_PPI)$p.value
	
	# Reorder gene factor levels by annd (high to low)
	HiG_PPI <- HiG_PPI %>%
		mutate(gene = factor(gene, levels = gene[order(-annd)]))

	plot_HiG_annd[[id]] = ggplot(HiG_PPI, aes(x = gene, y = annd, fill = isCTS)) +
	  geom_boxplot(aes(group = isCTS), width = 0.4, alpha = 0.3, outlier.shape = NA, color = NA) +  # Add boxplot first
	  geom_point(aes(color = isCTS), size = 3) +
	  scale_color_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#E6AB02")) +
	  geom_text_repel(aes(label = gene, color=isCTS), hjust = -0.1, vjust = 0) +
	  theme(
		legend.position = c(0, 0),
		legend.justification = c(1, 1),
		axis.text.x = element_blank(),  # Remove x-axis labels
		axis.ticks.x = element_blank()
	  ) +
	   annotate("text", x = (length(unique(HiG_PPI$gene)) + 1) / 2, 
			y = max(HiG_PPI$annd, na.rm=T) * 0.8, 
			label = paste0('wilcox p = ', signif(pval, 2), ' F,T: ', table(HiG_PPI$isCTS) %>% toString()), size = 4)  
	   #ggtitle(paste0('CTS genes of ' ,id2, ' Wilcox test ', pval_label))
  
}


pdf(file='gene_ranked_by_importance_dotBoxPlot.pdf', height=10.5)
x = ggarrange(plot_CTS_pr[[1]],plot_CTS_pr[[2]],
				plot_CTS_pr[[3]], plot_CTS_pr[[4]], 
          labels = c(paste0("CTS_", CTS[1], " network"), paste0("CTS_", CTS[2], " network" ),
					paste0("CTS_", CTS[3], " network"), paste0("CTS_", CTS[4], " network" )),
          ncol = 2, nrow = 3)  		  
print(x)
	  
x = ggarrange(plot_HiG_annd[[1]],plot_HiG_annd[[2]],
				plot_HiG_annd[[3]], plot_HiG_annd[[4]], 
          labels = c(paste0("HiG_", CTS[1], " network"), paste0("HiG_", CTS[2], " network" ),
					paste0("HiG_", CTS[3], " network"), paste0("HiG_", CTS[4], " network" )),
          ncol = 2, nrow = 3)  		  
print(x)	  

x = ggarrange(plot_HiG_annd[[1]],plot_HiG_annd[[2]],
				plot_HiG_annd[[3]], plot_HiG_annd[[4]], 
          labels = c(paste0("HiG_", CTS[1], " network"), paste0("HiG_", CTS[2], " network" ),
					paste0("HiG_", CTS[3], " network"), paste0("HiG_", CTS[4], " network" )),
          ncol = 2, nrow = 3)  
		  
print(x)	  

for(id in CTS){
	id2 = ifelse(grepl('.', id, fixed=TRUE),  unlist(strsplit(id, split='.',fixed=T))[1],  id)
	x = ggarrange(plot_CTS_annd[[id]], plot_HiG_annd[[id]],  
				plot_CTS_bc[[id]], plot_HiG_bc[[id]], 
				plot_CTS_pr[[id]], plot_HiG_pr[[id]], 
          labels = c(paste0("CTS_", id, " network"), paste0("HiG_", id2, " network" ),
					paste0("CTS_", id, " network"), paste0("HiG_", id2, " network" ),
					paste0("CTS_", id, " network"), paste0("HiG_", id2, " network" )),
          ncol = 2, nrow = 3)  
		  
	print(x)	  
    Sys.sleep(2)
}

dev.off()	
	
	
	
	
	
	
	
	
	
