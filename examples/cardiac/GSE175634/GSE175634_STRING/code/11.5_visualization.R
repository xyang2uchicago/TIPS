
## compare the cardiac CTS from hESC:C10, E2018 2018 cardiac.a, and E8.25 2019 C8 (not the in vivo and EB C4)
## compare HEP CTS from E8.25 2019 C13, E8.25 2018 endo.b, EB C11,
# and E8.25 2019 C17 (or C8)
require(dplyr)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)
library(scales)  # for color gradient
library(patchwork)  # to arrange plots
library(gridExtra)  # to arrange plots

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

###################
## assess CHD scores across three PPI_cat
#################	
CHD = readRDS( file='D:/projects/DS/result/CHD/CHD_Cilia_Genelist.rds')
names(CHD)
CHD = CHD$Griffin2023_PCGC_AllCurated

############################################
#### plot network for HiGCTS, each page is one HiG&CTS in 4 PPI threshiolds
# refer to 11.1-weighted_CTS_cardiac_network.R & 11.2.0_update_network_weights.R to  see the the graph is built
############################################

g = graph_list[["HiGCTS_CP.1"]]
range(E(g)$weight)
# [1] 0.0006951995 0.1414810923

p_listoflist = list()
for(int in c("HiGCTS_CP.1", "HiGCTS_CP", "HiGCTS_muscle", "HiGCTS_endoderm")){ 
#int = "HiGCTS_CP.1", 
	g = graph_list[[int]]
	
	p_listoflist[[int]]  = plot_weighted_PPIN(g, layout = "fr", 
		CHD = CHD, node_size_title = "|Wilcox score|", keep_isolated =TRUE)

}

names(p_listoflist)
#[1] "HiGCTS_CP.1"     "HiGCTS_CP"       "HiGCTS_muscle"   "HiGCTS_endoderm"

pdf(file='network_view_PPI_GSE175634_iPSC_CM.pdf', width=12)
 print( grid.arrange(grobs = p_listoflist[1:2] , ncol=2)) 
 print( grid.arrange(grobs = p_listoflist[3:4], ncol=2)) 
 
dev.off()


#########################################################
## OPTION: plot for Hig PPIN, too dense to show
######################################################### 
 
# score_thresholds = c('score200', 'score400', 'score700', 'weight_shrink')
# p_listoflist = list()

#for(score_threshold in score_thresholds){
score_threshold = 'weight_shrink'
input_path = paste0('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_', score_threshold)

# refer to 11.1_CTS_cardiac_network_degreeDistribution.R 
if(grepl('score',  score_threshold)) {
		graph_list <- readRDS( file= paste0(input_path,'/GSE175634_STRING_graph_perState_notsimplified.rds') )
		graph_list <- lapply(graph_list, simplify) #!!!!!!!!!!!!!!!!!!!
	} else {
		graph_list <- readRDS( file= paste0(input_path, '/GSE175634_STRING_graph_perState_simplified_',s,'weighted.rds') ) 
		}	
	 
p_list_HiG = list()
for(int in grep('^HiG_', names(graph_list), value=TRUE)){
		g = graph_list[[int]]
		
		p_list_HiG[[int]]  = plot_weighted_PPIN(g, layout = "fr", 
		CHD = CHD, node_size_title = "|Wilcox score|", keep_isolated =FALSE)
}
#	p_listoflist[[paste0('PPI',score_threshold)]] = p_list_HiG 
#}

names(p_list_HiG )
# [1] "HiG_0"        "HiG_1"        "HiG_2"        "HiG_CP"      
# [5] "HiG_4"        "HiG_5"        "HiG_6"        "HiG_7"       
# [9] "HiG_muscle"   "HiG_9"        "HiG_10"       "HiG_endoderm"
# [13] "HiG_12"       

CT_id = c(4, 9, 12)
noncardiac_id = c(1,3,8,11)
other_cardiac_id = setdiff(1:13, c(CT_id, noncardiac_id))  # 2  5  6  7 10 13

 
pdf(file= paste0('network_HiG_view_',score_threshold,'_GSE175634_iPSC_CM.pdf'), width=12, height=8)
print( grid.arrange(grobs = p_list_HiG[CT_id], ncol=4)) 
print( grid.arrange(grobs = p_list_HiG[other_cardiac_id], ncol=4)) 
print( grid.arrange(grobs = p_list_HiG[noncardiac_id], ncol=4)) 
 
dev.off()
 	
