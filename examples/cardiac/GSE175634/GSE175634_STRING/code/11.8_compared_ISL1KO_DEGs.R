# this script is used to compare the DEGs of ISL1KO with the DEGs of ISL1KO-connected genes in the HiGCTS_CP.1 module
# ISL1KO DEGs were extracted from D:/projects/DS/data/GSE195476_ISL1/doc/1-s2.0-S2213671123003739-mmc2.xlsx
# we calculated FRD based on the published marker genes' pvalues, verified all markers genes has |log2FC|>0.25, corresponding to FC>1.3
# code here
# we also extracted ISL1KO DEGs from Isl1 targets in D:/projects/DS/data/GSE80383_Isl1/GSE80383_Isl1/doc/Gao2019_Isl1_STable3_ChIPseq_embryo.xlsx
# code at D:\projects\DS\source\GSE80383_Isl1\Isl1_GS.R
# summary at D:\projects\DS\doc\GSE80383_Isl1\Gao2019_ISl1_data_processed.xlsx

#require(dplyr)
library(ggplot2)
library(igraph)
#library(tidygraph)
#library(ggraph)
library(scales)  # for color gradient
library(patchwork)  # to arrange plots
library(gridExtra)  # to arrange plots
#library(pracma)
#library(data.table)
library(ggpubr)
library(readxl)
library(gplots)
library(UpSetR)

 
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
## load published CHD and DEGs
#################	
CHD = readRDS( file='D:/projects/DS/result/CHD/CHD_Cilia_Genelist.rds')
names(CHD)
#CHD = CHD$Griffin2023_PCGC_AllCurated


ISL1_set = readRDS( file='D:/projects/DS/result/GSE195476_ISL1/ISL1_set.rds')
lengths(ISL1_set)
    # ISL1_NKO_d6CP    ISL1_WT_d18MNP      ISL1_WT_d6CP    NKX25_NKO_d6CP 
             # 1987              2917              1724               266 
    # NKX25_WT_d6CP         ISL_act_E         ISL_act_T         ISL_act_L 
            # 12141               121               155               946 
        # ISL_rep_E         ISL_rep_T         ISL_rep_L       ISL_act_MNP 
              # 449               242               586               144 
      # ISL_rep_MNP    ISL_DEG_d8CP_E    ISL_DEG_d8CP_L    ISL_DEG_d8CP_T 
              # 422               570              1532               397 
   # ISL_DEG_d18MNP ISL1_targets_CP_E ISL1_targets_CP_T ISL1_targets_CP_L 
              # 566               111               110               208  		
 # ISL1_targets_MNP 
              # 119 

# Prepare gene set list: 4 published + network nodes
gene_sets <- ISL1_set[c('ISL_act_E', 'ISL_act_T', 'ISL_act_L', 'ISL_rep_E', 'ISL_rep_T', 'ISL_rep_L', 'ISL_act_MNP', 'ISL_rep_MNP')]
gene_sets_2 <- ISL1_set[c('ISL_DEG_E', 'ISL_DEG_L', 'ISL_DEG_T', 'ISL_DEG_MNP')]


ISL1_set[['CHD_validated']] = CHD$Griffin2023_PCGC
ISL1_set[['CHD_curated']] = CHD$Griffin2023_PCGC_AllCurated
 
set_order <- c("CHD_validated", "CHD_curated",  "ISL1_targets_MNP" ,   "ISL1_targets_CP_L",  "ISL1_targets_CP_T", "ISL1_targets_CP_E")
UpSetR::upset(
  fromList(ISL1_set[set_order]),
  sets = set_order,        # this controls vertical order of sets
  keep.order = TRUE,
  nsets = 6,
  nintersects = NA,
  sets.bar.color = c("red", "pink","orange", "purple", "blue", "green" ),
  order.by = "freq",
  main.bar.color = "#636363",
  number.angles = 45,
  text.scale = c(1.3, 1.3, 1, 0.8, 0.8, 0.8),
  mainbar.y.label = "Intersection size",
  sets.x.label = "Total DEGs per stage"
)
intersect(ISL1_set[["CHD_curated"]], ISL1_set[['ISL1_targets_CP_T']])
# [1] "FOXP1" "GPC3"  "JAG1"  "MEF2C" "SMAD6"
intersect(ISL1_set[["CHD_curated"]], ISL1_set[['ISL1_targets_CP_E']])
# [1] "FBN2"  "FOXP1" "GPC3"  "JAG1"  "MEF2C" "MEIS2" "MSX2"  "PBX1" 
intersect(ISL1_set[["CHD_curated"]], ISL1_set[['ISL1_targets_CP_L']])
# [1] "ARHGAP31" "BMP10"    "CDH2"     "FBN2"     "GATA4"    "GPC3"    
# [7] "MEF2C"    "PRKAG2"   "SMAD2"   
 

lengths(gene_sets)
#     ISL_act_E     ISL_act_T     ISL_act_L     ISL_rep_E     ISL_rep_T     ISL_rep_L   ISL_act_MNP   ISL_rep_MNP 
#           121           155           946           449           242           586           144           422 
# CHD_validated   CHD_curated 
#            99           295 
lengths(gene_sets_2)
# ISL_DEG_E     ISL_DEG_L     ISL_DEG_T   ISL_DEG_MNP CHD_validated   CHD_curated 
#       570          1532           397           566            99           295 

################################################
## test if ISL1-connecgted genes are significantly enriched in the HiGCTS_CP.1 module
###############################################

## step 1: prepare the network and the gene sets
 g <- graph_list[["HiGCTS_CP.1"]]
 # direct neighboring genes of ISL1
 igraph::neighbors(g, 'ISL1', mode = "all") |> names()
# [1] "FGF10"  "MEIS1"  "NTRK2"  "HAS2"   "HAND2"  "BMP5"   "CITED2" "BMPER" 

 graph_name = 'HiGCTS_CP.1'
    ### step 2: enrichment test for each gene set
    res = PPIN_geneset_enrichment(g, gene_sets, key_gene='ISL1', graph_name=graph_name, GS_database='Gifford2023_iPSC') 
    #grid.arrange(res$p_module, res$p_network, res$p_dot, ncol = 3)
    grid.arrange(res$p_module, res$p_network, res$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE175634_ISL1/ISL1_UpDn_Enrichment_",graph_name,".pdf"), width = 7, height = 5)   

    write.table(res$enrichment_df, file= paste0("GSE175634_ISL1/ISL1_UpDn_Enrichment_",graph_name,".tsv"), sep='\t', row.names=F) 


    names(gene_sets_2)
    # [1] "ISL_DEG_E"     "ISL_DEG_L"     "ISL_DEG_T"     "CHD_validated" "CHD_curated"  
    res2 = PPIN_geneset_enrichment(g, gene_sets_2, key_gene='ISL1', graph_name=graph_name, GS_database='Gifford2023_iPSC') 
    #grid.arrange(res2$p_module, res2$p_component, res2$p_network, ncol = 3)
    grid.arrange(res2$p_module, res2$p_network, res2$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE175634_ISL1/ISL1_DEG_Enrichment_",graph_name,".pdf"), width = 7, height = 5)   
   
    write.table(res2$enrichment_df, file= paste0("GSE175634_ISL1/ISL1_DEG_Enrichment_",graph_name,".tsv") , sep='\t', row.names=F) 
 


################################################
## for each graph, test for the toppes degree gene as a control
# by setting key_gene=NULL
###############################################
for(i in 1:length(graph_list)){
    g = graph_list[[i]]
    graph_name = names(graph_list)[i]
 
    res = PPIN_geneset_enrichment(g, gene_sets, key_gene=NULL, graph_name=graph_name)  
    grid.arrange(res$p_module, res$p_component, res$p_network,  ncol = 3)
    dev.copy2pdf(file = paste0("GSE175634_ISL1/UpDn_Enrichment_",graph_name,"_",res$key_gene,".pdf"), width = 15, height = 5)   

    res2 = PPIN_geneset_enrichment(g, gene_sets_2, key_gene=NULL, graph_name=graph_name)  
    grid.arrange(res2$p_module, res2$p_component, res2$p_network, ncol = 3)
    dev.copy2pdf(file = paste0("GSE175634_ISL1/DEG_Enrichment_",graph_name,"_",res2$key_gene,".pdf"), width = 15, height = 5)   

}

################################################## 
################ ISL1 chip-seq peaks ###########
################################################## 
tmp = read.csv('D:/projects/DS/data/GSE195476_ISL1/Maven_ISL1_2022-main/analysis/bed/mapped_motifs/all_instances_curated_0based.bed',sep='\t', header=F)
dim(tmp)  #[1] 53353     6
head(tmp)
#     V1        V2        V3                V4 V5 V6
# 1 chr1 100629813 100629851 LHX-ISL1-28_1_558  0  -
# 2 chr1 100643263 100643280 LHX-ISL1-10_2_559  0  -
# 3 chr1 100643287 100643325 LHX-ISL1-28_3_559  0  -
# 4 chr1 100643336 100643346        EBF1_4_559  0  +
# 5 chr1 100643398 100643436 LHX-ISL1-28_5_559  0  -
# 6 chr1 101471599 101471637 LHX-ISL1-28_6_560  0  +
tmp$group = lapply(tmp$V4, function(x) unlist(strsplit(x, split="_"))[1]) %>% unlist
table(tmp$group)
#        EBF1        GATA        ISL1         LHX LHX-ISL1-10 LHX-ISL1-28  LHX-ISL1-9      NeuroD      NKX2.5 
#        2338        5185        2917        4832        2075       16759        2194        4010        4516 
#  NKX2.5-alt     Onecut2 
#        1827        6700 

########################################################
## GSE80383_Isl1/Isl1_GS.R generated genesets, mapped to  human symbols already
########################################################
load('D:/projects/DS/result/GSE80383_Isl1/ISL1_Mm2Hg.GS_uniqueSymbol_FC1.2_padj0.05.RData')
names(ISL1)
#  [1] "Isl1.embryo.bound"         "Isl1.iPSC.bound"           "Isl1KO.E8.75.up"           "Isl1KO.E8.75.dn"          
#  [5] "Isl1.E8.75.DEG"            "Isl1KO.E10.5RV.OFT.up"     "Isl1KO.E10.5RV.OFT.dn"     "Isl1.E10.5RV.OFT.DEG"     
#  [9] "Isl1.target"               "Isl1.target.E8.75"         "Isl1.target.E10.5"         "Isl1KO.CP.up"             
# [13] "Isl1KO.CP.dn"              "Isl1.CP.DEG"               "Isl1.target.CP"            "iPSC.open"                
# [17] "E8.75.open"                "Isl1.E8.75DE.pooled.bound" "Isl1.E10.5DE.pooled.bound"

g <- graph_list[["HiGCTS_CP.1"]]
 # direct neighboring genes of ISL1
 igraph::neighbors(g, 'ISL1', mode = "all") |> names()
# [1] "FGF10"  "MEIS1"  "NTRK2"  "HAS2"   "HAND2"  "BMP5"   "CITED2" "BMPER" 

 graph_name = 'HiGCTS_CP.1'
    ### step 2: enrichment test for  up/dn seperately
    ## ordered by Developmental timing → early (iPSC/E8.25–E8.75) → mid (E10.5) → late / in vitro (iPSC-derived CPC) → aggregate
    ## also by mechanistic category (binding → expression → chromatin → intersection)

    ISL1_gs = ISL1[c(    # In vitro (iPSC-derived CPCs)
    "Isl1.iPSC.bound"  ,   "iPSC.open" ,  
    # early, embryonic ISL1 activity
    "Isl1KO.CP.up"    ,       "Isl1KO.CP.dn"       ,           "Isl1.target.CP"    ,  #iSPC-derived CPC Gao et al., eLife 2019: strongly with E8.25–E9 embryonic peaks.
    # Gao said: '71% of genes bound by Isl1 in ESC-derived CPCs were also bound in E8.25-E9 embryos (Fig 3c)'
    "Isl1.embryo.bound" ,              
    "Isl1KO.E8.75.up"    ,       "Isl1KO.E8.75.dn"        ,  
    "E8.75.open" ,
    "Isl1.E8.75DE.pooled.bound" , # intersect(union(Isl1[['Isl1.embryo.bound']] ,Isl1[['Isl1.iPSC.bound']]), Isl1[['Isl1.E8.75.DEG']]  )
    # Mid-development (E10.5 embryo and integrative targets)
    "Isl1KO.E10.5RV.OFT.up"  ,   "Isl1KO.E10.5RV.OFT.dn" ,    
    "Isl1.target.E8.75"    ,     "Isl1.target.E10.5"     ,
    "Isl1.E10.5DE.pooled.bound", # intersect(union(Isl1[['Isl1.embryo.bound']] ,Isl1[['Isl1.iPSC.bound']]), Isl1[['Isl1KO.E10.5RV.OFT.DEG']] )
    # Cross-stage integrated
    "Isl1.target"        # intersect(union(Isl1[['Isl1.embryo.bound']] ,Isl1[['Isl1.iPSC.bound']]), union(Isl1[['Isl1.E8.75.DEG']] ,Isl1[['Isl1KO.E10.5RV.OFT.DEG']]) )        
 )]
    
    res = PPIN_geneset_enrichment(g, ISL1_gs, key_gene='ISL1', graph_name=graph_name, GS_database='Gao2019') 
    #grid.arrange(res$p_module, res$p_network, res$p_dot, ncol = 3)
    grid.arrange(res$p_module, res$p_network, res$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE80383_Isl1/ISL1_UpDn_Enrichment_",graph_name,".pdf"), width = 10, height = 7)   
   
    write.table(res2$enrichment_df, file= paste0("GSE80383_Isl1/Isl1_UpDn_Enrichment_",graph_name,".tsv") , sep='\t', row.names=F) 


## DEG = up | dn
 ISL1_gs = ISL1[c(    # In vitro (iPSC-derived CPCs)
    "Isl1.iPSC.bound"  ,   "iPSC.open" ,  
    # early, embryonic ISL1 activity
    "Isl1.CP.DEG"       ,           "Isl1.target.CP"    ,  #iSPC-derived CPC Gao et al., eLife 2019: strongly with E8.25–E9 embryonic peaks.
    # Gao said: '71% of genes bound by Isl1 in ESC-derived CPCs were also bound in E8.25-E9 embryos (Fig 3c)'
    "Isl1.embryo.bound" ,              
    "Isl1.E8.75.DEG"        ,  
    "E8.75.open" ,
    "Isl1.E8.75DE.pooled.bound" , # intersect(union(Isl1[['Isl1.embryo.bound']] ,Isl1[['Isl1.iPSC.bound']]), Isl1[['Isl1.E8.75.DEG']]  )
    # Mid-development (E10.5 embryo and integrative targets)
    "Isl1.E10.5RV.OFT.DEG" ,    
    "Isl1.target.E8.75"    ,     "Isl1.target.E10.5"     ,
    "Isl1.E10.5DE.pooled.bound", # intersect(union(Isl1[['Isl1.embryo.bound']] ,Isl1[['Isl1.iPSC.bound']]), Isl1[['Isl1KO.E10.5RV.OFT.DEG']] )
    # Cross-stage integrated
    "Isl1.target"        # intersect(union(Isl1[['Isl1.embryo.bound']] ,Isl1[['Isl1.iPSC.bound']]), union(Isl1[['Isl1.E8.75.DEG']] ,Isl1[['Isl1KO.E10.5RV.OFT.DEG']]) )        
 )]
    res2 = PPIN_geneset_enrichment(g, ISL1_gs, key_gene='ISL1', graph_name=graph_name, GS_database='Gao2019') 
    
    grid.arrange(res2$p_module, res2$p_network, res2$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE80383_Isl1/DEG_Enrichment_",graph_name,".pdf"), width = 10, height = 7)   

    write.table(res2$enrichment_df, file= paste0("GSE80383_Isl1/DEG_Enrichment_",graph_name,".tsv") , sep='\t', row.names=F) 
