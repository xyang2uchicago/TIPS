# this script is used to compare the DEGs of ISL1KO with the DEGs of ISL1KO-connected genes in the HiGCTS_CP.1 module
# ISL1KO DEGs were extracted from D:/projects/DS/data/GSE195476_ISL1/doc/1-s2.0-S2213671123003739-mmc2.xlsx
# we calculated FDR based on the published marker genes' pvalues, verified all markers genes has |log2FC|>0.25, corresponding to FC>1.3
# we also extracted consistent ISL1-binding peask at CP and MNPs, ampping peaks to to nearest proximal genes as ISL1-targets.
# code D:\projects\DS\source\GSE195476_ISL11\ISL1_GS.R
# summary at D:\projects\DS\doc\GSE195476_ISL11\Maven2023_ISl1_data_processed.xlsx
#
# Similarily, we also extracted ISL1KO DEGs from Isl1 targets in D:/projects/DS/data/GSE80383_Isl1/GSE80383_Isl1/doc/Gao2019_Isl1_STable3_ChIPseq_embryo.xlsx
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




########################################################
## GSE80383_Isl1/Isl1_GS.R generated genesets, mapped to  human symbols already
########################################################
load(here::here("examples", "cardiac", "data", "GSE80383_Isl1", "ISL1_Mm2Hg.GS_uniqueSymbol_FC1.2_padj0.05.RData"))
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
ISL1_neighbors <- igraph::neighbors(g, 'ISL1', mode = "all") |> names() %>% unique
ISL1_neighbors 
# [1] "FGF10"  "MEIS1"  "NTRK2"  "HAS2"   "HAND2"  "BMP5"   "CITED2" "BMPER"

# Marv2023 resutls: 
tmp = c('ID1', 'HAPLN1', 'RHOBTB3', 'FAM89A', 'HAS2', 'HAND2', 'PRTG', 'DUSP6', 'PTPN13', 'WLS', 'LAMA1', 'NKX3-1', 'BMPER', 'BMP5', 'LRRTM1')
intersect(tmp, ISL1_neighbors)
 # "HAS2"  "HAND2" "BMPER" "BMP5" 
tmp = c('ID1', 'HAPLN1', 'RHOBTB3', 'FAM89A', 'HAS2', 'HAND2', 'PRTG', 'DUSP6', 'PTPN13', 'WLS', 'LAMA1')
intersect(tmp, ISL1_neighbors)
# "HAS2"  "HAND2"
tmp = c('MEIS1', 'FAM89A', 'PRTG', 'PTPN13', 'BMPER', 'BMP5')
intersect(tmp, ISL1_neighbors)
# "MEIS1" "BMPER" "BMP5"

# Gao 2019 resutls:
tmp = c('DUSP6', 'MEIS1', 'ISL1', 'HAS2', 'LAMA1', 'ZNF608', 'HAND2', 'NKX3-1')
intersect(tmp, ISL1_neighbors) # "MEIS1" "HAS2"  "HAND2"
tmp = c('CITED2', 'DUSP6', 'MEIS1', 'IFI16', 'ISL1', 'FGF10', 'NTRK2', 'HAS2', 'LAMA1', 'ZNF608', 'ID1', 'WLS', 'PTPN13', 'LRRTM1', 'FAM89A', 'SLC7A2', 'HAND2', 'BMPER', 'PRTG', 'HAPLN1', 'SLC40A1', 'NKX3-1', 'CLGN', 'BMP5')
intersect(tmp, ISL1_neighbors) # "MEIS1" "HAS2"  "HAND2"
[1] "CITED2" "MEIS1"  "FGF10"  "NTRK2"  "HAS2"   "HAND2"  "BMPER"  "BMP5"  

	
graph_name = 'HiGCTS_CP.1'
    ### step 2: enrichment test for  up/dn seperately
    ## ordered by Developmental timing → early (iPSC/E8.25–E8.75) → mid (E10.5) → late / in vitro (iPSC-derived CPC) → aggregate
    ## also by mechanistic category (binding → expression → chromatin → intersection)

gene_sets = ISL1[c(    # In vitro (iPSC-derived CPCs)
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
    all(lengths(gene_sets)>0) #TRUE

    res = PPIN_geneset_enrichment(g, gene_sets, key_gene='ISL1', graph_name=graph_name, GS_database='Gao2019') 
    #grid.arrange(res$p_module, res$p_network, res$p_dot, ncol = 3)
    grid.arrange(res$p_module, res$p_network, res$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE80383_Isl1/ISL1_UpDn_Enrichment_",graph_name,".pdf"), width = 10, height = 7)   
   
    write.table(res$enrichment_df, file= paste0("GSE80383_Isl1/Isl1_UpDn_Enrichment_",graph_name,".tsv") , sep='\t', row.names=F) 


## DEG = up | dn
 gene_sets_2 = ISL1[c(    # In vitro (iPSC-derived CPCs)
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
    all(lengths(gene_sets_2)>0) #TRUE

	res2 = PPIN_geneset_enrichment(g, gene_sets_2, key_gene='ISL1', graph_name=graph_name, GS_database='Gao2019', sig_p_adjust=FALSE) 
    
    grid.arrange(res2$p_module, res2$p_network, res2$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE80383_Isl1/Isl1_Enrichment_",graph_name,".pdf"), width = 10, height = 7)   

    write.table(res2$enrichment_df, file= paste0("GSE80383_Isl1/Isl1_Enrichment_",graph_name,".tsv") , sep='\t', row.names=F) 

## Optional: for each graph, test for the toppes degree gene as a control
# by setting key_gene=NULL
###############################################
for(i in 1:length(graph_list)){
    g = graph_list[[i]]
    graph_name = names(graph_list)[i]
 
    res = PPIN_geneset_enrichment(g, gene_sets, key_gene=NULL, graph_name=graph_name, GS_database='Gao2019', sig_p_adjust=FALSE)  
    grid.arrange(res$p_module, res$p_network, res$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE80383_Isl1/UpDn_Enrichment_",graph_name,"_",res$key_gene,".pdf"), width = 10, height = 7)   

    res2 = PPIN_geneset_enrichment(g, gene_sets_2, key_gene=NULL, graph_name=graph_name, GS_database='Gao2019', sig_p_adjust=FALSE)  
    grid.arrange(res2$p_module, res2$p_network, res2$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE80383_Isl1/DEG_Enrichment_",graph_name,"_",res2$key_gene,".pdf"), width = 10, height = 7)   

}

########################################################
## GSE195476_ISL1/ISL1_GS.R generated genesets, iPSC thus human symbols already
########################################################
g <- graph_list[["HiGCTS_CP.1"]]
 # direct neighboring genes of ISL1
 igraph::neighbors(g, 'ISL1', mode = "all") |> names()
# [1] "FGF10"  "MEIS1"  "NTRK2"  "HAS2"   "HAND2"  "BMP5"   "CITED2" "BMPER" 

graph_name = 'HiGCTS_CP.1'

ISL1_set = readRDS( file=here::here("examples", "cardiac", "data", "GSE195476_ISL1", "ISL1_set.rds"))
ISL1_set[['ISL1.targets_CP']] = intersect(ISL1_set[['ISL1_WT_d6CP']], 
									unique(c(ISL1_set[['ISL_DEG_d8CP_E']], ISL1_set[['ISL_DEG_d8CP_L']], ISL1_set[['ISL_DEG_d8CP_T']]))
									)    # catdiac specific union 207     

lengths(ISL1_set)
    # ISL1_NKO_d6CP    ISL1_WT_d18MNP      ISL1_WT_d6CP    NKX25_NKO_d6CP 
             # 1987              2917              1724               266 
    # NKX25_WT_d6CP         ISL_up_E         ISL_up_T         ISL_up_L 
            # 12141               121               155               946 
        # ISL_dn_E         ISL_dn_T         ISL_dn_L       ISL_up_MNP 
              # 449               242               586               144 
      # ISL_dn_MNP    ISL_DEG_d8CP_E    ISL_DEG_d8CP_L    ISL_DEG_d8CP_T 
              # 422               570              1532               397 
   # ISL_DEG_d18MNP ISL1_targets_CP_E ISL1_targets_CP_T ISL1_targets_CP_L 
              # 566               111               110               208  		
 # ISL1_targets_MNP  ISL1.targets_CP
              # 119  			270

    gene_sets = ISL1_set[c(    # In vitro (iPSC day6-8)
    # early/transition/later CP (d8) ISL1 activity
    "ISL_up_E"    ,       "ISL_dn_E"       ,           "ISL1_targets_CP_E"    ,   
    "ISL_up_T"    ,       "ISL_dn_T"       ,           "ISL1_targets_CP_T"    , 
	"ISL_up_L"    ,       "ISL_dn_L"       ,           "ISL1_targets_CP_L"    ,  	
	"ISL1_WT_d6CP" ,
    # MNP (iPSC day 18)
    "ISL_up_MNP"  ,   "ISL_dn_MNP" ,    "ISL1_targets_MNP", 
	"ISL1_WT_d18MNP",
    # Cross-stage integrated
    "ISL1.targets_CP"        # cross-CP stage union
 )]
    all(lengths(gene_sets)>0) #TRUE
	
    res = PPIN_geneset_enrichment(g, gene_sets, key_gene='ISL1', graph_name=graph_name, GS_database='Maven2023', sig_p_adjust=FALSE) 
    #grid.arrange(res$p_module, res$p_network, res$p_dot, ncol = 3)
    grid.arrange(res$p_module, res$p_network, res$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE195476_ISL1/ISL1_UpDn_Enrichment_",graph_name,".pdf"), width = 10, height = 7)   
   
    write.table(res$enrichment_df, file= paste0("GSE195476_ISL1/ISL1_UpDn_Enrichment_",graph_name,".tsv") , sep='\t', row.names=F) 


## DEG = up | dn
 gene_sets_2 = ISL1_set[c(    # In vitro (iPSC day6-8)
    # early/transition/later CP (d8) ISL1 activity
    "ISL_DEG_d8CP_E"       ,           "ISL1_targets_CP_E"    ,   
    "ISL_DEG_d8CP_T"       ,           "ISL1_targets_CP_T"    , 
	"ISL_DEG_d8CP_L"       ,           "ISL1_targets_CP_L"    ,  	
	"ISL1_WT_d6CP" ,
    # MNP (iPSC day 18)
    "ISL_DEG_d18MNP" ,    "ISL1_targets_MNP", 
	"ISL1_WT_d18MNP",
    # Cross-stage integrated
    "ISL1.targets_CP"        # cross-CP stage union
	)]
	all(lengths(gene_sets_2)>0) #TRUE
	
    res2 = PPIN_geneset_enrichment(g, gene_sets_2, key_gene='ISL1', graph_name=graph_name, GS_database='Maven2023', sig_p_adjust=FALSE) 
    
    grid.arrange(res2$p_module, res2$p_network, res2$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE195476_ISL1/ISL1_Enrichment_",graph_name,".pdf"), width = 10, height = 7)   

    write.table(res2$enrichment_df, file= paste0("GSE195476_ISL1/ISL1_Enrichment_",graph_name,".tsv") , sep='\t', row.names=F) 


## Optional: for each graph, test for the toppes degree gene as a control
# by setting key_gene=NULL
###############################################
for(i in 1:length(graph_list)){
    g = graph_list[[i]]
    graph_name = names(graph_list)[i]
 
    res = PPIN_geneset_enrichment(g, gene_sets, key_gene=NULL, graph_name=graph_name, GS_database='Maven2023', sig_p_adjust=FALSE)  
    grid.arrange(res$p_module, res$p_network, res$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE195476_ISL1/UpDn_Enrichment_",graph_name,"_",res$key_gene,".pdf"), width = 10, height = 7)   

    res2 = PPIN_geneset_enrichment(g, gene_sets_2, key_gene=NULL, graph_name=graph_name, GS_database='Maven2023', sig_p_adjust=FALSE)  
    grid.arrange(res2$p_module, res2$p_network, res2$p_dot, layout_matrix = rbind(c(1, 2), c(3, 3)))
    dev.copy2pdf(file = paste0("GSE195476_ISL1/DEG_Enrichment_",graph_name,"_",res2$key_gene,".pdf"), width = 10, height = 7)   

}





# ################################################## 
# ################ ISL1 chip-seq motifs peaks ###########
# ################################################## 
# tmp = read.csv('D:/projects/DS/data/GSE195476_ISL1/Maven_ISL1_2022-main/analysis/bed/mapped_motifs/all_instances_curated_0based.bed',sep='\t', header=F)
# dim(tmp)  #[1] 53353     6
# head(tmp)
# #     V1        V2        V3                V4 V5 V6
# # 1 chr1 100629813 100629851 LHX-ISL1-28_1_558  0  -
# # 2 chr1 100643263 100643280 LHX-ISL1-10_2_559  0  -
# # 3 chr1 100643287 100643325 LHX-ISL1-28_3_559  0  -
# # 4 chr1 100643336 100643346        EBF1_4_559  0  +
# # 5 chr1 100643398 100643436 LHX-ISL1-28_5_559  0  -
# # 6 chr1 101471599 101471637 LHX-ISL1-28_6_560  0  +
# tmp$group = lapply(tmp$V4, function(x) unlist(strsplit(x, split="_"))[1]) %>% unlist
# table(tmp$group)
# #        EBF1        GATA        ISL1         LHX LHX-ISL1-10 LHX-ISL1-28  LHX-ISL1-9      NeuroD      NKX2.5 
# #        2338        5185        2917        4832        2075       16759        2194        4010        4516 
# #  NKX2.5-alt     Onecut2 
# #        1827        6700 