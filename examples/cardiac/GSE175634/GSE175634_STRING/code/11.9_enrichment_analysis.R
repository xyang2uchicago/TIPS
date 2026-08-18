library(ggplot2)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(seq2pathway.data)
library(seq2pathway)



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

#############################################################
## identify DEGs by both meta.p and mean.FC
########################################################  
## if working on entrez IDs
# library(EnrichmentBrowser)
# packageVersion("EnrichmentBrowser")
# #[1] '2.32.0'
# kegg.gs <- getGenesets(org="hsa", db="kegg")  # containing NCBI Entrez Gene IDs
# length(kegg.gs) # 367

# class(kegg.gs)
# #[1] "list"
# head(kegg.gs, 2)
# # $`hsa00010_Glycolysis_/_Gluconeogenesis`
# #  [1] "10327"  "124"    "125"    "126"    "127"    "128"    "130"    "130589" "131"    "160287" "1737"   "1738"  
# # [13] "2023"   "2026"   "2027"   "217"    "218"    "219"    "2203"   "221"    "222"    "223"    "224"    "226"   
# # [25] "229"    "230"    "2538"   "2597"   "26330"  "2645"   "2821"   "3098"   "3099"   "3101"   "387712" "3939"  
# # [37] "3945"   "3948"   "441531" "501"    "5105"   "5106"   "5160"   "5161"   "5162"   "5211"   "5213"   "5214"  
# # [49] "5223"   "5224"   "5230"   "5232"   "5236"   "5313"   "5315"   "55276"  "55902"  "57818"  "669"    "7167"  
# # [61] "80201"  "83440"  "84532"  "8789"   "92483"  "92579"  "9562"  

# # $`hsa00020_Citrate_cycle_(TCA_cycle)`
# #  [1] "1431"  "1737"  "1738"  "1743"  "2271"  "3417"  "3418"  "3419"  "3420"  "3421"  "4190"  "4191"  "47"    "48"   
# # [15] "4967"  "50"    "5091"  "5105"  "5106"  "5160"  "5161"  "5162"  "55753" "6389"  "6390"  "6391"  "6392"  "8801" 
# # [29] "8802"  "8803" 

# we worked on gene symbols for simplicity and consistency
library(GSA)

h.gs.gene2pathway <- GSA::GSA.read.gmt('F:/projects/Share/MsigDB/v2022.1_Hs/c3.tft.v2022.1.Hs.symbols.gmt') #  hallmark gene sets
names(h.gs.gene2pathway[[1]]) <- h.gs.gene2pathway[[2]]
names(h.gs.gene2pathway[[2]]) <- h.gs.gene2pathway[[2]] # 1127

c2.cp.gs.gene2pathway <- GSA::GSA.read.gmt('F:/projects/Share/MsigDB/v2022.1_Hs/c2.cp.v2022.1.Hs.symbols.gmt') # Canonical pathways
names(c2.cp.gs.gene2pathway[[1]]) <- c2.cp.gs.gene2pathway[[2]]
names(c2.cp.gs.gene2pathway[[2]]) <- c2.cp.gs.gene2pathway[[2]] # 3050

bp.gs.gene2pathway <- GSA::GSA.read.gmt('F:/projects/Share/MsigDB/v2022.1_Hs/c5.go.bp.v2022.1.Hs.symbols.gmt') # GO BP
names(bp.gs.gene2pathway[[1]]) <- bp.gs.gene2pathway[[2]]
names(bp.gs.gene2pathway[[2]]) <- bp.gs.gene2pathway[[2]]  # 7762

tft.gs.gene2pathway <- GSA::GSA.read.gmt('F:/projects/Share/MsigDB/v2022.1_Hs/c3.tft.v2022.1.Hs.symbols.gmt') #  transcription factor targets
names(tft.gs.gene2pathway[[1]]) <- tft.gs.gene2pathway[[2]]
names(tft.gs.gene2pathway[[2]]) <- tft.gs.gene2pathway[[2]] # 1126


################################################
## test if ISL1-connecgted genes are significantly enriched for any knwon GSs
###############################################

library(seq2pathway)
library(seq2pathway.data)
packageVersion("seq2pathway") #'1.28.0'
source('F:/projects/Seq2Path/seq2pathway/FisherTest_GO_BP_MF_CC_05242022.R')

## step 1: prepare the network and the gene sets
# g <- graph_list[["HiGCTS_CP.1"]]
# direct neighboring genes of ISL1
# direct_node = igraph::neighbors(g, 'ISL1', mode = "all") |> names()
# # [1] "FGF10"  "MEIS1"  "NTRK2"  "HAS2"   "HAND2"  "BMP5"   "CITED2" "BMPER" 
# nodes = V(g)$name
# ## step 2: enrichment test for each gene set


# generate the data.frame of background gene list required by FisherTest_GO_BP_MF_CC()
gene_names = read.table('F:/projects/scRNA/results/GSE175634_iPSC_CM/2024_3kHVG/geneNames_raw.tsv') 
dim(gene_names) # 38943     1
head(gene_names, 3)
#           V1
# 1     WASH7P
# 2 AL627309.1
# 3     CICP27

GO_GENCODE_df = data.frame(gene_name=gene_names[,1])
c2.cp_nodes <- h_nodes <- bp_nodes <- tft_nodes <- list()
savepath = 'seq2pathway_res/'

### enrichment test
for(i in 1:length(graph_list)){
    g <- graph_list[["HiGCTS_CP.1"]]
    # direct neighboring genes of ISL1
    direct_node = igraph::neighbors(g, 'ISL1', mode = "all") |> names()
    # [1] "FGF10"  "MEIS1"  "NTRK2"  "HAS2"   "HAND2"  "BMP5"   "CITED2" "BMPER" 
    nodes = V(g)$name
    graph_name = names(graph_list)[i]

    tmp <- FisherTest_GO_BP_MF_CC(gs= nodes, #gene.id[DEG[[j]]] %>% unique,
                                  genome='hg38',
                                  Ontology="newOntology", newOntology= c2.cp.gs.gene2pathway, 
                                  min_Intersect_Count=3,
                                  GO_GENCODE_df = GO_GENCODE_df)
    tmp <- tmp[[1]] 
    c2.cp_nodes[[graph_name]] <- subset(tmp, FDR<0.1 & Fisher_odds >2 & GO_gene_inBackground< 500)
      
    tmp <- FisherTest_GO_BP_MF_CC(gs= nodes, #gene.id[DEG[[j]]] %>% unique,
                                  genome='hg38',
                                  Ontology="newOntology", newOntology= h.gs.gene2pathway, 
                                  min_Intersect_Count=3,
                                  GO_GENCODE_df = GO_GENCODE_df)
    tmp <- tmp[[1]] 
    h_nodes[[graph_name]] <- subset(tmp, FDR<0.1 & Fisher_odds >2 & GO_gene_inBackground< 500)
    
    tmp <- FisherTest_GO_BP_MF_CC(gs= nodes,  #gene.id[DEG[[j]]] %>% unique,
                                  genome='hg38',
                                  Ontology="newOntology", newOntology= bp.gs.gene2pathway,
                                  min_Intersect_Count=3,
                                  GO_GENCODE_df = GO_GENCODE_df)
    tmp <- tmp[[1]]  
    bp_nodes[[graph_name]] <- subset(tmp, FDR<0.1 & Fisher_odds >2 & GO_gene_inBackground< 500)
        
    tmp <- FisherTest_GO_BP_MF_CC(gs= nodes, #gene.id[DEG[[j]]] %>% unique,
                                  genome='hg38',
                                  Ontology="newOntology", newOntology= tft.gs.gene2pathway, 
                                  min_Intersect_Count=3,
                                  GO_GENCODE_df = GO_GENCODE_df)
    tmp <- tmp[[1]] 
    tft_nodes[[graph_name]] <- subset(tmp, FDR<0.1 & Fisher_odds >2 & GO_gene_inBackground< 500)
}   
 
  # for(j in 1:length(DEG))
  names(c2.cp_nodes)  
  # [1] "HiGCTS_CP.1"     "HiG_0"           "HiG_1"           "HiG_2"          
#  [5] "HiG_CP"          "HiG_4"           "HiG_5"           "HiG_6"          
#  [9] "HiG_7"           "HiG_muscle"      "HiG_9"           "HiG_10"         
# [13] "HiG_endoderm"    "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm"
# [17] "HiGCTS_CP"       "CTS_muscle"      "CTS_endoderm"    "CTS_CP"         
# [21] "CTS_CP.1" 
  save(c2.cp_nodes, file= paste0(savepath, 'MsigDB.c2.cp_nodes_FISHER.sig.RData' ))
  save(h_nodes, file= paste0(savepath, 'MsigDB.h_nodes_FISHER.sig.RData' ))
  save(bp_nodes, file= paste0(savepath, 'bp_nodes_FISHER.sig.RData' ))
  save(tft_nodes, file= paste0(savepath, 'MsigDB.C3.tft_nodes_FISHER.sig.RData' ))
  
graph_name = 'HiGCTS_CP.1'
subset(bp_nodes[[graph_name]], FDR<0.01 & Fisher_odds > 5 
        & GO_gene_inBackground< 500 & Intersect_Count>=5)  

subset(tft_nodes[[graph_name]], FDR<0.01 & Fisher_odds > 5
        & GO_gene_inBackground< 500 & Intersect_Count>=4)  
                   GOID          Description Fisher_Pvalue Fisher_odds          FDR Intersect_Count GO_gene_inBackground GO_gene_raw_Count
# 1     GCTNWTTGK_UNKNOWN    GCTNWTTGK_UNKNOWN  1.377581e-06   21.327986 0.0001283312               6                  295               310
# 2               HNF3_Q6  (FOXA family) HNF3_Q6  2.566623e-06   28.068942 0.0001283312               5                  180               192
# 4  RTTTNNNYTGGM_UNKNOWN RTTTNNNYTGGM_UNKNOWN  3.218042e-05   26.086508 0.0008045105               4                  148               155
# 6     AAAYRNCTG_UNKNOWN    AAAYRNCTG_UNKNOWN  7.450782e-05   13.656957 0.0012417970               5                  362               375
# 7              FOXD3_01             FOXD3_01  9.413415e-05   19.639709 0.0013447735               4                  195               202
# 8               EVI1_04              EVI1_04  1.718693e-04   16.722844 0.0021094768               4                  228               241
# 9                 S8_01                S8_01  1.898529e-04   16.282639 0.0021094768               4                  234               247
# 10              LYF1_01              LYF1_01  3.095229e-04   14.277025 0.0025132971               4                  266               271
# 11             NKX25_02             NKX25_02  2.635945e-04   14.908642 0.0025132971               4                  255               267
# 12             HOXA4_Q2             HOXA4_Q2  2.838134e-04   14.614825 0.0025132971               4                  260               272
# 13     WTTGKCTG_UNKNOWN     WTTGKCTG_UNKNOWN  3.267286e-04    9.841441 0.0025132971               5                  497               519
# 15     CATTGTYY_SOX9_B1     CATTGTYY_SOX9_B1  9.451039e-04   10.518727 0.0063006927               4                  358               368
# 17     SMTTTTGT_UNKNOWN     SMTTTTGT_UNKNOWN  1.347075e-03    9.531041 0.0074837522               4                  394               408
# 19   AUTS2_TARGET_GENES   AUTS2_TARGET_GENES  1.587383e-03    9.099702 0.0076864743               4                  412               455
#                           Intersect_gene
# 1  CITED2 LAMA1 DUSP6 FGF10 MEIS1 NKX3-1
# 2            BMP5 HAPLN1 DUSP6 HAND2 ID1
# 4                 HAS2 NTRK2 BMP5 HAPLN1
# 6           CITED2 NTRK2 FGF10 MEIS1 ID1
# 7               CITED2 HAS2 LRRTM1 HAND2
# 8               HAS2 LRRTM1 SLC40A1 ISL1
# 9                 HAS2 BMP5 LRRTM1 MEIS1
# 10              HAPLN1 DUSP6 FGF10 MEIS1
# 11               CITED2 HAS2 NTRK2 FGF10
# 12               NTRK2 BMP5 LRRTM1 DUSP6
# 13           DUSP6 MEIS1 ISL1 HAS2 NTRK2
# 15                LAMA1 BMP5 FGF10 MEIS1
# 17               CITED2 DUSP6 ISL1 NTRK2
# 19              ISL1 DUSP6 CITED2 ZNF608

(df = subset(c2.cp_nodes[[graph_name]], FDR<0.01 )) #& Fisher_odds > 5 & GO_gene_inBackground< 500 & Intersect_Count>=3))  
#                              GOID                     Description Fisher_Pvalue Fisher_odds          FDR Intersect_Count
# 1            WP_HEART_DEVELOPMENT            WP_HEART_DEVELOPMENT  7.161285e-05    44.96409 0.0008593541               3
# 2 WP_NEURAL_CREST_DIFFERENTIATION WP_NEURAL_CREST_DIFFERENTIATION  7.934380e-04    19.12388 0.0047606283               3
# 3     REACTOME_SIGNALING_BY_NTRKS     REACTOME_SIGNALING_BY_NTRKS  1.860799e-03    14.09045 0.0065440215               3
#   GO_gene_inBackground GO_gene_raw_Count   Intersect_gene
# 1                   44                47 HAND2 ISL1 FGF10
# 2                   99               101    ISL1 ID1 PRTG
# 3                  133               134  ID1 DUSP6 NTRK2
GS_enrichment_Dotplot(df, GS_database= 'Msigdb.c2.cp', sig_p_adjust=TRUE, n_gene_characters=30)
dev.copy2pdf(file=paste0('seq2pathway_res/Msigdb.c2.cp_Dotplot_',graph_name,'.pdf'), height=3.5)