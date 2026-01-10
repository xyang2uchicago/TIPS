
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
# source('/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/celltype_specific_weight_v10.R')  # midway3 run
setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9/')

# setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9')
source('E:/Git_Holly/celltype_specific_weight_v10.R')



#########################################################
###### load published PPIN genes ############
#########################################################
PPI_color_palette = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")
PPI_size_platte = c("CTS" = 1, "HiGCTS" = 0.75, "HiG" = 0.25)
signauture_levels = c("CTS_muscle",  "CTS_endoderm",  "CTS_CP",   "CTS_CP.1" ,   "HiGCTS_muscle" , 
	"HiGCTS_endoderm" ,"HiGCTS_CP"  ,     "HiGCTS_CP.1" ,    "HiG_0" ,          "HiG_1",         
	"HiG_2" ,          "HiG_CP"   ,       "HiG_4"   ,        "HiG_5"  ,         "HiG_6" ,         
	"HiG_7" ,          "HiG_muscle" ,     "HiG_9" ,          "HiG_10" ,         "HiG_endoderm" ,  
	"HiG_12" )


# refer to 11.2.0_weighted_graph_attack_robustness.R
s = "combined"
file = paste0('GSE175634_STRING_graph_perState_simplified_',s,'weighted.rds')
graph_list <- readRDS( file)  
	
names(graph_list)
 # [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"           "HiG_5"          
 # [7] "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"   
# [13] "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"     
# [19] "CTS_endoderm"    "CTS_CP"          "CTS_CP.1"    

GS =  lapply(graph_list, function(x) V(x)$name)
lengths(GS)
    #       HiG_0           HiG_1           HiG_2          HiG_CP           HiG_4 
    #         325            1155             242             132             104 
    #       HiG_5           HiG_6           HiG_7      HiG_muscle           HiG_9 
    #         339             106             104             118            1467 
    #      HiG_10    HiG_endoderm          HiG_12   HiGCTS_muscle HiGCTS_endoderm 
    #          70              42              22              20               7 
    #   HiGCTS_CP     HiGCTS_CP.1      CTS_muscle    CTS_endoderm          CTS_CP 
    #          29              27              61              64              72 
    #    CTS_CP.1 
    #          77 
GS[[1]]

'TP53I11' %in% GS[['HiGCTS_CP.1']]  # FALSE
'TP53I11' %in% GS[['CTS_CP.1']]  # TRUE
'TP53I11' %in% GS[['HiG_CP']]  # FALSE
for(i in 1:length(GS)) if('TP53I11' %in% GS[[i]]) cat(names(GS)[i], '\t')
# HiG_1   HiG_5   HiG_9   CTS_CP.1

# refer to 3.3_motif_enrich_coding_promoters.R,  line 378-383
TF_targeting_TP53I11_CM = c("NFYB", "NFYA", "NFYC", "SP1", "SP2", "SP3", "SP4", "SP5", "RXRA", "HNF4A", "HNF4G", "NR2F1", "NR2F2", "PDX1", "LHX2", "EN2", "PRRX1")
intersect(TF_targeting_TP53I11_CM, GS[['HiG_1']])   # [1] "NR2F1" "NR2F2" "PRRX1"
intersect(TF_targeting_TP53I11_CM, GS[['HiG_5']])   # 0
intersect(TF_targeting_TP53I11_CM, GS[['HiG_muscle']])  # 0
## therefore, "ISL1" "NR2F1" "NR2F2" "PRRX1" co-target TP53I11 in CP to pull towards CF sublineage


 
 
#########################################################
###### load published gene sets ############
#########################################################

CHD = readRDS( file='D:/projects/DS/result/CHD/CHD_Cilia_Genelist.rds')
GS[['validated_CHD']] <- CHD[["Griffin2023_PCGC"]]
GS[['curated_CHD']] <- CHD[["Griffin2023_PCGC_AllCurated"]]
GS[['validated_CHD']]
######################################
## GSE245713 Fig 1F/STable 1
## Epigenetic delineation of the earliest cardiac lineage segregation by single-cell multi-omics
JCF_SHF_path = 'F:/projects/scRNA/data/GSE245713_JCF/doc/'
XieStb1 = read.delim(file=paste0(JCF_SHF_path,'STable1.txt'), sep='\t', header=T)
dim(XieStb1) #[1] 98  3
head(XieStb1, 3)
    # Gene Cluster  X
# 1  Mixl1  Common NA
# 2  Tdgf1  Common NA
# 3      T  Common NA
table(XieStb1$Cluster)
      # Common JCF specific SHF specific 
          # 29           39           30 
GS[['Xie2025_JCF_specific']]  <- XieStb1 %>% filter(Cluster == "JCF specific") %>% pull(Gene) %>% toupper()
GS[['Xie2025_SHF_specific']]  <- XieStb1 %>% filter(Cluster == "SHF specific") %>% pull(Gene) %>% toupper()
GS[['Xie2025_JCF_n_SHF']] <- XieStb1 %>% filter(Cluster == "Common") %>% pull(Gene) %>% toupper()


###################################
## Maven 2023 paper
## refer to D:\projects\DS\source\GSE195476_ISL1\ISL1_GS.R
ISL1_set = readRDS(file='D:/projects/DS/result/GSE195476_ISL1/ISL1_set.rds')
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
 # ISL1_targets_MNP 
              # 119
names(ISL1_set) = paste0('Maven2023_', names(ISL1_set))
GS = c(GS, ISL1_set)


#######################################################
## Gao 2019 paper: 
## refer toGSE80383_Isl1/Isl1_GS.R generated genesets, mapped to  human symbols already
load('D:/projects/DS/result/GSE80383_Isl1/ISL1_Mm2Hg.GS_uniqueSymbol_FC1.2_padj0.05.RData')
names(ISL1)
#  [1] "Isl1.embryo.bound"         "Isl1.iPSC.bound"           "Isl1KO.E8.75.up"           "Isl1KO.E8.75.dn"          
#  [5] "Isl1.E8.75.DEG"            "Isl1KO.E10.5RV.OFT.up"     "Isl1KO.E10.5RV.OFT.dn"     "Isl1.E10.5RV.OFT.DEG"     
#  [9] "Isl1.target"               "Isl1.target.E8.75"         "Isl1.target.E10.5"         "Isl1KO.CP.up"             
# [13] "Isl1KO.CP.dn"              "Isl1.CP.DEG"               "Isl1.target.CP"            "iPSC.open"                
# [17] "E8.75.open"                "Isl1.E8.75DE.pooled.bound" "Isl1.E10.5DE.pooled.bound"
unique_ID = c(1:4, 6:7, 10:13, 15:17)
names(ISL1) = paste0('Gao2019_', names(ISL1))
GS = c(GS, ISL1[unique_ID])



#########################################################
## Ameen2022 paper
## refer to F:\projects\scATAC\source\GSE181346_heart\scATACseq\ATAC_GS.R
ATAC_gene_1k = readRDS(file='F:/projects/scATAC/results/GSE181346_heart/scATACseq/ATAC_gene_1k.rds')
ATAC_intergenic_5k = readRDS(file='F:/projects/scATAC/results/GSE181346_heart/scATACseq/ATAC_intergenic_5k.rds')
names(ATAC_gene_1k)
# [1] "allCells"        "arteries"        "atrial"          "capillary"      
#  [5] "cfdev"           "cfmature"        "cfwk6"           "endocardium"    
#  [9] "epicardium"      "fiibroblast_all" "lymph"           "myocardium"     
# [13] "neuralcrest"     "pericytes"       "smc_all"         "smcdev"         
# [17] "smcwk19"         "smcwk6"          "valveearly"      "valvelate"      
# [21] "veins"           "ventricular"    
names(ATAC_gene_1k) = paste0('ATAC_gene_1k_', names(ATAC_gene_1k))
names(ATAC_intergenic_5k) = paste0('ATAC_intergenic_5k_', names(ATAC_intergenic_5k))
GS = c(GS, ATAC_gene_1k, ATAC_intergenic_5k)
 



#########################################################
###### load the identified CTS gene sets ############
#########################################################

subfolds <- c('hESC_Bargaje2017', 'lung_Treutlein2014', 
			  'gastrulationE8.25_Pijuan-Sala2019','gastrulationE8.25_Ibarra-Soria2018',
			  'EB_Zhao2019','simulated_EMT',
			  '2024_3kHVG/')
subsubfolds <- c('C_Consensus_929cells', 'C_Leiden_0.4', 
					 'C_SNNGraph_allcells', 'subcelltype',
					 'C_SNNGraph_k5', 'C_cell_type',
					 'BioTIP')

db=7
	inputDir = ifelse(db<7,
		"E:/Git_Holly/scRNAseq_examples/result", 
		"F:/projects/scRNA/results/GSE175634_iPSC_CM/")

load(file=paste0(inputDir, subfolds[db] ,'/',subsubfolds[db ],'/CTS_Lib_Scaledata.RData'))
lengths(CTS.Lib.Symbol)
	# muscle   endoderm      CP/CF         CM         CF         CP    mix pro
		# 63         66         77        113        109         76         69
# endoderm.1       CP.1
	   # 136         87
if(db==7) {
		sig.DNB = c(T, T, T, T, T, T, F, F, T)
		sig.IC = c(T, T, F, F, T, T, NA, NA, T)
		sig.deltaIC = c(T, T, NA, NA, F, T, NA, NA, T)
		CTS = CTS.Lib.Symbol[which(sig.DNB==T & sig.IC==T & sig.deltaIC==T)]
		}
lengths(CTS)
  # muscle endoderm       CP     CP.1 
      # 63       66       76       87 

'ARL4C' %in% GS[["Xie2025_JCF_specific"]] #[1] TRUE
'ARL4C' %in% CTS[['CP']]  #[1] TRUE

#########################################################
## for each CTS, report integrated table
#########################################################
annotate_genesets <- function(genes, genesets) {
  df <- as.data.frame(sapply(genesets, function(gs) genes %in% gs))
  rownames(df) <- genes
  df
}

for(i in 1:length(CTS)) {
       tmp = annotate_genesets(CTS[[i]], GS)
       write.table(tmp, file=paste0('CTS_', names(CTS)[i], '_overlap_with_published_genesets.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=FALSE)
      print(paste0('CTS_', names(CTS)[i], ' done'))
      print(apply(tmp, 1, function(x) sum(x)) %>% sort %>% tail(10))
}
# manually verify
# i=3
# tmp = annotate_genesets(CTS[[i]], GS)
# 'ARL4C' %in% rownames(tmp) # TRUE
# tmp['ARL4C', "Xie2025_JCF_specific"]

# int='CP'
# file = paste0('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9/CTS_',int,'_overlap_with_published_genesets.txt')
# CTS_table_report = read.delim(file=file, sep='\t', header=TRUE)
# CTS_table_report['ARL4C',"Xie2025_JCF_specific"] # TRUE
# df = CTS_table_report[,c( "curated_CHD" , "Xie2025_JCF_specific" , "Xie2025_SHF_specific" , "Xie2025_JCF_n_SHF")]              
# head(df)
# df$JCF = df$Xie2025_JCF_specific | df$Xie2025_JCF_n_SHF
# df['ARL4C',c("JCF","Xie2025_JCF_specific" )]
        # JCF Xie2025_JCF_specific
# ARL4C TRUE                 TRUE

#  [1] "CTS_muscle done"
#  TNNT2 PDLIM3 ATP2B4 SLC8A1 PDLIM5   WEE1   NPTN   CNN1 FILIP1  HSPB1 
#     43     43     45     45     47     47     47     51     53     57 
# [1] "CTS_endoderm done"
#   GSTK1  CAMK2D  CDKN1C  LGALS3     CLU   MGST2 FAM107B S100A10 TMEM141    UACA 
#      38      39      39      40      42      45      46      47      47      51 
# [1] "CTS_CP done"
# ARRDC3   CTSV   NID2  ARL4C   HEY1  MPZL1   BMP4   EMP2  SFRP1  KDM6B 
#     46     46     48     49     53     54     54     54     56     56 
# [1] "CTS_CP.1 done"
#   ISL1  DUSP6   NAV1 CITED2  BMPER  CGNL1    ID1  MEIS1  PLCE1 ZNF608 
#     45     45     46     47     48     49     51     54     54     55 
