wd = '/Users/felixyu/Documents/GitHub/TIPS/examples/IbarraSoria2018/'
setwd(paste0(wd, "results/GSE181346_heart_scATAC/"))

library(EnrichedHeatmap)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(purrr)
library(SingleCellExperiment)

heatmap_coding_target_only = TRUE
input_path = "../../data/"

species = 'mouse'   # we upcase gene signature symbols to match the gene in the annotation datasets (e.g. TF_human, pseudo scATAC-seq accessibility data)
# refer to code 12.0_rank_by_PageRank_BC.R
seed_TF = c("MEF2C", "GATA4", "MSX2")

########################################################
##  input 1   dataset sepeciric ---
 
db_specifc_result_path = paste0(wd, 'results/')
db_specifc_input_path = paste0(wd, 'data/')
shared_path <- paste0(wd, "../Shared_Data/")
# db_specifc_CTS_path = 'E:/Git_Holly/BioTIP/examples/result/gastrulationE8.25_Pijuan-Sala2019/C_SNNGraph_allcells/' # TODO
db_specifc_CTS_path = paste0(wd, 'data/')

### the gene expression dataset to calculate the delta edge for newly added edges between TF and CTS_target 
# sce_path= 'D:/projects/DS/result/GSE175634_iPSC_CM/2024_3kHVG/BioTIP/'
load(paste0(db_specifc_input_path, 'sce_16subtype.RData'))
sce # 4000 11039 
names(colData(sce))
# [1] "label"       "celltype"    "subcelltype"                   

celltype_col = "subcelltype"  
CP_cluster = 'cardiac.a'
CM_cluster = 'cardiac.c'
CF_cluster = 'extraembryonicMesoderm'
CMES_cluster = 'mixedMesoderm.a'

CTS_ID = 'cardiac.a'
CTS_name = paste0('CTS_', CTS_ID)

########################################################
##  input 2.1   hg38 coding gene symbols & TF annotations  -- shared ---
## the v38 gtf was doenloaded from on 3/12/2021
## https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz

coding_genes = readRDS(file = paste0(shared_path, 'coding_genes.rds'))  %>% unique
length(coding_genes) # 19930
names(coding_genes) = NULL
head(coding_genes, 3)
# [1]  "OR4F5"  "OR4F29" "OR4F16"

## get the TF database 
library(dorothea)
packageVersion('dorothea') # ‘1.8.0’
# human or mouse
data(dorothea_hs, package = "dorothea")
# data(dorothea_mm, package = "dorothea")
TF_human <- unique(dorothea_hs$tf)  # gene symbols
length(TF_human) # 1333


# CHD = readRDS( file='D:/projects/DS/result/CHD/CHD_Cilia_Genelist.rds')
CHD = readRDS( paste0(shared_path, 'CHD_Cilia_Genelist.rds'))
names(CHD)
CHD = CHD$Griffin2023_PCGC_AllCurated
length(CHD) #295

########################################################
##  input 2.1   maps -- shared ---
# bw_dir = 'F:/projects/scATAC/data/GSE181346_heart/GSE181346_processed_data/'
maps = read.delim(file=paste0(shared_path, 'readme_filename_map_xy.txt'), sep='\t', header=TRUE, comment.char='#', stringsAsFactors=FALSE)
dim(maps) #[1] 44  5
head(maps)
#           file_name                      identity NA. identity_simple category
# 1               all                    pseudobulk  NA      pseudobulk         
# 2 hft_allFibroblast pesudobulk_cardiac_fibroblast  NA   pesudobulk_CF         
# 3        hft_allSMC      pseudobulk_smooth_muscle  NA  pseudobulk_SMC         
# 4      hft_arteries                      arteries  NA        arteries         
# 5        hft_atrial         atrial_cardiomyocytes  NA       atrial_CM  PCW8_CM
# 6     hft_capillary                     capillary  NA       capillary  
# maps = maps[match(bw_files, paste0(maps$file_name, '.bw')),]
# dim(maps)
maps = maps[grepl('hft_', maps$file_name), ]
dim(maps)
unique(maps$identity)  # 21  5
#  [1] "pesudobulk_cardiac_fibroblast"  "pseudobulk_smooth_muscle"      
#  [3] "arteries"                       "atrial_cardiomyocytes"         
#  [5] "capillary"                      "pre_Cardiac_fibroblast"        
#  [7] "maturecardiac_fibroblast"       "cardiac_fibroblast_progenitors"
#  [9] "endocardium"                    "epicardium"                    
# [11] "lymphatic_endothelial"          "early_cardiomyocyte"           
# [13] "neural_crest"                   "pericytes"                     
# [15] "pre_smooth_muscle"              "smooth_muscle"                 
# [17] "Outflow_tract"                  "Fibroblast_like_1"             
# [19] "Fibroblast_like_2"              "Venal_endothelial"             
# [21] "Ventricular_cardiomyocytes"   
table(maps$category)
               # early_CP_pool           PCW19_CF PCW19_CF&PCW19_SMC 
                 # 8                  3                  1                  1 
          # PCW19_CM          PCW19_SMF            PCW6_CM            PCW8_CF 
                 # 1                  2                  1                  2 
           # PCW8_CM           PCW8_SMF 
                 # 1                  1 

########################################################
##  input 2.2 published gene openness score -- shared ---
## gene-level Marker genes (gene score) for fetal heart scATAC-seq clusters using differential test, published by authors 
# the table S3
library(openxlsx)
annotation_vitro = read.xlsx(xlsxFile = paste0(input_path, 'Ameen2022cell-supplement-10.xlsx'), sheet = 4)  

dim(annotation_vitro) #[1] 5262    8
head(annotation_vitro)
#  Chromosomes     start       end Gene.name   Log2FC      FDR  MeanDiff celltype
# 1        chr3  32818018  32897826    TRIM71 1.646251 3.25e-67 1.1276613     iPSC
# 2        chr4  92303622  93774556     GRID2 1.202665 1.59e-62 1.6245711     iPSC
# 3        chr4  92277075  92268767 LNCPRESS2 2.492502 1.30e-61 0.6258600     iPSC
# 4        chr3 109337572 109326141     DPPA4 2.435096 4.68e-58 0.5687305     iPSC
# 5        chr3  54639857  54632122      ESRG 3.116906 5.13e-58 0.6729837     iPSC
# 6        chr5 147870682 147882191   SCGB3A2 2.272151 2.51e-57 0.6195138     iPSC
table(annotation_vitro$celltype)  # No CFP accessibility data !!!! 
#    i-CF      i-CM      i-CP      i-EC     i-EPC     i-Mes i-Mes-End    i-MyoF     i-pCM     i-SMC      iPSC 
#       351       511       468       549       205       677       425        95       371       179       769 
#  iPSC-Mes 
#       662 

Ameen_annotation_list = list()
for (i in unique(annotation_vitro$celltype)) {
    Ameen_annotation_list[[i]] = annotation_vitro$Gene.name[which(annotation_vitro$celltype == i)]
}
lengths(Ameen_annotation_list)
#    iPSC  iPSC-Mes     i-Mes      i-CP    i-MyoF i-Mes-End     i-pCM      i-CM      i-EC     i-EPC     i-SMC 
#       769       662       677       468        95       425       371       511       549       205       179 
#      i-CF 
#       351 


# annotation_dir = 'F:/projects/scATAC/data/GSE181346_heart/doc/'
# annotation_fetal = read.xlsx(xlsxFile = paste0(annotation_dir, 'media-1.xlsx'), sheet = 4) 
# annotation_fetal = read.xlsx(xlsxFile = paste0(input_path, 'Ameen2022cell_media-1.xlsx'), sheet = 4)  
# annotation_fetal_0 = read.xlsx(xlsxFile = paste0(input_path, 'Ameen2022cell_media-1.xlsx'), sheet = 4)  
annotation_fetal = read.xlsx(xlsxFile = paste0(input_path, 'Ameen2022cell_Table1_sheet3_markergene_Laksshman2026update.xlsx'), sheet = 1)  
dim(annotation_fetal) #[1] 3536    9
head(annotation_fetal)
  # Chromosomes     start       end strand     name   Log2FC          FDR  MeanDiff Cluster
# 1       chr14  23376808  23379772      1    CMTM5 3.754208 1.021032e-34 3.3962760     aCM
# 2        chr1  11841017  11848079      1 NPPA-AS1 2.902697 6.194079e-30 2.2110155     aCM
# 3       chr14  23408277  23381990      2     MYH6 2.116923 1.484557e-23 1.1611431     aCM
# 4       chr10  86668449  86736068      1     LDB3 1.530406 4.418635e-22 1.6215322     aCM
# 5       chr15  73369264  73319859      2     HCN4 1.619006 6.473660e-22 0.9398254     aCM
# 6        chr2 178830802 178525989      2      TTN 1.202148 1.148231e-20 1.2172846     aCM
colnames(annotation_fetal)[9] = 'Cluster'
subset(annotation_fetal, name=='PRRX1')
     # Chromosomes     start       end strand  name   Log2FC          FDR  MeanDiff Cluster
# 1814        chr1 170662728 170739419      1 PRRX1 1.626475 9.563515e-16 0.7403641     FB1
# 2141        chr1 170662728 170739419      1 PRRX1 1.461361 6.203403e-10 0.5604366     FB2
subset(annotation_fetal, name=='ISL1')
     # Chromosomes    start      end strand name   Log2FC          FDR  MeanDiff Cluster
# 1918        chr5 51383391 51394738      1 ISL1 1.119516 1.207509e-05 0.3203454     FB1
# 2482        chr5 51383391 51394738      1 ISL1 1.159836 3.994425e-06 0.3741227     OFT


table(annotation_fetal$Cluster)  # No CFP accessibility data !!!! 
 # aCM    aEC    Cap     CF    CFP    eCM  Endo1  Endo2    EPC    FB1    FB2    lEC     NC 
 # 262    221    217    141    206    263    399     22     81    324     59     83    143 
 # OFT     PC  preCF preSMC    vCM    vEC   vSMC 
 # 262     52      2     75    476    176     72 

## rename Endo2 & Endo1, Endo1 to Endo, becasue we only have one Endo cluster in our data
# annotation_fetal$Cluster[which(annotation_fetal$Cluster=="Endo2")] = 'Endo'
# annotation_fetal$Cluster[which(annotation_fetal$Cluster=="Endo1")] = 'Endo'
annotation_fetal$Cluster[which(annotation_fetal$Cluster=="vSMC")] = 'SMC'

short_cluster_names = unique(annotation_fetal$Cluster)
length(unique(annotation_fetal$name)) #[1]  2823
 
for (i in short_cluster_names) {
    Ameen_annotation_list[[i]] = annotation_fetal$name[which(annotation_fetal$Cluster == i)]
}
lengths(Ameen_annotation_list)
#      iPSC  iPSC-Mes     i-Mes      i-CP    i-MyoF i-Mes-End     i-pCM      i-CM      i-EC     i-EPC     i-SMC 
    #   769       662       677       468        95       425       371       511       549       205       179 
    #  i-CF       aEC       aCM       Cap        CF       FB1       eCM      Endo       EPC       FB2       lEC 
    #   351       221       262       217       141       324       263       421        81        59        83 
    #    NC       OFT     preCF        PC    preSMC       vCM       vEC       SMC 
    #   143       262       208        52        75       476       176        72 

names(Ameen_annotation_list)

unique(maps$identity_simple)
#  [1] "pesudobulk_CF"         "pseudobulk_SMC"        "arteries"             
#  [4] "atrial_CM"             "capillary"             "pre_CF"               
#  [7] "matureCF"              "CF_progenitors"        "endocardium"          
# [10] "epicardium"            "lymphatic_endothelial" "early_CM"             
# [13] "neural_crest"          "pericytes"             "pre_SMC"              
# [16] "SMC"                   "Outflow_tract"         "Fibroblast_like_1"    
# [19] "Fibroblast_like_2"     "Venal_endothelial"     "Ventricular_CM"     
maps$cluster_abb = maps$identity_simple
maps$cluster_abb[which(maps$identity_simple=="arteries")] = 'aEC'
maps$cluster_abb[which(maps$identity_simple=="atrial_CM")] = 'aCM'
maps$cluster_abb[which(maps$identity_simple=="capillary")] = 'Cap'
maps$cluster_abb[which(maps$identity_simple=="pre_CF")] = 'preCF'
maps$cluster_abb[which(maps$identity_simple=="matureCF")] = 'CF'
maps$cluster_abb[which(maps$identity_simple=="CF_progenitors")] = 'CFP'
maps$cluster_abb[which(maps$identity_simple=="endocardium")] = 'Endo'
maps$cluster_abb[which(maps$identity_simple=="epicardium")] = 'EPC'
maps$cluster_abb[which(maps$identity_simple=="lymphatic_endothelial")] = 'lEC'
maps$cluster_abb[which(maps$identity_simple=="early_CM")] = 'eCM'
maps$cluster_abb[which(maps$identity_simple=="neural_crest")] = 'NC'
maps$cluster_abb[which(maps$identity_simple=="pericytes")] = 'PC'  
maps$cluster_abb[which(maps$identity_simple=="pre_SMC")] = 'preSMC'
maps$cluster_abb[which(maps$identity_simple=="SMC")] = 'SMC'
maps$cluster_abb[which(maps$identity_simple=="Outflow_tract")] = 'OFT'
maps$cluster_abb[which(maps$identity_simple=="Fibroblast_like_1")] = 'FB1'
maps$cluster_abb[which(maps$identity_simple=="Fibroblast_like_2")] = 'FB2'
maps$cluster_abb[which(maps$identity_simple=="Venal_endothelial")] = 'vEC'
maps$cluster_abb[which(maps$identity_simple=="Ventricular_CM")] = 'vCM'


setdiff( maps$cluster_abb , names(Ameen_annotation_list))  # [1] "pesudobulk_CF"  "pseudobulk_SMC" "Endo"    

names(Ameen_annotation_list)
 # [1] "iPSC"      "iPSC-Mes"  "i-Mes"     "i-CP"      "i-MyoF"    "i-Mes-End" "i-pCM"     "i-CM"     
 # [9] "i-EC"      "i-EPC"     "i-SMC"     "i-CF"      "aCM"       "aEC"       "Cap"       "CF"       
# [17] "CFP"       "eCM"       "Endo1"     "Endo2"     "EPC"       "FB1"       "FB2"       "lEC"      
# [25] "NC"        "OFT"       "PC"        "preCF"     "preSMC"    "vCM"       "vEC"       "SMC"     

colnames(annotation_vitro)[4] = 'name'
colnames(annotation_vitro)[8] = 'Cluster'

tmp = intersect(colnames(annotation_vitro), colnames(annotation_fetal))
ATAC_anno_df = rbind(as.data.frame(annotation_vitro[,tmp]), as.data.frame(annotation_fetal[,tmp]))
dim(ATAC_anno_df) #[1] 8798    8
head(ATAC_anno_df)
#   Chromosomes     start       end      name   Log2FC      FDR  MeanDiff Cluster
# 1        chr3  32818018  32897826    TRIM71 1.646251 3.25e-67 1.1276613    iPSC
# 2        chr4  92303622  93774556     GRID2 1.202665 1.59e-62 1.6245711    iPSC
# 3        chr4  92277075  92268767 LNCPRESS2 2.492502 1.30e-61 0.6258600    iPSC
# 4        chr3 109337572 109326141     DPPA4 2.435096 4.68e-58 0.5687305    iPSC
# 5        chr3  54639857  54632122      ESRG 3.116906 5.13e-58 0.6729837    iPSC
# 6        chr5 147870682 147882191   SCGB3A2 2.272151 2.51e-57 0.6195138    iPSC
rm(annotation_vitro, annotation_fetal)

###############################################################
# following has been run in 23_atac_0ver_GaoChIP_ATAC_CTS.R     NOT repeat !! 
#############################
# inputDir = "F:/projects/scRNA/results/GSE175634_iPSC_CM/'2024_3kHVG/BioTIP/"
load(file=paste0(db_specifc_CTS_path, 'BioTIP.res.RData'))
names(res)
#[1] "CTS.candidate" "CTS.score"     "Ic.shrink"     "significant"  
class(res$CTS.candidate)  # list
lengths(res$CTS.candidate)
#     endothelial.b            cardiac.a presomiticMesoderm.a        endothelial.b 
#                 33                   37                   71                   32 

## refer to code Ibarra2018_BioTIP_threshold30_replot.R for details
# Ic_shrink must reach its maximum value in the corresponding cluster, ensuring cluster-specific dominance

# CTS = res$CTS.candidate[which(res$significant==TRUE & Top_cluster== names(Top_cluster))]
# lengths(CTS)

CTS = res$CTS.candidate[2:1]
class(CTS) # list
lengths(CTS)
#   cardiac.a endothelial.b 
#          37            33 
 

########################################################
##  input 3 -- data-deriven --- DEGs 
## we used wilcox test for human iPSC dataset, but t-test for mouse dataset as previosuly published
## load DEGs refer to 3.5_Filter_Plot_Marker_diffBar.R
               # Create a new environment
#load(paste0(inputDir, "2024_3kHVG/DEG.wc_padj0.01_score40.RData"), envir = new_object)  # Load .RData into this environment
DEG = readRDS(paste0(db_specifc_input_path, "DEG_perState_min.prop0.25_lfc0.6_FDFR0.05.rds"))
names(DEG)
#  [1] "blood"                  "cardiac.b"              "cardiac.c"              "endothelial.a"         
#  [5] "endothelial.c"          "endothelial.d"          "extraembryonicMesoderm" "mesodermProgenitors"   
#  [9] "mixedMesoderm.a"        "mixedMesoderm.b"        "pharyngealMesoderm"     "presomiticMesoderm.a"  
# [13] "presomiticMesoderm.b"   "somiticMesoderm"        "endothelial.b"          "cardiac.a"
class(DEG[[1]])  #[1] ""character"
head(DEG[[1]], 4)
# [1] "Cited4" "Gypa"   "Smim1"  "Klf1"  

	

########################################################
##  input 4 -- data-deriven --- STRING v12 PPI 
## get the whole string v12 genes
library(igraph)
# refer to 11.2.0_weighted_graph_attack_robustness.R
file = paste0(db_specifc_result_path, '/PPI_weight/IbarraSoria2018_STRING_graph_perState_simplified_combinedweighted.rds')
graph_list <- readRDS( file)  	
names(graph_list)
 # [1] "HiG_blood"                  "HiG_cardiac.b"              "HiG_cardiac.c"             
#  [4] "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [7] "HiG_extraembryonicMesoderm" "HiG_mesodermProgenitors"    "HiG_mixedMesoderm.a"       
# [10] "HiG_mixedMesoderm.b"        "HiG_pharyngealMesoderm"     "HiG_presomiticMesoderm.a"  
# [13] "HiG_presomiticMesoderm.b"   "HiG_somiticMesoderm"        "HiG_endothelial.b"         
# [16] "HiG_cardiac.a"              "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"          
# [19] "CTS_endothelial.b"          "CTS_cardiac.a"  

############## upcase the mouse symbols to match the gene in the annotation datasets (e.g. TF_human, pseudo scATAC-seq accessibility data)

if(species == 'mouse') {
    for(i in seq_along(graph_list)) V(graph_list[[i]])$name = toupper(V(graph_list[[i]])$name)
    for(i in seq_along(DEG)) DEG[[i]] = toupper(DEG[[i]])
    for(i in seq_along(CTS)) CTS[[i]] = toupper(CTS[[i]])
    rownames(sce) = toupper(rownames(sce))
    cat('uppercase the mouse symbols to match the gene in the annotation datasets (e.g. TF_human, pseudo scATAC-seq accessibility data) done\n')
}
graph_list[[1]]
DEG[[1]]
CTS[[1]]



################# creat ATAC_anno_df table and save it ####################
### step 1.2) building the binary annotation matrix 
# rebuild_mat = TRUE # Decided by 24.1, 24.2, and 24.3
fileName = paste0('binary_annot_',CTS_name,'_scATAC_Maven2023_gene_ISL1_v3.tsv')

if(rebuild_mat) {

	annotation_sub = subset(ATAC_anno_df, name %in% CTS[[CTS_ID ]])
	length(unique(CTS[[CTS_ID]]))  # 37
	dim(annotation_sub) #[1]48  8
	length(unique(annotation_sub$name)) # 22
	colnames(annotation_sub)
# [1] "Chromosomes" "start"       "end"         "name"        "Log2FC"      "FDR"         "MeanDiff"   
# [8] "Cluster"                    

	unique(maps$category)
	# [1] ""              "PCW8_CM"       "PCW8_CF"       "PCW19_CF"      "early_CP_pool"
	# [6] "PCW6_CM"       "PCW19_SMF"     "PCW8_SMF"      "PCW19_CM"   
	subset(maps, category=="early_CP_pool")$cluster_abb
	subset(maps, category =="PCW8_CF")$cluster_abb # [1] "preCF" "FB2"  
	# [1] "preCF" "OFT"   "FB1"   
    ## NOte that    
    ## 1) --- age of CFP and preCF was inconsisdtent in Fig1 and table S1  -- waiting for response from authors ------------------------
    ## 2) --- we add PCW6_CM into the pool for flag gene accessibility at CP_pool
    ## 3) --- we add CFP accessibile gene was missing from Table S3 -- waiting for response from authors ------------------------


   ## here, we annotation the accessibility according to the published main figure 1
  	annotation_sub$PCW6CP_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c('CFP', 'OFT', 'FB1', 'eCM',  'i-CP' )]), '1', '0')  
	annotation_sub$PCW8_CM_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c('aCM',     'i-CM')]), '1', '0')
	annotation_sub$PCW8_CF_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("preCF", "FB2","EPC",      'i-CF')]), '1', '0')
	annotation_sub$PCW8_SMC_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("preSMC","EPC",   'i-SMC')]), '1', '0')
	annotation_sub$PCW19_CM_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("vCM")]), '1', '0')
	annotation_sub$PCW19_CF_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("CF","EPC")]), '1', '0')
	annotation_sub$PCW19_SMC_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("SMC","EPC")]), '1', '0')   
# FOr E8.25 snapshot, cell are to early to show CF accessibility, thus the accessibility of PCW6 CF clusters are included for ‘CF_accessibility’ (FB1, preCF) 
    annotation_sub$PCW6_CF_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("CFP","FB1")]), '1', '0')
	annotation_sub$PCW6_SMC_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("OFT")]), '1', '0')   
	annotation_sub$PCW6_CM_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("eCM")]), '1', '0')   
    annotation_sub$iEPC_access = ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("i-EPC")]), '1', '0')   

	tmp = apply(annotation_sub[,-(1:9)], 2, function(x) sum(as.numeric(x))) 
	which(tmp==0)
	# 0

	annotation_sub$CMES_hi = ifelse(annotation_sub$name %in% DEG[[CMES_cluster]], '1', '0') 
	annotation_sub$CP_hi = ifelse(annotation_sub$name %in% DEG[[CP_cluster]], '1', '0') 
	annotation_sub$CM_hi = ifelse(annotation_sub$name %in% DEG[[CM_cluster]], '1', '0') 
	annotation_sub$CF_hi = ifelse(annotation_sub$name %in% DEG[[CF_cluster]], '1', '0') 

	int = c("CMES_hi" ,"CP_hi", "CM_hi", "CF_hi", "PCW6CP_access", 
				"PCW8_CM_access", "PCW19_CM_access", 
				"PCW8_CF_access","PCW19_CF_access" 
				,"PCW8_SMC_access" , "PCW19_SMC_access"
                ,"PCW6_CM_access", "PCW6_CF_access", "PCW6_SMC_access", "iEPC_access"
				# , "Maven2023_gene_PRRX1_up_E"  ,    "Maven2023_gene_PRRX1_up_T"  ,    "Maven2023_gene_PRRX1_up_L"    
				#  , "Maven2023_gene_PRRX1_dn_E"  ,    "Maven2023_gene_PRRX1_dn_T"  ,    "Maven2023_gene_PRRX1_dn_L"      
				#  , "Maven2023_gene_PRRX1_WT_d6CP"
				#  , "Gao2019_gene_PRRX1_E825E9.bound"
				#  , "Gao2019_gene_PRRX1.iCPC_CPC.bound"
				)    
	mat <- as.matrix(sapply(annotation_sub[,int], as.numeric))
	rownames(mat) <- annotation_sub$name
	dim(mat) #[1] 48   15
	length(unique(rownames(mat))) # 22

	## collapse transcript of the same geens if the same pattern in mat
	# Sum across transcripts per gene
	mat <- rowsum(mat, group = rownames(mat))
	# Convert “>=1 transcript present” back to 0/1
	mat <- (mat > 0) * 1  # TRUE/FALSE → 1/0
	dim(mat) #[1] 22 15

	## add back the missing symbols
	(missing =setdiff(CTS[[CTS_ID]], rownames(mat)))
# 	[1] "RGS5"    "KLF6"    "ARL4C"   "GATA4"   "ST3GAL5" "TSPAN12" "PTH1R"   "DPYSL3"  "SPATS2L" "TUBB2A" 
# [11] "NTNG1"   "DBT"     "LIPA"    "EEF2K"   "FAM83D" 

	missing_mat = matrix(0, nrow=length(missing), ncol=ncol(mat))
	rownames(missing_mat) = missing
	colnames(missing_mat) = colnames(mat)
	mat = rbind(mat, missing_mat)
	dim(mat)  #[1] 37 15

	# table for CTS.CP ==  high/low in CP ; high/low in CM;  high/low in CF; 
	# open in PCW6; open in PCW8 CM; open in PCW19 CM; open in PCW8 CF; open in PCW19 CF
	mat = as.data.frame(mat)
	mat[, paste0('CTS_', CTS_ID)] = ifelse(rownames(mat) %in% CTS[[CTS_ID]], '1', '0')
	mat[, 'CP_hi'] = ifelse(rownames(mat) %in% DEG[[CP_cluster]], '1', '0')
	mat[, 'CM_hi'] = ifelse(rownames(mat) %in% DEG[[CM_cluster]], '1', '0')
	mat[, 'CF_hi'] = ifelse(rownames(mat) %in% DEG[[CF_cluster]], '1', '0')
	mat[, 'CMES_hi'] = ifelse(rownames(mat) %in% DEG[[CMES_cluster]], '1', '0')

########################################################
##  input 4 -- shared --- published ISL1-CP binding 
### maven2023 ISL1 set
# refer to D:\projects\DS\source\GSE195476_ISL1\ISL1_GS.R
# in which we did:
# 1) define consistent peaks (present in at least 2 replciates)
# 2) annotatePeak to the nearest proximal genes (within ±1 kb of transcription start sites, excluding “Distal Intergenic” peaks, hg19) to define ISL1-bound targets.
ISL1_set = readRDS(file=paste0(shared_path, 'GSE195476_ISL1/ISL1_set.rds'))
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
# ISL1_set = ISL1_set[c("ISL1_WT_d18MNP", "ISL1_targets_MNP", "ISL1_WT_d6CP",
#         "ISL1_NKO_d6CP", "NKX25_NKO_d6CP", "NKX25_WT_d6CP" )]
names(ISL1_set) = paste0('Maven2023_gene_', names(ISL1_set))

########################################################
##  input 5 -- shared --- published ISL1-CP binding 
## Gao2019 Isl1 chipseq set
# ISL1 targets were defined by assigning each ISL1 ChIP-seq peak to the nearest gene TSS, as previously reported.
load(paste0(shared_path, 'GSE80383_Isl1/ISL1_Mm2Hg.GS_uniqueSymbol_FC1.2_padj0.05.RData'))
names(ISL1)
#  [1] "Isl1.embryo.bound"         "Isl1.iPSC.bound"           "Isl1KO.E8.75.up"           "Isl1KO.E8.75.dn"          
#  [5] "Isl1.E8.75.DEG"            "Isl1KO.E10.5RV.OFT.up"     "Isl1KO.E10.5RV.OFT.dn"     "Isl1.E10.5RV.OFT.DEG"     
#  [9] "Isl1.target"               "Isl1.target.E8.75"         "Isl1.target.E10.5"         "Isl1KO.CP.up"             
# [13] "Isl1KO.CP.dn"              "Isl1.CP.DEG"               "Isl1.target.CP"            "iPSC.open"                
# [17] "E8.75.open"                "Isl1.E8.75DE.pooled.bound" "Isl1.E10.5DE.pooled.bound"
names(ISL1)[2] = 'Isl1.iCPC_CPC.bound'
names(ISL1)[16] = 'iPSC_CPC.open'

 	
	mat[, 'Maven2023_gene_ISL1_up_E'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_up_E']], '1', '0')
	mat[, 'Maven2023_gene_ISL1_up_T'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_up_T']], '1', '0')
	mat[, 'Maven2023_gene_ISL1_up_L'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_up_L']], '1', '0')

	mat[, 'Maven2023_gene_ISL1_dn_E'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_dn_E']], '1', '0')
	mat[, 'Maven2023_gene_ISL1_dn_T'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_dn_T']], '1', '0')
	mat[, 'Maven2023_gene_ISL1_dn_L'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_dn_L']], '1', '0')

	mat[, 'Maven2023_gene_ISL1_WT_d6CP'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL1_WT_d6CP']], '1', '0')
	mat[, 'Gao2019_gene_Isl1_E825E9.bound'] = ifelse(rownames(mat) %in% ISL1[['Isl1.embryo.bound']], '1', '0')
	mat[, 'Gao2019_gene_Isl1.iCPC_CPC.bound'] = ifelse(rownames(mat) %in% ISL1[['Isl1.iCPC_CPC.bound']], '1', '0')

	dim(mat)  #[1] 37 25
	colnames(mat)
 # [1] "CMES_hi"                          "CP_hi"                            "CM_hi"                           
 # [4] "CF_hi"                            "PCW6CP_access"                    "PCW8_CM_access"                  
 # [7] "PCW19_CM_access"                  "PCW8_CF_access"                   "PCW19_CF_access"                 
# [10] "PCW8_SMC_access"                  "PCW19_SMC_access"                 "PCW6_CM_access"                  
# [13] "PCW6_CF_access"                   "PCW6_SMC_access"                  "iEPC_access"                     
# [16] "CTS_cardiac.a"                    "Maven2023_gene_ISL1_up_E"         "Maven2023_gene_ISL1_up_T"        
# [19] "Maven2023_gene_ISL1_up_L"         "Maven2023_gene_ISL1_dn_E"         "Maven2023_gene_ISL1_dn_T"        
# [22] "Maven2023_gene_ISL1_dn_L"         "Maven2023_gene_ISL1_WT_d6CP"      "Gao2019_gene_Isl1_E825E9.bound"  
# [25] "Gao2019_gene_Isl1.iCPC_CPC.bound"


	 
####################################################
### extract subnetworks and add ISL1-bound links ### 
####################################################
	dim(mat) #[1] 37  25

	##  step 1) binary annotation  
	## ### block the genes instead of cluster (new_strategy: must have ISL1-CP binding while in this step regardless promoter accessibility
	# USED to be consistent with mouse dataset !! 
	mat = as.data.frame(mat)
	mat[,'ISL1_CP_bound'] = ifelse( (mat[,'Gao2019_gene_Isl1_E825E9.bound'] == 1 | 
								mat[,'Gao2019_gene_Isl1.iCPC_CPC.bound'] == 1 | 
								mat[,'Maven2023_gene_ISL1_WT_d6CP'] == 1) , 1, 0)
	table(mat$ISL1_CP_bound)
	#  0  1 
	#  2  35        

	# one label per CTS&HiG module gene to prioritize candidates CTS member genes
	# 1) must be CP_hi, CM_hi or CF_hi
	# 2) must be ISL1_CP_bound
	mat[,'ISL1_CP_candidate'] = ifelse(mat[,'CP_hi'] == 1 & mat[,'ISL1_CP_bound']  == 1 , 1, 0)   
	# 3) must be open at PCW8_CM or PCW19_CM for CM_hi
	mat[,'ISL1_CM_candidate'] = ifelse(mat[,'CM_hi'] == 1  & mat[,'ISL1_CP_bound']  , 1, 0) 
	# 4) must be open at PCW8_CF or PCW19_CF for CF_hi
	mat[,'ISL1_CF_candidate'] = ifelse( mat[,'CF_hi'] == 1 & mat[,'ISL1_CP_bound'] , 1, 0) 
	if('SMC_hi' %in% colnames(mat)) mat[,'ISL1_SMC_candidate'] = ifelse( mat[,'SMC_hi'] == 1 & mat[,'ISL1_CP_bound'] , 1, 0) 
	

	dim(mat) # dim(mat)[1] 37  29

	write.table(mat, file=fileName, sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)  #!!!!!!!!!!!!

}  else {
  mat = read.table(fileName, sep='\t', header=T, check.names=FALSE)
  cat('mat is loaded\n')
}


########################################################
##  input 6 -- data-driven --- RcisTarget predicted PRRX1 targets 
library(RcisTarget)
packageVersion("RcisTarget")  # ‘1.18.2’
library(data.table)

data(package = "RcisTarget")  # list which motif-annotation objects you actually have
# data(motifAnnotations_hgnc_v9)  ## mouse: data(motifAnnotations_mgi)
# dim(motifAnnotations_hgnc_v9)  # [1] 163192      7

data(motifAnnotations_hgnc)
dim(motifAnnotations) #[1] 253096      8

motifAnnot <- motifAnnotations
dim(motifAnnot)  #[1] 253096      8



# takes a while, do NOT repeat !!!
if(!file.exists(file= "cisTarget_targets_in_all_CTS.rds")) {
	# --- 1) Load rankings database (downloaded from cisTarget resources) ---
	# https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/

    dbFile <- "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

    feather_path = paste0(db_specifc_input_path, "cistarget/")

    if(!file.exists(file= paste0(feather_path, dbFile))){
        url <- "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

        # Local folder + destination file
        dest_dir <- feather_path
        dest_file <- file.path(
        dest_dir,
        basename(url)
        )

        # Create folder if it doesn't exist
        if (!dir.exists(dest_dir)) {
        dir.create(dest_dir, recursive = TRUE)
        }

        options(timeout = 600)

        # Download
        download.file(
        url,
        destfile = dest_file,
        mode = "wb",
        method = "libcurl"
        )

        cat("Downloaded to:", dest_file, "\n")
    }

	motifRankings <- importRankings(paste0(feather_path, dbFile))  # reads .feather into a rankingRcisTarget object

	# --- 1) calcualte the enriched motifs amogn linegae-specific HiGs
    ## haryngeal/cardiopharyngeal mesoderm is best known for producing cardiac muscle (especially SHF myocardium) and branchiomeric skeletal muscle, not epicardium-derived fibroblasts.
    ## from epiblast cells exiting through the primitive streak; extraembryonic mesoderm migrates proximally (vs embryonic mesoderm migrating anteriorly as “wings”).
## Pharyngeal / cardiopharyngeal mesoderm (SHF-like CP): SHF/CPM regulators (often ISL1/TBX1/FGF10 programs), early cardiac progenitor flavor.
## Proepicardial / epicardial: WT1/TBX18/TCF21-type epicardial identity, then EMT to EPDC-like fibroblast/smooth muscle programs.
## Extraembryonic mesoderm: proximal/ExM signatures, often PODXL + KRT8/18 plus strong ECM/adhesion/cytoskeleton programs

 	cisTarget.res_HiG <- cisTarget(
		geneSets = DEG[c(CM_cluster, CF_cluster)],
		motifRankings = motifRankings,
		motifAnnot = motifAnnot,
		nesThreshold = 3
	)
	# Genes in the gene sets NOT available in the dataset: 
        # cardiac.c:      33 (7% of 448)
        # extraembryonicMesoderm:         23 (6% of 369)


	saveRDS(cisTarget.res_HiG, file="cisTarget_targets_in_two_HiGs.rds") 
	
	# # "MEF2C" "GATA4" "MSX2" reported by code 12.0_rank_by_PageRank_BC.R
    # ## which DEG motif is the key TF for CTS.cardiac.a??
    # tmp = cisTarget.res_HiG
	# for (i in seed_TF) {
	  # hit <- grepl(i, tmp$motif) |
			 # grepl(i, tmp$TF_highConf) |
			 # grepl(i, tmp$TF_lowConf)

	  # if (any(hit)) {
		# cat("\n====", i, "====\n")
		# print(tmp[hit, , drop = FALSE])
	  # }
	# }


	# --- 2) calcualte the enriched motifs among CTS genes 
	
	cisTarget.res <- cisTarget(
		geneSets = CTS,
		motifRankings = motifRankings,
		motifAnnot = motifAnnot,
		nesThreshold = 3
	)
	# Genes in the gene sets NOT available in the dataset: 
    #     endothelial.b:  3 (9% of 33)
    #     presomiticMesoderm.a:   8 (11% of 71)
    #     endothelial.b:  3 (9% of 33)

	saveRDS(cisTarget.res, file="cisTarget_targets_in_all_CTS.rds")
		# # res columns include TF_highConf/TF_lowConf + enrichedGenes (comma-separated)
		# # (See RcisTarget vignette + addSignificantGenes docs)
	summary(cisTarget.res$NES)	
#        Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   3.010   3.180   3.310   3.571   3.772   5.840 
	summary(cisTarget.res[which(cisTarget.res$geneSet==CTS_ID),]$NES)   
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  3.030   3.150   3.290   3.476   3.598   4.770 

    #(NES_threshold = quantile(cisTarget.res[which(cisTarget.res$geneSet==CTS_ID),]$NES, seq(0,1,0.1))['80%'])  #3.814 

	# ## which CTS itself in in the motif, TF_highConf, or TF_lowConf
    # tmp = subset(cisTarget.res, geneSet %in% c(CP_cluster))
	# for (i in CTS[[CTS_ID]]) {
	  # hit <- grepl(i, tmp$motif) |
			 # grepl(i, tmp$TF_highConf) |
			 # grepl(i, tmp$TF_lowConf)

	  # if (any(hit)) {
		# cat("\n====", i, "====\n")
		# print(tmp[hit, , drop = FALSE])
	  # }
	# }

}

