setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9/GSE181346_heart_scATAC')

library(EnrichedHeatmap)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(purrr)


gtf_dir = 'F:/projects/Share/gtf/'
gtf_file = 'gencode.v38.annotation.gtf'
gtf = read.table(paste0(gtf_dir, gtf_file), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gtf) = c("Chromosomes", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
dim(gtf) #  3150424       9
head(gtf, 3)
#   Chromosomes source    feature start   end score strand frame
# 1        chr1 HAVANA       gene 11869 14409     .      +     .
# 2        chr1 HAVANA transcript 11869 14409     .      +     .
# 3        chr1 HAVANA       exon 11869 12227     .      +     .
#                                                                                                                                                                                                                                                                                                                                                                                 attribute
# 1                                                                                                                                                                                                                              gene_id ENSG00000223972.5; gene_type transcribed_unprocessed_pseudogene; gene_name DDX11L1; level 2; hgnc_id HGNC:37102; havana_gene OTTHUMG00000000961.2;
# 2                                           gene_id ENSG00000223972.5; transcript_id ENST00000456328.2; gene_type transcribed_unprocessed_pseudogene; gene_name DDX11L1; transcript_type processed_transcript; transcript_name DDX11L1-202; level 2; transcript_support_level 1; hgnc_id HGNC:37102; tag basic; havana_gene OTTHUMG00000000961.2; havana_transcript OTTHUMT00000362751.1;
# 3 gene_id ENSG00000223972.5; transcript_id ENST00000456328.2; gene_type transcribed_unprocessed_pseudogene; gene_name DDX11L1; transcript_type processed_transcript; transcript_name DDX11L1-202; exon_number 1; exon_id ENSE00002234944.1; level 2; transcript_support_level 1; hgnc_id HGNC:37102; tag basic; havana_gene OTTHUMG00000000961.2; havana_transcript OTTHUMT00000362751.1;
gtf = gtf[which(gtf$feature == "gene"),]
dim(gtf) #60649     9

# get$name from gtf$attribute when gene_type is protein_coding
gtf= gtf[which(grepl("gene_type protein_coding", gtf$attribute)),]
dim(gtf) # 19955     9

coding_genes = sapply(gtf$attribute, function(x)  strsplit(x, "; gene_name ")[[1]][2] %>% unlist %>% unique)
coding_genes = sapply(coding_genes, function(x)  strsplit(x, "; ")[[1]][1] %>% unlist %>% unique)
length(coding_genes %>% unique) # 19930
names(coding_genes) = NULL
head(coding_genes, 3)
# [1]  "OR4F5"  "OR4F29" "OR4F16"


########################################################
## 
bw_dir = 'F:/projects/scATAC/data/GSE181346_heart/GSE181346_processed_data/'
maps = read.delim(file=paste0(bw_dir, 'readme_filename_map_xy.txt'), sep='\t', header=TRUE, comment.char='#', stringsAsFactors=FALSE)
dim(maps) #[1] 44  5
head(maps)
#           file_name                      identity NA. identity_simple category
# 1               all                    pseudobulk  NA      pseudobulk         
# 2 hft_allFibroblast pesudobulk_cardiac_fibroblast  NA   pesudobulk_CF         
# 3        hft_allSMC      pseudobulk_smooth_muscle  NA  pseudobulk_SMC         
# 4      hft_arteries                      arteries  NA        arteries         
# 5        hft_atrial         atrial_cardiomyocytes  NA       atrial_CM  PCW8_CM
# 6     hft_capillary                     capillary  NA       capillary  
maps = maps[match(bw_files, paste0(maps$file_name, '.bw')),]
dim(maps)
maps = maps[grepl('hft_', maps$file_name), ]
unique(maps$identity)
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
            #  early_CP_pool      PCW19_CF      PCW19_CM     PCW19_SMF      PCW6_CM       PCW8_CF       PCW8_CM      PCW8_SMF 
            # 97             3             1             1             2             1             2             1             1


## gene-level Marker genes (gene score) for fetal heart scATAC-seq clusters using differential test, published by authors 
# the table S3
library(openxlsx)
annotation_dir = 'F:/projects/scATAC/data/GSE181346_heart/doc/'
annotation_df = read.xlsx(xlsxFile = paste0(annotation_dir, 'media-1.xlsx'), sheet = 4)
dim(annotation_df) #[1] 3536    9
head(annotation_df)
#  Chromosomes    start      end strand   name   Log2FC      FDR  MeanDiff Cluster
# 1        chr1 34792998 34795747      1   GJA4 1.489500 1.18e-16 0.9494266     aEC
# 2       chr15 40929340 40939072      1   DLL4 1.251478 1.18e-16 1.3185304     aEC
# 3        chr1 34859816 34712737      2 SMIM12 1.380869 7.89e-14 0.8370249     aEC
# 4       chr15 40894430 40903975      1  VPS18 1.060934 7.89e-14 0.8473623     aEC
# 5       chr17 64020634 64002594      2  ICAM2 1.941813 1.24e-13 0.7323670     aEC
# 6       chr10 69088106 69104811      1   SRGN 1.910139 3.02e-13 0.3881559     aEC

subset(annotation_df, name=='ISL1')
#      Chromosomes    start      end strand name   Log2FC      FDR  MeanDiff Cluster
# 947         chr5 51383391 51394738      1 ISL1 1.119516 1.21e-05 0.3203454     FB1
# 2276        chr5 51383391 51394738      1 ISL1 1.159836 3.99e-06 0.3741227     OFT


table(annotation_df$Cluster)
#    aCM    aEC    Cap     CF    eCM  Endo1  Endo2    EPC    FB1    FB2    lEC     NC    OFT 
#    262    221    217    141    263    399     22     81    324     59     83    143    262 
#     PC  preCF preSMC    SMC    vCM    vEC 
#     52    208     75     72    476    176 
## rename Endo2 & Endo1, Endo1 to Endo, becasue we only have one Endo cluster in our data
annotation_df$Cluster[which(annotation_df$Cluster=="Endo2")] = 'Endo'
annotation_df$Cluster[which(annotation_df$Cluster=="Endo1")] = 'Endo'

short_cluster_names = unique(annotation_df$Cluster)
length(unique(annotation_df$name)) #[1]  2823
annotation = list()
for (i in short_cluster_names) {
    annotation[[i]] = annotation_df$name[which(annotation_df$Cluster == i)]
}
lengths(annotation)
#    aEC    aCM    Cap     CF    FB1    eCM  Endo2   EPC    FB2    lEC     NC    OFT 
#    221    262    217    141    324    263     421   81     59     83    143    262 
#  preCF     PC preSMC    vCM    vEC    SMC 
#    208     52     75    476    176     72 

## rename names(annotation) to match maps$identity_simple
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
maps$cluster_abb[which(maps$identity_simple=="CF_progenitors")] = 'preCF'
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


all(names(annotation) %in% maps$cluster_abb)  # TRUE
# TRUE

###############################################################
# following has been run in 23_atac_0ver_GaoChIP_ATAC_CTS.R<  NOT repeat !! 
#############################
inputDir = "F:/projects/scRNA/results/GSE175634_iPSC_CM/"
load(file=paste0(inputDir, '2024_3kHVG/', 'BioTIP/', 'CTS_Lib_Scaledata.RData'))
lengths(CTS.Lib.Symbol)
	# muscle   endoderm      CP/CF         CM         CF         CP    mix pro
		# 63         66         77        113        109         76         69
# endoderm.1       CP.1
	   # 136         87

sig.DNB = c(T, T, T, T, T, T, F, F, T)
sig.IC = c(T, T, F, F, T, T, NA, NA, T)
sig.deltaIC = c(T, T, NA, NA, F, T, NA, NA, T)
CTS = CTS.Lib.Symbol[which(sig.DNB==T & sig.IC==T & sig.deltaIC==T)]

lengths(CTS)
  # muscle endoderm       CP     CP.1 
      # 63       66       76       87 

#######################
## load DEGs refer to 3.5_Filter_Plot_Marker_diffBar.R
new_object <- new.env()                # Create a new environment
load(paste0(inputDir, "2024_3kHVG/DEG.wc_padj0.01_score40.RData"), envir = new_object)  # Load .RData into this environment
markers.up <- new_object[[ls(new_object)[1]]]
names(markers.up)[which(names(markers.up)=='3')] = 'CP'
names(markers.up)[which(names(markers.up)=='8')] = 'muscle'
names(markers.up)[which(names(markers.up)=='11')] = 'endoderm'
class(markers.up[[1]])  #[1] "data.frame"
head(markers.up[[1]], 4)
      # names   scores logfoldchanges pvals pvals_adj
# 1     L1TD1 261.4430            NaN     0         0
# 2     TERF1 247.4386            NaN     0         0
# 3 MIR302CHG 228.6239            NaN     0         0
# 4      PTMA 222.5636            NaN     0         0
	
	
 load(paste0(inputDir,'2024_3kHVG/DEG.wc_padj0.01_score40.RData'))
 lapply(DEG.wc, nrow) %>% unlist
   # 0    1    2    3    4    5    6    7    8    9   10   11   12 
 # 344 1224  257  152  110  350  110  117  119 1633   81   44   25 
 DEG = list()
 for(i in 1:length(DEG.wc)) DEG[[i]] = DEG.wc[[i]]$names %>% unique
 lengths(DEG)
#  [1]  344 1224  257  152  110  350  110  117  119 1633   81   44   25
names(DEG) = 0:(length(DEG)-1)
names(DEG)[which(names(DEG)=='3')] = 'CP'
names(DEG)[which(names(DEG)=='8')] = 'muscle'
names(DEG)[which(names(DEG)=='11')] = 'endoderm'

### maven2023 ISL1 set
# refer to D:\projects\DS\source\GSE195476_ISL1\ISL1_GS.R
# in which we did:
# 1) define consistend peaks (present in at least 2 replciates)
# 2) annotatePeak to the nearest proximal genes (within ±1 kb of transcription start sites, excluding “Distal Intergenic” peaks, hg19) to define ISL1-bound targets.
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
# ISL1_set = ISL1_set[c("ISL1_WT_d18MNP", "ISL1_targets_MNP", "ISL1_WT_d6CP",
#         "ISL1_NKO_d6CP", "NKX25_NKO_d6CP", "NKX25_WT_d6CP" )]
names(ISL1_set) = paste0('Maven2023_gene_', names(ISL1_set))

## Gao2019 Isl1 chipseq set
# ISL1 targets were defined by assigning each ISL1 ChIP-seq peak to the nearest gene TSS, as previously reported.
load('D:/projects/DS/result/GSE80383_Isl1/ISL1_Mm2Hg.GS_uniqueSymbol_FC1.2_padj0.05.RData')
names(ISL1)
#  [1] "Isl1.embryo.bound"         "Isl1.iPSC.bound"           "Isl1KO.E8.75.up"           "Isl1KO.E8.75.dn"          
#  [5] "Isl1.E8.75.DEG"            "Isl1KO.E10.5RV.OFT.up"     "Isl1KO.E10.5RV.OFT.dn"     "Isl1.E10.5RV.OFT.DEG"     
#  [9] "Isl1.target"               "Isl1.target.E8.75"         "Isl1.target.E10.5"         "Isl1KO.CP.up"             
# [13] "Isl1KO.CP.dn"              "Isl1.CP.DEG"               "Isl1.target.CP"            "iPSC.open"                
# [17] "E8.75.open"                "Isl1.E8.75DE.pooled.bound" "Isl1.E10.5DE.pooled.bound"
names(ISL1)[2] = 'Isl1.iCPC_CPC.bound'
names(ISL1)[16] = 'iPSC_CPC.open'


### understand CTS.CP.1 ######
annotation_sub = subset(annotation_df, name %in% CTS[['CP.1']])
length(unique(CTS[['CP.1']]))  # 87
dim(annotation_sub) #[1] 25  9
length(unique(annotation_sub$name)) # 17


colnames(annotation_sub)
# [1] "Chromosomes"                      "start"                           
#  [3] "end"                              "strand"                          
#  [5] "name"                             "Log2FC"                          
#  [7] "FDR"                              "MeanDiff"                        
#  [9] "Cluster"                          " 

unique(maps$category)
# [1] ""              "PCW8_CM"       "PCW8_CF"       "PCW19_CF"      "early_CP_pool"
# [6] "PCW6_CM"       "PCW19_SMF"     "PCW8_SMF"      "PCW19_CM"   
subset(maps, category=="early_CP_pool")$cluster_abb
# [1] "preCF" "OFT"   "FB1" 
subset(maps, category=="PCW8_CF")$cluster_abb # [1] "preCF" "FB2"  

annotation_sub$PCW6CP_access = ifelse(annotation_sub$name %in% unlist(annotation[c('preCF', 'OFT', 'FB1', 'PCW6_CM')]), '1', '0')
annotation_sub$PCW8_CM_access = ifelse(annotation_sub$name %in% unlist(annotation[c('aCM')]), '1', '0')
annotation_sub$PCW8_CF_access = ifelse(annotation_sub$name %in% unlist(annotation[c("preCF", "FB2")]), '1', '0')
annotation_sub$PCW8_SMC_access = ifelse(annotation_sub$name %in% unlist(annotation[c("preSMC")]), '1', '0')
annotation_sub$PCW19_CM_access = ifelse(annotation_sub$name %in% unlist(annotation[c("vCM")]), '1', '0')
annotation_sub$PCW19_CF_access = ifelse(annotation_sub$name %in% unlist(annotation[c("CF")]), '1', '0')
annotation_sub$PCW19_SMC_access = ifelse(annotation_sub$name %in% unlist(annotation[c("SMC")]), '1', '0')  # !!!!??? "PC"

tmp = apply(annotation_sub[,-1:9], 2, function(x) sum(as.numeric(x))) 
which(tmp==0)


int = c("CMES_hi" ,"CP_hi", "CM_hi", "CF_hi", "PCW6CP_access", 
         "PCW8_CM_access", "PCW19_CM_access", 
         "PCW8_CF_access","PCW19_CF_access" 
        # ,"PCW8_SMC_access" , "PCW19_SMC_access"
        # , "Maven2023_gene_ISL1_up_E"  ,    "Maven2023_gene_ISL1_up_T"  ,    "Maven2023_gene_ISL1_up_L"    
        #  , "Maven2023_gene_ISL1_dn_E"  ,    "Maven2023_gene_ISL1_dn_T"  ,    "Maven2023_gene_ISL1_dn_L"      
        #  , "Maven2023_gene_ISL1_WT_d6CP"
        #  , "Gao2019_gene_Isl1_E825E9.bound"
        #  , "Gao2019_gene_Isl1.iCPC_CPC.bound"
          )    
mat <- as.matrix(sapply(annotation_sub[,int], as.numeric))
rownames(mat) <- annotation_sub$name
dim(mat) #[1] 25  18
length(unique(rownames(mat))) # 17

## collapse transcript of teh same geens if the same pattern in mat
# Sum across transcripts per gene
mat <- rowsum(mat, group = rownames(mat))
# Convert “>=1 transcript present” back to 0/1
mat <- (mat > 0) * 1  # TRUE/FALSE → 1/0
dim(mat) #[1] 17  18

## add back the missing symbols
(missing =setdiff(CTS[['CP.1']], rownames(mat)))
#  [1] "SAMD11"     "FYB2"       "DAB1"       "WLS"        "IGSF3"      "C1orf21"    "NAV1"      
#  [8] "LPGAT1"     "GPATCH11"   "LRRTM1"     "ZC3H6"      "ZEB2"       "LRP2"       "SLC40A1"   
# [15] "ITM2C"      "HACD2"      "TIPARP"     "LINC00888"  "CHIC2"      "PTPN13"     "AC108052.1"
# [22] "CLGN"       "GYPE"       "GYPB"       "HAND2"      "CPLANE1"    "ENC1"       "HAPLN1"    
# [29] "RHOBTB3"    "ZNF608"     "VTRNA1-3"   "LIN28B"     "PERP"       "CITED2"     "RPL22P16"  
# [36] "KIAA1324L"  "DYNC1I1"    "SLC7A2"     "ENTPD4"     "SLC25A37"   "NKX3-1"     "CA2"       
# [43] "NTRK2"      "PLCE1"      "ART5"       "CCKBR"      "MTRNR2L8"   "AL035078.1" "TP53I11"   
# [50] "SC5D"       "CCND2"      "PTMS"       "ATN1"       "DUSP6"      "GPC6"       "PRTG"      
# [57] "ORC6"       "C16orf87"   "LAMA1"      "CHST9"      "PSTPIP2"    "JAG1"       "NINL"      
# [64] "ID1"        "AC078993.1" "USP27X"     "MT-TF"      "MT-TV"      "MT-TL1"     "MT-TS2"    
missing_mat = matrix(0, nrow=length(missing), ncol=ncol(mat))
rownames(missing_mat) = missing
colnames(missing_mat) = colnames(mat)
mat = rbind(mat, missing_mat)
dim(mat)  #[1] 87 18

# table for CTS.CP ==  high/low in CP ; high/low in CM;  high/low in CF; 
# open in PCW6; open in PCW8 CM; open in PCW19 CM; open in PCW8 CF; open in PCW19 CF

mat[, 'CTS.CP.1'] = ifelse(rownames(mat) %in% CTS[['CP.1']], '1', '0')
mat[, 'CP_hi'] = ifelse(rownames(mat) %in% DEG[['CP']], '1', '0')
mat[, 'CM_hi'] = ifelse(rownames(mat) %in% DEG[['5']], '1', '0')
mat[, 'CF_hi'] = ifelse(rownames(mat) %in% DEG[['1']], '1', '0')
mat[, 'SMC_hi'] = ifelse(rownames(mat) %in% DEG[['muscle']], '1', '0')
mat[, 'CMES_hi'] = ifelse(rownames(mat) %in% DEG[['4']], '1', '0')

mat[, 'Maven2023_gene_ISL1_up_E'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_up_E']], '1', '0')
mat[, 'Maven2023_gene_ISL1_up_T'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_up_T']], '1', '0')
mat[, 'Maven2023_gene_ISL1_up_L'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_up_L']], '1', '0')

mat[, 'Maven2023_gene_ISL1_dn_E'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_dn_E']], '1', '0')
mat[, 'Maven2023_gene_ISL1_dn_T'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_dn_T']], '1', '0')
mat[, 'Maven2023_gene_ISL1_dn_L'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL_dn_L']], '1', '0')

mat[, 'Maven2023_gene_ISL1_WT_d6CP'] = ifelse(rownames(mat) %in% ISL1_set[['Maven2023_gene_ISL1_WT_d6CP']], '1', '0')
mat[, 'Gao2019_gene_Isl1_E825E9.bound'] = ifelse(rownames(mat) %in% ISL1[['Isl1.embryo.bound']], '1', '0')
mat[, 'Gao2019_gene_Isl1.iCPC_CPC.bound'] = ifelse(rownames(mat) %in% ISL1[['Isl1.iCPC_CPC.bound']], '1', '0')

tmp = apply(mat[,-1:7], 2, function(x) sum(as.numeric(x))) 
which(tmp==0)

Is_protein_coding = ifelse(rownames(mat) %in% coding_genes, TRUE, FALSE)

## get the TF database 
library(dorothea)
# human or mouse
data(dorothea_hs, package = "dorothea")
# data(dorothea_mm, package = "dorothea")
TF_human <- unique(dorothea_hs$tf)  # gene symbols
# TF_mouse <- unique(dorothea_mm$tf)
Is_TF = ifelse(rownames(mat) %in% TF_human, TRUE, FALSE)


## get the whole string v12 genes
library(igraph)
PPI_dir = 'F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9'
# refer to 11.2.0_weighted_graph_attack_robustness.R
file = paste0(PPI_dir, '/GSE175634_STRING_graph_perState_simplified_combinedweighted.rds')
graph_list <- readRDS( file)  	
names(graph_list)
 # [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"           "HiG_5"          
 # [7] "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"   
# [13] "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"     
# [19] "CTS_endoderm"    "CTS_CP"          "CTS_CP.1"    

PPI_gene = V(graph_list[['HiGCTS_CP.1']])$name 
In_PPI = ifelse(rownames(mat) %in% PPI_gene, TRUE, FALSE)

CHD = readRDS( file='D:/projects/DS/result/CHD/CHD_Cilia_Genelist.rds')
names(CHD)
CHD = CHD$Griffin2023_PCGC_AllCurated
IS_CHD = ifelse(rownames(mat) %in% CHD, TRUE, FALSE)

# Create a named character vector for RowSideColors based on Is_protein_coding (1 bar only)
# RowSideColors <- ifelse(Is_protein_coding, "blue", "darkgrey")
# names(RowSideColors) <- rownames(mat)

# heatmap(
#   mat, Rowv = TRUE, Colv = NA,
#   scale = 'none', col   = c("white", "orange"),
#   margins = c(4, 4),
#   cexCol =  1,
#   cexRow = 0.7,    # reduce to fit more if necessary, or adjust as you want
#   labRow = rownames(mat),  # show all gene names (y-axis labels)
#   distfun = function(x) dist(x, method = "binary"),
#   hclustfun = function(d) hclust(d, method = "ward.D2"),
#   RowSideColors = RowSideColors
# )


library(pheatmap)

## example: build row annotation with 3 yes/no columns
row_anno <- data.frame(
  is_TF = factor(ifelse(Is_TF, TRUE, FALSE), levels = c(FALSE, TRUE), labels = c("no", "yes")),
  coding = factor(ifelse(Is_protein_coding, TRUE, FALSE), levels = c(FALSE, TRUE), labels = c("no", "yes")),
  CHD= factor(ifelse(IS_CHD, TRUE, FALSE), levels = c(FALSE, TRUE), labels = c("no", "yes")),
  PPI = factor(ifelse(In_PPI, TRUE, FALSE), levels = c(FALSE, TRUE), labels = c("no", "yes"))
)
rownames(row_anno) <- rownames(mat)

## same color scheme for all three bars: no = lightgrey, yes = blue
ann_colors <- list(
  CHD = c(no = "lightgrey", yes = "blue"),
  PPI = c(no = "lightgrey", yes = "blue"),
  is_TF = c(no = "lightgrey", yes = "blue")
 )

ann_colors[['coding']] = c(no = "lightgrey", yes = "blue")



# ### block the genes instead of hcluster
# mat = as.data.frame(mat)
# # one label per gene
# # 1) must be  CP_hi, CM_hi or CF_hi
# # 2) must be open at CP_pool for CP_hi PPI
# mat[,'CP_candidate'] = ifelse(mat[,'CP_hi'] == 1 & (mat[,'PCW6CP_access'] == 1  ), 1, 0)   # | (mat[,'PCW8_CM_access'] == 1  trial 1 not used !
# # 3) must be open at PCW8_CM or PCW19_CM for CM_hi
# mat[,'CM_candidate'] = ifelse(mat[,'CM_hi'] == 1 & (mat[,'PCW8_CM_access'] == 1 | mat[,'PCW19_CM_access'] == 1 ), 1, 0) 
# # 4) must be open at PCW8_CF or PCW19_CF for CF_hi
# mat[,'CF_candidate'] = ifelse(mat[,'CF_hi'] == 1 & (mat[,'PCW8_CF_access'] == 1 | mat[,'PCW19_CF_access'] == 1 ), 1, 0) 

# mat[,'ISL1_CP_bound'] = ifelse( (mat[,'Gao2019_gene_Isl1_E825E9.bound'] == 1 | 
                            # mat[,'Gao2019_gene_Isl1.iCPC_CPC.bound'] == 1 | 
                            # mat[,'Maven2023_gene_ISL1_WT_d6CP'] == 1) , 1, 0)
                            

## ### version 2: block the genes instead of hcluster (new_strategy: must have ISL1-CP binding while regardless promoter accessibility
# USED to be consistent with mouse dataset !! 
mat = as.data.frame(mat)
mat[,'ISL1_CP_bound'] = ifelse( (mat[,'Gao2019_gene_Isl1_E825E9.bound'] == 1 | 
                            mat[,'Gao2019_gene_Isl1.iCPC_CPC.bound'] == 1 | 
                            mat[,'Maven2023_gene_ISL1_WT_d6CP'] == 1) , 1, 0)
table(mat$ISL1_CP_bound)
#  0  1 
# 19 68        

# one label per gene
# 1) must be  CP_hi, CM_hi or CF_hi
# 2) must be open at CP_pool for CP_hi PPI
mat[,'CP_candidate'] = ifelse(mat[,'CP_hi'] == 1 & mat[,'ISL1_CP_bound']  == 1 , 1, 0)   
# 3) must be open at PCW8_CM or PCW19_CM for CM_hi
mat[,'CM_candidate'] = ifelse(mat[,'CM_hi'] == 1  & mat[,'ISL1_CP_bound']  , 1, 0) 
# 4) must be open at PCW8_CF or PCW19_CF for CF_hi
mat[,'CF_candidate'] = ifelse( mat[,'CF_hi'] == 1 & mat[,'ISL1_CP_bound'] , 1, 0) 


# ### a plot with less rows: -- only ISL1-binding coding genes are shown
# coding_target_only =TRUE
if (coding_target_only) {
	dim(mat) # dim(mat)[1] 87 24
	x = which(row_anno$coding=='yes' & mat$ISL1_CP_bound==1)
	length(x)  # 65
	mat = mat[x,]
	row_anno = row_anno[x,]
    row_anno = row_anno[, -which(colnames(row_anno) == "coding")]
}	


pattern <- paste0(mat[, "CP_candidate"], mat[, "CM_candidate"], mat[, "CF_candidate"])
table(pattern)
# pattern
# 000 001 010 011 100 101 110 111 
# 28  15   2  15   2   3   5  17 



## map pattern → block name
pattern_to_block <- c(
  "100" = "CP_only",
  "110" = "CP_CM",
  "101" = "CP_CF",
  "111" = "CP_CM_CF",
  "010" = "CM_only",
  "001" = "CF_only",
  "011" = "CM_CF",
  "000" = "none"
)

block <- pattern_to_block[pattern]

## choose order you want to see in heatmap
block_levels <- c("CP_only","CP_CM","CM_only","CP_CM_CF","CM_CF","CF_only","CP_CF","none")
block <- factor(block, levels = block_levels)
names(block) <- rownames(mat)   # keep rownames aligned

ord <- order(block)
cols_to_show <- c("CP_candidate",   "CP_hi"   ,   "PCW6CP_access"    , #"PCW8_CM_access"  , 
                "CM_candidate"   ,  "CM_hi"  ,    "PCW8_CM_access"  , "PCW19_CM_access",
                "CF_candidate"   ,  "CF_hi"  ,    "PCW8_CF_access"  , "PCW19_CF_access"  ,
                "ISL1_CP_bound", "Maven2023_gene_ISL1_WT_d6CP", 
                "Gao2019_gene_Isl1_E825E9.bound", "Gao2019_gene_Isl1.iCPC_CPC.bound",
                "Maven2023_gene_ISL1_up_E" ,  
                "Maven2023_gene_ISL1_up_T"  ,  "Maven2023_gene_ISL1_up_L"  ,  "Maven2023_gene_ISL1_dn_E"   ,
                "Maven2023_gene_ISL1_dn_T" ,   "Maven2023_gene_ISL1_dn_L"   
)
row_label_col <- ifelse(mat[, "ISL1_CP_bound"] == 1, "red", "black")

# Ensure all values in mat are numeric before plotting
mat_numeric <- as.data.frame(apply(mat, 2, as.numeric), row.names = rownames(mat))
p = pheatmap(
  mat_numeric[ord, cols_to_show],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = c("white", "orange"),
  annotation_row = row_anno,
  annotation_colors = ann_colors,
  margins = c(4, 4),
  fontsize = 10,
  fontsize_col = 6,
  fontsize_row = 6,    # reduce to fit more if necessary, or adjust as you want
  labRow = rownames(mat)
)

# pdf(file='heatmap_blocked_CTS.CP.1_scATAC_Maven2023_gene_ISL1.pdf', height=9)
# print(p)
# dev.off()
# write.table(mat, file='heatmap_blocked_CTS.CP.1_scATAC_Maven2023_gene_ISL1.tsv', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)

if(coding_target_only) pdf(file='heatmap_blocked_CTS.CP.1_scATAC_Maven2023_gene_ISL1_v3_coding_target.pdf', height=9) else 
pdf(file='heatmap_blocked_CTS.CP.1_scATAC_Maven2023_gene_ISL1_v3.pdf', height=9)
print(p)
dev.off()

if(coding_target_only) fileName='heatmap_blocked_CTS.CP.1_scATAC_Maven2023_gene_ISL1_v3_coding_target.tsv' else fileName='heatmap_blocked_CTS.CP.1_scATAC_Maven2023_gene_ISL1_v3.tsv'
write.table(mat, file=fileName, sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)

## for Fig 5, only teh ISL1 CP bound genes are shown
x = which(mat_numeric[ord,]$ISL1_CP_bound == 1)
p = pheatmap(
  mat_numeric[ord[x], cols_to_show],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = c("white", "orange"),
  annotation_row = row_anno,
  annotation_colors = ann_colors,
  margins = c(4, 4),
  fontsize = 10,
  fontsize_col = 6,
  fontsize_row = 6,    # reduce to fit more if necessary, or adjust as you want
  labRow = rownames(mat)
)
pdf(file='heatmap_blocked_CTS.CP.1_scATAC_Maven2023_gene_ISL1_v3_ISL1_CP_bound.pdf', height=9)
print(p)
dev.off()

# pdf(file='heatmap_blocked_CTS.CP.1_scATAC_Maven2023_gene_ISL1.pdf', height=9)
# print(p)
# dev.off()
# write.table(mat, file='heatmap_blocked_CTS.CP.1_scATAC_Maven2023_gene_ISL1.tsv', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)



##########################################################
## -- PPI of CTS[['CP.1']] subset that are ISL1-target & exclusively high in CM or CF 
library(BioNet)
 packageVersion('BioNet') # '1.56.0'
 library(igraph)
 library("STRINGdb")
 packageVersion('STRINGdb') # '2.14.0'
 library(tibble)
 
 string_db <- STRINGdb$new( version="12.0", species=9606,   # species= 10090 for mouse;    ?? for human 2021 version
                            score_threshold = 200, #!!!!!!!!!!!default is 200
							network_type="full", 
                            input_directory="F:/projects/Share/PPIN/STRING/STRINGdb/human/full/") 
 string_db
 # version 12.0  / 11.5
 # proteins: 19699    / for mouse: 22048
 # interactions:  7533072   / for mouse 6859804

names(graph_list)
#  [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"          
#  [6] "HiG_5"           "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"          
# [11] "HiG_10"          "HiG_endoderm"    "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm"
# [16] "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"      "CTS_endoderm"    "CTS_CP"         
# [21] "CTS_CP.1"   

graph_ISL1_list = graph_list[c("HiGCTS_CP.1", "HiG_5", "HiG_1")]
names(graph_ISL1_list) = c("CTS_CP.1", "HiG_CM", "HiG_CF")
# [1] "HiGCTS_CP.1" "CTS_CP.1"   
## add the ISL1-targeting as dashlines
ISL1_target = rownames(mat)[which(mat[,"ISL1_CP_bound"] == 1)]
length(ISL1_target)  # 66

GRN_node_list = list( HiG_CM = rownames(mat)[which(mat[,"CM_candidate"] == 1)],
                      HiG_CF = rownames(mat)[which(mat[,"CF_candidate"] == 1)],
                      CTS_CP.1 = rownames(mat)[which(mat[,"CP_candidate"] == 1)])

for (name in names(graph_ISL1_list)) {
        graph = graph_ISL1_list[[name]]
        if(name != 'CTS_CP.1') {
            graph = induced_subgraph(graph, V(graph)$name %in% GRN_node_list[[name]])
        }
        ## add the GRN_node if not in the graph
        GRN_node_not_in_graph = setdiff(GRN_node_list[[name]], V(graph)$name)
        if (length(GRN_node_not_in_graph) > 0) {
            # 1) add the vertices
            graph <- add_vertices(graph,
                                nv   = length(GRN_node_not_in_graph),
                                name = GRN_node_not_in_graph)

            # 2) add PPI edges involving these new nodes
            all_nodes_now <- V(graph)$name  

            # PPI edges where at least one endpoint is new
            map_df = string_db$map(data.frame(gene=all_nodes_now), "gene", removeUnmappedRows = TRUE )
            hits <- map_df$STRING_id
            ppi_new <- string_db$get_subnetwork(hits) 
	 
            if (nrow(ppi_new) > 0) {
                ## Convert STRING IDs back to gene symbols
                # STRING uses internal ids — map back using map_df
                id2gene <- setNames(map_df$gene, map_df$STRING_id)

                ppi_genes_A <- id2gene[ppi_new$from]
                ppi_genes_B <- id2gene[ppi_new$to]

                # Remove self-loops and NA
                keep <- (!is.na(ppi_genes_A)) &
                        (!is.na(ppi_genes_B)) &
                        (ppi_genes_A != ppi_genes_B)

                ppi_genes_A <- ppi_genes_A[keep]
                ppi_genes_B <- ppi_genes_B[keep]

                # Build edge list
                ppi_edges <- c(rbind(ppi_genes_A, ppi_genes_B))

                # Add PPI edges
                m_before_ppi <- ecount(graph)
                if(length(ppi_edges) > 0) {
                    graph <- add_edges(graph, ppi_edges)  #!!!!
                    new_ppi_edges <- seq(m_before_ppi + 1, ecount(graph))
                }
            }
        }

        ## ensure node of CHD is red
        V(graph)$color = ifelse(V(graph)$name %in% CHD, "red", "lightgrey")
        E(graph)$color <- "grey70"
        E(graph)$lty   <- "solid"   
        ## add the edges between the ISL1-targeting and the (new and old) nodes in the graph
        ISL1_target_in_graph = intersect(ISL1_target, V(graph)$name)
        if(length(ISL1_target_in_graph) > 0) {
            # remember how many edges there were BEFORE adding
            m_before <- ecount(graph)

            edges_vec <- c(rbind('ISL1', ISL1_target_in_graph))
            graph <- add_edges(graph, edges_vec)

            # indices of the newly added edges
            new_edges <- seq(m_before + 1, ecount(graph))
            
            # style ONLY those edges
            new_edges_color = rep("grey70", length(new_edges))
            tmp = grep('ISL1_up', colnames(mat))
            new_edges_color[apply(mat[ISL1_target_in_graph,tmp, drop = FALSE], 1, sum) > 0] = "blue"
            tmp = grep('ISL1_dn', colnames(mat))
            new_edges_color[apply(mat[ISL1_target_in_graph,tmp, drop = FALSE], 1, sum) > 0] = "orange"

            E(graph)$color[new_edges] <- new_edges_color      # make them blue
            E(graph)$lty[new_edges]   <- "dashed"    # and dashed
     }


    graph_ISL1_list[[name]] = graph
     
}

# saveRDS(graph_ISL1_list, file='PPI_graph_ISL1_GRN_prediction.rds')
saveRDS(graph_ISL1_list, file='PPI_graph_ISL1_GRN_prediction_v3.rds')


node_color = c('black', 'darkred', 'blue')
names(node_color) = c('CTS_CP.1', 'HiG_CM', 'HiG_CF')

library(gplots)
# pdf(file='PPI_graph_ISL1_GRN_prediction.pdf')
pdf(file='PPI_graph_ISL1_GRN_prediction_v3.pdf')
venn(list(CTS_CP.1 = V(graph_ISL1_list[[1]])$name, 
    HiG_CM = V(graph_ISL1_list[[2]])$name, 
    HiG_CF = V(graph_ISL1_list[[3]])$name))

for(name in names(graph_ISL1_list)) {
    graph = graph_ISL1_list[[name]]
    set.seed(123)                      # for reproducibility
    center <- which(V(graph)$name == "ISL1")
    lay <- layout_as_star(graph, center = center)

    gene_color= rep(node_color[name], length(V(graph)$name))
    gene_color[V(graph)$name %in% setdiff(V(graph_ISL1_list[['HiG_CM']])$name, V(graph_ISL1_list[['HiG_CF']])$name)] = "darkred"
    gene_color[V(graph)$name %in% setdiff(V(graph_ISL1_list[['HiG_CF']])$name, V(graph_ISL1_list[['HiG_CM']])$name)] = "blue"
    if(name != 'CTS_CP.1') {
        gene_color[V(graph)$name %in% GRN_node_list[['CTS_CP.1']]] = "black"
    }
    plot(graph, layout = lay, 
        vertex.size = 5,
        vertex.label.color = gene_color,
        )
}
dev.off()

