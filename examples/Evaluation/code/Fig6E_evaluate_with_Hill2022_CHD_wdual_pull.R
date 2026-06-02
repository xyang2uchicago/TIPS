updir = "F:/projects/scRNA/results/cardiac_CTS_GRN/validation_CHD_C8vsC16_v3/"
setwd(updir)
 

library(igraph)
library(dplyr)
library(igraph)
library(openxlsx)
library(gplots)

library(sandwich)
library(lmtest)

library(metafor)
library(ggplot2)

library(SummarizedExperiment)
library(metafor)

pull_path = 'F:/projects/TIPS/doc/MSB2026_v2/'
pull_df = read.xlsx(paste0(pull_path, 'STable_final_2026May_edited.xlsx'), sheet=14)  #  
dim(pull_df)  # 151  16
head(pull_df, 3)
#   Database    CTS_ID HiG_method    PPI linkeage     x     y        w_CP   w_decendent       delta  abs_delta direction
# 1   Ibarra cardiac.a     1v1.25 STRING       CM GATA6  TBX5  0.14273740  0.0538731312 -0.08886427 0.08886427  decrease
# 2   Ibarra cardiac.a     1v1.25 STRING       CM  CBFB GATA6  0.03651142 -0.0056903109 -0.04220173 0.04220173  decrease
# 3   Ibarra cardiac.a     1v1.25 STRING       CM  CBFB  TBX5 -0.01286985 -0.0002101566  0.01265969 0.01265969  increase
#    status rank TF_highConf.(fr.cisTarget.or.ChIPseq)      motif  NES
# 1 changed    1                                  <NA>       <NA>   NA
# 2 changed    2             CBFB (directAnnotation).  hdpi__CBFB 4.67
# 3 changed    3             CBFB (directAnnotation).  hdpi__CBFB 4.67
colnames(pull_df)[grepl('TF_highConf.', colnames(pull_df))] = 'TF_highConf'

unique(pull_df$Database)
#"Ibarra"            "Pijuan"            "Elorbany"          

unique(pull_df$PPI)
# [1] "STRING" "IID"  


pal_Cluster <- c(
  CF1 = "#D55E00",
  CF2 = "#E69F00",
  CF3 = "#F4A582",
  CF4 = "#FDB863",
  CM1 = "#0072B2",
  CM2 = "#56B4E9",
  CM3 = "#74ADD1",
  EpiC = "#1B9E77",
  EpiL = "#66C2A5",
  PeriC = "#7570B3",
  SMC = "#984EA3",
  Other = "#999999",
  CF = "darkgreen",
  CM = "blue"
)
cardiac_clusters = c("CF1", "CF2", "CF3", "CF4", "CM1", "CM2", "CM3")


################################################################
# load Hill2022 gene expression matrix whith all genes (not only HVG)
# refer to F:\projects\scRNA\source\GSE203274_CHD\1_build_object_AllNuclei.R
## this is the independent iPSVC time-course dataset, a Seurat object , for validation !!!!
library(Seurat)
db_specifc_input_path = 'F:/projects/scRNA/results/GSE203274_CHD/'
sce = readRDS(paste0(db_specifc_input_path,"Hill_AllNuclei_GSE203274_CHD_object.rds"))
sce
# An object of class Seurat 
# 29266 features across 157293 samples within 1 assay 
# Active assay: RNA (29266 features, 2000 variable features)
#  4 layers present: counts, Hill_norm, data, scale.data
#  3 dimensional reductions calculated: pca, harmony, umap

colnames(sce@meta.data)
#   [1] "orig.ident"                                  
#  [2] "nCount_RNA"                                  
#  [3] "nFeature_RNA"                                
#  [4] "NAME"                                        
#  [5] "orig_ident"                                  
#  [6] "labID"                                       
#  [7] "procedure"                                   
#  [8] "age"                                         
#  [9] "gender"                                      
# [10] "echoEF"                                      
# [11] "vers10X"                                     
# [12] "diagnosis"                                   
# [13] "ageCont"                                     
# [14] "region"                                      
# [15] "batch_indices"                               
# [16] "percent_mt"                                  
# [17] "ClinicalRank"                                
# [18] "DEid"                                        
# [19] "MainCellType"                                
# [20] "Cluster"                                     
# [21] "labID2"                                      
# [22] "colors"                                      
# [23] "Diagnosis"                                   
# [24] "patientID"                                   
# [25] "cell_id"                                     
# [26] "biosample_id"                                
# [27] "donor_id"                                    
# [28] "species"                                     
# [29] "species__ontology_label"                     
# [30] "disease"                                     
# [31] "disease__ontology_label"                     
# [32] "organ"                                       
# [33] "organ__ontology_label"                       
# [34] "library_preparation_protocol"                
# [35] "sex"                                         
# [36] "library_preparation_protocol__ontology_label"
# [37] "RNA_snn_res.1.5"                             
# [38] "seurat_clusters"    


meta = read.delim(paste0(db_specifc_input_path, 'Hill_AllNuclei_metadata_CMCFsubclusters_combined.txt'), sep = "\t", header = TRUE, row.names=1)
# refer to 25.1_evaluate_with_Hill2022_CHD_wdual_pull.R
# meta = readRDS('meta_data_Hill2022_CHD_wiPSC_pullSore.rds')
dim(meta) #[1] 157,2943     9
meta[1:3,1:4]
#                         orig.ident nCount_RNA nFeature_RNA
# P8_1_AAGCCATGTCGGCTAC-1          P8      24524         7957
# P40_2_GTTGTGAAGCGCCATC-1        P40       8037         3749
# P26_1_AGTGCCGAGATAACGT-1        P26       5459         3252
#                                              NAME
# P8_1_AAGCCATGTCGGCTAC-1   P8_1_AAGCCATGTCGGCTAC-1
# P40_2_GTTGTGAAGCGCCATC-1 P40_2_GTTGTGAAGCGCCATC-1
# P26_1_AGTGCCGAGATAACGT-1 P26_1_AGTGCCGAGATAACGT-1
length(meta$labID %>% unique)  #[1] 13
length(meta$patientID %>% unique) #[1] 14
length(meta$labID2 %>% unique)#[1] 14
setdiff(meta$labID2, meta$labID)  #[1] "13_198_LV" "13_198_RV"
all(meta$patientID == meta$labID2)  # TRUE
length(unique(meta$biosample_id)) # 30
table(sce$MainCellType)
#   Adipo      CF      CM    Endo   ENDOC    EpiC    EpiL     LEC     Mac    Mast
#     454   21034   73296   35673    1436      41      51     658    7010      92
# Neurons   PeriC     SMC  Tcells
#    1985   12037    2405    1121

if('TYPE' %in% colnames(meta)) meta <- meta[meta$NAME != "TYPE", ]
all(rownames(meta) == rownames(sce@meta.data)) # TRUE
all(meta$cell == sce@meta.data$cell)  # TRUE
all(meta$Cluster == sce@meta.data$Cluster)  # FALSE
table(meta$Cluster, sce@meta.data$Cluster)
#           Adipo    CF    CM  Endo ENDOC  EpiC  EpiL   LEC   Mac  Mast Neurons PeriC   SMC Tcells
#   Adipo     454     0     0     0     0     0     0     0     0     0       0     0     0      0
#   CF1         0  8279     0     0     0     0     0     0     0     0       0     0     0      0
#   CF2         0  5191     0     0     0     0     0     0     0     0       0     0     0      0
#   CF3         0  4375     0     0     0     0     0     0     0     0       0     0     0      0
#   CF4         0  3189     0     0     0     0     0     0     0     0       0     0     0      0
#   CM          0     0   143     0     0     0     0     0     0     0       0     0     0      0
#   CM1         0     0 31139     0     0     0     0     0     0     0       0     0     0      0
#   CM2         0     0  6131     0     0     0     0     0     0     0       0     0     0      0
#   CM3         0     0 35883     0     0     0     0     0     0     0       0     0     0      0
#   Endo        0     0     0 35673     0     0     0     0     0     0       0     0     0      0
#   ENDOC       0     0     0     0  1436     0     0     0     0     0       0     0     0      0
#   EpiC        0     0     0     0     0    41     0     0     0     0       0     0     0      0
#   EpiL        0     0     0     0     0     0    51     0     0     0       0     0     0      0
#   LEC         0     0     0     0     0     0     0   658     0     0       0     0     0      0
#   Mac         0     0     0     0     0     0     0     0  7010     0       0     0     0      0
#   Mast        0     0     0     0     0     0     0     0     0    92       0     0     0      0
#   Neurons     0     0     0     0     0     0     0     0     0     0    1985     0     0      0
#   PeriC       0     0     0     0     0     0     0     0     0     0       0 12037     0      0
#   SMC         0     0     0     0     0     0     0     0     0     0       0     0  2405      0
#   Tcells      0     0     0     0     0     0     0     0     0     0       0     0     0   1121
sce@meta.data[,'Cluster'] = meta$Cluster  # TRUE
setdiff(colnames(meta), colnames(sce@meta.data)) #[1] "SubClusters"
sce <- AddMetaData(sce, metadata = meta[, "SubClusters", drop = FALSE])

library(SingleCellExperiment)
sce <- as.SingleCellExperiment(sce)
sce
# class: SingleCellExperiment 
# dim: 29266 157293 
# metadata(0):
# assays(2): counts logcounts
# rownames(29266): RP11-34P13.7 AL627309.1 ... RP11-352D3.2 BABAM1
# rowData names(0):
# colnames(157293): P8_1_AAGCCATGTCGGCTAC-1 P40_2_GTTGTGAAGCGCCATC-1 ... RV_198_1_GTCGTAACACAAGCCC-1
#   RV_198_1_ATCGTGAAGGTCATAA-1
# colData names(40): orig.ident nCount_RNA ... SubClusters ident
# reducedDimNames(3): PCA HARMONY UMAP
# mainExpName: RNA
# altExpNames(0):
range(logcounts(sce))  # [1] 0.00000 9.10805   this is a lognormalizaed and but NOT scaled counts !! 

# !!!!!!!!!  when the downloaded feature names don't have all CTS genes, manually check the aliase   !!!
# setdiff(CTS[[CTS_name]], rownames(sce))
#   "APELA" "GM266" "P3H2"    
# tmp = read.delim(paste0('F:/projects/scRNA/data/GSE203274_CHD/SCP1852/scalegenes.tsv.gz'), sep = "\t", header = F)
'C1orf168' %in% rownames(sce) #[1] TRUE
rownames(sce)[which(rownames(sce)=='C1orf168')] = 'FYB2'  #<<<<<<<<<<<<< Alas, this is an update !!!!
'PTPLB' %in% rownames(sce) #[1] TRUE
rownames(sce)[which(rownames(sce)=='PTPLB')] = 'HACD2'  #<<<<<<<<<<<<< Alas, this is an update !!!!
'LEPREL1' %in% rownames(sce) 
rownames(sce)[which(rownames(sce)=='LEPREL1')] = 'P3H2'  #<<<<<<<<<<<<< Alas, this is an update !!!!

# They appear to have both 'Diagnosis' and 'diagnosis'. Use the one that is populated.
DIAG_COL <- "Diagnosis"
CLUSTER_COL <- "Cluster"      # paper panel shows CF1-4; use "SubClusters" if you want CF1-6

CHD_levels = levels = c("Donor", "TOF", "DCM",  "HCM",  "HF_HLHS", "Neo_HLHS")
# Published legend looks like:
diag_cols <- c(
  "Donor"   = "#999999",  # grey
  "TOF"     = "#1F77B4",  # blue/teal
  "DCM"     = "#17BECF",  # cyan
  "HCM"     = "#F3E79B",  # pale yellow
  "HF_HLHS" = "#F4A582",  # peach
  "Neo_HLHS"= "#D73027"   # red
)

if(!"Diagnosis2" %in% colnames(colData(sce))) {
  colData(sce)$Diagnosis2 = colData(sce)$Diagnosis
  colData(sce)$Diagnosis2[which( colData(sce)$Diagnosis %in% c('HF_HLHS', 'Neo_HLHS'))] = 'HLHS'
}

################## load the pre-calcualted signatures ##############################################
meta  = readRDS('Hill2022_colData_sce_dualpull_colMeans.rds')
colnames(meta)

colData(sce) = meta 
################################################################

## becasue only CP.1 has both cf and cm pull #for(CTS_name in c( 'CP', 'CP.1')){
mat = assays(sce)[["logcounts"]]

## becasue only CP.1 has both cf and cm pull #for(CTS_name in c( 'CP', 'CP.1')){
table(pull_df$Database)
   #             Elorbany           Ibarra           Pijuan_Sala Pijuan_10celltype  Pijuan_9celltype 
   #                86                    49                    17                   260                    26 
sig_map <- data.frame(
  db = c( "Ibarra", "Ibarra", "Ibarra",  "Pijuan","Pijuan",   "Elorbany","Elorbany","Elorbany"),
  HiG_method = c("1v1.25", "1v1.25", "1v1.25",  "1v1.25", "1v1.25",    "1v1.25", "1v1.25","1v1.25"),
  PPI = c("STRING", "STRING", "IID",   "STRING", "IID",     "STRING", "STRING","IID"),
  CTS_name = c("cardiac.a", "cardiac.a","cardiac.a",   "8", "8",   "CP.1", "CP.1","CP.1"),
  cf_in_db = c("ExEM", "ExEM", "ExEM",   "C16", "C16",   "C1" , "C1","C1"),  #cm_col = c("CM_pullcardiac.a", "CM_pull8_wincreased_nodes", "CM_pullCP.1_wincreased_nodes"), ## when test C8vsC16
  cm_in_db = c("cardiac.c", "cardiac.c_", "cardiac.c",   "C17", "C17",   "C5", "C5","C5"), ## when test C8vsC18
  extend = c("", "noTBX5", "",   "", "",   "",  "wincreased","wincreased"),
  stringsAsFactors = FALSE
)

sig_map$ID = paste0(sig_map$db, '_', sig_map$HiG_method , '_', sig_map$PPI, '_', sig_map$CTS_name, sig_map$extend)
pull_df$ID = paste0(pull_df$Database, '_', '_', pull_df$PPI, '_', pull_df$CTS_ID)
print(sig_map)
        # db HiG_method    PPI  CTS_name cf_in_db   cm_in_db     extend                                    ID
# 1   Ibarra     1v1.25 STRING cardiac.a     ExEM  cardiac.c                   Ibarra_1v1.25_STRING_cardiac.a
# 2   Ibarra     1v1.25 STRING cardiac.a     ExEM cardiac.c_     noTBX5  Ibarra_1v1.25_STRING_cardiac.anoTBX5
# 3   Ibarra     1v1.25    IID cardiac.a     ExEM  cardiac.c                      Ibarra_1v1.25_IID_cardiac.a
# 4   Pijuan     1v1.25 STRING         8      C16        C17                           Pijuan_1v1.25_STRING_8
# 5   Pijuan     1v1.25    IID         8      C16        C17                              Pijuan_1v1.25_IID_8
# 6 Elorbany     1v1.25 STRING      CP.1       C1         C5                      Elorbany_1v1.25_STRING_CP.1
# 7 Elorbany     1v1.25 STRING      CP.1       C1         C5 wincreased Elorbany_1v1.25_STRING_CP.1wincreased
# 8 Elorbany     1v1.25    IID      CP.1       C1         C5 wincreased    Elorbany_1v1.25_IID_CP.1wincreased


setdiff(sig_map$db, pull_df$Database)  #character(0)
setdiff(sig_map$CTS_name, pull_df$CTS_ID) #character(0)


########################################################################################
################## first,  CF4 effect in Hill using the full regression model (Diagnosis + region), 
## as we did for individual signature, with the coefficient and CI for HLHS/Neo-HLHS vs Donor.
########################################################################################
source("E:/Git_Holly/TIPS/R/celltype_specific_weight_v10.R")
get_endpoint_effect  ## a customalized function to reuse 

# clean background to run base::qr()
if ("package:brainGraph" %in% search()) detach("package:brainGraph", unload = TRUE, character.only = TRUE)
if("brainGraph" %in% loadedNamespaces()) {
    unloadNamespace("brainGraph")
    .old_qr_array <- getS3method("qr", "array", optional = TRUE)

    registerS3method(
    "qr", "array",
    function(x, ...) {
        if (length(dim(x)) == 2L) {
        qr.default(unclass(x), ...)
        } else {
        apply(x, 3L, qr.default, ...)
        }
    }
    )
}




################################################
## note that get_cf4_effect() filters to Donor/Neo_HLHS only, while get_endpoint_effect() can keep other Diagnosis levels in the regression, shifting the coefficient.
## tehrefore this resutls similar to the full regression analysis, but with a different coefficient.
meta_tab_biosample <- bind_rows(lapply(seq_len(nrow(sig_map)), function(i) {
  id = sig_map$ID[i]
  out <- get_cf4_effect(sce = sce, 
    cf_col = paste0('CF_', id), cm_col = paste0('CM_', id),
    cluster = "CF4",  group_col = "Diagnosis",  case_level = "Neo_HLHS",
    control_level = "Donor", min_cells = 2, pseudobulk_level='biosample_id', 
    scale_outcome = TRUE, # => SD units, otherwise raw units
    adjust_region = TRUE  # if true fits imbalance_out ~ group + region, otherwise imbalance_out ~ group
  )
  if (!is.null(out)) out$signature <- id
  out
}))
## no reason to remove TBX5 for disease validation , neither for the intermediate contrl exeriment ##  update on 5/19/2026
x = which(meta_tab_biosample$signature %in% c('Ibarra_1v1.25_STRING_cardiac.anoTBX5', 'Elorbany_1v1.25_STRING_CP.1'))
if(length(x)>0) meta_tab_biosample = meta_tab_biosample[-x, ]

print(meta_tab_biosample)
est  <- as.numeric(meta_fit_biosample$b[1])
p <- meta_fit_biosample$pval
meta_fit_biosample <- metafor::rma.uni(yi = yi, sei = sei, data = meta_tab_biosample, method = "REML")
print(meta_fit_biosample)

meta_tab_biosample$lo <- meta_tab_biosample$yi - 1.96 * meta_tab_biosample$sei
meta_tab_biosample$hi <- meta_tab_biosample$yi + 1.96 * meta_tab_biosample$sei

forest(meta_fit_biosample, slab = paste0(meta_tab_biosample$signature, " (", meta_tab_biosample$n_case_biosamples, " vs ", meta_tab_biosample$n_ctrl_biosamples, ")"),
       xlab = "CF4 effect:  Neo_HLHS vs Donor (SD units)", main = paste0("Biosample-level summary, est=", round(est, 1), ", p=", signif(p, 2)))
dev.copy(pdf, file = 'CF4_effect_Biosample_Neo_HLHS.pdf')
dev.off()

saveRDS(meta_tab_biosample, file = 'meta_tab_biosample_Neo_HLHS.rds')

#########################################################################################
################## third ,  Cross-signature robustness summary, 
## id, Sensitivity of the CF4 Neo-HLHS endpoint to alternative dual-pull signature definitions.
########################################################################################

meta_tab_biosample <- bind_rows(lapply(seq_len(nrow(sig_map)), function(i) {
    id = sig_map$ID[i]
    out <- get_cf4_effect(sce = sce, cf_col = paste0('CF_', id), cm_col = paste0('CM_', id),
        cluster = "CF4",  group_col = "Diagnosis",  case_level = "HF_HLHS",
        control_level = "Donor", min_cells = 3, pseudobulk_level='biosample_id', 
        scale_outcome = TRUE, # => SD units, otherwise raw units
        adjust_region = TRUE  # if true fits imbalance_out ~ group + region, otherwise imbalance_out ~ group
    )
    if (!is.null(out)) out$signature <- id
    out
}))
## no reason to remove TBX5 for disease validation , neither for the intermediate contrl exeriment ##  update on 5/19/2026
x = which(meta_tab_biosample$signature %in% c('Ibarra_1v1.25_STRING_cardiac.anoTBX5', 'Elorbany_1v1.25_STRING_CP.1'))
if(length(x)>0) meta_tab_biosample = meta_tab_biosample[-x, ]


meta_tab_biosample
meta_fit_biosample <- metafor::rma.uni(yi = yi, sei = sei, data = meta_tab_biosample, method = "REML")
print(meta_fit_biosample)
est  <- as.numeric(meta_fit_biosample$b[1])
p <- meta_fit_biosample$pval
meta_tab_biosample$lo <- meta_tab_biosample$yi - 1.96 * meta_tab_biosample$sei
meta_tab_biosample$hi <- meta_tab_biosample$yi + 1.96 * meta_tab_biosample$sei

forest(meta_fit_biosample, slab = paste0(meta_tab_biosample$signature, " (", meta_tab_biosample$n_case_biosamples, " vs ", meta_tab_biosample$n_ctrl_biosamples, ")"),
       xlab = "CF4 effect:  HF_HLHS vs Donor (SD units)", main = paste0("Biosample-level summary, est=", round(est, 1), ", p=", signif(p, 2)))
dev.copy(pdf, file = 'CF4_effect_Biosample_HF_HLHS.pdf')
dev.off()

saveRDS(meta_tab_biosample, file = 'meta_tab_biosample_HF_HLHS.rds')  ## so good except the 10celltype_1vn !!!


meta_tab_biosample <- bind_rows(lapply(seq_len(nrow(sig_map)), function(i) {
    id = sig_map$ID[i]
	out <- get_cf4_effect(sce = sce, cf_col = paste0('CF_', id), cm_col = paste0('CM_', id),
    cluster = "CF4",  group_col = "Diagnosis2",  case_level = "HLHS",
    control_level = "Donor", min_cells = 3, pseudobulk_level='biosample_id', 
	scale_outcome = TRUE, # => SD units, otherwise raw units
    adjust_region = TRUE # if true fits imbalance_out ~ group + region, otherwise imbalance_out ~ group
  )
  if (!is.null(out)) out$signature <- id
  out
}))
## no reason to remove TBX5 for disease validation , neither for the intermediate contrl exeriment ##  update on 5/19/2026
x = which(meta_tab_biosample$signature %in% c('Ibarra_1v1.25_STRING_cardiac.anoTBX5', 'Elorbany_1v1.25_STRING_CP.1'))
if(length(x)>0) meta_tab_biosample = meta_tab_biosample[-x, ]


# meta_tab_biosample
# meta_fit_biosample <- metafor::rma.uni(yi = yi, sei = sei, data = meta_tab_biosample, method = "REML")
# print(meta_fit_biosample)
# est  <- as.numeric(meta_fit_biosample$b[1])
# pval <- as.numeric(meta_fit_biosample$pval[1]) 
# meta_tab_biosample$lo <- meta_tab_biosample$yi - 1.96 * meta_tab_biosample$sei
# meta_tab_biosample$hi <- meta_tab_biosample$yi + 1.96 * meta_tab_biosample$sei

# forest(meta_fit_biosample, slab = paste0(meta_tab_biosample$signature, " (", meta_tab_biosample$n_case_biosamples, " vs ", meta_tab_biosample$n_ctrl_biosamples, ")"),
       # xlab = "CF4 effect:  jointHLHS vs Donor (SD units)", 
	   # main = paste0("Biosample-level summary, est=", round(est, 1), ", p=",  format.pval(pval, digits = 2, eps = 1e-4)))
# dev.copy(pdf, file = 'CF4_effect_Biosample_jointHLHS.pdf')  
# dev.off()

# saveRDS(meta_tab_biosample, file = 'meta_tab_biosample_jointHLHS.rds')  # so good except the 10celltype_1vn !!!


 