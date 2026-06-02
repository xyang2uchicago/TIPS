
updir = "F:/projects/scRNA/results/cardiac_CTS_GRN/validation_TBX5_C8vsC16_v2/KO_SET2_SC_DF_processed/"  
# create the directory if not exists
if (!dir.exists(updir))  dir.create(updir, recursive = TRUE)
setwd(updir)


## refer to codes F:\projects\scRNA\source\GSE190475_TBX5_enhKO\
# 2_SC_process_fr_raw.py  # reanalyze from raw counts; QC; normalzation, PCA, leiden clustering
# 2.2_Large_screen_downstream.py; kepp 14 clusters with > 200 cells, kepp 11 of 14 clusters have TBX5KO or NC control ceells. 
# 3_h5ad_to_rds.R
 
files = list.files(path = 'F:/projects/scRNA/results/GSE190475_TBX5_enhKO/', pattern = '*.rds', full.names = TRUE)
files

(files = files[which(grepl('SET2_', files))])
# [1] "F:/projects/scRNA/results/GSE190475_TBX5_enhKO/GSE190475_KO_SET2_CM_DF.rds"
# [2] "F:/projects/scRNA/results/GSE190475_TBX5_enhKO/KO_SET2_SC_DF_processed.rds"
 

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
library('scater')
library(RColorBrewer)
library(Matrix)

library(grid)
library(ggrastr)
library("rlang")

pull_path = 'F:/projects/TIPS/doc/'
pull_df = read.xlsx(paste0(pull_path, 'STable_final_2026.xlsx'), sheet=16)  # for C8vsC18, extended for IID, Pijuan E6.5-E8.5, read this 
dim(pull_df)  # 437  17
head(pull_df, 3)
#   Database    CTS_ID HiG_method    PPI linkeage     x     y        w_CP   w_decendent       delta  abs_delta direction  status rank
# 1   Ibarra cardiac.a     1v1.25 STRING       CM GATA6  TBX5  0.14273740  0.0538731312 -0.08886427 0.08886427  decrease changed    1
# 2   Ibarra cardiac.a     1v1.25 STRING       CM  CBFB GATA6  0.03651142 -0.0056903109 -0.04220173 0.04220173  decrease changed    2
# 3   Ibarra cardiac.a     1v1.25 STRING       CM  CBFB  TBX5 -0.01286985 -0.0002101566  0.01265969 0.01265969  increase changed    3
#   TF_highConf.(fr.cisTarget.or.ChIPseq)      motif  NES
# 1                                  <NA>       <NA>   NA
# 2             CBFB (directAnnotation).  hdpi__CBFB 4.67
# 3             CBFB (directAnnotation).  hdpi__CBFB 4.67
  
# cols <- setNames(
#   c("orange", "blue", "lightgreen"),
#   c(CP_cluster, CM_cluster, CF_cluster)
# )
colnames(pull_df)[grepl('TF_highConf.', colnames(pull_df))] = 'TF_highConf'

unique(pull_df$Database)
#"Ibarra"            "Pijuan"            "Elorbany"          "Pijuan_9celltype"  "Pijuan_10celltype"
unique(pull_df$HiG_method)
# [1] "1v1.25" "1vn"  
unique(pull_df$PPI)
# [1] "STRING" "IID"  

##########################################################################################################################################
 
 
f = basename(files[2])
ID = gsub(".rds", "", f) %>% gsub("GSE190475_", "", .)
ID  # KO_SET2_SC_DF_processed

sample_levels = levels = c('WT', 'EXON', 'ENH3', 'ENH5')
 

#db_specifc_input_path1 = 'F:/projects/scRNA/data/GSE190475_TBX5_enhKO/'  # directly downloaded , for the objects of Figs5/6
db_specifc_input_path = 'F:/projects/scRNA/results/GSE190475_TBX5_enhKO/'  # processed to add the sample metadata, for the object of Fig1/3 

sce = readRDS(paste0(db_specifc_input_path, f))
sce
# class: SingleCellExperiment 
# dim: 26099 30913 
# metadata(17): Sample_colors hvg ... scrublet umap
# assays(4): X counts logcounts scaled
# rownames(26099): AL627309.1 AL627309.5 ... AC007325.4 AC007325.2
# rowData names(21): gene_ids feature_types ... dispersions_norm varm
# colnames(30913): AAACCCAAGAATAACC-1 AAACCCAAGAATACAC-13 ... TTTGTTGTCTGAGCAT-1
#   TTTGTTGTCTGTCAGA-7
# colData names(26): T_Reps B_Reps ... leiden_0.8 leiden_1.0
# reducedDimNames(4): X_pca X_pca_published X_umap X_umap_published
# mainExpName: NULL
# altExpNames(0):

rownames(sce)[which(rownames(sce) == 'C1orf168')] = 'FYB2'  # updating according to genecards.org

## check mitrocondro interaction
colnames(rowData(sce))
# [1] "gene_ids"              "feature_types"         "genome"               
#  [4] "n_counts"              "mean"                  "std"                  
#  [7] "mt"                    "ribo"                  "hb"                   
# [10] "n_cells_by_counts"     "mean_counts"           "log1p_mean_counts"    
# [13] "pct_dropout_by_counts" "total_counts"          "log1p_total_counts"   
# [16] "n_cells"               "highly_variable"       "means"                
# [19] "dispersions"           "dispersions_norm"      "varm"        

summary(rowData(sce)$n_counts)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#       3      36     635   11347    5518 5212929 

table(rowData(sce)$feature_types)
    #  Gene Expression  
    #            26099                      

## to verify mito%, we need raw counts
is.mito <- grepl("^MT-", rownames(sce))
table(is.mito)
# is.mito
# FALSE  TRUE 
# 26086    13

names(assays(sce))  #[1]   "counts"    "logcounts" "scaled"      
assays(sce)[["X"]][1:10] ## should be scaled data
#  [1] -0.012229370 -0.016430018 -0.007015169 -0.009826314 -0.046628002
#  [6] -0.007015169 -0.079207405 -0.011371258 -0.015601349 -0.044985756

colnames(colData(sce))
#  [1] "T_Reps"                  "B_Reps"                  "Sample"                  "n_counts_all"           
#  [5] "louvain"                 "n_genes_by_counts"       "log1p_n_genes_by_counts" "total_counts"           
#  [9] "log1p_total_counts"      "total_counts_mt"         "log1p_total_counts_mt"   "pct_counts_mt"          
# [13] "total_counts_ribo"       "log1p_total_counts_ribo" "pct_counts_ribo"         "total_counts_hb"        
# [17] "log1p_total_counts_hb"   "pct_counts_hb"           "n_genes"                 "T_batch"                
# [21] "doublet_score"           "predicted_doublet"       "leiden_0.2"              "leiden_0.5"             
# [25] "leiden_0.8"              "leiden_1.0"   
  
table(colData(sce)$louvain, colData(sce)$Sample)
#   ENH3 ENH5 EXON   WT
#   0 2715 4393 3179 3794
#   1  823  380 1522 5066
#   2   40  360   88 6084
#   3  745  359  419  421
#   4   73  171   53  228


# Clean up T_Reps field: remove Sample prefix in a vectorized way
# Prevents gsub(pattern, replacement, x) warning when 'pattern' is vectorized
# Use mapply to apply each sample as the prefix per element, or directly sub() per element:
sce$T_Reps <- mapply(function(trep, sample) sub(paste0('^', sample, '-'), '', trep), sce$T_Reps, sce$Sample)
table(sce$T_Reps)
# T1_B1 T1_B2 T2_B1 T2_B2 
#  6805  4401  4086  4727 

# Remove the sample prefix from B_Reps to get clean B_Rep labels
sce$B_Reps <- mapply(function(brep, sample) sub(paste0('^', sample, ':'), '', brep), sce$B_Reps, sce$Sample)
table(sce$B_Reps)
#  B1   B2 
#  133 4738 

## feature plot for the celltypes === somthing is inconsisten between the python h5ad and the rds object!! 
p1 = plotReducedDim(sce, dimred = "X_umap_published", colour_by = "louvain", point_size = 0.3)
p2 = plotReducedDim(sce, dimred = "X_umap", colour_by = "louvain", point_size = 0.3)
gridExtra::grid.arrange(p1, p2, ncol = 2)
# plotReducedDim(
#   sce,
#   "X_umap",
#   colour_by = I(sce$leiden_0.8 == 1),
#   point_size = 0.4
# )


reducedDim(sce, "UMAP") <- reducedDim(sce, "X_umap_published")


# 1) plotting, refer to 2evaluation_TBX5_v1/6.1_evaluate_with_Fig5_TBX5HOMOKOvsCTL_SET2_C8vsC18_clean.R 
################################################################


# ################################  2) scoring #############################################################
## Within Cluster CF4, does my curated CF/ECM gene signature score higher in HLHS than donors (i.e., program activation within the cluster)?
## imbalance=z(CM score)−z(CF score)
mat = assays(sce)[["logcounts"]]

## becasue only CP.1 has both cf and cm pull #for(CTS_name in c( 'CP', 'CP.1')){
table(pull_df$Database)
   #             Elorbany           Ibarra           Pijuan_Sala Pijuan_10celltype  Pijuan_9celltype 
   #                86                    49                    17                   260                    26 
sig_map <- data.frame(
  db = c( "Ibarra", "Ibarra", "Ibarra", "Pijuan","Pijuan",  "Elorbany","Elorbany"),
  HiG_method = c("1v1.25", "1v1.25", "1v1.25","1v1.25", "1v1.25",   "1v1.25","1v1.25"),
  PPI = c("STRING", "STRING", "IID", "STRING", "IID",    "STRING","IID"),
  CTS_name = c("cardiac.a", "cardiac.a","cardiac.a", "8", "8",   "CP.1","CP.1"),
  cf_in_db = c("ExEM", "ExEM", "ExEM", "C18", "C18",  "C1","C1"),  #cm_col = c("CM_pullcardiac.a", "CM_pull8_wincreased_nodes", "CM_pullCP.1_wincreased_nodes"), ## when test C8vsC16
  cm_in_db = c("cardiac.c", "cardiac.c_", "cardiac.c", "C17", "C17",    "C5","C5"), ## when test C8vsC18
  extend = c("", "noTBX5", "", "", "",     "wincreased","wincreased"),
  stringsAsFactors = FALSE
)

sig_map$ID = paste0(sig_map$db, '_', sig_map$HiG_method, '_', sig_map$PPI, '_', sig_map$CTS_name, sig_map$extend)
pull_df$ID = paste0(pull_df$Database, '_', pull_df$HiG_method, '_', pull_df$PPI, '_', pull_df$CTS_ID)
print(sig_map)
#                   db HiG_method    PPI  CTS_name cf_in_db   cm_in_db     extend                                    ID
# 1             Ibarra     1v1.25 STRING cardiac.a     ExEM  cardiac.c                   Ibarra_1v1.25_STRING_cardiac.a
# 2             Ibarra     1v1.25 STRING cardiac.a     ExEM cardiac.c_     noTBX5  Ibarra_1v1.25_STRING_cardiac.anoTBX5
# 3             Ibarra     1v1.25    IID cardiac.a     ExEM  cardiac.c                      Ibarra_1v1.25_IID_cardiac.a
# 4             Pijuan     1v1.25 STRING         8      C18        C17                           Pijuan_1v1.25_STRING_8
# 5             Pijuan     1v1.25    IID         8      C18        C17                              Pijuan_1v1.25_IID_8
# 6          Elorbany     1v1.25 STRING      CP.1       C1         C5 wincreased Elorbany_1v1.25_STRING_CP.1wincreased
# 7          Elorbany     1v1.25    IID      CP.1       C1         C5 wincreased    Elorbany_1v1.25_IID_CP.1wincreased


#############    testing on the SC dataset using author's louvain and Holly's leiden_0.8clustering results ###################################################

recalcualte = FALSE 
if(recalcualte) {
  
  for(i in 1:nrow(sig_map)){
	id = sig_map$ID[i]
	dbx = sig_map$db[i]
	CTS_namex = sig_map$CTS_name[i]
    HiG_methodx = sig_map$HiG_method[i]
	PPIx = sig_map$PPI[i]
	cf_pull_name = paste0('CF_', id)			#!!!ADDED
	cm_pull_name = paste0('CM_', id)			
	dual_pull_name = paste0('dual_', id)	
    
	cf_pull_genes_0 = cf_pull_genes = cm_pull_genes_0 = cm_pull_genes = dual_pull_genes = NULL
	sub_df = subset(pull_df, Database==dbx & HiG_method==HiG_methodx & PPI==PPIx  & CTS_ID==CTS_namex)
	#table(sub_df$TF_highConf == "ISL1 ChIP-seq", sub_df$linkeage)
   # prediction = ifelse(sub_df$TF_highConf == "ISL1 ChIP-seq", 'ChIP', 'cisTarget')
	
	# for(j in unique(prediction)) {
	# 	ID2 = paste0(ID, '_', prediction)
	# Step 1 — Compute per-cell CF/ECM signature score 
	cf_pull_genes_0 <- subset(sub_df, linkeage=='CF' )[, c('x','y')] %>%
					unlist() %>% unique()
	cm_pull_genes_0 <- subset(sub_df, linkeage=='CM' )[, c('x','y')] %>%
				unlist() %>% unique()
	(dual_pull_genes = intersect(cf_pull_genes_0, cm_pull_genes_0 )	)	
	if(length(dual_pull_genes) > 0) {
		cf_pull_genes = setdiff(cf_pull_genes_0, dual_pull_genes) 
		cm_pull_genes = setdiff(cm_pull_genes_0, dual_pull_genes)
	} else {
		cf_pull_genes = cf_pull_genes_0
		cm_pull_genes = cm_pull_genes_0
	}

    if(length(cf_pull_genes) > 0)  {
		colData(sce)[[cf_pull_name]] <- Matrix::colMeans(mat[cf_pull_genes, , drop = FALSE])
		cat(paste0('Calculating CF score for ', id, ': ', paste(cf_pull_genes, collapse = ', '), '\n'))
	} else {
		cat(paste0('No CF score for ', id), '\n')
		colData(sce)[[cf_pull_name]] <- NA
	}
	if(length(cm_pull_genes) > 0) {
		colData(sce)[[cm_pull_name]] <- Matrix::colMeans(mat[cm_pull_genes, , drop = FALSE])	
		cat(paste0('Calculating CM score for ', 	id, ': ', paste(cm_pull_genes, collapse = ', '), '\n'))
	} else {
		cat(paste0('No CM score for ', id), '\n')
		colData(sce)[[cm_pull_name]] <- NA
	}

    ## for the hiPSC dataset, add the score of CM_pull that adding the increased nodes (FGF10, LRRTM1) from dual-pull set
    if(id %in% c( 'Elorbany_1v1.25_STRING_CP.1wincreased', 'Elorbany_1v1.25_IID_CP.1wincreased')) {
          cm_pull_genes = c(cm_pull_genes, 'FGF10', 'LRRTM1')
          cat(paste0('Adding increased nodes for ', id, ': ', paste(cm_pull_genes, collapse = ', '), '\n'))
          colData(sce)[[cm_pull_name]] <- Matrix::colMeans(mat[cm_pull_genes, , drop = FALSE])	
        #  colData(sce)[[paste0(cm_pull_name,'_wincreased_nodes')]] <- Matrix::colMeans(mat[cm_pull_genes, , drop = FALSE])	
        }  
    ## for the Ibarra dataset, add the score without TBX5 for evaluation only
    if(id == 'Ibarra_1v1.25_STRING_cardiac.anoTBX5')  {
		cat(paste0('Removing TBX5 for ', id, ': ', paste(cm_pull_genes, collapse = ', '), '\n'))
        cm_pull_genes = setdiff(cm_pull_genes, 'TBX5')
		colData(sce)[[cm_pull_name]] <- Matrix::colMeans(mat[cm_pull_genes, , drop = FALSE])	
      #  colData(sce)[[paste0(cm_pull_name,'_noTBX5')]] <- Matrix::colMeans(mat[cm_pull_genes, , drop = FALSE])	
    }  
	 
	saveRDS(colData(sce), file = paste0(id, '_colData_sce_dualpull_colMeans.rds'))
	}
} else { ## END DO NOT recalculate !!!
		meta  = readRDS(paste0(id, '_colData_sce_dualpull_colMeans.rds'))
		colData(sce) = meta
}
# Adding CF score for Ibarra_1v1.25_STRING_cardiac.a: ANXA2, TWIST1
# Adding CM score for Ibarra_1v1.25_STRING_cardiac.a: GATA6, CBFB, ETV2, GATA4, MYOCD, PRDM6, RBFOX2, TBX5, SFRP5, RBM24
# Adding CF score for Ibarra_1v1.25_STRING_cardiac.anoTBX5: ANXA2, TWIST1
# Adding CM score for Ibarra_1v1.25_STRING_cardiac.anoTBX5: GATA6, CBFB, ETV2, GATA4, MYOCD, PRDM6, RBFOX2, TBX5, SFRP5, RBM24
# Adding CF score for Ibarra_1v1.25_IID_cardiac.a: ANXA2, ISL1, TWIST1
# Adding CM score for Ibarra_1v1.25_IID_cardiac.a: CBFB, GATA6, ETV2, GATA4, PRDM6, RBFOX2, TBX5, SFRP5, MYOCD, RBM24
# Adding CF score for Pijuan_1v1.25_STRING_8: DAB2, TPM27
# Adding CM score for Pijuan_1v1.25_STRING_8: FGF8, MPPED2, BMP7, RARB, NPM3
# Adding CF score for Pijuan_1v1.25_IID_8: DAB2 
# Adding CM score for Pijuan_1v1.25_IID_8:  FGF8, MPPED2, RARB
# Adding CF score for Elorbany_1v1.25_STRING_CP.1wincreased: FAM89A, PPFIBP1, SEMA3E
# Adding CM score for Elorbany_1v1.25_STRING_CP.1wincreased: H1F0, HAS2
# Adding increased nodes for Elorbany_1v1.25_STRING_CP.1wincreased: H1F0, HAS2, FGF10, LRRTM1
# Calculating CF score for Elorbany_1v1.25_IID_CP.1wincreased: FAM89A, PPFIBP1, SEMA3E
# Calculating CM score for Elorbany_1v1.25_IID_CP.1wincreased: HAS2
# Adding increased nodes for Elorbany_1v1.25_IID_CP.1wincreased: HAS2, FGF10, LRRTM1

colnames(colData(sce))
#  [1] "T_Reps"                                   "B_Reps"                                  
#  [3] "Sample"                                   "n_counts_all"                            
#  [5] "louvain"                                  "n_genes_by_counts"                       
#  [7] "log1p_n_genes_by_counts"                  "total_counts"                            
#  [9] "log1p_total_counts"                       "total_counts_mt"                         
# [11] "log1p_total_counts_mt"                    "pct_counts_mt"                           
# [13] "total_counts_ribo"                        "log1p_total_counts_ribo"                 
# [15] "pct_counts_ribo"                          "total_counts_hb"                         
# [17] "log1p_total_counts_hb"                    "pct_counts_hb"                           
# [19] "n_genes"                                  "T_batch"                                 
# [21] "doublet_score"                            "predicted_doublet"                       
# [23] "leiden_0.2"                               "leiden_0.5"                              
# [25] "leiden_0.8"                               "leiden_1.0"                              
# [27] "CF_Ibarra_1v1.25_STRING_cardiac.a"        "CM_Ibarra_1v1.25_STRING_cardiac.a"       
# [29] "CF_Ibarra_1v1.25_STRING_cardiac.anoTBX5"  "CM_Ibarra_1v1.25_STRING_cardiac.anoTBX5" 
# [31] "CF_Ibarra_1v1.25_IID_cardiac.a"           "CM_Ibarra_1v1.25_IID_cardiac.a"          
# [33] "CF_Pijuan_1v1.25_STRING_8"                "CM_Pijuan_1v1.25_STRING_8"               
# [35] "CF_Pijuan_1v1.25_IID_8"                   "CM_Pijuan_1v1.25_IID_8"                  
# [37] "CF_Elorbany_1v1.25_STRING_CP.1wincreased" "CM_Elorbany_1v1.25_STRING_CP.1wincreased"
# [39] "CF_Elorbany_1v1.25_IID_CP.1wincreased"    "CM_Elorbany_1v1.25_IID_CP.1wincreased"   

plot(colData(sce)[,"CM_Elorbany_1v1.25_STRING_CP.1wincreased"], colData(sce)[,"CM_Elorbany_1v1.25_IID_CP.1wincreased"])
############### 3) forest plor same CM or CF-like clusters endpoint across 3 datasets ##############


########################################################################################
source("E:/Git_Holly/TIPS/R/celltype_specific_weight_v10.R")
# get_endpoint_effect  ## a customalized function to reuse, which keeps other Diagnosis levels in the regression,  
## in contrast, get_cf4_effect() filters to Donor/Neo_HLHS only -- USED 

# get_cf4_effect ## a customalzied functiopn


## focus on the TBX5KO or NC control cells to do regression across library_names
table(sce$Sample) %>% sum() == ncol(sce) # [1] TRUE
unique(sce$Sample)  
# [1] WT   ENH3 EXON ENH5
# Levels: ENH3 ENH5 EXON WT

sce$Sample = factor(sce$Sample, levels = sample_levels)  ## for plotting
   
table(sce$Sample, sce$B_Reps)
#         R1   R2
#   WT   8862 6731
#   EXON 2794 2467
#   ENH3 1988 2408
#   ENH5 2826 2837




########################################################
### for simplicity, we reuse the function get_cf4_effect() to plot the forest plot
## for Author's louvin  
########################################################
# add dummy column:
sce$zero <- 0
sce$Cluster = sce$louvain  #!!!! used by get_cf4_effect()


# build meta_df:
# score_type	interpretation
# dual:CM-CF
# CM:	pure CM program
# CF:	pure CF program

scale_outcome = FALSE
df_plot = NULL
for(cl in c('0','1' ,'2', '3','4')){
	for(t in setdiff(unique(sce$Sample), "WT")){
 		meta_df <- bind_rows(lapply(1:nrow(sig_map), function(i) {

        id <- sig_map$ID[i]

		res_cm <- get_cf4_effect(
			sce = sce,
			cf_col = "zero",
			cm_col = paste0('CM_', id),
			cluster = cl,
			group_col = "Sample",
			case_level = t,
			control_level = "WT",
			pseudobulk_level = "T_Reps",
			scale_outcome = scale_outcome
		)
		if (!is.null(res_cm)) {
			res_cm$signature <- id  # nm
			res_cm$score_type <- "CM"
		}

		res_cf <- get_cf4_effect(
			sce = sce,
			cf_col = "zero",
			cm_col = paste0('CF_', id),  ## to reuse the function which output cm_col - cf-col
			cluster = cl,
			group_col = "Sample",
			case_level = t,
			control_level = "WT",
			pseudobulk_level = "T_Reps",
			scale_outcome = scale_outcome
		)
		if (!is.null(res_cf)) {
			res_cf$signature <- id		  # nm
			res_cf$score_type <- "CF"
		}

		bind_rows( res_cm, res_cf)  #bind_rows(res_dual, res_cm, res_cf)
		}))
		if(!is.null(meta_df) && nrow(meta_df) > 0){
			meta_df$cluster = cl
			meta_df$target = t	
			df_plot = bind_rows(df_plot, meta_df)
		} else {
			message(sprintf("No results for cluster=%s and target=%s", cl, t))
		}
	}
}

print(unique(df_plot$cluster)) #  "0" "1" "2" "3" "4"
df_plot$signature %>% unique()   
# [1] Ibarra_1v1.25_STRING_cardiac.a        Ibarra_1v1.25_STRING_cardiac.anoTBX5  Ibarra_1v1.25_IID_cardiac.a          
#  [4] Pijuan_1v1.25_STRING_8                Pijuan_1v1.25_IID_8                     
# [6] Elorbany_1v1.25_STRING_CP.1wincreased Elorbany_1v1.25_IID_CP.1wincreased  

#df_plot$signature = factor(df_plot$signature, levels = c("Elorbany" , "Pijuan_Sala","Ibarra"))K
   # Reverse the levels of signature so that they appear in the opposite order in plots
df_plot$signature = factor(df_plot$signature, levels = rev(unique(df_plot$signature)))

df_plot$target = factor(df_plot$target, levels = sample_levels)  ## for plotting

p = ggplot(df_plot, aes(x = signature, y = yi, color = score_type)) +
	geom_point(position = position_dodge(width = 0.5)) +
	scale_color_manual(values = c("CM-CF" = "grey", "CM" = "blue", "CF" = "green")) +
	geom_errorbar(aes(ymin = yi - 1.96*sei, ymax = yi + 1.96*sei),
					position = position_dodge(width = 0.5), width = 0.2) +
	facet_wrap(target~cluster, ncol=5) +
	theme_bw() +
	geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
	coord_flip() + 
	ggtitle(paste0("louvain clustering effect in Armendariz ", ID, ", scale_outcome=", scale_outcome)) 

if(scale_outcome == TRUE){
	p = p + ylab("Effect size (SD units)")	
} else {
	p = p + ylab("Effect size (raw units)")	
}

pdf(file=paste0('forest_TBX5ENHvsWT_louvain_scale_outcome_', scale_outcome, '.pdf'), width=10, height=7)
print(p)
dev.off()

saveRDS(df_plot, file = paste0('df_plot_TBX5ENHvsWT_louvain_scale_outcome_', scale_outcome, '.rds'))


############################  Fig for publication ####################
df_plot = readRDS( file = 'df_plot_TBX5ENHvsWT_louvain_scale_outcome_FALSE.rds')

## a simplified vesion for publication #######
## Raw units preserve the original magnitude and sign of the score difference.
## SD units standardize by each signature’s variance, so they are better for cross-signature comparability, 
## but they can make some effects look much larger or smaller just because that signature is more or less variable.
ID = 'KO_SET2_SC_DF_processed'
selective_signatures = c(
	'Ibarra_1v1.25_STRING_cardiac.a', 'Ibarra_1v1.25_STRING_cardiac.anoTBX5',  'Ibarra_1v1.25_IID_cardiac.a',          
 	'Pijuan_1v1.25_STRING_8',    'Pijuan_1v1.25_IID_8', 
    'Elorbany_1v1.25_STRING_CP.1wincreased','Elorbany_1v1.25_IID_CP.1wincreased')

scale_outcome = FALSE
df_plot = readRDS(paste0('df_plot_TBX5ENHvsWT_louvain_scale_outcome_', scale_outcome, '.rds'))
df_plot_simplified = df_plot %>% subset(signature %in% selective_signatures) %>%
	subset(target %in% c( 'ENH5'))  %>% 
	subset(cluster %in% c('1', '4'))

p = ggplot(df_plot_simplified, aes(x = signature, y = yi, color = score_type)) +
	geom_point(position = position_dodge(width = 0.5)) +
	scale_color_manual(values = c("CM-CF" = "grey", "CM" = "blue", "CF" = "green")) +
	geom_errorbar(aes(ymin = yi - 1.96*sei, ymax = yi + 1.96*sei),
					position = position_dodge(width = 0.5), width = 0.2) +
	facet_wrap(target~cluster, ncol=5) +
	theme_bw() +
	geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
	coord_flip() + ggtitle(paste0("louvain clustering effect in Armendariz ", ID, ", scale_outcome=", scale_outcome)) 
	
	if(scale_outcome == TRUE){
		p = p + ylab("Effect size (SD units)")	
	} else {
		p = p + ylab("Effect size (raw units)")	
	}

	pdf(file=paste0('selected_forest_TBX5ENHvsWT_louvain_scale_outcome_', scale_outcome, '.pdf'), width=7, height=2.5)
	print(p)
	dev.off()

