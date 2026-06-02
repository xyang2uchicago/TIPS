rebuild_mat = FALSE
source('F:/projects/scRNA/source/cardiac_CTS_GRN/IbarraSoria2018_E8.25_v9/24.0_acat_load_input_clean.R')
seed_TF  #"MEF2C", "GATA4", "MSX2"
names(graph_list)
#[1] "HiG_blood"                  "HiG_cardiac.b"              "HiG_cardiac.c"             
#  [4] "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"         
#  [7] "HiG_extraembryonicMesoderm" "HiG_mesodermProgenitors"    "HiG_mixedMesoderm.a"       
# [10] "HiG_mixedMesoderm.b"        "HiG_pharyngealMesoderm"     "HiG_presomiticMesoderm.a"  
# [13] "HiG_presomiticMesoderm.b"   "HiG_somiticMesoderm"        "HiG_endothelial.b"         
# [16] "HiG_cardiac.a"              "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"          
# [19] "CTS_endothelial.b"          "CTS_cardiac.a"     
names(DEG)
#   [1] "blood"                  "cardiac.b"              "cardiac.c"              "endothelial.a"         
#  [5] "endothelial.c"          "endothelial.d"          "extraembryonicMesoderm" "mesodermProgenitors"   
#  [9] "mixedMesoderm.a"        "mixedMesoderm.b"        "pharyngealMesoderm"     "presomiticMesoderm.a"  
# [13] "presomiticMesoderm.b"   "somiticMesoderm"        "endothelial.b"          "cardiac.a"     


celltype_col  #[1] "subcelltype"
CP_cluster # 'cardiac.a'
CM_cluster # 'cardiac.c'
CF_cluster # 'extraembryonicMesoderm'
CMES_cluster # 'mixedMesoderm.a'

CTS_name = 'cardiac.a'
CTS_ID = paste0('CTS_', CTS_name)

lengths(CTS)
#       endothelial.b            cardiac.a  
#                  33                   37            

################################################################
# load iPSC gene expression matrix as an independent evaluation matrix 
# db_specifc_input_path = 'F:/projects/TIPS/data/GSE175634_iPSC_CM_weighted_v9/'

## this is an singlecellExpression object 
# sce = readRDS(paste0(db_specifc_input_path, 'sce.rds'))
# dim(sce)  # 3000 230786
# names(colData(sce))
# meta = colData(sce0)
# saveRDS(meta, file = 'F:/projects/scRNA/results/GSE175634_iPSC_CM/meta_with_leiden_pseudotime.rds')

## this is a Seurat object !!
library(Seurat)
db_specifc_input_path = 'F:/projects/scRNA/results/GSE175634_iPSC_CM/'
sce = readRDS(paste0(db_specifc_input_path,"sc_SCT_normalized_count_only.rds"))
sce
# An object of class Seurat 
# 38847 features across 230786 samples within 1 assay 
# Active assay: SCT (38847 features, 3000 variable features)
# 1 layer present: data

colnames(sce@meta.data)
 # [1] "orig.ident"     "nCount_RNA"     "nFeature_RNA"   "cell"          
 # [5] "exp.grp"        "sample"         "diffday"        "individual"    
 # [9] "S.Score"        "G2M.Score"      "CC.Difference"  "demux.dbl.prb" 
# [13] "leiden"         "type"           "dpt_pseudotime" "percent.mt"    
# [17] "nCount_SCT"     "nFeature_SCT"  

meta = readRDS(paste0(db_specifc_input_path,"meta_with_leiden_pseudotime.rds"))
meta = as.data.frame(meta)
all(rownames(meta) == rownames(sce@meta.data)) # TRUE
all(meta$cell == sce@meta.data$cell)  # TRUE
sce <- AddMetaData(sce, metadata = meta[, setdiff(colnames(meta), colnames(sce@meta.data))])

colnames(sce@meta.data)
 # [1] "orig.ident"                       "nCount_RNA"                      
 # [3] "nFeature_RNA"                     "cell"                            
 # [5] "exp.grp"                          "sample"                          
 # [7] "diffday"                          "individual"                      
 # [9] "S.Score"                          "G2M.Score"                       
# [11] "CC.Difference"                    "demux.dbl.prb"                   
# [13] "leiden_published"                 "type"                            
# [15] "percent.mt"                       "nCount_SCT"                      
# [17] "nFeature_SCT"                     "PC1"                             
# [19] "PC2"                              "leiden_0.1"                      
# [21] "leiden_0.2"                       "leiden_0.3"                      
# [23] "leiden_0.4"                       "leiden_0.5"                      
# [25] "leiden_0.6"                       "leiden_0.7"                      
# [27] "leiden_0.8"                       "leiden_0.9"                      
# [29] "leiden_1"                         "leiden_0.5_type"                 
# [31] "dpt_pseudotime_published"         "dpt_pseudotime_leiden0.5_root382"

library(SingleCellExperiment)
sce <- as.SingleCellExperiment(sce)
sce
# class: SingleCellExperiment 
# dim: 38847 230786 
# metadata(0):
# assays(1): logcounts
# rownames(38847): WASH7P AL627309.1 ... AC007325.4 AC007325.2
# rowData names(0):
# colnames(230786): 1 2 ... 230785 230786
# colData names(35): orig.ident nCount_RNA ...
  # dpt_pseudotime_leiden0.5_root382 ident
# reducedDimNames(0):
# mainExpName: SCT
# altExpNames(0):
colData(sce)$diffday %>% unique
#[1] "day1"  "day3"  "day7"  "day5"  "day11" "day15" "day0" 
colData(sce)$diffday = factor(colData(sce)$diffday, levels = paste0('day', c('0', '1','3','5','7','11','15')))


celltype_col = "leiden_0.5"  
CP_cluster = '3'
CM_cluster = '5'
CF_cluster = '1'  # mesenchymal/ECM‑high (“CF‑like”).
CMES_cluster = '4'  # Mesp1+

cols <- setNames(
  c("orange", "blue", "lightgreen"),
  c(CP_cluster, CM_cluster, CF_cluster)
)


################################################################


# set local database-specific working directory
updir = "F:/projects/scRNA/results/cardiac_CTS_GRN/IbarraSoria2018_E8.25_v9/GSE181346_heart_scATAC"
setwd(updir)
setwd('../validation')

library(igraph)
library(dplyr)
library(igraph)
library(openxlsx)
library(gplots)
library(ggExtra)


pull_path = 'F:/projects/TIPS/doc/'
pull_df = read.xlsx(paste0(pull_path, 'STable_final_2026.xlsx'), sheet=15)
dim(pull_df)  # 122  15
head(pull_df, 3)
     # Database    CTS_ID linkeage     x     y                 w_CP
# 1 IbarraSoria cardiac.a       CM GATA6  TBX5  0.14273739867444099
# 2 IbarraSoria cardiac.a       CM  CBFB GATA6  3.65114162434605E-2
# 3 IbarraSoria cardiac.a       CM  CBFB  TBX5 -1.28698510207887E-2
             # w_decendent                  delta             abs_delta
# 1  5.3873131230046403E-2 -8.8864267444394796E-2 8.8864267444394796E-2
# 2 -5.6903109139440902E-3 -4.2201727157404598E-2 4.2201727157404598E-2
# 3 -2.1015662268983001E-4    1.26596943980988E-2   1.26596943980988E-2
  # direction  status rank               TF_highConf      motif  NES
# 1  decrease changed    1                      <NA>       <NA> <NA>
# 2  decrease changed    2 CBFB (directAnnotation).  hdpi__CBFB 4.67
# 3  increase changed    3 CBFB (directAnnotation).  hdpi__CBFB 4.67
table(pull_df$Database)

   # Elorbany IbarraSoria Pijuan_Sala 
         # 79          31          12 

db = 'IbarraSoria'  ## test the Elorbany CTS !!!!!!!!

## for validation, we narow to validate the linege-specific (not shared) nodese !!!!!!!!!!!!!!!!!!
pull_df = subset(pull_df, direction %in% c('decrease' ,'increase'))
pull_df = subset(pull_df, Database == db)
table(pull_df$direction, pull_df$linkeage)
          # CF CM
  # decrease  1 22
  # increase  4  4



############ validation # 1 (program usage vs pull score) ###############################
## Takeuchi 2025 sgRNA, defined 250 Gene Programs (GPs) responsived to TF perturbation. ###############
## Gene Program 5 (GP5) is highestly enriched in CHD genes  (padj = 1.6e-09), including 
## PRRX1, PBX3, HEY1, HOXB4, HOXB3 (padj=1e-4) that we identified as a CTS (n=76) in CP from iPSC dataset.
## GP5 also includse GATA6, ZEB2, IFI16, HAND2 (padj=5e-4) that we identified as another CTS (n=87) in CP from iPSC dataset.
GP_path = 'D:/projects/DS/data/TF_perturb_CM/doc/'

## reported CHD-associated gene programs
df = read.xlsx(xlsxFile = paste0(GP_path, 'media-2.xlsx'), sheet=11, startRow = 2)  # 250 programs
dim(df)  #[1] 17  5
(CHD_gp = df$Program) %>% as.character  
 # [1] "5"   "1"   "10"  "28"  "2"   "96"  "13"  "32"  "9"   "31"  "7"   "29" 
# [13] "16"  "12"  "6"   "45"  "241"


df = read.xlsx(xlsxFile = paste0(GP_path, 'media-2.xlsx'), sheet=8, startRow = 2)  # 250 programs
dim(df)  #[1] 300  251
df = df[, -1]
head(df[,1:10], 3)
      # 1       2       3       4     5     6       7      8     9      10
# 1  FBN2 CACNA1D  NCCRP1  MT-CO3 SFRP1 GRIK1 ZNF385D COL1A2  LMO2   H3F3B
# 2 LAMA2   RSPO3 ZNF385B  MT-CO2 RAB3C   CKB   PRR16  ITGA8 DPP10 GADD45G
# 3  MYL3     LBH   TBX18 MT-ATP8 PDE3A MFGE8   FABP5 COL1A1 CADM1  PHLDA1

length(unique(unlist(df))) # 25311
length(union(rownames(sce), unique(unlist(df)))) # 42032
universe = union(rownames(sce), unique(unlist(df)))
length(universe) # 42032 df and iPSC-shared universe, used to track PG gene members


library(dplyr)
library(tidyr)
source('E:/Git_Holly/TIPS/R/celltype_specific_weight_v10.R')

df_long <- df %>%
  pivot_longer(
    cols = everything(),
    names_to = "program",
    values_to = "gene_name"
  ) %>%
  filter(!is.na(gene_name) & gene_name != "") %>%   # remove empty cells
  distinct() %>% as.data.frame()

head(df_long)
# # A tibble: 6 × 2
  # program gene_name   
  # <chr>   <chr>  
# 1 1       FBN2   
# 2 2       CACNA1D
		

############# validation # 2 are the two subsets of GP5 tracks the CP→CM trajectory (or CP->CF trajectory) ?##########################################
#### turn GP5 into a quantitative score in the iPSC scRNA-seq and 
#### these overlapped genes behave differently across fates:
############################################################### 
# library(AUCell)  # AUCell becomes unstable for 4-5 genes
library(ggplot2)
library(Matrix)
library(dplyr)
library(ggpubr)
 
library(patchwork)

# -> GP5 and GP10 are reversed, PG10 may be regulated transcriptomic program


############# valdiation 3) temporal specificity / rewiring window  #############
## by quantifying lineage bias axis 
#ID = paste0('GP5_colMeans_', CTS_name)
## helper to format slope/angle text from lm()


####### loop over dot_type and gp #########################################################
#for(dot_type in c('cell','pseudobulk')){
dot_type = 'pseudobulk'

## against all GS genes not only asn CTS members 
ann_all = NULL
gp_names = c('CM-diff GP1','EMT-like GP28','Top CHD GP5'); names(gp_names) = c('1','28','5')

cf_pull_genes_0 <- subset(pull_df, Database==db & linkeage=='CF'  & CTS_ID==CTS_name)[, c('x','y')] %>%
			unlist() %>% unique()
cm_pull_genes_0 <- subset(pull_df, Database==db & linkeage=='CM'  & CTS_ID==CTS_name)[, c('x','y')] %>%
		unlist() %>% unique()
dual_pull_genes = intersect(cf_pull_genes_0, cm_pull_genes_0 )		
cf_pull_genes = setdiff(cf_pull_genes_0, dual_pull_genes) %>% intersect(., rownames(sce)) #!!!ADDED
cm_pull_genes = setdiff(cm_pull_genes_0, dual_pull_genes) %>% intersect(., rownames(sce))
	
cf_pull_name = paste0('CF_pull', CTS_name)			#!!!ADDED
cm_pull_name = paste0('CM_pull', CTS_name)			
dual_pull_name = paste0('dual_pull', CTS_name)	

if (db == "Elorbany" & CTS_name=='CP.1') {
	cm_pull_genes <- unique(c(cm_pull_genes, "FGF10", "LRRTM1"))
	cm_pull_name <- paste0(cm_pull_name, "_wincreased_nodes")
  }  
if (db == "Pijuan_Sala") {
	cm_pull_genes <- unique(c(cm_pull_genes, "HMGA2"))
	cm_pull_name <- paste0(cm_pull_name, "_wincreased_nodes")
  }          

colData(sce)[[cf_pull_name]] <- Matrix::colMeans(logcounts(sce)[cf_pull_genes, , drop = FALSE])
colData(sce)[[cm_pull_name]] <- Matrix::colMeans(logcounts(sce)[cm_pull_genes, , drop = FALSE])					
colData(sce)[[dual_pull_name]] <- Matrix::colMeans(logcounts(sce)[dual_pull_genes, , drop = FALSE])				


for(gp in c('1', '28','5')){ 

	print(CTS_name) # 'cardiac.a'   ## this is IbarraSoria CTS
	gp_present <- intersect(df[[as.character(gp)]], rownames(sce))		
	n_gp_present = length(gp_present)
			
	if(n_gp_present > 0){
		ID = paste0('GP',gp,'_in_sce')
		colData(sce)[[ID]] <- Matrix::colMeans(logcounts(sce)[gp_present, , drop = FALSE])		
		summary(colData(sce)[[ID]])
		# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
		# -0.7992 -0.2982 -0.1911  0.0000 -0.1103 13.8979 

		#### (pseudobulk, donor-aware)
		#  aggregate by individual × day × celltype.
		if(dot_type == 'pseudobulk') {
        df_score <- as.data.frame(colData(sce)) %>%
                dplyr::mutate(celltype = .data[[celltype_col]]) %>%
                dplyr::mutate(score = .data[[ID]]) %>%
                dplyr::mutate(cf_score = .data[[cf_pull_name]]) %>%
                dplyr::mutate(cm_score = .data[[cm_pull_name]]) %>%
                dplyr::mutate(dual_score = .data[[dual_pull_name]]) %>%
                dplyr::group_by(individual, diffday, celltype) %>%
                dplyr::summarise(
                    score_GPi = mean(score, na.rm = TRUE),
                    cf_pull_score =  mean(cf_score, na.rm = TRUE),
                    cm_pull_score =  mean(cm_score, na.rm = TRUE),
                    dual_pull_score =  mean(dual_score, na.rm = TRUE),
                    n_cells = n(),
                    .groups = "drop"
                    )  
                }
        if(dot_type == 'cell') {
            df_score <- as.data.frame(colData(sce)) %>%
                dplyr::mutate(celltype = .data[[celltype_col]]) %>%
                dplyr::mutate(score_GPi = .data[[ID]]) %>%
                dplyr::mutate(cf_pull_score = .data[[cf_pull_name]]) %>%
                dplyr::mutate(cm_pull_score = .data[[cm_pull_name]]) 
            }

		df_sub <- subset(df_score, celltype %in% c(CP_cluster, CM_cluster, CF_cluster))
        df_sub$individual <- factor(df_sub$individual)   # <<< ADDED
        df_sub$diffday    <- factor(df_sub$diffday)      # <<< ADDED
		df_sub$celltype <- factor(df_sub$celltype, levels =  c(CP_cluster, CM_cluster, CF_cluster))
		head(df_sub, 3)

		# fit <- lm(score_GPi ~ cm_pull_score + cf_pull_score + celltype + diffday + individual, data = df_sub)
		# summary(fit)
 	      	
		dot_size = ifelse(dot_type == 'pseudobulk', 1.5, 0.8)
		use_wt <- FALSE  # <<< ADDED TURE has bug unsolved !!
		ann_cf = ann_cm = ann_dual = NULL
    
		lineage = 'CF'
		if(length(cf_pull_genes)>0) {
            ann_cf <- make_ann_by_day_celltype(df_sub,
                xvar = "cf_pull_score",
                yvar = "score_GPi",
                p_cut = 0.05,
                use_weights = use_wt # ,weights_var = "n_cells"
                )  # <<< ADDED
				
			if(nrow(ann_cf)>0) {
			  ann_cf$lineage = 'CF'
              p_tmp1 =  ggplot(df_sub,  aes(x = cf_pull_score, y = score_GPi, color = celltype)) +
				geom_point(size = dot_size) +
				geom_smooth(method = "lm", se = FALSE) +
				scale_color_manual(values = cols, breaks = names(cols)) +
				theme_bw() +
				facet_wrap(~ diffday, nrow=1) + 
				labs(x = paste0( length(cf_pull_genes), "-gene CF-pull score: ", toString(cf_pull_genes) ),
                    y = paste0(n_gp_present, "-gene GP", gp, " signature score (cell-level)")) +
				ggtitle(paste0('CTS_',CTS_name, ": ", dot_type, ' level')) +
                coord_cartesian(clip = "off") +                                 # <<< ADDED
                theme(plot.margin = margin(5.5, 5.5, 20, 5.5)) +                # <<< ADDED (room for text)
                geom_text(                                                       # <<< ADDED
                    data = ann_cf,
                    aes(x = x, y = y, label = label, color = celltype),
                    inherit.aes = FALSE,
                    hjust = 0, vjust = 1, size = 2.5
                )
			  } else {
			  ann_cf = NULL
			  p_tmp1 =   ggplot() + theme_void() + labs(title = paste0('CTS_',CTS_name, ": ", lineage, '-bias -- No data'))
			  }
			  
			} else p_tmp1 =   ggplot() + theme_void() + labs(title = paste0('CTS_',CTS_name, ": ", lineage, '-bias -- No data'))
		
  	 
		lineage = 'CM'  # in fact this is redundant, just flip the x-axis compared to p_tmp1 !!
		if(length(cm_pull_genes)>0) {
            ann_cm <- make_ann_by_day_celltype(
                df_sub,
                xvar = "cm_pull_score",
                yvar = "score_GPi",
                p_cut = 0.05,
                use_weights = use_wt#, weights_var = "n_cells"
                )  # <<< ADDED
			ann_cm$lineage = 'CM'
            p_tmp2 = ggplot(df_sub, aes(x = cm_pull_score, y = score_GPi, color = celltype)) +
				geom_point(size = dot_size) +
				geom_smooth(method = "lm", se = FALSE) +
				scale_color_manual(values = cols, breaks = names(cols)) +
				theme_bw() +
				facet_wrap(~ diffday, nrow=1) +
				labs(x = paste0( length(unique(cm_pull_genes)), "-gene CM-pull score: ", toString(cm_pull_genes) ),
                    y = paste0(length(unique(gp_present)), "-gene GP", gp, " signature score (cell-level)")) +
				ggtitle(paste0('CTS_',CTS_name, ": ", dot_type, ' level')) +
                coord_cartesian(clip = "off") +                                 # <<< ADDED
                    theme(plot.margin = margin(5.5, 5.5, 20, 5.5)) +                # <<< ADDED
                    geom_text(                                                       # <<< ADDED
                        data = ann_cm,
                        aes(x = x, y = y, label = label, color = celltype),
                        inherit.aes = FALSE,
                        hjust = 0, vjust = 1, size = 2.5
                    )
			} else p_tmp2 =  ggplot() + theme_void() + labs(title = paste0('CTS_',CTS_name, ": ", lineage, '-bias -- No data'))

        lineage = 'dual'
        if(length(dual_pull_genes)>0) {
          ann_dual <- make_ann_by_day_celltype(df_sub,
            xvar = "dual_pull_score",
            yvar = "score_GPi",
            p_cut = 0.05,
            use_weights = use_wt#   ,weights_var = "n_cells"
			)  # <<< ADDED
			ann_dual$lineage = 'dual'
			p_tmp3 = ggplot(df_sub, aes(x = dual_pull_score, y = score_GPi, color = celltype)) +
				geom_point(size = dot_size) +
				geom_smooth(method = "lm", se = FALSE) +
				scale_color_manual(values = cols, breaks = names(cols)) +
				theme_bw() +
				facet_wrap(~ diffday, nrow=1) +
				labs(x = paste0( length(dual_pull_genes), "-gene dual-pull score: ", toString(dual_pull_genes) ),
					y = paste0(n_gp_present, "-gene GP", gp, " signature score (cell-level)")) +
				ggtitle(paste0('CTS_',CTS_name, ": ", dot_type, ' level')) +

				coord_cartesian(clip = "off") +                                 # <<< ADDED
					theme(plot.margin = margin(5.5, 5.5, 20, 5.5)) +                # <<< ADDED
					geom_text(                                                       # <<< ADDED
						data = ann_dual,
						aes(x = x, y = y, label = label, color = celltype),
						inherit.aes = FALSE,
						hjust = 0, vjust = 1, size = 2.5
					)
        } else p_tmp3 =  ggplot() + theme_void() + labs(title = paste0('CTS_',CTS_name, ": ", lineage, '-bias -- No data'))


		if(dot_type == 'pseudobulk') {
			Filename = paste0('GSE175634_iPSC_regression_GP',gp,'_vs_',CTS_name,'_perday_perindividual_',dot_type,'.pdf')
		} else {
			Filename = paste0('GSE175634_iPSC_regression_GP',gp,'_vs_',CTS_name,'_',dot_type,'.pdf')
		}

    pdf(file= Filename, height=12, width= 15)
 	gridExtra::grid.arrange(p_tmp1, p_tmp2, p_tmp3, ncol = 1)
	dev.off()
    
    ann = rbind(ann_cf, ann_cm, ann_dual)
    ann$gp_name = gp_names[gp]
    ann$CTS_name = CTS_name
    ann$GPid = paste0('GP',gp)
    ann$dot_type = dot_type

	}	# end of gp_present loo
  ann_all = rbind(ann_all, ann)
} # end of gp loop
#}	# end of CTS_name loop

dim(ann_all)  # 39  16
saveRDS(ann_all, file = 'ann_all_regression_GP_vs_CTS.rds')


#############
  ## Panel 2: iPSC time-course summary
  # in scatter plot:
##   slope (β) = regression coefficient estimate from GO/program score ~ CTS-pull score
# in summarizing heatmap:
##  standardized β = scale-normalized coupling effect, comparable across evaluation types

# Instead of showing all scatter facets, compute per-day × celltype:
# slope (β) for GO score ~ pull score
# show β as color (blue negative / red positive)
# show significance as a dot or outline
# optionally show n pseudobulks as dot size or transparency
# This compresses 14 scatterplots into one clean message:
# delivery:
# “coupling peaks at day5 and is celltype selective”
library("ggnewscale")
library(forcats)
day_levels <- c("day0","day1","day3","day5","day7","day11","day15")

ann_all = readRDS('ann_all_regression_GP_vs_CTS.rds')
dim(ann_all)  # 39  16

colnames(ann_all) 
#  [1] "diffday"  "celltype" "x"        "y"        "label"    "beta"     "beta_std" "t"        "slope"   
# [10] "p"        "r2"       "lineage"  "gp_name"  "CTS_name" "GPid"     "dot_type"


table(ann_all$CTS_name, ann_all$lineage)
            # CF CM dual
  # cardiac.a  6 14   19
 
table(ann_all$CTS_name, ann_all$GPid)
        # GP1 GP28 GP5
  # cardiac.a  14   15  10

library(dplyr)

ann_plot = ann_all %>%
  dplyr::filter(gp_name %in% c('CM-diff GP1','EMT-like GP28')) %>%
  dplyr::mutate(celltype = factor(celltype))
ann_plot$lineage = paste0(ann_plot$lineage , '_bias')

# summarize in case duplicates exist
ann_summary <- ann_plot %>%
  dplyr::group_by(diffday, celltype, gp_name, lineage) %>%
  dplyr::summarise(
    beta_std = mean(beta_std, na.rm = TRUE),
    p = min(p, na.rm = TRUE),
    # r2 = mean(r2, na.rm = TRUE),
	# t = mean(t, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(diffday = factor(diffday,levels = day_levels ) )

  
ann_summary <- ann_summary %>%
  dplyr::mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    )
  ) 
  
ann_line <- ann_summary %>%
  mutate(
    logp = -log10(p),
    celltype = factor(celltype,
                      levels = c(CP_cluster, CM_cluster, CF_cluster),
                      labels = c("CP","CM","CF/EMT")),
	lineage = factor(lineage, levels = c ("dual_bias", "CM_bias", "CF_bias"    ))
  )

cell_cols <- cols; names(cell_cols) = c("CP","CM","CF/EMT")

p_line <- ggplot(
  ann_line,
  aes(x = diffday, y = beta_std,
      color = celltype,
      group = celltype)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.8, alpha = 0.7) +
   # only point signficant ones
  geom_point(data = ann_line %>% dplyr::filter(p < 0.05), aes(size = logp), alpha = 0.9) +  
  facet_grid(gp_name ~ lineage) +
  scale_x_discrete(limits = day_levels, drop = FALSE) +
  scale_color_manual(values = cell_cols) +
  scale_size(range = c(2,6), name = "-log10(p)") +
  labs(
    x = "Differentiation day",
    y = "Standardized coupling β",
    color = "Cell type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

p_line

pdf(file='linedot_dualpull_GP_sig.pdf', width=7, height=5)
print(p_line)
dev.off()

