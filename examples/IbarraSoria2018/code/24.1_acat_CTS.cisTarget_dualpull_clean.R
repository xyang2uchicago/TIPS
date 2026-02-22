# Set True if running 24.0 the first time
rebuild_mat = FALSE
source('/Users/felixyu/Documents/GitHub/TIPS/examples/IbarraSoria2018/code/24.0_acat_load_input_clean.R')

source('../../../../R/celltype_specific_weight_v10.R')

## check the loaded objects =========================
seed_TF # 'ISL1'  # by top PageRank in code 12.xxxx
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

CTS_name # 'cardiac.a'
CTS_ID #  'CTS_cardiac.a'

lengths(CTS)
#       endothelial.b            cardiac.a  
#                  33                   37        
class(sce) # [1] "SingleCellExperiment"

dim(mat) # dim(mat)[1] 37  29

(updir = getwd())
# "/Users/felixyu/Documents/GitHub/TIPS/examples/IbarraSoria2018/results/GSE181346_heart_scATAC"
# create the directory if it doesn't exist?
dir.create(file.path(updir, paste0("cisTarget_predicted_", CTS_name)),
           showWarnings = FALSE, recursive = TRUE)
setwd(paste0(updir, '/cisTarget_predicted_',CTS_name))



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


NES_threshold = 3
cisTarget.res = readRDS( file= "../cisTarget_targets_in_all_CTS.rds")
# (NES_threshold = quantile(cisTarget.res[which(cisTarget.res$geneSet==CTS_name),]$NES, seq(0,1,0.1))['80%'])  #3.814 
write.table(subset(cisTarget.res, NES>= NES_threshold & geneSet %in% c(CP_cluster )) 
			, paste0("../cisTarget_targets_in_cardiac.a_NES", NES_threshold,".txt"),
			quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    




########################################################
##  --- binary glag CTS genes 
### understand CTS.CP.1 ######
#for(CTS_ID in c('CTS_CP.1','CTS_CP'))
dim(mat) #[1]  37  29
### step 1.2) building the binary annotation matrix 
files = list.files(pattern = '^heatmap_blocked_CTS_') %>% grep('_v3.tsv', ., value=TRUE)
if(length(files) > 0){
    fileName =  grep('_v3.tsv', files, value=TRUE)
	mat = read.table(fileName, sep='\t', header=T, check.names=FALSE)

    saved_variables = readRDS(file = "scATAC_cisTarget_variables.rds")

    # Recreate variables
    x = saved_variables$x
    key_TFs = saved_variables$key_TFs
    motifAnnot_sub = saved_variables$motifAnnot_sub

} else {
	#tmp = subset(cisTarget.res, geneSet == CTS_name & NES>= NES_threshold)
	#tmp = c(tmp$motif, tmp$TF_highConf)  #!!!!!!!!
	motifAnnot_sub = get_regulators_from_motifs(cisTarget.res, CTS_name, NES_threshold, motifAnnot = motifAnnot)
	keys = unique(motifAnnot_sub$regulators)
	if(any(is.na(keys))) keys = keys[!is.na(keys)]
	keys
# [1] "PRDM9"                                                                                                                            
 # [2] "CBFB"                                                                                                                             
 # [3] "ZNF582"                                                                                                                           
 # [4] "CRX;DBP;EP300;GATA1;GATA2;GATA3;GATA4;GATA5;GATA6;GTF2IRD1;HNF1A;HOXA13;IKZF2;MYB;MZF1;NFIA;NFIC;OTX1;OTX2;POU2F1;POU5F1;RAX;ZEB1"
 # [5] "ATF2;IKZF1"                                                                                                                       
 # [6] "PARP1"                                                                                                                            
 # [7] "ZNF408"                                                                                                                           
 # [8] "ZNF682"                                                                                                                           
 # [9] "ATF5"                                                                                                                             
# [10] "ZNF784"                                                                                                                           
# [11] "HMGA1;HMGA2"                                                                                                                      
# [12] "ZKSCAN7"                                                                                                                          
# [13] "ZNF134"                                                                                                                           
# [14] "CEBPA;CEBPB;CEBPD;CEBPE;CEBPG;DLX2;GCM1;HOXA2;HOXA3;TEAD4"                                                                        
# [15] "ZNF98"                                                                                                                            
# [16] "ZNF716"                                                                                                                           
# [17] "ZNF382"                                                                                                                           
# [18] "IRF5"                                                                                                                             
# [19] "ZFP62"                                                                                                                            
# [20] "ZBTB33"                                                                                                                           
# [21] "EVX1"                                                                                                                             
# [22] "ZNF343"                                                                                                                           
# [23] "ZFP92"                                                                                                                            
# [24] "HSF5"                                                                                                                             
# [25] "ZNF426"                                                                                                                           
# [26] "ZNF555"                                                                                                                           
# [27] "ZNF791"                                                                                                                           
# [28] "ZKSCAN2"                                                                                                                          
# [29] "NR4A1"                                                                                                                            
# [30] "ZNF34"                                                                                                                            
# [31] "ZNF77"                                                                                                                            
# [32] "ZNF429"                                                                                                                           
# [33] "ELK1;ERF;ERG;ETV2;ETV5;ETV7;FLI1;SPDEF"                                                                                           
# [34] "PRDM6"                                                                                                                            
# [35] "ZNF250"                                                                                                                           
# [36] "HIC1"                                                                                                                             
# [37] "GZF1"                                                                                                                             
# [38] "KLF3"                                                                                                                             
# [39] "ZFP64"                                                                                                                            
# [40] "GZF1;HSFY1;HSFY2"                                                                                                                 
# [41] "ZNF783"                                                                                                                           
# [42] "IRF4"                                                                                                                             
# [43] "ZSCAN9"                                                                                                                           
# [44] "ZNF671"                                                                                                                           
# [45] "RBFOX2"                                                                                                                           
# [46] "ZNF580"                                                                                                                           
# [47] "ZNF71"                                                                                                                            
# [48] "EP300;GTF2IRD1"                                                                                                                   
# [49] "ATF2"                                                                                                                             
# [50] "CEBPD;TEAD4"                                                                                                                      
# [51] "ERF;ETV7"                                                                                                                         
# [52] "GZF1;HSFY2"                                                                                                                       
# [53] "ETV2;ETV7"         
	length(keys) # 53

	
	## add the motif-based target for each key
	for(key in keys){
		motif_TF_highConf = motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf
		
		tmp = subset(cisTarget.res, geneSet == CTS_name & NES>=NES_threshold & (motif==motif_TF_highConf | TF_highConf==motif_TF_highConf))
		genes = unique(unlist(strsplit(tmp$enrichedGenes, ";")))
		mat[, paste0('cisTarget_', key, '.motif_target')] = ifelse(rownames(mat) %in% genes, '1', '0')
		}

	# 0
		
	####################################################
	### extract subnetworks and add seed_key-bound links ### 
	####################################################
	# literature review to identify topest TF here as a key to followup predictions
	

	dim(mat) #[1] 37  82  all CTS.cardiac.a genes !!  
	colnames(mat)
# [1] "CMES_hi"                                                                                                                                                 
 # [2] "CP_hi"                                                                                                                                                   
 # [3] "CM_hi"                                                                                                                                                   
 # [4] "CF_hi"                                                                                                                                                   
 # [5] "PCW6CP_access"                                                                                                                                           
 # [6] "PCW8_CM_access"                                                                                                                                          
 # [7] "PCW19_CM_access"                                                                                                                                         
 # [8] "PCW8_CF_access"                                                                                                                                          
 # [9] "PCW19_CF_access"                                                                                                                                         
# .....

	## from the enriched motifs find the potential regulators that are either HiG of CTS member of a cluster of itnerest
	x = colnames(mat)[grepl('cisTarget_', colnames(mat))]
	x = gsub('cisTarget_', '', x) %>% gsub('\\.motif_target', '', .)
	x = unlist(strsplit(x, ';')) %>% unique
	x
 # [1] "PRDM9"    "CBFB"     "ZNF582"   "CRX"      "DBP"      "EP300"    "GATA1"    "GATA2"    "GATA3"   
# [10] "GATA4"    "GATA5"    "GATA6"    "GTF2IRD1" "HNF1A"    "HOXA13"   "IKZF2"    "MYB"      "MZF1"    
# [19] "NFIA"     "NFIC"     "OTX1"     "OTX2"     "POU2F1"   "POU5F1"   "RAX"      "ZEB1"     "ATF2"    
# [28] "IKZF1"    "PARP1"    "ZNF408"   "ZNF682"   "ATF5"     "ZNF784"   "HMGA1"    "HMGA2"    "ZKSCAN7" 
# [37] "ZNF134"   "CEBPA"    "CEBPB"    "CEBPD"    "CEBPE"    "CEBPG"    "DLX2"     "GCM1"     "HOXA2"   
# [46] "HOXA3"    "TEAD4"    "ZNF98"    "ZNF716"   "ZNF382"   "IRF5"     "ZFP62"    "ZBTB33"   "EVX1"    
# [55] "ZNF343"   "ZFP92"    "HSF5"     "ZNF426"   "ZNF555"   "ZNF791"   "ZKSCAN2"  "NR4A1"    "ZNF34"   
# [64] "ZNF77"    "ZNF429"   "ELK1"     "ERF"      "ERG"      "ETV2"     "ETV5"     "ETV7"     "FLI1"    
# [73] "SPDEF"    "PRDM6"    "ZNF250"   "HIC1"     "GZF1"     "KLF3"     "ZFP64"    "HSFY1"    "HSFY2"   
# [82] "ZNF783"   "IRF4"     "ZSCAN9"   "ZNF671"   "RBFOX2"   "ZNF580"   "ZNF71"    

 key_TFs = seed_TF
	for(i in x){
		for(j in c(CM_cluster,  CF_cluster,  CP_cluster)){
		   if(i %in% DEG[[j]] & i %in% rownames(sce)) {
			cat(paste0(i, ' is in DEG of ', j, '\n'))
			key_TFs = c(key_TFs, i)
		   }
		}
	}
# CBFB is in DEG of cardiac.c
# CBFB is in DEG of extraembryonicMesoderm
# CBFB is in DEG of cardiac.a
# GATA4 is in DEG of cardiac.c
# GATA4 is in DEG of cardiac.a
# GATA5 is in DEG of cardiac.c
# GATA5 is in DEG of cardiac.a
# GATA6 is in DEG of cardiac.c
# GATA6 is in DEG of cardiac.a
# HMGA2 is in DEG of cardiac.c
# HMGA2 is in DEG of extraembryonicMesoderm
# HMGA2 is in DEG of cardiac.a
# ETV2 is in DEG of extraembryonicMesoderm
# PRDM6 is in DEG of extraembryonicMesoderm
# PRDM6 is in DEG of cardiac.a
# RBFOX2 is in DEG of cardiac.c
# RBFOX2 is in DEG of cardiac.a


for(i in x){
		if(i %in% CTS[[CTS_name]]) {
			cat(paste0(i, ' is in CTS of ', CTS_name, '\n'))
			key_TFs = c(key_TFs, i)
		   }
	}
# GATA4 is in CTS of cardiac.a
# GATA6 is in CTS of cardiac.a


	key_TFs= unique(key_TFs) 
	print(key_TFs) #  "MEF2C"  "GATA4"  "MSX2"   "CBFB"   "GATA5"  "GATA6"  "HMGA2"  "ETV2"   "PRDM6"  "RBFOX2"

	
	##  step 1) binary annotation  ; narrow seed_TF candidates to those having target genes highly expressed in CP, CM, or CF
		## ### block the genes instead of hcluster (new_strategy: must have PRRX1-CP binding while in this step regardless promoter accessibility
		# USED to be consistent with mouse dataset !! 
	mat = as.data.frame(mat)
		
		# one label per CTS&HiG module gene to prioritize candidates CTS member genes
		## refer to 12.0_rank_by_PageRank_BC.R !!!!!!!!!!!
		# 1) must be CP_hi, CM_hi or CF_hi
		# 2) must be targets of seed_TF

	# note that 'All Hox proteins bind this same consensus site with approximately equal affinity'
	# https://pmc.ncbi.nlm.nih.gov/articles/PMC10216783/#:~:text=All%20Hox%20proteins%20bind%20this,et%20al.%2C%202010).

	## step 1.2) filter out the key_TFs that come from one motif thus will share target genes for the downstream predictions	
	(x = intersect(which(grepl('cisTarget_', colnames(mat))), which(Reduce('|', lapply(key_TFs, function(p) grepl(p, colnames(mat), fixed = F)))) ))
	#[1] 31 33 40 62 63 74 82
    
	# to match who is whom
    x = NULL
 	for(j in key_TFs){
	  y = intersect(which(grepl('cisTarget_', colnames(mat))), which(Reduce('|', lapply(j, function(p) grepl(p, colnames(mat), fixed = F)))) )
      cat(j, '\t', y, '\t', colnames(mat)[y], '\n')
	  if(length(y)==0)   y =0
	  x = c(x,y[1])
	 }
	 names(x) = key_TFs
# MEF2C             
# GATA4    33      cisTarget_CRX;DBP;EP300;GATA1;GATA2;GATA3;GATA4;GATA5;GATA6;GTF2IRD1;HNF1A;HOXA13;IKZF2;MYB;MZF1;NFIA;NFIC;OTX1;OTX2;POU2F1;POU5F1;RAX;ZEB1.motif_target 
# MSX2              
# CBFB     31      cisTarget_CBFB.motif_target 
# GATA5    33      cisTarget_CRX;DBP;EP300;GATA1;GATA2;GATA3;GATA4;GATA5;GATA6;GTF2IRD1;HNF1A;HOXA13;IKZF2;MYB;MZF1;NFIA;NFIC;OTX1;OTX2;POU2F1;POU5F1;RAX;ZEB1.motif_target 
# GATA6    33      cisTarget_CRX;DBP;EP300;GATA1;GATA2;GATA3;GATA4;GATA5;GATA6;GTF2IRD1;HNF1A;HOXA13;IKZF2;MYB;MZF1;NFIA;NFIC;OTX1;OTX2;POU2F1;POU5F1;RAX;ZEB1.motif_target 
# HMGA2    40      cisTarget_HMGA1;HMGA2.motif_target 
# ETV2     62 82   cisTarget_ELK1;ERF;ERG;ETV2;ETV5;ETV7;FLI1;SPDEF.motif_target cisTarget_ETV2;ETV7.motif_target 
# PRDM6    63      cisTarget_PRDM6.motif_target 
# RBFOX2   74      cisTarget_RBFOX2.motif_target     (x = x[which(x>0)])
   x 
# MEF2C  GATA4   MSX2   CBFB  GATA5  GATA6  HMGA2   ETV2  PRDM6 RBFOX2 
     # 0     33      0     31     33     33     40     62     63     74 
  
    key_TFs = c("GATA4" ,  "CBFB" ,  "HMGA2" , "ETV2" ,  "PRDM6" , "RBFOX2")
	x = x[key_TFs]
    if(length(x)>0) {
        for(j in seq_along(x)){
		# if(grepl(";",colnames(mat)[x[j]])){
            # tmp = strsplit(colnames(mat)[x[j]], ";")[[1]] %>% gsub('\\.motif_target', '', .) %>% gsub('^cisTarget_', '', .)
            # names(x)[j] = intersect(seed_TF, tmp) %>% unique
            # } else {
        # names(x)[j] = intersect(seed_TF, colnames(mat)[x[j]])
        # }
     
        key = names(x)[j]
        mat[,paste0(key,'_CP_candidate')] = ifelse(mat[,'CP_hi'] == 1 & mat[,x[j]]  == 1 , 1, 0)   
        # 3) must be open at PCW8_CM or PCW19_CM for CM_hi
        mat[,paste0(key,'_CM_candidate')] = ifelse(mat[,'CM_hi'] == 1  & mat[,x[j]]  == 1 , 1, 0) 
        # 4) must be open at PCW8_CF or PCW19_CF for CF_hi
        mat[,paste0(key,'_CF_candidate')] = ifelse( mat[,'CF_hi'] == 1 & mat[,x[j]] == 1 , 1, 0) 
        } # end of for(j in 1:lengths(x))
        #key_TFs = c(key_TFs, key)
    } # end of if(length(x)>0)
	
    dim(mat) # 37  100

	cat(paste0('key_TFs: ', paste(key_TFs, collapse='_'), '\n'))  # key_TFs: 	GATA4_CBFB_HMGA2_ETV2_PRDM6_RBFOX2

		## record the whole mat for all CTS genes, vcan be subset later 
	if(length(key_TFs) > 0) {
		fileName=paste0('heatmap_blocked_',CTS_ID,'_scATAC_cisTarget_',paste(key_TFs, collapse='_'),'_v3.tsv')
		write.table(mat, file=fileName, sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)  #!!!!!!!!!!!!
        saveRDS(list(x = x, key_TFs = key_TFs, motifAnnot_sub = motifAnnot_sub), "scATAC_cisTarget_variables.rds")
	} else stop('No key TFs found for ', CTS_ID, '\n')
} 


# Step 1.3) heatmap confirming key TF’s self impact by checking its targets among the CTS_CP, candidates to be the highest pagerank TFs !!!!!!!!!!!!!1	


motif_TF_highConf = gsub('cisTarget_', '', colnames(mat)[x] ) %>% gsub('\\.motif_target', '', .)
for(key in motif_TF_highConf){  #!!!!!!!!
        if(grepl(';',key)) key_in_TFfamily = strsplit(key, ';', fixed=T) %>% unlist %>% intersect(key_TFs) %>% unique else key_in_TFfamily = key

        p = heatmap_pull_candidate(mat, graph_list, CTS_ID, CHD, key=key_in_TFfamily, coding_genes =coding_genes, TF = TF_human,
			show_SMC_access=TRUE)
            ## for Fig 5, only the PRRX1 CP bound genes are shown
        pdf(file=paste0('heatmap_blocked_',CTS_ID,'_scATAC_cisTarget_',key_in_TFfamily,'_v3_coding_target.pdf'), height=4)  
        print(p)
        dev.off()

}
#candidate genes:  4 
# simplified motif name:  cisTarget_CRX;DBP;EP300;GATA1;GATA2;GATA3;GATA4;GATA5;GATA6;GTF2IRD1;HNF1A;HOXA13;IKZF2;MYB;MZF1;NFIA;NFIC;OTX1;OTX2;POU2F1;POU5F1;RAX;ZEB1.motif_target 
# candidate genes:  6 
# direct motif:  cisTarget_CBFB.motif_target  used from multiple potentoal matches
# candidate genes:  12 
# simplified motif name:  cisTarget_HMGA1;HMGA2.motif_target 
# candidate genes:  2 
# simplified motif name:  cisTarget_ELK1;ERF;ERG;ETV2;ETV5;ETV7;FLI1;SPDEF.motif_target cisTarget_ETV2;ETV7.motif_target 
# candidate genes:  3 
# direct motif:  cisTarget_PRDM6.motif_target  used from multiple potentoal matches
# candidate genes:  4 
# direct motif:  cisTarget_RBFOX2.motif_target  used from multiple potentoal matches

##########################################################
## -- subset of CTS[['CP']] that are 
## PRRX1-target 
## exclusively highly expressed (HiG) in CM or CF 
## further narrow the subset to be open at CP, but the bone PPIN of CTS at CP does NOT require CP accessibility, considering PRRX1's pioneering role 
library(BioNet)
packageVersion('BioNet') # '1.56.0'
library(igraph)
library(tibble)
 

mat = read.table(paste0('heatmap_blocked_',CTS_ID,'_scATAC_cisTarget_',paste(key_TFs, collapse='_'),'_v3.tsv'), sep='\t', header=T, check.names=FALSE)	
dim(mat)  #[1] 54  71
 
for(key in key_TFs){
	key_column = which(grepl(key, colnames(mat)) & grepl('cisTarget_', colnames(mat)))
    if(key=='HOX') key_in_TFfamily = 'HOXB2' else key_in_TFfamily = key
    
    ## further narrow the subset to be open at CP  or descendents
  	graph_TF_list = identify_TF_targeted_pull_candidate(mat, graph_list, CTS_ID, CHD, key= key,
                                keep_selfloop=TRUE, # whether to keep the self-loop of the key
								TF_bound_column_name = key_column,
								TF_appendix = key,
                                edge_colored_by_Maven2023_ISL1KO = FALSE ,
								key_in_TFfamily = key_in_TFfamily) 
									# saveRDS(graph_TF_list, file='PPI_graph_PRRX1_GRN_prediction.rds')
	saveRDS(graph_TF_list, file=paste0('PPI_graph_',key_in_TFfamily,'_GRN_prediction_',CTS_ID,'_v3.rds'))
    }

names(graph_TF_list)
	#[1] "CTSHiG.CP_TF.target"       "CTS.CP_TF.target_HiGCM"   "CTS.CP_TF.target_HiGCF"
	#[4] "CTSHiG.CP_PRRX1.CP.bound_CPopen"
 

for(key in key_TFs){
	if(key=='HOX') key_in_TFfamily = 'HOXB2' else key_in_TFfamily = key

    graph_TF_list = readRDS(file=paste0('PPI_graph_',key_in_TFfamily,'_GRN_prediction_',CTS_ID,'_v3.rds'))
    plot_TF_targeted_pull_candidate(graph_TF_list, key_in_TFfamily, CTS_ID, saveFigure = TRUE)

}
	# => PPI_graph_GATA4_GRN_prediction_CTS_CP_v3.pdf;  
 




##########################################################
#### step 2.1) quantify the edge weight changes and validate by Maven Isk1_KO in CP results ===================

# ## TO get teh cluster ID, remove the .1, .2, etc. from the cluster names which indicating additional CTS identified from the same cluster
CM_ID = paste0('HiG_',CM_cluster) %>% gsub('\\.[1-9]','',.)
CF_ID = paste0('HiG_',CF_cluster) %>% gsub('\\.[1-9]','',.)

# ## heler to paly around the layout to be beautful 
make_layout <- function(g, seed) {
  set.seed(seed)
  layout_with_fr(g)
}

for(key in key_TFs){
	if(key=='HOX') key_in_TFfamily = 'HOXB2' else key_in_TFfamily = key  

    graph_TF_list = readRDS(file=paste0('PPI_graph_',key_in_TFfamily,'_GRN_prediction_',CTS_ID,'_v3.rds'))
    names(graph_TF_list)
    # [1] "CTSHiG.CP_TF.target"        "CTS.CP_TF.target_HiGCM"     "CTS.CP_TF.target_HiGCF"    
    # [4] "CTSHiG.CP_TF.target_CPopen"
    edge_attr_names(graph_TF_list[["CTS.CP_TF.target_HiGCM"]])
    #[1] "weight"         "corexp_sign"    "coexp_target"   "norm_PPI_score" "color"          "lty"   	

    res = fill_TF_targeting_predicted_edges(graph_TF_list, linkeage_name = 'CM', graph_list, 
                                            sce, celltype_col = celltype_col , CT_cluster_id = CP_cluster ,
                                            descendent_cluster_id = CM_cluster, TF_symbol=key_in_TFfamily,
                                            HVG=rownames(sce) )
    cat(paste0(key,' CM vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), '\n'))  # 2
    if(vcount(res[["g_CT_sub"]]) > 0)  saveRDS(res, file=paste0('PPI_graph_',key_in_TFfamily,'_GRN_prediction_',CTS_ID,'_CM_final.rds'))
    

    res = fill_TF_targeting_predicted_edges(graph_TF_list, linkeage_name = 'CF', graph_list, 
                                            sce, celltype_col = celltype_col ,  CT_cluster_id = CP_cluster ,
                                            descendent_cluster_id = CF_cluster,  TF_symbol=key_in_TFfamily,
                                            HVG=rownames(sce) )
    cat(paste0(key, ' CF vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), '\n'))  #  0  
    if(vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file=paste0('PPI_graph_',key_in_TFfamily,'_GRN_prediction_',CTS_ID,'_CF_final.rds'))
    names(res)
    #[1] "g_CT_sub"              "g_descendent_sub"
  
}
# GATA4 CM vcount(res[["g_CT_sub"]]): 2
# GATA4 CF vcount(res[["g_CT_sub"]]): 0
# CBFB CM vcount(res[["g_CT_sub"]]): 3
# CBFB CF vcount(res[["g_CT_sub"]]): 0
# HMGA2 CM vcount(res[["g_CT_sub"]]): 5
# HMGA2 CF vcount(res[["g_CT_sub"]]): 4
# ETV2 CM vcount(res[["g_CT_sub"]]): 2
# ETV2 CF vcount(res[["g_CT_sub"]]): 0
# PRDM6 CM vcount(res[["g_CT_sub"]]): 3
# PRDM6 CF vcount(res[["g_CT_sub"]]): 0
# RBFOX2 CM vcount(res[["g_CT_sub"]]): 2
# RBFOX2 CF vcount(res[["g_CT_sub"]]): 0


### step 2.2, reporting the results	and visualization #################
(files = list.files(pattern='PPI_graph_.*_GRN_prediction_.*_final.rds')) 	
# [1] "PPI_graph_CBFB_GRN_prediction_CTS_cardiac.a_CM_final.rds"  
# [2] "PPI_graph_ETV2_GRN_prediction_CTS_cardiac.a_CM_final.rds"  
# [3] "PPI_graph_GATA4_GRN_prediction_CTS_cardiac.a_CM_final.rds" 
# [4] "PPI_graph_HMGA2_GRN_prediction_CTS_cardiac.a_CF_final.rds" 
# [5] "PPI_graph_HMGA2_GRN_prediction_CTS_cardiac.a_CM_final.rds" 
# [6] "PPI_graph_PRDM6_GRN_prediction_CTS_cardiac.a_CM_final.rds" 
# [7] "PPI_graph_RBFOX2_GRN_prediction_CTS_cardiac.a_CM_final.rds"
 

# for CF-pull with results, reporting and plotting

final_table = NULL
# fwe only got CF-pull predictions for the CTS_CP signature genes,  no for CM pull predictions
for(f in files){
	pull = ifelse(grepl('CF', f),'CF', 'CM')
    pattern <- paste(key_TFs, collapse = "|")
    key = key_in_TFfamily = regmatches(f, regexpr(pattern, f))
	
    res = readRDS(file= f)
    g1 = res[["g_CT_sub"]]
    g2 = res[["g_descendent_sub"]]

    if(vcount(g1) > 0 & vcount(g2) > 0  & vcount(g1) == vcount(g2)) {
		change_df = edge_change_table(g1=g1, g2=g2, weight_attr = "weight", missing_as = 0, undirected = TRUE)  
		predict = prioritize_edge_change(g1, edge_change_df=change_df, top_n = 5, title= paste0(pull,'-pull subnetwork_',key_in_TFfamily))
	# =>TIPS_delta_edge_reweighting_CF-pull subnetwork_ETS2.pdf
	} else cat('No edges in ',pull,'-pull subnetwork\n')
	# No edges in CF-pull subnetwork


     change_df
    #    from    to         w1         w2       delta  abs_delta direction  status rank
    # 1 GATA4 MYOCD 0.03826759 0.01032437 -0.02794322 0.02794322  decrease changed    1

	x = paste0(key, '_motif_target')
	(y = motifAnnot_sub[which(motifAnnot_sub$regulators %in% key | grepl(key, motifAnnot_sub$regulators)),])
# 	           motif_TF_highConf regulators parsed_from_id mapped_from_annot
# 10 ETS2 (directAnnotation).        ETS2           ETS2              <NA>
    motif_TF_highConf = y$motif_TF_highConf
	tmp = subset(cisTarget.res, geneSet == CTS_name & NES>=NES_threshold & (motif==motif_TF_highConf | TF_highConf==motif_TF_highConf))
    all_regulators = strsplit(y$regulators, ';', fixed=T) %>% unlist %>% unique
    
    change_df = cbind(linkeage = pull, change_df, TF_highConf = tmp$TF_highConf, motif = tmp$motif, NES = tmp$NES)
    change_df$TF_highConf[which(change_df$from !=key_in_TFfamily & change_df$to !=key_in_TFfamily)] = ''
	change_df$motif[which(change_df$from !=key_in_TFfamily & change_df$to !=key_in_TFfamily)] = ''
	change_df$NES[which(change_df$from !=key_in_TFfamily & change_df$to !=key_in_TFfamily)] = ''
	
	final_table = rbind(final_table, change_df)

}

dim(final_table) # 21  13
####### generateing final table of the predicted subnetwork #################
 
write.table(final_table, file=paste0('PPI_graph_GRN_prediction_',CTS_ID,'_dualpull_final_table.tsv'),
			quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
