rebuild_mat = TRUE
CF_cluster = '1'
source('F:/projects/scRNA/source/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9/IID/24.0_acat_load_input_clean.R')

names(graph_list)
#  [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"          
#  [6] "HiG_5"           "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"          
# [11] "HiG_10"          "HiG_endoderm"    "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm"
# [16] "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"      "CTS_endoderm"    "CTS_CP"         
# [21] "CTS_CP.1"   
names(DEG)
#  [1] "0"        "1"        "2"        "CP"       "4"        "5"        "6"        "7"        "muscle"   "9"        "10"       "endoderm"
# [13] "12" 
colnames(mat)
#  [1] "CMES_hi"                          "CP_hi"                           
#  [3] "CM_hi"                            "CF_hi"                           
#  [5] "PCW6CP_access"                    "PCW8_CM_access"                  
#  [7] "PCW19_CM_access"                  "PCW8_CF_access"                  
#  [9] "PCW19_CF_access"                  "PCW8_SMC_access"                 
# [11] "PCW19_SMC_access"                 "PCW6_CM_access"                  
# [13] "PCW6_CF_access"                   "PCW6_SMC_access"                 
# [15] "iEPC_access"                      "CTS_CP.1"                        
# [17] "Maven2023_gene_ISL1_up_E"         "Maven2023_gene_ISL1_up_T"        
# [19] "Maven2023_gene_ISL1_up_L"         "Maven2023_gene_ISL1_dn_E"        
# [21] "Maven2023_gene_ISL1_dn_T"         "Maven2023_gene_ISL1_dn_L"        
# [23] "Maven2023_gene_ISL1_WT_d6CP"      "Gao2019_gene_Isl1_E825E9.bound"  
# [25] "Gao2019_gene_Isl1.iCPC_CPC.bound" "ISL1_CP_bound"                   
# [27] "ISL1_CP_candidate"                "ISL1_CM_candidate"               
# [29] "ISL1_CF_candidate"     

CTS_name = 'CP.1'
CTS_ID = paste0('CTS_', CTS_name)

lengths(CTS)
#   muscle endoderm       CP     CP.1 
#       63       66       76       87 

class(sce) # [1] "SingleCellExperiment"

updir = "F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_IID/GSE181346_heart_scATAC_C3vs1"
path = paste0(updir, '/cisTarget_predicted_',CTS_name)
if(!dir.exists(path)) {dir.create(path)}
setwd(path)


NES_threshold = 4.5
## top 3 pageRankd
df = read.table('../../table_top5_PageRank_perPPI.tsv', header=T, sep='\t') %>% filter(signature == CTS_ID) %>% filter(rank_by_PR <=3)  
(seed_TF = df$gene)
# "FGF10"  "KDR"    "PTPN13"

########################################################
##  input 6 -- data-driven --- RcisTarget enriched motifs among CTS genes
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

mat = read.table(paste0('../binary_annot_',CTS_ID,'_scATAC_Maven2023_gene_ISL1_v3.tsv'), header=T, sep='\t', check.names=F)
dim(mat) #[1] 76 29

cisTarget.res = readRDS( file= "../cisTarget_targets_in_all_CTS.rds")
write.table(subset(cisTarget.res, NES>=NES_threshold & geneSet %in% c('CP','CP.1')) 
				,file=paste0("../cisTarget_targets_in_CTS_CP_CP.1_NES",NES_threshold,".txt"),
				quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')




########################################################
##   data-driven --- binary glag CTS genes 
### understand CTS.CP.1 ######
#for(CTS_ID in c('CTS_CP.1','CTS_CP'))
source('e:/Git_Holly/TIPS/R/celltype_specific_weight_v10.R')

### step 1.2) building the binary annotation matrix 
files = list.files(pattern = '^heatmap_blocked_CTS_') %>% grep('_v3.tsv', ., value=TRUE)

if(length(files) > 0)  {

   
    ## strategy 2) used: teh saem as CTS_CP , simply narrow down to NEs>4
	#tmp = subset(cisTarget.res, geneSet == CTS_name & NES>= NES_threshold)
	#tmp = c(tmp$motif, tmp$TF_highConf)
	motifAnnot_sub = get_regulators_from_motifs(subset(cisTarget.res, geneSet == CTS_name & NES>=NES_threshold), 
											CTS_name, motifAnnot = motifAnnot)
		
	keys = unique(motifAnnot_sub$regulators)
	if(any(is.na(keys))) keys = keys[!is.na(keys)]
	keys
     keys
# "EBF1"  "BCL6B"

	length(keys)  # 2
    
	motifAnnot_sub[which(motifAnnot_sub$regulators %in% keys),] %>% head(10)
	
	hit = lapply(keys, function(x) unlist(stringr::str_split(x, ';'))) %>% unlist %>% unique
    length(hit) #[1] 50
 
    intersect(hit, CTS[[CTS_name]])  # 0, no overlap, all are predicted upstream regulators

	
		# table for CTS.CP ==  high/low in CP ; high/low in CM;  high/low in CF; 
		# open in PCW6; open in PCW8 CM; open in PCW19 CM; open in PCW8 CF; open in PCW19 CF
	for(key in keys){
		motif_TF_highConf = motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf
		
		tmp = subset(cisTarget.res, geneSet == CTS_name & NES>=NES_threshold & (motif==motif_TF_highConf | TF_highConf==motif_TF_highConf))
		genes = unique(unlist(strsplit(tmp$enrichedGenes, ";")))
        mat[, paste0('cisTarget_', key, '.motif_target')] = ifelse(rownames(mat) %in% genes, '1', '0')

	}



		tmp = apply(mat, 2, function(x) sum(as.numeric(x))) 
		which(tmp==0)
		# PCW19_CF_access  PCW8_SMC_access PCW19_SMC_access      iEPC_access 
        #        9               10               11               15 


		
	####################################################
	### extract subnetworks and add PRRX1-bound links ### 
	####################################################
	# literature review to identify topest TF here as a key to followup predictions
	

	dim(mat) #[1] 87  31 all CTS.CP genes !!  for CTS.CP, we should study key TFs such as PRRX1 !!!!!!
	colnames(mat)
#  [1] "CMES_hi"                          "CP_hi"                           
#  [3] "CM_hi"                            "CF_hi"                           
#  [5] "PCW6CP_access"                    "PCW8_CM_access"                  
#  [7] "PCW19_CM_access"                  "PCW8_CF_access"                  
#  [9] "PCW19_CF_access"                  "PCW8_SMC_access"                 
# [11] "PCW19_SMC_access"                 "PCW6_CM_access"                  
# [13] "PCW6_CF_access"                   "PCW6_SMC_access"                 
# [15] "iEPC_access"                      "CTS_CP.1"                        
# [17] "Maven2023_gene_ISL1_up_E"         "Maven2023_gene_ISL1_up_T"        
# [19] "Maven2023_gene_ISL1_up_L"         "Maven2023_gene_ISL1_dn_E"        
# [21] "Maven2023_gene_ISL1_dn_T"         "Maven2023_gene_ISL1_dn_L"        
# [23] "Maven2023_gene_ISL1_WT_d6CP"      "Gao2019_gene_Isl1_E825E9.bound"  
# [25] "Gao2019_gene_Isl1.iCPC_CPC.bound" "ISL1_CP_bound"                   
# [27] "ISL1_CP_candidate"                "ISL1_CM_candidate"               
# [29] "ISL1_CF_candidate"                "cisTarget_EBF1.motif_target"     
# [31] "cisTarget_BCL6B.motif_target"                                                                                                                                                                                              
 
		##  step 1) binary annotation  
		## ### block the genes instead of hcluster (new_strategy: must have PRRX1-CP binding while in this step regardless promoter accessibility
		# USED to be consistent with mouse dataset !! 
		mat = as.data.frame(mat)
		
	## from the enriched motifs find the potential regulators that are either HiG of CTS member of a cluster of itnerest

    x = colnames(mat)[grepl('cisTarget_', colnames(mat))]
    x = gsub('cisTarget_', '', x) %>% gsub('\\.motif_target', '', .)
    x = unlist(strsplit(x, ';')) %>% unique
    x
# [1]"EBF1"  "BCL6B"


   key_TFs = seed_TF
	for(i in x){
		for(j in c(CM_cluster,  CF_cluster,  CP_cluster)){
		   if(i %in% DEG[[j]] & i %in% rownames(sce)) {
			cat(paste0(i, ' is in DEG of ', j, '\n'))
			key_TFs = c(key_TFs, i)
		   }
		}
	}
  # Nothing found

   for(i in x){
		if(i %in% CTS[[CTS_name]]) {
			cat(paste0(i, ' is in CTS of ', CTS_name, '\n'))
			key_TFs = c(key_TFs, i)
		   }
	}
    # Nothing found

	key_TFs= unique(key_TFs) 
	print(key_TFs) #  "FGF10"  "KDR"    "PTPN13"


		## step 1.2) filter out the key_TFs that come from one motif thus will share target genes for the downstream predictions	

	 # to match who is whom
    x = NULL
 	for(j in key_TFs){
	  y = intersect(which(grepl('cisTarget_', colnames(mat))), which(Reduce('|', lapply(j, function(p) grepl(p, colnames(mat), fixed = F)))) )
      cat(j, '\t', y, '\t', colnames(mat)[y], '\n')
	  if(length(y)==0)   y =0
	  x = c(x,y[1])
	 }
	 names(x) = key_TFs
   
# FGF10             
# KDR               
# PTPN13  
   x
# FGF10    KDR PTPN13 
#      0      0      0 

if(any(x>0)) {

 key_TFs  = names(x)[x>0] else int_key_TFs = NULL
	#key_TFs = c('TEAD1','ETS2') #!!!!!!!!!!!!!!,  ETV1  was excluded because it shares the same motif with ETS2/ELF/HOXA2/HOXA32
    x = x[key_TFs]
	for(i in 1:length(x)){
				key = names(x)[i]
				mat[,paste0(key,'_CP_candidate')] = ifelse(mat[,'CP_hi'] == 1 & mat[,x[i]]  == 1 , 1, 0)   
				# 3) must be open at PCW8_CM or PCW19_CM for CM_hi
				mat[,paste0(key,'_CM_candidate')] = ifelse(mat[,'CM_hi'] == 1  & mat[,x[i]]  == 1 , 1, 0) 
				# 4) must be open at PCW8_CF or PCW19_CF for CF_hi
				mat[,paste0(key,'_CF_candidate')] = ifelse( mat[,'CF_hi'] == 1 & mat[,x[i]] == 1 , 1, 0) 
			}
		 
		dim(mat) # dim(mat)[1] 76  56
		
        dim(mat)  #[1] 76 44
colnames(mat)
# [38] "cisTarget_ETV2;ONECUT2.motif_target"                                                                                                                                                                                    
# [39] "TEAD1_CP_candidate"                                                                                                                                                                                                     
# [40] "TEAD1_CM_candidate"                                                                                                                                                                                                     
# [41] "TEAD1_CF_candidate"                                                                                                                                                                                                     
# [42] "ETS2_CP_candidate"                                                                                                                                                                                                      
# [43] "ETS2_CM_candidate"                                                                                                                                                                                                      
# [44] "ETS2_CF_candidate"                                                                                                                                                                                             


		## record the whole mat for all CTS genes, vcan be subset later 
		fileName=paste0('heatmap_blocked_',CTS_ID,'_scATAC_cisTarget_v3.tsv')
		write.table(mat, file=fileName, sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)  #!!!!!!!!!!!!

}   # endo of if(rebuild_mat)

mat = read.table(paste0('heatmap_blocked_',CTS_ID,'_scATAC_cisTarget_v3.tsv'), sep='\t', header=T, check.names=FALSE)

# Step 1.3) heatmap confirming key TF’s self impact by checking its targets among the CTS_CP, candidates to be the highest pagerank TFs !!!!!!!!!!!!!1	
dim(mat)  #[1] 76 44
colnames(mat)
# [38] "cisTarget_ETV2;ONECUT2.motif_target"                                                                                                                                                                                    
# [39] "TEAD1_CP_candidate"                                                                                                                                                                                                     
# [40] "TEAD1_CM_candidate"                                                                                                                                                                                                     
# [41] "TEAD1_CF_candidate"                                                                                                                                                                                                     
# [42] "ETS2_CP_candidate"                                                                                                                                                                                                      
# [43] "ETS2_CM_candidate"                                                                                                                                                                                                      
# [44] "ETS2_CF_candidate"                                                                                                                                                                                    
      

# key_TFs= c('TEAD1','ETS2') 
# fileName=paste0('heatmap_blocked_',CTS_ID,'_scATAC_cisTarget_v3.tsv')
# mat = read.table(fileName, sep='\t',header=T, check.names=F)
# x = c(  30,  35)
# names(x) = key_TFs

for(key in key_TFs){ 
	if(grepl(';',key)) key_in_TFfamily = strsplit(key, ';', fixed=T) %>% unlist %>% intersect(seed_TF) %>% unique else key_in_TFfamily = key

    p = heatmap_pull_candidate(mat, graph_list, CTS_ID, CHD, key=key_in_TFfamily[1], coding_genes =coding_genes, TF = TF_human,
	    show_SMC_access=TRUE)
            ## for Fig 5, only the PRRX1 CP bound genes are shown
    pdf(file=paste0('heatmap_blocked_',CTS_ID,'_scATAC_cisTarget_',key,'_v3_coding_target.pdf'), height=7)  
    print(p)
    dev.off()
 
}
# candidate genes:  8 
# simplified motif name:  cisTarget_CEBPA;CEBPB;CEBPD;CEBPE;CEBPG;PDX1;SOX9;SRY;TEAD1.motif_target 
# candidate genes:  9 
# simplified motif name:  cisTarget_DLX2;DLX3;DRGX;ELF1;ELF2;ELF3;ELF4;ELK1;ELK3;ELK4;ERF;ERG;ETS1;ETS2;ETV1;ETV2;ETV3;ETV4;ETV5;ETV7;FEV;FLI1;FOXO1;GABPA;GABPB1;GATAD1;HOXA2;HOXA3;KMT2B;ONECUT2;SOX15;TFAP2C;ZBTB11;ZBTB25;ZNF580.motif_target 
 
 
 
##########################################################
## -- subset of CTS[['CP']] that are 
## PRRX1-target 
## exclusively highly expressed (HiG) in CM or CF 
## further narrow the subset to be open at CP, but the bone PPIN of CTS at CP does NOT require CP accessibility, considering PRRX1's pioneering role 
 library(BioNet)
 packageVersion('BioNet') # '1.56.0'
 library(igraph)
 library(tibble)
 


mat = read.table(paste0('heatmap_blocked_',CTS_ID,'_scATAC_cisTarget_v3.tsv'), sep='\t', header=T, check.names=FALSE)	
dim(mat)  #[1] 76  44
 
for(key in key_TFs){ 
     if(key=='HOX') key_in_TFfamily = 'HOXB2' else {
	  if(grepl(';',key)) key_in_TFfamily = strsplit(key, ';', fixed=T) %>% unlist %>% intersect(seed_TF) %>% unique else key_in_TFfamily = key
      }

	key_column = which(grepl(key, colnames(mat)) & grepl('cisTarget_', colnames(mat)))

  	graph_TF_list = identify_TF_targeted_pull_candidate(mat, graph_list, CTS_ID, CHD, key= key,
                                keep_selfloop=TRUE, # whether to keep the self-loop of the key
								TF_bound_column_name = key_column,
								TF_appendix = key,
                                edge_colored_by_Maven2023_ISL1KO = FALSE ,
								key_in_TFfamily = key_in_TFfamily) 
									# saveRDS(graph_TF_list, file='PPI_graph_PRRX1_GRN_prediction.rds')
	saveRDS(graph_TF_list, file=paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_v3.rds'))
	names(graph_TF_list)
	#[1] "CTSHiG.CP_TF.target"       "CTS.CP_TF.target_HiGCM"   "CTS.CP_TF.target_HiGCF"
	#[4] "CTSHiG.CP_PRRX1.CP.bound_CPopen"
}
 

for(key in key_TFs){ 
    if(key=='HOX') key_in_TFfamily = 'HOXB2' else {
	  if(grepl(';',key)) key_in_TFfamily = strsplit(key, ';', fixed=T) %>% unlist %>% intersect(seed_TF) %>% unique else key_in_TFfamily = key
      }
    key_column = which(grepl(key, colnames(mat)) & grepl('cisTarget_', colnames(mat)))
 
    graph_TF_list = readRDS(file=paste0('PPI_graph_',key_in_TFfamily,'_GRN_prediction_',CTS_ID,'_v3.rds'))

	plot_TF_targeted_pull_candidate(graph_TF_list, key_in_TFfamily, CTS_ID, saveFigure = TRUE, 
		saveFileName = paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_v3.pdf'))
	# => PPI_graph_ETS2_GRN_prediction_CTS_CP.1_v3.pdf; PPI_graph_EHF;ETS2;SPI1;ZNF639_GRN_prediction_CTS_CP.1_v3.pdf
}


##########################################################
#### step 2.1) quantify the edge weight changes and validate by Maven Isk1_KO in CP results

source('E:/Git_Holly/TIPS/R/celltype_specific_weight_v10.R')


celltype_col  #[1] "leiden_0.5"
CP_cluster # '3'
CM_cluster # '5'
CF_cluster # '1'


# ## heler to paly around the layout to be beautful 
make_layout <- function(g, seed) {
  set.seed(seed)
  layout_with_fr(g)
}

for(key in key_TFs){
    if(key=='HOX') key_in_TFfamily = 'HOXB2' else {
	  if(grepl(';',key)) key_in_TFfamily = strsplit(key, ';', fixed=T) %>% unlist %>% intersect(seed_TF) %>% unique else key_in_TFfamily = key
      }
	key_column = which(grepl(key, colnames(mat)) & grepl('cisTarget_', colnames(mat)))
 

	graph_TF_list = readRDS(file=paste0('PPI_graph_',key_in_TFfamily,'_GRN_prediction_',CTS_ID,'_v3.rds'))
	names(graph_TF_list)
	# [1] "CTSHiG.CP_TF.target"        "CTS.CP_TF.target_HiGCM"     "CTS.CP_TF.target_HiGCF"    
	# [4] "CTSHiG.CP_TF.target_CPopen"
	edge_attr_names(graph_TF_list[["CTS.CP_TF.target_HiGCM"]])
	#[1] "weight"         "corexp_sign"    "coexp_target"   "norm_PPI_score" "color"          "lty"   	

	res = fill_TF_targeting_predicted_edges(graph_TF_list, linkeage_name = 'CM', graph_list, 
											sce, celltype_col = celltype_col , CT_cluster_id = CP_cluster ,
											decendent_cluster_id = CM_cluster, TF_symbol=key_in_TFfamily )
	vcount(res[["g_CT_sub"]])  # 3 !!! for both EZH2 and HOX # NOTHING found !!
    if(vcount(res[["g_CT_sub"]]) > 0)  saveRDS(res, file=paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_CM_final.rds'))
	 

	res = fill_TF_targeting_predicted_edges(graph_TF_list, linkeage_name = 'CF', graph_list, 
											sce, celltype_col = celltype_col ,  CT_cluster_id = CP_cluster ,
											decendent_cluster_id = CF_cluster,  TF_symbol=key_in_TFfamily )
	vcount(res[["g_CT_sub"]]) #  3 for HOXB2
	if(vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file=paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_CF_final.rds'))
	names(res)
	#[1] "g_CT_sub"              "g_decendent_sub"

}



### step 2.2, reporting the results	and visualization #################
(files = list.files(pattern='PPI_graph_.*_GRN_prediction_.*_final.rds')) 	
#[1] "PPI_graph_ETS2_GRN_prediction_CTS_CP_CF_final.rds" 
#[2] "PPI_graph_TEAD1_GRN_prediction_CTS_CP_CF_final.rds"


## NOtice taht the two ETS2 motifs ended with the saem [predictions !!!!!!!!!!!!] 

final_table = NULL

# fwe only got CF-pull predictions for the CTS_CP signature genes,  no for CM pull predictions
#for(pull in c('CF','CM')){
pull = 'CF'
for(key in key_TFs){
    key_in_TFfamily = key
	if(pull=='CF') f = grep(key, files, value = TRUE)  
  
  	res = readRDS(file=f) 
	g1 = res[["g_CT_sub"]]
	g2 = res[["g_decendent_sub"]]

	if(vcount(g1) > 0 & vcount(g2) > 0  & vcount(g1) == vcount(g2)) {
		change_df = edge_change_table(g1=g1, g2=g2, weight_attr = "weight", missing_as = 0, undirected = TRUE)  
		predict = prioritize_edge_change(g1, edge_change_df=change_df, top_n = 5, title= paste0(pull,'-pull subnetwork_',key_in_TFfamily))
	# =>TIPS_delta_edge_reweighting_CF-pull subnetwork_ETS2.pdf
	} else cat('No edges in ',pull,'-pull subnetwork\n')
	# No edges in CF-pull subnetwork

	change_df
	# from    to          w1           w2      delta  abs_delta direction  status rank
	# 1 EZH2 PRRX1 -0.05310376 -0.007506509 0.04559726 0.04559726  increase changed    1
	# 2 EZH2 HOXB4 -0.02289344 -0.003131056 0.01976238 0.01976238  increase changed    2
	#x = colnames(mat)[grepl(key, colnames(mat)) & grepl('cisTarget_', colnames(mat))] %>% gsub('_target','',.) %>% gsub('cisTarget_','',.)
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

dim(final_table) # 3  13
####### generateing final table of the predicted subnetwork #################
 
write.table(final_table, file=paste0('PPI_graph_GRN_prediction_',CTS_ID,'_duallpull_final_table.tsv'),
			quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')

} # end of if(any(x>0))