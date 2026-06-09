rebuild_mat = FALSE
source('F:/projects/scRNA/source/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_IID/24.0_acat_load_input_clean.R')

names(graph_list)
#  [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"          
#  [6] "HiG_5"           "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"          
# [11] "HiG_10"          "HiG_endoderm"    "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm"
# [16] "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"      "CTS_endoderm"    "CTS_CP"         
# [21] "CTS_CP.1"   

names(DEG)
#   [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12"

lengths(CTS)
#   muscle endoderm       CP     CP.1 
#       63       66       76       87 

	 
####################################################
### extract subnetworks and add ISL1-bound links ### 
####################################################

up_dir ="F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_IID/GSE181346_heart_scATAC_C3vs1"
#for(CTS_name in c('CP', 'CP.1')) {
CTS_name = 'CP.1'
	CTS_ID = paste0('CTS_',CTS_name)
	# refer to code 12.0_rank_by_PageRank_BC.R
	if(CTS_name == 'CP.1') seed_TF = c("ISL1" , "HAND2") else seed_TF = c('HOXB2') ## ETS2 was not identified but always cobinding with ISL1 and HAND2

	dir.create(file.path(up_dir, paste0("ChIPseq_predicted_", CTS_name)),
           showWarnings = FALSE, recursive = TRUE)
	setwd(paste0(up_dir, '/ChIPseq_predicted_',CTS_name))
	
	fileName = paste0('../binary_annot_',CTS_ID,'_scATAC_Maven2023_gene_ISL1_v3.tsv')
	mat = read.table(fileName, header=TRUE, sep='\t', check.names=FALSE)
	dim(mat) #[1] 87 29


    key = 'ISL1'
	##  step 1) binary annotation  
	## ### block the genes instead of hcluster (new_strategy: must have ISL1-CP binding while in this step regardless promoter accessibility
	# USED to be consistent with mouse dataset !! 
	mat = as.data.frame(mat)
	dim(mat)  #[1] 87 29
	colnames(mat)
# [1] "CMES_hi"                          "CP_hi"                            "CM_hi"                           
#  [4] "CF_hi"                            "PCW6CP_access"                    "PCW8_CM_access"                  
#  [7] "PCW19_CM_access"                  "PCW8_CF_access"                   "PCW19_CF_access"                 
# [10] "PCW8_SMC_access"                  "PCW19_SMC_access"                 "PCW6_CM_access"                  
# [13] "PCW6_CF_access"                   "PCW6_SMC_access"                  "iEPC_access"                     
# [16] "CTS_CP"                           "Maven2023_gene_ISL1_up_E"         "Maven2023_gene_ISL1_up_T"        
# [19] "Maven2023_gene_ISL1_up_L"         "Maven2023_gene_ISL1_dn_E"         "Maven2023_gene_ISL1_dn_T"        
# [22] "Maven2023_gene_ISL1_dn_L"         "Maven2023_gene_ISL1_WT_d6CP"      "Gao2019_gene_Isl1_E825E9.bound"  
# [25] "Gao2019_gene_Isl1.iCPC_CPC.bound" "ISL1_CP_bound"                    "ISL1_CP_candidate"               
# [28] "ISL1_CM_candidate"                "ISL1_CF_candidate"                      

   source('e:/Git_Holly/TIPS/R/celltype_specific_weight_v10.R')

# Step 1.3) heatmap confirming key TF’s self impact by checking its targets among the CTS_CP, candidates to be the highest pagerank TFs !!!!!!!!!!!!!1	
	p = heatmap_pull_candidate(mat, graph_list, CTS_ID, CHD, key=key, coding_genes =coding_genes, TF = TF_human,
         chip_targets =TRUE, show_SMC_access=TRUE)
# candidate genes:  57 

	pdf(file=paste0('heatmap_blocked_',CTS_ID,'_scATAC_chIPseq_',key,'_v3_coding_target.pdf'), height=7)  
	print(p)
	dev.off()



##########################################################
## -- subset of CTS[['CP.1']] that are 
## ISL1-target 
## exclusively highly expressed (HiG) in CM or CF 
## further narrow the subset to be open at CP, but the bone PPIN of CTS at CP does NOT require CP accessibility, considering Isl1's pioneering role 
library(BioNet)
packageVersion('BioNet') # '1.56.0'
library(igraph)
library(tibble)
library(gplots) 
#for(CTS_ID in c('CTS_CP.1','CTS_CP'))


source('E:/Git_Holly/TIPS/R/celltype_specific_weight_v10.R')

{
	# Step 1.4) subPPI for targets of the key_TF
	dim(mat)  #[1] 87 229
 
    key_target_column = 'ISL1_CP_bound'
	graph_TF_list = identify_TF_targeted_pull_candidate(mat, graph_list, CTS_ID, CHD, key= key,
                                keep_selfloop=TRUE, # whether to keep the self-loop of the key
								TF_bound_column_name = key_target_column,
								TF_appendix = key,
                                edge_colored_by_Maven2023_ISL1KO = TRUE,
								key_in_TFfamily = key) 
									 
	saveRDS(graph_TF_list, file=paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_v3.rds'))
	lengths(graph_TF_list)
    #      CTSHiG.CP_TF.target     CTS.CP_TF.target_HiGCM     CTS.CP_TF.target_HiGCF CTSHiG.CP_TF.target_CPopen 
    #                    27                          13                          14                         10 



	# Step 1.5) venn diagram comparing the keyTF targets at CP stage and two terminal lineage-specific stages
    graph_TF_list = readRDS(file=paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_v3.rds'))
	
	plot_TF_targeted_pull_candidate(graph_TF_list, key, CTS_ID, saveFigure = TRUE)
	# => PPI_graph_ISL1_GRN_prediction_CTS_CP.1_v3.pdf
}



##########################################################
#### step 2.1) quantify the edge weight changes and validate by Maven Isk1_KO in CP results

source('E:/Git_Holly/TIPS/R/celltype_specific_weight_v10.R')


celltype_col  #[1] "leiden_0.5"
CP_cluster #'3'
CM_cluster # '5'
CF_cluster # '1'

# if(grepl('.', key, fixed=T)) warnning(paste('pls assign gene symbol to ', key))
# if(key== "SOX9.TEAD1") key_symbol = 'TEAD1' else key_symbol = key

# ## heler to paly around the layout to be beautful 
make_layout <- function(g, seed) {
  set.seed(seed)
  layout_with_fr(g)
}


graph_TF_list = readRDS(file=paste0('PPI_graph_ISL1_GRN_prediction_',CTS_ID,'_v3.rds'))
names(graph_TF_list)
# [1] "CTSHiG.CP_TF.target"        "CTS.CP_TF.target_HiGCM"    
# [3] "CTS.CP_TF.target_HiGCF"     "CTSHiG.CP_TF.target_CPopen"
edge_attr_names(graph_TF_list[["CTS.CP_TF.target_HiGCM"]])
#[1] "weight"         "corexp_sign"    "coexp_target"   "norm_PPI_score" "color"          "lty"   	

res = fill_TF_targeting_predicted_edges(graph_TF_list, linkeage_name = 'CM', graph_list, 
                                            sce, celltype_col = celltype_col , CT_cluster_id = CP_cluster ,
                                            decendent_cluster_id = CM_cluster, TF_symbol=key,
                                            HVG=rownames(sce) )
cat(paste0(key,' CM vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), '\n'))  
# ISL1 CM vcount(res[["g_CT_sub"]]): 12

if(vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file=paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_CM_final.rds'))
names(res)
#[1] "g_CT_sub"              "g_decendent_sub"


res = fill_TF_targeting_predicted_edges(graph_TF_list, linkeage_name = 'CF', graph_list, 
										sce, celltype_col = celltype_col , CT_cluster_id = CP_cluster ,
                                        decendent_cluster_id = CF_cluster, TF_symbol=key,
                                        HVG=rownames(sce) )
cat(paste0(key,' CF vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), '\n'))  
# ISL1 CF vcount(res[["g_CT_sub"]]): 114

if(vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file=paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_CF_final.rds'))


(files = list.files(pattern='PPI_graph_.*_GRN_prediction_.*_final.rds')) 	
#[1] "PPI_graph_ISL1_GRN_prediction_CTS_CP.1_CF_final.rds"
# [2] "PPI_graph_ISL1_GRN_prediction_CTS_CP.1_CM_final.rds"

### step 2.2, reporting the results	and visualization #################

final_table = NULL
for( pull in c('CM', 'CF')) {
    
    res = readRDS(file= paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_',pull,'_final.rds')) 
    g1 = res[["g_CT_sub"]]
    g2 = res[["g_decendent_sub"]]

    if(vcount(g1) > 0 & vcount(g2) > 0  & vcount(g1) == vcount(g2)) {
         change_df = edge_change_table(g1=g1, g2=g2, weight_attr = "weight", missing_as = 0, undirected = TRUE)  
        } else cat(paste0('No edges in ', pull, '-pull subnetwork\n'))
      

    predict = prioritize_edge_change(g1, edge_change_df = change_df, top_n = 5, title= paste0(pull, '_pull_subnetwork_',key) )
    # => TIPS_delta_edge_reweighting_CF_pull_subnetwork.pdf
    
    tmp = cbind(linkeage = pull, change_df, ChIP = key)
     
    final_table = rbind(final_table, tmp)

}

dim(final_table) # 32 11
####### generateing final table of the predicted subnetwork #################
 
write.table(final_table, file=paste0('PPI_graph_',key,'_GRN_prediction_',CTS_ID,'_dualpull_final_table.tsv'),
			quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
