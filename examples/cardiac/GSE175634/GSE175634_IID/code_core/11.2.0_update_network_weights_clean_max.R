##############################
## TIPS ANALYSIS on midway3
##############################
# scp -p -r E:\Git_Holly\TIPS-main\BioTIP_update_06162025.R   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/.
# scp -p -r F:\projects\scRNA\source\cardiac_CTS_GRN\GSE175634_iPSC_CM_weighted_v9\IID\11.2.0_update_network_weights_clean_max.R   xyang2@midway3.rcc.uchicago.edu:/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/.
# scp -p -r F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_IID/* xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weighted_IID/.
# scp -p -r xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weighted_IID/* F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_IID/.
# ssh xyang2@midway3.rcc.uchicago.edu

# rcchelp balance  
# rcchelp usage
# squeue -u xyang2
# squeue -p bigmem --state=PD | wc -l   # amd-hm  amd gpu
# squeue -p gpu --state=PD | wc -l

# sinteractive -p amd  --mem=180G  --account=pi-xyang2 --time=3:00:00
# sinteractive -p bigmem  --mem=500G  --account=pi-xyang2 --time=1:00:00  -c 1  
# sinteractive -p caslake  --cpus-per-task=1 --mem=180G --time=6:00:00  --account=pi-xyang2 
# sinteractive -p gpu --account=pi-xyang2 --gres=gpu:1 --mem=180GB  --time=6:00:00 -c 4  (used)
# cd /project/imoskowitz/xyang2/heart_dev/Atlas2_results
# R
## To quit your interactive job:
# exit or Ctrl-D

# module load python/anaconda-2022.05
# source activate /project/xyang2/software-packages/env/velocity_2025Feb_xy

code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_IID", "code_core")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(ppi_path)
library('SingleCellExperiment')
library(Seurat)
library(dplyr)
library(BioTIP)
packageVersion('BioTIP')
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)

sce <- readRDS(paste0(data_dir, "sce.rds"))
sce
dim(sce) # 3000 230786
colnames(colData(sce))
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

### SCT-normalized counts and SCT-normalized residules were calcualted
input = 'Scaledata'  
assayName = ifelse(input == 'Scaledata', 'logcounts', 'counts')

##  using the max score for each edge of IID.db
#########################################################################

graph_list <- readRDS( file= paste0(results_dir,'GSE175634_IID_graph_perState_notsimplified.rds') )
N0 = sapply(graph_list, vcount)


diagnoze=FALSE
if(diagnoze){
	graph_list0 = graph_list

	range(E(graph_list0[[1]])$weight)  #[1] 0.01 0.71
	any(which_multiple(graph_list0[[1]]))
	table(which_multiple(graph_list0[[1]]))
	E(graph_list[[1]])[which_multiple(graph_list0[[1]])]
	edge_attr_names(graph_list0[[1]])

	graph_list <- lapply(graph_list0, simplify) # !!!!!!!!!!!!!!!!!!!
	range(E(graph_list[[1]])$weight)  #[1]  0.01 0.71

	graph_list <- lapply(graph_list0, function(g) simplify(g, edge.attr.comb =  "max"))  
	range(E(graph_list[[1]])$weight)  #[1]  0.01 0.71
}

graph_list <- lapply(graph_list, function(g) simplify(g, edge.attr.comb = list(weight = "max")))
range(E(graph_list[[1]])$weight)  #[1]  0.01 0.71


########## clean 1) remove name-duplicated Vertex due to inconsistence in IID.db , not appl;icable for IID ###########
source(here::here("R", paste0("celltype_specific_weight_v", celltype_specific_weight_version, ".R")))

step1 = TRUE
if(step1){
	## calculating co-expression-based speciricity scores, runs 1 hour when cores=1, DO NOT repeat !
	## first, add a meta column to match the graph_list names
	cluster = as.vector(colData(sce)[,'leiden_0.5'])
	cluster[which(cluster=='3')] = 'CP'
	cluster[which(cluster=='8')] = 'muscle'
	cluster[which(cluster=='11')] = 'endoderm'
	colData(sce)$cluster = factor(cluster, levels = c('0' ,'1','2', 'CP','4', '5', '6', '7','muscle','9', '10','endoderm', '12' ))
	rm(cluster)

	network_specificity_list = calculate_network_specificity(sce,
                                          graph_list,
                                          assayName = assayName, # "logcounts"
                                          celltype_col = "cluster",
                                          method = "pearson",
                                          cores = core_count,
                                          shrink = TRUE,
                                          min_n_Eg = 5)

	saveRDS(network_specificity_list, 'network_specificity_list.rds')
	names(network_specificity_list)
	# [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"
	 # [5] "HiG_4"           "HiG_5"           "HiG_6"           "HiG_7"
	 # [9] "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"
	# [13] "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"
	# [17] "HiGCTS_CP.1"     "CTS_muscle"      "CTS_endoderm"    "CTS_CP"
	# [21] "CTS_CP.1"
	names(network_specificity_list[[1]])
	# [1] "scores"       "genes"        "coexp_target" "corexp_sign"
	names(network_specificity_list[[1]][['scores']])
	# [1] "ratio"    "zscore"   "diff"     "combined"
}

step2 = TRUE
if(step2){
	# list of celltype-xpecific (coexpression) scores
	network_specificity_list <- readRDS("network_specificity_list.rds")
    library(data.table)
    for (net in names(network_specificity_list)) {
        cat("Analyzing: ", net)
        spec_data <- network_specificity_list[[net]]
        corexp_named <- spec_data$corexp_sign

        stopifnot(
            is.matrix(corexp_named),
            all(dim(corexp_named) == dim(spec_data$coexp_target))
        )
        network_specificity_list[[net]]$corexp_sign <- corexp_named
    }
    names(network_specificity_list[[1]])
    # [1] "scores"       "genes"        "coexp_target" "corexp_sign"
    table(network_specificity_list[[1]]$corexp_sign)
    # negative positive
    # 19744    85881   STRING
	# 17828    76421   IID  

	#for(s in c("combined","ratio", "zscore", "diff"))
    s = "combined"
	{										 
		weighted_graph_list <- update_network_weights(graph_list,
									  network_specificity_list,
									  specificity_method = s, #"combined", #"ratio", "zscore", "diff",
									  verbose = FALSE,
									  cores = core_count)
		saveRDS(weighted_graph_list, file = paste0('GSE175634_IID_graph_perState_simplified_',s,'weighted.rds'))
	}
  
	# double check outputs
    g <- weighted_graph_list[[1]]
    vertex_attr_names(g) # [1]  "name"   "weight" "FDR"   
    edge_attr_names(g) #  "weight"         "norm_PPI_score" "corexp_sign"    "coexp_target"


}								  
## copy resutls locally
# scp -p -r  xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPImax_weight/* F:\projects\scRNA\results\cardiac_CTS_GRN\GSE175634_iPSC_CM_weight_PPImax\.
		



## One-time debug block removed here: a manual comparison of this run's
## output against a local backup file (_bk.rds) to confirm reproducibility
## on a prior pass — not part of the pipeline.
