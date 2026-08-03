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

#setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_IID')
setwd('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weighted_IID/')  #!!!!!!!!
library('SingleCellExperiment')
library(Seurat)
library(dplyr) 
# devtools::install_github("xyang2uchicago/BioTIP")
library(BioTIP)
packageVersion('BioTIP')  # 1.16.0
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)

# source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/BioTIP_update_06162025.R')
PPI_color_platte = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

sce = readRDS('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/sce.rds')
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

inputdir = '/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weighted_IID/'   ## shared input !!!!!!!!
graph_list <- readRDS( file= paste0(inputdir,'GSE175634_IID_graph_perState_notsimplified.rds') ) 
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
# source('/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/celltype_specific_weight_v9_xy.R')
source('/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/celltype_specific_weight_v10.R')

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
                                          cores = 4,
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
									  cores = 4) 	
		saveRDS(weighted_graph_list, file = paste0('GSE175634_IID_graph_perState_simplified_',s,'weighted.rds'))
	}
  
	# double check outputs
    g <- weighted_graph_list[[1]]
    vertex_attr_names(g) # [1]  "name"   "weight" "FDR"   
    edge_attr_names(g) #  "weight"         "norm_PPI_score" "corexp_sign"    "coexp_target"


}								  
## copy resutls locally
# scp -p -r  xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPImax_weight/* F:\projects\scRNA\results\cardiac_CTS_GRN\GSE175634_iPSC_CM_weight_PPImax\.
		

###################################################
## compare weights methoids
## check new weights for ISL1 in "HiGCTS_CP.1"	
#####################################################
step3 = TRUE
library(ggplot2)
library(hexbin)


if(step3){
	pdf(file=paste0('compare_specificity_method_vs_PPIscores.pdf'))

	for(net_name in names(graph_list)){
		#net_name = "HiGCTS_CP.1"	
		plot_data = NULL
		par(mfrow=c(2,2))
		#for(s in c("combined","ratio", "zscore", "diff")) {
        s = "combined"
			graph_list = readRDS( file = paste0('GSE175634_IID_graph_perState_simplified_',s,'weighted.rds'))
			g = graph_list[[net_name]]
			
			# Check for ISL1 vertex and mark edges
			isl1_vertex <- which(V(g)$name == 'ISL1')
			is_isl1_edge <- rep(FALSE, length(E(g)))
			if(length(isl1_vertex) > 0) {
				isl1_edges <- incident(g, isl1_vertex, mode = "all")
				is_isl1_edge[isl1_edges] <- TRUE
			}
			
			temp_data <- data.frame(
				net_name=net_name,
				original_weight = E(g)$norm_PPI_score,
				log_weight = log10(E(g)$weight),
				method = s,
				is_isl1 = is_isl1_edge
			)

			plot_data <- rbind(plot_data, temp_data)
		}
		# plot(E(g)$original_weight, log10(E(g)$weight), 
			# main = paste(s, net_name),
			# col = edge_colors ,
			# ptch = )
		p2 <- ggplot(plot_data, aes(x = original_weight, y = log_weight)) +
			stat_binhex(bins = 50, alpha=0.7) +
			geom_point(data = subset(plot_data, is_isl1 == TRUE), 
					   aes(x = original_weight, y = log_weight), 
					   color = "red", size = 1.5) +
			scale_fill_viridis_c() +
			facet_wrap(~method, scales = "free") +
			theme_minimal() +
			labs(title = paste("Edge Weight Density -", net_name),
				 x = "E(g)$norm_PPI_score", 
				 y = "log10(E(g)$weight)",
				 fill = "Count") +
			theme(legend.position = "bottom")

		print(p2)
#	}
	dev.off()
}								  

net_name = "HiGCTS_CP.1"	
plot_data = NULL

#for(s in c("combined","ratio", "zscore", "diff")) {
s = "combined"
	graph_list = readRDS( file = paste0('GSE175634_IID_graph_perState_simplified_',s,'weighted.rds'))
	g = graph_list[[net_name]]
	
	# Check for ISL1 vertex and mark edges
	isl1_vertex <- which(V(g)$name == 'ISL1')
	is_isl1_edge <- rep(FALSE, length(E(g)))
	if(length(isl1_vertex) > 0) {
		isl1_edges <- incident(g, isl1_vertex, mode = "all")
		is_isl1_edge[isl1_edges] <- TRUE
	}
	
	temp_data <- data.frame(
		net_name=net_name,
		original_weight = E(g)$norm_PPI_score,
		log_weight = log10(E(g)$weight),
		method = s,
		is_isl1 = is_isl1_edge
	)

	plot_data <- rbind(plot_data, temp_data)
#}
table(plot_data$is_isl1)
# FALSE
#     6
if(any(plot_data$is_isl1==TRUE)){
    temp = subset(plot_data,is_isl1==TRUE)
    pdf(file='compare_specificity_method_HiGCTS_CP.1_vs_PPIscores.pdf')	
    ggplot(temp, aes(x = original_weight, y = log_weight, color=as.factor(original_weight))) +
                geom_point() +
                scale_fill_viridis_c() +
                facet_wrap(~method, scales = "free") +
                theme_minimal() +
                labs(title = paste("Isl1 linkages -", net_name),
                    x = "E(g)$norm_PPI_score", 
                    y = "log10(E(g)$weight)") +
                theme(legend.position = "bottom")
    dev.off()		
}


###############################
# on local desktop, check if the max-weighted PPI is the same as Holly's 1st run resutls
##################################
library(igraph)
setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_IID/')

graph_list1 = readRDS('GSE175634_IID_graph_perState_simplified_combinedweighted_bk.rds')
graph_list2 = readRDS('GSE175634_IID_graph_perState_simplified_combinedweighted.rds')

for(i in names(graph_list1)){
## all.equal() handles floating-point differences
   if (isFALSE(all.equal(E(graph_list1[[i]])$weight,
                     E(graph_list2[[i]])$weight))) {
   cat(i, "\t")
   }
}

graphs_identical <- function(g1, g2, tol = 1e-12) {
  library(igraph)
  
  # 1. Same number of vertices & edges
  if (vcount(g1) != vcount(g2)) return(FALSE)
  if (ecount(g1) != ecount(g2)) return(FALSE)
  
  # 2. Same vertex names (if present)
  if (!is.null(V(g1)$name) && !is.null(V(g2)$name)) {
    if (!setequal(V(g1)$name, V(g2)$name)) return(FALSE)
  }
  
  # 3. Extract edge lists
  el1 <- ends(g1, E(g1))
  el2 <- ends(g2, E(g2))
  
  # 4. Create canonical edge IDs (order-independent)
  edge_id <- function(el) {
    apply(el, 1, function(x) paste(sort(x), collapse = "_"))
  }
  
  id1 <- edge_id(el1)
  id2 <- edge_id(el2)
  
  # 5. Match edges
  if (!setequal(id1, id2)) return(FALSE)
  
  # 6. Align weights by edge identity
  w1 <- E(g1)$weight
  w2 <- E(g2)$weight
  
  names(w1) <- id1
  names(w2) <- id2
  
  w1 <- w1[order(names(w1))]
  w2 <- w2[order(names(w2))]
  
  # 7. Compare weights with tolerance
  if (!isTRUE(all.equal(w1, w2, tolerance = tol))) return(FALSE)
  
  return(TRUE)
}

for (i in seq_along(graph_list1)) {
  if (graphs_identical(graph_list1[[i]], graph_list2[[i]])) {
    cat(i, "\n")
  }
}
# 1 
# 2 
# 3 
# 4 
# 5 
# 6 
# 7 
# 8 
# 9 
# 10 
# 11 
# 12 
# 13 
# 14 
# 15 
# 16 
# 17 
# 18 
# 19 
# 20 
# 21 


## confirmed they are the same !! good to go ahead 