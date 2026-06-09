### In this version, we teste the v9 code for GSE175634_iPSC_CM
### Additiobnally, the max(STRING original scores) in the simplifiy()  


# keep the 6/6/2025 transfer the infor needed to remove duplicated vertex by setting inpuitdir to PPI_weight
# scp -p -r F:\projects\scRNA\results\cardiac_CTS_GRN\GSE175634_iPSC_CM_weight\correct_n_edges_HiG_STRING2.14.0.rds xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weight/.



##############################
## TIPS ANALYSIS on midway3
##############################
# scp -p -r E:\Git_Holly\TIPS-main\BioTIP_update_06162025.R   xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/.
# scp -p -r xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weighted_v9/* F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9/.
# ssh xyang2@midway3.rcc.uchicago.edu

# rcchelp balance  
# rcchelp usage
# squeue -u xyang2
# squeue -p bigmem --state=PD | wc -l   # amd-hm  amd gpu
# squeue -p gpu --state=PD | wc -l

# sinteractive -p amd-hm  --mem=180G  --account=pi-xyang2 --time=3:00:00
# sinteractive -p bigmem  --mem=500G  --account=pi-xyang2 --time=1:00:00  -c 1  
# sinteractive -p caslake  --cpus-per-task=1 --mem=180G --time=6:00:00  --account=pi-xyang2 
# sinteractive -p gpu --account=pi-xyang2 --gres=gpu:1 --mem=180GB  --time=6:00:00 -c 4  (used)
# cd /project/imoskowitz/xyang2/heart_dev/Atlas2_results
# R
## To quit your interactive job:
# exit or Ctrl-D

# module load python/anaconda-2022.05
# source activate /project/xyang2/software-packages/env/velocity_2025Feb_xy

#setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weight_shrink')
setwd('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weighted_v9/')  #!!!!!!!!
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
library(stringr)

# source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/BioTIP_update_06162025.R')
PPI_color_platte = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

sce = readRDS('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP/sce.rds')
sce
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

##  using the max score for each edge of STRING.db
#########################################################################
diagnoze=FALSE
if(diagnoze){
	graph_list0 = graph_list

	range(E(graph_list0[[1]])$weight)  #[1] 0.200 0.999
	any(which_multiple(graph_list0[[1]]))
	table(which_multiple(graph_list0[[1]]))
	E(graph_list[[1]])[which_multiple(graph_list0[[1]])]
	edge_attr_names(graph_list0[[1]])

	graph_list <- lapply(graph_list0, simplify) # !!!!!!!!!!!!!!!!!!!
	range(E(graph_list[[1]])$weight)  #[1] 0.400 1.998

	graph_list <- lapply(graph_list0, simplify, edge.attr.comb = 'max')  # by default is 'sum'
	range(E(graph_list[[1]])$weight)  #[1]  0.200 0.999
}



# scp -p -r F:\projects\scRNA\results\cardiac_CTS_GRN\GSE175634_iPSC_CM_weight\GSE175634_STRING_graph_perState_notsimplified.rds xyang2@midway3.rcc.uchicago.edu:/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weight/.
inputdir = '/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/PPI_weight/'   ## shared input !!!!!!!!
graph_list <- readRDS( file= paste0(inputdir,'GSE175634_STRING_graph_perState_notsimplified.rds') ) 
N0 = sapply(graph_list, vcount)

########## clean 1) remove name-duplicated Vertex due to inconsistence in STRING.db ###########
## refer to 11.1.0_correct_vertex_duplication.R 
correct_n_edges = readRDS(paste0(inputdir,'correct_n_edges_HiG_STRING2.14.0.rds'))
for(g_name in unique(correct_n_edges$graph_id)){
	vertices_to_remove = subset(correct_n_edges, graph_id == g_name)$vetex_index_to_remove
	if(any(is.na(vertices_to_remove))) vertices_to_remove = vertices_to_remove[!is.na(vertices_to_remove)]
	graph_list[[g_name]] = delete_vertices(graph_list[[g_name]], vertices_to_remove)
}
N = sapply(graph_list, vcount)
(N0-N)[which(N0-N>0)]
# HiG_1 HiG_2 HiG_4 HiG_6 HiG_7 HiG_9 
    # 1     4     1     2     4     1 
range(E(graph_list[[1]])$weight)  #[1] 0.200 0.999

graph_list <- lapply(graph_list, simplify, edge.attr.comb = 'max')  # by default is 'sum'
range(E(graph_list[[1]])$weight)  #[1]  0.200 0.999

names(graph_list)
 # [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"          "HiG_4"           "HiG_5"          
 # [7] "HiG_6"           "HiG_7"           "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"   
# [13] "HiG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"       "HiGCTS_CP.1"     "CTS_muscle"     
# [19] "CTS_endoderm"    "CTS_CP"          "CTS_CP.1"    
edge_counts <- sapply(graph_list, ecount)
edge_counts
          # HiG_0           HiG_1           		HiG_2          	HiG_CP           HiG_4 
          # 10066           31459->31169            7394-> 7129           1306            2330->2269 
          # HiG_5           HiG_6           	HiG_7      			HiG_muscle           HiG_9 
           # 4638             987->894            1712->1587            1989           43172 -> 42812
         # HiG_10    HiG_endoderm          HiG_12   HiGCTS_muscle HiGCTS_endoderm 
           # 1213             141              70              68              11 
      # HiGCTS_CP     HiGCTS_CP.1      CTS_muscle    CTS_endoderm          CTS_CP 
             # 38              41             316             390             128 
       # CTS_CP.1 
            # 159


# scp -p -r E:/Git_Holly/TIPS-main/celltype_specific_weight_v9_xy.R  xyang2@midway3.rcc.uchicago.edu:/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/.
# source('F:/projects/scRNA/source/cardiac_CTS_GRN/celltype_specific_weight_v5.R')
# source('/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/celltype_specific_weight_v5.R')
source('/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/celltype_specific_weight_v9_xy.R')

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
											 assayName,
											 celltype_col = "cluster",
											 method = "pearson", 
											 cores = 4,
											 shrink = TRUE,    
											 verbose = FALSE)
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
    # 19744    85881
	   
	for(s in c("combined","ratio", "zscore", "diff"))
	# s = "combined"
	{										 
		weighted_graph_list <- update_network_weights(graph_list,
									  network_specificity_list,
									  specificity_method = s, #"combined", #"ratio", "zscore", "diff",
									  verbose = FALSE,
									  cores = 4) 	
		saveRDS(weighted_graph_list, file = paste0('GSE175634_STRING_graph_perState_simplified_',s,'weighted.rds'))
	}
  
	# double check outputs
    g <- weighted_graph_list[[1]]
    vertex_attr_names(g) # [1]  "name"   "weight" "FDR"    "size"   "color"
    edge_attr_names(g) #  "combined_score"  "weight"          "width"           "original_weight"
#[5] "corexp_sign"     "coexp_target"


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

	for(net_name in names(graph_list))
		{
		#net_name = "HiGCTS_CP.1"	
		plot_data = NULL
		par(mfrow=c(2,2))
		for(s in c("combined","ratio", "zscore", "diff")) {
			graph_list = readRDS( file = paste0('GSE175634_STRING_graph_perState_simplified_',s,'weighted.rds'))
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
				original_weight = E(g)$original_weight,
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
				 x = "E(g)$original_weight", 
				 y = "log10(E(g)$weight)",
				 fill = "Count") +
			theme(legend.position = "bottom")

		print(p2)
	}
	dev.off()
}								  

net_name = "HiGCTS_CP.1"	
plot_data = NULL

for(s in c("combined","ratio", "zscore", "diff")) {
	graph_list = readRDS( file = paste0('GSE175634_STRING_graph_perState_simplified_',s,'weighted.rds'))
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
		original_weight = E(g)$original_weight,
		log_weight = log10(E(g)$weight),
		method = s,
		is_isl1 = is_isl1_edge
	)

	plot_data <- rbind(plot_data, temp_data)
}

temp = subset(plot_data,is_isl1==TRUE)
pdf(file='compare_specificity_method_HiGCTS_CP.1_vs_PPIscores.pdf')	
ggplot(temp, aes(x = original_weight, y = log_weight, color=as.factor(original_weight))) +
			geom_point() +
			scale_fill_viridis_c() +
			facet_wrap(~method, scales = "free") +
			theme_minimal() +
			labs(title = paste("Isl1 linkages -", net_name),
				 x = "E(g)$original_weight", 
				 y = "log10(E(g)$weight)") +
			theme(legend.position = "bottom")
dev.off()		