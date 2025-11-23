library('SingleCellExperiment')
library(Seurat)
library(dplyr) 
library(scuttle)
library(ggplot2)
library(hexbin)

## dependence to run BioTIP
library(igraph)
require(psych)
library(stringr)

########## BEGINNING OF USER INPUT ##########

wd = '/Users/felixyu/Documents/IbarraSoria2018/'
setwd(paste0(wd, "results/PPI_weight/"))

db <- "IbarraSoria2018"

celltype_specific_weight_version <- '10'
BioTIP_version <- '06232025'

specificity_methods <- c("combined", "ratio", "zscore", "diff")

isl1_cluster <- "HiGCTS_cardiac.a" # cluster containing ISL1 gene

core_count <- 8 # number of cores used for parallel processing in steps 1 and 2. Use core_count = 1 if on Windows.

step1 <- FALSE # calculate gene correlations and specificity
step2 <- FALSE # update network edge weights
step3 <- TRUE # graph comparing specificity methods for all clusters and isl1_cluster

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_', BioTIP_version, '.R'))

load(paste0(wd, 'data/sce_16subtype.RData'))

########## END OF USER INPUT ##########

sce
colnames(colData(sce)) # label celltype subcelltype

unique(colData(sce)$celltype)
# [1] extraembryonicMesoderm endothelial            blood                  mesodermProgenitors    presomiticMesoderm    
# [6] somiticMesoderm        mixedMesoderm          pharyngealMesoderm     cardiac               
# 9 Levels: blood cardiac endothelial extraembryonicMesoderm mesodermProgenitors mixedMesoderm ... somiticMesoderm

unique(colData(sce)$subcelltype)
# [1] extraembryonicMesoderm endothelial.a          endothelial.b          endothelial.c          endothelial.d         
# [6] blood                  mesodermProgenitors    presomiticMesoderm.b   presomiticMesoderm.a   somiticMesoderm       
# [11] mixedMesoderm.a        pharyngealMesoderm     mixedMesoderm.b        cardiac.a              cardiac.b             
# [16] cardiac.c
# 16 Levels: blood cardiac.a cardiac.b cardiac.c endothelial.a endothelial.b endothelial.c ... somiticMesoderm

# Calculate log-normalized counts
if (!"logcounts" %in% assayNames(sce)) {
    sce <- scater::logNormCounts(sce)
}
assayName = 'logcounts'



# Read in the graph_list file
graph_list <- readRDS( file= paste0('../',db, '_STRING_graph_perState_notsimplified.rds'))  
N0 = sapply(graph_list, vcount)

# Simplifying graph list (removing duplicate edges)
graph_list <- lapply(graph_list, simplify, edge.attr.comb = 'max')  # by default is 'sum'

edge_counts <- sapply(graph_list, ecount)


graphs_with_duplicates <- sapply(graph_list, function(g) {
  vertex_names <- V(g)$name
  if(is.null(vertex_names)) {
    # If no names, use vertex indices
    vertex_names <- V(g)
  }
  any(duplicated(vertex_names))
})

# See which graphs have duplicates
(which(graphs_with_duplicates))
# named integer(0)


# STEP 1

# Run the below part only once
if (step1){
print('entering into step 1')
# Calculate co-expression based specificity scores
# Get the current order of subcelltypes
subtype_order <- unique(colData(sce)$subcelltype)
colData(sce)$subcelltype <- factor(colData(sce)$subcelltype, levels = subtype_order)
	
network_specificity_list = calculate_network_specificity(sce, 
											 graph_list, 
											 assayName,
											 celltype_col = "subcelltype", #"cluster",
											 method = "pearson", 
											 cores = core_count,
											 shrink = TRUE,    
											 verbose = FALSE)
											 
saveRDS(network_specificity_list, file = paste0('network_specificity_list.rds'))

}


# STEP 2
if(step2){
    print('entering into step 2')
	network_specificity_list <- readRDS(paste0('network_specificity_list.rds'))

	(names(network_specificity_list[[1]]))
    # [1] "scores"       "genes"        "coexp_target" "corexp_sign" 
	(table(network_specificity_list[[1]]$corexp_sign))
    # negative positive 
    #    51792    82897

	library(data.table)

	for(net in names(network_specificity_list)) {
		cat('Analyzing: ', net, '\n')
		spec_data <- network_specificity_list[[net]]
		corexp_named <- spec_data$corexp_sign

		stopifnot(
			  is.matrix(corexp_named),
			  all(dim(corexp_named) == dim(spec_data$coexp_target))
			  )
        
		network_specificity_list[[net]]$corexp_sign <- corexp_named
	}


	(names(network_specificity_list[[1]]))
	#[1] "scores"      "genes"       "corexp_sign"  "coexp_target"
	(table(network_specificity_list[[1]]$corexp_sign))
    # negative positive 
    #    51792    82897 
	   
	for(s in specificity_methods)
	{		
		weighted_graph_list <- update_network_weights(graph_list,
									  network_specificity_list,
									  specificity_method = s,
									  verbose = FALSE,
                                      cores = core_count) 	
		saveRDS(weighted_graph_list, file = paste0(db, '_STRING_graph_perState_simplified_',s,'weighted.rds'))
	}

}

# STEP 3		
# Compare weight methods
# Check new weights for ISL1 in "cardiac.a"
if(step3){
    print('entering into step 3')
	pdf(file=paste0('compare_specificity_method_vs_PPIscores.pdf')) # create this file for comparison

	for(net_name in names(graph_list))
		{
		plot_data = NULL
		par(mfrow=c(2,2))
		for(s in specificity_methods) {
			weighted_graph_list <- readRDS( file = paste0(db, '_STRING_graph_perState_simplified_',s,'weighted.rds'))
			g = weighted_graph_list[[net_name]]
			
			# Check for ISL1 vertex and mark edges
                isl1_vertex <- which(tolower(V(g)$name) == "isl1")
                is_isl1_edge <- rep(FALSE, length(E(g)))
                if(length(isl1_vertex) > 0) {
                    isl1_edges <- incident(g, isl1_vertex, mode = "all")
                    is_isl1_edge[isl1_edges] <- TRUE
                }
                
                temp_data <- data.frame(
                    net_name=net_name,
                    norm_PPI_score = E(g)$norm_PPI_score,
                    log_weight = log10(E(g)$weight),
                    method = s,
                    is_isl1 = is_isl1_edge
                )

                plot_data <- rbind(plot_data, temp_data)
		}
		p2 <- ggplot(plot_data, aes(x = norm_PPI_score, y = log_weight)) +
			stat_binhex(bins = 50, alpha=0.7) +
			geom_point(data = subset(plot_data, is_isl1 == TRUE), 
					   aes(x = norm_PPI_score, y = log_weight), 
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
	}
	dev.off()
}

net_name = isl1_cluster	# Specific to the 2018 dataset
plot_data = NULL

for(s in specificity_methods) {
	graph_list <- readRDS( file = paste0(db, '_STRING_graph_perState_simplified_',s,'weighted.rds'))
	g = graph_list[[net_name]]
	
	# Check for ISL1 vertex and mark edges
	isl1_vertex <- which(tolower(V(g)$name) == 'isl1')
	is_isl1_edge <- rep(FALSE, length(E(g)))
	if(length(isl1_vertex) > 0) {
		isl1_edges <- incident(g, isl1_vertex, mode = "all")
		is_isl1_edge[isl1_edges] <- TRUE
	}
    else{
        print(paste0("No ISL1 found in ", isl1_cluster))
        break
    }
	
	temp_data <- data.frame(
		net_name=net_name,
		norm_PPI_score = E(g)$norm_PPI_score,
		log_weight = log10(E(g)$weight),
		method = s,
		is_isl1 = is_isl1_edge
	)

	plot_data <- rbind(plot_data, temp_data)
}

if(length(isl1_vertex) > 0){
    temp = subset(plot_data,is_isl1==TRUE)
    pdf(file=paste0('compare_specificity_method_', isl1_cluster,'_vs_PPIscores.pdf'))	
        p3 <- ggplot(temp, aes(x = norm_PPI_score, y = log_weight, color=as.factor(norm_PPI_score))) +
                    geom_point() +
                    scale_fill_viridis_c() +
                    facet_wrap(~method, scales = "free") +
                    theme_minimal() +
                    labs(title = paste("Isl1 linkages -", net_name),
                        x = "E(g)$norm_PPI_score", 
                        y = "log10(E(g)$weight)") +
                    theme(legend.position = "bottom")
        print(p3)
    dev.off()	
}





















