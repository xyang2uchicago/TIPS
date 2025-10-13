library('SingleCellExperiment')
library(Seurat)
library(dplyr) 
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)
library(stringr)

wd = '/Users/felixyu/Documents/IbarraSoria2018/'
setwd(paste0(wd, "results/"))
output_dir <- 'PPI_weight/'

PPI_color_platte = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

load(paste0(wd, 'results/sce_16subtype.RData'))

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
graph_list <- readRDS( file= 'IbarraSoria2018_STRING_graph_perState_notsimplified.rds')  
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
source(paste0(wd, 'code/celltype_specific_weight_v9.R'))
source(paste0(wd, 'code/BioTIP_update_06162025.R'))

# Run the below part only once
print('entering into step 1')
step1 = TRUE

if (step1){
# Calculate co-expression based specificity scores (1 hour run)

               
# Get the current order of subcelltypes
subtype_order <- unique(colData(sce)$subcelltype)
colData(sce)$subcelltype <- factor(colData(sce)$subcelltype, levels = subtype_order)
	
network_specificity_list = calculate_network_specificity(sce, 
											 graph_list, 
											 assayName,
											 celltype_col = "subcelltype", #"cluster",
											 method = "pearson", 
											 cores = 4,
											 shrink = TRUE,    
											 verbose = FALSE)
											 
saveRDS(network_specificity_list, file = paste0(output_dir,'network_specificity_list.rds'))

}


# STEP 2
print('entering into step 2')

step2 = TRUE
if(step2){

	network_specificity_list <- readRDS(paste0(output_dir,'network_specificity_list.rds'))

	(names(network_specificity_list[[1]]))
    # [1] "scores"       "genes"        "coexp_target" "corexp_sign" 
	(table(network_specificity_list[[1]]$corexp_sign))
    # negative positive 
    #    51792    82897

	library(data.table)

	for(net in names(network_specificity_list)) {
		cat('Analyzing: ', net)
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
	   
	for(s in c("combined", "ratio", "zscore", "diff"))
	{		
		weighted_graph_list <- update_network_weights(graph_list,
									  network_specificity_list,
									  specificity_method = s,
									  verbose = FALSE,
                                      cores = 1) 	
		saveRDS(weighted_graph_list, file = paste0(output_dir, '2018_STRING_graph_perState_simplified_',s,'weighted.rds'))
	}

}

# STEP 3		
# Compare weight methods
# Check new weights for ISL1 in "cardiac.a"
print('entering into step 3')
step3 = TRUE
library(ggplot2)
library(hexbin)

if(step3){
	pdf(file=paste0(output_dir, 'compare_specificity_method_vs_PPIscores.pdf')) # create this file for comparison

	for(net_name in names(graph_list))
		{
		plot_data = NULL
		par(mfrow=c(2,2))
		for(s in c("combined", "ratio", "zscore", "diff")) {
			weighted_graph_list <- readRDS( file = paste0(output_dir, '2018_STRING_graph_perState_simplified_',s,'weighted.rds'))
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
                    original_weight = E(g)$original_weight,
                    log_weight = log10(E(g)$weight),
                    method = s,
                    is_isl1 = is_isl1_edge
                )

                plot_data <- rbind(plot_data, temp_data)
		}
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

net_name = "HiGCTS_cardiac.a"	# Specific to the 2018 dataset
plot_data = NULL

for(s in c("combined", "ratio", "zscore", "diff")) {
	graph_list <- readRDS( file = paste0(output_dir, '2018_STRING_graph_perState_simplified_',s,'weighted.rds'))
	g = graph_list[[net_name]]
	
	# Check for ISL1 vertex and mark edges
	isl1_vertex <- which(tolower(V(g)$name) == 'isl1')
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
pdf(file=paste0(output_dir,'compare_specificity_method_HiGCTS_cardiac.a_vs_PPIscores.pdf'))	
    p3 <- ggplot(temp, aes(x = original_weight, y = log_weight, color=as.factor(original_weight))) +
                geom_point() +
                scale_fill_viridis_c() +
                facet_wrap(~method, scales = "free") +
                theme_minimal() +
                labs(title = paste("Isl1 linkages -", net_name),
                    x = "E(g)$original_weight", 
                    y = "log10(E(g)$weight)") +
                theme(legend.position = "bottom")
    print(p3)
dev.off()	





















