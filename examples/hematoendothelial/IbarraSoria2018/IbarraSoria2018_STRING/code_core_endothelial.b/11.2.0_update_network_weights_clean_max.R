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

code_dir <- here::here("examples", "hematoendothelial", "IbarraSoria2018", "IbarraSoria2018_STRING", "code_core_endothelial.b")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(ppi_path)

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_', BioTIP_version, '.R'))

load(paste0(data_dir, 'sce_16subtype.RData'))

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
    assayName <- assayNames(sce)[1]
    x <- assay(sce, assayName)
    if (max(x) > 20) {
        sce <- scater::logNormCounts(sce, assay.type = assayName)
    } else {
        assayNames(sce)[1] <- 'logcounts'
    }
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
    # [1] "scores"       "genes"        "coexp_focal" "corexp_sign" 
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
			  all(dim(corexp_named) == dim(spec_data$coexp_focal))
			  )
        
		network_specificity_list[[net]]$corexp_sign <- corexp_named
	}


	(names(network_specificity_list[[1]]))
	#[1] "scores"      "genes"       "corexp_sign"  "coexp_focal"
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






















