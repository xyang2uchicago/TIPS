 
 #source('F:/projects/scRNA/source/cardiac_CTS_GRN/celltype_specific_weight_v6.R')
#  source('/project/xyang2/heart_dev/source_midway3/GSE175634_iPSC_CM/celltype_specific_weight_v6.R')
# library(BioTIP)  # or load Felix's last version of BioTIP_xxxx.R 
#  source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP_update_06162025.R')

source('/Users/felixyu/Documents/GSE87038_weighted/code/celltype_specific_weight_v7.R')


############ calculate_specificity ############ 
 # Create example correlation matrices
 genes <- c("GENE1", "GENE2", "GENE3", "GENE4")
 target_matrix <- matrix(runif(16), 4, 4)
 rownames(target_matrix) <- colnames(target_matrix) <- genes
 
 other_matrix1 <- matrix(runif(16), 4, 4)
 rownames(other_matrix1) <- colnames(other_matrix1) <- genes
 
 other_matrix2 <- matrix(runif(16), 4, 4)
 rownames(other_matrix2) <- colnames(other_matrix2) <- genes
 
 coexp_other_list = list(other_matrix1, other_matrix2)
 for(other_coexp in coexp_other_list) {
    cat(dim(other_coexp), '\n')
 }
 
 # Calculate specificity
 spec_scores <- calculate_specificity(
   target_matrix, 
   coexp_other_list
 )
 
 ############# compute_cluster_correlation ############
 library(SingleCellExperiment)
 library(scuttle)
 set.seed(123)
 counts <- matrix(rpois(5000, lambda = 10), nrow = 100, ncol = 50, 
                 dimnames = list(paste0("GENE", 1:100), paste0("Cell", 1:50)))
 coldata <- data.frame(cluster = sample(c("A", "B", "C"), 50, replace = TRUE))
 toy_sce <- SingleCellExperiment(assays = list(counts = counts), colData = coldata)
 toy_sce <- logNormCounts(toy_sce)

 # Get co-expression matrix for a set of genes in cluster "1"
 genes_of_interest <- c("GENE1", "GENE2", "GENE3", "GENE4")
 
 coexp_matrix <- compute_cluster_correlation(
   sce = toy_sce,
   cluster_id = "A",
   genes = genes_of_interest,
   celltype_col = "cluster",
   assayName = "logcounts"
 )
 
 
 ############ calculate_network_specificity ############
 set.seed(123)
 counts <- matrix(rpois(5000, lambda = 10), nrow = 100, ncol = 50, 
                 dimnames = list(paste0("GENE", 1:100), paste0("Cell", 1:50)))
 coldata <- data.frame(cluster = sample(c("A", "B", "C"), 50, replace = TRUE))
 toy_sce <- SingleCellExperiment(assays = list(counts = counts), colData = coldata)
 toy_sce <- logNormCounts(toy_sce)

 gene_names <- paste0("GENE", 1:10)
 toy_network <- erdos.renyi.game(10, p = 0.3) %>% set_vertex_attr("name", value = gene_names)

 genes_of_interest <- c("GENE1", "GENE2", "GENE3", "GENE4")

 specificity_scores_shrink <- calculate_network_specificity(
   sce = toy_sce,
   graph_list = list('A' = toy_network),
   celltype_col = "cluster",
   cores = 1
 )
  
  specificity_scores <- calculate_network_specificity(
   sce = toy_sce,
   graph_list = list('A' = toy_network),
   celltype_col = "cluster",
   shrink = FALSE,
   cores = 1
 )
  
  plot(specificity_scores_shrink[['A']]$scores$combined,
       specificity_scores[['A']]$scores$combined)
	   
	
############### update_network_weights() ##########
# generate a toy network with edge weights, with a cell type naemd 'A'
 set.seed(123)
 counts <- matrix(rpois(5000, lambda = 10), nrow = 100, ncol = 50, 
                 dimnames = list(paste0("GENE", 1:100), paste0("Cell", 1:50)))
 coldata <- data.frame(cluster = sample(c("A", "B", "C"), 50, replace = TRUE))
 toy_sce <- SingleCellExperiment(assays = list(counts = counts), colData = coldata)
 toy_sce <- logNormCounts(toy_sce)
 gene_names <- paste0("GENE", 1:10)
 toy_network <- erdos.renyi.game(10, p = 0.3) %>% set_vertex_attr("name", value = gene_names)
 E(toy_network)$weight = sample(seq(0.2, 1, 0.01), ecount(toy_network))
 edge_attr_names(toy_network) 
# [1] "weight"

 #  calculate specific_scores for celltype A
 specificity_scores_shrink <- calculate_network_specificity(
   sce = toy_sce,
   graph_list = list('A' = toy_network),
   celltype_col = "cluster",
   cores = 1
)
 # update teh network edge feature by adding the specificity_scores 
  updated_networks <- update_network_weights(
    list('A' = toy_network), 
    specificity_scores_shrink,
    specificity_method = "combined"
  )
  edge_attr_names(updated_networks[['A']] ) 
# [1] "weight"          "original_weight" "corexp_sign"     "coexp_target"  

  plot(E(updated_networks[['A']])$coexp_target, E(updated_networks[['A']])$weight,
	xlab = 'Pearson cor',  ylab ='Coexp&PPI-combined Weights' , main= 'Cell cluster A-specific' )

 
 
 ####### robustness_MonteCarlo ###### 
 
 gene_names <- c("ISL1", "FGF10", "MEIS1", "HAPLN1", "NTRK2", "DUSP6", "HAS2", 
                "H1F0", "HAND2", "BMP5", "ID1", "CITED2", "BMPER", "WLS", 
                "NKX3-1", "LAMA1", "LRRTM1", "PTPN13", "IFI16", "SLC7A2",
                "GENE21", "GENE22", "GENE23", "GENE24", "GENE25", "GENE26", "GENE27")

# Generate random edges (41 edges among 27 genes)
set.seed(123)  # For reproducible results
edge_list <- t(replicate(41, sample(gene_names, 2)))

# Create network and add weights
toy_network <- graph_from_edgelist(edge_list, directed = FALSE)
E(toy_network)$weight <- runif(ecount(toy_network), 0.2, 1.0)

 res1 = robustness_MonteCarlo(toy_network, "vertex", "random", N=10)
 res2 = robustness_MonteCarlo(toy_network, "edge", "random", N=10)
 res3 = robustness_MonteCarlo(toy_network, "vertex", "btwn.cent")
 res4 = robustness_MonteCarlo(toy_network, "edge", "btwn.cent")
 dim(res1) #[1] 27  5
 dim(res2) #[1] 42  5
 dim(res3) #[1] 27  5
 dim(res4)  #[1] 42  5
 