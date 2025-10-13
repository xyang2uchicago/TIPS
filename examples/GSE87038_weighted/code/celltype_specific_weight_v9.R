############ load libraary for testing #########
library(SingleCellExperiment)
library(scater)
library(Matrix)
library(qlcMatrix)
library(igraph)
library(parallel)
# Auto-detect based on file system
if (file.exists("/project/imoskowitz/")) {
    # We're on midway3
    #   source('/project/imoskowitz/xyang2/heart_dev/GSE175634_iPSC_CM/BioTIP_update_06162025.R')
    source("/project/xyang2/felixy/GSE87038_weighted/code/BioTIP_update_06162025.R")
} else {
    # We're on local machine
    source("/Users/felixyu/Documents/GSE87038_weighted/code/BioTIP_update_06162025.R")
}

################  The above lines should be removed when submission to Bioconductor !!!###############################

## version 9 removed invert_weights in update_network_weights()
## version 8 replace NA with 0 in compute_cluster_correlation()
## version 8 added min_n_Eg to calculate_network_specificity()
## version 7 removed sparse matrix handling and added parallelization and vectorization of update_network_weights()
## version 6 call the new Biotip::cor.shrink()  function which can handle sparser matrix
## version 5 considers BioTIP::cor.shrink() into compute_cluster_correlation() by adding a parameter <shrink>
## version 4 add <corexp_sign> as an output to calculate_network_specificity(), calculate_network_specificity()

#' Calculate Cell-Type Specificity Scores for Gene Co-expression
#'
#' @description
#' Calculates several measures of cell-type specificity by comparing
#' co-expression of genes in a target cell type/cluster to their co-expression
#' in other cell types/clusters.
#'
#' @param coexp_target A numeric matrix of gene-gene correlations (co-expression)
#'   in the target cell type. Row and column names should be gene identifiers.
#' @param coexp_other_list A list of numeric matrices with gene-gene correlations
#'   in other cell types/clusters. Each matrix should have the same dimensions
#'   and gene ordering as \code{coexp_target}.
#'
#' @return A list containing four different specificity score matrices:
#'   \describe{
#'     \item{ratio}{Ratio of target co-expression to mean co-expression in other cell types}
#'     \item{zscore}{Z-score measure comparing target to other cell types}
#'     \item{diff}{Simple difference between target and mean of other cell types}
#'     \item{combined}{Product of target co-expression and (1 - other co-expression)}
#'   }
#'
#' @details
#' This function calculates multiple specificity metrics that highlight gene-gene
#' connections that are specific to a target cell type:
#'
#' 1. Ratio: Simple ratio between target and other cell types, high values indicate
#'    more specific connections.
#' 2. Z-score: How many standard deviations the target co-expression is above the
#'    mean co-expression in other cell types.
#' 3. Difference: Simple subtraction to find connections stronger in target.
#' 4. Combined: A metric that rewards both high co-expression in target and low
#'    co-expression in other cell types.
#'
#' Absolute values of correlations are used to focus on connection strength rather
#' than direction (positive or negative correlation).
#'
#' @examples
#' \dontrun{
#' # Create example correlation matrices
#' genes <- c("GENE1", "GENE2", "GENE3", "GENE4")
#' target_matrix <- matrix(runif(16), 4, 4)
#' rownames(target_matrix) <- colnames(target_matrix) <- genes
#'
#' other_matrix1 <- matrix(runif(16), 4, 4)
#' rownames(other_matrix1) <- colnames(other_matrix1) <- genes
#'
#' other_matrix2 <- matrix(runif(16), 4, 4)
#' rownames(other_matrix2) <- colnames(other_matrix2) <- genes
#'
#' # Calculate specificity
#' spec_scores <- calculate_specificity(
#'     target_matrix,
#'     list(other_matrix1, other_matrix2)
#' )
#' }
#'
#' @export
#' @author Holly Yang \email{xyang2.at.uchicago.edu}

calculate_specificity <- function(coexp_target, coexp_other_list) {
    # Convert correlation to absolute value
    coexp_target_abs <- abs(coexp_target)

    # Calculate mean and standard deviation of co-expression in other cell types
    other_coexp_mean <- matrix(0, nrow = nrow(coexp_target), ncol = ncol(coexp_target))
    other_coexp_sd <- matrix(0, nrow = nrow(coexp_target), ncol = ncol(coexp_target))
    dimnames(other_coexp_mean) <- dimnames(coexp_target)
    dimnames(other_coexp_sd) <- dimnames(coexp_target)

    for (other_coexp in coexp_other_list) {
        other_coexp_abs <- abs(other_coexp)
        other_coexp_mean <- other_coexp_mean + other_coexp_abs
    }

    # Check if there are other cell types to compare with
    if (length(coexp_other_list) > 0) {
        other_coexp_mean <- other_coexp_mean / length(coexp_other_list)

        # Calculate standard deviation if multiple cell types
        if (length(coexp_other_list) > 1) {
            for (other_coexp in coexp_other_list) {
                other_coexp_abs <- abs(other_coexp)
                other_coexp_sd <- other_coexp_sd + (other_coexp_abs - other_coexp_mean)^2
            }
            other_coexp_sd <- sqrt(other_coexp_sd / (length(coexp_other_list) - 1))
        } else {
            # Default small value to avoid division by zero
            other_coexp_sd[] <- 0.1
        }
    } else {
        # If no other cell types, set to zeros (maximum specificity)
        other_coexp_mean[] <- 0
        other_coexp_sd[] <- 0.1
    }

    # Calculate various specificity measures

    # Option 1: Simple ratio (avoid division by zero)
    specificity_ratio <- coexp_target_abs / (other_coexp_mean + 0.01)

    # Option 2: Z-score like measure
    specificity_zscore <- (coexp_target_abs - other_coexp_mean) / (other_coexp_sd + 0.01)

    # Option 3: Difference-based measure
    specificity_diff <- coexp_target_abs - other_coexp_mean

    # Option 4: Combined measure (rewards high co-expression in target and low in others)
    specificity_combined <- coexp_target_abs * (1 - other_coexp_mean)

    return(list(
        ratio = specificity_ratio,
        zscore = specificity_zscore,
        diff = specificity_diff,
        combined = specificity_combined
    ))
}

#' Compute Gene-Gene Correlation Matrix for a Specific Cell Type/Cluster
#'
#' @description
#' Computes a gene-gene correlation matrix (co-expression matrix) for a set of
#' genes in a specific cell type or cluster from a SingleCellExperiment object.
#'
#' @param sce A SingleCellExperiment object containing gene expression data.
#' @param cluster_id Character or numeric identifier for the target cluster
#'   or cell type.
#' @param genes Character vector of gene identifiers to include in the
#'   correlation calculation.
#' @param celltype_col Character string specifying the column in sce's colData that
#'   contains the cell type/cluster identifiers.
#' @param assayName Character string specifying which assay to use from the SCE object.
#'   Options are "logcounts" (default), "counts", or other assay names.
#' @param method Character string specifying the correlation method to use.
#'   Options are "pearson" (default), "spearman", or "kendall".
#'   When <shrink> is TRUE, we apply "pearson" and this parameter is ignored.
#' @param shrink A flag specifying whether to shrink the matrix of gene-gene correlation or not.
#'   If TRUE, call BioTIP::cor.shrink() else call cor() function.
#'   This appraoch uses the method outlined by Schafer and Strimmer in
#'   "A Shrinkage Approach to Large-Scale Covariance Matrix Estimation
#'   and Implications for Functional Genomics" (2005). Here, we shrink between-gene correlations
#'   towards 0 due to the low global gene expressional dependence in a stable state
#'   Comparing to fun = 'cor', the 'BioTIP' method without shinkage is modified
#'   to ignore missing values, analogous to how \code{cor(X, use = "pairwise.complete.obs")} works.
#'   For between-sample correlation matrix, we shrink
#'   towards the average correlation to reflect the similar gene-expression profiles in a stable state.
#' @param min_cells Integer specifying the minimum number of cells required
#'   for a reliable correlation. A warning is issued if fewer cells are available.
#'
#' @return A symmetric numeric matrix of gene-gene correlations with row and
#'   column names matching the input gene identifiers.
#'
#' @details
#' This function extracts expression data for specified genes from a particular
#' cell cluster, then calculates pairwise correlations between all genes.
#' The resulting matrix can be used for co-expression network analysis, gene
#' module detection, or calculation of cell-type specificity scores.
#'
#' @examples
#' { #
#'     library(SingleCellExperiment)
#'     library(scuttle)
#'     set.seed(123)
#'     counts <- matrix(rpois(5000, lambda = 10),
#'         nrow = 100, ncol = 50,
#'         dimnames = list(paste0("GENE", 1:100), paste0("Cell", 1:50))
#'     )
#'     coldata <- data.frame(cluster = sample(c("A", "B", "C"), 50, replace = TRUE))
#'     toy_sce <- SingleCellExperiment(assays = list(counts = counts), colData = coldata)
#'     toy_sce <- logNormCounts(toy_sce)
#'
#'     # Get co-expression matrix for a set of genes in cluster "1"
#'     genes_of_interest <- c("GENE1", "GENE2", "GENE3", "GENE4")
#'     coexp_matrix <- compute_cluster_correlation(
#'         sce = sce_object,
#'         cluster_id = "A",
#'         genes = genes_of_interest,
#'         celltype_col = "cluster",
#'         assayName = "logcounts"
#'     )
#' }
#'
#' @importFrom SingleCellExperiment logcounts counts
#' @export
#' @seealso \code{\link{cor.shrink}}
#' @author Holly Yang \email{xyang2.at.uchicago.edu}
#'
compute_cluster_correlation <- function(sce,
                                        cluster_id,
                                        genes,
                                        celltype_col = "cluster",
                                        assayName = c("logcounts", "counts"),
                                        method = "pearson",
                                        min_cells = 10,
                                        shrink = TRUE) { #**  new
    require(Matrix)
    require(qlcMatrix)
    # Match assay name argument
    assayName <- match.arg(assayName)

    # Subset SCE to specific cluster
    sce_cluster <- sce[, sce[[celltype_col]] == cluster_id]

    # Check if we have enough cells
    if (ncol(sce_cluster) < min_cells) {
        warning("Cluster ", cluster_id, " has only ", ncol(sce_cluster), " cells. Correlations may be unreliable.")
    }

    # Get normalized expression matrix for these genes
    if (assayName == "logcounts") {
        expr_mat <- logcounts(sce_cluster[genes, ])
    } else if (assayName == "counts") {
        expr_mat <- counts(sce_cluster[genes, ])
    } else {
        expr_mat <- log1p(counts(sce_cluster[genes, ]))
    }

    # Calculate correlation matrix
    # Added v7 convert any sparse matrix into dense matrix for processing
    if (shrink) {
        expr_mat <- as.matrix(expr_mat)
        cor_mat <- cor.shrink(expr_mat, shrink = TRUE, target = 0, MARGIN = 1)
    } else {
        expr_mat <- t(as.matrix(expr_mat))
        cor_mat <- stats::cor(expr_mat, method = method)
    }

    cor_mat[is.na(cor_mat)] <- 0  # v8 added
    return(cor_mat)
}


#' Calculate Cell-Type Specificity Scores for Gene Networks
#'
#' @description
#' This function computes cell-type specificity scores for genes in each network.
#' It calculates co-expression within the target cell type and compares it to
#' co-expression in other cell types to identify connections that are specific
#' to the target cell type.
#'
#' @param sce A SingleCellExperiment object containing gene expression data.
#' @param graph_list A named list of igraph objects, where each network corresponds
#'   to a specific cell type indicated in its name (e.g., 'HiG_0', 'CTS_CP').
#' @param assayName Character string specifying which assay to use from the SCE object.
#'   Options are "logcounts" (default), "counts", or other assay names.
#' @param celltype_col Character string specifying the column in sce's colData that
#'   contains the cell type/cluster identifiers.
#' @param method Character string specifying the correlation method to use.
#'   Options are "pearson" (default), "spearman", or "kendall".
#'    When <shrink> was TRUE, "pearson"was used and this parameter is ignored.
#' @param cores Integer specifying the number of cores to use for parallel processing.
#' @param verbose Logical indicating whether to print detailed progress messages.
#' @param shrink A flag specifying whether to shrink the matrix of gene-gene correlation or not.
#'   If TRUE, call BioTIP::cor.shrink() else call cor() function.
#' @param min_n_Vg An integer to set the minimal number of nodes in the graph to calculate the pecificity score.
#' @return A nested list with specificity scores for each network. Each network entry contains
#'   several specificity measures (ratio, zscore, diff, combined).
#'
#' @details
#' The function extracts the appropriate cell type identifier from each network name,
#' computes gene co-expression within that cell type, and compares it to co-expression
#' in all other cell types to generate several specificity metrics.
#'
#' @examples
#' \dontrun{
#' # Calculate specificity scores
#' specificity_scores <- calculate_network_specificity(
#'     sce = sce_object,
#'     graph_list = my_networks,
#'     celltype_col = "cluster"
#' )
#' }
#'
#' @importFrom SingleCellExperiment logcounts counts
#' @importFrom igraph V ecount as_edgelist
#' @importFrom parallel mclapply
#' @export
#' @seealso \code{\link{compute_cluster_correlation}}
#' @author Holly Yang \email{xyang2.at.uchicago.edu}
#'
#'

calculate_network_specificity <- function(sce,
                                          graph_list,
                                          assayName = c("logcounts", "counts"),
                                          celltype_col = "cluster",
                                          method = "pearson",
                                          cores = 4,
                                          shrink = TRUE,
                                          min_n_Vg = 5,
                                          min_n_Eg = 10,
                                          verbose = FALSE) {
    # Get all cluster IDs present in the SCE
    all_clusters <- unique(sce[[celltype_col]])
    cat("Available clusters in SCE:", paste(all_clusters, collapse = ", "), "\n")

    # Match assay name argument
    assayName <- match.arg(assayName)

    # Extract cluster ID from network name
    get_cluster_id <- function(network_name) {
        # Extract cluster ID after the underscore
        cluster_id <- sub("^[^_]*_", "", network_name)
        return(cluster_id)
    }

    # Initialize results list
    specificity_scores <- list()

    # Process each network
    for (net_name in names(graph_list)) {
        cat("\nProcessing network:", net_name, "\n")

        # Get the graph
        g <- graph_list[[net_name]]

        # Skip very small networks if needed
      if(ecount(g) < min_n_Eg) {  # v8 changed
            cat("Skipping small network:", net_name, "with only", ecount(g), "edges\n")
            next
        }

        # Extract cluster ID from network name
        target_cluster <- get_cluster_id(net_name)
        target_cluster <- gsub("\\.\\d+$", "", target_cluster)
        if (verbose) cat("Target cluster for this network:", target_cluster, "\n")

        # Check if this cluster exists in the SCE
        if (!target_cluster %in% all_clusters) {
            warning("Cluster '", target_cluster, "' not found in SCE for network '", net_name, "'. Skipping.")
            next
        }

        # Define other clusters as all except the target
        other_clusters <- setdiff(all_clusters, target_cluster)

        # Get gene names from this network
        net_genes <- V(g)$name
        if (verbose) cat("Network has", length(net_genes), "genes\n")

        # Check which genes exist in the SCE
        genes_in_sce <- intersect(rownames(sce), net_genes)
        if (verbose) cat("Found", length(genes_in_sce), "genes in SCE\n")

        if (length(genes_in_sce) < min_n_Vg) {
            warning("Too few genes found in SCE for network '", net_name, "'. Skipping.")
            next
        }

        # Calculate co-expression for target cluster using the independent function
        if (verbose) cat("Computing co-expression for target cluster:", target_cluster, "\n")
        coexp_target <- compute_cluster_correlation(
            sce = sce,
            cluster_id = target_cluster,
            genes = genes_in_sce,
            celltype_col = celltype_col,
            assayName = assayName,
            method = method,
            shrink = shrink
        )

        # Calculate co-expression for other clusters
        if (verbose) cat("Computing co-expression for", length(other_clusters), "other clusters\n")

        # Use parallel processing if multiple other clusters
        if (cores > 1 && length(other_clusters) > 1) {
            coexp_other_list <- mclapply(other_clusters, function(clust) {
                compute_cluster_correlation(
                    sce = sce,
                    cluster_id = clust,
                    genes = genes_in_sce,
                    celltype_col = celltype_col,
                    assayName = assayName,
                    method = method,
                    shrink = shrink # added v7
                )
            }, mc.cores = min(cores, length(other_clusters)))
        } else {
            coexp_other_list <- lapply(other_clusters, function(clust) {
                compute_cluster_correlation(
                    sce = sce,
                    cluster_id = clust,
                    genes = genes_in_sce,
                    celltype_col = celltype_col,
                    assayName = assayName,
                    method = method,
                    shrink = shrink # added v7
                )
            })
        }

        names(coexp_other_list) <- other_clusters

        # Calculate specificity scores using the independent function

        scores <- calculate_specificity(coexp_target, coexp_other_list)

        # Store results
        specificity_scores[[net_name]] <- list(
            scores = scores,
            genes = genes_in_sce, # Store genes for which scores were calculated
            coexp_target = coexp_target,
            corexp_sign = matrix(
                ifelse(coexp_target > 0, "positive", "negative"),
                nrow = nrow(coexp_target),
                ncol = ncol(coexp_target),
                dimnames = dimnames(coexp_target)
            ) # added v7
        )

        cat("Completed calculating specificity for network:", net_name, "\n")
    }


    return(specificity_scores = specificity_scores)
}

#' Update Graph Edge Weights Based on Cell-Type Specificity
#'
#' @description
#' This function updates edge weights by combining STRING protein interaction
#' confidence scores with cell-type specific co-expression data. The resulting
#' weights reflect both interaction confidence and cell-type specificity.
#' Requires specificity scores from calculate_network_specificity() that include
#' both numeric scores and co-expression sign information.
#'
#' @param graph_list A named list of PPIN objects to update.
#'   Within each PPIN g , the E(g)$weight was the STRING_combined_score / 1000.
#' @param specificity_scores A list of specificity scores as returned by
#'   \code{calculate_network_specificity}.
#' @param specificity_method Character string specifying which specificity measure to use.
#'   Options are "combined" (default), "ratio", "zscore", or "diff".
#' @param max_min_norm Logical indicating whether to apply max-min normalization
#'   to both STRING and specificity matrices before combining. Defaults to FALSE.
#' @param verbose Logical indicating whether to print detailed progress messages.
#'
#' @param cores Integer specifying the number of cores to use for parallel processing.
#'
#' @return A list of igraph objects with the following updated edge features, where n_e is the number of edges.
#' \itemize{
#'   \item \code{...}: Besides 'weight', other edge features initiated from the input igraph object.
#'   \item \code{weight}: Final combined edge weights (PPI * specificity).
#'   \item \code{original_weight}:  Original PPI scores (STRING combined_score/1000)
#'   \item \code{corexp_sign}: Sign of co-expression correlation ('positive' or 'negative')
#'   \item \code{coexp_target}: Cell-type specific co-expression scores used for weighting
#' }
#'
#' @details
#' The function combines the original edge weights with cell-type specificity scores
#' to produce weights that favor edges with both high confidence protein interactions
#' and cell-type specific co-expression patterns. By default, weights are inverted
#' (1/weight) so that higher combined scores result in shorter distances for path-based
#' algorithms.
#'
#' @examples
#' {
#'     set.seed(123)
#'     counts <- matrix(rpois(5000, lambda = 10),
#'         nrow = 100, ncol = 50,
#'         dimnames = list(paste0("GENE", 1:100), paste0("Cell", 1:50))
#'     )
#'     coldata <- data.frame(cluster = sample(c("A", "B", "C"), 50, replace = TRUE))
#'     toy_sce <- SingleCellExperiment(assays = list(counts = counts), colData = coldata)
#'     toy_sce <- logNormCounts(toy_sce)
#'     gene_names <- paste0("GENE", 1:10)
#'     toy_network <- erdos.renyi.game(10, p = 0.3) %>% set_vertex_attr("name", value = gene_names)
#'     E(toy_network)$weight <- sample(seq(0.2, 1, 0.01), ecount(toy_network))
#'     edge_attr_names(toy_network)
#'     # [1] "weight"
#'
#'     #  calculate specific_scores for celltype A
#'     specificity_scores_shrink <- calculate_network_specificity(
#'         sce = toy_sce,
#'         graph_list = list("A" = toy_network),
#'         celltype_col = "cluster",
#'         cores = 1
#'     )
#'     # update the network edge feature by adding the specificity_scores
#'     updated_networks <- update_network_weights(
#'         list("A" = toy_network),
#'         specificity_scores_shrink,
#'         specificity_method = "combined"
#'     )
#'     edge_attr_names(updated_networks[["A"]])
#'     # [1] "weight"          "original_weight" "corexp_sign"     "coexp_target"
#'
#'     plot(E(updated_networks[["A"]])$coexp_target, E(updated_networks[["A"]])$weight,
#'         xlab = "Pearson cor", ylab = "Coexp&PPI-combined Weights", main = "Cell cluster A-specific"
#'     )
#' }
#'
#' @importFrom igraph V E ecount ends
#' @export
#' @author Holly Yang \email{xyang2.at.uchicago.edu}
#'
update_network_weights <- function(graph_list,
                                   specificity_scores,
                                   specificity_method = c("combined", "ratio", "zscore", "diff"),
                                   max_min_norm = FALSE,
                                   verbose = FALSE,
                                   cores = 1) {
    # Match specificity method argument
    specificity_method <- match.arg(specificity_method)

    # max-min normalization function
    normalize <- function(x) {
        if (all(x == 0)) {
            return(x)
        } # Handle all-zero matrix
        x_norm <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        x_norm[is.na(x_norm)] <- 0
        return(x_norm)
    }

    update_graph_weights <- function(net_name) {
        if (verbose) cat("\nUpdating weights for network:", net_name, "\n")

        if (!net_name %in% names(graph_list)) {
            warning("Network '", net_name, "' not found in graph_list. Skipping.")
            return(NULL)
        }

        # Get the graph
        g <- graph_list[[net_name]]
        E(g)$original_weight <- E(g)$weight

        # Get specificity data
        spec_data <- specificity_scores[[net_name]]
        genes_in_sce <- spec_data$genes
        chosen_specificity <- spec_data$scores[[specificity_method]]

        if (is.null(chosen_specificity)) {
            warning(
                "Specificity method '", specificity_method,
                "' not found for network '", net_name, "'. Skipping."
            )
            return(g)
        }

        # Get STRING scores
        edge_list <- igraph::as_edgelist(g, names = TRUE)
        string_scores <- E(g)$weight

        # Create STRING score matrix
        string_matrix <- matrix(0,
            nrow = length(genes_in_sce), ncol = length(genes_in_sce),
            dimnames = list(genes_in_sce, genes_in_sce)
        )
        for (i in seq_len(nrow(edge_list))) {
            gene1 <- edge_list[i, 1]
            gene2 <- edge_list[i, 2]
            if (gene1 %in% genes_in_sce && gene2 %in% genes_in_sce) {
                string_matrix[gene1, gene2] <- string_scores[i]
                string_matrix[gene2, gene1] <- string_scores[i]
            }
        }

        # Normalize if requested
        if (max_min_norm) {
            string_norm <- normalize(string_matrix)
            specificity_norm <- normalize(chosen_specificity)
        } else {
            string_norm <- string_matrix
            specificity_norm <- chosen_specificity
        }

        # Edge endpoints
        EE <- igraph::ends(g, es = E(g), names = TRUE)

        # Valid edges
        valid <- EE[, 1] %in% genes_in_sce & EE[, 2] %in% genes_in_sce
        idx <- cbind(EE[valid, 1], EE[valid, 2])

        # Preallocate (all NA)
        nE <- igraph::ecount(g)
        corexp_vals <- rep(NA_character_, nE)
        coexp_target_vals <- rep(NA_real_, nE)
        new_weights <- E(g)$original_weight # default: keep old weights

        # Fill for valid edges
        corexp_vals[valid] <- as.character(spec_data$corexp_sign[idx])
        coexp_target_vals[valid] <- specificity_norm[idx]

        # Compute combined weights
        candidate_weights <- (string_norm[idx] * specificity_norm[idx])

        # Only accept positive weights
        keep <- which(valid)[candidate_weights > 0 & !is.na(candidate_weights)]

        new_weights[keep] <- candidate_weights[candidate_weights > 0 & !is.na(candidate_weights)] # v9 changed

        # Assign in bulk
        g <- set_edge_attr(g, "corexp_sign", value = corexp_vals)
        g <- set_edge_attr(g, "coexp_target", value = coexp_target_vals)
        g <- set_edge_attr(g, "weight", value = new_weights)

        if (verbose) cat("Completed updating weights for network:", net_name, "\n")
        return(g)
    }

    results <- list()

    # Use parallel processing if multiple other clusters
    if (cores > 1 && length(names(specificity_scores)) > 1) {
        results <- mclapply(names(specificity_scores), function(net_name) {
            update_graph_weights(net_name)
        }, mc.cores = min(cores, length(names(specificity_scores))))
    } else {
        results <- lapply(names(specificity_scores), function(net_name) {
            update_graph_weights(net_name)
        })
    }

    results <- setNames(results, names(specificity_scores))

    return(results)
}


#' Calculate Weighted Betweenness Centrality with Centralization Metrics
#'
#' @description
#' A wrapper function that calculates betweenness centrality with proper edge
#' weight handling and returns centralization metrics in the same format as
#' \code{igraph::centr_betw()}.
#'
#' @param g An igraph object. Edge weights will be used if present.
#'
#' @return
#' A list with the same structure as \code{igraph::centr_betw()}:
#' \itemize{
#'   \item \code{res}: Vector of betweenness centrality scores for each vertex
#'   \item \code{centralization}: Normalized centralization score
#'   \item \code{theoretical_max}: Theoretical maximum centralization value
#' }
#'
#' @details
#' This function extends \code{igraph::centr_betw()} by properly handling edge
#' weights when present. It calculates weighted betweenness centrality and
#' provides centralization metrics that indicate how centralized the network is
#' around its most central vertex.
#'
#' @examples
#' library(igraph)
#'
#' # Create example network
#' g <- make_star(10)
#' result <- new_centr_betw(g)
#' print(result$res) # Betweenness scores
#' print(result$centralization) # Centralization metric
#'
#' @seealso
#' \code{\link[igraph]{betweenness}}, \code{\link[igraph]{centr_betw}}
#' @importFrom igraph vcount betweenness edge_attr_names E
#' @export
#' @author Holly Yang \email{xyang2.at.uchicago.edu}
#'
new_centr_betw <- function(g) {
    library(igraph)
    if (!igraph::is_igraph(g)) {
        stop("Input 'g' must be an igraph object")
    }
    if (vcount(g) == 0) {
        stop("Graph must contain at least one vertex")
    }
    # Calculate weighted betweenness scores
    if ("weight" %in% igraph::edge_attr_names(g)) {
        betw_scores <- igraph::betweenness(g, weights = E(g)$weight)
    } else {
        # Fallback to unweighted if no weights
        betw_scores <- igraph::betweenness(g)
    }

    # Calculate centralization (to match centr_betw structure)
    n <- igraph::vcount(g)
    max_betw <- max(betw_scores)
    centralization_score <- sum(max_betw - betw_scores)

    # Theoretical maximum (for undirected graphs)
    theoretical_max <- (n - 1) * (n - 2) / 2

    # Normalize
    centralization_normalized <- if (theoretical_max > 0) {
        centralization_score / theoretical_max
    } else {
        0
    }

    # Return in same format as centr_betw()
    return(list(
        res = betw_scores,
        centralization = centralization_normalized,
        theoretical_max = theoretical_max
    ))
}


#' Network Robustness Analysis considering edge weights
#'
#' @description
#' This is an update of brainGraph::robustness function by
#' 	1) fixing issues for random attack
#' 	2) fully considering the edge weight when applicable
#' Analyzes network robustness by simulating targeted attacks or random failures,
#' measuring how the largest connected component changes as elements are removed.
#'
#' @param g An igraph object
#' @param type Attack target: "vertex" or "edge"
#' @param measure Character string specifying the attack strategy:
#'   \itemize{
#'     \item "btwn.cent": Strategic attack targeting highest betweenness centrality
#'     \item "degree": Strategic attack targeting highest degree (or strength for edge-weighted igraph) nodes (vertex only)
#'     \item "random": Random failure simulation
#'   }
#' @param N Integer. Number of Monte Carlo simulations for random attacks.
#'   Only used when \code{measure = "random"}. Default is 1000.
#'
#' @return
#' A data.table with the following columns:
#' \itemize{
#'   \item \code{type}: Type of attack performed
#'   \item \code{measure}: Attack strategy used
#'   \item \code{comp.size}: Size of largest connected component
#'   \item \code{comp.pct}: Percentage of original largest component remaining
#'   \item \code{removed.pct}: Percentage of vertices/edges removed
#' }
#' @details
#' This function simulates network robustness under different attack scenarios:
#' When edge weights are present, the function uses weighted centrality measures
#' and preserves weights throughout the analysis.
#'
#' @examples
#' library(igraph)
#' library(data.table)
#'
#' # Create example gene network
#' gene_names <- c(
#'     "ISL1", "FGF10", "MEIS1", "HAPLN1", "NTRK2", "DUSP6",
#'     "HAS2", "H1F0", "HAND2", "BMP5", "ID1", "CITED2",
#'     "BMPER", "WLS", "NKX3-1", "LAMA1", "LRRTM1", "PTPN13",
#'     "IFI16", "SLC7A2", "GENE21", "GENE22", "GENE23",
#'     "GENE24", "GENE25", "GENE26", "GENE27"
#' )
#'
#' # Generate random edges (41 edges among 27 genes)
#' set.seed(123)
#' edge_list <- t(replicate(41, sample(gene_names, 2)))
#'
#' # Create network and add weights
#' toy_network <- graph_from_edgelist(edge_list, directed = FALSE)
#' E(toy_network)$weight <- runif(ecount(toy_network), 0.2, 1.0)
#'
#' # Analyze vertex robustness under random attacks
#' result1 <- robustness_MonteCarlo(toy_network, "vertex", "random", N = 10)
#'
#' # Analyze edge robustness under betweenness centrality attacks
#' result2 <- robustness_MonteCarlo(toy_network, "edge", "btwn.cent")
#'
#' # Analyze vertex robustness under degree-based attacks
#' result3 <- robustness_MonteCarlo(toy_network, "vertex", "degree")
#'
#' @seealso
#' \code{\link[igraph]{betweenness}}, \code{\link[igraph]{degree}},
#' \code{\link[igraph]{strength}}, \code{\link[igraph]{edge_betweenness}},
#' \code{\link{new_centr_betw}}
#'
#' @importFrom igraph is_igraph vcount ecount components V E delete_vertices
#'   delete_edges betweenness degree strength edge_betweenness edge_attr_names
#'   graph_attr_names as_edgelist graph_from_edgelist
#' @importFrom data.table data.table :=
#' @export
#'
#' @author Holly Yang (xyang2_at_uchicago.edu)
#'
robustness_MonteCarlo <- function(g, type = c("vertex", "edge"), measure = c(
                                      "btwn.cent",
                                      "degree", "random"
                                  ), N = 1000) {
    library(igraph)
    library(data.table)
    # Input validation
    if (!igraph::is_igraph(g)) stop("Input 'g' must be an igraph object")
    if (!is.numeric(N) || N <= 0 || N != as.integer(N)) stop("Parameter 'N' must be a positive integer")
    if (vcount(g) == 0) stop("Graph must contain at least one vertex")
    if (ecount(g) == 0) stop("Graph must contain at least one edge")

    type <- match.arg(type)
    measure <- match.arg(measure)

    # Force sequential processing
    i <- NULL
    stopifnot(is_igraph(g))
    type <- match.arg(type)
    measure <- match.arg(measure)
    orig_max <- max(igraph::components(g)$csize)
    n <- switch(type,
        vertex = vcount(g),
        edge = ecount(g)
    )
    removed.pct <- seq.int(0, 1, length.out = n + 1L)

    if (measure == "random") {
        otype <- paste("Random", type, "removal")
        rand <- matrix(rep.int(seq_len(n), N), nrow = n, ncol = N)
        index <- apply(rand, 2L, sample)
    } else {
        otype <- paste("Targeted", type, "attack")
        max.comp.removed <- rep.int(orig_max, n)
    }

    if (type == "vertex") {
        if (measure == "random") {
            # Sequential vertex random
            max.comp <- matrix(0, nrow = n, ncol = N)
            for (i in seq_len(N)) {
                ord <- igraph::V(g)$name[index[, i]]
                tmp <- rep.int(orig_max, n)
                g.new <- g
                for (j in seq_len(n - 1L)) {
                    g.new <- igraph::delete_vertices(g.new, ord[j])
                    tmp[j + 1L] <- max(igraph::components(g.new)$csize)
                }
                max.comp[, i] <- tmp
            }
            max.comp.removed <- rowMeans(max.comp)
        } else {
            if (measure == "btwn.cent") {
                val <- new_centr_betw(g)$res
            } else {
                val <- if ("weight" %in% edge_attr_names(g)) {
                    strength(g, weights = E(g)$weight) # weighted degree
                } else {
                    degree(g) # unweighted degree
                }
            }
            ord <- V(g)$name[order(val, decreasing = TRUE)]
            for (j in seq_len(n - 1L)) {
                g <- delete_vertices(g, ord[j])
                max.comp.removed[j + 1L] <- max(igraph::components(g)$csize)
            }
        }
    } else { # edge attacks
        if (measure == "degree") {
            stop("For edge attacks, must choose \"btwn.cent\" or \"random\"!")
        } else if (measure == "random") {
            max.comp <- matrix(0, nrow = n, ncol = N)
            for (i in seq_len(N)) {
                el <- igraph::as_edgelist(g, names = FALSE)[index[, i], ]
                tmp <- rep.int(orig_max, n)
                for (j in seq_len(n - 1L)) {
                    # g.rand <- igraph::graph_from_edgelist(el[-seq_len(j), , drop = FALSE], directed = FALSE)  # lossing edge weights
                    # Instead of recreating from edgelist, use delete_edges:
                    remaining_edges <- seq_len(ecount(g))[-seq_len(j)]
                    g.rand <- delete_edges(g, seq_len(j))
                    tmp[j + 1L] <- max(igraph::components(g.rand)$csize)
                }
                max.comp[, i] <- tmp
            }
            max.comp.removed <- rowMeans(max.comp)
        } else {
            # Edge betweenness attack
            if ("weight" %in% edge_attr_names(g)) {
                edge_betw <- edge_betweenness(g, weights = E(g)$weight)
            } else {
                edge_betw <- edge_betweenness(g)
            }
            ord <- order(edge_betw, decreasing = TRUE)
            el <- as_edgelist(g, names = FALSE)[ord, ]
            for (j in seq_len(n - 1L)) {
                g <- graph_from_edgelist(el[-seq_len(j), , drop = FALSE], directed = FALSE)
                max.comp.removed[j + 1L] <- max(igraph::components(g)$csize)
            }
        }
    }

    max.comp.removed <- c(max.comp.removed, 0)
    comp.pct <- max.comp.removed / orig_max
    out <- data.table(
        type = otype, measure = measure, comp.size = max.comp.removed,
        comp.pct = comp.pct, removed.pct = removed.pct
    )

    # Add graph name if available
    if ("name" %in% graph_attr_names(g)) {
        out[, `:=`(eval(getOption("bg.group")), g$name)]
    }

    return(out)
}



#' Calculate Strength Distribution of Network Vertices
#'
#' @description
#' Computes the strength distribution (or degree distribution if unweighted)
#' of vertices in an igraph network object. This function extends
#' \code{igraph::degree_distribution()} by providing options for weighted
#' networks and different normalization strategies.
#'
#' @param g An igraph object representing the network.
#' @param use_weights Logical. Whether to use edge weights for strength
#'   calculation. If \code{FALSE} or if no weights are present, degree
#'   distribution is calculated instead. Default is \code{TRUE}.
#' @param normalized Logical. Whether to normalize the strength/degree values.
#'   For degree: normalized by (n-1) where n is number of vertices.
#'   For strength: normalized by (n-1). Default is \code{TRUE}.
#' @param cumulative Logical. Whether to return the cumulative distribution.
#'   Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{igraph::degree()}.
#'
#' @return
#' A numeric vector representing the probability density of the strength
#' (or degree) distribution. If \code{cumulative = TRUE}, returns the
#' cumulative distribution function values.
#'
#' @details
#' This function calculates the distribution of vertex strengths in weighted
#' networks or vertex degrees in unweighted networks. The strength of a vertex
#' is the sum of weights of all edges incident to that vertex.
#'
#' When \code{use_weights = FALSE} or when the graph has no edge weights,
#' the function falls back to calculating the degree distribution.
#'
#' For normalized distributions, degree values are divided by (n-1) where n
#' is the number of vertices, representing the maximum possible degree.
#'
#' @examples
#' library(igraph)
#'
#' # Create example gene network
#' gene_names <- c(
#'     "ISL1", "FGF10", "MEIS1", "HAPLN1", "NTRK2", "DUSP6",
#'     "HAS2", "H1F0", "HAND2", "BMP5", "ID1", "CITED2",
#'     "BMPER", "WLS", "NKX3-1", "LAMA1", "LRRTM1", "PTPN13",
#'     "IFI16", "SLC7A2", "GENE21", "GENE22", "GENE23",
#'     "GENE24", "GENE25", "GENE26", "GENE27"
#' )
#'
#' # Generate random edges (41 edges among 27 genes)
#' set.seed(123)
#' edge_list <- t(replicate(41, sample(gene_names, 2)))
#'
#' # Create network and add weights
#' toy_network <- graph_from_edgelist(edge_list, directed = FALSE)
#' E(toy_network)$weight <- runif(ecount(toy_network), 0.2, 1.0)
#'
#' # Calculate strength distribution
#' strength_dist <- strength_distribution(toy_network)
#'
#' # Calculate normalized degree distribution (ignoring weights)
#' degree_dist <- strength_distribution(toy_network, use_weights = FALSE)
#'
#' # Calculate cumulative strength distribution
#' cum_strength_dist <- strength_distribution(toy_network, cumulative = TRUE)
#'
#' @seealso
#' \code{\link[igraph]{degree_distribution}}, \code{\link[igraph]{strength}},
#' \code{\link[igraph]{degree}}
#'
#' @importFrom igraph is_igraph vcount ecount edge_attr_names degree strength E
#' @importFrom graphics hist
#'
#' @export
#'
#' @author: Holly Yang (xyang2_at_uchicago.edu)
#'
strength_distribution <- function(g, use_weights = TRUE, normalized = TRUE, cumulative = FALSE, ...) {
    if (!igraph::is_igraph(g)) {
        stop("Input must be an igraph object")
    }
    # Handle empty graphs
    if (vcount(g) == 0) {
        return(numeric(0))
    }

    if ((!"weight" %in% edge_attr_names(g)) | (!use_weights)) {
        cs <- degree(g, normalized = normalized, ...)
        # Handle empty degree sequence
        if (length(cs) == 0) {
            return(numeric(0))
        }
        if (!normalized) {
            hi <- hist(cs, -1:max(cs), plot = FALSE)$density
        } else {
            n_bins <- min(100, length(unique(cs)) * 2) # Adaptive number of bins
            breaks <- seq(0, max(cs) + 0.001, length.out = n_bins + 1)
            hi <- hist(cs, breaks = breaks, plot = FALSE)$density
        }
    } else {
        cs <- strength(g, weights = E(g)$weight)
        # Handle empty graphs
        if (length(cs) == 0) {
            return(numeric(0))
        }
        if (normalized) cs <- cs / (vcount(g) - 1)

        # Handle the case where cs might have non-integer values (from strength)
        # Create proper breaks that span the full range
        # Create appropriate bins
        if (normalized) {
            # For normalized: create many small bins between 0 and max
            n_bins <- min(100, length(unique(cs)) * 2) # Adaptive number of bins
            breaks <- seq(0, max(cs) + 0.001, length.out = n_bins + 1)
        } else {
            min_cs <- floor(min(cs))
            max_cs <- ceiling(max(cs))
            breaks <- seq(min_cs - 0.5, max_cs + 0.5, by = 1)
        }
        hi <- hist(cs, breaks = breaks, plot = FALSE)$density
    }

    if (!cumulative) {
        res <- hi
    } else {
        res <- rev(cumsum(rev(hi)))
    }
    res
}
