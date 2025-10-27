############ load libraary for testing #########
library(SingleCellExperiment)
library(scater)
library(Matrix)
library(qlcMatrix)
library(igraph)
library(parallel)

BioTIP_version <- '06232025'
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_', BioTIP_version, '.R'))

################  The above lines should be removed when submission to Bioconductor !!!###############################
## version 10 add synthetic_simulation() and more visualization functions

## version 10 update robustness_MonteCarlo(..., measure = "btwn.cent") to use 1/weight when calling new_centr_betw() and betweenness()  
## version 10 add parameter graph_list into update_graph_weights()
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
#' 1. **Ratio** – The ratio of absolute co-expression in the target cell type to
#'    the mean absolute co-expression across other cell types:
#'    \deqn{ ratio = |r_{target}| / (mean(|r_{other}|) + 0.01) }
#'    Higher values indicate that a gene pair is more strongly co-expressed in
#'    the target cell type than elsewhere. This score penalizes high co-expression 
#     elsewhereis via the denominator, therefore, it is very sensitive to 
#'    very low background co-expression (division amplifies small denominators).
#'
#' 2. **Z-score** – A standardized measure showing how many standard deviations
#'    the target co-expression differs from the mean across other cell types:
#'    \deqn{ zscore = (|r_{target}| - mean(|r_{other}|)) / (sd(|r_{other}|) + 0.01) }
#'    The z-score measures how far the target’s co-expression deviates from 
#'    the cross-type average in units of standard deviation. Therefore, it  
#'    highlights edges whose co-expression is unusually strong (or weak) in the target.
#'
#' 3. **Difference** – The simple difference between target and mean co-expression:
#'    \deqn{ diff = |r_{target}| - mean(|r_{other}|) }
#'    Positive values indicate stronger co-expression in the target cell type.
#'	  It gives credit as long as the target is stronger than the mean, good for 
#'    identifying edges that are stronger in the target than others (quantitative gain).
#'
#' 4. **Combined** – A multiplicative metric rewarding high co-expression in the target
#'    and low co-expression in other cell types:
#'    \deqn{ combined = |r_{target}| * (1 - mean(|r_{other}|)) }
#'    This score Rewards high co-expression in target and explicitly penalizes high 
#'    co-expression elsewhere through (1 - mean).
#'    Emphasizes edges that are both strong and specific to the target cell type.
#'
#' All correlations are converted to absolute values before comparison to capture
#' connection strength regardless of direction (positive or negative correlation).
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
#' @param specificity_method Character string specifying which specificity measure
#'   to use when combining STRING confidence with co-expression specificity scores.
#'   Options are "combined" (default), "ratio", "zscore", or "diff".
#'
#'   See \code{\link{calculate_network_specificity}} for detailed definitions
#'   of these four measures.
#'
#' @param max_min_norm Logical indicating whether to apply max-min normalization
#'   to both STRING and specificity matrices before combining. Defaults to FALSE.
#' @param verbose Logical indicating whether to print detailed progress messages.
#'
#' @param cores Integer specifying the number of cores to use for parallel processing.
#'
#' @return A list of igraph objects with the following updated edge features, where n_e is the number of edges.
#' \itemize{
#'   \item \code{...}: Besides 'weight', other edge features initiated from the input igraph object.
#'   \item \code{weight}: Final combined edge weights (nornalized PPI * specificity). A positive numeric vector. Larger edge weights correspond to stronger connections.
#'   \item \code{norm_PPI_score}:  Original PPI scores (STRING combined_score/1000)
#'   \item \code{corexp_sign}: Sign of co-expression correlation ('positive' or 'negative')
#'   \item \code{coexp_target}: Cell-type specific co-expression scores used for weighting
#' }
#'
#' @details
#' The function combines the original edge weights with cell-type specificity scores
#' to produce weights that favor edges with both high confidence protein interactions
#' and cell-type specific co-expression patterns. By default, larger edge weights 
#' correspond to stronger connections.
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
#'     # [1] "weight"          "norm_PPI_score" "corexp_sign"     "coexp_target"
#'
#'     plot(E(updated_networks[["A"]])$coexp_target, E(updated_networks[["A"]])$weight,
#'         xlab = "Pearson cor", ylab = "Coexp&PPI-combined Weights", main = "Cell cluster A-specific"
#'     )
#' }
#'
#' @importFrom igraph V E ecount ends
#' @export
#' @author Holly Yang \email{xyang2.at.uchicago.edu}
#' @author Felix Yu \email{felixy.at.uchicago.edu}
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

    # an internal function to update each igraph object of the inputting graph_list
    update_graph_weights <- function(net_name, graph_list) {
        if (verbose) cat("\nUpdating weights for network:", net_name, "\n")

        if (!net_name %in% names(graph_list)) {
            warning("Network '", net_name, "' not found in graph_list. Skipping.")
            return(NULL)
        }

        # Get the graph
        g <- graph_list[[net_name]]
        E(g)$norm_PPI_score <- E(g)$weight # Store original normalized PPI weights

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
        new_weights <- E(g)$norm_PPI_score # default: keep old weights

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
            update_graph_weights(net_name, graph_list)
        }, mc.cores = min(cores, length(names(specificity_scores))))
    } else {
        results <- lapply(names(specificity_scores), function(net_name) {
            update_graph_weights(net_name, graph_list)
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
#' weights when present, being interpreted as distances. 
#' This function calculates weighted betweenness centrality and
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
        betw_scores <- igraph::betweenness(g, weights = E(g)$weight)  ## weight should be distance-based weight !!
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
#' When edge weight present as a positive numeric vector, larger edge weights correspond to stronger connections. 
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
#' Edge weight is an optional positive weight vector. When edge weights are present, 
#' the function uses weighted centrality measures
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
	
	if (measure == "btwn.cent" & "weight" %in% edge_attr_names(g)) { # v10
	# since v9, graph_list contains connection-based E(g)$weight, but betweenness requires distance-based weight
        E(g)$weight = 1/E(g)$weight
	}			

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
                edge_betw <- edge_betweenness(g, weights = E(g)$weight) # new_centr_betw takes distance weight as input
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
       # out[, `:=`(eval(getOption("bg.group")), g$name)]
       out[, graph_name := graph_attr(g, "name")]  # v11
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


#### additional visualization fuctions in v10 #######

#' Extract edge weights and assign PPI categories from a list of graphs
#'
#' This function iterates over a list of igraph objects, extracts edge weights,
#' and classifies each network into CTS, HiGCTS, or HiG categories based on
#' graph names. For each network, it returns a tidy data frame of edge weights
#' with associated metadata (network name, PPI category, cluster ID, and edge count).
#' Optionally, clusters listed in `unstable_cluster_ID` are labeled as 'unstable'.
#'
#' @param graph_list A named list of igraph objects containing PPI networks.
#' @param PPI_color_palette A named color vector for PPI categories (not directly used here but for consistency).
#' @param unstable_cluster_ID A character vector of cluster IDs to mark as 'unstable'.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{sample}{Graph name.}
#'     \item{PPI_cat}{Category of PPIN: 'CTS', 'HiGCTS', or 'HiG'.}
#'     \item{edge_weight}{Numeric weight of each edge.}
#'     \item{num_edges}{Number of edges in the network.}
#'     \item{cluster_ID}{Cluster identifier parsed from the graph name.}
#'     \item{cluster_cat}{Cluster classification: 'unstable' or 'stable'.}
#'   }
#'
#' @details
#' If a network lacks explicit edge weights, all edges are assigned a weight of 1.
#' Graph names are assumed to follow the pattern "CTS_", "HiGCTS_", or "HiG_".
#'
#' @examples
#' edge_df <- extract_edge_weights_by_category(graph_list, PPI_color_palette, unstable_cluster_ID)
#'
extract_edge_weights_by_category <- function(graph_list, PPI_color_palette, unstable_cluster_ID) {  
  # Initialize storage for edge weights
  all_edge_data <- data.frame()  
  # Extract edge weights from each graph
  for (graph_name in names(graph_list)) {
    graph <- graph_list[[graph_name]]    
    # Get edge weights (assumes edges have 'weight' attribute)
    edge_weights <- E(graph)$weight    
    # If no weights, use default of 1
    if (is.null(edge_weights)) {
      edge_weights <- rep(1, ecount(graph))
    }
    # Assign PPI category
    PPI_cat <- case_when(
      grepl("^HiGCTS_", graph_name) ~ "HiGCTS",
      grepl("^HiG_", graph_name) ~ "HiG", 
      grepl("^CTS_", graph_name) ~ "CTS",
      TRUE ~ "Other"
    )    
    # Store data
    edge_data <- data.frame(
      sample = graph_name,
      PPI_cat = PPI_cat,
      edge_weight = edge_weights,
      num_edges = length(edge_weights),
	  cluster_ID = unlist(strsplit(graph_name, split='_'))[2]
    )    
    all_edge_data <- rbind(all_edge_data, edge_data)
  }  
  # Filter to main categories
  all_edge_data <- all_edge_data[all_edge_data$PPI_cat %in% c("CTS", "HiGCTS", "HiG"), ]
  all_edge_data$PPI_cat <- factor(all_edge_data$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))
  all_edge_data$cluster_cat <- ifelse(all_edge_data$cluster_ID %in% unstable_cluster_ID, 'unstable', 'stable')
  
  return(all_edge_data)
}

#' Plot and summarize PPI edge-weight distributions by category
#'
#' This function visualizes and statistically compares the distribution of
#' edge weights across PPI categories. It generates (1) density plots for each
#' category and (2) combined violin/box plots with pairwise Wilcoxon tests.
#'
#' @param edge_data A data frame produced by `extract_edge_weights_by_category()`
#'   containing columns `PPI_cat` and `edge_weight`.
#' @param PPI_color_palette A named vector of colors corresponding to PPI categories.
#'
#' @return A list containing:
#'   \describe{
#'     \item{density_plot}{Faceted density plots showing edge-weight distributions per PPI category.}
#'     \item{boxplot}{Combined box-violin plots with Wilcoxon significance comparisons.}
#'     \item{summary_stats}{Data frame summarizing mean, median, and total number of edges per category.}
#'   }
#'
#' @details
#' Pairwise comparisons (CTS–HiGCTS, CTS–HiG, HiGCTS–HiG) are performed using
#' `wilcox.test`, with Bonferroni correction and significance labels displayed
#' directly on the plot. The density and box-violin plots share a common
#' color scheme for category consistency.
#'
#' @examples
#' plots <- plot_edge_weight_distributions(edge_df, PPI_color_palette)
#' plots$density_plot
#' plots$boxplot
#' 
plot_edge_weight_distributions <- function(edge_data, PPI_color_palette) {  
  # Calculate summary statistics
  summary_stats <- edge_data %>%
    group_by(PPI_cat) %>%
    summarise(
      mean_weight = round(mean(edge_weight), 3),
      median_weight = round(median(edge_weight), 3),
      total_edges = n(),
      .groups = 'drop'
    )  
  # 1. Density plot for HiG PPIs, colored by cluster_categories (stable or instable)
 p1 <- ggplot(edge_data, aes(
    x = log10(edge_weight + 1),         # log10 transform inside aes()
    fill = PPI_cat, 
    color = PPI_cat)) + 
  geom_density(alpha = 0.3, size = 1.2, adjust = 1.2) +  # adjust controls smoothness 
  #facet_wrap(~PPI_cat, ncol = 1, scales = "free_y") +
  labs(
    title = "hiPSC",
    x = expression(log[10]~"(PPI edge weight + 1)"),
    y = "Density"
  ) +  
  scale_fill_manual(values = PPI_color_palette) +
  scale_color_manual(values = PPI_color_palette) + 
  # Optional: label ticks in original (non-log) units
  scale_x_continuous(
    breaks = seq(0, 2, by = 0.5),
    labels = c("1", "3", "10", "32", "100")
  ) +  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
   # strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold"),
   # plot.subtitle = element_text(hjust = 0.5, color = "gray60"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

  
  # 2. violin comparison
  p2 <- ggplot(edge_data, aes(x = PPI_cat, y = edge_weight, fill = PPI_cat)) +
    # geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_violin(alpha = 0.3, width = 0.8, trim = TRUE) +
    stat_compare_means(comparisons = list(c("HiG", "HiGCTS"), c("HiG", "CTS"), c("HiGCTS", "CTS")), 
                      method = "wilcox.test",
					  label = "p.signif",
                      p.adjust.method = "bonferroni"
                      ) +
    labs(
      title = "wilcox.test",
      subtitle = paste("Total edges - CTS:", summary_stats$total_edges[1],
                       " | HiGCTS:", summary_stats$total_edges[2],
                       " | HiG:", summary_stats$total_edges[3]),
      x = "PPI Category",
      y = "Edge weight"
    ) +
    scale_fill_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray60")
    )
    
  return(list(
    density_plot = p1,
    boxplot = p2,
    summary_stats = summary_stats
  ))
}

#' Plot network-wise edge betweenness metrics across PPIN categories
#'
#' This function visualizes the number of network communities identified
#' by edge-betweenness clustering across three categories of
#' protein–protein interaction networks (PPINs): CTS, HiGCTS, and HiG.
#' It produces (1) a line plot connecting HiG clusters and overlaying
#' corresponding CTS and HiGCTS transition states, and (2) a boxplot
#' comparing the three categories with statistical tests.
#'
#' @param nEB_data A named numeric vector containing the number of
#'   communities per PPIN (e.g., from edge-betweenness clustering). 
#'   Names must encode sample categories with prefixes such as
#'   "CTS_", "HiGCTS_", or "HiG_".
#' @param PPI_color_palette A named vector of colors for each PPIN
#'   category, e.g. `c(CTS = "blue", HiGCTS = "pink", HiG = "gold")`.
#' @param method Character string specifying the statistical test to use
#'   for pairwise comparisons in the boxplot. Either `"t.test"` or `"wilcox.test"`.
#'   Defaults to `"t.test"`.
#'
#' @return A list with two ggplot objects:
#'   \describe{
#'     \item{line_plot}{A line plot showing the HiG cluster trend
#'       (ordered by increasing number of communities) with CTS and
#'       HiGCTS points overlaid at corresponding positions.}
#'     \item{boxplot}{A boxplot with optional jitter displaying the
#'       distribution of community counts across PPIN categories,
#'       annotated with pairwise significance results.}
#'   }
#'
#' @details
#' - The function first classifies each sample name into one of the
#'   PPIN categories (`CTS`, `HiGCTS`, or `HiG`) and extracts the
#'   associated cluster identifier.
#' - HiG networks are sorted by their community counts to define
#'   the x-axis order, and corresponding CTS/HiGCTS samples with
#'   matching cluster IDs are aligned above these positions.
#' - The boxplot compares distributions among the three PPIN categories,
#'   using either a t-test or Wilcoxon test for significance annotation.
#'
#' @examples
#' plots <- plot_nEB_ggplot(nEB_data, PPI_color_palette, method = "wilcox.test")
#' plots$line_plot
#' plots$boxplot
#'
plot_nEB_ggplot <- function(nEB_data, PPI_color_palette, method=c('t.test', 'wilcox.test')) {
  method = match.arg(method)
  # Convert to data frame
  plot_data <- data.frame(
    sample = factor(names(nEB_data), levels = names(nEB_data)),
    value = as.numeric(nEB_data),
    index = 1:length(nEB_data),
	stringsAsFactors = FALSE  #not to automatically convert character strings into factors when creating a data frame.
  )
  # Assign categories based on sample names
  plot_data$PPI_cat <- case_when(
    grepl("^HiGCTS_", plot_data$sample) ~ "HiGCTS",
    grepl("^HiG_", plot_data$sample) ~ "HiG", 
    grepl("^CTS_", plot_data$sample) ~ "CTS",
    TRUE ~ "Other"
  )
  # Extract sample cluster (remove prefix)
  plot_data$sample_type <- case_when(
    plot_data$PPI_cat == "HiG" ~ gsub("^HiG_", "", plot_data$sample),
    plot_data$PPI_cat == "HiGCTS" ~ gsub("^HiGCTS_", "", plot_data$sample),
    plot_data$PPI_cat == "CTS" ~ gsub("^CTS_", "", plot_data$sample),
    TRUE ~ plot_data$sample  #catches anything that didn't match the previous conditions
  )
  # Get HiG samples as the base x-axis (13 positions)
  hig_data <- plot_data[plot_data$PPI_cat == "HiG", ]
  hig_data <- hig_data[order(hig_data$value), ]  # Order by value (smallest to largest)  # Maintain original order
  hig_data$x_position <- 1:nrow(hig_data)
  # Map other categories to corresponding HiG positions
  other_data <- plot_data[plot_data$PPI_cat != "HiG", ]
  # Create mapping for overlaid points
  overlay_data <- data.frame()  
  for (i in 1:nrow(other_data)) {
    sample_type <- other_data$sample_type[i]   
    # Find corresponding HiG position
    matching_hig <- which(hig_data$sample_type == sample_type)    
    if (length(matching_hig) > 0) {
      # If there's a matching HiG sample type, use that position
      x_pos <- hig_data$x_position[matching_hig[1]]
    } else {
      # If no exact match, skip this point or handle differently
      next
    }    
    overlay_data <- rbind(overlay_data, data.frame(
      sample = other_data$sample[i],
      value = other_data$value[i],
      PPI_cat = other_data$PPI_cat[i],
      sample_type = sample_type,
      x_position = x_pos
    ))
  }  
  
  # Create the line plot
  p1 <- ggplot() +
    # First, draw line connecting only HiG samples
    geom_line(data = hig_data, aes(x = x_position, y = value), 
              color = PPI_color_palette["HiG"], size = 1.2, alpha = 0.8) +
    # Add HiG points
    geom_point(data = hig_data, aes(x = x_position, y = value, color = PPI_cat), 
               size = 3, alpha = 0.9) +
    # Add overlaid points for other categories
    geom_point(data = overlay_data, aes(x = x_position, y = value, color = PPI_cat), 
               size = 3, alpha = 0.9) +
    labs(
      title = "cluster_edge_betweenness",
      subtitle = paste("HiG line (n=", sum(plot_data$PPI_cat == "HiG"), 
                       ") | HiGCTS dots (n=", sum(plot_data$PPI_cat == "HiGCTS"),
                       ") | CTS dots (n=", sum(plot_data$PPI_cat == "CTS"), ")", sep = ""),
      x = "cell clusters, transition states are named",
      y = "number of communities in PPI",
      color = "PPI_cat"
    ) +
    scale_x_continuous(
      breaks = 1:nrow(hig_data),
      labels = hig_data$sample
    ) +
    scale_color_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray60"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  # Create boxplot comparing the three PPI_cat categories**
  # Filter data to include only the three main categories and only the transition clusters !!
  boxplot_data <- plot_data[plot_data$PPI_cat %in% c("CTS", "HiGCTS","HiG" ) & 
				plot_data$sample_type %in% overlay_data$sample_type, ]
  # Set factor levels to control order
  boxplot_data$PPI_cat <- factor(boxplot_data$PPI_cat, levels = c("CTS","HiGCTS", "HiG" ))
  
  comparisons <- list(c("HiG", "HiGCTS"), c("HiG", "CTS"), c("HiGCTS", "CTS"))
  # Calculate summary statistics for subtitle
  summary_stats <- boxplot_data %>%
      group_by(PPI_cat) %>%
      summarise(
        mean_val = round(mean(value),1),
        median_val = round(median(value), 1),
        .groups = 'drop'
    )

  p2 <- ggplot(boxplot_data, aes(x = PPI_cat, y = value, fill = PPI_cat)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
    geom_jitter(aes(color = PPI_cat), width = 0.2, size = 2, alpha = 0.8) +
    labs(
      title = paste(method, "test, transition clusters"),
      subtitle = paste("Mean/Median - CTS:", summary_stats$mean_val[1], "/", summary_stats$median_val[1],
                       " | HiGCTS:", summary_stats$mean_val[2], "/", summary_stats$median_val[2],
                       " | HiG:", summary_stats$mean_val[3], "/", summary_stats$median_val[3]),
      x = "PPI Category",
      y = "number of communities in PPI"
    ) +
    scale_fill_manual(values = PPI_color_palette) +
    scale_color_manual(values = PPI_color_palette) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray60"),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 11),
      legend.position = "none",  # Remove legend since colors match categories
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
	stat_compare_means(comparisons = comparisons, 
                      method = method,
                      p.adjust.method = 'none',#"bonferroni",
                      label = "p.signif",
                      bracket.size = 0.6,
                      step.increase = 0.1) 
  
   return(list(line_plot = p1, boxplot = p2))
}


#' @title Synthetic network simulation and robustness validation
#'
#' @description
#' This function generates synthetic weighted networks (random, scale-free, small-world, and
#' degree-preserving rewired) matched to an empirical PPIN, assigns edge weights drawn
#' from the empirical distribution, and evaluates fragmentation resilience under targeted-node
#' attack using the \code{robustness_MonteCarlo()} function.
#' It returns area-under-curve (AUC) statistics and publication-ready ggplot objects
#' summarizing edge-weight distributions, fragmentation curves, and AUC bar plots.
#'
#' @param g_real An \code{igraph} object representing the real PPIN.
#'               Must contain a numeric edge attribute \code{"weight"}.
#' @param main  Optional character string used as the plot title and ID in the
#'              output data frame. Defaults to "Network".
#'
#' @details
#' **Workflow**
#' \enumerate{
#'   \item Extract the empirical edge-weight distribution from \code{g_real}.
#'   \item Generate four synthetic topologies with the same node and edge counts:
#'         \itemize{
#'           \item Erdős–Rényi random graph via \code{sample_gnp()}.
#'           \item Scale-free graph via \code{sample_pa()}.
#'           \item Small-world graph via \code{sample_smallworld()}.
#'           \item Degree-preserving rewired graph via \code{rewire(..., with = keeping_degseq())}.
#'         }
#'   \item Assign edge weights to all synthetic networks by resampling from the empirical
#'         distribution of \code{E(g_real)$weight}, maintaining overall weight heterogeneity.
#'   \item Evaluate network fragmentation for each network type using
#'         \code{robustness_MonteCarlo(type = "vertex", measure = "btwn.cent")}.
#'   \item Compute the area under the fragmentation curve (AUC) as a measure of
#'         network resilience using \code{pracma::trapz()}.
#'   \item Return an AUC summary table and three ggplot objects:
#'         \enumerate{
#'           \item \code{p_weights} — histogram/density of edge-weight distributions.
#'           \item \code{p_line} — fragmentation curves (fraction of vertices removed vs.
#'                 fraction of largest component remaining).
#'           \item \code{p_AUC} — bar plot of AUC (higher = more robust).
#'         }
#' }
#'
#' @return
#' A named list with the following elements:
#' \describe{
#'   \item{\code{auc_summary}}{Data frame of AUC and normalized AUC for each network type.}
#'   \item{\code{p_weights}}{ggplot histogram/density comparing edge-weight distributions.}
#'   \item{\code{p_line}}{ggplot line plot of fragmentation curves.}
#'   \item{\code{p_AUC}}{ggplot bar plot of AUC values.}
#'	 \item{\code{network_colors}}{named color platte used here.}
#' }
#'
#' @note
#' The degree-preserving rewired network retains each node's degree but randomizes edge
#' connections, thereby disrupting modular organization.  If this network collapses even
#' faster than the real PPIN, it indicates that modular structure confers partial buffering
#' within the real network.
#'
#' @examples
#' data(graph_CTS_CP.1) 
#' g_real <- graph_CTS_CP.1  # graph_list[["CTS_CP.1"]]
#' result <- synthetic_simulation(g_real, main = "CTS_CP.1")
#' result$auc_summary
#' result$p_weights + result$p_line + result$p_AUC
#'
#' @import igraph data.table ggplot2 pracma
#' @export

synthetic_simulation = function(g_real, main=NULL){
	if(is.null(main)) main = 'Network'
	
	# g_real = graph_list[["CTS_CP.1"]]
	################################################
	# Step 1. Extract empirical weight distribution
	emp_weights <- E(g_real)$weight
	summary(emp_weights)
	#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	# 0.00008 0.01101 0.01947 0.02692 0.03209 0.41399 

	################################################
	# step 2: Generate synthetic networks 
	n_nodes <- vcount(g_real)
	# approximate density of the real PPIN
	p_edge  <- ecount(g_real) / choose(n_nodes, 2)

	# 1. random (Erdős–Rényi)
	g_random <- sample_gnp(n = n_nodes, p = p_edge)
	# 2. scale-free (preferential attachment)
	m_links  <- round(ecount(g_real) / n_nodes)  # average edges per new node
	g_scale_free <- sample_pa(n = n_nodes, m = m_links, directed = FALSE)
	# 3. small-world
	g_small_world <- sample_smallworld(dim = 1, size = n_nodes, nei = 3, p = 0.05)
    # 4. Add degree-preserving randomization control (optional)
	# This “degree-preserving” null keeps local degree statistics but removes biological wiring patterns and topological coherence
	g_deg_preserved <- rewire(g_real, with = keeping_degseq(niter = 1e5))


	################################################
	# Step 3 – assign synthetic edge weights using the empirical distribution
	# Add continuous weights so that your synthetic graphs resemble the real PPIN in weight heterogeneity.
	# This step ensures your robustness tests and DNB/Ic calculations remain comparable.

	# resample empirical weights
	E(g_random)$weight     <- sample(emp_weights, ecount(g_random), replace = TRUE)
	E(g_scale_free)$weight <- sample(emp_weights, ecount(g_scale_free), replace = TRUE)
	E(g_small_world)$weight <- sample(emp_weights, ecount(g_small_world), replace = TRUE)
	E(g_deg_preserved)$weight <- sample(emp_weights, ecount(g_deg_preserved), replace = TRUE)

	################################################
	# Step 4 – verify weight distributions

	# hist(emp_weights, breaks = 100, col = "gray", freq = FALSE,
		 # main = "Edge-weight distributions", xlab = "Weight")
	# hist(E(g_random)$weight, breaks = 100, col = rgb(1,0,0,0.3), add = TRUE)
	# legend("topright", legend = c("Real PPIN", "Random network"), 
		   # fill = c("gray", rgb(1,0,0,0.3)), bty = "n")
	# # Combine real and random weights into one data frame
	df_weights <- data.frame(
	  weight = c(emp_weights, E(g_random)$weight),
	  network = rep(c("Real PPIN", "Random network"),
					times = c(length(emp_weights), ecount(g_random)))
	)

	# ggplot histogram / density overlay
	p_weights <- ggplot(df_weights, aes(x = weight, fill = network, color = network)) +
	  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.4, bins = 80) +
	  geom_density(alpha = 0.6) +
	  scale_fill_manual(values = c("gray50", "#E64B35")) +
	  scale_color_manual(values = c("gray20", "#E64B35")) +
	  labs(
		title = "Edge-weight distribution",
		x = "Edge weight",
		y = "Density"
	  ) +
	  theme_classic(base_size = 14) +
	  theme(
		legend.title = element_blank(),
		plot.title = element_text(hjust = 0.5, face = "bold")
	  )		   
		   
	################################################
	# Step 5 – run your robustness_MonteCarlo()
	graph_attr(g_real, "name") <- "real_PPI"
	graph_attr(g_random, "name") <- "random_network"
	graph_attr(g_scale_free, "name") <- "scale_free"
	graph_attr(g_small_world, "name") <- "small_world"
	graph_attr(g_deg_preserved, "name") <- "degree_preserving"

	sim_real   <- robustness_MonteCarlo(g_real, type = "vertex", measure = "btwn.cent", N = 100)
	sim_random <- robustness_MonteCarlo(g_random, type = "vertex", measure = "btwn.cent", N = 100)
	sim_scale  <- robustness_MonteCarlo(g_scale_free, type = "vertex", measure = "btwn.cent", N = 100)
	sim_small  <- robustness_MonteCarlo(g_small_world, type = "vertex", measure = "btwn.cent", N = 100)
	sim_deg <- robustness_MonteCarlo(g_deg_preserved, type="vertex", measure="btwn.cent", N=100)
	
	################################################
	# Step 6 -  compare AUCs 

	auc_real   <- trapz(sim_real$removed.pct,   sim_real$comp.pct)
	auc_random <- trapz(sim_random$removed.pct, sim_random$comp.pct)
	auc_scale   <- trapz(sim_scale$removed.pct,   sim_scale$comp.pct)
	auc_small   <- trapz(sim_small$removed.pct,   sim_small$comp.pct)
	auc_deg <- trapz(sim_deg$removed.pct, sim_deg$comp.pct)

	#Combine and summarize AUC results
	auc_summary <- data.frame(
	  network = c("real_PPIN", "random", "scale_free", "small_world", "deg_preserving"),
	  AUC = c(auc_real, auc_random, auc_scale, auc_small, auc_deg )
	)

	auc_summary$normalized_AUC <- auc_summary$AUC / max(auc_summary$AUC)
	auc_summary
		  # network       AUC normalized_AUC
	# 1   real_PPIN 0.2327409      0.2342621
	# 2      random 0.9935065      1.0000000
	# 3  scale_free 0.9935065      1.0000000
	# 4 small_world 0.9935065      1.0000000
	# 5 deg_preserving 0.1541933      0.1552011
	auc_summary$ID = main
	
	################################################
	# Step 7 — Visualize fragmentation curves

	# Combine all simulation outputs
	sim_all <- rbindlist(list(
	  cbind(sim_real,       network = "real_PPIN"),
	  cbind(sim_random,     network = "random"),
	  cbind(sim_scale,      network = "scale_free"),
	  cbind(sim_small,      network = "small_world"),
	  cbind(sim_deg,      network = "deg_preserving")
	), fill = TRUE)

	# Plot
	network_colors <- c(
	  "real_PPIN" = "#E64B35FF",
	  "deg_preserving" = "#D39200FF",
	  "random" = "#4DBBD5FF",
	  "scale_free" = "#00A087FF",
	  "small_world" = "#3C5488FF"
	)
	
	network_levels <- c("real_PPIN", "deg_preserving", "random", "scale_free", "small_world")
	sim_all$network    <- factor(sim_all$network, levels = network_levels)
	auc_summary$network <- factor(auc_summary$network, levels = network_levels)

	p_line = ggplot(sim_all, aes(x = removed.pct, y = comp.pct, color = network)) +
	  geom_line(size = 1.2) +
	  theme_classic(base_size = 12) +
	  labs(
		title = paste(main, "fragmentation under targeted-node attack"),
		x = "Fraction of vertices removed",
		y = "Fraction of largest component remaining"
	  ) +
	  scale_color_manual(values = network_colors) +
	  scale_fill_manual(values = network_colors)	

	# Plot AUC summary
	p_AUC = ggplot(auc_summary, aes(x = network, y = AUC, fill = network)) +
	  geom_col(width = 0.7) +
	  theme_minimal(base_size = 12) +
	  labs(
		title = paste(main, "resilience (AUC of fragmentation curve)"),
		x = "Network type",
		y = "AUC (higher = more robust)"
	  ) +  
	  scale_color_manual(values = network_colors) +
	  scale_fill_manual(values = network_colors)  +	
	  theme(legend.position = "none", 
		axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9))
  

	################################################
	# Step 8 — return
	outputs = list(auc_summary=auc_summary, 
					p_weights = p_weights,
					p_line = p_line,
					p_AUC = p_AUC,
					network_colors = network_colors)

	return(outputs)
}


#' @title Plot weighted protein–protein interaction (PPIN) network
#'
#' @description
#' Generate a weighted PPIN visualization using ggraph, where
#' edge width and node size reflect interaction strength.
#' Positive and negative co-expression edges are colored orange and blue,
#' and optionally, disease-related (e.g., CHD) genes are highlighted in red.
#'
#' @param g An \code{igraph} object containing both vertex and edge attributes.
#'   Expected attributes include:
#'   \itemize{
#'     \item \strong{Vertex attributes:}
#'       \code{name} (character) — gene identifier for each node;
#'       \code{weight} (numeric) — node-level score, representing statistical strength (e.g., |Wilcox score|).
#'     \item \strong{Edge attributes:}
#'       \code{weight} (numeric) — edge strength, where higher values indicate stronger interactions;
#'       \code{corexp_sign} (character) — co-expression direction with expected values
#'       \code{"positive"} and/or \code{"negative"}.
#'   }
#' @param layout Character, graph layout algorithm to use.
#'   Defaults to \code{"fr"} (Fruchterman–Reingold).
#' @param CHD Optional character vector of gene names to highlight.
#'
#' @details
#' The function removes isolated nodes (degree = 0) before plotting.
#' Node size corresponds to \code{V(g)$weight}, and edge width corresponds to
#' \code{E(g)$weight}. Edges are colored by \code{E(g)$corexp_sign}.
#' Curated CHD genes (if provided) are highlighted in red.
#'
#' @return A \code{ggplot} object showing the weighted PPIN.
#' @examples
#' plot_weighted_PPIN(graph_list[["HiGCTS_CP.1"]], CHD = CHD_genes)
#' @import ggraph ggplot2 igraph
#' @export
#'
plot_weighted_PPIN = function(g, layout = "fr", 
		CHD = NULL, node_size_title = "|Wilcox score|") {
	# ---- Attribute checks ----
	if (!"weight" %in% vertex_attr_names(g))
		stop("Please assign V(g)$weight for plotting node size")

	if (!"weight" %in% edge_attr_names(g))
		stop("Please assign E(g)$weight for plotting edge width")

	if (!"corexp_sign" %in% edge_attr_names(g))
		stop("Please assign E(g)$corexp_sign for plotting edge color")
    # ---- Validate that edge signs include 'positive' and 'negative' ----
    expected_signs <- c("positive", "negative")
    missing_signs <- setdiff(expected_signs, unique(E(g)$corexp_sign))
	if (length(missing_signs) > 1)  
		stop("Edge attribute 'corexp_sign' does not include expected values: 'positive', 'negative'")
	
	# ---- Disease-gene highlight ----
	if (is.null(CHD)) {
		warning("No disease genes supplied for highlighting (CHD = NULL)")
		V(g)$is_CHD <- FALSE
	  } else {
		V(g)$is_CHD <- V(g)$name %in% CHD
	  }

	# ---- Remove isolated nodes ----
	g_connected <- delete_vertices(g, which(degree(g) == 0))
		
	set.seed(1234)  # Use to ensure layout repeatable  
	layout_coords <- create_layout(g_connected, layout = layout)  
	
	p  = ggraph(layout_coords) +  
	  geom_edge_link(aes(
		width = weight,                   
		color = corexp_sign           
	  ), alpha = 0.7) +	  
	  geom_node_point(aes(
		size = weight,    
		color = is_CHD   
		#shape = is_CHD
	  )) +  
	  geom_node_text(aes(label = name), repel = TRUE, size = 3) +  # <--- Add labels
	  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "gray70")) +
	  scale_edge_color_manual(values = c("positive" = "orange", "negative" = "blue")) +
	  scale_size_continuous(range = c(2, 5), name = node_size_title) +
	  scale_edge_width_continuous(range = c(0.1, 1), 
					limits = range(E(g)$weight, na.rm = TRUE), 
					name = "E weights") +	 
	  theme_void() +
	  labs(
		 edge_color = "Specific co-exp",
		 color = "Curated CHD genes"
	  ) +
	  theme(legend.position = "right") +
	  ggtitle(paste0(db, ': ', int, ' ', vcount(g_connected), '/', vcount(g), ' PPI genes'))

	return(p)
}	