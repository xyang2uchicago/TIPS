# 04/27/2026 — fix bug coexp_target -> coexp_focal in calculate_network_specificity();
#              add cts_overlap_by_cluster(), fmt_lm_term(), make_lineage_pseudobulk_plot(),
#              fmt_lm(), get_endpoint_effect(), get_cf4_effect(), get_signature_genes(),
#              sample_matched_genes(), run_random_null_test(), pseudobulk_volcano_plot()
############ load library for testing #########
library(SingleCellExperiment)
library(scater)
library(Matrix)
library(qlcMatrix)
library(igraph)
library(parallel)

BioTIP_version <- "06232025"
source(paste0("https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_", BioTIP_version, ".R"))

################  The above lines should be removed when submission to Bioconductor !!!###############################
## version 2026, add functions to prioritize lineage-pull regulator from CTS
## (can extended to one key upstream regulators from cisTarget motif target-enrichment analysis of CTS)
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
#' co-expression of genes in a focal cell type/cluster to their co-expression
#' in other cell types/clusters.
#'
#' @param coexp_focal A numeric matrix of gene-gene correlations (co-expression)
#'   in the focal cell type. Row and column names should be gene identifiers.
#' @param coexp_other_list A list of numeric matrices with gene-gene correlations
#'   in other cell types/clusters. Each matrix should have the same dimensions
#'   and gene ordering as \code{coexp_focal}.
#'
#' @return A list containing four Fdifferent specificity score matrices:
#'   \describe{
#'     \item{ratio}{Ratio of focal co-expression to mean co-expression in other cell types}
#'     \item{zscore}{Z-score measure comparing focal to other cell types}
#'     \item{diff}{Simple difference between focal and mean of other cell types}
#'     \item{combined}{Product of focal co-expression and (1 - other co-expression)}
#'   }
#'
#' @details
#' This function calculates multiple specificity metrics that highlight gene-gene
#' connections that are specific to a focal cell type:
#'
#' 1. **Ratio** – The ratio of absolute co-expression in the focal cell type to
#'    the mean absolute co-expression across other cell types:
#'    \deqn{ ratio = |r_{focal}| / (mean(|r_{other}|) + 0.01) }
#'    Higher values indicate that a gene pair is more strongly co-expressed in
#'    the focal cell type than elsewhere. This score penalizes high co-expression
#     elsewhereis via the denominator, therefore, it is very sensitive to
#'    very low background co-expression (division amplifies small denominators).
#'
#' 2. **Z-score** – A standardized measure showing how many standard deviations
#'    the focal co-expression differs from the mean across other cell types:
#'    \deqn{ zscore = (|r_{focal}| - mean(|r_{other}|)) / (sd(|r_{other}|) + 0.01) }
#'    The z-score measures how far the focal’s co-expression deviates from
#'    the cross-type average in units of standard deviation. Therefore, it
#'    highlights edges whose co-expression is unusually strong (or weak) in the focal.
#'
#' 3. **Difference** – The simple difference between focal and mean co-expression:
#'    \deqn{ diff = |r_{focal}| - mean(|r_{other}|) }
#'    Positive values indicate stronger co-expression in the focal cell type.
#' 	  It gives credit as long as the focal is stronger than the mean, good for
#'    identifying edges that are stronger in the focal than others (quantitative gain).
#'
#' 4. **Combined** – A multiplicative metric rewarding high co-expression in the focal
#'    and low co-expression in other cell types:
#'    \deqn{ combined = |r_{focal}| * (1 - mean(|r_{other}|)) }
#'    This score Rewards high co-expression in focal and explicitly penalizes high
#'    co-expression elsewhere through (1 - mean).
#'    Emphasizes edges that are both strong and specific to the focal cell type.
#'
#' All correlations are converted to absolute values before comparison to capture
#' connection strength regardless of direction (positive or negative correlation).
#'
#' @examples
#' \dontrun{
#' # Create example correlation matrices
#' genes <- c("GENE1", "GENE2", "GENE3", "GENE4")
#' focal_matrix <- matrix(runif(16), 4, 4)
#' rownames(focal_matrix) <- colnames(focal_matrix) <- genes
#'
#' other_matrix1 <- matrix(runif(16), 4, 4)
#' rownames(other_matrix1) <- colnames(other_matrix1) <- genes
#'
#' other_matrix2 <- matrix(runif(16), 4, 4)
#' rownames(other_matrix2) <- colnames(other_matrix2) <- genes
#'
#' # Calculate specificity
#' spec_scores <- calculate_specificity(
#'     focal_matrix,
#'     list(other_matrix1, other_matrix2)
#' )
#' }
#'
#' @export
#' @author Holly Yang \email{xyang2.at.uchicago.edu}

calculate_specificity <- function(coexp_focal, coexp_other_list) {
    # Convert correlation to absolute value
    coexp_focal_abs <- abs(coexp_focal)

    # Calculate mean and standard deviation of co-expression in other cell types
    other_coexp_mean <- matrix(0, nrow = nrow(coexp_focal), ncol = ncol(coexp_focal))
    other_coexp_sd <- matrix(0, nrow = nrow(coexp_focal), ncol = ncol(coexp_focal))
    dimnames(other_coexp_mean) <- dimnames(coexp_focal)
    dimnames(other_coexp_sd) <- dimnames(coexp_focal)

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
    specificity_ratio <- coexp_focal_abs / (other_coexp_mean + 0.01)

    # Option 2: Z-score like measure
    specificity_zscore <- (coexp_focal_abs - other_coexp_mean) / (other_coexp_sd + 0.01)

    # Option 3: Difference-based measure
    specificity_diff <- coexp_focal_abs - other_coexp_mean

    # Option 4: Combined measure (rewards high co-expression in focal and low in others)
    specificity_combined <- coexp_focal_abs * (1 - other_coexp_mean)

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
#' @param cluster_id Character or numeric identifier for the focal cluster
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

    cor_mat[is.na(cor_mat)] <- 0 # v8 added
    return(cor_mat)
}


#' Calculate Cell-Type Specificity Scores for Gene Networks
#'
#' @description
#' This function computes cell-type specificity scores for genes in each network.
#' It calculates co-expression within the focal cell type and compares it to
#' co-expression in other cell types to identify connections that are specific
#' to the focal cell type.
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
    # Rename 'strand' in mcols(rowRanges(sce)) to avoid collision with GRanges reserved slot
    if (!is.null(rowRanges(sce)) && ncol(mcols(rowRanges(sce))) > 0) {
        mc <- mcols(rowRanges(sce))
        if ('strand' %in% colnames(mc)) {
            colnames(mc)[colnames(mc) == 'strand'] <- 'gene_strand'
            mcols(rowRanges(sce)) <- mc
        }
    }

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
        if (ecount(g) < min_n_Eg) { # v8 changed
            cat("Skipping small network:", net_name, "with only", ecount(g), "edges\n")
            next
        }

        # Extract cluster ID from network name
        focal_cluster <- get_cluster_id(net_name)
        focal_cluster <- gsub("\\.\\d+$", "", focal_cluster)
        if (verbose) cat("Focal cluster for this network:", focal_cluster, "\n")

        # Check if this cluster exists in the SCE
        if (!focal_cluster %in% all_clusters) {
            warning("Cluster '", focal_cluster, "' not found in SCE for network '", net_name, "'. Skipping.")
            next
        }

        # Define other clusters as all except the focal
        other_clusters <- setdiff(all_clusters, focal_cluster)

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

        # Calculate co-expression for focal cluster using the independent function
        if (verbose) cat("Computing co-expression for focal cluster:", focal_cluster, "\n")
        coexp_focal <- compute_cluster_correlation(
            sce = sce,
            cluster_id = focal_cluster,
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

        scores <- calculate_specificity(coexp_focal, coexp_other_list)

        # Store results
        specificity_scores[[net_name]] <- list(
            scores = scores,
            genes = genes_in_sce, # Store genes for which scores were calculated
            coexp_focal = coexp_focal,
            corexp_sign = matrix(
                ifelse(coexp_focal > 0, "positive", "negative"),
                nrow = nrow(coexp_focal),
                ncol = ncol(coexp_focal),
                dimnames = dimnames(coexp_focal)
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
#'   \item \code{coexp_focal}: Cell-type specific co-expression scores used for weighting
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
#'     # [1] "weight"          "norm_PPI_score" "corexp_sign"     "coexp_focal"
#'
#'     plot(E(updated_networks[["A"]])$coexp_focal, E(updated_networks[["A"]])$weight,
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
        coexp_focal_vals <- rep(NA_real_, nE)
        new_weights <- E(g)$norm_PPI_score # default: keep old weights

        # Fill for valid edges
        corexp_vals[valid] <- as.character(spec_data$corexp_sign[idx])
        coexp_focal_vals[valid] <- specificity_norm[idx]

        # Compute combined weights
        candidate_weights <- (string_norm[idx] * specificity_norm[idx])

        # Only accept positive weights
        keep <- which(valid)[candidate_weights > 0 & !is.na(candidate_weights)]

        new_weights[keep] <- candidate_weights[candidate_weights > 0 & !is.na(candidate_weights)] # v9 changed

        # Assign in bulk
        g <- set_edge_attr(g, "corexp_sign", value = corexp_vals)
        g <- set_edge_attr(g, "coexp_focal", value = coexp_focal_vals)
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
        betw_scores <- igraph::betweenness(g, weights = E(g)$weight) ## weight should be distance-based weight !!
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
        E(g)$weight <- 1 / E(g)$weight
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
        out[, graph_name := graph_attr(g, "name")] # v11
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
#' @author X. Yang
#' @export
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
            cluster_ID = unlist(strsplit(graph_name, split = "_"))[2]
        )
        all_edge_data <- rbind(all_edge_data, edge_data)
    }
    # Filter to main categories
    all_edge_data <- all_edge_data[all_edge_data$PPI_cat %in% c("CTS", "HiGCTS", "HiG"), ]
    all_edge_data$PPI_cat <- factor(all_edge_data$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))
    all_edge_data$cluster_cat <- ifelse(all_edge_data$cluster_ID %in% unstable_cluster_ID, "unstable", "stable")

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
#' @author X. Yang
#' @export
#'
plot_edge_weight_distributions <- function(edge_data, PPI_color_palette) {
    # Calculate summary statistics
    summary_stats <- edge_data %>%
        group_by(PPI_cat) %>%
        summarise(
            mean_weight = round(mean(edge_weight), 3),
            median_weight = round(median(edge_weight), 3),
            total_edges = n(),
            .groups = "drop"
        )
    # 1. Density plot for HiG PPIs, colored by cluster_categories (stable or instable)
    p1 <- ggplot(edge_data, aes(
        x = log10(edge_weight + 1), # log10 transform inside aes()
        fill = PPI_cat,
        color = PPI_cat
    )) +
        geom_density(alpha = 0.3, size = 1.2, adjust = 1.2) + # adjust controls smoothness
        # facet_wrap(~PPI_cat, ncol = 1, scales = "free_y") +
        labs(
            title = "hiPSC",
            x = expression(log[10] ~ "(PPI edge weight + 1)"),
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
        stat_compare_means(
            comparisons = list(c("HiG", "HiGCTS"), c("HiG", "CTS"), c("HiGCTS", "CTS")),
            method = "wilcox.test",
            label = "p.signif",
            p.adjust.method = "bonferroni"
        ) +
        labs(
            title = "wilcox.test",
            subtitle = paste(
                "Total edges - CTS:", summary_stats$total_edges[1],
                " | HiGCTS:", summary_stats$total_edges[2],
                " | HiG:", summary_stats$total_edges[3]
            ),
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
#' @author X. Yang
#' @export
#'
plot_nEB_ggplot <- function(nEB_data, PPI_color_palette, method = c("t.test", "wilcox.test")) {
    method <- match.arg(method)
    # Convert to data frame
    plot_data <- data.frame(
        sample = factor(names(nEB_data), levels = names(nEB_data)),
        value = as.numeric(nEB_data),
        index = 1:length(nEB_data),
        stringsAsFactors = FALSE # not to automatically convert character strings into factors when creating a data frame.
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
        TRUE ~ plot_data$sample # catches anything that didn't match the previous conditions
    )
    # Get HiG samples as the base x-axis (13 positions)
    hig_data <- plot_data[plot_data$PPI_cat == "HiG", ]
    hig_data <- hig_data[order(hig_data$value), ] # Order by value (smallest to largest)  # Maintain original order
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
        geom_line(
            data = hig_data, aes(x = x_position, y = value),
            color = PPI_color_palette["HiG"], size = 1.2, alpha = 0.8
        ) +
        # Add HiG points
        geom_point(
            data = hig_data, aes(x = x_position, y = value, color = PPI_cat),
            size = 3, alpha = 0.9
        ) +
        # Add overlaid points for other categories
        geom_point(
            data = overlay_data, aes(x = x_position, y = value, color = PPI_cat),
            size = 3, alpha = 0.9
        ) +
        labs(
            title = "cluster_edge_betweenness",
            subtitle = paste("HiG line (n=", sum(plot_data$PPI_cat == "HiG"),
                ") | HiGCTS dots (n=", sum(plot_data$PPI_cat == "HiGCTS"),
                ") | CTS dots (n=", sum(plot_data$PPI_cat == "CTS"), ")",
                sep = ""
            ),
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
    boxplot_data <- plot_data[plot_data$PPI_cat %in% c("CTS", "HiGCTS", "HiG") &
        plot_data$sample_type %in% overlay_data$sample_type, ]
    # Set factor levels to control order
    boxplot_data$PPI_cat <- factor(boxplot_data$PPI_cat, levels = c("CTS", "HiGCTS", "HiG"))

    comparisons <- list(c("HiG", "HiGCTS"), c("HiG", "CTS"), c("HiGCTS", "CTS"))
    # Calculate summary statistics for subtitle
    summary_stats <- boxplot_data %>%
        group_by(PPI_cat) %>%
        summarise(
            mean_val = round(mean(value), 1),
            median_val = round(median(value), 1),
            .groups = "drop"
        )

    p2 <- ggplot(boxplot_data, aes(x = PPI_cat, y = value, fill = PPI_cat)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
        geom_jitter(aes(color = PPI_cat), width = 0.2, size = 2, alpha = 0.8) +
        labs(
            title = paste(method, "test, transition clusters"),
            subtitle = paste(
                "Mean/Median - CTS:", summary_stats$mean_val[1], "/", summary_stats$median_val[1],
                " | HiGCTS:", summary_stats$mean_val[2], "/", summary_stats$median_val[2],
                " | HiG:", summary_stats$mean_val[3], "/", summary_stats$median_val[3]
            ),
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
            legend.position = "none", # Remove legend since colors match categories
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        stat_compare_means(
            comparisons = comparisons,
            method = method,
            p.adjust.method = "none", # "bonferroni",
            label = "p.signif",
            bracket.size = 0.6,
            step.increase = 0.1
        )

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
#' 	 \item{\code{network_colors}}{named color palette used here.}
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
#' g_real <- graph_CTS_CP.1 # graph_list[["CTS_CP.1"]]
#' result <- synthetic_simulation(g_real, main = "CTS_CP.1")
#' result$auc_summary
#' result$p_weights + result$p_line + result$p_AUC
#'
#' @import igraph data.table ggplot2 pracma
#' @author X. Yang
#' @export
#'
synthetic_simulation <- function(g_real, main = NULL) {
    if (is.null(main)) main <- "Network"

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
    p_edge <- ecount(g_real) / choose(n_nodes, 2)

    # 1. random (Erdős–Rényi)
    g_random <- sample_gnp(n = n_nodes, p = p_edge)
    # 2. scale-free (preferential attachment)
    m_links <- round(ecount(g_real) / n_nodes) # average edges per new node
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
    E(g_random)$weight <- sample(emp_weights, ecount(g_random), replace = TRUE)
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
            times = c(length(emp_weights), ecount(g_random))
        )
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

    sim_real <- robustness_MonteCarlo(g_real, type = "vertex", measure = "btwn.cent", N = 100)
    sim_random <- robustness_MonteCarlo(g_random, type = "vertex", measure = "btwn.cent", N = 100)
    sim_scale <- robustness_MonteCarlo(g_scale_free, type = "vertex", measure = "btwn.cent", N = 100)
    sim_small <- robustness_MonteCarlo(g_small_world, type = "vertex", measure = "btwn.cent", N = 100)
    sim_deg <- robustness_MonteCarlo(g_deg_preserved, type = "vertex", measure = "btwn.cent", N = 100)

    ################################################
    # Step 6 -  compare AUCs

    auc_real <- trapz(sim_real$removed.pct, sim_real$comp.pct)
    auc_random <- trapz(sim_random$removed.pct, sim_random$comp.pct)
    auc_scale <- trapz(sim_scale$removed.pct, sim_scale$comp.pct)
    auc_small <- trapz(sim_small$removed.pct, sim_small$comp.pct)
    auc_deg <- trapz(sim_deg$removed.pct, sim_deg$comp.pct)

    # Combine and summarize AUC results
    auc_summary <- data.frame(
        network = c("real_PPIN", "random", "scale_free", "small_world", "deg_preserving"),
        AUC = c(auc_real, auc_random, auc_scale, auc_small, auc_deg)
    )

    auc_summary$normalized_AUC <- auc_summary$AUC / max(auc_summary$AUC)
    auc_summary
    # network       AUC normalized_AUC
    # 1   real_PPIN 0.2327409      0.2342621
    # 2      random 0.9935065      1.0000000
    # 3  scale_free 0.9935065      1.0000000
    # 4 small_world 0.9935065      1.0000000
    # 5 deg_preserving 0.1541933      0.1552011
    auc_summary$ID <- main

    ################################################
    # Step 7 — Visualize fragmentation curves

    # Combine all simulation outputs
    sim_all <- rbindlist(list(
        cbind(sim_real, network = "real_PPIN"),
        cbind(sim_random, network = "random"),
        cbind(sim_scale, network = "scale_free"),
        cbind(sim_small, network = "small_world"),
        cbind(sim_deg, network = "deg_preserving")
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
    sim_all$network <- factor(sim_all$network, levels = network_levels)
    auc_summary$network <- factor(auc_summary$network, levels = network_levels)

    p_line <- ggplot(sim_all, aes(x = removed.pct, y = comp.pct, color = network)) +
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
    p_AUC <- ggplot(auc_summary, aes(x = network, y = AUC, fill = network)) +
        geom_col(width = 0.7) +
        theme_minimal(base_size = 12) +
        labs(
            title = paste(main, "resilience (AUC of fragmentation curve)"),
            x = "Network type",
            y = "AUC (higher = more robust)"
        ) +
        scale_color_manual(values = network_colors) +
        scale_fill_manual(values = network_colors) +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)
        )


    ################################################
    # Step 8 — return
    outputs <- list(
        auc_summary = auc_summary,
        p_weights = p_weights,
        p_line = p_line,
        p_AUC = p_AUC,
        network_colors = network_colors
    )

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
#' @author X. Yang; Felix Yu
#' @export
#'
plot_weighted_PPIN <- function(g, layout = "fr", CHD = NULL, node_size_title = "|Wilcox score|") {
    # ---- Attribute checks ----
    if (!"weight" %in% vertex_attr_names(g)) {
        stop("Please assign V(g)$weight for plotting node size")
    }

    if (!"weight" %in% edge_attr_names(g)) {
        stop("Please assign E(g)$weight for plotting edge width")
    }

    if (!"corexp_sign" %in% edge_attr_names(g)) {
        stop("Please assign E(g)$corexp_sign for plotting edge color")
    }

    # ---- Validate edge signs ----
    expected_signs <- c("positive", "negative")
    missing_signs <- setdiff(expected_signs, unique(E(g)$corexp_sign))
    if (length(missing_signs) > 1) {
        stop("Edge attribute 'corexp_sign' does not include expected values: 'positive', 'negative'")
    }

    # ---- Disease-gene highlight ----
    if (is.null(CHD)) {
        warning("No disease genes supplied for highlighting (CHD = NULL)")
        V(g)$is_CHD <- FALSE
    } else {
        V(g)$is_CHD <- toupper(V(g)$name) %in% CHD
    }

    # ---- Determine if this is a CTS network ----
    # Only CTS networks will use HiG for shape mapping
    graph_name <- graph_attr(g, "name") # if you use graph_attr, or pass the network ID externally
    use_HiG_shape <- !is.null(graph_name) && grepl("^CTS_", graph_name)

    # ---- Layout ----
    set.seed(1234)
    layout_coords <- create_layout(g, layout = layout)

    # ---- Start plotting ----
    p <- ggraph(layout_coords) +
        # edges
        geom_edge_link(aes(
            width = weight,
            color = corexp_sign
        ), alpha = 0.7) +

        # nodes
        geom_node_point(aes(
            size = weight,
            color = is_CHD,
            shape = if (use_HiG_shape) is_HiG else NULL # only map shape for CTS
        )) +

        # node labels
        geom_node_text(aes(label = name), repel = TRUE, size = 3) +

        # color scales
        scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "gray70")) +
        scale_edge_color_manual(values = c("positive" = "orange", "negative" = "blue")) +

        # node size / edge width
        scale_size_continuous(range = c(2, 5), name = node_size_title) +
        scale_edge_width_continuous(
            range = c(0.1, 1),
            limits = range(E(g)$weight, na.rm = TRUE),
            name = "E weights"
        ) +

        # theme
        theme_void() +
        labs(edge_color = "Specific co-exp", color = "Curated CHD genes") +
        theme(legend.position = "right") +
        ggtitle(paste0(vcount(g), " PPI genes"))

    # ---- Conditionally add shape scale for CTS only ----
    if (use_HiG_shape) {
        p <- p + scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), name = "is_HiG")
    }

    return(p)
}

#' PPIN-based gene set enrichment analysis for a key gene
#'
#' @description
#' Performs enrichment analysis of genes connected to a given key transcription factor
#' or regulator (e.g., ISL1, or the gene with the top degree) within a protein-protein interaction (PPIN) or co-expression
#' graph, testing overlap with a list of predefined gene sets (e.g., CHD, DEGs).
#'
#' The function tests two levels:
#' 1. **Direct module** – genes directly connected to the key gene.
#' 2. **Connected component** – all genes in the same connected component as the key gene.
#'
#' For each level, Fisher’s exact tests are performed to assess enrichment of overlap
#' with each provided gene set. Both raw and Benjamini–Hochberg–adjusted p-values
#' are returned, along with visualization plots.
#'
#' @param g An `igraph` object representing the protein–protein interaction or co-expression network.
#' @param gene_sets A named list of gene sets (each a character vector of gene symbols).
#'   The order of gene sets in this list will be preserved when plot enrichment dot plots.
#'   Example: `list(ISL_act_E = c("MYL7","TNNI1"), CHD = c("NKX2-5","TBX5"))`.
#' @param key_gene Character string specifying the hub or query gene of interest.
#'   If `NULL` (default) or not found in the network, the function automatically
#'   selects the gene with the highest degree (i.e., the most connected node) in `g`
#'   as the key gene. This ensures the analysis can proceed even if the specified
#'   `key_gene` is absent or unspecified.
#' @param graph_name Optional character string specifying the network name for figure titles
#'   (default: `"HiGCTS_CP.1"`).
#' @param GS_database Optional character string specifying the source for input ('gene_sets')
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item `enrichment_df`: Data frame of enrichment results for directly connected genes, and the all network member genes.
#'   \item `key_gene_module`: Character vector of directly connected genes (including the key gene).
#'   \item `key_gene_network`: Character vector of PPIN member genes of the input ('g').
#'   \item `key_gene`: Character vector of the analyzed gene of interest, eith the input ('key_gene') or the hub with highest degree.
#'   \item `p_module`: `ggplot` object showing enrichment (-log10 raw p-values) for the direct module.
#'   \item `p_network`: `ggplot` object showing enrichment (-log10 raw p-values) for the igraph member genes.
#'   \item `p_dot`: `ggplot` object showing dot plots for borh -log10 p-values but als odds ratios.
#' }
#'
#' @details
#' - Each enrichment test uses a one-sided Fisher’s exact test (`alternative = "greater"`).
#' - The background universe includes all unique genes from the network and all provided gene sets.
#' - Adjusted p-values (`padj`) are calculated using the Benjamini–Hochberg (BH) method.
#' - Red dashed line in plots indicates p = 0.05.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_ring(10)
#' V(g)$name <- paste0("Gene", 1:10)
#' gene_sets <- list(
#'     SetA = c("Gene1", "Gene2", "Gene3"),
#'     SetB = c("Gene4", "Gene5", "Gene6")
#' )
#' res <- PPIN_geneset_enrichment(g, gene_sets, key_gene = "Gene1", graph_name = "ExampleGraph")
#'
#' # View results
#' res$enrichment_df
#'
#' # Plot directly connected gene enrichment
#' print(res$p_module)
#'
#' # Plot network-level enrichment
#' print(res$p_network)
#'
#' # Plot for both with odds ratios
#' print(res$p_dot)
#' }
#' @author X. Yang
#' @export
#'
PPIN_geneset_enrichment <- function(g, gene_sets, key_gene = "ISL1", graph_name = "HiGCTS_CP.1", GS_database = "genesets") {
    library(ggplot2)

    if (is.null(key_gene)) {
        d <- degree(g, mode = "all")
        key_gene <- names(d)[which.max(d)]
        message("testing for the top-degree gene: ", key_gene)
    } else if (!key_gene %in% igraph::V(g)$name) {
        warning("Given Key gene is not in the network, using NULL to use the top degree gene")
    }

    universe <- unique(unlist(c(gene_sets, list(PPIN = igraph::V(g)$name))))

    # level 1: direct connected genes of key_gene
    key_gene_neighbors <- igraph::neighbors(g, key_gene, mode = "all") |> names()
    key_gene_module <- unique(c(key_gene, key_gene_neighbors))
    length(key_gene_module) # 9
    module_enrichment <- lapply(names(gene_sets), function(set_name) {
        gs <- gene_sets[[set_name]]

        in_module <- universe %in% key_gene_module
        in_set <- universe %in% gs

        Intersect_gene <- universe[in_module & in_set]

        tbl <- table(in_module, in_set)
        ft <- fisher.test(tbl, alternative = "greater")

        data.frame(
            geneset = set_name,
            Intersect_Count = length(Intersect_gene),
            Intersect_gene = paste(Intersect_gene, collapse = ", "),
            module_size = sum(in_module),
            set_size = sum(in_set),
            p_value = ft$p.value,
            Fisher_odds = ft$estimate
        )
    })

    module_enrichment_df <- do.call(rbind, module_enrichment)
    module_enrichment_df$padj <- p.adjust(module_enrichment_df$p_value, method = "BH")
    module_enrichment_df <- module_enrichment_df[order(module_enrichment_df$padj), ]

    ## level 3: all gene members in the network
    key_gene_network <- igraph::V(g)$name
    network_enrichment <- lapply(names(gene_sets), function(set_name) {
        gs <- gene_sets[[set_name]]

        in_network <- universe %in% key_gene_network
        in_set <- universe %in% gs

        Intersect_gene <- universe[in_network & in_set]

        tbl <- table(in_network, in_set)
        ft <- fisher.test(tbl, alternative = "greater")

        data.frame(
            geneset = set_name,
            Intersect_Count = length(Intersect_gene),
            Intersect_gene = paste(Intersect_gene, collapse = ", "),
            module_size = sum(in_network),
            set_size = sum(in_set),
            p_value = ft$p.value,
            Fisher_odds = ft$estimate
        )
    })

    network_enrichment_df <- do.call(rbind, network_enrichment)
    network_enrichment_df$padj <- p.adjust(network_enrichment_df$p_value, method = "BH")
    network_enrichment_df <- network_enrichment_df[order(network_enrichment_df$padj), ]

    # step 3: visualization
    module_enrichment_df$GS_category <- "key_neighbor"
    network_enrichment_df$GS_category <- "network_gene"
    enrichment_df <- rbind(module_enrichment_df, network_enrichment_df) %>%
        mutate(
            GS_category = factor(GS_category, levels = c("key_neighbor", "network_gene")), # preserve row order
            geneset = factor(geneset, levels = names(gene_sets)), # preserve column order
            log10p = -log10(p_value), # compute -log10(p)
            sig = ifelse(p_value < 0.05 & Fisher_odds > 2, "significant", "nonsignificant"),
        )

    # Dot plot
    p_dot <- GS_enrichment_Dotplot(enrichment_df,
        GS_database = GS_database,
        sig_p_adjust = FALSE, reverse_order = FALSE, n_gene_characters = 0
    )

    p_module <- ggplot(module_enrichment_df, aes(x = reorder(geneset, -log10(p_value)), y = -log10(p_value))) +
        geom_bar(stat = "identity", fill = "#407bc2") +
        geom_text(
            aes(label = Intersect_gene),
            vjust = -0.4, size = 3, lineheight = 0.9
        ) +
        labs(
            x = "Gene Set",
            y = expression(-log[10](p_value)),
            title = paste0("Enrichment of ", length(key_gene_module), " directly ", key_gene, "-interactive genes in ", graph_name)
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8)

    p_network <- ggplot(network_enrichment_df, aes(x = reorder(geneset, -log10(p_value)), y = -log10(p_value))) +
        geom_bar(stat = "identity", fill = "#407bc2") +
        geom_text(
            aes(label = Intersect_gene),
            vjust = -0.4, size = 3, lineheight = 0.9
        ) +
        labs(
            x = "Gene Set",
            y = expression(-log[10](p_value)),
            title = paste0("Enrichment of ", length(key_gene_network), " genes in ", graph_name)
        ) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8)

    return(list(
        enrichment_df = enrichment_df,
        key_gene_module = key_gene_module,
        key_gene_network = key_gene_network,
        p_module = p_module,
        p_network = p_network,
        p_dot = p_dot,
        key_gene = key_gene
    ))
}


#' @title GS_enrichment_Dotplot
#' @description
#' Generate a dot plot visualization of gene set enrichment results, where each point represents
#' one enriched gene set. The plot displays odds ratio on the x-axis, gene set name (and intersecting genes)
#' on the y-axis, with point size indicating the number of overlapping genes and color representing
#' statistical significance (-log10 FDR). Only significant dot (FDR<0.05 & odds ratio>2) are colored.
#'
#' @details
#' This function assumes that the **first column** of the input dataframe (`df`) contains
#' the gene set identifier or name. The dataframe must also include the following columns:
#' \itemize{
#'   \item \strong{Fisher_odds}: Numeric value of the enrichment odds ratio.
#'   \item \strong{FDR}: Adjusted p-value (False Discovery Rate).
#'   \item \strong{Intersect_Count}: Number of overlapping (intersecting) genes.
#'   \item \strong{Intersect_gene}: Character string listing intersecting genes,
#'         separated by spaces or commas.
#' }
#'
#' Each y-axis label includes both the gene set ID and the corresponding intersecting genes
#' on a second line for clarity.
#'
#' @param df A dataframe containing gene set enrichment results.
#'   The first column must be the gene set identifier. Required columns:
#'   \code{Fisher_odds}, \code{FDR}, \code{Intersect_Count}, and \code{Intersect_gene}.
#' @param GS_database Character string specifying the name of the gene set database
#'   (e.g. "Msigdb.c2.cp", "Reactome", "WikiPathways"). Default is \code{"Msigdb.c2.cp"}.
#' @param sig_p_adjust logic value, default is TRUE to call significance based on the multi-test adjusted p-value
#' @param reverse_order logic value, default is TRUE to reverse for top-down order in the dot plot
#' @param n_gene_characters  Integer, the character number of intersect gene symbols per gene_set to be shown
#' @param GS_category an optinal character indicate the column of the input dataframe (`df`)
#'
#' @return
#' A \code{ggplot} object representing the enrichment dot plot, where:
#' \itemize{
#'   \item x-axis: Odds ratio (\code{Fisher_odds})
#'   \item y-axis: Gene set name and intersecting genes
#'   \item Point color: -log10(FDR)
#'   \item Point size: Number of intersecting genes
#' }
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'     GSID = c("WP_HEART_DEVELOPMENT", "WP_NEURAL_CREST_DIFFERENTIATION"),
#'     Fisher_odds = c(44.96, 19.12),
#'     FDR = c(0.00086, 0.00476),
#'     Intersect_Count = c(3, 3),
#'     Intersect_gene = c("HAND2 ISL1 FGF10", "ISL1 ID1 PRTG")
#' )
#' p <- GS_enrichment_Dotplot(df, GS_database = "WikiPathways")
#' print(p)
#' }
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[dplyr]{mutate}}
#'
#' @author X. Yang
#' @export
#'
GS_enrichment_Dotplot <- function(df, GS_database = "Msigdb.c2.cp", sig_p_adjust = TRUE, reverse_order = TRUE,
                                  n_gene_characters = 30) {
    library(ggplot2)
    library(dplyr)
    colnames(df)[1] <- "GSID"

    if (is.factor(df$GSID)) geneset_levels <- levels(df$GSID) else geneset_levels <- unique(df$GSID)

    if (sig_p_adjust) {
        if (!all(c("Fisher_odds", "FDR", "Intersect_Count", "Intersect_gene") %in% colnames(df))) {
            stop("Required columns 'Fisher_odds','FDR','Intersect_Count' and/or 'Intersect_gene' are missing.")
        }
    } else {
        if (!all(c("Fisher_odds", "p_value", "Intersect_Count", "Intersect_gene") %in% colnames(df))) {
            stop("Required columns 'Fisher_odds','p_value','Intersect_Count' and/or 'Intersect_gene' are missing.")
        }
    }

    if (n_gene_characters > 0) {
        df$y_label <- paste0(df$GSID, "\n(", substr(df$Intersect_gene, 1, n_gene_characters), ")")
    } else {
        df$y_label <- df$GSID
    }
    if (reverse_order) df$y_label <- factor(df$y_label, levels = rev(df$y_label)) # reverse for top-down order

    if (sig_p_adjust) {
        df$sig <- ifelse(df$FDR < 0.05 & df$Fisher_odds > 2, "significant", "nonsignificant")
    } else {
        df$sig <- ifelse(df$p_value < 0.05 & df$Fisher_odds > 2, "significant", "nonsignificant")
    }

    if (is.null(GS_database)) GS_database <- "geneset"
    if (sig_p_adjust) {
        p <- ggplot(df, aes(x = Fisher_odds, y = y_label)) +
            geom_point(aes(
                size = Intersect_Count, color = -log10(FDR),
                shape = sig, fill = ifelse(sig == "significant", -log10(FDR), NA)
            ), stroke = 1.1) +
            scale_color_gradient(low = "blue", high = "red", name = "-log10FDR")
    } else {
        p <- ggplot(df, aes(x = Fisher_odds, y = y_label)) +
            geom_point(aes(
                size = Intersect_Count, color = -log10(p_value),
                shape = sig, fill = ifelse(sig == "significant", -log10(p_value), NA)
            ), stroke = 1.1) +
            scale_color_gradient(low = "blue", high = "red", name = "-log10P")
    }

    p <- p +
        scale_shape_manual(values = c("significant" = 21, "nonsignificant" = 21)) +
        scale_fill_gradient(low = "white", high = "red", na.value = "white") +
        scale_size_continuous(name = "Intersect count", range = c(4, 10)) +
        guides(fill = "none", shape = "none") +
        labs(
            x = "Odds ratio", y = NULL,
            title = paste(GS_database, "enrichment")
        ) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.y = element_text(size = 10, lineheight = 0.9),
            axis.text.x = element_text(size = 11),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank(),
            legend.title = element_text(size = 9), # smaller legend title
            legend.text = element_text(size = 8) # smaller legend labels
        )

    if ("GS_category" %in% colnames(df)) p <- p + facet_wrap(~GS_category, scales = "free", ncol = 2)

    return(p)
}


#' @title Volcano plot highlighting specific genes
#'
#' @description
#' Generates a volcano plot from differential expression data (typically RNA-seq or scRNA-seq results)
#' and highlights user-specified genes of interest.
#' The function plots log2 fold changes (x-axis) versus –log10 adjusted p-values (y-axis),
#' categorizes genes as upregulated, downregulated, or non-significant,
#' and visually emphasizes selected genes (e.g., known regulators or markers).
#'
#' @param highlight_genes Character vector of gene symbols to highlight on the volcano plot.
#'   Default includes a panel of cardiac developmental regulators:
#'   \code{c('ISL1', 'GATA6', 'FGF10', 'HAND2', 'CITED2', 'HAS2',
#'   'GATA4', 'MEF2C', 'MSX2', 'TWIST1', 'FGF8', 'IRX3', 'ALX1', 'IGFBPL1')}.
#' @param df A data frame containing differential expression results.
#'   The first column should contain gene symbols, and the table must include
#'   columns named \code{"avg_log2FC"} (log2 fold change) and \code{"p_val_adj"} (adjusted p-value or FDR).
#' @param logFC_cut Numeric value defining the absolute log2 fold-change threshold for significance.
#'   Genes with \code{|avg_log2FC| > logFC_cut} and \code{p_val_adj < p_adj.cut}
#'   are considered significantly differentially expressed. Default = 0.25.
#' @param p_adj.cut Numeric value defining the adjusted p-value (FDR) cutoff for significance. Default = 0.05.
#' @param unique_symbol logical value, if TRUE, select the one with largest(abs(pi)) for transcript of the same gene symbol to highligh
#' @param main Character string specifying the main plot title. If \code{NULL}, defaults to "Volcano plot".
#'
#' @details
#' The function automatically categorizes genes into "Up", "Down", and "NotSig" in Isl1KO reported by Maven2023 table S1, filtered by
#' abs(log₂FC) > 0.25 & FDR < 0.05
#' based on user-defined thresholds, and visually distinguishes them using color:
#' \itemize{
#'   \item Upregulated genes in Isl1KO — red (\code{"#E64B35"})
#''   \item Downregulated genes in Isl1KO — blue (\code{"#4DBBD5"})
#'   \item Non-significant genes in Isl1KO — light gray (\code{"grey80"})
#' }
#' Highlighted genes are plotted in orange and labeled using \code{ggrepel}.
#'
#' @return
#' A \code{ggplot2} object representing the volcano plot. The plot can be
#' further customized, saved, or integrated into multi-panel figures.
#'
#' @examples
#' \dontrun{
#' # Example using DEGs_E (with columns: gene, avg_log2FC, p_val_adj)
#' p <- Volcano_for_highlight_genes(
#'     df = DEGs_E,
#'     highlight_genes = c("ISL1", "GATA6", "HAND2", "MEF2C"),
#'     logFC_cut = 0.25,
#'     p_adj.cut = 0.05,
#'     main = "ISL1-dependent DEGs in Early CPs"
#' )
#'
#' # Display or save
#' print(p)
#' ggsave("Volcano_ISL1_EarlyCP.pdf", p, width = 5, height = 4)
#' }
#'
#' @seealso
#' \code{\link[ggplot2]{ggplot}}, \code{\link[ggrepel]{geom_text_repel}}
#'
#' @author X. Yang
#' @export
#'
Volcano_for_highlight_genes <- function(highlight_genes = c(
                                            "ISL1", "GATA6", "FGF10", "HAND2", "CITED2",
                                            "HAS2", "GATA4", "MEF2C", "MSX2", "TWIST1",
                                            "FGF8", "IRX3", "ALX1", "IGFBPL1"
                                        ), df,
                                        logFC_cut = 0.25, p_adj.cut = 0.05, unique_symbol = FALSE, main = NA) {
    if (is.null(main)) main <- "Volcano plot"
    # Ensure column names match
    if (colnames(df)[1] != "gene") {
        colnames(df)[1] <- "gene"
        warning("Using the 1st column of df as gene symbols")
    }
    if (!all(c("avg_log2FC", "p_val_adj") %in% colnames(df))) stop("df must have columns named with'avg_log2FC','p_val_adj'")


    if (unique_symbol) { #  select the one with largest(abs(pi)) for transcript of the same gene symbol to highligh
        df$pi <- -log10(df$p_val_adj) * df$avg_log2FC
        df$abs_pi <- abs(df$pi)
        # For genes with multiple transcripts, keep the one with largest |π|
        df <- df %>%
            group_by(gene) %>%
            slice_max(order_by = abs_pi, n = 1, with_ties = FALSE) %>%
            ungroup()
    }

    # Add significance & highlight flags
    volc_data <- df %>%
        mutate(
            negLog10P = -log10(p_val_adj),
            Significance = case_when(
                avg_log2FC > logFC_cut & p_val_adj < p_adj.cut ~ "Up",
                avg_log2FC < -logFC_cut & p_val_adj < p_adj.cut ~ "Down",
                TRUE ~ "NotSig"
            ),
            Highlight = ifelse(gene %in% highlight_genes, "Yes", "No")
        )

    # Set colors for categories
    volcano_colors <- c("Up" = "#E64B35", "Down" = "#4DBBD5", "NotSig" = "grey80")

    # Volcano plot
    p <- ggplot(volc_data, aes(x = avg_log2FC, y = negLog10P)) +
        geom_point(aes(color = Significance), alpha = 0.7, size = 1.8) +
        scale_color_manual(values = volcano_colors) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.6) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.6) +
        geom_point(
            data = subset(volc_data, Highlight == "Yes"),
            color = "orange", size = 2.5
        ) +
        geom_text_repel(
            data = subset(volc_data, Highlight == "Yes"),
            aes(label = gene),
            size = 3,
            color = "black",
            box.padding = 0.4,
            segment.color = "grey60",
            max.overlaps = 100
        ) +
        theme_classic(base_size = 12) +
        labs(
            x = expression(Log[2] ~ Fold ~ Change),
            y = expression(-Log[10] ~ FDR),
            title = main
        ) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none"
        )

    return(p)
}


################################################ 2026 new functions #############################################
# This function plots gene symbols ranked by pagerank, highlight the CHD in red dot and the TF in red triangle #########################

#' Rank and plot genes across PPIN signatures
#'
#' For each signature, rank genes by betweenness and plot `key` vs rank.
#' CHD genes are highlighted in red; TFs are triangles.
#'
#' @param df Data frame with columns `gene`, `signature`, `BetweennessCentrality`, and `key`.
#' @param CHD,TF_human Character vectors of genes to highlight.
#' @param signatures Signatures to plot (must exist in `df$signature`).
#' @param key Column name to plot (e.g. "PageRank").
#' @param top_TF_rank Numbers of top ranked TFs to label.
#' @param gene_top_n  Numbers of gene_of_interested (eg. TF, CHD) to label
#' @param saveFigure If TRUE, save `CP_rank_gene_by_<key>.pdf`.
#'
#' @return Named list of subset of df for the top TFs (one per signature).
#' @export plot_TF_CHD_in_PPIN, Triangle = TF; circle = non-TF.
#' Red points indicate CHD genes.
#' Labels are shown for TF/CHD genes within the top gene_top_n ranks, plus the top top_TF_rank TFs;
#' label text is red for TFs and black otherwise.
#' @examples
#' \dontrun{
#' set.seed(1)
#' df <- data.frame(
#'     gene = paste0("G", 1:200),
#'     signature = rep(c("HiG_CP", "CTS_CP", "CTS_CP.1", "HiGCTS_CP", "HiGCTS_CP.1"), each = 40),
#'     BetweennessCentrality = rexp(200),
#'     PageRank = runif(200)
#' )
#' CHD <- c("G2", "G50", "G120")
#' TF <- c("G10", "G11", "G12")
#' p_list <- rank_TF_CHD_in_PPIN(
#'     df,
#'     CHD = CHD, TF_human = TF,
#'     key = "PageRank",
#'     saveFigure = FALSE
#' )
#' p_list[["HiG_CP"]]
#' }
#' @export
#'
rank_TF_CHD_in_PPIN <- function(df, CHD, TF_human,
                                signatures = c("HiG_CP", "CTS_CP", "CTS_CP.1", "HiGCTS_CP", "HiGCTS_CP.1"),
                                key = "PageRank", #' BetweennessCentrality'
                                top_TF_rank = 3, gene_top_n = 20, saveFigure = TRUE) {
    # tmp = try(df$signature, silent = TRUE)
    if (is.null(df$signature)) stop("signature column not found in df")
    if (any(!signatures %in% df$signature)) stop("signature not found in df$signature")
    if (!(key %in% colnames(df))) stop("key not found in df")

    int <- union(CHD, TF_human)
    p <- df_TF_top <- list()

    for (signatureID in signatures) {
        df_plot <- subset(df, signature == signatureID) %>%
            # arrange(desc(.data[[key]])) %>%
            arrange(desc(!!sym(key))) %>% ## advanced way
            mutate(rank = row_number()) %>%
            mutate(is_CHD = gene %in% CHD) %>%
            mutate(is_TF = gene %in% TF_human)

        labeled_genes <- union(intersect(df_plot$gene[1:gene_top_n], int), subset(df_plot, is_TF)$gene[1:top_TF_rank])

        p[[signatureID]] <- ggplot(df_plot, aes(x = rank, y = .data[[key]])) +
            # base layer: all genes
            geom_point(aes(shape = is_TF, size = is_TF), color = "black", alpha = 0.9) +
            # overlay: CHD genes in red (keeps TF shape)
            geom_point(
                data = subset(df_plot, is_CHD),
                aes(shape = is_TF, size = is_TF),
                color = "red"
            ) +
            scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
            scale_size_manual(values = c(`FALSE` = 1.0, `TRUE` = 2.2)) +
            labs(x = paste("Rank in", signatureID), y = key, title = signatureID) +
            theme_classic(base_size = 12) +
            ggrepel::geom_text_repel(
                data = subset(df_plot, gene %in% labeled_genes),
                aes(label = gene, color = is_TF),
                size = 3, box.padding = 0.3,
                segment.color = "grey50",
                max.overlaps = 30
            ) +
            scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black")) +
            theme(legend.position = "none")

        df_TF_top[[signatureID]] <- filter(df_plot, is_TF & rank <= gene_top_n)
        if (nrow(df_TF_top[[signatureID]]) > top_TF_rank) df_TF_top[[signatureID]] <- df_TF_top[[signatureID]][1:top_TF_rank, ] # !!!!!!!!!!$gene
    }

    # ---- Correct legend panel (matches plot semantics exactly) ----
    legend_text <- paste0(
        "Point shape/size: TF (triangle, larger) vs non-TF (circle, smaller)\n",
        "Point color: CHD genes are overlaid in red\n",
        "Label color: TF labels in red; non-TF labels in black\n",
        "Labels shown: TF/CHD within top ", gene_top_n, " ranks + top ", top_TF_rank, " TFs"
    )

    p_legend <- ggplot() +
        annotate("point", x = 0.08, y = 0.78, shape = 17, colour = "black", size = 3) +
        annotate("text", x = 0.15, y = 0.78, label = "TF", hjust = 0, size = 3) +
        annotate("point", x = 0.08, y = 0.62, shape = 16, colour = "black", size = 2) +
        annotate("text", x = 0.15, y = 0.62, label = "non-TF", hjust = 0, size = 3) +
        annotate("point", x = 0.08, y = 0.46, shape = 16, colour = "red", size = 2) +
        annotate("text", x = 0.15, y = 0.46, label = "CHD (red points)", colour = "red", hjust = 0, size = 3) +
        annotate("text", x = 0.02, y = 0.26, label = legend_text, hjust = 0, size = 3) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        theme_void()

    if (saveFigure) {
        pdf(file = paste0("CP_rank_gene_by_", key, ".pdf"), width = 7, height = 7)
        print((p[[signatures[1]]] | p[[signatures[2]]] | p[[signatures[3]]]) /
            (p_legend | p[[signatures[4]]] | p[[signatures[5]]]))
        dev.off()
    }
    return(df_TF_top)
}


#############################################################################################
## subset a igraph object to a subset that directly connected to a given vetrex of interest ###########

#' Subgraph around seed vertices
#'
#' Induce the subgraph of vertices within `order` steps of `seed`.
#' Returns `NULL` if none of the seeds are present.
#'
#' @param g An igraph object.
#' @param seed Character vector of vertex names.
#' @param order Neighborhood distance (1 = direct neighbors).
#' @param mode Neighborhood mode for directed graphs ("all", "out", "in").
#'
#' @return An igraph subgraph or `NULL`.
#' @export

direct_connected_to_seed <- function(g, seed, order = 1, mode = "all") {
    seeds <- intersect(seed, V(g)$name)
    if (length(seeds) > 0) {
        v_1hop <- unique(unlist(neighborhood(g, order = order, nodes = seeds, mode = mode)))
        g_sub <- induced_subgraph(g, v_1hop)
    } else {
        g_sub <- NULL
    }
    g_sub
}

lign_pair <- function(g1, g2, nodes = NULL) {
    library(igraph)

    if (is.null(nodes)) nodes <- sort(intersect(V(g1)$name, V(g2)$name))

    g1a <- induced_subgraph(g1, nodes)
    g2a <- induced_subgraph(g2, nodes)

    # reorder vertices to EXACTLY the same order in both graphs
    g1a <- permute(g1a, match(nodes, V(g1a)$name))
    g2a <- permute(g2a, match(nodes, V(g2a)$name))

    stopifnot(identical(V(g1a)$name, V(g2a)$name))
    list(g1 = g1a, g2 = g2a, nodes = nodes)
}


#############################################################################################
## identify the TF-targeted pull candidate genes #############################################################################################
# a function used by identify_TF_targeted_pull_candidate to add missing vertices to the graph
add_vertex_if_missing <- function(g, key) {
    vname <- unique(key)
    missing <- setdiff(vname, V(g)$name)
    if (length(missing) > 0) {
        g <- add_vertices(g, nv = length(missing), name = missing)
        # optional: set default attributes for new nodes
        if ("weight" %in% vertex_attr_names(g)) V(g)$weight[match(missing, V(g)$name)] <- NA_real_
        if ("FDR" %in% vertex_attr_names(g)) V(g)$FDR[match(missing, V(g)$name)] <- NA_real_
        if ("color" %in% vertex_attr_names(g)) V(g)$color[match(missing, V(g)$name)] <- "grey80"
    }
    g
}


#' Build TF-targeted CTS/HiG subnetworks
#'
#' Builds four lineage-biased subnetworks (CP, HiGCM, HiGCF, CPopen) and adds
#' dashed TF→target edges from `key` to TF-bound targets present in each subgraph.
#' Newly added edges are marked with `E(g)$lty == "dashed"`.
#'
#' @param mat Gene-by-annotation table (rownames are genes).
#' @param graph_list Named list of igraph objects containing `CTS_ID` and `paste0("HiG", CTS_ID)`.
#' @param CTS_ID Name of the CTS network in `graph_list`.
#' @param CHD CHD genes to highlight (colored red).
#' @param key TF node name to add/connect.
#' @param keep_selfloop Whether to keep `key`→`key`.
#' @param TF_bound_column_name Column in `mat` (0/1) marking TF-bound targets.
#' @param TF_appendix If not NULL, use `{key}_CP_candidate/{key}_CM_candidate/{key}_CF_candidate`.
#' @param edge_colored_by_Maven2023_ISL1KO If TRUE, color new edges using `Maven2023_gene_ISL1_(up|dn)*` columns.
#' @param key_in_TFfamily default is the same as {key}, eitherwise a string  of the predicted TF in a TF motif family
#'
#' @return Named list of 4 igraph subnetworks with updated `color` and `lty` attributes.
#' @importFrom igraph V E induced_subgraph is_directed as.undirected get_edge_ids
#'   add_edges ecount
#'
#' @examples
#' \dontrun{
#' library(igraph)
#'
#' # Minimal toy PPIN (note lower-case to demonstrate internal toupper())
#' g <- make_ring(5)
#' V(g)$name <- c("isl1", "gata4", "gata6", "meis1", "fgf10")
#'
#' graph_list <- list(
#'     HiGCTS_CP.1 = g,
#'     CTS_CP.1    = g
#' )
#'
#' # Minimal annotation matrix
#' mat <- data.frame(
#'     ISL1_CP_bound = c(1, 1, 0, 1, 1),
#'     CP_candidate = c(1, 1, 0, 0, 0),
#'     CM_candidate = c(0, 1, 0, 1, 0),
#'     CF_candidate = c(0, 0, 1, 0, 1),
#'     PCW6CP_access = c(1, 1, 1, 1, 1),
#'     PCW8_CM_access = c(0, 1, 0, 1, 0),
#'     PCW19_CM_access = c(0, 0, 0, 0, 0),
#'     PCW8_CF_access = c(0, 0, 1, 0, 1),
#'     PCW19_CF_access = c(0, 0, 0, 0, 0),
#'     PCW8_SMC_access = c(0, 0, 0, 0, 0),
#'     PCW19_SMC_access = c(0, 0, 0, 0, 0),
#'     Maven2023_gene_ISL1_up_E = c(0, 1, 0, 0, 1),
#'     Maven2023_gene_ISL1_dn_E = c(0, 0, 0, 1, 0)
#' )
#' rownames(mat) <- toupper(V(g)$name)
#'
#' keep_selfloop <- TRUE
#' out <- identify_TF_targeted_pull_candidate(
#'     mat = mat,
#'     graph_list = graph_list,
#'     CTS_ID = "CTS_CP.1",
#'     CHD = c("GATA4"),
#'     key = "ISL1"
#' )
#'
#' names(out)
#' plot(out[[1]])
#' }
#'
#' @export
#'
identify_TF_targeted_pull_candidate <- function(mat, graph_list, CTS_ID, CHD, key = "ISL1",
                                                keep_selfloop = TRUE, # whether to keep the self-loop of the key
                                                TF_bound_column_name = "ISL1_CP_bound", TF_appendix = NULL,
                                                edge_colored_by_Maven2023_ISL1KO = TRUE, key_in_TFfamily = key) {
    four_names <- c("CTSHiG.CP_TF.target", "CTS.CP_TF.target_HiGCM", "CTS.CP_TF.target_HiGCF", "CTSHiG.CP_TF.target_CPopen")

    graph_TF_list <- graph_list[c(paste0("HiG", CTS_ID), CTS_ID, CTS_ID, CTS_ID)] # note that the follwoing three PPINs be replaced later in this sections !!!
    names(graph_TF_list) <- four_names

    ## add the ISL1-targeting as dashlines
    TF_target <- rownames(mat)[which(mat[, TF_bound_column_name] == 1)]
    length(TF_target) # 35

    ## get the subset of CTS that are Isl1-bound in CP
    if (is.null(TF_appendix)) {
        GRN_node_list <- list(
            rownames(mat)[which(mat[, "CP_candidate"] == 1)],
            rownames(mat)[which(mat[, "CM_candidate"] == 1)],
            rownames(mat)[which(mat[, "CF_candidate"] == 1)],
            rownames(mat)[which(mat[, "CP_candidate"] == 1)]
        )
    } else {
        GRN_node_list <- list(
            "CTSHiG.CP_TF.target" = rownames(mat)[which(mat[, paste0(key, "_CP_candidate")] == 1)],
            "CTS.CP_TF.target_HiGCM" = rownames(mat)[which(mat[, paste0(key, "_CM_candidate")] == 1)],
            "CTS.CP_TF.target_HiGCF" = rownames(mat)[which(mat[, paste0(key, "_CF_candidate")] == 1)],
            "CTSHiG.CP_TF.target_CPopen" = rownames(mat)[which(mat[, paste0(key, "_CP_candidate")] == 1)]
        )
    }
    names(GRN_node_list) <- four_names
    lengths(GRN_node_list)
    #     CTSHiG.CP_Isl1.CP.bound     CTS.CP_Isl1.CP.bound_HiGCM
    # 12                             11
    # CTS.CP_Isl1.CP.bound_HiGCF CTSHiG.CP_Isl1.CP.bound_CPopen
    # 6                             12


    ## further narrow the subset to be open at CP, but the bone PPIN of CTS at CP (the 1st object in GRN_node_list) does NOT require CP accessibility,
    # considering Isl1's pioneering role in CP cell fate determination
    for (i in 2:length(GRN_node_list)) {
        # extend the lineage-biased subset by adding back the lieage-psecifically open genes
        if (grepl("CM", names(GRN_node_list)[i])) {
            GRN_node_list[[i]] <- intersect(GRN_node_list[[i]], rownames(mat)[which(mat[, "PCW6CP_access"] == 1 |
                mat[, "PCW8_CM_access"] == 1 | mat[, "PCW19_CM_access"] == 1)])
        } else if (grepl("CF", names(GRN_node_list)[i]) | grepl("SMC", names(GRN_node_list)[i])) {
            GRN_node_list[[i]] <- intersect(GRN_node_list[[i]], rownames(mat)[which(mat[, "PCW6CP_access"] == 1 |
                mat[, "PCW8_CF_access"] == 1 | mat[, "PCW19_CF_access"] == 1 |
                mat[, "PCW8_SMC_access"] == 1 | mat[, "PCW19_SMC_access"] == 1 |
                mat[, "iEPC_access"] == 1)])
        } else {
            GRN_node_list[[i]] <- intersect(GRN_node_list[[i]], rownames(mat)[which(mat[, "PCW6CP_access"] == 1)])
        }
    } # end for(i in 2:length(GRN_node_list)){

    lengths(GRN_node_list)
    #        CTSHiG.CP_Isl1.CP.bound     CTS.CP_Isl1.CP.bound_HiGCM
    # 12                              5
    # CTS.CP_Isl1.CP.bound_HiGCF CTSHiG.CP_Isl1.CP.bound_CPopen
    # 1                              1

    mat <- as.data.frame(apply(mat, 2, as.numeric), row.names = rownames(mat))

    for (name in names(graph_TF_list)) {
        graph <- graph_TF_list[[name]]
        V(graph)$name <- toupper(V(graph)$name) ## for mosue dataset to reuse the mat annotaiton

        graph <- induced_subgraph(graph, intersect(GRN_node_list[[name]], V(graph)$name))
        
        ## ensure node of CHD is red
        V(graph)$color <- ifelse(V(graph)$name %in% CHD, "red", "lightgrey")
        E(graph)$color <- "grey70"
        E(graph)$lty <- "solid"
        ## add the edges between the ISL1-targeting and the (new and old) nodes in the graph

        TF_target_in_graph <- intersect(TF_target, V(graph)$name)
        if (!keep_selfloop) TF_target_in_graph <- setdiff(TF_target_in_graph, key)

        if (length(TF_target_in_graph) > 0) {
            ## check adjacency against an undirected view of the current PPIN
            ## (so an existing edge counts regardless of direction)
            g_check <- if (is_directed(graph)) as.undirected(graph, mode = "collapse") else graph

            flag <- !(key_in_TFfamily %in% V(graph)$name)
            g_check <- add_vertex_if_missing(g_check, key_in_TFfamily)


            if (is.null(vertex_attr(g_check, "color"))) V(g_check)$color <- "lightgrey"
            if (flag) V(g_check)[name == key_in_TFfamily]$color <- "yellow"

            ## get.edge.ids returns 0 for pairs with no edge
            vp_pairs <- as.vector(rbind(key_in_TFfamily, TF_target_in_graph))
            eid <- get_edge_ids(g_check, vp = vp_pairs)

            ## only add TF->target links for missing pairs
            targets_to_add <- TF_target_in_graph[eid == 0]

            if (length(targets_to_add) > 0) {
                m_before <- ecount(g_check)

                edges_vec <- as.vector(rbind(key_in_TFfamily, targets_to_add))
                g_check <- add_edges(g_check, edges_vec)

                new_edges <- seq.int(m_before + 1, ecount(g_check))
                # ---- ensure edge attributes exist for all edges ----
                if (is.null(edge_attr(g_check, "color"))) E(g_check)$color <- rep(NA_character_, ecount(g_check))
                if (is.null(edge_attr(g_check, "lty"))) E(g_check)$lty <- rep(NA_character_, ecount(g_check))

                ## style ONLY those newly added TF->target edges (order matches targets_to_add)
                new_edges_color <- rep("grey70", length(targets_to_add))

                if (edge_colored_by_Maven2023_ISL1KO) { ## style ONLY those newly added TF->target edges (order matches targets_to_add)
                    up_cols <- grep("^Maven2023_gene_ISL1_up", colnames(mat), value = TRUE)
                    if (length(up_cols) > 0) {
                        new_edges_color[rowSums(mat[targets_to_add, up_cols, drop = FALSE]) > 0] <- "blue"
                    }
                    dn_cols <- grep("^Maven2023_gene_ISL1_dn", colnames(mat), value = TRUE)
                    if (length(dn_cols) > 0) {
                        new_edges_color[rowSums(mat[targets_to_add, dn_cols, drop = FALSE]) > 0] <- "orange"
                    }
                } # end  if(edge_colored_by_Maven2023_ISL1KO)

                E(g_check)$color[new_edges] <- new_edges_color
                E(g_check)$lty[new_edges] <- "dashed"
            } # end if (length(targets_to_add) > 0)
            graph <- g_check
        } # end if (length(targets_to_add) > 0)

        graph_TF_list[[name]] <- graph
    } # end for (name in names(graph_TF_list))
    return(graph_TF_list)
}


##########################################################
## plot venn diagram of CTS and its predicted dual pull candidates, which are
## and co-expressed with the key TF (as CTS identification) in transition state, lineage-based express highly, targets of the key TF, and
## open either at CP pool (PCW6) or at its specificed lineages at PCW8 or PCW19 (as dual pull candidates)
# also plotthe graphs of predicted subnetwork from CTS PPIN
# graph_TF_list: a list of igraph object, the output of identify_TF_targeted_pull_candidate()
# TF_key: a string of the key TF symbol
# CTS_ID: a string of the CTS ID
# saveFigure: a logical value to indicate whether to save the figure

plot_TF_targeted_pull_candidate <- function(graph_TF_list, TF_key, CTS_ID, saveFigure = TRUE, saveFileName = NULL) {
    require(gplots)
    # key = gsub(';','\\.',key)
    # graph_TF_list = readRDS(file=paste0('PPI_graph_',TF_key,'_GRN_prediction_',CTS_ID,'_v3.rds'))
    node_color <- c("black", "blue", "darkgreen")
    names(node_color) <- names(graph_TF_list)[1:3]

    # pdf(file='PPI_graph_PRRX1_GRN_prediction.pdf')
    tmp <- list(V(graph_TF_list[[1]])$name, V(graph_TF_list[[2]])$name, V(graph_TF_list[[3]])$name)
    names(tmp) <- names(graph_TF_list)[1:3]
    names(tmp) <- gsub("_TF", paste0("_", TF_key), names(tmp))
    x <- venn(tmp)
    ints <- attr(x, "intersections")
    ints <- ints[-which(names(ints) == paste0("CTSHiG.CP_", TF_key, ".target"))] # drop the first element, regardless of name

    venn_text <- vapply(
        names(ints),
        function(nm) {
            sprintf(
                "%s (%d): %s",
                nm,
                length(ints[[nm]]),
                paste(ints[[nm]], collapse = ", ")
            )
        },
        character(1)
    )
    if (is.null(saveFileName)) saveFileName <- paste0("PPI_graph_", TF_key, "_GRN_prediction_", CTS_ID, "_v3.pdf")
    if (saveFigure) pdf(file = saveFileName) # !!!!!!!
    plot(x)
    # place text above the Venn
    text(
        x = par("usr")[1],
        y = par("usr")[4],
        labels = paste(venn_text, collapse = "\n"),
        adj = c(0, 1),
        cex = 0.6
    )


    for (name in names(graph_TF_list)) {
        graph <- graph_TF_list[[name]]
        set.seed(123) # for reproducibility
        center <- which(V(graph)$name == TF_key)
        lay <- layout_as_star(graph, center = center)

        gene_color <- rep(node_color[name], length(V(graph)$name))
        ## those in CTS_CP & in HiG_CM but not in HiG_CF  is blue !!!!!
        gene_color[V(graph)$name %in% setdiff(V(graph_TF_list[["CTS.CP_TF.target_HiGCM"]])$name, V(graph_TF_list[["CTS.CP_TF.target_HiGCF"]])$name)] <- "blue"
        gene_color[V(graph)$name %in% setdiff(V(graph_TF_list[["CTS.CP_TF.target_HiGCF"]])$name, V(graph_TF_list[["CTS.CP_TF.target_HiGCM"]])$name)] <- "darkgreen"
        if (name != "CTSHiG.CP_TF.target") {
            gene_color[V(graph)$name %in% V(graph_TF_list[["CTSHiG.CP_TF.target_CPopen"]])$name] <- "black"
        }
        plot(graph,
            layout = lay,
            vertex.size = 5,
            vertex.label.color = gene_color,
            main = name
        )
    }
    if (saveFigure) dev.off()
}


##### calculate within cluster gene-gene co-expression for graph_pair, to be used in prioritize_edge_change() ##############################################################

## Helper function:
fill_coexp_for_edges <- function(g, sce, celltype_col = "cluster", cluster_id,
                                 assayName = "logcounts",
                                 ppi_weight = 1,
                                 min_cells = 10,
                                 shrink = TRUE) {
    if (igraph::ecount(g) == 0) {
        return(g)
    }

    # cells in this cluster
    cl <- as.character(SummarizedExperiment::colData(sce)[, celltype_col])
    cells <- which(cl == cluster_id)
    if (length(cells) < min_cells) {
        warning("Too few cells in cluster ", cluster_id, " (", length(cells), "); leaving edges unchanged.")
        # return(g)
        return(NULL)
    }

    # Helper: get newly added edges
    new_edge_ids <- function(g, dashed_value = "dashed") {
        if (!"lty" %in% igraph::edge_attr_names(g)) {
            return(integer(0))
        }
        which(igraph::E(g)$lty == dashed_value)
    }

    g_id_toupdate <- new_edge_ids(g)
    if (length(g_id_toupdate) > 0) {
        EE <- igraph::ends(g, es = igraph::E(g)[g_id_toupdate], names = TRUE)
        genes_needed <- unique(as.vector(EE))
        genes_in_sce <- intersect(genes_needed, rownames(sce))

        g_cor_mat <- compute_cluster_correlation(sce,
            cluster_id = cluster_id,
            genes = genes_in_sce,
            celltype_col = celltype_col,
            assayName = assayName,
            min_cells = min_cells,
            shrink = shrink
        )
        # fetch correlation for each edge
        cor_vec <- mapply(function(a, b) g_cor_mat[a, b], EE[, 1], EE[, 2])
        E(g)$coexp_target[g_id_toupdate] <- abs(cor_vec)
        E(g)$corexp_sign[g_id_toupdate] <- ifelse(cor_vec >= 0, "positive", "negative")
        E(g)$weight[g_id_toupdate] <- abs(cor_vec) * ppi_weight
        E(g)$norm_PPI_score[g_id_toupdate] <- ppi_weight
    }
    return(g)
}


#############################################################################################
## Edge weight change table between two igraph networks #############################################################################################
#' Edge weight change table between two graphs
#'
#' Compares two igraph networks by merging edges (by endpoint names) and computing
#' `delta = w2 - w1` plus change labels (gained/lost/changed/unchanged), if PPI(i,j)>0
#' otehrwise delta = r2 - r1 where r() is the shriked coexprssion value of two genes in a cell clsuters
#'
#' @param g1,g2 igraph objects.
#' @param weight_attr Edge attribute to compare (default: "weight").
#' @param missing_as Weight to use when an edge is absent (default: 0).
#' @param undirected If TRUE, treat A--B the same as B--A.
#'
#' @return A data.frame with `from`, `to`, `w1`, `w2`, `delta`, `abs_delta`, `direction`, `status`, `rank`.
#' @export
edge_change_table <- function(g1, g2, weight_attr = "weight", missing_as = 0, undirected = TRUE) {
    # require vertex names
    if (is.null(V(g1)$name) || is.null(V(g2)$name)) stop("Both graphs must have V(g)$name")

    # require corexp_sign if you want signed weights
    if (is.null(edge_attr(g1, "corexp_sign")) || is.null(edge_attr(g2, "corexp_sign"))) {
        stop("'corexp_sign' is not an edge attribute of one or both graphs")
    }

    # helper to extract a clean edge table with signed weights
    extract_edges <- function(g, w_attr, undirected) {
        e <- igraph::as_data_frame(g, what = "edges")

        e$from <- as.character(e$from)
        e$to <- as.character(e$to)

        if (undirected) {
            a <- pmin(e$from, e$to)
            b <- pmax(e$from, e$to)
            e$from <- a
            e$to <- b
        }

        # raw weight (fallback to 1 if missing attribute)
        w <- if (w_attr %in% names(e)) e[[w_attr]] else rep(1, nrow(e))
        # fix NA weights (e.g. newly added dashed TF->target edges)
        w[is.na(w)] <- 0

        # signed weight using corexp_sign
        sign_vec <- e$corexp_sign
        sign_vec[is.na(sign_vec)] <- "positive"
        signed_w <- ifelse(sign_vec == "positive", w, -w)

        data.frame(from = e$from, to = e$to, w = signed_w, stringsAsFactors = FALSE)
    }

    e1 <- extract_edges(g1, weight_attr, undirected)
    e2 <- extract_edges(g2, weight_attr, undirected)

    # collapse duplicates (keep max abs weight; adjust if you prefer sum/mean)
    collapse_fun <- function(df) {
        key <- paste(df$from, df$to, sep = "||")
        # pick the entry with largest |w|
        idx <- tapply(seq_len(nrow(df)), key, function(ii) ii[which.max(abs(df$w[ii]))])
        df[unlist(idx), , drop = FALSE]
    }
    e1 <- collapse_fun(e1)
    names(e1)[3] <- "w1"
    e2 <- collapse_fun(e2)
    names(e2)[3] <- "w2"

    m <- merge(e1, e2, by = c("from", "to"), all = TRUE)

    m$w1[is.na(m$w1)] <- missing_as
    m$w2[is.na(m$w2)] <- missing_as

    m$delta <- m$w2 - m$w1
    m$abs_delta <- abs(m$delta)

    m$direction <- ifelse(m$delta > 1e-10, "increase",
        ifelse(m$delta < -1e-10, "decrease", "unchanged")
    )

    m$status <- ifelse(m$w1 == missing_as & m$w2 != missing_as, "gained",
        ifelse(m$w1 != missing_as & m$w2 == missing_as, "lost",
            ifelse(m$w1 != m$w2, "changed", "unchanged")
        )
    )

    m <- m[order(m$abs_delta, decreasing = TRUE), ]
    m$rank <- seq_len(nrow(m))
    rownames(m) <- NULL
    m
}


#############################################################################################
#' Annotate and plot edge-weight changes on a graph
#'
#' Adds per-edge change metrics (e.g., `delta`) from `edge_change_df` to `g1` and
#' optionally plots the network with edge color/width driven by the change.
#'
#' @param g1 An igraph object.
#' @param edge_change_df Data frame with at least `regulator`, `target`, and `delta`.
#' @param top_n Number of edges with largest |delta| to label.
#' @param title Plot title (also used in default PDF name).
#' @param savepdf If TRUE, save a PDF.
#' @param ... Passed to `plot.igraph()`.
#'
#' @return `g1` with added edge attributes prefixed `TIPS_`.
#' @export
#'
prioritize_edge_change <- function(g1, edge_change_df, top_n = 5, title = "CM-pull subnetwork",
                                   savepdf = TRUE) { # , layout = NULL
    # Initialize edge attributes
    E(g1)$TIPS_w1 <- E(g1)$weight
    E(g1)$TIPS_w2 <- E(g1)$weight
    E(g1)$TIPS_delta <- 0
    E(g1)$TIPS_direction <- "unchanged"

    # Match dataframe edges to graph edges (robust way)
    # edge endpoints in g1
    edge_ends <- ends(g1, E(g1), names = TRUE)

    # helper to find edge index
    find_edge <- function(a, b, ends_mat) {
        which(
            (ends_mat[, 1] == a & ends_mat[, 2] == b) |
                (ends_mat[, 1] == b & ends_mat[, 2] == a)
        )
    }
    # Update edges using edge_change_df
    for (i in seq_len(nrow(edge_change_df))) {
        a <- edge_change_df$from[i]
        b <- edge_change_df$to[i]

        eid <- find_edge(a, b, edge_ends)

        if (length(eid) == 1) {
            E(g1)$TIPS_w1[eid] <- edge_change_df$w1[i]
            E(g1)$TIPS_w2[eid] <- edge_change_df$w2[i]
            E(g1)$TIPS_delta[eid] <- edge_change_df$delta[i]
            E(g1)$TIPS_direction[eid] <- edge_change_df$direction[i]
            # 	E(g1)$TIPS_status[eid]    <- edge_change_df$status[i]
        }
    }

    # Edge color by directionE(g1)$color <- "grey80"
    E(g1)$color[E(g1)$TIPS_direction == "increase"] <- "#D73027" # red
    E(g1)$color[E(g1)$TIPS_direction == "decrease"] <- "#4575B4" # blue
    # Edge width by |delta|
    E(g1)$width <- 1 + 4 * abs(E(g1)$TIPS_delta) / max(abs(E(g1)$TIPS_delta), na.rm = TRUE)
    # make unchanged edges thinner:
    E(g1)$width[E(g1)$TIPS_status == "unchanged"] <- 0.5

    # label only top-ranked changes
    top_edges <- edge_change_df$rank <= top_n
    label_edges <- unlist(
        mapply(find_edge,
            edge_change_df$from[top_edges],
            edge_change_df$to[top_edges],
            MoreArgs = list(ends_mat = edge_ends)
        )
    )

    E(g1)$label <- NA
    E(g1)$label[label_edges] <- round(E(g1)$TIPS_delta[label_edges], 3)
    E(g1)$label.cex <- 0.8

    if (grepl("CM", title)) label <- "CMvsCP" else label <- "CFvsCP"
    if (savepdf) pdf(file = paste0("TIPS_delta_edge_reweighting_", title, ".pdf"), width = 7, height = 7)
    # if(is.null(layout)) layout = layout_with_fr(g1)
    # avoide the error that layour reuires positive edge weight, but we may have zero
    plot(
        g1,
        layout = layout_with_fr(g1, weights = NA),
        edge.curved = 0.15,
        vertex.size = 22,
        vertex.label.cex = 0.9,
        main = paste0("TIPS delta-edge reweighting (", title, ")")
    )
    mtext(paste0(label, " Top ", top_n, " edges labeled by delta"), side = 1, line = -1, cex = 2)

    if (savepdf) dev.off()

    g1
}

get_regulators_from_motifs <- function(cisTarget.res, CTS_name, NES_threshold = 3, motifAnnot = NULL,
                                       sep = ";", toupper_out = FALSE) {
    motifs <- subset(cisTarget.res, geneSet == CTS_name & NES >= NES_threshold)
    motifs <- c(motifs$motif, motifs$TF_highConf)

    motifs_in <- motifs
    motifs <- trimws(motifs)

    # ---- helpers ----
    extract_symbols <- function(x) {
        x <- as.character(x)
        x <- x[!is.na(x)]
        x <- trimws(x)
        x <- x[nzchar(x)]
        if (length(x) == 0) {
            return(character(0))
        }

        # remove "(directAnnotation)" etc.
        x <- gsub("\\([^)]*\\)", "", x)

        # "." often separates blocks like "A ... . B ..."
        x <- gsub("\\.+", ";", x)

        # unify separators
        x <- gsub(",", ";", x)

        parts <- unlist(strsplit(paste(x, collapse = ";"), ";", fixed = TRUE), use.names = FALSE)
        parts <- trimws(parts)
        parts <- parts[nzchar(parts)]

        # strip trailing punctuation, remove internal whitespace
        parts <- gsub("[[:punct:]]+$", "", parts)
        parts <- gsub("\\s+", "", parts)

        # keep gene-like tokens only
        parts <- parts[grepl("^[A-Za-z][A-Za-z0-9-]*$", parts)]
        if (toupper_out) parts <- toupper(parts)

        parts
    }

    collapse_symbols <- function(x) {
        x <- unique(x)
        x <- sort(x)
        if (length(x) == 0) {
            return(NA_character_)
        }
        paste(x, collapse = sep)
    }

    clean_to_key <- function(x) collapse_symbols(extract_symbols(x))

    # ---- 1) parse motif strings / IDs (fallback) ----
    parsed <- rep(NA_character_, length(motifs))

    # (a) annotation-label strings like "ETS2 (directAnnotation). " or long multi-TF strings
    i <- grepl("directAnnotation|inferredBy", motifs, ignore.case = TRUE) | grepl("\\(", motifs)
    parsed[i] <- vapply(motifs[i], clean_to_key, character(1))

    # (b) plain gene symbol (e.g., "SPIB")
    i <- is.na(parsed) & grepl("^[A-Za-z][A-Za-z0-9-]*$", motifs)
    parsed[i] <- if (toupper_out) toupper(motifs[i]) else motifs[i]

    # (c) motif ID patterns that embed TF names
    i <- is.na(parsed) & grepl("^hdpi__", motifs)
    parsed[i] <- vapply(sub("^hdpi__([^_]+).*$", "\\1", motifs[i]), clean_to_key, character(1))

    i <- is.na(parsed) & grepl("^hocomoco__", motifs)
    parsed[i] <- vapply(sub("^hocomoco__([^_]+)_.*$", "\\1", motifs[i]), clean_to_key, character(1))

    i <- is.na(parsed) & grepl("^swissregulon__", motifs)
    parsed[i] <- vapply(sub("^swissregulon__[^_]+__([^_]+).*$", "\\1", motifs[i]), clean_to_key, character(1))

    i <- is.na(parsed) & grepl("^taipale_tf_pairs__", motifs)
    parsed[i] <- vapply(sub("^taipale_tf_pairs__([^_]+)_([^_]+)_.*$", "\\1;\\2", motifs[i]), clean_to_key, character(1))

    i <- is.na(parsed) & grepl("^kznf__", motifs)
    parsed[i] <- vapply(sub("^kznf__([^_]+)_.*$", "\\1", motifs[i]), clean_to_key, character(1))

    i <- is.na(parsed) & grepl("^dbtfbs__", motifs)
    parsed[i] <- vapply(sub("^dbtfbs__([^_]+)_.*$", "\\1", motifs[i]), clean_to_key, character(1))

    # ---- 2) motifAnnot mapping (preferred when available) ----
    ann_map <- rep(NA_character_, length(motifs))
    if (!is.null(motifAnnot)) {
        if (!all(c("motif", "TF") %in% colnames(motifAnnot))) {
            stop("motifAnnot must have at least columns: motif, TF")
        }

        mot <- trimws(as.character(motifAnnot$motif))
        tf <- as.character(motifAnnot$TF)
        ok <- !is.na(mot) & nzchar(mot) & !is.na(tf) & nzchar(tf)
        mot <- mot[ok]
        tf <- tf[ok]

        # collapse TFs per motif into a single cleaned key string
        idx_by_motif <- split(seq_along(mot), mot)
        tf_by_motif <- vapply(
            idx_by_motif,
            function(ii) collapse_symbols(extract_symbols(tf[ii])),
            character(1)
        )

        ann_map <- unname(tf_by_motif[motifs])
    }

    # ---- 3) final regulators: annotation mapping when available, otherwise parsed ----
    regulators <- ifelse(!is.na(ann_map) & nzchar(ann_map), ann_map, parsed)

    motifAnnot_sub <- data.frame(
        motif_TF_highConf = motifs_in,
        regulators = regulators,
        parsed_from_id = parsed,
        mapped_from_annot = ann_map,
        stringsAsFactors = FALSE
    )

    # --- 4 check redundent regulator and pick the one with highest NES ------
    reg_dup <- motifAnnot_sub$regulators[duplicated(motifAnnot_sub$regulators)]
    if (length(reg_dup) > 0) {
        motifAnnot_unique <- motifAnnot_sub[!duplicated(motifAnnot_sub$regulators), ]
        for (reg in reg_dup) {
            motif_found <- motifAnnot_sub[which(motifAnnot_sub$regulators %in% reg), ]$motif
            cis_found <- subset(cisTarget.res, geneSet == CTS_name & motif %in% motif_found)
            if (nrow(cis_found) > 0) {
                cis_found <- subset(cis_found, NES == max(cis_found$NES))
                ## in case equal NES, pick the one randomly
                if (nrow(cis_found) > 1) cis_found <- cis_found[sample(1:nrow(cis_found), 1), ]
                x <- which(motifAnnot_unique$regulators == reg)
                motifAnnot_unique[x, ] <- motifAnnot_sub[which(motifAnnot_sub$motif == cis_found$motif), ]
            }
        }
    } else {
        motifAnnot_unique <- motifAnnot_sub
    }
    return(motifAnnot_unique)
}

# get_regulators_from_motifs_v0 <- function(motifs, motifAnnot = NULL) {
#    stopifnot(is.character(motifs))

# #  # --- 1) Regex-based parsing for common motif name styles ---
#   parsed <- rep(NA_character_, length(motifs))

#   # hocomoco__NFKB1_HUMAN.H11MO.0.A  -> NFKB1
#   i <- grepl("^hocomoco__", motifs)
#   parsed[i] <- sub("^hocomoco__([^_]+)_.*$", "\\1", motifs[i])

#   # swissregulon__hs__EZH2  -> EZH2 ; swissregulon__mm__Atf2 -> Atf2
#   i <- grepl("^swissregulon__", motifs)
#   parsed[i] <- sub("^swissregulon__[^_]+__([^_]+).*$", "\\1", motifs[i])

#   # taipale_tf_pairs__ETV2_ONECUT2_... -> ETV2;ONECUT2 (keep both)
#   i <- grepl("^taipale_tf_pairs__", motifs)
#   parsed[i] <- sub("^taipale_tf_pairs__([^_]+)_([^_]+)_.*$", "\\1;\\2", motifs[i])

#   # kznf__ZNF263_Transfac... -> ZNF263
#   i <- grepl("^kznf__", motifs)
#   parsed[i] <- sub("^kznf__([^_]+)_.*$", "\\1", motifs[i])


#   # --- 2) Annotation-based mapping (recommended for transfac_pro__, metacluster__, tfdimers__, etc.) ---
#   ann_map <- rep(NA_character_, length(motifs))
#   if (!is.null(motifAnnot)) {
#     # motifAnnot should contain columns: motif, TF, (often) directAnnotation, etc.
#     if (!all(c("motif", "TF") %in% colnames(motifAnnot))) {
#       stop("motifAnnot must have at least columns: motif, TF")
#     }
#     # For each motif, collect unique TFs
#     tf_by_motif <- tapply(motifAnnot$TF, motifAnnot$motif, function(x) paste(unique(x), collapse = ";"))
#     ann_map <- unname(tf_by_motif[motifs])
#   }

#   # choose annotation mapping when available, otherwise parsed
#   regulators <- ifelse(!is.na(ann_map) & ann_map != "", ann_map, parsed)

#   data.frame(
#     motif = motifs,
#     regulators = regulators,
#     parsed_from_id = parsed,
#     mapped_from_annot = ann_map,
#     stringsAsFactors = FALSE
#   )
# }


heatmap_pull_candidate <- function(
  mat, graph_list, CTS_ID, CHD, key = "ISL1", coding_genes = NULL, TF = NULL,
  chip_targets = FALSE, show_SMC_access = FALSE
) {
    library(pheatmap)
    key <- gsub(";", "\\.", key)

    PPI_gene <- V(graph_list[[CTS_ID]])$name
    In_PPI <- ifelse(rownames(mat) %in% PPI_gene, TRUE, FALSE)
    IS_CHD <- ifelse(rownames(mat) %in% CHD, TRUE, FALSE)
    Is_protein_coding <- ifelse(rownames(mat) %in% coding_genes, TRUE, FALSE)
    # TF_mouse <- unique(dorothea_mm$tf)
    Is_TF <- ifelse(rownames(mat) %in% TF, TRUE, FALSE)

    ## example: build row annotation with 3 yes/no columns
    row_anno <- data.frame(
        is_TF = factor(ifelse(Is_TF, TRUE, FALSE), levels = c(FALSE, TRUE), labels = c("no", "yes")),
        coding = factor(ifelse(Is_protein_coding, TRUE, FALSE), levels = c(FALSE, TRUE), labels = c("no", "yes")),
        CHD = factor(ifelse(IS_CHD, TRUE, FALSE), levels = c(FALSE, TRUE), labels = c("no", "yes")),
        PPI = factor(ifelse(In_PPI, TRUE, FALSE), levels = c(FALSE, TRUE), labels = c("no", "yes"))
    )
    rownames(row_anno) <- rownames(mat)

    ## same color scheme for all three bars: no = lightgrey, yes = blue
    ann_colors <- list(
        CHD = c(no = "lightgrey", yes = "blue"),
        PPI = c(no = "lightgrey", yes = "blue"),
        is_TF = c(no = "lightgrey", yes = "blue")
    )

    ann_colors[["coding"]] <- c(no = "lightgrey", yes = "blue")

    ########## update column names to non-key-specific names ##########
    {
        # replace to non-key-specific column names
        colnames(mat) <- gsub(paste0(key, "_CP_candidate"), "CP_candidate", colnames(mat))
        colnames(mat) <- gsub(paste0(key, "_CM_candidate"), "CM_candidate", colnames(mat))
        colnames(mat) <- gsub(paste0(key, "_CF_candidate"), "CF_candidate", colnames(mat))
    }

    # ### a plot with less rows: -- only PRRX1-binding coding genes are shown
    x <- which(rownames(mat) %in% coding_genes & (mat[, "CP_candidate"] == 1 | mat[, "CM_candidate"] == 1 | mat[, "CF_candidate"] == 1))
    length(x) # 3
    mat_sub <- mat[x, ]
    row_anno_sub <- row_anno[x, ]
    row_anno_sub <- row_anno_sub[, -which(colnames(row_anno_sub) == "coding")]
    cat("candidate genes: ", nrow(mat_sub), "\n")

    pattern <- paste0(mat_sub[, "CP_candidate"], mat_sub[, "CM_candidate"], mat_sub[, "CF_candidate"])
    (y <- table(pattern))
    # 	pattern
    # 101 111
    #   1   2
    ## map pattern → block name
    pattern_to_block <- c(
        "100" = "CP_only",
        "110" = "CP_CM",
        "101" = "CP_CF",
        "111" = "CP_CM_CF",
        "010" = "CM_only",
        "001" = "CF_only",
        "011" = "CM_CF",
        "000" = "none"
    )

    block <- pattern_to_block[pattern]

    ## choose order you want to see in heatmap
    # block_levels <- c("CP_only","CP_CM","CM_only","CP_CM_CF","CM_CF","CF_only","CP_CF","none")
    block_levels <- pattern_to_block[which(names(pattern_to_block) %in% unique(pattern))]
    block <- factor(block, levels = block_levels)
    names(block) <- rownames(mat_sub) # keep rownames aligned

    ord <- order(block)
    cols_to_show <- c(
        "CP_candidate", "CP_hi", "PCW6CP_access", # "PCW8_CM_access"  ,
        "CM_candidate", "CM_hi", "PCW8_CM_access", "PCW19_CM_access",
        "CF_candidate", "CF_hi", "iEPC_access", "PCW8_CF_access", "PCW19_CF_access"
    )
    if (show_SMC_access) {
        cols_to_show <- c(cols_to_show, "PCW8_SMC_access", "PCW19_SMC_access")
    }
    if (chip_targets) {
        cols_to_show <- c(
            cols_to_show, "ISL1_CP_bound", "Maven2023_gene_ISL1_WT_d6CP",
            "Gao2019_gene_Isl1_E825E9.bound", "Gao2019_gene_Isl1.iCPC_CPC.bound",
            "Maven2023_gene_ISL1_up_E",
            "Maven2023_gene_ISL1_up_T", "Maven2023_gene_ISL1_up_L", "Maven2023_gene_ISL1_dn_E",
            "Maven2023_gene_ISL1_dn_T", "Maven2023_gene_ISL1_dn_L"
        )

        title <- paste0("candidate: ", key, " CP ChIP-seq targets that are highly expressed in a state")
    } else {
        # cols_to_show = c(cols_to_show, paste0('cisTarget_',key,'.motif_target'))

        y <- which(grepl(key, colnames(mat_sub)) & grepl("cisTarget_", colnames(mat_sub)))
        muti_match_flag <- grepl(";", colnames(mat_sub)[y])
        if (length(y) > 0 && any(!muti_match_flag)) {
            y <- y[!muti_match_flag]
            muti_match_flag <- FALSE
            cat("direct motif: ", colnames(mat_sub)[y], " used from multiple potentoal matches\n")
        }
        if (any(muti_match_flag)) {
            cat("simplified motif name: ", colnames(mat_sub)[y], "\n")
            colnames(mat_sub)[y] <- paste0("cisTarget_", key, ".motif_target") ## in case motif name is too long, simplify it
        }
        cols_to_show <- c(cols_to_show, colnames(mat_sub)[y])

        title <- paste0("candidate: ", key, " cisTarget targets that are highly expressed in a state")
    }

    row_label_col <- ifelse(mat_sub[, "CP_candidate"] == 1, "red", "black")

    # Ensure all values in mat are numeric before plotting
    # mat_sub_numeric <- as.data.frame(apply(mat_sub, 2, as.numeric), row.names = rownames(mat_sub)) %>% as.matrix
    mat_sub_numeric <- data.matrix(mat_sub)
    cols_to_show <- intersect(cols_to_show, colnames(mat_sub_numeric))
    if (length(cols_to_show) == 0) stop("No columns found in mat_sub for plotting.")
    p <- pheatmap(
        mat_sub_numeric[ord, cols_to_show, drop = FALSE],
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        color = c("white", "orange"),
        annotation_row = row_anno_sub,
        annotation_colors = ann_colors,
        margins = c(4, 4),
        fontsize = 10,
        fontsize_col = 6,
        fontsize_row = 6, # reduce to fit more if necessary, or adjust as you want
        main = title
    )
    return(p)
}


#############################################################################################

# graph_TF_list: the dual-pull predicted subnetwork of PPIN (CTS at bifurcation)
# graph_list: the edge-precalcualted PPIN list
# sce: the "SingleCellExperiment" object of the datasets
# celltype_col: a string to specify the columns of cell cluster
# descendant_cluster_id: a string of the cluster ID of the descendant of interest
fill_TF_targeting_predicted_edges <- function(graph_TF_list, linkage_name = "CM", graph_list,
                                              sce, celltype_col = "cluster", CT_cluster_id = CP_cluster,
                                              descendant_cluster_id = CM_cluster, TF_symbol = "ISL1", HVG = NULL,
                                              shrink = TRUE) {
    if (!paste0("CTS.CP_TF.target_HiG", linkage_name) %in% names(graph_TF_list)) stop(paste0("the parent PPIN 'CTS.CP_TF.target_HiG", linkage_name, " is missing from graph_TF_list"))

    CM_ID <- paste0("HiG_", descendant_cluster_id)

    # helper function to find the missing edges
    edge_keys <- function(g) {
        el <- igraph::as_edgelist(g, names = TRUE)
        a <- pmin(el[, 1], el[, 2])
        b <- pmax(el[, 1], el[, 2])
        paste(a, b, sep = "--")
    }

    # for CM-pull
    g1 <- graph_TF_list[[paste0("CTS.CP_TF.target_HiG", linkage_name)]]
    g_descendant_sub <- graph_list[[CM_ID]]

    ## check if all nodes are in the HVG
    if (!is.null(HVG)) {
        if (!all(V(g1)$name %in% HVG)) stop("all nodes are not in the HVG")
        if (!all(V(g_descendant_sub)$name %in% HVG)) stop("all nodes are not in the HVG")
    }
    ## add back the coexpression between keyTF and targets that was missing from PPI database but suggested by ChIP-seq or cisTarget
    if (vcount(g1) > 0) g1 <- fill_coexp_for_edges(g1, sce, celltype_col = celltype_col, cluster_id = CT_cluster_id)

    if (vcount(g1) > 0) {
        ## the subset of CTS genes that also highly expressed in CM thus in the edge-calculated PPIN of HiG
        g_descendant_sub <- induced_subgraph(g_descendant_sub, vids = V(g_descendant_sub)$name %in% V(g1)$name)
        # add key TF back to the subset of g_descendant_sub if it is within HiG but not in PPINs thus was missing
        if (!(TF_symbol %in% V(g_descendant_sub)$name) & (TF_symbol %in% DEG[[descendant_cluster_id]])) {
            g_descendant_sub <- igraph::add_vertices(g_descendant_sub, 1, name = TF_symbol)
            V(graph)[name == key_in_TFfamily]$color <- "yellow"
            g_descendant_sub <- fill_coexp_for_edges(g_descendant_sub, sce, celltype_col = celltype_col, cluster_id = descendant_cluster_id)
        }
    }

    cat("---- DEBUG: Edge presence ----\n")
    cat("In g1:", any(edge_keys(g1) == "DAB2--HMGA2"), "\n")
    cat("In g_descendant_sub:", any(edge_keys(g_descendant_sub) == "DAB2--HMGA2"), "\n")

    ## fill back coexpression in the descendant cluster for missing linkages
    if (vcount(g1) > 0 & vcount(g_descendant_sub) > 0) {
        # find the missing edges
        missing_in_desc <- setdiff(edge_keys(g1), edge_keys(g_descendant_sub)) # edges present in g1 but absent in g_descendant_sub

        if (length(missing_in_desc) > 0) {
            # add back missing_in_desc to g_descendant_sub, assign lty = dashed
            pairs <- do.call(rbind, strsplit(missing_in_desc, "--", fixed = TRUE))
            from <- pairs[, 1]
            to <- pairs[, 2]

            # make sure vertices exist (should, but safe)
            missing_v <- setdiff(unique(c(from, to)), V(g_descendant_sub)$name)
            if (length(missing_v) > 0) {
                g_descendant_sub <- add_vertices(g_descendant_sub, nv = length(missing_v), name = missing_v)
            }

            # build edge vector: from1,to1, from2,to2,...
            edges_vec <- as.vector(rbind(from, to))

            m_before <- ecount(g_descendant_sub)
            g_descendant_sub <- add_edges(g_descendant_sub, edges_vec)
            new_edges <- seq.int(m_before + 1, ecount(g_descendant_sub))

            # ensure attrs exist, then style new edges
            if (!"lty" %in% edge_attr_names(g_descendant_sub)) E(g_descendant_sub)$lty <- "solid"
            E(g_descendant_sub)$lty[new_edges] <- "dashed"

            cat("---- DEBUG: After edge insertion ----\n")
            print(
                igraph::as_data_frame(g_descendant_sub, what = "edges")[
                    (from == "DAB2" & to == "HMGA2") |
                        (from == "HMGA2" & to == "DAB2"),
                ]
            )

            g_descendant_sub <- fill_coexp_for_edges(g_descendant_sub, sce, celltype_col = celltype_col, cluster_id = descendant_cluster_id, shrink = shrink)


            cat("---- DEBUG: After coexpression fill ----\n")

            edges_df <- igraph::as_data_frame(g_descendant_sub, what = "edges")

            subset_edge <- edges_df[
                (edges_df$from == "DAB2" & edges_df$to == "HMGA2") |
                    (edges_df$from == "HMGA2" & edges_df$to == "DAB2"),
            ]

            print(subset_edge)
        }
    } else {
        g_descendant_sub <- NULL
    }

    return(list(g_CT_sub = g1, g_descendant_sub = g_descendant_sub))
}

########################################################
#' Merge cisTarget/ChIP-derived edge-delta tables into a single styled igraph network
#'
#' This function merges many TF-specific “Top edges labeled by delta” results (stored in a
#' single table) into one undirected graph, de-duplicates repeated edges by keeping the
#' row with the largest absolute change (\code{abs_delta}), and applies a consistent
#' visualization style:
#' \itemize{
#'   \item Node color: CHD genes in red; “added TF” nodes in yellow; all others in grey.
#'   \item Edge color: red for \code{direction == "increase"}, blue for \code{"decrease"},
#'         grey otherwise.
#'   \item Edge width: proportional to \code{abs_delta}.
#'   \item Edge line type: dashed if an edge is not present in a reference STRING PPIN
#'         graph (\code{g_string}); otherwise solid.
#'   \item Edge labels: \code{top_n_label} edges (by \code{abs_delta}) labeled with \code{delta}.
#' }
#'
#' The function is designed for outputs from your TIPS delta-edge reweighting pipeline
#' where \code{final_table} contains columns such as \code{from}, \code{to}, \code{delta},
#' \code{abs_delta}, \code{direction}, \code{status}, and \code{rank}.
#'
#' @param final_table A data.frame of edge changes. Must contain at least columns:
#'   \code{linkage}, \code{from}, \code{to}, \code{delta}, \code{abs_delta},
#'   \code{direction}, \code{status}, \code{rank}. Additional columns are allowed.
#'
#' @param descendant Character scalar. Which lineage to keep (matched against the
#'   \code{linkage} column), e.g. \code{"CM"} or \code{"CF"}. Default \code{"CM"}.
#'
#' @param CHD Character vector of CHD genes to highlight as red nodes. Default \code{character(0)}.
#'
#' @param added_TF Character vector of TF symbols that were added as nodes because they
#'   were missing from the STRING PPIN. These are colored yellow unless they overlap with
#'   \code{CHD} (CHD overrides). Default \code{character(0)}.
#'
#' @param top_n_label Integer. Number of edges to label (largest \code{abs_delta}).
#'   Set \code{0} to disable labeling. Default \code{5}.
#'
#' @param g_string Optional igraph object representing the reference STRING PPIN
#'   subnetwork (same naming convention/casing as \code{final_table}). If provided,
#'   dashed edges indicate edges absent from \code{g_string}. If NULL (default), line
#'   type defaults to solid for all edges (or you can set dashed by status upstream).
#'
#' @param normalize_case Logical. If TRUE (default), vertex names in both graphs are
#'   uppercased prior to edge comparison against \code{g_string}. Set FALSE if your
#'   naming is already consistent and case-sensitive.
#'
#' @return An undirected igraph object with vertex attributes:
#'   \code{name}, \code{color}, and edge attributes including:
#'   \code{delta}, \code{abs_delta}, \code{direction}, \code{status}, \code{rank},
#'   \code{weight} (=\code{abs_delta}), \code{color}, \code{width}, \code{lty},
#'   and \code{label} (for the top edges).
#'
#' @details
#' ## De-duplication
#' Because multiple TF panels may contain the same edge (or the same TF may contribute
#' repeated rows), the function creates an undirected edge key by sorting endpoints
#' (\code{pmin(from,to)}, \code{pmax(from,to)}) and keeps only the row with maximum
#' \code{abs_delta} for each unique pair.
#'
#' ## Line type semantics (recommended)
#' If you want dashed lines to strictly mean “not in STRING PPIN”, pass the appropriate
#' STRING subgraph as \code{g_string}. This is more robust than using \code{status}
#' because \code{status} can reflect gained/lost under the delta comparison rather than
#' “absent from STRING”.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' library(dplyr)
#'
#' # final_table: output table read from file
#' final_table <- read.table("final_table.tsv",
#'     header = TRUE, sep = "\t",
#'     stringsAsFactors = FALSE
#' )
#'
#' # Optional reference STRING PPIN subgraph for the same lineage
#' # g_string <- graph_list[["HiG_CM_cluster"]]  # example
#'
#' g <- make_merged_TIPS_graph(
#'     final_table,
#'     descendant = "CM",
#'     CHD = c("GATA6", "TBX5"),
#'     added_TF = c("ETV2", "PRDM6"),
#'     top_n_label = 10,
#'     g_string = NULL
#' )
#'
#' set.seed(1)
#' plot(g,
#'     layout = layout_with_fr(g, weights = NA),
#'     edge.curved = 0.15, vertex.size = 22, vertex.label.cex = 0.9
#' )
#'
#' # With STRING reference for dashed edges:
#' # g <- make_merged_TIPS_graph(final_table, "CM", CHD, added_TF, 10, g_string = g_string)
#' }
#'
#' @export

make_merged_TIPS_graph <- function(final_table,
                                   descendant = "CM",
                                   CHD = character(0),
                                   added_TF = character(0), # TFs that were added because missing in STRING PPIN
                                   top_n_label = 5,
                                   g_string = NULL,
                                   normalize_case = TRUE) {
    library(dplyr)
    req <- c("linkage", "from", "to", "delta", "abs_delta", "direction", "status", "rank")
    stopifnot(all(req %in% colnames(final_table)))

    ## two helper functions
    library(igraph)

    edge_keys <- function(g, undirected = TRUE) {
        el <- igraph::as_edgelist(g, names = TRUE)
        if (undirected && !igraph::is_directed(g)) {
            a <- pmin(el[, 1], el[, 2])
            b <- pmax(el[, 1], el[, 2])
            paste(a, b, sep = "||")
        } else {
            paste(el[, 1], el[, 2], sep = "||")
        }
    }

    set_lty_by_reference <- function(g, g_string, undirected = TRUE) {
        keys_string <- edge_keys(g_string, undirected = undirected)
        keys_g <- edge_keys(g, undirected = undirected)

        E(g)$lty <- "solid"
        E(g)$lty[!(keys_g %in% keys_string)] <- "dashed"
        g
    }

    # Build graph with all edge attributes carried in
    g <- graph_from_data_frame(
        final_table %>% dplyr::select(from, to, delta, abs_delta, direction, status, rank),
        directed = FALSE
    )
    E(g)$weight <- E(g)$abs_delta
    g <- igraph::simplify(g, edge.attr.comb = "max")

    # ---- Node colors: CHD red; added TF yellow; else grey ----
    V(g)$color <- "grey80"
    if (length(CHD) > 0) V(g)$color[V(g)$name %in% CHD] <- "red"
    if (length(added_TF) > 0) V(g)$color[V(g)$name %in% added_TF] <- "yellow" # overrides red

    # ---- Edge style ----
    # color by 'direction' of changes from CP to descedent
    E(g)$color <- "grey80"
    E(g)$color[E(g)$direction == "increase"] <- "#D73027" # red
    E(g)$color[E(g)$direction == "decrease"] <- "#4575B4" # blue

    # # lty: default from status (optional, dashed for gain )
    E(g)$lty <- "solid"
    E(g)$lty[E(g)$status == "gained"] <- "dashed"

    # OVERRIDE: dashed if not in STRING PPIN
    if (!is.null(g_string)) {
        if (normalize_case) {
            igraph::V(g)$name <- toupper(igraph::V(g)$name)
            igraph::V(g_string)$name <- toupper(igraph::V(g_string)$name)
        }
        g <- set_lty_by_reference(g, g_string)
    }

    # width by abs_delta
    max_abs <- max(E(g)$abs_delta, na.rm = TRUE)
    if (is.finite(max_abs) && max_abs > 0) {
        E(g)$width <- 1 + 4 * (E(g)$abs_delta / max_abs)
    } else {
        E(g)$width <- 2
    }

    # ---- Labels: only top_n_label edges by abs_delta ----
    E(g)$label <- NA
    if (top_n_label > 0) {
        top_e <- order(E(g)$abs_delta, decreasing = TRUE)[seq_len(min(top_n_label, ecount(g)))]
        E(g)$label[top_e] <- round(E(g)$delta[top_e], 3)
        E(g)$label.cex <- 0.8
    }

    g
}


#' Cluster-wise CTS enrichment using Fisher's exact test
#'
#' @description
#' Tests over-representation of a CTS gene set in each cluster (eg GEne Program) 
#' using a 2×2 contingency table and Fisher's exact test
#' (equivalent to a hypergeometric test).
#'
#' @param df Data frame containing cluster and gene columns.
#' @param cts_genes Character vector of CTS genes.
#' @param cluster_col Column name in `df` defining clusters. Default `"Cluster"`.
#' @param gene_col Column name in `df` defining gene IDs. Default `"gene_name"`.
#' @param universe Optional background gene vector. If `NULL`,
#'   all unique genes in `df[[gene_col]]` are used.
#' @param adjust_method Multiple testing correction method passed to
#'   `p.adjust()`. Default `"BH"`.
#'
#' @details
#' The universe defines the background for enrichment.
#' Both cluster genes and CTS genes are restricted to the universe
#' before computing the 2×2 table.
#'
#' @return Data frame with one row per cluster containing:
#' \itemize{
#'   \item \code{GS}: cluster name
#'   \item \code{n_cluster_genes}: number of cluster genes (in universe)
#'   \item \code{n_CTS}: number of CTS genes (in universe)
#'   \item \code{n_overlap}: number of overlapping genes
#'   \item \code{overlap}: comma-separated overlapping genes
#'   \item \code{odds_ratio}: Fisher odds ratio
#'   \item \code{p_value}: Fisher p-value
#'   \item \code{adj_p}: adjusted p-value
#'   \item \code{neglog10_p}, \code{neglog10_adj_p}
#' }
#'
#' @examples
#' # Example
#' df <- data.frame(Cluster = c("A","A","B"),
#'                  gene_name = c("G1","G2","G2"))
#' cts_overlap_by_cluster(df, c("G2"))
#'
#' @export
cts_overlap_by_cluster <- function(df, cts_genes, cluster_col = "Cluster",
                                   gene_col = "gene_name", universe = NULL,
                                   adjust_method = "BH") {
  stopifnot(cluster_col %in% colnames(df), gene_col %in% colnames(df))

  # Clean inputs
  cts_genes <- unique(as.character(cts_genes))
  df[[gene_col]] <- as.character(df[[gene_col]])
  df[[cluster_col]] <- as.character(df[[cluster_col]])  # updated 4/27/2026

  # Define universe (background)
  # If you don't provide universe, use all unique genes observed in df
  if (is.null(universe)) {
    universe <- unique(df[[gene_col]])
  } else {
    universe <- unique(as.character(universe))
  }

  # Intersect CTS with universe #!!!!!!!!!!
  cts_in_universe <- intersect(cts_genes, universe)

  clusters <- sort(unique(df[[cluster_col]]))

  res <- lapply(clusters, function(cl) {
    genes_cl <- unique(df[[gene_col]][df[[cluster_col]] == cl])   # genes_cl <- unique(df[df[[cluster_col]] == cl, gene_col])
    genes_cl <- intersect(genes_cl, universe)  #!!!!!!!!!!
	
	A = intersect(genes_cl, cts_in_universe)
    a <- length(A)               # in cluster AND CTS
    b <- length(setdiff(genes_cl, cts_in_universe))                 # in cluster NOT CTS
    c <- length(setdiff(cts_in_universe, genes_cl))                 # CTS NOT in cluster
    d <- length(setdiff(universe, union(genes_cl, cts_in_universe)))# neither

    tab <- matrix(c(a, b, c, d), nrow = 2,
                  dimnames = list(Cluster = c("in_cluster", "not_in_cluster"),
                                  CTS = c("in_CTS", "not_in_CTS")))

    ft <- fisher.test(tab)

    data.frame(
      GS = cl,
      n_cluster_genes = length(genes_cl),
      n_CTS = length(cts_in_universe),
      n_overlap = a,
	  overlap = toString(A),
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, res)
  out$adj_p <- p.adjust(out$p_value, method = adjust_method)

  # Optional: add -log10 columns for convenience
  out$neglog10_p <- -log10(pmax(out$p_value, .Machine$double.xmin))
  out$neglog10_adj_p <- -log10(pmax(out$adj_p, .Machine$double.xmin))

  out[order(out$adj_p, out$p_value), ]
}



#############################################################################################
## internal functions, reused  for validation of CTS dual-pull in independent datasets
## ---------- helpers for annotation ----------

fmt_lm_term <- function(fit, term){
  ## ADDED
  weighted_sd <- function(x, w = NULL) {
    ok <- is.finite(x)
    if (!is.null(w)) ok <- ok & is.finite(w) & w > 0
    x <- x[ok]

    if (length(x) < 2) return(NA_real_)

    if (is.null(w)) {
      return(stats::sd(x))
    } else {
      w <- w[ok]
      wm <- sum(w * x) / sum(w)
      sqrt(sum(w * (x - wm)^2) / sum(w))
    }
  }
  co <- summary(fit)$coefficients
  if (!term %in% rownames(co)) return(NULL)

  b  <- co[term, "Estimate"] # effect size (coupling strength)
  t  <- co[term, "t value"]  # effect size weighted by precision, comparable within the same dataset and model structure.
  p  <- co[term, "Pr(>|t|)"]
  r2 <- summary(fit)$r.squared
 
 ## ADDED: compute standardized beta
  mf <- model.frame(fit)
  x  <- mf[[term]]
  y  <- model.response(mf)
  w  <- fit$weights

  sx <- weighted_sd(x, w)
  sy <- weighted_sd(y, w)

  beta_std <- if (is.finite(sx) && is.finite(sy) && sy > 0) { # comparable 
    b * sx / sy
  } else {
    NA_real_
  }

  list(
    beta = b,                 ## ADDED
    t = t,
    beta_std = beta_std,      ## ADDED
    slope = b,                ## keep for backward compatibility
    p = p,
    r2 = r2,
    label = sprintf("beta=%.2f, R²=%.2f, p=%s",  #t
                    b, r2, format.pval(p, digits=2, eps=1e-16))
  )
}



# Build annotation table: per diffday × celltype, fit y ~ x, keep significant
# precision_w was used by voom [PMID: 24485249]and Dreamlet [PMID: 37205331]
# where w=1/SEM², where SEM = sd(signal)/sqrt(n) and n is number of cells

  make_ann_by_day_celltype <- function(df_sub, xvar, yvar = "score_GPi",
                                      p_cut = 0.05, use_weights = FALSE,
                                      split_by_col = "diffday",
                                      color_by_col = "celltype"#,
                                      #weights_var = "n_cells"
                                      ) {

    # ensure factors are droplevelled so split works cleanly
    df_sub[[split_by_col]] = droplevels(factor(df_sub[[split_by_col]]))
    df_sub[[color_by_col]] = droplevels(factor(df_sub[[color_by_col]]))

    split_by_col_sym <- rlang::sym(split_by_col)
    color_by_col_sym <- rlang::sym(color_by_col)
    # per facet (day) and per celltype regression
    ann <- df_sub %>%
      dplyr::group_by(!!split_by_col_sym, !!color_by_col_sym) %>%
      dplyr::group_modify(function(d, ...) {
        # need at least 2 points and some x variation
        if (nrow(d) < 2) return(data.frame())
        if (sd(d[[xvar]], na.rm = TRUE) == 0) return(data.frame())

        # fit <- if (use_weights && weights_var %in% names(d)) {
          #lm(reformulate(xvar, response = yvar), data = d, weights = d[[weights_var]])
        fit <- if(use_weights) {
          w = 1 / (d[[paste0(yvar, "_sem")]]^2)
          w[!is.finite(w)] <- NA
          d$w = w
          lm(reformulate(xvar, response = yvar), data = d, weights = w)
        } else {
          lm(reformulate(xvar, response = yvar), data = d)
        }

        out <- fmt_lm_term(fit, term = xvar)
        if (is.null(out)) return(data.frame())
        if (is.na(out$p) || out$p >= p_cut) return(data.frame())

        # place label near top-left of each facet, and stack by celltype
        x_left <- quantile(d[[xvar]], 0.02, na.rm = TRUE)
        y_top  <- quantile(d[[yvar]], 0.98, na.rm = TRUE)

        data.frame(
          x = x_left,
          y = y_top,
          label = out$label,
          beta = out$beta,
          beta_std = out$beta_std,  # standardize by the standard deviation of the response variable,
          t = out$t,
          slope = out$slope,          ## keep for compatibility, the same as beta
          p = out$p,
          r2 = out$r2
        )
      }) %>%
      dplyr::ungroup()

    # stack multiple celltype labels vertically within each facet
    if (nrow(ann) > 0) {
      ann <- ann %>%
      dplyr::group_by(!!split_by_col_sym) %>%
      dplyr::arrange(!!color_by_col_sym, .by_group = TRUE) %>%
      dplyr::mutate(y = y - (dplyr::row_number() - 1) * 0.08 * diff(range(df_sub[[yvar]], na.rm=TRUE))) %>%
      dplyr::ungroup()
  }

  ann
}


# Create a helper function to reduce repetition
# for simple pseudobulk association. However,  HOT used 
make_lineage_pseudobulk_plot <- function(df_sub, lineage, 
    pull_genes, pull_score_col, 
    CTS_name, 
    yvar = "score_GPi", 
    gp_present ,
    use_wt = FALSE, 
    dot_size = 1.5, 
    split_by_col = 'diffday', 
    color_by_col = 'celltype', 
    cols = NULL, 
    x_name = NULL, 
    y_name = NULL ) { 

     # pull_score_sem_col = paste0(pull_score_col, "_sem") 
      lab_prefix = lineage 
      if(length(pull_genes) > 0) { 
        ann <- make_ann_by_day_celltype( 
          df_sub, 
          xvar = pull_score_col, 
          yvar = yvar, 
          p_cut = 0.05, 
          use_weights = use_wt, 
          split_by_col = split_by_col ) 
          
          # For dual, the SEM is on y-axis; others, on x. # << tidy-eval symbols (replaces aes_string) 
          x_sym <- rlang::sym(pull_score_col) 
          y_sym <- rlang::sym("score_GPi") 
          color_by_sym <- rlang::sym(color_by_col) 
          p_tmp = ggplot(df_sub, ggplot2::aes(x = !!x_sym, y = !!y_sym, color = !!color_by_sym)) + 
               geom_point(size = dot_size) + 
               geom_smooth(method = "lm", 
                se = FALSE, aes(weight = if (use_wt) 1/(score_GPi_sem^2) else NULL)) + 
               scale_color_manual(values = cols, breaks = names(cols)) + 
               theme_bw() + theme(legend.position = "none") + 
               facet_wrap(as.formula(paste("~", split_by_col)), nrow = 1 ) + # facet_wrap(~ diffday, nrow=1) + updated for Hill2022 CHD dataset !!!!!!!!!!!!!!!!!!!!!! 
               labs( x = if (!is.null(x_name)) x_name else paste0(length(unique(pull_genes)), "-gene ", lab_prefix, "-pull score: ", toString(pull_genes)), 
                     y = if (!is.null(y_name)) y_name else paste0(length(unique(gp_present)), "-gene", id, " signature score ") ) + 
                ggtitle(paste0('CTS_', CTS_name, ": ", dot_type, ' level', id)) + 
                coord_cartesian(clip = "off") + 
                theme(plot.margin = margin(5.5, 5.5, 20, 5.5)) 
                
          if (is.null(ann) || nrow(ann) == 0 || !all(c("x","y","label") %in% colnames(ann))) { 
                 p_tmp <- p_tmp # do nothing (no annotations) 
          } else { 
                  p_tmp <- p_tmp + geom_text( data = ann, aes(x = x, y = y, label = label, color = split), # color = celltype), updated for Hill2022 CHD dataset !!!!!!!!!!!!!!!!!!!!!! 
                  inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.5 ) 
                  } 
        } else p_tmp = ggplot() + theme_void() + labs(title = paste0('CTS_',CTS_name, ": ", lineage, '-bias -- No data')) 
        
        p_tmp 
}


## helper to format beta/angle text from lm()
fmt_lm <- function(fit, termlabels = "lineage_bias"){
  co <- summary(fit)$coefficients
  if (!termlabels %in% rownames(co)) return(NA_character_)

  p  <- co[termlabels, "Pr(>|t|)"]
  r2 <- summary(fit)$r.squared

 sprintf("β=%.3f, p=%s",
        co[termlabels, "Estimate"],
        format.pval(p, digits=2))

#   sprintf("R²=%.2f, p=%s",
#           r2, format.pval(p, digits=2, eps=1e-16))
}

 
# fmt_lm_by <- function(d, termlabels = "lineage_bias", response = "score_GPi", weights=NULL){
#   fit <- lm(reformulate(termlabels = termlabels, response = response, weights = weights), data = d )
#   fmt_lm(fit, termlabels = termlabels)
# }



## testing across all clusters/main cell types, similar resutls, 
## called by xxx_evaluate_with_Hill2022_CHD_wdual_pull_PatientID_clean.R
get_endpoint_effect <- function(
  sce,
  cf_col,
  cm_col,
  cell_clustering_col = c('Cluster', 'MainCellType'),
  cluster = "CF4",
  group_col = "Diagnosis",
  disease = "Neo_HLHS",
  control_level = "Donor",
  pseudobulk_level = c("biosample_id", "patientID"),
  min_cells = 3,
  scale_outcome = FALSE,
  adjust_region = TRUE
) {
  pseudobulk_level <- match.arg(pseudobulk_level)
  cell_clustering_col = match.arg(cell_clustering_col)
  if(cell_clustering_col == 'Cluster') {
    dat <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
        dplyr::filter(.data$Cluster == cluster, !is.na(.data[[group_col]])) %>%
        dplyr::transmute(
        patientID   = factor(patientID),
        biosample_id = factor(biosample_id),
        region      = factor(region),
        group       = factor(.data[[group_col]]),
        imbalance   = .data[[cm_col]] - .data[[cf_col]]
        ) %>%
        droplevels()
  } else {
  dat <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
        dplyr::filter(.data$MainCellType == cluster, !is.na(.data[[group_col]])) %>%
        dplyr::transmute(
        patientID   = factor(patientID),
        biosample_id = factor(biosample_id),
        region      = factor(region),
        group       = factor(.data[[group_col]]),
        imbalance   = .data[[cm_col]] - .data[[cf_col]]
        ) %>%
        droplevels()
  }

  if (!control_level %in% levels(dat$group)) {
    stop("control_level not found in ", group_col)
  }

  dat$group <- stats::relevel(dat$group, ref = control_level)

  # pseudobulk
  if (pseudobulk_level == "biosample_id") {
    pb <- dat %>%
      dplyr::group_by(biosample_id, patientID, group, region) %>%
      dplyr::summarise(
        imbalance = mean(imbalance, na.rm = TRUE),
        n_cells = dplyr::n(),
        .groups = "drop"
      )
  } else {
    pb <- dat %>%
      dplyr::group_by(patientID, group, region) %>%
      dplyr::summarise(
        imbalance = mean(imbalance, na.rm = TRUE),
        n_cells = dplyr::n(),
        .groups = "drop"
      )
  }

  pb <- pb %>%
    dplyr::filter(n_cells >= min_cells) %>%
    droplevels()

  if ( dplyr::n_distinct(pb$group) < 2) return(NULL)
  if (!disease %in% levels(pb$group)) {
    return(data.frame(
      yi = NA_real_, sei = NA_real_, p = NA_real_,
      term = paste0("group", disease),
      disease = disease,
      control = control_level,
      group_col = group_col,
      cluster = cluster,
      n_units = nrow(pb),
      n_patients = dplyr::n_distinct(pb$patientID),
      n_case_units = 0,
      n_ctrl_units = sum(pb$group == control_level),
      pseudobulk_level = pseudobulk_level,
      adjust_region = adjust_region,
      scale_outcome = scale_outcome,
      note = "disease level absent after filtering"
    ))
  }

  # Optional scaling for cross-signature comparability
  # This changes coefficient units, but not t/p if rows and model are identical.
  pb$y <- pb$imbalance
  if (scale_outcome) {
    sdy <- sd(pb$y, na.rm = TRUE)
    if (!is.finite(sdy) || sdy == 0) return(NULL)
    pb$y <- pb$y / sdy
  }

  # same old-style full model
  if (adjust_region && dplyr::n_distinct(pb$region) > 1) {
    fit <- stats::lm(y ~ group + region, data = pb)
    model_used <- "y ~ group + region"
  } else {
    fit <- stats::lm(y ~ group, data = pb)
    model_used <- "y ~ group"
  }

  # same variance logic as before
  if (pseudobulk_level == "biosample_id") {
    ct <- lmtest::coeftest(
      fit,
      vcov. = sandwich::vcovCL(fit, cluster = pb$patientID)
    )
  } else {
    ct <- lmtest::coeftest(
      fit,
      vcov. = sandwich::vcovHC(fit, type = "HC1")
    )
  }

  co <- unclass(ct)
  pcol <- grep("^Pr\\(", colnames(co), value = TRUE)[1]
  term <- paste0("group", disease)

  if (!term %in% rownames(co)) {
    return(data.frame(
      yi = NA_real_, sei = NA_real_, p = NA_real_,
      term = term,
      disease = disease,
      control = control_level,
      group_col = group_col,
      cluster = cluster,
      n_units = nrow(pb),
      n_patients = dplyr::n_distinct(pb$patientID),
      n_case_units = sum(pb$group == disease),
      n_ctrl_units = sum(pb$group == control_level),
      pseudobulk_level = pseudobulk_level,
      adjust_region = adjust_region,
      scale_outcome = scale_outcome,
      model = model_used,
      note = "term not estimable / aliased"
    ))
  }

  data.frame(
    yi = co[term, "Estimate"],
    sei = co[term, "Std. Error"],
    p = co[term, pcol],
    term = term,
    disease = disease,
    control = control_level,
    group_col = group_col,
    cluster = cluster,
    n_units = nrow(pb),
    n_patients = dplyr::n_distinct(pb$patientID),
    # n_case_units = sum(pb$group == disease),
    # n_ctrl_units = sum(pb$group == control_level),
    # n_case_patients = dplyr::n_distinct(pb$patientID[pb$group == disease]),
    # n_ctrl_patients = dplyr::n_distinct(pb$patientID[pb$group == control_level]),
    n_case_biosamples = if (pseudobulk_level == "biosample_id")
      dplyr::n_distinct(pb$biosample_id[pb$group == disease]) else NA_integer_,
    n_ctrl_biosamples = if (pseudobulk_level == "biosample_id")
      dplyr::n_distinct(pb$biosample_id[pb$group == control_level]) else NA_integer_,
    n_case_units = sum(pb$group == disease),
    n_ctrl_units = sum(pb$group == control_level),
    n_case_patients = dplyr::n_distinct(pb$patientID[pb$group == disease]),
    n_ctrl_patients = dplyr::n_distinct(pb$patientID[pb$group == control_level]),
    pseudobulk_level = pseudobulk_level,
    adjust_region = adjust_region,
    scale_outcome = scale_outcome,
    model = model_used,
    note = NA_character_
  )
}


## testing for a particular cluster /main cell type, used in the paper
get_cf4_effect <- function(
  sce,
  cf_col,
  cm_col,
  cluster = "CF4",
  group_col = "Diagnosis",
  case_level = "Neo_HLHS",
  control_level = "Donor",
  min_cells = 3,
  pseudobulk_level = c("patientID", "biosample_id", "library_name", 'T_Reps', 'B_Reps'),
  scale_outcome = FALSE,
  adjust_region = TRUE
) {
  pseudobulk_level <- match.arg(pseudobulk_level)

  dat <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::filter(.data$Cluster == cluster, !is.na(.data[[group_col]])) %>%
    dplyr::mutate(
      group = factor(.data[[group_col]], levels = c(control_level, case_level)),
      imbalance = .data[[cm_col]] - .data[[cf_col]]
    ) %>%
    dplyr::filter(group %in% c(control_level, case_level))

  # pseudobulk
  if (pseudobulk_level == "patientID") {
    dat <- dat %>%
      dplyr::group_by(patientID, group, region) %>%
      dplyr::summarise(
        imbalance = mean(imbalance, na.rm = TRUE),
        n_cells = dplyr::n(),
        .groups = "drop"
      )
  } else if (pseudobulk_level == "biosample_id") {
    dat <- dat %>%
      dplyr::group_by(biosample_id, patientID, group, region) %>%
      dplyr::summarise(
        imbalance = mean(imbalance, na.rm = TRUE),
        n_cells = dplyr::n(),
        .groups = "drop"
      )
  } else if (pseudobulk_level == "library_name") {
    dat <- dat %>%
      dplyr::group_by(library_name, Sample, group) %>%
      dplyr::summarise(
        imbalance = mean(imbalance, na.rm = TRUE),
        n_cells = dplyr::n(),
        .groups = "drop"
      )
  }
  else if (pseudobulk_level == "T_Reps") {
    dat <- dat %>%
      dplyr::group_by(T_Reps, B_Reps, Sample, group) %>%
      dplyr::summarise(
        imbalance = mean(imbalance, na.rm = TRUE),
        n_cells = dplyr::n(),
        .groups = "drop"
      )
  } else if (pseudobulk_level == "B_Reps") {
      dat <- dat %>%
      dplyr::group_by(B_Reps, Sample, group) %>%
      dplyr::summarise(
        imbalance = mean(imbalance, na.rm = TRUE),
        n_cells = dplyr::n(),
        .groups = "drop"
      )
  }  

  dat <- dat %>% dplyr::filter(n_cells >= min_cells)
  if (dplyr::n_distinct(dat$group) < 2) return(NULL)

  # optional scaling for meta-analysis comparability
  if (scale_outcome) {
    sdx <- sd(dat$imbalance, na.rm = TRUE)
    if (!is.finite(sdx) || sdx == 0) return(NULL)
    dat <- dat %>% dplyr::mutate(imbalance_out = as.numeric(scale(imbalance)))
  } else {
    dat <- dat %>% dplyr::mutate(imbalance_out = imbalance)
  }

  if('region' %in% colnames(dat))  dat$region <- factor(dat$region)

    # full model if region is estimable
  form <- if (adjust_region && 'region' %in% colnames(dat) && dplyr::n_distinct(dat$region) > 1) {
      imbalance_out ~ group + region
    } else {
      imbalance_out ~ group
  }

  X <- model.matrix(form, data = dat)
  # check nrow(X) to make sure the design matrix is full rank, i.e. that the coefficients in lm(form, data = dat) are actually estimable and not aliased by collinearity or an overparameterized design.
  if (nrow(X) <= ncol(X) || Matrix::rankMatrix(X)[1] < ncol(X)) {
    form <- imbalance_out ~ group
    X <- model.matrix(form, data = dat)
    if (nrow(X) <= ncol(X) || Matrix::rankMatrix(X)[1] < ncol(X)) return(NULL)
  }

  fit <- stats::lm(form, data = dat)

  # patient-level: ordinary lm summary
  # biosample-level: cluster-robust SE by patientID
  if (pseudobulk_level %in% c("patientID", "library_name","B_Reps")) {
    co <- summary(fit)$coefficients
  } else if (pseudobulk_level == "biosample_id") {
    ct <- lmtest::coeftest(
      fit,
      vcov. = sandwich::vcovCL(fit, cluster = dat$patientID)
    )
    co <- unclass(ct)
  } else if (pseudobulk_level == "T_Reps") {
    ct <- lmtest::coeftest(
      fit,
      vcov. = sandwich::vcovCL(fit, cluster = dat$B_Reps)
    )
    co <- unclass(ct)
  }

  term <- paste0("group", case_level)
  if (!term %in% rownames(co)) return(NULL)

  p_col <- grep("^Pr\\(", colnames(co), value = TRUE)[1]

 if (pseudobulk_level == "patientID") {
   out = data.frame( 
      yi = co[term, "Estimate"],
      sei = co[term, "Std. Error"],
      p = co[term, p_col],
      n_case_patientID = dplyr::n_distinct(dat$patientID[dat$group == case_level]) , 
      n_ctrl_patientID = dplyr::n_distinct(dat$patientID[dat$group == control_level])  ,
      n_case_units = sum(dat$group == case_level),
      n_ctrl_units = sum(dat$group == control_level),
      pseudobulk_level = pseudobulk_level,
      scale_outcome = scale_outcome
      )
  } else if (pseudobulk_level == "library_name") {
    out = data.frame( 
      yi = co[term, "Estimate"],
      sei = co[term, "Std. Error"],
      p = co[term, p_col],
      n_case_library_name = dplyr::n_distinct(dat$library_name[dat$group == case_level]) ,
      n_ctrl_library_name = dplyr::n_distinct(dat$library_name[dat$group == control_level]) ,
      n_case_units = sum(dat$group == case_level),
      n_ctrl_units = sum(dat$group == control_level),
      pseudobulk_level = pseudobulk_level,
      scale_outcome = scale_outcome
      )
  } else if(pseudobulk_level == "biosample_id") {
    out = data.frame(
      yi = co[term, "Estimate"],
      sei = co[term, "Std. Error"],
      p = co[term, p_col],
      n_case_biosamples = dplyr::n_distinct(dat$biosample_id[dat$group == case_level]),
      n_ctrl_biosamples =  dplyr::n_distinct(dat$biosample_id[dat$group == control_level]) ,
      n_case_patients = dplyr::n_distinct(dat$patientID[dat$group == case_level]),
      n_ctrl_patients = dplyr::n_distinct(dat$patientID[dat$group == control_level]),
      n_case_units = sum(dat$group == case_level),
      n_ctrl_units = sum(dat$group == control_level),
      # n_case_units = sum(dat$group == case_level),
      # n_ctrl_units = sum(dat$group == control_level),
      # n_case = dplyr::n_distinct(dat$patientID[dat$group == case_level]),
      # n_ctrl = dplyr::n_distinct(dat$patientID[dat$group == control_level]),
      pseudobulk_level = pseudobulk_level,
      scale_outcome = scale_outcome
    )
  } else if(pseudobulk_level %in% c("T_Reps", "B_Reps")) {
out = data.frame( 
      yi = co[term, "Estimate"],
      sei = co[term, "Std. Error"],
      p = co[term, p_col],
      n_case_T_Reps = if (pseudobulk_level == "T_Reps")
        dplyr::n_distinct(dat$T_Reps[dat$group == case_level]) else NA_integer_,
      n_ctrl_T_Reps = if (pseudobulk_level == "T_Reps")
        dplyr::n_distinct(dat$T_Reps[dat$group == control_level]) else NA_integer_,
      n_case_B_Reps = if (pseudobulk_level == "B_Reps")
        dplyr::n_distinct(dat$B_Reps[dat$group == case_level]) else NA_integer_,
      n_ctrl_B_Reps = if (pseudobulk_level == "B_Reps")
        dplyr::n_distinct(dat$B_Reps[dat$group == control_level]) else NA_integer_,
      n_case_units = sum(dat$group == case_level),
      n_ctrl_units = sum(dat$group == control_level),
      pseudobulk_level = pseudobulk_level,
      scale_outcome = scale_outcome
      )
  } 
  return(out)
}


## Helper function to get the signature genes
get_signature_genes <- function(pull_df, id) {
  add_lineage_specifically_edge_increased_cm = grepl("wincreased", id) 
  cf0 <- pull_df %>%
    dplyr::filter(ID == id,  linkeage == "CF" ) %>%
    dplyr::select(x, y) %>% unlist() %>% unique()

  cm0 <- pull_df %>%
    dplyr::filter(ID == id,  linkeage == "CM" ) %>%
    dplyr::select(x, y) %>% unlist() %>% unique()

  dual <- intersect(cf0, cm0)
  cf <- setdiff(cf0, dual)
  cm <- setdiff(cm0, dual)

  increased =  pull_df %>%
    dplyr::filter(ID == id,  direction == "increase" ) %>%
    dplyr::select(x, y) %>% unlist() %>% unique() %>% intersect(dual,.)

  if (add_lineage_specifically_edge_increased_cm) {
  ## for the hiPSC dataset, add the score of CM_pull that adding the increased nodes (FGF10, LRRTM1) from dual-pull set
    if(id %in% c( 'Elorbany_1v1.25_STRING_CP.1wincreased', 'Elorbany_1v1.25_IID_CP.1wincreased')) {
          cm  = c(cm , 'FGF10', 'LRRTM1')
             }  
    # ## for the Ibarra dataset, add the score without TBX5 for evaluation only
    # if(id == 'Ibarra_1v1.25_STRING_cardiac.anoTBX5')  {
		#   cm  = setdiff(cm , 'TBX5')
	  #    }  
  }
  list(cf = cf, cm = cm)
}

# Helper to sample matched random genes
# 'bin' here refers to a grouping of genes based on their average expression levels. 
# Genes are binned into discrete intervals (bins) according to quantiles of their mean expression.
# This allows sampling of control genes matched to the target gene set by expression level,
# thus controlling for gene expression differences.
# when n_bins <= 1 means plain random sampling
sample_matched_genes <- function(target_genes, pool, bins = NULL, exclude = character()) {
  target_genes <- unique(target_genes)
  if (length(target_genes) == 0) stop("No target genes available for sampling.")

  pool_use <- setdiff(pool, exclude)

  ## plain random sampling
  if (is.null(bins)) {
    if (length(pool_use) < length(target_genes)) {
      stop("Not enough genes available for random sampling.")
    }
    return(sample(pool_use, length(target_genes), replace = FALSE))
  }

  ## expression-matched sampling
  target_genes <- intersect(target_genes, names(bins))
  if (length(target_genes) == 0) stop("No target genes found in bins.")

  need <- table(bins[target_genes])
  out <- character(0)

  for (b in names(need)) {
    k <- as.integer(need[[b]])
    cand <- setdiff(pool[bins[pool] == as.integer(b)], c(exclude, out))

    ## fallback if that bin is too small
    if (length(cand) < k) cand <- setdiff(pool, c(exclude, out))
    if (length(cand) < k) stop("Not enough genes available for matched sampling.")

    out <- c(out, sample(cand, k, replace = FALSE))
  }

  out
}

############################################################
## - use get_cf4_effect() for BOTH observed and null effects
## - keeps missing handling and counting consistent
## - returns biosample counts, patient counts, and observation counts
## - adds emp_p so your old accessor still works
############################################################
run_random_null_test <- function(
  sce,
  cf_genes,
  cm_genes,
  obs_cf_col,
  obs_cm_col,
  cluster = "CF",
  case_level = "HLHS",
  group_col = "Diagnosis2",
  control_level = "Donor",
  pseudobulk_level = "biosample_id",
  min_cells = 3,
  scale_outcome = FALSE,
  adjust_region = TRUE,
  n_perm = 200,
  n_bins = 10,
  seed = 1
) {
  set.seed(seed)

  lc <- SingleCellExperiment::logcounts(sce)
  all_genes <- rownames(sce)

  cf_genes <- intersect(cf_genes, all_genes)
  cm_genes <- intersect(cm_genes, all_genes)

  ## exclude the real signature genes from the null pool
  pool <- setdiff(all_genes, union(cf_genes, cm_genes))

  ## CHANGED: only build bins when n_bins > 1
  bins <- NULL
  if (!is.null(n_bins) && n_bins > 1) {
    gene_means <- Matrix::rowMeans(lc)
    brks <- unique(stats::quantile(
      gene_means,
      probs = seq(0, 1, length.out = n_bins + 1),
      na.rm = TRUE
    ))
    if (length(brks) > 2) {
      bins <- cut(gene_means, breaks = brks, include.lowest = TRUE, labels = FALSE)
      names(bins) <- all_genes
    }
  }

  ## observed effect
  obs <- get_cf4_effect(
    sce = sce,
    cf_col = obs_cf_col,
    cm_col = obs_cm_col,
    cluster = cluster,
    group_col = group_col,
    case_level = case_level,
    control_level = control_level,
    pseudobulk_level = pseudobulk_level,
    min_cells = min_cells,
    scale_outcome = scale_outcome,
    adjust_region = adjust_region
  )

  if (is.null(obs)) stop("Observed effect is not estimable under current settings.")

  null_res <- vector("list", n_perm)

  for (i in seq_len(n_perm)) {
    cf_rand <- sample_matched_genes(
      target_genes = cf_genes,
      pool = pool,
      bins = bins
    )
    cm_rand <- sample_matched_genes(
      target_genes = cm_genes,
      pool = pool,
      bins = bins,
      exclude = cf_rand
    )

    SummarizedExperiment::colData(sce)$CF_rand <- Matrix::colMeans(lc[cf_rand, , drop = FALSE])
    SummarizedExperiment::colData(sce)$CM_rand <- Matrix::colMeans(lc[cm_rand, , drop = FALSE])

    tmp <- get_cf4_effect(
      sce = sce,
      cf_col = "CF_rand",
      cm_col = "CM_rand",
      cluster = cluster,
      group_col = group_col,
      case_level = case_level,
      control_level = control_level,
      pseudobulk_level = pseudobulk_level,
      min_cells = min_cells,
      scale_outcome = scale_outcome,
      adjust_region = adjust_region
    )
    if (!is.null(tmp)) {                                        # <<< CHANGED
      if (!"n_case_obs" %in% names(tmp) &&  "n_case_units" %in% names(tmp))  tmp$n_case_obs <- tmp$n_case_units                     # <<< CHANGED
      if (!"n_ctrl_obs" %in% names(tmp) &&  "n_ctrl_units" %in% names(tmp))   tmp$n_ctrl_obs <- tmp$n_ctrl_units      
    }

    null_res[[i]] <- if (is.null(tmp)) {
      data.frame(
        iter = i,
        yi = NA_real_,
        sei = NA_real_,
        p = NA_real_,
        n_case_biosamples = NA_integer_,
        n_ctrl_biosamples = NA_integer_,
        n_case_obs = NA_integer_,
        n_ctrl_obs = NA_integer_,
        n_case_patients = NA_integer_,
        n_ctrl_patients = NA_integer_
      )
    } else {
      cbind(
        data.frame(iter = i),
        tmp[, c(
          "yi", "sei", "p",
          "n_case_biosamples", "n_ctrl_biosamples",
          "n_case_obs", "n_ctrl_obs",
          "n_case_patients", "n_ctrl_patients"
        ), drop = FALSE]
      )
    }
  }

  null_df <- dplyr::bind_rows(null_res) %>%
    dplyr::filter(is.finite(yi))

  if (nrow(null_df) == 0) stop("All null iterations failed.")

  ## CHANGED: add +1 correction
  emp_p_left  <- (1 + sum(null_df$yi <= obs$yi, na.rm = TRUE)) / (nrow(null_df) + 1)
  emp_p_right <- (1 + sum(null_df$yi >= obs$yi, na.rm = TRUE)) / (nrow(null_df) + 1)
  emp_p_two   <- (1 + sum(abs(null_df$yi) >= abs(obs$yi), na.rm = TRUE)) / (nrow(null_df) + 1)

  list(
    obs = obs,
    null = null_df,
    emp_p_left = emp_p_left,
    emp_p_right = emp_p_right,
    emp_p_two = emp_p_two,

    ## CHANGED: keeps your old style working
    emp_p = if (obs$yi < 0) emp_p_left else emp_p_right
  )
}


## plot the null distribution and the observed effect ##########
  plot_null_test <- function(res, title = "") {
  obs_y <- res$obs$yi
  emp_p <- if (obs_y < 0) res$emp_p_left else res$emp_p_right

  ggplot(res$null, aes(x = yi)) +
    geom_histogram(bins = 40, fill = "grey80", color = "white") +
    geom_vline(xintercept = obs_y, color = "red", linewidth = 1) +
    annotate(
      "text",
      x = obs_y,
      y = Inf,
      label = paste0("Observed = ", round(obs_y, 2),
                     "\nEmpirical p = ", signif(emp_p, 3)),
      vjust = 1.2,
      hjust = ifelse(obs_y < median(res$null$yi, na.rm = TRUE), 0, 1),
      color = "red",
      size = 4
    ) +
    labs(
      x = "Random-signature effect (SD units)",
      y = "Count",
      title = title
    ) +
    theme_classic(base_size = 14)
}

# # a helper function to plot the regression of the GO term score on the CF and CM pull scores
# # regression on cell-level only (excluding pseudobulk level) when db == 'Elorbany' 
# celltype_regression_scatterplot = function(sce, pull_df, db, CTS_name, id,
#       CP_cluster, CM_cluster, CF_cluster, celltype_colors, go_names, celltype_col, data_logcounts=NULL) {
#     if(is.null(data_logcounts)) data_logcounts = logcounts(sce)
#     dot_size = ifelse(db=='Elorbany', 2, 0.6)

#     ID = make.names(paste0(id,'_', CTS_name))

#     id_genes = go_gene_lists[[id]]
# 		id_present <- intersect(id_genes, rownames(sce))
#         n_id_present = length(id_present)
# 		colData(sce)[[ID]] <- Matrix::colMeans(data_logcounts[id_present, , drop = FALSE])
# 		summary(colData(sce)[[ID]])
# 		# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 		# -0.7992 -0.2982 -0.1911  0.0000 -0.1103 13.8979 

# 		cf_pull_genes_0 <- subset(pull_df, Database==db & linkeage=='CF'  & CTS_ID==CTS_name)[, c('x','y')] %>%
# 					unlist() %>% unique()
# 		cm_pull_genes_0 <- subset(pull_df, Database==db & linkeage=='CM'  & CTS_ID==CTS_name)[, c('x','y')] %>%
# 				unlist() %>% unique()
# 		dual_pull_genes = intersect(cf_pull_genes_0, cm_pull_genes_0 )		
# 		cf_pull_genes = setdiff(cf_pull_genes_0, dual_pull_genes)  #!!!ADDED
# 		cm_pull_genes = setdiff(cm_pull_genes_0, dual_pull_genes)
		
#         if(db == 'Elorbany') cm_pull_genes = c(cm_pull_genes, 'FGF10', 'LRRTM1') 
#         if(db == 'Pijuan_Sala') cm_pull_genes = c(cm_pull_genes, 'HMGA2')

# 		cf_pull_name = paste0('CF_pull', CTS_name)			#!!!ADDED
# 		if(db %in% c('Elorbany', 'Pijuan_Sala')) cm_pull_name = paste0('CM_pull_wincreased_nodes', CTS_name)		else cm_pull_name = paste0('CM_pull', CTS_name)		

# 		colData(sce)[[cf_pull_name]] <- Matrix::colMeans(data_logcounts[cf_pull_genes, , drop = FALSE])
# 		colData(sce)[[cm_pull_name]] <- Matrix::colMeans(data_logcounts[cm_pull_genes, , drop = FALSE])
			
				
# 		#### (pseudobulk, donor-aware)
# 		#  aggregate by individual × day × celltype.
# 		df_score <- as.data.frame(colData(sce)) %>%
# 		  dplyr::mutate(celltype = .data[[celltype_col]]) %>%
# 		  dplyr::mutate(score_GOi = .data[[ID]]) %>%
# 		  dplyr::mutate(cf_pull_score = .data[[cf_pull_name]]) %>%
# 		  dplyr::mutate(cm_pull_score = .data[[cm_pull_name]]) 
      
#     if(db =='Pijuan_Sala') {
#         df_score = df_score %>% dplyr::group_by(sample, celltype) 	 
#         df_score$individual <- factor(df_score$sample)

# 		    df_sub <- subset(df_score, celltype %in% c(CP_cluster, CM_cluster, CF_cluster))
# 		    head(df_sub, 3)

# 		    fit <- lm(score_GOi ~ cm_pull_score + cf_pull_score + celltype + individual , data = df_sub)
# 		    summary(fit)
#       } else if(db == 'IbarraSoria') {
#         df_sub <- subset(df_score, celltype %in% c(CP_cluster, CM_cluster, CF_cluster))
#         fit <- lm(score_GOi ~ celltype , data = df_sub)
#         summary(fit)
#       } else {
#         df_score = df_score %>% dplyr::group_by(individual, diffday, celltype)  %>%
#             dplyr::summarise(
#               n_cells = dplyr::n(),
#               score_GOi     = mean(score_GOi, na.rm = TRUE),
#               cf_pull_score = mean(cf_pull_score, na.rm = TRUE),
#               cm_pull_score = mean(cm_pull_score, na.rm = TRUE),
#               .groups = "drop"
#             ) 
#         df_score$individual <- factor(df_score$individual)
#         df_score$diffday <- factor(df_score$diffday)
# 		    df_sub <- subset(df_score, celltype %in% c(CP_cluster, CM_cluster, CF_cluster))
#         df_sub$celltype <- factor(df_sub$celltype, levels =  c(CP_cluster, CM_cluster, CF_cluster))
# 		    head(df_sub, 3)

# 		    fit <- lm(score_GOi ~ cm_pull_score + cf_pull_score + celltype + individual + diffday , data = df_sub)
# 		    summary(fit)
#       }
        
#     if(length(cf_pull_genes)>0) {
# 		fit2 <- lm(score_GOi ~ cf_pull_score ,data = df_sub)
#         lab_reg <- paste0("Overall: ", fmt_lm(fit2, 'cf_pull_score'))
#        # per-celltype regressions (often useful given your colored lines)
#         lab_reg_by <- df_sub %>%
#             dplyr::mutate(celltype = droplevels(celltype)) %>%   # <-- key
#             split(.$celltype) %>%
#             lapply(function(d) {
#                 paste0(unique(d$celltype), ": ", fmt_lm_by(d, 'cf_pull_score'))
#             }) %>%
#             unlist() %>%
#             paste(collapse = "\n")

# 		lineage = 'CF'
#             # plot limits for placing annotation
#             x_lim <- quantile(df_sub$cf_pull_score, c(0.01, 0.99), na.rm = TRUE)
#             y_lim <- range(df_sub$score_GOi,    na.rm = TRUE)
#             y_top <- y_lim[2] + 0.18 * diff(y_lim)      ## CHANGED: extra headroom for text
#             y_txt1 <- y_lim[2] + 0.16 * diff(y_lim)     ## CHANGED: 1st line higher
#             y_txt2 <- y_lim[2] + 0.08 * diff(y_lim)     ## CHANGED: 2nd line higher

       
#             tmp <- ggplot(df_sub,  aes(x = cf_pull_score, y = score_GOi, color = celltype)) +
#                 geom_point(size = dot_size, alpha = 0.8) +
#                 geom_smooth(method = "lm", se = FALSE) +   # , aes(weight = if (use_wt) 1/(score_GPi_sem^2) else NULL)   use_wt = FALSE here
#                 scale_color_manual(values = celltype_colors, breaks = names(celltype_colors), labels=c('CP','CM','EMT')) +
#                 theme_bw() +
#                 labs(# x = paste0( length(cf_pull_genes), "-gene CF-pull over ",   length(cm_pull_genes), "-gene CM-pull score (cell-level)" ),
# 					          x = paste0( length(cf_pull_genes), "-gene CF-pull score "),
#                     y = paste0(n_id_present, "-gene ",  id, " score ")) +
#                 ggtitle(paste0("CTS_", CTS_name, '\n',go_names[id])) +
#                 coord_cartesian(ylim = c(y_lim[1], y_top), clip = "off") +
#                 annotate("text", 
#                         x = x_lim[1], y = y_txt1,
#                         label = lab_reg, hjust = 0, vjust = 1, size = 2.5) +
#                 annotate("text",
#                         x = x_lim[1], y = y_txt2,
#                         label = lab_reg_by, hjust = 0, vjust = 1, size = 2.5) +
#                 geom_abline(intercept = 0, beta = 1, linetype = "dashed", color = "black")
#             if(db=='Elorbany')  {
#               cat('pls modify 24.4.2_overlap_with_Takeuchi2026_v4_wdual_pull.R, which calls pseudobulk_regression_scatterplot()')
#               p_tmp1 =   ggplot() + theme_void() + labs(title = paste0('CTS_',CTS_name, ": ", lineage, '-bias -- No data'))
#             }
# 			}  else p_tmp1 =   ggplot() + theme_void() + labs(title = paste0('CTS_',CTS_name, ": ", lineage, '-bias -- No data'))
			 
# 	if(length(cm_pull_genes)>0) {
# 		fit3 <- lm(score_GOi ~ cm_pull_score , data = df_sub)
#         lab_reg <- paste0("Overall: ", fmt_lm(fit3, 'cm_pull_score'))
# 		# per-celltype regressions (often useful given your colored lines)
#         lab_reg_by <- df_sub %>%
#             dplyr::mutate(celltype = droplevels(celltype)) %>%   # <-- key
#             split(.$celltype) %>%
#             lapply(function(d) {
#                 paste0(unique(d$celltype), ": ", fmt_lm_by(d, 'cm_pull_score'))
#             }) %>%
#             unlist() %>%
#             paste(collapse = "\n")

# 		lineage = 'CM'
#             # plot limits for placing annotation
#             x_lim <- quantile(df_sub$cm_pull_score, c(0.01, 0.99), na.rm = TRUE)
#             y_lim <- range(df_sub$score_GOi,    na.rm = TRUE)
#             y_top <- y_lim[2] + 0.18 * diff(y_lim)      ## CHANGED: extra headroom for text
#             y_txt1 <- y_lim[2] + 0.16 * diff(y_lim)     ## CHANGED: 1st line higher
#             y_txt2 <- y_lim[2] + 0.08 * diff(y_lim)     ## CHANGED: 2nd line higher
 
#             tmp <- ggplot(df_sub,  aes(x = cm_pull_score, y = score_GOi, color = celltype)) +
#                 geom_point(size = dot_size, alpha = 0.8) +
#                 geom_smooth(method = "lm", se = FALSE) +  #, aes(weight = if (use_wt) 1/(score_GPi_sem^2) else NULL)  use_wt = FALSE here
#                 scale_color_manual(values = celltype_colors, breaks = names(celltype_colors), labels=c('CP','CM','EMT')) +
#                 theme_bw() +
#                 labs(#x = paste0( length(cm_pull_genes), "-gene CM-pull over ",   length(cf_pull_genes), "-gene CF-pull score (cell-level)" ),
#                     x = paste0( length(cm_pull_genes), "-gene CM-pull score"),
# 					          y = paste0(n_id_present, "-gene ", id, " score")) +
#                 ggtitle(paste0("CTS_", CTS_name, '\n',go_names[id])) +
#                 coord_cartesian(ylim = c(y_lim[1], y_top), clip = "off") +
#                 annotate("text", 
#                         x = x_lim[1], y = y_txt1,
#                         label = lab_reg, hjust = 0, vjust = 1, size = 2.5) +
#                 annotate("text",
#                         x = x_lim[1], y = y_txt2,
#                         label = lab_reg_by, hjust = 0, vjust = 1, size = 2.5) +
#                 geom_abline(intercept = 0, beta = 1, linetype = "dashed", color = "black")
#             if(db=='Elorbany')  {
#                cat('pls modify 24.4.2_overlap_with_Takeuchi2026_v4_wdual_pull.R, which calls pseudobulk_regression_scatterplot()')
#                p_tmp2 =  ggplot() + theme_void() + labs(title = paste0('CTS_',CTS_name, ": ", lineage, '-bias -- No data'))
#             } else {
#               p_tmp2 = ggExtra::ggMarginal(
#                 tmp,
#                 type = "density",
#                 margins = "both",          # top + right
#                 groupColour = TRUE,        # color densities by celltype
#                 groupFill = FALSE,         # lines only (set TRUE if you want filled)
#                 size = 6                   # thickness of marginal panels
#                 )
#             }

# 			} else p_tmp2 =  ggplot() + theme_void() + labs(title = paste0('CTS_',CTS_name, ": ", lineage, '-bias -- No data'))
	
#   return(list(p_cf = p_tmp1, p_cm = p_tmp2))
# }


# ##########################################################################################################################################
# ## ENH3 vs WT within louvain = (C4) at pseudobulk level  ############
# # use RAW counts to build pseudobulk counts
# selected_genes = which(rowData(sce)$n_cells_by_counts >10 & rowData(sce)$feature_types == 'Gene Expression' &
#     #rowData(sce)$pct_dropout_by_counts <0.95 &  # In how many cells is this gene completely absent? -- too sparse
#     #rowData(sce)$pct_dropout_by_counts >0.05 &   # too uniform
#     rowData(sce)$total_counts >100
#     )  
#     #rowData(sce)$highly_variable     )
# counts <- assay(sce, "counts")[selected_genes, ]
# dim(counts)
# # [1]15303  1833
     
# meta <- as.data.frame(colData(sce))

# # highlight genes of interest
# genes_highlight <- c("TBX5" ,  "MYOCD" , "GATA6" , "RBM24", "SFRP5" , "HMGA2" , "CBFB"   ,"PRDM6" 
#   , "ANXA2",  "TWIST1", "ISL1" , "GATA6", "RBFOX2", "ETV2"
#   ,"PECAM1" ,"TNNT2" ,  "COL3A1", "FN1"  , "NKX2-5"
# )  #  "COL3A1", "VIM",  "MEF2C", "COL1A2", "GATA4",  'TBX1',  "NKX2-5", "FN1", "TNNT2",


pseudobulk_volcano_plot = function(sce, counts, meta, clustering_column, genes_highlight,
        n_background = 1000 , Rep_column = c('T_Reps', 'B_Reps'), reference_name = 'WT'){
  genotypes = setdiff(unique(colData(sce)$Sample), reference_name)
  res_list_all = list()
  for(genotype in genotypes){
    g_list = list()
    res_list = list()

    for(Cluster_of_interest in unique(colData(sce)[[clustering_column]]) ){
    #Cluster_of_interest = '4'

        # keep only cluster 4 (C3)
        cells_use <- as.character(meta[[clustering_column]]) == as.character(Cluster_of_interest)
        counts_sub <- counts[, cells_use]
        meta_sub <- meta[cells_use, ]

        # create pseudobulk sample ID
        # meta_sub$pb_id <- paste(meta_sub$Sample, meta_sub$B_Reps, sep = "_")
        meta_sub$pb_id <- paste(meta_sub$Sample, meta_sub[[Rep_column]], sep = "_")
        # aggregate counts
        pb_counts <- rowsum(t(counts_sub), group = meta_sub$pb_id)
        pb_counts <- t(pb_counts)   # genes x samples
        if(Rep_column == 'T_Reps')  pb_meta <- meta_sub %>% distinct(pb_id, Sample, T_Reps) else {
		    pb_meta <- meta_sub %>% distinct(pb_id, Sample, B_Reps)      
		}

        # match order
        pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$pb_id), ]
        # keep only WT vs ENH3
        keep <- pb_meta$Sample %in% c(reference_name, genotype)
		# --- NEW: ensure both groups present ---
		if (length(unique(pb_meta$Sample)) < 2) {
		  message("Skipping cluster ", Cluster_of_interest, " (only one condition present)")
		  next
		}
		if (any(table(pb_meta$Sample) < 2)) {
		  message("Skipping cluster ", Cluster_of_interest, " (not enough replicates)")
		  next
		}

        pb_counts <- pb_counts[, keep]
        pb_meta <- pb_meta[keep, ]
        pb_meta$Sample <- factor(pb_meta$Sample, levels = c(reference_name, genotype))

        # differential expression (edgeR)
        library(edgeR)

        dge <- DGEList(counts = pb_counts)
        dge <- calcNormFactors(dge)

        design <- model.matrix(~ Sample, data = pb_meta)
        dge <- estimateDisp(dge, design)
        fit <- glmQLFit(dge, design)
        qlf <- glmQLFTest(fit, coef = 2)  # ENH3 vs WT

        res <- topTags(qlf, n = Inf)$table
        res$gene <- rownames(res)
        table(res$logCPM>3)
    # FALSE  TRUE 
    #  3077 12982 
  
		# prepare for volcano plot
        res$logP <- -log10(res$PValue)
        res$signif <- res$FDR < 0.05 & abs(res$logFC) > 1 & res$logCPM > 3 & !is.na(res$logFC)

        genes_to_label <- unique(c(genes_highlight, res$gene[which(res$signif)]) )
        res$label <- ifelse(res$gene %in% genes_to_label, res$gene, NA)

        res$direction <- "ns"
        res$direction[(res$signif | res$gene %in% genes_highlight) & res$logFC < 0] <- "down"
        res$direction[(res$signif | res$gene %in% genes_highlight) & res$logFC > 0] <- "up"
        
    #    print(res[which(res$gene=='TBX5'),])
    #           logFC   logCPM        F    PValue       FDR   gene     logP
    # COL3A1 2.295694 4.871047 5.833796 0.0157223 0.7086984 COL3A1 1.803484
    #             padj signif  label direction
    # COL3A1 0.7086984  FALSE COL3A1        up

        # --- NEW: build reduced plotting dataset ---
		res_base <- subset(res, logCPM > 3)

		# background (non-significant & not highlighted)
		bg_df <- res_base %>% dplyr::filter(direction == "ns", !gene %in% genes_to_label)
		bg_n <- min(n_background, nrow(bg_df))
		# sample background
		set.seed(1)
		bg_df <- dplyr::slice_sample(bg_df, n = bg_n)

		# foreground (keep all signal)
		fg_df <- res_base %>%
		  dplyr::filter(direction != "ns" | gene %in% genes_highlight)

		# combine
		res_plot <- dplyr::bind_rows(bg_df, fg_df) %>%
		  dplyr::distinct(gene, .keep_all = TRUE)

		g_list[[Cluster_of_interest]] = ggplot(res_plot, aes(x = logFC, y = logP)) +
			#geom_point(color = "grey70", alpha = 0.6, size = 1) +
			# sampled background only
			geom_point(data = subset(res_plot, direction == "ns" & !gene %in% genes_highlight),
				color = "grey80", alpha = 0.6, size = 0.8) +
			geom_point(data = subset(res_plot, direction == "down"), color = "blue", size = 1.5) +
			geom_point(data = subset(res_plot, direction == "up"), color = "red", size = 1.5) +
			geom_point(data = subset(res_plot, gene %in% genes_highlight), color = "black", size = 1.5) +
			# labels for highlighted genes: always keep
			ggrepel::geom_text_repel(
			  data = subset(res_plot, gene %in% genes_highlight),
			  aes(label = gene),
			  size = 3,
			  max.overlaps = Inf
			) +
			# labels for other selected genes
			ggrepel::geom_text_repel(
			  data = subset(res_plot, !is.na(label) & !gene %in% genes_highlight),
			  aes(label = label),
			  size = 3
			) +
			geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
			geom_hline(yintercept = -log10(0.05), linetype = "dashed") +

			theme_classic() +
			labs(
				title = paste0(genotype, " vs ",reference_name," (C", Cluster_of_interest, " pseudobulk), colored: FDR<0.05, logCPM>3, |logFC|>1"),
				x = paste0("log2 Fold Change (", genotype, "-/- vs ",reference_name,")"),
				y = "-log10(p-value)"
			)  
    # g_list[[Cluster_of_interest]] 
    res_list[[Cluster_of_interest]] = res
    }
    res_list_all[[genotype]] = res_list
    pdf(paste0("volcano_plot_",genotype, "_vs_",reference_name,"_", clustering_column, "_pseudobulk_", Rep_column, ".pdf" ))
    for(i in 1:length(g_list)){print(g_list[[i]])}
    dev.off()
   }
  
  return(res_list_all)
}
