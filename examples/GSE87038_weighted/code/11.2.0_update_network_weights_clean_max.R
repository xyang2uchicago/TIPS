library("SingleCellExperiment")
library(Seurat)
library(dplyr)
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)
library(stringr)

wd = "/Users/felixyu/Documents/GSE87038_weighted/"
setwd(paste0(wd, "results/PPI_weight/"))

PPI_color_platte <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

load(paste0(wd, "data/sce_E8.25_uncorrected.RData"))

celltype_specific_weight_version <- '10'
BioTIP_version <- '06232025'

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_', BioTIP_version, '.R'))

sce
colnames(colData(sce))
#  [1] "cell"             "barcode"          "sample"           "pool"
#  [5] "stage"            "sequencing.batch" "theiler"          "doub.density"
#  [9] "doublet"          "cluster"          "cluster.sub"      "cluster.stage"
# [13] "cluster.theiler"  "stripped"         "celltype"         "colour"
# [17] "sizeFactor"       "batch"            "label"

# Calculate log-normalized counts
# FIXED
if (!"logcounts" %in% assayNames(sce)) {
    sce <- scater::logNormCounts(sce)
}
assayName = 'logcounts'

## Problem detection and clean 2:
## the doubled range of weights after using simplify() on your STRING graph
## clean 2: using the max score for each edge of STRING.db
#########################################################################
diagnose <- FALSE
if (diagnose) {
    graph_list0 <- graph_list

    range(E(graph_list0[[1]])$weight) # [1] 0.200 0.999
    any(which_multiple(graph_list0[[1]]))
    table(which_multiple(graph_list0[[1]]))
    E(graph_list[[1]])[which_multiple(graph_list0[[1]])]
    edge_attr_names(graph_list0[[1]])

    graph_list <- lapply(graph_list0, simplify) # !!!!!!!!!!!!!!!!!!!
    range(E(graph_list[[1]])$weight) # [1] 0.400 1.998

    graph_list <- lapply(graph_list0, simplify, edge.attr.comb = "max") # by default is 'sum'
    range(E(graph_list[[1]])$weight) # [1]  0.200 0.999
}

inputdir <- paste0(wd, "results/")
graph_list <- readRDS(file = paste0(inputdir, "GSE87038_STRING_graph_perState_notsimplified.rds"))
N0 <- sapply(graph_list, vcount)

########## clean 1) remove name-duplicated Vertex due to inconsistency in STRING.db ###########
## refer to 11.1.0_correct_vertex_duplication.R
correct_n_edges <- readRDS(paste0(inputdir, "correct_n_edges_HiG_STRING2.14.0.rds"))
for (g_name in unique(correct_n_edges$graph_id)) {
    # Subset the rows for this graph
    rows <- subset(correct_n_edges, graph_id == g_name)

    # Extract the vertex_index_to_remove as a character vector
    vertices_str <- rows$vertex_index_to_remove
    vertices_str <- vertices_str[!is.na(vertices_str)]

    # Split each string by "," and convert to numeric, then combine all
    vertices_to_remove <- unlist(
        lapply(vertices_str, function(s) as.numeric(strsplit(s, ",")[[1]]))
    )

    # Remove duplicates and sort decreasingly to avoid index shifting when deleting
    vertices_to_remove <- sort(unique(vertices_to_remove), decreasing = TRUE)
    cat("Removing from", g_name, ":", paste(vertices_to_remove, collapse = ", "), "\n")
    if (length(vertices_to_remove) > 0) {
        graph_list[[g_name]] <- delete_vertices(graph_list[[g_name]], vertices_to_remove)
    }
}

N <- sapply(graph_list, vcount)
(N0 - N)[which(N0 - N > 0)]
#  HiG_5 HiG_17 HiG_18 HiG_19
#      2      1      2      2
range(E(graph_list[[1]])$weight) # [1] 0.200 0.999

graph_list <- lapply(graph_list, simplify, edge.attr.comb = "max") # by default is 'sum'
range(E(graph_list[[1]])$weight) # [1]  0.200 0.999

names(graph_list)
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_13"
# [26] "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"
# [31] "CTS_16.1"    "CTS_13"      "CTS_8"
edge_counts <- sapply(graph_list, ecount)
edge_counts
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6
#        4978        9544        7120        5362        5808       10993
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#        6937        8655        5828        7004       11719        8543
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#        7827        7572        9301        8632       10581        8508
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#        5880          17          15          54          10          31
#   HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#           9          10          52          54         207          70
#    CTS_16.1      CTS_13       CTS_8
#         210         147         165

graphs_with_duplicates <- sapply(graph_list, function(g) {
    vertex_names <- V(g)$name
    if (is.null(vertex_names)) {
        # If no names, use vertex indices
        vertex_names <- V(g)
    }
    any(duplicated(vertex_names))
})

# See which graphs have duplicates
which(graphs_with_duplicates)

step1 <- FALSE
if (step1) {
    ## calculating co-expression-based specificity scores, runs 1 hour when cores=1, DO NOT repeat !
    ## first, add a meta column to match the graph_list names
    colData(sce)$cluster <- colData(sce)$label

    network_specificity_list <- calculate_network_specificity(sce,
        graph_list,
        assayName,
        celltype_col = "cluster",
        method = "pearson",
        cores = 4,
        shrink = TRUE,
        verbose = TRUE
    )
    saveRDS(network_specificity_list, "network_specificity_list.rds")
    names(network_specificity_list)
    #  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
    #  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
    # [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
    # [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
    # [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_8"
    # [26] "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"      "CTS_16.1"
    # [31] "CTS_13"      "CTS_8"

    names(network_specificity_list[[1]])
    # [1] "scores"      "genes"       "corexp_sign"

    names(network_specificity_list[[1]][["scores"]])
    # [1] "ratio"    "zscore"   "diff"     "combined"
}

step2 <- FALSE
if (step2) {
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
    # [1] "scores"      "genes"       "corexp_sign"

    table(network_specificity_list[[1]]$corexp_sign)
    # negative positive
    # 40106    51703

    for (s in c("combined", "ratio", "zscore", "diff")) {
        weighted_graph_list <- update_network_weights(graph_list,
            network_specificity_list,
            specificity_method = s,
            verbose = TRUE,
            cores = 1
        )
        saveRDS(weighted_graph_list, file = paste0("GSE87038_STRING_graph_perState_simplified_", s, "weighted.rds"))
    }

    # double check outputs
    g <- weighted_graph_list[[1]]
    vertex_attr_names(g) # "name"   "weight" "FDR"  
    edge_attr_names(g) # "weight"         "norm_PPI_score" "corexp_sign"    "coexp_target"
}

###################################################
## compare weights methods
## check new weights for ISL1 in "HiGCTS_8"
#####################################################
step3 <- TRUE
library(ggplot2)
library(hexbin)


if (step3) {
    pdf(file = paste0("compare_specificity_method_vs_PPIscores.pdf"))

    for (net_name in names(network_specificity_list))
    {
        # net_name = "HiGCTS_8"
        plot_data <- NULL
        par(mfrow = c(2, 2))

        for (s in c("combined", "ratio", "zscore", "diff")) {
            graph_list <- readRDS(file = paste0("GSE87038_STRING_graph_perState_simplified_", s, "weighted.rds"))
            g <- graph_list[[net_name]]
            # Safeguard to skip graphs with missing original weights
            if (is.null(E(g)$norm_PPI_score) || length(E(g)$norm_PPI_score) == 0) {
                cat("⚠️  Skipping", net_name, "- no 'norm_PPI_score' on edges.\n")
                next
            }

            # Check for ISL1 vertex and mark edges
            isl1_vertex <- which(tolower(V(g)$name) == "isl1")
            is_isl1_edge <- rep(FALSE, length(E(g)))
            if (length(isl1_vertex) > 0) {
                isl1_edges <- incident(g, isl1_vertex, mode = "all")
                is_isl1_edge[isl1_edges] <- TRUE
            }

            temp_data <- data.frame(
                net_name = net_name,
                norm_PPI_score = E(g)$norm_PPI_score,
                log_weight = log10(E(g)$weight),
                method = s,
                is_isl1 = is_isl1_edge
            )

            plot_data <- rbind(plot_data, temp_data)
        }
        if (!is.null(plot_data) && "is_isl1" %in% colnames(plot_data)) {
            p2 <- ggplot(plot_data, aes(x = norm_PPI_score, y = log_weight)) +
                stat_binhex(bins = 50, alpha = 0.7) +
                geom_point(
                    data = subset(plot_data, is_isl1 == TRUE),
                    aes(x = norm_PPI_score, y = log_weight),
                    color = "red", size = 1.5
                ) +
                scale_fill_viridis_c() +
                facet_wrap(~method, scales = "free") +
                theme_minimal() +
                labs(
                    title = paste("Edge Weight Density -", net_name),
                    x = "E(g)$norm_PPI_score",
                    y = "log10(E(g)$weight)",
                    fill = "Count"
                ) +
                theme(legend.position = "bottom")

            print(p2)
        } else {
            message("⚠️  Skipping plot for ", net_name, " — plot_data is empty or malformed.")
        }
    }
    dev.off()
}

step4 <- TRUE
if (step4) {
    net_name <- "HiGCTS_8"
    plot_data <- NULL

    for (s in c("combined", "ratio", "zscore", "diff")) {
        graph_list <- readRDS(file = paste0("GSE87038_STRING_graph_perState_simplified_", s, "weighted.rds"))
        g <- graph_list[[net_name]]

        graphs_with_duplicates <- sapply(graph_list, function(g) {
            vertex_names <- V(g)$name
            if (is.null(vertex_names)) {
                # If no names, use vertex indices
                vertex_names <- V(g)
            }
            any(duplicated(vertex_names))
        })

        which(graphs_with_duplicates)

        # Check for ISL1 vertex and mark edges
        isl1_vertex <- which(tolower(V(g)$name) == "isl1")
        is_isl1_edge <- rep(FALSE, length(E(g)))
        if (length(isl1_vertex) > 0) {
            isl1_edges <- incident(g, isl1_vertex, mode = "all")
            is_isl1_edge[isl1_edges] <- TRUE
        }

        temp_data <- data.frame(
            net_name = net_name,
            norm_PPI_score = E(g)$norm_PPI_score,
            log_weight = log10(E(g)$weight),
            method = s,
            is_isl1 = is_isl1_edge
        )

        plot_data <- rbind(plot_data, temp_data)
    }

    temp <- subset(plot_data, is_isl1 == TRUE)
    pdf(file = paste0("compare_specificity_method_HiGCTS_8_vs_PPIscores.pdf"))
    print(
        ggplot(temp, aes(x = norm_PPI_score, y = log_weight, color = as.factor(norm_PPI_score))) +
            geom_point() +
            scale_fill_viridis_c() +
            facet_wrap(~method, scales = "free") +
            theme_minimal() +
            labs(
                title = paste("Isl1 linkages -", net_name),
                x = "E(g)$norm_PPI_score",
                y = "log10(E(g)$weight)"
            ) +
            theme(legend.position = "bottom")
    )

    dev.off()
}
