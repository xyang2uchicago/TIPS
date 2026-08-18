library("SingleCellExperiment")
library(dplyr)
library(scuttle)

library(igraph)

########## BEGINNING OF USER INPUT ##########

source(here::here("examples", "config.R"))
wd = tips_path("examples", "hematoendothelial", "GSE87038", "GSE87038_STRING/")
setwd(paste0(wd, "results_core_13/PPI_weight/"))

db <- "GSE87038"

specificity_methods <- c("combined") # Other methods: "ratio", "zscore", "diff"

cluster_labels <- "label"

core_count <- 1 # number of cores used for parallel processing in steps 1 and 2. Use core_count = 1 if on Windows.

step1 <- TRUE # calculate gene correlations and specificity
step2 <- TRUE # update network edge weights


celltype_specific_weight_version <- '10'
BioTIP_version <- '06232025'

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_', BioTIP_version, '.R'))

load(paste0(wd, "../data/sce_E8.25_uncorrected.RData"))

########## END OF USER INPUT ##########

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
    assayName <- assayNames(sce)[1]
    x <- assay(sce, assayName)
    if (max(x) > 20) {
        sce <- scater::logNormCounts(sce, assay.type = assayName)
    } else {
        assayNames(sce)[1] <- 'logcounts'
    }
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

graph_list <- readRDS(file = paste0(paste0("../", db, "_STRING_graph_perState_notsimplified.rds")))
N0 <- sapply(graph_list, vcount)

########## clean 1) remove name-duplicated Vertex due to inconsistency in STRING.db ###########
## refer to 11.1.0_correct_vertex_duplication.R
correct_n_edges <- readRDS(paste0("../correct_n_edges_HiG_STRING2.14.0.rds"))
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
#  HiG_1  HiG_2  HiG_3  HiG_4  HiG_5  HiG_6  HiG_9 HiG_17 HiG_18 HiG_19 HiG_16
#      2      1      2      2      2      1      2      1      2      3      2
#  HiG_8
#      2
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
#        9420       15155       12319        9506        9657       15204
#       HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#       12042       15430       10466       11456       16115       12819
#      HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#       12444       10034       15491       13989       15111       13730
#       HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#       10272          20          18         110          12          40
#   HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#          16          14          52          54         207          70
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

if (step1) {
    ## first, add a meta column to match the graph_list names
    colData(sce)$cluster <- colData(sce)[[cluster_labels]]

    network_specificity_list <- calculate_network_specificity(sce,
        graph_list,
        assayName,
        celltype_col = "cluster",
        method = "pearson",
        cores = core_count,
        shrink = TRUE,
        verbose = FALSE
    )
    saveRDS(network_specificity_list, "network_specificity_list.rds")
    names(network_specificity_list)
    # All 33 networks computed (all 19 clusters present via cluster_labels="label")
    # HiGCTS_13 included: 16 edges at lfc_HiGCTS=0.5 (was 9 edges at lfc=0.6)

    names(network_specificity_list[[1]])
    # [1] "scores"      "genes"       "corexp_sign"

    names(network_specificity_list[[1]][["scores"]])
    # [1] "ratio"    "zscore"   "diff"     "combined"
}

if (step2) {
    network_specificity_list <- readRDS("network_specificity_list.rds")

    library(data.table)

    for (net in names(network_specificity_list)) {
        cat("Analyzing: ", net, "\n")
        spec_data <- network_specificity_list[[net]]
        corexp_named <- spec_data$corexp_sign

        stopifnot(
            is.matrix(corexp_named),
            all(dim(corexp_named) == dim(spec_data$coexp_focal))
        )

        network_specificity_list[[net]]$corexp_sign <- corexp_named
    }

    names(network_specificity_list[[1]])
    # [1] "scores"      "genes"       "corexp_sign"

    table(network_specificity_list[[1]]$corexp_sign)
    # network_specificity_list[[1]] is HiG_6 (first computed network)

    for (s in specificity_methods) {
        weighted_graph_list <- update_network_weights(graph_list,
            network_specificity_list,
            specificity_method = s,
            verbose = FALSE,
            cores = 1
        )
        saveRDS(weighted_graph_list, file = paste0(db, "_STRING_graph_perState_simplified_", s, "weighted.rds"))
    }

    # double check outputs
    g <- weighted_graph_list[[1]]
    vertex_attr_names(g) # "name"   "weight" "FDR"  
    edge_attr_names(g) # "weight"         "norm_PPI_score" "corexp_sign"    "coexp_focal"
}
