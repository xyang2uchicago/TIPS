library("SingleCellExperiment")
library(Seurat)
library(dplyr)
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)
library(stringr)

########## BEGINNING OF USER INPUT ##########

wd = "C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/"
setwd(paste0(wd, "results/PPI_weight/"))

db <- "GSE116893"

specificity_methods <- c("combined") # Other methods: "ratio", "zscore", "diff"

cluster_labels <- "leiden_0.6"


core_count <- 1 # number of cores used for parallel processing in steps 1 and 2. Use core_count = 1 if on Windows.

step1 <- FALSE # calculate gene correlations and specificity
step2 <- TRUE # update network edge weights

celltype_specific_weight_version <- '10'
BioTIP_version <- '06232025'

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))
source(paste0('https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_', BioTIP_version, '.R'))

sce <- readRDS(paste0(wd, "data/sce.rds"))

########## END OF USER INPUT ##########

sce
colnames(colData(sce))
#  [1] "orig.ident"             "nCount_RNA"             "nFeature_RNA"          
#  [4] "Stage_Code"             "Tissue"                 "Patient_No"            
#  [7] "Risk_Category"          "First_Avail_TP"         "MYCN_Status"           
# [10] "ALK_Status"             "TP53_Status"            "Response"              
# [13] "Vital_Status"           "Age_at_IDX_in_months"   "Gender"                
# [16] "Treatment"              "First_Avail_Time_Point" "sample_name"           
# [19] "biospecimen_id"         "percent.mt"             "seurat_clusters"       
# [22] "sample_label_wo_prefix" "cell_type"              "S.Score"               
# [25] "G2M.Score"              "Phase"                  "malignancy"            
# [28] "cell_state"             "MES_Score"              "ADRN_Score"            
# [31] "MES_ADRN_diff"          "Event"                  "leiden_0.2"            
# [34] "leiden_0.4"             "leiden_0.6"             "leiden_0.8"            
# [37] "leiden_1"               "leiden_1.2"             "leiden_0.3"            
# [40] "leiden_0.5"             "leiden_0.7"             "RNA_snn_res.0.2"       
# [43] "ident"

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
# HiG_17 
#      1
range(E(graph_list[[1]])$weight) # [1] 0.200 0.999

graph_list <- lapply(graph_list, simplify, edge.attr.comb = "max") # by default is 'sum'
range(E(graph_list[[1]])$weight) # [1]  0.200 0.999

names(graph_list)
#  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"    
#  [7] "HiG_7"     "HiG_8"     "HiG_10"    "HiG_11"    "HiG_12"    "HiG_13"   
# [13] "HiG_14"    "HiG_16"    "HiG_17"    "HiG_15"    "HiG_9"     "HiGCTS_15"
# [19] "HiGCTS_16" "HiGCTS_9"  "CTS_15"    "CTS_16"    "CTS_9"

edge_counts <- sapply(graph_list, ecount)
edge_counts
#     HiG_1     HiG_2     HiG_3     HiG_4     HiG_5     HiG_6     HiG_7     HiG_8 
#     13101      4272      3520      4372      4623      9872      4363      3946 
#    HiG_10    HiG_11    HiG_12    HiG_13    HiG_14    HiG_16    HiG_17    HiG_15 
#      5037     10588     10459     16636     24134     10239     37844     23162 
#     HiG_9 HiGCTS_15 HiGCTS_16  HiGCTS_9    CTS_15    CTS_16     CTS_9 
#      2407       381       413         0       931       493       454

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
        verbose = TRUE
    )
    saveRDS(network_specificity_list, "network_specificity_list.rds")
    names(network_specificity_list)
    #  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"    
    #  [7] "HiG_7"     "HiG_8"     "HiG_10"    "HiG_11"    "HiG_12"    "HiG_13"   
    # [13] "HiG_14"    "HiG_16"    "HiG_17"    "HiG_15"    "HiG_9"     "HiGCTS_15"
    # [19] "HiGCTS_16" "CTS_15"    "CTS_16"    "CTS_9" 

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
    # negative positive 
    #    66468    71916

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