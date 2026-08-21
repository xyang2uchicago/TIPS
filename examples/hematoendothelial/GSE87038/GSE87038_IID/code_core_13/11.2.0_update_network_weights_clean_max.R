library("SingleCellExperiment")
library(Seurat)
library(dplyr)
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)
library(stringr)

########## BEGINNING OF USER INPUT ##########

code_dir <- here::here("examples", "hematoendothelial", "GSE87038", "GSE87038_IID", "code_core_13")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
outdir <- ppi_path
setwd(outdir)

source(paste0("https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v", celltype_specific_weight_version, ".R"))
source(paste0("https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_", BioTIP_version, ".R"))

load(paste0(data_dir, "sce_E8.25_uncorrected.RData"))
rownames(sce) <- toupper(rownames(sce)) # !!!!!!!


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
assayName <- "logcounts"

graph_list <- readRDS(file.path(wd, "results_core_13", paste0(db, "_IID_graph_perState_notsimplified.rds")))


graph_list <- lapply(graph_list, function(g) {
  V(g)$name <- toupper(V(g)$name)
  g
}) # makes sure the sce and the graph list work together without case sensitivity issues.


N0 <- sapply(graph_list, vcount)

graphs_with_duplicates <- sapply(graph_list, function(g) any(duplicated(V(g)$name)))
stopifnot(!any(graphs_with_duplicates))


N <- sapply(graph_list, vcount)
(N0 - N)[which(N0 - N > 0)]
# named numeric(0)
range(E(graph_list[[1]])$weight) # [1] 0 0.05

graph_list <- lapply(graph_list, simplify, edge.attr.comb = "max")
range(E(graph_list[[1]])$weight) # [1] 0 0.05
names(graph_list)
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_13"
# [26] "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"
# [31] "CTS_16.1"    "CTS_13"      "CTS_8"
# (33 networks; HiGCTS_7/11/16/16.1/13/8 and CTS_7/11/16 have <10 edges and are excluded by NSL step)
edge_counts <- sapply(graph_list, ecount)
edge_counts
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6
#      2921        4256        3585        2810        2566        4089
#     HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#      3428        4705        2968        2902        3970        3840
#    HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#      3234        2429        4661        3905        4610        3661
#     HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#      3095           0           1          12           1           2
# HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#         1           0           3           7          17           6
#  CTS_16.1      CTS_13       CTS_8
#        29          19          23

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
# named integer(0)

if (step1) {
  # Hardcode "label" (1-19 fine-grained) — auto-detection can pick the coarse "cluster"
  # column (6 values) instead, which causes calculate_network_specificity to skip CTS_13.
  colData(sce)$cluster <- as.character(colData(sce)$label)
  cat("Using cluster column: label (hardcoded)\n")

  network_specificity_list <- calculate_network_specificity(sce,
    graph_list,
    assayName,
    celltype_col = "cluster",
    method = "pearson",
    cores = core_count,
    shrink = TRUE,
    verbose = FALSE
  )
  saveRDS(network_specificity_list, "IID_network_specificity_list.rds")

  names(network_specificity_list)
  #  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"
  #  [7] "HiG_9"     "HiG_10"    "HiG_12"    "HiG_14"    "HiG_17"    "HiG_18"
  # [13] "HiG_19"    "HiG_7"     "HiG_11"    "HiG_15"    "HiG_16"    "HiG_13"
  # [19] "HiG_8"     "HiGCTS_15" "CTS_15"    "CTS_16.1"  "CTS_13"    "CTS_8"
  # (24 networks: HiGCTS_7/11/16/16.1/13/8 and CTS_7/11/16 excluded by <10 edge threshold)

  names(network_specificity_list[[1]])
  # [1] "scores"      "genes"       "coexp_focal" "corexp_sign"

  names(network_specificity_list[[1]][["scores"]])
  # [1] "ratio"    "zscore"   "diff"     "combined"
}

if (step2) {
  network_specificity_list <- readRDS("IID_network_specificity_list.rds")

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
  # [1] "scores"      "genes"       "coexp_focal" "corexp_sign"

  table(network_specificity_list[[1]]$corexp_sign)
  # negative positive
  #  32800    41729

  for (s in specificity_methods) {
    weighted_graph_list <- update_network_weights(graph_list,
      network_specificity_list,
      specificity_method = s,
      verbose = FALSE,
      cores = 1
    )
    saveRDS(weighted_graph_list, file = paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds"))
  }

  # double check outputs
  g <- weighted_graph_list[[1]]
  vertex_attr_names(g) # "name"   "weight" "FDR"
  edge_attr_names(g)   # "weight"  "norm_PPI_score"  "corexp_sign"  "coexp_focal"
  length(weighted_graph_list) # 24
}

