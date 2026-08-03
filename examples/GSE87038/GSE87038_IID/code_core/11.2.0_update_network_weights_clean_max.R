library("SingleCellExperiment")
library(Seurat)
library(dplyr)
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)
library(stringr)

########## BEGINNING OF USER INPUT ##########

wd <- "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/GSE87038_IID/"
outdir <- paste0(wd, "results_core/PPI_weight/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

db <- "GSE87038"

specificity_methods <- c("combined") # Other methods: "ratio", "zscore", "diff"

isl1_cluster <- "HiGCTS_8" # cluster containing ISL1 gene

core_count <- 1 # number of cores used for parallel processing in steps 1 and 2. Use core_count = 1 if on Windows.

step1 <- TRUE # calculate gene correlations and specificity
step2 <- TRUE # update network edge weights
step3 <- TRUE # graph comparing specificity methods for all clusters
step4 <- TRUE # graph comparing specificity methods for isl1_cluster

celltype_specific_weight_version <- "10"
BioTIP_version <- "06232025"

source(paste0("https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v", celltype_specific_weight_version, ".R"))
source(paste0("https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_", BioTIP_version, ".R"))

load(paste0(wd, "../data/sce_E8.25_uncorrected.RData"))
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

graph_list <- readRDS(file.path(wd, "results_core", paste0(db, "_IID_graph_perState_notsimplified.rds")))

graph_list <- graph_list[names(graph_list) != "CTS_13"]

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
range(E(graph_list[[1]])$weight) # [1] 0.01, 0.41

graph_list <- lapply(graph_list, simplify, edge.attr.comb = "max")
range(E(graph_list[[1]])$weight) # [1] 0.01 0.41
names(graph_list)
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_13"
# [26] "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"
# [31] "CTS_16.1"    "CTS_8"
edge_counts <- sapply(graph_list, ecount)
edge_counts
#       HiG_1       HiG_2       HiG_3       HiG_4       HiG_5       HiG_6
#      2254        3987        2959        2249        2056        4706
#     HiG_9      HiG_10      HiG_12      HiG_14      HiG_17      HiG_18
#      2758        3608        2430        2891        4503        3948
#    HiG_19       HiG_7      HiG_11      HiG_15      HiG_16      HiG_13
#      3399        2523        4158        3660        5008        3490
#     HiG_8    HiGCTS_7   HiGCTS_11   HiGCTS_15   HiGCTS_16 HiGCTS_16.1
#      2279           0           2          10           6           6
# HiGCTS_13    HiGCTS_8       CTS_7      CTS_11      CTS_15      CTS_16
#         3           1          10          18          50          28
#  CTS_16.1       CTS_8
#        57          65

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
  ## first, add a meta column to match the graph_list names
  # IID fix, to make sure sce$cluster matches the graph_list names as closely as possible, to avoid issues with calculating specificity and updating weights later on
  net_ids <- names(graph_list)
  focal_ids <- sub("^[^_]*_", "", net_ids)
  focal_ids <- sub("\\.\\d+$", "", focal_ids)
  focal_ids <- unique(focal_ids)

  cd <- as.data.frame(colData(sce))
  scores <- sapply(names(cd), function(col) {
    vals <- unique(as.character(cd[[col]]))
    sum(focal_ids %in% vals)
  })
  best_col <- names(scores)[which.max(scores)]

  colData(sce)$cluster <- as.character(cd[[best_col]])
  cat("Using cluster column:", best_col, "\n")


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
  # [19] "HiG_8"     "HiGCTS_15" "CTS_7"     "CTS_11"    "CTS_15"    "CTS_16"
  # [25] "CTS_16.1"  "CTS_8"

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
  edge_attr_names(g) # "weight"         "norm_PPI_score" "corexp_sign"    "coexp_focal"
}

###################################################
## compare weights methods
## check new weights for ISL1 in "HiGCTS_8"
#####################################################
library(ggplot2)
library(hexbin)

if (step3) {
  pdf(file = "compare_specificity_method_edgeweights.pdf")

  for (net_name in names(network_specificity_list)) {
    plot_data <- NULL

    for (s in specificity_methods) {
      weighted_graph_list <- readRDS(
        paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds")
      )

      g <- weighted_graph_list[[net_name]]

      # ---- SAFETY CHECKS ----
      if (is.null(g) || !igraph::is_igraph(g) || igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
        message("⚠️ Skipping ", net_name, " (empty/invalid graph) for method ", s)
        next
      }
      vn <- igraph::V(g)$name
      if (is.null(vn)) {
        message("⚠️ Skipping ", net_name, " (missing vertex names) for method ", s)
        next
      }
      w <- igraph::E(g)$weight
      if (is.null(w) || length(w) == 0) {
        message("⚠️ Skipping ", net_name, " (missing edge weights) for method ", s)
        next
      }
      # make log10 safe
      w[w <= 0] <- .Machine$double.xmin
      log_w <- log10(w)
      # -----------------------

      # mark edges incident to Isl1 (if present)
      isl1_vertex <- which(tolower(vn) == "isl1")
      is_isl1_edge <- rep(FALSE, igraph::ecount(g))
      if (length(isl1_vertex) > 0) {
        isl1_edges <- igraph::incident(g, isl1_vertex, mode = "all")
        is_isl1_edge[isl1_edges] <- TRUE
      }

      temp_data <- data.frame(
        net_name = net_name,
        log_weight = log_w,
        method = s,
        is_isl1 = is_isl1_edge
      )

      plot_data <- rbind(plot_data, temp_data)
    }

    if (!is.null(plot_data) && nrow(plot_data) > 0) {
      p2 <- ggplot(plot_data, aes(x = log_weight)) +
        geom_histogram(bins = 80) +
        geom_vline(
          data = subset(plot_data, is_isl1 == TRUE),
          aes(xintercept = log_weight),
          color = "red",
          alpha = 0.15
        ) +
        facet_wrap(~method, scales = "free") +
        theme_minimal() +
        labs(
          title = paste("Edge weight distribution -", net_name),
          x = "log10(E(g)$weight)",
          y = "Edge count"
        )

      print(p2)
    } else {
      message("⚠️ Skipping plot for ", net_name, " — no data collected.")
    }
  }

  dev.off()
}


if (step4) {
  net_name <- isl1_cluster
  plot_data <- NULL

  for (s in specificity_methods) {
    weighted_graph_list <- readRDS(
      paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds")
    )

    g <- weighted_graph_list[[net_name]]

    # ---- SAFETY CHECKS ----
    if (is.null(g) || !igraph::is_igraph(g) || igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
      message("⚠️ Skipping ", net_name, " (empty/invalid) for method ", s)
      next
    }
    vn <- igraph::V(g)$name
    if (is.null(vn)) {
      message("⚠️ Skipping ", net_name, " (missing vertex names) for method ", s)
      next
    }
    w <- igraph::E(g)$weight
    if (is.null(w) || length(w) == 0) {
      message("⚠️ Skipping ", net_name, " (missing edge weights) for method ", s)
      next
    }
    w[w <= 0] <- .Machine$double.xmin
    log_w <- log10(w)
    # -----------------------

    isl1_vertex <- which(tolower(vn) == "isl1")
    if (length(isl1_vertex) == 0) {
      message("No ISL1 found in ", isl1_cluster, " for method ", s)
      next
    }

    is_isl1_edge <- rep(FALSE, igraph::ecount(g))
    isl1_edges <- igraph::incident(g, isl1_vertex, mode = "all")
    is_isl1_edge[isl1_edges] <- TRUE

    temp_data <- data.frame(
      net_name = net_name,
      log_weight = log_w,
      method = s,
      is_isl1 = is_isl1_edge
    )

    plot_data <- rbind(plot_data, temp_data)
  }

  if (!is.null(plot_data) && any(plot_data$is_isl1)) {
    temp <- subset(plot_data, is_isl1 == TRUE)

    pdf(file = paste0("compare_specificity_method_", isl1_cluster, "_isl1_edgeweights.pdf"))

    print(
      ggplot(temp, aes(x = log_weight)) +
        geom_histogram(bins = 60) +
        facet_wrap(~method, scales = "free") +
        theme_minimal() +
        labs(
          title = paste("ISL1-incident edge weights -", net_name),
          x = "log10(E(g)$weight)",
          y = "Edge count"
        )
    )

    dev.off()
  } else {
    message("⚠️ Step4: no ISL1 edges collected for ", isl1_cluster)
  }
  # ⚠️ Skipping HiGCTS_8 (empty/invalid) for method combined
  # ⚠️ Step4: no ISL1 edges collected for HiGCTS_8
}
