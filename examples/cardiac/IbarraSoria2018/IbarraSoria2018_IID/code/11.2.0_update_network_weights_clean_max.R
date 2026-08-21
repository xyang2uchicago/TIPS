library("SingleCellExperiment")
library(Seurat)
library(dplyr)
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)
library(stringr)

########## BEGINNING OF USER INPUT ##########

source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID/")
outdir <- paste0(wd, "results/PPI_weight/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

db <- "IbarraSoria2018"

celltype_specific_weight_version <- "10"
BioTIP_version <- "06232025"

specificity_methods <- c("combined")

# cluster containing ISL1 gene (in this dataset)
isl1_cluster <- "HiGCTS_cardiac.a"

core_count <- 1 # use 1 on Windows / low-memory machines

step1 <- TRUE # calculate gene correlations and specificity
step2 <- TRUE # update network edge weights
step3 <- TRUE # compare edgeweight distributions across networks
step4 <- TRUE # ISL1-network focused distribution

source(here::here("R", "celltype_specific_weight_v10.R"))
source(paste0("https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_", BioTIP_version, ".R"))

load(file.path(wd, "../data", "sce_16subtype.RData"))

########## END OF USER INPUT ##########

# Normalize gene names in SCE to uppercase (important for IID vs expression matching)
rownames(sce) <- toupper(rownames(sce))

sce
# class: SingleCellExperiment
# dim: 4000 11039

colnames(colData(sce)) # label celltype subcelltype

unique(colData(sce)$celltype)
# [1] extraembryonicMesoderm endothelial            blood
# [4] mesodermProgenitors    presomiticMesoderm     somiticMesoderm
# [7] mixedMesoderm          pharyngealMesoderm     cardiac
# 9 Levels: blood cardiac endothelial ... somiticMesoderm

unique(colData(sce)$subcelltype)
#  [1] extraembryonicMesoderm endothelial.a          endothelial.b
#  [4] endothelial.c          endothelial.d          blood
#  [7] mesodermProgenitors    presomiticMesoderm.b   presomiticMesoderm.a
# [10] somiticMesoderm        mixedMesoderm.a        pharyngealMesoderm
# [13] mixedMesoderm.b        cardiac.a              cardiac.b
# [16] cardiac.c
# 16 Levels: blood cardiac.a cardiac.b cardiac.c endothelial.a ... somiticMesoderm

# Ensure logcounts exists
if (!"logcounts" %in% assayNames(sce)) {
  sce <- scater::logNormCounts(sce)
}
assayName <- "logcounts"

# -----------------------------
# 1) Read IID graphs
# -----------------------------
graph_list <- readRDS(file.path(wd, "results", paste0(db, "_IID_graph_perState_notsimplified.rds")))

graph_list <- lapply(graph_list, function(g) {
  V(g)$name <- toupper(V(g)$name)
  g
})

# Simplify edges using max (avoid doubling weights when multi-edges exist)
graph_list <- lapply(graph_list, simplify, edge.attr.comb = "max")

split_by_prefix <- function(graph_list, prefix) {
  keep <- startsWith(names(graph_list), paste0(prefix, "_"))
  gl <- graph_list[keep]
  # strip prefix so names match subcelltype
  names(gl) <- sub(paste0("^", prefix, "_"), "", names(gl))
  names(gl) <- sub("\\.\\d+$", "", names(gl))
  gl
}

graph_families <- list(
  HiG = split_by_prefix(graph_list, "HiG"),
  HiGCTS = split_by_prefix(graph_list, "HiGCTS"),
  CTS = split_by_prefix(graph_list, "CTS")
)

subtype_order <- unique(as.character(colData(sce)$subcelltype))
colData(sce)$subcelltype <- factor(as.character(colData(sce)$subcelltype), levels = subtype_order)


isl1_subtype <- sub("^[^_]*_", "", isl1_cluster) # -> "cardiac.a"
cat("ISL1 focal subtype:", isl1_subtype, "\n")

lapply(graph_families, length) # HiG: 16 HiGCTS: 2 CTS: 2

names(graph_list)
#  [1] "HiG_blood"                  "HiG_cardiac.b"
#  [3] "HiG_cardiac.c"              "HiG_endothelial.a"
#  [5] "HiG_endothelial.c"          "HiG_endothelial.d"
#  [7] "HiG_extraembryonicMesoderm" "HiG_mesodermProgenitors"
#  [9] "HiG_mixedMesoderm.a"        "HiG_mixedMesoderm.b"
# [11] "HiG_pharyngealMesoderm"     "HiG_presomiticMesoderm.a"
# [13] "HiG_presomiticMesoderm.b"   "HiG_somiticMesoderm"
# [15] "HiG_endothelial.b"          "HiG_cardiac.a"
# [17] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"
# [19] "CTS_endothelial.b"          "CTS_cardiac.a"
edge_counts <- sapply(graph_list, ecount)
edge_counts
#                  HiG_blood              HiG_cardiac.b
#                       3118                       3498
#              HiG_cardiac.c          HiG_endothelial.a
#                       3764                       4382
#          HiG_endothelial.c          HiG_endothelial.d
#                       3966                       3247
# HiG_extraembryonicMesoderm    HiG_mesodermProgenitors
#                       2654                       2187
#        HiG_mixedMesoderm.a        HiG_mixedMesoderm.b
#                       2134                       2750
#     HiG_pharyngealMesoderm   HiG_presomiticMesoderm.a
#                       2717                       2129
#   HiG_presomiticMesoderm.b        HiG_somiticMesoderm
#                       3462                       2162
#          HiG_endothelial.b              HiG_cardiac.a
#                       3156                       3072
#       HiGCTS_endothelial.b           HiGCTS_cardiac.a
#                          1                          5
#          CTS_endothelial.b              CTS_cardiac.a
#                          5                         17

names(graph_families$HiG)
#  [1] "blood"                  "cardiac.b"              "cardiac.c"
#  [4] "endothelial.a"          "endothelial.c"          "endothelial.d"
#  [7] "extraembryonicMesoderm" "mesodermProgenitors"    "mixedMesoderm.a"
# [10] "mixedMesoderm.b"        "pharyngealMesoderm"     "presomiticMesoderm.a"
# [13] "presomiticMesoderm.b"   "somiticMesoderm"        "endothelial.b"
# [16] "cardiac.a"
names(graph_families$HiGCTS)
# [1] "endothelial.b" "cardiac.a"
names(graph_families$CTS)
# [1] "endothelial.b" "cardiac.a"

for (fam in names(graph_families)) {
  gl <- graph_families[[fam]]
  graph_names <- names(gl)
  sce_levels <- levels(colData(sce)$subcelltype)

  missing_in_sce <- setdiff(graph_names, sce_levels)
  missing_in_graphs <- setdiff(sce_levels, graph_names)

  cat("\n--- Name matching check:", fam, "---\n")
  cat("Graphs not in sce subcelltype levels:", length(missing_in_sce), "\n")
  if (length(missing_in_sce) > 0) print(missing_in_sce)

  cat("sce subcelltype levels missing graphs:", length(missing_in_graphs), "\n")
  if (length(missing_in_graphs) > 0) print(missing_in_graphs)
  cat("-----------------------------------\n")
}

# --- Name matching check: HiG ---
# Graphs not in sce subcelltype levels: 0
# sce subcelltype levels missing graphs: 0
# -----------------------------------

# --- Name matching check: HiGCTS ---
# Graphs not in sce subcelltype levels: 0
# sce subcelltype levels missing graphs: 14
#  [1] "extraembryonicMesoderm" "endothelial.a"          "endothelial.c"
#  [4] "endothelial.d"          "blood"                  "mesodermProgenitors"
#  [7] "presomiticMesoderm.b"   "presomiticMesoderm.a"   "somiticMesoderm"
# [10] "mixedMesoderm.a"        "pharyngealMesoderm"     "mixedMesoderm.b"
# [13] "cardiac.b"              "cardiac.c"
# -----------------------------------

# --- Name matching check: CTS ---
# Graphs not in sce subcelltype levels: 0
# sce subcelltype levels missing graphs: 14
#  [1] "extraembryonicMesoderm" "endothelial.a"          "endothelial.c"
#  [4] "endothelial.d"          "blood"                  "mesodermProgenitors"
#  [7] "presomiticMesoderm.b"   "presomiticMesoderm.a"   "somiticMesoderm"
# [10] "mixedMesoderm.a"        "pharyngealMesoderm"     "mixedMesoderm.b"
# [13] "cardiac.b"              "cardiac.c"
# -----------------------------------


############################################################
library(data.table)
library(ggplot2)
library(hexbin)

subtype_order <- unique(colData(sce)$subcelltype)
colData(sce)$subcelltype <- factor(colData(sce)$subcelltype,
  levels = subtype_order
)

# Helper: "HiGCTS_cardiac.a" -> "cardiac.a"
net_to_subtype <- function(net_name) sub("^[^_]*_", "", net_name)

############################################################
# STEP 1
############################################################
if (step1) {
  ref_graph_list <- graph_families$HiG
  stopifnot(length(ref_graph_list) > 0)

  network_specificity_ref <- calculate_network_specificity(
    sce,
    ref_graph_list,
    assayName,
    celltype_col = "subcelltype",
    method = "pearson",
    cores = core_count,
    shrink = TRUE,
    verbose = FALSE
  )

  # Expand specificity to ALL prefixed IID graphs
  network_specificity_list <- vector("list", length(graph_list))
  names(network_specificity_list) <- names(graph_list)

  for (net in names(graph_list)) {
    subt <- net_to_subtype(net)
    network_specificity_list[[net]] <- network_specificity_ref[[subt]]
  }

  saveRDS(network_specificity_list,
    file = "network_specificity_list.rds"
  )
}

############################################################
# STEP 2
############################################################

if (step2) {
  network_specificity_list <- readRDS("network_specificity_list.rds")

  for (net in names(network_specificity_list)) {
    spec_data <- network_specificity_list[[net]]
    corexp_named <- spec_data$corexp_sign

    stopifnot(
      is.matrix(corexp_named),
      all(dim(corexp_named) == dim(spec_data$coexp_focal))
    )

    network_specificity_list[[net]]$corexp_sign <- corexp_named
  }

  for (s in specificity_methods) {
    weighted_graph_list <- update_network_weights(
      graph_list,
      network_specificity_list,
      specificity_method = s,
      verbose = FALSE,
      cores = core_count
    )

    saveRDS(weighted_graph_list,
      file = paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds")
    )
  }
}

############################################################
# STEP 3
############################################################
if (step3) {
  pdf(file = "2018_compare_specificity_method_vs_PPIscores.pdf")

  for (net_name in names(graph_list)) {
    plot_data <- NULL

    for (s in specificity_methods) {
      weighted_graph_list <- readRDS(
        paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds")
      )

      g <- weighted_graph_list[[net_name]]
      if (!is_igraph(g) || ecount(g) == 0) next

      isl1_vertex <- which(tolower(V(g)$name) == "isl1")
      is_isl1_edge <- rep(FALSE, ecount(g))
      if (length(isl1_vertex) > 0) {
        isl1_edges <- incident(g, isl1_vertex, mode = "all")
        is_isl1_edge[isl1_edges] <- TRUE
      }

      temp_data <- data.frame(
        net_name = net_name,
        norm_PPI_score = E(g)$norm_PPI_score,
        log_weight = log10(pmax(E(g)$weight, .Machine$double.xmin)),
        method = s,
        is_isl1 = is_isl1_edge
      )

      plot_data <- rbind(plot_data, temp_data)
    }

    p2 <- ggplot(
      plot_data,
      aes(
        x = norm_PPI_score,
        y = log_weight
      )
    ) +
      stat_binhex(bins = 50, alpha = 0.7) +
      geom_point(
        data = subset(plot_data, is_isl1 == TRUE),
        aes(
          x = norm_PPI_score,
          y = log_weight
        ),
        color = "red",
        size = 1.5
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
  }

  dev.off()
}

############################################################
# STEP 4
############################################################

if (step4) {
  net_name <- isl1_cluster
  plot_data <- NULL

  for (s in specificity_methods) {
    weighted_graph_list <- readRDS(
      paste0(db, "_IID_graph_perState_simplified_", s, "weighted.rds")
    )

    g <- weighted_graph_list[[net_name]]
    if (!is_igraph(g) || ecount(g) == 0) next

    isl1_vertex <- which(tolower(V(g)$name) == "isl1")
    if (length(isl1_vertex) == 0) next

    is_isl1_edge <- rep(FALSE, ecount(g))
    isl1_edges <- incident(g, isl1_vertex, mode = "all")
    is_isl1_edge[isl1_edges] <- TRUE

    temp_data <- data.frame(
      net_name = net_name,
      norm_PPI_score = E(g)$norm_PPI_score,
      log_weight = log10(pmax(E(g)$weight, .Machine$double.xmin)),
      method = s,
      is_isl1 = is_isl1_edge
    )

    plot_data <- rbind(plot_data, temp_data)
  }

  if (!is.null(plot_data) && any(plot_data$is_isl1)) {
    pdf(file = paste0(
      "2018_compare_specificity_method_",
      isl1_cluster,
      "_vs_PPIscores.pdf"
    ))

    p3 <- ggplot(
      subset(plot_data, is_isl1 == TRUE),
      aes(
        x = norm_PPI_score,
        y = log_weight,
        color = as.factor(norm_PPI_score)
      )
    ) +
      geom_point() +
      facet_wrap(~method, scales = "free") +
      theme_minimal() +
      labs(
        title = paste("Isl1 linkages -", net_name),
        x = "E(g)$norm_PPI_score",
        y = "log10(E(g)$weight)"
      ) +
      theme(legend.position = "bottom")

    print(p3)
    dev.off()
  } else {
    message("⚠️ Step4: no ISL1 edges collected for ", isl1_cluster)
  }
}
