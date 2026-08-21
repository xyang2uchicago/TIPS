library("SingleCellExperiment")
library(Seurat)
library(dplyr)
library(scuttle)

## dependence to run BioTIP
library(igraph)
require(psych)
library(stringr)

########## BEGINNING OF USER INPUT ##########

code_dir <- here::here("examples", "hematoendothelial", "IbarraSoria2018", "IbarraSoria2018_IID", "code_core_endothelial.b")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
outdir <- ppi_path
setwd(outdir)

source(paste0("https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v", celltype_specific_weight_version, ".R"))
source(paste0("https://raw.githubusercontent.com/xyang2uchicago/BioTIP/refs/heads/master/R/BioTIP_update_", BioTIP_version, ".R"))

load(paste0(data_dir, "sce_16subtype.RData"))

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
    assayName <- assayNames(sce)[1]
    x <- assay(sce, assayName)
    if (max(x) > 20) {
        sce <- scater::logNormCounts(sce, assay.type = assayName)
    } else {
        assayNames(sce)[1] <- 'logcounts'
    }
}
assayName <- "logcounts"

# -----------------------------
# 1) Read IID graphs
# -----------------------------
graph_list <- readRDS(file.path(wd, "results_core_endothelial.b", paste0(db, "_IID_graph_perState_notsimplified.rds")))

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

