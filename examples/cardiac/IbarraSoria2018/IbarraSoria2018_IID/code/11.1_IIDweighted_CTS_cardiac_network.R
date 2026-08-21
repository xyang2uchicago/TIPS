library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(scran)
packageVersion("scran") # 1.37.0
library("SingleCellExperiment")

########## BEGINNING OF USER INPUT ##########

source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID/")
outdir <- file.path(wd, "data")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

db <- "IbarraSoria2018"

db_species <- 10090 # 10090 for mouse, 9606 for human

load(file = "../../data/BioTIP.res.RData")
load("../../data/sce_16subtype.RData")

CTS <- res$CTS.candidate[which(res$significant)]
names(CTS)

########## END OF USER INPUT ##########
# downloaded from GitHub
DEG <- readRDS(file = paste0("../../data/DEG_perState_min.prop0.25_lfc", logFC.cut <- 0.6, "_FDFR0.01.rds"))
lengths(DEG)

markers.up <- readRDS("../../data/markers.up_ttest_min.prop0.25.rds")
#                  blood              cardiac.b              cardiac.c
#                    422                    414                    448
#          endothelial.a          endothelial.c          endothelial.d
#                    456                    481                    427
# extraembryonicMesoderm    mesodermProgenitors        mixedMesoderm.a
#                    369                    362                    335
#        mixedMesoderm.b     pharyngealMesoderm   presomiticMesoderm.a
#                    372                    376                    344
#   presomiticMesoderm.b        somiticMesoderm          endothelial.b
#                    420                    333                    420
#              cardiac.a
#                    409

######################################################
# 3) load IID db and build GRN
################################################################
library(igraph)
library(data.table)

# This file downloaded from:
iid_file <- Sys.glob(paste0(shared_data_path, "/PPIN/human_annotated_PPIs.txt.gz"))
(length(iid_file))
stopifnot(length(iid_file) == 1)

iid <- fread(iid_file[1], sep = "\t", header = TRUE)
iid <- iid[mouse == 1]

iid[, n_exp_pmids := as.integer(n_exp_pmids)]
iid[is.na(n_exp_pmids), n_exp_pmids := 0L]

iid[, w := fifelse(n_exp_pmids >= 100L, 1000, n_exp_pmids * 10) / 1000]

edges <- iid[, .(
  from   = as.character(symbol1),
  to     = as.character(symbol2),
  weight = as.numeric(w)
)]

edges <- edges[!is.na(from) & from != "" & !is.na(to) & to != ""]
edges <- edges[from != to]

edges[, `:=`(
  u = pmin(from, to),
  v = pmax(from, to)
)]

edges <- edges[, .(weight = max(weight, na.rm = TRUE)), by = .(u, v)]

g_iid_global <- graph_from_data_frame(
  edges[, .(from = u, to = v, weight)],
  directed = FALSE
)

(cat("Built g_iid_global:", igraph::vcount(g_iid_global), "nodes,", igraph::ecount(g_iid_global), "edges\n"))
# Built g_iid_global: 17182 nodes, 1517715 edges

graph_list <- list()

get_iid_subnetwork <- function(g_global, hits) {
  hits <- toupper(hits)
  hits <- unique(hits[!is.na(hits) & hits != ""])
  hits <- intersect(hits, V(g_global)$name)
  if (length(hits) < 2) {
    return(make_empty_graph())
  }
  induced_subgraph(g_global, vids = hits)
}


# HiG
for (i in names(DEG)) {
  # Filter differentially expressed genes based on logFC and FDR thresholds
  diff_exp <- markers.up[[i]]
  diff_exp$symbol <- rownames(diff_exp)
  diff_exp <- subset(diff_exp, summary.logFC > logFC.cut & FDR < 0.01)

  hits <- diff_exp$symbol

  # Get the subnetwork
  graph <- get_iid_subnetwork(g_iid_global, hits)

  # Add node attributes (only for nodes that exist in the graph)
  if (vcount(graph) > 0) {
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), "summary.logFC"]
    V(graph)$FDR <- diff_exp[match(V(graph)$name, diff_exp$symbol), "FDR"]
  }

  if (vcount(graph) > 0 && all(toupper(V(graph)$name) %in% toupper(DEG[[i]]))) {
    # map vertex names to the exact form used in DEG[[i]]
    V(graph)$name <- DEG[[i]][match(toupper(V(graph)$name), toupper(DEG[[i]]))]
  }

  graph_list[[paste0("HiG_", i)]] <- graph
}

# instead build for transitory PPI_cats only with the significantly up-regulated CTS genes !!!!!
# build for (up-regulated_marker intersecting CTS)
# HiGCTS
for (i in names(CTS)) {
  # Get unique cluster ids for clusters containing subclusters labeled by numerical id
  j <- if (grepl("\\.", i) && grepl("^[0-9]+$", sub("^[^.]*\\.", "", i))) sub("\\..*$", "", i) else i

  # Get the full marker table for that cluster
  deg_table <- markers.up[[j]]
  deg_table$symbol <- rownames(deg_table)

  # Intersect with CTS
  deg_table <- deg_table[deg_table$symbol %in% CTS[[i]], ]

  # Filter on significance
  diff_exp <- subset(deg_table, summary.logFC > logFC.cut & FDR < 0.01)

  hits <- diff_exp$symbol
  graph <- get_iid_subnetwork(g_iid_global, hits)


  if (vcount(graph) > 0) {
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), "summary.logFC"]
    V(graph)$FDR <- diff_exp[match(V(graph)$name, diff_exp$symbol), "FDR"]
  }

  graph_list[[paste0("HiGCTS_", i)]] <- graph
}

## lastly, build for CTS (no significance filter — raw CTS[[i]] membership,
## matching the STRING arm's CTS loop, which also takes CTS[[i]] as-is)
# refer to 6.3_DE.statistics_CTS.R

# CTS
for (i in names(CTS)) {
  j <- if (grepl("\\.", i) && grepl("^[0-9]+$", sub("^[^.]*\\.", "", i))) sub("\\..*$", "", i) else i
  diff_exp <- markers.up[[j]][CTS[[i]], ]
  diff_exp$symbol <- rownames(diff_exp)

  hits <- diff_exp$symbol
  graph <- get_iid_subnetwork(g_iid_global, hits)


  if (vcount(graph) > 0) {
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), "summary.logFC"]
    V(graph)$FDR <- diff_exp[match(V(graph)$name, diff_exp$symbol), "FDR"]
  }

  graph_list[[paste0("CTS_", i)]] <- graph
}


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

df_graph_info <- data.frame(
  name = names(graph_list),
  vcount = sapply(graph_list, igraph::vcount),
  ecount = sapply(graph_list, igraph::ecount),
  stringsAsFactors = FALSE
)
# (df_graph_info)
#                                                  name vcount ecount
# HiG_blood                                   HiG_blood    374   3118
# HiG_cardiac.b                           HiG_cardiac.b    391   3498
# HiG_cardiac.c                           HiG_cardiac.c    421   3764
# HiG_endothelial.a                   HiG_endothelial.a    424   4382
# HiG_endothelial.c                   HiG_endothelial.c    447   3966
# HiG_endothelial.d                   HiG_endothelial.d    392   3247
# HiG_extraembryonicMesoderm HiG_extraembryonicMesoderm    345   2654
# HiG_mesodermProgenitors       HiG_mesodermProgenitors    330   2187
# HiG_mixedMesoderm.a               HiG_mixedMesoderm.a    309   2134
# HiG_mixedMesoderm.b               HiG_mixedMesoderm.b    343   2750
# HiG_pharyngealMesoderm         HiG_pharyngealMesoderm    345   2717
# HiG_presomiticMesoderm.a     HiG_presomiticMesoderm.a    315   2129
# HiG_presomiticMesoderm.b     HiG_presomiticMesoderm.b    379   3462
# HiG_somiticMesoderm               HiG_somiticMesoderm    305   2162
# HiG_endothelial.b                   HiG_endothelial.b    389   3156
# HiG_cardiac.a                           HiG_cardiac.a    384   3072
# HiGCTS_endothelial.b             HiGCTS_endothelial.b     12      1
# HiGCTS_cardiac.a                     HiGCTS_cardiac.a     13      5
# CTS_endothelial.b                   CTS_endothelial.b     27      5
# CTS_cardiac.a                           CTS_cardiac.a     37     17


saveRDS(graph_list, file = file.path(wd, "results", paste0(db, "_IID_graph_perState_notsimplified.rds"))) # !!!!!!!!!!!!!!!!!!!

graph_list <- readRDS(file = file.path(wd, "results", paste0(db, "_IID_graph_perState_notsimplified.rds")))
graph_list <- lapply(graph_list, simplify) # !!!!!!!!!!!!!!!!!!! # FIXED

# Check which graphs have duplicate vertex names
graphs_with_duplicates <- sapply(graph_list, function(g) {
  vertex_names <- toupper(V(g)$name)
  if (is.null(vertex_names)) {
    # If no names, use vertex indices
    vertex_names <- V(g)
  }
  any(duplicated(vertex_names))
})

# See which graphs have duplicates
which(graphs_with_duplicates)
# named integer(0)

# Show actual edges for duplicated vertices
g1 <- graph_list[[names(graph_list)[1]]]
vertex_names <- toupper(V(g1)$name)
(duplicated_names <- unique(vertex_names[duplicated(vertex_names)]))
# character(0)
if (length(duplicated_names) > 0) {
  for (dup_name in duplicated_names) {
    dup_indices <- which(vertex_names == dup_name)
    all(incident(g1, dup_indices[1], mode = "all") == incident(g1, dup_indices[2], mode = "all"))

    edges <- incident(g1, dup_indices[1], mode = "all")
    edge_list1 <- get.edgelist(g1)[edges, ]
    edge_list1 %>% dim() # [1] 290   2
    weights1 <- E(g1)[edges]$weight

    edges <- incident(g1, dup_indices[2], mode = "all")
    edge_list2 <- get.edgelist(g1)[edges, ]
    edge_list2 %>% dim() # [1] 64  2
    weights2 <- E(g1)[edges]$weight

    all(edge_list2[, 2] %in% edge_list1[, 2]) # FALSE !!
  }
}
cat("end\n")


# w <- unlist(lapply(graph_list, function(g) E(g)$weight))
# table(w)
# #  0.01  0.02  0.03  0.04  0.05  0.06  0.07  0.08  0.09   0.1  0.11  0.12  0.13
# # 56743  3697  1242   561   268   150   128    96    51    33    22    41    39
# #  0.14  0.15  0.17  0.18  0.19  0.21  0.22  0.23  0.24  0.27  0.32  0.34  0.41
# #     9     9    24     7     1    12     2     6     8     3     2     9     6
# #  0.52   0.6
# #     2     1

# sum(table(w))
# #[1] 63172
