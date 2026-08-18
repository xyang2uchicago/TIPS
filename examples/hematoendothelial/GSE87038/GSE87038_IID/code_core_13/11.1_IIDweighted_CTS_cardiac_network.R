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
wd <- tips_path("examples", "hematoendothelial", "GSE87038", "GSE87038_IID/")
outdir <- file.path(wd, "data")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

db <- "GSE87038"

db_species <- 10090 # 10090 for mouse, 9606 for human

load(file = "../../data/BioTIP.res.RData")
load("../../data/sce_E8.25_uncorrected.RData")

CTS <- res$CTS.candidate[which(res$significant)]

# Find subclusters
if (any(duplicated(names(CTS)))) cat('renamed duplicated CTS by extended with ".x" and reorder!')
## manually rename
names(CTS)[7] <- "16.1"
CTS <- c(CTS[1:4], CTS[7], CTS[5:6])
names(CTS)
# [1] "7"    "11"   "15"   "16"   "16.1" "13"   "8"

markers.up <- readRDS("../../data/markers.up_ttest_min.prop0.25.rds")

########## END OF USER INPUT ##########
# downloaded from GitHub
DEG <- readRDS(file = paste0("../../data/DEG_perState_min.prop0.25_lfc", logFC.cut <- 0.6, "_FDFR0.01.rds"))
lengths(DEG)

######################################################
# 3) load IID db and build GRN
################################################################
library(igraph)
library(data.table)

# Mouse PPI file from IID (https://iid.ophid.utoronto.ca/) — mouse_annotated_PPIs.txt.gz
iid_file <- Sys.glob("../data/mouse_annotated_PPIs.txt.gz")
(length(iid_file))
stopifnot(length(iid_file) == 1)

iid <- fread(iid_file[1], sep = "\t", header = TRUE)
iid[, n_exp_pmids := as.integer(n_exp_pmids)]
iid[is.na(n_exp_pmids), n_exp_pmids := 0L]

iid[, w := fifelse(n_exp_pmids >= 100L, 1000, n_exp_pmids * 10) / 1000]

edges <- iid[, .(
  from   = toupper(as.character(symbol1)),
  to     = toupper(as.character(symbol2)),
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
# Built g_iid_global: 17025 nodes, 943097 edges

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

## lastly, build for CTS (node weights not assigned — consistent with STRING 11.1)
# CTS
for (i in names(CTS)) {
  hits <- toupper(CTS[[i]])
  graph <- get_iid_subnetwork(g_iid_global, hits)
  graph_list[[paste0("CTS_", i)]] <- graph
}


names(graph_list)
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_13"
# [26] "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"
# [31] "CTS_16.1"    "CTS_13"      "CTS_8"

df_graph_info <- data.frame(
  name = names(graph_list),
  vcount = sapply(graph_list, igraph::vcount),
  ecount = sapply(graph_list, igraph::ecount),
  stringsAsFactors = FALSE
)

(df_graph_info)
#                    name vcount ecount
# HiG_1             HiG_1    387   2921
# HiG_2             HiG_2    524   4256
# HiG_3             HiG_3    509   3585
# HiG_4             HiG_4    385   2810
# HiG_5             HiG_5    388   2566
# HiG_6             HiG_6    584   4089
# HiG_9             HiG_9    442   3428
# HiG_10           HiG_10    535   4705
# HiG_12           HiG_12    418   2968
# HiG_14           HiG_14    443   2902
# HiG_17           HiG_17    595   3970
# HiG_18           HiG_18    537   3840
# HiG_19           HiG_19    621   3234
# HiG_7             HiG_7    374   2429
# HiG_11           HiG_11    537   4661
# HiG_15           HiG_15    496   3905
# HiG_16           HiG_16    610   4610
# HiG_13           HiG_13    495   3661
# HiG_8             HiG_8    411   3095
# HiGCTS_7       HiGCTS_7     14      0
# HiGCTS_11     HiGCTS_11     23      1
# HiGCTS_15     HiGCTS_15     41     12
# HiGCTS_16     HiGCTS_16     15      1
# HiGCTS_16.1 HiGCTS_16.1     33      2
# HiGCTS_13     HiGCTS_13     14      1
# HiGCTS_8       HiGCTS_8     11      0
# CTS_7             CTS_7     27      3
# CTS_11           CTS_11     47      7
# CTS_15           CTS_15     61     17
# CTS_16           CTS_16     33      6
# CTS_16.1       CTS_16.1     74     29
# CTS_13           CTS_13     55     19
# CTS_8             CTS_8     52     23


saveRDS(graph_list, file = file.path(wd, "results_core_13", paste0(db, "_IID_graph_perState_notsimplified.rds"))) # !!!!!!!!!!!!!!!!!!!

graph_list <- readRDS(file = file.path(wd, "results_core_13", paste0(db, "_IID_graph_perState_notsimplified.rds")))
graph_list <- lapply(graph_list, simplify, edge.attr.comb = "max") # !!!!!!!!!!!!!!!!!!! # FIXED

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
g1 <- graph_list[["HiG_1"]]
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
