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

wd = "/Users/felixyu/Documents/GSE87038_weighted/"
setwd(paste0(wd, "results/"))
score_threshold <- "weight"
PPI_color_platte <- c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

load("../data/sce_E8.25_uncorrected.RData")

load(file = "../data/BioTIP.res.RData")

CTS <- res$CTS.candidate[which(res$significant)]

if (any(duplicated(names(CTS)))) cat('renamed duplicated CTS by extended with ".x" and reorder!')
## manually rename
names(CTS)[7] <- "16.1"
CTS <- c(CTS[1:4], CTS[7], CTS[5:6])
names(CTS)
# [1] "7"    "11"   "15"   "16"   "16.1" "13"   "8"

########################
#### DEG generation taken from score400 workflow ####
## load DEGs (or marker genes)
# refer to BioTIP publication,  GSE87038_markers.R;  evaluate_CTS_GSE130146.R
# we used the findMarkers function in scran package, using pairwise Welch t-tests for genes that are detected in a minimum 25% per cluster, with fold changes larger than 2.
# refer to CTS_cardiac_network_robustness_notsimplified.R
########################
logFC.cut <- 0.6 # 1

load("../data/sce_E8.25_uncorrected.RData")
sce
# class: SingleCellExperiment
# dim: 10938 7240
# metadata(0):
#   assays(2): counts logcounts
table(sce$label, sce$celltype)
#       Somitic mesoderm Intermediate mesoderm ExE mesoderm Paraxial mesoderm Allantois Pharyngeal mesoderm Cardiomyocytes Mesenchyme
# 1               156                    36            0               650         0                   3              0          0
# 2                 0                     0            0                 0         0                   0            334          3
# 3                 0                     0            0                 0         0                   7              2        369
# 4                11                   598          348                31         0                  37              0          0
# 5               551                   192            0                 0         0                   0              0          0
# 6                 0                     0            0                 0         0                   0              0          0
# 7                 0                     0            0                 0         0                   0              0          0
# 8                 0                     5            0                21         0                 742             21          8
# 9                 0                     0            0               519         0                  12              0          0
# 10                0                     0            0                 3         0                   1              0          0
# 11                0                     0            0                 0       179                   0              0          4
# 12                0                     7          276                 0        19                  42              0        122
# 13                0                     0            0                 0         0                   0              0          0
# 14                0                     0           27                 0       328                   0              0         29
# 15                0                     0            0                 0         0                   0              0          0
# 16                0                     0            0                 0         0                   0              0        175
# 17                0                     0            0                 0         0                   0            107          0
# 18                0                     0            0                 0         0                   0              0        203
# 19                0                     0            0                 0         0                   0              0         63
#
#             Haematoendothelial progenitors Endothelium Blood progenitors 1 Blood progenitors 2
# 1                               0           0                   0                   0
# 2                               0           0                   0                   0
# 3                               3           0                   0                   0
# 4                               0           0                   0                   0
# 5                               0           0                   0                   0
# 6                               0         283                   0                   0
# 7                               0           0                   0                 137
# 8                               3           0                   0                   0
# 9                               3           0                   0                   0
# 10                            204          70                   0                   0
# 11                              0           0                   0                   0
# 12                              0           0                   0                   0
# 13                            223           0                   0                   0
# 14                             13           0                   0                   0
# 15                              2           5                  34                  19
# 16                              0           0                   0                   0
# 17                              0           0                   0                   0
# 18                              0           0                   0                   0
# 19                              0           0                   0                   0

# find markers for every cluster compared to all remaining cells, report only the positive one
## begin do not repeat !!!!!!!!!!

markers.up <- findMarkers(sce,
    test = "t", # if wilcox test rather than t-test, get AUC rather than lfc
    groups = sce$label, # lfc=logFC.cut ,
    min.prop = 0.25,
    direction = "up"
) # , block=mnn.all$sample)
DEG <- list()
unique_CTS_ID <- names(CTS)
if (any(grepl(".", names(CTS), fixed = T))) unique_CTS_ID <- unique_CTS_ID[-which(grepl(".", names(CTS), fixed = T))]

for (i in c(setdiff(names(markers.up), names(CTS)), unique_CTS_ID)) {
    interesting.up <- markers.up[[i]]
    DEG[[i]] <- subset(interesting.up, summary.logFC > logFC.cut & FDR < 0.01) %>% rownames()
}

saveRDS(DEG, file = paste0("../data/DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDFR0.05.rds"))
saveRDS(markers.up, file = "../data/markers.up_ttest_min.prop0.25.rds")
dim(markers.up[[1]])


markers.up_all <- findMarkers(sce,
    test = "t", # if wilcox test rather than t-test, get AUC rather than lfc
    groups = sce$label,
    min.prop = NULL,
    direction = "up"
) # , block=mnn.all$sample)
saveRDS(markers.up_all, file = "../data/markers.up_all_ttest.rds")
dim(markers.up[[1]])


DEG <- readRDS(file = paste0("../data/DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDFR0.05.rds"))
lengths(DEG)

######################################################
# 3) load STRING db and build GRN
## second trial, using STRING, IFT57/74/88 are all included
# https://www.bioconductor.org/packages/release/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf
# refer to CTS_cardiac_ntwork_robustness_notsimplified.R
################################################################
library(BioNet)
packageVersion("BioNet") # '1.69.0'
library(igraph)
library("STRINGdb")
packageVersion("STRINGdb") # '2.21.0'
library(tibble)

string_db <- STRINGdb$new(
    version = "12.0", species = 10090, # species= 10090 for mouse
    score_threshold = 200, # !!!!!!!!!!!default is 200
    network_type = "full",
    input_directory = "../data/PPIN"
)
string_db
# version: 12.0
# species: 10090
# proteins: 21840
# interactions: 6859804>

graph_list <- list()
# rather than build for steady PPI_cats
# for(i in setdiff(names(DEG), names(CTS))){
# build for traditional up-regulated markers

# HiG
for (i in names(DEG)) {
    # Filter differentially expressed genes based on logFC and FDR thresholds
    diff_exp <- markers.up[[i]]
    diff_exp$symbol <- rownames(diff_exp)
    diff_exp <- subset(diff_exp, summary.logFC > logFC.cut & FDR < 0.01)

    # Map to STRING
    mapped <- string_db$map(diff_exp, "symbol", removeUnmappedRows = TRUE)
    hits <- mapped$STRING_id

    # Get the subnetwork
    graph <- string_db$get_subnetwork(hits)

    ## Handle missing nodes (if any)  ## e.g., due to graph missing the [1] "GM1673" in DEG[['5']]
    if (vcount(graph) < length(hits)) {
        missing_node <- setdiff(hits, V(graph)$name)
        if (length(missing_node) > 0) mapped <- mapped[-which(mapped$STRING_id %in% missing_node), ]
    }
    # Verify mapping integrity
    stopifnot(all(mapped[match(V(graph)$name, mapped$STRING_id), ]$STRING_id == V(graph)$name))

    # Add node attributes: replace STRING IDs with gene symbols
    V(graph)$name <- mapped[match(V(graph)$name, mapped$STRING_id), ]$symbol

    # Add node weights from logFC values and additioal node attributes
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$summary.logFC # summary.logFC if scanpy using t.test

    V(graph)$FDR <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$FDR # FDR of scanpy_wilcox marker identification
    # E(graph)

    # Add edge weights from STRING interaction scores #!!!!!!!!!!!!!!!!!
    # STRING scores are between 0-1000, normalize to 0-1 for better visualization
    E(graph)$weight <- E(graph)$combined_score / 1000

    # Fix case issues if needed (for mouse genes)
    if (all(mapped$symbol %in% toupper(DEG[[i]]))) V(graph)$name <- DEG[[i]][match(V(graph)$name, toupper(DEG[[i]]))]

    # Optional: Set visual attributes based on weights
    V(graph)$size <- scales::rescale(V(graph)$weight, to = c(3, 15))
    E(graph)$width <- scales::rescale(E(graph)$weight, to = c(0.5, 3)) # useless, using weight directly otherwise can't compare across clusters
    if (!is.null(V(graph)$weight) && length(V(graph)$weight) > 0 && any(is.finite(V(graph)$weight))) {
        V(graph)$color <- colorRampPalette(c("lightblue", "red"))(100)[
            cut(V(graph)$weight, breaks = 100)
        ]
    } else {
        V(graph)$color <- "grey" # fallback for empty or NA weights
    }

    graph_list[[i]] <- graph
}
names(graph_list) <- paste0("HiG_", names(DEG))

# instead build for transitory PPI_cats only with the significantly up-regulated CTS genes !!!!!
# build for (up-regulated_marker intersecting CTS)
# HiGCTS
for (i in names(CTS)) {
    if (grepl(".", i, fixed = T)) j <- unlist(strsplit(i, split = ".", fixed = T))[1] else j <- i

    # Get the full marker table for that cluster
    deg_table <- markers.up[[j]]
    deg_table$symbol <- rownames(deg_table)

    # Intersect with CTS
    deg_table <- deg_table[deg_table$symbol %in% CTS[[i]], ]

    # Filter on significance
    diff_exp <- subset(deg_table, summary.logFC > logFC.cut & FDR < 0.01)

    mapped <- string_db$map(diff_exp, "symbol", removeUnmappedRows = TRUE)

    hits <- mapped$STRING_id
    # string_db$plot_network( hits )

    graph <- string_db$get_subnetwork(hits) # t
    # translate STRING_id to symbol
    all(mapped[match(V(graph)$name, mapped$STRING_id), ]$STRING_id == V(graph)$name) # TRUE
    V(graph)$name <- mapped[match(V(graph)$name, mapped$STRING_id), ]$symbol
    # The 'scores' column reflects the strength of the differential expression for each gene based on the Wilcoxon rank-sum test. A high score suggests that the gene's expression is significantly different between the groups under comparison.
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$summary.logFC
    V(graph)$FDR <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$FDR

    # Add edge weights from STRING interaction scores #!!!!!!!!!!!!!!!!!
    # STRING scores are between 0-1000, normalize to 0-1 for better visualization
    E(graph)$weight <- E(graph)$combined_score / 1000

    # Fix case issues if needed (for mouse genes)
    if (all(mapped$symbol %in% toupper(DEG[[i]]))) V(graph)$name <- DEG[[i]][match(V(graph)$name, toupper(DEG[[i]]))]

    # Optional: Set visual attributes based on weights
    V(graph)$size <- scales::rescale(V(graph)$weight, to = c(3, 15))
    E(graph)$width <- scales::rescale(E(graph)$weight, to = c(0.5, 3))
    if (!is.null(V(graph)$weight) && length(V(graph)$weight) > 0 && any(is.finite(V(graph)$weight))) {
        V(graph)$color <- colorRampPalette(c("lightblue", "red"))(100)[
            cut(V(graph)$weight, breaks = 100)
        ]
    } else {
        V(graph)$color <- "grey" # fallback for empty or NA weights
    }

    graph_list[[paste0("HiGCTS_", i)]] <- graph
}

## lastly, build for CTS
# refer to 6.3_DE.statistics_CTS.R
markers.up_all <- readRDS("../data/markers.up_all_ttest.rds")

# CTS
for (i in names(CTS)) {
    if (grepl(".", i, fixed = T)) j <- unlist(strsplit(i, split = ".", fixed = T))[1] else j <- i
    diff_exp <- markers.up_all[[j]][CTS[[i]], ]
    diff_exp$symbol <- rownames(diff_exp)
    mapped <- string_db$map(diff_exp, "symbol", removeUnmappedRows = TRUE)
    mapped %>% dim()
    # [1] 84  21
    length(unique(mapped$symbol))
    # [1] 84

    hits <- mapped$STRING_id
    # string_db$plot_network( hits )

    graph <- string_db$get_subnetwork(hits) # t
    # translate STRING_id to symbol
    all(mapped[match(V(graph)$name, mapped$STRING_id), ]$STRING_id == V(graph)$name) # TRUE
    V(graph)$name <- mapped[match(V(graph)$name, mapped$STRING_id), ]$symbol
    # The 'scores' column reflects the strength of the differential expression for each gene based on the Wilcoxon rank-sum test. A high score suggests that the gene's expression is significantly different between the groups under comparison.
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$summary.logFC
    V(graph)$FDR <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$FDR
    E(graph)$weight <- E(graph)$combined_score / 1000

    # Optional: Set visual attributes based on weights
    V(graph)$size <- scales::rescale(V(graph)$weight, to = c(3, 15))
    E(graph)$width <- scales::rescale(E(graph)$weight, to = c(0.5, 3))
    if (!is.null(V(graph)$weight) && length(V(graph)$weight) > 0 && any(is.finite(V(graph)$weight))) {
        V(graph)$color <- colorRampPalette(c("lightblue", "red"))(100)[
            cut(V(graph)$weight, breaks = 100)
        ]
    } else {
        V(graph)$color <- "grey" # fallback for empty or NA weights
    }

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
# HiG_1             HiG_1    303   9956
# HiG_2             HiG_2    435  19088
# HiG_3             HiG_3    411  14240
# HiG_4             HiG_4    304  10724
# HiG_5             HiG_5    322  11622
# HiG_6             HiG_6    512  21986
# HiG_9             HiG_9    358  13874
# HiG_10           HiG_10    422  17310
# HiG_12           HiG_12    341  11656
# HiG_14           HiG_14    364  14008
# HiG_17           HiG_17    524  23448
# HiG_18           HiG_18    457  17098
# HiG_19           HiG_19    529  15668
# HiG_7             HiG_7    336  15144
# HiG_11           HiG_11    441  18602
# HiG_15           HiG_15    406  17264
# HiG_16           HiG_16    524  21162
# HiG_13           HiG_13    403  17016
# HiG_8             HiG_8    332  11760
# HiGCTS_7       HiGCTS_7     14     34
# HiGCTS_11     HiGCTS_11     19     30
# HiGCTS_15     HiGCTS_15     30    108
# HiGCTS_16     HiGCTS_16     13     20
# HiGCTS_16.1 HiGCTS_16.1     30     62
# HiGCTS_13     HiGCTS_13     13     18
# HiGCTS_8       HiGCTS_8     10     20
# CTS_7             CTS_7     31    104
# CTS_11           CTS_11     51    108
# CTS_15           CTS_15     66    414
# CTS_16           CTS_16     39    140
# CTS_16.1       CTS_16.1     79    420
# CTS_13           CTS_13     60    294
# CTS_8             CTS_8     54    330

saveRDS(graph_list, file = "GSE87038_STRING_graph_perState_notsimplified.rds") # !!!!!!!!!!!!!!!!!!!

graph_list <- readRDS(file = "GSE87038_STRING_graph_perState_notsimplified.rds")
graph_list <- lapply(graph_list, simplify, edge.attr.comb ='max') # !!!!!!!!!!!!!!!!!!! # FIXED

# Check which graphs have duplicate vertex names
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
#       HiG_5      HiG_17      HiG_18      HiG_19    HiGCTS_7   HiGCTS_11
#           5          11          12          13          20          21
#   HiGCTS_15   HiGCTS_16 HiGCTS_16.1   HiGCTS_13    HiGCTS_8       CTS_7
#          22          23          24          25          26          27
#      CTS_11      CTS_15      CTS_16    CTS_16.1      CTS_13       CTS_8
#          28          29          30          31          32          33

# Show actual edges for duplicated vertices
g1 <- graph_list[["HiG_1"]]
vertex_names <- V(g1)$name
(duplicated_names <- unique(vertex_names[duplicated(vertex_names)]))
#  "H3F3B"
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
