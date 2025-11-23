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

# Set working directory
wd <- "/Users/felixyu/Documents/IbarraSoria2018/"
setwd(paste0(wd, "results/"))

# Database settings
db <- "IbarraSoria2018"
db_species <- 10090 # 10090 for mouse, 9606 for human

# Load BioTIP results and sce object
load(file = "../data/BioTIP.res.RData")
load("../data/sce_16subtype.RData") 

########## END OF USER INPUT ##########

# Extract significant CTS clusters
CTS <- res$CTS.candidate[which(res$significant)]
names(CTS)
# [1] "endothelial.b" "cardiac.a"

logFC.cut <- 0.6

sce
# class: SingleCellExperiment 
# dim: 4000 11039 

table(sce$label, sce$celltype)

# Generate markers for each cluster (positive only)
markers.up <- findMarkers(
    sce,
    test = "t",            # Welch t-test
    groups = sce$subcelltype,    # use subcelltype labels from sce
    min.prop = 0.25,       # only consider genes expressed in >=25% of cells
    direction = "up"
)

DEG <- list()

# Handle potential "." in cluster names
unique_CTS_ID <- names(CTS)

# Subclusters are labeled numerically as <cluster_name.subcluster_number>
# remove only if the part AFTER the '.' is numeric
unique_CTS_ID <- unique_CTS_ID[ !(
    grepl("\\.", unique_CTS_ID) &
    grepl("^[0-9]+$", sub("^[^.]*\\.", "", unique_CTS_ID))
) ]


# Assign DEGs per cluster
for (i in c(setdiff(names(markers.up), names(CTS)), unique_CTS_ID)) {
    interesting.up <- markers.up[[i]]
    DEG[[i]] <- rownames(interesting.up[interesting.up$summary.logFC > logFC.cut & interesting.up$FDR < 0.01, ])
}

saveRDS(DEG, file = paste0("../data/DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDFR0.05.rds"))
saveRDS(markers.up, file = "../data/markers.up_ttest_min.prop0.25.rds")

dim(markers.up[[1]]) # just to check dimensions

# Compute all markers without min.prop filter ---
markers.up_all <- findMarkers(
    sce,
    test = "t",
    groups = sce$subcelltype,
    min.prop = NULL,
    direction = "up"
)
saveRDS(markers.up_all, file = "../data/markers.up_all_ttest.rds")
dim(markers.up_all[[1]])

DEG <- readRDS(file = paste0("../data/DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDFR0.05.rds"))
lengths(DEG)

######################################################
# 3) load STRING db and build GRN 
## second trial, using STRING, IFT57/74/88 are all included
# https://www.bioconductor.org/packages/release/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf
# refer to CTS_cardiac_ntwork_robustness_notsimplified.R
################################################################
library(BioNet)
packageVersion('BioNet') # '1.56.0'
library(igraph)
library("STRINGdb")
packageVersion('STRINGdb') # '2.14.0'
library(tibble)

string_db <- STRINGdb$new( version="12.0", species=10090,   # species= 10090 for mouse;    ?? for human 2021 version
                        score_threshold = 200, #!!!!!!!!!!!default is 200
                        network_type="full", 
                        input_directory="../data/PPIN") 
string_db
# version 12.0  / 11.5
# proteins: 19699    / for mouse: 22048
# interactions:  7533072   / for mouse 6859804

graph_list <- list()

# HiG
for (i in names(DEG)) {
    # Filter differentially expressed genes based on logFC and FDR thresholds
    diff_exp <- markers.up[[i]]
    diff_exp$symbol <- rownames(diff_exp)
    diff_exp <- diff_exp[diff_exp$summary.logFC > logFC.cut & diff_exp$FDR < 0.01, ]

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
    graph <- delete_edge_attr(graph, "combined_score") # Remove combined_score as there is no use for it.

    # Fix case issues if needed (for mouse genes)
    if (all(mapped$symbol %in% toupper(DEG[[i]]))) V(graph)$name <- DEG[[i]][match(V(graph)$name, toupper(DEG[[i]]))]

    graph_list[[i]] <- graph
}
names(graph_list) <- paste0("HiG_", names(DEG))

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
    diff_exp <- deg_table[deg_table$summary.logFC > logFC.cut & deg_table$FDR < 0.01, ]

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
    graph <- delete_edge_attr(graph, "combined_score") # Remove combined_score as there is no use for it.

    # Fix case issues if needed (for mouse genes)
    if (all(mapped$symbol %in% toupper(DEG[[i]]))) V(graph)$name <- DEG[[i]][match(V(graph)$name, toupper(DEG[[i]]))]

    graph_list[[paste0("HiGCTS_", i)]] <- graph
}

markers.up_all <- readRDS("../data/markers.up_all_ttest.rds")

# CTS
for (i in names(CTS)) {
    # Get unique cluster ids for clusters containing subclusters labeled by numerical id
    j <- if (grepl("\\.", i) && grepl("^[0-9]+$", sub("^[^.]*\\.", "", i))) sub("\\..*$", "", i) else i

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
    graph <- delete_edge_attr(graph, "combined_score") # Remove combined_score as there is no use for it.

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
# [15] "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"          
# [17] "CTS_endothelial.b"          "CTS_cardiac.a"    
 

 saveRDS(graph_list, file= paste0(db, '_STRING_graph_perState_notsimplified.rds'))  #!!!!!!!!!!!!!!!!!!!
 
 graph_list <- readRDS( file= paste0(db, '_STRING_graph_perState_notsimplified.rds'))  
 graph_list <- lapply(graph_list, simplify) #!!!!!!!!!!!!!!!!!!!

# Check which graphs have duplicate vertex names
graphs_with_duplicates <- sapply(graph_list, function(g) {
  vertex_names <- V(g)$name
  if(is.null(vertex_names)) {
    # If no names, use vertex indices
    vertex_names <- V(g)
  }
  any(duplicated(vertex_names))
})

# See which graphs have duplicates
which(graphs_with_duplicates)
# named integer(0)