library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(readxl)
library(scran)
packageVersion("scran") # 1.37.0
library("SingleCellExperiment")

########## BEGINNING OF USER INPUT ##########

wd = "C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/"
setwd(paste0(wd, "results/"))

db <- "GSE116893"

db_species <- 9606 # 10090 for mouse, 9606 for human

download_files <- FALSE

# Download data files
if(download_files){
    if (!dir.exists("../data/PPIN")) {
        dir.create("../data/PPIN")
    }

    # STRING db
    base_url <- "https://stringdb-static.org/download"
    files <- c(
        paste0("protein.aliases.v12.0/", db_species, ".protein.aliases.v12.0.txt.gz"),
        paste0("protein.info.v12.0/", db_species, ".protein.info.v12.0.txt.gz"),
        paste0("protein.links.v12.0/", db_species, ".protein.links.v12.0.txt.gz")
    )

    for (f in files) {
        download.file(paste0(base_url, "/", f), paste0("../data/PPIN/", basename(f)))
    }
}

load(file = "../data/BioTIP.res.RData")
sce <- readRDS("../data/sce.rds")

CTS <- res$CTS.candidate[which(res$significant)]

# Find subclusters
if (any(duplicated(names(CTS)))) cat('renamed duplicated CTS by extended with ".x" and reorder!')

df <- read_excel("../data/STable4_MES_ADRN_TF.xlsx", skip = 2)

########## END OF USER INPUT ##########

########################
#### DEG generation taken from score400 workflow ####
## load DEGs (or marker genes)
# we used the findMarkers function in scran package, using pairwise Welch t-tests for genes that are detected in a minimum 25% per cluster, with fold changes larger than 2.
# refer to CTS_cardiac_network_robustness_notsimplified.R
########################
logFC.cut <- 0.6

# Reference gene sets for downstream intersection checks
NOTCH_loop     <- c("NOTCH1", "NOTCH2", "NOTCH3", "JAG1", "MAML2", "HES1")
NCC_Cheung2005 <- c("SOX9", "FOXD3", "SNAI2")
MES_TF  <- subset(df, group == "MES")$gene   # 20 TFs
ADRN_TF <- subset(df, group == "ADRN")$gene  # 18 TFs

# Find marker genes
if (!file.exists(paste0("../data/DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDR0.05.rds"))) {
    markers.up <- findMarkers(sce,
        test      = "t",
        groups    = sce$leiden_0.6,
        min.prop  = 0.25,
        direction = "up"
    )
    DEG <- list()

    unique_CTS_ID <- names(CTS)
    unique_CTS_ID <- unique_CTS_ID[!(
        grepl("\\.", unique_CTS_ID) &
        grepl("^[0-9]+$", sub("^[^.]*\\.", "", unique_CTS_ID))
    )]

    for (i in c(setdiff(names(markers.up), names(CTS)), unique_CTS_ID)) {
        interesting.up <- markers.up[[i]]
        DEG[[i]] <- subset(interesting.up, summary.logFC > logFC.cut & FDR < 0.01) %>% rownames()
    }

    saveRDS(DEG,        file = paste0("../data/DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDR0.05.rds"))
    saveRDS(markers.up, file = "../data/markers.up_ttest_min.prop0.25.rds")
}

DEG        <- readRDS(file = paste0("../data/DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDR0.05.rds"))
markers.up <- readRDS(file = "../data/markers.up_ttest_min.prop0.25.rds")
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
    version = "12.0", species = db_species, # species= 10090 for mouse
    score_threshold = 200, # !!!!!!!!!!!default is 200
    network_type = "full",
    input_directory = "../data/PPIN"
)
string_db
# version: 12.0
# species: 9606
# proteins: 19699
# interactions: 7533072

graph_list <- list()
# rather than build for steady PPI_cats
# for(i in setdiff(names(DEG), names(CTS))){
# build for traditional up-regulated markers

# HiG
for (i in names(DEG)) {
    if (length(DEG[[i]]) == 0) next

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

    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$summary.logFC
    V(graph)$FDR    <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$FDR

    # Add edge weights from STRING interaction scores
    # STRING scores are between 0-1000, normalize to 0-1 for better visualization
    E(graph)$weight <- E(graph)$combined_score / 1000
    graph <- delete_edge_attr(graph, "combined_score")

    if (all(mapped$symbol %in% toupper(DEG[[i]]))) V(graph)$name <- DEG[[i]][match(V(graph)$name, toupper(DEG[[i]]))]

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

    mapped <- string_db$map(diff_exp, "symbol", removeUnmappedRows = TRUE)

    hits <- mapped$STRING_id

    graph <- string_db$get_subnetwork(hits)
    all(mapped[match(V(graph)$name, mapped$STRING_id), ]$STRING_id == V(graph)$name) # TRUE
    V(graph)$name <- mapped[match(V(graph)$name, mapped$STRING_id), ]$symbol
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$summary.logFC
    V(graph)$FDR    <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$FDR

    E(graph)$weight <- E(graph)$combined_score / 1000
    graph <- delete_edge_attr(graph, "combined_score")

    if (all(mapped$symbol %in% toupper(DEG[[i]]))) V(graph)$name <- DEG[[i]][match(V(graph)$name, toupper(DEG[[i]]))]

    graph_list[[paste0("HiGCTS_", i)]] <- graph
}

## lastly, build for CTS

# CTS
for (i in names(CTS)) {

    CTS_genes <- data.frame(symbol = CTS[[i]])

    mapped <- string_db$map(CTS_genes, "symbol", removeUnmappedRows = TRUE)
    mapped %>% dim()
    # [1] 110  2
    length(unique(mapped$symbol))
    # [1] 110

    hits <- mapped$STRING_id

    graph <- string_db$get_subnetwork(hits) # t
    # translate STRING_id to symbol
    all(mapped[match(V(graph)$name, mapped$STRING_id), ]$STRING_id == V(graph)$name) # TRUE
    V(graph)$name <- mapped[match(V(graph)$name, mapped$STRING_id), ]$symbol
    # The 'scores' column reflects the strength of the differential expression for each gene based on the Wilcoxon rank-sum test. A high score suggests that the gene's expression is significantly different between the groups under comparison.
    # V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$summary.logFC
    # V(graph)$FDR <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$FDR
    E(graph)$weight <- E(graph)$combined_score / 1000
    graph <- delete_edge_attr(graph, "combined_score") # Remove combined_score as there is no use for it.

    # Fix case issues if needed (for mouse genes)
    if (all(mapped$symbol %in% toupper(CTS[[i]]))) V(graph)$name <- CTS[[i]][match(V(graph)$name, toupper(CTS[[i]]))]

    graph_list[[paste0("CTS_", i)]] <- graph
}


names(graph_list)
#  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"    
#  [7] "HiG_7"     "HiG_8"     "HiG_10"    "HiG_11"    "HiG_12"    "HiG_13"   
# [13] "HiG_14"    "HiG_16"    "HiG_17"    "HiG_15"    "HiG_9"     "HiGCTS_15"
# [19] "HiGCTS_16" "HiGCTS_9"  "CTS_15"    "CTS_16"    "CTS_9"

df_graph_info <- data.frame(
    name = names(graph_list),
    vcount = sapply(graph_list, igraph::vcount),
    ecount = sapply(graph_list, igraph::ecount),
    stringsAsFactors = FALSE
)

(df_graph_info)

#                name vcount ecount
# HiG_1         HiG_1    372  26202
# HiG_2         HiG_2    354   8544
# HiG_3         HiG_3    284   7040
# HiG_4         HiG_4    173   8744
# HiG_5         HiG_5    353   9246
# HiG_6         HiG_6    579  19744
# HiG_7         HiG_7    335   8726
# HiG_8         HiG_8    311   7892
# HiG_10       HiG_10    278  10074
# HiG_11       HiG_11    520  21176
# HiG_12       HiG_12    548  20918
# HiG_13       HiG_13    378  33272
# HiG_14       HiG_14   1044  48268
# HiG_16       HiG_16    555  20478
# HiG_17       HiG_17   1386  75900
# HiG_15       HiG_15    728  46324
# HiG_9         HiG_9    221   4814
# HiGCTS_15 HiGCTS_15     50    762
# HiGCTS_16 HiGCTS_16    124    826
# HiGCTS_9   HiGCTS_9      5      0
# CTS_15       CTS_15     75   1862
# CTS_16       CTS_16    138    986
# CTS_9         CTS_9    110    908

saveRDS(graph_list, file = paste0(db, "_STRING_graph_perState_notsimplified.rds")) # !!!!!!!!!!!!!!!!!!!

graph_list <- readRDS(file = paste0(db, "_STRING_graph_perState_notsimplified.rds"))
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
# HiG_17 
#     15
