library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(igraph)

########## BEGINNING OF USER INPUT ##########

wd = "C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/"
setwd(paste0(wd, "results/"))
score_threshold <- "weight"

db <- "GSE116893"
db_species <- 9606 # 10090 for mouse, 9606 for human

# ATTENTION: MANUAL INPUT REQUIRED AFTER LINE 154

########## END OF USER INPUT ##########

graph_list <- readRDS(file = paste0(db, "_STRING_graph_perState_notsimplified.rds"))
graph_list <- lapply(graph_list, simplify, edge.attr.comb ='max') # !!!!!!!!!!!!!!!!!!! # FIXED

N <- sapply(graph_list, vcount)


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

# Show actual edges for duplicated vertices
g1 <- graph_list[["HiG_17"]]
vertex_names <- V(g1)$name
(duplicated_names <- unique(vertex_names[duplicated(vertex_names)]))

duplicated_names
#  "LINC02210-CRHR1"

# refer to ../xxx_score200/11.1.1_CTS_cardiac_network_degreeDistribution.R
library("STRINGdb")
packageVersion("STRINGdb") # '2.21.0'
string_db <- STRINGdb$new(
    version = "12.0", species = db_species,
    score_threshold = 200,
    network_type = "full",
    input_directory = "../data/PPIN"
)
string_db

DEG <- readRDS(paste0("../data/DEG_perState_min.prop0.25_lfc0.6_FDFR0.05.rds"))

DEG <- lapply(DEG, function(x) data.frame(names = x, stringsAsFactors = FALSE))

lapply(DEG, nrow) %>% unlist()

DEG_list <- lapply(DEG, unique)

## build PPIN again and track back the correct number of edges for the duplicated gene ##

markers.up <- readRDS(paste0("../data/DEG_perState_min.prop0.25_lfc0.6_FDFR0.05.rds"))

graph_list <- readRDS(file = paste0(db, "_STRING_graph_perState_notsimplified.rds"))
graph_list <- lapply(graph_list, simplify, edge.attr.comb ='max') # !!!!!!!!!!!!!!!!!!! # FIXED

N <- sapply(graph_list, vcount)

unmapped_genes <- c() # initialize vector

# BioMart database used to map STRING protein id to official gene id

library(biomaRt)

# Set up biomart connection
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to map a STRING protein ID to Ensembl gene and symbol
map_protein_to_gene <- function(string_id) {
    peptide_id <- sub("9606\\.", "", string_id)

    result <- getBM(
        attributes = c("ensembl_peptide_id", "ensembl_gene_id", "external_gene_name"),
        filters = "ensembl_peptide_id",
        values = peptide_id,
        mart = mart
    )

    if (nrow(result) == 0) {
        return(data.frame(
            STRING_id = string_id,
            gene_id = NA,
            gene_symbol = NA,
            stringsAsFactors = FALSE
        ))
    } else {
        return(data.frame(
            STRING_id = string_id,
            gene_id = result$ensembl_gene_id[1],
            gene_symbol = result$external_gene_name[1],
            stringsAsFactors = FALSE
        ))
    }
}

correct_n_edges <- NULL

# MANUAL INPUT REQUIRED
for (i in sub("^HiG_", "", names(which(graphs_with_duplicates)))) {
    cat("\n")
    cat("Analyzing HiG_", i, "\n")

    diff_exp <- DEG[[i]]  # data frame with 'names' column

    mapped <- string_db$map(diff_exp, "names", removeUnmappedRows = TRUE)

    # For cluster 17: current STRING aliases no longer map LINC02210-CRHR1 (a read-through
    # lncRNA overlapping CRHR1 on chr17). The saved graph was built when STRING returned
    # both ENSP00000488912 (correct lncRNA entry) and ENSP00000381333 (CRHR1 protein,
    # locus-overlap artifact). Inject both so the duplicate is detected and the correct
    # edge count is stored. The removal block below then discards ENSP00000381333.
    if (i == "17") {
        existing_linc <- mapped$STRING_id[mapped$names == "LINC02210-CRHR1"]
        if (!"9606.ENSP00000488912" %in% existing_linc)
            mapped <- rbind(mapped, data.frame(names = "LINC02210-CRHR1",
                                               STRING_id = "9606.ENSP00000488912",
                                               stringsAsFactors = FALSE))
        if (!"9606.ENSP00000381333" %in% existing_linc)
            mapped <- rbind(mapped, data.frame(names = "LINC02210-CRHR1",
                                               STRING_id = "9606.ENSP00000381333",
                                               stringsAsFactors = FALSE))
    }

    dup_name <- unique(mapped$names[duplicated(mapped$names)])
    cat("Duplicated gene names:", paste(dup_name, collapse = ", "), "\n")

    dup_entries <- mapped[mapped$names %in% dup_name, ]

    # Apply biomart function to get gene ID and symbol info
    annotation <- do.call(rbind, lapply(dup_entries$STRING_id, map_protein_to_gene))

    # Combine for biological validation
    annotated_dup <- cbind(dup_entries, annotation[, c("gene_id", "gene_symbol")])

    print(annotated_dup)

    # In the manual code below, we want to delete all duplicate gene names that are associated with
    # the wrong Ensemble ID (those that are not in the official database corresponding to our gene).

    # Find official human gene symbols: https://www.genenames.org/
    # Find Ensemble ID at NCBI: https://www.ncbi.nlm.nih.gov
    # Map protein id to gene id using BioMart: https://useast.ensembl.org/info/data/biomart/index.html


    #               names            STRING_id         gene_id     gene_symbol
    # 96            H3F3A 9606.ENSP00000355780 ENSG00000163041           H3-3A ✅
    # 97            H3F3A 9606.ENSP00000254810 ENSG00000132475           H3-3B
    # 388             HN1 9606.ENSP00000498587 ENSG00000148400          NOTCH1
    # 389             HN1 9606.ENSP00000439228            <NA>            <NA>
    # 390             HN1 9606.ENSP00000348316 ENSG00000189159            JPT1 ✅
    # 513 LINC02210-CRHR1 9606.ENSP00000488912 ENSG00000263715 LINC02210-CRHR1 ✅ (correct lncRNA entry)
    # 514 LINC02210-CRHR1 9606.ENSP00000381333 ENSG00000120088           CRHR1    (chr17 locus overlap, remove)

    # Official Symbol: H3-3A / Official Ensemble ID: ENSG00000163041
    # Official Symbol: JPT1 / Official Ensemble ID: ENSG00000189159
    # Official Symbol: LINC02210-CRHR1 / Official Ensemble ID: ENSG00000263715

    if (i == "17") {
        # Use STRING ID-based removal to avoid positional index shift bugs
        mapped <- mapped[mapped$STRING_id != "9606.ENSP00000254810", ]  # H3F3A: remove H3-3B, keep H3-3A
        mapped <- mapped[mapped$STRING_id != "9606.ENSP00000498587", ]  # HN1: remove NOTCH1
        mapped <- mapped[mapped$STRING_id != "9606.ENSP00000439228", ]  # HN1: remove unresolved entry
        mapped <- mapped[mapped$STRING_id != "9606.ENSP00000381333", ]  # LINC02210-CRHR1: remove CRHR1 locus overlap
    }


    # Rebuild graph after removal, translate STRING IDs to gene symbols
    hits <- mapped$STRING_id
    graph <- string_db$get_subnetwork(hits)
    all(mapped[match(V(graph)$name, mapped$STRING_id), ]$STRING_id == V(graph)$name) # TRUE
    V(graph)$name <- mapped[match(V(graph)$name, mapped$STRING_id), ]$names

    for (j in seq_along(dup_name)) {
        edges <- incident(graph, dup_name[j], mode = "all")
        n <- get.edgelist(graph)[edges, ]

        res <- data.frame(
            "graph_id" = paste0("HiG_", i),
            "names"    = dup_name[j],
            "n_edge"   = nrow(n),
            "STRING_id" = subset(mapped, names == dup_name[j])$STRING_id
        )

        if (is.null(correct_n_edges)) {
            correct_n_edges <- res
        } else {
            correct_n_edges <- rbind(correct_n_edges, res)
        }
    }
}

(correct_n_edges)
#   graph_id           names n_edge            STRING_id vertex_index_to_remove
# 1   HiG_17           H3F3A     68 9606.ENSP00000355780                       
# 2   HiG_17             HN1     12 9606.ENSP00000348316                       
# 3   HiG_17 LINC02210-CRHR1     20 9606.ENSP00000488912                    921

################################################################
## remove duplicated vertex directly from un-simplified graph ##
################################################################
graph_list <- readRDS(file = paste0(db, "_STRING_graph_perState_notsimplified.rds"))

N <- sapply(graph_list, vcount)
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

# initialization
correct_n_edges$vertex_index_to_remove <- vector("list", nrow(correct_n_edges))

for (i in names(which(graphs_with_duplicates))) {
    g_name <- i
    cat("Processing graph:", g_name, "\n")

    g <- graph_list[[g_name]]
    vertex_names <- V(g)$name

    duplicated_names <- unique(vertex_names[duplicated(vertex_names)])
    cat("Duplicated gene names found:", paste(duplicated_names, collapse = ", "), "\n")

    if (length(duplicated_names) == 0) {
        cat("No duplicates in this graph. Skipping.\n")
        next
    }

    for (dup_name in duplicated_names) {
        dup_indices <- which(vertex_names == dup_name)

        n <- array(dim = length(dup_indices))
        for (j in seq_along(dup_indices)) {
            idx <- dup_indices[j]
            n[j] <- length(incident(g, idx, mode = "all"))
        }

        x <- which(correct_n_edges$graph_id == g_name &
            toupper(correct_n_edges$names) %in% toupper(dup_name))

        if (length(x) == 0) {
            cat("    WARNING: No match found in correct_n_edges for", dup_name, "in graph", g_name, "\n")
            next
        }

        # Make sure both vectors are numeric and non-empty
        observed <- as.numeric(n)
        expected <- as.numeric(correct_n_edges$n_edge[x])

        # If expected has more entries than observed, just compare to the first few
        if (length(expected) > length(observed)) {
            expected <- expected[seq_along(observed)]
        }

        # Compute diffs safely
        diffs <- observed - expected
        keep_index <- which.min(abs(diffs))

        cat("    Keep index in `dup_indices`:", keep_index, "\n")
        cat("    Removing indices:", paste(dup_indices[-keep_index], collapse = ", "), "\n")

        correct_n_edges$vertex_index_to_remove[x] <- list(dup_indices[-keep_index])
    }
}

# Final check
# Identify entries with zero (or containing zero)
zero_indices <- sapply(correct_n_edges$vertex_index_to_remove, function(x) {
    length(x) == 1 && x == 0
})

if (any(zero_indices)) {
    cat("\nZeroes detected in `vertex_index_to_remove`, setting to NA...\n")
    correct_n_edges$vertex_index_to_remove[zero_indices] <- NA
}



correct_n_edges_copy <- correct_n_edges

# Convert the list column to character, collapsing each vector by commas
correct_n_edges_copy$vertex_index_to_remove <- sapply(
    correct_n_edges_copy$vertex_index_to_remove,
    function(x) {
        if (is.null(x) || all(is.na(x))) {
            return(NA_character_)
        }
        paste(x, collapse = ",")
    }
)

write.table(correct_n_edges_copy, file = "correct_n_edges_HiG_STRING2.14.0.txt", sep = "\t", row.names = FALSE)

saveRDS(correct_n_edges_copy, file = "correct_n_edges_HiG_STRING2.14.0.rds")
