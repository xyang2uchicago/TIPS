library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr)
library(igraph)

########## BEGINNING OF USER INPUT ##########

code_dir <- here::here("examples", "hematoendothelial", "GSE87038", "GSE87038_STRING", "code_core_13")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)

# ATTENTION: MANUAL INPUT REQUIRED AFTER LINE ~150 (search "MANUAL INPUT
# REQUIRED" below). This filters by verified-bad STRING protein IDs (content
# keyed, not row-position keyed like some sibling folders' 11.1.1), so it is
# reasonably robust to DEG size changing between runs — but if BioMart/STRING
# ever surfaces a NEW duplicate gene not in the current bad_ids list, this
# will need to be re-derived by hand.

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
#  HiG_1  HiG_2  HiG_3  HiG_4  HiG_5  HiG_6  HiG_9 HiG_17 HiG_18 HiG_19 HiG_16  HiG_8
#      1      2      3      4      5      6      7     11     12     13     17     19
# (12 graphs at lfc=0.5 vs 4 at lfc=0.6; Hn1l gains HiG_1/3/4/9/16/8, Cxx1a gains HiG_2/6, Hist1h4i new in HiG_19)

# Show actual edges for duplicated vertices
g1 <- graph_list[["HiG_5"]]
vertex_names <- V(g1)$name
(duplicated_names <- unique(vertex_names[duplicated(vertex_names)]))

duplicated_names
#  "Hn1l"

if (length(duplicated_names) > 0) {
    for (dup_name in duplicated_names) {
        dup_indices <- which(vertex_names == dup_name)
        # cat("Edges for duplicated vertex '", dup_name, "':\n", sep = "")
        all(incident(g1, dup_indices[1], mode = "all") == incident(g1, dup_indices[2], mode = "all"))

        # First instance
        edges1 <- incident(g1, dup_indices[1], mode = "all")
        edge_list1 <- get.edgelist(g1)[edges1, , drop = FALSE]
        if (is.null(dim(edge_list1))) {
            edge_list1 <- matrix(edge_list1, ncol = 2)
        }
        dim(edge_list1) # 3 2
        weights1 <- E(g1)[edges1]$weight

        # Second instance
        edges2 <- incident(g1, dup_indices[2], mode = "all")
        edge_list2 <- get.edgelist(g1)[edges2, , drop = FALSE]
        if (is.null(dim(edge_list2))) {
            edge_list2 <- matrix(edge_list2, ncol = 2)
        }
        dim(edge_list2) # 1 2
        weights2 <- E(g1)[edges2]$weight

        all(edge_list2[, 2] %in% edge_list1[, 2]) # FALSE !!
    }
}


### why "Hn1l" is duplicated in 'HiG_5'?
# refer to ../xxx_score200/11.1.1_CTS_cardiac_network_degreeDistribution.R
library("STRINGdb")
packageVersion("STRINGdb") # '2.21.0'
string_db <- STRINGdb$new(
    version = "12.0", species = db_species,
    score_threshold = score_threshold,
    network_type = "full",
    input_directory = paste0(shared_data_path, "/PPIN")
)
string_db

DEG <- readRDS(paste0(data_dir, "DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDFR0.01.rds"))

DEG <- lapply(DEG, function(x) data.frame(names = x, stringsAsFactors = FALSE))

lapply(DEG, nrow) %>% unlist()

DEG_list <- lapply(DEG, unique)

(lengths(DEG[["5"]])) # 1
any(duplicated(DEG[["5"]])) # [1] FALSE

## build PPIN again and track back the correct number of edges for the duplicated gene ##

markers.up <- readRDS(paste0(data_dir, "DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDFR0.01.rds"))

graph_list <- readRDS(file = paste0(db, "_STRING_graph_perState_notsimplified.rds"))
graph_list <- lapply(graph_list, simplify, edge.attr.comb ='max') # !!!!!!!!!!!!!!!!!!! # FIXED

N <- sapply(graph_list, vcount)

unmapped_genes <- c() # initialize vector

# BioMart database used to map STRING protein id to official gene id

library(biomaRt)

# Set up biomart connection
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Function to map a STRING protein ID to Ensembl gene and symbol
map_protein_to_gene <- function(string_id) {
    peptide_id <- sub("10090\\.", "", string_id)

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

    # Find official mouse gene symbols: https://www.informatics.jax.org/
    # Find Ensemble ID at NCBI: https://www.ncbi.nlm.nih.gov
    # Map protein id to gene id using BioMart: https://useast.ensembl.org/info/data/biomart/index.html

    # Biomart-verified incorrect STRING protein IDs — filter by STRING_id (robust to any DEG size).
    # Applies across all affected clusters; no-op for clusters that don't contain these genes.
    #
    # HN1L  → keep ENSMUSP00000024981 (Jpt2, ENSMUSG00000024165) ✅
    #          remove ENSMUSP00000078129 (AY358078) and ENSMUSP00000156906 (Cramp1)
    #          affected: HiG_1, HiG_3, HiG_4, HiG_5, HiG_9, HiG_16, HiG_18, HiG_19, HiG_8
    #
    # CXX1A → keep ENSMUSP00000086158 (Rtl8a, ENSMUSG00000067925) ✅
    #          remove ENSMUSP00000086157 (Rtl8b)
    #          affected: HiG_2, HiG_6, HiG_17
    #
    # HIST1H4I → keep ENSMUSP00000100042 (H4c9/Hist1h4i, ENSMUSG00000060639) ✅
    #             remove ENSMUSP00000085006 (H4c11/Hist1h4k, ENSMUSG00000067455)
    #             affected: HiG_19
    bad_ids <- c(
        "10090.ENSMUSP00000078129",
        "10090.ENSMUSP00000156906",
        "10090.ENSMUSP00000086157",
        "10090.ENSMUSP00000085006"
    )
    mapped <- mapped[!(mapped$STRING_id %in% bad_ids), ]

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
#    graph_id    names n_edge                STRING_id
# 1     HiG_1     HN1L      8 10090.ENSMUSP00000024981
# 2     HiG_2    CXX1A     12 10090.ENSMUSP00000086158
# 3     HiG_3     HN1L      4 10090.ENSMUSP00000024981
# 4     HiG_4     HN1L      8 10090.ENSMUSP00000024981
# 5     HiG_5     HN1L      6 10090.ENSMUSP00000024981
# 6     HiG_6    CXX1A      6 10090.ENSMUSP00000086158
# 7     HiG_9     HN1L      6 10090.ENSMUSP00000024981
# 8    HiG_17    CXX1A      8 10090.ENSMUSP00000086158
# 9    HiG_18     HN1L      4 10090.ENSMUSP00000024981
# 10   HiG_19 HIST1H4I     62 10090.ENSMUSP00000100042
# 11   HiG_19     HN1L      4 10090.ENSMUSP00000024981
# 12   HiG_16     HN1L      8 10090.ENSMUSP00000024981
# 13    HiG_8     HN1L      8 10090.ENSMUSP00000024981

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
#  HiG_1  HiG_2  HiG_3  HiG_4  HiG_5  HiG_6  HiG_9 HiG_17 HiG_18 HiG_19 HiG_16  HiG_8
#      1      2      3      4      5      6      7     11     12     13     17     19

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


# FIXES:

# What needs to be done: Re-run 11.1.1 with the fix to regenerate correct_n_edges_HiG_STRING2.14.0.rds with correct values, then re-run all scripts that load it (11.1.2, 11.2.0, 11.2.1, 11.3, 11.6, Fig3). The CTS-level analyses don't need to be touched.