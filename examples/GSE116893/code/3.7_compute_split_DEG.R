########## BEGINNING OF USER INPUT ##########
wd <- "C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/"

rebuild_deg <- FALSE  # set TRUE to rerun findMarkers (~5 min); FALSE loads saved markers
########## END OF USER INPUT ##########

# Computes DEG[["9_hi"]] / DEG[["9_lo"]] and adds HiG_9_hi / HiG_9_lo to the
# notsimplified graph list, following 11.1_STRINGweighted_CTS_network.R exactly.
# Outputs:
#   data/DEG_perState_min.prop0.25_lfc0.6_FDR0.05.rds       (adds 9_hi / 9_lo entries)
#   data/markers.up_C9split_ttest_min.prop0.25.rds           (saved for 3.8 reuse)
#   results/GSE116893_STRING_graph_perState_notsimplified.rds (adds HiG_9_hi / HiG_9_lo)
# Run 3.8 afterward to simplify + add specificity layer → simplified_combinedweighted.rds

library(SingleCellExperiment)
library(scran)
library(dplyr)
library(igraph)
library(STRINGdb)

data_dir  <- paste0(wd, "data/")
res_path  <- paste0(wd, "results/")
markers_path <- paste0(data_dir, "markers.up_C9split_ttest_min.prop0.25.rds")

# ---- load ----
sce <- readRDS(paste0(data_dir, "sce_C9split.rds"))
cat("leiden_C9split groups:", paste(sort(unique(colData(sce)[["leiden_C9split"]])), collapse=", "), "\n")

DEG        <- readRDS(paste0(data_dir, "DEG_perState_min.prop0.25_lfc0.6_FDR0.05.rds"))
graph_list <- readRDS(paste0(res_path, "GSE116893_STRING_graph_perState_notsimplified.rds"))
cat("Existing graph_list entries:", paste(names(graph_list), collapse=", "), "\n")

logFC.cut <- 0.6

# ---- STRING db (same settings as 11.1) ----
string_db <- STRINGdb$new(
    version         = "12.0",
    species         = 9606,
    score_threshold = 200,
    network_type    = "full",
    input_directory = paste0(data_dir, "PPIN")
)

# ---- findMarkers on leiden_C9split ----
if (rebuild_deg || !file.exists(markers_path)) {
    cat("\nRunning findMarkers on leiden_C9split...\n")
    markers_split <- findMarkers(
        sce,
        test      = "t",
        groups    = colData(sce)[["leiden_C9split"]],
        min.prop  = 0.25,
        direction = "up"
    )
    saveRDS(markers_split, markers_path)
    cat("Saved markers to:", markers_path, "\n")

    for (target in c("9_hi", "9_lo")) {
        DEG[[target]] <- subset(markers_split[[target]], summary.logFC > logFC.cut & FDR < 0.01) %>% rownames()
        cat(target, "DEG n =", length(DEG[[target]]), "\n")
    }

    deg_path <- paste0(data_dir, "DEG_perState_min.prop0.25_lfc0.6_FDR0.05.rds")
    saveRDS(DEG, file = deg_path)
    cat("Saved DEG to:", deg_path, "\n")
} else {
    cat("\nLoading saved markers from:", markers_path, "\n")
    markers_split <- readRDS(markers_path)
    cat("9_hi DEG n =", length(DEG[["9_hi"]]), " | 9_lo DEG n =", length(DEG[["9_lo"]]), "\n")
}

# ---- build raw (unsimplified) HiG subgraphs and add to notsimplified list ----
for (target in c("9_hi", "9_lo")) {
    cat("\n--- Processing", target, "---\n")

    interesting.up <- markers_split[[target]]
    diff_exp <- as.data.frame(interesting.up)
    diff_exp$symbol <- rownames(diff_exp)
    diff_exp <- subset(diff_exp, summary.logFC > logFC.cut & FDR < 0.01)

    mapped <- string_db$map(diff_exp, "symbol", removeUnmappedRows = TRUE)
    hits   <- mapped$STRING_id

    graph  <- string_db$get_subnetwork(hits)

    if (vcount(graph) < length(hits)) {
        missing_node <- setdiff(hits, V(graph)$name)
        if (length(missing_node) > 0) mapped <- mapped[-which(mapped$STRING_id %in% missing_node), ]
    }

    V(graph)$name   <- mapped[match(V(graph)$name, mapped$STRING_id), ]$symbol
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$summary.logFC
    V(graph)$FDR    <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$FDR
    E(graph)$weight <- E(graph)$combined_score / 1000
    graph           <- delete_edge_attr(graph, "combined_score")
    # NOTE: do NOT simplify here — that happens in 3.8 alongside specificity computation

    if (any(duplicated(V(graph)$name))) {
        dup <- unique(V(graph)$name[duplicated(V(graph)$name)])
        warning("Duplicate vertex names in HiG_", target, ": ", paste(dup, collapse=", "),
                " — resolve manually via biomaRt as in 11.1.1 before running 3.8")
    }

    graph_list[[paste0("HiG_", target)]] <- graph
    cat("HiG_", target, "vcount:", vcount(graph), " ecount:", ecount(graph),
        "| multi-edges:", any(which_multiple(graph)), "\n")
}

# ---- save to notsimplified ----
ns_path <- paste0(res_path, "GSE116893_STRING_graph_perState_notsimplified.rds")
saveRDS(graph_list, file = ns_path)
cat("\nSaved notsimplified graph_list to:", ns_path, "\n")
cat("Entries now:", paste(names(graph_list), collapse=", "), "\n")
cat("\nNext step: run 3.8_compute_split_specificity.R\n")
