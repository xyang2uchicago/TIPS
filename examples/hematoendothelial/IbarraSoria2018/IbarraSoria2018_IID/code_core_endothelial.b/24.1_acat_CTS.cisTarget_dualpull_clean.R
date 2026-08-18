# Set True if running 24.0 the first time
rebuild_mat <- TRUE
source(here::here("examples", "hematoendothelial", "IbarraSoria2018", "IbarraSoria2018_IID", "code_core_endothelial.b", "24.0_acat_load_input_clean.R"))

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

## check the loaded objects =========================
seed_TF  # "GATA2" "TAL1"  -- from CTS_endothelial.b (code 12.0)
CTS_ID   # "endothelial.b"
CTS_name # "CTS_endothelial.b"
CP_cluster # "endothelial.b"
CM_cluster # "blood"
CF_cluster # "endothelial.c"
celltype_col # "subcelltype"

dim(mat) # [1] 33  4
colnames(mat)
# [1] "CP_hi"  "CM_hi"  "CF_hi"  "CTS_endothelial.b"

## set subfolder =======================
(updir <- getwd())
dir.create(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)),
    showWarnings = FALSE, recursive = TRUE
)
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

NES_threshold <- 3

########################################################
##  input 6 -- RcisTarget
library(RcisTarget)
packageVersion("RcisTarget") # '1.29.0'
library(data.table)

data(motifAnnotations_hgnc)
motifAnnot <- motifAnnotations
dim(motifAnnot) # [1] 253096      8

cisTarget.res <- readRDS(file = "../cisTarget_targets_in_all_CTS.rds")
write.table(subset(cisTarget.res, NES >= NES_threshold & geneSet %in% c(CP_cluster)),
    paste0("../cisTarget_targets_in_", CP_cluster, "_NES", NES_threshold, ".txt"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

########################################################
##  binary flag CTS genes -- load or build
dim(mat) # [1] 37  4
files <- list.files(pattern = "^heatmap_blocked_CTS_") %>% grep("_v3.tsv", ., value = TRUE)
if (length(files) > 0) {
    mat <- read.table(grep("_v3.tsv", files, value = TRUE), sep = "\t", header = T, check.names = FALSE)

    saved_variables <- readRDS("cisTarget_variables.rds")

    x              <- saved_variables$x
    key_TFs        <- saved_variables$key_TFs
    motifAnnot_sub <- saved_variables$motifAnnot_sub
} else {
    motifAnnot_sub <- get_regulators_from_motifs(cisTarget.res, CTS_ID, NES_threshold, motifAnnot = motifAnnot)
    keys <- unique(motifAnnot_sub$regulators)
    if (any(is.na(keys))) keys <- keys[!is.na(keys)]
    keys
    # [1] "PRDM9"
    # [2] "CBFB"
    # [3] "ZNF582"
    # [4] "CRX;DBP;EP300;GATA1;GATA2;GATA3;GATA4;GATA5;GATA6;GTF2IRD1;HNF1A;HOXA13;IKZF2;MYB;MZF1;NFIA;NFIC;OTX1;OTX2;POU2F1;POU5F1;RAX;ZEB1"
    # [5] "ATF2;IKZF1"
    # ...
    # [11] "HMGA1;HMGA2"
    # ...
    # [33] "ELK1;ERF;ERG;ETV2;ETV5;ETV7;FLI1;SPDEF"
    # [34] "PRDM6"
    # ...
    # [45] "RBFOX2"
    # length(keys) # 53

    for (key in keys) {
        motif_TF_highConf <- motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf
        tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold &
                      (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf))
        genes <- unique(unlist(strsplit(tmp$enrichedGenes, ";")))
        mat[, paste0("cisTarget_", key, ".motif_target")] <- ifelse(rownames(mat) %in% genes, "1", "0")
    }

    dim(mat) # [1] 37 57

    x_all <- colnames(mat)[grepl("cisTarget_", colnames(mat))]
    x_all <- gsub("cisTarget_", "", x_all) %>% gsub("\\.motif_target", "", .)
    x_all <- unlist(strsplit(x_all, ";")) %>% unique()

    key_TFs <- seed_TF
    for (i in x_all) {
        for (j in c(CM_cluster, CF_cluster, CP_cluster)) {
            if (i %in% DEG[[j]] & i %in% rownames(sce)) {
                cat(paste0(i, " is in DEG of ", j, "\n"))
                key_TFs <- c(key_TFs, i)
            }
        }
    }
    # GATA1 is in DEG of blood
    # GATA2 is in DEG of endothelial.c
    # GATA2 is in DEG of endothelial.b
    # GFI1B is in DEG of blood
    # TAL1 is in DEG of blood
    # TAL1 is in DEG of endothelial.c
    # TAL1 is in DEG of endothelial.b
    # ETV2 is in DEG of endothelial.c
    # ETV2 is in DEG of endothelial.b
    # FLI1 is in DEG of endothelial.c
    # FLI1 is in DEG of endothelial.b
    # HOXB2 is in DEG of endothelial.c
    # HOXB2 is in DEG of endothelial.b
    # HMGA2 is in DEG of endothelial.c
    # HMGA2 is in DEG of endothelial.b
    # ELK3 is in DEG of endothelial.c
    # ELK3 is in DEG of endothelial.b
    # ETS2 is in DEG of endothelial.c
    # ETS2 is in DEG of endothelial.b
    # SOX17 is in DEG of endothelial.c
    # SOX17 is in DEG of endothelial.b
    # SOX4 is in DEG of endothelial.c
    # SOX4 is in DEG of endothelial.b
    # FEV is in DEG of endothelial.b
    # ATF3 is in DEG of endothelial.c
    # ATF4 is in DEG of blood
    # MYC is in DEG of blood

    for (i in x_all) {
        if (i %in% CTS[[CTS_ID]]) {
            cat(paste0(i, " is in CTS of ", CTS_ID, "\n"))
            key_TFs <- c(key_TFs, i)
        }
    }
    # GATA2 is in CTS of endothelial.b
    # TAL1 is in CTS of endothelial.b
    # ETV2 is in CTS of endothelial.b
    # SOX17 is in CTS of endothelial.b
    # FEV is in CTS of endothelial.b

    key_TFs <- unique(key_TFs)
    key_TFs
    #  [1] "GATA2" "TAL1"  "GATA1" "GFI1B" "ETV2"  "FLI1"  "HOXB2" "HMGA2" "ELK3"  "ETS2"  "SOX17" "SOX4"  "FEV"   "ATF3"  "ATF4"  "MYC"

    mat <- as.data.frame(mat)

    ## per-key cisTarget column lookup (0 if no matching column)
    x <- NULL
    for (j in key_TFs) {
        y <- intersect(
            which(grepl("cisTarget_", colnames(mat))),
            which(Reduce("|", lapply(j, function(p) grepl(p, colnames(mat), fixed = FALSE))))
        )
        cat(j, "\t", y, "\t", if (length(y) > 0) colnames(mat)[y] else "NO MATCH", "\n")
        if (length(y) == 0) y <- 0
        x <- c(x, y[1])
    }
    names(x) <- key_TFs
    # GATA2  -> col 14 (GATA family; CTS + DEG member)
    # TAL1   -> col 14 (same col as GATA2)
    # GATA1  -> col 14 (same col as GATA2)
    # GFI1B  -> col 14 (same col as GATA family)
    # ETV2   -> col 22 (DLX2;DLX3;ERF;ETV2;... motif_target; CTS member)
    # FLI1   -> col 22 (same col as ETV2)
    # HOXB2  -> col 22 (same col as ETV2)
    # HMGA2  -> col 32 (HMGA1;HMGA2 motif_target)
    # ELK3   -> col 33 (ELF/ELK/SOX broad family)
    # ETS2   -> col 33 (same col as ELK3)
    # SOX17  -> col 33 (same col as ELK3; CTS member)
    # SOX4   -> col 33 (same col as ELK3)
    # FEV    -> col 36 (broad ETS; also specific col 45, 60; CTS member)
    # ATF3   -> col 46 (broad bHLH motif)
    # ATF4   -> col 46 (same col as ATF3)
    # MYC    -> col 46 (same col as ATF3)

    ## MANUAL SELECTION: one per cisTarget col group; prefer CTS members (ETV2/col22, GATA2/col14, HMGA2/col32, SOX17/col33, FEV/col36)
    key_TFs <- c("ETV2", "GATA2", "HMGA2", "SOX17", "FEV")
    x <- x[key_TFs]

    if (length(x) > 0) {
        for (j in seq_along(x)) {
            key <- names(x)[j]
            mat[, paste0(key, "_CP_candidate")] <- ifelse(mat[, "CP_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CM_candidate")] <- ifelse(mat[, "CM_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CF_candidate")] <- ifelse(mat[, "CF_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
        }
    }

    cat(paste0("key_TFs: ", paste(key_TFs, collapse = "_"), "\n"))
    # key_TFs: ETV2_GATA2_HMGA2_SOX17_FEV

    if (length(key_TFs) > 0) {
        fileName <- paste0("heatmap_blocked_", CTS_name, "_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv")
        write.table(mat, file = fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
        saveRDS(list(x = x, key_TFs = key_TFs, motifAnnot_sub = motifAnnot_sub), "cisTarget_variables.rds")
        cat("Saved heatmap to:", fileName, "\n")
    } else {
        stop("No key TFs found for ", CTS_name, "\n")
    }
}

motif_TF_highConf <- gsub("cisTarget_", "", colnames(mat)[x]) %>% gsub("\\.motif_target", "", .)
print(motif_TF_highConf)

###  heatmaps confirming key TF's self impact ================
for (key in motif_TF_highConf) {
    if (grepl(";", key)) {
        key_in_TFfamily <- strsplit(key, ";", fixed = TRUE) %>%
            unlist() %>% intersect(key_TFs) %>% unique()
    } else {
        key_in_TFfamily <- key
    }

    p <- tryCatch(
        heatmap_pull_candidate(mat, graph_list, CTS_name, CHD,
            key = key_in_TFfamily, coding_genes = coding_genes, TF = TF_mouse,
            show_SMC_access = FALSE
        ),
        error = function(e) {
            message("heatmap skipped for '", key_in_TFfamily, "' (0 candidate genes): ", e$message)
            NULL
        }
    )
    if (!is.null(p)) {
        pdf(file = paste0("heatmap_blocked_", CTS_name, "_cisTarget_", key_in_TFfamily, "_v3_coding_target.pdf"), height = 4)
        print(p)
        dev.off()
        cat("  => saved heatmap PDF for", key_in_TFfamily, "\n")
    }
}

##########################################################
## identify_TF_targeted_pull_candidate
## Redefined locally to use HiG_cardiac.a / CTS_cardiac.a naming (IID convention)
library(BioNet)
packageVersion("BioNet") # '1.56.0'
library(igraph)
library(tibble)

mat <- read.table(paste0("heatmap_blocked_", CTS_name, "_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv"),
    sep = "\t", header = TRUE, check.names = FALSE)
dim(mat)

# No ATAC data: dummy access columns so downstream filtering passes everything
for (col in c("PCW6CP_access", "PCW8_CM_access", "PCW19_CM_access",
              "PCW8_CF_access", "PCW19_CF_access",
              "PCW8_SMC_access", "PCW19_SMC_access", "iEPC_access")) {
    mat[[col]] <- 1L
}

identify_TF_targeted_pull_candidate <- function(
    mat, graph_list, CTS_name, CHD,
    key, keep_selfloop = TRUE,
    TF_bound_column_name,
    TF_appendix = NULL,
    edge_colored_by_Maven2023_ISL1KO = FALSE,
    key_in_TFfamily = key
) {
    CTS_ID <- sub("^CTS_", "", CTS_name)
    four_names <- c(
        "CTSHiG.CP_TF.target",
        "CTS.CP_TF.target_HiGCM",
        "CTS.CP_TF.target_HiGCF",
        "CTSHiG.CP_TF.target_CPopen"
    )
    graph_TF_list <- list(
        graph_list[[paste0("HiG_", CTS_ID)]],
        graph_list[[paste0("CTS_", CTS_ID)]],
        graph_list[[paste0("CTS_", CTS_ID)]],
        graph_list[[paste0("CTS_", CTS_ID)]]
    )
    names(graph_TF_list) <- four_names
    if (any(vapply(graph_TF_list, is.null, logical(1))))
        stop("Missing graph(s): ", paste(names(graph_TF_list)[vapply(graph_TF_list, is.null, logical(1))], collapse = ", "))

    TF_target <- rownames(mat)[which(mat[, TF_bound_column_name] == 1)]

    if (is.null(TF_appendix)) {
        GRN_node_list <- list(
            rownames(mat)[which(mat[, "CP_candidate"] == 1)],
            rownames(mat)[which(mat[, "CM_candidate"] == 1)],
            rownames(mat)[which(mat[, "CF_candidate"] == 1)],
            rownames(mat)[which(mat[, "CP_candidate"] == 1)]
        )
    } else {
        GRN_node_list <- list(
            CTSHiG.CP_TF.target        = rownames(mat)[which(mat[, paste0(key, "_CP_candidate")] == 1)],
            CTS.CP_TF.target_HiGCM    = rownames(mat)[which(mat[, paste0(key, "_CM_candidate")] == 1)],
            CTS.CP_TF.target_HiGCF    = rownames(mat)[which(mat[, paste0(key, "_CF_candidate")] == 1)],
            CTSHiG.CP_TF.target_CPopen = rownames(mat)[which(mat[, paste0(key, "_CP_candidate")] == 1)]
        )
    }
    names(GRN_node_list) <- four_names

    for (i in 2:length(GRN_node_list)) {
        if (grepl("CM", names(GRN_node_list)[i])) {
            GRN_node_list[[i]] <- intersect(GRN_node_list[[i]],
                rownames(mat)[which(mat[, "PCW6CP_access"] == 1 | mat[, "PCW8_CM_access"] == 1 | mat[, "PCW19_CM_access"] == 1)])
        } else if (grepl("CF", names(GRN_node_list)[i]) | grepl("SMC", names(GRN_node_list)[i])) {
            GRN_node_list[[i]] <- intersect(GRN_node_list[[i]],
                rownames(mat)[which(
                    mat[, "PCW6CP_access"] == 1 | mat[, "PCW8_CF_access"] == 1 | mat[, "PCW19_CF_access"] == 1 |
                    mat[, "PCW8_SMC_access"] == 1 | mat[, "PCW19_SMC_access"] == 1 | mat[, "iEPC_access"] == 1)])
        } else {
            GRN_node_list[[i]] <- intersect(GRN_node_list[[i]],
                rownames(mat)[which(mat[, "PCW6CP_access"] == 1)])
        }
    }

    mat <- as.data.frame(apply(mat, 2, as.numeric), row.names = rownames(mat))

    for (name in names(graph_TF_list)) {
        graph <- graph_TF_list[[name]]
        V(graph)$name <- toupper(V(graph)$name)
        graph <- induced_subgraph(graph, intersect(GRN_node_list[[name]], V(graph)$name))
        V(graph)$color <- ifelse(V(graph)$name %in% CHD, "red", "lightgrey")
        E(graph)$color <- "grey70"
        E(graph)$lty   <- "solid"

        TF_target_in_graph <- intersect(TF_target, V(graph)$name)
        if (!keep_selfloop) TF_target_in_graph <- setdiff(TF_target_in_graph, key_in_TFfamily)

        if (length(TF_target_in_graph) > 0) {
            g_check <- if (is_directed(graph)) as.undirected(graph, mode = "collapse") else graph
            flag  <- !(key_in_TFfamily %in% V(graph)$name)
            g_check <- add_vertex_if_missing(g_check, key_in_TFfamily)
            if (is.null(vertex_attr(g_check, "color"))) V(g_check)$color <- rep("lightgrey", vcount(g_check))
            if (flag) V(g_check)[name == key_in_TFfamily]$color <- "yellow"

            vp_pairs <- as.vector(rbind(key_in_TFfamily, TF_target_in_graph))
            eid      <- get_edge_ids(g_check, vp = vp_pairs)
            targets_to_add <- TF_target_in_graph[eid == 0]

            if (length(targets_to_add) > 0) {
                m_before   <- ecount(g_check)
                edges_vec  <- as.vector(rbind(key_in_TFfamily, targets_to_add))
                g_check    <- add_edges(g_check, edges_vec)
                new_edges  <- seq.int(m_before + 1, ecount(g_check))
                if (is.null(edge_attr(g_check, "color"))) E(g_check)$color <- rep("grey70", ecount(g_check))
                if (is.null(edge_attr(g_check, "lty")))   E(g_check)$lty   <- rep("solid",  ecount(g_check))
                E(g_check)$color[new_edges] <- "grey70"
                E(g_check)$lty[new_edges]   <- "dashed"
            }
            graph <- g_check
        }
        graph_TF_list[[name]] <- graph
    }
    return(graph_TF_list)
}

for (key in key_TFs) {
    key_column <- which(grepl(key, colnames(mat)) & grepl("cisTarget_", colnames(mat)))
    key_in_TFfamily <- key
    cat("Processing TF:", key, "| cisTarget column index:", key_column, "\n")

    graph_TF_list <- identify_TF_targeted_pull_candidate(mat, graph_list, CTS_name, CHD,
        key = key, keep_selfloop = TRUE,
        TF_bound_column_name = key_column, TF_appendix = key,
        edge_colored_by_Maven2023_ISL1KO = FALSE, key_in_TFfamily = key_in_TFfamily)
    saveRDS(graph_TF_list, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
    cat("  saved:", paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"), "\n")
}

for (key in key_TFs) {
    key_in_TFfamily <- key
    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
    plot_TF_targeted_pull_candidate(graph_TF_list, key_in_TFfamily, CTS_name, saveFigure = TRUE)
    cat("  plotted:", key_in_TFfamily, "\n")
}
# => PPI_graph_<key>_GRN_prediction_CTS_endothelial.b_v3.pdf

##########################################################
make_layout <- function(g, seed) { set.seed(seed); layout_with_fr(g) }

for (key in key_TFs) {
    key_in_TFfamily <- key
    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))

    res <- fill_TF_targeting_predicted_edges(graph_TF_list,
        linkage_name = "CM", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CM_cluster, TF_symbol = key_in_TFfamily,
        HVG = rownames(sce))
    cat(key, "CM vcount:", vcount(res[["g_CT_sub"]]), "\n")
    if (vcount(res[["g_CT_sub"]]) > 0)
        saveRDS(res, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CM_final.rds"))

    res_cf <- fill_TF_targeting_predicted_edges(graph_TF_list,
        linkage_name = "CF", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
        HVG = rownames(sce))
    cat(key, "CF vcount:", vcount(res_cf[["g_CT_sub"]]), "\n")
    if (vcount(res_cf[["g_CT_sub"]]) > 0)
        saveRDS(res_cf, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CF_final.rds"))
}
# ETV2  CM vcount: 0  (not saved)
# ETV2  CF vcount: 0  (not saved)
# GATA2 CM vcount: 0  (not saved)
# GATA2 CF vcount: 1
# HMGA2 CM vcount: 2
# HMGA2 CF vcount: 6
# SOX17 CM vcount: 0  (not saved)
# SOX17 CF vcount: 2
# FEV   CM vcount: 2
# FEV   CF vcount: 4

### reporting ==========
print((files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds")))

final_table <- NULL
for (f in files) {
    pull <- ifelse(grepl("CF", f), "CF", "CM")
    pattern <- paste(key_TFs, collapse = "|")
    key <- key_in_TFfamily <- regmatches(f, regexpr(pattern, f))
    cat("Processing:", f, "| pull:", pull, "| key:", key, "\n")

    res <- readRDS(file = f)
    g1  <- res[["g_CT_sub"]]
    g2  <- res[["g_descendant_sub"]]

    if (vcount(g1) > 0 & vcount(g2) > 0 & vcount(g1) == vcount(g2)) {
        change_df <- edge_change_table(g1 = g1, g2 = g2, weight_attr = "weight", missing_as = 0, undirected = TRUE)
        predict   <- prioritize_edge_change(g1, edge_change_df = change_df, top_n = 5,
                                            title = paste0(pull, "-pull subnetwork_", key_in_TFfamily))
    } else {
        cat("No edges in", pull, "pull for", key_in_TFfamily, "\n"); next
    }

    y <- motifAnnot_sub[which(motifAnnot_sub$regulators %in% key | grepl(key, motifAnnot_sub$regulators)), ]
    motif_TF_highConf_val <- y$motif_TF_highConf
    tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold &
                      (motif == motif_TF_highConf_val | TF_highConf == motif_TF_highConf_val))
    tmp1 <- if (nrow(tmp) > 0) tmp[1, , drop = FALSE] else data.frame(TF_highConf = NA, motif = NA, NES = NA)
    change_df <- cbind(linkage = pull, change_df,
                       TF_highConf = tmp1$TF_highConf, motif = tmp1$motif, NES = tmp1$NES)
    change_df$TF_highConf[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
    change_df$motif[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)]       <- ""
    change_df$NES[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)]         <- ""
    final_table <- rbind(final_table, change_df)
}

dim(final_table) # 13  13
print(final_table)

write.table(final_table,
    file = paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)
