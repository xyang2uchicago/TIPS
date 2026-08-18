# Set True if running 24.0 the first time
rebuild_mat <- TRUE
source(here::here("examples", "cardiac", "GSE87038", "GSE87038_IID", "code_core", "24.0_acat_load_input_clean.R"))

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

## check the loaded objects =========================
seed_TF  # "NR2F1" "IRX5" "ALX1"  -- from CTS_8 (no HiGCTS_8 in IID)
CTS_ID   # "8"
CTS_name # "CTS_8"
CP_cluster # "8"
CM_cluster # "17"
CF_cluster # "18"
celltype_col # "label"

dim(mat) # [1] 54  4
colnames(mat)
# [1] "CP_hi"  "CM_hi"  "CF_hi"  "CTS_8"

## set subfolder =======================
(updir <- getwd())
dir.create(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)),
    showWarnings = FALSE, recursive = TRUE)
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

NES_threshold <- 3

########################################################
##  input 6 -- RcisTarget
library(RcisTarget)
library(data.table)
data(motifAnnotations_hgnc)
motifAnnot <- motifAnnotations
dim(motifAnnot) # [1] 253096      8

cisTarget.res <- readRDS(file = "../cisTarget_targets_in_all_CTS.rds")
write.table(subset(cisTarget.res, NES >= NES_threshold & geneSet %in% c(CP_cluster)),
    paste0("../cisTarget_targets_in_", CP_cluster, "_NES", NES_threshold, ".txt"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

########################################################
##  binary flag CTS genes -- load or build
files <- list.files(pattern = "^heatmap_blocked_CTS_") %>% grep("_v3.tsv", ., value = TRUE)
if (length(files) > 0) {
    mat <- read.table(grep("_v3.tsv", files, value = TRUE), sep = "\t", header = T, check.names = FALSE)
    saved_variables <- readRDS("cisTarget_variables.rds")
    x            <- saved_variables$x
    key_TFs      <- saved_variables$key_TFs
    motifAnnot_sub <- saved_variables$motifAnnot_sub
} else {
    motifAnnot_sub <- get_regulators_from_motifs(cisTarget.res, CTS_ID, NES_threshold, motifAnnot = motifAnnot)
    keys <- unique(motifAnnot_sub$regulators)
    if (any(is.na(keys))) keys <- keys[!is.na(keys)]
    keys
    #  [1] "PARP1"
    #  [2] "HOXD3"
    #  [3] "ZSCAN9"
    #  [4] "CEBPB;PDX1;STAT6"
    #  [5] "TCF3"
    #  [6] "POU3F3"
    #  [7] "AHCTF1"
    #  [8] "ISX;LHX2;LHX9;SHOX2"
    #  [9] "FOXM1;IKZF1;TBL1XR1"
    # [10] "FOXM1"
    # [11] "HMGA1;HMGA2"
    # [12] "STAT5A"
    # [13] "POU4F3"
    # [14] "RREB1"
    # [15] "PRDM5;ZNF324;ZNF324B;ZNF341;ZNF580"
    # [16] "SF1"
    # [17] "HMG20A"
    # [18] "KLF2"
    # [19] "BCL11A;BCL6;CHURC1;DDX20;DPF2;EGR1;EGR2;EP300;FLI1;GTF2F1;HINFP;IKZF1;IKZF3;IRF4;KLF10;KLF13;KLF15;KLF16;KLF6;KLF9;MAZ;MTF1;MZF1;NFXL1;OVOL2;PATZ1;PAX4;PML;PRDM9;PURA;RARA;RARB;RARG;RBAK;RCOR1;RELA;RXRA;RXRB;RXRG;SIN3A;SP1;SP2;SP3;SP4;SP5;SPI1;SREBF1;SREBF2;STAT5A;TAF1;THRA;THRB;VEZF1;WRNIP1;WT1;ZBED1;ZBTB14;ZBTB17;ZBTB5;ZBTB7B;ZNF148;ZNF202;ZNF212;ZNF214;ZNF263;ZNF341;ZNF398;ZNF432;ZNF444;ZNF467;ZNF496;ZNF529;ZNF543;ZNF596;ZNF615;ZNF701;ZNF737;ZNF740;ZNF875;ZNF880"
    # ...

    for (key in keys) {
        motif_TF_highConf <- motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf
        tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold &
                      (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf))
        genes <- unique(unlist(strsplit(tmp$enrichedGenes, ";")))
        mat[, paste0("cisTarget_", key, ".motif_target")] <- ifelse(rownames(mat) %in% genes, "1", "0")
    }

    dim(mat) # [1] 54 37 (4 base + 33 cisTarget columns)
    colnames(mat)

    x_all <- colnames(mat)[grepl("cisTarget_", colnames(mat))]
    x_all <- gsub("cisTarget_", "", x_all) %>% gsub("\\.motif_target", "", .)
    x_all <- unlist(strsplit(x_all, ";")) %>% unique()

    key_TFs <- seed_TF
    for (i in x_all) {
        for (j in c(CM_cluster, CF_cluster, CP_cluster)) {
            if (i %in% DEG[[j]] & i %in% rownames(sce)) {
                cat(i, "is in DEG of", j, "\n")
                key_TFs <- c(key_TFs, i)
            }
        }
    }
    # TCF3 is in DEG of 8
    # HMGA2 is in DEG of 17
    # HMGA2 is in DEG of 18
    # HMGA2 is in DEG of 8
    # KLF6 is in DEG of 17
    # KLF6 is in DEG of 18

    for (i in x_all) {
        if (i %in% CTS[[CTS_ID]]) {
            cat(i, "is in CTS of", CTS_ID, "\n")
            key_TFs <- c(key_TFs, i)
        }
    }
    # RARB is in CTS of 8

    key_TFs <- unique(key_TFs)
    key_TFs # [1] "NR2F1" "IRX5"  "ALX1"  "TCF3"  "HMGA2" "KLF6"  "RARB"

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
    # NR2F1  -> 0 (no cisTarget column)
    # IRX5   -> 0
    # ALX1   -> 0
    # TCF3   -> col for TCF3.motif_target
    # HMGA2  -> col for HMGA1;HMGA2.motif_target
    # KLF6   -> col for big BCL11A;...;KLF6;... motif_target
    # RARB   -> same col as KLF6 (same big motif)

    ## MANUAL SELECTION: drop NR2F1/IRX5/ALX1 (no cisTarget column)
    ## keep one from KLF6/RARB family (both map to same motif — keep both)
    key_TFs <- c("TCF3", "HMGA2", "KLF6", "RARB")
    x <- x[key_TFs]
    # TCF3 HMGA2 KLF6 RARB

    if (length(x) > 0) {
        for (j in seq_along(x)) {
            key <- names(x)[j]
            mat[, paste0(key, "_CP_candidate")] <- ifelse(mat[, "CP_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CM_candidate")] <- ifelse(mat[, "CM_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CF_candidate")] <- ifelse(mat[, "CF_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
        }
    }

    cat(paste0("key_TFs: ", paste(key_TFs, collapse = "_"), "\n"))
    # key_TFs: TCF3_HMGA2_KLF6_RARB

    fileName <- paste0("heatmap_blocked_", CTS_name, "_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv")
    write.table(mat, file = fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    saveRDS(list(x = x, key_TFs = key_TFs, motifAnnot_sub = motifAnnot_sub), "cisTarget_variables.rds")
    cat("Saved heatmap to:", fileName, "\n")
}

motif_TF_highConf <- gsub("cisTarget_", "", colnames(mat)[x]) %>% gsub("\\.motif_target", "", .)
print(motif_TF_highConf)
# [1] "TCF3"
# [2] "HMGA1;HMGA2"
# [3] "BCL11A;BCL6;...;KLF6;...;RARB;...ZNF880"
# [4] same big motif

###  heatmaps
for (key in motif_TF_highConf) {
    if (grepl(";", key)) {
        key_in_TFfamily <- strsplit(key, ";", fixed = TRUE) %>%
            unlist() %>% intersect(key_TFs) %>% unique()
    } else {
        key_in_TFfamily <- key
    }
    p <- tryCatch(
        heatmap_pull_candidate(mat, graph_list, CTS_name, CHD,
            key = key_in_TFfamily, coding_genes = coding_genes, TF = TF_human,
            show_SMC_access = FALSE),
        error = function(e) { message("heatmap skipped for '", key_in_TFfamily, "': ", e$message); NULL }
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
## Redefined locally to use HiG_8/CTS_8 naming (IID graph lacks legacy aliases)
library(BioNet)
library(igraph)
library(tibble)

mat <- read.table(paste0("heatmap_blocked_", CTS_name, "_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv"),
    sep = "\t", header = TRUE, check.names = FALSE)
dim(mat) # [1] 54 ... (4 + cisTarget + candidate columns)

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

##########################################################
make_layout <- function(g, seed) { set.seed(seed); layout_with_fr(g) }

for (key in key_TFs) {
    key_in_TFfamily <- key
    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))

    res <- fill_TF_targeting_predicted_edges(graph_TF_list, linkage_name = "CM", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CM_cluster, TF_symbol = key_in_TFfamily, HVG = rownames(sce))
    cat(key, "CM vcount:", vcount(res[["g_CT_sub"]]), "\n")
    if (vcount(res[["g_CT_sub"]]) > 0)
        saveRDS(res, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CM_final.rds"))

    res_cf <- fill_TF_targeting_predicted_edges(graph_TF_list, linkage_name = "CF", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily, HVG = rownames(sce))
    cat(key, "CF vcount:", vcount(res_cf[["g_CT_sub"]]), "\n")
    if (vcount(res_cf[["g_CT_sub"]]) > 0)
        saveRDS(res_cf, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CF_final.rds"))
}
# TCF3  CM vcount: 0   (not saved)
# TCF3  CF vcount: 0   (not saved)
# HMGA2 CM vcount: 4
# HMGA2 CF vcount: 2
# KLF6  CM vcount: 2
# KLF6  CF vcount: 0   (not saved)
# RARB  CM vcount: 2
# RARB  CF vcount: 0   (not saved)

### reporting ==========
print((files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds")))
# [1] "PPI_graph_HMGA2_GRN_prediction_CTS_8_CF_final.rds"
# [2] "PPI_graph_HMGA2_GRN_prediction_CTS_8_CM_final.rds"
# [3] "PPI_graph_KLF6_GRN_prediction_CTS_8_CM_final.rds"
# [4] "PPI_graph_RARB_GRN_prediction_CTS_8_CM_final.rds"

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

dim(final_table) # 6  13
print(final_table)

write.table(final_table,
    file = paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
cat("Saved final_table.tsv\n")
