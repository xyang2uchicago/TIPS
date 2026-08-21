code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_STRING", "code_core")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
source(file.path(code_dir, "24.0_acat_load_input_clean_CP.1.R"))

source(here::here("R", paste0("celltype_specific_weight_v", celltype_specific_weight_version, ".R")))

## check the loaded objects =========================
seed_TF  # character(0)
CTS_ID   # "CP.1"
CTS_name # "CTS_CP.1"

## set subfolder =======================
(updir <- getwd())
dir.create(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)),
    showWarnings = FALSE, recursive = TRUE
)
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

########################################################
##  input 6 -- RcisTarget predicted TF enriched among CTS genes
library(RcisTarget)
packageVersion("RcisTarget")
library(data.table)

data(motifAnnotations_hgnc)
motifAnnot <- motifAnnotations

cisTarget.res <- readRDS(file = "../cisTarget_targets_in_all_CTS.rds")
write.table(subset(cisTarget.res, NES >= NES_threshold & geneSet %in% c(CTS_ID)),
    paste0("../cisTarget_targets_in_CP.1_NES", NES_threshold, ".txt"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

########################################################
##   --- binary flag CTS genes
dim(mat)
### = load the binary annotation matrix built in 24.0 ============================
files <- list.files(pattern = "^heatmap_blocked_CTS_") %>% grep("_v3.tsv", ., value = TRUE)
if (length(files) > 0) {
    fileName <- grep("_v3.tsv", files, value = TRUE)
    mat <- read.table(fileName, sep = "\t", header = T, check.names = FALSE)

    saved_variables <- readRDS(file = "cisTarget_variables.rds")

    x <- saved_variables$x
    key_TFs <- saved_variables$key_TFs
    motifAnnot_sub <- saved_variables$motifAnnot_sub
} else {
    motifAnnot_sub <- get_regulators_from_motifs(cisTarget.res, CTS_ID, NES_threshold, motifAnnot = motifAnnot)
    keys <- unique(motifAnnot_sub$regulators)
    if (any(is.na(keys))) keys <- keys[!is.na(keys)]

    for (key in keys) {
        motif_TF_highConf <- motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf
        tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold & (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf))
        genes <- unique(unlist(strsplit(tmp$enrichedGenes, ";")))
        mat[, paste0("cisTarget_", key, ".motif_target")] <- ifelse(rownames(mat) %in% genes, "1", "0")
    }

    dim(mat)
    colnames(mat)

    x <- colnames(mat)[grepl("cisTarget_", colnames(mat))]
    x <- gsub("cisTarget_", "", x) %>% gsub("\\.motif_target", "", .)
    x <- unlist(strsplit(x, ";")) %>% unique()

    key_TFs <- seed_TF
    for (i in x) {
        for (j in c(CM_cluster, CF_cluster, CP_cluster)) {
            if (i %in% DEG[[j]] & i %in% rownames(sce)) {
                cat(paste0(i, " is in DEG of ", j, "\n"))
                key_TFs <- c(key_TFs, i)
            }
        }
    }
    for (i in x) {
        if (i %in% CTS[[CTS_ID]]) {
            cat(paste0(i, " is in CTS of ", CTS_ID, "\n"))
            key_TFs <- c(key_TFs, i)
        }
    }
    key_TFs <- unique(key_TFs)
    key_TFs
    # At NES>=4.5 for CTS_CP.1: EBF1 and BCL6B are enriched motif TFs but neither is in DEG of any cluster
    # → key_TFs = character(0) (matches original 24.1.2 behavior: no GRN predictions at NES=4.5)

    mat <- as.data.frame(mat)

    ## one column index per TF (take first match)
    x <- sapply(key_TFs, function(tf) {
        idx <- which(grepl("cisTarget_", colnames(mat)) & grepl(tf, colnames(mat), fixed = FALSE))
        if (length(idx) > 0) idx[1] else NA_integer_
    })
    x <- x[!is.na(x)]
    if (length(x) > 0) {
        print(x)
        print(colnames(mat)[x])
    }

    names(x) <- names(x)
    if (length(x) > 0) {
        for (j in seq_along(x)) {
            key <- names(x)[j]
            mat[, paste0(key, "_CP_candidate")] <- ifelse(mat[, "CP_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CM_candidate")] <- ifelse(mat[, "CM_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CF_candidate")] <- ifelse(mat[, "CF_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
        }
    }

    dim(mat)
    cat(paste0("key_TFs: ", paste(key_TFs, collapse = "_"), "\n"))
    # key_TFs: (empty at NES=4.5 — EBF1/BCL6B not in DEG of any cluster)

    if (length(key_TFs) > 0) {
        fileName <- paste0("heatmap_blocked_", CTS_name, "_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv")
        write.table(mat, file = fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
        saveRDS(list(x = x, key_TFs = key_TFs, motifAnnot_sub = motifAnnot_sub), "cisTarget_variables.rds")
    } else {
        stop("No key TFs found for ", CTS_name, " at NES=", NES_threshold,
             " — EBF1/BCL6B enriched at NES>=4.5 for CTS_CP.1 but not in DEG of any cluster (matches original 24.1.2)")
    }
}

motif_TF_highConf <- gsub("cisTarget_", "", colnames(mat)[x]) %>% gsub("\\.motif_target", "", .)
print(motif_TF_highConf)

###  heatmap confirming key TF self impact among CTS_CP.1 ================
for (key in motif_TF_highConf) {
    if (grepl(";", key)) {
        key_in_TFfamily <- strsplit(key, ";", fixed = T) %>%
            unlist() %>%
            intersect(key_TFs) %>%
            unique()
    } else {
        key_in_TFfamily <- key
    }

    p <- tryCatch(
        heatmap_pull_candidate(mat, graph_list, CTS_name, CHD,
            key = key_in_TFfamily, coding_genes = coding_genes, TF = TF_human,
            show_SMC_access = FALSE
        ),
        error = function(e) {
            message("heatmap skipped for '", key_in_TFfamily, "': ", e$message)
            NULL
        }
    )
    if (!is.null(p)) {
        pdf(file = paste0("heatmap_blocked_", CTS_name, "_cisTarget_", key_in_TFfamily, "_v3_coding_target.pdf"), height = 4)
        print(p)
        dev.off()
    }
}

##########################################################
library(BioNet)
library(igraph)
library(tibble)

mat <- read.table(paste0("heatmap_blocked_", CTS_name, "_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv"), sep = "\t", header = T, check.names = FALSE)
print(dim(mat))

# No ATAC data: add dummy access columns so downstream functions run without errors
for (col in c("PCW6CP_access", "PCW8_CM_access", "PCW19_CM_access",
              "PCW8_CF_access", "PCW19_CF_access",
              "PCW8_SMC_access", "PCW19_SMC_access", "iEPC_access")) {
    mat[[col]] <- 1L
}

for (key in key_TFs) {
    key_column <- which(grepl(key, colnames(mat)) & grepl("cisTarget_", colnames(mat)))
    key_in_TFfamily <- key

    graph_TF_list <- identify_TF_targeted_pull_candidate(mat, graph_list, CTS_name, CHD,
        key = key,
        keep_selfloop = TRUE,
        TF_bound_column_name = key_column,
        TF_appendix = key,
        edge_colored_by_Maven2023_ISL1KO = FALSE,
        key_in_TFfamily = key_in_TFfamily
    )
    saveRDS(graph_TF_list, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
}

names(graph_TF_list)
# [1] "CTSHiG.CP.1_TF.target"        "CTS.CP.1_TF.target_HiGCM"
# [3] "CTS.CP.1_TF.target_HiGCF"     "CTSHiG.CP.1_TF.target_CPopen"

for (key in key_TFs) {
    key_in_TFfamily <- key
    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
    tryCatch(
        plot_TF_targeted_pull_candidate(graph_TF_list, key_in_TFfamily, CTS_name, saveFigure = TRUE),
        error = function(e) message("plot skipped for '", key_in_TFfamily, "' (empty subgraph): ", e$message)
    )
}

##########################################################
make_layout <- function(g, seed) {
    set.seed(seed)
    layout_with_fr(g)
}

for (key in key_TFs) {
    key_in_TFfamily <- key

    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))

    res <- fill_TF_targeting_predicted_edges(graph_TF_list,
        linkage_name = "CM", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CM_cluster, TF_symbol = key_in_TFfamily,
        HVG = rownames(sce)
    )
    cat(paste0(key, ' CM vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), "\n"))
    # KLF6 CM: 25 | GATA3 CM: 26
    if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CM_final.rds"))

    res_cf <- fill_TF_targeting_predicted_edges(graph_TF_list,
        linkage_name = "CF", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
        HVG = rownames(sce)
    )
    cat(paste0(key, ' CF vcount(res_cf[["g_CT_sub"]]): ', vcount(res_cf[["g_CT_sub"]]), "\n"))
    # KLF6 CF: 31 | GATA3 CF: 32
    if (vcount(res_cf[["g_CT_sub"]]) > 0) saveRDS(res_cf, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CF_final.rds"))
}

### reporting ==========
print((files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds")))
# [1] "PPI_graph_GATA3_GRN_prediction_CTS_CP.1_CF_final.rds"
# [2] "PPI_graph_GATA3_GRN_prediction_CTS_CP.1_CM_final.rds"
# [3] "PPI_graph_KLF6_GRN_prediction_CTS_CP.1_CF_final.rds"
# [4] "PPI_graph_KLF6_GRN_prediction_CTS_CP.1_CM_final.rds"

final_table <- NULL
for (f in files) {
    pull <- ifelse(grepl("CF", f), "CF", "CM")
    pattern <- paste(key_TFs, collapse = "|")
    key <- key_in_TFfamily <- regmatches(f, regexpr(pattern, f))

    res <- readRDS(file = f)
    g1 <- res[["g_CT_sub"]]
    g2 <- res[["g_descendant_sub"]]

    if (vcount(g1) > 0 & vcount(g2) > 0 & vcount(g1) == vcount(g2)) {
        change_df <- edge_change_table(g1 = g1, g2 = g2, weight_attr = "weight", missing_as = 0, undirected = TRUE)
        predict <- prioritize_edge_change(g1, edge_change_df = change_df, top_n = 5, title = paste0(pull, "-pull subnetwork_", key_in_TFfamily))
    } else {
        cat("No edges in", pull, "-pull subnetwork for", key_in_TFfamily, "\n")
        next
    }

    x_key <- paste0(key, "_motif_target")
    y <- motifAnnot_sub[which(motifAnnot_sub$regulators %in% key | grepl(key, motifAnnot_sub$regulators)), ]
    motif_TF_highConf_val <- y$motif_TF_highConf
    tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold &
                      (motif == motif_TF_highConf_val | TF_highConf == motif_TF_highConf_val))

    tmp1 <- if (nrow(tmp) > 0) tmp[1, , drop = FALSE] else data.frame(TF_highConf = NA, motif = NA, NES = NA)
    change_df <- cbind(linkage = pull, change_df,
                       TF_highConf = tmp1$TF_highConf, motif = tmp1$motif, NES = tmp1$NES)
    change_df$TF_highConf[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
    change_df$motif[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
    change_df$NES[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""

    final_table <- rbind(final_table, change_df)
}

if (!is.null(final_table)) {
    dim(final_table)
    write.table(final_table,
        file = paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
        quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
    )
}
