# Set True if running 24.0 the first time
rebuild_mat <- TRUE
source(here::here("examples", "hematoendothelial", "GSE87038", "GSE87038_STRING", "code_core_13", "24.0_acat_load_input_clean.R"))

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

## check the loaded objects =========================
seed_TF # NULL  (set NULL in 24.0 USER INPUT; 12.0 gave ETV2/TAL1 for HiGCTS_13 in original code; key_TFs overridden manually below)
names(graph_list)
# [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"       "HiG_6"
# [7] "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"      "HiG_17"      "HiG_18"
# [13] "HiG_19"      "HiG_7"       "HiG_11"      "HiG_15"      "HiG_16"      "HiG_13"
# [19] "HiG_8"       "HiGCTS_7"    "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1"
# [25] "HiGCTS_13"   "HiGCTS_8"    "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"
# [31] "CTS_16.1"    "CTS_13"      "CTS_8"
# (33 networks: HiGCTS_13 present via lfc_HiGCTS=0.5 in 11.1)
names(DEG)
#  [1] "1"  "2"  "3"  "4"  "5"  "6"  "9"  "10" "12" "14" "17" "18" "19" "7"  "11" "15" "16" "13" "8"

celltype_col # [1] "label"
CP_cluster # '13'
CM_cluster # '7'
CF_cluster # '6'

lengths(CTS) # a list
#    7   11   15   16    8 16.1   13
#   32   52   67   40   54   79   60

class(sce) # [1] "SingleCellExperiment"

dim(mat) # [1] 60  4
colnames(mat)
# [1] "CP_hi"  "CM_hi"  "CF_hi"  "CTS_13"

seed_TF # NULL
CTS_ID # '13'
CTS_name # 'CTS_13'

## set subfold =======================
(updir <- getwd())
dir.create(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)),
    showWarnings = FALSE, recursive = TRUE
)
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

NES_threshold <- 3

########################################################
##  input 6 -- data-driven --- RcisTarget predicted TF that are enriched among CTS genes
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
##   --- binary flag CTS genes

dim(mat) # [1] 54  4
### = load the binary annotation matrix build in code 24.0xxx ============================
files <- list.files(pattern = "^heatmap_blocked_CTS_") %>% grep("_v3.tsv", ., value = TRUE)
if (length(files) > 0) {
    fileName <- grep("_v3.tsv", files, value = TRUE)
    mat <- read.table(fileName, sep = "\t", header = T, check.names = FALSE)

    saved_variables <- readRDS(file = "cisTarget_variables.rds")

    x <- saved_variables$x
    key_TFs <- saved_variables$key_TFs
    motifAnnot_sub <- saved_variables$motifAnnot_sub
} else {
    #### focus on the this CTS member gene-enriched motifs, with NES_threshold control =========================
    motifAnnot_sub <- get_regulators_from_motifs(cisTarget.res, CTS_ID, NES_threshold, motifAnnot = motifAnnot)
    keys <- unique(motifAnnot_sub$regulators)
    if (any(is.na(keys))) keys <- keys[!is.na(keys)]
    keys
    # At NES>=4.5: 3 motif groups
    # [1] "HAND1;HAND2..."   "PLAG1..."   "FOXD2;FOXA2..."
    # (NES=3 had ~19 groups; NES=4.5 focuses on HAND1/2, PLAG1, FOXD2/FOXA2)

    ## add to the binary annotation matrix the motif-based TF-target annotations
    for (key in keys) {
        motif_TF_highConf <- motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf
        tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold & (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf))
        genes <- unique(unlist(strsplit(tmp$enrichedGenes, ";")))
        mat[, paste0("cisTarget_", key, ".motif_target")] <- ifelse(rownames(mat) %in% genes, "1", "0")
    }

    dim(mat) # [1] 60 78  (4 base cols + 74 cisTarget cols at NES=3)
    colnames(mat)

    ##  extend seed_TF candidates to those 1) being CTS itself, or 2) CTS_enriched motifs having target genes highly expressed in CP, CM, or CF ===============
    x <- colnames(mat)[grepl("cisTarget_", colnames(mat))]
    x <- gsub("cisTarget_", "", x) %>% gsub("\\.motif_target", "", .)
    x <- unlist(strsplit(x, ";")) %>% unique()
    x
    # At NES>=4.5: individual TFs = HAND1, HAND2 (from HAND1;HAND2 motif),
    #              PLAG1, FOXD2, FOXA2 (from respective motifs); stripped of annotation text

    key_TFs <- seed_TF
    for (i in x) {
        for (j in c(CM_cluster, CF_cluster, CP_cluster)) {
            if (i %in% DEG[[j]] & i %in% rownames(sce)) {
                cat(paste0(i, " is in DEG of ", j, "\n"))
                key_TFs <- c(key_TFs, i)
            }
        }
    }
    # (no output: HAND1, HAND2, PLAG1, FOXD2, FOXA2 are NOT in DEG_13, DEG_7, or DEG_6)

    for (i in x) {
        if (i %in% CTS[[CTS_ID]]) {
            cat(paste0(i, " is in CTS of ", CTS_ID, "\n"))
            key_TFs <- c(key_TFs, i)
        }
    }
    # (no output: none of the NES>=4.5 TF families are in CTS_13)

    key_TFs <- unique(key_TFs)
    key_TFs # NULL  (seed_TF=NULL; motif family names don't match CTS_13 gene names; override manually below)

    mat <- as.data.frame(mat)

    # Manually choose TFs that are present in cisTarget. If two TFs share the same cisTarget motif, choose one.
    ## MANUAL SELECTION at NES=3: HAND2 (lateral plate mesoderm / vascular progenitor; motif targets
    ## TIE1, RHOJ, PLXND1, ETS1 — all in CTS_13 and enriched in DEG_6/endothelium), FOXA2 (pioneer
    ## TF; targets ETV2, GNGT2, GRAP, DUSP2 in DEG_6 and CTS_13, driving hemogenic endothelium),
    ## GATA3 (GATA family; hematopoietic lineage specification, uniquely matched column at NES=3),
    ## SOX17 (master regulator of hemogenic endothelium identity and EHT transition, uniquely matched).
    ## Drop PLAG1/ZNFs/CLOCK (non-specific). Skip ETV2 and RUNX1 (each match multiple cisTarget
    ## columns at NES=3 — ERG;ETV2;FLI1 + DLX/ETV family, and RUNX1 + CBFB;RUNX1 family —
    ## would break names(x) <- key_TFs. GATA3/SOX17 subsume biological coverage uniquely.)
    ## Order must match column order in mat (HAND2→FOXA2→GATA3→SOX17).
    key_TFs <- c("HAND2", "FOXA2", "GATA3", "SOX17")
    # Use first-match-per-TF to avoid NA columns when a key_TF name appears in multiple cisTarget columns
    # (e.g., "GATA3" appears in both cisTarget_GATA3.motif_target and several super-family columns).
    x <- sapply(key_TFs, function(p) {
        matches <- which(grepl("cisTarget_", colnames(mat)) & grepl(p, colnames(mat)))
        if (length(matches) == 0) NA_integer_ else matches[1]
    })
    x <- x[!is.na(x)]
    if (length(x) > 0) {
        for (j in seq_along(x)) {
            key <- names(x)[j]
            mat[, paste0(key, "_CP_candidate")] <- ifelse(mat[, "CP_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CM_candidate")] <- ifelse(mat[, "CM_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CF_candidate")] <- ifelse(mat[, "CF_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
        }
    }

    cat(paste0("key_TFs: ", paste(key_TFs, collapse = "_"), "\n")) # key_TFs: HAND2_FOXA2_GATA3_SOX17

    if (length(key_TFs) > 0) {
        fileName <- paste0("heatmap_blocked_", CTS_name, "_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv")
        write.table(mat, file = fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
        saveRDS(list(x = x, key_TFs = key_TFs, motifAnnot_sub = motifAnnot_sub), "cisTarget_variables.rds")
    } else {
        stop("No key TFs found for ", CTS_name, "\n")
    }
}

motif_TF_highConf <- gsub("cisTarget_", "", colnames(mat)[x]) %>% gsub("\\.motif_target", "", .)
print(motif_TF_highConf)

###  heatmap confirming key TF's self impact by checking its targets among the CTS_CP ================
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
    }
}

##########################################################
## --  identify_TF_targeted_pull_candidate -- subset of CTS[['CP']] that are
## exclusively highly expressed (HiG) in CM (or CF)
library(BioNet)
packageVersion("BioNet") # '1.56.0'
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

# code_core has no ChIP-seq; HiGCTS slot in identify_TF_targeted_pull_candidate is redundant
# since CP candidates are already selected by DEG + cisTarget. Use CTS as the backbone.
graph_list[[paste0("HiG", CTS_name)]] <- graph_list[[CTS_name]]

for (key in key_TFs) {
    key_column <- which(grepl(key, colnames(mat)) & grepl("cisTarget_", colnames(mat)))
    if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key

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
# [1] "CTSHiG.CP_TF.target"        "CTS.CP_TF.target_HiGCM"
# [3] "CTS.CP_TF.target_HiGCF"     "CTSHiG.CP_TF.target_CPopen"

##########################################################
make_layout <- function(g, seed) {
    set.seed(seed)
    layout_with_fr(g)
}

for (key in key_TFs) {
    if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key

    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))

    res <- fill_TF_targeting_predicted_edges(graph_TF_list,
        linkage_name = "CM", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CM_cluster, TF_symbol = key_in_TFfamily,
        HVG = rownames(sce)
    )
    cat(paste0(key, ' CM vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), "\n"))
    if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CM_final.rds"))

    res_cf <- fill_TF_targeting_predicted_edges(graph_TF_list,
        linkage_name = "CF", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
        HVG = rownames(sce)
    )
    cat(paste0(key, ' CF vcount(res_cf[["g_CT_sub"]]): ', vcount(res_cf[["g_CT_sub"]]), "\n"))
    if (vcount(res_cf[["g_CT_sub"]]) > 0) saveRDS(res_cf, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CF_final.rds"))
}
# HAND2 CM vcount(res[["g_CT_sub"]]): 2
# HAND2 CF vcount(res_cf[["g_CT_sub"]]): 9
# FOXA2 CM vcount(res[["g_CT_sub"]]): 2
# FOXA2 CF vcount(res_cf[["g_CT_sub"]]): 8
# GATA3 CM vcount(res[["g_CT_sub"]]): 0
# GATA3 CF vcount(res_cf[["g_CT_sub"]]): 5
# SOX17 CM vcount(res[["g_CT_sub"]]): 2
# SOX17 CF vcount(res_cf[["g_CT_sub"]]): 4

### reporting ==========
print((files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds")))
# [1] "PPI_graph_FOXA2_GRN_prediction_CTS_13_CF_final.rds"
# [2] "PPI_graph_FOXA2_GRN_prediction_CTS_13_CM_final.rds"
# [3] "PPI_graph_GATA3_GRN_prediction_CTS_13_CF_final.rds"
# [4] "PPI_graph_HAND2_GRN_prediction_CTS_13_CF_final.rds"
# [5] "PPI_graph_HAND2_GRN_prediction_CTS_13_CM_final.rds"
# [6] "PPI_graph_SOX17_GRN_prediction_CTS_13_CF_final.rds"
# [7] "PPI_graph_SOX17_GRN_prediction_CTS_13_CM_final.rds"

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
        if (nrow(change_df) == 0) {
            cat("No edge changes for", key_in_TFfamily, "-", pull, "\n")
            next
        }
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

dim(final_table) # [1] 36 13  (CF: 33 rows across FOXA2/HAND2/GATA3/SOX17 subnetworks; CM: 3 rows)
# CF FOXA2 (NES=4.68): DUSP2→FOXA2(-), FOXA2→GRAP(+), FOXA2→GNGT2(+), ETV2→FOXA2(-),
#     FOXA2→UPP1(+), FOXA2→N4BP3(+), FOXA2→PLXND1(+)
# CF HAND2 (NES=5.17): ETS1→HAND2(-), HAND2→TIE1(+), PLXND1→RHOJ(-), HAND2→PLXND1(+),
#     HAND2→KCNE3(+), HAND2→TSPAN8(+), HAND2→RHOJ(+), PLXND1→TIE1(-),
#     ETS1→TIE1(-), CD38→HAND2(+), KCNE3→PLXND1(-), CD38→ETS1(-),
#     RHOJ→TIE1(-), HAND2→UBASH3B(+) [+ CDH5/ETV2/ETS1/PECAM1 edges]
# CF GATA3 (NES=3.66): GATA3→PECAM1(+), ETV2→GATA3(-), CDH5→GATA3(+), ETS1→GATA3(+)
# CF SOX17 (NES=3.63): MYZAP→SOX17(-), DUSP2→SOX17(-), SOX17→TIE1(-)
# CM: DUSP2→FOXA2(-) (NES=4.68); HAND2→UBASH3B(-) (NES=5.17); DUSP2→SOX17(-) (NES=3.63)

write.table(final_table,
    file = paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)
