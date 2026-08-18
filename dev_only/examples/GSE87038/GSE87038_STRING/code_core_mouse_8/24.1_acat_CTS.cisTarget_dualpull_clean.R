## 24.1_acat_CTS.cisTarget_dualpull_clean.R — code_core_mouse_8 (CORE: no ATAC / no ChIP)
##
## Faithful to GSE87038_STRING/code_core_13/24.1_acat_CTS.cisTarget_dualpull_clean.R, but:
##   - CP/CM/CF clusters set for the C8 arm (CP=8, CM=17, CF=16)
##   - key_TFs selected automatically (PageRank + cisTarget), porting the selector from
##     PI's GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/code_core_4_9vs11/24.1...R,
##     instead of code_core_13's hand-picked key_TFs <- c("HAND2","FOXA2","GATA3","SOX17")
##   - Native mouse gene case kept throughout (matches 24.0's mm10 cisTarget + native TF_mouse)
##
## identify_TF_targeted_pull_candidate()/fill_TF_targeting_predicted_edges()
## (R/celltype_specific_weight_v10.R) hard-uppercase graph vertex names internally and
## reference the global DEG/sce directly -- native case from 24.0 is deliberately switched
## to uppercase right before this section starts (see "From here on..." comment below) to
## match. Several grepl()/intersect() calls further down also need ignore.case=TRUE /
## toupper() for the same reason -- see inline comments.
##
## Automatic selector (use_automatic_key_TFs=TRUE) picked ALX1, ZFP207, FOS
## (NES_threshold=3); the top PageRank TFs for HiGCTS_8 (ISL1, IRX3 -- ISL1 is the
## literature CP marker, confirmed by 12.0) were skipped because neither has a cisTarget
## motif column at this NES threshold. ALX1/FOS come from 100-/105-member cisTarget
## megaclusters (flag 39%/30% of all CTS_8 genes -- not specific), and ZFP207 is not a
## classical sequence-specific TF. Pending PI review (2026-07-29), default is manual:
## Prrx2, Sox4, Ezh2 -- solo (non-megacluster) cisTarget columns that are also
## independently DEG in the CP cluster (8) itself; EZH2 in particular has published roles
## in cardiac progenitor proliferation/differentiation. Provisional, not vetted by domain
## expert -- compare against code_mouse_8's manual TCF3/HMGA2/RARB (different data
## modality: ATAC/ChIP-informed, hg38) once a call is made.
##
## (NB: when rerunning after changing manual_key_TFs/toggle, delete
## cisTarget_predicted_8/heatmap_blocked_CTS_*_v3.tsv + cisTarget_variables.rds first --
## the "load cached mat" branch above will otherwise silently reuse the OLD key_TFs.)

# Set True if running 24.0 the first time
rebuild_mat <- TRUE
source(here::here("examples", "cardiac", "GSE87038", "GSE87038_STRING", "code_core_mouse_8", "24.0_acat_load_input_clean.R"))

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

## check the loaded objects =========================
seed_TF # NULL (set NULL in 24.0 USER INPUT; key_TFs resolved automatically below)
names(graph_list)
names(DEG)

celltype_col # [1] "label"
CP_cluster # '8'
CM_cluster # '17'
CF_cluster # '16'

lengths(CTS) # a list
#    7   11   15   16   13    8 16.1
#   32   52   67   40   60   54   79

class(sce) # [1] "SingleCellExperiment"

dim(mat) # 54  4
colnames(mat) # "CP_hi" "CM_hi" "CF_hi" "CTS_8"

seed_TF
CTS_ID # '8'
CTS_name # 'CTS_8'

## set subfold =======================
(updir <- getwd())
dir.create(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)),
    showWarnings = FALSE, recursive = TRUE
)
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

NES_threshold <- 3

## Felix 12.0: top N TFs by PageRank in HiGCTS_<CP> (among TFs, within top gene_top_n genes)
top_TF_rank <- 3
gene_top_n  <- 20

## Set TRUE to select key_TFs automatically (PageRank + cisTarget). Default FALSE uses
## the manual pick below pending PI review -- see header note for rationale/tradeoffs.
use_automatic_key_TFs <- FALSE
manual_key_TFs <- c("Prrx2", "Sox4", "Ezh2")

########################################################
##  input 6 -- data-driven --- RcisTarget predicted TF that are enriched among CTS genes (mm10)
library(RcisTarget)
packageVersion("RcisTarget") # '1.29.0'
library(data.table)

## motifAnnot (mgi) already loaded by 24.0 above — reused here, no re-load needed.
dim(motifAnnot) # 252126      7

cisTarget.res <- readRDS(file = "../cisTarget_targets_in_all_CTS.rds")
write.table(subset(cisTarget.res, NES >= NES_threshold & geneSet %in% c(CP_cluster)),
    paste0("../cisTarget_targets_in_", CP_cluster, "_NES", NES_threshold, ".txt"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

########################################################
##   --- binary flag CTS genes

dim(mat) # 54  4
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

    ## add to the binary annotation matrix the motif-based TF-target annotations
    for (key in keys) {
        motif_TF_highConf <- motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf
        tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold & (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf))
        genes <- unique(unlist(strsplit(tmp$enrichedGenes, ";")))
        mat[, paste0("cisTarget_", key, ".motif_target")] <- ifelse(rownames(mat) %in% genes, "1", "0")
    }

    dim(mat) # 54 74  (4 base cols + 70 cisTarget motif cols at NES=3; +9 _CP/CM/CF_candidate cols -> 83 after key_TFs below)
    colnames(mat)

    ## ====================================================================
    ## key_TFs selected automatically: PageRank (12.0) + cisTarget motif review.
    ## Ported from PI's GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/
    ## code_core_4_9vs11/24.1_acat_CTS.cisTarget_dualpull_clean.R (replaces
    ## code_core_13's hand-picked key_TFs <- c("HAND2","FOXA2","GATA3","SOX17")).
    ## ====================================================================
    split_motif_regulators <- function(keys) {
      keys <- as.character(keys)
      keys <- keys[!is.na(keys) & nzchar(keys)]
      unique(unlist(strsplit(keys, ";", fixed = TRUE)))
    }

    cistarget_cols_for_tf <- function(tf, mat) {
      tf <- toupper(as.character(tf))
      motif_cols <- which(
        grepl("^cisTarget_.+\\.motif_target$", colnames(mat)) &
          !grepl("\\.merged_motif_target$", colnames(mat))
      )
      idx <- integer()
      col_solo <- match(paste0("cisTarget_", tf, ".motif_target"), colnames(mat))
      if (!is.na(col_solo)) idx <- c(idx, col_solo)
      for (j in motif_cols) {
        if (j %in% idx) next
        regs <- sub("^cisTarget_", "", colnames(mat)[j])
        regs <- sub("\\.motif_target$", "", regs)
        parts <- toupper(strsplit(regs, ";", fixed = TRUE)[[1]])
        if (tf %in% parts) idx <- c(idx, j)
      }
      unique(idx)
    }

    tf_has_motif_columns <- function(tf, mat) {
      length(cistarget_cols_for_tf(tf, mat)) > 0L
    }

    ## All TFs in PageRank order within HiGCTS_<CP>, else CTS_<CP>.
    pagerank_tfs_in_higcts <- function(df_PageRank, TF_mouse, CP_cluster, gene_top_n = 20L) {
      higcts_sig <- paste0("HiGCTS_", CP_cluster)
      cts_sig <- paste0("CTS_", CP_cluster)
      if (higcts_sig %in% df_PageRank$signature) {
        pr_sig <- higcts_sig
      } else if (cts_sig %in% df_PageRank$signature) {
        pr_sig <- cts_sig
        message("[24.1] No ", higcts_sig, " in df_PageRank — using ", cts_sig)
      } else {
        stop(
          "No ", higcts_sig, " or ", cts_sig,
          " in df_PageRank — run 11.3 then 12.0 first",
          call. = FALSE
        )
      }
      tf_mouse <- toupper(as.character(TF_mouse))
      df_plot <- df_PageRank[df_PageRank$signature == pr_sig, , drop = FALSE]
      df_plot$gene <- toupper(as.character(df_plot$gene))
      df_plot <- df_plot[order(-df_plot$PageRank), , drop = FALSE]
      df_plot$rank <- seq_len(nrow(df_plot))
      is_tf <- df_plot$gene %in% tf_mouse
      unique(df_plot$gene[is_tf & df_plot$rank <= gene_top_n])
    }

    ## Felix 24.1 auto-extend: individual motif regulators in DEG and/or CTS
    extend_key_tfs_from_cistarget <- function(
        key_TFs, motifAnnot_sub, TF_mouse, DEG, CTS, CTS_ID, CP_cluster, CM_cluster, CF_cluster, sce
    ) {
      regulators <- split_motif_regulators(motifAnnot_sub$regulators)
      regulators <- intersect(toupper(regulators), toupper(as.character(TF_mouse)))
      deg_clusters <- as.character(c(CP_cluster, CM_cluster, CF_cluster))
      cts_genes <- toupper(as.character(CTS[[as.character(CTS_ID)]]))
      sce_genes <- toupper(rownames(sce))
      for (i in regulators) {
        in_deg <- any(vapply(deg_clusters, function(j) i %in% toupper(DEG[[j]]), logical(1L)))
        in_cts <- i %in% cts_genes
        if ((in_deg && i %in% sce_genes) || in_cts) {
          key_TFs <- c(key_TFs, i)
        }
      }
      unique(key_TFs)
    }

    ## One TF per cisTarget column (Felix manual: pick one TF when motifs are shared).
    dedupe_key_tfs_by_cistarget_col <- function(candidate_tfs, mat, top_TF_rank = 3L) {
      candidate_tfs <- toupper(as.character(candidate_tfs))
      candidate_tfs <- candidate_tfs[!is.na(candidate_tfs) & nzchar(candidate_tfs)]
      candidate_tfs <- unique(candidate_tfs)
      col_idx <- integer()
      tf_names <- character()
      for (tf in candidate_tfs) {
        cols <- cistarget_cols_for_tf(tf, mat)
        if (!length(cols)) next
        j <- cols[1L]
        if (j %in% col_idx) next
        col_idx <- c(col_idx, j)
        tf_names <- c(tf_names, tf)
        if (length(col_idx) >= top_TF_rank) break
      }
      x <- stats::setNames(col_idx, tf_names)
      if (length(x)) {
        for (k in seq_along(x)) message(names(x)[k], "\t", x[k], "\t", colnames(mat)[x[k]])
      }
      x
    }

    key_tfs_pageRank_and_cistarget <- function(
        df_PageRank, TF_mouse, CP_cluster, mat, motifAnnot_sub,
        DEG, CTS, CTS_ID, CM_cluster, CF_cluster, sce,
        top_TF_rank = 3L, gene_top_n = 20L
    ) {
      tf_mouse <- toupper(as.character(TF_mouse))
      pr_tfs <- pagerank_tfs_in_higcts(df_PageRank, TF_mouse, CP_cluster, gene_top_n)
      pr_with_col <- pr_tfs[vapply(pr_tfs, tf_has_motif_columns, logical(1L), mat = mat)]
      pr_skipped <- setdiff(pr_tfs, pr_with_col)
      if (length(pr_skipped)) {
        message(
          "[24.1] PageRank TF(s) skipped (no cisTarget column at NES>=", NES_threshold, "): ",
          paste(pr_skipped, collapse = ", ")
        )
      }

      solo_cols <- grep("^cisTarget_[^;]+\\.motif_target$", colnames(mat), value = TRUE)
      solo_tfs <- toupper(gsub("^cisTarget_|\\.motif_target$", "", solo_cols))
      solo_tfs <- intersect(solo_tfs, tf_mouse)

      motif_tfs <- split_motif_regulators(motifAnnot_sub$regulators)
      motif_tfs <- intersect(toupper(motif_tfs), tf_mouse)
      motif_tfs <- motif_tfs[vapply(motif_tfs, tf_has_motif_columns, logical(1L), mat = mat)]

      extended <- extend_key_tfs_from_cistarget(
        character(), motifAnnot_sub, TF_mouse, DEG, CTS, CTS_ID,
        CP_cluster, CM_cluster, CF_cluster, sce
      )

      candidates <- unique(c(pr_with_col, extended, solo_tfs, motif_tfs))
      x <- dedupe_key_tfs_by_cistarget_col(candidates, mat, top_TF_rank = top_TF_rank)
      list(x = x, pr_with_col = pr_with_col)
    }

    if (use_automatic_key_TFs) {
      ppi_path <- paste0(wd, "results_core/PPI_weight/")
      df_PageRank <- readRDS(file.path(ppi_path, "df_PAGERANK_strength_ANND.rewiring.P.rds"))

      kt_res <- key_tfs_pageRank_and_cistarget(
        df_PageRank, TF_mouse, CP_cluster, mat, motifAnnot_sub,
        DEG, CTS, CTS_ID, CM_cluster, CF_cluster, sce,
        top_TF_rank = top_TF_rank, gene_top_n = gene_top_n
      )
      x <- kt_res$x
      key_TFs <- names(x)
      key_TFs <- key_TFs[!is.na(key_TFs) & nzchar(key_TFs)]
      message(
        "[24.1] key_TFs (automatic, PageRank + cisTarget motifs, top ", top_TF_rank, "): ",
        if (length(key_TFs)) paste(key_TFs, collapse = ", ") else "(none)"
      )
    } else {
      ## Manual pick pending PI review -- see header note. Reuses the same
      ## dedupe_key_tfs_by_cistarget_col() helper the automatic path uses, so each TF
      ## still needs a cisTarget column in mat (structural requirement of the dual-pull
      ## method) -- it's just picked by human judgment instead of PageRank+column-order.
      x <- dedupe_key_tfs_by_cistarget_col(manual_key_TFs, mat, top_TF_rank = length(manual_key_TFs))
      key_TFs <- names(x)
      key_TFs <- key_TFs[!is.na(key_TFs) & nzchar(key_TFs)]
      message(
        "[24.1] key_TFs (manual): ",
        if (length(key_TFs)) paste(key_TFs, collapse = ", ") else "(none)"
      )
      missing_manual <- setdiff(toupper(manual_key_TFs), key_TFs)
      if (length(missing_manual)) {
        message("[24.1] manual key_TFs with no cisTarget column (dropped): ", paste(missing_manual, collapse = ", "))
      }
    }

    if (!length(key_TFs)) {
      stop("No key_TFs — run 11.3/12.0, check NES_threshold, or set manual_key_TFs", call. = FALSE)
    }

    if (length(x) > 0) {
        for (j in seq_along(x)) {
            key <- names(x)[j]
            mat[, paste0(key, "_CP_candidate")] <- ifelse(mat[, "CP_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CM_candidate")] <- ifelse(mat[, "CM_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
            mat[, paste0(key, "_CF_candidate")] <- ifelse(mat[, "CF_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
        }
    }

    cat(paste0("key_TFs: ", paste(key_TFs, collapse = "_"), "\n"))

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
        ## key_TFs is uppercase; family members here are native mouse case -- match
        ## case-insensitively and keep the native-case spelling (matches mat/CTS/DEG).
        family_members <- strsplit(key, ";", fixed = T) %>% unlist()
        key_in_TFfamily <- family_members[toupper(family_members) %in% toupper(key_TFs)] %>%
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

## From here on, the shared GRN-extraction helpers (R/celltype_specific_weight_v10.R:
## identify_TF_targeted_pull_candidate, fill_TF_targeting_predicted_edges) internally
## uppercase graph vertex names ("for mouse dataset to reuse the mat annotation") and
## reference the global DEG/sce directly. Uppercase to match -- native case was only
## needed above for mm10 cisTarget gene-set matching, which is already complete.
rownames(mat) <- toupper(rownames(mat))
for (i in seq_along(graph_list)) V(graph_list[[i]])$name <- toupper(V(graph_list[[i]])$name)
rownames(sce) <- toupper(rownames(sce))
DEG <- lapply(DEG, toupper)

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
    ## key_TFs is uppercase (from dedupe_key_tfs_by_cistarget_col); mat's cisTarget_
    ## column names are native mouse case (from motifAnnot_sub$regulators) -- match
    ## case-insensitively.
    key_column <- which(grepl(key, colnames(mat), ignore.case = TRUE) & grepl("cisTarget_", colnames(mat)))
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

for (key in key_TFs) {
    if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key

    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
    plot_TF_targeted_pull_candidate(graph_TF_list, key_in_TFfamily, CTS_name, saveFigure = TRUE)
}

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

### reporting ==========
print((files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds")))

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

    ## key/key_in_TFfamily are uppercase (from filenames); motifAnnot_sub$regulators
    ## and change_df$from/to are native mouse case -- match case-insensitively.
    x_key <- paste0(key, "_motif_target")
    y <- motifAnnot_sub[which(
        toupper(motifAnnot_sub$regulators) %in% toupper(key) |
            grepl(key, motifAnnot_sub$regulators, ignore.case = TRUE)
    ), ]
    ## unique(): case-insensitive regulator match above can return >1 row (e.g. a
    ## solo column and a family column both containing key) -- avoid a recycling
    ## warning against cisTarget.res$motif/TF_highConf below.
    motif_TF_highConf_val <- unique(y$motif_TF_highConf)
    tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold &
                      (motif %in% motif_TF_highConf_val | TF_highConf %in% motif_TF_highConf_val))

    tmp1 <- if (nrow(tmp) > 0) tmp[1, , drop = FALSE] else data.frame(TF_highConf = NA, motif = NA, NES = NA)
    change_df <- cbind(linkage = pull, change_df,
                       TF_highConf = tmp1$TF_highConf, motif = tmp1$motif, NES = tmp1$NES)
    is_key_edge <- toupper(change_df$from) == toupper(key_in_TFfamily) | toupper(change_df$to) == toupper(key_in_TFfamily)
    change_df$TF_highConf[which(!is_key_edge)] <- ""
    change_df$motif[which(!is_key_edge)] <- ""
    change_df$NES[which(!is_key_edge)] <- ""

    final_table <- rbind(final_table, change_df)
}

dim(final_table) # [1] 11 13  (CF: ALX1=1, FOS=2, ZFP207=1; CM: ALX1=3, FOS=3, ZFP207=1)
# CF ALX1: ALX1->DAB2(increase); CM ALX1: ALX1->BMP7, ALX1->DAB2, ALX1->FAM213A (all increase)
# CF FOS: FOS->MPPED2, FOS->TPM2 (decrease); CM FOS: FOS->MPPED2(lost), FOS->TPM2(lost), FGF8->FOS(decrease)
# CF/CM ZFP207: MPPED2->ZFP207 (decrease/lost)
# NB: automatic selector's key_TFs (ALX1/ZFP207/FOS) differ from code_mouse_8's manual
# TCF3/HMGA2/RARB and from ISL1 (the literature CP marker, confirmed by 12.0 but skipped
# here for lacking a cisTarget motif column at NES>=3) -- review before treating as final.

write.table(final_table,
    file = paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)
