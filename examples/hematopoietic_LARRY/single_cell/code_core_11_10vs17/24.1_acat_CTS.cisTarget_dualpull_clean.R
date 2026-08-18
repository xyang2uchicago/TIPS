## 24.1_acat_CTS.cisTarget_dualpull_clean.R — Felix code_core_13 (CORE: no ATAC / no ChIP)
##
## Faithful to GSE87038_STRING/code_core_13/24.1_acat_CTS.cisTarget_dualpull_clean.R
## Larry-only: paths, mm10 cisTarget, load_larry_sce().
## key_TFs = HiGCTS_<CP> PageRank TFs with cisTarget columns, then motif-enriched regulators
## (Felix: 12.0 PageRank prior + 24.1 cisTarget motif review; dual-pull needs motif columns).
##
## Manuscript terms (Table 1):
##   Reweighted edge        = context-specific progenitor-to-descendant edge (dualpull_final_table rows)
##   Lineage-biased edge    = reweighted edge preferentially supporting one descendant branch
##   Lineage-lean gene      = gene on arm-specific biased edges, excluding shared/bridge genes
## Heatmap QC shows both CM and CF arms (CP/CM/CF_hi); shared nodes are allowed in subnetworks.

########## BEGINNING OF USER INPUT ##########
code_dir <- get0("code_dir", ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/code_core_4_9vs11")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)

rebuild_mat_24.0 <- FALSE
rebuild_heatmap  <- TRUE    ## FALSE to reload heatmap_blocked_*_v3.tsv (Felix: file exists -> reload)

NES_threshold <- 3

## Felix 12.0: top N TFs by PageRank in HiGCTS_<CP> (among TFs, within top gene_top_n genes)
top_TF_rank <- 3
gene_top_n  <- 20

## Optional override (NULL = use PageRank only). Dedupe shared-motif TFs after cisTarget column map.
key_TFs_override <- NULL
# key_TFs_override <- c("GATA2")

## Motif target gene selection per key_TF (set in run_TIPS_core.R; override here if sourcing 24.1 alone):
##   solo  = cisTarget_<TF>.motif_target only (Felix default)
##   fam   = tightest family motif column containing TF (e.g. GATA1;GATA2;...;GATA6)
##   merge = union of all motif columns containing TF
if (!exists("motif_target_strategy", inherits = FALSE)) {
  motif_target_strategy <- Sys.getenv("TIPS_MOTIF_TARGET_STRATEGY", "solo")
}
# motif_target_strategy <- "fam"  # uncomment to override run_TIPS_core when sourcing 24.1 alone

save_heatmap_qc <- TRUE

## TIPS_extend: Felix-first dual-pull; HiG fallback + PageRank re-pick when CM_hi/CF_hi empty (1_5vs7-like arms)
if (!exists("tips_extend_auto_key_tfs", inherits = FALSE)) {
  tips_extend_auto_key_tfs <- TRUE
}
########## END OF USER INPUT ##########

motif_target_strategy <- match.arg(motif_target_strategy, c("solo", "fam", "merge"))
message("[24.1] motif_target_strategy=", motif_target_strategy)

rebuild_mat <- rebuild_mat_24.0
source(file.path(code_dir, "24.0_acat_load_input_clean.R"))

source(file.path(tips_r_dir, paste0("celltype_specific_weight_v", celltype_specific_weight_version, ".R")))
tips_extend_r <- file.path(
  normalizePath(dirname(get0("wd", ifnotfound = dirname(code_dir))), winslash = "/"),
  "TIPS_extend_v2.R"
)
if (file.exists(tips_extend_r)) {
  source(tips_extend_r)
} else {
  warning("[24.1] TIPS_extend_v2.R not found: ", tips_extend_r, call. = FALSE)
}

library(BioNet)
library(igraph)
library(tibble)
library(purrr)
library(RcisTarget)
library(data.table)
library(SingleCellExperiment)

if (!exists("motifAnnot")) {
  motifAnnot <- load_motif_annotations()
}

sce <- load_larry_sce()
sce <- sync_sce_cluster_labels(sce, group_col = group_col)
sce <- sce_upper_rownames(sce)
assign("DEG", DEG, envir = .GlobalEnv)

## check loaded objects (Felix code_core_13) =========================
message("[24.1] seed_TF: ", if (length(seed_TF)) paste(seed_TF, collapse = ", ") else "(none)")
names(graph_list)
names(DEG)
colnames(mat)
lengths(CTS)
stopifnot(all(c("CP_hi", "CM_hi", "CF_hi", CTS_name) %in% colnames(mat)))

(updir <- getwd())
dir.create(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)), showWarnings = FALSE, recursive = TRUE)
setwd(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)))

########################################################
## input 6 — cisTarget enriched motifs among CTS genes
cisTarget.res <- readRDS(file = "../cisTarget_targets_in_all_CTS.rds")
write.table(
  subset(cisTarget.res, NES >= NES_threshold & geneSet %in% CP_cluster),
  paste0("../cisTarget_targets_in_", CP_cluster, "_NES", NES_threshold, ".txt"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

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

motif_col_regulator_count <- function(col_idx, mat) {
  regs <- sub("^cisTarget_", "", colnames(mat)[col_idx])
  regs <- sub("\\.motif_target$", "", regs)
  length(strsplit(regs, ";", fixed = TRUE)[[1]])
}

fam_col_for_tf <- function(tf, mat) {
  tf <- toupper(as.character(tf))
  cols <- cistarget_cols_for_tf(tf, mat)
  if (!length(cols)) return(NA_integer_)
  solo_j <- match(paste0("cisTarget_", tf, ".motif_target"), colnames(mat))
  fam_cols <- cols[if (is.na(solo_j)) TRUE else cols != solo_j]
  if (!length(fam_cols)) return(solo_j)
  fam_cols[which.min(vapply(fam_cols, motif_col_regulator_count, integer(1L), mat = mat))]
}

add_merged_motif_columns <- function(mat, tfs, strategy = "solo") {
  strategy <- match.arg(strategy, c("solo", "fam", "merge"))
  if (!identical(strategy, "merge")) return(mat)
  tfs <- unique(toupper(as.character(tfs)))
  tfs <- tfs[!is.na(tfs) & nzchar(tfs)]
  for (tf in tfs) {
    cols <- cistarget_cols_for_tf(tf, mat)
    if (length(cols) <= 1L) next
    merged_name <- paste0("cisTarget_", tf, ".merged_motif_target")
    if (merged_name %in% colnames(mat)) next
    vals <- mat[, cols, drop = FALSE]
    merged <- apply(vals, 1L, function(row) any(row == "1" | row == 1L))
    mat[[merged_name]] <- ifelse(merged, "1", "0")
    message(
      "[24.1] merge motif targets ", tf, ": ", length(cols),
      " column(s) -> ", sum(merged), " genes (", merged_name, ")"
    )
  }
  mat
}

cistarget_col_for_tf <- function(tf, mat, strategy = "solo") {
  strategy <- match.arg(strategy, c("solo", "fam", "merge"))
  cols <- cistarget_cols_for_tf(tf, mat)
  if (!length(cols)) return(NA_integer_)
  if (identical(strategy, "solo")) return(cols[1L])
  if (identical(strategy, "fam")) return(fam_col_for_tf(tf, mat))
  tf <- toupper(as.character(tf))
  merged_name <- paste0("cisTarget_", tf, ".merged_motif_target")
  if (merged_name %in% colnames(mat)) return(match(merged_name, colnames(mat)))
  cols[1L]
}

## One TF per cisTarget column (Felix manual: pick one TF when motifs are shared).
dedupe_key_tfs_by_cistarget_col <- function(
    candidate_tfs, mat, top_TF_rank = 3L, motif_target_strategy = "solo"
) {
  motif_target_strategy <- match.arg(motif_target_strategy, c("solo", "fam", "merge"))
  candidate_tfs <- toupper(as.character(candidate_tfs))
  candidate_tfs <- candidate_tfs[!is.na(candidate_tfs) & nzchar(candidate_tfs)]
  candidate_tfs <- unique(candidate_tfs)
  col_idx <- integer()
  tf_names <- character()
  for (tf in candidate_tfs) {
    j <- cistarget_col_for_tf(tf, mat, strategy = motif_target_strategy)
    if (is.na(j) || j %in% col_idx) next
    col_idx <- c(col_idx, j)
    tf_names <- c(tf_names, tf)
    if (length(col_idx) >= top_TF_rank) break
  }
  x <- stats::setNames(col_idx, tf_names)
  if (length(x)) {
    for (k in seq_along(x)) {
      message(names(x)[k], "\t", x[k], "\t", colnames(mat)[x[k]])
    }
  }
  x
}

map_key_tfs_to_cistarget_cols <- function(
    key_tfs, mat, top_TF_rank = Inf, motif_target_strategy = "solo"
) {
  dedupe_key_tfs_by_cistarget_col(
    key_tfs, mat, top_TF_rank = top_TF_rank, motif_target_strategy = motif_target_strategy
  )
}

motif_regulators_from_colname <- function(col_name) {
  regs <- sub("^cisTarget_", "", col_name)
  regs <- sub("\\.(merged_)?motif_target$", "", regs)
  toupper(strsplit(regs, ";", fixed = TRUE)[[1]])
}

tips_hig_genes <- function(graph_list, CTS_name, CP_cluster) {
  hig_names <- c(paste0("HiGCTS_", CP_cluster), paste0("HiG", CTS_name), CTS_name)
  genes <- character()
  for (nm in hig_names) {
    g <- graph_list[[nm]]
    if (!is.null(g) && inherits(g, "igraph")) {
      genes <- c(genes, V(g)$name)
    }
  }
  unique(toupper(genes))
}

## Pick one TF symbol from a motif family: prefer CTS, then HiG, then key_TFs.
resolve_key_tf_symbol <- function(
    tf_key, mat = NULL, mat_col_idx = NA_integer_,
    key_TFs = character(), CTS = NULL, CTS_ID = NULL,
    graph_list = NULL, CTS_name = NULL, CP_cluster = NULL,
    TF_mouse = NULL
) {
  tf_key <- toupper(as.character(tf_key))
  tf_key <- tf_key[!is.na(tf_key) & nzchar(tf_key)]
  if (!length(tf_key)) {
    stop("resolve_key_tf_symbol: empty tf_key", call. = FALSE)
  }

  candidates <- unique(tf_key)
  if (!is.null(mat) && !is.na(mat_col_idx) && mat_col_idx >= 1L) {
    candidates <- unique(c(candidates, motif_regulators_from_colname(colnames(mat)[mat_col_idx])))
  }
  if (any(grepl(";", candidates, fixed = TRUE))) {
    candidates <- unique(unlist(strsplit(candidates, ";", fixed = TRUE)))
  }
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (length(candidates) == 1L) {
    if (identical(candidates, "HOX")) return("HOXB2")
    return(candidates)
  }

  cts_genes <- if (!is.null(CTS) && !is.null(CTS_ID) && as.character(CTS_ID) %in% names(CTS)) {
    toupper(as.character(CTS[[as.character(CTS_ID)]]))
  } else {
    character()
  }
  hig_genes <- if (!is.null(graph_list)) {
    tips_hig_genes(graph_list, CTS_name, CP_cluster)
  } else {
    character()
  }

  pick_first <- function(pool, from) {
    hit <- intersect(from, pool)
    if (length(hit)) hit[1L] else character()
  }

  resolved <- pick_first(candidates, cts_genes)
  if (!length(resolved)) resolved <- pick_first(candidates, hig_genes)
  if (!length(resolved)) resolved <- pick_first(toupper(key_TFs), candidates)
  if (!length(resolved) && !is.null(TF_mouse)) {
    resolved <- pick_first(candidates, toupper(as.character(TF_mouse)))
  }
  if (!length(resolved)) resolved <- candidates[1L]
  if (identical(resolved, "HOX")) resolved <- "HOXB2"
  resolved
}

resolve_key_tf_map <- function(
    x, mat, key_TFs, CTS, CTS_ID, graph_list, CTS_name, CP_cluster, TF_mouse
) {
  key_TFs <- toupper(as.character(key_TFs))
  key_TFs <- key_TFs[!is.na(key_TFs) & nzchar(key_TFs)]
  stats::setNames(
    vapply(names(x), function(k) {
      resolve_key_tf_symbol(
        k, mat = mat, mat_col_idx = unname(x[[k]]),
        key_TFs = key_TFs, CTS = CTS, CTS_ID = CTS_ID,
        graph_list = graph_list, CTS_name = CTS_name, CP_cluster = CP_cluster,
        TF_mouse = TF_mouse
      )
    }, character(1L)),
    names(x)
  )
}

tf_symbol_from_ppi_filename <- function(f) {
  ## [Caesar 2026-08-12] regmatches(regexpr()) returns the WHOLE match, so this
  ## returned the filename instead of the TF symbol -> key_in_TFfamily was wrong
  ## -> ns blanked TF_highConf/motif/NES on every row. sub() takes capture group 1.
  pat <- "^PPI_graph_([^_]+)_GRN_prediction_.*_final\\.rds$"
  if (!grepl(pat, f, perl = TRUE)) return(character())
  sub(pat, "\\1", f, perl = TRUE)
}

## All TFs in PageRank order within HiGCTS_<CP>, else CTS_<CP> (Felix 12.0 / rank_TF_CHD_in_PPIN)
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

## Felix 24.1 auto-extend: individual motif regulators in DEG and/or CTS (lung code_core pattern)
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

key_tfs_pageRank_and_cistarget <- function(
    df_PageRank, TF_mouse, CP_cluster, mat, motifAnnot_sub,
    DEG, CTS, CTS_ID, CM_cluster, CF_cluster, sce,
    top_TF_rank = 3L, gene_top_n = 20L, motif_target_strategy = "solo"
) {
  motif_target_strategy <- match.arg(motif_target_strategy, c("solo", "fam", "merge"))
  tf_mouse <- toupper(as.character(TF_mouse))
  pr_tfs <- pagerank_tfs_in_higcts(df_PageRank, TF_mouse, CP_cluster, gene_top_n)
  pr_with_col <- pr_tfs[vapply(pr_tfs, tf_has_motif_columns, logical(1L), mat = mat)]
  pr_skipped <- setdiff(pr_tfs, pr_with_col)
  if (length(pr_skipped)) {
    message(
      "[24.1] PageRank TF(s) skipped (no cisTarget column at NES>=",
      NES_threshold, "): ",
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
  mat <- add_merged_motif_columns(mat, candidates, motif_target_strategy)
  x <- dedupe_key_tfs_by_cistarget_col(
    candidates, mat, top_TF_rank = top_TF_rank, motif_target_strategy = motif_target_strategy
  )
  list(x = x, mat = mat, pr_with_col = pr_with_col)
}

########################################################
## step 1.2) binary flag CTS genes + cisTarget motif columns
files <- list.files(pattern = "^heatmap_blocked_CTS_") %>% grep("_v3.tsv", ., value = TRUE)
resolved_key_TFs <- NULL

if (length(files) > 0 && !rebuild_heatmap && file.exists("cisTarget_variables.rds")) {
  saved_variables <- readRDS("cisTarget_variables.rds")
  saved_strategy <- if (!is.null(saved_variables$motif_target_strategy)) {
    saved_variables$motif_target_strategy
  } else {
    "solo"
  }
  if (!identical(saved_strategy, motif_target_strategy)) {
    message(
      "[24.1] motif_target_strategy changed (", saved_strategy, " -> ", motif_target_strategy,
      "); set rebuild_heatmap <- TRUE"
    )
    rebuild_heatmap <- TRUE
  }
}

if (length(files) > 0 && !rebuild_heatmap) {
  message("[24.1] reload cached heatmap: ", files[1])
  fileName <- grep("_v3.tsv", files, value = TRUE)[1]
  mat <- read.table(fileName, sep = "\t", header = TRUE, check.names = FALSE)
  saved_variables <- readRDS("cisTarget_variables.rds")
  x <- saved_variables$x
  key_TFs <- saved_variables$key_TFs
  motifAnnot_sub <- saved_variables$motifAnnot_sub
  resolved_key_TFs <- saved_variables$resolved_key_TFs
  if (exists("tips_maybe_extend_key_tfs", mode = "function")) {
    df_PageRank <- readRDS(file.path(ppi_path, "df_PAGERANK_strength_ANND.rewiring.P.rds"))
  pr_with_col <- pagerank_tfs_in_higcts(df_PageRank, TF_mouse, CP_cluster, gene_top_n)
    pr_with_col <- pr_with_col[vapply(pr_with_col, tf_has_motif_columns, logical(1L), mat = mat)]
    ext <- tips_maybe_extend_key_tfs(
      x, mat, pr_with_col,
      function(tf) cistarget_col_for_tf(tf, mat, strategy = motif_target_strategy),
      graph_list, CM_cluster, CF_cluster,
      top_TF_rank = top_TF_rank,
      auto = tips_extend_auto_key_tfs,
      key_TFs_override = NULL
    )
    x <- ext$x
    mat <- ext$mat
    key_TFs <- ext$key_TFs
  }
} else {
  motifAnnot_sub <- get_regulators_from_motifs(
    cisTarget.res, CTS_ID, NES_threshold, motifAnnot = motifAnnot
  )
  keys <- unique(motifAnnot_sub$regulators)
  keys <- keys[!is.na(keys)]
  if (identical(species, "mouse")) {
    motifAnnot_sub$regulators <- toupper(motifAnnot_sub$regulators)
    keys <- toupper(keys)
  }

  for (key in keys) {
    motif_TF_highConf <- motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf
    tmp <- subset(
      cisTarget.res,
      geneSet == CTS_ID & NES >= NES_threshold &
        (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf)
    )
    genes <- unique(unlist(strsplit(tmp$enrichedGenes, ";")))
    if (identical(species, "mouse")) genes <- toupper(genes)
    mat[, paste0("cisTarget_", key, ".motif_target")] <- ifelse(rownames(mat) %in% genes, "1", "0")
  }

  mat <- as.data.frame(mat)

  df_PageRank <- readRDS(file.path(ppi_path, "df_PAGERANK_strength_ANND.rewiring.P.rds"))
  if (!exists("TF_mouse") || !length(TF_mouse)) {
    if (requireNamespace("dorothea", quietly = TRUE)) {
      data(dorothea_mm, package = "dorothea", envir = environment())
      TF_mouse <- unique(dorothea_mm$tf)
    } else {
      TF_mouse <- character()
    }
  }
  if (identical(species, "mouse")) TF_mouse <- toupper(TF_mouse)

  message("[24.1] --- cisTarget column map (key_TFs) ---")
  kt_res <- key_tfs_pageRank_and_cistarget(
    df_PageRank, TF_mouse, CP_cluster, mat, motifAnnot_sub,
    DEG, CTS, CTS_ID, CM_cluster, CF_cluster, sce,
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n,
    motif_target_strategy = motif_target_strategy
  )
  x <- kt_res$x
  mat <- kt_res$mat
  pr_with_col <- kt_res$pr_with_col

  if (!is.null(key_TFs_override) && length(key_TFs_override)) {
    key_TFs_override <- toupper(key_TFs_override[!is.na(key_TFs_override) & nzchar(key_TFs_override)])
    message("[24.1] key_TFs_override: ", paste(key_TFs_override, collapse = ", "))
    message("[24.1] --- cisTarget column map (key_TFs_override) ---")
    mat <- add_merged_motif_columns(mat, key_TFs_override, motif_target_strategy)
    x <- map_key_tfs_to_cistarget_cols(
      key_TFs_override, mat, top_TF_rank = top_TF_rank,
      motif_target_strategy = motif_target_strategy
    )
    if (exists("tips_refresh_dualpull_candidates", mode = "function")) {
      mat <- tips_refresh_dualpull_candidates(mat, x, graph_list, CM_cluster, CF_cluster)
    }
  } else if (exists("tips_maybe_extend_key_tfs", mode = "function")) {
    ext <- tips_maybe_extend_key_tfs(
      x, mat, pr_with_col,
      function(tf) cistarget_col_for_tf(tf, mat, strategy = motif_target_strategy),
      graph_list, CM_cluster, CF_cluster,
      top_TF_rank = top_TF_rank,
      auto = tips_extend_auto_key_tfs,
      key_TFs_override = NULL
    )
    x <- ext$x
    mat <- ext$mat
  }
  key_TFs <- names(x)
  key_TFs <- key_TFs[!is.na(key_TFs) & nzchar(key_TFs)]
  message("[24.1] key_TFs (PageRank + cisTarget motifs, top ", top_TF_rank, "): ",
          if (length(key_TFs)) paste(key_TFs, collapse = ", ") else "(none)")

  if (!length(key_TFs)) {
    stop("No key_TFs — run 11.3/12.0, check NES_threshold, or set key_TFs_override", call. = FALSE)
  }
  if (!length(x)) {
    stop(
      "No cisTarget columns matched key_TFs — check NES_threshold / motif enrichment for: ",
      paste(key_TFs, collapse = ", "),
      call. = FALSE
    )
  }

  resolved_key_TFs <- resolve_key_tf_map(
    x, mat, key_TFs, CTS, CTS_ID, graph_list, CTS_name, CP_cluster, TF_mouse
  )
  message(
    "[24.1] resolved PPI TF symbols (CTS > HiG > key_TFs): ",
    paste(names(resolved_key_TFs), resolved_key_TFs, sep = "=", collapse = ", ")
  )

  fileName <- paste0("heatmap_blocked_", CTS_name, "_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv")
  write.table(mat, file = fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  saveRDS(
    list(
      x = x, key_TFs = key_TFs, resolved_key_TFs = resolved_key_TFs,
      motifAnnot_sub = motifAnnot_sub,
      motif_target_strategy = motif_target_strategy
    ),
    "cisTarget_variables.rds"
  )
  message("[24.1] saved: ", fileName)
}

key_TFs <- key_TFs[!is.na(key_TFs) & nzchar(key_TFs)]
if (!exists("TF_mouse") || !length(TF_mouse)) {
  if (requireNamespace("dorothea", quietly = TRUE)) {
    data(dorothea_mm, package = "dorothea", envir = environment())
    TF_mouse <- unique(dorothea_mm$tf)
  } else {
    TF_mouse <- character()
  }
}
if (identical(species, "mouse")) TF_mouse <- toupper(TF_mouse)
if (is.null(resolved_key_TFs) || !length(resolved_key_TFs)) {
  resolved_key_TFs <- resolve_key_tf_map(
    x, mat, key_TFs, CTS, CTS_ID, graph_list, CTS_name, CP_cluster, TF_mouse
  )
  message(
    "[24.1] resolved PPI TF symbols (CTS > HiG > key_TFs): ",
    paste(names(resolved_key_TFs), resolved_key_TFs, sep = "=", collapse = ", ")
  )
}

print(resolved_key_TFs)

### heatmap QC — both CM and CF arms visible (Felix heatmap_pull_candidate) ================
if (isTRUE(save_heatmap_qc)) {
  for (key in names(x)) {
    key_in_TFfamily <- unname(resolved_key_TFs[[key]])
    if (!nzchar(key_in_TFfamily)) next
    p <- tryCatch(
      heatmap_pull_candidate(
        mat, graph_list, CTS_name, CHD,
        key = key, coding_genes = coding_genes, TF = TF_mouse,
        show_SMC_access = FALSE
      ),
      error = function(e) {
        message("heatmap skipped for '", key_in_TFfamily, "': ", e$message)
        NULL
      }
    )
    if (!is.null(p)) {
      pdf(
        file = paste0("heatmap_blocked_", CTS_name, "_cisTarget_", key_in_TFfamily, "_v3_coding_target.pdf"),
        height = 4
      )
      print(p)
      dev.off()
    }
  }
}

##########################################################
## identify_TF_targeted_pull_candidate — dual-pull subset
##########################################################
mat <- read.table(
  paste0("heatmap_blocked_", CTS_name, "_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE
)
if (exists("tips_refresh_dualpull_candidates", mode = "function")) {
  mat <- tips_refresh_dualpull_candidates(mat, x, graph_list, CM_cluster, CF_cluster)
}

## No ATAC data: dummy access columns (Felix code_core_13)
for (col in c("PCW6CP_access", "PCW8_CM_access", "PCW19_CM_access",
              "PCW8_CF_access", "PCW19_CF_access",
              "PCW8_SMC_access", "PCW19_SMC_access", "iEPC_access")) {
  mat[[col]] <- 1L
}

graph_list[[paste0("HiG", CTS_name)]] <- graph_list[[CTS_name]]

for (key in key_TFs) {
  key_column <- unname(x[[key]])
  if (is.na(key_column)) {
    key_column <- which(grepl(key, colnames(mat)) & grepl("cisTarget_", colnames(mat)))[1]
  }
  key_in_TFfamily <- unname(resolved_key_TFs[[key]])

  graph_TF_list <- identify_TF_targeted_pull_candidate(
    mat, graph_list, CTS_name, CHD,
    key = key,
    keep_selfloop = TRUE,
    TF_bound_column_name = key_column,
    TF_appendix = key,
    edge_colored_by_Maven2023_ISL1KO = FALSE,
    key_in_TFfamily = key_in_TFfamily
  )
  saveRDS(
    graph_TF_list,
    paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds")
  )
}

names(graph_TF_list)

##########################################################
## step 2.1) quantify edge weight changes (CM / CF pull) — Felix linkage_name "CM"/"CF"
##########################################################
for (key in key_TFs) {
  key_in_TFfamily <- unname(resolved_key_TFs[[key]])

  graph_TF_list <- readRDS(paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))

  res <- fill_TF_targeting_predicted_edges(
    graph_TF_list, linkage_name = "CM", graph_list,
    sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
    descendant_cluster_id = CM_cluster, TF_symbol = key_in_TFfamily,
    HVG = rownames(sce)
  )
  message(key, " (", key_in_TFfamily, ") CM vcount(g_CT_sub): ", vcount(res[["g_CT_sub"]]))
  if (vcount(res[["g_CT_sub"]]) > 0) {
    saveRDS(res, paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CM_final.rds"))
  }

  res_cf <- fill_TF_targeting_predicted_edges(
    graph_TF_list, linkage_name = "CF", graph_list,
    sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
    descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
    HVG = rownames(sce)
  )
  message(key, " (", key_in_TFfamily, ") CF vcount(g_CT_sub): ", vcount(res_cf[["g_CT_sub"]]))
  if (vcount(res_cf[["g_CT_sub"]]) > 0) {
    saveRDS(res_cf, paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CF_final.rds"))
  }
}

##########################################################
## step 2.2) reporting — lineage-biased edges (dualpull_final_table)
##########################################################
files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final\\.rds$")
print(files)

final_table <- NULL
for (f in files) {
  pull <- ifelse(grepl("_CF_final", f), "CF", "CM")
  key_raw <- tf_symbol_from_ppi_filename(f)
  if (!length(key_raw) || !nzchar(key_raw)) {
    message("Skip (cannot parse TF from filename): ", f)
    next
  }
  key_in_TFfamily <- resolve_key_tf_symbol(
    key_raw, key_TFs = key_TFs, CTS = CTS, CTS_ID = CTS_ID,
    graph_list = graph_list, CTS_name = CTS_name, CP_cluster = CP_cluster,
    TF_mouse = TF_mouse
  )

  res <- readRDS(f)
  g1 <- res[["g_CT_sub"]]
  g2 <- res[["g_descendant_sub"]]

  if (vcount(g1) > 0 && vcount(g2) > 0 && vcount(g1) == vcount(g2)) {
    change_df <- edge_change_table(g1, g2, weight_attr = "weight", missing_as = 0, undirected = TRUE)
    if (nrow(change_df) == 0) {
      message("No edge changes for ", key_in_TFfamily, " - ", pull)
      next
    }
    prioritize_edge_change(
      g1, edge_change_df = change_df, top_n = 5,
      title = paste0(pull, "-pull subnetwork_", key_in_TFfamily)
    )
  } else {
    message("No edges in ", pull, "-pull subnetwork for ", key_in_TFfamily)
    next
  }

  y <- motifAnnot_sub[
    motifAnnot_sub$regulators %in% key_in_TFfamily |
      grepl(key_in_TFfamily, motifAnnot_sub$regulators, fixed = TRUE),
  ]
  motif_TF_highConf_val <- if (nrow(y)) y$motif_TF_highConf[1] else NA
  tmp <- subset(
    cisTarget.res,
    geneSet == CTS_ID & NES >= NES_threshold &
      (motif == motif_TF_highConf_val | TF_highConf == motif_TF_highConf_val)
  )
  tmp1 <- if (nrow(tmp)) tmp[1, , drop = FALSE] else data.frame(TF_highConf = NA, motif = NA, NES = NA)
  change_df <- cbind(
    linkage = pull, change_df,
    TF_highConf = tmp1$TF_highConf, motif = tmp1$motif, NES = tmp1$NES
  )
  ns <- change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily
  change_df$TF_highConf[ns] <- ""
  change_df$motif[ns] <- ""
  change_df$NES[ns] <- ""
  final_table <- rbind(final_table, change_df)
}

if (is.null(final_table) || !nrow(final_table)) {
  if (exists("tips_diagnose_dualpull_failure", mode = "function")) {
    tips_diagnose_dualpull_failure(mat, key_TFs, x, graph_list, CM_cluster, CF_cluster)
  }
  stop("empty final_table — review key_TFs and CM/CF candidate colSums", call. = FALSE)
}

write.table(
  final_table,
  paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)
message("[24.1] saved: PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv")
message("[24.1] next: source('24.3_acat_CTS.cisTarget_merge_GRN.R')")
