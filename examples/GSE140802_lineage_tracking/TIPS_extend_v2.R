## =============================================================================
## TIPS_extend_v2.R — full TIPS helper library (DEFAULT — source this file)
## =============================================================================
##
## Superset of TIPS_extend.R:
##   Part 1 — identical dual-pull helpers (24.1 setup)
##   Part 2 — NEW lineage-lean gene modules (post dualpull_final_table.tsv)
##
## Versioning:
##   source(".../TIPS_extend_v2.R")  <- all functions; new downstream scripts
##   source(".../TIPS_extend.R")     <- Part 1 only (Larry 24.1)
##   TIPS_extend_v3.R (future)
##
## Override: Sys.setenv(TIPS_EXTEND_R = "/path/to/TIPS_extend_v2.R")
## =============================================================================

TIPS_EXTEND_VERSION <- "v2"


# =============================================================================
# PART 1 — Dual-pull candidate assignment and key-TF selection
# (same functions as TIPS_extend.R)
# =============================================================================
#
# PURPOSE: Run inside Felix 24.1 BEFORE dualpull_final_table.tsv exists.
# Assigns CP/CM/CF candidate genes per key TF from cisTarget heatmaps, with
# HiG-subgraph fallback when CM_hi/CF_hi rows are empty (e.g. code_core_1_5vs7).
# Used by: 24.1_acat_CTS.cisTarget_dualpull_clean.R (GSE140802 Larry arms).
# =============================================================================

#' Motif-hit mask for one cisTarget column
#'
#' @description Returns a logical vector marking genes that are motif targets of
#'   the key TF in column \code{col_idx}.
#' @param mat cisTarget heatmap matrix (genes x motifs/TFs).
#' @param col_idx Column index for one key TF.
#' @return Logical vector of length \code{nrow(mat)}.
tips_col_is_hit <- function(mat, col_idx) {
  if (is.na(col_idx) || col_idx < 1L) return(rep(FALSE, nrow(mat)))
  v <- mat[, col_idx]
  v == 1 | v == 1L | v == "1"
}

#' Vertex names from a hierarchical interaction graph (HiG)
#'
#' @description Genes present in \code{HiG_<cluster_id>} for arm-local fallback.
#' @param graph_list Named list of igraph objects (\code{HiG_<cluster>}).
#' @param cluster_id Cluster id without the \code{HiG_} prefix.
#' @return Uppercase gene symbols in that HiG.
tips_hig_vertices <- function(graph_list, cluster_id) {
  g <- graph_list[[paste0("HiG_", cluster_id)]]
  if (is.null(g) || !inherits(g, "igraph")) return(character())
  toupper(as.character(igraph::V(g)$name))
}

#' Assign CP/CM/CF dual-pull candidates for one key TF
#'
#' @description Felix-first assignment using \code{CP_hi}, \code{CM_hi},
#'   \code{CF_hi}. When CM or CF candidates are empty, intersects motif targets
#'   with genes in the corresponding HiG arm subgraph.
#' @param mat cisTarget matrix with \code{CP_hi}, \code{CM_hi}, \code{CF_hi}.
#' @param key Key TF name (column prefix for candidate columns).
#' @param col_idx cisTarget column index for \code{key}.
#' @param graph_list List of HiG igraph objects.
#' @param CM_cluster,CF_cluster Cluster ids for blood (CM) and alternate (CF) arms.
#' @return \code{mat} with \code{<key>_CP_candidate}, \code{_CM_candidate},
#'   \code{_CF_candidate} columns added or updated.
tips_assign_dualpull_candidates <- function(
    mat, key, col_idx, graph_list, CM_cluster, CF_cluster
) {
  key <- as.character(key)
  cp_col <- paste0(key, "_CP_candidate")
  cm_col <- paste0(key, "_CM_candidate")
  cf_col <- paste0(key, "_CF_candidate")
  hit <- tips_col_is_hit(mat, col_idx)

  mat[[cp_col]] <- ifelse((mat[, "CP_hi"] == 1 | mat[, "CP_hi"] == 1L) & hit, 1, 0)
  mat[[cm_col]] <- ifelse((mat[, "CM_hi"] == 1 | mat[, "CM_hi"] == 1L) & hit, 1, 0)
  mat[[cf_col]] <- ifelse((mat[, "CF_hi"] == 1 | mat[, "CF_hi"] == 1L) & hit, 1, 0)

  motif_genes <- toupper(rownames(mat)[hit])
  hig_cm <- tips_hig_vertices(graph_list, CM_cluster)
  hig_cf <- tips_hig_vertices(graph_list, CF_cluster)

  if (!sum(mat[[cm_col]]) && length(motif_genes) && length(hig_cm)) {
    cm_hit <- intersect(motif_genes, hig_cm)
    if (length(cm_hit)) {
      mat[[cm_col]] <- ifelse(toupper(rownames(mat)) %in% cm_hit, 1, 0)
      message(
        "[TIPS_extend] ", key, ": CM_hi empty — CM candidates from motif targets in HiG_",
        CM_cluster, " (n=", length(cm_hit), ")"
      )
    }
  }
  if (!sum(mat[[cf_col]]) && length(motif_genes) && length(hig_cf)) {
    cf_hit <- intersect(motif_genes, hig_cf)
    if (length(cf_hit)) {
      mat[[cf_col]] <- ifelse(toupper(rownames(mat)) %in% cf_hit, 1, 0)
      message(
        "[TIPS_extend] ", key, ": CF_hi empty — CF candidates from motif targets in HiG_",
        CF_cluster, " (n=", length(cf_hit), ")"
      )
    }
  }
  mat
}

#' Refresh dual-pull candidates for all key TFs
#'
#' @description Loops over named key TFs and calls
#'   \code{\link{tips_assign_dualpull_candidates}} for each.
#' @inheritParams tips_assign_dualpull_candidates
#' @param x Named integer vector mapping TF name to cisTarget column index.
#' @return Updated cisTarget matrix \code{mat}.
tips_refresh_dualpull_candidates <- function(mat, x, graph_list, CM_cluster, CF_cluster) {
  for (j in seq_along(x)) {
    mat <- tips_assign_dualpull_candidates(
      mat, names(x)[j], unname(x[[j]]), graph_list, CM_cluster, CF_cluster
    )
  }
  mat
}

#' Count CM/CF dual-pull candidates for one key TF
#'
#' @param mat cisTarget matrix with \code{<key>_CM_candidate} columns.
#' @param key Key TF name.
#' @return Named integer vector \code{c(CM = n, CF = n)}.
tips_dualpull_candidate_sums <- function(mat, key) {
  key <- as.character(key)
  cm_col <- paste0(key, "_CM_candidate")
  cf_col <- paste0(key, "_CF_candidate")
  c(
    CM = if (cm_col %in% colnames(mat)) sum(mat[[cm_col]]) else 0L,
    CF = if (cf_col %in% colnames(mat)) sum(mat[[cf_col]]) else 0L
  )
}

#' Test whether a key TF has dual-pull candidate support
#'
#' @description TRUE when the TF has at least one CM or CF candidate gene.
#' @inheritParams tips_dualpull_candidate_sums
#' @return Logical scalar.
tips_key_tf_has_dualpull_support <- function(mat, key) {
  s <- tips_dualpull_candidate_sums(mat, key)
  s["CM"] > 0 || s["CF"] > 0
}

#' Summary table of dual-pull candidate counts per key TF
#'
#' @param mat cisTarget matrix with candidate columns.
#' @param keys Character vector of key TF names.
#' @return Data frame with columns \code{TF}, \code{CM_candidates},
#'   \code{CF_candidates}.
tips_dualpull_candidate_report <- function(mat, keys) {
  keys <- as.character(keys)
  out <- lapply(keys, function(k) {
    s <- tips_dualpull_candidate_sums(mat, k)
    data.frame(TF = k, CM_candidates = unname(s["CM"]), CF_candidates = unname(s["CF"]),
               row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' PageRank walk to pick key TFs with dual-pull support
#'
#' @description Walks PageRank-ordered TFs and keeps the first
#'   \code{top_TF_rank} TFs that have CM or CF candidates after assignment.
#' @param candidate_tfs Ordered TF names (typically PageRank).
#' @param mat cisTarget matrix.
#' @param col_idx_for_tf Function mapping TF name to cisTarget column index.
#' @param graph_list HiG graph list.
#' @param CM_cluster,CF_cluster Arm cluster ids.
#' @param top_TF_rank Maximum number of key TFs to return.
#' @return List with \code{x} (named column indices) and updated \code{mat}.
tips_extend_pick_key_tfs <- function(
    candidate_tfs, mat, col_idx_for_tf, graph_list, CM_cluster, CF_cluster,
    top_TF_rank = 3L
) {
  candidate_tfs <- toupper(as.character(candidate_tfs))
  candidate_tfs <- candidate_tfs[!is.na(candidate_tfs) & nzchar(candidate_tfs)]
  candidate_tfs <- unique(candidate_tfs)

  col_idx <- integer()
  tf_names <- character()
  for (tf in candidate_tfs) {
    j <- col_idx_for_tf(tf)
    if (is.na(j) || j %in% col_idx) next
    tmp_mat <- tips_assign_dualpull_candidates(
      mat, tf, j, graph_list, CM_cluster, CF_cluster
    )
    if (!tips_key_tf_has_dualpull_support(tmp_mat, tf)) {
      s <- tips_dualpull_candidate_sums(tmp_mat, tf)
      message(
        "[TIPS_extend] skip ", tf, " (CM_candidates=", s["CM"], ", CF_candidates=", s["CF"], ")"
      )
      next
    }
    col_idx <- c(col_idx, j)
    tf_names <- c(tf_names, tf)
    mat <- tmp_mat
    if (length(col_idx) >= top_TF_rank) break
  }
  list(x = stats::setNames(col_idx, tf_names), mat = mat)
}

#' Print dual-pull failure diagnostics
#'
#' @description Reports candidate counts and HiG sizes when 24.1 cannot build
#'   a dual-pull subnetwork.
#' @inheritParams tips_extend_pick_key_tfs
#' @param key_TFs Current key TF names under consideration.
#' @param x Named column-index vector for key TFs.
#' @return Invisible \code{NULL}.
tips_diagnose_dualpull_failure <- function(
    mat, key_TFs, x, graph_list, CM_cluster, CF_cluster
) {
  message("[TIPS_extend] dual-pull candidate report:")
  print(tips_dualpull_candidate_report(mat, key_TFs))
  message(
    "[TIPS_extend] HiG_", CM_cluster, " vertices: ", length(tips_hig_vertices(graph_list, CM_cluster)),
    "; HiG_", CF_cluster, " vertices: ", length(tips_hig_vertices(graph_list, CF_cluster))
  )
  message(
    "[TIPS_extend] Fix: set key_TFs_override in 24.1 USER INPUT, rebuild_heatmap <- TRUE,",
    " or tips_extend_auto_key_tfs <- TRUE"
  )
  invisible(NULL)
}

#' Felix-first key-TF assignment with optional PageRank re-pick
#'
#' @description Refreshes candidates for current key TFs; if none have support
#'   and \code{auto = TRUE}, re-picks from PageRank-ordered \code{pr_with_col}.
#' @param x Named column indices for current key TFs.
#' @param pr_with_col PageRank-ordered TF names for auto re-pick.
#' @param col_idx_for_tf Function mapping TF name to column index.
#' @param graph_list HiG graph list.
#' @param CM_cluster,CF_cluster Arm cluster ids.
#' @param top_TF_rank Max key TFs when re-picking.
#' @param auto Enable PageRank re-pick when current TFs lack support.
#' @param key_TFs_override Manual TF list; disables auto re-pick when set.
#' @return List \code{list(x, mat, key_TFs)}.
tips_maybe_extend_key_tfs <- function(
    x, mat, pr_with_col, col_idx_for_tf, graph_list, CM_cluster, CF_cluster,
    top_TF_rank = 3L, auto = TRUE, key_TFs_override = NULL
) {
  mat <- tips_refresh_dualpull_candidates(mat, x, graph_list, CM_cluster, CF_cluster)
  key_TFs <- names(x)
  need_pick <- !length(key_TFs) ||
    all(!vapply(key_TFs, tips_key_tf_has_dualpull_support, logical(1L), mat = mat))
  if (isTRUE(auto) && is.null(key_TFs_override) && need_pick && length(pr_with_col)) {
    message("[TIPS_extend] re-picking key_TFs with dual-pull support")
    pick <- tips_extend_pick_key_tfs(
      pr_with_col, mat, col_idx_for_tf, graph_list, CM_cluster, CF_cluster,
      top_TF_rank = top_TF_rank
    )
    if (length(pick$x)) {
      x <- pick$x
      mat <- pick$mat
      key_TFs <- names(x)
      message("[TIPS_extend] key_TFs: ", paste(key_TFs, collapse = ", "))
    }
  }
  list(x = x, mat = mat, key_TFs = names(x))
}


# =============================================================================
# PART 2 — Lineage-lean gene modules from dual-pull final tables  [NEW in v2]
# =============================================================================
#
# PURPOSE: Run AFTER PPI_graph_GRN_prediction_*_dualpull_final_table.tsv exists.
# Builds CM/CF gene modules for UCell scoring, cross-atlas heatmaps, and AUROC.
#
# Manuscript terms:
#   lineage-biased edge  = any gene on CM or CF arm (module_set = "full")
#   lineage-lean gene    = arm-exclusive genes plus bridge reclaim:
#     (a) tips_reclaim_bridge — opposite Delta-w on bridge gene or shared edge
#     (b) tips_reclaim_cf_increase_only — CF increase majority, no CM increase
#         (e.g. HAND2, FOXA2; excludes TAL1 with one CF increase edge)
#
# Used by: HEP_cross_atlas_transfer_Ucell.R, cross-atlas heatmaps.
# Entry points: tips_lineage_lean_genes(), tips_lineage_lean_modules()
# =============================================================================


# --- 2a. Path and table I/O ---------------------------------------------------

#' Normalize gene symbols to uppercase unique values
#'
#' @param x Character vector of gene symbols (any case; may contain NA/blank).
#' @return Character vector of unique uppercase non-empty symbols.
tips_norm_sym <- function(x) {
  toupper(unique(na.omit(as.character(x))[nzchar(as.character(x))]))
}

#' Resolve path to a dual-pull final table
#'
#' @description Constructs the standard TIPS path to
#'   \code{PPI_graph_GRN_prediction_CTS_<id>_dualpull_final_table.tsv}.
#' @param results_wd Root results directory (e.g. \code{results_core_13}).
#' @param cts_id CTS cluster id (e.g. \code{"13"}, \code{"endothelial.b"}).
#' @param nes_cut NES threshold used in the run (e.g. \code{3}).
#' @param has_nes45 If \code{TRUE} and \code{nes_cut >= 4.5}, use the
#'   \code{_NES4.5} cisTarget subdirectory.
#' @return Character file path.
tips_dualpull_path <- function(results_wd, cts_id, nes_cut, has_nes45 = FALSE) {
  sub <- if (isTRUE(has_nes45) && nes_cut >= 4.5) {
    paste0("cisTarget_predicted_", cts_id, "_NES4.5")
  } else {
    paste0("cisTarget_predicted_", cts_id)
  }
  file.path(results_wd, sub, paste0("PPI_graph_GRN_prediction_CTS_", cts_id, "_dualpull_final_table.tsv"))
}

#' Read and split a dual-pull final table into CM/CF arms
#'
#' @param path Path to \code{*_dualpull_final_table.tsv}.
#' @return List with \code{cm_df}, \code{cf_df}, \code{cm_all}, \code{cf_all},
#'   and \code{bridge} (genes on both arms). Empty components if file missing.
tips_read_dualpull <- function(path) {
  empty <- list(
    cm_df = data.frame(), cf_df = data.frame(),
    cm_all = character(), cf_all = character(), bridge = character()
  )
  if (!file.exists(path)) return(empty)
  df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("linkage", "from", "to", "direction")
  if (!all(req %in% colnames(df))) {
    stop("Unexpected dualpull columns in ", path, ": need ", paste(req, collapse = ", "))
  }
  df <- df[df$direction %in% c("increase", "decrease"), , drop = FALSE]
  cm_df <- df[df$linkage == "CM", , drop = FALSE]
  cf_df <- df[df$linkage == "CF", , drop = FALSE]
  cm_all <- tips_norm_sym(c(cm_df$from, cm_df$to))
  cf_all <- tips_norm_sym(c(cf_df$from, cf_df$to))
  list(cm_df = cm_df, cf_df = cf_df, cm_all = cm_all, cf_all = cf_all, bridge = intersect(cm_all, cf_all))
}


# --- 2b. Edge and gene direction helpers --------------------------------------

#' Canonical undirected edge key for one PPI/GRN edge
#'
#' @param from,to Gene symbols at the two ends of an edge.
#' @return Character \code{"GENE1|GENE2"} with genes sorted and uppercased.
tips_dualpull_edge_key <- function(from, to) {
  paste(sort(tips_norm_sym(c(from, to))), collapse = "|")
}

#' Direction classes for one gene on one dual-pull arm
#'
#' @param df Rows from dualpull table filtered to \code{linkage == "CM"} or CF.
#' @param gene Gene symbol.
#' @return Character vector of \code{"increase"} and/or \code{"decrease"}, or empty.
tips_gene_dirs <- function(df, gene) {
  g <- toupper(gene)
  i <- which(toupper(df$from) == g | toupper(df$to) == g)
  if (!length(i)) character() else unique(df$direction[i])
}

#' Collapse duplicate edges to a single direction per arm
#'
#' @description Edges with conflicting directions on the same arm are dropped.
#' @param df Dual-pull rows for one arm.
#' @return Data frame with columns \code{edge_key} and \code{direction}.
tips_collapse_edges <- function(df) {
  if (!nrow(df)) return(data.frame(edge_key = character(), direction = character()))
  df$edge_key <- mapply(tips_dualpull_edge_key, df$from, df$to, USE.NAMES = FALSE)
  out <- aggregate(direction ~ edge_key, df, function(x) {
    x <- unique(x[x %in% c("increase", "decrease")])
    if (length(x) == 1L) x else NA_character_
  })
  out[!is.na(out$direction), , drop = FALSE]
}

#' Test whether a gene participates in a given direction on an arm
#'
#' @inheritParams tips_gene_dirs
#' @param dir Direction to test: \code{"increase"} or \code{"decrease"}.
#' @return Logical scalar.
tips_has_gene_direction <- function(df, gene, dir) {
  g <- toupper(gene)
  i <- which(toupper(df$from) == g | toupper(df$to) == g)
  if (!length(i)) FALSE else any(df$direction[i] == dir)
}

#' Count increase and decrease edges touching one gene on one arm
#'
#' @inheritParams tips_gene_dirs
#' @return Named integer vector \code{c(increase = n, decrease = n)}.
tips_gene_edge_dir_counts <- function(df, gene) {
  g <- toupper(gene)
  i <- which(toupper(df$from) == g | toupper(df$to) == g)
  if (!length(i)) return(c(increase = 0L, decrease = 0L))
  c(
    increase = as.integer(sum(df$direction[i] == "increase")),
    decrease = as.integer(sum(df$direction[i] == "decrease"))
  )
}


# --- 2c. Bridge reclaim rules -------------------------------------------------

#' Reclaim bridge genes with opposite lineage directions
#'
#' @description Bridge genes appear on both CM and CF arms. Reclaim to CM when
#'   CM shows increase and CF decrease (gene-level unidirectional pattern), or
#'   when a shared edge is CM-increase / CF-decrease (edge-level). Symmetric
#'   rules apply for CF reclaim.
#' @param cm_df,cf_df Dual-pull rows for CM and CF linkages.
#' @param bridge Character vector of bridge gene symbols.
#' @return Named list \code{list(CM = <genes>, CF = <genes>)} of reclaimed symbols.
#' @seealso \code{\link{tips_reclaim_cf_increase_only}}
tips_reclaim_bridge <- function(cm_df, cf_df, bridge) {
  out <- list(CM = character(), CF = character())
  if (!length(bridge)) return(out)
  for (g in bridge) {
    cm_d <- tips_gene_dirs(cm_df, g)
    cf_d <- tips_gene_dirs(cf_df, g)
    if (length(cm_d) == 1L && cm_d == "increase" && length(cf_d) == 1L && cf_d == "decrease") {
      out$CM <- c(out$CM, toupper(g))
    }
    if (length(cm_d) == 1L && cm_d == "decrease" && length(cf_d) == 1L && cf_d == "increase") {
      out$CF <- c(out$CF, toupper(g))
    }
  }
  cm_e <- tips_collapse_edges(cm_df)
  cf_e <- tips_collapse_edges(cf_df)
  if (nrow(cm_e) && nrow(cf_e)) {
    m <- merge(cm_e, cf_e, by = "edge_key", suffixes = c("_cm", "_cf"))
    for (i in seq_len(nrow(m))) {
      genes <- intersect(tips_norm_sym(strsplit(m$edge_key[i], "|", fixed = TRUE)[[1]]), bridge)
      if (!length(genes)) next
      if (m$direction_cm[i] == "increase" && m$direction_cf[i] == "decrease") {
        out$CM <- union(out$CM, genes)
      }
      if (m$direction_cm[i] == "decrease" && m$direction_cf[i] == "increase") {
        out$CF <- union(out$CF, genes)
      }
    }
  }
  out
}

#' Reclaim CF-biased bridge TFs with CF increase edge majority
#'
#' @description For bridge genes with mixed CF directions but no CM increase,
#'   reclaim to CF when CF increase edges outnumber CF decrease edges.
#'   Captures endothelial-biased regulators (HAND2, FOXA2) while excluding
#'   genes like TAL1 that have only one CF increase among mostly decrease edges.
#' @inheritParams tips_reclaim_bridge
#' @return Character vector of reclaimed CF gene symbols.
#' @seealso \code{\link{tips_reclaim_bridge}}
tips_reclaim_cf_increase_only <- function(cm_df, cf_df, bridge) {
  if (!length(bridge)) return(character())
  out <- character()
  for (g in bridge) {
    if (tips_has_gene_direction(cm_df, g, "increase")) next
    cf_n <- tips_gene_edge_dir_counts(cf_df, g)
    if (cf_n["increase"] > cf_n["decrease"]) out <- c(out, toupper(g))
  }
  unique(out)
}


# --- 2d. Public API: lineage-lean module builders -----------------------------

#' Build lineage-lean or lineage-biased gene sets from a dual-pull table
#'
#' @description Main entry point for Part 2. Supports three module modes:
#' \describe{
#'   \item{\code{full}}{All genes on each arm (lineage-biased-edge modules).}
#'   \item{\code{lean_exclusive_only}}{Arm-exclusive genes only; no bridge reclaim.}
#'   \item{\code{lean}}{Arm-exclusive plus \code{\link{tips_reclaim_bridge}} and
#'     \code{\link{tips_reclaim_cf_increase_only}}.}
#' }
#' @param path Path to \code{*_dualpull_final_table.tsv}.
#' @param module_set One of \code{"lean"}, \code{"lean_exclusive_only"},
#'   \code{"full"}, or alias \code{"biased"} (\code{-> lean}).
#' @return List with:
#' \describe{
#'   \item{\code{CM}, \code{CF}}{Gene symbol vectors per arm.}
#'   \item{\code{bridge_genes}}{Genes present on both arms.}
#'   \item{\code{gene_set_type}}{Label for the module type.}
#'   \item{\code{gene_table}}{Per-gene audit: branch and membership reason.}
#'   \item{\code{n_*}}{Counts of arm-exclusive and reclaimed genes.}
#' }
#' @seealso \code{\link{tips_lineage_lean_modules}}
tips_lineage_lean_genes <- function(path, module_set = "lean") {
  empty <- list(
    CM = character(), CF = character(), bridge_genes = character(),
    gene_set_type = module_set, gene_table = NULL,
    n_cm_arm_exclusive = 0L, n_cf_arm_exclusive = 0L,
    n_cm_reclaim = 0L, n_cf_reclaim = 0L
  )
  dp <- tips_read_dualpull(path)
  if (!nrow(dp$cm_df) && !nrow(dp$cf_df)) return(empty)

  cm_all <- dp$cm_all
  cf_all <- dp$cf_all
  bridge <- dp$bridge
  cm_df <- dp$cm_df
  cf_df <- dp$cf_df

  if (identical(module_set, "full")) {
    return(modifyList(empty, list(
      CM = cm_all, CF = cf_all, bridge_genes = bridge,
      gene_set_type = "all_lineage_biased_edge_genes",
      n_cm_arm_exclusive = length(cm_all), n_cf_arm_exclusive = length(cf_all)
    )))
  }

  module_set <- if (identical(module_set, "biased")) "lean" else module_set
  if (!module_set %in% c("lean", "lean_exclusive_only")) {
    stop("module_set must be lean, lean_exclusive_only, full, or biased")
  }

  cm_ex <- setdiff(cm_all, bridge)
  cf_ex <- setdiff(cf_all, bridge)
  cm_rc <- cf_rc <- cf_inc_rc <- character()
  gene_set_type <- if (module_set == "lean") "lineage_lean" else "lean_exclusive_only"

  if (module_set == "lean") {
    rc <- tips_reclaim_bridge(cm_df, cf_df, bridge)
    cm_rc <- rc$CM
    cf_rc <- rc$CF
    cf_inc_rc <- tips_reclaim_cf_increase_only(cm_df, cf_df, bridge)
  }

  cf_all_rc <- unique(c(cf_rc, cf_inc_rc))
  mk <- function(genes, branch, membership) {
    if (length(genes)) data.frame(gene = genes, branch = branch, membership = membership) else NULL
  }
  gene_table <- do.call(rbind, list(
    mk(cm_ex, "CM", "arm_exclusive"),
    mk(cf_ex, "CF", "arm_exclusive"),
    mk(cm_rc, "CM", "opposite_direction_reclaim"),
    mk(cf_rc, "CF", "opposite_direction_reclaim"),
    mk(cf_inc_rc, "CF", "cf_increase_only_bridge")
  ))

  list(
    CM = unique(c(cm_ex, cm_rc)),
    CF = unique(c(cf_ex, cf_all_rc)),
    bridge_genes = bridge,
    gene_set_type = gene_set_type,
    gene_table = gene_table,
    n_cm_arm_exclusive = length(cm_ex),
    n_cf_arm_exclusive = length(cf_ex),
    n_cm_reclaim = length(cm_rc),
    n_cf_reclaim = length(cf_all_rc)
  )
}

#' Lineage-lean modules with Felix HEP arm names
#'
#' @description Thin wrapper around \code{\link{tips_lineage_lean_genes}} that
#'   maps CM to Blood and CF to Endothelium for cross-atlas UCell scripts.
#' @inheritParams tips_lineage_lean_genes
#' @return List \code{list(Blood = <CM genes>, Endothelium = <CF genes>)}.
tips_lineage_lean_modules <- function(path, module_set = "lean") {
  g <- tips_lineage_lean_genes(path, module_set = module_set)
  list(Blood = g$CM, Endothelium = g$CF)
}
