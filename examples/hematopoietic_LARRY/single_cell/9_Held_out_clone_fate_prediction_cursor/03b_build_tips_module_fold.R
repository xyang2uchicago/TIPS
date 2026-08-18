## 03b_build_tips_module_fold.R
## Re-invoke unmodified code_core_11_10vs17 (11.1 → 24.3) on training-clone cells.
## Frozen: leiden_r0_8 C11 and BioTIP CTS. Refit: DEG, STRING reweight, PageRank, cisTarget dual-pull.
## A script error or missing/stale dual-pull/graph stops that fold (status=error).
## An empty module after a successful run is status=empty_module (scientific outcome).

heldout_code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor"
)
if (!exists("heldout_cfg")) {
  source(file.path(heldout_code_dir, "00_configuration.R"))
  source(file.path(heldout_code_dir, "00_helpers.R"))
}
cfg <- heldout_ensure()

empty_module <- function() {
  list(genes = character(), genes_lean = character(), edge_table = data.frame())
}

extract_fold_module <- function(results_dir, cts_id, linkage) {
  f <- list.files(
    file.path(results_dir, paste0("cisTarget_predicted_", cts_id)),
    pattern = "dualpull_final_table\\.tsv$", full.names = TRUE
  )
  if (!length(f)) return(empty_module())
  tab <- utils::read.delim(f[1], stringsAsFactors = FALSE)
  other <- if (identical(linkage, "CF")) "CM" else "CF"
  list(
    genes = norm_genes(with(tab[tab$linkage == linkage, ], c(from, to))),
    genes_lean = linkage_set(tab, linkage, other),
    edge_table = tab[tab$linkage == linkage, , drop = FALSE]
  )
}

extract_graph_genes <- function(graph_rds, cluster_id) {
  if (!file.exists(graph_rds)) return(character())
  gl <- readRDS(graph_rds)
  key <- intersect(
    c(paste0("HiG_", cluster_id), paste0("HiGCTS_", cluster_id), as.character(cluster_id)),
    names(gl)
  )
  if (!length(key)) return(character())
  norm_genes(igraph::V(gl[[key[1]]])$name)
}

redirect_fold_paths <- function(fold_root, cfg) {
  results_dir <- file.path(fold_root, paste0("results_core_", cfg$tips_arm_tag))
  data_dir <- file.path(fold_root, "data")
  ppi_path <- file.path(results_dir, "PPI_weight")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ppi_path, recursive = TRUE, showWarnings = FALSE)
  cache <- cfg$string_id_cache
  dest <- file.path(data_dir, "unique_STRING_mapping_correction.txt")
  if (file.exists(cache) && !file.exists(dest)) file.copy(cache, dest)
  list(results_dir = results_dir, data_dir = data_dir, ppi_path = ppi_path)
}

assert_fold_globals <- function(fold_seurat, paths, cfg) {
  g <- .GlobalEnv
  Sys.setenv(SEURAT_RDS = fold_seurat)
  assign("seurat_rds", fold_seurat, envir = g)
  assign("code_dir", paste0(cfg$tips_arm_code_dir, "/"), envir = g)
  assign("results_dir", paths$results_dir, envir = g)
  assign("data_dir", paths$data_dir, envir = g)
  assign("ppi_path", paths$ppi_path, envir = g)
  deg <- tips_deg_paths(get("logFC.cut", g), get("fdr.cut", g), get("min_prop", g), paths$data_dir)
  assign("deg_rdata", deg$deg_rdata, envir = g)
  assign("deg_file", deg$deg_rdata, envir = g)
  assign("markers_rdata", deg$markers_rdata, envir = g)
  assign("markers_file", deg$markers_rdata, envir = g)
  db <- get("db", g)
  assign("graph_notsimp_file", file.path(paths$results_dir, paste0(db, "_STRING_graph_perState_notsimplified.rds")), envir = g)
  assign("graph_weighted_file", file.path(paths$ppi_path, paste0(db, "_STRING_graph_perState_simplified_combinedweighted.rds")), envir = g)
  assign("correct_edges_file", file.path(paths$results_dir, "correct_n_edges_HiG_STRING2.14.0.rds"), envir = g)
}

run_tips_core_on_fold <- function(fold_seurat, fold_root, cfg, t0) {
  if (!dir.exists(cfg$tips_arm_code_dir)) stop("Missing TIPS arm: ", cfg$tips_arm_code_dir)
  old <- Sys.getenv(c("SEURAT_RDS", "BIOTIP_CTS_RDATA", "BIOTIP_CTS_SUMMARY_CSV",
                      "TIPS_GROUP_COL", "TIPS_CP_CLUSTER", "TIPS_CM_CLUSTER", "TIPS_CF_CLUSTER",
                      "TIPS_WD", "TIPS_ALLOW_CELL_LEVEL"))
  oldwd <- getwd()
  on.exit({
    do.call(Sys.setenv, as.list(old))
    setwd(oldwd)
  }, add = TRUE)
  Sys.setenv(
    SEURAT_RDS = fold_seurat,
    BIOTIP_CTS_RDATA = cfg$cell_cts_rdata,
    BIOTIP_CTS_SUMMARY_CSV = cfg$cell_cts_summary,
    TIPS_GROUP_COL = cfg$state_col,
    TIPS_CP_CLUSTER = cfg$state_id,
    TIPS_CM_CLUSTER = cfg$cm_cluster,
    TIPS_CF_CLUSTER = cfg$cf_cluster,
    TIPS_WD = paste0(cfg$tips_arm_wd, "/"),
    TIPS_ALLOW_CELL_LEVEL = "TRUE"
  )
  assign("TIPS_CONFIGURED", FALSE, envir = .GlobalEnv)
  source(file.path(cfg$tips_arm_code_dir, "00_configuration.R"), local = FALSE)
  tips_configure(
    TAG = cfg$tips_arm_tag, CP_cluster = cfg$state_id,
    CM_cluster = cfg$cm_cluster, CF_cluster = cfg$cf_cluster,
    group_col = cfg$state_col, db_species = 10090L,
    Real_names = c(CP = "progenitor", CM = "Baso", CF = "Meg"),
    motif_target_strategy = "merge", wd = paste0(cfg$tips_arm_wd, "/")
  )
  paths <- redirect_fold_paths(fold_root, cfg)
  assert_fold_globals(fold_seurat, paths, cfg)
  assign("rebuild_mat", TRUE, envir = .GlobalEnv)
  assign("rebuild_heatmap", TRUE, envir = .GlobalEnv)
  steps <- c(
    "mutrans_sce_xy.R", "11.1_STRINGweighted_CTS_network.R",
    "11.1.1_check_vertex_duplication.R", "11.2.0_update_network_weights_clean_max.R",
    "11.3_CTS_network_ANND_pagerank.R", "12.0_rank_by_PageRank_BC.R",
    "24.1_acat_CTS.cisTarget_dualpull_clean.R", "24.3_acat_CTS.cisTarget_merge_GRN.R"
  )
  for (s in steps) {
    f <- file.path(cfg$tips_arm_code_dir, s)
    if (!file.exists(f)) stop("Missing essential TIPS step: ", s)
    message("[03b] sourcing ", s)
    assert_fold_globals(fold_seurat, paths, cfg)
    source(f, local = FALSE)
  }
  dual <- list.files(
    file.path(paths$results_dir, paste0("cisTarget_predicted_", cfg$state_id)),
    pattern = "dualpull_final_table\\.tsv$", full.names = TRUE
  )
  g_rds <- file.path(paths$results_dir, "GSE140802_STRING_graph_perState_notsimplified.rds")
  floor_t <- t0 - 2
  if (!length(dual) || !file.exists(dual[1])) stop("dual-pull table missing after TIPS")
  if (file.info(dual[1])$mtime < floor_t) stop("dual-pull table is stale (not generated this fold)")
  if (!file.exists(g_rds)) stop("STRING graph RDS missing after TIPS: ", g_rds)
  if (file.info(g_rds)$mtime < floor_t) stop("STRING graph RDS is stale (not generated this fold)")
  paths$dual_table <- dual[1]
  paths$graph_rds <- g_rds
  paths
}

write_fold_status <- function(man, cfg) {
  f <- file.path(cfg$results_dir, "tables", "03b_fold_status.tsv")
  row <- data.frame(
    fold_tag = man$fold_tag, fold_id = man$fold_id, repeat_id = man$repeat_id,
    status = man$status, error = man$error,
    n_meg_genes = length(man$module$CF$genes_lean),
    n_baso_genes = length(man$module$CM$genes_lean),
    dual_table = man$dual_table, graph_rds = man$graph_rds,
    elapsed_mins = man$elapsed_mins,
    stringsAsFactors = FALSE
  )
  if (file.exists(f)) {
    old <- utils::read.delim(f, stringsAsFactors = FALSE)
    old <- old[as.character(old$fold_tag) != man$fold_tag, , drop = FALSE]
    row <- rbind(old, row)
  }
  write_tsv(row, f)
}

build_tips_module_for_fold <- function(fold_id, repeat_id, test_clone_ids, cfg, work_meta, clone_mat) {
  tag <- sprintf("fold%d_rep%d", fold_id, repeat_id)
  fold_root <- file.path(cfg$fold_modules_dir, tag)
  man_f <- file.path(fold_root, "fold_manifest.rds")
  if (file.exists(man_f) && !isTRUE(cfg$overwrite_fold_modules)) {
    man <- readRDS(man_f)
    if (fold_cache_ok(man, test_clone_ids, fold_id, repeat_id, cfg)) {
      message("[03b] reuse ", man_f, " (cache identity matched)")
      return(man)
    }
    message("[03b] ignoring stale cache ", man_f)
  }
  dir.create(fold_root, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  man <- list(
    fold_tag = tag, fold_id = fold_id, repeat_id = repeat_id,
    test_clone_ids = sort(as.integer(test_clone_ids)),
    tips_mode = cfg$tips_mode,
    state_id = cfg$state_id, cf_cluster = cfg$cf_cluster, cm_cluster = cfg$cm_cluster,
    lock_fingerprint = heldout_lock_fingerprint(cfg),
    code_hash = heldout_code_hash(cfg),
    module = list(CF = empty_module(), CM = empty_module()),
    graph_notsimp = list(CF = character(), CM = character()),
    deg = list(), status = "error", error = "",
    dual_table = "", graph_rds = "", elapsed_mins = NA_real_
  )
  man <- tryCatch({
    seu_path <- if (file.exists(cfg$seurat_rds_tips_arm)) cfg$seurat_rds_tips_arm else cfg$seurat_rds
    suppressPackageStartupMessages(library(Seurat))
    seu <- readRDS(seu_path)
    if (!"cell_index0" %in% names(seu@meta.data)) {
      seu$cell_index0 <- work_meta$cell_index0[match(colnames(seu), work_meta$cell_id)]
    }
    idx0 <- as.integer(seu$cell_index0)
    if (mean(is.na(idx0)) > 0.5) stop("cell_index0 missing on training Seurat: ", seu_path)
    is_test <- cells_in_clones(clone_mat, idx0, test_clone_ids)
    is_multi <- clone_hits(clone_mat, idx0)$multi
    n_before <- ncol(seu)
    n_test <- sum(is_test, na.rm = TRUE)
    n_multi <- sum(is_multi, na.rm = TRUE)
    keep <- colnames(seu)[!(is_test | is_multi)]
    seu_train <- subset(seu, cells = keep)
    idx_train <- as.integer(seu_train$cell_index0)
    n_remain <- sum(cells_in_clones(clone_mat, idx_train, test_clone_ids))
    n_multi_remain <- sum(clone_hits(clone_mat, idx_train)$multi)
    audit <- data.frame(
      fold_tag = tag, n_cells_before = n_before,
      n_test_clone_cells_identified = n_test,
      n_multi_clone_cells_removed = n_multi,
      n_removed = n_before - ncol(seu_train),
      n_cells_after = ncol(seu_train),
      n_test_clone_cells_remaining_in_training_object = n_remain,
      n_multi_clone_cells_remaining = n_multi_remain,
      stringsAsFactors = FALSE
    )
    write_tsv(audit, file.path(fold_root, "exclusion_audit.tsv"))
    af <- file.path(cfg$results_dir, "tables", "03b_exclusion_audit.tsv")
    if (file.exists(af)) {
      olda <- utils::read.delim(af, stringsAsFactors = FALSE)
      olda <- olda[as.character(olda$fold_tag) != tag, , drop = FALSE]
      audit <- rbind(olda, audit)
    }
    write_tsv(audit, af)
    if (n_remain != 0L) stop("test-clone cells remaining in training object: ", n_remain)
    message("[03b] ", tag, ": keep ", ncol(seu_train), "/", n_before,
            " (removed test=", n_test, " multi=", n_multi, ")")
    fold_seurat <- file.path(fold_root, "seurat_train.rds")
    saveRDS(seu_train, fold_seurat)
    rm(seu, seu_train); gc()
    paths <- run_tips_core_on_fold(fold_seurat, fold_root, cfg, t0)
    cf <- extract_fold_module(paths$results_dir, cfg$state_id, "CF")
    cm <- extract_fold_module(paths$results_dir, cfg$state_id, "CM")
    deg_f <- get0("deg_rdata", envir = .GlobalEnv, ifnotfound = "")
    man$seurat_train_rds <- fold_seurat
    man$results_dir <- paths$results_dir
    man$dual_table <- paths$dual_table
    man$graph_rds <- paths$graph_rds
    man$module <- list(CF = cf, CM = cm)
    man$graph_notsimp <- list(
      CF = extract_graph_genes(paths$graph_rds, cfg$cf_cluster),
      CM = extract_graph_genes(paths$graph_rds, cfg$cm_cluster)
    )
    man$deg <- if (nzchar(deg_f) && file.exists(deg_f)) lapply(readRDS(deg_f), function(x) norm_genes(x)) else list()
    man$n_excluded <- n_test + n_multi
    man$elapsed_mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    man$status <- if (!length(cf$genes_lean) && !length(cm$genes_lean)) "empty_module" else "success"
    man$error <- ""
    man
  }, error = function(e) {
    man$elapsed_mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    man$status <- "error"
    man$error <- conditionMessage(e)
    man
  })
  saveRDS(man, man_f)
  write_fold_status(man, cfg)
  message("[03b] ", tag, " status=", man$status,
          " CF n=", length(man$module$CF$genes_lean),
          " CM n=", length(man$module$CM$genes_lean),
          if (nzchar(man$error)) paste0(" error: ", man$error) else "")
  man
}

build_all_fold_modules <- function(cfg, lock, work_meta) {
  clone_mat <- read_clone_matrix(cfg$public_clone)
  primary <- lock$primary
  fold_list <- lock$fold_list
  if (is.null(fold_list) || !length(fold_list)) fold_list <- list(primary$fold)
  n_repeats <- if (!is.null(lock$n_repeats)) lock$n_repeats else 1L
  out <- list()
  for (r in seq_len(n_repeats)) {
    primary$fold <- fold_list[[r]]
    for (f in sort(unique(primary$fold))) {
      key <- paste(r, f, sep = "_")
      test_clones <- primary$clone_id[primary$fold == f]
      out[[key]] <- build_tips_module_for_fold(f, r, test_clones, cfg, work_meta, clone_mat)
    }
  }
  saveRDS(out, file.path(cfg$results_dir, "rds", "03b_fold_manifests.rds"))
  out
}

if (!isTRUE(get0("heldout_source_03b_quiet", ifnotfound = FALSE))) {
  lock_rds <- file.path(cfg$results_dir, "rds", "02_locked_clones_folds.rds")
  ds_rds <- file.path(cfg$results_dir, "rds", "01_clone_dataset.rds")
  if (file.exists(lock_rds) && file.exists(ds_rds)) {
    message("[03b] building TIPS modules for all folds (real code_core_11_10vs17)")
    build_all_fold_modules(cfg, readRDS(lock_rds), readRDS(ds_rds)$work_meta)
  } else {
    message("[03b] functions loaded; run 01+02 first to build folds")
  }
}
