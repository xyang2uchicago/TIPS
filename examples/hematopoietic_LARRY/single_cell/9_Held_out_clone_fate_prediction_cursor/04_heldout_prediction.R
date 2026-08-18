## 04_heldout_prediction.R
## Outer clone-grouped CV:
##   - drop every cell of testing clones from TIPS / DEG construction
##   - freeze genes + scoring rule
##   - score only day-2 C11 cells of testing clones
##   - pool out-of-fold clone scores
##
## Later day-6 outcomes of test clones are not used to build modules.

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor"
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_ensure()

lock_rds <- file.path(cfg$results_dir, "rds", "02_locked_clones_folds.rds")
sets_rds <- file.path(cfg$results_dir, "rds", "03_gene_sets.rds")
ds_rds <- file.path(cfg$results_dir, "rds", "01_clone_dataset.rds")
if (!file.exists(lock_rds)) stop("Run 02 first")
if (!file.exists(sets_rds)) stop("Run 03 first")
if (!file.exists(cfg$seurat_rds)) stop("Missing Seurat: ", cfg$seurat_rds)

lock <- readRDS(lock_rds)
sets <- readRDS(sets_rds)
ds <- readRDS(ds_rds)
primary <- lock$primary
fold_list <- lock$fold_list
if (is.null(fold_list) || !length(fold_list)) fold_list <- list(primary$fold)
n_repeats <- if (!is.null(lock$n_repeats)) lock$n_repeats else length(fold_list)
if (n_repeats > 1L) {
  cfg$n_random <- min(cfg$n_random, 20L)
  message("[04] n_random=", cfg$n_random, " (capped because n_repeats=", n_repeats, ")")
}
work <- ds$work_meta
d2_cells <- ds$d2_c11_cells
d2_cells <- d2_cells[d2_cells$clone_id %in% primary$clone_id, ]

heldout_source_03b_quiet <- TRUE
fold_mans <- NULL
if (identical(cfg$tips_mode, "run_tips_pipeline")) {
  assign("heldout_source_03b_quiet", TRUE, envir = .GlobalEnv)
  source(file.path(code_dir, "03b_build_tips_module_fold.R"), local = TRUE)
  man_all <- file.path(cfg$results_dir, "rds", "03b_fold_manifests.rds")
  if (file.exists(man_all) && identical(Sys.getenv("HELDOUT_RUN"), "all")) {
    cfg$overwrite_fold_modules <- FALSE
  }
  message("[04] ensuring per-fold TIPS modules (reuse only if cache identity matches)")
  fold_mans <- build_all_fold_modules(cfg, lock, work)
}

candidate <- unique(c(
  sets$tips_meg_frozen, sets$tips_baso_frozen,
  sets$ppi_unweighted_meg, sets$cts_c11,
  sets$mutrans_td_meg, sets$mutrans_td_baso,
  sets$metacell_tips_meg, sets$hvg,
  unlist(sets$tips_cf_edges[, c("from", "to")], use.names = FALSE),
  unlist(sets$tips_cm_edges[, c("from", "to")], use.names = FALSE),
  if (!is.null(fold_mans)) unlist(lapply(fold_mans, function(m) {
    c(m$module$CF$genes, m$module$CF$genes_lean, m$module$CM$genes, m$module$CM$genes_lean,
      m$graph_notsimp$CF, m$graph_notsimp$CM, unlist(m$deg))
  }))
))

message("[04] loading Seurat (large)")
suppressPackageStartupMessages(library(Seurat))
seu <- readRDS(cfg$seurat_rds)
expr <- extract_seurat_expr(seu, candidate)
## keep a larger HVG panel for DEG / random matching
hvg_use <- match_genes_to_rownames(sets$hvg, rownames(seu))
if (length(setdiff(hvg_use, rownames(expr)))) {
  extra <- extract_seurat_expr(seu, setdiff(hvg_use, rownames(expr)))
  expr <- rbind(expr, extra)
}
meta <- seu@meta.data
meta$cell_id <- rownames(meta)
rm(seu)
gc()
message("[04] expression matrix ", nrow(expr), " genes x ", ncol(expr), " cells")

## map working clone_id onto Seurat cells
work_map <- setNames(work$clone_id, work$cell_id)
clone_of <- work_map[colnames(expr)]
state_of <- as.character(meta[[cfg$state_col]][match(colnames(expr), rownames(meta))])
names(state_of) <- colnames(expr)
type_of <- as.character(meta$Cell.type.annotation[match(colnames(expr), rownames(meta))])
names(type_of) <- colnames(expr)
time_of <- as_int_like(meta$Time.point[match(colnames(expr), rownames(meta))])
names(time_of) <- colnames(expr)

map_edge_genes <- function(edges, feature_names) {
  mp <- setNames(as.character(feature_names), toupper(as.character(feature_names)))
  edges$from_r <- unname(mp[toupper(edges$from)])
  edges$to_r <- unname(mp[toupper(edges$to)])
  edges[!is.na(edges$from_r) & !is.na(edges$to_r), ]
}

cf_edges <- map_edge_genes(sets$tips_cf_edges, rownames(expr))
cm_edges <- map_edge_genes(sets$tips_cm_edges, rownames(expr))

mean_expr <- function(genes, cells) {
  use <- match_genes_to_rownames(genes, rownames(expr))
  cells <- intersect(cells, colnames(expr))
  if (!length(use) || !length(cells)) return(setNames(numeric(), character()))
  rowMeans(as.matrix(expr[use, cells, drop = FALSE]), na.rm = TRUE)
}

wilcox_up_genes <- function(pos_cells, neg_cells, universe, n_keep) {
  pos_cells <- intersect(pos_cells, colnames(expr))
  neg_cells <- intersect(neg_cells, colnames(expr))
  universe <- match_genes_to_rownames(universe, rownames(expr))
  if (length(pos_cells) < 10L || length(neg_cells) < 10L || length(universe) < 5L) {
    return(character())
  }
  m <- as.matrix(expr[universe, c(pos_cells, neg_cells), drop = FALSE])
  npos <- length(pos_cells)
  p <- apply(m, 1, function(v) {
    a <- v[seq_len(npos)]; b <- v[(npos + 1):length(v)]
    if (!is.finite(stats::sd(c(a, b)))) return(1)
    suppressWarnings(stats::wilcox.test(a, b, alternative = "greater", exact = FALSE)$p.value)
  })
  lfc <- rowMeans(m[, seq_len(npos), drop = FALSE]) - rowMeans(m[, (npos + 1):ncol(m), drop = FALSE])
  keep <- names(sort(p[lfc > 0 & is.finite(p)]))
  head(keep, n_keep)
}

matched_random <- function(target_genes, universe, train_c11, n_draw, seed) {
  target_genes <- match_genes_to_rownames(target_genes, rownames(expr))
  universe <- setdiff(match_genes_to_rownames(universe, rownames(expr)), target_genes)
  if (!length(target_genes)) {
    return(replicate(n_draw, character(), simplify = FALSE))
  }
  if (length(universe) < length(target_genes)) {
    return(replicate(n_draw, sample(universe, length(universe)), simplify = FALSE))
  }
  mu_t <- mean_expr(target_genes, train_c11)
  mu_u <- mean_expr(universe, train_c11)
  set.seed(seed)
  lapply(seq_len(n_draw), function(i) {
    picked <- character()
    pool <- names(mu_u)
    for (g in names(mu_t)) {
      if (!length(pool)) break
      d <- abs(mu_u[pool] - mu_t[g])
      w <- which.min(d)
      picked <- c(picked, pool[w])
      pool <- pool[-w]
    }
    picked
  })
}

reweight_module <- function(edges, train_cells, positive_delta = TRUE) {
  cp <- intersect(train_cells, names(state_of)[state_of == as.character(cfg$state_id)])
  cf <- intersect(train_cells, names(state_of)[state_of == as.character(cfg$cf_cluster)])
  if (!nrow(edges)) return(character())
  w1 <- pair_pearson(expr, edges$from_r, edges$to_r, cp)
  w2 <- pair_pearson(expr, edges$from_r, edges$to_r, cf)
  delta <- w2 - w1
  keep <- if (positive_delta) which(is.finite(delta) & delta > 0) else which(is.finite(delta))
  if (!length(keep)) {
    message("[04] no positive-delta edges; using all CF-linkage nodes")
    keep <- which(is.finite(delta) | TRUE)
  }
  unique(c(edges$from[keep], edges$to[keep]))
}

score_clone_mean <- function(cell_scores, cell_ids, clone_ids, fun = cfg$score_fun) {
  df <- data.frame(clone_id = clone_ids, score = unname(cell_scores[cell_ids]), stringsAsFactors = FALSE)
  df <- df[is.finite(df$score), ]
  if (!nrow(df)) return(setNames(numeric(), character()))
  agg <- stats::aggregate(score ~ clone_id, df, if (identical(fun, "median")) stats::median else mean)
  setNames(agg$score, as.character(agg$clone_id))
}

fold_modules <- list()
pred_rows <- list()
rand_rows <- list()
success_rows <- list()

for (r in seq_len(n_repeats)) {
  primary$fold <- fold_list[[r]]
  message("[04] repeat ", r, "/", n_repeats)
  for (f in sort(unique(primary$fold))) {
  message("[04] repeat ", r, " fold ", f)
  test_clones <- primary$clone_id[primary$fold == f]
  train_clones <- primary$clone_id[primary$fold != f]
  test_cell <- d2_cells$cell_id[d2_cells$clone_id %in% test_clones]
  test_cell <- intersect(test_cell, colnames(expr))
  is_test_clone <- !is.na(clone_of) & clone_of %in% test_clones
  train_cells <- colnames(expr)[!is_test_clone]
  train_c11 <- intersect(train_cells, names(state_of)[state_of == as.character(cfg$state_id)])
  train_d6 <- intersect(train_cells, names(time_of)[time_of == cfg$day_late])
  train_meg <- intersect(train_d6, names(type_of)[type_of == cfg$meg_label])
  train_not_meg <- setdiff(train_d6, train_meg)

  ppi_meg <- sets$ppi_unweighted_meg
  ppi_baso <- sets$ppi_unweighted_baso
  fold_status <- "success"
  fold_error <- ""
  if (identical(cfg$tips_mode, "run_tips_pipeline")) {
    man <- fold_mans[[paste(r, f, sep = "_")]]
    if (is.null(man)) man <- build_tips_module_for_fold(f, r, test_clones, cfg, work, read_clone_matrix(cfg$public_clone))
    fold_status <- if (is.null(man$status)) "success" else man$status
    fold_error <- if (is.null(man$error)) "" else man$error
    tips_genes <- man$module$CF$genes_lean
    if (!length(tips_genes)) tips_genes <- man$module$CF$genes
    tips_baso_genes <- man$module$CM$genes_lean
    if (!length(tips_baso_genes)) tips_baso_genes <- man$module$CM$genes
    if (length(man$graph_notsimp$CF)) ppi_meg <- man$graph_notsimp$CF
    if (length(man$graph_notsimp$CM)) ppi_baso <- man$graph_notsimp$CM
  } else if (identical(cfg$tips_mode, "frozen_genes")) {
    tips_genes <- sets$tips_meg_frozen
    tips_baso_genes <- sets$tips_baso_frozen
  } else {
    message("[04] HELDOUT_TIPS_MODE=", cfg$tips_mode, " — Pearson-delta proxy, not code_core_11_10vs17")
    tips_genes <- reweight_module(cf_edges, train_cells, positive_delta = TRUE)
    reweight_baso <- function(edges) {
      cp <- intersect(train_cells, names(state_of)[state_of == as.character(cfg$state_id)])
      cm <- intersect(train_cells, names(state_of)[state_of == as.character(cfg$cm_cluster)])
      if (!nrow(edges)) return(character())
      w1 <- pair_pearson(expr, edges$from_r, edges$to_r, cp)
      w2 <- pair_pearson(expr, edges$from_r, edges$to_r, cm)
      delta <- w2 - w1
      keep <- which(is.finite(delta) & delta > 0)
      if (!length(keep)) keep <- seq_len(nrow(edges))
      unique(c(edges$from[keep], edges$to[keep]))
    }
    tips_baso_genes <- reweight_baso(cm_edges)
  }
  n_keep <- max(length(tips_genes), 5L)
  meg_markers <- wilcox_up_genes(train_meg, train_not_meg, hvg_use, n_keep)
  node_de_megup <- wilcox_up_genes(train_meg, train_c11, hvg_use, n_keep)
  rand_sets <- matched_random(tips_genes, hvg_use, train_c11, cfg$n_random, cfg$seed + 1000L * r + f)

  fold_modules[[paste(r, f, sep = "_")]] <- list(
    repeat_id = r, fold = f, status = fold_status, error = fold_error,
    tips_meg = norm_genes(tips_genes), tips_baso = norm_genes(tips_baso_genes),
    cts_c11 = sets$cts_c11, mutrans_td_meg = sets$mutrans_td_meg,
    ppi_unweighted_meg = norm_genes(ppi_meg),
    tips_meg_frozen = sets$tips_meg_frozen,
    metacell_tips_meg = sets$metacell_tips_meg,
    meg_markers_train = norm_genes(meg_markers),
    node_de_megup = norm_genes(node_de_megup)
  )

  score_one <- function(genes, method) {
    empty <- !length(match_genes_to_rownames(genes, rownames(expr)))
    out <- data.frame(
      repeat_id = r, fold = f, method = method,
      clone_id = as.integer(test_clones),
      score = NA_real_, n_genes = 0L, module_empty = empty,
      stringsAsFactors = FALSE
    )
    if (empty) return(out)
    sc <- fit_z_then_score(expr, genes, train_c11, test_cell, cfg$score_fun)
    cl <- score_clone_mean(sc$score, test_cell, d2_cells$clone_id[match(test_cell, d2_cells$cell_id)])
    out$n_genes <- length(sc$genes_used)
    m <- match(as.integer(names(cl)), out$clone_id)
    ok <- !is.na(m)
    out$score[m[ok]] <- unname(cl[ok])
    out
  }

  tips_now <- score_one(tips_genes, "tips_meg")
  baso_now <- score_one(tips_baso_genes, "tips_baso")
  pred_rows[[length(pred_rows) + 1L]] <- tips_now
  pred_rows[[length(pred_rows) + 1L]] <- score_one(sets$cts_c11, "cts_c11")
  pred_rows[[length(pred_rows) + 1L]] <- score_one(sets$mutrans_td_meg, "mutrans_td_meg")
  pred_rows[[length(pred_rows) + 1L]] <- score_one(ppi_meg, "ppi_unweighted_meg")
  pred_rows[[length(pred_rows) + 1L]] <- score_one(sets$tips_meg_frozen, "tips_meg_frozen")
  if (length(sets$metacell_tips_meg)) {
    pred_rows[[length(pred_rows) + 1L]] <- score_one(sets$metacell_tips_meg, "metacell_tips_meg")
  }
  pred_rows[[length(pred_rows) + 1L]] <- score_one(meg_markers, "meg_markers_train")
  pred_rows[[length(pred_rows) + 1L]] <- score_one(node_de_megup, "node_de_megup")
  pred_rows[[length(pred_rows) + 1L]] <- baso_now
  pred_rows[[length(pred_rows) + 1L]] <- score_one(sets$mutrans_td_baso, "mutrans_td_baso")
  pred_rows[[length(pred_rows) + 1L]] <- score_one(ppi_baso, "ppi_unweighted_baso")

  success_rows[[length(success_rows) + 1L]] <- data.frame(
    fold = f, rep = r, status = fold_status, error = fold_error,
    meg_n_genes = length(norm_genes(tips_genes)),
    baso_n_genes = length(norm_genes(tips_baso_genes)),
    meg_recovered = as.integer(length(norm_genes(tips_genes)) > 0L && !identical(fold_status, "error")),
    baso_recovered = as.integer(length(norm_genes(tips_baso_genes)) > 0L && !identical(fold_status, "error")),
    n_test_clones = length(test_clones),
    n_test_clones_scored_meg = sum(is.finite(tips_now$score)),
    n_test_clones_scored_baso = sum(is.finite(baso_now$score)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(rand_sets)) {
    sc <- fit_z_then_score(expr, rand_sets[[i]], train_c11, test_cell, cfg$score_fun)
    cl <- score_clone_mean(sc$score, test_cell, d2_cells$clone_id[match(test_cell, d2_cells$cell_id)])
    rr <- data.frame(
      repeat_id = r, fold = f, draw = i, method = "random_matched",
      clone_id = as.integer(test_clones), score = NA_real_,
      n_genes = length(sc$genes_used), stringsAsFactors = FALSE
    )
    m <- match(as.integer(names(cl)), rr$clone_id)
    ok <- !is.na(m)
    rr$score[m[ok]] <- unname(cl[ok])
    rand_rows[[length(rand_rows) + 1L]] <- rr
  }
}

} ## end repeats

pred_long <- do.call(rbind, pred_rows)
covar_cols <- intersect(c(
  "clone_id", "n_d2_c11", "n_d6", "n_meg_d6", "n_baso_d6",
  "n_undiff_d6", "n_mature_d6", "frac_meg_d6", "frac_baso_d6",
  "frac_meg_mature_d6", "frac_baso_mature_d6",
  "library", "starting_population",
  "n_d6_well1", "n_meg_d6_well1", "n_baso_d6_well1", "frac_meg_d6_well1",
  "n_d6_well2", "n_meg_d6_well2", "n_baso_d6_well2", "frac_meg_d6_well2",
  "in_both_later_wells", "meg_positive", "baso_positive"
), names(primary))
pred_long <- merge(pred_long, primary[, covar_cols], by = "clone_id", all.x = TRUE)
write_tsv(pred_long, file.path(cfg$results_dir, "tables", "04_oof_clone_predictions_long.tsv"))

pred <- stats::aggregate(
  cbind(score, n_genes) ~ clone_id + method, pred_long,
  function(x) mean(x, na.rm = TRUE)
)
pred$score[!is.finite(pred$score)] <- NA_real_
empty_flag <- stats::aggregate(module_empty ~ clone_id + method, pred_long, function(x) any(x))
names(empty_flag)[3] <- "module_empty"
pred <- merge(pred, empty_flag, by = c("clone_id", "method"), all.x = TRUE)
pred <- merge(pred, primary[, covar_cols], by = "clone_id", all.x = TRUE)
pred <- merge(pred, primary[, c("clone_id", "fold")], by = "clone_id", all.x = TRUE)
write_tsv(pred, file.path(cfg$results_dir, "tables", "04_oof_clone_predictions.tsv"))

succ <- do.call(rbind, success_rows)
write_tsv(succ, file.path(cfg$results_dir, "tables", "04_fold_module_success.tsv"))
n_eligible <- nrow(primary)
tips_pred <- pred[pred$method == "tips_meg", ]
n_scored <- sum(is.finite(tips_pred$score))
msg_cov <- sprintf(
  "A non-empty Meg module was recovered in %d/%d training folds, providing out-of-fold scores for %d/%d eligible clones.",
  sum(succ$meg_recovered), nrow(succ), n_scored, n_eligible
)
writeLines(msg_cov, file.path(cfg$results_dir, "04_coverage.txt"))
message("[04] ", msg_cov)
msg_baso <- sprintf(
  "A non-empty Baso module was recovered in %d/%d training folds, providing out-of-fold scores for %d/%d eligible clones.",
  sum(succ$baso_recovered), nrow(succ),
  sum(is.finite(pred$score[pred$method == "tips_baso"])), n_eligible
)
message("[04] ", msg_baso)

rand_long <- if (length(rand_rows)) do.call(rbind, rand_rows) else data.frame()
if (nrow(rand_long)) {
  rand_long <- merge(rand_long, primary[, covar_cols], by = "clone_id", all.x = TRUE)
  write_tsv(rand_long, file.path(cfg$results_dir, "tables", "04_random_set_scores.tsv"))
}

saveRDS(
  list(pred = pred, pred_long = pred_long, rand_long = rand_long,
       fold_modules = fold_modules, success = succ,
       tips_mode = cfg$tips_mode, n_repeats = n_repeats,
       n_eligible = n_eligible, coverage = msg_cov),
  file.path(cfg$results_dir, "rds", "04_oof_predictions.rds")
)

mod_sizes <- do.call(rbind, lapply(names(fold_modules), function(k) {
  m <- fold_modules[[k]]
  chr <- names(m)[vapply(m, is.character, logical(1))]
  data.frame(
    repeat_fold = k,
    fold = m$fold,
    repeat_id = m$repeat_id,
    set = chr,
    n = vapply(m[chr], length, integer(1)),
    genes = vapply(m[chr], function(x) paste(sort(x), collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
}))
write_tsv(mod_sizes, file.path(cfg$results_dir, "tables", "04_fold_module_genes.tsv"))
message("[04] pooled OOF predictions: ", nrow(pred), " clone-method rows across ",
        n_repeats, " repeats; methods: ", paste(unique(pred$method), collapse = ", "))
message("[04] test-clone day-6 labels were not used to choose genes.")
