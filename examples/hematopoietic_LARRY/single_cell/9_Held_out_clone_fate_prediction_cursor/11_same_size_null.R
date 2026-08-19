## 11_same_size_null.R
##
## Does TIPS pick an unusually fate-informative compact network?
##
## The paired TIPS−CTS Δρ test in step 10 compares a 5-gene module to the
## full 71-gene CTS. That confounds size with selection. This step holds size
## (and, for PPI, degree / connectedness) fixed:
##   1. cts_ksubset — random k-subsets of C11 CTS
##   2. ppi_degree  — degree-matched k nodes from the C11 STRING graph
##   3. ppi_connected — connected k-node subgraphs, then pick the candidate
##      whose sorted log-degree vector is closest to TIPS (the "subnetwork" null)
##   4. mutrans_ksubset — random k-subsets of MuTrans Meg TD genes
##   5. meg_marker_ksubset — random k-subsets of training Meg-up HVGs (top pool)
##   6. node_de_ksubset — random k-subsets of training C11-vs-Meg up HVGs (top pool)
##   7. weinreb_ksubset — random k-subsets of Weinreb Table S3 Meg progenitor DGE
##
## Training Meg markers and node DE in step 04 are already truncated to k genes.
## Those specific top-k lists are the step-10 paired comparators. Here the null
## is random compact draws from a larger training-only DE pool (default top 50).
## Weinreb S3 is still not a held-out predictor; it is only a published Meg
## gene universe for the same-size compactness null.
##
## Each draw is scored with the same OOF rule as TIPS (train-C11 z, test day-2
## cells, clone mean). Gene lists are not refit. TIPS genes themselves are not
## a CTS subset (GATA2 is outside CTS); the CTS null is a competing compact
## CTS module, not "TIPS as one CTS subset among many."
##
##   Sys.setenv(HELDOUT_RUN = "11")
##   source(".../run_heldout_pipeline.R")

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = here::here("examples", "hematopoietic_LARRY", "single_cell", "9_Held_out_clone_fate_prediction_cursor")
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_ensure()

n_draw <- as.integer(Sys.getenv("HELDOUT_N_NULL", "500"))
if (!is.finite(n_draw) || n_draw < 50L) n_draw <- 500L
n_cand <- as.integer(Sys.getenv("HELDOUT_N_NULL_CAND", "30"))
if (!is.finite(n_cand) || n_cand < 5L) n_cand <- 30L
n_pool <- as.integer(Sys.getenv("HELDOUT_N_NULL_POOL", "50"))
if (!is.finite(n_pool) || n_pool < 10L) n_pool <- 50L
seed_null <- cfg$seed + 1111L

need_files <- c(
  lock = file.path(cfg$results_dir, "rds", "02_locked_clones_folds.rds"),
  sets = file.path(cfg$results_dir, "rds", "03_gene_sets.rds"),
  ds = file.path(cfg$results_dir, "rds", "01_clone_dataset.rds"),
  pred = file.path(cfg$results_dir, "rds", "04_oof_predictions.rds")
)
miss <- names(need_files)[!file.exists(need_files)]
if (length(miss)) stop("Run earlier steps first; missing: ", paste(need_files[miss], collapse = ", "))
if (!file.exists(cfg$seurat_rds)) stop("Missing Seurat: ", cfg$seurat_rds)
if (!requireNamespace("igraph", quietly = TRUE)) stop("11_same_size_null.R requires igraph")

lock <- readRDS(need_files[["lock"]])
sets <- readRDS(need_files[["sets"]])
ds <- readRDS(need_files[["ds"]])
obj04 <- readRDS(need_files[["pred"]])
primary <- lock$primary
fold_modules <- obj04$fold_modules
if (is.null(fold_modules) || !length(fold_modules)) stop("04 fold_modules missing")
work <- ds$work_meta
d2_cells <- ds$d2_c11_cells
d2_cells <- d2_cells[d2_cells$clone_id %in% primary$clone_id, ]

fold_graph_path <- function(fold_id, repeat_id = 1L) {
  file.path(
    cfg$fold_modules_dir,
    paste0("fold", fold_id, "_rep", repeat_id),
    paste0("results_core_", cfg$tips_arm_tag),
    "GSE140802_STRING_graph_perState_notsimplified.rds"
  )
}

load_c11_graph <- function(path) {
  if (!file.exists(path)) path <- cfg$cell_tips_graph
  gl <- load_tips_state_graphs(path)
  g <- graph_for_cluster(gl, cfg$state_id)
  if (is.null(g)) return(NULL)
  igraph::V(g)$name <- toupper(as.character(igraph::V(g)$name))
  g
}

deg_vec <- function(g, genes) {
  d <- igraph::degree(g)
  names(d) <- toupper(as.character(igraph::V(g)$name))
  unname(d[toupper(as.character(genes))])
}

deg_l1 <- function(d_tips, d_null) {
  a <- sort(log1p(pmax(as.numeric(d_tips), 0)))
  b <- sort(log1p(pmax(as.numeric(d_null), 0)))
  k <- min(length(a), length(b))
  if (!k || length(a) != length(b)) return(Inf)
  sum(abs(a[seq_len(k)] - b[seq_len(k)]))
}

sample_cts_k <- function(pool, k, n, seed) {
  set.seed(seed)
  pool <- unique(pool)
  if (length(pool) < k) {
    return(replicate(n, pool, simplify = FALSE))
  }
  replicate(n, sample(pool, k), simplify = FALSE)
}

sample_degree_matched <- function(pool, pool_deg, target_deg, n, seed, k = NULL) {
  set.seed(seed)
  pool <- unique(toupper(as.character(pool)))
  names(pool_deg) <- toupper(as.character(names(pool_deg)))
  pool_deg <- pool_deg[pool]
  pool_deg <- pool_deg[is.finite(pool_deg)]
  pool <- names(pool_deg)
  tgt <- as.numeric(target_deg)
  if (is.null(k)) k <- length(tgt)
  tgt <- tgt[is.finite(tgt)]
  if (!length(tgt) && length(pool_deg)) tgt <- rep(stats::median(pool_deg), k)
  if (length(tgt) < k && length(tgt)) {
    tgt <- c(tgt, rep(stats::median(tgt), k - length(tgt)))
  }
  k <- length(tgt)
  if (!k || length(pool) < k) {
    return(replicate(n, head(pool, k), simplify = FALSE))
  }
        lapply(seq_len(n), function(i) {
            picked <- character()
            left <- pool
            left_d <- pool_deg
            for (td in sample(tgt)) {
                if (!length(left)) break
                ad <- abs(left_d - td)
                band <- max(1, 0.25 * td)
                in_band <- which(ad <= band)
                if (!length(in_band)) in_band <- which.min(ad)
                w <- sample(in_band, 1L)
                picked <- c(picked, left[w])
                left <- left[-w]
                left_d <- left_d[-w]
            }
            picked
        })
}

grow_connected <- function(g, k, n_try = 40L) {
  vn <- toupper(as.character(igraph::V(g)$name))
  igraph::V(g)$name <- vn
  d <- igraph::degree(g)
  names(d) <- vn
  seeds <- vn[d > 0]
  if (!length(seeds)) seeds <- vn
  if (length(vn) < k) return(vn)
  for (t in seq_len(n_try)) {
    s <- sample(seeds, 1L)
    while (length(s) < k) {
      nbr <- unique(unlist(lapply(s, function(v) {
        as.character(igraph::V(g)$name[igraph::neighbors(g, v)])
      })))
      nbr <- setdiff(nbr, s)
      if (!length(nbr)) break
      s <- c(s, sample(nbr, 1L))
    }
    if (length(s) >= k) return(s[seq_len(k)])
  }
  NULL
}

sample_connected_degree <- function(g, target_deg, n, n_cand, seed) {
  set.seed(seed)
  k <- length(target_deg)
  lapply(seq_len(n), function(i) {
    best <- NULL
    best_d <- Inf
    for (j in seq_len(n_cand)) {
      cand <- grow_connected(g, k)
      if (is.null(cand) || length(cand) < k) next
      dd <- deg_l1(target_deg, deg_vec(g, cand))
      if (dd < best_d) {
        best <- cand
        best_d <- dd
      }
    }
    if (is.null(best)) {
      vn <- toupper(as.character(igraph::V(g)$name))
      return(sample(vn, min(k, length(vn))))
    }
    best
  })
}

induced_connected <- function(g, genes) {
  genes <- intersect(toupper(as.character(genes)), toupper(as.character(igraph::V(g)$name)))
  if (length(genes) < 2L) return(NA)
  sg <- igraph::induced_subgraph(g, genes)
  igraph::is_connected(sg)
}

tips_obs_rho <- NA_real_
tips_row <- obj04$pred
if (!is.null(tips_row) && nrow(tips_row)) {
  tt <- tips_row[tips_row$method == "tips_meg", ]
  if (nrow(tt)) tips_obs_rho <- spearman_safe(tt$score, tt$frac_meg_d6)$rho
}

## --- expression ---
tips_all <- unique(unlist(lapply(fold_modules, function(m) m$tips_meg)))
cts_all <- norm_genes(sets$cts_c11)
mutrans_all <- norm_genes(sets$mutrans_td_meg)
weinreb_all <- norm_genes(sets$weinreb_meg)
hvg_all <- norm_genes(sets$hvg)
ppi_nodes_frozen <- character()
g_frozen <- load_c11_graph(cfg$cell_tips_graph)
if (!is.null(g_frozen)) ppi_nodes_frozen <- toupper(as.character(igraph::V(g_frozen)$name))
candidate <- unique(c(tips_all, cts_all, mutrans_all, weinreb_all, hvg_all, ppi_nodes_frozen))
message("[11] loading Seurat for ", length(candidate), " candidate genes")
suppressPackageStartupMessages(library(Seurat))
seu <- readRDS(cfg$seurat_rds)
expr <- extract_seurat_expr(seu, candidate)
meta <- seu@meta.data
meta$cell_id <- rownames(meta)
rm(seu)
gc()
message("[11] expression ", nrow(expr), " genes x ", ncol(expr), " cells")

work_map <- setNames(work$clone_id, work$cell_id)
clone_of <- work_map[colnames(expr)]
state_of <- as.character(meta[[cfg$state_col]][match(colnames(expr), rownames(meta))])
names(state_of) <- colnames(expr)
type_of <- as.character(meta$Cell.type.annotation[match(colnames(expr), rownames(meta))])
names(type_of) <- colnames(expr)
time_of <- as_int_like(meta$Time.point[match(colnames(expr), rownames(meta))])
names(time_of) <- colnames(expr)

feat_map <- setNames(as.character(rownames(expr)), toupper(as.character(rownames(expr))))
in_expr <- function(genes) {
  hit <- unname(feat_map[toupper(as.character(genes))])
  unique(hit[!is.na(hit)])
}

cts_expr <- in_expr(cts_all)
mutrans_expr <- in_expr(mutrans_all)
weinreb_expr <- in_expr(weinreb_all)
hvg_expr <- in_expr(hvg_all)
message("[11] CTS in expression: ", length(cts_expr), "/", length(cts_all),
        "  MuTrans TD: ", length(mutrans_expr), "/", length(mutrans_all),
        "  Weinreb S3 Meg: ", length(weinreb_expr), "/", length(weinreb_all),
        "  HVG: ", length(hvg_expr), "/", length(hvg_all))

## Training-only upregulated pool (mean difference). Larger than the step-04
## top-k module so we can draw same-size random compact sets.
up_pool <- function(pos_cells, neg_cells, universe, n_keep) {
  pos_cells <- intersect(pos_cells, colnames(expr))
  neg_cells <- intersect(neg_cells, colnames(expr))
  universe <- intersect(universe, rownames(expr))
  if (length(pos_cells) < 10L || length(neg_cells) < 10L || length(universe) < 5L) {
    return(character())
  }
  m <- as.matrix(expr[universe, c(pos_cells, neg_cells), drop = FALSE])
  npos <- length(pos_cells)
  lfc <- rowMeans(m[, seq_len(npos), drop = FALSE], na.rm = TRUE) -
    rowMeans(m[, (npos + 1):ncol(m), drop = FALSE], na.rm = TRUE)
  keep <- names(sort(lfc[is.finite(lfc) & lfc > 0], decreasing = TRUE))
  head(keep, n_keep)
}

precompute_clone_z <- function(genes, train_cells, test_cells, test_clones_by_cell) {
  genes <- intersect(genes, rownames(expr))
  train_cells <- intersect(train_cells, colnames(expr))
  test_cells <- intersect(test_cells, colnames(expr))
  clones <- unique(as.integer(test_clones_by_cell[test_cells]))
  clones <- clones[!is.na(clones)]
  out <- matrix(NA_real_, nrow = length(genes), ncol = length(clones),
                dimnames = list(genes, as.character(clones)))
  if (!length(genes) || !length(train_cells) || !length(test_cells) || !length(clones)) {
    return(out)
  }
  m_train <- as.matrix(expr[genes, train_cells, drop = FALSE])
  mu <- rowMeans(m_train, na.rm = TRUE)
  sdv <- apply(m_train, 1, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  z <- sweep(as.matrix(expr[genes, test_cells, drop = FALSE]), 1, mu, "-")
  z <- sweep(z, 1, sdv, "/")
  z[!is.finite(z)] <- 0
  for (cid in clones) {
    cells <- names(test_clones_by_cell)[test_clones_by_cell == cid]
    cells <- intersect(cells, test_cells)
    if (!length(cells)) next
    out[, as.character(cid)] <- rowMeans(z[, cells, drop = FALSE], na.rm = TRUE)
  }
  out
}

score_sets <- function(clone_z, gene_sets) {
  cn <- colnames(clone_z)
  mat <- matrix(NA_real_, nrow = length(gene_sets), ncol = length(cn),
                dimnames = list(NULL, cn))
  ng <- integer(length(gene_sets))
  for (i in seq_along(gene_sets)) {
    use <- intersect(in_expr(gene_sets[[i]]), rownames(clone_z))
    ng[i] <- length(use)
    if (!length(use)) next
    mat[i, ] <- colMeans(clone_z[use, , drop = FALSE], na.rm = TRUE)
  }
  list(score = mat, n_genes = ng)
}

families <- c(
  "cts_ksubset", "ppi_degree", "ppi_connected",
  "mutrans_ksubset", "meg_marker_ksubset", "node_de_ksubset",
  "weinreb_ksubset"
)
fold_scores <- list(
  cts_ksubset = list(),
  ppi_degree = list(),
  ppi_connected = list(),
  mutrans_ksubset = list(),
  meg_marker_ksubset = list(),
  node_de_ksubset = list(),
  weinreb_ksubset = list(),
  tips = list()
)
audit_rows <- list()
set.seed(seed_null)

fold_keys <- names(fold_modules)
for (fk in fold_keys) {
  m <- fold_modules[[fk]]
  f <- as.integer(m$fold)
  r <- as.integer(m$repeat_id)
  tips_genes <- in_expr(m$tips_meg)
  k <- length(tips_genes)
  if (k < 2L) {
    message("[11] skip fold ", f, ": empty TIPS module")
    next
  }
  test_clones <- primary$clone_id[primary$fold == f]
  test_cell <- intersect(d2_cells$cell_id[d2_cells$clone_id %in% test_clones], colnames(expr))
  is_test <- !is.na(clone_of) & clone_of %in% test_clones
  train_cells <- colnames(expr)[!is_test]
  train_c11 <- intersect(train_cells, names(state_of)[state_of == as.character(cfg$state_id)])
  train_d6 <- intersect(train_cells, names(time_of)[time_of == cfg$day_late])
  train_meg <- intersect(train_d6, names(type_of)[type_of == cfg$meg_label])
  train_not_meg <- setdiff(train_d6, train_meg)
  test_clone_by_cell <- setNames(as.integer(clone_of[test_cell]), test_cell)

  gpath <- fold_graph_path(f, r)
  g <- load_c11_graph(gpath)
  if (is.null(g)) g <- g_frozen
  ppi_pool <- character()
  tips_deg <- rep(NA_real_, k)
  names(tips_deg) <- tips_genes
  connected_tips <- NA
  if (!is.null(g)) {
    ppi_pool <- in_expr(igraph::V(g)$name)
    gd <- igraph::degree(g)
    names(gd) <- toupper(as.character(igraph::V(g)$name))
    tips_in_g <- intersect(toupper(tips_genes), names(gd))
    tips_deg[match(tips_in_g, toupper(tips_genes))] <- unname(gd[tips_in_g])
    connected_tips <- induced_connected(g, tips_genes)
  }
  meg_pool <- up_pool(train_meg, train_not_meg, hvg_expr, n_pool)
  node_pool <- up_pool(train_meg, train_c11, hvg_expr, n_pool)
  overlap_cts <- length(intersect(toupper(tips_genes), toupper(cts_expr)))
  overlap_td <- length(intersect(toupper(tips_genes), toupper(mutrans_expr)))
  overlap_s3 <- length(intersect(toupper(tips_genes), toupper(weinreb_expr)))

  pool_genes <- unique(c(cts_expr, ppi_pool, tips_genes, mutrans_expr, meg_pool, node_pool, weinreb_expr))
  clone_z <- precompute_clone_z(pool_genes, train_c11, test_cell, test_clone_by_cell)
  message("[11] fold ", f, " k=", k, " test_clones=", length(test_clones),
          " CTS=", length(cts_expr), " PPI=", length(ppi_pool),
          " MuTrans=", length(mutrans_expr), " MegPool=", length(meg_pool),
          " NodeDE=", length(node_pool), " WeinrebS3=", length(weinreb_expr),
          " TIPS∩CTS=", overlap_cts, " TIPS∩TD=", overlap_td,
          " TIPS∩S3=", overlap_s3, " TIPS_connected=", connected_tips)

  cts_sets <- sample_cts_k(cts_expr, k, n_draw, seed_null + 10L * f)
  td_sets <- sample_cts_k(mutrans_expr, k, n_draw, seed_null + 40L * f)
  meg_sets <- sample_cts_k(meg_pool, k, n_draw, seed_null + 50L * f)
  node_sets <- sample_cts_k(node_pool, k, n_draw, seed_null + 60L * f)
  s3_sets <- sample_cts_k(weinreb_expr, k, n_draw, seed_null + 70L * f)
  if (!is.null(g) && length(ppi_pool) >= k && sum(is.finite(tips_deg)) >= 2L) {
    pool_deg <- gd[toupper(ppi_pool)]
    names(pool_deg) <- ppi_pool
    deg_sets <- sample_degree_matched(ppi_pool, pool_deg, tips_deg, n_draw, seed_null + 20L * f, k = k)
    conn_sets <- sample_connected_degree(g, ifelse(is.finite(tips_deg), tips_deg, stats::median(tips_deg, na.rm = TRUE)), n_draw, n_cand, seed_null + 30L * f)
  } else {
    deg_sets <- replicate(n_draw, character(), simplify = FALSE)
    conn_sets <- replicate(n_draw, character(), simplify = FALSE)
    message("[11] fold ", f, ": PPI graph too small; PPI nulls will be empty")
  }

  sc_tips <- score_sets(clone_z, list(tips_genes))
  sc_cts <- score_sets(clone_z, cts_sets)
  sc_deg <- score_sets(clone_z, deg_sets)
  sc_conn <- score_sets(clone_z, conn_sets)
  sc_td <- score_sets(clone_z, td_sets)
  sc_meg <- score_sets(clone_z, meg_sets)
  sc_node <- score_sets(clone_z, node_sets)
  sc_s3 <- score_sets(clone_z, s3_sets)

  fold_scores$tips[[fk]] <- sc_tips$score[1, ]
  fold_scores$cts_ksubset[[fk]] <- sc_cts$score
  fold_scores$ppi_degree[[fk]] <- sc_deg$score
  fold_scores$ppi_connected[[fk]] <- sc_conn$score
  fold_scores$mutrans_ksubset[[fk]] <- sc_td$score
  fold_scores$meg_marker_ksubset[[fk]] <- sc_meg$score
  fold_scores$node_de_ksubset[[fk]] <- sc_node$score
  fold_scores$weinreb_ksubset[[fk]] <- sc_s3$score

  audit_rows[[length(audit_rows) + 1L]] <- data.frame(
    fold = f, repeat_id = r, k = k,
    tips_genes = paste(sort(tips_genes), collapse = ","),
    n_cts_pool = length(cts_expr),
    n_ppi_pool = length(ppi_pool),
    n_mutrans_pool = length(mutrans_expr),
    n_meg_marker_pool = length(meg_pool),
    n_node_de_pool = length(node_pool),
    n_weinreb_pool = length(weinreb_expr),
    n_tips_in_cts = overlap_cts,
    n_tips_in_mutrans = overlap_td,
    n_tips_in_weinreb = overlap_s3,
    n_tips_in_ppi = sum(toupper(tips_genes) %in% toupper(ppi_pool)),
    tips_connected = connected_tips,
    tips_degrees = paste(sprintf("%s=%s", names(tips_deg), ifelse(is.finite(tips_deg), as.character(tips_deg), "NA")), collapse = ";"),
    mean_cts_n_genes = mean(sc_cts$n_genes),
    mean_ppi_degree_n_genes = mean(sc_deg$n_genes),
    mean_ppi_connected_n_genes = mean(sc_conn$n_genes),
    mean_mutrans_n_genes = mean(sc_td$n_genes),
    mean_meg_marker_n_genes = mean(sc_meg$n_genes),
    mean_node_de_n_genes = mean(sc_node$n_genes),
    mean_weinreb_n_genes = mean(sc_s3$n_genes),
    stringsAsFactors = FALSE
  )
}

pool_draws <- function(fold_mat_list) {
  ## each element is n_draw x clones OR a named clone vector for TIPS
  if (!length(fold_mat_list)) return(NULL)
  first <- fold_mat_list[[1]]
  if (is.null(dim(first))) {
    sc <- do.call(c, unname(fold_mat_list))
    nm <- names(sc)
    cid <- suppressWarnings(as.integer(nm))
    names(sc) <- as.character(cid)
    return(list(score = sc, clone_id = cid))
  }
  clones <- unique(unlist(lapply(fold_mat_list, colnames), use.names = FALSE))
  out <- matrix(NA_real_, nrow = n_draw, ncol = length(clones),
                dimnames = list(NULL, clones))
  for (mat in fold_mat_list) {
    if (is.null(mat) || !nrow(mat)) next
    cn <- intersect(colnames(mat), clones)
    out[, cn] <- mat[, cn, drop = FALSE]
  }
  list(score = out, clone_id = as.integer(clones))
}

frac <- setNames(primary$frac_meg_d6, as.character(primary$clone_id))
n_meg <- setNames(primary$n_meg_d6, as.character(primary$clone_id))
n_tot <- setNames(primary$n_d6, as.character(primary$clone_id))
lab <- setNames(as.integer(primary$meg_positive %in% c(TRUE, "TRUE", 1, "1")), as.character(primary$clone_id))

rho_of_scores <- function(sc, clone_id) {
  y <- frac[as.character(clone_id)]
  if (is.null(dim(sc))) {
    return(spearman_safe(as.numeric(sc), y)$rho)
  }
  apply(sc, 1, function(v) spearman_safe(v, y)$rho)
}

auc_of_scores <- function(sc, clone_id) {
  y <- lab[as.character(clone_id)]
  if (is.null(dim(sc))) return(roc_auc(y, as.numeric(sc)))
  apply(sc, 1, function(v) roc_auc(y, v))
}

mc_p <- function(null_stat, obs) {
  x <- null_stat[is.finite(null_stat)]
  if (!length(x) || !is.finite(obs)) return(NA_real_)
  (1 + sum(x >= obs)) / (1 + length(x))
}

tips_pool <- pool_draws(fold_scores$tips)
tips_rho_check <- rho_of_scores(tips_pool$score, tips_pool$clone_id)
message("[11] TIPS Spearman from this scoring=", signif(tips_rho_check, 3),
        " (step 04 reported ", signif(tips_obs_rho, 3), ")")

sum_family <- function(name, fold_list) {
  pooled <- pool_draws(fold_list)
  if (is.null(pooled) || is.null(dim(pooled$score))) {
    return(list(tab = NULL, rho = NULL, auc = NULL))
  }
  rr <- rho_of_scores(pooled$score, pooled$clone_id)
  aa <- auc_of_scores(pooled$score, pooled$clone_id)
  tab <- data.frame(
    family = name,
    n_draw = length(rr),
    n_clones = length(pooled$clone_id),
    k_tips = length(in_expr(fold_modules[[1]]$tips_meg)),
    tips_rho = unname(tips_rho_check),
    tips_auroc = auc_of_scores(tips_pool$score, tips_pool$clone_id),
    null_rho_median = stats::median(rr, na.rm = TRUE),
    null_rho_q025 = unname(stats::quantile(rr, 0.025, na.rm = TRUE, names = FALSE)),
    null_rho_q975 = unname(stats::quantile(rr, 0.975, na.rm = TRUE, names = FALSE)),
    p_tips_ge_null_rho = mc_p(rr, tips_rho_check),
    null_auroc_median = stats::median(aa, na.rm = TRUE),
    null_auroc_q025 = unname(stats::quantile(aa, 0.025, na.rm = TRUE, names = FALSE)),
    null_auroc_q975 = unname(stats::quantile(aa, 0.975, na.rm = TRUE, names = FALSE)),
    p_tips_ge_null_auroc = mc_p(aa, auc_of_scores(tips_pool$score, tips_pool$clone_id)),
    n_finite_rho = sum(is.finite(rr)),
    stringsAsFactors = FALSE
  )
  list(tab = tab, rho = rr, auc = aa, pooled = pooled)
}

res <- lapply(setNames(families, families), function(fam) sum_family(fam, fold_scores[[fam]]))
sum_tab <- do.call(rbind, lapply(res, `[[`, "tab"))
rownames(sum_tab) <- NULL
audit <- do.call(rbind, audit_rows)

draw_long <- do.call(rbind, lapply(families, function(fam) {
  rr <- res[[fam]]$rho
  if (is.null(rr)) return(NULL)
  data.frame(family = fam, draw = seq_along(rr), spearman_rho = rr,
             auroc = res[[fam]]$auc, stringsAsFactors = FALSE)
}))

write_tsv(sum_tab, file.path(cfg$results_dir, "tables", "11_same_size_null_summary.tsv"))
write_tsv(audit, file.path(cfg$results_dir, "tables", "11_same_size_null_fold_audit.tsv"))
write_tsv(draw_long, file.path(cfg$results_dir, "tables", "11_same_size_null_draws.tsv"))
saveRDS(
  list(
    n_draw = n_draw, seed = seed_null, summary = sum_tab, audit = audit,
    draws = draw_long, tips_rho = tips_rho_check, tips_rho_04 = tips_obs_rho
  ),
  file.path(cfg$results_dir, "rds", "11_same_size_null.rds")
)

fam_lab <- c(
  cts_ksubset = "Same-size CTS subsets",
  ppi_degree = "Degree-matched PPI nodes",
  ppi_connected = "Degree-matched connected PPI",
  mutrans_ksubset = "Same-size MuTrans TD subsets",
  meg_marker_ksubset = "Same-size training Meg-marker subsets",
  node_de_ksubset = "Same-size node-DE subsets",
  weinreb_ksubset = "Same-size Weinreb S3 Meg subsets"
)
message("[11] TIPS is not a CTS subset (GATA2 is outside CTS). CTS/MuTrans/S3 nulls = competing k-gene modules from those lists.")
message("[11] Meg-marker and node-DE nulls draw k genes from a training-only top-", n_pool, " upregulated HVG pool.")
message("[11] Weinreb S3 null draws k genes from published Meg progenitor DGE; S3 is not used as a held-out predictor.")
for (i in seq_len(nrow(sum_tab))) {
  message("[11] ", fam_lab[[sum_tab$family[i]]], ": TIPS ρ=", signif(sum_tab$tips_rho[i], 3),
          " vs null median ", signif(sum_tab$null_rho_median[i], 3),
          " [", signif(sum_tab$null_rho_q025[i], 3), ", ", signif(sum_tab$null_rho_q975[i], 3), "]",
          "  P(TIPS ≥ null)=", signif(sum_tab$p_tips_ge_null_rho[i], 3))
}

## ---- figures ----
figdir <- file.path(cfg$results_dir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("[11] ggplot2 not installed; tables written, figures skipped")
} else {
  suppressPackageStartupMessages(library(ggplot2))
  col_tips <- "#C00000"
  col_anno <- "grey45"
  theme_held <- function() {
    theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10.5, face = "bold"),
        plot.subtitle = element_text(size = 8, lineheight = 1.15),
        plot.caption = element_text(size = 7.5, colour = col_anno, hjust = 0)
      )
  }
  cairo_ok <- isTRUE(capabilities("cairo"))
  save_pdf <- function(file, plot, width, height) {
    tmp <- paste0(file, ".tmp.pdf")
    if (file.exists(tmp)) unlink(tmp)
    wrote <- FALSE
    if (cairo_ok) {
      wrote <- isTRUE(tryCatch({
        ggsave(tmp, plot, width = width, height = height, bg = "white",
               device = grDevices::cairo_pdf)
        TRUE
      }, error = function(e) FALSE))
    }
    if (!wrote) {
      wrote <- isTRUE(tryCatch({
        ggsave(tmp, plot, width = width, height = height, bg = "white", device = "pdf")
        TRUE
      }, error = function(e) FALSE))
    }
    if (!wrote || !file.exists(tmp)) stop("Could not write PDF: ", file)
    if (file.exists(file)) unlink(file)
    if (file.exists(file)) {
      alt <- sub("\\.pdf$", paste0("_", format(Sys.time(), "%H%M%S"), ".pdf"), file)
      file.rename(tmp, alt)
      message("[11] ", basename(file), " locked; wrote ", basename(alt))
    } else {
      file.rename(tmp, file)
    }
  }

  draw_long$family_lab <- unname(fam_lab[draw_long$family])
  draw_long$family_lab <- factor(draw_long$family_lab, levels = unname(fam_lab[families]))
  ann <- sum_tab
  ann$family_lab <- factor(unname(fam_lab[ann$family]), levels = levels(draw_long$family_lab))
  ann$lab <- sprintf(
    "TIPS ρ=%.3f\nnull median=%.3f\nP=%.3f",
    ann$tips_rho, ann$null_rho_median, ann$p_tips_ge_null_rho
  )
  p <- ggplot(draw_long, aes(x = spearman_rho)) +
    geom_histogram(bins = 40, fill = "grey85", colour = "white", linewidth = 0.2) +
    geom_vline(data = ann, aes(xintercept = tips_rho), colour = col_tips, linewidth = 0.7) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab),
              hjust = -0.05, vjust = 1.1, size = 2.8, colour = "grey20", inherit.aes = FALSE) +
    facet_wrap(~family_lab, ncol = 4) +
    labs(
      title = "Same-size null: is the TIPS module unusually fate-informative?",
      subtitle = sprintf(
        "k=%d genes, %d OOF clones, %d draws per family. Red line = TIPS. P = (1 + #{null ρ ≥ TIPS}) / (1 + n).",
        length(in_expr(fold_modules[[1]]$tips_meg)), length(tips_pool$clone_id), n_draw
      ),
      x = "Out-of-fold Spearman vs later Meg fraction",
      y = "Count",
      caption = paste(
        "CTS / MuTrans / Weinreb S3: random k-subsets of those lists (TIPS is not a CTS subset; GATA2 is outside CTS).",
        "PPI nulls: C11 STRING, size k, degree-matched; connected = snowball subgraphs then closest degree sequence.",
        "Meg-marker / node-DE: random k-subsets of a training-only top-50 upregulated HVG pool.",
        "This tests compact selection, not TIPS vs the full list, and does not use S3 as a held-out predictor."
      )
    ) +
    theme_held() +
    theme(strip.background = element_rect(fill = "grey95", colour = NA))
  outp <- file.path(figdir, "11_same_size_null.pdf")
  save_pdf(outp, p, width = 14.4, height = 7.4)
  message("[11] figure -> ", outp)
}
