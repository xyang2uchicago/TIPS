## 00_helpers.R — shared utilities for held-out clone-fate prediction

norm_genes <- function(x) {
  x <- unique(toupper(as.character(x)))
  x[!is.na(x) & nzchar(x) & x != "NA"]
}

as_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

as_int_like <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

match_genes_to_rownames <- function(genes, feature_names) {
  genes <- norm_genes(genes)
  map <- setNames(as.character(feature_names), toupper(as.character(feature_names)))
  hit <- unname(map[genes])
  unique(hit[!is.na(hit)])
}

linkage_set <- function(df, linkage, other) {
  df$from <- toupper(as.character(df$from))
  df$to <- toupper(as.character(df$to))
  nodes <- norm_genes(with(df[df$linkage == linkage, ], c(from, to)))
  shared <- intersect(nodes, norm_genes(with(df[df$linkage == other, ], c(from, to))))
  pos <- norm_genes(with(df[df$linkage == linkage & df$delta > 0, ], c(from, to)))
  setdiff(nodes, setdiff(shared, pos))
}

linkage_edges <- function(df, linkage) {
  sub <- df[df$linkage == linkage, c("from", "to", "w1", "w2", "delta"), drop = FALSE]
  sub$from <- toupper(as.character(sub$from))
  sub$to <- toupper(as.character(sub$to))
  sub
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[write] ", path)
  invisible(path)
}

read_public_meta <- function(path) {
  meta <- utils::read.delim(
    gzfile(path), header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  names(meta) <- make.names(names(meta))
  meta$cell_index0 <- seq_len(nrow(meta)) - 1L
  meta$cell_id <- paste0("cell_", meta$cell_index0)
  meta$Time.point <- as_int_like(meta$Time.point)
  meta$Well <- as_int_like(meta$Well)
  meta
}

read_working_meta <- function(cfg) {
  if (file.exists(cfg$obs_metadata_csv)) {
    meta <- utils::read.csv(cfg$obs_metadata_csv, stringsAsFactors = FALSE, check.names = FALSE)
    if ("X" %in% names(meta) && !"cell_id" %in% names(meta)) {
      names(meta)[names(meta) == "X"] <- "cell_id"
    }
  } else if (file.exists(cfg$seurat_rds)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Need Seurat or obs_metadata.csv")
    }
    seu <- readRDS(cfg$seurat_rds)
    meta <- seu@meta.data
    meta$cell_id <- rownames(meta)
    rm(seu)
    gc()
  } else {
    stop("No working metadata: ", cfg$obs_metadata_csv)
  }
  meta$Time.point <- as_int_like(meta$Time.point)
  meta$Well <- as_int_like(meta$Well)
  meta$cell_index0 <- as_int_like(meta$cell_index0)
  if (!"cell_id" %in% names(meta) || !nzchar(meta$cell_id[1])) {
    meta$cell_id <- paste0("cell_", meta$cell_index0)
  }
  meta[[cfg$state_col]] <- as.character(meta[[cfg$state_col]])
  meta
}

read_clone_matrix <- function(path) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install Matrix")
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  mat <- Matrix::readMM(con)
  mat <- methods::as(mat, "CsparseMatrix")
  mat
}

clone_hits <- function(clone_mat, cell_index0) {
  n <- length(cell_index0)
  n_clones <- integer(n)
  clone_id <- rep(NA_integer_, n)
  n_pub <- nrow(clone_mat)
  ok <- !is.na(cell_index0) & cell_index0 >= 0L & cell_index0 < n_pub
  if (!any(ok)) {
    return(list(clone_id = clone_id, n_clones = n_clones, multi = n_clones > 1L))
  }
  rows <- as.integer(cell_index0[ok]) + 1L
  sub <- clone_mat[rows, , drop = FALSE]
  rs <- as.integer(Matrix::rowSums(sub > 0))
  n_clones[ok] <- rs
  idx <- Matrix::summary(sub)
  if (nrow(idx)) {
    idx <- idx[rs[idx$i] == 1L, , drop = FALSE]
    if (nrow(idx)) clone_id[ok][idx$i] <- as.integer(idx$j)
  }
  list(clone_id = clone_id, n_clones = n_clones, multi = n_clones > 1L)
}

## Cells with >1 clone get clone_id=NA (not the first clone).
assign_clone_ids <- function(clone_mat, cell_index0) {
  h <- clone_hits(clone_mat, cell_index0)
  out <- h$clone_id
  attr(out, "n_multi_clone_cells") <- sum(h$multi)
  attr(out, "n_clones") <- h$n_clones
  attr(out, "multi_clone") <- h$multi
  out
}

## TRUE if cell_index0 has a 1 in any of the given clone columns (1-based).
cells_in_clones <- function(clone_mat, cell_index0, clone_ids) {
  hit <- rep(FALSE, length(cell_index0))
  clone_ids <- unique(as.integer(clone_ids))
  clone_ids <- clone_ids[!is.na(clone_ids) & clone_ids >= 1L & clone_ids <= ncol(clone_mat)]
  n_pub <- nrow(clone_mat)
  ok <- !is.na(cell_index0) & cell_index0 >= 0L & cell_index0 < n_pub
  if (!any(ok) || !length(clone_ids)) return(hit)
  rows <- as.integer(cell_index0[ok]) + 1L
  rs <- Matrix::rowSums(clone_mat[rows, clone_ids, drop = FALSE] > 0)
  hit[ok] <- as.numeric(rs) > 0
  hit
}

is_undiff_label <- function(x) {
  grepl("undiff", as.character(x), ignore.case = TRUE)
}

module_score_matrix <- function(expr, genes, score_fun = "mean") {
  ## expr: genes x cells, rownames = features
  use <- match_genes_to_rownames(genes, rownames(expr))
  if (!length(use)) {
    return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))
  }
  m <- expr[use, , drop = FALSE]
  ## gene-wise z-score across the supplied cells (caller must pass training or
  ## a frozen training transform). Then average.
  zs <- t(scale(t(as.matrix(m))))
  zs[!is.finite(zs)] <- 0
  sc <- if (identical(score_fun, "median")) {
    apply(zs, 2, stats::median, na.rm = TRUE)
  } else {
    colMeans(zs, na.rm = TRUE)
  }
  sc[!is.finite(sc)] <- NA_real_
  sc
}

fit_z_then_score <- function(expr, genes, train_cells, test_cells, score_fun = "mean") {
  use <- match_genes_to_rownames(genes, rownames(expr))
  out <- setNames(rep(NA_real_, length(test_cells)), test_cells)
  if (!length(use) || !length(train_cells)) return(list(score = out, genes_used = use, mu = NULL, sd = NULL))
  train_cells <- intersect(train_cells, colnames(expr))
  test_cells <- intersect(test_cells, colnames(expr))
  m_train <- as.matrix(expr[use, train_cells, drop = FALSE])
  mu <- rowMeans(m_train, na.rm = TRUE)
  sdv <- apply(m_train, 1, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  z_test <- sweep(as.matrix(expr[use, test_cells, drop = FALSE]), 1, mu, "-")
  z_test <- sweep(z_test, 1, sdv, "/")
  z_test[!is.finite(z_test)] <- 0
  sc <- if (identical(score_fun, "median")) {
    apply(z_test, 2, stats::median, na.rm = TRUE)
  } else {
    colMeans(z_test, na.rm = TRUE)
  }
  out[names(sc)] <- sc
  list(score = out, genes_used = use, mu = mu, sd = sdv)
}

pair_pearson <- function(expr, from, to, cells) {
  cells <- intersect(cells, colnames(expr))
  map <- setNames(as.character(rownames(expr)), toupper(as.character(rownames(expr))))
  from_m <- unname(map[toupper(as.character(from))])
  to_m <- unname(map[toupper(as.character(to))])
  g <- unique(c(from_m, to_m))
  g <- g[!is.na(g)]
  if (length(cells) < 5L || length(g) < 2L) {
    return(rep(NA_real_, length(from)))
  }
  m <- as.matrix(expr[g, cells, drop = FALSE])
  vapply(seq_along(from_m), function(i) {
    a <- from_m[i]; b <- to_m[i]
    if (is.na(a) || is.na(b) || !a %in% rownames(m) || !b %in% rownames(m)) return(NA_real_)
    xa <- m[a, ]; xb <- m[b, ]
    if (stats::sd(xa, na.rm = TRUE) == 0 || stats::sd(xb, na.rm = TRUE) == 0) return(NA_real_)
    suppressWarnings(stats::cor(xa, xb, method = "pearson", use = "pairwise.complete.obs"))
  }, numeric(1))
}

extract_seurat_expr <- function(seu, genes, assay = "RNA", layer = "data") {
  feats <- rownames(seu)
  use <- match_genes_to_rownames(genes, feats)
  if (!length(use)) stop("None of the requested genes are in the Seurat object")
  mat <- tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = layer)[use, , drop = FALSE],
    error = function(e) {
      Seurat::GetAssayData(seu, assay = assay, slot = layer)[use, , drop = FALSE]
    }
  )
  mat
}

stratified_clone_folds <- function(clones, n_folds = 5L, seed = 1L, min_meg_pos_per_fold = 8L) {
  stopifnot(all(c("clone_id", "library", "starting_population", "n_meg_d6") %in% names(clones)))
  set.seed(seed)
  clones$meg_pos <- as.integer(clones$n_meg_d6 > 0)
  clones$library <- as_chr(clones$library)
  clones$starting_population <- as_chr(clones$starting_population)
  strata <- paste(clones$library, clones$starting_population, clones$meg_pos, sep = "|")
  fold <- rep(NA_integer_, nrow(clones))
  for (s in unique(strata)) {
    ii <- which(strata == s)
    if (length(ii) == 1L) {
      fold[ii] <- sample.int(n_folds, 1L)
    } else {
      fold[ii] <- sample(rep(seq_len(n_folds), length.out = length(ii)))
    }
  }
  tab <- table(fold, clones$meg_pos)
  meg_per_fold <- if ("1" %in% colnames(tab)) tab[, "1"] else rep(0, n_folds)
  mode <- "grouped_5fold"
  if (any(meg_per_fold < min_meg_pos_per_fold)) {
    ## fall back: stratify only by meg_pos (+ library if large enough)
    message("[folds] Meg-positive clones per fold too small (",
            paste(meg_per_fold, collapse = ","), "); stratifying by Meg-pos only")
    fold <- rep(NA_integer_, nrow(clones))
    for (mp in c(0L, 1L)) {
      ii <- which(clones$meg_pos == mp)
      fold[ii] <- sample(rep(seq_len(n_folds), length.out = length(ii)))
    }
    mode <- "grouped_5fold_megpos"
  }
  clones$fold <- fold
  attr(clones, "fold_mode") <- mode
  attr(clones, "meg_per_fold") <- as.integer(table(factor(clones$fold, levels = seq_len(n_folds)), clones$meg_pos)[, "1"])
  clones
}

quasibinomial_or_per_sd <- function(n_success, n_total, score,
                                    library = NULL, starting_population = NULL) {
  empty <- list(
    or = NA_real_, p = NA_real_, deviance = NA_real_, n = 0L,
    log_or = NA_real_, se = NA_real_, null_deviance = NA_real_,
    or_score_only = NA_real_, p_score_only = NA_real_,
    formula = NA_character_
  )
  ok <- is.finite(score) & is.finite(n_success) & is.finite(n_total) & n_total > 0
  empty$n <- sum(ok)
  if (sum(ok) < 8L) return(empty)
  z <- as.numeric(scale(score[ok]))
  dat <- data.frame(
    y1 = n_success[ok],
    y0 = pmax(n_total[ok] - n_success[ok], 0),
    z = z,
    stringsAsFactors = FALSE
  )
  unpack <- function(fit) {
    if (is.null(fit) || !("z" %in% names(stats::coef(fit)))) {
      return(list(or = NA_real_, log_or = NA_real_, se = NA_real_, p = NA_real_,
                  deviance = NA_real_, null_deviance = NA_real_))
    }
    sm <- summary(fit)
    est <- unname(stats::coef(fit)[["z"]])
    se <- unname(sm$coefficients["z", "Std. Error"])
    p <- tryCatch(unname(sm$coefficients["z", "Pr(>|t|)"]), error = function(e) NA_real_)
    list(or = exp(est), log_or = est, se = se, p = p,
         deviance = stats::deviance(fit), null_deviance = fit$null.deviance)
  }
  fit_z <- function(form) {
    tryCatch(stats::glm(form, data = dat, family = stats::quasibinomial()), error = function(e) NULL)
  }
  so <- unpack(fit_z(cbind(y1, y0) ~ z))
  rhs <- "z"
  if (!is.null(library)) {
    dat$lib <- droplevels(factor(as_chr(library[ok])))
    if (nlevels(dat$lib) >= 2L) rhs <- paste(rhs, "+ lib")
  }
  if (!is.null(starting_population)) {
    dat$sp <- droplevels(factor(as_chr(starting_population[ok])))
    if (nlevels(dat$sp) >= 2L) rhs <- paste(rhs, "+ sp")
  }
  form_txt <- paste("cbind(y1, y0) ~", rhs)
  prim <- if (identical(rhs, "z")) so else unpack(fit_z(stats::as.formula(form_txt)))
  if (!identical(rhs, "z") && !is.finite(prim$or)) prim <- so
  list(
    or = prim$or, log_or = prim$log_or, se = prim$se, p = prim$p,
    deviance = prim$deviance, null_deviance = prim$null_deviance,
    n = sum(ok),
    or_score_only = so$or, p_score_only = so$p,
    formula = form_txt
  )
}

spearman_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5L) return(list(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = unname(ct$p.value))
}

roc_auc <- function(label, score) {
  ok <- is.finite(score) & !is.na(label)
  y <- as.integer(label[ok]); s <- score[ok]
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 < 1L || n0 < 1L) return(NA_real_)
  r <- rank(s, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

pr_auc <- function(label, score) {
  ok <- is.finite(score) & !is.na(label)
  y <- as.integer(label[ok]); s <- score[ok]
  n1 <- sum(y == 1L)
  if (n1 < 1L || sum(y == 0L) < 1L) return(NA_real_)
  o <- order(s, decreasing = TRUE)
  y <- y[o]
  tp <- cumsum(y); fp <- cumsum(1L - y)
  prec <- tp / pmax(tp + fp, 1)
  rec <- tp / n1
  rec <- c(0, rec); prec <- c(1, prec)
  sum((rec[-1] - rec[-length(rec)]) * (prec[-1] + prec[-length(prec)]) / 2)
}

log_loss_binom <- function(n_success, n_total, score) {
  ok <- is.finite(score) & n_total > 0
  if (!sum(ok)) return(NA_real_)
  z <- as.numeric(scale(score[ok]))
  y <- cbind(n_success[ok], pmax(n_total[ok] - n_success[ok], 0))
  fit <- tryCatch(stats::glm(y ~ z, family = stats::binomial()), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  p <- stats::fitted(fit)
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  k <- n_success[ok]; n <- n_total[ok]
  -mean(k * log(p) + (n - k) * log(1 - p)) / mean(n)
}

heldout_code_hash <- function(cfg) {
  files <- sort(c(
    file.path(cfg$code_dir, "03b_build_tips_module_fold.R"),
    list.files(cfg$tips_arm_code_dir, pattern = "\\.[Rr]$", full.names = TRUE)
  ))
  files <- files[file.exists(files)]
  if (!length(files)) return(NA_character_)
  paste(unname(tools::md5sum(files)), collapse = "")
}

heldout_lock_fingerprint <- function(cfg) {
  if (!file.exists(cfg$lock_file)) return(NA_character_)
  unname(tools::md5sum(cfg$lock_file))
}

fold_cache_ok <- function(man, test_clone_ids, fold_id, repeat_id, cfg) {
  if (!is.list(man) || is.null(man$test_clone_ids) || is.null(man$code_hash)) return(FALSE)
  if (!identical(man$status, "success") && !identical(man$status, "empty_module")) return(FALSE)
  identical(sort(as.integer(man$test_clone_ids)), sort(as.integer(test_clone_ids))) &&
    identical(as.integer(man$fold_id), as.integer(fold_id)) &&
    identical(as.integer(man$repeat_id), as.integer(repeat_id)) &&
    identical(man$tips_mode, cfg$tips_mode) &&
    identical(as.character(man$state_id), as.character(cfg$state_id)) &&
    identical(as.character(man$cf_cluster), as.character(cfg$cf_cluster)) &&
    identical(as.character(man$cm_cluster), as.character(cfg$cm_cluster)) &&
    identical(man$lock_fingerprint, heldout_lock_fingerprint(cfg)) &&
    identical(man$code_hash, heldout_code_hash(cfg))
}

session_stamp <- function(cfg, extra = list()) {
  info <- c(
    list(
      time = as.character(Sys.time()),
      R = R.version.string,
      tips_mode = cfg$tips_mode,
      state = paste(cfg$state_col, cfg$state_id, sep = "="),
      seed = cfg$seed
    ),
    extra
  )
  write_tsv(as.data.frame(info, stringsAsFactors = FALSE), file.path(cfg$results_dir, "session_stamp.tsv"))
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (!length(a) && !length(b)) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

fisher_overlap <- function(a, b, universe) {
  u <- norm_genes(universe)
  a <- intersect(norm_genes(a), u)
  b <- intersect(norm_genes(b), u)
  ft <- tryCatch(stats::fisher.test(table(u %in% a, u %in% b)), error = function(e) NULL)
  list(
    k = length(intersect(a, b)), n_a = length(a), n_b = length(b),
    p = if (is.null(ft)) NA_real_ else unname(ft$p.value),
    or = if (is.null(ft) || is.null(ft$estimate)) NA_real_ else unname(ft$estimate)
  )
}

## Same merge/delta logic as celltype_specific_weight_v10.R::edge_change_table.
## corexp_sign is used when present; unsigned weights otherwise.
edge_change_table <- function(g1, g2, weight_attr = "weight", missing_as = 0, undirected = TRUE) {
  if (is.null(igraph::V(g1)$name) || is.null(igraph::V(g2)$name)) {
    stop("Both graphs must have V(g)$name")
  }
  extract_edges <- function(g) {
    e <- igraph::as_data_frame(g, what = "edges")
    e$from <- as.character(e$from)
    e$to <- as.character(e$to)
    if (undirected) {
      a <- pmin(e$from, e$to)
      b <- pmax(e$from, e$to)
      e$from <- a
      e$to <- b
    }
    w <- if (weight_attr %in% names(e)) e[[weight_attr]] else rep(1, nrow(e))
    w[is.na(w)] <- 0
    if ("corexp_sign" %in% names(e)) {
      sg <- e$corexp_sign
      sg[is.na(sg)] <- "positive"
      w <- ifelse(sg == "positive", w, -w)
    }
    data.frame(from = e$from, to = e$to, w = w, stringsAsFactors = FALSE)
  }
  collapse_fun <- function(df) {
    key <- paste(df$from, df$to, sep = "||")
    idx <- tapply(seq_len(nrow(df)), key, function(ii) ii[which.max(abs(df$w[ii]))])
    df[unlist(idx), , drop = FALSE]
  }
  e1 <- collapse_fun(extract_edges(g1)); names(e1)[3] <- "w1"
  e2 <- collapse_fun(extract_edges(g2)); names(e2)[3] <- "w2"
  m <- merge(e1, e2, by = c("from", "to"), all = TRUE)
  m$w1[is.na(m$w1)] <- missing_as
  m$w2[is.na(m$w2)] <- missing_as
  m$delta <- m$w2 - m$w1
  m$abs_delta <- abs(m$delta)
  m$direction <- ifelse(m$delta > 1e-10, "increase",
    ifelse(m$delta < -1e-10, "decrease", "unchanged"))
  m$status <- ifelse(m$w1 == missing_as & m$w2 != missing_as, "gained",
    ifelse(m$w1 != missing_as & m$w2 == missing_as, "lost",
      ifelse(m$w1 != m$w2, "changed", "unchanged")))
  m <- m[order(m$abs_delta, decreasing = TRUE), ]
  m$rank <- seq_len(nrow(m))
  rownames(m) <- NULL
  m
}

load_tips_state_graphs <- function(path) {
  if (!file.exists(path) || !requireNamespace("igraph", quietly = TRUE)) return(NULL)
  gl <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(gl) || !length(gl)) return(NULL)
  gl
}

graph_for_cluster <- function(gl, cluster) {
  if (is.null(gl)) return(NULL)
  keys <- c(paste0("HiG_", cluster), as.character(cluster), paste0("CTS_", cluster))
  hit <- intersect(keys, names(gl))
  if (!length(hit)) return(NULL)
  g <- gl[[hit[1]]]
  igraph::V(g)$name <- toupper(as.character(igraph::V(g)$name))
  g
}
