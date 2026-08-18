## 06_branch_resolvability.R
## Design 8-row Meg vs Baso table. Rows 1–7 are outcome-blind. Row 8 is joined
## last and does not change any Baso threshold.
##
## Row 4 p-value and row 6 bootstrap use Pearson Δ-correlation on existing
## TIPS/STRING edges. That is not a full STRING-reweight per shuffle/bootstrap.
## Observed row 4 also reports edge_change_table() on the real reweighted graphs.

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor"
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_ensure()

n_perm <- as.integer(Sys.getenv("HELDOUT_EDGE_N_PERM", "200"))
n_boot <- as.integer(Sys.getenv("HELDOUT_STAB_N_BOOT", "50"))

sets_rds <- file.path(cfg$results_dir, "rds", "03_gene_sets.rds")
if (file.exists(sets_rds)) {
  sets <- readRDS(sets_rds)
} else {
  sz_f <- file.path(cfg$results_dir, "tables", "03_gene_set_sizes.tsv")
  if (!file.exists(sz_f)) stop("Run 03 first")
  sz <- utils::read.delim(sz_f, stringsAsFactors = FALSE)
  sets <- lapply(setNames(strsplit(as.character(sz$genes), ",", fixed = TRUE), sz$set), norm_genes)
}
if (!file.exists(cfg$cell_tips_table)) stop("Missing TIPS dual-pull table")
tips_df <- utils::read.delim(cfg$cell_tips_table, stringsAsFactors = FALSE)
cf_edges <- linkage_edges(tips_df, "CF")
cm_edges <- linkage_edges(tips_df, "CM")

meg_nodes <- norm_genes(sets$tips_meg_frozen)
baso_nodes <- norm_genes(sets$tips_baso_frozen)
cts <- norm_genes(sets$cts_c11)
td_meg <- norm_genes(sets$mutrans_td_meg)
td_baso <- norm_genes(sets$mutrans_td_baso)
hvg <- if (length(sets$hvg)) norm_genes(sets$hvg) else unique(c(cts, td_meg, td_baso))

cts_sum <- if (file.exists(cfg$cell_cts_summary)) {
  utils::read.csv(cfg$cell_cts_summary, stringsAsFactors = FALSE)
} else NULL
c11_row <- if (!is.null(cts_sum) && "CTS_ID" %in% names(cts_sum)) {
  cts_sum[as.character(cts_sum$CTS_ID) == as.character(cfg$state_id), , drop = FALSE]
} else NULL
biotip_ok <- function(row) {
  if (is.null(row) || !nrow(row)) return(NA)
  mci <- if ("MCI_P" %in% names(row)) row$MCI_P[1] else NA
  icg <- if ("IC_g_local_p" %in% names(row)) row$IC_g_local_p[1] else NA
  ics <- if ("IC_s_local_p" %in% names(row)) row$IC_s_local_p[1] else NA
  loc <- if ("localmax" %in% names(row)) tolower(as.character(row$localmax[1])) == "yes" else NA
  isTRUE(mci < 0.05 && icg < 0.05 && (is.na(ics) || ics < 0.05) && (is.na(loc) || loc))
}
biotip_txt <- function(row) {
  if (is.null(row) || !nrow(row)) return("C11 not in CTS summary")
  sprintf("MCI_P=%s IC_g=%s IC_s=%s localmax=%s pass=%s",
          row$MCI_P[1],
          if ("IC_g_local_p" %in% names(row)) row$IC_g_local_p[1] else NA,
          if ("IC_s_local_p" %in% names(row)) row$IC_s_local_p[1] else NA,
          if ("localmax" %in% names(row)) row$localmax[1] else NA,
          biotip_ok(row))
}

ek <- function(a, b) paste(pmin(toupper(as.character(a)), toupper(as.character(b))),
                           pmax(toupper(as.character(a)), toupper(as.character(b))), sep = "|")

g_notsimp <- load_tips_state_graphs(cfg$cell_tips_graph)
g_w <- load_tips_state_graphs(cfg$cell_tips_graph_weighted)
g_src <- if (!is.null(g_w)) g_w else g_notsimp
g_cp <- graph_for_cluster(if (!is.null(g_notsimp)) g_notsimp else g_src, cfg$state_id)
g_cf <- graph_for_cluster(g_src, cfg$cf_cluster)
g_cm <- graph_for_cluster(g_src, cfg$cm_cluster)
g_cp_w <- graph_for_cluster(g_src, cfg$state_id)

backbone_txt <- function(genes, g) {
  genes <- norm_genes(genes)
  if (is.null(g) || !length(genes)) return("graph missing")
  hit <- intersect(genes, igraph::V(g)$name)
  if (!length(hit)) return(sprintf("0/%d STRING nodes", length(genes)))
  sg <- igraph::induced_subgraph(g, hit)
  n_conn <- sum(igraph::degree(sg) > 0)
  sprintf("%d/%d STRING nodes, %d/%d connected in induced subgraph",
          length(hit), length(genes), n_conn, length(genes))
}

integrity_txt <- function(g, src) {
  if (is.null(g)) return("graph RDS not loaded")
  cc <- igraph::components(g)
  sprintf("%s v=%d e=%d connected=%s n_comp=%d giant=%d",
          src, igraph::vcount(g), igraph::ecount(g), igraph::is_connected(g),
          cc$no, max(cc$csize))
}
g_src_lab <- if (!is.null(g_w)) "weighted HiG" else "notsimplified"

obs_edge_txt <- function(g_from, g_to, edges) {
  dual <- sprintf("dual-pull n=%d n_pos=%d mean(delta>0)=%.3f",
                  nrow(edges), sum(edges$delta > 0),
                  mean(edges$delta[edges$delta > 0], na.rm = TRUE))
  if (is.null(g_from) || is.null(g_to)) return(dual)
  dw <- tryCatch(edge_change_table(g_from, g_to), error = function(e) NULL)
  if (is.null(dw) || !nrow(dw) || !nrow(edges)) return(dual)
  sub <- dw[ek(dw$from, dw$to) %in% ek(edges$from, edges$to), , drop = FALSE]
  if (!nrow(sub)) return(paste0(dual, "; TIPS edges not in reweighted graph"))
  sprintf("edge_change_table on TIPS edges n=%d n_pos=%d mean(delta>0)=%.3f",
          nrow(sub), sum(sub$delta > 0), mean(sub$delta[sub$delta > 0], na.rm = TRUE))
}

fmt_fet <- function(fo) {
  sprintf("k=%d FET p=%.3g OR=%.2f", fo$k, fo$p, fo$or)
}
holly_cts_td <- function(pattern) {
  d <- file.path(cfg$cell_tips_root, "Holly")
  f <- list.files(d, pattern = pattern, full.names = TRUE)
  if (!length(f)) return("")
  x <- utils::read.delim(f[1], stringsAsFactors = FALSE)
  hit <- x[grepl("CTS", x$set_a, ignore.case = TRUE) & grepl("TD|MuTrans", x$set_b, ignore.case = TRUE), ]
  if (!nrow(hit)) hit <- x[grepl("CTS", x$set_b, ignore.case = TRUE) & grepl("TD|MuTrans", x$set_a, ignore.case = TRUE), ]
  if (!nrow(hit)) return("")
  sprintf("; Holly %s p=%.3g k=%d", basename(f[1]), hit$p_value[1], hit$overlap_k[1])
}
row2_meg <- paste0(fmt_fet(fisher_overlap(cts, td_meg, hvg)), holly_cts_td("Mega_cellCTS"))
row2_baso <- paste0(fmt_fet(fisher_overlap(cts, td_baso, hvg)), holly_cts_td("Baso_cellCTS"))

align_edges <- function(edges, rn) {
  map <- setNames(rn, toupper(rn))
  out <- edges
  out$from <- unname(map[toupper(as.character(edges$from))])
  out$to <- unname(map[toupper(as.character(edges$to))])
  out[!is.na(out$from) & !is.na(out$to) & out$from != out$to, , drop = FALSE]
}

mean_abs_dcor <- function(mat, edges, cp, tg) {
  mean(abs(pair_pearson(mat, edges$from, edges$to, tg) -
             pair_pearson(mat, edges$from, edges$to, cp)), na.rm = TRUE)
}

shuffle_p <- function(mat, edges, cp, tg, n_perm, seed) {
  edges <- align_edges(edges, rownames(mat))
  if (!nrow(edges) || length(cp) < 8L || length(tg) < 8L) {
    return(list(text = "too few TIPS edges/cells for permutation", p = NA_real_))
  }
  real <- mean_abs_dcor(mat, edges, cp, tg)
  pool <- unique(c(cp, tg))
  set.seed(seed)
  perm <- replicate(n_perm, {
    tg2 <- sample(pool, length(tg))
    mean_abs_dcor(mat, edges, setdiff(pool, tg2), tg2)
  })
  p <- (1 + sum(perm >= real, na.rm = TRUE)) / (1 + sum(is.finite(perm)))
  list(p = p, text = sprintf("mean|d-cor| on TIPS edges=%.3f perm p=%.3g (n=%d)", real, p, n_perm))
}

string_pairs <- function(g, genes) {
  hit <- intersect(norm_genes(genes), if (is.null(g)) character() else igraph::V(g)$name)
  if (length(hit) < 2L) return(NULL)
  e <- igraph::as_data_frame(igraph::induced_subgraph(g, hit), what = "edges")
  if (!nrow(e)) return(NULL)
  cbind(as.character(e$from), as.character(e$to))
}

boot_stab <- function(mat, pairs, tips_genes, cp, tg, top_n, n_boot, seed) {
  if (is.null(pairs) || !nrow(pairs) || length(cp) < 8L || length(tg) < 8L) {
    return(list(text = "too few STRING pairs/cells for bootstrap", jacc = NA_real_))
  }
  rn <- rownames(mat)
  map <- setNames(rn, toupper(rn))
  pairs[, 1] <- unname(map[toupper(pairs[, 1])])
  pairs[, 2] <- unname(map[toupper(pairs[, 2])])
  pairs <- pairs[complete.cases(pairs) & pairs[, 1] != pairs[, 2], , drop = FALSE]
  if (!nrow(pairs)) return(list(text = "STRING pairs not in expression matrix", jacc = NA_real_))
  pick <- function(cp_ids, tg_ids) {
    dlt <- abs(pair_pearson(mat, pairs[, 1], pairs[, 2], tg_ids) -
                 pair_pearson(mat, pairs[, 1], pairs[, 2], cp_ids))
    top <- pairs[utils::head(order(dlt, decreasing = TRUE, na.last = NA), top_n), , drop = FALSE]
    list(genes = unique(as.vector(top)), edges = ek(top[, 1], top[, 2]))
  }
  set.seed(seed)
  boots <- lapply(seq_len(n_boot), function(i) pick(sample(cp, replace = TRUE), sample(tg, replace = TRUE)))
  jg <- mean(utils::combn(n_boot, 2, function(ij) jaccard(boots[[ij[1]]]$genes, boots[[ij[2]]]$genes)), na.rm = TRUE)
  je <- mean(utils::combn(n_boot, 2, function(ij) jaccard(boots[[ij[1]]]$edges, boots[[ij[2]]]$edges)), na.rm = TRUE)
  freq <- sort(table(unlist(lapply(boots, `[[`, "genes"))) / n_boot, decreasing = TRUE)
  tips_hit <- intersect(norm_genes(tips_genes), toupper(names(freq)))
  tips_txt <- if (!length(tips_hit)) "TIPS genes never selected" else {
    paste(sprintf("%s(%.0f%%)", tips_hit, 100 * as.numeric(freq[tips_hit])), collapse = ",")
  }
  list(
    jacc = jg,
    text = sprintf("cell-boot gene Jaccard=%.2f edge Jaccard=%.2f (B=%d); TIPS freq %s",
                   jg, je, n_boot, tips_txt)
  )
}

jacc_across <- function(md, set_name) {
  tips <- md[md$set == set_name, ]
  sets_l <- strsplit(as.character(tips$genes), ",", fixed = TRUE)
  if (length(sets_l) < 2L) return(NA_real_)
  mean(utils::combn(length(sets_l), 2, function(ij) {
    jaccard(norm_genes(sets_l[[ij[1]]]), norm_genes(sets_l[[ij[2]]]))
  }), na.rm = TRUE)
}
jacc_fold <- jacc_fold_baso <- NA_real_
mod_f <- file.path(cfg$results_dir, "tables", "04_fold_module_genes.tsv")
if (file.exists(mod_f)) {
  md <- utils::read.delim(mod_f, stringsAsFactors = FALSE)
  jacc_fold <- jacc_across(md, "tips_meg")
  jacc_fold_baso <- jacc_across(md, "tips_baso")
}

sh_cf <- sh_cm <- list(text = "Seurat not loaded", p = NA)
st_cf <- st_cm <- list(text = "Seurat not loaded", jacc = NA)
seu_path <- if (file.exists(cfg$seurat_rds_tips_arm)) cfg$seurat_rds_tips_arm else cfg$seurat_rds
if (file.exists(seu_path) && requireNamespace("Seurat", quietly = TRUE) &&
    requireNamespace("igraph", quietly = TRUE)) {
  message("[06] loading Seurat for rows 4 and 6")
  suppressPackageStartupMessages(library(Seurat))
  seu <- readRDS(seu_path)
  cand <- unique(c(cts, meg_nodes, baso_nodes, cf_edges$from, cf_edges$to, cm_edges$from, cm_edges$to))
  mat <- as.matrix(extract_seurat_expr(seu, cand))
  rownames(mat) <- as.character(rownames(mat))
  st <- as.character(seu[[cfg$state_col]][, 1])
  names(st) <- colnames(seu)
  rm(seu); gc()
  cp <- names(st)[st == as.character(cfg$state_id)]
  cf <- names(st)[st == as.character(cfg$cf_cluster)]
  cm <- names(st)[st == as.character(cfg$cm_cluster)]
  sh_cf <- shuffle_p(mat, cf_edges, cp, cf, n_perm, cfg$seed)
  sh_cm <- shuffle_p(mat, cm_edges, cp, cm, n_perm, cfg$seed + 1L)
        sp <- string_pairs(g_cp, cts)
        if (is.null(sp)) {
          g2 <- match_genes_to_rownames(cts, rownames(mat))
          sp <- if (length(g2) >= 2L) t(utils::combn(g2, 2)) else NULL
        }
  st_cf <- boot_stab(mat, sp, meg_nodes, cp, cf, top_n = max(10L, nrow(cf_edges)), n_boot, cfg$seed)
  st_cm <- boot_stab(mat, sp, baso_nodes, cp, cm, top_n = max(10L, nrow(cm_edges)), n_boot, cfg$seed + 2L)
}

jacc <- if (file.exists(cfg$string_iid_jaccard)) {
  utils::read.delim(cfg$string_iid_jaccard, stringsAsFactors = FALSE)
} else NULL
j_cf <- if (!is.null(jacc) && "CF-lean" %in% jacc$Group) jacc$jaccard[jacc$Group == "CF-lean"][1] else NA
j_cm <- if (!is.null(jacc) && "CM-lean" %in% jacc$Group) jacc$jaccard[jacc$Group == "CM-lean"][1] else NA
mc_string <- file.path(cfg$metacell_tips_string, "PPI_weight",
                       "GSE140802_STRING_graph_perState_simplified_combinedweighted.rds")
mc_iid <- file.path(cfg$metacell_tips_iid, "PPI_weight",
                    "GSE140802_STRING_graph_perState_simplified_combinedweighted.rds")
mc_txt <- "no IID arm at leiden_r0_8"
if (file.exists(mc_string) && file.exists(mc_iid)) {
  gs <- load_tips_state_graphs(mc_string)
  gi <- load_tips_state_graphs(mc_iid)
  jv <- function(cl) {
    a <- graph_for_cluster(gs, cl); b <- graph_for_cluster(gi, cl)
    if (is.null(a) || is.null(b)) return(NA_real_)
    jaccard(igraph::V(a)$name, igraph::V(b)$name)
  }
  mc_txt <- sprintf("metacell vertex Jaccard STRING–IID CF11=%.2f CM9=%.2f", jv("11"), jv("9"))
}
row7_meg <- sprintf("CF-lean Jaccard=%.3f; %s", j_cf, mc_txt)
row7_baso <- sprintf("CM-lean Jaccard=%s; %s",
                     if (is.na(j_cm)) "NA" else sprintf("%.3f", j_cm), mc_txt)

wm <- fisher_overlap(meg_nodes, sets$weinreb_meg, hvg)
wb <- fisher_overlap(baso_nodes, sets$weinreb_baso, hvg)
row8_meg <- sprintf("Weinreb FET %s", fmt_fet(wm))
row8_baso <- sprintf("Weinreb FET %s", fmt_fet(wb))
perf_path <- file.path(cfg$results_dir, "tables", "05_method_performance.tsv")
if (file.exists(perf_path)) {
  perf <- utils::read.delim(perf_path, stringsAsFactors = FALSE)
  fmt <- function(arm, method) {
    hit <- perf[perf$endpoint == arm & perf$method == method, ]
    if ("denominator" %in% names(perf)) {
      hit <- hit[hit$denominator == "all_day6_progeny" | is.na(hit$denominator), ]
    }
    if (!nrow(hit)) return("05 not in table")
    sprintf("Spearman=%.3f OR_p=%.3g AUROC=%.3f", hit$spearman_rho[1], hit$or_p[1], hit$auroc[1])
  }
  row8_meg <- paste(fmt("Meg", "tips_meg"), row8_meg, sep = "; ")
  row8_baso <- paste(fmt("Baso", "tips_baso"), row8_baso, sep = "; ")
}

resolv <- data.frame(
  level = c(
    "Transition state",
    "Independent transition support",
    "Backbone coverage",
    "Edge evidence",
    "Network integrity",
    "Resampling stability",
    "Backbone robustness",
    "External validation"
  ),
  diagnostic = c(
    "BioTIP empirical significance (C11 summary table)",
    "CTS–MuTrans overlap (Fisher; Holly FET appended when present)",
    "CTS ∪ TIPS-arm genes represented and connected in C11 STRING graph",
    "edge_change_table on TIPS edges + shuffled descendant labels on those edges",
    "Retained connected component, node and edge counts (weighted HiG, else notsimplified)",
    "Module gene/edge Jaccard across cell bootstraps (STRING-pair |delta-cor|)",
    "STRING–IID Jaccard (metacell 4_9vs11 proxy; no leiden_r0_8 IID arm)",
    "Held-out clone-fate + Weinreb fate-gene overlap"
  ),
  Meg = c(
    biotip_txt(c11_row),
    row2_meg,
    backbone_txt(unique(c(cts, meg_nodes)), g_cp),
    paste(obs_edge_txt(g_cp_w, g_cf, cf_edges), sh_cf$text, sep = "; "),
    integrity_txt(g_cf, g_src_lab),
    paste0(st_cf$text, sprintf("; CV-fold Jaccard=%.3f", jacc_fold)),
    row7_meg,
    row8_meg
  ),
  Baso = c(
    paste0("same C11 BioTIP; pass=", biotip_ok(c11_row)),
    row2_baso,
    backbone_txt(unique(c(cts, baso_nodes)), g_cp),
    paste(obs_edge_txt(g_cp_w, g_cm, cm_edges), sh_cm$text, sep = "; "),
    integrity_txt(g_cm, g_src_lab),
    paste0(st_cm$text, sprintf("; CV-fold Jaccard=%.3f", jacc_fold_baso)),
    row7_baso,
    row8_baso
  ),
  stringsAsFactors = FALSE
)
write_tsv(resolv, file.path(cfg$results_dir, "tables", "06_internal_resolvability.tsv"))
write_tsv(resolv, file.path(cfg$results_dir, "tables", "06_branch_resolvability.tsv"))
if (file.exists(perf_path)) {
  write_tsv(resolv, file.path(cfg$results_dir, "tables", "06_resolvability_with_heldout.tsv"))
  message("[06] eight design rows; held-out joined only on External validation.")
} else {
  message("[06] eight design rows; 05 not run yet so External validation has Weinreb only.")
}

baso_genes <- unique(c(sets$weinreb_baso, td_baso))
audit <- data.frame(
  gene = baso_genes,
  in_Weinreb_Baso = baso_genes %in% sets$weinreb_baso,
  in_MuTrans_TD_4to9 = baso_genes %in% td_baso,
  in_C11_CTS = baso_genes %in% cts,
  in_TIPS_Baso_module = baso_genes %in% baso_nodes,
  in_unweighted_PPI_CM = baso_genes %in% norm_genes(sets$ppi_unweighted_baso),
  stringsAsFactors = FALSE
)
audit$lost_at <- ifelse(
  audit$in_TIPS_Baso_module, "retained_in_TIPS",
  ifelse(audit$in_unweighted_PPI_CM, "TIPS_delta_filter",
  ifelse(audit$in_C11_CTS, "PPI_or_dualpull",
  ifelse(audit$in_MuTrans_TD_4to9 | audit$in_Weinreb_Baso, "not_in_CTS", "other")))
)
write_tsv(audit, file.path(cfg$results_dir, "tables", "06_baso_gene_audit.tsv"))
funnel <- as.data.frame(table(audit$lost_at), stringsAsFactors = FALSE)
names(funnel) <- c("stage", "n_genes")
write_tsv(funnel, file.path(cfg$results_dir, "tables", "06_baso_evidence_funnel.tsv"))
saveRDS(list(resolv = resolv, audit = audit), file.path(cfg$results_dir, "rds", "06_resolvability.rds"))
print(resolv[, c("level", "diagnostic")])
message("[06] Baso genes tracked: ", nrow(audit))
