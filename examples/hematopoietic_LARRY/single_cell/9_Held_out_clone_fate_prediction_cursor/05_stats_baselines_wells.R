## 05_stats_baselines_wells.R
## Primary quasibinomial: cbind(n_meg, n_d6) ~ scale(score) + library + starting_population.
## Score-only OR is a sensitivity column. Log-loss is not reported (it was in-sample).
## Primary denominator: all recovered day-6 progeny. Sensitivity: mature (non-undiff) progeny.

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor"
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_ensure()

pred_rds <- file.path(cfg$results_dir, "rds", "04_oof_predictions.rds")
if (!file.exists(pred_rds)) stop("Run 04 first")
obj <- readRDS(pred_rds)
pred <- obj$pred
rand_long <- if (!is.null(obj$rand_long)) obj$rand_long else {
  rf <- file.path(cfg$results_dir, "tables", "04_random_set_scores.tsv")
  if (file.exists(rf)) utils::read.delim(rf, stringsAsFactors = FALSE) else data.frame()
}
succ <- obj$success
n_eligible <- if (!is.null(obj$n_eligible)) obj$n_eligible else length(unique(pred$clone_id))

eval_endpoint <- function(dd, n_s, n_t, frac, label_pos) {
  if (!nrow(dd) || !n_s %in% names(dd) || !n_t %in% names(dd)) {
    return(data.frame(
      n_eligible = n_eligible, n_clones = 0L, n_scored = 0L, n_positive = 0L,
      or_per_sd = NA_real_, or_p = NA_real_, or_score_only = NA_real_, p_score_only = NA_real_,
      glm_formula = NA_character_, quasibinomial_deviance = NA_real_,
      spearman_rho = NA_real_, spearman_p = NA_real_, auroc = NA_real_, auprc = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  lib <- if ("library" %in% names(dd)) dd$library else NULL
  sp <- if ("starting_population" %in% names(dd)) dd$starting_population else NULL
  qb <- quasibinomial_or_per_sd(dd[[n_s]], dd[[n_t]], dd$score, library = lib, starting_population = sp)
  spm <- spearman_safe(dd$score, dd[[frac]])
  lab <- if (is.logical(dd[[label_pos]]) || is.numeric(dd[[label_pos]])) dd[[label_pos]] else dd[[label_pos]] %in% c(TRUE, "TRUE", 1, "1")
  data.frame(
    n_eligible = n_eligible,
    n_clones = nrow(dd),
    n_scored = sum(is.finite(dd$score)),
    n_positive = sum(lab, na.rm = TRUE),
    or_per_sd = qb$or, or_p = qb$p,
    or_score_only = qb$or_score_only, p_score_only = qb$p_score_only,
    glm_formula = qb$formula,
    quasibinomial_deviance = qb$deviance,
    spearman_rho = spm$rho, spearman_p = spm$p,
    auroc = roc_auc(lab, dd$score),
    auprc = pr_auc(lab, dd$score),
    stringsAsFactors = FALSE
  )
}

meg_methods <- c(
  "tips_meg", "cts_c11", "mutrans_td_meg", "ppi_unweighted_meg",
  "tips_meg_frozen", "metacell_tips_meg", "meg_markers_train", "node_de_megup"
)
baso_methods <- c("tips_baso", "mutrans_td_baso", "ppi_unweighted_baso", "cts_c11")

rows <- list()
for (m in intersect(meg_methods, unique(pred$method))) {
  dd <- pred[pred$method == m, ]
  r <- eval_endpoint(dd, "n_meg_d6", "n_d6", "frac_meg_d6", "meg_positive")
  r$method <- m
  r$endpoint <- "Meg"
  r$denominator <- "all_day6_progeny"
  rows[[length(rows) + 1L]] <- r
  if ("n_mature_d6" %in% names(dd) && "frac_meg_mature_d6" %in% names(dd)) {
    rmature <- eval_endpoint(dd, "n_meg_d6", "n_mature_d6", "frac_meg_mature_d6", "meg_positive")
    rmature$method <- m
    rmature$endpoint <- "Meg"
    rmature$denominator <- "mature_day6_progeny"
    rows[[length(rows) + 1L]] <- rmature
  }
}
for (m in intersect(baso_methods, unique(pred$method))) {
  dd <- pred[pred$method == m, ]
  r <- eval_endpoint(dd, "n_baso_d6", "n_d6", "frac_baso_d6", "baso_positive")
  r$method <- m
  r$endpoint <- "Baso"
  r$denominator <- "all_day6_progeny"
  rows[[length(rows) + 1L]] <- r
  if ("n_mature_d6" %in% names(dd) && "frac_baso_mature_d6" %in% names(dd)) {
    rmature <- eval_endpoint(dd, "n_baso_d6", "n_mature_d6", "frac_baso_mature_d6", "baso_positive")
    rmature$method <- m
    rmature$endpoint <- "Baso"
    rmature$denominator <- "mature_day6_progeny"
    rows[[length(rows) + 1L]] <- rmature
  }
}
perf <- do.call(rbind, rows)
perf <- perf[, c("endpoint", "denominator", "method", setdiff(names(perf), c("endpoint", "denominator", "method")))]
write_tsv(perf, file.path(cfg$results_dir, "tables", "05_method_performance.tsv"))
print(perf[perf$denominator == "all_day6_progeny", ])

## --- random-set null (each draw evaluated separately) ---
if (nrow(rand_long) && "draw" %in% names(rand_long)) {
  tips_row <- perf[perf$method == "tips_meg" & perf$endpoint == "Meg" & perf$denominator == "all_day6_progeny", ][1, ]
  draw_ids <- sort(unique(rand_long$draw))
  rand_metrics <- do.call(rbind, lapply(draw_ids, function(d) {
    dd <- rand_long[rand_long$draw == d, ]
    r <- eval_endpoint(dd, "n_meg_d6", "n_d6", "frac_meg_d6", "meg_positive")
    r$draw <- d
    r
  }))
  write_tsv(rand_metrics, file.path(cfg$results_dir, "tables", "05_random_set_metrics.tsv"))
  p_rho <- (1 + sum(rand_metrics$spearman_rho >= tips_row$spearman_rho, na.rm = TRUE)) / (1 + sum(is.finite(rand_metrics$spearman_rho)))
  p_or <- (1 + sum(rand_metrics$or_per_sd >= tips_row$or_per_sd, na.rm = TRUE)) / (1 + sum(is.finite(rand_metrics$or_per_sd)))
  rand_sum <- data.frame(
    endpoint = "Meg", method = "random_matched",
    n_draws = length(draw_ids),
    spearman_median = stats::median(rand_metrics$spearman_rho, na.rm = TRUE),
    spearman_q025 = stats::quantile(rand_metrics$spearman_rho, 0.025, na.rm = TRUE, names = FALSE),
    spearman_q975 = stats::quantile(rand_metrics$spearman_rho, 0.975, na.rm = TRUE, names = FALSE),
    or_median = stats::median(rand_metrics$or_per_sd, na.rm = TRUE),
    or_q025 = stats::quantile(rand_metrics$or_per_sd, 0.025, na.rm = TRUE, names = FALSE),
    or_q975 = stats::quantile(rand_metrics$or_per_sd, 0.975, na.rm = TRUE, names = FALSE),
    tips_spearman = tips_row$spearman_rho,
    tips_or = tips_row$or_per_sd,
    p_tips_ge_random_spearman = p_rho,
    p_tips_ge_random_or = p_or,
    stringsAsFactors = FALSE
  )
  write_tsv(rand_sum, file.path(cfg$results_dir, "tables", "05_tips_vs_random.tsv"))
}

## --- bootstrap CI for primary TIPS Meg ---
tips <- pred[pred$method == "tips_meg", ]
set.seed(cfg$seed)
boot_or <- replicate(cfg$n_bootstrap, {
  ii <- sample.int(nrow(tips), replace = TRUE)
  quasibinomial_or_per_sd(tips$n_meg_d6[ii], tips$n_d6[ii], tips$score[ii],
                          library = tips$library[ii], starting_population = tips$starting_population[ii])$or
})
boot_rho <- replicate(cfg$n_bootstrap, {
  ii <- sample.int(nrow(tips), replace = TRUE)
  spearman_safe(tips$score[ii], tips$frac_meg_d6[ii])$rho
})
obs_or <- perf$or_per_sd[perf$method == "tips_meg" & perf$endpoint == "Meg" & perf$denominator == "all_day6_progeny"][1]
obs_rho <- perf$spearman_rho[perf$method == "tips_meg" & perf$endpoint == "Meg" & perf$denominator == "all_day6_progeny"][1]
boot <- data.frame(
  stat = c("or_per_sd", "spearman_rho"),
  point = c(obs_or, obs_rho),
  ci95_lo = c(stats::quantile(boot_or, 0.025, na.rm = TRUE), stats::quantile(boot_rho, 0.025, na.rm = TRUE)),
  ci95_hi = c(stats::quantile(boot_or, 0.975, na.rm = TRUE), stats::quantile(boot_rho, 0.975, na.rm = TRUE)),
  n_boot = cfg$n_bootstrap,
  stringsAsFactors = FALSE
)
write_tsv(boot, file.path(cfg$results_dir, "tables", "05_tips_meg_bootstrap.tsv"))

## --- paired clone-level bootstrap of TIPS − baseline ---
## Simple percentile CIs (shared seed per baseline). Manuscript-facing version
## with shared indices, BCa, bootstrap p, AUROC, and figures: 10_clone_bootstrap.R
paired_rows <- list()
base_m <- setdiff(intersect(meg_methods, unique(pred$method)), "tips_meg")
wide <- tips[, c("clone_id", "score", "n_meg_d6", "n_d6", "frac_meg_d6", "library", "starting_population")]
names(wide)[2] <- "score_tips"
for (m in base_m) {
  b <- pred[pred$method == m, c("clone_id", "score")]
  names(b)[2] <- "score_base"
  w <- merge(wide, b, by = "clone_id")
  if (nrow(w) < 8L) next
  set.seed(cfg$seed + 13L)
  d_rho <- replicate(cfg$n_bootstrap, {
    ii <- sample.int(nrow(w), replace = TRUE)
    spearman_safe(w$score_tips[ii], w$frac_meg_d6[ii])$rho -
      spearman_safe(w$score_base[ii], w$frac_meg_d6[ii])$rho
  })
  d_or <- replicate(cfg$n_bootstrap, {
    ii <- sample.int(nrow(w), replace = TRUE)
    quasibinomial_or_per_sd(w$n_meg_d6[ii], w$n_d6[ii], w$score_tips[ii],
                            library = w$library[ii], starting_population = w$starting_population[ii])$or -
      quasibinomial_or_per_sd(w$n_meg_d6[ii], w$n_d6[ii], w$score_base[ii],
                            library = w$library[ii], starting_population = w$starting_population[ii])$or
  })
  paired_rows[[length(paired_rows) + 1L]] <- data.frame(
    baseline = m, n_clones = nrow(w),
    delta_spearman = obs_rho - perf$spearman_rho[perf$method == m & perf$endpoint == "Meg" & perf$denominator == "all_day6_progeny"][1],
    delta_spearman_ci_lo = stats::quantile(d_rho, 0.025, na.rm = TRUE),
    delta_spearman_ci_hi = stats::quantile(d_rho, 0.975, na.rm = TRUE),
    delta_or = obs_or - perf$or_per_sd[perf$method == m & perf$endpoint == "Meg" & perf$denominator == "all_day6_progeny"][1],
    delta_or_ci_lo = stats::quantile(d_or, 0.025, na.rm = TRUE),
    delta_or_ci_hi = stats::quantile(d_or, 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
if (length(paired_rows)) {
  write_tsv(do.call(rbind, paired_rows), file.path(cfg$results_dir, "tables", "05_tips_minus_baseline_bootstrap.tsv"))
}

## --- permutation of fate labels within library ---
perm_one <- function() {
  out <- tips
  for (lib in unique(out$library)) {
    ii <- which(out$library == lib)
    jj <- sample(ii)
    out$n_meg_d6[ii] <- tips$n_meg_d6[jj]
    out$n_d6[ii] <- tips$n_d6[jj]
    out$frac_meg_d6[ii] <- tips$frac_meg_d6[jj]
    out$meg_positive[ii] <- tips$meg_positive[jj]
  }
  list(
    or = quasibinomial_or_per_sd(out$n_meg_d6, out$n_d6, out$score,
                                 library = out$library, starting_population = out$starting_population)$or,
    rho = spearman_safe(out$score, out$frac_meg_d6)$rho
  )
}
set.seed(cfg$seed + 7L)
perm <- replicate(cfg$n_perm, perm_one(), simplify = FALSE)
perm_or <- vapply(perm, `[[`, numeric(1), "or")
perm_rho <- vapply(perm, `[[`, numeric(1), "rho")
perm_tab <- data.frame(
  stat = c("or_per_sd", "spearman_rho"),
  observed = c(obs_or, obs_rho),
  n_perm = cfg$n_perm,
  p_ge_observed = c(
    mean(perm_or >= obs_or, na.rm = TRUE),
    mean(perm_rho >= obs_rho, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write_tsv(perm_tab, file.path(cfg$results_dir, "tables", "05_tips_meg_permutation.tsv"))

## --- two-well replication (well-specific positivity for AUROC/AUPRC) ---
well_eval <- function(dd, n_s, n_t, frac, pos_col, tag) {
  r <- eval_endpoint(dd, n_s, n_t, frac, pos_col)
  r$well <- tag
  r
}
tips$meg_pos_well1 <- tips$n_meg_d6_well1 >= 1L
tips$meg_pos_well2 <- tips$n_meg_d6_well2 >= 1L
w1 <- tips[!is.na(tips$frac_meg_d6_well1) & tips$n_d6_well1 > 0, ]
w2 <- tips[!is.na(tips$frac_meg_d6_well2) & tips$n_d6_well2 > 0, ]
both <- tips[tips$in_both_later_wells %in% c(TRUE, "TRUE", "true", 1, "1"), ]
well_rows <- rbind(
  well_eval(w1, "n_meg_d6_well1", "n_d6_well1", "frac_meg_d6_well1", "meg_pos_well1", "well1"),
  well_eval(w2, "n_meg_d6_well2", "n_d6_well2", "frac_meg_d6_well2", "meg_pos_well2", "well2")
)
conc <- if (nrow(both) >= 5L) {
  spearman_safe(both$frac_meg_d6_well1, both$frac_meg_d6_well2)
} else list(rho = NA_real_, p = NA_real_)
well_rows$well1_well2_spearman <- conc$rho
well_rows$well1_well2_spearman_p <- conc$p
write_tsv(well_rows, file.path(cfg$results_dir, "tables", "05_two_well_replication.tsv"))

if (requireNamespace("aod", quietly = TRUE) && nrow(tips) >= 8L) {
  z <- as.numeric(scale(tips$score))
  dat <- data.frame(y = tips$n_meg_d6, n = tips$n_d6, z = z)
  bb <- tryCatch(aod::betabin(cbind(y, n - y) ~ z, random = ~1, data = dat), error = function(e) NULL)
  if (!is.null(bb)) {
    cf <- tryCatch(summary(bb)@Coef, error = function(e) NULL)
    write_tsv(
      data.frame(note = "aod::betabin score-only", or_per_sd = if (!is.null(cf)) exp(cf[2, 1]) else NA_real_),
      file.path(cfg$results_dir, "tables", "05_tips_meg_betabinomial.tsv")
    )
  }
}

if (!is.null(succ)) {
  write_tsv(succ, file.path(cfg$results_dir, "tables", "05_fold_module_success.tsv"))
}

saveRDS(list(perf = perf, boot = boot, perm = perm_tab, wells = well_rows),
        file.path(cfg$results_dir, "rds", "05_stats.rds"))
message("[05] Primary endpoint is Meg / all recovered day-6 progeny. Mature-only is a sensitivity.")
message("[05] TIPS Meg OR=", signif(obs_or, 3),
        " Spearman=", signif(obs_rho, 3),
        " perm p(rho)=", signif(perm_tab$p_ge_observed[2], 3),
        " scored=", sum(is.finite(tips$score)), "/", n_eligible)
