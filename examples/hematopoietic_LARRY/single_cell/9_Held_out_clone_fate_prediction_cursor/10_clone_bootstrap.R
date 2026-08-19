## 10_clone_bootstrap.R
##
## Manuscript question: is held-out TIPS Spearman (ρ ≈ 0.42) meaningfully
## larger than a baseline ρ (CTS ≈ 0.33, unweighted PPI, training Meg markers,
## node DE, MuTrans TD)?
##
## Method: paired clone bootstrap. Out-of-fold scores are treated as fixed;
## clones are the resampling unit. The same bootstrap index vector is used for
## TIPS and every baseline, so Δ = stat(TIPS) − stat(baseline) is paired.
## This is sampling variability of the pooled OOF association, not a refit of
## TIPS / DEG / MuTrans gene lists.
##
## Primary: Spearman ρ vs later Meg fraction (all recovered day-6 progeny).
## Secondary: AUROC for meg_positive (n_meg_d6 >= 1).
## Sensitivity: quasibinomial OR per SD (library + starting population).
##
## If the TIPS−CTS 95% CI for Δρ includes 0, do not claim predictive
## superiority; keep the selling point on branch-specific resolution.
##
##   Sys.setenv(HELDOUT_RUN = "10")
##   source(".../run_heldout_pipeline.R")
## or:
##   source(".../10_clone_bootstrap.R")

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = here::here("examples", "hematopoietic_LARRY", "single_cell", "9_Held_out_clone_fate_prediction_cursor")
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_ensure()

pred_rds <- file.path(cfg$results_dir, "rds", "04_oof_predictions.rds")
if (!file.exists(pred_rds)) stop("Run 04 first: missing ", pred_rds)
obj <- readRDS(pred_rds)
pred <- obj$pred

n_boot <- as.integer(Sys.getenv("HELDOUT_N_BOOT_DELTA", "5000"))
if (!is.finite(n_boot) || n_boot < 200L) n_boot <- 5000L
do_or <- !(Sys.getenv("HELDOUT_BOOT_OR", "1") %in% c("0", "false", "FALSE"))
seed_delta <- cfg$seed + 1010L

focal <- c(
  cts_c11 = "Fixed C11 CTS",
  ppi_unweighted_meg = "Unweighted PPI",
  meg_markers_train = "Training Meg markers",
  node_de_megup = "Node-level DE",
  mutrans_td_meg = "MuTrans TD"
)
method_order <- c("tips_meg", names(focal))
method_lab <- c(tips_meg = "TIPS", focal)

need <- c("tips_meg", names(focal))
have <- unique(pred$method)
miss <- setdiff(need, have)
if (length(miss)) stop("Missing methods in 04 predictions: ", paste(miss, collapse = ", "))

as_pos <- function(x) {
  as.integer(x %in% c(TRUE, "TRUE", "true", 1, "1"))
}

wide_from_pred <- function(pred, methods) {
  tips <- pred[pred$method == "tips_meg", , drop = FALSE]
  keep <- c(
    "clone_id", "score", "n_meg_d6", "n_d6", "frac_meg_d6", "meg_positive",
    "library", "starting_population"
  )
  keep <- intersect(keep, names(tips))
  w <- tips[, keep, drop = FALSE]
  names(w)[names(w) == "score"] <- "tips_meg"
  for (m in setdiff(methods, "tips_meg")) {
    b <- pred[pred$method == m, c("clone_id", "score"), drop = FALSE]
    names(b)[2] <- m
    w <- merge(w, b, by = "clone_id")
  }
  w$meg_positive <- as_pos(w$meg_positive)
  ok <- is.finite(w$frac_meg_d6) & is.finite(w$n_d6) & w$n_d6 > 0
  for (m in methods) ok <- ok & is.finite(w[[m]])
  w[ok, , drop = FALSE]
}

percentile_ci <- function(x, probs = c(0.025, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) < 20L) return(c(NA_real_, NA_real_))
  unname(stats::quantile(x, probs = probs, names = FALSE, type = 7, na.rm = TRUE))
}

## Davison & Hinkley (1997) one-sided bootstrap p for H1: Δ > 0.
## Recorded in tables only. A one-sided superiority test was not prespecified;
## inference is the two-sided 95% CI (whether it includes 0).
boot_p_gt0 <- function(d) {
  d <- d[is.finite(d)]
  if (!length(d)) return(NA_real_)
  (1 + sum(d <= 0)) / (1 + length(d))
}

boot_p_two_sided <- function(d) {
  d <- d[is.finite(d)]
  if (!length(d)) return(NA_real_)
  n <- length(d)
  p_le <- (1 + sum(d <= 0)) / (1 + n)
  p_ge <- (1 + sum(d >= 0)) / (1 + n)
  min(1, 2 * min(p_le, p_ge))
}

## Efron BCa using jackknife acceleration. Falls back to percentile if z0/a
## are not finite (e.g. all bootstrap draws on one side of the point estimate).
bca_ci <- function(theta_obs, theta_boot, theta_jack, probs = c(0.025, 0.975)) {
  perc <- percentile_ci(theta_boot, probs)
  tb <- theta_boot[is.finite(theta_boot)]
  tj <- theta_jack[is.finite(theta_jack)]
  if (length(tb) < 50L || length(tj) < 8L || !is.finite(theta_obs)) return(perc)
  prop_less <- mean(tb < theta_obs)
  if (prop_less <= 0 || prop_less >= 1) return(perc)
  z0 <- stats::qnorm(prop_less)
  tdot <- mean(tj)
  d <- tdot - tj
  den <- sum(d^2)
  if (den <= 0) return(perc)
  a <- sum(d^3) / (6 * den^1.5)
  if (!is.finite(a) || !is.finite(z0)) return(perc)
  za <- stats::qnorm(probs)
  zba <- z0 + za
  adj <- z0 + zba / (1 - a * zba)
  p <- stats::pnorm(adj)
  p <- pmin(pmax(p, 1 / (length(tb) + 1)), length(tb) / (length(tb) + 1))
  unname(stats::quantile(tb, probs = p, names = FALSE, type = 7))
}

## Williams / Steiger test of two dependent Spearman correlations that share Y.
## Applied to ranks so the null is Δρ_Spearman = 0.
steiger_spearman <- function(x, z, y) {
  empty <- list(r_xz = NA_real_, t = NA_real_, df = NA_real_, p = NA_real_)
  ok <- is.finite(x) & is.finite(z) & is.finite(y)
  if (sum(ok) < 8L) return(empty)
  x <- rank(x[ok]); z <- rank(z[ok]); y <- rank(y[ok])
  n <- length(x)
  rxy <- suppressWarnings(stats::cor(x, y))
  rzy <- suppressWarnings(stats::cor(z, y))
  rxz <- suppressWarnings(stats::cor(x, z))
  if (!all(is.finite(c(rxy, rzy, rxz)))) return(empty)
  num <- (rxy - rzy) * sqrt((n - 3) * (1 + rxz))
  den <- 2 * (1 - rxy^2 - rzy^2 - rxz^2 + 2 * rxy * rzy * rxz)
  if (!is.finite(den) || den <= 0) return(empty)
  tstat <- num / sqrt(den)
  df <- n - 3
  list(
    r_xz = unname(rxz),
    t = unname(tstat),
    df = df,
    p = 2 * stats::pt(-abs(tstat), df = df)
  )
}

or_on_idx <- function(w, ii, score) {
  quasibinomial_or_per_sd(
    w$n_meg_d6[ii], w$n_d6[ii], score[ii],
    library = w$library[ii],
    starting_population = w$starting_population[ii]
  )$or
}

rds10 <- file.path(cfg$results_dir, "rds", "10_clone_bootstrap.rds")
fig_only <- Sys.getenv("HELDOUT_FIGURES_ONLY", "0") %in% c("1", "true", "TRUE")
if (fig_only) {
  if (!file.exists(rds10)) stop("HELDOUT_FIGURES_ONLY=1 but missing ", rds10)
  obj10 <- readRDS(rds10)
  w <- obj10$wide
  n <- nrow(w)
  n_pos <- sum(w$meg_positive == 1L)
  n_boot <- obj10$n_boot
  seed_delta <- obj10$seed
  boot <- obj10$boot
  levels_tab <- obj10$levels
  delta_tab <- obj10$delta
  verdict <- obj10$verdict
  cts_rho <- delta_tab[delta_tab$baseline == "cts_c11" & delta_tab$statistic == "spearman_rho", ][1, ]
  message("[10] reused bootstrap draws from ", rds10)
} else {

w <- wide_from_pred(pred, need)
n <- nrow(w)
if (n < 20L) stop("Too few complete clones for bootstrap: ", n)
score_mat <- as.matrix(w[, need, drop = FALSE])
colnames(score_mat) <- need
y <- w$frac_meg_d6
lab <- w$meg_positive
n_pos <- sum(lab == 1L)
pos_idx <- which(lab == 1L)
neg_idx <- which(lab == 0L)

rho_on_idx <- function(ii) {
  vapply(seq_len(ncol(score_mat)), function(j) {
    x <- score_mat[ii, j]
    yy <- y[ii]
    ok <- is.finite(x) & is.finite(yy)
    if (sum(ok) < 5L) return(NA_real_)
    suppressWarnings(unname(stats::cor(x[ok], yy[ok], method = "spearman")))
  }, numeric(1))
}
auc_on_idx <- function(ii) {
  vapply(seq_len(ncol(score_mat)), function(j) {
    roc_auc(lab[ii], score_mat[ii, j])
  }, numeric(1))
}

obs_rho <- rho_on_idx(seq_len(n))
obs_auc <- auc_on_idx(seq_len(n))
names(obs_rho) <- need
names(obs_auc) <- need
obs_or <- setNames(rep(NA_real_, length(need)), need)
if (do_or) {
  for (m in need) {
    obs_or[[m]] <- or_on_idx(w, seq_len(n), w[[m]])
  }
}

message("[10] n_clones=", n, " n_meg_positive=", n_pos,
        " n_boot=", n_boot, " seed=", seed_delta)
message("[10] observed Spearman: ",
        paste(sprintf("%s=%.3f", method_lab[need], obs_rho), collapse = "; "))

## One shared index matrix so every method and every statistic sees the same clones.
set.seed(seed_delta)
boot_idx <- replicate(n_boot, sample.int(n, replace = TRUE))
set.seed(seed_delta + 1L)
boot_idx_strat <- replicate(n_boot, c(
  sample(pos_idx, replace = TRUE),
  sample(neg_idx, replace = TRUE)
))

run_boot <- function(idx_mat, with_or = FALSE, tag = "unstratified") {
  B <- ncol(idx_mat)
  rho <- matrix(NA_real_, B, length(need), dimnames = list(NULL, need))
  auc <- matrix(NA_real_, B, length(need), dimnames = list(NULL, need))
  ors <- if (with_or) matrix(NA_real_, B, length(need), dimnames = list(NULL, need)) else NULL
  tick <- max(1L, floor(B / 5))
  for (b in seq_len(B)) {
    if (b %% tick == 0L) message("[10] ", tag, " bootstrap ", b, "/", B)
    ii <- idx_mat[, b]
    rho[b, ] <- rho_on_idx(ii)
    auc[b, ] <- auc_on_idx(ii)
    if (with_or) {
      for (j in seq_along(need)) {
        ors[b, j] <- or_on_idx(w, ii, score_mat[, j])
      }
    }
  }
  list(rho = rho, auc = auc, or = ors)
}

boot <- run_boot(boot_idx, with_or = do_or, tag = "unstratified")
boot_s <- run_boot(boot_idx_strat, with_or = FALSE, tag = "stratified")

message("[10] jackknife for BCa (n=", n, ")")
jack_rho <- matrix(NA_real_, n, length(need), dimnames = list(NULL, need))
jack_auc <- matrix(NA_real_, n, length(need), dimnames = list(NULL, need))
for (i in seq_len(n)) {
  ii <- seq_len(n)[-i]
  jack_rho[i, ] <- rho_on_idx(ii)
  jack_auc[i, ] <- auc_on_idx(ii)
}

sum_one <- function(obs, boot_col, jack_col) {
  ci <- percentile_ci(boot_col)
  bca <- bca_ci(obs, boot_col, jack_col)
  data.frame(
    point = unname(obs),
    boot_mean = mean(boot_col, na.rm = TRUE),
    ci95_lo = ci[1], ci95_hi = ci[2],
    bca_lo = bca[1], bca_hi = bca[2],
    n_boot_finite = sum(is.finite(boot_col)),
    stringsAsFactors = FALSE
  )
}

sum_delta <- function(obs_t, obs_b, boot_t, boot_b, jack_t, jack_b, strat_t, strat_b) {
  d_obs <- obs_t - obs_b
  d_boot <- boot_t - boot_b
  d_jack <- jack_t - jack_b
  d_strat <- strat_t - strat_b
  ci <- percentile_ci(d_boot)
  bca <- bca_ci(d_obs, d_boot, d_jack)
  sci <- percentile_ci(d_strat)
  includes0 <- function(lo, hi) {
    if (!is.finite(lo) || !is.finite(hi)) return(NA)
    lo <= 0 && hi >= 0
  }
  data.frame(
    delta = unname(d_obs),
    delta_boot_mean = mean(d_boot, na.rm = TRUE),
    ci95_lo = ci[1], ci95_hi = ci[2],
    bca_lo = bca[1], bca_hi = bca[2],
    ci_includes_0 = includes0(ci[1], ci[2]),
    bca_includes_0 = includes0(bca[1], bca[2]),
    p_boot_gt0_onesided = boot_p_gt0(d_boot),
    p_boot_two_sided = boot_p_two_sided(d_boot),
    prop_boot_gt0 = mean(d_boot > 0, na.rm = TRUE),
    strat_ci95_lo = sci[1], strat_ci95_hi = sci[2],
    strat_ci_includes_0 = includes0(sci[1], sci[2]),
    n_boot_finite = sum(is.finite(d_boot)),
    stringsAsFactors = FALSE
  )
}

level_rows <- lapply(need, function(m) {
  cbind(
    method = m,
    label = unname(method_lab[[m]]),
    n_clones = n,
    n_meg_positive = n_pos,
    n_boot = n_boot,
    statistic = "spearman_rho",
    sum_one(obs_rho[[m]], boot$rho[, m], jack_rho[, m]),
    stringsAsFactors = FALSE
  )
})
auc_rows <- lapply(need, function(m) {
  cbind(
    method = m,
    label = unname(method_lab[[m]]),
    n_clones = n,
    n_meg_positive = n_pos,
    n_boot = n_boot,
    statistic = "auroc",
    sum_one(obs_auc[[m]], boot$auc[, m], jack_auc[, m]),
    stringsAsFactors = FALSE
  )
})
levels_tab <- rbind(do.call(rbind, level_rows), do.call(rbind, auc_rows))
if (do_or) {
  or_rows <- lapply(need, function(m) {
    cbind(
      method = m,
      label = unname(method_lab[[m]]),
      n_clones = n,
      n_meg_positive = n_pos,
      n_boot = n_boot,
      statistic = "or_per_sd",
      sum_one(obs_or[[m]], boot$or[, m], rep(NA_real_, n)),
      stringsAsFactors = FALSE
    )
  })
  levels_tab <- rbind(levels_tab, do.call(rbind, or_rows))
}

delta_rows <- lapply(names(focal), function(m) {
  st <- steiger_spearman(w$tips_meg, w[[m]], y)
  d_rho <- sum_delta(
    obs_rho[["tips_meg"]], obs_rho[[m]],
    boot$rho[, "tips_meg"], boot$rho[, m],
    jack_rho[, "tips_meg"], jack_rho[, m],
    boot_s$rho[, "tips_meg"], boot_s$rho[, m]
  )
  d_auc <- sum_delta(
    obs_auc[["tips_meg"]], obs_auc[[m]],
    boot$auc[, "tips_meg"], boot$auc[, m],
    jack_auc[, "tips_meg"], jack_auc[, m],
    boot_s$auc[, "tips_meg"], boot_s$auc[, m]
  )
  rho_row <- data.frame(
    baseline = m,
    label = unname(focal[[m]]),
    n_clones = n,
    n_meg_positive = n_pos,
    n_boot = n_boot,
    statistic = "spearman_rho",
    tips = unname(obs_rho[["tips_meg"]]),
    baseline_point = unname(obs_rho[[m]]),
    d_rho,
    score_spearman = unname(st$r_xz),
    steiger_t = unname(st$t),
    steiger_df = unname(st$df),
    steiger_p = unname(st$p),
    stringsAsFactors = FALSE
  )
  auc_row <- data.frame(
    baseline = m,
    label = unname(focal[[m]]),
    n_clones = n,
    n_meg_positive = n_pos,
    n_boot = n_boot,
    statistic = "auroc",
    tips = unname(obs_auc[["tips_meg"]]),
    baseline_point = unname(obs_auc[[m]]),
    d_auc,
    score_spearman = NA_real_,
    steiger_t = NA_real_,
    steiger_df = NA_real_,
    steiger_p = NA_real_,
    stringsAsFactors = FALSE
  )
  out <- rbind(rho_row, auc_row)
  if (do_or) {
    d_or <- sum_delta(
      obs_or[["tips_meg"]], obs_or[[m]],
      boot$or[, "tips_meg"], boot$or[, m],
      rep(NA_real_, n), rep(NA_real_, n),
      rep(NA_real_, n_boot), rep(NA_real_, n_boot)
    )
    or_row <- data.frame(
      baseline = m,
      label = unname(focal[[m]]),
      n_clones = n,
      n_meg_positive = n_pos,
      n_boot = n_boot,
      statistic = "or_per_sd",
      tips = unname(obs_or[["tips_meg"]]),
      baseline_point = unname(obs_or[[m]]),
      d_or,
      score_spearman = NA_real_,
      steiger_t = NA_real_,
      steiger_df = NA_real_,
      steiger_p = NA_real_,
      stringsAsFactors = FALSE
    )
    out <- rbind(out, or_row)
  }
  out
})
delta_tab <- do.call(rbind, delta_rows)
rownames(delta_tab) <- NULL
delta_tab$label <- factor(delta_tab$label, levels = unname(focal))

## Headline CTS comparison drives the manuscript claim language.
cts_rho <- delta_tab[delta_tab$baseline == "cts_c11" & delta_tab$statistic == "spearman_rho", ][1, ]
claim_ok <- isTRUE(!cts_rho$ci_includes_0) && isTRUE(!cts_rho$bca_includes_0)
verdict <- data.frame(
  question = "Is TIPS Spearman meaningfully larger than fixed C11 CTS?",
  tips_rho = cts_rho$tips,
  cts_rho = cts_rho$baseline_point,
  delta_rho = cts_rho$delta,
  percentile_ci95_lo = cts_rho$ci95_lo,
  percentile_ci95_hi = cts_rho$ci95_hi,
  bca_ci95_lo = cts_rho$bca_lo,
  bca_ci95_hi = cts_rho$bca_hi,
  percentile_includes_0 = cts_rho$ci_includes_0,
  bca_includes_0 = cts_rho$bca_includes_0,
  p_boot_two_sided = cts_rho$p_boot_two_sided,
  steiger_p = cts_rho$steiger_p,
  claim_predictive_superiority = ifelse(claim_ok, "supported", "not_supported"),
  recommended_selling_point = ifelse(
    claim_ok,
    "Predictive gain vs CTS is supported by the two-sided clone-bootstrap CI; still report branch-specific resolution.",
    "Do not claim predictive superiority vs CTS. Keep the selling point on branch-specific regulatory resolution and interpretability."
  ),
  n_clones = n,
  n_meg_positive = n_pos,
  n_boot = n_boot,
  seed = seed_delta,
  note = paste(
    "Paired clone bootstrap of OOF scores; gene lists are not refit.",
    "Claim uses the two-sided 95% CI (includes 0 or not).",
    "A one-sided superiority test was not prespecified; p_boot_gt0_onesided in the delta table is diagnostic only."
  ),
  stringsAsFactors = FALSE
)

write_tsv(levels_tab, file.path(cfg$results_dir, "tables", "10_method_bootstrap_ci.tsv"))
write_tsv(delta_tab, file.path(cfg$results_dir, "tables", "10_tips_minus_baseline_bootstrap.tsv"))
write_tsv(verdict, file.path(cfg$results_dir, "tables", "10_manuscript_verdict.tsv"))

saveRDS(
  list(
    wide = w, n_boot = n_boot, seed = seed_delta,
    boot_idx = boot_idx, boot_idx_strat = boot_idx_strat,
    boot = boot, boot_strat = boot_s,
    jack_rho = jack_rho, jack_auc = jack_auc,
    obs_rho = obs_rho, obs_auc = obs_auc, obs_or = obs_or,
    levels = levels_tab, delta = delta_tab, verdict = verdict
  ),
  file.path(cfg$results_dir, "rds", "10_clone_bootstrap.rds")
)

}  # end !fig_only bootstrap

## ---- figures ----
figdir <- file.path(cfg$results_dir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
has_gg <- requireNamespace("ggplot2", quietly = TRUE)
has_pw <- requireNamespace("patchwork", quietly = TRUE)
if (!has_gg) {
  message("[10] ggplot2 not installed; tables were written, figures skipped")
} else {
  suppressPackageStartupMessages(library(ggplot2))
  col_tips <- "#C00000"
  col_base <- "grey55"
  col_sig <- "#1F4E79"
  col_ns <- "grey55"
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
      message("[10] ", basename(file), " locked; wrote ", basename(alt))
    } else {
      file.rename(tmp, file)
    }
  }

  rho_lev <- levels_tab[levels_tab$statistic == "spearman_rho", , drop = FALSE]
  rho_lev$label <- factor(rho_lev$label, levels = rev(unname(method_lab[method_order])))
  rho_lev$is_tips <- rho_lev$method == "tips_meg"
  rho_d <- delta_tab[delta_tab$statistic == "spearman_rho", , drop = FALSE]
  rho_d$label <- factor(rho_d$label, levels = rev(unname(focal)))
  rho_d$sig <- ifelse(rho_d$ci_includes_0, "CI includes 0", "CI excludes 0")
  rho_d$sig <- factor(rho_d$sig, levels = c("CI excludes 0", "CI includes 0"))
  rho_d$ann <- sprintf(
    "Δ=%.3f  [%0.3f, %0.3f]%s",
    rho_d$delta, rho_d$ci95_lo, rho_d$ci95_hi,
    ifelse(rho_d$ci_includes_0, "  n.s.", "*")
  )

  p_levels <- ggplot(rho_lev, aes(x = point, y = label)) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi, colour = is_tips),
                  width = 0.18, linewidth = 0.65, orientation = "y") +
    geom_point(aes(colour = is_tips, size = is_tips)) +
    geom_text(aes(label = sprintf("%.3f", point)), hjust = -0.25, size = 3, colour = "grey20") +
    scale_colour_manual(values = c("FALSE" = col_base, "TRUE" = col_tips), guide = "none") +
    scale_size_manual(values = c("FALSE" = 2.4, "TRUE" = 3.4), guide = "none") +
    coord_cartesian(xlim = c(
      min(0, min(rho_lev$ci95_lo, na.rm = TRUE) - 0.02),
      max(rho_lev$ci95_hi, na.rm = TRUE) + 0.12
    ), clip = "off") +
    labs(
      title = "Clone-bootstrap Spearman vs later Meg fraction",
      subtitle = sprintf("n=%d clones (%d Meg-positive); %d paired resamples", n, n_pos, n_boot),
      x = "Spearman rho (percentile 95% CI)",
      y = NULL
    ) +
    theme_held() +
    theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank())

  p_delta <- ggplot(rho_d, aes(x = delta, y = label, colour = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.45) +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi),
                  width = 0.18, linewidth = 0.7, orientation = "y") +
    geom_point(size = 2.8) +
    geom_text(aes(x = pmax(ci95_hi, 0), y = label, label = ann),
              hjust = -0.05, size = 2.6, colour = "grey20", inherit.aes = FALSE) +
    scale_colour_manual(values = c("CI excludes 0" = col_sig, "CI includes 0" = col_ns), guide = "none") +
    coord_cartesian(xlim = c(
      min(rho_d$ci95_lo, na.rm = TRUE) - 0.02,
      max(rho_d$ci95_hi, na.rm = TRUE) + 0.22
    ), clip = "off") +
    labs(
      title = "TIPS minus baseline (paired Δρ)",
      subtitle = "Same clones resampled for TIPS and the baseline",
      x = "Δ Spearman rho (percentile 95% CI)",
      y = NULL
    ) +
    theme_held() +
    theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank())

  d_cts <- boot$rho[, "tips_meg"] - boot$rho[, "cts_c11"]
  dens <- data.frame(delta = d_cts[is.finite(d_cts)])
  cts_sub <- sprintf(
    "Observed Δρ = %.3f (0.42 vs 0.33 in the text).\nTwo-sided 95%% CI [%.3f, %.3f]%s; BCa [%.3f, %.3f]%s.\nA one-sided superiority test was not prespecified; inference is the CI.",
    cts_rho$delta, cts_rho$ci95_lo, cts_rho$ci95_hi,
    ifelse(cts_rho$ci_includes_0, " includes 0", " excludes 0"),
    cts_rho$bca_lo, cts_rho$bca_hi,
    ifelse(cts_rho$bca_includes_0, " includes 0", " excludes 0")
  )
  p_dens <- ggplot(dens, aes(x = delta)) +
    geom_histogram(bins = 40, fill = "grey85", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = cts_rho$delta, colour = col_tips, linewidth = 0.6) +
    geom_vline(xintercept = c(cts_rho$ci95_lo, cts_rho$ci95_hi),
               colour = col_tips, linetype = "dotted", linewidth = 0.45) +
    labs(
      title = "Headline: TIPS − fixed C11 CTS",
      subtitle = cts_sub,
      x = "Bootstrap Δ Spearman rho",
      y = "Count",
      caption = verdict$recommended_selling_point[1]
    ) +
    theme_held()

  if (has_pw) {
    combined <- p_levels + p_delta + p_dens +
      patchwork::plot_layout(widths = c(1, 1.15, 1.15)) +
      patchwork::plot_annotation(
        caption = paste(
          "Out-of-fold scores are fixed; uncertainty is clone-sampling variability.",
          "Primary endpoint: Spearman vs Meg / all recovered day-6 progeny."
        ),
        theme = theme(plot.caption = element_text(size = 7.5, colour = col_anno, hjust = 0))
      )
    out_fig <- file.path(figdir, "10_clone_bootstrap_delta_rho.pdf")
    save_pdf(out_fig, combined, width = 13.2, height = 4.8)
    message("[10] figure -> ", out_fig)
  } else {
    save_pdf(file.path(figdir, "10a_method_rho_ci.pdf"), p_levels, 5.2, 4.2)
    save_pdf(file.path(figdir, "10b_delta_rho_forest.pdf"), p_delta, 5.6, 4.2)
    save_pdf(file.path(figdir, "10c_tips_minus_cts_density.pdf"), p_dens, 5.6, 4.4)
  }

  auc_d <- delta_tab[delta_tab$statistic == "auroc", , drop = FALSE]
  auc_d$label <- factor(auc_d$label, levels = rev(unname(focal)))
  auc_d$sig <- ifelse(auc_d$ci_includes_0, "CI includes 0", "CI excludes 0")
  auc_d$ann <- sprintf("Δ=%.3f  [%0.3f, %0.3f]", auc_d$delta, auc_d$ci95_lo, auc_d$ci95_hi)
  p_auc <- ggplot(auc_d, aes(x = delta, y = label, colour = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.45) +
    geom_errorbar(aes(xmin = ci95_lo, xmax = ci95_hi),
                  width = 0.18, linewidth = 0.7, orientation = "y") +
    geom_point(size = 2.8) +
    geom_text(aes(x = pmax(ci95_hi, 0), y = label, label = ann),
              hjust = -0.05, size = 2.6, colour = "grey20", inherit.aes = FALSE) +
    scale_colour_manual(
      values = c("CI excludes 0" = col_sig, "CI includes 0" = col_ns),
      guide = "none"
    ) +
    coord_cartesian(xlim = c(
      min(auc_d$ci95_lo, na.rm = TRUE) - 0.02,
      max(auc_d$ci95_hi, na.rm = TRUE) + 0.18
    ), clip = "off") +
    labs(
      title = "TIPS minus baseline (paired Δ AUROC)",
      subtitle = sprintf("Meg-positive = n_meg_d6 >= 1; n=%d clones", n),
      x = "Δ AUROC (percentile 95% CI)",
      y = NULL,
      caption = "Secondary to Spearman. Same paired clone resamples as the ρ figure."
    ) +
    theme_held() +
    theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank())
  save_pdf(file.path(figdir, "10_clone_bootstrap_delta_auroc.pdf"), p_auc, 6.4, 4.0)
}

fmt_ci <- function(lo, hi) sprintf("[%0.3f, %0.3f]", lo, hi)
message("[10] verdict: ", verdict$claim_predictive_superiority)
message("[10] TIPS ρ=", signif(cts_rho$tips, 3),
        "  CTS ρ=", signif(cts_rho$baseline_point, 3),
        "  Δρ=", signif(cts_rho$delta, 3),
        "  percentile CI ", fmt_ci(cts_rho$ci95_lo, cts_rho$ci95_hi),
        if (isTRUE(cts_rho$ci_includes_0)) " includes 0" else " excludes 0",
        "  BCa ", fmt_ci(cts_rho$bca_lo, cts_rho$bca_hi),
        "  two-sided boot P=", signif(cts_rho$p_boot_two_sided, 3),
        "  Steiger P=", signif(cts_rho$steiger_p, 3),
        "  (one-sided P not used for inference)")
for (m in names(focal)) {
  r <- delta_tab[delta_tab$baseline == m & delta_tab$statistic == "spearman_rho", ][1, ]
  message("[10] TIPS − ", r$label, ": Δρ=", signif(r$delta, 3),
          " CI ", fmt_ci(r$ci95_lo, r$ci95_hi),
          if (isTRUE(r$ci_includes_0)) " includes 0" else " excludes 0")
}
message("[10] ", verdict$recommended_selling_point)
