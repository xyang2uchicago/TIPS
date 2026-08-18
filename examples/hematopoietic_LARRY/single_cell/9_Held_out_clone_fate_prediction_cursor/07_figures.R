## 07_figures.R — main LARRY held-out figure panels a–g

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor"
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_ensure()
figdir <- file.path(cfg$results_dir, "figures")

has_gg <- requireNamespace("ggplot2", quietly = TRUE)
if (has_gg) suppressPackageStartupMessages(library(ggplot2))

theme_held <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

## a: design counts
counts_f <- file.path(cfg$results_dir, "01_feasibility_counts.tsv")
if (file.exists(counts_f) && has_gg) {
  ct <- utils::read.delim(counts_f, stringsAsFactors = FALSE)
  ct$item <- factor(ct$item, levels = rev(ct$item))
  p <- ggplot(ct, aes(x = n, y = item)) +
    geom_col(fill = "#4A86C8") +
    geom_text(aes(label = n), hjust = -0.1, size = 3) +
    labs(title = "a. LARRY clone-level held-out design", x = "Count", y = NULL) +
    expand_limits(x = max(ct$n) * 1.15) + theme_held()
  ggsave(file.path(figdir, "07a_design_counts.pdf"), p, width = 8, height = 5, bg = "white")
}

## b: TIPS module genes across folds
mod_f <- file.path(cfg$results_dir, "tables", "04_fold_module_genes.tsv")
if (file.exists(mod_f) && has_gg) {
  md <- utils::read.delim(mod_f, stringsAsFactors = FALSE)
  tips <- md[md$set == "tips_meg", ]
  p <- ggplot(tips, aes(x = factor(fold), y = n)) +
    geom_col(fill = "#C0504D") +
    labs(title = "b. Held-out TIPS Meg module size per fold", x = "Fold", y = "Genes") +
    theme_held()
  ggsave(file.path(figdir, "07b_module_size.pdf"), p, width = 5, height = 4, bg = "white")
}

## c: day-2 TIPS vs later Meg fraction
pred_f <- file.path(cfg$results_dir, "tables", "04_oof_clone_predictions.tsv")
if (file.exists(pred_f) && has_gg) {
  pr <- utils::read.delim(pred_f, stringsAsFactors = FALSE)
  tips <- pr[pr$method == "tips_meg", ]
  p <- ggplot(tips, aes(x = score, y = frac_meg_d6, size = n_d6)) +
    geom_point(alpha = 0.7, colour = "#C0504D") +
    geom_smooth(method = "lm", se = TRUE, colour = "grey20", linewidth = 0.6, show.legend = FALSE) +
    labs(
      title = "c. Out-of-fold day-2 TIPS score vs later Meg fraction",
      subtitle = "Point size = day-6 progeny; barcode-matched clonal relatives, not literal descendants",
      x = "Day-2 TIPS Meg score (OOF)", y = "Day-6 Meg fraction"
    ) +
    theme_held()
  ggsave(file.path(figdir, "07c_score_vs_meg_fraction.pdf"), p, width = 6.5, height = 5, bg = "white")
}

## d: method comparison
perf_f <- file.path(cfg$results_dir, "tables", "05_method_performance.tsv")
if (file.exists(perf_f) && has_gg) {
  pf <- utils::read.delim(perf_f, stringsAsFactors = FALSE)
  if ("denominator" %in% names(pf)) pf <- pf[pf$denominator == "all_day6_progeny" | is.na(pf$denominator), ]
  meg <- pf[pf$endpoint == "Meg", ]
  meg$method <- factor(meg$method, levels = rev(meg$method[order(meg$spearman_rho)]))
  p <- ggplot(meg, aes(x = spearman_rho, y = method)) +
    geom_col(fill = "#4A86C8") +
    labs(title = "d. Out-of-fold Spearman vs later Meg fraction", x = "Spearman rho", y = NULL) +
    theme_held()
  ggsave(file.path(figdir, "07d_baselines_spearman.pdf"), p, width = 7, height = 5, bg = "white")
}

## e: two wells
well_f <- file.path(cfg$results_dir, "tables", "05_two_well_replication.tsv")
if (file.exists(pred_f) && has_gg) {
  pr <- utils::read.delim(pred_f, stringsAsFactors = FALSE)
  tips <- pr[pr$method == "tips_meg" & pr$in_both_later_wells %in% c(TRUE, "TRUE", "true"), ]
  if (nrow(tips)) {
    p1 <- ggplot(tips, aes(x = score, y = frac_meg_d6_well1, size = n_d6_well1)) +
      geom_point(alpha = 0.7, colour = "#548235") +
      labs(title = "e. Well 1", x = "Day-2 TIPS", y = "Meg fraction well 1") + theme_held()
    p2 <- ggplot(tips, aes(x = score, y = frac_meg_d6_well2, size = n_d6_well2)) +
      geom_point(alpha = 0.7, colour = "#833C0C") +
      labs(title = "Well 2", x = "Day-2 TIPS", y = "Meg fraction well 2") + theme_held()
    p3 <- ggplot(tips, aes(x = frac_meg_d6_well1, y = frac_meg_d6_well2, size = n_d6)) +
      geom_point(alpha = 0.7) +
      labs(title = "Well concordance", x = "Meg fraction well 1", y = "Meg fraction well 2") + theme_held()
    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined <- patchwork::wrap_plots(p1, p2, p3, nrow = 1)
      ggsave(file.path(figdir, "07e_two_wells.pdf"), combined, width = 12, height = 4, bg = "white")
    } else {
      ggsave(file.path(figdir, "07e_well1.pdf"), p1, width = 4.5, height = 4, bg = "white")
      ggsave(file.path(figdir, "07e_well2.pdf"), p2, width = 4.5, height = 4, bg = "white")
      ggsave(file.path(figdir, "07e_concordance.pdf"), p3, width = 4.5, height = 4, bg = "white")
    }
  }
}

## f: 8-row Meg/Baso resolvability
res_f <- file.path(cfg$results_dir, "tables", "06_resolvability_with_heldout.tsv")
if (!file.exists(res_f)) res_f <- file.path(cfg$results_dir, "tables", "06_branch_resolvability.tsv")
if (!file.exists(res_f)) res_f <- file.path(cfg$results_dir, "tables", "06_internal_resolvability.tsv")
if (file.exists(res_f)) {
  rs <- utils::read.delim(res_f, stringsAsFactors = FALSE)
  out_f <- file.path(figdir, "07f_resolvability_table.pdf")
  if (has_gg && all(c("level", "Meg", "Baso") %in% names(rs))) {
    rs$level <- factor(rs$level, levels = rev(rs$level))
    long <- rbind(
      data.frame(level = rs$level, arm = "Meg", text = rs$Meg, stringsAsFactors = FALSE),
      data.frame(level = rs$level, arm = "Baso", text = rs$Baso, stringsAsFactors = FALSE)
    )
    long$text <- vapply(as.character(long$text), function(t) {
      paste(strwrap(t, width = 44), collapse = "\n")
    }, character(1))
    p <- ggplot(long, aes(x = arm, y = level)) +
      geom_tile(fill = "grey95", color = "white") +
      geom_text(aes(label = text), size = 2.4, lineheight = 0.9) +
      labs(title = "f. Matched Meg/Baso branch resolvability", x = NULL, y = NULL) +
      theme_held() + theme(panel.grid = element_blank())
    ggsave(out_f, p, width = 12, height = 8, bg = "white")
  } else {
    grDevices::pdf(out_f, width = 12, height = 8)
    grid::grid.newpage()
    show <- if (all(c("level", "Meg", "Baso") %in% names(rs))) {
      rs[, intersect(c("level", "diagnostic", "Meg", "Baso"), names(rs)), drop = FALSE]
    } else rs
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.table(show, rows = NULL)
    } else {
      plot.new()
      text(0, 1, paste(utils::capture.output(print(show)), collapse = "\n"),
           adj = c(0, 1), family = "mono", cex = 0.45)
    }
    grDevices::dev.off()
  }
}

## g: Baso audit funnel
fun_f <- file.path(cfg$results_dir, "tables", "06_baso_evidence_funnel.tsv")
if (file.exists(fun_f) && has_gg) {
  fn <- utils::read.delim(fun_f, stringsAsFactors = FALSE)
  p <- ggplot(fn, aes(x = reorder(stage, -n_genes), y = n_genes)) +
    geom_col(fill = "#7F6000") +
    labs(title = "g. Baso evidence loss (Weinreb Baso + MuTrans TD)", x = NULL, y = "Genes") +
    theme_held() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(figdir, "07g_baso_funnel.pdf"), p, width = 7, height = 4.5, bg = "white")
}

message("[07] figures -> ", figdir)
