## 08_summary_figure.R — held-out summary figures
##
## One-row three-panel design/reproducibility/overlap figure:
##   a. Clone eligibility funnel (replaces the raw 07a count bar)
##   b. Gene-by-fold TIPS presence matrix
##   c. Three-row forest: held-out TIPS vs Weinreb, MuTrans, metacell TIPS
##      (Fig. 2c already reports metacell/MuTrans vs Weinreb and the gene matrix)
##
## Also still writes the two-panel Meg prediction figure (score vs fraction
## and baseline Spearman) as 08_summary_meg_heldout.pdf.

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = here::here("examples", "hematopoietic_LARRY", "single_cell", "9_Held_out_clone_fate_prediction_cursor")
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_ensure()
figdir <- file.path(cfg$results_dir, "figures")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("08_summary_figure.R requires ggplot2")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("08_summary_figure.R requires patchwork")
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

col_tips <- "#C00000"
col_base <- "grey55"
col_anno <- "grey45"
col_bar <- "#4A86C8"

theme_held <- function() {
  theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10.5, face = "bold"),
      plot.subtitle = element_text(size = 8, lineheight = 1.15),
      plot.caption = element_text(size = 7.5, colour = col_anno, hjust = 0)
    )
}

fmt_p <- function(p, digits = 1) {
  if (!is.finite(p) || p <= 0) return(as.character(p))
  if (p >= 0.001) return(formatC(p, format = "f", digits = 3))
  e <- floor(log10(p))
  m <- p / (10^e)
  sprintf("%.*fe%d", digits, m, e)
}

fmt_rho <- function(x) sprintf("%.3f", x)
fmt_n <- function(n) format(as.integer(n), big.mark = ",", scientific = FALSE, trim = TRUE)

split_genes <- function(x) {
  norm_genes(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
}

cairo_ok <- isTRUE(capabilities("cairo"))
save_pdf <- function(file, plot, width, height) {
  if (cairo_ok) {
    ggsave(file, plot, width = width, height = height, bg = "white",
           device = grDevices::cairo_pdf)
  } else {
    ggsave(file, plot, width = width, height = height, bg = "white")
  }
}

counts_f <- file.path(cfg$results_dir, "01_feasibility_counts.tsv")
mod_f <- file.path(cfg$results_dir, "tables", "04_fold_module_genes.tsv")
excl_f <- file.path(cfg$results_dir, "tables", "03b_exclusion_audit.tsv")
pred_f <- file.path(cfg$results_dir, "tables", "04_oof_clone_predictions.tsv")
perf_f <- file.path(cfg$results_dir, "tables", "05_method_performance.tsv")
sets_rds <- file.path(cfg$results_dir, "rds", "03_gene_sets.rds")
if (!file.exists(counts_f)) stop("Run 01 first: missing ", counts_f)
if (!file.exists(mod_f)) stop("Run 04 first: missing ", mod_f)
if (!file.exists(pred_f)) stop("Run 04 first: missing ", pred_f)
if (!file.exists(perf_f)) stop("Run 05 first: missing ", perf_f)
if (!file.exists(sets_rds)) stop("Run 03 first: missing ", sets_rds)

ct <- utils::read.delim(counts_f, stringsAsFactors = FALSE)
n_of <- function(item) {
  hit <- ct$n[ct$item == item]
  if (!length(hit)) stop("Missing feasibility item: ", item)
  as.integer(hit[1])
}
n_c11 <- n_of("C11_day2_cells")
n_assigned <- n_of("C11_day2_cells_with_unique_clone")
n_clones <- n_of("distinct_C11_clones")
n_eval <- n_of("C11_clones_with_annotated_day6_progeny")
n_meg <- n_of("C11_clones_with_at_least_one_Meg")
n_baso <- n_of("C11_clones_with_at_least_one_Baso")

excl_line <- "All test-clone cells were removed before TIPS refitting."
if (file.exists(excl_f)) {
  ex <- utils::read.delim(excl_f, stringsAsFactors = FALSE)
  rem <- as.integer(ex$n_test_clone_cells_remaining_in_training_object)
  n_folds_excl <- nrow(ex)
  if (length(rem) && all(rem == 0L)) {
    excl_line <- sprintf(
      "All test-clone cells were removed before TIPS refitting (0 remaining in %d/%d training objects).",
      n_folds_excl, n_folds_excl
    )
  } else {
    excl_line <- sprintf(
      "Test-clone cells remaining in training objects: %s.",
      paste(rem, collapse = ",")
    )
  }
}

funnel <- data.frame(
  stage = factor(
    c("C11 day-2 cells", "Uniquely assigned cells", "Distinct C11 clones", "Evaluable clones"),
    levels = rev(c("C11 day-2 cells", "Uniquely assigned cells", "Distinct C11 clones", "Evaluable clones"))
  ),
  n = c(n_c11, n_assigned, n_clones, n_eval),
  stringsAsFactors = FALSE
)
funnel$lab <- sprintf("%s  %s", fmt_n(funnel$n), as.character(funnel$stage))
flow_line <- sprintf(
  "%s C11 cells -> %s uniquely assigned cells -> %s clones -> %s evaluable clones -> %d Meg-positive/%d Baso-positive.",
  fmt_n(n_c11), fmt_n(n_assigned), fmt_n(n_clones), fmt_n(n_eval), n_meg, n_baso
)

p_funnel <- ggplot(funnel, aes(x = n, y = stage)) +
  geom_col(fill = col_bar, width = 0.72) +
  geom_text(aes(label = lab), hjust = 0, nudge_x = max(funnel$n) * 0.03, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.55)), labels = fmt_n) +
  labs(
    title = "Clone eligibility and held-out design",
    subtitle = flow_line,
    caption = excl_line,
    x = "Count",
    y = NULL
  ) +
  theme_held() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

md <- utils::read.delim(mod_f, stringsAsFactors = FALSE)
arm_rows <- md[md$set %in% c("tips_meg", "tips_baso"), , drop = FALSE]
if (!nrow(arm_rows)) stop("No tips_meg/tips_baso rows in ", mod_f)
presence <- do.call(rbind, lapply(seq_len(nrow(arm_rows)), function(i) {
  gs <- split_genes(arm_rows$genes[i])
  if (!length(gs)) {
    return(data.frame(gene = character(), fold = integer(), arm = character(), present = logical()))
  }
  data.frame(
    gene = gs,
    fold = as.integer(arm_rows$fold[i]),
    arm = ifelse(arm_rows$set[i] == "tips_meg", "Meg TIPS", "Baso TIPS"),
    present = TRUE,
    stringsAsFactors = FALSE
  )
}))
folds <- sort(unique(as.integer(arm_rows$fold)))
grid <- do.call(rbind, lapply(split(presence, presence$arm), function(dd) {
  genes <- unique(dd$gene)
  g <- expand.grid(gene = genes, fold = folds, stringsAsFactors = FALSE)
  g$arm <- dd$arm[1]
  key <- paste(dd$gene, dd$fold, sep = "\t")
  g$present <- paste(g$gene, g$fold, sep = "\t") %in% key
  g
}))
grid$arm <- factor(grid$arm, levels = c("Meg TIPS", "Baso TIPS"))
grid$gene <- factor(grid$gene, levels = rev(unique(c(
  sort(unique(grid$gene[grid$arm == "Meg TIPS"])),
  sort(unique(grid$gene[grid$arm == "Baso TIPS"]))
))))
n_fold_ok <- length(unique(arm_rows$fold[arm_rows$set == "tips_meg" & arm_rows$n > 0]))
n_fold_tot <- length(folds)

pf <- utils::read.delim(perf_f, stringsAsFactors = FALSE)
if ("denominator" %in% names(pf)) {
  pf <- pf[pf$denominator == "all_day6_progeny" | is.na(pf$denominator), , drop = FALSE]
}
meg_row <- pf[pf$method == "tips_meg" & pf$endpoint == "Meg", , drop = FALSE]
if (!nrow(meg_row)) stop("No tips_meg Meg row in ", perf_f)
meg_row <- meg_row[1, ]
n_scored <- as.integer(meg_row$n_scored)
n_eligible <- as.integer(meg_row$n_eligible)
fold_note <- sprintf(
  "Non-empty Meg module in %d/%d folds; scored %d/%d clones.",
  n_fold_ok, n_fold_tot, n_scored, n_eligible
)

p_fold <- ggplot(grid, aes(x = factor(fold), y = gene, fill = present)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(values = c(`TRUE` = col_tips, `FALSE` = "grey90"), guide = "none") +
  facet_grid(arm ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Fold-specific module recovery",
    subtitle = fold_note,
    x = "Held-out fold",
    y = NULL
  ) +
  theme_held() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(size = 8, face = "bold")
  )

sets <- readRDS(sets_rds)
heldout_tips <- unique(presence$gene[presence$arm == "Meg TIPS" & presence$present])
weinreb <- norm_genes(sets$weinreb_meg)
mutrans <- norm_genes(sets$mutrans_td_meg)
metacell <- norm_genes(sets$metacell_tips_meg)
hvg <- norm_genes(sets$hvg)
universe <- unique(c(hvg, weinreb, norm_genes(sets$weinreb_baso)))

fisher_overlap <- function(a, b, U) {
  A <- intersect(a, U)
  B <- intersect(b, U)
  N <- length(U)
  k <- length(intersect(A, B))
  na <- length(A)
  nb <- length(B)
  tab <- matrix(c(k, na - k, nb - k, N - na - nb + k), nrow = 2, byrow = TRUE)
  ft <- stats::fisher.test(tab, alternative = "two.sided")
  data.frame(
    n_a = na, n_b = nb, overlap_k = k, universe_n = N,
    p_value = ft$p.value, odds_ratio = unname(ft$estimate),
    ci_low = unname(ft$conf.int[1]), ci_high = unname(ft$conf.int[2]),
    overlap_genes = paste(sort(intersect(A, B)), collapse = ","),
    stringsAsFactors = FALSE
  )
}

comparisons <- list(
  list(label = "vs Weinreb fate genes", a = heldout_tips, b = weinreb),
  list(label = "vs MuTrans drivers",    a = heldout_tips, b = mutrans),
  list(label = "vs metacell TIPS",      a = heldout_tips, b = metacell)
)
fet <- do.call(rbind, lapply(comparisons, function(cmp) {
  cbind(fisher_overlap(cmp$a, cmp$b, universe), label = cmp$label)
}))
fet$label <- factor(fet$label, levels = rev(vapply(comparisons, `[[`, "", "label")))
write_tsv(fet, file.path(cfg$results_dir, "tables", "08_heldout_overlap_FET.tsv"))

clip_or <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x) & x > 0] <- 2000
  x[!is.finite(x)] <- 0.05
  pmin(pmax(x, 0.05), 2000)
}
fet$or_plot <- clip_or(fet$odds_ratio)
fet$ci_low_plot <- clip_or(fet$ci_low)
fet$ci_high_plot <- clip_or(fet$ci_high)
fet$lab <- sprintf("k=%d, P=%s", fet$overlap_k, vapply(fet$p_value, fmt_p, character(1)))

p_forest <- ggplot(fet, aes(x = or_plot, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low_plot, xmax = ci_high_plot), width = 0.18, linewidth = 0.65, colour = col_tips, orientation = "y") +
  geom_point(size = 2.8, colour = col_tips) +
  geom_text(aes(x = ci_high_plot, label = lab), hjust = -0.08, size = 2.7, colour = "grey20") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.72))) +
  labs(
    title = "Held-out TIPS overlap",
    x = "Odds ratio (log scale, 95% CI)",
    y = NULL
  ) +
  theme_held() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank()
  )

three <- p_funnel + p_fold + p_forest +
  plot_layout(nrow = 1, widths = c(1.15, 0.85, 1.15)) +
  plot_annotation(
    caption = paste(
      "c: held-out TIPS vs Weinreb fate genes, MuTrans drivers, and metacell TIPS.",
      "Metacell/MuTrans vs Weinreb and gene-level overlap are in Fig. 2c."
    ),
    theme = theme(plot.caption = element_text(size = 7.5, colour = col_anno, hjust = 0))
  )
out_three <- file.path(figdir, "08_heldout_design_three_panel.pdf")
save_pdf(out_three, three, width = 12.4, height = 3.9)
message("[08] ", flow_line)
message("[08] ", fold_note)
message("[08] ", excl_line)
message("[08] three-panel -> ", out_three)

## ---- existing two-panel Meg prediction figure (Figure 2d) ----
pr <- utils::read.delim(pred_f, stringsAsFactors = FALSE)
tips <- pr[pr$method == "tips_meg", , drop = FALSE]
if (!nrow(tips)) stop("No tips_meg rows in ", pred_f)
baso_row <- pf[pf$method == "tips_baso" & pf$endpoint == "Baso", , drop = FALSE]
if (!nrow(baso_row)) stop("No tips_baso Baso row in ", perf_f)
baso_row <- baso_row[1, ]
stats_line <- sprintf(
  "n=%d clones; %d Meg-positive; rho=%s; quasibinomial P=%s; AUROC=%.3f.",
  as.integer(meg_row$n_scored), as.integer(meg_row$n_positive),
  fmt_rho(meg_row$spearman_rho), fmt_p(meg_row$or_p), meg_row$auroc
)
baso_line <- sprintf("Baso arm: rho=%s, AUROC=%.3f.", fmt_rho(baso_row$spearman_rho), baso_row$auroc)
p_note <- paste(
  "P was obtained from a quasibinomial model adjusted for",
  "library and starting population."
)

p_left <- ggplot(tips, aes(x = score, y = frac_meg_d6)) +
  geom_point(aes(size = n_d6), alpha = 0.7, colour = col_tips, stroke = 0) +
  scale_size_area(name = "Day-6 n", max_size = 6) +
  labs(
    title = "Clone-held-out TIPS vs later Meg fraction",
    subtitle = stats_line,
    x = "Out-of-fold TIPS score",
    y = "Later Meg fraction"
  ) +
  theme_held() +
  theme(legend.position = "bottom", legend.margin = margin(t = -4))

method_map <- c(
  tips_meg = "TIPS",
  cts_c11 = "C11 CTS",
  ppi_unweighted_meg = "Unweighted PPI",
  meg_markers_train = "Training Meg markers",
  node_de_megup = "Node-level DE",
  mutrans_td_meg = "MuTrans TD"
)
dot <- pf[pf$endpoint == "Meg" & pf$method %in% names(method_map), , drop = FALSE]
dot$label <- unname(method_map[dot$method])
dot$label <- factor(dot$label, levels = rev(unname(method_map)))
dot$is_tips <- dot$method == "tips_meg"
dot$rho_lab <- fmt_rho(dot$spearman_rho)
x_max <- max(dot$spearman_rho, na.rm = TRUE)

p_right <- ggplot(dot, aes(x = spearman_rho, y = label)) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = spearman_rho, y = label, yend = label, colour = is_tips), linewidth = 0.6) +
  geom_point(aes(colour = is_tips, size = is_tips)) +
  geom_text(aes(label = rho_lab), hjust = -0.25, size = 3, colour = "grey20") +
  scale_colour_manual(values = c("FALSE" = col_base, "TRUE" = col_tips), guide = "none") +
  scale_size_manual(values = c("FALSE" = 2.6, "TRUE" = 3.6), guide = "none") +
  coord_cartesian(
    xlim = c(min(0, min(dot$spearman_rho, na.rm = TRUE)), x_max + 0.12),
    clip = "off"
  ) +
  labs(
    title = "Comparison with baselines",
    x = "Spearman rho",
    y = NULL,
    caption = baso_line
  ) +
  theme_held() +
  theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank())

pred_fig <- p_left + p_right +
  plot_layout(widths = c(1.2, 1), guides = "keep") +
  plot_annotation(
    caption = p_note,
    theme = theme(plot.caption = element_text(size = 8, colour = col_anno, hjust = 0))
  )
out_pred <- file.path(figdir, "08_summary_meg_heldout.pdf")
save_pdf(out_pred, pred_fig, width = 10.2, height = 4.6)
message("[08] ", stats_line)
message("[08] prediction panel -> ", out_pred)
