## 8.2.2_Caesar_forestMatrix_baso_TIPScll_TIPSmeta.R — forest plot + compact overlap figure (STRING),
## Basophil lineage. Replaces the venn3 diagrams in 8.2.2_Caesar_venn3_baso_TIPScll_TIPSmeta.R with
## one figure: a forest plot of pairwise Fisher-exact odds ratios (log scale, 95% CI), precision vs
## recall against Weinreb Basophil, named tiles for genes in >=2 sources, and an UpSet of
## intersection sizes (unique genes counted; symbols only when the intersection is small).
##
## Only 3 comparisons (not 5, cf. the Megakaryocyte version): the metacell TIPS pipeline has no
## CM/dualpull output for the 4->9 (CP->Baso) arm (its dualpull table has zero linkage=="CM" rows),
## so there is no "metacell TIPS" gene set for Basophil, and no matrix column for it either. MuTrans
## TD genes are their own evidence source, not a stand-in labeled "metacell TIPS".
##
## Same 3 gene sets and `style` convention as 8.2.2_Caesar_venn3_baso_TIPScll_TIPSmeta.R; companion
## only, does not source/modify it.
##   source("F:/projects/TIPS/source/GSE140802_lineage_tracking/8.2.2_Caesar_forestMatrix_baso_TIPScll_TIPSmeta.R")

for (pkg in c("readxl", "ggplot2", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, call. = FALSE)
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

## ---- config ----
GSE140802_cell_ROOT <- Sys.getenv("GSE140802_cell_ROOT", "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory")
pub_xlsx    <- Sys.getenv("WEINREB_TABLE_S3", "F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/aaw3381-Weinreb-Table-S3.xlsx")
pub_lineage <- "Basophil"; pub_sheet <- "DGE of progenitors in vitro"
TD.cut <- 0.7
td_csv_dir <- Sys.getenv("MUTRANS_TD_CSV_DIR", here::here("examples", "hematopoietic_LARRY", "data", "larry_figures"))
CP_cluster <- 4; CM_cluster <- 9  # larry_transition_4_to_9
fisher_alt <- Sys.getenv("TIPS_FET_ALTERNATIVE", "two.sided")
tingjun_tips  <- Sys.getenv("TINGJUN_TIPS_ROOT", "F:/projects/TIPS/results/GSE140802_lineage_tracking/Tingjunl/TIPS/7_scaledata_leiden_r0_8_TIPS_STRING")
holly_out_dir <- file.path(tingjun_tips, "Holly")
norm_genes <- function(x) { x <- unique(toupper(as.character(x))); x[!is.na(x) & nzchar(x)] }

## style: 'biased' (default) = full/unfiltered TIPS node set (Weinreb's fate lists aren't fate-exclusive
## either, so an exclusive TIPS set would be apples-to-oranges). 'leanning' = old strict filter: keep a
## gene shared with `other`'s edges only if it also has a delta>0 edge on the `linkage` side.
style <- "biased"  # "biased" | "leanning"
linkage_set <- function(df, linkage, other) {
  nodes <- norm_genes(with(df[df$linkage == linkage, ], c(from, to)))
  if (style != "leanning") return(nodes)
  shared <- intersect(nodes, norm_genes(with(df[df$linkage == other, ], c(from, to))))
  pos    <- norm_genes(with(df[df$linkage == linkage & df$delta > 0, ], c(from, to)))
  setdiff(nodes, setdiff(shared, pos))
}

## ---- gene sets: Weinreb Basophil, cell-level TIPS (C11->C10 CM arm), MuTrans TD (A4->A9) ----
pub_df <- readxl::read_excel(pub_xlsx, sheet = pub_sheet)
pub_df <- pub_df[!is.na(pub_df$`Gene symbol`) & nzchar(as.character(pub_df$`Gene symbol`)), ]
baso <- norm_genes(pub_df$`Gene symbol`[pub_df$Lineage == pub_lineage])
tips_cll_df <- utils::read.delim(file.path(tingjun_tips, "results_core_11_10vs17", "cisTarget_predicted_11",
                                            "PPI_graph_GRN_prediction_CTS_11_dualpull_final_table.tsv"), stringsAsFactors = FALSE)
tips_cll <- linkage_set(tips_cll_df, "CM", "CF")
td_df <- utils::read.csv(file.path(td_csv_dir, paste0("td_genes_scores_", CP_cluster, "_to_", CM_cluster, "_seacell.csv")), stringsAsFactors = FALSE)
td_mutrans <- norm_genes(td_df$Genes[abs(td_df$corr) > TD.cut])

## ---- background universe (cell-level HVG union Weinreb Table S3) ----
hvg <- norm_genes(names(readRDS(file.path(GSE140802_cell_ROOT, "cell_cycle_signature_cor_3khvg.rds"))))
universe <- unique(c(hvg, norm_genes(pub_df$`Gene symbol`)))
message("[8.2.2 forestMatrix] universe n=", length(universe), "; Basophil n=", length(baso),
        "; TIPS_cll(C11->C10) n=", length(tips_cll), "; TD_MuTrans(A4->A9) n=", length(td_mutrans))

## ---- pairwise Fisher's exact test: odds ratio + 95% CI ----
fisher_overlap <- function(a, b, U) {
  A <- intersect(a, U); B <- intersect(b, U); N <- length(U)
  k <- length(intersect(A, B)); na <- length(A); nb <- length(B)
  tab <- matrix(c(k, na - k, nb - k, N - na - nb + k), nrow = 2, byrow = TRUE)
  ft <- stats::fisher.test(tab, alternative = fisher_alt)
  data.frame(n_a = na, n_b = nb, overlap_k = k, p_value = ft$p.value, odds_ratio = unname(ft$estimate),
             ci_low = unname(ft$conf.int[1]), ci_high = unname(ft$conf.int[2]),
             overlap_genes = paste(sort(intersect(A, B)), collapse = ", "), stringsAsFactors = FALSE)
}

## ---- the 3 available comparisons (see header); red = vs Weinreb, blue = method agreement ----
comparisons <- list(
  list(label = "Cell-level TIPS vs Weinreb fate genes", a = tips_cll, b = baso,        group = "vs_Weinreb"),
  list(label = "MuTrans drivers vs Weinreb fate genes",  a = td_mutrans, b = baso,     group = "vs_Weinreb"),
  list(label = "Cell-level TIPS vs MuTrans drivers",     a = tips_cll, b = td_mutrans, group = "method_agreement")
)
fet <- do.call(rbind, lapply(comparisons, function(cmp) cbind(fisher_overlap(cmp$a, cmp$b, universe), label = cmp$label, group = cmp$group)))
fet$label <- factor(fet$label, levels = rev(vapply(comparisons, `[[`, "", "label")))  # top-to-bottom = requested order
dir.create(holly_out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(fet, file.path(holly_out_dir, "forestMatrix_Baso_FET.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("[8.2.2 forestMatrix] ", paste(sprintf("%s: k=%d OR=%.3g [%.3g,%.3g] p=%.4g", fet$label, fet$overlap_k,
        fet$odds_ratio, fet$ci_low, fet$ci_high, fet$p_value), collapse = " | "))

## ---- forest plot: log-scale OR with 95% CI, labeled with k and P; clip only the plotted position ----
## (a real k=0 row can occur here -- OR/CI can be 0 or Inf; the reported values above are never clipped)
clip <- function(x) pmin(pmax(x, 0.05), 2000)
fet$or_plot <- clip(fet$odds_ratio); fet$ci_low_plot <- clip(fet$ci_low); fet$ci_high_plot <- clip(fet$ci_high)
fet$lab <- sprintf("k=%d, P=%.2g", fet$overlap_k, fet$p_value)

forest_p <- ggplot(fet, aes(x = or_plot, y = label, color = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = ci_low_plot, xmax = ci_high_plot), width = 0.15, linewidth = 0.7, orientation = "y") +
  geom_point(size = 3) +
  geom_text(aes(x = ci_high_plot, label = lab), hjust = -0.08, size = 3.2, color = "black", show.legend = FALSE) +
  scale_x_log10(expand = expansion(mult = c(0.08, 0.55))) +
  scale_color_manual(values = c(vs_Weinreb = "#D62728", method_agreement = "#4A86C8"),
                      labels = c(vs_Weinreb = "vs Weinreb lineage-tracing", method_agreement = "TIPS/MuTrans agreement"), name = NULL) +
  labs(x = "Odds ratio (log scale, 95% CI)", y = NULL, title = paste0(pub_lineage, ": convergence with lineage-traced fate genes")) +
  theme_minimal(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank())

## ---- compact figure: precision/recall + named shared genes + UpSet (unique genes counted) ----
evidence_sets <- list(`Weinreb fate` = baso, `Cell TIPS` = tips_cll, `MuTrans` = td_mutrans)
genes <- unique(unlist(evidence_sets))
pres <- sapply(evidence_sets, function(s) genes %in% s)
rownames(pres) <- genes
n_src <- rowSums(pres)
ord_ge2 <- rownames(pres)[n_src >= 2L]
ord_ge2 <- ord_ge2[order(-n_src[ord_ge2], ord_ge2)]
if (!length(ord_ge2)) stop("no gene has >=2 evidence sources")
message("[8.2.2 forestMatrix] genes with >=2 evidence sources (n=", length(ord_ge2), "): ",
        paste(ord_ge2, collapse = ", "))
mat <- expand.grid(gene = ord_ge2, source = names(evidence_sets), stringsAsFactors = FALSE)
mat$present <- pres[cbind(mat$gene, mat$source)]
mat$gene <- factor(mat$gene, levels = rev(ord_ge2))
mat$source <- factor(mat$source, levels = names(evidence_sets))
mat_ge2 <- ggplot(mat, aes(x = source, y = gene, fill = present)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_manual(values = c(`TRUE` = "#4A86C8", `FALSE` = "grey92"), guide = "none") +
  labs(x = NULL, y = NULL, title = "Shared genes (\u22652 sources)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(),
        plot.title = element_text(size = 11))

n_weinreb <- length(baso)
pr <- do.call(rbind, lapply(
  list(list(method = "Cell TIPS", genes = tips_cll),
       list(method = "MuTrans", genes = td_mutrans)),
  function(x) {
    g <- unique(x$genes)
    k <- length(intersect(g, baso))
    data.frame(method = x$method, n_module = length(g), k = k,
               precision = if (length(g)) k / length(g) else NA_real_,
               recall = if (n_weinreb) k / n_weinreb else NA_real_,
               stringsAsFactors = FALSE)
  }
))
pr$method <- factor(pr$method, levels = rev(pr$method))
pr_long <- rbind(
  data.frame(method = pr$method, metric = "Precision (module \u2229 Weinreb / module)",
             value = pr$precision, lab = sprintf("%d/%d", pr$k, pr$n_module), stringsAsFactors = FALSE),
  data.frame(method = pr$method, metric = sprintf("Recall (module \u2229 Weinreb / %d Weinreb)", n_weinreb),
             value = pr$recall, lab = sprintf("%d/%d", pr$k, n_weinreb), stringsAsFactors = FALSE)
)
pr_long$metric <- factor(pr_long$metric, levels = unique(pr_long$metric))
message("[8.2.2 forestMatrix] vs Weinreb Basophil: ",
        paste(sprintf("%s precision=%d/%d recall=%d/%d", pr$method, pr$k, pr$n_module, pr$k, n_weinreb),
              collapse = " | "))

pr_p <- ggplot(pr_long, aes(x = value, y = method, fill = metric)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68) +
  geom_text(aes(label = lab), position = position_dodge(width = 0.72),
            hjust = -0.08, size = 3.1, color = "grey15") +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     limits = c(0, 1), expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c("#4A86C8", "#D62728"), name = NULL) +
  labs(x = NULL, y = NULL, title = "Precision vs recall of the Weinreb fate list") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11))

n_label_max <- 8L
set_cols <- names(evidence_sets)
mem <- data.frame(gene = genes, stringsAsFactors = FALSE)
for (nm in set_cols) mem[[nm]] <- as.integer(mem$gene %in% evidence_sets[[nm]])
mem$key <- apply(mem[, set_cols, drop = FALSE], 1, function(r) paste(set_cols[as.integer(r) == 1L], collapse = " + "))
inter <- do.call(rbind, lapply(split(mem, mem$key), function(d) {
  present <- as.integer(colSums(d[, set_cols, drop = FALSE]) > 0)
  row <- data.frame(key = d$key[1], n = nrow(d), n_sets = sum(present),
                    genes = paste(sort(d$gene), collapse = ", "),
                    gene_lab = paste(sort(d$gene), collapse = "\n"),
                    stringsAsFactors = FALSE)
  for (j in seq_along(set_cols)) row[[set_cols[j]]] <- present[j]
  row
}))
inter <- inter[order(-inter$n_sets, -inter$n, inter$key), ]
rownames(inter) <- NULL
inter$idx <- seq_len(nrow(inter))
inter$fill <- ifelse(inter$n_sets >= 2L, "#1F4E79", "grey62")
inter$show_genes <- inter$n <= n_label_max
n_lab_lines <- if (any(inter$show_genes)) {
  max(vapply(strsplit(inter$gene_lab[inter$show_genes], "\n"), length, integer(1)))
} else 1L
y_top <- max(inter$n) + 1.2 + 0.38 * min(n_lab_lines, 8)

p_bar <- ggplot(inter, aes(x = idx, y = n)) +
  geom_col(aes(fill = fill), width = 0.65, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(aes(label = n), vjust = -0.25, size = 3, colour = "grey15") +
  geom_text(data = inter[inter$show_genes, , drop = FALSE], aes(label = gene_lab),
            vjust = -1.3, size = 2.2, lineheight = 0.9, colour = "grey10") +
  scale_x_continuous(breaks = inter$idx, expand = expansion(add = 0.55)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, y_top), clip = "off") +
  labs(y = "Intersection size", x = NULL,
       title = paste0("Set overlap (unique genes counted, not listed when n > ", n_label_max, ")")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(size = 11), plot.margin = margin(8, 4, 2, 8))

mat_long <- do.call(rbind, lapply(seq_len(nrow(inter)), function(i) {
  data.frame(idx = inter$idx[i], set = factor(set_cols, levels = rev(set_cols)),
             present = as.integer(unlist(inter[i, set_cols])), stringsAsFactors = FALSE)
}))
seg <- do.call(rbind, lapply(split(mat_long, mat_long$idx), function(d) {
  on <- d[d$present == 1L, ]
  if (nrow(on) < 2L) return(NULL)
  data.frame(idx = on$idx[1], y1 = on$set[which.min(as.integer(on$set))],
             y2 = on$set[which.max(as.integer(on$set))], stringsAsFactors = FALSE)
}))
p_mat <- ggplot(mat_long, aes(x = idx, y = set)) +
  geom_segment(data = seg, aes(x = idx, xend = idx, y = y1, yend = y2),
               inherit.aes = FALSE, colour = "#1F4E79", linewidth = 0.7) +
  geom_point(aes(colour = factor(present), size = factor(present))) +
  scale_colour_manual(values = c("0" = "grey80", "1" = "#1F4E79"), guide = "none") +
  scale_size_manual(values = c("0" = 2.0, "1" = 3.2), guide = "none") +
  scale_x_continuous(breaks = inter$idx, expand = expansion(add = 0.55)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(2, 4, 6, 8))

p_set <- ggplot(
  data.frame(set = factor(set_cols, levels = rev(set_cols)),
             n = vapply(evidence_sets, length, integer(1))),
  aes(x = n, y = set)
) +
  geom_col(fill = "#4A86C8", width = 0.65) +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = "Set size", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(2, 10, 6, 2))

p_up <- (p_bar + plot_spacer() + plot_layout(widths = c(4.2, 1.05))) /
  (p_mat + p_set + plot_layout(widths = c(4.2, 1.05))) +
  plot_layout(heights = c(1.45, 1))

out_pdf <- file.path(holly_out_dir, "forestMatrix_Baso_STRING.pdf")
p_out <- forest_p /
  (pr_p + mat_ge2 + plot_layout(widths = c(1.15, 1))) /
  p_up +
  plot_layout(heights = c(1.85, 1.85, 2.55))
ggsave(out_pdf, p_out, width = 10.4, height = 10.6, bg = "white")
message("[8.2.2 forestMatrix] wrote ", out_pdf)
utils::write.table(inter[, c("key", "n", "n_sets", "genes")],
                   file.path(holly_out_dir, "forestMatrix_Baso_STRING_intersections.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
