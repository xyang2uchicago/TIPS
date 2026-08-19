## 8.2.2_Caesar_forestMatrix_baso_TIPScll_TIPSmeta.R — forest plot + gene-evidence matrix (STRING),
## Basophil lineage. Replaces the venn3 diagrams in 8.2.2_Caesar_venn3_baso_TIPScll_TIPSmeta.R with
## one figure: a forest plot of pairwise Fisher-exact odds ratios (log scale, 95% CI), and a gene x
## evidence-source matrix for genes supported by >=2 sources.
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

## ---- gene-evidence matrix: genes supported by >=2 of {Weinreb, Cell TIPS, MuTrans} ----
evidence_sets <- list(`Weinreb fate` = baso, `Cell TIPS` = tips_cll, `MuTrans` = td_mutrans)
genes <- unique(unlist(evidence_sets))
pres <- sapply(evidence_sets, function(s) genes %in% s)
rownames(pres) <- genes
gene_keep <- rownames(pres)[rowSums(pres) >= 2]
if (!length(gene_keep)) stop("no gene has >=2 evidence sources")
message("[8.2.2 forestMatrix] genes with >=2 evidence sources (n=", length(gene_keep), "): ", paste(sort(gene_keep), collapse = ", "))

mat <- expand.grid(gene = gene_keep, source = names(evidence_sets), stringsAsFactors = FALSE)
mat$present <- pres[cbind(mat$gene, mat$source)]
mat$gene <- factor(mat$gene, levels = rev(sort(gene_keep)))
mat$source <- factor(mat$source, levels = names(evidence_sets))

matrix_p <- ggplot(mat, aes(x = source, y = gene, fill = present)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_manual(values = c(`TRUE` = "#4A86C8", `FALSE` = "grey92"), guide = "none") +
  labs(x = NULL, y = NULL, title = paste0(pub_lineage, ": genes supported by \u22652 evidence sources")) +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())

combined <- forest_p / matrix_p + plot_layout(heights = c(length(comparisons), max(3, length(gene_keep)) * 0.5))
out_pdf <- file.path(holly_out_dir, "forestMatrix_Baso_STRING.pdf")
ggsave(out_pdf, combined, width = 9, height = 5 + 0.3 * max(3, length(gene_keep)), bg = "white")
message("[8.2.2 forestMatrix] wrote ", out_pdf)
