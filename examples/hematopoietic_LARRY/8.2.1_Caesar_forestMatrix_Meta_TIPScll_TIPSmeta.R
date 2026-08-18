## 8.2.1_Caesar_forestMatrix_Meta_TIPScll_TIPSmeta.R — forest plot + gene-evidence matrix (STRING),
## Megakaryocyte lineage. Replaces the venn3 diagrams in 8.2.1_Caesar_venn3_Meta_TIPScll_TIPSmeta.R
## with one figure: a forest plot of 5 pairwise Fisher-exact odds ratios (log scale, 95% CI), and a
## gene x evidence-source matrix for genes supported by >=2 sources. Same 4 gene sets and `style`
## convention as that script; companion only, does not source/modify it.
##   source("F:/projects/TIPS/source/GSE140802_lineage_tracking/8.2.1_Caesar_forestMatrix_Meta_TIPScll_TIPSmeta.R")

for (pkg in c("readxl", "ggplot2", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, call. = FALSE)
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

## ---- config ----
GSE140802_cell_ROOT <- Sys.getenv("GSE140802_cell_ROOT", "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory")
GSE140802_ROOT <- Sys.getenv("GSE140802_ROOT", "F:/projects/TIPS/source/GSE140802_lineage_tracking")
TAG <- "4_9vs11"; TD.cut <- 0.7
pub_xlsx   <- Sys.getenv("WEINREB_TABLE_S3", "F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/aaw3381-Weinreb-Table-S3.xlsx")
td_csv_dir <- Sys.getenv("MUTRANS_TD_CSV_DIR", "F:/projects/TIPS/results/GSE140802_lineage_tracking/larry/figures")
pub_lineage <- "Megakaryocyte"; pub_sheet <- "DGE of progenitors in vitro"
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
load_cell_tips_c11 <- function() {
  f <- file.path(tingjun_tips, "results_core_11_10vs17", "cisTarget_predicted_11", "PPI_graph_GRN_prediction_CTS_11_dualpull_final_table.tsv")
  linkage_set(utils::read.delim(f, stringsAsFactors = FALSE), "CF", "CM")
}

## ---- gene sets: Weinreb Mega, cell-level TIPS (C17vC11), metacell TIPS (S11vS4), MuTrans TD (S4->S11) ----
tips_wd <- paste0(normalizePath(file.path(GSE140802_ROOT, "7_data_MuTrans_TIPS_STRING"), winslash = "/", mustWork = FALSE), "/")
source(file.path(tips_wd, paste0("code_core_", TAG), "00_configuration.R"))
tips_configure(TAG = TAG, wd = tips_wd)
pub_df <- readxl::read_excel(pub_xlsx, sheet = pub_sheet)
pub_df <- pub_df[!is.na(pub_df$`Gene symbol`) & nzchar(as.character(pub_df$`Gene symbol`)), ]
mega <- norm_genes(pub_df$`Gene symbol`[pub_df$Lineage == pub_lineage])
tips_df <- utils::read.delim(file.path(results_dir, paste0("cisTarget_predicted_", CTS_ID), paste0("PPI_graph_GRN_prediction_CTS_", CTS_ID, "_dualpull_final_table.tsv")), stringsAsFactors = FALSE)
tips_meta <- linkage_set(tips_df, "CF", "CM")
cell_tips <- load_cell_tips_c11()
td_df <- utils::read.csv(file.path(td_csv_dir, paste0("td_genes_scores_", CP_cluster, "_to_", CF_cluster, "_seacell.csv")), stringsAsFactors = FALSE)
td <- norm_genes(td_df$Genes[abs(td_df$corr) > TD.cut])

## ---- background universe (cell-level HVG union Weinreb Table S3) ----
hvg <- norm_genes(names(readRDS(file.path(GSE140802_cell_ROOT, "cell_cycle_signature_cor_3khvg.rds"))))
universe <- unique(c(hvg, norm_genes(pub_df$`Gene symbol`)))
message("[8.2.1 forestMatrix] universe n=", length(universe), "; Mega n=", length(mega), "; cell_TIPS n=",
        length(cell_tips), "; metacell_TIPS n=", length(tips_meta), "; MuTrans_TD n=", length(td))

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

## ---- the 5 comparisons, in the requested order; red = vs Weinreb, blue = method agreement ----
comparisons <- list(
  list(label = "Cell-level TIPS vs Weinreb fate genes", a = cell_tips, b = mega,      group = "vs_Weinreb"),
  list(label = "Metacell TIPS vs Weinreb fate genes",   a = tips_meta, b = mega,      group = "vs_Weinreb"),
  list(label = "MuTrans drivers vs Weinreb fate genes", a = td,        b = mega,      group = "vs_Weinreb"),
  list(label = "Cell-level TIPS vs MuTrans drivers",    a = cell_tips, b = td,        group = "method_agreement"),
  list(label = "Cell-level vs metacell TIPS",           a = cell_tips, b = tips_meta, group = "method_agreement")
)
fet <- do.call(rbind, lapply(comparisons, function(cmp) cbind(fisher_overlap(cmp$a, cmp$b, universe), label = cmp$label, group = cmp$group)))
fet$label <- factor(fet$label, levels = rev(vapply(comparisons, `[[`, "", "label")))  # top-to-bottom = requested order
dir.create(holly_out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(fet, file.path(holly_out_dir, "forestMatrix_Mega_FET.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("[8.2.1 forestMatrix] ", paste(sprintf("%s: k=%d OR=%.3g [%.3g,%.3g] p=%.4g", fet$label, fet$overlap_k,
        fet$odds_ratio, fet$ci_low, fet$ci_high, fet$p_value), collapse = " | "))

## ---- forest plot: log-scale OR with 95% CI, labeled with k and P; clip only the plotted position ----
## (very large/undefined ORs can arise from small modules -- the reported values above are never clipped)
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

## ---- gene-evidence matrix: genes supported by >=2 of {Weinreb, cell TIPS, metacell TIPS, MuTrans} ----
evidence_sets <- list(`Weinreb fate` = mega, `Cell TIPS` = cell_tips, `Metacell TIPS` = tips_meta, `MuTrans` = td)
genes <- unique(unlist(evidence_sets))
pres <- sapply(evidence_sets, function(s) genes %in% s)
rownames(pres) <- genes
gene_keep <- rownames(pres)[rowSums(pres) >= 2]
if (!length(gene_keep)) stop("no gene has >=2 evidence sources")
message("[8.2.1 forestMatrix] genes with >=2 evidence sources (n=", length(gene_keep), "): ", paste(sort(gene_keep), collapse = ", "))

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
out_pdf <- file.path(holly_out_dir, "forestMatrix_Mega_STRING.pdf")
ggsave(out_pdf, combined, width = 9, height = 5 + 0.3 * max(3, length(gene_keep)), bg = "white")
message("[8.2.1 forestMatrix] wrote ", out_pdf)
