## 8.2.2_Caesar_venn3_baso_TIPScll_TIPSmeta.R — venn3 diagrams + Fisher exact tests, Basophil lineage (STRING)
## Sets: Basophil (Weinreb) x {TIPS_cll (C11->C10, cell-level TIPS CM arm) | cell_CTS (C11 CTS, BioTIP module)} x TD_MuTrans (larry_transition_4_to_9).
## Companion to 8.2.1_Caesar_FET_venn3_Mega_TD_TIPS_simplified.R (does NOT modify or source it).
## Metacell TIPS pipeline has no CM/dualpull output for the 4->9 arm, so the TD-scored genes stand in for it.
##   source("F:/projects/TIPS/source/GSE140802_lineage_tracking/8.2.2_Caesar_venn3_baso_TIPScll_TIPSmeta.R")

style = 'biased' # 'leanning'  because the fate-associted gene liste have shared genes between different fates, we assess its overlapping with TIPS outcome using the full TIPS model genes.

for (pkg in c("readxl", "ggVennDiagram")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, call. = FALSE)
}
suppressPackageStartupMessages({ library(ggplot2); library(ggVennDiagram) })
## ---- config ----
GSE140802_cell_ROOT <- Sys.getenv("GSE140802_cell_ROOT", "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory")
pub_xlsx    <- Sys.getenv("WEINREB_TABLE_S3", "F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/aaw3381-Weinreb-Table-S3.xlsx")
pub_lineage <- "Basophil"; pub_sheet <- "DGE of progenitors in vitro"
TD.cut <- 0.7
td_csv_dir <- Sys.getenv("MUTRANS_TD_CSV_DIR", here::here("examples", "hematopoietic_LARRY", "data", "larry_figures"))
CP_cluster <- 4; CM_cluster <- 9  # larry_transition_4_to_9: TAG=4_9vs11's progenitor (4) -> Basophil-fate cluster (9)
fisher_alt <- Sys.getenv("TIPS_FET_ALTERNATIVE", "two.sided")
tingjun_biotip <- Sys.getenv("TINGJUN_BIOTIP_ROOT", "F:/projects/TIPS/results/GSE140802_lineage_tracking/Tingjunl/BioTIP")
tingjun_tips   <- Sys.getenv("TINGJUN_TIPS_ROOT", "F:/projects/TIPS/results/GSE140802_lineage_tracking/Tingjunl/TIPS/7_scaledata_leiden_r0_8_TIPS_STRING")
holly_out_dir  <- file.path(tingjun_tips, "Holly")
norm_genes <- function(x) { x <- unique(toupper(as.character(x))); x[!is.na(x) & nzchar(x)] }

## ---- Basophil (Weinreb) ----
pub_df <- readxl::read_excel(pub_xlsx, sheet = pub_sheet)
pub_df <- pub_df[!is.na(pub_df$`Gene symbol`) & nzchar(as.character(pub_df$`Gene symbol`)), ]
baso <- norm_genes(pub_df$`Gene symbol`[pub_df$Lineage == pub_lineage])

## ---- TIPS_cll (C11->C10): CM-linkage nodes of the cell-level TIPS run ----
## A gene shared with CF's nodes is kept here only if it has a delta>0 edge on the CM side (a
## decreasing-only appearance isn't a real driving signal) - e.g. MYC: delta>0 in CM, <0 in CF.
tips_cll_df <- utils::read.delim(file.path(tingjun_tips, "results_core_11_10vs17", "cisTarget_predicted_11",
                                            "PPI_graph_GRN_prediction_CTS_11_dualpull_final_table.tsv"), stringsAsFactors = FALSE)
cf_nodes <- norm_genes(with(tips_cll_df[tips_cll_df$linkage == "CF", ], c(from, to)))
cm_nodes <- norm_genes(with(tips_cll_df[tips_cll_df$linkage == "CM", ], c(from, to)))

if(style =='leanning') {
	cm_pos   <- norm_genes(with(tips_cll_df[tips_cll_df$linkage == "CM" & tips_cll_df$delta > 0, ], c(from, to)))
	tips_cll <- setdiff(cm_nodes, setdiff(intersect(cf_nodes, cm_nodes), cm_pos))
} else {
	tips_cll <- cm_nodes
}

## ---- cell_CTS = "C11 CTS" BioTIP module intersected with the cell-level tested gene universe ----
## (reproduces panel E of make_Sx_panels.R: CTS.Lib.Symbol[["11"]] ∩ rownames(testres); same set as 8.2.1's cell_CTS)
Bdir <- file.path(tingjun_biotip, "results/BioTIP_leiden_r0_8")
e <- new.env(); load(file.path(Bdir, "BioTIP_leiden_r0_8CTS_Lib_Scaledata.RData"), e)
o <- new.env(); load(file.path(Bdir, "BioTIP_leiden_r0_8optimized_test_sd_selection_Scaledata.RData"), o)
cts_uni <- unique(unlist(lapply(get("testres", o), rownames)))
cell_CTS <- norm_genes(intersect(get("CTS.Lib.Symbol", e)[["11"]], cts_uni))

## ---- TD_MuTrans (larry_transition_4_to_9): MuTrans TD score, NOT a TIPS/dualpull output (none exists for the CM arm) ----
td_df <- utils::read.csv(file.path(td_csv_dir, paste0("td_genes_scores_", CP_cluster, "_to_", CM_cluster, "_seacell.csv")), stringsAsFactors = FALSE)
td_mutrans <- norm_genes(td_df$Genes[abs(td_df$corr) > TD.cut])

## ---- background universe: cell-level HVG (3k, cell-cycle-cor RDS) union Weinreb Table S3 ----
hvg <- norm_genes(names(readRDS(file.path(GSE140802_cell_ROOT, "cell_cycle_signature_cor_3khvg.rds"))))
universe <- unique(c(hvg, norm_genes(pub_df$`Gene symbol`)))
message("[8.2.2 Caesar] universe n=", length(universe), "; Basophil n=", length(baso), "; TIPS_cll(C11->C10) n=", length(tips_cll),
        "; cell_CTS(C11 CTS) n=", length(cell_CTS), "; TD_MuTrans(A4->A9) n=", length(td_mutrans))

## ---- pairwise Fisher's exact test on a 2x2 overlap table ----
fisher_overlap <- function(a, b, la, lb, U) {
  A <- intersect(a, U); B <- intersect(b, U); N <- length(U)
  k <- length(intersect(A, B)); na <- length(A); nb <- length(B)
  tab <- matrix(c(k, na - k, nb - k, N - na - nb + k), nrow = 2, byrow = TRUE)
  ft <- stats::fisher.test(tab, alternative = fisher_alt)
  data.frame(set_a = la, set_b = lb, n_a = na, n_b = nb, overlap_k = k, universe_n = N,
             p_value = ft$p.value, odds_ratio = unname(ft$estimate),
             overlap_genes = paste(sort(intersect(A, B)), collapse = ", "), stringsAsFactors = FALSE)
}
pairwise_fet <- function(sets, U) { combos <- combn(names(sets), 2, simplify = FALSE); do.call(rbind, lapply(combos, function(p) fisher_overlap(sets[[p[1]]], sets[[p[2]]], p[1], p[2], U))) }

## ---- 3-circle venn diagram ----
plot_venn3 <- function(sets, title, outfile) {
  p <- ggVennDiagram(sets, label = "count", set_color = "black") +
    scale_fill_gradient(low = "#F8FBFF", high = "#4A86C8") +
    labs(title = title) + theme(legend.position = "none")
  ggsave(outfile, p, width = 9, height = 8, bg = "white")
  message("[8.2.2 Caesar] wrote ", outfile)
}

## ---- run venn3 + FET (STRING); called once per middle circle (TIPS_cll, cell_CTS) ----
run_venn_fet <- function(sets, tag) {
  plot_venn3(sets, paste0(paste(names(sets), collapse = " x "), " (STRING)"),
             file.path(holly_out_dir, paste0("venn3_", tag, "_STRING.pdf")))
  fet <- pairwise_fet(sets, universe)
  fet$ppi_arm <- "STRING"
  utils::write.table(fet, file.path(holly_out_dir, paste0("FET_venn3_", tag, "_STRING.tsv")),
                      sep = "\t", quote = FALSE, row.names = FALSE)
  message("[8.2.2 Caesar] ", tag, ": ", paste(sprintf("%s vs %s k=%d p=%.4g", fet$set_a, fet$set_b, fet$overlap_k, fet$p_value), collapse = " | "))
  fet
}
dir.create(holly_out_dir, recursive = TRUE, showWarnings = FALSE)
fet_tipscll <- run_venn_fet(list(`Basophil (Weinreb)` = baso, `TIPS_cll (C11->C10)` = tips_cll, `TD_MuTrans (A4->A9)` = td_mutrans), "Baso_TIPScll_TDMuTrans")
fet_cellcts <- run_venn_fet(list(`Basophil (Weinreb)` = baso, `cell_CTS (C11 CTS)` = cell_CTS, `TD_MuTrans (A4->A9)` = td_mutrans), "Baso_cellCTS_TDMuTrans")

# [8.2.2 Caesar] universe n=3239; Basophil n=67; TIPS_cll(C11->C10) n=3; cell_CTS(C11 CTS) n=71; TD_MuTrans(A4->A9) n=108
# [8.2.2 Caesar] wrote F:/projects/TIPS/results/GSE140802_lineage_tracking/Tingjunl/TIPS/7_scaledata_leiden_r0_8_TIPS_STRING/Holly/venn3_Baso_TIPScll_TDMuTrans_STRING.pdf
# [8.2.2 Caesar] Baso_TIPScll_TDMuTrans: Basophil (Weinreb) vs TIPS_cll (C11->C10) k=1 p=0.0608 | Basophil (Weinreb) vs TD_MuTrans (A4->A9) k=7 p=0.001052 | TIPS_cll (C11->C10) vs TD_MuTrans (A4->A9) k=0 p=1
# [8.2.2 Caesar] wrote F:/projects/TIPS/results/GSE140802_lineage_tracking/Tingjunl/TIPS/7_scaledata_leiden_r0_8_TIPS_STRING/Holly/venn3_Baso_cellCTS_TDMuTrans_STRING.pdf
# [8.2.2 Caesar] Baso_cellCTS_TDMuTrans: Basophil (Weinreb) vs cell_CTS (C11 CTS) k=2 p=0.658 | Basophil (Weinreb) vs TD_MuTrans (A4->A9) k=7 p=0.001052 | cell_CTS (C11 CTS) vs TD_MuTrans (A4->A9) k=8 p=0.0002689
