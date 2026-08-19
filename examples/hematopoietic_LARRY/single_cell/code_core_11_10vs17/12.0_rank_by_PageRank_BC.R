## 12.0_rank_by_PageRank_BC.R — Felix TIPS_core (LARRY MuTrans; native mouse symbols)
##
## CODING FIX (not strategy): Felix rank_TF_CHD_in_PPIN() does pdf() then print() then dev.off().
## If print() errors, dev.off() is skipped → corrupt CP_rank_gene_by_*.pdf. finalize_pdf_devices()
## closes any leaked device. Pad signatures to 5 panels (Felix patchwork layout uses [1:5]).

########## BEGINNING OF USER INPUT ##########
code_dir <- get0("code_dir", ifnotfound = file.path(here::here("examples", "hematopoietic_LARRY", "single_cell"), "code_core_11_10vs17"))
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)
top_TF_rank <- 3
gene_top_n  <- 20
########## END OF USER INPUT ##########

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

if (requireNamespace("dorothea", quietly = TRUE)) {
  data(dorothea_mm, package = "dorothea", envir = environment())
  TF_mouse <- toupper(unique(dorothea_mm$tf))
} else {
  TF_mouse <- character(0)
}

source(file.path(tips_r_dir, paste0("celltype_specific_weight_v", celltype_specific_weight_version, ".R")))

finalize_pdf_devices <- function() {
  while (grDevices::dev.cur() > 1L) grDevices::dev.off()
}

## Felix layout: (p[[sig1]]|p[[sig2]]|p[[sig3]]) / (legend|p[[sig4]]|p[[sig5]])
pad_signatures_for_rank_plot <- function(sigs, n = 5L) {
  sigs <- unique(as.character(sigs))
  if (!length(sigs)) return(sigs)
  while (length(sigs) < n) sigs <- c(sigs, tail(sigs, 1L))
  sigs[seq_len(n)]
}

db_specifc_output_path <- ppi_path
df_PageRank <- readRDS(file.path(db_specifc_output_path, "df_PAGERANK_strength_ANND.rewiring.P.rds"))
df_PageRank$gene <- toupper(as.character(df_PageRank$gene))

signatures_found <- intersect(
  c(paste0("HiG_", CP_cluster), paste0("CTS_", CP_cluster), paste0("HiGCTS_", CP_cluster)),
  unique(df_PageRank$signature)
)
signatures_plot <- pad_signatures_for_rank_plot(signatures_found)
HiGCTS_sig <- paste0("HiGCTS_", CP_cluster)
CTS_sig <- paste0("CTS_", CP_cluster)

res_pr <- tryCatch(
  rank_TF_CHD_in_PPIN(
    df_PageRank, CHD, TF_mouse,
    signatures = signatures_plot,
    key = "PageRank", top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
  ),
  error = function(e) {
    finalize_pdf_devices()
    stop(e)
  },
  finally = finalize_pdf_devices()
)

df_betweenness <- read.table(
  file.path(db_specifc_output_path, "df_betweeness.tsv"),
  sep = "\t", header = TRUE
)
df_betweenness$gene <- toupper(as.character(df_betweenness$gene))

res_bw <- tryCatch(
  rank_TF_CHD_in_PPIN(
    df_betweenness, CHD, TF_mouse,
    signatures = signatures_plot,
    key = "BetweennessCentrality",
    top_TF_rank = top_TF_rank, gene_top_n = gene_top_n, saveFigure = TRUE
  ),
  error = function(e) {
    finalize_pdf_devices()
    stop(e)
  },
  finally = finalize_pdf_devices()
)

pick_seed_tfs_from_sig <- function(res_pr, res_bw, sig) {
  if (!sig %in% names(res_pr) || is.null(res_pr[[sig]])) return(character())
  pr_genes <- as.character(subset(res_pr[[sig]])$gene)
  pr_genes <- pr_genes[nzchar(pr_genes)]
  if (!length(pr_genes)) return(character())
  bw_sub <- res_bw[[sig]]
  if (is.null(bw_sub) || !nrow(bw_sub)) return(pr_genes)
  bw_genes <- as.character(subset(bw_sub, BetweennessCentrality > 0)$gene)
  seed <- intersect(pr_genes, bw_genes)
  if (!length(seed)) seed <- pr_genes
  seed
}

seed_TF_auto <- character()
seed_source <- NULL
if (HiGCTS_sig %in% signatures_found) {
  seed_TF_auto <- pick_seed_tfs_from_sig(res_pr, res_bw, HiGCTS_sig)
  seed_source <- HiGCTS_sig
} else if (CTS_sig %in% signatures_found) {
  seed_TF_auto <- pick_seed_tfs_from_sig(res_pr, res_bw, CTS_sig)
  seed_source <- CTS_sig
  message(
    "[12.0] No ", HiGCTS_sig, " — using ", CTS_sig,
    " (top ", top_TF_rank, " TF in top ", gene_top_n, " PageRank)"
  )
}

if (length(seed_TF_auto) && !length(seed_TF)) {
  seed_TF <- seed_TF_auto
  message("[12.0] seed_TF from ", seed_source, ": ", paste(seed_TF, collapse = ", "))
} else if (!length(seed_TF_auto)) {
  message(
    "[12.0] No ", HiGCTS_sig, " or ", CTS_sig,
    " — set TIPS_SEED_TF manually for 24.1 if needed"
  )
}

message("[12.0] next: rebuild_mat <- TRUE; source('24.1_acat_CTS.cisTarget_dualpull_clean.R')")
