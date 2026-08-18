## 8.2_venn_TIPS_two_figures.R — TIPS set overlap (Larry 4_9vs11)
## STRING: Fig 1–3 (venn3, venn4, venn5)
## IID:    Fig 1 only (venn3)
## Cross-PPI: venn2 TIPS_A4 (STRING vs IID)
## Fig 1 (3 sets): Mega + TD + TIPS  -> ggVenn
## Fig 2 (4 sets): Mega + CTS + TD + TIPS  -> ggVenn
## Fig 3 (5 sets): TIPS + DEG A4 + DEG A11 + Mega + TD -> ggVenn (5 circles)
## Intersection gene tables exported as TSV
##
## Needs: readxl, ggVennDiagram
##
##   source("F:/projects/TIPS/source/GSE140802_lineage_tracking/8.2_venn_TIPS_two_figures.R")

########## USER INPUT (shared) ##########
GSE140802_ROOT <- Sys.getenv(
  "GSE140802_ROOT",
  "F:/projects/TIPS/source/GSE140802_lineage_tracking"
)
TAG <- "4_9vs11"
TD.cut <- 0.7
td_csv_dir <- Sys.getenv(
  "MUTRANS_TD_CSV_DIR",
  "F:/projects/TIPS/results/GSE140802_lineage_tracking/larry/figures"
)
pub_xlsx <- Sys.getenv(
  "WEINREB_TABLE_S3",
  "F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/aaw3381-Weinreb-Table-S3.xlsx"
)
pub_lineage <- "Megakaryocyte"
biotip_cts_rdata <- Sys.getenv(
  "BIOTIP_CTS_RDATA",
  "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/larry_BioTIP/BioTIP_attractor_07252026_used/CTS_Lib_data.RData"
)
########## END USER INPUT ##########

for (pkg in c("readxl", "ggVennDiagram")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install ", pkg, ": install.packages('", pkg, "')", call. = FALSE)
  }
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggVennDiagram)
})

norm_genes <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  toupper(x)
}

pairwise_intersection_table <- function(gene_sets) {
  labs <- names(gene_sets)
  rows <- vector("list", length(labs) * (length(labs) - 1L) / 2L)
  k <- 0L
  for (i in seq_len(length(labs) - 1L)) {
    for (j in (i + 1L):length(labs)) {
      inter <- intersect(gene_sets[[i]], gene_sets[[j]])
      k <- k + 1L
      rows[[k]] <- data.frame(
        set_a = labs[i],
        set_b = labs[j],
        n = length(inter),
        genes = paste(sort(inter), collapse = ", "),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

intersection_row <- function(gene_sets, set_names, label) {
  genes <- Reduce(intersect, gene_sets[set_names])
  data.frame(
    region = label,
    n = length(genes),
    genes = paste(sort(genes), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

VENN_SCRIPT_VERSION <- "20260728_5set"  # check: print(VENN_SCRIPT_VERSION)

plot_venn_gg <- function(gene_sets, title, outfile, width = NULL, height = NULL) {
  n_sets <- length(gene_sets)
  if (n_sets < 2L) {
    stop("Need at least 2 gene sets for a Venn diagram.")
  }
  cat_names <- names(gene_sets)
  if (is.null(cat_names) || length(cat_names) != n_sets || !all(nzchar(cat_names))) {
    stop("Each gene set needs a non-empty name for Venn labels.")
  }
  if (is.null(width)) {
    width <- switch(as.character(n_sets), `2` = 8, `3` = 9, `4` = 10, `5` = 12, `6` = 13, `7` = 14, 11)
  }
  if (is.null(height)) {
    height <- switch(as.character(n_sets), `2` = 7, `3` = 8, `4` = 9, `5` = 11, `6` = 12, `7` = 13, 10)
  }
  p <- ggVennDiagram(
    gene_sets,
    category.names = cat_names,
    label = "count",
    label_size = if (n_sets <= 3L) 3.4 else if (n_sets == 4L) 3.0 else 2.5,
    set_size = if (n_sets <= 3L) 4.8 else if (n_sets == 4L) 4.2 else 3.6,
    set_color = "black"
  ) +
    scale_fill_gradient(low = "#F8FBFF", high = "#4A86C8") +
    labs(title = title) +
    coord_cartesian(clip = "off") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.margin = margin(14, 24, 14, 24),
      legend.position = "none"
    )
  ggsave(outfile, p, width = width, height = height, bg = "white")
  message("[8.2] wrote ", outfile)
  invisible(p)
}

run_venn_arm <- function(ppi_arm, which_figs = c(1L, 2L, 3L)) {
  if (!1L %in% which_figs) {
    stop("which_figs must include 1 (venn3)", call. = FALSE)
  }
  message("[8.2] ===== ", ppi_arm, " ", TAG, " (figs: ", paste(which_figs, collapse = ","), ") =====")
  tips_wd <- paste0(
    normalizePath(
      file.path(GSE140802_ROOT, paste0("7_data_MuTrans_TIPS_", ppi_arm)),
      winslash = "/", mustWork = FALSE
    ),
    "/"
  )
  code_dir <- file.path(tips_wd, paste0("code_core_", TAG))
  source(file.path(code_dir, "00_configuration.R"))
  tips_configure(TAG = TAG, wd = tips_wd)

  venn_dir <- file.path(results_dir, "venn")
  dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(venn_dir)
  message("[8.2] results_dir = ", results_dir)
  message("[8.2] all figures/tables -> ", venn_dir)

  if (!file.exists(deg_rdata)) stop("Run 11.1 first — missing: ", deg_rdata)
  DEG <- lapply(readRDS(deg_rdata), norm_genes)

  if (2L %in% which_figs) {
    if (!file.exists(biotip_cts_rdata)) stop("Missing BioTIP CTS file: ", biotip_cts_rdata)
    load(biotip_cts_rdata)
    stopifnot(exists("CTS.Lib.Symbol"))
    cts_genes <- norm_genes(CTS.Lib.Symbol[[CTS_ID]])
  }

  pub_df <- readxl::read_excel(pub_xlsx, sheet = "DGE of progenitors in vitro")
  pub_df <- pub_df[!is.na(pub_df[["Gene symbol"]]) & nzchar(as.character(pub_df[["Gene symbol"]])), , drop = FALSE]
  pub_genes <- split(norm_genes(pub_df[["Gene symbol"]]), pub_df$Lineage)
  if (!pub_lineage %in% names(pub_genes)) {
    stop("Lineage not in Weinreb table: ", pub_lineage)
  }
  pub_mega <- pub_genes[[pub_lineage]]

  td_file <- file.path(
    td_csv_dir,
    paste0("td_genes_scores_", CP_cluster, "_to_", CF_cluster, "_seacell.csv")
  )
  if (!file.exists(td_file)) stop("Missing TD file: ", td_file)
  td_df <- read.csv(td_file, stringsAsFactors = FALSE)
  td_genes <- norm_genes(td_df$Genes[abs(td_df$corr) > TD.cut])

  tips_table <- file.path(
    results_dir,
    paste0("cisTarget_predicted_", CTS_ID),
    paste0("PPI_graph_GRN_prediction_CTS_", CTS_ID, "_dualpull_final_table.tsv")
  )
  if (!file.exists(tips_table)) stop("Run 24.1 first — missing: ", tips_table)
  TIPS <- read.delim(tips_table, stringsAsFactors = FALSE)
  tips_cf <- subset(TIPS, linkage == "CF", select = c(from, to))
  TIPS_A4 <- norm_genes(c(tips_cf$from, tips_cf$to))

  if (any(c(2L, 3L) %in% which_figs)) {
    deg_a4 <- DEG[[CP_cluster]]
    deg_a11 <- DEG[[CF_cluster]]
  }

  lab_tips <- "TIPS"
  lab_deg4 <- paste0("DEG A", CP_cluster)
  lab_deg11 <- paste0("DEG A", CF_cluster)

  ## Fig 1 (venn3): Mega + TD + TIPS
  venn_pub <- list(
    `Mega (Weinreb)` = pub_mega,
    `TD (A4->A11)` = td_genes,
    TIPS = TIPS_A4
  )
  plot_venn_gg(
    venn_pub,
    paste0("Mega x TD x TIPS (", pub_lineage, "; ", ppi_arm, ")"),
    file.path(venn_dir, "venn3_Mega_TD_TIPS_ggplot.pdf")
  )
  write.table(
    pairwise_intersection_table(venn_pub),
    file.path(venn_dir, "intersections_venn3_pairwise.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    intersection_row(venn_pub, names(venn_pub), "Mega & TD & TIPS"),
    file.path(venn_dir, "intersections_venn3_all_Mega_TD_TIPS.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # ## Fig 2 (venn4): Mega + CTS + TD + TIPS
  # if (!2L %in% which_figs) {
  #   message("[8.2] skipping venn4/venn5 for ", ppi_arm)
  #   return(invisible(list(venn_dir = venn_dir, TIPS_A4 = TIPS_A4)))
  # }
  # venn_pub <- list(
  #   `Mega (Weinreb)` = pub_mega,
  #   `CTS (A11 vs A4)` = cts_genes,
  #   `TD (A4->A11)` = td_genes,
  #   TIPS = TIPS_A4
  # )
  # plot_venn_gg(
  #   venn_pub,
  #   paste0("Mega x CTS x TD x TIPS (", pub_lineage, "; ", ppi_arm, ")"),
  #   file.path(venn_dir, "venn4_Mega_CTS_TD_TIPS_ggplot.pdf")
  # )
  # write.table(
  #   pairwise_intersection_table(venn_pub),
  #   file.path(venn_dir, "intersections_venn4_pairwise.tsv"),
  #   sep = "\t", quote = FALSE, row.names = FALSE
  # )
  # write.table(
  #   intersection_row(venn_pub, names(venn_pub), "Mega & CTS & TD & TIPS"),
  #   file.path(venn_dir, "intersections_venn4_all_Mega_CTS_TD_TIPS.tsv"),
  #   sep = "\t", quote = FALSE, row.names = FALSE
  # )

#   ## Fig 3 (venn5): TIPS + DEG A4 + DEG A11 + Mega + TD (5 circles)
#   if (!3L %in% which_figs) {
#     return(invisible(list(venn_dir = venn_dir, TIPS_A4 = TIPS_A4)))
#   }
#   venn5 <- setNames(
#     list(TIPS_A4, deg_a4, deg_a11, pub_mega, td_genes),
#     c(lab_tips, lab_deg4, lab_deg11, "Mega (Weinreb)",
#       paste0("TD (A", CP_cluster, "->A", CF_cluster, ")"))
#   )
#   plot_venn_gg(
#     venn5,
#     paste0("TIPS x DEG A", CP_cluster, " x DEG A", CF_cluster,
#            " x Mega x TD (", pub_lineage, "; ", ppi_arm, ")"),
#     file.path(venn_dir, "venn5_TIPS_DEG4_DEG11_Mega_TD_ggplot.pdf")
#   )
#   write.table(
#     pairwise_intersection_table(venn5),
#     file.path(venn_dir, "intersections_venn5_pairwise.tsv"),
#     sep = "\t", quote = FALSE, row.names = FALSE
#   )
#   write.table(
#     intersection_row(
#       venn5, names(venn5),
#       paste0("TIPS & DEG A", CP_cluster, " & A", CF_cluster, " & Mega & TD")
#     ),
#     file.path(venn_dir, "intersections_venn5_all5.tsv"),
#     sep = "\t", quote = FALSE, row.names = FALSE
#   )
#   write.table(
#     intersection_row(
#       venn5,
#       c(lab_tips, lab_deg4, lab_deg11),
#       paste0("TIPS & DEG A", CP_cluster, " & A", CF_cluster)
#     ),
#     file.path(venn_dir, "intersections_venn5_triple_TIPS_DEG4_DEG11.tsv"),
#     sep = "\t", quote = FALSE, row.names = FALSE
#   )

#   message("[8.2] TIPS ∩ DEG A", CP_cluster, ": ",
#           paste(intersect(TIPS_A4, deg_a4), collapse = ", "))
#   message("[8.2] TIPS ∩ DEG A", CF_cluster, ": ",
#           paste(intersect(TIPS_A4, deg_a11), collapse = ", "))
#   message(
#     "[8.2] TIPS ∩ DEG A", CP_cluster, " ∩ A", CF_cluster, ": ",
#     paste(Reduce(intersect, list(TIPS_A4, deg_a4, deg_a11)), collapse = ", ")
#   )
#   invisible(list(venn_dir = venn_dir, TIPS_A4 = TIPS_A4))
# }

jaccard_score <- function(x, y) {
  x <- norm_genes(x)
  y <- norm_genes(y)
  u <- union(x, y)
  if (!length(u)) return(NA_real_)
  length(intersect(x, y)) / length(u)
}

message("[8.2] plot_venn_gg loaded (", VENN_SCRIPT_VERSION, ")")
res_string <- run_venn_arm("STRING", which_figs = c(1L, 2L, 3L))
res_iid <- run_venn_arm("IID", which_figs = 1L)

## Fig cross-PPI: TIPS_A4 STRING vs IID (2 circles)
venn_tips2 <- list(
  `TIPS (STRING)` = res_string$TIPS_A4,
  `TIPS (IID)` = res_iid$TIPS_A4
)
compare_dir <- file.path(
  GSE140802_ROOT,
  "7_data_MuTrans_TIPS_STRING",
  paste0("results_core_", TAG),
  "venn"
)
dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)

plot_venn_gg(
  venn_tips2,
  paste0("TIPS_A4: STRING vs IID (", TAG, ")"),
  file.path(compare_dir, "venn2_TIPS_A4_STRING_vs_IID_ggplot.pdf"),
  width = 8,
  height = 7
)
tips2_tab <- pairwise_intersection_table(venn_tips2)
tips2_tab$jaccard <- jaccard_score(venn_tips2[[1]], venn_tips2[[2]])
write.table(
  tips2_tab,
  file.path(compare_dir, "intersections_venn2_TIPS_A4_STRING_vs_IID.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    region = "TIPS_A4 STRING & IID",
    n = length(intersect(venn_tips2[[1]], venn_tips2[[2]])),
    genes = paste(sort(intersect(venn_tips2[[1]], venn_tips2[[2]])), collapse = ", "),
    stringsAsFactors = FALSE
  ),
  file.path(compare_dir, "intersections_venn2_TIPS_A4_STRING_vs_IID_shared.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message(
  "[8.2] TIPS_A4 STRING (n=", length(venn_tips2[[1]]), ") vs IID (n=", length(venn_tips2[[2]]),
  ") Jaccard=", round(tips2_tab$jaccard[1], 3),
  " | shared: ", paste(sort(intersect(venn_tips2[[1]], venn_tips2[[2]])), collapse = ", ")
)


# Shared (6): BHLHE41, COL1A2, GATA1, GATA2, PCYT1B, SPTA1 → Jaccard ≈ 0.667
# STRING-only: HBB-BS
# IID-only: IKZF2, SPTB