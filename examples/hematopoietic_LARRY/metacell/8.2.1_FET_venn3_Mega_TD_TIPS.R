## 8.2.1_FET_venn3_Mega_TD_TIPS.R — Fisher exact tests for venn3 (Panel F stats)
##
## Companion to 8.2_venn_TIPS_two_figures.R (does NOT modify that script).
##
## Confirmed settings:
##   1) Background = union(~3000 HVG in SCE, all Weinreb Table S3 in vitro symbols)
##   2) TIPS = CF signature genes from dualpull_final_table (per PPI backbone):
##        STRING n=7 (BHLHE41/GATA1/GATA2 CF edges; see tips_cf_genes())
##        IID    n=8 (all CF edges)
##   3) SCE from seu_attractor_MuTrans_HVG.rds via load_larry_sce()
##   4) fisher.test two-sided (R default)
##
##   source("F:/projects/TIPS/source/GSE140802_lineage_tracking/8.2.1_FET_venn3_Mega_TD_TIPS.R")

########## USER INPUT (keep in sync with 8.2) ##########
GSE140802_ROOT <- Sys.getenv(
  "GSE140802_ROOT",
  "F:/projects/TIPS/source/GSE140802_lineage_tracking"
)
TAG <- "4_9vs11"
TD.cut <- 0.7
td_csv_dir <- Sys.getenv(
  "MUTRANS_TD_CSV_DIR",
  here::here("examples", "hematopoietic_LARRY", "data", "larry_figures")
)
pub_xlsx <- Sys.getenv(
  "WEINREB_TABLE_S3",
  "F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/aaw3381-Weinreb-Table-S3.xlsx"
)
pub_lineage <- "Megakaryocyte"
pub_sheet <- "DGE of progenitors in vitro"

## PPI arms to test (same as 8.2 venn3 runs)
ppi_arms <- c("STRING", "IID")

## Background modes written to output TSV
background_modes <- c("hvg3000_union_s3")  # or also "hvg3000"

## fisher.test alternative: "two.sided" (R default) or "greater" (enrichment one-sided)
fisher_alternative <- Sys.getenv("TIPS_FET_ALTERNATIVE", "two.sided")
########## END USER INPUT ##########

norm_genes <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  toupper(x)
}

filter_to_universe <- function(genes, universe) {
  intersect(norm_genes(genes), norm_genes(universe))
}

## TIPS CF signature genes (Panel G / venn3 / FET):
##   Source file per arm:
##     .../cisTarget_predicted_4/PPI_graph_GRN_prediction_CTS_4_dualpull_final_table.tsv
##   IID: all unique from/to on linkage == "CF" (n=8; Ikzf2, Sptb + 6 shared).
##   STRING: CF edges from BHLHE41/GATA1/GATA2 dual-pull only (n=7).
##     The combined dualpull table also contains MXD1 cisTarget sub-run edges
##     (rows where from or to is MXD1); those are NOT part of the STRING CF
##     signature module (Hbb-bs + 6 shared) — not an arbitrary subtraction.
tips_cf_genes <- function(tips_cf_df, ppi_arm) {
  cf <- tips_cf_df
  if (toupper(ppi_arm) == "STRING") {
    is_mxd1 <- toupper(cf$from) == "MXD1" | toupper(cf$to) == "MXD1"
    cf <- cf[!is_mxd1, , drop = FALSE]
  }
  norm_genes(c(cf$from, cf$to))
}

## 2×2 Fisher's exact test for overlap enrichment.
## Rows = set_a, columns = set_b (matches manual cbind):
##   fisher.test(cbind(c(k, n_b-k), c(n_a-k, N-n_a-n_b+k)))
## where set_a is first label (e.g. TD), set_b is second (e.g. TIPS).
fisher_overlap <- function(
    set_a, set_b, universe, label_a, label_b,
    alternative = fisher_alternative
) {
  U <- unique(norm_genes(universe))
  A <- filter_to_universe(set_a, U)
  B <- filter_to_universe(set_b, U)
  N <- length(U)
  a <- length(A)
  b <- length(B)
  k <- length(intersect(A, B))
  a_only <- a - k
  b_only <- b - k
  neither <- N - a - b + k

  empty_row <- data.frame(
    set_a = label_a, set_b = label_b,
    universe_n = N, n_a = a, n_b = b, overlap_k = k,
    cell_in_both = k, cell_a_only = a_only, cell_b_only = b_only, cell_neither = neither,
    overlap_genes = "",
    p_value = NA_real_, odds_ratio = NA_real_,
    fisher_alternative = alternative,
    stringsAsFactors = FALSE
  )

  if (N == 0L) return(empty_row)

  if (a == 0L || b == 0L) {
    empty_row$overlap_genes <- paste(sort(intersect(A, B)), collapse = ", ")
    return(empty_row)
  }

  tab <- matrix(
    c(k, a_only, b_only, neither),
    nrow = 2L, byrow = TRUE,
    dimnames = list(c("in_A", "not_in_A"), c("in_B", "not_in_B"))
  )
  ft <- stats::fisher.test(tab, alternative = alternative)

  data.frame(
    set_a = label_a,
    set_b = label_b,
    universe_n = N,
    n_a = a,
    n_b = b,
    overlap_k = k,
    cell_in_both = k,
    cell_a_only = a_only,
    cell_b_only = b_only,
    cell_neither = neither,
    fisher_cbind_col1 = paste(c(k, b_only), collapse = ","),
    fisher_cbind_col2 = paste(c(a_only, neither), collapse = ","),
    overlap_genes = paste(sort(intersect(A, B)), collapse = ", "),
    p_value = unname(ft$p.value),
    odds_ratio = unname(ft$estimate),
    fisher_alternative = alternative,
    stringsAsFactors = FALSE
  )
}

pairwise_fet_table <- function(gene_sets, universe) {
  labs <- names(gene_sets)
  rows <- vector("list", length(labs) * (length(labs) - 1L) / 2L)
  k <- 0L
  for (i in seq_len(length(labs) - 1L)) {
    for (j in (i + 1L):length(labs)) {
      k <- k + 1L
      rows[[k]] <- fisher_overlap(
        gene_sets[[i]], gene_sets[[j]], universe,
        labs[i], labs[j]
      )
    }
  }
  do.call(rbind, rows)
}

load_weinreb_table_s3_genes <- function(pub_xlsx, sheet = pub_sheet) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Install readxl: install.packages('readxl')", call. = FALSE)
  }
  if (!file.exists(pub_xlsx)) {
    stop("Missing Weinreb Table S3: ", pub_xlsx, call. = FALSE)
  }
  pub_df <- readxl::read_excel(pub_xlsx, sheet = sheet)
  sym_col <- "Gene symbol"
  if (!sym_col %in% colnames(pub_df)) {
    stop("Expected column '", sym_col, "' in ", sheet, call. = FALSE)
  }
  pub_df <- pub_df[
    !is.na(pub_df[[sym_col]]) & nzchar(as.character(pub_df[[sym_col]])),
    ,
    drop = FALSE
  ]
  norm_genes(pub_df[[sym_col]])
}

load_sce_hvg_genes <- function(code_dir) {
  mutrans_helper <- file.path(code_dir, "mutrans_sce_xy.R")
  if (!file.exists(mutrans_helper)) {
    stop("Missing ", mutrans_helper, " — needed for SCE background", call. = FALSE)
  }
  source(mutrans_helper, local = TRUE)
  if (!exists("load_larry_sce", mode = "function")) {
    stop("mutrans_sce_xy.R did not define load_larry_sce()", call. = FALSE)
  }
  sce <- load_larry_sce()
  hvg <- norm_genes(rownames(sce))
  if (length(hvg) < 500L) {
    warning(
      "[8.21] SCE has only ", length(hvg), " genes — expected ~3000 HVG",
      call. = FALSE
    )
  }
  list(sce = sce, hvg = hvg)
}

build_backgrounds <- function(hvg, s3_all) {
  list(
    hvg3000 = unique(hvg),
    hvg3000_union_s3 = unique(c(hvg, s3_all))
  )
}

load_venn3_gene_sets <- function(ppi_arm) {
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

  pub_df <- readxl::read_excel(pub_xlsx, sheet = pub_sheet)
  pub_df <- pub_df[
    !is.na(pub_df[["Gene symbol"]]) & nzchar(as.character(pub_df[["Gene symbol"]])),
    ,
    drop = FALSE
  ]
  pub_by_lineage <- split(as.character(pub_df[["Gene symbol"]]), pub_df$Lineage)
  pub_by_lineage <- lapply(pub_by_lineage, norm_genes)
  if (!pub_lineage %in% names(pub_by_lineage)) {
    stop("Lineage not in Weinreb table: ", pub_lineage, call. = FALSE)
  }
  pub_mega <- pub_by_lineage[[pub_lineage]]

  td_file <- file.path(
    td_csv_dir,
    paste0("td_genes_scores_", CP_cluster, "_to_", CF_cluster, "_seacell.csv")
  )
  if (!file.exists(td_file)) stop("Missing TD file: ", td_file, call. = FALSE)
  td_df <- utils::read.csv(td_file, stringsAsFactors = FALSE)
  td_genes <- norm_genes(td_df$Genes[abs(td_df$corr) > TD.cut])

  tips_table <- file.path(
    results_dir,
    paste0("cisTarget_predicted_", CTS_ID),
    paste0("PPI_graph_GRN_prediction_CTS_", CTS_ID, "_dualpull_final_table.tsv")
  )
  if (!file.exists(tips_table)) {
    stop("Run 24.1 first — missing: ", tips_table, call. = FALSE)
  }
  tips_df <- utils::read.delim(tips_table, stringsAsFactors = FALSE)
  tips_cf <- tips_df[tips_df$linkage == "CF", c("from", "to"), drop = FALSE]
  tips_genes <- tips_cf_genes(tips_cf, ppi_arm)

  gene_sets <- list(
    `Mega (Weinreb)` = pub_mega,
    `TD (A4->A11)` = td_genes,
    TIPS = tips_genes
  )

  list(
    code_dir = code_dir,
    results_dir = results_dir,
    gene_sets = gene_sets,
    ppi_arm = ppi_arm,
    cp = CP_cluster,
    cf = CF_cluster,
    cts_id = CTS_ID
  )
}

run_fet_arm <- function(ppi_arm, backgrounds) {
  message("[8.2.1] ===== ", ppi_arm, " ", TAG, " =====")
  loaded <- load_venn3_gene_sets(ppi_arm)
  venn_dir <- file.path(loaded$results_dir, "venn")
  dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

  all_rows <- list()
  for (bg_name in intersect(names(backgrounds), background_modes)) {
    universe <- backgrounds[[bg_name]]
    tab <- pairwise_fet_table(loaded$gene_sets, universe)
    tab$ppi_arm <- ppi_arm
    tab$background <- bg_name
    tab$tag <- TAG
    tab$pub_lineage <- pub_lineage
    tab$td_cut <- TD.cut
    all_rows[[bg_name]] <- tab
  }
  out <- do.call(rbind, all_rows)
  rownames(out) <- NULL

  gene_list_file <- file.path(
    venn_dir,
    paste0("FET_venn3_gene_lists_", ppi_arm, "_", TAG, ".tsv")
  )
  gl_rows <- lapply(names(loaded$gene_sets), function(nm) {
    data.frame(
      ppi_arm = loaded$ppi_arm,
      set = nm,
      n = length(loaded$gene_sets[[nm]]),
      genes = paste(sort(loaded$gene_sets[[nm]]), collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  utils::write.table(
    do.call(rbind, gl_rows), gene_list_file,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  message("[8.2.1] wrote gene lists: ", gene_list_file)

  out_file <- file.path(
    venn_dir,
    paste0("FET_venn3_Mega_TD_TIPS_", ppi_arm, "_", TAG, ".tsv")
  )
  utils::write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[8.2.1] wrote ", out_file)

  ## Annotation-friendly one-liners (for figure legends)
  primary <- out[out$background == "hvg3000_union_s3", , drop = FALSE]
  for (i in seq_len(nrow(primary))) {
    message(sprintf(
      "[8.2.1] %s | %s vs %s: k=%d, OR=%.2f, p=%.4g (N=%d; n_a=%d, n_b=%d)",
      ppi_arm,
      primary$set_a[i], primary$set_b[i],
      primary$overlap_k[i], primary$odds_ratio[i],
      primary$p_value[i], primary$universe_n[i],
      primary$n_a[i], primary$n_b[i]
    ))
  }

  invisible(list(table = out, venn_dir = venn_dir, gene_sets = loaded$gene_sets))
}

## --- main ---
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install readxl: install.packages('readxl')", call. = FALSE)
}

## SCE background shared across PPI arms (same Larry metacells / HVG)
ref_code_dir <- file.path(
  GSE140802_ROOT,
  "7_data_MuTrans_TIPS_IID",
  paste0("code_core_", TAG)
)
source(file.path(ref_code_dir, "00_configuration.R"))
tips_configure(TAG = TAG, wd = paste0(
  normalizePath(
    file.path(GSE140802_ROOT, "7_data_MuTrans_TIPS_IID"),
    winslash = "/", mustWork = FALSE
  ),
  "/"
))
sce_info <- load_sce_hvg_genes(ref_code_dir)
s3_all <- load_weinreb_table_s3_genes(pub_xlsx, sheet = pub_sheet)
backgrounds <- build_backgrounds(sce_info$hvg, s3_all)

bg_out_dir <- file.path(
  GSE140802_ROOT,
  "7_data_MuTrans_TIPS_STRING",
  paste0("results_core_", TAG),
  "venn"
)
dir.create(bg_out_dir, recursive = TRUE, showWarnings = FALSE)
bg_union <- backgrounds$hvg3000_union_s3
bg_df <- data.frame(
  gene = sort(bg_union),
  in_hvg3000 = bg_union %in% backgrounds$hvg3000,
  in_table_s3 = bg_union %in% s3_all,
  stringsAsFactors = FALSE
)
utils::write.table(
  bg_df,
  file.path(bg_out_dir, paste0("FET_background_hvg3000_union_s3_", TAG, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message(
  "[8.2.1] fisher_alternative=", fisher_alternative
)
message(
  "[8.2.1] background: HVG n=", length(backgrounds$hvg3000),
  "; Table S3 n=", length(s3_all),
  "; union n=", length(bg_union),
  " (S3-only adds ", sum(!bg_df$in_hvg3000), " genes)"
)

res <- lapply(ppi_arms, run_fet_arm, backgrounds = backgrounds)

## Combined summary across arms / backgrounds
combined <- do.call(rbind, lapply(res, `[[`, "table"))
combined_file <- file.path(
  GSE140802_ROOT,
  "7_data_MuTrans_TIPS_STRING",
  paste0("results_core_", TAG),
  "venn",
  paste0("FET_venn3_Mega_TD_TIPS_ALL_", TAG, ".tsv")
)
dir.create(dirname(combined_file), recursive = TRUE, showWarnings = FALSE)
utils::write.table(combined, combined_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("[8.2.1] wrote combined summary: ", combined_file)

# [tips_configure] TAG=4_9vs11 db_species=10090 (mouse) group_col=attractor CP=4 CM=9 CF=11 motif_target_strategy=solo -> F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_IID/results_core_4_9vs11/
# [mutrans_sce] Seurat: F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/larry_BioTIP/BioTIP_attractor/seu_attractor_MuTrans_HVG.rds
# [assert_sce_gene_overlap] Seurat build vs DEG: 1467 / 1467
# [mutrans_sce] assay=logcounts (Seurat RNA data layer) | metacells: 1200 | attractors: 14 | genes: 3000


# -
                                                                                                                          
# [8.2.1] fisher_alternative=two.sided
# [8.2.1] background: HVG n=3000; Table S3 n=391; union n=3204 (S3-only adds 204 genes)
# [8.2.1] ===== STRING 4_9vs11 =====
# [tips_configure] TAG=4_9vs11 db_species=10090 (mouse) group_col=attractor CP=4 CM=9 CF=11 motif_target_strategy=solo -> F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/results_core_4_9vs11/

                                                                                                                         

# -
                                                                                                                          
# [8.2.1] wrote gene lists: F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/results_core_4_9vs11//venn/FET_venn3_gene_lists_STRING_4_9vs11.tsv
# [8.2.1] wrote F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/results_core_4_9vs11//venn/FET_venn3_Mega_TD_TIPS_STRING_4_9vs11.tsv
# [8.2.1] STRING | Mega (Weinreb) vs TD (A4->A11): k=3, OR=4.46, p=0.03691 (N=3204; n_a=58, n_b=41)
# [8.2.1] STRING | Mega (Weinreb) vs TIPS: k=2, OR=22.33, p=0.006381 (N=3204; n_a=58, n_b=7)
# [8.2.1] STRING | TD (A4->A11) vs TIPS: k=1, OR=13.11, p=0.08629 (N=3204; n_a=41, n_b=7)
# [8.2.1] ===== IID 4_9vs11 =====
# [tips_configure] TAG=4_9vs11 db_species=10090 (mouse) group_col=attractor CP=4 CM=9 CF=11 motif_target_strategy=solo -> F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_IID/results_core_4_9vs11/

                                                                                                                        

# -
                                                                                                                          
# [8.2.1] wrote gene lists: F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_IID/results_core_4_9vs11//venn/FET_venn3_gene_lists_IID_4_9vs11.tsv
# [8.2.1] wrote F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_IID/results_core_4_9vs11//venn/FET_venn3_Mega_TD_TIPS_IID_4_9vs11.tsv
# [8.2.1] IID | Mega (Weinreb) vs TD (A4->A11): k=3, OR=4.46, p=0.03691 (N=3204; n_a=58, n_b=41)
# [8.2.1] IID | Mega (Weinreb) vs TIPS: k=2, OR=18.61, p=0.008409 (N=3204; n_a=58, n_b=8)
# [8.2.1] IID | TD (A4->A11) vs TIPS: k=1, OR=11.24, p=0.098 (N=3204; n_a=41, n_b=8)
# [8.2.1] wrote combined summary: F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/results_core_4_9vs11/venn/FET_venn3_Mega_TD_TIPS_ALL_4_9vs11.tsv
# > 

