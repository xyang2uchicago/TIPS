## 25_venn_CTS_MuTrans_pub_markers.R — CTS × HiG × MuTrans TD × Weinreb lineage (3×3 Venn grid)
##
## Run after 11.1 (DEG + CTS). Optional: any time for QC / figure.
##   source("25_venn_CTS_MuTrans_pub_markers.R")

########## BEGINNING OF USER INPUT ##########
code_dir <- get0("code_dir", ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/code_core_4_9vs11")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)

TD.cut <- 0.7  # used by MuTrans author,the coefficient to transition path to define candidate TD genes 

## MuTrans TD gene scores (seacell CSV from larry figures)
td_csv_dir <- Sys.getenv(
  "MUTRANS_TD_CSV_DIR",
  "F:/projects/TIPS/results/GSE140802_lineage_tracking/larry/figures"
)

## Weinreb et al. Table S3 — in vitro progenitor DEG by lineage
pub_xlsx <- Sys.getenv(
  "WEINREB_TABLE_S3",
  "F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/aaw3381-Weinreb-Table-S3.xlsx"
)

## Third attractor arm for exploratory Basophil panel (A4→A9); CM/CF from 00_configuration
# third_arm_cluster <- Sys.getenv("TIPS_THIRD_ARM_CLUSTER", "9")
# third_pub_lineage <- "Basophil"
third_arm_cluster <- Sys.getenv("TIPS_THIRD_ARM_CLUSTER", "10")
third_pub_lineage <- "Mast cell"

venn_pdf <- file.path(results_dir, "venn_TD_CTS_HiG.pdf")
########## END OF USER INPUT ##########

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install readxl: install.packages('readxl')", call. = FALSE)
}
if (!requireNamespace("gplots", quietly = TRUE)) {
  stop("Install gplots: install.packages('gplots')", call. = FALSE)
}

library(gplots)

## --- inputs: DEG, CTS, SCE ---
if (!file.exists(deg_rdata)) stop("Run 11.1 first — missing: ", deg_rdata)
DEG <- readRDS(deg_rdata)
DEG <- lapply(DEG, as.character)

load(biotip_cts_rdata)
stopifnot(exists("CTS.Lib.Symbol"))
CTS <- CTS.Lib.Symbol

sce <- load_larry_sce()
N <- ncol(sce)

## --- Weinreb published lineage markers ---
if (!file.exists(pub_xlsx)) stop("Missing Weinreb Table S3: ", pub_xlsx)
pub_df <- readxl::read_excel(pub_xlsx, sheet = "DGE of progenitors in vitro")
pub_genes <- split(as.character(pub_df[["Gene symbol"]]), pub_df$Lineage)
message("[25] Weinreb lineages: ", paste(names(pub_genes), collapse = ", "))

## --- MuTrans TD genes (|corr| > TD.cut) ---
read_td <- function(from_cl, to_cl) {
  f <- file.path(td_csv_dir, paste0("td_genes_scores_", from_cl, "_to_", to_cl, "_seacell.csv"))
  if (!file.exists(f)) stop("Missing TD file: ", f)
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$Genes[which(abs(df$corr) > TD.cut)]
}

TD_CP_CM <- read_td(CP_cluster, CM_cluster)
TD_CP_CF <- read_td(CP_cluster, CF_cluster)
TD_CP_third <- read_td(CP_cluster, third_arm_cluster)

message(
  "[25] TD genes |corr|>", TD.cut, ": ",
  "A", CP_cluster, "->A", CM_cluster, "=", length(TD_CP_CM), "; ",
  "A", CP_cluster, "->A", CF_cluster, "=", length(TD_CP_CF), "; ",
  "A", CP_cluster, "->A", third_arm_cluster, "=", length(TD_CP_third)
)

cts_genes <- as.character(CTS[[CTS_ID]])
deg_cm <- DEG[[CM_cluster]]
deg_cf <- DEG[[CF_cluster]]
deg_third <- DEG[[third_arm_cluster]]
pub_baso <- pub_genes[["Basophil"]]
pub_mega <- pub_genes[["Megakaryocyte"]]
pub_third <- pub_genes[[third_pub_lineage]]

## --- 3×3 Venn grid ---
pdf(venn_pdf, width = 12, height = 12)
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))

## row 1: CTS × HiG × Weinreb
venn(list(CTS = cts_genes, HiG = deg_cm, Baso = pub_baso))
venn(list(CTS = cts_genes, HiG = deg_cf, Mega = pub_mega))
venn(list(CTS = cts_genes, HiG = deg_third, Third = pub_third))

## row 2: TD × HiG × Weinreb
venn(list(TD = TD_CP_CM, HiG = deg_cm, Baso = pub_baso))
venn(list(TD = TD_CP_CF, HiG = deg_cf, Mega = pub_mega))
venn(list(TD = TD_CP_third, HiG = deg_third, Third = pub_third))

## row 3: CTS × TD × Weinreb
venn(list(CTS = cts_genes, TD = TD_CP_CM, Baso = pub_baso))
venn(list(CTS = cts_genes, TD = TD_CP_CF, Mega = pub_mega))
venn(list(CTS = cts_genes, TD = TD_CP_third, Third = pub_third))

dev.off()
message("[25] wrote ", venn_pdf)

## three-way intersections (CTS ∩ TD ∩ Weinreb lineage)
venn_row3_intersect <- list(
  Baso = Reduce(intersect, list(cts_genes, TD_CP_CM, pub_baso)),
  Mega = Reduce(intersect, list(cts_genes, TD_CP_CF, pub_mega)),
  Third = Reduce(intersect, list(cts_genes, TD_CP_third, pub_third))
)
message("[25] row-3 intersections:")
message("  Basophil: ", toString(venn_row3_intersect$Baso))
message("  Mega: ", toString(venn_row3_intersect$Mega))
message("  ", third_pub_lineage, ": ", toString(venn_row3_intersect$Third))

## CTS ∩ TD (per arm)
message("[25] CTS ∩ TD:")
message("  A", CM_cluster, ": ", toString(intersect(cts_genes, TD_CP_CM)))
message("  A", CF_cluster, ": ", toString(intersect(cts_genes, TD_CP_CF)))
message("  A", third_arm_cluster, ": ", toString(intersect(cts_genes, TD_CP_third)))

saveRDS(
  list(
    TD_cut = TD.cut,
    TD = list(CM = TD_CP_CM, CF = TD_CP_CF, third = TD_CP_third),
    venn_row3_intersect = venn_row3_intersect,
    clusters = list(CP = CP_cluster, CM = CM_cluster, CF = CF_cluster, third = third_arm_cluster)
  ),
  file.path(results_dir, "venn_CTS_MuTrans_pub_markers_latest.rds")
)
