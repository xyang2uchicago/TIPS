## hematoendothelial GSE87038 (E8.25 mouse gastrulation) — STRING arm, focal cluster 13
## Full docs: README.md in this folder
## USER PARAMETERS below → tips_configure() → 11.1 … 24.3
##
## NOTE on 11.1.1: this step performs a live BioMart lookup to resolve which
## duplicate STRING vertex to keep for a handful of genes, then filters by a
## fixed set of verified-bad STRING protein IDs (bad_ids) — robust to DEG
## size/order, unlike the row-index approach some sibling folders use. Safe
## to chain here since STRING db version is pinned ("12.0").
##
## NOTE: logFC.cut=0.5 and NES_threshold=3 match the manuscript and the
## committed data — see 00_configuration.R for details.

species       <- "mouse"
celltype_col  <- "label"
CP_cluster    <- "13"
CM_cluster    <- "7"
CF_cluster    <- "6"
CTS_ID        <- "13"
NES_threshold <- 3

code_dir <- here::here("examples", "hematoendothelial", "GSE87038", "GSE87038_STRING", "code_core_13")
source(file.path(code_dir, "00_configuration.R"))
tips_configure(
  species = species, celltype_col = celltype_col,
  CP_cluster = CP_cluster, CM_cluster = CM_cluster, CF_cluster = CF_cluster,
  CTS_ID = CTS_ID, NES_threshold = NES_threshold
)

## Log console output + messages to results folder (after tips_configure)
log_file <- file.path(results_dir, "running_message.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)

.old_message <- getOption("message")
options(message = function(...) {
  m <- if (exists(".makeMessage", mode = "function", where = asNamespace("base"))) {
    get(".makeMessage", asNamespace("base"))(..., domain = NULL)
  } else {
    paste(..., collapse = "")
  }
  cat(m, file = log_con, append = TRUE)
  if (is.function(.old_message)) {
    .old_message(...)
  } else {
    cat(m, file = stderr())
  }
})

on.exit({
  options(message = .old_message)
  sink()
  close(log_con)
}, add = TRUE)
message("[log] writing to ", log_file)

source(file.path(code_dir, "11.1_STRINGweighted_CTS_cardiac_network.R"))
source(file.path(code_dir, "11.1.1_check_vertex_duplication.R"))
source(file.path(code_dir, "11.2.0_update_network_weights_clean_max.R"))
source(file.path(code_dir, "11.3_CTS_cardiac_network_ANND_pagerank.R"))
source(file.path(code_dir, "12.0_rank_by_PageRank_BC.R"))
source(file.path(code_dir, "24.1_acat_CTS.cisTarget_dualpull_clean.R"))
source(file.path(code_dir, "24.3_acat_CTS.cisTarget_merge_GRN.R"))

message("[run_TIPS_core] done -> ", results_dir)
