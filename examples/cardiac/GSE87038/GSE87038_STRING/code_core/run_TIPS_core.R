## cardiac GSE87038 (E8.25 mouse gastrulation) — STRING arm
## Full docs: README.md in this folder
## USER PARAMETERS below → tips_configure() → 11.1 … 24.3
##
## NOTE on 11.1.1: this step performs a live BioMart lookup to resolve which
## duplicate STRING vertex to keep for a handful of genes, then applies
## hardcoded row-index removals that were manually determined by inspecting
## that lookup's output. It is safe to chain here because the STRING db
## version is pinned ("12.0") and the duplicate list is regenerated fresh
## within this same run. If STRING/BioMart results ever drift, the manual
## indices in 11.1.1 will need to be re-derived by hand — see the comments
## inside that file before trusting a rerun's output blindly.

species       <- "mouse"
celltype_col  <- "label"
CP_cluster    <- "8"
CM_cluster    <- "17"
CF_cluster    <- "16"
CTS_ID        <- "8"
seed_TF       <- "ISL1"
NES_threshold <- 3

code_dir <- here::here("examples", "cardiac", "GSE87038", "GSE87038_STRING", "code_core")
source(file.path(code_dir, "00_configuration.R"))
tips_configure(
  species = species, celltype_col = celltype_col,
  CP_cluster = CP_cluster, CM_cluster = CM_cluster, CF_cluster = CF_cluster,
  CTS_ID = CTS_ID, seed_TF = seed_TF, NES_threshold = NES_threshold
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
