## hematoendothelial IbarraSoria2018 (E8.5 mouse gastrulation) — IID arm, focal cluster endothelial.b
## Full docs: README.md in this folder (if present)
## USER PARAMETERS below → tips_configure() → 11.1 … 24.3

species       <- "mouse"
celltype_col  <- "subcelltype"
CP_cluster    <- "endothelial.b"
CM_cluster    <- "blood"
CF_cluster    <- "endothelial.c"
CTS_ID        <- "endothelial.b"
seed_TF       <- c("GATA2", "TAL1")
NES_threshold <- 3

code_dir <- here::here("examples", "hematoendothelial", "IbarraSoria2018", "IbarraSoria2018_IID", "code_core_endothelial.b")
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

source(file.path(code_dir, "11.1_IIDweighted_CTS_cardiac_network.R"))
source(file.path(code_dir, "11.2.0_update_network_weights_clean_max.R"))
source(file.path(code_dir, "11.3_CTS_cardiac_network_ANND_pagerank.R"))
source(file.path(code_dir, "12.0_rank_by_PageRank_BC.R"))
source(file.path(code_dir, "24.1_acat_CTS.cisTarget_dualpull_clean.R"))
source(file.path(code_dir, "24.3_acat_CTS.cisTarget_merge_GRN.R"))

message("[run_TIPS_core] done -> ", results_dir)
