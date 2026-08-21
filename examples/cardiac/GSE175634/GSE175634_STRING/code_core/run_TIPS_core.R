## cardiac GSE175634 (human iPSC-CM differentiation) — STRING arm
## Full docs: README.md in this folder
## USER PARAMETERS below → tips_configure() → 11.1 → 11.1.1 → 11.2.0 → 11.3 → 12.0
##   → {24.1, 24.3} for CP, then {24.1, 24.3}_CP.1
##
## CP and CP.1 are NOT interchangeable via a parameter loop: the CP.1-suffixed
## files contain genuinely different logic (no "MANUAL SELECTION" key_TFs
## override) — verified by diffing against the CP files. Each pair is
## dispatched to its own file rather than collapsed into a loop.
##
## Human STRING db (species=9606) is not shipped locally — set download_files
## below to TRUE once to fetch it from stringdb-static.org (~GB scale).

species        <- "human"
celltype_col   <- "leiden_0.5"
CP_cluster     <- "CP"
CM_cluster     <- "5"
CF_cluster     <- "1"
NES_threshold  <- 4.5
download_files <- FALSE

## core_count for calculate_network_specificity/update_network_weights (11.2.0).
## Overridable via TIPS_CORE_COUNT so the same driver can be submitted at
## cores=4 (gpu partition, PI-validated) or cores=1 (bigmem, also PI-validated
## -- their own comment already documents "runs 1 hour when cores=1" as the
## reference runtime) without editing this file per partition.
core_count <- as.integer(Sys.getenv("TIPS_CORE_COUNT", "4"))

code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_STRING", "code_core")
source(file.path(code_dir, "00_configuration.R"))
tips_configure(
  species = species, celltype_col = celltype_col,
  CP_cluster = CP_cluster, CM_cluster = CM_cluster, CF_cluster = CF_cluster,
  NES_threshold = NES_threshold, download_files = download_files, core_count = core_count
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

## Build the full graph_list once — it already covers both CP and CP.1 clusters.
source(file.path(code_dir, "11.1_STRINGweighted_CTS_cardiac_network.R"))
source(file.path(code_dir, "11.1.1_check_vertex_duplication.R"))
source(file.path(code_dir, "11.2.0_update_network_weights_clean_max.R"))
source(file.path(code_dir, "11.3_CTS_cardiac_network_ANND_pagerank.R"))
source(file.path(code_dir, "12.0_rank_by_PageRank_BC.R"))

## CP branch — expected to fully succeed.
setwd(results_dir)
source(file.path(code_dir, "24.1_acat_CTS.cisTarget_dualpull_clean.R"))
source(file.path(code_dir, "24.3_acat_CTS.cisTarget_merge_GRN.R"))

## CP.1 branch — 24.0_..._CP.1.R locally sets CTS_ID="CP.1", seed_TF=character(0)
## (no HiGCTS_CP.1 seed TFs identified by 12.0). Expected/documented to hit
## 24.1's stop("No key TFs found...") at NES=4.5 — tryCatch lets that halt
## skip past cleanly instead of aborting the driver after CP's results are
## already saved.
setwd(results_dir)
tryCatch({
  source(file.path(code_dir, "24.1_acat_CTS.cisTarget_dualpull_clean_CP.1.R"))
  source(file.path(code_dir, "24.3_acat_CTS.cisTarget_merge_GRN_CP.1.R"))
}, error = function(e) {
  message("[run_TIPS_core] CTS_ID=CP.1 stopped: ", conditionMessage(e))
})

message("[run_TIPS_core] done -> ", results_dir)
