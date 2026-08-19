## run_heldout_pipeline.R
##
##   source("F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor/run_heldout_pipeline.R")
##
## HELDOUT_RUN:
##   step1     (default)  01 only — inspect counts before any prediction
##   all                  01–08, 10, and 11, including 03b (real TIPS per fold; expensive)
##   from_lock            02–08, 10, and 11 using an existing LOCKED file
##   01 .. 08, 03b, 10, 11 a single step

code_dir <- here::here("examples", "hematopoietic_LARRY", "single_cell", "9_Held_out_clone_fate_prediction_cursor")
heldout_code_dir <- code_dir

run <- Sys.getenv("HELDOUT_RUN", "step1")
if (run == "all") {
  Sys.setenv(HELDOUT_AUTOLOCK = "1")
  if (is.na(Sys.getenv("HELDOUT_OVERWRITE_FOLDS", unset = NA))) {
    Sys.setenv(HELDOUT_OVERWRITE_FOLDS = "1")
  }
}

source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_configure()
scripts <- c(
  "01" = "01_construct_clone_dataset.R",
  "02" = "02_lock_criteria_and_folds.R",
  "03" = "03_gene_sets.R",
  "03b" = "03b_build_tips_module_fold.R",
  "04" = "04_heldout_prediction.R",
  "05" = "05_stats_baselines_wells.R",
  "06" = "06_branch_resolvability.R",
  "07" = "07_figures.R",
  "08" = "08_summary_figure.R",
  "10" = "10_clone_bootstrap.R",
  "11" = "11_same_size_null.R"
)

wanted <- switch(
  run,
  step1 = "01",
  all = names(scripts),
  from_lock = names(scripts)[-1],
  {
    if (run %in% names(scripts)) run else stop("Unknown HELDOUT_RUN=", run)
  }
)

session_stamp(cfg, list(HELDOUT_RUN = run, overwrite_folds = cfg$overwrite_fold_modules))
log_file <- file.path(cfg$results_dir, "running_message.txt")
con <- file(log_file, open = "at")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

message("[run] HELDOUT_RUN=", run, " steps=", paste(wanted, collapse = ","))
for (s in wanted) {
  message("\n========== ", scripts[[s]], " ==========")
  sys.source(file.path(code_dir, scripts[[s]]), envir = new.env(parent = globalenv()))
}
message("[run] done -> ", cfg$results_dir)
