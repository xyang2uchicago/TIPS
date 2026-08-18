## 02_lock_criteria_and_folds.R
## Lock clone-size criteria from size distributions only (never from predictive
## performance), then assign clone-grouped CV folds.
##
## Auto-lock rule (size-only; chosen after 01 histograms, never from performance):
##   min_n_d2 = 1
##   min_n_d6 = 1
##     Keep every C11 clone with any annotated day-6 progeny. The primary
##     quasibinomial already weights by n_d6, which is why the design asked
##     to retain numerator and denominator. A min_n_d6 of 5 would drop most
##     Meg-positive C11 clones in this dataset (see 01).
##   min_n_d6_per_well = 1
##   meg_positive_rule = n_meg_d6 >= 1
## Sensitivity thresholds min_n_d6 in {3,5,10} are stored, not used for the primary set.
## If Meg-positive clones are few, 02 stores repeated grouped folds.

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor"
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_ensure()

ds_path <- file.path(cfg$results_dir, "rds", "01_clone_dataset.rds")
if (!file.exists(ds_path)) stop("Run 01 first: ", ds_path)
ds <- readRDS(ds_path)
clones <- ds$clones

n_d6_pos <- clones$n_d6[clones$n_d6 >= 1]
q25 <- if (length(n_d6_pos)) as.numeric(stats::quantile(n_d6_pos, 0.25, names = FALSE)) else 1
n_meg_any <- sum(clones$n_meg_d6 > 0 & clones$n_d6 >= 1)
n_meg_ge5 <- sum(clones$n_meg_d6 > 0 & clones$n_d6 >= 5)
proposed <- data.frame(
  parameter = c("min_n_d2", "min_n_d6", "min_n_d6_per_well", "meg_positive_min_n_meg", "primary_label"),
  value = c("1", "1", "1", "1", "keep_all_d6_binomial_weight"),
  rationale = c(
    "Most C11 clones have a single day-2 cell (01 median n_d2_c11=1)",
    paste0(
      "Keep n_d6>=1. Quasibinomial weights by n_d6 (design: retain counts, not only fractions). ",
      "Meg-positive C11 clones: ", n_meg_any, " at n_d6>=1 vs ", n_meg_ge5,
      " at n_d6>=5 (q25 n_d6 among d6+ clones=", q25, "). Do not drop the rare Meg clones."
    ),
    "Well-specific tests use n_d6_well>=1; concordance subset is clones in both wells",
    "Meg-biased clone for AUROC/AUPRC: at least one annotated day-6 Meg. Locked before scoring.",
    "Do not retune this label after seeing predictive performance"
  ),
  stringsAsFactors = FALSE
)
write_tsv(proposed, cfg$propose_file)

autolock <- Sys.getenv("HELDOUT_AUTOLOCK", "0") %in% c("1", "true", "TRUE")
if (!file.exists(cfg$lock_file)) {
  if (autolock || Sys.getenv("HELDOUT_RUN", "") %in% c("all", "from_lock")) {
    file.copy(cfg$propose_file, cfg$lock_file, overwrite = FALSE)
    message("[02] wrote LOCKED criteria from size-only proposed rule")
  } else {
    stop(
      "No locked criteria yet.\n",
      "Inspect ", cfg$propose_file, " and 01 histograms, then either:\n",
      "  - copy it to ", cfg$lock_file, " (edit thresholds if you wish), or\n",
      "  - Sys.setenv(HELDOUT_AUTOLOCK='1') to accept the size-only rule.\n",
      "Do not choose thresholds by looking at predictive performance."
    )
  }
}

lock <- utils::read.delim(cfg$lock_file, stringsAsFactors = FALSE)
get_lock <- function(name, default) {
  v <- lock$value[lock$parameter == name]
  if (!length(v) || is.na(v[1])) return(default)
  v[1]
}
min_n_d2 <- as.integer(get_lock("min_n_d2", 1))
min_n_d6 <- as.integer(get_lock("min_n_d6", 5))
min_n_d6_well <- as.integer(get_lock("min_n_d6_per_well", 3))
meg_pos_min <- as.integer(get_lock("meg_positive_min_n_meg", 1))

clones$pass_primary <- clones$n_d2_c11 >= min_n_d2 & clones$n_d6 >= min_n_d6
clones$pass_well <- clones$n_d2_c11 >= min_n_d2 &
  clones$n_d6_well1 >= min_n_d6_well & clones$n_d6_well2 >= min_n_d6_well
clones$meg_positive <- clones$n_meg_d6 >= meg_pos_min
clones$baso_positive <- clones$n_baso_d6 >= 1L

primary <- clones[clones$pass_primary, , drop = FALSE]
if (!nrow(primary)) stop("No clones pass locked size criteria")

n_meg_pos <- sum(primary$meg_positive)
message("[02] primary clones=", nrow(primary),
        " Meg-positive=", n_meg_pos,
        " Baso-positive=", sum(primary$baso_positive),
        " both-wells=", sum(clones$pass_well))
n_repeats <- cfg$n_repeats
## Real TIPS per fold is expensive; do not auto-expand to 10 repeats.
if (!identical(cfg$tips_mode, "run_tips_pipeline") &&
    n_meg_pos < cfg$n_folds * cfg$min_meg_pos_per_fold) {
  n_repeats <- max(n_repeats, 10L)
  message("[02] only ", n_meg_pos,
          " Meg-positive clones — using repeated grouped CV (n_repeats=", n_repeats, ")")
} else if (identical(cfg$tips_mode, "run_tips_pipeline") && n_repeats > 1L) {
  message("[02] run_tips_pipeline: n_repeats=", n_repeats,
          " (each repeat is a full 11.1–24.3 run per fold)")
}

fold_list <- lapply(seq_len(n_repeats), function(r) {
  stratified_clone_folds(
    primary,
    n_folds = cfg$n_folds,
    seed = cfg$seed + 100L * r,
    min_meg_pos_per_fold = 1L
  )$fold
})
primary$fold <- fold_list[[1]]
fold_mode <- paste0("repeated_grouped_", cfg$n_folds, "fold_x", n_repeats)
meg_per_fold <- as.integer(table(factor(primary$fold, levels = seq_len(cfg$n_folds)), primary$meg_positive)[, "TRUE"])
if (!length(meg_per_fold)) {
  meg_per_fold <- as.integer(table(factor(primary$fold, levels = seq_len(cfg$n_folds)), as.integer(primary$meg_positive))[, "1"])
}
message("[02] fold mode=", fold_mode, " Meg-pos in repeat1 folds=", paste(meg_per_fold, collapse = ","))

## sensitivity membership (size only)
sens <- data.frame(
  min_n_d6 = c(1L, 3L, 5L, 10L),
  n_clones = NA_integer_,
  n_meg_pos = NA_integer_,
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(sens))) {
  keep <- clones$n_d2_c11 >= min_n_d2 & clones$n_d6 >= sens$min_n_d6[i]
  sens$n_clones[i] <- sum(keep)
  sens$n_meg_pos[i] <- sum(keep & clones$meg_positive)
}

saveRDS(
  list(
    clones = clones, primary = primary, lock = lock, fold_mode = fold_mode,
    sensitivity = sens, fold_list = fold_list, n_repeats = n_repeats, n_folds = cfg$n_folds
  ),
  file.path(cfg$results_dir, "rds", "02_locked_clones_folds.rds")
)
write_tsv(primary, file.path(cfg$results_dir, "tables", "02_primary_clones_with_folds.tsv"))
write_tsv(sens, file.path(cfg$results_dir, "tables", "02_sensitivity_size_thresholds.tsv"))
write_tsv(
  data.frame(fold = seq_along(meg_per_fold), n_meg_positive = meg_per_fold, mode = fold_mode),
  file.path(cfg$results_dir, "tables", "02_fold_balance.tsv")
)

message("[02] locked. Folds assigned. Still no scores.")
