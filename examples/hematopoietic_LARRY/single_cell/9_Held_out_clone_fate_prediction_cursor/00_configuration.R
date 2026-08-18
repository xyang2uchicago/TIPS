## 00_configuration.R — held-out clone-fate prediction (LARRY / GSE140802)
##
## Design: F:/projects/TIPS/doc/NatureComm/Larry_prediction.docx
## This is held-out prediction of later clonal fate, not prospective prediction.
##
## C11 = leiden_r0_8 cluster 11 in the cell-level NM-trajectory object (transition-
## prone progenitor). It is NOT MuTrans attractor 11 (that attractor is Meg).
## Later Meg/Baso fate uses Weinreb Cell.type.annotation, never cluster 17 / 10.

heldout_configure <- function(envir = .GlobalEnv) {
  cfg <- list(
    code_dir = normalizePath(
      "F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor",
      winslash = "/", mustWork = FALSE
    ),
    results_dir = normalizePath(
      Sys.getenv(
        "HELDOUT_RESULTS_DIR",
        "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/9_heldout_clone_fate"
      ),
      winslash = "/", mustWork = FALSE
    ),
    public_data_dir = normalizePath(
      Sys.getenv(
        "LARRY_PUBLIC_DIR",
        "F:/projects/TIPS/data/GSE140802_lineage_tracking/Github"
      ),
      winslash = "/", mustWork = FALSE
    ),
    seurat_rds = normalizePath(
      Sys.getenv(
        "HELDOUT_SEURAT_RDS",
        "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/Reanalyzed_NMtrajectory_Seurat5_noCellCycle.PCA_UMAP.rds"
      ),
      winslash = "/", mustWork = FALSE
    ),
    obs_metadata_csv = normalizePath(
      Sys.getenv(
        "HELDOUT_OBS_CSV",
        "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/h5ad_export/obs_metadata.csv"
      ),
      winslash = "/", mustWork = FALSE
    ),
    weinreb_xlsx = normalizePath(
      Sys.getenv(
        "WEINREB_TABLE_S3",
        "F:/projects/TIPS/data/GSE140802_lineage_tracking/doc/aaw3381-Weinreb-Table-S3.xlsx"
      ),
      winslash = "/", mustWork = FALSE
    ),
    mutrans_td_dir = normalizePath(
      Sys.getenv(
        "MUTRANS_TD_CSV_DIR",
        "F:/projects/TIPS/results/GSE140802_lineage_tracking/larry/figures"
      ),
      winslash = "/", mustWork = FALSE
    ),
    cell_tips_root = normalizePath(
      Sys.getenv(
        "TINGJUN_TIPS_ROOT",
        "F:/projects/TIPS/results/GSE140802_lineage_tracking/Tingjunl/TIPS/7_scaledata_leiden_r0_8_TIPS_STRING"
      ),
      winslash = "/", mustWork = FALSE
    ),
    cell_biotip_dir = normalizePath(
      Sys.getenv(
        "TINGJUN_BIOTIP_C11",
        "F:/projects/TIPS/results/GSE140802_lineage_tracking/Tingjunl/BioTIP/results/BioTIP_leiden_r0_8"
      ),
      winslash = "/", mustWork = FALSE
    ),
    metacell_tips_string = normalizePath(
      "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/results_core_4_9vs11",
      winslash = "/", mustWork = FALSE
    ),
    metacell_tips_iid = normalizePath(
      "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_IID/results_core_4_9vs11",
      winslash = "/", mustWork = FALSE
    ),
    string_iid_jaccard = normalizePath(
      "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/larry_BioTIP/jaccard_STRING_vs_IID_4_9vs11.tsv",
      winslash = "/", mustWork = FALSE
    ),
    hvg_rds = normalizePath(
      "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/cell_cycle_signature_cor_3khvg.rds",
      winslash = "/", mustWork = FALSE
    ),
    ## Frozen C11 state (fate-blind leiden). Sensitivity: set HELDOUT_REMAP_STATE=1.
    state_col = Sys.getenv("HELDOUT_STATE_COL", "leiden_r0_8"),
    state_id = Sys.getenv("HELDOUT_STATE_ID", "11"),
    cf_cluster = Sys.getenv("HELDOUT_CF_CLUSTER", "17"),  # cell-level Meg cluster (module construction only)
    cm_cluster = Sys.getenv("HELDOUT_CM_CLUSTER", "10"),  # cell-level Baso cluster
    day_early = 2L,
    day_late = 6L,
    meg_label = Sys.getenv("HELDOUT_MEG_LABEL", "Meg"),
    baso_label = Sys.getenv("HELDOUT_BASO_LABEL", "Baso"),
    n_folds = as.integer(Sys.getenv("HELDOUT_N_FOLDS", "5")),
    n_repeats = as.integer(Sys.getenv("HELDOUT_N_REPEATS", "1")),
    n_random = as.integer(Sys.getenv("HELDOUT_N_RANDOM", "100")),
    n_bootstrap = as.integer(Sys.getenv("HELDOUT_N_BOOT", "1000")),
    n_perm = as.integer(Sys.getenv("HELDOUT_N_PERM", "1000")),
    seed = as.integer(Sys.getenv("HELDOUT_SEED", "140802")),
    td_cut = as.numeric(Sys.getenv("HELDOUT_TD_CUT", "0.7")),
    ## Primary TIPS construction: re-invoke code_core_11_10vs17 (11.1–24.3)
    ## on training-clone cells. Not a Pearson-delta proxy.
    ##   run_tips_pipeline — actual algorithm (default; report this)
    ##   proxy             — Pearson-delta on existing dual-pull edges (smoke test)
    ##   frozen_genes      — fate-blind TIPS gene list (leakage diagnostic)
    tips_mode = Sys.getenv("HELDOUT_TIPS_MODE", "run_tips_pipeline"),
    score_fun = Sys.getenv("HELDOUT_SCORE_FUN", "mean"),
    min_meg_pos_per_fold = as.integer(Sys.getenv("HELDOUT_MIN_MEG_POS_PER_FOLD", "8")),
    tips_arm_wd = normalizePath(
      Sys.getenv(
        "TIPS_ARM_WD",
        "F:/projects/TIPS/results/GSE140802_lineage_tracking/Tingjunl/TIPS/7_scaledata_leiden_r0_8_TIPS_STRING"
      ),
      winslash = "/", mustWork = FALSE
    ),
    tips_arm_tag = "11_10vs17",
    seurat_rds_tips_arm = normalizePath(
      Sys.getenv(
        "HELDOUT_TIPS_SEURAT_RDS",
        "F:/projects/TIPS/results/GSE140802_lineage_tracking/Tingjunl/BioTIP/results/seu_HVG3000_leiden_r0_8.rds"
      ),
      winslash = "/", mustWork = FALSE
    ),
    ## HELDOUT_RUN=all sets overwrite unless the user exported 0.
    overwrite_fold_modules = Sys.getenv("HELDOUT_OVERWRITE_FOLDS", "0") %in% c("1", "true", "TRUE")
  )

  cfg$public_meta <- file.path(cfg$public_data_dir, "stateFate_inVitro_metadata.txt.gz")
  cfg$public_clone <- file.path(cfg$public_data_dir, "stateFate_inVitro_clone_matrix.mtx.gz")
  cfg$cell_tips_table <- file.path(
    cfg$cell_tips_root, "results_core_11_10vs17", "cisTarget_predicted_11",
    "PPI_graph_GRN_prediction_CTS_11_dualpull_final_table.tsv"
  )
  cfg$cell_tips_graph <- file.path(
    cfg$cell_tips_root, "results_core_11_10vs17",
    "GSE140802_STRING_graph_perState_notsimplified.rds"
  )
  cfg$cell_tips_graph_weighted <- file.path(
    cfg$cell_tips_root, "results_core_11_10vs17", "PPI_weight",
    "GSE140802_STRING_graph_perState_simplified_combinedweighted.rds"
  )
  cfg$cell_cts_rdata <- file.path(
    cfg$cell_biotip_dir, "BioTIP_leiden_r0_8CTS_Lib_Scaledata.RData"
  )
  cfg$cell_cts_testres <- file.path(
    cfg$cell_biotip_dir, "BioTIP_leiden_r0_8optimized_test_sd_selection_Scaledata.RData"
  )
  cfg$cell_cts_summary <- file.path(
    cfg$cell_biotip_dir, "BioTIP_leiden_r0_8CTS_summary_Scaledata.csv"
  )
  cfg$lock_file <- file.path(cfg$results_dir, "02_clone_size_criteria_LOCKED.tsv")
  cfg$propose_file <- file.path(cfg$results_dir, "02_clone_size_criteria_PROPOSED.tsv")
  cfg$tips_arm_code_dir <- file.path(cfg$tips_arm_wd, paste0("code_core_", cfg$tips_arm_tag))
  cfg$tips_arm_results_dir <- file.path(cfg$tips_arm_wd, paste0("results_core_", cfg$tips_arm_tag))
  cfg$fold_modules_dir <- file.path(cfg$results_dir, "03_fold_modules")
  cfg$string_id_cache <- file.path(cfg$tips_arm_wd, "data", "unique_STRING_mapping_correction.txt")
  if (!file.exists(cfg$string_id_cache)) {
    cfg$string_id_cache <- "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/data/unique_STRING_mapping_correction.txt"
  }

  dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$results_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$results_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$results_dir, "rds"), recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$fold_modules_dir, recursive = TRUE, showWarnings = FALSE)

  assign("heldout_cfg", cfg, envir = envir)
  message("[heldout] results -> ", cfg$results_dir)
  invisible(cfg)
}

heldout_ensure <- function() {
  if (!exists("heldout_cfg", envir = .GlobalEnv)) heldout_configure()
  heldout_cfg
}
