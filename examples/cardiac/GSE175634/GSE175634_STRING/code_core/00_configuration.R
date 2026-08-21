## 00_configuration.R — cardiac GSE175634 (human iPSC-CM differentiation), STRING arm
##
## Call tips_configure() from run_TIPS_core.R (or let ensure_tips_configured()
## do it automatically when a numbered script is run standalone).
##
## See GSE175634_IID/code_core/00_configuration.R for the shared dataset-level
## notes (Wilcoxon DEG source, unfiltered CTS loop, CP/CP.1 second arm, the
## CP_cluster="3"->"CP" bug fix). STRING-specific notes:
## - Human STRING db (species=9606) is NOT present locally by default — only
##   mouse (10090) PPIN files ship in examples/Shared_Data/PPIN/. Set
##   download_files=TRUE once to fetch the 9606.protein.*.v12.0.txt.gz files
##   from stringdb-static.org (same pattern as GSE87038_STRING's 11.1).
## - 11.2.0 computes only the "combined" specificity method (step3/step4's
##   ISL1 multi-method comparison plots and the ratio/zscore/diff variants
##   were removed as unused; matches every other code_core folder now).
## - 11.1.1 does row-position-based manual duplicate-vertex fixes (STRING db
##   version-pinned "12.0", same fragile-but-accepted pattern as GSE87038_STRING's
##   11.1.1) — preserved as-is.
## - 11.3's producer and 12.0's consumer both use the "rewring" typo
##   consistently (verified directly), matching the IID arm.

source(here::here("examples", "tips_core_shared_config.R"))

tips_configure <- function(
    species       = "human",
    celltype_col  = "leiden_0.5",
    CP_cluster    = "CP",
    CM_cluster    = "5",
    CF_cluster    = "1",
    CTS_ID        = CP_cluster,
    seed_TF       = "HOXB2",  # from 12.0: top PageRank in HiGCTS_CP with BetweennessCentrality > 0; overwritten when 12.0 runs
    NES_threshold = 4.5,
    core_count    = 4,      # original per-file value; calculate_network_specificity "runs 1 hour" per PI's own comment
    celltype_specific_weight_version = "10",
    heatmap_coding_target_only = TRUE,
    specificity_methods = c("combined"),
    score_threshold = 200,   # STRINGdb interaction confidence cutoff
    download_files  = FALSE, # set TRUE once to fetch human (9606) STRING PPIN files
    wd = here::here("examples", "cardiac", "GSE175634", "GSE175634_STRING/"),
    envir = .GlobalEnv
) {
  wd <- normalizePath(wd, winslash = "/", mustWork = FALSE)
  if (!grepl("/$", wd)) wd <- paste0(wd, "/")
  assign("wd", wd, envir = envir)

  assign("db", "GSE175634", envir = envir)
  assign("db_species", 9606L, envir = envir)
  assign("species", species, envir = envir)
  assign("celltype_col", celltype_col, envir = envir)
  assign("CP_cluster", CP_cluster, envir = envir)
  assign("CM_cluster", CM_cluster, envir = envir)
  assign("CF_cluster", CF_cluster, envir = envir)
  assign("CTS_ID", CTS_ID, envir = envir)
  assign("CTS_name", paste0("CTS_", CTS_ID), envir = envir)
  assign("seed_TF", seed_TF, envir = envir)
  assign("NES_threshold", NES_threshold, envir = envir)
  assign("core_count", core_count, envir = envir)
  assign("celltype_specific_weight_version", celltype_specific_weight_version, envir = envir)
  assign("heatmap_coding_target_only", heatmap_coding_target_only, envir = envir)
  assign("specificity_methods", specificity_methods, envir = envir)
  assign("score_threshold", score_threshold, envir = envir)
  assign("download_files", download_files, envir = envir)
  assign("s", "combined", envir = envir)
  ## TIPS_SKIP_SPECIFICITY=1 skips step1 (calculate_network_specificity, the
  ## ~1hr step) and step2 in 11.2.0, reusing network_specificity_list.rds and
  ## the *_combinedweighted.rds graph if they already exist from a prior run
  ## that got this far (e.g. a resubmit after a later step failed). step2
  ## always reads network_specificity_list.rds from disk regardless of step1,
  ## so setting both FALSE together is safe once step1 has actually saved it.
  skip_specificity <- Sys.getenv("TIPS_SKIP_SPECIFICITY", "FALSE") %in% c("1", "TRUE", "true", "T")
  assign("step1", !skip_specificity, envir = envir)
  assign("step2", !skip_specificity, envir = envir)
  assign("PPI_color_platte", c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02"), envir = envir)
  assign("top_TF_rank", 3L, envir = envir)
  assign("gene_top_n", 20L, envir = envir)

  ## TIPS_RESULTS_TAG lets two concurrent runs (e.g. a gpu-partition run at
  ## core_count=4 and a bigmem-partition run at core_count=1, submitted to
  ## race the HPC queue) write to separate folders instead of clobbering each
  ## other's results_core/ mid-run. Empty by default -- unset, this is a no-op.
  results_tag <- Sys.getenv("TIPS_RESULTS_TAG", "")
  results_dir <- paste0(wd, "results_core", results_tag, "/")
  data_dir    <- paste0(wd, "../data/")   # GSE175634/data/, shared with the IID arm
  ppi_path    <- results_dir              # this dataset keeps PPI outputs directly in results_core/, not a PPI_weight/ subfolder
  assign("results_dir", results_dir, envir = envir)
  assign("data_dir", data_dir, envir = envir)
  assign("ppi_path", ppi_path, envir = envir)
  assign("db_specifc_output_path", ppi_path, envir = envir)

  shared_path <- paste0(shared_data_path, "/")
  assign("shared_path", shared_path, envir = envir)
  assign("coding_genes", tips_core_load_coding_genes(shared_data_path, source = "human"), envir = envir)
  assign("CHD", NULL, envir = envir)

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  assign("TIPS_CONFIGURED", TRUE, envir = envir)
  message(
    "[tips_configure] GSE175634_STRING species=", species,
    " CP=", CP_cluster, " CM=", CM_cluster, " CF=", CF_cluster,
    " CTS_ID=", CTS_ID, " -> ", results_dir
  )
  invisible(TRUE)
}
