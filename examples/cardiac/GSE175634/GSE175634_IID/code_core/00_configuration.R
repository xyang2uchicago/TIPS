## 00_configuration.R — cardiac GSE175634 (human iPSC-CM differentiation), IID arm
##
## Call tips_configure() from run_TIPS_core.R (or let ensure_tips_configured()
## do it automatically when a numbered script is run standalone).
##
## This dataset is PI's own HPC working code (not a repo-native script like
## the other 9 code_core folders): every numbered script hardcoded a Windows
## drive path or a University of Chicago midway3 cluster path, with comments
## documenting a manual "run locally -> scp to cluster -> run there -> scp
## back" workflow. That workflow is collapsed here into a single local
## directory per the established results_core/ convention.

source(here::here("examples", "tips_core_shared_config.R"))

tips_configure <- function(
    species       = "human",
    celltype_col  = "leiden_0.5",
    CP_cluster    = "CP",
    CM_cluster    = "5",
    CF_cluster    = "1",
    CTS_ID        = CP_cluster,
    seed_TF       = c("PRRX1", "HOXB2"),  # from 12.0: top PageRank in HiGCTS_CP with BetweennessCentrality > 0; overwritten when 12.0 runs
    NES_threshold = 4.5,
    core_count    = 1,      # parallel cores for 11.2.0 steps 1-2; default to 1, not everyone has cores to spare
    celltype_specific_weight_version = "10",
    heatmap_coding_target_only = TRUE,
    wd = here::here("examples", "cardiac", "GSE175634", "GSE175634_IID/"),
    envir = .GlobalEnv
) {
  wd <- normalizePath(wd, winslash = "/", mustWork = FALSE)
  if (!grepl("/$", wd)) wd <- paste0(wd, "/")
  assign("wd", wd, envir = envir)

  assign("db", "GSE175634", envir = envir)
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
  assign("s", "combined", envir = envir)
  assign("step1", TRUE, envir = envir)
  assign("step2", TRUE, envir = envir)
  assign("PPI_color_platte", c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02"), envir = envir)
  assign("top_TF_rank", 3L, envir = envir)
  assign("gene_top_n", 20L, envir = envir)

  results_dir <- paste0(wd, "results_core/")
  data_dir    <- paste0(wd, "../data/")   # GSE175634/data/, shared with the STRING arm
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
    "[tips_configure] GSE175634_IID species=", species,
    " CP=", CP_cluster, " CM=", CM_cluster, " CF=", CF_cluster,
    " CTS_ID=", CTS_ID, " -> ", results_dir
  )
  invisible(TRUE)
}
