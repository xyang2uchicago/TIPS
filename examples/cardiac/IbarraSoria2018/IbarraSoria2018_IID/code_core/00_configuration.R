## 00_configuration.R — cardiac IbarraSoria2018 (E8.5 mouse gastrulation), IID arm
##
## Call tips_configure() from run_TIPS_core.R (or let ensure_tips_configured()
## do it automatically when a numbered script is run standalone).

source(here::here("examples", "tips_core_shared_config.R"))

tips_configure <- function(
    species       = "mouse",
    celltype_col  = "subcelltype",
    CP_cluster    = "cardiac.a",
    CM_cluster    = "cardiac.c",
    CF_cluster    = "extraembryonicMesoderm",
    CTS_ID        = CP_cluster,
    seed_TF       = c("GATA4", "TBX5", "MSX2"),   # from CTS_cardiac.a (code 12.0 output)
    NES_threshold = 3,
    logFC.cut     = 0.6,
    core_count    = 1,      # parallel cores for 11.2.0 steps 1-2; keep 1 on Windows
    celltype_specific_weight_version = "10",
    BioTIP_version = "06232025",
    heatmap_coding_target_only = TRUE,
    wd = here::here("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID/"),
    envir = .GlobalEnv
) {
  wd <- normalizePath(wd, winslash = "/", mustWork = FALSE)
  if (!grepl("/$", wd)) wd <- paste0(wd, "/")
  assign("wd", wd, envir = envir)

  assign("db", "IbarraSoria2018", envir = envir)
  assign("db_species", 10090L, envir = envir)
  assign("species", species, envir = envir)
  assign("celltype_col", celltype_col, envir = envir)
  assign("CP_cluster", CP_cluster, envir = envir)
  assign("CM_cluster", CM_cluster, envir = envir)
  assign("CF_cluster", CF_cluster, envir = envir)
  assign("CTS_ID", CTS_ID, envir = envir)
  assign("CTS_name", paste0("CTS_", CTS_ID), envir = envir)
  assign("seed_TF", seed_TF, envir = envir)
  assign("NES_threshold", NES_threshold, envir = envir)
  assign("logFC.cut", logFC.cut, envir = envir)
  assign("core_count", core_count, envir = envir)
  assign("celltype_specific_weight_version", celltype_specific_weight_version, envir = envir)
  assign("BioTIP_version", BioTIP_version, envir = envir)
  assign("heatmap_coding_target_only", heatmap_coding_target_only, envir = envir)
  assign("specificity_methods", c("combined"), envir = envir)
  assign("s", "combined", envir = envir)
  assign("step1", TRUE, envir = envir)
  assign("step2", TRUE, envir = envir)
  assign("PPI_color_palette", c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02"), envir = envir)
  assign("top_TF_rank", 3L, envir = envir)
  assign("gene_top_n", 20L, envir = envir)

  results_dir <- paste0(wd, "results_core/")
  data_dir    <- paste0(wd, "../data/")   # IbarraSoria2018/data/, shared with the STRING arm
  ppi_path    <- paste0(results_dir, "PPI_weight/")
  assign("results_dir", results_dir, envir = envir)
  assign("data_dir", data_dir, envir = envir)
  assign("ppi_path", ppi_path, envir = envir)
  assign("db_specifc_output_path", ppi_path, envir = envir)
  ## Note: "outdir"/"inputdir" are deliberately NOT set here — they mean a
  ## different directory in different numbered scripts (11.1's own working
  ## dir vs. 11.2.0/11.3's PPI_weight dir). Each script sets its own from
  ## data_dir/ppi_path directly to avoid ambiguity.

  shared_path <- paste0(shared_data_path, "/")
  assign("shared_path", shared_path, envir = envir)
  assign("coding_genes", tips_core_load_coding_genes(shared_data_path, source = "human"), envir = envir)
  assign("CHD", NULL, envir = envir)

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ppi_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  assign("TIPS_CONFIGURED", TRUE, envir = envir)
  message(
    "[tips_configure] IbarraSoria2018_IID species=", species,
    " CP=", CP_cluster, " CM=", CM_cluster, " CF=", CF_cluster,
    " CTS_ID=", CTS_ID, " -> ", results_dir
  )
  invisible(TRUE)
}
