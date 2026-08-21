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
##
## Notable dataset-specific facts (see individual numbered scripts for detail):
## - DEG source uses Wilcoxon scores/pvals_adj (score>40, padj<0.01), not the
##   t-test summary.logFC/FDR convention other datasets use. The score/padj
##   thresholds are already baked into the source .RData/.rds filenames
##   (DEG.wc_padj0.01_score40.*), so 11.1's own subset() re-filter is
##   redundant/defensive, not a separately tunable parameter.
## - The CTS loop intentionally uses the UNFILTERED DEG.wc_nofiltering.rds
##   (loaded as markers.up_all) — this is NOT the "stale intern artifact"
##   pattern fixed elsewhere this session; here it's the deliberate,
##   consistent "CTS gets no significance filter" convention, just sourced
##   from a genuinely separate (unfiltered) file rather than the same table.
## - GSE175634 has a second CTS "arm" (CP vs CP.1) inside a single dataset,
##   unlike other folders' IID/STRING split. run_TIPS_core.R loops CTS_ID
##   (and CP_cluster, which must track it) over c("CP","CP.1") instead of
##   having separate _CP.1-suffixed files.
## - Confirmed bug, fixed here: 24.0 hardcoded CP_cluster="3", but the DEG
##   list's cluster "3" was already renamed to "CP" upstream — DEG[["3"]] is
##   NULL, so CP_hi was always computed as all-zero (verified directly
##   against the committed binary_annot_CTS_CP_DEG.tsv, which shows CP_hi
##   summing to 0). Fixed to CP_cluster="CP", matching 12.0's already-correct
##   usage of the same variable.
## - BioTIP here is used as an ordinary installed package (library(BioTIP)) —
##   unlike other datasets, there is no source(BioTIP_update_*.R) call. The
##   celltype_specific_weight_v10.R script (sourced from GitHub, same as
##   everywhere else) supplies calculate_network_specificity/
##   update_network_weights directly; the installed BioTIP package is not
##   the source of those two functions and is not required to be.

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
    core_count    = 4,      # original per-file value; calculate_network_specificity "runs 1 hour" per PI's own comment
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
