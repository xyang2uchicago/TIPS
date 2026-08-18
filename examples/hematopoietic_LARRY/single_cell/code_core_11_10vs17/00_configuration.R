## 00_configuration.R — Felix TIPS_core paths + Larry MuTrans BioTIP inputs
##
## Call tips_configure() from run_TIPS_core.R (TAG + CP/CM/CF roles).
## Step scripts: source this file, then ensure_tips_configured(code_dir).
##
## Two sources (no other pipelines):
##   1) Felix:  F:/projects/TIPS/results/felix/TIPS/examples/Shared_Data/
##              F:/projects/TIPS/results/felix/TIPS/R/
##              F:/projects/Share/PPIN/STRING/STRINGdb/mouse/full/  (STRING)
##   2) Larry:  .../larry_BioTIP/BioTIP_attractor/  (Seurat + BioTIP CTS from your pipeline)

tips_deg_paths <- function(logFC.cut, fdr.cut, min_prop, data_dir) {
  tag <- paste0("min.prop", min_prop, "_lfc", logFC.cut, "_FDFR", fdr.cut)
  list(
    deg_tag = tag,
    deg_rdata = file.path(data_dir, paste0("DEG_perState_", tag, ".rds")),
    markers_rdata = file.path(data_dir, paste0("markers.up_ttest_", tag, ".rds"))
  )
}

tips_extend_r_path <- function(wd = NULL, code_dir = NULL) {
  env <- Sys.getenv("TIPS_EXTEND_R", "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  if (!is.null(wd) && nzchar(wd)) {
    return(normalizePath(file.path(dirname(wd), "TIPS_extend_v2.R"), winslash = "/", mustWork = FALSE))
  }
  if (!is.null(code_dir) && nzchar(code_dir)) {
    return(normalizePath(
      file.path(dirname(dirname(code_dir)), "TIPS_extend_v2.R"),
      winslash = "/", mustWork = FALSE
    ))
  }
  normalizePath("TIPS_extend_v2.R", winslash = "/", mustWork = FALSE)
}

tips_source_extend <- function(wd = NULL, code_dir = NULL) {
  path <- tips_extend_r_path(wd = wd, code_dir = code_dir)
  if (!file.exists(path)) {
    stop("TIPS_extend_v2.R not found: ", path, call. = FALSE)
  }
  if (!exists(".tips_extend_sourced", envir = .GlobalEnv) ||
      !identical(get(".tips_extend_sourced", envir = .GlobalEnv), path)) {
    source(path, local = FALSE)
    assign(".tips_extend_sourced", path, envir = .GlobalEnv)
  }
  invisible(path)
}

resolve_deg_paths <- function(results_dir, min_prop, data_dir) {
  latest <- file.path(results_dir, "deg_cutoffs_latest.rds")
  if (file.exists(latest)) {
    dc <- readRDS(latest)
    list(
      logFC.cut = dc$logFC.cut,
      fdr.cut = dc$fdr.cut,
      min.prop = if (!is.null(dc$min.prop)) dc$min.prop else min_prop,
      deg_tag = dc$deg_tag,
      deg_rdata = dc$deg_file,
      markers_rdata = dc$markers_file
    )
  } else {
    lfc <- as.numeric(Sys.getenv("DEG_LOGFC_CUT", "0.5"))
    fdr <- as.numeric(Sys.getenv("DEG_FDR_CUT", "0.05"))
    paths <- tips_deg_paths(lfc, fdr, min_prop, data_dir)
    c(list(logFC.cut = lfc, fdr.cut = fdr, min.prop = min_prop), paths)
  }
}

## Keep CTS modules passing BioTIP summary table (Felix: res$CTS.candidate[which(res$significant)]).
filter_cts_by_summary <- function(
    CTS,
    summary_path,
    mci_p_cut = 0.05,
    ic_g_local_p_cut = 0.05,
    ic_s_local_p_cut = 0.05,
    ic_s_delta_p_cut = 0.05,
    require_localmax = TRUE,
    always_keep = character()
) {
  if (!file.exists(summary_path)) {
    stop("Missing CTS summary: ", summary_path, call. = FALSE)
  }
  sm <- read.csv(summary_path, stringsAsFactors = FALSE)
  sm$CTS_ID <- as.character(sm$CTS_ID)
  if ("IC_g_p" %in% names(sm) && !"IC_g_local_p" %in% names(sm)) {
    sm$IC_g_local_p <- sm$IC_g_p
  }
  if ("IC_s_p" %in% names(sm) && !"IC_s_local_p" %in% names(sm)) {
    sm$IC_s_local_p <- sm$IC_s_p
  }
  if (!"IC_g_local_p" %in% names(sm)) {
    stop("CTS summary missing IC_g_local_p (or IC_g_p): ", summary_path, call. = FALSE)
  }
  sig <- !is.na(sm$MCI_P) & sm$MCI_P < mci_p_cut &
    !is.na(sm$IC_g_local_p) & sm$IC_g_local_p < ic_g_local_p_cut
  if ("IC_s_local_p" %in% names(sm)) {
    ic_s_ok <- is.na(sm$IC_s_local_p) | sm$IC_s_local_p < ic_s_local_p_cut
    sig <- sig & ic_s_ok
  }
  if ("IC_s_delta_p" %in% names(sm)) {
    ic_s_delta_ok <- is.na(sm$IC_s_delta_p) | sm$IC_s_delta_p < ic_s_delta_p_cut
    sig <- sig & ic_s_delta_ok
  }
  ## Lastly: BioTIP local maximum (CTS_summary_data.csv column localmax)
  if (require_localmax && "localmax" %in% names(sm)) {
    sig <- sig & tolower(sm$localmax) == "yes"
  }
  sig_ids <- unique(c(sm$CTS_ID[sig], as.character(always_keep)))
  sig_ids <- intersect(sig_ids, names(CTS))
  dropped <- setdiff(names(CTS), sig_ids)
  msg_ic_s <- if ("IC_s_local_p" %in% names(sm)) {
    paste0(" IC_s_local_p<", ic_s_local_p_cut, " (NA passes)")
  } else {
    ""
  }
  msg_ic_s_delta <- if ("IC_s_delta_p" %in% names(sm)) {
    paste0(" IC_s_delta_p<", ic_s_delta_p_cut, " (NA passes)")
  } else {
    ""
  }
  message(
    "[filter_cts_by_summary] ", basename(summary_path),
    " MCI_P<", mci_p_cut, " IC_g_local_p<", ic_g_local_p_cut, msg_ic_s,
    msg_ic_s_delta,
    if (require_localmax) " localmax=yes (last)" else "",
    " -> keep ", length(sig_ids), "/", length(CTS),
    ": ", paste(sig_ids, collapse = ", ")
  )
  if (length(dropped)) {
    message("[filter_cts_by_summary] dropped: ", paste(dropped, collapse = ", "))
  }
  if (!length(sig_ids)) {
    stop("No significant CTS modules after summary filter.", call. = FALSE)
  }
  CTS[sig_ids]
}

## STRING + cisTarget resources keyed by NCBI taxon (Felix: db_species at top of run_TIPS_core.R).
resolve_species_resources <- function(db_species, shared_path) {
  db_species <- as.integer(db_species)
  if (db_species == 10090L) {
    list(
      species = "mouse",
      string_ppin_default = "F:/projects/Share/PPIN/STRING/STRINGdb/mouse/full",
      cistarget_db = "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
      cistarget_dir_default = "F:/projects/Share/TF/cistarget/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based",
      cistarget_dir_env = "CISTARGET_MM10_DIR",
      cistarget_url_base = paste0(
        "https://resources.aertslab.org/cistarget/databases/mus_musculus/",
        "mm10/refseq_r80/mc_v10_clust/gene_based/"
      ),
      coding_genes_rds = file.path(shared_path, "coding_genes_mouse.rds"),
      motif_annot_type = "mgi"
    )
  } else if (db_species == 9606L) {
    list(
      species = "human",
      string_ppin_default = "F:/projects/Share/PPIN/STRING/STRINGdb/human/full",
      cistarget_db = "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
      cistarget_dir_default = file.path(shared_path, "cistarget"),
      cistarget_dir_env = "CISTARGET_HG38_DIR",
      cistarget_url_base = paste0(
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/",
        "hg38/refseq_r80/mc_v10_clust/gene_based/"
      ),
      coding_genes_rds = file.path(shared_path, "coding_genes_human.rds"),
      motif_annot_type = "hgnc"
    )
  } else {
    stop(
      "Unsupported db_species: ", db_species, " (use 10090 mouse or 9606 human)",
      call. = FALSE
    )
  }
}

default_real_names <- function(TAG = "4_9vs11") {
  switch(
    TAG,
    "4_9vs11" = c(CP = "progenitor", CM = "Baso", CF = "Meg"),
    "4_10vs11" = c(CP = "progenitor", CM = "Mast", CF = "Meg"),
    "1_5vs7" = c(CP = "progenitor", CM = "Mono", CF = "Neut"),
    c(CP = "progenitor", CM = "CM", CF = "CF")
  )
}

apply_cts_summary_filter <- function(CTS, always_keep = CTS_ID) {
  latest <- file.path(results_dir, "cts_filter_latest.rds")
  if (file.exists(latest)) {
    cf <- readRDS(latest)
    return(filter_cts_by_summary(
      CTS,
      summary_path = cf$summary_path,
      mci_p_cut = cf$mci_p_cut,
      ic_g_local_p_cut = cf$ic_g_local_p_cut,
      ic_s_local_p_cut = if (!is.null(cf$ic_s_local_p_cut)) cf$ic_s_local_p_cut else cts_ic_s_p_cut,
      ic_s_delta_p_cut = if (!is.null(cf$ic_s_delta_p_cut)) cf$ic_s_delta_p_cut else cts_ic_s_delta_p_cut,
      require_localmax = cf$require_localmax,
      always_keep = always_keep
    ))
  }
  filter_cts_by_summary(
    CTS, biotip_cts_summary_csv,
    always_keep = always_keep
  )
}

## TAG + attractor role assignment (CP/CM/CF) and all TIPS paths.
tips_configure <- function(
    TAG = Sys.getenv("TIPS_CORE_TAG", "4_9vs11"),
    CP_cluster = Sys.getenv("TIPS_CP_CLUSTER", "4"),
    CM_cluster = Sys.getenv("TIPS_CM_CLUSTER", "9"),
    CF_cluster = Sys.getenv("TIPS_CF_CLUSTER", "11"),
    group_col = Sys.getenv("TIPS_GROUP_COL", "attractor"),
    db_species = as.integer(Sys.getenv("TIPS_DB_SPECIES", "10090")),
    Real_names = NULL,
    motif_target_strategy = NULL,
    wd = Sys.getenv(
      "TIPS_WD",
      "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/"
    ),
    envir = .GlobalEnv
) {
  wd <- normalizePath(wd, winslash = "/", mustWork = FALSE)
  if (!grepl("/$", wd)) wd <- paste0(wd, "/")

  TIPS_CORE_TAG <- TAG
  CTS_CORE_ID <- Sys.getenv("TIPS_CTS_ID", CP_cluster)

  assign("TIPS_CORE_TAG", TIPS_CORE_TAG, envir = envir)
  assign("CTS_CORE_ID", CTS_CORE_ID, envir = envir)
  assign("code_dir", paste0(wd, "code_core_", TIPS_CORE_TAG, "/"), envir = envir)
  assign("results_dir", paste0(wd, "results_core_", TIPS_CORE_TAG, "/"), envir = envir)
  assign("data_dir", paste0(wd, "data/"), envir = envir)
  assign("ppi_path", paste0(get("results_dir", envir = envir), "PPI_weight/"), envir = envir)

  felix_examples <- normalizePath(
    Sys.getenv("FELIX_EXAMPLES", "F:/projects/TIPS/results/felix/TIPS/examples"),
    winslash = "/", mustWork = FALSE
  )
  shared_path <- normalizePath(
    Sys.getenv("TIPS_SHARED_DATA", file.path(felix_examples, "Shared_Data")),
    winslash = "/", mustWork = FALSE
  )
  tips_r_dir <- normalizePath(
    Sys.getenv("TIPS_R_DIR", "F:/projects/TIPS/results/felix/TIPS/R"),
    winslash = "/", mustWork = FALSE
  )
  larry_biotip_dir <- normalizePath(
    Sys.getenv(
      "LARRY_BIOTIP_DIR",
      "F:/projects/TIPS/results/GSE140802_lineage_tracking/inVitro_NMtrajectory/larry_BioTIP/BioTIP_attractor"
    ),
    winslash = "/", mustWork = FALSE
  )

  assign("felix_examples", felix_examples, envir = envir)
  assign("shared_path", shared_path, envir = envir)
  assign("tips_r_dir", tips_r_dir, envir = envir)
  assign("larry_biotip_dir", larry_biotip_dir, envir = envir)

  assign("db", "GSE140802", envir = envir)
  assign("db_species", as.integer(db_species), envir = envir)

  sp_res <- resolve_species_resources(db_species, shared_path)

  ## Larry MuTrans attractors (group column = attractor in Seurat meta.data)
  # CP_cluster <- Sys.getenv("TIPS_CP_CLUSTER", "4")   # source
  # CM_cluster <- Sys.getenv("TIPS_CM_CLUSTER", "9")  # Baso
  # CF_cluster <- Sys.getenv("TIPS_CF_CLUSTER", "11")  # Meg
  assign("CP_cluster", as.character(CP_cluster), envir = envir)
  assign("CM_cluster", as.character(CM_cluster), envir = envir)
  assign("CF_cluster", as.character(CF_cluster), envir = envir)
  if (is.null(Real_names)) {
    Real_names <- default_real_names(TAG)
  } else if (is.null(names(Real_names)) && length(Real_names) == 3L) {
    names(Real_names) <- c("CP", "CM", "CF")
  }
  assign("Real_names", Real_names, envir = envir)
  if (is.null(motif_target_strategy) || !nzchar(as.character(motif_target_strategy)[1L])) {
    motif_target_strategy <- Sys.getenv("TIPS_MOTIF_TARGET_STRATEGY", "solo")
  }
  motif_target_strategy <- match.arg(motif_target_strategy, c("solo", "fam", "merge"))
  assign("motif_target_strategy", motif_target_strategy, envir = envir)
  assign("group_col", as.character(group_col), envir = envir)
  assign("tips_sce_assay", "logcounts", envir = envir)
  assign("CTS_ID", Sys.getenv("TIPS_CTS_ID", CP_cluster), envir = envir)
  assign("CTS_name", paste0("CTS_", get("CTS_ID", envir = envir)), envir = envir)
  seed_TF <- strsplit(Sys.getenv("TIPS_SEED_TF", ""), ",", fixed = TRUE)[[1]]
  seed_TF <- seed_TF[nzchar(seed_TF)]
  assign("seed_TF", seed_TF, envir = envir)

  assign("celltype_col", "label", envir = envir)
  assign("cluster_labels", "label", envir = envir)
  assign("species", sp_res$species, envir = envir)
  assign("celltype_specific_weight_version", "10", envir = envir)
  assign("NES_threshold", as.numeric(Sys.getenv("TIPS_NES_THRESHOLD", "3")), envir = envir)
  assign("min_prop", 0.25, envir = envir)

  deg_info <- resolve_deg_paths(
    get("results_dir", envir = envir),
    get("min_prop", envir = envir),
    get("data_dir", envir = envir)
  )
  assign("logFC.cut", deg_info$logFC.cut, envir = envir)
  assign("fdr.cut", deg_info$fdr.cut, envir = envir)
  assign("deg_tag", deg_info$deg_tag, envir = envir)
  assign("deg_rdata", deg_info$deg_rdata, envir = envir)
  assign("markers_rdata", deg_info$markers_rdata, envir = envir)
  assign("deg_file", deg_info$deg_rdata, envir = envir)
  assign("markers_file", deg_info$markers_rdata, envir = envir)

  assign(
    "biotip_cts_rdata",
    Sys.getenv("BIOTIP_CTS_RDATA", file.path(larry_biotip_dir, "CTS_Lib_data.RData")),
    envir = envir
  )
  assign(
    "biotip_cts_summary_csv",
    Sys.getenv("BIOTIP_CTS_SUMMARY_CSV", file.path(larry_biotip_dir, "CTS_summary_data.csv")),
    envir = envir
  )

  assign("cts_mci_p_cut", 0.05, envir = envir)
  assign("cts_ic_g_p_cut", 0.05, envir = envir)
  assign("cts_ic_s_p_cut", 0.05, envir = envir)
  assign("cts_ic_s_delta_p_cut", 0.05, envir = envir)
  assign("cts_require_localmax", TRUE, envir = envir)

  seurat_rds <- {
    env <- Sys.getenv("SEURAT_RDS", "")
    default <- file.path(larry_biotip_dir, "seu_attractor_MuTrans_HVG.rds")
    if (nzchar(env) && file.exists(env)) env else default
  }
  LARRY_SEURAT_RDS <- paste0(
    "F:/projects/TIPS/results/GSE140802_lineage_tracking/",
    "inVitro_NMtrajectory/larry_BioTIP/BioTIP_attractor/",
    "seu_attractor_MuTrans_HVG.rds"
  )
  if (!file.exists(seurat_rds) && file.exists(LARRY_SEURAT_RDS)) {
    seurat_rds <- LARRY_SEURAT_RDS
  }
  assign("seurat_rds", seurat_rds, envir = envir)
  assign("LARRY_SEURAT_RDS", LARRY_SEURAT_RDS, envir = envir)

  string_ppin_dir <- normalizePath(
    Sys.getenv("STRING_PPIN_DIR", sp_res$string_ppin_default),
    winslash = "/", mustWork = FALSE
  )
  if (!dir.exists(string_ppin_dir)) {
    fallback <- file.path(get("data_dir", envir = envir), "PPIN")
    if (dir.exists(fallback)) {
      string_ppin_dir <- fallback
    } else {
      string_ppin_dir <- file.path(felix_examples, "GSE87038/data/PPIN")
    }
  }
  assign("string_ppin_dir", string_ppin_dir, envir = envir)

  cistarget_db <- sp_res$cistarget_db
  cistarget_share_dir <- normalizePath(
    Sys.getenv(sp_res$cistarget_dir_env, sp_res$cistarget_dir_default),
    winslash = "/", mustWork = FALSE
  )
  cistarget_url <- paste0(sp_res$cistarget_url_base, cistarget_db)
  cistarget_feather <- file.path(cistarget_share_dir, cistarget_db)
  assign("cistarget_db", cistarget_db, envir = envir)
  assign("cistarget_share_dir", cistarget_share_dir, envir = envir)
  assign("cistarget_url", cistarget_url, envir = envir)
  assign("cistarget_feather", cistarget_feather, envir = envir)
  assign("motif_annot_type", sp_res$motif_annot_type, envir = envir)
  assign("coding_genes_rds", sp_res$coding_genes_rds, envir = envir)
  assign("coding_genes_mouse_rds", sp_res$coding_genes_rds, envir = envir)
  ## legacy aliases (older scripts referenced mm10_* names)
  assign("cistarget_mm10_db", cistarget_db, envir = envir)
  assign("cistarget_mm10_share_dir", cistarget_share_dir, envir = envir)
  assign("cistarget_mm10_url", cistarget_url, envir = envir)
  assign("cistarget_mm10_feather", cistarget_feather, envir = envir)

  assign("CHD", NULL, envir = envir)
  assign("CM_linkage", "CM", envir = envir)
  assign("CF_linkage", "CF", envir = envir)
  assign("heatmap_coding_target_only", TRUE, envir = envir)

  assign("ppi_tag", "STRING", envir = envir)

  db_val <- get("db", envir = envir)
  results_dir_val <- get("results_dir", envir = envir)
  ppi_path_val <- get("ppi_path", envir = envir)
  assign(
    "graph_notsimp_file",
    file.path(results_dir_val, paste0(db_val, "_STRING_graph_perState_notsimplified.rds")),
    envir = envir
  )
  assign(
    "graph_weighted_file",
    file.path(ppi_path_val, paste0(db_val, "_STRING_graph_perState_simplified_combinedweighted.rds")),
    envir = envir
  )
  assign(
    "correct_edges_file",
    file.path(results_dir_val, "correct_n_edges_HiG_STRING2.14.0.rds"),
    envir = envir
  )
  assign(
    "correct_edges_review_file",
    sub("\\.rds$", "_auto_pick_review.txt", file.path(results_dir_val, "correct_n_edges_HiG_STRING2.14.0.rds")),
    envir = envir
  )

  dir.create(get("results_dir", envir = envir), recursive = TRUE, showWarnings = FALSE)
  dir.create(get("ppi_path", envir = envir), recursive = TRUE, showWarnings = FALSE)
  dir.create(get("data_dir", envir = envir), recursive = TRUE, showWarnings = FALSE)

  assign("TIPS_CONFIGURED", TRUE, envir = envir)
  message(
    "[tips_configure] TAG=", TIPS_CORE_TAG,
    " db_species=", db_species, " (", sp_res$species, ")",
    " group_col=", group_col,
    " CP=", CP_cluster, " CM=", CM_cluster, " CF=", CF_cluster,
    " motif_target_strategy=", motif_target_strategy,
    " -> ", get("results_dir", envir = envir)
  )
  tips_source_extend(wd = wd, code_dir = get("code_dir", envir = envir))
  invisible(TRUE)
}

infer_tips_wd <- function(code_dir = "") {
  env_wd <- Sys.getenv("TIPS_WD", "")
  if (nzchar(env_wd)) return(env_wd)
  if (nzchar(code_dir)) {
    if (grepl("MuTrans_TIPS_IID", code_dir, fixed = TRUE)) {
      return("F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_IID/")
    }
    if (grepl("MuTrans_TIPS_STRING", code_dir, fixed = TRUE)) {
      return("F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/")
    }
    parent <- sub("code_core_.*/$", "", code_dir)
    if (nzchar(parent)) return(parent)
  }
  NULL
}

ensure_tips_configured <- function(code_dir = get0("code_dir", ifnotfound = "")) {
  needs_configure <- !isTRUE(get0("TIPS_CONFIGURED", ifnotfound = FALSE)) ||
    !exists("coding_genes_rds", envir = .GlobalEnv, inherits = FALSE) ||
    !exists("Real_names", envir = .GlobalEnv, inherits = FALSE)
  if (!needs_configure) {
    return(invisible(TRUE))
  }
  if (!exists("tips_configure", mode = "function")) {
    if (!nzchar(code_dir)) {
      stop("code_dir required for ensure_tips_configured()", call. = FALSE)
    }
    source(file.path(code_dir, "00_configuration.R"), local = FALSE)
  }
  tag <- get0("TIPS_CORE_TAG", ifnotfound = "")
  if (!nzchar(tag) && nzchar(code_dir)) {
    tag <- sub("^.*code_core_", "", sub("/$", "", code_dir))
  }
  if (!nzchar(tag)) tag <- Sys.getenv("TIPS_CORE_TAG", "4_9vs11")
  default_cp <- switch(tag, "4_9vs11" = "4", "4_10vs11" = "4", "1_5vs7" = "1", "4")
  default_cm <- switch(tag, "4_9vs11" = "9", "4_10vs11" = "10", "1_5vs7" = "5", "9")
  default_cf <- switch(tag, "4_9vs11" = "11", "4_10vs11" = "11", "1_5vs7" = "7", "11")
  db_sp <- get0("db_species", ifnotfound = as.integer(Sys.getenv("TIPS_DB_SPECIES", "10090")))
  wd_arg <- infer_tips_wd(code_dir)
  default_wd <- if (!is.null(wd_arg)) wd_arg else Sys.getenv(
    "TIPS_WD",
    "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/"
  )
  tips_configure(
    TAG = tag,
    CP_cluster = get0("CP_cluster", ifnotfound = Sys.getenv("TIPS_CP_CLUSTER", default_cp)),
    CM_cluster = get0("CM_cluster", ifnotfound = Sys.getenv("TIPS_CM_CLUSTER", default_cm)),
    CF_cluster = get0("CF_cluster", ifnotfound = Sys.getenv("TIPS_CF_CLUSTER", default_cf)),
    group_col = get0("group_col", ifnotfound = Sys.getenv("TIPS_GROUP_COL", "attractor")),
    db_species = db_sp,
    Real_names = get0("Real_names", ifnotfound = NULL),
    motif_target_strategy = get0("motif_target_strategy", ifnotfound = NULL),
    wd = default_wd
  )
}

## Load SCE in memory from Larry Seurat (sources mutrans_sce_xy.R in this code_dir)
load_larry_sce <- function() {
  .code_dir <- if (exists("code_dir", inherits = TRUE)) {
    get("code_dir", inherits = TRUE)
  } else {
    "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/code_core_4_9vs11"
  }
  source(file.path(.code_dir, "mutrans_sce_xy.R"), local = FALSE)

  .seurat <- if (exists("seurat_rds", inherits = TRUE)) {
    get("seurat_rds", inherits = TRUE)
  } else {
    LARRY_SEURAT_RDS
  }
  if (!file.exists(.seurat)) .seurat <- LARRY_SEURAT_RDS

  .sce_cache <- if (exists("sce_rdata", inherits = TRUE)) {
    get("sce_rdata", inherits = TRUE)
  } else {
    NULL
  }
  use_cache <- isTRUE(as.logical(Sys.getenv("USE_SCE_CACHE", "FALSE")))
  sce_obj <- NULL
  if (use_cache && !is.null(.sce_cache) && file.exists(.sce_cache)) {
    message("[load_larry_sce] cache: ", .sce_cache)
    sce_try <- readRDS(.sce_cache)
    assert_tips_sce_assay(sce_try)
    cache_ok <- tryCatch({
      validate_tips_sce(sce_try, context = "cache")
      TRUE
    }, error = function(e) {
      message("[load_larry_sce] ignoring bad cache: ", conditionMessage(e))
      FALSE
    })
    if (isTRUE(cache_ok)) {
      sce_obj <- sce_try
    }
  }
  if (is.null(sce_obj)) {
    sce_obj <- build_mutrans_sce(
      save_path = if (use_cache) .sce_cache else NULL,
      seurat_path = .seurat
    )
  }
  sync_sce_cluster_labels(sce_obj, group_col = group_col)
}
