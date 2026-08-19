## 24.0_acat_load_input_clean.R — Felix code_core_13 (Larry MuTrans; cisTarget per db_species)
##
## CODING FIX: mouse toupper block below matches Felix GSE87038 24.0 (not a strategy change).

library(ggplot2)
library(dplyr)
library(purrr)
library(igraph)
library(RcisTarget)
library(data.table)

## Motif annotations — Felix / RcisTarget idiom (see IbarraSoria2018_STRING/code_core_endothelial.b/24.0):
##   data(motifAnnotations_*); motifAnnot <- motifAnnotations   # NOT motifAnnotations_mgi
## Ibarra (mouse STRING): toupper symbols + hg38 + motifAnnotations_hgnc
## Larry MuTrans: db_species 10090 -> mm10 + motifAnnotations_mgi; 9606 -> hg38 + motifAnnotations_hgnc

load_motif_annotations_mgi <- function() {
  for (d in c("motifAnnotations_mgi", "motifAnnotations_mgi_v9")) {
    env <- new.env(parent = emptyenv())
    ok <- tryCatch({
      utils::data(list = d, package = "RcisTarget", envir = env)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) next
    if (exists("motifAnnotations", envir = env, inherits = FALSE)) {
      return(get("motifAnnotations", envir = env))
    }
    if (exists(d, envir = env, inherits = FALSE)) {
      return(get(d, envir = env))
    }
  }
  stop(
    "Could not load mouse motif annotations. Felix pattern:\n",
    "  data(motifAnnotations_mgi); motifAnnot <- motifAnnotations\n",
    "Update RcisTarget: BiocManager::install('RcisTarget')",
    call. = FALSE
  )
}

load_motif_annotations_hgnc <- function() {
  for (d in c("motifAnnotations_hgnc", "motifAnnotations_hgnc_v9")) {
    env <- new.env(parent = emptyenv())
    ok <- tryCatch({
      utils::data(list = d, package = "RcisTarget", envir = env)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) next
    if (exists("motifAnnotations", envir = env, inherits = FALSE)) {
      return(get("motifAnnotations", envir = env))
    }
    if (exists(d, envir = env, inherits = FALSE)) {
      return(get(d, envir = env))
    }
  }
  stop(
    "Could not load human motif annotations. Felix pattern:\n",
    "  data(motifAnnotations_hgnc); motifAnnot <- motifAnnotations\n",
    "Update RcisTarget: BiocManager::install('RcisTarget')",
    call. = FALSE
  )
}

load_motif_annotations <- function(annot_type = NULL) {
  if (is.null(annot_type)) {
    annot_type <- get0("motif_annot_type", ifnotfound = "mgi")
  }
  if (identical(annot_type, "hgnc")) {
    load_motif_annotations_hgnc()
  } else {
    load_motif_annotations_mgi()
  }
}

########## BEGINNING OF USER INPUT ##########
code_dir <- get0("code_dir", ifnotfound = file.path(here::here("examples", "hematopoietic_LARRY", "metacell", "7_data_MuTrans_TIPS_IID"), "code_core_4_9vs11"))
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)
if (!exists("rebuild_mat")) rebuild_mat <- TRUE
########## END OF USER INPUT ##########

if (!file.exists(deg_rdata)) stop("Missing DEG: ", deg_rdata, " — run 11.1 first")
if (!file.exists(graph_weighted_file)) stop("Missing weighted graphs: ", graph_weighted_file, " — run 11.2.0 first")

DEG <- readRDS(deg_rdata)
DEG <- lapply(DEG, as.character)
graph_list <- readRDS(graph_weighted_file)

load(biotip_cts_rdata)
stopifnot(exists("CTS.Lib.Symbol"))
CTS <- apply_cts_summary_filter(CTS.Lib.Symbol, always_keep = CTS_ID)
if (!CTS_ID %in% names(CTS)) stop("CTS_ID not in significant CTS: ", CTS_ID)

coding_rds <- get0("coding_genes_rds", ifnotfound = get0("coding_genes_mouse_rds"))
if (!is.null(coding_rds) && file.exists(coding_rds)) {
  coding_genes <- unique(readRDS(coding_rds))
} else {
  coding_genes <- unique(unlist(CTS, use.names = FALSE))
}

if (requireNamespace("dorothea", quietly = TRUE)) {
  if (identical(species, "human")) {
    data(dorothea_hs, package = "dorothea", envir = environment())
    TF_mouse <- unique(dorothea_hs$tf)
  } else {
    data(dorothea_mm, package = "dorothea", envir = environment())
    TF_mouse <- unique(dorothea_mm$tf)
  }
} else {
  TF_mouse <- character()
}

fileName <- paste0("binary_annot_", CTS_name, "_DEG.tsv")

if (rebuild_mat || !file.exists(fileName)) {
  genes <- as.character(CTS[[CTS_ID]])
  mat <- data.frame(
    CP_hi = as.integer(genes %in% DEG[[CP_cluster]]),
    CM_hi = as.integer(genes %in% DEG[[CM_cluster]]),
    CF_hi = as.integer(genes %in% DEG[[CF_cluster]]),
    row.names = genes, stringsAsFactors = FALSE
  )
  mat[[CTS_name]] <- 1L
  write.table(mat, fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
} else {
  mat <- read.table(fileName, header = TRUE, sep = "\t", check.names = FALSE)
}

if (!file.exists("cisTarget_targets_in_all_CTS.rds")) {
  if (.Platform$OS.type == "windows" && requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(BiocParallel::SerialParam())
  }
  message(
    "[24.0] cisTarget db_species=", db_species, " (", species, ") ",
    motif_annot_type, " -> ", cistarget_feather
  )
  feather_path <- cistarget_share_dir
  if (!dir.exists(feather_path)) dir.create(feather_path, recursive = TRUE)
  if (!file.exists(cistarget_feather)) {
    options(timeout = 600)
    download.file(cistarget_url, destfile = cistarget_feather, mode = "wb", method = "libcurl")
  }
  motifAnnot <- load_motif_annotations()
  motifRankings <- importRankings(cistarget_feather)

  gene_sets_HiG <- DEG[intersect(c(CM_cluster, CF_cluster), names(DEG))]
  cisTarget.res_HiG <- cisTarget(
    geneSets = gene_sets_HiG,
    motifRankings = motifRankings,
    motifAnnot = motifAnnot,
    nesThreshold = 3
  )
  saveRDS(cisTarget.res_HiG, file = "cisTarget_targets_in_two_HiGs.rds")

  cisTarget.res <- cisTarget(
    geneSets = CTS,
    motifRankings = motifRankings,
    motifAnnot = motifAnnot,
    nesThreshold = 3
  )
  saveRDS(cisTarget.res, file = "cisTarget_targets_in_all_CTS.rds")
}

if (file.exists("cisTarget_targets_in_all_CTS.rds")) {
  cisTarget.res <- readRDS("cisTarget_targets_in_all_CTS.rds")
}

if (!exists("motifAnnot")) {
  motifAnnot <- load_motif_annotations()
}

## Felix code_core_13: uppercase mouse symbols before 24.1 dual-pull (identify_TF uses toupper internally)
if (identical(species, "mouse")) {
  for (i in seq_along(graph_list)) {
    g <- graph_list[[i]]
    if (vcount(g) > 0 && !is.null(V(g)$name)) {
      V(g)$name <- toupper(V(g)$name)
      graph_list[[i]] <- g
    }
  }
  DEG <- lapply(DEG, toupper)
  CTS <- lapply(CTS, toupper)
  coding_genes <- toupper(coding_genes)
  TF_mouse <- toupper(TF_mouse)
  if (!is.null(rownames(mat))) rownames(mat) <- toupper(rownames(mat))
  message("[24.0] uppercase mouse symbols done")
}
