## tips_core_shared_config.R — shared engine for cardiac/hematoendothelial
## code_core `tips_configure()` implementations.
##
## Each code_core* folder still has its own 00_configuration.R (mirrors
## hematopoietic_LARRY's per-folder pattern) that sources this file and
## defines a folder-local tips_configure() wrapper whose default argument
## values are that folder's known-good literals. This file holds only the
## logic that's genuinely identical across all folders.
##
## Unlike hematopoietic_LARRY's resolve_species_resources(), cisTarget
## resources here are NOT species-driven: every cardiac/hematoendothelial
## script uppercases gene symbols to human orthologs and queries the human
## hg38/hgnc cisTarget database regardless of species (confirmed across all
## 10 code_core folders). The coding-gene reference list DOES vary, but by
## tree (cardiac always uses the human list; hematoendothelial uses the
## native-species list + toupper()), not by species — so it's an explicit
## per-folder choice via tips_core_load_coding_genes(), never auto-derived.

source(here::here("examples", "config.R"))

ensure_tips_configured <- function(code_dir = get0("code_dir", ifnotfound = "")) {
  if (isTRUE(get0("TIPS_CONFIGURED", ifnotfound = FALSE))) {
    return(invisible(TRUE))
  }
  if (!exists("tips_configure", mode = "function")) {
    if (!nzchar(code_dir)) stop("code_dir required for ensure_tips_configured()", call. = FALSE)
    source(file.path(code_dir, "00_configuration.R"), local = FALSE)
  }
  tips_configure()
}

## Fixed cisTarget resource — always human, regardless of species (see note above).
tips_core_cistarget_db     <- "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
tips_core_cistarget_url    <- paste0(
  "https://resources.aertslab.org/cistarget/databases/homo_sapiens/",
  "hg38/refseq_r80/mc_v10_clust/gene_based/", tips_core_cistarget_db
)
tips_core_motif_annot_type <- "hgnc"

## Loads coding_genes.rds (human) or coding_genes_mouse.rds per the folder's
## own explicit choice — cardiac always uses "human" even for mouse datasets;
## hematoendothelial uses "mouse" + uppercase = TRUE for its mouse datasets.
tips_core_load_coding_genes <- function(shared_path, source = c("human", "mouse"), uppercase = FALSE) {
  source <- match.arg(source)
  file <- if (source == "human") "coding_genes.rds" else "coding_genes_mouse.rds"
  genes <- unique(readRDS(file.path(shared_path, file)))
  if (uppercase) genes <- toupper(genes)
  names(genes) <- NULL
  genes
}

## Downloads/verifies the cisTarget feather DB into <shared_path>/cistarget/.
tips_core_ensure_cistarget_feather <- function(shared_path) {
  feather_path <- file.path(shared_path, "cistarget", "")
  if (!dir.exists(feather_path)) dir.create(feather_path, recursive = TRUE)
  dest <- file.path(feather_path, tips_core_cistarget_db)
  if (!file.exists(dest)) {
    options(timeout = 600)
    download.file(tips_core_cistarget_url, destfile = dest, mode = "wb", method = "libcurl")
  }
  dest
}
