## 03_gene_sets.R
## Load comparator gene sets. Weinreb Table S3 is NOT a held-out predictor.
## Training-only DEG sets are built in 04 from training clones.

code_dir <- get0(
  "heldout_code_dir",
  ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/9_Held_out_clone_fate_prediction_cursor"
)
source(file.path(code_dir, "00_configuration.R"))
source(file.path(code_dir, "00_helpers.R"))
cfg <- heldout_ensure()

sets <- list()

## --- cell-level TIPS dual-pull table (fate-blind construction; genes frozen as candidates) ---
if (!file.exists(cfg$cell_tips_table)) stop("Missing cell-level TIPS table: ", cfg$cell_tips_table)
tips_df <- utils::read.delim(cfg$cell_tips_table, stringsAsFactors = FALSE)
sets$tips_meg_frozen <- linkage_set(tips_df, "CF", "CM")
sets$tips_baso_frozen <- linkage_set(tips_df, "CM", "CF")
sets$ppi_unweighted_meg <- norm_genes(with(tips_df[tips_df$linkage == "CF", ], c(from, to)))
sets$ppi_unweighted_baso <- norm_genes(with(tips_df[tips_df$linkage == "CM", ], c(from, to)))
sets$tips_cf_edges <- linkage_edges(tips_df, "CF")
sets$tips_cm_edges <- linkage_edges(tips_df, "CM")

## --- C11 CTS (cell-level BioTIP) ---
if (!file.exists(cfg$cell_cts_rdata)) stop("Missing C11 CTS RData: ", cfg$cell_cts_rdata)
e <- new.env()
load(cfg$cell_cts_rdata, envir = e)
cts_lib <- NULL
for (nm in c("CTS.Lib.Symbol", "CTS.Lib", "CTS")) {
  if (exists(nm, envir = e)) { cts_lib <- get(nm, envir = e); break }
}
if (is.null(cts_lib)) stop("CTS list not found in ", cfg$cell_cts_rdata)
c11_key <- intersect(c(cfg$state_id, paste0("CTS_", cfg$state_id), as.character(cfg$state_id)), names(cts_lib))
if (!length(c11_key)) {
  ## try numeric-name match
  c11_key <- names(cts_lib)[as.character(names(cts_lib)) == cfg$state_id]
}
if (!length(c11_key)) stop("No CTS module for state ", cfg$state_id, "; names: ", paste(names(cts_lib), collapse = ", "))
cts_c11 <- norm_genes(cts_lib[[c11_key[1]]])
if (file.exists(cfg$cell_cts_testres)) {
  o <- new.env()
  load(cfg$cell_cts_testres, envir = o)
  if (exists("testres", envir = o)) {
    uni <- unique(unlist(lapply(get("testres", o), rownames)))
    cts_c11 <- intersect(cts_c11, norm_genes(uni))
  }
}
sets$cts_c11 <- cts_c11

## --- MuTrans transition drivers ---
td_meg <- file.path(cfg$mutrans_td_dir, "td_genes_scores_4_to_11_seacell.csv")
td_baso <- file.path(cfg$mutrans_td_dir, "td_genes_scores_4_to_9_seacell.csv")
read_td <- function(path, cut) {
  if (!file.exists(path)) return(character())
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  gcol <- intersect(c("Genes", "genes", "gene"), names(d))[1]
  ccol <- intersect(c("corr", "Corr"), names(d))[1]
  if (is.na(gcol) || is.na(ccol)) return(character())
  norm_genes(d[[gcol]][abs(d[[ccol]]) > cut])
}
sets$mutrans_td_meg <- read_td(td_meg, cfg$td_cut)
sets$mutrans_td_baso <- read_td(td_baso, cfg$td_cut)

## --- HVG universe for random sets ---
if (file.exists(cfg$hvg_rds)) {
  hvg <- readRDS(cfg$hvg_rds)
  sets$hvg <- norm_genes(if (is.null(names(hvg))) hvg else names(hvg))
} else {
  sets$hvg <- unique(c(sets$cts_c11, sets$mutrans_td_meg))
}

## --- Weinreb Table S3: enrichment / audit only, never a held-out predictor ---
if (file.exists(cfg$weinreb_xlsx) && requireNamespace("readxl", quietly = TRUE)) {
  pub <- readxl::read_excel(cfg$weinreb_xlsx, sheet = "DGE of progenitors in vitro")
  gcol <- intersect(c("Gene symbol", "Gene.symbol"), names(pub))[1]
  lcol <- intersect(c("Lineage"), names(pub))[1]
  sets$weinreb_meg <- norm_genes(pub[[gcol]][pub[[lcol]] == "Megakaryocyte"])
  sets$weinreb_baso <- norm_genes(pub[[gcol]][pub[[lcol]] == "Basophil"])
} else {
  sets$weinreb_meg <- character()
  sets$weinreb_baso <- character()
  message("[03] Weinreb Table S3 not loaded (optional for prediction)")
}

## --- metacell TIPS (concordance, not primary scorer) ---
mc_tab <- file.path(cfg$metacell_tips_string, "cisTarget_predicted_4", "PPI_graph_GRN_prediction_CTS_4_dualpull_final_table.tsv")
if (file.exists(mc_tab)) {
  mdf <- utils::read.delim(mc_tab, stringsAsFactors = FALSE)
  sets$metacell_tips_meg <- linkage_set(mdf, "CF", "CM")
  sets$metacell_tips_baso <- linkage_set(mdf, "CM", "CF")
} else {
  sets$metacell_tips_meg <- character()
  sets$metacell_tips_baso <- character()
}

sizes <- data.frame(
  set = names(sets)[vapply(sets, is.character, logical(1))],
  n = vapply(sets[vapply(sets, is.character, logical(1))], length, integer(1)),
  genes = vapply(sets[vapply(sets, is.character, logical(1))], function(x) paste(sort(x), collapse = ","), character(1)),
  stringsAsFactors = FALSE
)
write_tsv(sizes, file.path(cfg$results_dir, "tables", "03_gene_set_sizes.tsv"))
saveRDS(sets, file.path(cfg$results_dir, "rds", "03_gene_sets.rds"))
message("[03] gene sets:\n", paste(sprintf("  %s n=%d", sizes$set, sizes$n), collapse = "\n"))
message("[03] Weinreb Table S3 is stored for audit only — not used as a held-out comparator.")
