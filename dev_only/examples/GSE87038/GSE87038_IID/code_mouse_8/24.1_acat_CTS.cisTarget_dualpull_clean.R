## 24.1_acat_CTS.cisTarget_dualpull_clean.R — code_mouse_8 (IID, C8 arm, full ATAC/ChIP)
##
## Rewritten from GSE87038_IID/code/24.1...R to use the standard shared
## identify_TF_targeted_pull_candidate()/fill_TF_targeting_predicted_edges()
## (R/celltype_specific_weight_v10.R) instead of the original's locally-redefined
## override + a "linkeage_name" typo that doesn't match the real "linkage_name" param
## (would error "unused argument" if actually run). Manual key_TFs unchanged: TCF3,
## HMGA2, KLF6, RARB (IID's own established pick -- differs from STRING's
## TCF3/HMGA2/RARB; IID additionally finds KLF6).
##
## IID's graph_list has no HiGCTS_8 (didn't pass the edge-count threshold in 11.x, which
## is out of scope here) -- the original script's own custom override substituted HiG_8
## in its place. Replicated the same substitution via an alias so the standard shared
## function (which looks up "HiG"+CTS_name = "HiGCTS_8") finds it.

# Set True if running 24.0 the first time
rebuild_mat <- TRUE
source(here::here("examples", "cardiac", "GSE87038", "GSE87038_IID", "code_mouse_8", "24.0_acat_load_input_clean.R"))

source(here::here("R", "celltype_specific_weight_v10.R"))

## check the loaded objects =========================
seed_TF # "NR2F1" "IRX5"  "ALX1"
names(graph_list)
names(DEG)

celltype_col
CP_cluster
CM_cluster
CF_cluster

dim(mat)
colnames(mat)

seed_TF
CTS_ID
CTS_name

## set subfold =======================
(updir <- getwd())
dir.create(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)),
  showWarnings = FALSE, recursive = TRUE
)
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

NES_threshold <- 3

########################################################
##  input 6 -- data-driven --- RcisTarget predicted TF that are enriched among CTS genes
library(RcisTarget)
packageVersion("RcisTarget") # '1.29.0'
library(data.table)

data(motifAnnotations_hgnc)
motifAnnot <- motifAnnotations
dim(motifAnnot)

cisTarget.res <- readRDS(file = "../cisTarget_targets_in_all_CTS.rds")
write.table(subset(cisTarget.res, NES >= NES_threshold & geneSet %in% c(CP_cluster)),
  paste0("../cisTarget_targets_in_", CP_cluster, "_NES", NES_threshold, ".txt"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

########################################################
##   --- binary flag CTS genes

dim(mat)
files <- list.files(pattern = "^heatmap_blocked_CTS_") %>% grep("_v3.tsv", ., value = TRUE)
if (length(files) > 0) {
  fileName <- grep("_v3.tsv", files, value = TRUE)
  mat <- read.table(fileName, sep = "\t", header = T, check.names = FALSE)
  saved_variables <- readRDS(file = "scATAC_cisTarget_variables.rds")
  x <- saved_variables$x
  key_TFs <- saved_variables$key_TFs
  motifAnnot_sub <- saved_variables$motifAnnot_sub
} else {
  motifAnnot_sub <- get_regulators_from_motifs(cisTarget.res, CTS_ID, NES_threshold, motifAnnot = motifAnnot)
  keys <- unique(motifAnnot_sub$regulators)
  if (any(is.na(keys))) keys <- keys[!is.na(keys)]
  keys

  for (key in keys) {
    motif_TF_highConf <- motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf
    tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold & (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf))
    genes <- unique(unlist(strsplit(tmp$enrichedGenes, ";")))
    mat[, paste0("cisTarget_", key, ".motif_target")] <- ifelse(rownames(mat) %in% genes, "1", "0")
  }

  dim(mat)
  colnames(mat)

  ## from the enriched motifs find the potential regulators that are either HiG of CTS member of a cluster of interest
  x <- colnames(mat)[grepl("cisTarget_", colnames(mat))]
  x <- gsub("cisTarget_", "", x) %>% gsub("\\.motif_target", "", .)
  x <- unlist(strsplit(x, ";")) %>% unique()
  x

  ##  extend seed_TF candidates
  key_TFs <- seed_TF
  for (i in x) {
    for (j in c(CM_cluster, CF_cluster, CP_cluster)) {
      if (i %in% DEG[[j]] & i %in% rownames(sce)) {
        cat(paste0(i, " is in DEG of ", j, "\n"))
        key_TFs <- c(key_TFs, i)
      }
    }
  }
  key_TFs

  for (i in x) {
    if (i %in% CTS[[CTS_ID]]) {
      cat(paste0(i, " is in CTS of ", CTS_ID, "\n"))
      key_TFs <- c(key_TFs, i)
    }
  }

  key_TFs <- unique(key_TFs)
  key_TFs

  mat <- as.data.frame(mat)

  ##  filter out the key_TFs that come from one motif thus will share target genes for the downstream predictions
  ## ====== note that this step is manual: selecting TFs of the same motif closest to the seed TF ======
  x <- NULL
  for (j in key_TFs) {
    y <- intersect(which(grepl("cisTarget_", colnames(mat))), which(Reduce("|", lapply(j, function(p) grepl(p, colnames(mat), fixed = F)))))
    cat(j, "\t", y, "\t", colnames(mat)[y], "\n")
    if (length(y) == 0) y <- 0
    x <- c(x, y[1])
  }
  names(x) <- key_TFs
  x

  # Manually choose TFs that are present in cisTarget. If two TFs are part of the same cisTarget, choose one.
  # Rebuild x fresh, one column per final key_TF (first match) -- x above is sized to the
  # EXTENDED (pre-override) key_TFs list, which is a different length than this final
  # 4-TF list; blindly reapplying names(x) <- key_TFs here would silently misalign or
  # leave a stray NA name.
  key_TFs <- c("TCF3", "HMGA2", "KLF6", "RARB")
  x <- NULL
  for (j in key_TFs) {
    y <- intersect(which(grepl("cisTarget_", colnames(mat))), which(grepl(j, colnames(mat), fixed = FALSE)))
    if (length(y) == 0) y <- NA_integer_
    x <- c(x, y[1])
  }
  names(x) <- key_TFs
  if (anyNA(x)) {
    message("[24.1] manual key_TF(s) with no cisTarget column under mm10 (dropped): ", paste(names(x)[is.na(x)], collapse = ", "))
    key_TFs <- key_TFs[!is.na(x)]
    x <- x[!is.na(x)]
  }
  x
  if (length(x) > 0) {
    for (j in seq_along(x)) {
      key <- names(x)[j]
      mat[, paste0(key, "_CP_candidate")] <- ifelse(mat[, "CP_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
      mat[, paste0(key, "_CM_candidate")] <- ifelse(mat[, "CM_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
      mat[, paste0(key, "_CF_candidate")] <- ifelse(mat[, "CF_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
    }
  }

  dim(mat)
  cat(paste0("key_TFs: ", paste(key_TFs, collapse = "_"), "\n")) # key_TFs: TCF3_HMGA2_KLF6_RARB

  if (length(key_TFs) > 0) {
    fileName <- paste0("heatmap_blocked_", CTS_name, "_scATAC_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv")
    write.table(mat, file = fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    saveRDS(list(x = x, key_TFs = key_TFs, motifAnnot_sub = motifAnnot_sub), "scATAC_cisTarget_variables.rds")
  } else {
    stop("No key TFs found for ", CTS_name, "\n")
  }
}

motif_TF_highConf <- gsub("cisTarget_", "", colnames(mat)[x]) %>% gsub("\\.motif_target", "", .)
print(motif_TF_highConf)

###  heatmap confirming key TF's self impact by checking its targets among the CTS_CP ================
for (key in motif_TF_highConf) {
  if (grepl(";", key)) {
    key_in_TFfamily <- strsplit(key, ";", fixed = T) %>%
      unlist() %>%
      intersect(key_TFs) %>%
      unique()
  } else {
    key_in_TFfamily <- key
  }

  p <- tryCatch(
    heatmap_pull_candidate(mat, graph_list, CTS_name, CHD,
      key = key_in_TFfamily, coding_genes = coding_genes, TF = TF_human,
      show_SMC_access = TRUE
    ),
    error = function(e) {
      message("heatmap skipped for '", key_in_TFfamily, "' (0 candidate genes): ", e$message)
      NULL
    }
  )
  if (!is.null(p)) {
    pdf(file = paste0("heatmap_blocked_", CTS_name, "_scATAC_cisTarget_", key_in_TFfamily, "_v3_coding_target.pdf"), height = 4)
    print(p)
    dev.off()
  }
}

##########################################################
## --  identify_TF_targeted_pull_candidate -- subset of CTS[['CP']] that are
## exclusively highly expressed (HiG) in CM (or CF)
library(BioNet)
packageVersion("BioNet")
library(igraph)
library(tibble)

mat <- read.table(paste0("heatmap_blocked_", CTS_name, "_scATAC_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv"), sep = "\t", header = T, check.names = FALSE)
dim(mat)

## IID has no HiGCTS_8 (see header note) -- alias HiG_8 in its place, matching the
## original script's own substitution.
graph_list[[paste0("HiG", CTS_name)]] <- graph_list[[paste0("HiG_", CTS_ID)]]

for (key in key_TFs) {
  key_column <- which(grepl(key, colnames(mat)) & grepl("cisTarget_", colnames(mat)))
  if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key

  graph_TF_list <- identify_TF_targeted_pull_candidate(mat, graph_list, CTS_name, CHD,
    key = key,
    keep_selfloop = TRUE,
    TF_bound_column_name = key_column,
    TF_appendix = key,
    edge_colored_by_Maven2023_ISL1KO = FALSE,
    key_in_TFfamily = key_in_TFfamily
  )
  saveRDS(graph_TF_list, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
}

names(graph_TF_list)

for (key in key_TFs) {
  if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key
  graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
  plot_TF_targeted_pull_candidate(graph_TF_list, key_in_TFfamily, CTS_name, saveFigure = TRUE)
}

##########################################################
make_layout <- function(g, seed) {
  set.seed(seed)
  layout_with_fr(g)
}

for (key in key_TFs) {
  if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key
  graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))

  res <- fill_TF_targeting_predicted_edges(graph_TF_list,
    linkage_name = "CM", graph_list,
    sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
    descendant_cluster_id = CM_cluster, TF_symbol = key_in_TFfamily,
    HVG = rownames(sce)
  )
  cat(paste0(key, ' CM vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), "\n"))
  if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CM_final.rds"))

  res_cf <- fill_TF_targeting_predicted_edges(graph_TF_list,
    linkage_name = "CF", graph_list,
    sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
    descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
    HVG = rownames(sce)
  )
  cat(paste0(key, ' CF vcount(res_cf[["g_CT_sub"]]): ', vcount(res_cf[["g_CT_sub"]]), "\n"))
  if (vcount(res_cf[["g_CT_sub"]]) > 0) saveRDS(res_cf, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CF_final.rds"))
}

### reporting ==========
print((files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds")))

final_table <- NULL
for (f in files) {
  pull <- ifelse(grepl("CF", f), "CF", "CM")
  pattern <- paste(key_TFs, collapse = "|")
  key <- key_in_TFfamily <- regmatches(f, regexpr(pattern, f))

  res <- readRDS(file = f)
  g1 <- res[["g_CT_sub"]]
  g2 <- res[["g_descendant_sub"]]

  if (vcount(g1) > 0 & vcount(g2) > 0 & vcount(g1) == vcount(g2)) {
    change_df <- edge_change_table(g1 = g1, g2 = g2, weight_attr = "weight", missing_as = 0, undirected = TRUE)
    if (nrow(change_df) == 0) {
      cat("No edge changes for", key_in_TFfamily, "-", pull, "\n")
      next
    }
    predict <- prioritize_edge_change(g1, edge_change_df = change_df, top_n = 5, title = paste0(pull, "-pull subnetwork_", key_in_TFfamily))
  } else {
    cat("No edges in", pull, "-pull subnetwork for", key_in_TFfamily, "\n")
    next
  }

  y <- motifAnnot_sub[which(motifAnnot_sub$regulators %in% key | grepl(key, motifAnnot_sub$regulators)), ]
  motif_TF_highConf_val <- unique(y$motif_TF_highConf)
  tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold &
                  (motif %in% motif_TF_highConf_val | TF_highConf %in% motif_TF_highConf_val))

  tmp1 <- if (nrow(tmp) > 0) tmp[1, , drop = FALSE] else data.frame(TF_highConf = NA, motif = NA, NES = NA)
  change_df <- cbind(linkage = pull, change_df,
                     TF_highConf = tmp1$TF_highConf, motif = tmp1$motif, NES = tmp1$NES)
  change_df$TF_highConf[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
  change_df$motif[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
  change_df$NES[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""

  final_table <- rbind(final_table, change_df)
}

dim(final_table)

write.table(final_table,
  file = paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)
