library(EnrichedHeatmap)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(purrr)
library(SingleCellExperiment)

########## BEGINNING OF USER INPUT ##########
source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "GSE87038", "GSE87038_IID/")
dir.create(paste0(wd, "results_mouse_8/GSE181346_heart_scATAC/"), showWarnings = FALSE, recursive = TRUE)
setwd(paste0(wd, "results_mouse_8/GSE181346_heart_scATAC/"))
db_specifc_input_path <- paste0(wd, "../data/")
db_specifc_CTS_path   <- paste0(wd, "../data/")
ppi_path              <- paste0(wd, "results/PPI_weight/")
shared_path            <- paste0(shared_data_path, "/")
species      <- "mouse"
celltype_col <- "label"
CP_cluster   <- "8"
CM_cluster   <- "17"
CF_cluster   <- "18"
CMES_cluster <- "4"
CTS_ID       <- "8"
seed_TF      <- c("NR2F1", "IRX5", "ALX1")
heatmap_coding_target_only <- TRUE
########## END OF USER INPUT ##########
CTS_name <- paste0("CTS_", CTS_ID)

### the gene expression dataset to calculate the delta edge for newly added edges between TF and CTS_target
load(paste0(db_specifc_input_path, "sce_E8.25_uncorrected.RData"))
sce

########################################################
##  input 2.1   hg38 coding gene symbols & TF annotations  -- shared ---
coding_genes <- readRDS(file = paste0(shared_path, "coding_genes.rds")) %>% unique()
length(coding_genes) # 19930
names(coding_genes) <- NULL

## get the TF database
library(dorothea)
packageVersion("dorothea") # '1.22.0'
data(dorothea_hs, package = "dorothea")
TF_human <- unique(dorothea_hs$tf) # gene symbols
length(TF_human) # 1333

CHD <- readRDS(paste0(shared_path, "CHD_Cilia_Genelist.rds"))
CHD <- CHD$Griffin2023_PCGC_AllCurated
length(CHD) # 295

########################################################
##  input 2.2   maps -- shared ---
maps <- read.delim(file = paste0(shared_path, "readme_filename_map_xy.txt"), sep = "\t", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
maps <- maps[grepl("hft_", maps$file_name), ]

########################################################
##  input 2.3 published gene openness score -- shared ---
library(openxlsx)
annotation_vitro <- read.xlsx(xlsxFile = paste0(db_specifc_input_path, "Ameen2022cell-supplement-10.xlsx"), sheet = 4)

Ameen_annotation_list <- list()
for (i in unique(annotation_vitro$celltype)) {
  Ameen_annotation_list[[i]] <- annotation_vitro$Gene.name[which(annotation_vitro$celltype == i)]
}

annotation_fetal <- read.xlsx(xlsxFile = paste0(db_specifc_input_path, "Ameen2022cell_Table1_sheet3_markergene_Laksshman2026update.xlsx"), sheet = 1)
colnames(annotation_fetal)[9] <- "Cluster"
annotation_fetal$Cluster[which(annotation_fetal$Cluster == "vSMC")] <- "SMC"
short_cluster_names <- unique(annotation_fetal$Cluster)

for (i in short_cluster_names) {
  Ameen_annotation_list[[i]] <- annotation_fetal$name[which(annotation_fetal$Cluster == i)]
}

maps$cluster_abb <- maps$identity_simple
maps$cluster_abb[which(maps$identity_simple == "arteries")] <- "aEC"
maps$cluster_abb[which(maps$identity_simple == "atrial_CM")] <- "aCM"
maps$cluster_abb[which(maps$identity_simple == "capillary")] <- "Cap"
maps$cluster_abb[which(maps$identity_simple == "pre_CF")] <- "preCF"
maps$cluster_abb[which(maps$identity_simple == "matureCF")] <- "CF"
maps$cluster_abb[which(maps$identity_simple == "CF_progenitors")] <- "CFP"
maps$cluster_abb[which(maps$identity_simple == "endocardium")] <- "Endo"
maps$cluster_abb[which(maps$identity_simple == "epicardium")] <- "EPC"
maps$cluster_abb[which(maps$identity_simple == "lymphatic_endothelial")] <- "lEC"
maps$cluster_abb[which(maps$identity_simple == "early_CM")] <- "eCM"
maps$cluster_abb[which(maps$identity_simple == "neural_crest")] <- "NC"
maps$cluster_abb[which(maps$identity_simple == "pericytes")] <- "PC"
maps$cluster_abb[which(maps$identity_simple == "pre_SMC")] <- "preSMC"
maps$cluster_abb[which(maps$identity_simple == "SMC")] <- "SMC"
maps$cluster_abb[which(maps$identity_simple == "Outflow_tract")] <- "OFT"
maps$cluster_abb[which(maps$identity_simple == "Fibroblast_like_1")] <- "FB1"
maps$cluster_abb[which(maps$identity_simple == "Fibroblast_like_2")] <- "FB2"
maps$cluster_abb[which(maps$identity_simple == "Venal_endothelial")] <- "vEC"
maps$cluster_abb[which(maps$identity_simple == "Ventricular_CM")] <- "vCM"

colnames(annotation_vitro)[4] <- "name"
colnames(annotation_vitro)[8] <- "Cluster"
tmp <- intersect(colnames(annotation_vitro), colnames(annotation_fetal))
ATAC_anno_df <- rbind(as.data.frame(annotation_vitro[, tmp]), as.data.frame(annotation_fetal[, tmp]))
rm(annotation_vitro, annotation_fetal)

###############################################################
## BioTIP CTS identification (run previously -- NOT repeat)
###############################################################
load(file = paste0(db_specifc_CTS_path, "BioTIP.res.RData"))
class(res$CTS.candidate) # list
lengths(res$CTS.candidate)

load(paste0(db_specifc_CTS_path, "BioTIP_top0.1FDR0.2_IC_sim.Permutation.RData"))

x <- which(duplicated(names(BioTIP_scores)))
if (length(x) > 0) {
  if (all(names(BioTIP_scores) == names(res$CTS.candidate))) {
    x <- which(duplicated(names(res$CTS.candidate)))
    names(res$CTS.candidate)[x] <- paste(names(res$CTS.candidate)[x], "1", sep = ".")
    names(res$CTS.score[x]) <- names(res$CTS.candidate)[x]
    names(res$Ic.shrink[x]) <- names(res$CTS.candidate)[x]
    names(res$significant[x]) <- names(res$CTS.candidate)[x]
    names(SimResults_g)[x] <- names(res$CTS.candidate)[x]
    names(BioTIP_scores)[x] <- names(res$CTS.candidate)[x]
  } else {
    cat("pls let names(BioTIP_scores) == names(res$CTS.candidate)")
  }
  rm(x)
}

(Top_cluster <- lapply(BioTIP_scores, function(x) which.max(x) %>% names()) %>% unlist())
x <- grep(".", names(Top_cluster), fixed = T, value = T)
if (length(x) > 0) names(Top_cluster)[which(names(Top_cluster) == x)] <- unlist(strsplit(x, split = ".", fixed = T))[1]

CTS <- res$CTS.candidate[which(res$significant == TRUE & Top_cluster == names(Top_cluster))]
lengths(CTS)
stopifnot(CTS_ID %in% names(CTS))

########################################################
##  input 3 -- data-driven --- DEGs
DEG <- readRDS(paste0(db_specifc_input_path, "DEG_perState_min.prop0.25_lfc0.6_FDFR0.01.rds"))
lengths(DEG)
head(DEG[[CTS_ID]], 4)

########################################################
##  input 4 -- data-driven --- IID PPI
library(igraph)
file <- paste0(ppi_path, "GSE87038_IID_graph_perState_simplified_combinedweighted.rds")
graph_list <- readRDS(file)
if ("CTS_13" %in% names(graph_list)) graph_list <- graph_list[-which(names(graph_list) == "CTS_13")]
names(graph_list)

############## uppercase the mouse symbols to match the gene in the annotation datasets
if (species == "mouse") {
  for (i in seq_along(graph_list)) V(graph_list[[i]])$name <- toupper(V(graph_list[[i]])$name)
  for (i in seq_along(DEG)) DEG[[i]] <- toupper(DEG[[i]])
  for (i in seq_along(CTS)) CTS[[i]] <- toupper(CTS[[i]])
  rownames(sce) <- toupper(rownames(sce))
  cat("uppercase the mouse symbols to match the gene in the annotation datasets done\n")
}
graph_list[[1]]
DEG[[1]]
CTS[[1]]

################# create ATAC_anno_df table and save it ####################
fileName <- paste0("binary_annot_", CTS_name, "_scATAC_Maven2023_gene_ISL1_v3.tsv")

if (rebuild_mat) {
  annotation_sub <- subset(ATAC_anno_df, name %in% CTS[[CTS_ID]])

  annotation_sub$PCW6CP_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("CFP", "OFT", "FB1", "eCM", "i-CP")]), "1", "0")
  annotation_sub$PCW8_CM_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("aCM", "i-CM")]), "1", "0")
  annotation_sub$PCW8_CF_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("preCF", "FB2", "EPC", "i-CF")]), "1", "0")
  annotation_sub$PCW8_SMC_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("preSMC", "EPC", "i-SMC")]), "1", "0")
  annotation_sub$PCW19_CM_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("vCM")]), "1", "0")
  annotation_sub$PCW19_CF_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("CF", "EPC")]), "1", "0")
  annotation_sub$PCW19_SMC_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("SMC", "EPC")]), "1", "0")
  annotation_sub$PCW6_CF_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("CFP", "FB1")]), "1", "0")
  annotation_sub$PCW6_SMC_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("OFT")]), "1", "0")
  annotation_sub$PCW6_CM_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("eCM")]), "1", "0")
  annotation_sub$iEPC_access <- ifelse(annotation_sub$name %in% unlist(Ameen_annotation_list[c("i-EPC")]), "1", "0")

  annotation_sub$CMES_hi <- ifelse(annotation_sub$name %in% DEG[[CMES_cluster]], "1", "0")
  annotation_sub$CP_hi <- ifelse(annotation_sub$name %in% DEG[[CP_cluster]], "1", "0")
  annotation_sub$CM_hi <- ifelse(annotation_sub$name %in% DEG[[CM_cluster]], "1", "0")
  annotation_sub$CF_hi <- ifelse(annotation_sub$name %in% DEG[[CF_cluster]], "1", "0")

  int <- c(
    "CMES_hi", "CP_hi", "CM_hi", "CF_hi", "PCW6CP_access",
    "PCW8_CM_access", "PCW19_CM_access",
    "PCW8_CF_access", "PCW19_CF_access",
    "PCW8_SMC_access", "PCW19_SMC_access",
    "PCW6_CM_access", "PCW6_CF_access", "PCW6_SMC_access", "iEPC_access"
  )
  mat <- as.matrix(sapply(annotation_sub[, int], as.numeric))
  rownames(mat) <- annotation_sub$name
  mat <- rowsum(mat, group = rownames(mat))
  mat <- (mat > 0) * 1

  (missing <- setdiff(CTS[[CTS_ID]], rownames(mat)))
  missing_mat <- matrix(0, nrow = length(missing), ncol = ncol(mat))
  rownames(missing_mat) <- missing
  colnames(missing_mat) <- colnames(mat)
  mat <- rbind(mat, missing_mat)

  mat <- as.data.frame(mat)
  mat[, paste0("CTS_", CTS_ID)] <- ifelse(rownames(mat) %in% CTS[[CTS_ID]], "1", "0")
  mat[, "CP_hi"] <- ifelse(rownames(mat) %in% DEG[[CP_cluster]], "1", "0")
  mat[, "CM_hi"] <- ifelse(rownames(mat) %in% DEG[[CM_cluster]], "1", "0")
  mat[, "CF_hi"] <- ifelse(rownames(mat) %in% DEG[[CF_cluster]], "1", "0")
  mat[, "CMES_hi"] <- ifelse(rownames(mat) %in% DEG[[CMES_cluster]], "1", "0")

  ########################################################
  ##  input 5 -- shared --- published ISL1-CP binding
  ISL1_set <- readRDS(file = paste0(shared_path, "GSE195476_ISL1/ISL1_set.rds"))
  names(ISL1_set) <- paste0("Maven2023_gene_", names(ISL1_set))

  load(paste0(shared_path, "GSE80383_Isl1/ISL1_Mm2Hg.GS_uniqueSymbol_FC1.2_padj0.05.RData"))
  names(ISL1)[2] <- "Isl1.iCPC_CPC.bound"
  names(ISL1)[16] <- "iPSC_CPC.open"

  mat[, "Maven2023_gene_ISL1_up_E"] <- ifelse(rownames(mat) %in% ISL1_set[["Maven2023_gene_ISL_up_E"]], "1", "0")
  mat[, "Maven2023_gene_ISL1_up_T"] <- ifelse(rownames(mat) %in% ISL1_set[["Maven2023_gene_ISL_up_T"]], "1", "0")
  mat[, "Maven2023_gene_ISL1_up_L"] <- ifelse(rownames(mat) %in% ISL1_set[["Maven2023_gene_ISL_up_L"]], "1", "0")
  mat[, "Maven2023_gene_ISL1_dn_E"] <- ifelse(rownames(mat) %in% ISL1_set[["Maven2023_gene_ISL_dn_E"]], "1", "0")
  mat[, "Maven2023_gene_ISL1_dn_T"] <- ifelse(rownames(mat) %in% ISL1_set[["Maven2023_gene_ISL_dn_T"]], "1", "0")
  mat[, "Maven2023_gene_ISL1_dn_L"] <- ifelse(rownames(mat) %in% ISL1_set[["Maven2023_gene_ISL_dn_L"]], "1", "0")
  mat[, "Maven2023_gene_ISL1_WT_d6CP"] <- ifelse(rownames(mat) %in% ISL1_set[["Maven2023_gene_ISL1_WT_d6CP"]], "1", "0")
  mat[, "Gao2019_gene_Isl1_E825E9.bound"] <- ifelse(rownames(mat) %in% ISL1[["Isl1.embryo.bound"]], "1", "0")
  mat[, "Gao2019_gene_Isl1.iCPC_CPC.bound"] <- ifelse(rownames(mat) %in% ISL1[["Isl1.iCPC_CPC.bound"]], "1", "0")

  ####################################################
  ### extract subnetworks and add ISL1-bound links ###
  ####################################################
  mat <- as.data.frame(mat)
  mat[, "ISL1_CP_bound"] <- ifelse((mat[, "Gao2019_gene_Isl1_E825E9.bound"] == 1 |
    mat[, "Gao2019_gene_Isl1.iCPC_CPC.bound"] == 1 |
    mat[, "Maven2023_gene_ISL1_WT_d6CP"] == 1), 1, 0)

  mat[, "ISL1_CP_candidate"] <- ifelse(mat[, "CP_hi"] == 1 & mat[, "ISL1_CP_bound"] == 1, 1, 0)
  mat[, "ISL1_CM_candidate"] <- ifelse(mat[, "CM_hi"] == 1 & mat[, "ISL1_CP_bound"], 1, 0)
  mat[, "ISL1_CF_candidate"] <- ifelse(mat[, "CF_hi"] == 1 & mat[, "ISL1_CP_bound"], 1, 0)
  if ("SMC_hi" %in% colnames(mat)) mat[, "ISL1_SMC_candidate"] <- ifelse(mat[, "SMC_hi"] == 1 & mat[, "ISL1_CP_bound"], 1, 0)

  dim(mat)
  write.table(mat, file = fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
} else {
  mat <- read.table(fileName, header = TRUE, sep = "\t", check.names = FALSE)
}

########################################################
##  input 6 -- data-driven --- RcisTarget predicted TF that are enriched among CTS genes
library(RcisTarget)
packageVersion("RcisTarget") # '1.29.0'
library(data.table)

## Mouse cistarget (mm10), uppercased to match this pipeline's existing uppercase-mouse
## convention -- see GSE87038_STRING/code_mouse_8/24.0...R for the full rationale
## (verified ~98% CTS/DEG gene overlap against the uppercased mm10 rankings, 2026-07-28).
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

motifAnnot <- load_motif_annotations_mgi()
dim(motifAnnot)

if (!file.exists(file = "cisTarget_targets_in_all_CTS.rds")) {
  dbFile <- "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
  feather_path <- paste0(shared_path, "cistarget/")

  if (!file.exists(file = paste0(feather_path, dbFile))) {
    url <- paste0(
      "https://resources.aertslab.org/cistarget/databases/mus_musculus/",
      "mm10/refseq_r80/mc_v10_clust/gene_based/", dbFile
    )
    if (!dir.exists(feather_path)) dir.create(feather_path, recursive = TRUE)
    options(timeout = 600)
    download.file(url, destfile = paste0(feather_path, dbFile), mode = "wb", method = "libcurl")
  }

  motifRankings <- importRankings(paste0(feather_path, dbFile))
  colnames(motifRankings@rankings) <- toupper(colnames(motifRankings@rankings))

  cisTarget.res_HiG <- cisTarget(
    geneSets = DEG[c(CM_cluster, CF_cluster)],
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
