library(ggplot2)
library(dplyr)
library(purrr)
library(SingleCellExperiment)

########## BEGINNING OF USER INPUT ##########
wd <- "C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/"
setwd(paste0(wd, "results/results_15/"))

shared_path <- paste0(wd, "../Shared_Data/")
data_dir    <- paste0(wd, "data/")
biotip_dir  <- paste0(wd, "BioTIP/")
ppi_path    <- paste0(wd, "results/PPI_weight/")

celltype_specific_weight_version <- '10'
species      <- "human"
CHD          <- NULL               # no disease gene labeling for this dataset
celltype_col <- "leiden_0.6"        # original leiden clustering; cluster "9" intact
CTS_ID       <- "15"
CP_cluster   <- "15"               # CTS cluster = CP (same as GSE87038 pattern)
CM_cluster   <- "7"                # terminal ADRN-Proliferating (large, same state as CTS_15)
CF_cluster   <- "9"                # terminal MES cluster (alternative NCC fate)

seed_TF      <- NULL               # fill from 12.0 output once available

rebuild_mat                <- TRUE
heatmap_coding_target_only <- TRUE
########## END OF USER INPUT ##########

CTS_name <- paste0("CTS_", CTS_ID)

########################################################
##  input 1 -- SCE (original leiden_0.6 clustering; use sce_C9split.rds only for code_9)
sce <- readRDS(paste0(data_dir, "sce.rds"))
dim(sce)
table(colData(sce)[[celltype_col]])

########################################################
##  input 2.1 -- TF database (human)
library(dorothea)
data(dorothea_hs, package = "dorothea")
TF_human <- unique(dorothea_hs$tf)
length(TF_human) # 1333

########################################################
##  input 2.2 -- BioTIP CTS identification (computed in 3.5_compile_inputs.R)
load(file = paste0(data_dir, "BioTIP.res.RData"))
names(res)
# [1] "CTS.candidate" "CTS.score" "Ic.shrink" "significant"
lengths(res$CTS.candidate)

CTS <- res$CTS.candidate[res$significant]
lengths(CTS)
# 15  16   9

########################################################
##  input 3 -- DEGs per leiden_0.6 cluster
##  Note: DEG[["9_hi"]] and DEG[["9_lo"]] are NOT in this file yet —
##  they require the split DEG computation (add to 3.6 or a separate 3.7).
##  For the CTS_15 run, DEG[["15"]] and DEG[["16"]] are sufficient.
DEG <- readRDS(paste0(data_dir, "DEG_perState_min.prop0.25_lfc0.6_FDR0.05.rds"))
lengths(DEG)
head(DEG[["15"]], 4)
class(DEG[[1]])

########################################################
##  input 4 -- STRING v12 PPI graph list
library(igraph)
graph_list <- readRDS(paste0(ppi_path, "GSE116893_STRING_graph_perState_simplified_combinedweighted.rds"))
names(graph_list)

graph_list[[1]]
DEG[[1]]
CTS[[1]]

########################################################
##  build binary annotation matrix
##  Columns: CP_hi (DEG in CP_cluster), CF_hi (DEG in CF_cluster), CTS_X membership
##  No ATAC-seq or ChIP-seq data available for this dataset.
fileName <- paste0("binary_annot_", CTS_name, "_NB_DEG.tsv")

if (rebuild_mat) {
    genes <- CTS[[CTS_ID]]

    mat <- data.frame(
        CP_hi = as.integer(genes %in% DEG[[CP_cluster]]),
        CM_hi = as.integer(genes %in% DEG[[CM_cluster]]),
        CF_hi = as.integer(genes %in% DEG[[CF_cluster]]),
        row.names = genes
    )
    mat[[paste0("CTS_", CTS_ID)]] <- 1L

    dim(mat)
    colnames(mat)
    colSums(mat)

    write.table(mat, file = fileName, sep = "\t", quote = FALSE,
                row.names = TRUE, col.names = TRUE)
} else {
    mat <- read.table(fileName, header = TRUE, sep = "\t", check.names = FALSE)
}

########################################################
##  input 5 -- RcisTarget: TF motif enrichment among CTS genes
library(RcisTarget)
packageVersion("RcisTarget")
library(data.table)

data(motifAnnotations_hgnc)
motifAnnot <- motifAnnotations
dim(motifAnnot)

if (!file.exists("cisTarget_targets_in_all_CTS.rds")) {
    # hg38 gene-based rankings — download once if missing
    # https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/
    dbFile       <- "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
    feather_path <- paste0(data_dir, "cistarget/")

    if (!file.exists(paste0(feather_path, dbFile))) {
        url <- paste0(
            "https://resources.aertslab.org/cistarget/databases/homo_sapiens/",
            "hg38/refseq_r80/mc_v10_clust/gene_based/", dbFile
        )
        if (!dir.exists(feather_path)) dir.create(feather_path, recursive = TRUE)
        options(timeout = 600)
        download.file(url, destfile = paste0(feather_path, dbFile),
                      mode = "wb", method = "libcurl")
        cat("Downloaded to:", paste0(feather_path, dbFile), "\n")
    }

    motifRankings <- importRankings(paste0(feather_path, dbFile))

    # enriched motifs among CP and CF terminal DEGs
    gene_sets_HiG <- DEG[intersect(c(CM_cluster, CF_cluster), names(DEG))]
    cisTarget.res_HiG <- cisTarget(
        geneSets      = gene_sets_HiG,
        motifRankings = motifRankings,
        motifAnnot    = motifAnnot,
        nesThreshold  = 3
    )
    saveRDS(cisTarget.res_HiG, file = "cisTarget_targets_in_two_HiGs.rds")

    # enriched motifs among all significant CTS gene sets
    cisTarget.res <- cisTarget(
        geneSets      = CTS,
        motifRankings = motifRankings,
        motifAnnot    = motifAnnot,
        nesThreshold  = 3
    )
    saveRDS(cisTarget.res, file = "cisTarget_targets_in_all_CTS.rds")

    summary(cisTarget.res$NES)
    summary(cisTarget.res[cisTarget.res$geneSet == CTS_ID, ]$NES)
}
