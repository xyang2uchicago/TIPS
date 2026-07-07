library(ggplot2)
library(dplyr)
library(purrr)
library(SingleCellExperiment)

########## BEGINNING OF USER INPUT ##########
wd <- "C:/Users/felix/Documents/GitHub/TIPS/examples/Yu2025/"
setwd(paste0(wd, "results/results_9/"))

shared_path  <- paste0(wd, "../Shared_Data/")
data_dir     <- paste0(wd, "data/")
biotip_dir   <- paste0(wd, "BioTIP/")
ppi_path     <- paste0(wd, "results/PPI_weight/")
coding_genes <- readRDS(paste0(shared_path, "coding_genes.rds")) %>% unique()

celltype_specific_weight_version <- '10'
species      <- "human"
CHD          <- NULL               # no disease gene labeling for this dataset
celltype_col <- "leiden_C9split"   # required: 9_hi and 9_lo only exist in this column
CTS_ID       <- "9"
CP_cluster   <- "9_hi"             # NCC-high progenitor cells (migratory NCC, more plastic)
CM_cluster   <- "7"                # terminal ADRN-Proliferating (main ADRN fate)
CF_cluster   <- "9_lo"             # NCC-low committed cells (alternative MES/NCC fate)

seed_TF      <- NULL

rebuild_mat                <- TRUE
heatmap_coding_target_only <- TRUE
########## END OF USER INPUT ##########

CTS_name <- paste0("CTS_", CTS_ID)

########################################################
##  input 1 -- SCE with C9 NCC split (required for leiden_C9split)
sce <- readRDS(paste0(data_dir, "sce_C9split.rds"))
dim(sce)
table(colData(sce)[[celltype_col]])
#     1    10    11    12    13    14    15    16    17     2     3     4     5
# 28394 10039  2911  1862  1834   642   348   189    73 28380 27401 25381 24098
#     6     7     8  9_hi  9_lo
# 19850 19743 18215  1001 10867

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

CTS <- res$CTS.candidate[res$significant]
lengths(CTS)
#  15  16   9
#  76 139 110

########################################################
##  input 3 -- DEGs per leiden_C9split (includes 9_hi and 9_lo from 3.7_compute_split_DEG.R)
DEG <- readRDS(paste0(data_dir, "DEG_perState_min.prop0.25_lfc0.6_FDR0.05.rds"))
lengths(DEG)
# DEG[["9_hi"]]: 613 genes (NCC-high markers)
# DEG[["9_lo"]]: 231 genes (NCC-low markers)

########################################################
##  input 4 -- STRING v12 PPI graph list
library(igraph)
graph_list <- readRDS(paste0(ppi_path, "Yu2025_STRING_graph_perState_simplified_combinedweighted.rds"))
# HiGCTS_9 (5v, 0e) is excluded from _combinedweighted (too small for specificity); load locally
ns_list <- readRDS(paste0(wd, "results/Yu2025_STRING_graph_perState_notsimplified.rds"))
graph_list[["HiGCTS_9"]] <- simplify(ns_list[["HiGCTS_9"]], edge.attr.comb = "max")
rm(ns_list)
names(graph_list)
# includes HiG_9_hi (611v), HiG_9_lo (230v), HiGCTS_9 (5v, 0e)

########################################################
##  build binary annotation matrix
##  Columns: CP_hi (DEG in 9_hi), CM_hi (DEG in 7), CF_hi (DEG in 9_lo), CTS_9 membership
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

# cisTarget RDS was copied from results_15 (covers all 3 significant CTS: 15, 16, 9)
if (!file.exists("cisTarget_targets_in_all_CTS.rds")) {
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
    }

    motifRankings <- importRankings(paste0(feather_path, dbFile))

    gene_sets_HiG <- DEG[intersect(c(CM_cluster, CF_cluster), names(DEG))]
    cisTarget.res_HiG <- cisTarget(
        geneSets      = gene_sets_HiG,
        motifRankings = motifRankings,
        motifAnnot    = motifAnnot,
        nesThreshold  = 3
    )
    saveRDS(cisTarget.res_HiG, file = "cisTarget_targets_in_two_HiGs.rds")

    cisTarget.res <- cisTarget(
        geneSets      = CTS,
        motifRankings = motifRankings,
        motifAnnot    = motifAnnot,
        nesThreshold  = 3
    )
    saveRDS(cisTarget.res, file = "cisTarget_targets_in_all_CTS.rds")
}
