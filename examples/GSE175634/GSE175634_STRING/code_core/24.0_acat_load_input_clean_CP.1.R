library(ggplot2)
library(dplyr)
library(purrr)
library(SingleCellExperiment)

########## BEGINNING OF USER INPUT ##########
wd <- "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE175634/GSE175634_STRING/"
setwd(paste0(wd, "results_core/"))

data_dir    <- paste0(wd, "../data/")  # GSE175634/data/ shared across PPI workflows
ppi_path    <- paste0(wd, "results_core/")
shared_path <- paste0(wd, "../../Shared_Data/")

coding_genes <- readRDS(file = paste0(shared_path, "coding_genes.rds")) %>% unique()
length(coding_genes) # 19930
names(coding_genes) <- NULL

celltype_specific_weight_version <- '10'
species      <- "human"
# leiden_0.5 raw cluster numbers match the original code (3=CP, 5=CM, 1=CF)
celltype_col <- "leiden_0.5"
CP_cluster   <- "3"   # cardiogenic progenitor (transitional state)
CM_cluster   <- "5"   # cardiomyocyte (terminal CM fate)
CF_cluster   <- "1"   # alternative fate; original code cluster assignment

CTS_ID  <- "CP.1"
seed_TF <- character(0)  # HiGCTS_CP.1 contains no non-cardiac TF with BC > 0

rebuild_mat                <- TRUE
heatmap_coding_target_only <- TRUE
########## END OF USER INPUT ##########

CTS_name <- paste0("CTS_", CTS_ID)

########################################################
##  input 1 -- dataset specific
sce <- readRDS(paste0(data_dir, "sce.rds"))
sce
# dim: 3000 x 230786

########################################################
##  input 2.1 -- TF annotations -- shared
library(dorothea)
data(dorothea_hs, package = "dorothea")
TF_human <- unique(dorothea_hs$tf)
length(TF_human) # 1333

CHD <- NULL

########################################################
##  input 2.2 -- BioTIP CTS identification -- dataset specific
load(file = paste0(data_dir, "BioTIP.res.RData"))
names(res)
# [1] "CTS.candidate" "CTS.score" "Ic.shrink" "significant"

CTS <- res$CTS.candidate[which(res$significant)]
lengths(CTS)
# muscle endoderm       CP     CP.1
#     63       66       76       87

########################################################
##  input 3 -- DEGs per cluster
DEG <- readRDS(paste0(data_dir, "DEG_perState_wilcox_padj0.01_score40.rds"))
lengths(DEG)
#    0    1    2   CP    4    5    6    7 muscle    9   10 endoderm   12
#  344 1224  257  152  110  350  110  117    119 1633   81       44   25

########################################################
##  input 4 -- STRING v12 PPI graph list
library(igraph)
graph_list <- readRDS(paste0(ppi_path, "GSE175634_STRING_graph_perState_simplified_combinedweighted.rds"))
names(graph_list)
# [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"
# [5] "HiG_4"           "HiG_5"           "HiG_6"           "HiG_7"
# [9] "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"
# [13] "HiG_12"         "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"
# [17] "HiGCTS_CP.1"    "CTS_muscle"      "CTS_endoderm"    "CTS_CP"
# [21] "CTS_CP.1"

########################################################
##  build binary annotation matrix
##  Columns: CP_hi (DEG in CP), CM_hi (DEG in muscle), CF_hi (DEG in endoderm), CTS_CP.1 membership
fileName <- paste0("binary_annot_", CTS_name, "_DEG.tsv")

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
    # [1] 87  4
    colnames(mat)
    # [1] "CP_hi"    "CM_hi"    "CF_hi"    "CTS_CP.1"
    colSums(mat)
    # CP_hi  CM_hi  CF_hi CTS_CP.1
    #    29      1      4       87

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
