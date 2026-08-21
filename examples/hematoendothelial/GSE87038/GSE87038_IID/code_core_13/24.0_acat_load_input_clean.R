library(ggplot2)
library(dplyr)
library(purrr)
library(SingleCellExperiment)

########## BEGINNING OF USER INPUT ##########
code_dir <- here::here("examples", "hematoendothelial", "GSE87038", "GSE87038_IID", "code_core_13")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)
rebuild_mat <- TRUE
########## END OF USER INPUT ##########

CTS_name <- paste0("CTS_", CTS_ID)

########################################################
##  input 1   dataset specific ---
load(paste0(data_dir, 'sce_E8.25_uncorrected.RData'))
sce # 10938 7240
names(colData(sce))
# [1] "cell"             "barcode"          "sample"           "pool"             "stage"
#  [6] "sequencing.batch" "theiler"          "doub.density"     "doublet"          "cluster"
# [11] "cluster.sub"      "cluster.stage"    "cluster.theiler"  "stripped"         "celltype"
# [16] "colour"           "sizeFactor"       "batch"            "label"

########################################################
##  input 2.1   TF annotations  -- shared ---
library(dorothea)
packageVersion("dorothea") # '1.22.0'
data(dorothea_mm, package = "dorothea")
TF_mouse <- unique(dorothea_mm$tf)
length(TF_mouse) # 1113

CHD <- NULL

########################################################
##  input 2.2   BioTIP CTS identification  -- dataset specific ---
load(file = paste0(data_dir, 'BioTIP.res.RData'))
names(res)
# [1] "CTS.candidate" "CTS.score"     "Ic.shrink"     "significant"
class(res$CTS.candidate) # list
lengths(res$CTS.candidate)
# 17  7 11 15 14  6 16 13  8 16
# 47 32 52 67 45 90 40 60 54 79

load(paste0(data_dir, "BioTIP_top0.1FDR0.2_IC_sim.Permutation.RData"))
names(BioTIP_scores)
# [1] "17" "7"  "11" "15" "14" "6"  "16" "13" "8"  "16"
dim(SimResults_g[[1]]) # 19  100

x <- which(duplicated(names(BioTIP_scores)))
if (length(x) > 0) {
    if (all(names(BioTIP_scores) == names(res$CTS.candidate))) {
        x <- which(duplicated(names(res$CTS.candidate)))
        names(res$CTS.candidate)[x] <- paste(names(res$CTS.candidate)[x], "1", sep = ".")
        names(res$CTS.score[x])     <- names(res$CTS.candidate)[x]
        names(res$Ic.shrink[x])     <- names(res$CTS.candidate)[x]
        names(res$significant[x])   <- names(res$CTS.candidate)[x]
        names(SimResults_g)[x]      <- names(res$CTS.candidate)[x]
        names(BioTIP_scores)[x]     <- names(res$CTS.candidate)[x]
    }
    rm(x)
}
lengths(res$CTS.candidate)
# 17    7   11   15   14    6   16   13    8 16.1
# 47   32   52   67   45   90   40   60   54   79

(Top_cluster <- lapply(BioTIP_scores, function(x) which.max(x) %>% names()) %>% unlist())
# 17    7   11   15   14    6   16   13    8 16.1
# "2"  "7" "11" "15" "12" "10" "16" "10"  "8" "16"

x <- grep(".", names(Top_cluster), fixed = T, value = T)
if (length(x) > 0) names(Top_cluster)[which(names(Top_cluster) == x)] <- unlist(strsplit(x, split = ".", fixed = T))[1]

CTS <- res$CTS.candidate[which(res$significant)]
lengths(CTS)
#    7   11   15   16 16.1   13    8
#   32   52   67   40   79   60   54

########################################################
##  input 3 -- data-driven --- DEGs
DEG <- readRDS(paste0(data_dir, "DEG_perState_min.prop0.25_lfc", logFC.cut, "_FDFR0.01.rds"))
lengths(DEG)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
# 423 566 549 425 431 637 412 448 480 577 572 452 531 475 537 657 653 578 677
head(DEG[["8"]], 4)
# [1] "Arg1"    "Prrx2"   "Igfbpl1" "Id2"
class(DEG[[1]]) # [1] "character"

########################################################
##  input 4 -- data-driven --- IID PPI
library(igraph)
file <- paste0(ppi_path, "GSE87038_IID_graph_perState_simplified_combinedweighted.rds")
graph_list <- readRDS(file)
names(graph_list)
# [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"
# [7] "HiG_9"     "HiG_10"    "HiG_12"    "HiG_14"    "HiG_17"    "HiG_18"
# [13] "HiG_19"    "HiG_7"     "HiG_11"    "HiG_15"    "HiG_16"    "HiG_13"
# [19] "HiG_8"     "HiGCTS_15" "CTS_15"    "CTS_16.1"  "CTS_13"    "CTS_8"
# (24 networks; CTS_7/11/16 and HiGCTS_7/11/16/16.1/13/8 excluded by <10 edges in IID PPI)

if (species == "mouse") {
    for (i in seq_along(graph_list)) V(graph_list[[i]])$name <- toupper(V(graph_list[[i]])$name)
    for (i in seq_along(DEG)) DEG[[i]] <- toupper(DEG[[i]])
    for (i in seq_along(CTS)) CTS[[i]] <- toupper(CTS[[i]])
    rownames(sce) <- toupper(rownames(sce))
    cat("uppercase mouse symbols done\n")
}

graph_list[[1]]
DEG[[1]]
CTS[[1]]

########################################################
##  build binary annotation matrix
##  Columns: CP_hi (DEG in CP), CM_hi (DEG in CM), CF_hi (DEG in CF), CTS_X membership
##  No ATAC-seq or ChIP-seq data used.
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

    dim(mat) # [1] 60  4
    colnames(mat)
    # [1] "CP_hi"  "CM_hi"  "CF_hi"  "CTS_13"
    colSums(mat)
    # CP_hi CM_hi CF_hi CTS_13
    #    15     5    38     60
    # (same as STRING code_core_13: mat uses BioTIP CTS_13 genes, same DEG lfc=0.5)

    write.table(mat, file = fileName, sep = "\t", quote = FALSE,
                row.names = TRUE, col.names = TRUE)
} else {
    mat <- read.table(fileName, header = TRUE, sep = "\t", check.names = FALSE)
}

########################################################
##  input 5 -- RcisTarget: TF motif enrichment among CTS genes
library(RcisTarget)
packageVersion("RcisTarget") # '1.29.0'
library(data.table)

data(motifAnnotations_hgnc)
motifAnnot <- motifAnnotations
dim(motifAnnot) # [1] 253096      8

ct_file <- "cisTarget_targets_in_all_CTS.rds"

if (!file.exists(ct_file)) {
    motifRankings <- importRankings(tips_core_ensure_cistarget_feather(shared_data_path))

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
    saveRDS(cisTarget.res, file = ct_file)
    cat("Regenerated cisTarget_targets_in_all_CTS.rds with CTS_13\n")
}
