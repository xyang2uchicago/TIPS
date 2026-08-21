library(ggplot2)
library(dplyr)
library(purrr)
library(SingleCellExperiment)

########## BEGINNING OF USER INPUT ##########
code_dir <- here::here("examples", "hematoendothelial", "IbarraSoria2018", "IbarraSoria2018_STRING", "code_core_endothelial.b")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)
########## END OF USER INPUT ##########

CTS_name <- paste0("CTS_", CTS_ID)

########################################################
##  input 1   dataset specific ---
load(paste0(data_dir, 'sce_16subtype.RData'))
sce # 4000 11039
names(colData(sce))
# [1] "label"       "celltype"    "subcelltype"

########################################################
##  input 2.1   TF annotations  -- shared ---
library(dorothea)
packageVersion("dorothea") # '1.22.0'
# mouse (IbarraSoria2018 uses toupper() conversion to match GSE87038 C13 style)
data(dorothea_mm, package = "dorothea")
TF_mouse <- unique(dorothea_mm$tf) # gene symbols
length(TF_mouse) # 1113

CHD <- NULL

########################################################
##  input 2.2   BioTIP CTS identification  -- dataset specific ---
## CTSs were required to contain at least 30 genes (n >= 30); significance by delta-Ic with p < 0.05
load(file = paste0(data_dir, 'BioTIP.res.RData'))
names(res)
# [1] "CTS.candidate" "CTS.score"     "Ic.shrink"     "significant"
class(res$CTS.candidate) # list
lengths(res$CTS.candidate)
#     endothelial.b            cardiac.a presomiticMesoderm.a        endothelial.b
#                 33                   37                   71                   32

CTS <- res$CTS.candidate[2:1]
class(CTS) # list
lengths(CTS)
#   cardiac.a endothelial.b
#          37            33

########################################################
##  input 3 -- data-driven --- DEGs
DEG <- readRDS(paste0(data_dir, "DEG_perState_min.prop0.25_lfc0.6_FDFR0.01.rds"))
names(DEG)
#  [1] "blood"                  "cardiac.b"              "cardiac.c"
#  [4] "endothelial.a"          "endothelial.c"          "endothelial.d"
#  [7] "extraembryonicMesoderm" "mesodermProgenitors"     "mixedMesoderm.a"
# [10] "mixedMesoderm.b"        "pharyngealMesoderm"      "presomiticMesoderm.a"
# [13] "presomiticMesoderm.b"   "somiticMesoderm"         "endothelial.b"
# [16] "cardiac.a"
class(DEG[[1]]) # [1] "character"
head(DEG[[1]], 4)
# [1] "Cited4" "Gypa"   "Smim1"  "Klf1"

########################################################
##  input 4 -- data-driven --- STRING v12 PPI
library(igraph)
file <- paste0(ppi_path, "IbarraSoria2018_STRING_graph_perState_simplified_combinedweighted.rds")
graph_list <- readRDS(file)
names(graph_list)
# [1] "HiG_blood"                  "HiG_cardiac.b"              "HiG_cardiac.c"
# [4] "HiG_endothelial.a"          "HiG_endothelial.c"          "HiG_endothelial.d"
# [7] "HiG_extraembryonicMesoderm" "HiG_mesodermProgenitors"    "HiG_mixedMesoderm.a"
# [10] "HiG_mixedMesoderm.b"       "HiG_pharyngealMesoderm"     "HiG_presomiticMesoderm.a"
# [13] "HiG_presomiticMesoderm.b"  "HiG_somiticMesoderm"        "HiG_endothelial.b"
# [16] "HiG_cardiac.a"             "HiGCTS_endothelial.b"       "HiGCTS_cardiac.a"
# [19] "CTS_endothelial.b"         "CTS_cardiac.a"

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

    dim(mat) # [1] 33  4
    colnames(mat)
    # [1] "CP_hi"  "CM_hi"  "CF_hi"  "CTS_endothelial.b"
    colSums(mat)

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

if (!file.exists("cisTarget_targets_in_all_CTS.rds")) {
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
    saveRDS(cisTarget.res, file = "cisTarget_targets_in_all_CTS.rds")
}
