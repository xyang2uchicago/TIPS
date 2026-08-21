library(ggplot2)
library(dplyr)
library(purrr)
library(SingleCellExperiment)

########## BEGINNING OF USER INPUT ##########
code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_IID", "code_core")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)
rebuild_mat <- TRUE
CP_cluster  <- "3"   # cardiogenic progenitor (transitional state)
########## END OF USER INPUT ##########

CTS_name <- paste0("CTS_", CTS_ID)

########################################################
##  input 1 -- dataset specific
sce <- readRDS(paste0(data_dir, "sce.rds"))
sce
# dim: 3000 x 230786
names(colData(sce))
# [1] "orig.ident"                       "nCount_RNA"
# [3] "nFeature_RNA"                     "cell"
# [5] "exp.grp"                          "sample"
# [7] "diffday"                          "individual"
# [9] "S.Score"                          "G2M.Score"
# [11] "CC.Difference"                   "demux.dbl.prb"
# [13] "leiden_published"                "type"
# [15] "percent.mt"                      "nCount_SCT"
# [17] "nFeature_SCT"                    "PC1"
# [19] "PC2"                             "leiden_0.1"
# [21] "leiden_0.2"                      "leiden_0.3"
# [23] "leiden_0.4"                      "leiden_0.5"
# [25] "leiden_0.6"                      "leiden_0.7"
# [27] "leiden_0.8"                      "leiden_0.9"
# [29] "leiden_1"                        "leiden_0.5_type"
# [31] "dpt_pseudotime_published"        "dpt_pseudotime_leiden0.5_root382"

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
##  input 4 -- IID PPI graph list
library(igraph)
graph_list <- readRDS(paste0(ppi_path, "GSE175634_IID_graph_perState_simplified_combinedweighted.rds"))
names(graph_list)
# [1] "HiG_0"           "HiG_1"           "HiG_2"           "HiG_CP"
# [5] "HiG_4"           "HiG_5"           "HiG_6"           "HiG_7"
# [9] "HiG_muscle"      "HiG_9"           "HiG_10"          "HiG_endoderm"
# [13] "HiG_12"         "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"
# [17] "HiGCTS_CP.1"    "CTS_muscle"      "CTS_endoderm"    "CTS_CP"
# [21] "CTS_CP.1"

## IID gene symbols may be mixed case — upcase for consistency with DEG/CTS
for (i in seq_along(graph_list)) V(graph_list[[i]])$name <- toupper(V(graph_list[[i]])$name)

########################################################
##  build binary annotation matrix
##  Columns: CP_hi (DEG in CP), CM_hi (DEG in muscle), CF_hi (DEG in endoderm), CTS_CP membership
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
    # [1] 76  4
    colnames(mat)
    # [1] "CP_hi"  "CM_hi"  "CF_hi"  "CTS_CP"
    colSums(mat)
    # CP_hi  CM_hi  CF_hi CTS_CP
    #     0      0     18     57

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
