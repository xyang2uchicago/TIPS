## 24.0_acat_load_input_clean.R — code_core_mouse_8 (CORE: no ATAC / no ChIP)
##
## Faithful to GSE87038_STRING/code_core_13/24.0_acat_load_input_clean.R, but:
##   - CP/CM/CF clusters set for the C8 arm (CP=8, CM=17, CF=16)
##   - cisTarget switched to mm10 + motifAnnotations_mgi (PI's GSE140802 pattern:
##     GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/code_core_4_9vs11/24.0...R)
##   - Native mouse gene case kept throughout (NOT toupper'd) to match the mm10 database.
##     (code_core_13 uppercases mouse symbols to match the HUMAN hg38 database instead —
##     that trick does not apply once the database itself is mouse-native.)

library(ggplot2)
library(dplyr)
library(purrr)
library(SingleCellExperiment)
library(RcisTarget)

## Motif annotations — mouse (Felix / RcisTarget idiom; see GSE140802_lineage_tracking
## 7_data_MuTrans_TIPS_STRING/code_core_4_9vs11/24.0_acat_load_input_clean.R):
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

########## BEGINNING OF USER INPUT ##########
source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "GSE87038", "GSE87038_STRING/")
setwd(paste0(wd, "results_core_mouse_8/"))

data_dir    <- paste0(wd, "../data/")
## 11.x has not been rerun for the C8 arm; reuse the existing native-case weighted graph.
ppi_path    <- paste0(wd, "results/PPI_weight/")
shared_path <- paste0(shared_data_path, "/")

coding_genes <- readRDS(file = paste0(shared_path, "coding_genes_mouse.rds")) %>% unique()
length(coding_genes) # 26310
names(coding_genes) <- NULL

celltype_specific_weight_version <- '10'
celltype_col <- "label"
CP_cluster   <- "8"    # cardiac progenitor cluster
CM_cluster   <- "17"   # cardiomyocyte cluster
CF_cluster   <- "16"   # cardiac fibroblast cluster

CTS_ID  <- "8"
seed_TF <- NULL  # determined by 12.0; key_TFs resolved automatically in 24.1

heatmap_coding_target_only <- TRUE
rebuild_mat <- TRUE
########## END OF USER INPUT ##########

CTS_name <- paste0("CTS_", CTS_ID)

########################################################
##  input 1   dataset specific ---
load(paste0(data_dir, "sce_E8.25_uncorrected.RData"))
sce
names(colData(sce))

########################################################
##  input 2.1   TF annotations  -- shared ---
library(dorothea)
packageVersion("dorothea") # '1.22.0'
data(dorothea_mm, package = "dorothea")
TF_mouse <- unique(dorothea_mm$tf) # gene symbols, native case
length(TF_mouse) # 1113

CHD <- NULL

########################################################
##  input 2.2   BioTIP CTS identification  -- dataset specific ---
## CTSs were required to contain at least 30 genes (n >= 30); significance by delta-Ic with p < 0.05
## Additional filter: Ic_shrink must reach its maximum value in the corresponding cluster.
load(file = paste0(data_dir, "BioTIP.res.RData"))
names(res)
class(res$CTS.candidate) # list
lengths(res$CTS.candidate)
#   17    7   11   15   14    6   16   13    8   16
#   47   32   52   67   45   90   40   60   54   79

load(paste0(data_dir, "BioTIP_top0.1FDR0.2_IC_sim.Permutation.RData"))
names(BioTIP_scores)
dim(SimResults_g[[1]])

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
lengths(res$CTS.candidate)
#   17    7   11   15   14    6   16   13    8 16.1
#   47   32   52   67   45   90   40   60   54   79

CTS <- res$CTS.candidate[which(res$significant)]
lengths(CTS)
#    7   11   15   16   13    8 16.1
#   32   52   67   40   60   54   79
stopifnot(CTS_ID %in% names(CTS))

########################################################
##  input 3 -- data-driven --- DEGs
## 11.1 for the C8 arm was not rerun at lfc=0.5 (unlike C13); keep the lfc=0.6 DEG
## set that the existing graph_list was actually built from.
DEG <- readRDS(paste0(data_dir, "DEG_perState_min.prop0.25_lfc0.6_FDFR0.01.rds"))
lengths(DEG)
head(DEG[[CTS_ID]], 4)
class(DEG[[1]]) # [1] "character"

########################################################
##  input 4 -- data-driven --- STRING v12 PPI
library(igraph)
file <- paste0(ppi_path, "GSE87038_STRING_graph_perState_simplified_combinedweighted.rds")
graph_list <- readRDS(file)
names(graph_list)

## Native mouse case kept throughout (no toupper): matches TF_mouse, mm10 cisTarget.
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

    dim(mat) # 54  4
    colnames(mat) # "CP_hi" "CM_hi" "CF_hi" "CTS_8"
    colSums(mat)
    # CP_hi CM_hi CF_hi CTS_8
    #    10     9     3    54

    write.table(mat, file = fileName, sep = "\t", quote = FALSE,
                row.names = TRUE, col.names = TRUE)
} else {
    mat <- read.table(fileName, header = TRUE, sep = "\t", check.names = FALSE)
}

########################################################
##  input 5 -- RcisTarget: TF motif enrichment among CTS genes (mm10)
packageVersion("RcisTarget") # '1.29.0'
library(data.table)

motifAnnot <- load_motif_annotations_mgi()
dim(motifAnnot) # 252126      7

if (!file.exists("cisTarget_targets_in_all_CTS.rds")) {
    dbFile       <- "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
    feather_path <- paste0(shared_path, "cistarget/")

    if (!file.exists(paste0(feather_path, dbFile))) {
        url <- paste0(
            "https://resources.aertslab.org/cistarget/databases/mus_musculus/",
            "mm10/refseq_r80/mc_v10_clust/gene_based/", dbFile
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
