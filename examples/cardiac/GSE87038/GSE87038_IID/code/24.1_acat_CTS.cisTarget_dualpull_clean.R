# Set True if running 24.0 the first time
rebuild_mat <- TRUE
source(here::here("examples", "config.R"))
wd <- tips_path("examples", "cardiac", "GSE87038", "GSE87038_IID/")
source(here::here("examples", "cardiac", "GSE87038", "GSE87038_IID", "code", "24.0_acat_load_input_clean.R"))

source(here::here("R", "celltype_specific_weight_v10.R"))

## check the loaded objects =========================
seed_TF # "NR2F1" "IRX5"  "ALX1"  # by top PageRank in code 12.xxxx
names(graph_list)
#  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"
#  [7] "HiG_9"     "HiG_10"    "HiG_12"    "HiG_14"    "HiG_17"    "HiG_18"
# [13] "HiG_19"    "HiG_7"     "HiG_11"    "HiG_15"    "HiG_16"    "HiG_13"
# [19] "HiG_8"     "HiGCTS_15" "CTS_7"     "CTS_11"    "CTS_15"    "CTS_16"
# [25] "CTS_16.1"  "CTS_8"
names(DEG)
#  [1] "1"  "2"  "3"  "4"  "5"  "6"  "9"  "10" "12" "14" "17" "18" "19" "7"  "11" "15" "16" "13" "8"

celltype_col # [1] "label"
CP_cluster # '8'
CM_cluster # '17'
CF_cluster # '18'
CMES_cluster # '4'

lengths(CTS) # a lsit
#   7 11 15 16  8  16.1
#  32 52 67 40  54  79

class(sce) # [1] "SingleCellExperiment"

dim(mat) # [1] 54  29
colnames(mat)
#  [1] "CMES_hi"                          "CP_hi"
#  [3] "CM_hi"                            "CF_hi"
#  [5] "PCW6CP_access"                    "PCW8_CM_access"
#  [7] "PCW19_CM_access"                  "PCW8_CF_access"
#  [9] "PCW19_CF_access"                  "PCW8_SMC_access"
# [11] "PCW19_SMC_access"                 "PCW6_CM_access"
# [13] "PCW6_CF_access"                   "PCW6_SMC_access"
# [15] "iEPC_access"                      "CTS_8"
# [17] "Maven2023_gene_ISL1_up_E"         "Maven2023_gene_ISL1_up_T"
# [19] "Maven2023_gene_ISL1_up_L"         "Maven2023_gene_ISL1_dn_E"
# [21] "Maven2023_gene_ISL1_dn_T"         "Maven2023_gene_ISL1_dn_L"
# [23] "Maven2023_gene_ISL1_WT_d6CP"      "Gao2019_gene_Isl1_E825E9.bound"
# [25] "Gao2019_gene_Isl1.iCPC_CPC.bound" "ISL1_CP_bound"
# [27] "ISL1_CP_candidate"                "ISL1_CM_candidate"
# [29] "ISL1_CF_candidate"


seed_TF # "NR2F1" "IRX5"  "ALX1"
CTS_ID # '8'
CTS_name #' CTS_8'

## set subfold =======================
(updir <- getwd())
# "/Users/felixyu/Documents/GitHub/TIPS/examples/cardiac/GSE87038/results/GSE181346_heart_scATAC"
# create the directory if it doesn't exist?
dir.create(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)),
  showWarnings = FALSE, recursive = TRUE
)
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

NES_threshold <- 3

########################################################
##  input 6 -- data-driven --- RcisTarget predicted TF that are enriched among CTS genes
library(RcisTarget)
packageVersion("RcisTarget") # ‘1.29.0’
library(data.table)

data(package = "RcisTarget") # list which motif-annotation objects you actually have
# data(motifAnnotations_hgnc_v9)  ## mouse: data(motifAnnotations_mgi)
# dim(motifAnnotations_hgnc_v9)  # [1] 163192      7

data(motifAnnotations_hgnc)
dim(motifAnnotations) # [1] 253096      8

motifAnnot <- motifAnnotations
dim(motifAnnot) # [1] 253096      8


cisTarget.res <- readRDS(file = "../cisTarget_targets_in_all_CTS.rds")
# (NES_threshold = quantile(cisTarget.res[which(cisTarget.res$geneSet==CTS_ID),]$NES, seq(0,1,0.1))['80%'])  #3.642
write.table(subset(cisTarget.res, NES >= NES_threshold & geneSet %in% c(CP_cluster)),
  paste0("../cisTarget_targets_in_", CP_cluster, "_NES", NES_threshold, ".txt"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)


########################################################
##   --- binary flag CTS genes

dim(mat) # 54 29
### = load the binary annotation matrix build in code 24.0xxx ============================
files <- list.files(pattern = "^heatmap_blocked_CTS_") %>% grep("_v3.tsv", ., value = TRUE)
if (length(files) > 0) {
  fileName <- grep("_v3.tsv", files, value = TRUE)
  mat <- read.table(fileName, sep = "\t", header = T, check.names = FALSE)

  saved_variables <- readRDS(file = "scATAC_cistTarget_variables.rds")

  # Recreate variables
  x <- saved_variables$x
  key_TFs <- saved_variables$key_TFs
  motifAnnot_sub <- saved_variables$motifAnnot_sub
} else {
  #### focus on the this CTS member gene-enriched motifs, with NES_threshold control =========================
  motifAnnot_sub <- get_regulators_from_motifs(cisTarget.res, CTS_ID, NES_threshold, motifAnnot = motifAnnot)
  keys <- unique(motifAnnot_sub$regulators)
  if (any(is.na(keys))) keys <- keys[!is.na(keys)]
  keys
  #  [1] "PARP1"
  #  [2] "HOXD3"
  #  [3] "ZSCAN9"
  #  [4] "CEBPB;PDX1;STAT6"
  #  [5] "TCF3"
  #  [6] "POU3F3"
  #  [7] "AHCTF1"

  #  [8] "ISX;LHX2;LHX9;SHOX2"

  #  [9] "FOXM1;IKZF1;TBL1XR1"

  # [10] "FOXM1"

  # [11] "HMGA1;HMGA2"

  # [12] "STAT5A"

  # [13] "POU4F3"

  # [14] "RREB1"

  # [15] "PRDM5;ZNF324;ZNF324B;ZNF341;ZNF580"

  # [16] "SF1"

  # [17] "HMG20A"

  # [18] "KLF2"

  # [19] "BCL11A;BCL6;CHURC1;DDX20;DPF2;EGR1;EGR2;EP300;FLI1;GTF2F1;HINFP;IKZF1;IKZF3;IRF4;KLF10;KLF13;KLF15;KLF16;KLF6;KLF9;MAZ;MTF1;MZF1;NFXL1;OVOL2;PATZ1;PAX4;PML;PRDM9;PURA;RARA;RARB;RARG;RBAK;RCOR1;RELA;RXRA;RXRB;RXRG;SIN3A;SP1;SP2;SP3;SP4;SP5;SPI1;SREBF1;SREBF2;STAT5A;TAF1;THRA;THRB;VEZF1;WRNIP1;WT1;ZBED1;ZBTB14;ZBTB17;ZBTB5;ZBTB7B;ZNF148;ZNF202;ZNF212;ZNF214;ZNF263;ZNF341;ZNF398;ZNF432;ZNF444;ZNF467;ZNF496;ZNF529;ZNF543;ZNF596;ZNF615;ZNF701;ZNF737;ZNF740;ZNF875;ZNF880"
  # [20] "ZNF148;ZNF281"

  # [21] "ATF2;IKZF1"

  # [22] "GTF2H3"

  # [23] "ZNF282"

  # [24] "BARX1"

  # [25] "IRF4"

  # [26] "ZKSCAN2"

  # [27] "HMGA1"

  # [28] "STAT6"

  # [29] "LHX9"

  # [30] "PRDM5;ZNF324;ZNF341;ZNF580"

  # [31] "DPF2;GTF2F1;IKZF3;MAZ;NFXL1;PRDM9;RBAK;SIN3A;TAF1;VEZF1;WRNIP1;ZBTB5;ZNF263;ZNF341;ZNF444;ZNF467;ZNF496;ZNF596;ZNF701;ZNF875"

  # [32] "ZNF281"

  # [33] "ATF2"
  length(tmp) # 10


  ## add to the binary annotation matrix the motif-based TF-target annotations
  for (key in keys) {
    motif_TF_highConf <- motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf

    tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold & (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf))
    genes <- unique(unlist(strsplit(tmp$enrichedGenes, ";")))
    mat[, paste0("cisTarget_", key, ".motif_target")] <- ifelse(rownames(mat) %in% genes, "1", "0")
  }

  # 0

  ####################################################
  ### extract subnetworks and add seed_key-bound links ###
  ####################################################
  # literature review to identify topest TF here as a key to followup predictions


  dim(mat) # [1] 54 65  all CTS.8 genes !!
  colnames(mat)

  ## from the enriched motifs find the potential regulators that are either HiG of CTS member of a cluster of interest  ===================
  x <- colnames(mat)[grepl("cisTarget_", colnames(mat))]
  x <- gsub("cisTarget_", "", x) %>% gsub("\\.motif_target", "", .)
  x <- unlist(strsplit(x, ";")) %>% unique()
  x
  #   [1] "PARP1"   "HOXD3"   "ZSCAN9"  "CEBPB"   "PDX1"    "STAT6"   "TCF3"
  #   [8] "POU3F3"  "AHCTF1"  "ISX"     "LHX2"    "LHX9"    "SHOX2"   "FOXM1"
  #  [15] "IKZF1"   "TBL1XR1" "HMGA1"   "HMGA2"   "STAT5A"  "POU4F3"  "RREB1"
  #  [22] "PRDM5"   "ZNF324"  "ZNF324B" "ZNF341"  "ZNF580"  "SF1"     "HMG20A"
  #  [29] "KLF2"    "BCL11A"  "BCL6"    "CHURC1"  "DDX20"   "DPF2"    "EGR1"
  #  [36] "EGR2"    "EP300"   "FLI1"    "GTF2F1"  "HINFP"   "IKZF3"   "IRF4"
  #  [43] "KLF10"   "KLF13"   "KLF15"   "KLF16"   "KLF6"    "KLF9"    "MAZ"
  #  [50] "MTF1"    "MZF1"    "NFXL1"   "OVOL2"   "PATZ1"   "PAX4"    "PML"
  #  [57] "PRDM9"   "PURA"    "RARA"    "RARB"    "RARG"    "RBAK"    "RCOR1"
  #  [64] "RELA"    "RXRA"    "RXRB"    "RXRG"    "SIN3A"   "SP1"     "SP2"
  #  [71] "SP3"     "SP4"     "SP5"     "SPI1"    "SREBF1"  "SREBF2"  "TAF1"
  #  [78] "THRA"    "THRB"    "VEZF1"   "WRNIP1"  "WT1"     "ZBED1"   "ZBTB14"
  #  [85] "ZBTB17"  "ZBTB5"   "ZBTB7B"  "ZNF148"  "ZNF202"  "ZNF212"  "ZNF214"
  #  [92] "ZNF263"  "ZNF398"  "ZNF432"  "ZNF444"  "ZNF467"  "ZNF496"  "ZNF529"
  #  [99] "ZNF543"  "ZNF596"  "ZNF615"  "ZNF701"  "ZNF737"  "ZNF740"  "ZNF875"
  # [106] "ZNF880"  "ZNF281"  "ATF2"    "GTF2H3"  "ZNF282"  "BARX1"   "ZKSCAN2"

  ##  extend seed_TF candidates to those 1) being CTS itself, or 2) CTS_enriched motifs having target genes highly expressed in CP, CM, or CF ===============
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
  # TCF3 is in DEG of 8
  # HMGA2 is in DEG of 17
  # HMGA2 is in DEG of 18
  # HMGA2 is in DEG of 8
  # KLF6 is in DEG of 17
  # KLF6 is in DEG of 18
  # [1] "NR2F1" "IRX5"  "ALX1"  "TCF3"  "HMGA2" "HMGA2" "HMGA2" "KLF6"  "KLF6"


  for (i in x) {
    if (i %in% CTS[[CTS_ID]]) {
      cat(paste0(i, " is in CTS of ", CTS_ID, "\n"))
      key_TFs <- c(key_TFs, i)
    }
  }
  # RARB is in CTS of 8

  key_TFs <- unique(key_TFs)
  key_TFs # [1] "NR2F1" "IRX5"  "ALX1"  "TCF3"  "HMGA2" "KLF6"  "RARB"


  mat <- as.data.frame(mat)

  ##  filter out the key_TFs that come from one motif thus will share target genes for the downstream predictions	=======================
  ## ====== note taht this step need manually watch, such as selecting the TFs of the same motif that are most close to the seed TF ======
  # (x = which(grepl(seed_TF[3], colnames(mat)) | grepl(seed_TF[2], colnames(mat)) | grepl(seed_TF[3], colnames(mat)) ) )  # TF candidateds for CTS.CP
  (x <- intersect(which(grepl("cisTarget_", colnames(mat))), which(Reduce("|", lapply(key_TFs, function(p) grepl(p, colnames(mat), fixed = F))))))
  x
  # [1] 34 40 48
  colnames(mat)[x]
  # [1] "cisTarget_TCF3.motif_target"
  # [2] "cisTarget_HMGA1;HMGA2.motif_target"
  # [3] "cisTarget_BCL11A;BCL6;CHURC1;DDX20;DPF2;EGR1;EGR2;EP300;FLI1;GTF2F1;HINFP;IKZF1;IKZF3;IRF4;KLF10;KLF13;KLF15;KLF16;KLF6;KLF9;MAZ;MTF1;MZF1;NFXL1;OVOL2;PATZ1;PAX4;PML;PRDM9;PURA;RARA;RARB;RARG;RBAK;RCOR1;RELA;RXRA;RXRB;RXRG;SIN3A;SP1;SP2;SP3;SP4;SP5;SPI1;SREBF1;SREBF2;STAT5A;TAF1;THRA;THRB;VEZF1;WRNIP1;WT1;ZBED1;ZBTB14;ZBTB17;ZBTB5;ZBTB7B;ZNF148;ZNF202;ZNF212;ZNF214;ZNF263;ZNF341;ZNF398;ZNF432;ZNF444;ZNF467;ZNF496;ZNF529;ZNF543;ZNF596;ZNF615;ZNF701;ZNF737;ZNF740;ZNF875;ZNF880.motif_target"

  # Manually choose TFs that are present in cisTarget. If two TFs are part of the same cisTarget, choose one.
  # Optionally use automated code from IbarraSoria2018 code.


  x <- NULL
  for (j in key_TFs) {
    y <- intersect(which(grepl("cisTarget_", colnames(mat))), which(Reduce("|", lapply(j, function(p) grepl(p, colnames(mat), fixed = F)))))
    cat(j, "\t", y, "\t", colnames(mat)[y], "\n")
    if (length(y) == 0) y <- 0
    x <- c(x, y[1])
  }
  names(x) <- key_TFs
  x
  # NR2F1  IRX5  ALX1  TCF3 HMGA2  KLF6  RARB
  #     0     0     0    34    40    48    48


  key_TFs <- c("TCF3", "HMGA2", "KLF6", "RARB")
  names(x) <- key_TFs
  if (length(x) > 0) {
    for (j in seq_along(x)) {
      # if(grepl(";",colnames(mat)[x[j]])){
      #     tmp = strsplit(colnames(mat)[x[j]], ";")[[1]] %>% gsub('\\.motif_target', '', .) %>% gsub('^cisTarget_', '', .)
      #     names(x)[j] = intersect(seed_TF, tmp) %>% unique
      #     } else {
      # names(x)[j] = intersect(seed_TF, colnames(mat)[x[j]])
      # }

      key <- names(x)[j]
      mat[, paste0(key, "_CP_candidate")] <- ifelse(mat[, "CP_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
      # 3) must be open at PCW8_CM or PCW19_CM for CM_hi
      mat[, paste0(key, "_CM_candidate")] <- ifelse(mat[, "CM_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
      # 4) must be open at PCW8_CF or PCW19_CF for CF_hi
      mat[, paste0(key, "_CF_candidate")] <- ifelse(mat[, "CF_hi"] == 1 & mat[, x[j]] == 1, 1, 0)
    } # end of for(j in seq_alongs(x))
  } # end of if(length(x)>0)

  dim(mat) # 54 77

  cat(paste0("key_TFs: ", paste(key_TFs, collapse = "_"), "\n")) # key_TFs: TCF3_HMGA2_KLF6_RARB

  ## record the whole mat for all CTS genes, which we can load later  ===================
  if (length(key_TFs) > 0) {
    fileName <- paste0("heatmap_blocked_", CTS_name, "_scATAC_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv")
    write.table(mat, file = fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) # !!!!!!!!!!!!

    saveRDS(list(x = x, key_TFs = key_TFs, motifAnnot_sub = motifAnnot_sub), "scATAC_cistTarget_variables.rds")
  } else {
    stop("No key TFs found for ", CTS_name, "\n")
  }
}

###  heatmap confirming key TF’s self impact by checking its targets among the CTS_CP ================
# key_TFs =  c('TCF3','HMGA2','RARB')
# fileName=paste0('heatmap_blocked_',CTS_name,'_scATAC_cisTarget_',paste(key_TFs, collapse='_'),'_v3.tsv')
# mat = read.table(fileName, sep='\t',header=T, check.names=F)
# x = c(  34 ,40, 48)
# names(x) = key_TFs


dim(mat) # [1] 54 77
colnames(mat)


motif_TF_highConf <- gsub("cisTarget_", "", colnames(mat)[x]) %>% gsub("\\.motif_target", "", .)
print(motif_TF_highConf)
# [1] "TCF3"
# [2] "HMGA1;HMGA2"
# [3] "BCL11A;BCL6;CHURC1;DDX20;DPF2;EGR1;EGR2;EP300;FLI1;GTF2F1;HINFP;IKZF1;IKZF3;IRF4;KLF10;KLF13;KLF15;KLF16;KLF6;KLF9;MAZ;MTF1;MZF1;NFXL1;OVOL2;PATZ1;PAX4;PML;PRDM9;PURA;RARA;RARB;RARG;RBAK;RCOR1;RELA;RXRA;RXRB;RXRG;SIN3A;SP1;SP2;SP3;SP4;SP5;SPI1;SREBF1;SREBF2;STAT5A;TAF1;THRA;THRB;VEZF1;WRNIP1;WT1;ZBED1;ZBTB14;ZBTB17;ZBTB5;ZBTB7B;ZNF148;ZNF202;ZNF212;ZNF214;ZNF263;ZNF341;ZNF398;ZNF432;ZNF444;ZNF467;ZNF496;ZNF529;ZNF543;ZNF596;ZNF615;ZNF701;ZNF737;ZNF740;ZNF875;ZNF880"
# [4] "BCL11A;BCL6;CHURC1;DDX20;DPF2;EGR1;EGR2;EP300;FLI1;GTF2F1;HINFP;IKZF1;IKZF3;IRF4;KLF10;KLF13;KLF15;KLF16;KLF6;KLF9;MAZ;MTF1;MZF1;NFXL1;OVOL2;PATZ1;PAX4;PML;PRDM9;PURA;RARA;RARB;RARG;RBAK;RCOR1;RELA;RXRA;RXRB;RXRG;SIN3A;SP1;SP2;SP3;SP4;SP5;SPI1;SREBF1;SREBF2;STAT5A;TAF1;THRA;THRB;VEZF1;WRNIP1;WT1;ZBED1;ZBTB14;ZBTB17;ZBTB5;ZBTB7B;ZNF148;ZNF202;ZNF212;ZNF214;ZNF263;ZNF341;ZNF398;ZNF432;ZNF444;ZNF467;ZNF496;ZNF529;ZNF543;ZNF596;ZNF615;ZNF701;ZNF737;ZNF740;ZNF875;ZNF880"

for (key in motif_TF_highConf) { # !!!!!!!!
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
    error = function(e) { message("heatmap skipped for '", paste(key_in_TFfamily, collapse=";"), "': ", e$message); NULL }
  )
  ## for Fig 5, only the PRRX1 CP bound genes are shown
  if (!is.null(p)) {
    pdf(file = paste0("heatmap_blocked_", CTS_name, "_scATAC_cisTarget_", key_in_TFfamily, "_v3_coding_target.pdf"), height = 4)
    print(p)
    dev.off()
  }
}

# candidate genes:  1
# direct motif:  cisTarget_TCF3.motif_target  used from multiple potentoal matches
# candidate genes:  5
# simplified motif name:  cisTarget_HMGA1;HMGA2.motif_target
# candidate genes:  1
# simplified motif name:  cisTarget_BCL11A;BCL6;CHURC1;DDX20;DPF2;EGR1;EGR2;EP300;FLI1;GTF2F1;HINFP;IKZF1;IKZF3;IRF4;KLF10;KLF13;KLF15;KLF16;KLF6;KLF9;MAZ;MTF1;MZF1;NFXL1;OVOL2;PATZ1;PAX4;PML;PRDM9;PURA;RARA;RARB;RARG;RBAK;RCOR1;RELA;RXRA;RXRB;RXRG;SIN3A;SP1;SP2;SP3;SP4;SP5;SPI1;SREBF1;SREBF2;STAT5A;TAF1;THRA;THRB;VEZF1;WRNIP1;WT1;ZBED1;ZBTB14;ZBTB17;ZBTB5;ZBTB7B;ZNF148;ZNF202;ZNF212;ZNF214;ZNF263;ZNF341;ZNF398;ZNF432;ZNF444;ZNF467;ZNF496;ZNF529;ZNF543;ZNF596;ZNF615;ZNF701;ZNF737;ZNF740;ZNF875;ZNF880.motif_target

##########################################################
## --  identify_TF_targeted_pull_candidat -- subset of CTS[['CP']] that are
## exclusively highly expressed (HiG) in CM (or CF)
## further narrow the subset to be open at CP or CM (or CF) ======================
library(BioNet)
packageVersion("BioNet") # '1.70.0'
library(igraph)
library(tibble)


mat <- read.table(paste0("heatmap_blocked_", CTS_name, "_scATAC_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv"), sep = "\t", header = T, check.names = FALSE)
dim(mat) # [1] 54  77

# make graph_list compatible with older helper code
if (!"HiG8" %in% names(graph_list) && "HiG_8" %in% names(graph_list)) {
  graph_list[["HiG8"]] <- graph_list[["HiG_8"]]
}
if (!"8" %in% names(graph_list) && "CTS_8" %in% names(graph_list)) {
  graph_list[["8"]] <- graph_list[["CTS_8"]]
}

graph_list_fixed <- graph_list

graph_list_fixed[["HiG8"]] <- graph_list_fixed[["HiG_8"]]
graph_list_fixed[["8"]] <- graph_list_fixed[["CTS_8"]]

# Had to redeine the indentiy_TF_targeted_pull_candidate function to make it compatible with the graph_list_fixed
identify_TF_targeted_pull_candidate <- function(
  mat, graph_list, CTS_name, CHD,
  key, keep_selfloop = TRUE,
  TF_bound_column_name,
  TF_appendix = NULL,
  edge_colored_by_Maven2023_ISL1KO = FALSE,
  key_in_TFfamily = key
) {
  CTS_ID <- sub("^CTS_", "", CTS_name)

  four_names <- c(
    "CTSHiG.CP_TF.target",
    "CTS.CP_TF.target_HiGCM",
    "CTS.CP_TF.target_HiGCF",
    "CTSHiG.CP_TF.target_CPopen"
  )

  graph_TF_list <- list(
    graph_list[[paste0("HiG_", CTS_ID)]],
    graph_list[[paste0("CTS_", CTS_ID)]],
    graph_list[[paste0("CTS_", CTS_ID)]],
    graph_list[[paste0("CTS_", CTS_ID)]]
  )
  names(graph_TF_list) <- four_names

  if (any(vapply(graph_TF_list, is.null, logical(1)))) {
    stop("Missing graph(s): ", paste(names(graph_TF_list)[vapply(graph_TF_list, is.null, logical(1))], collapse = ", "))
  }

  TF_target <- rownames(mat)[which(mat[, TF_bound_column_name] == 1)]
  length(TF_target)

  if (is.null(TF_appendix)) {
    GRN_node_list <- list(
      rownames(mat)[which(mat[, "CP_candidate"] == 1)],
      rownames(mat)[which(mat[, "CM_candidate"] == 1)],
      rownames(mat)[which(mat[, "CF_candidate"] == 1)],
      rownames(mat)[which(mat[, "CP_candidate"] == 1)]
    )
  } else {
    GRN_node_list <- list(
      CTSHiG.CP_TF.target = rownames(mat)[which(mat[, paste0(key, "_CP_candidate")] == 1)],
      CTS.CP_TF.target_HiGCM = rownames(mat)[which(mat[, paste0(key, "_CM_candidate")] == 1)],
      CTS.CP_TF.target_HiGCF = rownames(mat)[which(mat[, paste0(key, "_CF_candidate")] == 1)],
      CTSHiG.CP_TF.target_CPopen = rownames(mat)[which(mat[, paste0(key, "_CP_candidate")] == 1)]
    )
  }
  names(GRN_node_list) <- four_names

  for (i in 2:length(GRN_node_list)) {
    if (grepl("CM", names(GRN_node_list)[i])) {
      GRN_node_list[[i]] <- intersect(
        GRN_node_list[[i]],
        rownames(mat)[which(mat[, "PCW6CP_access"] == 1 | mat[, "PCW8_CM_access"] == 1 | mat[, "PCW19_CM_access"] == 1)]
      )
    } else if (grepl("CF", names(GRN_node_list)[i]) | grepl("SMC", names(GRN_node_list)[i])) {
      GRN_node_list[[i]] <- intersect(
        GRN_node_list[[i]],
        rownames(mat)[which(
          mat[, "PCW6CP_access"] == 1 |
            mat[, "PCW8_CF_access"] == 1 |
            mat[, "PCW19_CF_access"] == 1 |
            mat[, "PCW8_SMC_access"] == 1 |
            mat[, "PCW19_SMC_access"] == 1 |
            mat[, "iEPC_access"] == 1
        )]
      )
    } else {
      GRN_node_list[[i]] <- intersect(
        GRN_node_list[[i]],
        rownames(mat)[which(mat[, "PCW6CP_access"] == 1)]
      )
    }
  }

  mat <- as.data.frame(apply(mat, 2, as.numeric), row.names = rownames(mat))

  for (name in names(graph_TF_list)) {
    graph <- graph_TF_list[[name]]

    V(graph)$name <- toupper(V(graph)$name)
    graph <- induced_subgraph(graph, intersect(GRN_node_list[[name]], V(graph)$name))

    V(graph)$color <- ifelse(V(graph)$name %in% CHD, "red", "lightgrey")
    E(graph)$color <- "grey70"
    E(graph)$lty <- "solid"

    TF_target_in_graph <- intersect(TF_target, V(graph)$name)
    if (!keep_selfloop) TF_target_in_graph <- setdiff(TF_target_in_graph, key_in_TFfamily)

    if (length(TF_target_in_graph) > 0) {
      g_check <- if (is_directed(graph)) as.undirected(graph, mode = "collapse") else graph

      flag <- !(key_in_TFfamily %in% V(graph)$name)
      g_check <- add_vertex_if_missing(g_check, key_in_TFfamily)

      if (is.null(vertex_attr(g_check, "color"))) V(g_check)$color <- rep("lightgrey", vcount(g_check))
      if (flag) V(g_check)[name == key_in_TFfamily]$color <- "yellow"

      vp_pairs <- as.vector(rbind(key_in_TFfamily, TF_target_in_graph))
      eid <- get_edge_ids(g_check, vp = vp_pairs)
      targets_to_add <- TF_target_in_graph[eid == 0]

      if (length(targets_to_add) > 0) {
        m_before <- ecount(g_check)
        edges_vec <- as.vector(rbind(key_in_TFfamily, targets_to_add))
        g_check <- add_edges(g_check, edges_vec)
        new_edges <- seq.int(m_before + 1, ecount(g_check))

        if (is.null(edge_attr(g_check, "color"))) E(g_check)$color <- rep("grey70", ecount(g_check))
        if (is.null(edge_attr(g_check, "lty"))) E(g_check)$lty <- rep("solid", ecount(g_check))

        new_edges_color <- rep("grey70", length(targets_to_add))
        if (edge_colored_by_Maven2023_ISL1KO) {
          up_cols <- grep("^Maven2023_gene_ISL1_up", colnames(mat))
          if (length(up_cols) > 0) {
            new_edges_color[rowSums(mat[targets_to_add, up_cols, drop = FALSE]) > 0] <- "red"
          }
          dn_cols <- grep("^Maven2023_gene_ISL1_dn", colnames(mat))
          if (length(dn_cols) > 0) {
            new_edges_color[rowSums(mat[targets_to_add, dn_cols, drop = FALSE]) > 0] <- "blue"
          }
        }

        E(g_check)$color[new_edges] <- new_edges_color
        E(g_check)$lty[new_edges] <- "dashed"
      }

      graph <- g_check
    }

    graph_TF_list[[name]] <- graph
  }

  return(graph_TF_list)
}

key <- "TCF3"
key_column <- which(grepl(key, colnames(mat), fixed = TRUE) & grepl("cisTarget_", colnames(mat)))

graph_TF_list <- identify_TF_targeted_pull_candidate(
  mat, graph_list, CTS_name, CHD,
  key = key,
  keep_selfloop = TRUE,
  TF_bound_column_name = key_column,
  TF_appendix = key,
  edge_colored_by_Maven2023_ISL1KO = FALSE,
  key_in_TFfamily = key
)

for (key in key_TFs) {
  key_column <- which(grepl(key, colnames(mat)) & grepl("cisTarget_", colnames(mat)))
  if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key

  ## further narrow the subset to be open at CP  or descendants
  graph_TF_list <- identify_TF_targeted_pull_candidate(mat, graph_list_fixed, CTS_name, CHD,
    key = key,
    keep_selfloop = TRUE, # whether to keep the self-loop of the key
    TF_bound_column_name = key_column,
    TF_appendix = key,
    edge_colored_by_Maven2023_ISL1KO = FALSE,
    key_in_TFfamily = key_in_TFfamily
  )
  # saveRDS(graph_TF_list, file='PPI_graph_PRRX1_GRN_prediction.rds')

  saveRDS(graph_TF_list, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
}

names(graph_TF_list)
# [1] "CTSHiG.CP_TF.target"        "CTS.CP_TF.target_HiGCM"
# [3] "CTS.CP_TF.target_HiGCF"     "CTSHiG.CP_TF.target_CPopen"


for (key in key_TFs) {
  if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key

  graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
  plot_TF_targeted_pull_candidate(graph_TF_list, key_in_TFfamily, CTS_name, saveFigure = TRUE)
}
# => PPI_graph_GATA4_GRN_prediction_CTS_CP_v3.pdf;


##########################################################
#### =quantify the edge weight changes and validate by Maven Isk1_KO in CP results ===================

# ## TO get teh cluster ID, remove the .1, .2, etc. from the cluster names which indicating additional CTS identified from the same cluster
# CM_ID = paste0('HiG_',CM_cluster) %>% gsub('\\.[1-9]','',.)
# CF_ID = paste0('HiG_',CF_cluster) %>% gsub('\\.[1-9]','',.)

# ## heler to paly around the layout to be beautful
make_layout <- function(g, seed) {
  set.seed(seed)
  layout_with_fr(g)
}

for (key in key_TFs) {
  if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key

  graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
  names(graph_TF_list)
  # [1] "CTSHiG.CP_TF.target"        "CTS.CP_TF.target_HiGCM"     "CTS.CP_TF.target_HiGCF"
  # [4] "CTSHiG.CP_TF.target_CPopen"
  edge_attr_names(graph_TF_list[["CTS.CP_TF.target_HiGCM"]])
  # [1] "weight"         "corexp_sign"    "coexp_target"   "norm_PPI_score" "color"          "lty"

  res <- fill_TF_targeting_predicted_edges(graph_TF_list,
    linkage_name = "CM", graph_list,
    sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
    descendant_cluster_id = CM_cluster, TF_symbol = key_in_TFfamily,
    HVG = rownames(sce)
  )
  cat(paste0(key, ' CM vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), "\n")) # 2
  if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CM_final.rds"))


  res <- fill_TF_targeting_predicted_edges(graph_TF_list,
    linkage_name = "CF", graph_list,
    sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
    descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
    HVG = rownames(sce)
  )
  cat(paste0(key, ' CF vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), "\n")) #  0
  if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CF_final.rds"))
  names(res)
  # [1] "g_CT_sub"              "g_descendant_sub"
}
# TCF3 CM vcount(res[["g_CT_sub"]]): 0
# TCF3 CF vcount(res[["g_CT_sub"]]): 0
# HMGA2 CM vcount(res[["g_CT_sub"]]): 3
# HMGA2 CF vcount(res[["g_CT_sub"]]): 2
# KLF6 CM vcount(res[["g_CT_sub"]]): 2
# KLF6 CF vcount(res[["g_CT_sub"]]): 0
# RARB CM vcount(res[["g_CT_sub"]]): 2
# RARB CF vcount(res[["g_CT_sub"]]): 0


### =reporting the results	and visualization #==========
(files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds"))
# [1] "PPI_graph_HMGA2_GRN_prediction_CTS_8_CF_final.rds"
# [2] "PPI_graph_HMGA2_GRN_prediction_CTS_8_CM_final.rds"
# [3] "PPI_graph_KLF6_GRN_prediction_CTS_8_CM_final.rds"
# [4] "PPI_graph_RARB_GRN_prediction_CTS_8_CM_final.rds"


# for CF-pull with results, reporting and plotting

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
    predict <- prioritize_edge_change(g1, edge_change_df = change_df, top_n = 5, title = paste0(pull, "-pull subnetwork_", key_in_TFfamily))
    # =>TIPS_delta_edge_reweighting_CF-pull subnetwork_ETS2.pdf
  } else {
    cat("No edges in ", pull, "-pull subnetwork\n")
  }
  # No edges in CF-pull subnetwork


  x <- paste0(key, "_motif_target")
  (y <- motifAnnot_sub[which(motifAnnot_sub$regulators %in% key | grepl(key, motifAnnot_sub$regulators)), ])
  motif_TF_highConf <- y$motif_TF_highConf
  tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold & (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf))
  all_regulators <- strsplit(y$regulators, ";", fixed = T) %>%
    unlist() %>%
    unique()

  change_df <- cbind(linkage = pull, change_df, TF_highConf = tmp$TF_highConf, motif = tmp$motif, NES = tmp$NES)
  change_df$TF_highConf[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
  change_df$motif[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
  change_df$NES[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""

  final_table <- rbind(final_table, change_df)
}

dim(final_table) # 5  13
####### generateing final table of the predicted subnetwork #################

write.table(final_table,
  file = paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

## verify the unchanged edges in CF is due to 'shrinkage' correlation effect which is more robust
key_in_TFfamily <- "HMGA2"
graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
res <- fill_TF_targeting_predicted_edges(graph_TF_list,
  linkage_name = "CF", graph_list,
  sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
  descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
  HVG = rownames(sce),
  shrink = FALSE
)
E(res[[2]])$weight # [1] 0.02714859
E(res[[1]])$weight # [1] 0

res <- fill_TF_targeting_predicted_edges(graph_TF_list,
  linkage_name = "CF", graph_list,
  sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
  descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
  HVG = rownames(sce),
  shrink = TRUE
)
E(res[[2]])$weight # [1] 0
E(res[[1]])$weight # [1] 0
