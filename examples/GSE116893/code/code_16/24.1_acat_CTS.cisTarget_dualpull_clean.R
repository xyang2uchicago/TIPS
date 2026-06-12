# Set True if running 24.0 the first time
rebuild_mat <- TRUE
source("C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/code/code_16/24.0_acat_load_input_clean.R")

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

## check the loaded objects =========================
seed_TF
names(graph_list)
#  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"    
#  [7] "HiG_7"     "HiG_8"     "HiG_10"    "HiG_11"    "HiG_12"    "HiG_13"   
# [13] "HiG_14"    "HiG_16"    "HiG_17"    "HiG_15"    "HiG_9"     "HiGCTS_15"
# [19] "HiGCTS_16" "CTS_15"    "CTS_16"    "CTS_9" 
names(DEG)
#  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "10" "11" "12" "13" "14" "16" "17"
# [16] "15" "9" 

celltype_col # [1] "leiden_0.6"
CP_cluster # '16'
CM_cluster # '7'
CF_cluster # '9'

lengths(CTS) # a list
#  15  16   9
#  76 139 110

class(sce) # [1] "SingleCellExperiment"

dim(mat) # [1] 139 4
colnames(mat)
# [1] "CP_hi"  "CM_hi"  "CF_hi"  "CTS_16"

seed_TF
CTS_ID # '16'
CTS_name # 'CTS_16'

## set subfold =======================
(updir <- getwd())
# "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/results/GSE181346_heart_scATAC"
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

(dim(mat)) # [1] 76  4
### = load the binary annotation matrix build in code 24.0xxx ============================
files <- list.files(pattern = "^heatmap_blocked_CTS_") %>% grep("_v3.tsv", ., value = TRUE)
if (length(files) > 0) {
    fileName <- grep("_v3.tsv", files, value = TRUE)
    mat <- read.table(fileName, sep = "\t", header = T, check.names = FALSE)

    saved_variables <- readRDS(file = "NB_cisTarget_variables.rds")

    # Recreate variables
    x <- saved_variables$x
    key_TFs <- saved_variables$key_TFs
    motifAnnot_sub <- saved_variables$motifAnnot_sub
} else {
    #### focus on the this CTS member gene-enriched motifs, with NES_threshold control =========================
    motifAnnot_sub <- get_regulators_from_motifs(cisTarget.res, CTS_ID, NES_threshold, motifAnnot = motifAnnot)
    keys <- unique(motifAnnot_sub$regulators)
    if (any(is.na(keys))) keys <- keys[!is.na(keys)]

    ## add to the binary annotation matrix the motif-based TF-target annotations
    for (key in keys) {
        motif_TF_highConf <- motifAnnot_sub[match(key, motifAnnot_sub$regulators), ]$motif_TF_highConf

        tmp <- subset(cisTarget.res, geneSet == CTS_ID & NES >= NES_threshold & (motif == motif_TF_highConf | TF_highConf == motif_TF_highConf))
        genes <- unique(unlist(strsplit(tmp$enrichedGenes, ";")))
        mat[, paste0("cisTarget_", key, ".motif_target")] <- ifelse(rownames(mat) %in% genes, "1", "0")
    }

    ####################################################
    ### extract subnetworks and add seed_key-bound links ###
    ####################################################
    # literature review to identify top TF here as a key to followup predictions


    (dim(mat)) # [1] 139 154
    (colnames(mat))
    # ...
    #   [9] "cisTarget_ZNF730.motif_target"
    #  [31] "cisTarget_E2F3.motif_target"
    #  [44] "cisTarget_PLAGL1.motif_target"
    #  [48] "cisTarget_E2F1;E2F3;HES1;HES5;HES7;HEY1;HEY2;MAX;MYCN;TFAP2A;TFAP2B;TFAP2C;TFAP2E.motif_target"
    #  [67] "cisTarget_MTF2.motif_target"
    #  [74] "cisTarget_EBF1.motif_target"
    #  [97] "cisTarget_ZNF704.motif_target"
    # [113] "cisTarget_GATA1;GATA2;GATA3;GATA4;GATA5;GATA6;MZF1;NEUROD1;NFIA;NFIC;NR3C1;PAX7;TAF1;TAF1L;THAP1;YY1;YY2;ZFP42;ZNF597.motif_target"
    # [119] "cisTarget_DNMT1;IRF4;IRF5;IRF6;JUND;NEUROD1;PHF1;PRDM15;PRDM5;SP2;SPZ1;WT1;ZBTB20;ZFP64;ZNF17;ZNF296;ZNF350;ZNF431;ZNF486;ZNF574;ZNF589;ZNF891.motif_target"
    # [140] "cisTarget_CHD1.motif_target"
    # [154] "cisTarget_RARA.motif_target"


    ## from the enriched motifs find the potential regulators that are either HiG of CTS member of a cluster of interest  ===================
    x <- colnames(mat)[grepl("cisTarget_", colnames(mat))]
    x <- gsub("cisTarget_", "", x) %>% gsub("\\.motif_target", "", .)
    x <- unlist(strsplit(x, ";")) %>% unique()
    print(x)
    #   [1] "SREBF2"  "WT1"     "YY1"     "E2F2"    "ZBTB20"  "ZFP64"   "CHURC1"
    #   [8] "ZNF730"  "EGR1"    "KLF1"    "KLF10"   "KLF12"   "KLF13"   "KLF14"
    #  [43] "EGR4"    "ZNF311"  "GTF3C2"  "SP7"     "ZNF142"  "BACH1"   "BACH2"
    # [85] "BCL6"    "CHD1"    "CPSF4"   "E2F3"    "E2F4"    "EGR2"    "FOXO3"
    # [141] "TFAP2C"  "ZIC4"    "ZIC5"    "HES1"    "HES5"    "HES7"    "HEY1"
    # [148] "HEY2"    "MAX"     "MYCN"    "TFAP2E"  "CNOT3"   "TFDP1"   "CTCF"
    # [225] "GATA2"   "GATA3"   "GATA4"   "GATA5"   "GATA6"   "NFIA"    "NR3C1"
    # [232] "PAX7"    ...  "JUND" ... "EBF1" ... "NFIB" ... "ZNF704"
    # (270 total unique TF names)

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
    # SREBF2 is in DEG of 16
    # ZBTB20 is in DEG of 9
    # ZBTB20 is in DEG of 16
    # KLF12 is in DEG of 7
    # KLF12 is in DEG of 9
    # KLF7 is in DEG of 9
    # BACH2 is in DEG of 16
    # CHD1 is in DEG of 16
    # PLAGL1 is in DEG of 16
    # MYCN is in DEG of 7
    # EBF1 is in DEG of 7
    # NFIB is in DEG of 7
    # NFIC is in DEG of 16
    # ZBTB21 is in DEG of 16
    # ZNF704 is in DEG of 7
    # ZNF141 is in DEG of 16
    # GATA6 is in DEG of 16
    # NFIA is in DEG of 16
    # JUND is in DEG of 16
    # CHD1 is in CTS of 16


    for (i in x) {
        if (i %in% CTS[[CTS_ID]]) {
            cat(paste0(i, " is in CTS of ", CTS_ID, "\n"))
            key_TFs <- c(key_TFs, i)
        }
    }

    key_TFs <- unique(key_TFs)
    (key_TFs)
    # [1] "SREBF2" "ZBTB20" "KLF12"  "KLF7"   "BACH2"  "CHD1"   "PLAGL1" "MYCN"
    # [9] "EBF1"   "NFIB"   "NFIC"   "ZBTB21" "ZNF704" "ZNF141" "GATA6"  "NFIA"
    # [17] "JUND"


    mat <- as.data.frame(mat)

    ##  filter out the key_TFs that come from one motif thus will share target genes for the downstream predictions	=======================
    ## ====== note that this step is manual, such as selecting the TFs of the same motif that are most close to the seed TF ======
    # (x = which(grepl(seed_TF[3], colnames(mat)) | grepl(seed_TF[2], colnames(mat)) | grepl(seed_TF[3], colnames(mat)) ) )  # TF candidateds for CTS.CP
    (x <- intersect(which(grepl("cisTarget_", colnames(mat))), which(Reduce("|", lapply(key_TFs, function(p) grepl(p, colnames(mat), fixed = F))))))
    (x)
    #  [1]   5   7  10  16  21  26  34  40  44  48  49  54  56  72  74  76  77  79  81
    # [20]  86  96  97 103 108 113 118 119 124 127 132 140 142 149 150 153
    (colnames(mat)[x])
    #  [1] "cisTarget_SREBF2;WT1;YY1.motif_target"
    #  [2] "cisTarget_ZBTB20;ZFP64.motif_target"
    #  [3] "cisTarget_EGR1;KLF1;KLF10;KLF12;KLF13;...;ZNF740.motif_target"   (shared, col 10 = KLF12 & KLF7)
    #  [4] "cisTarget_BACH1;BACH2;...;ZNF740.motif_target"                     (shared, col 16 = BACH2)
    #  [5] "cisTarget_BCL6;CHD1;...;ZNF880.motif_target"                       (shared, col 21 = CHD1 fallback)
    #  [6] "cisTarget_CHURC1;E2F4;...;ZXDC.motif_target"                       (shared, col 26 = KLF12/KLF7 overlap)
    #  [9] "cisTarget_PLAGL1.motif_target"                                      (dedicated, col 44)
    # [10] "cisTarget_E2F1;E2F3;HES1;HES5;HES7;HEY1;HEY2;MAX;MYCN;...motif_target" (shared, col 48 = MYCN)
    # [14] "cisTarget_EBF1;MZF1;STAT1.motif_target"                            (col 72, EBF1 partial)
    # [15] "cisTarget_EBF1.motif_target"                                        (dedicated, col 74)
    # [16] "cisTarget_CHD1;ZGPAT;ZNF202.motif_target"                          (col 76, CHD1 partial)
    # [19] "cisTarget_DIDO1;HINFP;ILF3;MAZ;NFIB;NFIC;WT1;ZBTB21;...motif_target" (shared, col 81 = NFIB/NFIC/ZBTB21)
    # [22] "cisTarget_ZNF704.motif_target"                                      (dedicated, col 97)
    # [24] "cisTarget_E2F1;E2F4;E2F6;E2F7;EGR2;...;ZNF398.motif_target"        (shared, col 108 = ZNF141)
    # [25] "cisTarget_GATA1;...;GATA6;...;NFIA;NFIC;...;PAX7;...motif_target"  (shared, col 113 = GATA6/NFIA)
    # [26] "cisTarget_CLOCK;DLX3;ERF;ETV2;...;MYCN;...motif_target"             (shared, col 118 = MYCN alt)
    # [27] "cisTarget_DNMT1;IRF4;IRF5;IRF6;JUND;...;ZBTB20;ZFP64;...motif_target" (shared, col 119 = JUND)
    # [31] "cisTarget_CHD1.motif_target"                                        (dedicated, col 140)

    # Manually choose TFs that are present in cisTarget. If two TFs are part of the same cisTarget, choose one.
    # Per-TF loop (IbarraSoria2018 pattern): prefers dedicated single-TF motif; falls back to first partial match.
    x <- NULL
    for (j in key_TFs) {
        y_exact <- which(colnames(mat) == paste0("cisTarget_", j, ".motif_target"))
        if (length(y_exact) > 0) {
            y <- y_exact
        } else {
            y <- intersect(which(grepl("cisTarget_", colnames(mat))),
                           which(grepl(j, colnames(mat), fixed = FALSE)))
        }
        cat(j, "\t", y, "\t", colnames(mat)[y], "\n")
        if (length(y) == 0) y <- 0
        x <- c(x, y[1])
    }
    names(x) <- key_TFs
    x <- x[which(x > 0)]
    key_TFs <- names(x)

    # Per-TF lookup output (j | col_index | col_name):
    # SREBF2    5   cisTarget_SREBF2;WT1;YY1.motif_target
    # ZBTB20    7   cisTarget_ZBTB20;ZFP64.motif_target
    # KLF12    10   cisTarget_EGR1;KLF1;KLF10;KLF12;...;ZNF740.motif_target   (same col as KLF7)
    # KLF7     10   cisTarget_EGR1;KLF1;KLF10;KLF12;...;ZNF740.motif_target   (redundant with KLF12)
    # BACH2    16   cisTarget_BACH1;BACH2;...;ZNF740.motif_target
    # CHD1    140   cisTarget_CHD1.motif_target                                (dedicated; CHD1 in CTS_16 gene set)
    # PLAGL1   44   cisTarget_PLAGL1.motif_target                              (dedicated)
    # MYCN     48   cisTarget_E2F1;E2F3;HES1;...;MYCN;...motif_target
    # EBF1     74   cisTarget_EBF1.motif_target                                (dedicated; DEG of CM cluster 7)
    # NFIB     81   cisTarget_DIDO1;HINFP;ILF3;MAZ;NFIB;NFIC;WT1;ZBTB21;...  (same col as NFIC, ZBTB21)
    # NFIC     81   (redundant with NFIB)
    # ZBTB21   81   (redundant with NFIB)
    # ZNF704   97   cisTarget_ZNF704.motif_target                              (dedicated; DEG of CM cluster 7)
    # ZNF141  103   cisTarget_CDC5L;EGR1;EGR2;...;ZNF141;...motif_target
    # GATA6   113   cisTarget_GATA1;...;GATA6;...;NFIA;NFIC;...;PAX7;...      (same col as NFIA; PAX7 in this motif)
    # NFIA    113   (redundant with GATA6; pick NFIA — more NCC-relevant than GATA6)
    # JUND    119   cisTarget_DNMT1;IRF4;IRF5;IRF6;JUND;...;ZBTB20;...        (AP-1/MES; closest to JUNB from paper)
    #
    # Chosen key_TFs (paper Aim 2 Fig 6A: progenitor=PAX, proliferation=MYC, MES=JUNB):
    #   CHD1  : dedicated motif; in CTS_16 gene set; chromatin remodeler — CTS-intrinsic
    #   MYCN  : MYC-family proliferation module; DEG of CM (cluster 7, ADRN)
    #   EBF1  : dedicated motif; DEG of CM (cluster 7, ADRN); sympathetic lineage marker
    #   JUND  : AP-1/MES module; DEG of CP (cluster 16); closest available proxy for JUNB
    key_TFs <- c("CHD1", "MYCN", "EBF1", "JUND")
    x <- x[key_TFs]
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
        } # end of for(j in seq_along(x))
    } # end of if(length(x)>0)

    dim(mat) # 76 99

    cat(paste0("key_TFs: ", paste(key_TFs, collapse = "_"), "\n")) # key_TFs: TCF3_HMGA2_RARB

    ## record the whole mat for all CTS genes, which we can load later  ===================
    if (length(key_TFs) > 0) {
        fileName <- paste0("heatmap_blocked_", CTS_name, "_NB_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv")
        write.table(mat, file = fileName, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) # !!!!!!!!!!!!

        saveRDS(list(x = x, key_TFs = key_TFs, motifAnnot_sub = motifAnnot_sub), "NB_cisTarget_variables.rds")
    } else {
        stop("No key TFs found for ", CTS_name, "\n")
    }
}

###  heatmap confirming key TF’s self impact by checking its targets among the CTS_CP ================
# key_TFs =  c(‘CHD1’,’MYCN’,’EBF1’,’JUND’)
# fileName=paste0(‘heatmap_blocked_’,CTS_name,’_NB_cisTarget_’,paste(key_TFs, collapse=’_’),’_v3.tsv’)
# mat = read.table(fileName, sep=’\t’,header=T, check.names=F)
# x = c(140, 48, 74, 119)
# names(x) = key_TFs


print(dim(mat)) # [1] 139 166
print(colnames(mat))
#   [1] "CP_hi"
#   [2] "CM_hi"
#   [3] "CF_hi"
#   [4] "CTS_16"
#   ... cisTarget motif columns [5..154] ...
#   "CHD1_CP_candidate"
#   "CHD1_CM_candidate"
#   "CHD1_CF_candidate"
#   "MYCN_CP_candidate"
#   "MYCN_CM_candidate"
#   "MYCN_CF_candidate"
#   "EBF1_CP_candidate"
#   "EBF1_CM_candidate"
#   "EBF1_CF_candidate"
#   "JUND_CP_candidate"
#   "JUND_CM_candidate"
#   "JUND_CF_candidate"

motif_TF_highConf <- gsub("cisTarget_", "", colnames(mat)[x]) %>% gsub("\\.motif_target", "", .)
print(motif_TF_highConf)
# x = c(CHD1=140, MYCN=48, EBF1=74, JUND=119)
# [1] "CHD1"
# [2] "E2F1;E2F3;HES1;HES5;HES7;HEY1;HEY2;MAX;MYCN;TFAP2A;TFAP2B;TFAP2C;TFAP2E"
# [3] "EBF1"
# [4] "DNMT1;IRF4;IRF5;IRF6;JUND;NEUROD1;PHF1;PRDM15;PRDM5;SP2;SPZ1;WT1;ZBTB20;ZFP64;ZNF17;ZNF296;ZNF350;ZNF431;ZNF486;ZNF574;ZNF589;ZNF891"

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
            show_SMC_access = FALSE
        ),
        error = function(e) {
            message("heatmap skipped for '", key_in_TFfamily, "' (0 candidate genes): ", e$message)
            NULL
        }
    )
    if (!is.null(p)) {
        pdf(file = paste0("heatmap_blocked_", CTS_name, "_NB_cisTarget_", key_in_TFfamily, "_v3_coding_target.pdf"), height = 4)
        print(p)
        dev.off()
    }
}
# => heatmap PDFs saved for any key_TFs with >0 CP candidate genes

##########################################################
## --  identify_TF_targeted_pull_candidat -- subset of CTS[['CP']] that are
## exclusively highly expressed (HiG) in CM (or CF)
## further narrow the subset to be open at CP or CM (or CF) ======================
library(BioNet)
packageVersion("BioNet") # '1.56.0'
library(igraph)
library(tibble)


mat <- read.table(paste0("heatmap_blocked_", CTS_name, "_NB_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv"), sep = "\t", header = T, check.names = FALSE)
print(dim(mat)) # [1] 139 166

# No ATAC data for NB: add dummy access columns (all 1) so identify_TF_targeted_pull_candidate
# accessibility filter becomes a no-op
for (col in c("PCW6CP_access", "PCW8_CM_access", "PCW19_CM_access",
              "PCW8_CF_access", "PCW19_CF_access",
              "PCW8_SMC_access", "PCW19_SMC_access", "iEPC_access")) {
    mat[[col]] <- 1L
}

for (key in key_TFs) {
    key_column <- which(grepl(key, colnames(mat)) & grepl("cisTarget_", colnames(mat)))
    if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key

    ## further narrow the subset to be open at CP or descendants
    graph_TF_list <- identify_TF_targeted_pull_candidate(mat, graph_list, CTS_name, CHD,
        key = key,
        keep_selfloop = TRUE, # whether to keep the self-loop of the key
        TF_bound_column_name = key_column,
        TF_appendix = key,
        edge_colored_by_Maven2023_ISL1KO = FALSE,
        key_in_TFfamily = key_in_TFfamily
    )
    saveRDS(graph_TF_list, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
}

print(names(graph_TF_list))
# [1] "CTSHiG.CP_TF.target"        "CTS.CP_TF.target_HiGCM"    "CTS.CP_TF.target_HiGCF"
# [4] "CTSHiG.CP_TF.target_CPopen"


for (key in key_TFs) {
    if (key == "HOX") key_in_TFfamily <- "HOXB2" else key_in_TFfamily <- key

    graph_TF_list <- readRDS(file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_v3.rds"))
    plot_TF_targeted_pull_candidate(graph_TF_list, key_in_TFfamily, CTS_name, saveFigure = TRUE)
}
# => PPI_graph_<key_TF>_GRN_prediction_CTS_8_v3.pdf


##########################################################
#### =quantify the edge weight changes and validate by Maven Isk1_KO in CP results ===================

# ## TO get teh cluster ID, remove the .1, .2, etc. from the cluster names which indicating additional CTS identified from the same cluster
# CM_ID = paste0('HiG_',CM_cluster) %>% gsub('\\.[1-9]','',.)
# CF_ID = paste0('HiG_',CF_cluster) %>% gsub('\\.[1-9]','',.)

# ## helper to play around with the layout
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


    res_cf <- fill_TF_targeting_predicted_edges(graph_TF_list,
        linkage_name = "CF", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
        HVG = rownames(sce)
    )
    cat(paste0(key, ' CF vcount(res_cf[["g_CT_sub"]]): ', vcount(res_cf[["g_CT_sub"]]), "\n"))
    if (vcount(res_cf[["g_CT_sub"]]) > 0) saveRDS(res_cf, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CF_final.rds"))
    names(res_cf)
    # [1] "g_CT_sub"              "g_descendant_sub"
}
# CHD1 CM vcount(res[["g_CT_sub"]]): 3
# MYCN CM vcount(res[["g_CT_sub"]]): 3
# EBF1 CM vcount(res[["g_CT_sub"]]): 4
# JUND CM vcount(res[["g_CT_sub"]]): 3
# CHD1 CF vcount(res_cf[["g_CT_sub"]]): 4
# MYCN CF vcount(res_cf[["g_CT_sub"]]): 6
# EBF1 CF vcount(res_cf[["g_CT_sub"]]): 7
# JUND CF vcount(res_cf[["g_CT_sub"]]): 6


### =reporting the results	and visualization #==========
print((files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds")))
# [1] "PPI_graph_CHD1_GRN_prediction_CTS_16_CF_final.rds"
# [2] "PPI_graph_CHD1_GRN_prediction_CTS_16_CM_final.rds"
# [3] "PPI_graph_EBF1_GRN_prediction_CTS_16_CF_final.rds"
# [4] "PPI_graph_EBF1_GRN_prediction_CTS_16_CM_final.rds"
# [5] "PPI_graph_JUND_GRN_prediction_CTS_16_CF_final.rds"
# [6] "PPI_graph_JUND_GRN_prediction_CTS_16_CM_final.rds"
# [7] "PPI_graph_MYCN_GRN_prediction_CTS_16_CF_final.rds"
# [8] "PPI_graph_MYCN_GRN_prediction_CTS_16_CM_final.rds"

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

    # tmp may have 0 or >1 rows if the motif matches multiple entries; take first row only.
    tmp1 <- if (nrow(tmp) > 0) tmp[1, , drop = FALSE] else data.frame(TF_highConf = NA, motif = NA, NES = NA)
    change_df <- cbind(linkage = pull, change_df,
                       TF_highConf = tmp1$TF_highConf, motif = tmp1$motif, NES = tmp1$NES)
    change_df$TF_highConf[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
    change_df$motif[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
    change_df$NES[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""

    final_table <- rbind(final_table, change_df)
}

print(dim(final_table)) # [1] 31 13  (8 files: 4 CM + 4 CF)
####### generating final table of the predicted subnetwork #################

write.table(final_table,
    file = paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

