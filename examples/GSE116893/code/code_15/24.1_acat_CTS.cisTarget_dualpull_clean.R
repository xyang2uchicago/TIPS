# Set True if running 24.0 the first time
rebuild_mat <- TRUE
source("C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/code/code_15/24.0_acat_load_input_clean.R")

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
CP_cluster # '15'
CM_cluster # '7'
CF_cluster # '9'

lengths(CTS) # a list
#  15  16   9 
#  76 139 110

class(sce) # [1] "SingleCellExperiment"

dim(mat) # [1] 76 4
colnames(mat)
# [1] "CP_hi"  "CM_hi"  "CF_hi"  "CTS_15"

seed_TF
CTS_ID # '15'
CTS_name # 'CTS_15'

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


    (dim(mat)) # [1] 76 87
    (colnames(mat))
    # ...
    # [70] "cisTarget_NFYA;NFYB;NFYC.motif_target"                                                                                                                                                                                                                     
    # [71] "cisTarget_YY2.motif_target"                                                                                                                                                                                                                                
    # [72] "cisTarget_E2F1;E2F3;E2F4;E2F6;E2F7;TFDP1;TFDP2;ZNF566;ZNF574.motif_target"                                                                                                                                                                                 
    # [73] "cisTarget_E2F1;E2F2;E2F3;E2F4;E2F5;E2F6;E2F7;RB1;TFDP1;TFDP2.motif_target"                                                                                                                                                                                 
    # [74] "cisTarget_E2F4;E2F7;ZFP69B.motif_target"                                                                                                                                                                                                                   
    # [75] "cisTarget_E2F6.motif_target"                                                                                                                                                                                                                               
    # [76] "cisTarget_E2F8;TFAP2C.motif_target"                                                                                                                                                                                                                        
    # [77] "cisTarget_E2F1;E2F4;KLF1;KLF10;KLF11;KLF12;KLF13;KLF14;KLF15;KLF17;KLF2;KLF5;KLF7;KLF8;MAZ;PATZ1;SP1;SP2;SP3;SP4;SP5;SP6;SP7;SP8;SP9;ZNF281.motif_target"                                                                                                  
    # [78] "cisTarget_WT1.motif_target"                                                                                                                                                                                                                                
    # [79] "cisTarget_RELA.motif_target"                                                                                                                                                                                                                               
    # [80] "cisTarget_MYBL2.motif_target"                                                                                                                                                                                                                              
    # [81] "cisTarget_SP1;SP2;SP3;SP4.motif_target"                                                                                                                                                                                                                    
    # [82] "cisTarget_SMAD3.motif_target"                                                                                                                                                                                                                              
    # [83] "cisTarget_ZNF548.motif_target"                                                                                                                                                                                                                             
    # [84] "cisTarget_ZFX.motif_target"                                                                                                                                                                                                                                
    # [85] "cisTarget_ACAA1;E2F2;E2F3;FOXN4;ZNF501.motif_target"                                                                                                                                                                                                       
    # [86] "cisTarget_SP1.motif_target"                                                                                                                                                                                                                                
    # [87] "cisTarget_CHD1.motif_target"   


    ## from the enriched motifs find the potential regulators that are either HiG of CTS member of a cluster of interest  ===================
    x <- colnames(mat)[grepl("cisTarget_", colnames(mat))]
    x <- gsub("cisTarget_", "", x) %>% gsub("\\.motif_target", "", .)
    x <- unlist(strsplit(x, ";")) %>% unique()
    print(x)
    #   [1] "CEBPZ"   "DRAP1"   "FOS"     "FOXI1"   "HMGXB3"  "IRF3"    "MBD2"   
    #   [8] "NFYA"    "NFYB"    "NFYC"    "PBX1"    "PBX3"    "SP2"     "YBX1"   
    #  [15] "E2F8"    "E2F1"    "E2F4"    "E2F6"    "TFDP1"   "E2F3"    "E2F7"   
    #  [22] "TFDP2"   "TAF1"    "YY1"     "YY2"     "SP4"     "ZNF566"  "ZNF574" 
    #  [29] "E2F2"    "E2F5"    "RB1"     "ZFP69B"  "ZNF142"  "ZNF141"  "KLF1"   
    #  [36] "ZBTB33"  "ZNF846"  "ZNF692"  "TFAP2B"  "TFAP2C"  "TFAP2E"  "EGR2"   
    #  [43] "EGR4"    "KLF10"   "KLF11"   "KLF12"   "KLF13"   "KLF14"   "KLF15"  
    #  [50] "KLF16"   "KLF17"   "KLF18"   "KLF2"    "KLF3"    "KLF4"    "KLF5"   
    #  [57] "KLF6"    "KLF7"    "KLF8"    "KLF9"    "MAZ"     "PATZ1"   "SP1"    
    #  [64] "SP3"     "SP5"     "SP6"     "SP7"     "SP8"     "SP9"     "VEZF1"  
    #  [71] "ZBTB14"  "ZBTB17"  "ZNF148"  "ZNF281"  "ZNF398"  "CHURC1"  "CTCF"   
    #  [78] "ETS1"    "MAFA"    "MZF1"    "NR3C1"   "OVOL2"   "SPZ1"    "SREBF2" 
    #  [85] "TFAP4"   "WT1"     "ZBTB7A"  "RELA"    "TRIM28"  "ZNF555"  "ZNF91"  
    #  [92] "ZNF730"  "MYBL1"   "MYBL2"   "EGR1"    "SMAD3"   "PML"     "ZNF548" 
    #  [99] "ZFY"     "ZNF571"  "MTF2"    "ZBTB7B"  "MYB"     "TFAP2A"  "SMAD9"  
    # [106] "ZNF449"  "ZNF530"  "ZNF534"  "ZNF662"  "ZNF383"  "ZNF682"  "NRL"    
    # [113] "TIGD1"   "ZNF445"  "ZNF425"  "ZFX"     "ACAA1"   "FOXN1"   "FOXN2"  
    # [120] "FOXN4"   "HES7"    "ONECUT2" "ZNF501"  "ZNF350"  "ZNF674"  "MAX"    
    # [127] "MYC"     "CHD1"    "ZGPAT"   "ZNF202" 

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
    # FOS is in DEG of 9
    # PBX1 is in DEG of 9
    # PBX3 is in DEG of 7
    # PBX3 is in DEG of 15
    # E2F3 is in DEG of 15
    # E2F7 is in DEG of 15
    # TFDP2 is in DEG of 15
    # E2F2 is in DEG of 15
    # RB1 is in DEG of 15
    # KLF12 is in DEG of 7
    # KLF12 is in DEG of 9
    # KLF12 is in DEG of 15
    # KLF13 is in DEG of 15
    # KLF7 is in DEG of 9
    # SP3 is in DEG of 15
    # ZNF148 is in DEG of 15
    # CTCF is in DEG of 15
    # ETS1 is in DEG of 15
    # NR3C1 is in DEG of 15
    # SREBF2 is in DEG of 15
    # ZNF730 is in DEG of 15
    # MTF2 is in DEG of 15
    # MYB is in DEG of 15
    # CHD1 is in DEG of 15


    for (i in x) {
        if (i %in% CTS[[CTS_ID]]) {
            cat(paste0(i, " is in CTS of ", CTS_ID, "\n"))
            key_TFs <- c(key_TFs, i)
        }
    }

    key_TFs <- unique(key_TFs)
    (key_TFs)
    # [1] "FOS"    "PBX1"   "PBX3"   "E2F3"   "E2F7"   "TFDP2"  "E2F2"   "RB1"   
    # [9] "KLF12"  "KLF13"  "KLF7"   "SP3"    "ZNF148" "CTCF"   "ETS1"   "NR3C1" 
    # [17] "SREBF2" "ZNF730" "MTF2"   "MYB"    "CHD1"  


    mat <- as.data.frame(mat)

    ##  filter out the key_TFs that come from one motif thus will share target genes for the downstream predictions	=======================
    ## ====== note that this step is manual, such as selecting the TFs of the same motif that are most close to the seed TF ======
    # (x = which(grepl(seed_TF[3], colnames(mat)) | grepl(seed_TF[2], colnames(mat)) | grepl(seed_TF[3], colnames(mat)) ) )  # TF candidateds for CTS.CP
    (x <- intersect(which(grepl("cisTarget_", colnames(mat))), which(Reduce("|", lapply(key_TFs, function(p) grepl(p, colnames(mat), fixed = F))))))
    (x)
    #  [1]  5  8  9 10 13 14 15 24 27 28 30 34 35 36 37 42 44 60 61 62 65 66 68 69 72
    # [26] 73 74 77 80 81 85 87
    (colnames(mat)[x])
    #  [1] "cisTarget_CEBPZ;DRAP1;FOS;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX1;PBX3;SP2;YBX1.motif_target"                                                                                                                                                           
    #  [2] "cisTarget_E2F1;E2F3;E2F4;E2F6;E2F7;E2F8;TFDP1;TFDP2.motif_target"                                                                                                                                                                                          
    #  [3] "cisTarget_CEBPZ;DRAP1;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX1;PBX3;SP2.motif_target"                                                                                                                                                                    
    #  [4] "cisTarget_CEBPZ;DRAP1;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX3;SP2.motif_target"                                                                                                                                                                         
    #  [5] "cisTarget_E2F1;E2F3;E2F4;E2F6;E2F7;E2F8;SP4;TFDP1;TFDP2;ZNF566;ZNF574.motif_target"                                                                                                                                                                        
    #  [6] "cisTarget_E2F1;E2F2;E2F3;E2F4;E2F5;E2F6;E2F7;E2F8;RB1;TFDP1;TFDP2.motif_target"                                                                                                                                                                            
    #  [7] "cisTarget_E2F1;E2F3;E2F4;E2F6;E2F7;TFDP1;TFDP2;ZFP69B.motif_target"                                                                                                                                                                                        
    #  [8] "cisTarget_E2F1;E2F2;E2F3;E2F4;E2F5;E2F6;E2F7.motif_target"                                                                                                                                                                                                 
    #  [9] "cisTarget_E2F1;E2F4;E2F6;E2F7;EGR2;EGR4;KLF1;KLF10;KLF11;KLF12;KLF13;KLF14;KLF15;KLF16;KLF17;KLF18;KLF2;KLF3;KLF4;KLF5;KLF6;KLF7;KLF8;KLF9;MAZ;PATZ1;SP1;SP2;SP3;SP4;SP5;SP6;SP7;SP8;SP9;TFDP1;TFDP2;VEZF1;ZBTB14;ZBTB17;ZNF148;ZNF281;ZNF398.motif_target"
    # [10] "cisTarget_CHURC1;CTCF;ETS1;MAFA;MZF1;NR3C1;OVOL2;SPZ1;SREBF2;TFAP4;WT1;YY1;ZBTB7A.motif_target"                                                                                                                                                            
    # [11] "cisTarget_E2F3.motif_target"                                                                                                                                                                                                                               
    # [12] "cisTarget_E2F1;E2F3;E2F4;TFDP1.motif_target"                                                                                                                                                                                                               
    # [13] "cisTarget_ZNF730.motif_target"                                                                                                                                                                                                                             
    # [14] "cisTarget_MYBL1;MYBL2.motif_target"                                                                                                                                                                                                                        
    # [15] "cisTarget_EGR1;KLF1;KLF10;KLF12;KLF13;KLF14;KLF15;KLF18;KLF2;KLF3;KLF5;KLF6;KLF7;KLF9;MAZ;SP1;SP2;SP3;SP4;SP5;SP9;WT1.motif_target"                                                                                                                        
    # [16] "cisTarget_MTF2.motif_target"                                                                                                                                                                                                                               
    # [17] "cisTarget_MYB.motif_target"                                                                                                                                                                                                                                
    # [18] "cisTarget_ACAA1;E2F2;E2F3;FOXN1;FOXN2;FOXN4;HES7;ONECUT2;ZBTB14;ZNF501.motif_target"                                                                                                                                                                       
    # [19] "cisTarget_SP1;SP3.motif_target"                                                                                                                                                                                                                            
    # [20] "cisTarget_EGR2;KLF12;KLF7;KLF8;SP4;ZNF682.motif_target"                                                                                                                                                                                                    
    # [21] "cisTarget_CHD1;ZGPAT;ZNF202.motif_target"                                                                                                                                                                                                                  
    # [22] "cisTarget_CEBPZ;DRAP1;FOS;FOXI1;HMGXB3;IRF3;NFYA;NFYB;NFYC;PBX1;PBX3;SP2.motif_target"                                                                                                                                                                     
    # [23] "cisTarget_E2F1;E2F3;E2F6;E2F8;TFDP1.motif_target"                                                                                                                                                                                                          
    # [24] "cisTarget_MBD2;NFYA;NFYB;NFYC;PBX3.motif_target"                                                                                                                                                                                                           
    # [25] "cisTarget_E2F1;E2F3;E2F4;E2F6;E2F7;TFDP1;TFDP2;ZNF566;ZNF574.motif_target"                                                                                                                                                                                 
    # [26] "cisTarget_E2F1;E2F2;E2F3;E2F4;E2F5;E2F6;E2F7;RB1;TFDP1;TFDP2.motif_target"                                                                                                                                                                                 
    # [27] "cisTarget_E2F4;E2F7;ZFP69B.motif_target"                                                                                                                                                                                                                   
    # [28] "cisTarget_E2F1;E2F4;KLF1;KLF10;KLF11;KLF12;KLF13;KLF14;KLF15;KLF17;KLF2;KLF5;KLF7;KLF8;MAZ;PATZ1;SP1;SP2;SP3;SP4;SP5;SP6;SP7;SP8;SP9;ZNF281.motif_target"                                                                                                  
    # [29] "cisTarget_MYBL2.motif_target"                                                                                                                                                                                                                              
    # [30] "cisTarget_SP1;SP2;SP3;SP4.motif_target"                                                                                                                                                                                                                    
    # [31] "cisTarget_ACAA1;E2F2;E2F3;FOXN4;ZNF501.motif_target"                                                                                                                                                                                                       
    # [32] "cisTarget_CHD1.motif_target"     

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

    key_TFs <- c("E2F3", "MYB", "ZNF730", "MTF2")
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
# key_TFs =  c('TCF3','HMGA2','RARB')
# fileName=paste0('heatmap_blocked_',CTS_name,'_NB_cisTarget_',paste(key_TFs, collapse='_'),'_v3.tsv')
# mat = read.table(fileName, sep='\t',header=T, check.names=F)
# x = c(  34 ,40, 48)
# names(x) = key_TFs


print(dim(mat)) # [1] 76 105
print(colnames(mat))
#   [1] "CP_hi"                                                                                                                                                                                                                                                     
#   [2] "CM_hi"                                                                                                                                                                                                                                                     
#   [3] "CF_hi"                                                                                                                                                                                                                                                     
#   [4] "CTS_15"                                                                                                                                                                                                                                                    
#   [5] "cisTarget_CEBPZ;DRAP1;FOS;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX1;PBX3;SP2;YBX1.motif_target"                                                                                                                                                           
#   [6] "cisTarget_E2F8.motif_target"                                                                                                                                                                                                                               
#   [7] "cisTarget_E2F1;E2F4;E2F6;TFDP1.motif_target"                                                                                                                                                                                                               
#   [8] "cisTarget_E2F1;E2F3;E2F4;E2F6;E2F7;E2F8;TFDP1;TFDP2.motif_target"                                                                                                                                                                                          
#   [9] "cisTarget_CEBPZ;DRAP1;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX1;PBX3;SP2.motif_target"                                                                                                                                                                    
#  [10] "cisTarget_CEBPZ;DRAP1;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX3;SP2.motif_target"                                                                                                                                                                         
#  [11] "cisTarget_TAF1;YY1;YY2.motif_target"                                                                                                                                                                                                                       
#  [12] "cisTarget_E2F4.motif_target"                                                                                                                                                                                                                               
#  [13] "cisTarget_E2F1;E2F3;E2F4;E2F6;E2F7;E2F8;SP4;TFDP1;TFDP2;ZNF566;ZNF574.motif_target"                                                                                                                                                                        
#  [14] "cisTarget_E2F1;E2F2;E2F3;E2F4;E2F5;E2F6;E2F7;E2F8;RB1;TFDP1;TFDP2.motif_target"                                                                                                                                                                            
#  [15] "cisTarget_E2F1;E2F3;E2F4;E2F6;E2F7;TFDP1;TFDP2;ZFP69B.motif_target"                                                                                                                                                                                        
#  [16] "cisTarget_SP2.motif_target"                                                                                                                                                                                                                                
#  [17] "cisTarget_ZNF142.motif_target"                                                                                                                                                                                                                             
#  [18] "cisTarget_ZNF141.motif_target"                                                                                                                                                                                                                             
#  [19] "cisTarget_E2F4;E2F6.motif_target"                                                                                                                                                                                                                          
#  [20] "cisTarget_KLF1.motif_target"                                                                                                                                                                                                                               
#  [21] "cisTarget_ZBTB33.motif_target"                                                                                                                                                                                                                             
#  [22] "cisTarget_E2F1.motif_target"                                                                                                                                                                                                                               
#  [23] "cisTarget_ZNF846.motif_target"                                                                                                                                                                                                                             
#  [24] "cisTarget_E2F1;E2F2;E2F3;E2F4;E2F5;E2F6;E2F7.motif_target"                                                                                                                                                                                                 
#  [25] "cisTarget_ZFP69B;ZNF692.motif_target"                                                                                                                                                                                                                      
#  [26] "cisTarget_E2F8;TFAP2B;TFAP2C;TFAP2E.motif_target"
#  ...
#  [87] "cisTarget_CHD1.motif_target"                                                                                                                                                                                                                               
#  [88] "E2F3_CP_candidate"                                                                                                                                                                                                                                         
#  [89] "E2F3_CM_candidate"                                                                                                                                                                                                                                         
#  [90] "E2F3_CF_candidate"                                                                                                                                                                                                                                         
#  [91] "MYB_CP_candidate"                                                                                                                                                                                                                                          
#  [92] "MYB_CM_candidate"                                                                                                                                                                                                                                          
#  [93] "MYB_CF_candidate"                                                                                                                                                                                                                                          
#  [94] "PBX3_CP_candidate"                                                                                                                                                                                                                                         
#  [95] "PBX3_CM_candidate"                                                                                                                                                                                                                                         
#  [96] "PBX3_CF_candidate"                                                                                                                                                                                                                                         
#  [97] "ZNF730_CP_candidate"                                                                                                                                                                                                                                       
#  [98] "ZNF730_CM_candidate"                                                                                                                                                                                                                                       
#  [99] "ZNF730_CF_candidate"                                                                                                                                                                                                                                       
# [100] "MTF2_CP_candidate"                                                                                                                                                                                                                                         
# [101] "MTF2_CM_candidate"                                                                                                                                                                                                                                         
# [102] "MTF2_CF_candidate"                                                                                                                                                                                                                                         
# [103] "NA_CP_candidate"                                                                                                                                                                                                                                           
# [104] "NA_CM_candidate"                                                                                                                                                                                                                                           
# [105] "NA_CF_candidate"   

motif_TF_highConf <- gsub("cisTarget_", "", colnames(mat)[x]) %>% gsub("\\.motif_target", "", .)
print(motif_TF_highConf)
#  [1] "CEBPZ;DRAP1;FOS;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX1;PBX3;SP2;YBX1"                                                                                                                                                           
#  [2] "E2F1;E2F3;E2F4;E2F6;E2F7;E2F8;TFDP1;TFDP2"                                                                                                                                                                                          
#  [3] "CEBPZ;DRAP1;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX1;PBX3;SP2"                                                                                                                                                                    
#  [4] "CEBPZ;DRAP1;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX3;SP2"                                                                                                                                                                         
#  [5] "E2F1;E2F3;E2F4;E2F6;E2F7;E2F8;SP4;TFDP1;TFDP2;ZNF566;ZNF574"                                                                                                                                                                        
#  [6] "E2F1;E2F2;E2F3;E2F4;E2F5;E2F6;E2F7;E2F8;RB1;TFDP1;TFDP2"                                                                                                                                                                            
#  [7] "E2F1;E2F3;E2F4;E2F6;E2F7;TFDP1;TFDP2;ZFP69B"                                                                                                                                                                                        
#  [8] "E2F1;E2F2;E2F3;E2F4;E2F5;E2F6;E2F7"                                                                                                                                                                                                 
#  [9] "E2F1;E2F4;E2F6;E2F7;EGR2;EGR4;KLF1;KLF10;KLF11;KLF12;KLF13;KLF14;KLF15;KLF16;KLF17;KLF18;KLF2;KLF3;KLF4;KLF5;KLF6;KLF7;KLF8;KLF9;MAZ;PATZ1;SP1;SP2;SP3;SP4;SP5;SP6;SP7;SP8;SP9;TFDP1;TFDP2;VEZF1;ZBTB14;ZBTB17;ZNF148;ZNF281;ZNF398"
# [10] "CHURC1;CTCF;ETS1;MAFA;MZF1;NR3C1;OVOL2;SPZ1;SREBF2;TFAP4;WT1;YY1;ZBTB7A"                                                                                                                                                            
# [11] "E2F3"                                                                                                                                                                                                                               
# [12] "E2F1;E2F3;E2F4;TFDP1"                                                                                                                                                                                                               
# [13] "ZNF730"                                                                                                                                                                                                                             
# [14] "MYBL1;MYBL2"                                                                                                                                                                                                                        
# [15] "EGR1;KLF1;KLF10;KLF12;KLF13;KLF14;KLF15;KLF18;KLF2;KLF3;KLF5;KLF6;KLF7;KLF9;MAZ;SP1;SP2;SP3;SP4;SP5;SP9;WT1"                                                                                                                        
# [16] "MTF2"                                                                                                                                                                                                                               
# [17] "MYB"                                                                                                                                                                                                                                
# [18] "ACAA1;E2F2;E2F3;FOXN1;FOXN2;FOXN4;HES7;ONECUT2;ZBTB14;ZNF501"                                                                                                                                                                       
# [19] "SP1;SP3"                                                                                                                                                                                                                            
# [20] "EGR2;KLF12;KLF7;KLF8;SP4;ZNF682"                                                                                                                                                                                                    
# [21] "CHD1;ZGPAT;ZNF202"                                                                                                                                                                                                                  
# [22] "CEBPZ;DRAP1;FOS;FOXI1;HMGXB3;IRF3;NFYA;NFYB;NFYC;PBX1;PBX3;SP2"                                                                                                                                                                     
# [23] "E2F1;E2F3;E2F6;E2F8;TFDP1"                                                                                                                                                                                                          
# [24] "MBD2;NFYA;NFYB;NFYC;PBX3"                                                                                                                                                                                                           
# [25] "E2F1;E2F3;E2F4;E2F6;E2F7;TFDP1;TFDP2;ZNF566;ZNF574"                                                                                                                                                                                 
# [26] "E2F1;E2F2;E2F3;E2F4;E2F5;E2F6;E2F7;RB1;TFDP1;TFDP2"                                                                                                                                                                                 
# [27] "E2F4;E2F7;ZFP69B"                                                                                                                                                                                                                   
# [28] "E2F1;E2F4;KLF1;KLF10;KLF11;KLF12;KLF13;KLF14;KLF15;KLF17;KLF2;KLF5;KLF7;KLF8;MAZ;PATZ1;SP1;SP2;SP3;SP4;SP5;SP6;SP7;SP8;SP9;ZNF281"                                                                                                  
# [29] "MYBL2"                                                                                                                                                                                                                              
# [30] "SP1;SP2;SP3;SP4"                                                                                                                                                                                                                    
# [31] "ACAA1;E2F2;E2F3;FOXN4;ZNF501"                                                                                                                                                                                                       
# [32] "CHD1"    

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

# candidate genes:  0 
# simplified motif name:  cisTarget_CEBPZ;DRAP1;FOS;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX1;PBX3;SP2;YBX1.motif_target cisTarget_CEBPZ;DRAP1;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX1;PBX3;SP2.motif_target cisTarget_CEBPZ;DRAP1;FOXI1;HMGXB3;IRF3;MBD2;NFYA;NFYB;NFYC;PBX3;SP2.motif_target cisTarget_CEBPZ;DRAP1;FOS;FOXI1;HMGXB3;IRF3;NFYA;NFYB;NFYC;PBX1;PBX3;SP2.motif_target cisTarget_MBD2;NFYA;NFYB;NFYC;PBX3.motif_target 

##########################################################
## --  identify_TF_targeted_pull_candidat -- subset of CTS[['CP']] that are
## exclusively highly expressed (HiG) in CM (or CF)
## further narrow the subset to be open at CP or CM (or CF) ======================
library(BioNet)
packageVersion("BioNet") # '1.56.0'
library(igraph)
library(tibble)


mat <- read.table(paste0("heatmap_blocked_", CTS_name, "_NB_cisTarget_", paste(key_TFs, collapse = "_"), "_v3.tsv"), sep = "\t", header = T, check.names = FALSE)
print(dim(mat)) # [1] 76 99

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
    cat(paste0(key, ' CM vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), "\n"))
    if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CM_final.rds"))


    res <- fill_TF_targeting_predicted_edges(graph_TF_list,
        linkage_name = "CF", graph_list,
        sce, celltype_col = celltype_col, CT_cluster_id = CP_cluster,
        descendant_cluster_id = CF_cluster, TF_symbol = key_in_TFfamily,
        HVG = rownames(sce)
    )
    cat(paste0(key, ' CF vcount(res[["g_CT_sub"]]): ', vcount(res[["g_CT_sub"]]), "\n"))
    if (vcount(res[["g_CT_sub"]]) > 0) saveRDS(res, file = paste0("PPI_graph_", key_in_TFfamily, "_GRN_prediction_", CTS_name, "_CF_final.rds"))
    names(res)
    # [1] "g_CT_sub"              "g_descendant_sub"
}
# E2F3  CM vcount: 4   CF vcount: 0
# MYB   CM vcount: 6   CF vcount: 0
# ZNF730 CM vcount: 3  CF vcount: 0
# MTF2  CM vcount: 3   CF vcount: 0
# CF direction absent: CF_hi=0 in mat (no CTS_15 gene overlaps DEG[["9"]]) — biological result, not fixable by SCE/column choice.


### =reporting the results	and visualization #==========
print((files <- list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds")))
# [1] "PPI_graph_E2F3_GRN_prediction_CTS_15_CM_final.rds"
# [2] "PPI_graph_MTF2_GRN_prediction_CTS_15_CM_final.rds"
# [3] "PPI_graph_MYB_GRN_prediction_CTS_15_CM_final.rds"
# [4] "PPI_graph_ZNF730_GRN_prediction_CTS_15_CM_final.rds"


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

    # tmp may have 0 or >1 rows; take first row only to avoid cbind row-count mismatch.
    tmp1 <- if (nrow(tmp) > 0) tmp[1, , drop = FALSE] else data.frame(TF_highConf = NA, motif = NA, NES = NA)
    change_df <- cbind(linkage = pull, change_df, TF_highConf = tmp1$TF_highConf, motif = tmp1$motif, NES = tmp1$NES)
    change_df$TF_highConf[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
    change_df$motif[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""
    change_df$NES[which(change_df$from != key_in_TFfamily & change_df$to != key_in_TFfamily)] <- ""

    final_table <- rbind(final_table, change_df)
}

print(dim(final_table)) # [1] 18 13  (4 CM-only files; CF direction absent)
####### generateing final table of the predicted subnetwork #################

write.table(final_table,
    file = paste0("PPI_graph_GRN_prediction_", CTS_name, "_dualpull_final_table.tsv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)

