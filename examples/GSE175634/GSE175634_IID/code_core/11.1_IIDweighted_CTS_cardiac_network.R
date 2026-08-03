

## compare the cardiac CTS from hESC:C10, E2018 2018 cardiac.a, and E8.25 2019 C8 (not the in vivo and EB C4)
## compare HEP CTS from E8.25 2019 C13, E8.25 2018 endo.b, EB C11,
# and E8.25 2019 C17 (or C8)
library(gplots)
require(dplyr)
library(data.table)
library(ggplot2)
library("gridExtra")
library(ggrepel)
library(ggpubr) # resuired by stat_compare_means()

setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_IID')
score_threshold = 'weight'
PPI_color_platte = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

HVG <- list()

load('F:/projects/scRNA/results/GSE87038_gastrulation/uncorrected/E8.25_mesoderm/sce_E8.25_uncorrected.RData')
HVG[[1]] <- rownames(sce)
load('E:/Git_Holly/scRNAseq_examples/data/gastrulationE8.25_Ibarra-Soria2018/sce_16subtype.RData')
HVG[[2]] <- rownames(sce)
load('E:/Git_Holly/scRNAseq_examples/data/EB_Zhao2019/sce.GSE130146_noenderdormPgerm.RData')
HVG[[3]] <- rownames(sce)
load('E:/Git_Holly/scRNAseq_examples/data/hESC_Bargaje2017/sce_hESC.RData')
HVG[[4]] <- rownames(sce)

tmp = read.csv('F:/projects/scRNA/results/GSE175634_iPSC_CM/2024_3kHVG/3khvg.tsv', header=F)
HVG[[5]] <- tmp[,1]

names(HVG) <- c('Pijuan-Sala2019','Ibarra-Soria2018', 'EB', 'hESC2017', 'hiPSC_CM2022')
lengths(HVG)
# Pijuan-Sala2019 Ibarra-Soria2018               EB         hESC2017     hiPSC_CM2022 
#           10938             4000             4000               96             3000 
rm(sce)
################################################################
############testing the GSE175634 dataset ###############
## TF-target using PROGENy, becasue 'PROGENy with 500 or 1000 footprint genes per pathway yields the best performance' 
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1949-z
# PROGENy is a tool that infers pathway activity for 14 signaling pathways (Androgen, Estrogen, EGFR, Hypoxia, JAK-STAT, MAPK, NFkB, PI3K, p53, TGFb, TNFa, Trail, VEGF, and WNT) from gene expression data
# DoRothEA is a gene set resource containing signed transcription factor (TF)-target interactions [13]. Those interactions were curated and collected from different types of evidence such as literature curated resources, ChIP-seq peaks, TF binding site motifs, and interactions inferred directly from gene expression
## 1) gent the CTS genes for transitory PPI_cats   
################################################################
## manually record the Ic predictions

subfolds <- c('hESC_Bargaje2017', 'lung_Treutlein2014', 
			  'gastrulationE8.25_Pijuan-Sala2019','gastrulationE8.25_Ibarra-Soria2018',
			  'EB_Zhao2019','simulated_EMT',
			  '2024_3kHVG/')
subsubfolds <- c('C_Consensus_929cells', 'C_Leiden_0.4', 
					 'C_SNNGraph_allcells', 'subcelltype',
					 'C_SNNGraph_k5', 'C_cell_type',
					 'BioTIP')


db=7
	inputDir = ifelse(db<7,
		"E:/Git_Holly/scRNAseq_examples/result", 
		"F:/projects/scRNA/results/GSE175634_iPSC_CM/")

  load(file=paste0(inputDir, subfolds[db] ,'/',subsubfolds[db ],'/CTS_Lib_Scaledata.RData'))
  lengths(CTS.Lib.Symbol)
	# muscle   endoderm      CP/CF         CM         CF         CP    mix pro
		# 63         66         77        113        109         76         69
# endoderm.1       CP.1
	   # 136         87
	if(db==7) {
		sig.DNB = c(T, T, T, T, T, T, F, F, T)
		sig.IC = c(T, T, F, F, T, T, NA, NA, T)
		sig.deltaIC = c(T, T, NA, NA, F, T, NA, NA, T)
		CTS = CTS.Lib.Symbol[which(sig.DNB==T & sig.IC==T & sig.deltaIC==T)]
		}
	lengths(CTS)
  # muscle endoderm       CP     CP.1 
      # 63       66       76       87 


########################
## load DEGs refer to 3.5_Filter_Plot_Marker_diffBar.R
new_object <- new.env()                # Create a new environment
load(paste0(inputDir, "2024_3kHVG/DEG.wc_padj0.01_score40.RData"), envir = new_object)  # Load .RData into this environment
markers.up <- new_object[[ls(new_object)[1]]]
names(markers.up)[which(names(markers.up)=='3')] = 'CP'
names(markers.up)[which(names(markers.up)=='8')] = 'muscle'
names(markers.up)[which(names(markers.up)=='11')] = 'endoderm'
class(markers.up[[1]])  #[1] "data.frame"
head(markers.up[[1]], 4)
      # names   scores logfoldchanges pvals pvals_adj
# 1     L1TD1 261.4430            NaN     0         0
# 2     TERF1 247.4386            NaN     0         0
# 3 MIR302CHG 228.6239            NaN     0         0
# 4      PTMA 222.5636            NaN     0         0
range(markers.up[[1]]$scores)
#[1]  40.06131 261.44299	This is a well filtered DEG set
	
 load(paste0(inputDir,'2024_3kHVG/DEG.wc_padj0.01_score40.RData'))
 lapply(DEG.wc, nrow) %>% unlist
   # 0    1    2    3    4    5    6    7    8    9   10   11   12 
 # 344 1224  257  152  110  350  110  117  119 1633   81   44   25 
 DEG = list()
 for(i in 1:length(DEG.wc)) DEG[[i]] = DEG.wc[[i]]$names %>% unique
 lengths(DEG)
#  [1]  344 1224  257  152  110  350  110  117  119 1633   81   44   25
names(DEG) = 0:(length(DEG)-1)
names(DEG)[which(names(DEG)=='3')] = 'CP'
names(DEG)[which(names(DEG)=='8')] = 'muscle'
names(DEG)[which(names(DEG)=='11')] = 'endoderm'

######################################################
# 3) load IID db and build GRN 
 ## sdownloaded on 2026 March by Horatio Ai
 # https://iid.ophid.utoronto.ca/
 ################################################################
  
library(igraph)
library(data.table)

# refer to 0_IIDweighted_PPIN_build.R
g_iid_global <- readRDS(file='iid2025_human_mouse_conserved_global.rds')
vcount(g_iid_global) # 17182
ecount(g_iid_global) # 1517715


get_iid_subnetwork <- function(g_global, hits) {
    hits <- toupper(hits)
    hits <- unique(hits[!is.na(hits) & hits != ""])
    hits <- intersect(hits, V(g_global)$name)
    if (length(hits) < 2) return(make_empty_graph())
    induced_subgraph(g_global, vids = hits)
}



graph_list <- list()
# build for traditional up-regulated markers using IID db
for (i in names(DEG)) {
    # Filter differentially expressed genes based on logFC and FDR thresholds
    diff_exp <- markers.up[[i]]
    diff_exp$symbol <- diff_exp$names
   # diff_exp <- subset(diff_exp, summary.logFC > logFC.cut & FDR < 0.01) # we used summary.logFC < 0.6 in mouse datasets which used t-test historically to define HiGs
    diff_exp <- subset(diff_exp, scores> 40 & pvals_adj < 0.01) ## which was already filtered before saving the DEG object, jsut ensure about the cutoffs
    hits <- diff_exp$symbol

    # Get the subnetwork
    graph <- get_iid_subnetwork(g_iid_global, hits)

  # Add node attributes (only for nodes that exist in the graph)
  if (vcount(graph) > 0) {
    # The 'scores' column reflects the strength of the differential expression for each gene based on the Wilcoxon rank-sum test. A high score suggests that the gene's expression is significantly different between the groups under comparison.
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), "scores"]  # $scores   #summary.logFC if scanpy using t.test
    V(graph)$FDR    <- diff_exp[match(V(graph)$name, diff_exp$symbol), "pvals_adj"]  ## FDR of scanpy_wilcox marker identification
  }

  if (vcount(graph) > 0 && all(toupper(V(graph)$name) %in% toupper(DEG[[i]]))) {
    # map vertex names to the exact form used in DEG[[i]]
    V(graph)$name <- DEG[[i]][match(toupper(V(graph)$name), toupper(DEG[[i]]))]
  }

  graph_list[[paste0("HiG_", i)]] <- graph
}

 
# instead build for transitory PPI_cats only with the significantly up-regulated CTS genes !!!!!
# build for (up-regulated_marker intersecting CTS)
# HiGCTS
for (i in names(CTS)) {
    # Get unique cluster ids for clusters containing subclusters labeled by numerical id
    j <- if (grepl("\\.", i) && grepl("^[0-9]+$", sub("^[^.]*\\.", "", i))) sub("\\..*$", "", i) else i

    # Get the full marker table for that cluster
    deg_table <- markers.up[[j]]
 
    # Intersect with CTS
    deg_table <- deg_table[deg_table$names %in% CTS[[i]], ]

    # Filter on significance
    diff_exp <- subset(deg_table, scores> 40 & pvals_adj < 0.01)
    graph <- get_iid_subnetwork(g_iid_global, diff_exp$names)

    if (vcount(graph) > 0) {
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), "scores"]
    V(graph)$FDR    <- diff_exp[match(V(graph)$name, diff_exp$symbol), "pvals_adj"]
    }

    graph_list[[paste0("HiGCTS_", i)]] <- graph
}

## lastly, build for CTS
 # refer to 6.3_DE.statistics_CTS.R 
 markers.up_all = readRDS(paste0(inputDir,subfolds[db] ,'/DEG.wc_nofiltering.rds'))
 names(markers.up_all)[which(names(markers.up_all)=='3')] = 'CP'
 names(markers.up_all)[which(names(markers.up_all)=='8')] = 'muscle'
 names(markers.up_all)[which(names(markers.up_all)=='11')] = 'endoderm'
 for (i in names(CTS)) {
    j <- if (grepl("\\.", i) && grepl("^[0-9]+$", sub("^[^.]*\\.", "", i))) sub("\\..*$", "", i) else i
    rownames(markers.up_all[[i]]) = markers.up_all[[i]]$names
    diff_exp <- markers.up_all[[j]][CTS[[i]], ]
 
    hits <- diff_exp$names
    graph <- get_iid_subnetwork(g_iid_global, hits)

    if (vcount(graph) > 0) {
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), "scores"]
    V(graph)$FDR    <- diff_exp[match(V(graph)$name, diff_exp$symbol), "pvals_adj"]
    }

    graph_list[[paste0("CTS_", i)]] <- graph
}


lengths(graph_list)
    #        HiG_0           HiG_1           HiG_2          HiG_CP 
    #         307            1090             213             117 
    #       HiG_4           HiG_5           HiG_6           HiG_7 
    #          87             326              92              90 
    #  HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
    #         114            1352              60              39 
    #      HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP 
    #          20              19               7              27 
    # HiGCTS_CP.1      CTS_muscle    CTS_endoderm          CTS_CP 
    #          26              58              61              67 
    #    CTS_CP.1 
    #          70 
 
saveRDS(graph_list, file= 'GSE175634_IID_graph_perState_notsimplified.rds')  #!!!!!!!!!!!!!!!!!!!
 
 graph_list <- readRDS( file= 'GSE175634_IID_graph_perState_notsimplified.rds')  
 graph_list <- lapply(graph_list, simplify, edge.attr.comb ='max') #!!!!!!!!!!!!!!!!!!!
 lengths(graph_list)
    #       HiG_0           HiG_1           HiG_2          HiG_CP 
    #         307            1090             213             117 
    #       HiG_4           HiG_5           HiG_6           HiG_7 
    #          87             326              92              90 
    #  HiG_muscle           HiG_9          HiG_10    HiG_endoderm 
    #         114            1352              60              39 
    #      HiG_12   HiGCTS_muscle HiGCTS_endoderm       HiGCTS_CP 
    #          20              19               7              27 
    # HiGCTS_CP.1      CTS_muscle    CTS_endoderm          CTS_CP 
    #          26              58              61              67 
    #    CTS_CP.1 
    #          70 

## degree distribution was not discussed again because it is not in the manuscript  

# Check which graphs have duplicate vertex names
graphs_with_duplicates <- sapply(graph_list, function(g) {
  vertex_names <- V(g)$name
  if(is.null(vertex_names)) {
    # If no names, use vertex indices
    vertex_names <- V(g)
  }
  any(duplicated(vertex_names))
})

# See which graphs have duplicates
which(graphs_with_duplicates)
named integer(0)
