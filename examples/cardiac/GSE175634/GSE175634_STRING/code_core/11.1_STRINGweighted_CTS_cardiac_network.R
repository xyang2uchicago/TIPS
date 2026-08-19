# this version added E(graph)$weight;  formated the V(graph)$weight to be previous V(graph)$score (for iPSC) or V(graph)$avg_logFC
# Notice in igraph, E(g)$weight is the distance, ie. the higher weight the far two nodes;  while in PPI, the high score, the short distance
# therfore, we process the E(g)$weight to be 1/score or corefficience
#
# this version simplify duplicated edges between any two vertices 
# Q1: Does the GRN of 37 CTS genes in E8.25_2018 dataset fragile?
# # to build PPI_cat-specific GRN for each PPI_cat 
# identify PPI_cat-specially highly expressed genes for steady PPI_cat but using the CTS for transitory PPI_cats

# Q2: Does the GRN of 19 CTS genes in hESC dataset fragile? 
# -- using string physicial links
# -- using HumanBase 'cardiac muscle' links
# to build PPI_cat-specific GRN for C5, C9 (CT), C7, C8, and C10 (CT) 
# 1-3) identify PPI_cat-specially highly expressed genes for C5, C7, and C8, and using the identified CTS for C9, C10
# 4) evaluate the GRN robustness (stability) of each PPI_cat
# 4.1) estimate the significance of degree distribution
# 4.2) estimate the significane of robustness
# 5) Does the most individual contributors matter 
#  ('specialized processes (e.g., cell differentiation) are mainly related to regulators with low Knn,
#    whereas essential processes are mainly related to regulators with high page rank or degree')
# 6) Does the GRN of  transitory PPI_cat exhibit bisignatures but not the steady PPI_cat?
# 7) Will the GRN of tranistory PPI_cat (smalle number of nodes) be compareable to the Subset the GRN of steady PPI_cat with the same number of nodes?


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

setwd('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weight')
score_threshold = 'weight'
PPI_color_platte = c("CTS" = "#7570B3", "HiGCTS" = "#E7298A", "HiG" = "#E6AB02")

HVG <- list()
## refer to F:\projects\scRNA\source\GSE87038_gastrulation\SC3_GSE87038_E8.25_HEP.R, 
# Note the published is a subset of the original 7k cell object 
# load('../data/gastrulationE8.25_Pijuan-Sala2019/sce_E8.25_HEP.RData')
# sce
# # dim: 3073 1362
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

  load(file=here::here("examples", "cardiac", "GSE175634", "data", "CTS_Lib_Scaledata.RData"))
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
load(here::here("examples", "cardiac", "GSE175634", "data", "DEG.wc_padj0.01_score40.RData"), envir = new_object)  # Load .RData into this environment
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
	
	
 load(here::here("examples", "cardiac", "GSE175634", "data", "DEG.wc_padj0.01_score40.RData"))
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
# 3) load STRING db and build GRN 
 ## second trial, using STRING, IFT57/74/88 are all included
 # https://www.bioconductor.org/packages/release/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf
 # refer to CTS_cardiac_ntwork_robustness_notsimplified.R
 ################################################################
 library(BioNet)
 packageVersion('BioNet') # '1.56.0'
 library(igraph)
 library("STRINGdb")
 packageVersion('STRINGdb') # '2.14.0'
 library(tibble)
 
 string_db <- STRINGdb$new( version="12.0", species=9606,   # species= 10090 for mouse;    ?? for human 2021 version
                            score_threshold = 200, #!!!!!!!!!!!default is 200
							network_type="full", 
                            input_directory="F:/projects/Share/PPIN/STRING/STRINGdb/human/full/") 
 string_db
 # version 12.0  / 11.5
 # proteins: 19699    / for mouse: 22048
 # interactions:  7533072   / for mouse 6859804
 
 graph_list <- list()
 # rather than build for steady PPI_cats
 # for(i in setdiff(names(DEG), names(CTS))){
 # build for traditional up-regulated markers
 # build for traditional up-regulated markers
for(i in names(DEG)){
	# Filter differentially expressed genes based on logFC and FDR thresholds
	diff_exp = markers.up[[i]] 
	## Add gene symbols as a column and map to STRING IDs
	diff_exp$symbol = diff_exp$names
    mapped <- string_db$map(diff_exp, "symbol", removeUnmappedRows = TRUE )
    mapped %>% dim
    # [1] 84  21
    length(unique(mapped$symbol))
    #[1] 84
    # Get the subnetwork for the mapped genes
	hits <- mapped$STRING_id
    graph <- string_db$get_subnetwork(hits) 
	edge_attr_names(graph) # "combined_score"
	
    ## Handle missing nodes (if any)  ## e.g., due to graph missing the [1] "GM1673" in DEG[['5']]
    if(vcount(graph) < length(hits)) {
		missing_node = setdiff(hits, V(graph)$name)
		if(length(missing_node)>0) mapped = mapped[-which(mapped$STRING_id %in% missing_node),]
    }
    # Verify mapping integrity
	stopifnot(all(mapped[match(V(graph)$name, mapped$STRING_id), ]$STRING_id == V(graph)$name))
   
	# Add node attributes: replace STRING IDs with gene symbols	
    V(graph)$name =  mapped[match(V(graph)$name, mapped$STRING_id),]$symbol
	
	# Add node weights from logFC values and additioal node attributes
    V(graph)$weight = diff_exp[match(V(graph)$name, diff_exp$symbol),]$scores   #summary.logFC if scanpy using t.test
    V(graph)$FDR = diff_exp[match(V(graph)$name, diff_exp$symbol),]$pvals_adj    # FDR of scanpy_wilcox marker identification
    # E(graph)   
   
    # Add edge weights from STRING interaction scores #!!!!!!!!!!!!!!!!!
	# STRING scores are between 0-1000, normalize to 0-1 for better visualization
	E(graph)$weight <- E(graph)$combined_score / 1000

	# Fix case issues if needed (for mouse genes)
    if(all(mapped$symbol %in% toupper(DEG[[i]]))) V(graph)$name = DEG[[i]][match(V(graph)$name, toupper(DEG[[i]]))] 
    
    # Optional: Set visual attributes based on weights
	V(graph)$size <- scales::rescale(V(graph)$weight, to = c(3, 15))
	E(graph)$width <- scales::rescale(E(graph)$weight, to = c(0.5, 3))  # useless, using weight directly otherwise can't compare across clusters
	V(graph)$color <- colorRampPalette(c("lightblue", "red"))(100)[
	  cut(V(graph)$weight, breaks = 100)
	]

	graph_list[[i]] = graph
 }
  names(graph_list) = paste0('HiG_', names(DEG))
 
 # instead build for transitory PPI_cats only with the significantly up-regulated CTS genes !!!!! 
 # build for (up-regulated_marker intersecting CTS)
 for(i in names(CTS)){
   if(grepl('.',i, fixed=T)) j = unlist(strsplit(i, split='.', fixed=T))[1]  else j = i
   diff_exp = markers.up[[j]][which(markers.up[[j]]$names %in% CTS[[i]]), ]
   
   diff_exp$symbol = diff_exp$names
   mapped <- string_db$map(diff_exp , "symbol", removeUnmappedRows = TRUE )
   mapped %>% dim
   # [1] 84  21
   length(unique(mapped$symbol))
   #[1] 84
   
   hits <- mapped$STRING_id
   # string_db$plot_network( hits )
   
   graph <- string_db$get_subnetwork(hits) #t
   # translate STRING_id to symbol
   all(mapped[match(V(graph)$name, mapped$STRING_id),]$STRING_id == V(graph)$name) # TRUE
   V(graph)$name =  mapped[match(V(graph)$name, mapped$STRING_id),]$symbol
   # The 'scores' column reflects the strength of the differential expression for each gene based on the Wilcoxon rank-sum test. A high score suggests that the gene's expression is significantly different between the groups under comparison.
   V(graph)$weight = diff_exp[match(V(graph)$name, diff_exp$symbol),]$scores  
   V(graph)$FDR = diff_exp[match(V(graph)$name, diff_exp$symbol),]$pvals_adj
    
	# Add edge weights from STRING interaction scores #!!!!!!!!!!!!!!!!!
	# STRING scores are between 0-1000, normalize to 0-1 for better visualization
	E(graph)$weight <- E(graph)$combined_score / 1000
   
 	# Fix case issues if needed (for mouse genes)
    if(all(mapped$symbol %in% toupper(DEG[[i]]))) V(graph)$name = DEG[[i]][match(V(graph)$name, toupper(DEG[[i]]))] 
    
    # Optional: Set visual attributes based on weights
	V(graph)$size <- scales::rescale(V(graph)$weight, to = c(3, 15))
	E(graph)$width <- scales::rescale(E(graph)$weight, to = c(0.5, 3))
	V(graph)$color <- colorRampPalette(c("lightblue", "red"))(100)[
	  cut(V(graph)$weight, breaks = 100)
	]
  
   graph_list[[paste0('HiGCTS_',i)]] = graph
   
 }
 
 ## lastly, build for CTS
 # refer to 6.3_DE.statistics_CTS.R 
 markers.up_all = readRDS(here::here("examples", "cardiac", "GSE175634", "data", "DEG.wc_nofiltering.rds"))
 names(markers.up_all)[which(names(markers.up_all)=='3')] = 'CP'
 names(markers.up_all)[which(names(markers.up_all)=='8')] = 'muscle'
 names(markers.up_all)[which(names(markers.up_all)=='11')] = 'endoderm'


 for(i in names(CTS)){
   if(grepl('.',i, fixed=T)) j = unlist(strsplit(i, split='.', fixed=T))[1]  else j = i
   diff_exp = markers.up_all[[j]][which(markers.up_all[[j]]$names %in% CTS[[i]]), ]
  
   diff_exp$symbol = diff_exp$names
   mapped <- string_db$map(diff_exp , "names", removeUnmappedRows = TRUE )
   mapped %>% dim
   # [1] 84  21
   length(unique(mapped$names))
   #[1] 84
   
   hits <- mapped$STRING_id
   # string_db$plot_network( hits )
   
   graph <- string_db$get_subnetwork(hits) #t
   # translate STRING_id to symbol
   all(mapped[match(V(graph)$name, mapped$STRING_id),]$STRING_id == V(graph)$name) # TRUE
   V(graph)$name =  mapped[match(V(graph)$name, mapped$STRING_id),]$symbol
   # The 'scores' column reflects the strength of the differential expression for each gene based on the Wilcoxon rank-sum test. A high score suggests that the gene's expression is significantly different between the groups under comparison.
   V(graph)$weight = diff_exp[match(V(graph)$name, diff_exp$symbol),]$scores  
   V(graph)$FDR = diff_exp[match(V(graph)$name, diff_exp$symbol),]$pvals_adj
   
	# Add edge weights from STRING interaction scores #!!!!!!!!!!!!!!!!!
	# STRING scores are between 0-1000, normalize to 0-1 for better visualization
	E(graph)$weight <- E(graph)$combined_score / 1000
   
    # Optional: Set visual attributes based on weights
	V(graph)$size <- scales::rescale(V(graph)$weight, to = c(3, 15))
	E(graph)$width <- scales::rescale(E(graph)$weight, to = c(0.5, 3))
	V(graph)$color <- colorRampPalette(c("lightblue", "red"))(100)[
	  cut(V(graph)$weight, breaks = 100)
	]
	
	graph_list[[paste0('CTS_',i)]] = graph
   
 }
 
 names(graph_list)
 # [1] "UpG_0"           "UpG_1"           "UpG_2"           "UpG_CP"         
 # [5] "UpG_4"           "UpG_5"           "UpG_6"           "UpG_7"          
 # [9] "UpG_muscle"      "UpG_9"           "UpG_10"          "UpG_endoderm"   
# [13] "UpG_12"          "HiGCTS_muscle"   "HiGCTS_endoderm" "HiGCTS_CP"      
# [17] "HiGCTS_CP.1"     "CTS_muscle"      "CTS_endoderm"    "CTS_CP"         
# [21] "CTS_CP.1"      
 
  # for(i in 1:length(graph_list)){
    # tmp = vertex_attr_names(graph_list[[i]])
    # if(tmp[2] == "logFC_ave") {
		# old_name <- "logFC_ave"
		# vertex_attr(graph_list[[i]], "scores") <- vertex_attr(graph_list[[i]], old_name)
		# graph_list[[i]] <- delete_vertex_attr(graph_list[[i]], old_name)
	# }   
 # }

  # new_names <- gsub("UpGCTS", "HiGCTS", names(graph_list))
 saveRDS(graph_list, file= 'GSE175634_STRING_graph_perState_notsimplified.rds')  #!!!!!!!!!!!!!!!!!!!
 
 graph_list <- readRDS( file= 'GSE175634_STRING_graph_perState_notsimplified.rds')  
 graph_list <- lapply(graph_list, simplify) #!!!!!!!!!!!!!!!!!!!
 
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
# HiG_1 HiG_2 HiG_4 HiG_6 HiG_7 HiG_9 
    # 2     3     5     7     8    10 

# Show actual edges for duplicated vertices
g1 <- graph_list[['HiG_1']]
vertex_names <- V(g1)$name
(duplicated_names <- unique(vertex_names[duplicated(vertex_names)]))
#  "H3F3B"
if(length(duplicated_names) > 0) {
  for(dup_name in duplicated_names) {
    dup_indices <- which(vertex_names == dup_name)
    cat("Edges for duplicated vertex '", dup_name, "':\n", sep = "")
    all(incident(g1, dup_indices[1], mode = "all") == incident(g1, dup_indices[2], mode = "all"))
	
	edges <- incident(g1, dup_indices[1], mode = "all")
    edge_list1 <- get.edgelist(g1)[edges, ]
    edge_list1 %>% dim  #[1] 290   2
	weights1 <- E(g1)[edges]$weight
	
	edges <- incident(g1, dup_indices[2], mode = "all")
    edge_list2 <- get.edgelist(g1)[edges, ]
    edge_list2 %>% dim# [1] 64  2
	weights2 <- E(g1)[edges]$weight
  
    all(edge_list2[,2] %in% edge_list1[,2])  # FALSE !!
   }
}