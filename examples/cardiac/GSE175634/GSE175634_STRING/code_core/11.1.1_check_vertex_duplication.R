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

code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_STRING", "code_core")
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)

## WARNING: the duplicate-vertex fixes below (`if (i=='1') mapped = mapped[-29,]`
## etc.) are hardcoded ROW POSITIONS, not gene-name-keyed — they are only valid
## for this exact STRINGdb version ("12.0") + score_threshold (200) + gene
## symbol set. Re-running against a different STRING db version or DEG source
## would silently remove the wrong rows. Preserved as-is (same fragile-but-
## accepted pattern as GSE87038_STRING's own 11.1.1) since the db version is
## pinned and this has been verified against the current data.

graph_list <- readRDS( file= 'GSE175634_STRING_graph_perState_notsimplified.rds')
graph_list <- lapply(graph_list, simplify) #!!!!!!!!!!!!!!!!!!!

N = sapply(graph_list, vcount)

## degree distribution was not discussed again because it is not in the manuscript
#####################  #################
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


### why "H3F3B" is duplicated in 'HiG_1'?
# refer to ../xxx_score200/11.1.1_CTS_cardiac_network_degreeDistribution.R
library("STRINGdb")
packageVersion('STRINGdb') # '2.14.0'
string_db <- STRINGdb$new( version="12.0", species=db_species,
                            score_threshold = score_threshold,
							network_type="full",
                            input_directory=paste0(shared_data_path, "/PPIN"))
string_db

load(here::here("examples", "cardiac", "GSE175634", "data", "DEG.wc_padj0.01_score40.RData"))
lapply(DEG.wc, nrow) %>% unlist
   # 0    1    2    3    4    5    6    7    8    9   10   11   12 
 # 344 1224  257  152  110  350  110  117  119 1633   81   44   25 
DEG = list()
for(i in 1:length(DEG.wc)) DEG[[i]] = DEG.wc[[i]]$names %>% unique
lengths(DEG)
#  [1]  344 1224  257  152  110  350  110  117  119 1633   81   44   25
names(DEG) = names(DEG.wc) #  0:(length(DEG)-1)

length(DEG[['1']]) # 1224
any(duplicated(DEG[['1']])) #[1] FALSE
N['HiG_1']  # 1156

## build PPIN again and track back the correct number of edges for the duplciated gene ####

new_object <- new.env()                # Create a new environment
load(here::here("examples", "cardiac", "GSE175634", "data", "DEG.wc_padj0.01_score40.RData"), envir = new_object)  # Load .RData into this environment
markers.up <- new_object[[ls(new_object)[1]]]

graph_list <- readRDS( file= 'GSE175634_STRING_graph_perState_notsimplified.rds')  
graph_list <- lapply(graph_list, simplify) #!!!!!!!!!!!!!!!!!!!

N = sapply(graph_list, vcount)
graphs_with_duplicates <- sapply(graph_list, function(g) {
  vertex_names <- V(g)$name
  if(is.null(vertex_names)) {
    # If no names, use vertex indices
    vertex_names <- V(g)
  }
  any(duplicated(vertex_names))
})
which(graphs_with_duplicates)
# HiG_1 HiG_2 HiG_4 HiG_6 HiG_7 HiG_9 
    # 2     3     5     7     8    10 

correct_n_edges = NULL
# colnames(correct_n_edges) = names(which(graphs_with_duplicates))
# rownames(correct_n_edges) = c('graph_id','names', 'n_edge', 'STRING_id')

for(i in substr(names(which(graphs_with_duplicates)), 5,10)){
	diff_exp = markers.up[[i]]
	mapped <- string_db$map(diff_exp, "names", removeUnmappedRows = TRUE )
	mapped %>% dim # [1] 1156  6
	length(unique(mapped$names))  #[1] 1155
	length(unique(mapped$STRING_id))  # 1156
	(dup_name = mapped[which(duplicated(mapped$names)),]$names ) #  "H3F3B"
	cat('analyzing',i,'\n')
	(mapped[which(mapped$names %in% dup_name),] )  #  "H3F3B"
   # names   scores logfoldchanges pvals pvals_adj            STRING_id  genecard(manually fetch)
# 29 H3F3B 121.2121            NaN     0         0 9606.ENSP00000355780	H3-3A
# 30 H3F3B 121.2121            NaN     0         0 9606.ENSP00000254810  H3-3B OK
	if(i=='1') mapped = mapped[-29,]
       # names   scores logfoldchanges pvals pvals_adj            STRING_id
# 101     GARS 67.71815            NaN     0         0 9606.ENSP00000373918	GARS1	OK
# 102     GARS 67.71815            NaN     0         0 9606.ENSP00000371253	GART
# 113 HIST1H4C 65.46398            NaN     0         0 9606.ENSP00000367034  H4C3  higher ranked alias
# 114 HIST1H4C 65.46398            NaN     0         0 9606.ENSP00000244537 H4C6
# 142     SARS 60.43468            NaN     0         0 9606.ENSP00000472465	ENSG00000269547
# 143     SARS 60.43468            NaN     0         0 9606.ENSP00000472847	SARS2
# 144     SARS 60.43468            NaN     0         0 9606.ENSP00000358939	SARS1	OK
	if(i=='2') mapped = mapped[-c(102, 114, 142, 143),]
       # names   scores logfoldchanges pvals pvals_adj            STRING_id
# 102 HIST1H4C 42.83387            NaN     0         0 9606.ENSP00000367034	H4C3  higher ranked alias
# 103 HIST1H4C 42.83387            NaN     0         0 9606.ENSP00000244537 H4C6
	if(i=='4') mapped = mapped[-c(103),]
      # names   scores logfoldchanges pvals pvals_adj            STRING_id
# 18 HIST1H4C 79.42637            NaN     0         0 9606.ENSP00000244537 H4C6
# 19 HIST1H4C 79.42637            NaN     0         0 9606.ENSP00000367034  H4C3  higher ranked alias
# 45    H3F3B 57.50134            NaN     0         0 9606.ENSP00000355780	H3-3A
# 46    H3F3B 57.50134            NaN     0         0 9606.ENSP00000254810  H3-3B OK
	if(i=='6') mapped = mapped[-c(18,45),]
       # names   scores logfoldchanges pvals pvals_adj            STRING_id
# 22  HIST1H4C 69.50114            NaN     0         0 9606.ENSP00000244537 H4C6
# 23  HIST1H4C 69.50114            NaN     0         0 9606.ENSP00000367034	H4C3  higher ranked alias
# 44      GARS 60.83009            NaN     0         0 9606.ENSP00000373918	GARS1	OK
# 45      GARS 60.83009            NaN     0         0 9606.ENSP00000371253	GART
# 103     SARS 40.59730            NaN     0         0 9606.ENSP00000472465	ENSG00000269547
# 104     SARS 40.59730            NaN     0         0 9606.ENSP00000358939	SARS1	OK
# 105     SARS 40.59730            NaN     0         0 9606.ENSP00000472847	SARS2	
	if(i=='7') mapped = mapped[-c(22,45, 103, 105),]
     # names   scores logfoldchanges pvals pvals_adj            STRING_id
# 833  H3F3A 54.18238            NaN     0         0 9606.ENSP00000355780	H3-3A OK
# 834  H3F3A 54.18238            NaN     0         0 9606.ENSP00000254810  H3-3B
# 1255 H3F3B 44.90171            NaN     0         0 9606.ENSP00000254810  H3-3B OK
# 1256 H3F3B 44.90171            NaN     0         0 9606.ENSP00000355780	H3-3A
	if(i=='9') mapped = mapped[-c(834, 1256),]
	
   hits <- mapped$STRING_id
   graph <- string_db$get_subnetwork(hits) #t
   all(mapped[match(V(graph)$name, mapped$STRING_id),]$STRING_id == V(graph)$name) # TRUE  
   V(graph)$name =  mapped[match(V(graph)$name, mapped$STRING_id),]$names   

   for(j in 1:length(dup_name)){
	edges = incident(graph, dup_name[j], mode = "all")
	n = get.edgelist(graph)[edges, ]
	
	res = data.frame('graph_id'=paste0('HiG_',i),
				'names' = dup_name[j], 
				'n_edge' = nrow(n) , 
				'STRING_id' = subset(mapped, names==dup_name[j])$STRING_id)
	 
	if(is.null(correct_n_edges)) correct_n_edges = res else {
		correct_n_edges = rbind(correct_n_edges, res) 
	}
   }
}   
correct_n_edges   
   # graph_id    names n_edge            STRING_id
# 1     HiG_1    H3F3B    578 9606.ENSP00000254810
# 2     HiG_2     GARS    100 9606.ENSP00000373918
# 3     HiG_2 HIST1H4C     78 9606.ENSP00000367034
# 4     HiG_2     SARS     84 9606.ENSP00000358939
# 5     HiG_2     SARS     84 9606.ENSP00000358939
# 6     HiG_4 HIST1H4C     36 9606.ENSP00000367034
# 7     HiG_6 HIST1H4C     42 9606.ENSP00000367034
# 8     HiG_6    H3F3B    104 9606.ENSP00000254810
# 9     HiG_7 HIST1H4C     28 9606.ENSP00000367034
# 10    HiG_7     GARS     56 9606.ENSP00000373918
# 11    HiG_7     SARS     46 9606.ENSP00000358939
# 12    HiG_7     SARS     46 9606.ENSP00000358939
# 13    HiG_9    H3F3A    154 9606.ENSP00000355780
# 14    HiG_9    H3F3B    720 9606.ENSP00000254810

# saveRDS(correct_n_edges, file='correct_n_edges_HiG_STRING2.14.0.rds')
# getwd()
# #[1] "F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weight"

################################################################
## remove duplciated vertex directly from un-simplified graph ##
################################################################
graph_list <- readRDS( file= 'GSE175634_STRING_graph_perState_notsimplified.rds')  
#graph_list <- lapply(graph_list, simplify) #!!!!!!!!!!!!!!!!!!!

N = sapply(graph_list, vcount)
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
	
# correct_n_edges = readRDS('correct_n_edges_HiG_STRING2.14.0.rds')
correct_n_edges$vetex_index_to_remove = 0  # initilization

for(i in names(which(graphs_with_duplicates) )){
	g_name = i
	g = graph_list[[g_name]]
	vertex_names <- V(g)$name
    (duplicated_names <- unique(vertex_names[duplicated(vertex_names)]))
	
	for(dup_name in duplicated_names) {
		(dup_indices <- which(vertex_names == dup_name))
		n = array(dim = length(dup_indices))
		for( j in 1:length(dup_indices)){
			edges <- incident(g, dup_indices[j], mode = "all")
			edge_list <- get.edgelist(g)[edges, ]
			n[j] = edge_list %>% nrow  #[1] 290   2
		}
		x = which(correct_n_edges$graph_id==g_name & correct_n_edges$names %in% dup_name)
		keep_index = which.min(n - correct_n_edges$n_edge[x]) 	    
		correct_n_edges$vetex_index_to_remove[x] = dup_indices[-keep_index]
	}		
}
if(any(correct_n_edges$vetex_index_to_remove == 0)) {
	correct_n_edges$vetex_index_to_remove[which(correct_n_edges$vetex_index_to_remove == 0)]=NA
	}

write.table(correct_n_edges, file='correct_n_edges_HiG_STRING2.14.0.txt', sep='\t')
saveRDS(correct_n_edges, file='correct_n_edges_HiG_STRING2.14.0.rds')
getwd()
#[1] "F:/projects/scRNA/results/cardiac_CTS_GRN/GSE175634_iPSC_CM_weighted_v9"

 