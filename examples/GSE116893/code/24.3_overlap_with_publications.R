# source('F:/projects/scRNA/source/cardiac_CTS_GRN/GSE87038_Pijuan2019_v9/C8_vs_C16/24.3_overlap_with_publications.R')
CF_cluster = '16'  # mesenchymal/ECM‑high (“CF‑like”).
if(length(CF_cluster)==1) EMT_ID  = CF_cluster else EMT_ID = 'C3C16C18C19'


rebuild_mat = FALSE
source('F:/projects/scRNA/source/cardiac_CTS_GRN/GSE87038_Pijuan2019_v9/24.0_acat_load_input_clean_v2.R')

celltype_col  # 'leiden_0.5'

names(graph_list)
#  [1] "HiG_1"       "HiG_2"       "HiG_3"       "HiG_4"       "HiG_5"      
#  [6] "HiG_6"       "HiG_9"       "HiG_10"      "HiG_12"      "HiG_14"     
# [11] "HiG_17"      "HiG_18"      "HiG_19"      "HiG_7"       "HiG_11"     
# [16] "HiG_15"      "HiG_16"      "HiG_13"      "HiG_8"       "HiGCTS_7"   
# [21] "HiGCTS_11"   "HiGCTS_15"   "HiGCTS_16"   "HiGCTS_16.1" "HiGCTS_8"   
# [26] "CTS_7"       "CTS_11"      "CTS_15"      "CTS_16"      "CTS_16.1"   
# [31] "CTS_8"   

names(DEG)
#   [1] "1"  "2"  "3"  "4"  "5"  "6"  "9"  "10" "12" "14" "17" "18" "19" "7" 
#  [15] "11" "15" "16" "13" "8" 


lengths(CTS)
#   7   11   15   16    8 16.1 
#  32   52   67   40   54   79 

(updir = getwd())
#[1] "F:/projects/scRNA/results/cardiac_CTS_GRN/GSE87038_Pijuan2019_v9/GSE181346_heart_scATAC"
fpath = paste0('F:/projects/scRNA/results/cardiac_CTS_GRN/GSE87038_Pijuan2019_v9/validation_C8vs', EMT_ID)
if(!file.exists(fpath)) dir.create(fpath)
setwd(fpath)

library(igraph)
library(dplyr)
library(igraph)
library(openxlsx)
library(gplots)
library(ggExtra)




####################################################
###########overlap with published gene sets
## refer to F:\projects\scRNA\source\cardiac_CTS_GRN\GSE175634_iPSC_CM_weighted_v9\24.4.0_GS_publications_collection.R
library(dplyr)
library(tidyr)
source('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v10.R')

GP_path = 'D:/projects/DS/data/TF_perturb_CM/doc/'


df = read.xlsx(xlsxFile = paste0(GP_path, 'media-2.xlsx'), sheet=8, startRow = 2)  # 250 programs
dim(df)  #[1] 300  251
df = df[, -1]
head(df[,1:10], 3)
      # 1       2       3       4     5     6       7      8     9      10
# 1  FBN2 CACNA1D  NCCRP1  MT-CO3 SFRP1 GRIK1 ZNF385D COL1A2  LMO2   H3F3B
# 2 LAMA2   RSPO3 ZNF385B  MT-CO2 RAB3C   CKB   PRR16  ITGA8 DPP10 GADD45G
# 3  MYL3     LBH   TBX18 MT-ATP8 PDE3A MFGE8   FABP5 COL1A1 CADM1  PHLDA1

length(unique(unlist(df))) # 25311

GS = readRDS(paste0(shared_input_path, 'published_GS_collection.rds'))
lengths(GS)
 
length(union(rownames(sce), unique(unlist(df)))) # 28549
universe = union(rownames(sce), unique(unlist(df)))
length(universe) # 28549
universe = union(unlist(GS), universe)
length(universe) # 29225


  
GS_long = data.frame()
for(i in names(GS)){
	for(j in GS[[i]]){
		GS_long = rbind(GS_long, data.frame(scRNA_signature = i, gene = j))
	}
}
dim(GS_long) #[1] 11748     2

 
enrich_res = list()	
for(i in CTS_name){  
	enrich_res[[paste0('CTS_',i)]] = cts_overlap_by_cluster(GS_long, CTS[[i]], gene_col = 'gene',
					universe = universe, cluster_col = 'scRNA_signature')
}

for(i in c(CP_cluster, CF_cluster, CM_cluster)){
 enrich_res[[paste0('HiG_',i)]] <- cts_overlap_by_cluster(GS_long, DEG[[i]], gene_col = 'gene',
					universe = universe, cluster_col = 'scRNA_signature')
}

## how about Mutrans TF genes?
for(i in c("Elorbany_TD_A10_to_9" , "Elorbany_TD_A10_to_4")){
	enrich_res[[i]] <- cts_overlap_by_cluster(GS_long, GS[[i]], gene_col = 'gene',
					universe = universe, cluster_col = 'scRNA_signature')
}

## how about Mutrans Regulon target genes?
for (i in grep('Regulon', names(GS), value=T)){
	enrich_res[[i]] <- cts_overlap_by_cluster(GS_long, GS[[i]], gene_col = 'gene',
					universe = universe, cluster_col = 'scRNA_signature')
}


tmp <- enrich_res %>%
  dplyr::bind_rows(.id = "scRNA_signature") %>%   # keeps list names if present
  dplyr::filter(odds_ratio > 2, adj_p < 0.05, n_overlap>=3)
dim(tmp)  # [1]  228 11  


write.table(tmp, file= 	'GS_CTS_or_DEG_padj0.05_OR2_minOverlap3.tsv', sep='\t', row.names=F)		


######### plot 
library(dplyr)
library(ggplot2)
library(forcats)
library(tidytext)

tmp = read.table( file= 'GS_CTS_or_DEG_padj0.05_OR2_minOverlap3.tsv', sep='\t', header=T)
## simplify the GP to be compared
# GP5 is the most strongly enriched for CHD genes, 
# GP1 epresent cardiac differentiation and heart construction, 
# GP28 represents CHD-related extracellular matrix.
tmp = tmp[-which(tmp$GS %in% c('Takeuchi_gp10', 'Takeuchi_gp2')),]

### order by odds ratio
my_dot_plot = function(tmp2, ntop=10, ncol=2, CTS_name, CP_cluster, CM_cluster, Levels=NULL){
    if(is.null(Levels)) {
		Levels = c(paste0('CTS_', CTS_name), paste0('HiG_',CP_cluster), paste0('HiG_',CM_cluster), paste0('HiG_',CF_cluster))
		if(any(grep('Elorbany_TD', tmp2$scRNA_signature))) Levels = c(Levels, "Elorbany_TD_A10_to_4", "Elorbany_TD_A10_to_9")
		if(any(grep('Regulon', tmp2$scRNA_signature))) Levels = c(Levels, grep('Regulon', tmp2$scRNA_signature, value=T) %>% unique)
	}
	
	df_top <- tmp2 %>% 
	filter(!is.infinite(odds_ratio)) %>%
	group_by(scRNA_signature) %>%
	arrange(adj_p, .by_group = TRUE) %>%
	slice_head(n = ntop) %>%      # top 5 smallest adj_p
	ungroup() %>%
	mutate(
		neglog10 = -log10(adj_p)
	)%>%
	mutate(
		GS = reorder_within(GS, neglog10, scRNA_signature)
	)
	
	df_plot <- df_top %>%
	group_by(scRNA_signature) %>%
	mutate(GS = fct_reorder(GS, neglog10)) %>%
	ungroup() %>%
	mutate(scRNA_signature = factor(scRNA_signature, 
			levels = Levels))

	p <- ggplot(df_plot, aes(x = odds_ratio, y = GS)) +
	geom_point(aes(size = n_overlap,
					color = neglog10)) +
	facet_wrap(~ scRNA_signature, scales = "free", ncol=ncol) +
	scale_y_reordered() +
	scale_color_gradient(
		low = "#3B4CC0",
		high = "#B40426",
		name = expression(-log[10](FDR))
	) +
	scale_size(range = c(2, 5), name = "Intersect count") +
	labs(
		x = paste("Odds ratio (N=", length(universe), ", OR>2, padj<0.05, overlap>=3)"),
		y = NULL
	) +
	theme_bw(base_size = 12) +
	theme(
		strip.background = element_blank(),
		panel.grid.major.y = element_line(color = "grey90"),
		panel.grid.minor = element_blank()
	)

	return(p)
}

# order by GS source
my_dot_plot2 = function(tmp2, ntop=10, ncol=2, CTS_name, CP_cluster, CM_cluster, Levels=NULL){ 

    if(is.null(Levels)) {
        Levels = c(paste0('CTS_', CTS_name),
                   paste0('HiG_', CP_cluster),
                   paste0('HiG_', CM_cluster),
                   paste0('HiG_', CF_cluster))
        if(any(grep('Elorbany_TD', tmp2$scRNA_signature))) 
            Levels = c(Levels, "Elorbany_TD_A10_to_4", "Elorbany_TD_A10_to_9")
        if(any(grep('Regulon', tmp2$scRNA_signature))) 
            Levels = c(Levels, grep('Regulon', tmp2$scRNA_signature, value=TRUE) %>% unique)
    }

    df_top <- tmp2 %>% 
        filter(!is.infinite(odds_ratio)) %>%
        group_by(scRNA_signature) %>%
        arrange(adj_p, .by_group = TRUE) %>%
        slice_head(n = ntop) %>%
        ungroup() %>%
        mutate(
            neglog10 = -log10(adj_p)
        )

    ## reorder GS rows within each facet by source group
    gs_levels <- df_top %>%
        distinct(GS) %>%
        mutate(
            source_group = case_when(
                grepl("^Maven2023", GS) ~ 1,
                grepl("^Elorbany", GS) ~ 2,
                grepl("^Xie", GS) ~ 3,
                grepl("^Takeuchi", GS) ~ 4,
                grepl("^PCGC", GS) ~ 5,
                TRUE ~ 6
            )
        ) %>%
        arrange(source_group, GS) %>%
        pull(GS)

    df_plot <- df_top %>%
        mutate(
            GS = factor(GS, levels = rev(gs_levels)),
            scRNA_signature = factor(scRNA_signature, levels = Levels)
        )

    p <- ggplot(df_plot, aes(x = odds_ratio, y = GS)) +
        geom_point(aes(size = n_overlap, color = neglog10)) +
        facet_wrap(~ scRNA_signature, scales = "free", ncol=ncol) +
        scale_color_gradient(
            low = "#3B4CC0",
            high = "#B40426",
            name = expression(-log[10](FDR))
        ) +
        scale_size(range = c(2, 5), name = "Intersect count") +
        labs(
            x = paste("Odds ratio (N=", length(universe), 
                      ", OR>2, padj<0.05, overlap>=3)"),
            y = NULL
        ) +
        theme_bw(base_size = 12) +
        theme(
            strip.background = element_blank(),
            panel.grid.major.y = element_line(color = "grey90"),
            panel.grid.minor = element_blank()
        )

    return(p)
}


# excluding up/downregulated genes , Regulon target genes, or Elorbany_TD
tmp2 = subset(tmp, !(grepl('up', tmp$GS) | grepl('dn', tmp$GS) | 
	grepl('Regulon', tmp$scRNA_signature) | grepl('Elorbany_TD', tmp$scRNA_signature))) 
dim(tmp2) # 64 11

 
pdf(file='GS_CTS_or_DEG_padj0.05_top15.pdf', height=7, width=10)  #!!!!!!!!!!!!

p = my_dot_plot(tmp2, ntop=15, ncol=2, CTS_name=CTS_name, CP_cluster=CP_cluster, CM_cluster=CM_cluster)  # excluding up/downregulated genes and Regulon target genes
print(p)
dev.off()

### replot CTS specifically 
tmp2_CTS = subset(tmp2, scRNA_signature %in% c(paste0('CTS_',CTS_name)))
dim(tmp2_CTS) # 11   11

pdf(file='GS_CTS_padj0.05_top15_ordered_by_GS_source.pdf', height=7, width=10)  #!!!!!!!!!!!!
p = my_dot_plot2(tmp2_CTS, ntop=15, ncol=2, CTS_name=CTS_name, CP_cluster=CP_cluster, CM_cluster=CM_cluster)  # excluding up/downregulated genes and Regulon target genes
print(p)
dev.off()



pdf(file='GS_CTS_padj0.05_top15.pdf', height=7, width=10)  #!!!!!!!!!!!!
p = my_dot_plot(tmp2, ntop=15, ncol=1, CTS_name=CTS_name, CP_cluster=CP_cluster, CM_cluster=CM_cluster)  # excluding up/downregulated genes and Regulon target genes
print(p)
dev.off()



################################
## new bubble plot 
tmp = read.table( file= 'GS_CTS_or_DEG_padj0.05_OR2_minOverlap3.tsv', sep='\t', header=T)
dim(tmp)  #[1]  228   11
colnames(tmp)
 # [1] "scRNA_signature" "GS"              "n_cluster_genes" "n_CTS"           "n_overlap"       "overlap"        
 # [7] "odds_ratio"      "p_value"         "adj_p"           "neglog10_p"      "neglog10_adj_p" 

unique(tmp$GS)
int = c(  "Elorbany_C1_RegulonWT1_target" ,  "Takeuchi_gp28"     ,  
 "Elorbany_C5_RegulonMEF2A_target" ,  "Takeuchi_gp1"  ,
"Elorbany_C3_RegulonPOU5F1_target",   "Takeuchi_gp5"    )
unique(tmp$scRNA_signature)

  
tmp = subset(tmp, GS %in% int & scRNA_signature %in% c("CTS_8", "HiG_8"  , paste0("HiG_",CF_cluster) , "HiG_17" ))
dim(tmp) # 15 11
tmp$scRNA_signature = factor(tmp$scRNA_signature, levels =  c("CTS_8", "HiG_8"  , paste0("HiG_",CF_cluster)  , "HiG_17" ))
 
tmp$GS[which(tmp$GS == "Takeuchi_gp1")] = 'Takeuchi_CM-differentiation-like CHD GP1'
tmp$GS[which(tmp$GS == "Takeuchi_gp28")] = 'Takeuchi_Mesenchymal-like CHD GP28'
tmp$GS[which(tmp$GS == "Takeuchi_gp5")] = 'Takeuchi_Unknown leading CHD GP5'
tmp$GS = factor(tmp$GS, levels = c("Elorbany_C3_RegulonPOU5F1_target",   "Takeuchi_Unknown leading CHD GP5"  ,
  "Elorbany_C1_RegulonWT1_target" ,  "Takeuchi_Mesenchymal-like CHD GP28"     ,  
 "Elorbany_C5_RegulonMEF2A_target" ,  "Takeuchi_CM-differentiation-like CHD GP1"  
  ) )

g = ggplot(tmp, aes(x = scRNA_signature, y =GS)) +
          #   y = reorder(GS,odds_ratio))) +
  geom_point(aes(size = n_overlap, color = neglog10_adj_p)) +
  scale_color_viridis_c() +
  theme_classic() +
  ylab('Enriched GS') +
  ggtitle('Significance:  OR>2, padj<0.05, overlap>=3')
  
pdf(file='bubble_GS_CTS_padj0.05_top15.pdf',width=7, height=3.5)
print(g)
dev.off()