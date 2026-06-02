setwd("F:/projects/scRNA/results/cardiac_CTS_GRN/CHD/")
library(dplyr)
library(igraph)
library(gplots)
library(ggplot2)
library(openxlsx)

## TIPS identification

pull_path = 'F:/projects/TIPS/doc/'
pull_df = read.xlsx(paste0(pull_path, 'STable_final_2026.xlsx'), sheet=16)  # for C8vsC18 read this 
dim(pull_df)  # 95  15
head(pull_df, 3)
     # Database    CTS_ID linkeage     x     y                 w_CP
# 1 IbarraSoria cardiac.a       CM GATA6  TBX5  0.14273739867444099
# 2 IbarraSoria cardiac.a       CM  CBFB GATA6  3.65114162434605E-2
# 3 IbarraSoria cardiac.a       CM  CBFB  TBX5 -1.28698510207887E-2
             # w_decendent                  delta             abs_delta
# 1  5.3873131230046403E-2 -8.8864267444394796E-2 8.8864267444394796E-2
# 2 -5.6903109139440902E-3 -4.2201727157404598E-2 4.2201727157404598E-2
# 3 -2.1015662268983001E-4    1.26596943980988E-2   1.26596943980988E-2
  # direction  status rank               TF_highConf      motif  NES
# 1  decrease changed    1                      <NA>       <NA> <NA>
# 2  decrease changed    2 CBFB (directAnnotation).  hdpi__CBFB 4.67
# 3  increase changed    3 CBFB (directAnnotation).  hdpi__CBFB 4.67
         

#############    testing on the SC dataset using author's louvain and Holly's leiden_0.8clustering results ###################################################
sig_map <- data.frame(
  db = c( "IbarraSoria", "IbarraSoria", "Pijuan_Sala", "Elorbany"),
  signature = c("cardiac.a", "cardiac.a_noTBX5", "8", "CP.1"),
  cf_col = c("CF_pullcardiac.a", "CF_pullcardiac.a", "CF_pull8", "CF_pullCP.1"),
  #cm_col = c("CM_pullcardiac.a", "CM_pull8_wincreased_nodes", "CM_pullCP.1_wincreased_nodes"), ## when test C8vsC16
  cm_col = c("CM_pullcardiac.a", "CM_pullcardiac.a_noTBX5", "CM_pull8", "CM_pullCP.1_wincreased_nodes"), ## when test C8vsC16
  stringsAsFactors = FALSE
)

sig_map = sig_map[-2,]

PCGC_path = 'D:/projects/DS/data/CHD_NDD/PCGC/'
CHD1 = read.csv(paste0(PCGC_path, 'Table1_99CHDgene.txt'), sep='\t')
CHD2 = read.csv(paste0(PCGC_path, 'STable1_295CHDgene.txt'), sep='\t')
dim(CHD1)  # 99  1
dim(CHD2)  # 295  1

overap_test = function(CHD, pull_genes, N=22000){
    x = length(intersect(CHD, pull_genes))
    y = length(setdiff(pull_genes, CHD))
    z = length(setdiff(CHD, pull_genes))
    contingency_table = matrix(c(x, y, z, N-x-y-z), nrow=2, byrow=TRUE)
    res = fisher.test(contingency_table)
    return(list(p_value = res$p.value, OR = res$estimate, CI = res$conf.int))
}

par(mfrow=c(2,3))
#for(l in c('CM', 'CF', 'all')){
l = 'all'
  g_list <- list()
  for (i in seq_len(nrow(sig_map))) {
    name <- sig_map$db[i]
    sig_df <- subset(pull_df, Database == name )
    if(l != 'all') sig_df = sig_df[which(sig_df$linkeage == l),]
    pull_genes <- sig_df[,c('x','y')] %>% unlist %>% unique
    venn(list(PCGC = CHD1[,1], pull_genes = pull_genes))
    int = intersect(CHD1[,1], pull_genes)
    res = overap_test(CHD1[,1], pull_genes)
    mtext(paste0(name, '\n', toString(int), "\nOR=", round(res$OR, 1), 
       " p=", format(res$p_value, scientific = TRUE, digits = 1)), side = 3, line = 0)
  }
#}
#for(l in c('CM', 'CF', 'all')){
l = 'all'   
  g_list <- list()
  for (i in seq_len(nrow(sig_map))) {
    name <- sig_map$db[i]
    sig_df <- subset(pull_df, Database == name )
    if(l != 'all') sig_df = sig_df[which(sig_df$linkeage == l),]
    pull_genes <- sig_df[,c('x','y')] %>% unlist %>% unique
    venn(list(PCGC = CHD2[,1], pull_genes = pull_genes))
    int = intersect(CHD2[,1], pull_genes)
    res = overap_test(CHD2[,1], pull_genes)
    mtext(paste0(name, '\n', toString(int), "\nOR=", round(res$OR, 1), 
       " p=", format(res$p_value, scientific = TRUE, digits = 1)), side = 3, line = 0)
  }
#}

dev.copy2pdf(file = "pullgene_CHD_enrichment.pdf")