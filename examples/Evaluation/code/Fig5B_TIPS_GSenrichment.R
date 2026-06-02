# source('F:/projects/scRNA/source/cardiac_CTS_GRN/IbarraSoria2018_E8.25_v9/24.3_overlap_with_publications.R')
setwd("F:/projects/scRNA/results/cardiac_CTS_GRN/evaluation_reproducibility/enrichment/")


library(forcats)
library(tidytext)
library(tidyr)
library(ggplot2)
library(ggrepel)


## ================================================================================
################ enrichment test with published gene sets ################################
## ================================================================================ 
## refer to 
# 1) F:\projects\scRNA\source\cardiac_CTS_GRN\IbarraSoria2018_E8.25_v9\24.4.1_overlap_with_publications.R
# 2) F:\projects\scRNA\source\cardiac_CTS_GRN\GSE87038_Pijuan2019_v9\C8_vs_C16\24.3_overlap_with_publications.R
# 3) F:\projects\scRNA\source\cardiac_CTS_GRN\GSE175634_iPSC_CM_weighted_v9\24.4.1_overlap_with_publications.R


## ================================================================================
################ plot ################################
## ================================================================================ 

CTS_res = read.table( file= '/data/Evaluation/GS_CTS_padj0.05_OR2_minOverlap3.tsv', sep='\t', header=T)
dim(CTS_res) # [1] 71  11
## simplify the GP to be compared
# GP5 is the most strongly enriched for CHD genes, 
# GP1 epresent cardiac differentiation and heart construction, 
# GP28 represents CHD-related extracellular matrix.
CTS_res  = CTS_res[-which(CTS_res$GS %in% c('Takeuchi_gp10', 'Takeuchi_gp2')),]
# leave the CTS vs TD overlap for unique topic
CTS_res = CTS_res[-which(grepl('Elorbany_TD', CTS_res$GS)),]
table( CTS_res$scRNA_signature)
#   CTS_Elorbany.CP CTS_Elorbany.CP.1        CTS_Ibarra        CTS_Pijuan 
#             16                16                15                16 
CTS_res = CTS_res[which(CTS_res$scRNA_signature  != 'CTS_Elorbany.CP'),]
CTS_res$scRNA_signature[which(CTS_res$scRNA_signature == 'CTS_Elorbany.CP.1')] = 'CTS_Elorbany'


TIPS_res = read.table( file= '/data/Evaluation/GS_TIPS_padj0.05_OR2_minOverlap2.tsv', sep='\t', header=T)
dim(TIPS_res) # [1] 62  11
TIPS_res  = TIPS_res[-which(TIPS_res$GS %in% c('Takeuchi_gp10', 'Takeuchi_gp2')),]
# leave the CTS vs TD overlap for unique topic
TIPS_res = TIPS_res[-which(grepl('Elorbany_TD', TIPS_res$GS)),]
table( TIPS_res$scRNA_signature)
#    TIPS_Elorbany_IID TIPS_Elorbany_STRING      TIPS_Ibarra_IID   TIPS_Ibarra_STRING 
#                   13                   14                   11                   11 
#      TIPS_Pijuan_IID   TIPS_Pijuan_STRING 
#                    3                    4 
TIPS_res = TIPS_res[which(grepl('_STRING', TIPS_res$scRNA_signature )),]


### combine two resutls
CTS_res$category = unlist(lapply(strsplit(CTS_res$scRNA_signature, '_'), '[', 1))
CTS_res$database = unlist(lapply(strsplit(CTS_res$scRNA_signature, '_'), '[', 2))
TIPS_res$category = unlist(lapply(strsplit(TIPS_res$scRNA_signature, '_'), '[', 1))
TIPS_res$database = unlist(lapply(strsplit(TIPS_res$scRNA_signature, '_'), '[', 2))

enrich_df <-  rbind(CTS_res, TIPS_res)
dim(enrich_df) # [1] 76 13
## simplify the GS names
grep('_gene', enrich_df$GS, value=T) %>% unique
#  [1] "Maven2023_gene_ISL1_targets_CP_E" "Maven2023_gene_ISL1_targets_CP_L"
#  [3] "Maven2023_gene_ISL_dn_E"          "Maven2023_gene_ISL_dn_L"         
#  [5] "Maven2023_gene_ISL_dn_T"          "Maven2023_gene_ISL1_WT_d6CP"     
#  [7] "Maven2023_gene_ISL1_targets_CP_T" "Maven2023_gene_ISL_up_T"         
#  [9] "Maven2023_gene_ISL_up_L"          "Maven2023_gene_ISL_up_E"     
enrich_df$GS <- gsub("_gene", "", enrich_df$GS)

# remove Maven2023_gene_ISL_upxxx and Maven2023_gene_ISL_dnxxxx to be simplify
enrich_df = enrich_df[-which(grepl('_ISL_up', enrich_df$GS) | grepl('_ISL_dn', enrich_df$GS)),]
dim(enrich_df)  # 49  13

write.table(enrich_df, 'PublishedGS_TIPS_CTS_Fig5b.txt', sep='\t', row.names=F)  #!!!!!!


# make CTS and TIPS per each dataset togetehr
col_levels <- expand.grid(
  dataset = dataset_order,
  group = group_order
) %>%
  mutate(col_id = paste(dataset, group, sep = "\n")) %>%
  pull(col_id)
 
## significance cutoff
padj_cutoff <- 0.05
dataset_order <- c("Ibarra", "Pijuan","Elorbany")
group_order <- c("CTS", "TIPS")

heat_df <- enrich_df %>%
  mutate(
    dataset = factor(database, levels = dataset_order),
    group = factor(category, levels = group_order),
    value = ifelse(adj_p < padj_cutoff, -log10(adj_p), NA_real_)
  ) %>%
  complete(GS, dataset, group, fill = list(value = NA_real_)) %>%
  mutate(col_id = factor(
    paste(dataset, group, sep = "\n"),
    levels = col_levels    #as.vector(outer(dataset_order, group_order, paste, sep = "\n"))
  ))

## reorder GS rows within each facet by source group
gs_levels <- heat_df %>%
  distinct(GS) %>%
  mutate(source_group = case_when(
    grepl("^Maven2023", GS) ~ 1,
    grepl("^Elorbany", GS) ~ 2,
    grepl("^Xie", GS) ~ 3,
    grepl("^Takeuchi", GS) ~ 4,
    grepl("^PCGC", GS) ~ 5,
    TRUE ~ 6
  )) %>%
  arrange(source_group, GS) %>%
  pull(GS)

heat_df <- heat_df %>%
  mutate(GS = factor(GS, levels = rev(gs_levels)))

ggplot(heat_df, aes(x = col_id, y = GS, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white", high = "darkblue", na.value = "lightgrey",
    name = expression(-log[10]("padj"))
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.ticks = element_blank()
  )

  dev.copy2pdf(file= 'CTS_vsTIPS_heatmap_published_GS.pdf')


  ################### scatter plot ######### 


padj_cutoff <- 0.05

concord_df <- enrich_df %>%
  mutate(
    dataset = database,
    group = category,
    neglog10_padj = -log10(adj_p)
  ) %>%
  select(dataset, group, GS, neglog10_padj, adj_p) %>%
  pivot_wider(
    names_from = group,
    values_from = c(neglog10_padj, adj_p)
  ) %>%
  mutate(
    neglog10_padj_CTS  = replace_na(neglog10_padj_CTS, 0),
    neglog10_padj_TIPS = replace_na(neglog10_padj_TIPS, 0),
    adj_p_CTS  = replace_na(adj_p_CTS, 1),
    adj_p_TIPS = replace_na(adj_p_TIPS, 1),
    CTS_sig  = adj_p_CTS  < padj_cutoff,
    TIPS_sig = adj_p_TIPS < padj_cutoff,
    status = case_when(
      CTS_sig & TIPS_sig ~ "Both",
      CTS_sig & !TIPS_sig ~ "CTS only",
      !CTS_sig & TIPS_sig ~ "TIPS only",
      TRUE ~ "NS"
    )
  )

label_gs <- unique(concord_df$GS)
# c(
#   "Maven2023_ISL1_WT_d6CP",
#   "Maven2023_ISL1_targets_CP_T",
#   "Elorbany_C5_RegulonMEF2A_target",
#   "Elorbany_C1_RegulonWT1_target",
#   "PCGC_AllCurated",
#   "Takeuchi_gp5",
#   "Xie_JCF specific",
#   "Xie_SHF specific"
# )

concord_df$dataset <- factor(concord_df$dataset, levels = dataset_order)

g = ggplot(concord_df, aes(x = neglog10_padj_CTS, y = neglog10_padj_TIPS)) +
  geom_point(aes(color = status), size = 2.8, alpha = 0.9) +
  geom_label_repel(
    data = subset(concord_df, GS %in% label_gs),
    aes(label = GS),
    size = 2.8,
    label.size = 0.2,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    segment.alpha = 0.8,
    force = 1.5,
    nudge_x = 0.2,
    nudge_y = 0.2,
    max.overlaps = Inf
  ) +
  facet_wrap(~ dataset, scales = "free") +
  expand_limits(x = 0, y = 0) +
  scale_color_manual(
    values = c(
      "Both" = "blue",
      "TIPS only" = "darkred",
      "CTS only" = "darkorange3",
      "NS" = "grey70"
    ),
    drop = TRUE
  ) +
  labs(
    x = expression(CTS ~ -log[10]("padj")),
    y = expression(TIPS ~ -log[10]("padj")),
    color = NULL
  ) +
  theme_classic(base_size = 12)

  pdf(file= 'CTS_vsTIPS_scatterplot_published_GS.pdf', width=10, height=6)
  print(g)
  dev.off()
