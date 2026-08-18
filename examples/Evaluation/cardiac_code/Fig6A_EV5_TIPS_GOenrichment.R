# source('F:/projects/scRNA/source/cardiac_CTS_GRN/IbarraSoria2018_E8.25_v9/24.3_overlap_with_publications.R')
setwd("F:/projects/scRNA/results/cardiac_CTS_GRN/evaluation_reproducibility/enrichment/")

library(dplyr)

#library(openxlsx)

library(forcats)
library(tidytext)
 
library(tidyr)
library(ggplot2)
library(ggrepel)


## ================================================================================
################ enrichment with published gene sets ################################
## ================================================================================ 

CTS_res = read.table( file= 'GO_BP_CTS_padj0.05_OR2_minOverlap3.tsv', sep='\t', header=T)
dim(CTS_res) # [1] 80  11
table( CTS_res$GO_BP)
#  CTS_Elorbany.CP CTS_Elorbany.CP.1        CTS_Ibarra        CTS_Pijuan 
#                3                13                17                47 

CTS_res = CTS_res[which(CTS_res$GO_BP  != 'CTS_Elorbany.CP'),]
CTS_res$GO_BP[which(CTS_res$GO_BP == 'CTS_Elorbany.CP.1')] = 'CTS_Elorbany'


TIPS_res = read.table( file= 'GO_BP_TIPS_padj0.05_OR2_minOverlap2.tsv', sep='\t', header=T)
dim(TIPS_res) # [1] 93  11
table( TIPS_res$GO_BP)
  #  TIPS_Elorbany_IID TIPS_Elorbany_STRING      TIPS_Ibarra_IID   TIPS_Ibarra_STRING 
  #                 18                   18                   22                   22 
  # TIPS_Pijuan_STRING 
  #                 13   
TIPS_res = TIPS_res[which(grepl('_STRING', TIPS_res$GO_BP)),]



### combine two resutls
CTS_res$category = unlist(lapply(strsplit(CTS_res$GO_BP, '_'), '[', 1))
CTS_res$database = unlist(lapply(strsplit(CTS_res$GO_BP, '_'), '[', 2))
TIPS_res$category = unlist(lapply(strsplit(TIPS_res$GO_BP, '_'), '[', 1))
TIPS_res$database = unlist(lapply(strsplit(TIPS_res$GO_BP, '_'), '[', 2))

enrich_df <-   rbind(CTS_res, TIPS_res)
dim(enrich_df) # [1] 130 13
 
  
## significance cutoff
padj_cutoff <- 0.05
dataset_order <- c("Ibarra", "Pijuan","Elorbany")
group_order <- c("CTS", "TIPS")


# make CTS and TIPS per each dataset togetehr
col_levels <- expand.grid(
  dataset = dataset_order,
  group = group_order
) %>%
  mutate(col_id = paste(dataset, group, sep = "\n")) %>%
  pull(col_id)
 

heat_df <- enrich_df %>%
  mutate(
    dataset = factor(database, levels = dataset_order),
    group = factor(category, levels = group_order),
    value = ifelse(adj_p < padj_cutoff, -log10(adj_p), NA_real_)
  ) %>%
  complete(GO_BP, dataset, group, fill = list(value = NA_real_)) %>%
  mutate(col_id = factor(
    paste(dataset, group, sep = "\n"),
    levels = col_levels    #as.vector(outer(dataset_order, group_order, paste, sep = "\n"))
  ))

## reorder GS rows within each facet by source group
gs_levels <- heat_df %>%
  distinct(GO_BP) %>%
  mutate(source_group = case_when(
    grepl("^CTS_Ibarra", GO_BP) ~ 1,
    grepl("^CTS_Pijuan", GO_BP) ~ 2,
    grepl("^CTS_Elorbany", GO_BP) ~ 3,
    grepl("^TIPS_Ibarra", GO_BP) ~ 4,
    grepl("^TIPS_Pijuan", GO_BP) ~ 5,
    grepl("^TIPS_Elorbany", GO_BP) ~ 6,
  )) %>%
  arrange(source_group, GO_BP) %>%
  pull(GO_BP)

heat_df <- heat_df %>% filter(!is.na(GS)) %>%
  mutate(GO_BP = factor(GO_BP, levels = rev(gs_levels)))
dim(heat_df) # 130  7 

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

  dev.copy2pdf(file= 'CTS_vsTIPS_heatmap_GO_BP.pdf')

## narrowdown to those have 3 ore more calls
x = table(heat_df$GS)
heat_df2 = heat_df[which(heat_df$GS %in% names(x[which(x >= 3)])),]
dim(heat_df2) # [1] 25 13

write.table(heat_df2, 'GO_BP_3ormorecall_Fig5a.txt', sep='\t', row.names=F)  #!!!!!!

ggplot(heat_df2, aes(x = col_id, y = GS, fill = value)) +
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
  ) +
  ggtitle("CTS vs TIPS enrichment with GO_BP (3 or more calls)")

  dev.copy2pdf(file= 'CTS_vsTIPS_heatmap_GO_BP_overlapped_calls.pdf')  ## Fig 5A


  ################### scatter plot ######### 
enrich_df$GO_BP <- gsub("_STRING$", "", enrich_df$GO_BP)

padj_cutoff <- 0.05

concord_df <- enrich_df %>%
  dplyr::mutate(
    dataset = database,
    group = category,
    neglog10_padj = -log10(adj_p)
  ) %>%
  dplyr::select(dataset, group, GS, neglog10_padj, adj_p) %>%
  dplyr::group_by(dataset, GS, group) %>%
  dplyr::summarise(
    neglog10_padj = max(neglog10_padj, na.rm = TRUE),
    adj_p = min(adj_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = c(neglog10_padj, adj_p)
  ) %>%
  dplyr::mutate(
    neglog10_padj_CTS  = tidyr::replace_na(neglog10_padj_CTS, 0),
    neglog10_padj_TIPS = tidyr::replace_na(neglog10_padj_TIPS, 0),
    adj_p_CTS  = tidyr::replace_na(adj_p_CTS, 1),
    adj_p_TIPS = tidyr::replace_na(adj_p_TIPS, 1),
    CTS_sig  = adj_p_CTS < padj_cutoff,
    TIPS_sig = adj_p_TIPS < padj_cutoff,
    status = dplyr::case_when(
      CTS_sig & TIPS_sig ~ "Both",
      CTS_sig & !TIPS_sig ~ "CTS only",
      !CTS_sig & TIPS_sig ~ "TIPS only",
      TRUE ~ "NS"
    )
  )
  table(concord_df$status)
    # Both  CTS only TIPS only 
    #    22        55        31 
label_gs <- unique(concord_df$GS[which(concord_df$status %in% c("Both", "TIPS only"))])
# c(
#   "pharyngeal system development",
#   "outflow tract septum morphogenesis",
#   "epithelial to mesenchymal transition",
#   "mesenchymal cell differentiationet"
 
# )

concord_df$dataset <- factor(concord_df$dataset, levels = dataset_order)

## remove the GOprefix  from the GS names
concord_df$GS <- sub("^GO:[0-9]+_", "", concord_df$GS)

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

  pdf(file= 'CTS_vsTIPS_scatterplot_GO_BP.pdf', width=10, height=6)  ## Fig EV5!!
  print(g)
  dev.off()



