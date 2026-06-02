setwd("F:/projects/scRNA/results/cardiac_CTS_GRN/evaluation_TF_perturb/")
library(dplyr)
library(igraph)
library(gplots)
library(ggplot2)
library(openxlsx)
library(ggpubr)



#### TF perturb data , published GO-enrichment reuslts
data_path = 'D:/projects/DS/data/TF_perturb_CM/doc/'

df_CHD = read.xlsx(paste0(data_path, 'media-2.xlsx'), sheet = 11, colNames = TRUE, startRow = 2)
dim(df_CHD)  # [1] 17  5

head(df_CHD,3)
#   Program CHD_Gene_Overlap      P-value
# 1       5               14 1.643646e-09
# 2       1               12 1.362282e-07
# 3      10                9 5.161361e-05
#                                                                                                             Overlapping.Genes
# 1 ['GATA4', 'JAG1', 'GPC3', 'HAND2', 'FBN1', 'SMARCB1', 'GATA5', 'TBX5', 'ZFPM2', 'GATA6', 'MEIS2', 'ZEB2', 'SMAD6', 'PRKD1']
# 2                  ['CITED2', 'FOXP1', 'KDR', 'TAB2', 'KDM6A', 'ACTC1', 'MYH7', 'MYBPC3', 'MEIS2', 'MYH6', 'NKX2-5', 'TBX20']
# 3                                           ['GATA4', 'CITED2', 'SF3B4', 'JAG1', 'GATA5', 'GATA6', 'MED13L', 'RAF1', 'RAD21']
#    adj_p_value
# 1 3.780387e-08
# 2 1.566625e-06
# 3 2.967783e-04

## count  the n_overlap_genes per row
df_CHD$n_overlap = sapply(strsplit(df_CHD$Overlapping.Genes, ', '), length)
# calcuaoted the odds ratio given 142 genes from the CHDgene database with the 300 genes defining each program, and tghe 22000 human genes
df_CHD$OR <- sapply(df_CHD$n_overlap, function(k) {
  unname(fisher.test(matrix(c(k, 142 - k, 300 - k, 22000 - 300 - 142 + k), nrow = 2))$estimate)
})

df_GO = read.xlsx(paste0(data_path, 'media-2.xlsx'), sheet = 10, colNames = TRUE, startRow = 2)
dim(df_GO)  # [1] 5934  10
head(df_GO,3)
#   Gene.program                                           Term Overlap
# 1            1   Cardiac Muscle Cell Development (GO:0055013)   15/26
# 2            1                Myofibril Assembly (GO:0030239)   17/46
# 3            1 Actomyosin Structure Organization (GO:0031032)   20/77
#        P-value Adjusted.P-value Old.P-value Old.Adjusted.P-value Odds.Ratio
# 1 2.052729e-21     3.247418e-18           0                    0   94.20574
# 2 7.393806e-20     5.848501e-17           0                    0   40.74656
# 3 1.201208e-19     6.334372e-17           0                    0   24.61529
#   Combined.Score
# 1       4487.501
# 2       1794.929
# 3       1072.384
#                                                                                                                    Genes
# 1                            FHOD3;ADPRHL1;SORBS2;LRRC10;MYLK3;SLC8A1;TTN;CSRP3;HEY2;TCAP;ALPK3;ALPK2;NKX2-5;PDLIM5;MYH6
# 2                     FHOD3;TMOD1;ADPRHL1;ACTN2;MYLK3;LMOD2;TTN;CSRP3;OBSCN;SYNPO2L;TNNT2;TCAP;PGM5;FLNC;MYOZ2;MYH6;MYH7
# 3 FHOD3;TMOD1;MYH7B;ACTN2;MYLK3;LMOD2;TTN;CSRP3;FRMD5;OBSCN;SYNPO2L;TNNT2;EPB41L3;TCAP;PGM5;PHACTR1;FLNC;MYOZ2;MYH6;MYH7
df_GO$n_overlap <- as.numeric(sub("/.*", "", df_GO$Overlap))


#tmp = df_GO[which(df_GO$'Adjusted.P-value' < 1e-5 & df_GO$'Odds.Ratio' > 20 & df_GO$n_overlap >=8  ## if want to make story for both GP1, GP28, and GP5
tmp = df_GO[which(df_GO$'Adjusted.P-value' < 1e-7 & df_GO$'Odds.Ratio' > 20 & df_GO$n_overlap >10   # more strigent cutoff to focus on the GP1 and GP28
          # & df_GO$Gene.program %in% df_CHD$Program[which(df_CHD$OR>4 & df_CHD$adj_p_value < 0.005 & df_CHD$n_overlap >=8)]
           ), ]
dim(tmp) #[1] 179  11  // 299 10
unique(tmp$Gene.program)
#  1   2   4   7  11  14  19  20  25  28  29  30  32  38  47 194
# ## otherwise GP5 is included :  [1]   1   2   4   5   7  11  13  14  18  19  20  21  22  25  28  29  30  32  38  43  47 113 194 227 236
unique(tmp$Gene.program) %>% length  # 16  // 25

df_CHD$Program[which(df_CHD$OR>4 & df_CHD$adj_p_value < 0.01 & df_CHD$n_overlap > 5)]
#[1]  5  1 10 28

selected = subset(tmp, Gene.program %in% df_CHD$Program[which(df_CHD$OR>4 & df_CHD$adj_p_value < 0.01 & df_CHD$n_overlap > 5)])
dim(selected) #[1] 10  11 // 18  11
# the top 3 enriched (lowest adjusted p-value) terms for each unique Gene.program
selected = selected[order(selected$Gene.program, selected$'Adjusted.P-value'), ]
top2_selected = do.call(rbind, lapply(split(selected, selected$Gene.program), function(df) head(df, 3)))
top2_selected
#     Gene.program                                           Term Overlap      P-value Adjusted.P-value
# 1.1            1   Cardiac Muscle Cell Development (GO:0055013)   15/26 2.052729e-21     3.247418e-18
# 1.2            1                Myofibril Assembly (GO:0030239)   17/46 7.393806e-20     5.848501e-17
# 1.3            1 Actomyosin Structure Organization (GO:0031032)   20/77 1.201208e-19     6.334372e-17
# 28            28      Collagen Fibril Organization (GO:0030199)   11/42 2.039477e-11     1.628522e-08
#     Old.P-value Old.Adjusted.P-value Odds.Ratio Combined.Score
# 1.1           0                    0   94.20574      4487.5015
# 1.2           0                    0   40.74656      1794.9290
# 1.3           0                    0   24.61529      1072.3844
# 28            0                    0   24.14991       594.4678
#                                                                                                                      Genes
# 1.1                            FHOD3;ADPRHL1;SORBS2;LRRC10;MYLK3;SLC8A1;TTN;CSRP3;HEY2;TCAP;ALPK3;ALPK2;NKX2-5;PDLIM5;MYH6
# 1.2                     FHOD3;TMOD1;ADPRHL1;ACTN2;MYLK3;LMOD2;TTN;CSRP3;OBSCN;SYNPO2L;TNNT2;TCAP;PGM5;FLNC;MYOZ2;MYH6;MYH7
# 1.3 FHOD3;TMOD1;MYH7B;ACTN2;MYLK3;LMOD2;TTN;CSRP3;FRMD5;OBSCN;SYNPO2L;TNNT2;EPB41L3;TCAP;PGM5;PHACTR1;FLNC;MYOZ2;MYH6;MYH7
# 28                                                  COL1A1;COL3A1;COL1A2;COL5A1;LOX;LUM;COL5A2;SERPINH1;LOXL1;TGFBR1;LOXL2
#     n_overlap
# 1.1        15
# 1.2        17
# 1.3        20
# 28         11
### 687             5      Aortic Valve Morphogenesis (GO:0003180)    9/37 2.934784e-09     2.037678e-06
### 687            0                    0   21.72901       426.9019
#### 687                                          JAG1;HEY1;CDH11;SNAI1;GATA5;TWIST1;GATA4;SOX9;SMAD6
### 687          9

## select TFs with adjusted p-value < 0.05
# venn(list(CHD = df_CHD$Program[which(df_CHD$OR>4 & df_CHD$adj_p_value < 0.005)], 
#             GO = df_GO$Gene.program[which(df_GO$'Adjusted.P-value' < 0.005 & df_GO$'Odds.Ratio' > 4)]))

##################### plot ###############################
library(ggplot2)
library(patchwork)

## panel A: counts
df_count <- data.frame(
  set = factor(c("Neither", "GO only", "GO + CHD"),
               levels = c("Neither", "GO only", "GO + CHD")),
  n = c(234, 14, 2)
)
df_count$label <- paste0(df_count$n, " (", round(100 * df_count$n / 250, 1), "%)")

p1 <- ggplot(df_count, aes(x = "Gene programs", y = n, fill = set)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), size = 4) +
  coord_flip() +
  labs(
    title = "250 gene programs screened",
    subtitle = "GO-enriched: OR > 20, FDR < 1e-7, overlap >10\nCHD-enriched: OR > 4, FDR < 0.01, overlap > 5",
    x = NULL, y = "Number of programs"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

## panel B: the 3 overlapping programs and their leading GO terms
df_hit <- data.frame(
  program = c("GP1", "GP1", "GP28"),
  term = c(
    "Cardiac Muscle Cell Development",
    "Myofibril Assembly",
    "Collagen Fibril Organization"
  ),
  GO_ID = c("GO:0055013", "GO:0030239", "GO:0030199")
)

p2 <- ggplot(df_hit, aes(x = 1, y = GO_ID)) +
  geom_point(size = 3) +
  geom_segment(aes(x = 1.03, xend = 1.14, yend =GO_ID), linewidth = 0.4) +
  geom_text(aes(x = 1.16, label = paste0(term, "\n (", GO_ID, ")")),
            hjust = 0, size = 4) +
  xlim(0.98, 2.4) +
  labs(
    title = "Programs enriched in both GO and CHD",
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 12) +
  theme(plot.title.position = "plot")

p1 + p2 + plot_layout(widths = c(1.8, 1.1))
dev.copy2pdf(file='GO_CHD_enriched_programs_selection.pdf')

