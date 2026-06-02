setwd("F:/projects/scRNA/results/cardiac_CTS_GRN/evaluation_reproducibility/PijuanE825_AlternativeBranch/")

## =========================
## 1) Packages
## =========================
library("ggplot2")
library("pheatmap")
library("openxlsx")
library("dplyr")

## =========================
## 2) File locations
## =========================
##  HiG_method: "1v1.25"
##  PPI: "STRING"

pull_path = 'F:/projects/TIPS/doc/'
pull_df = read.xlsx(paste0(pull_path, 'STable_final_2026.xlsx'), sheet=15)  # for PijuanSampleE825 alternative branch results only 
dim(pull_df)  # 53  17
head(pull_df, 3)
#       Database CTS_ID counter_cluster linkeage     x      y       w_CP  w_decendent       delta
# 1 Pijuan_Sala      8              16       CF  DAB2  HMGA2  0.0000000 -0.005129130 -0.00512913
# 2 Pijuan_Sala      8              16       CM HMGA2 MPPED2 -0.1528846  0.205327212  0.35821185
# 3 Pijuan_Sala      8              16       CM  FGF8  HMGA2 -0.1243642  0.002014388  0.12637861
#    abs_delta direction  status rank                       TF_highConf                motif  NES
# 1 0.00512913  decrease  gained    1 HMGA1; HMGA2 (directAnnotation).  transfac_pro__M07320 3.47
# 2 0.35821185  increase changed    1 HMGA1; HMGA2 (directAnnotation).  transfac_pro__M07320 3.47
# 3 0.12637861  increase changed    2 HMGA1; HMGA2 (directAnnotation).  transfac_pro__M07320 3.47

colnames(pull_df)[grepl('TF_highConf.', colnames(pull_df))] = 'TF_highConf'
unique(pull_df$Database)   # [1] "Pijuan_Sala"
table(pull_df$counter_cluster)
#  3 16 18 19 
# 12 14 13 14 
## =========================
## 4) gele lists & Helper functions
## =========================

df = pull_df

CM_bias_lists = list(
  "C8vsC16" = df[which(df$counter_cluster == "16" &  df$linkeage == "CM"), c('x','y')] %>% unlist() %>% unique(),
  "C8vsC18" = df[which(df$counter_cluster == "18" &  df$linkeage == "CM"), c('x','y')] %>% unlist() %>% unique(),
  "C8vsC19" = df[which(df$counter_cluster == "19" &  df$linkeage == "CM"), c('x','y')] %>% unlist() %>% unique(),
  "C8vsC3" = df[which(df$counter_cluster == "3" &  df$linkeage == "CM"), c('x','y')] %>% unlist() %>% unique()
)
lengths(CM_bias_lists)
# C8vsC16 C8vsC18 C8vsC19  C8vsC3 
 #     7       7       7       7 


CF_bias_lists = list(
  "C8vsC16" = df[which(df$counter_cluster == "16" &  df$linkeage == "CF"), c('x','y')] %>% unlist() %>% unique(),
  "C8vsC18" = df[which(df$counter_cluster == "18" &  df$linkeage == "CF"), c('x','y')] %>% unlist() %>% unique(),
  "C8vsC19" = df[which(df$counter_cluster == "19" &  df$linkeage == "CF"), c('x','y')] %>% unlist() %>% unique(),
  "C8vsC3" = df[which(df$counter_cluster == "3" &  df$linkeage == "CF"), c('x','y')] %>% unlist() %>% unique()
)
lengths(CF_bias_lists)
# C8vsC16 C8vsC18 C8vsC19  C8vsC3 
#        4       5       6       4 

dual_bias_lists = list()
for(i in 1:length(CM_bias_lists)) dual_bias_lists[[i]] = intersect(CM_bias_lists[[i]], CF_bias_lists[[i]])
lengths(dual_bias_lists)
# [1]  3 2 3 2
names(dual_bias_lists) = names(CM_bias_lists)



#########  extract leangeage-leaning genes
CM_lean_list = CF_lean_list = list()
for(i in 1:length(CM_bias_lists)) {
  CM_lean_list[[i]] = setdiff(CM_bias_lists[[i]], dual_bias_lists[[i]])
  CF_lean_list[[i]] = setdiff(CF_bias_lists[[i]], dual_bias_lists[[i]])
}
names(CM_lean_list) = names(CF_lean_list) = names(CM_bias_lists) 


print(CM_lean_list)
# $C8vsC16
# "FGF8"   "MPPED2" "BMP7"   "RARB"   "NPM3"  

# $C8vsC18
# [1] "FGF8"   "MPPED2" "BMP7"   "RARB"   "NPM3"  

# $C8vsC19
# [1] "FGF8"   "BMP7"   "MPPED2" "RARB"  

# $C8vsC3
# [1] "FGF8"   "MPPED2" "BMP7"   "RARB"   "NPM3"  


print(CF_lean_list)
# $C8vsC16
# [1] "DAB2" "TPM2" 

# $C8vsC18
# [1] "DAB2" "TPM2" "P3H2"

# $C8vsC19
# [1] "DAB2" "TPM2" "P3H2"

# $C8vsC3
# [1] "DAB2" "TPM2"



jaccard <- function(x, y) {
  x <- unique(x)
  y <- unique(y)
  u <- union(x, y)
  if (length(x) == 0 && length(y) == 0) return(NA_real_)
  if (length(u) == 0) return(NA_real_)
  length(intersect(x, y)) / length(u)
}

## =========================
## 5) Read leange_lean gene lists
## =========================
names(CM_lean_list) = paste0('CM_', names(CM_lean_list))
names(CF_lean_list) = paste0('CF_', names(CF_lean_list))
gene_lists = c(CM_lean_list, CF_lean_list)


## Quick check: number of genes in each list
gene_counts <- data.frame(
  List = names(gene_lists),
  N_Genes = sapply(gene_lists, length),
  row.names = NULL
)
print(gene_counts)
write.csv(gene_counts,  "gene_list_sizes_lineage_lean.csv", row.names = FALSE)

## =========================
## 5.2 ) Read leange_biased gene lists
## =========================
names(CM_bias_lists) = paste0('CM_', names(CM_bias_lists))
names(CF_bias_lists) = paste0('CF_', names(CF_bias_lists))

gene_biased_lists = c(CM_bias_lists, CF_bias_lists)
names(gene_biased_lists)

jac_mat_biased <- outer(
  names(gene_biased_lists), names(gene_biased_lists),
  Vectorize(function(a, b) jaccard(gene_biased_lists[[a]], gene_biased_lists[[b]]))
)

rownames(jac_mat_biased) <- names(gene_biased_lists)
colnames(jac_mat_biased) <- names(gene_biased_lists)

print(round(jac_mat_biased, 3))

write.csv(round(jac_mat_biased, 3), "pairwise_jaccard_lineage_bias_matrix.csv", row.names = TRUE)


## =========================
## 5.3 ) dual_biased gene list overlaps
## =========================
names(dual_bias_lists)
# [1] "Ibarra_STRING"   "Ibarra_IID"      "Pijuan_STRING"   "Pijuan_IID"     
# [5] "Elorbany_STRING" "Elorbany_IID" 

jac_mat_dual <- matrix(
  NA_real_,
  nrow = length(names(dual_bias_lists)),
  ncol = length(names(dual_bias_lists)),
  dimnames = list(names(dual_bias_lists), names(dual_bias_lists))
)

## Fill matrix
for (i in seq_along(names(dual_bias_lists))) {
  for (j in seq_along(names(dual_bias_lists))) {
    jac_mat_dual[i, j] <- jaccard(
      dual_bias_lists[[names(dual_bias_lists)[i]]],
      dual_bias_lists[[names(dual_bias_lists)[j]]]
    )
  }
}
print(round(jac_mat_dual, 3))

write.csv(round(jac_mat_dual, 3), "pairwise_jaccard_dual_biased_matrix.csv", row.names = TRUE)


###########################################
## summary plot 
## =========================
library(ggplot2)
## =========================
## Summary Jaccard plot:
## dual-pull vs lineage-lean
## Uses existing jac_mat and jac_mat_dual
jac_mat = read.csv("pairwise_jaccard_lineage_lean_matrix.csv", row.names = 1)
jac_mat_dual = read.csv("pairwise_jaccard_dual_biased_matrix.csv", row.names = 1)

## Helper: convert Jaccard matrix to long format
mat_to_pairs <- function(mat, category) {   
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)  
  colnames(df) <- c("list1", "list2", "jaccard")
  
  ## remove self-comparisons
  df <- subset(df, list1 != list2)
  
  ## keep upper triangle only
  ord1 <- match(df$list1, rownames(mat))
  ord2 <- match(df$list2, colnames(mat))
  df <- df[ord1 < ord2, ]
  
  df$Category <- category
  df
}

## Parse names for lineage-lean names like:
## CM_Ibarra_STRING, CF_Pijuan_IID
parse_lean_name <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  data.frame(
    lineage = parts[1],
    dataset = parts[2],
    network = parts[3],
    stringsAsFactors = FALSE
  )
}

## Parse names for dual names like:
## Ibarra_STRING, Pijuan_IID
parse_dual_name <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  data.frame(
    lineage = "Dual-pull",
    dataset = parts[1],
    network = parts[2],
    stringsAsFactors = FALSE
  )
}

## =========================
## 1) Lineage-lean summary from jac_mat
## =========================
jac_mat <- as.matrix(jac_mat)
mode(jac_mat) <- "numeric"
lean_df <- mat_to_pairs(jac_mat, "Lineage-lean")

info1 <- do.call(rbind, lapply(lean_df$list1, parse_lean_name))
info2 <- do.call(rbind, lapply(lean_df$list2, parse_lean_name))

lean_df$lineage1 <- info1$lineage
lean_df$dataset1 <- info1$dataset
lean_df$network1 <- info1$network

lean_df$lineage2 <- info2$lineage
lean_df$dataset2 <- info2$dataset
lean_df$network2 <- info2$network

## keep only CM-vs-CM or CF-vs-CF comparisons
lean_df <- subset(lean_df, lineage1 == lineage2)

lean_df$Group <- paste0(lean_df$lineage1, "-lean")

## =========================
## 2) Dual-pull summary from jac_mat_dual
## =========================
jac_mat_dual <- as.matrix(jac_mat_dual)
mode(jac_mat_dual) <- "numeric"
dual_df <- mat_to_pairs(jac_mat_dual, "Dual-pull")

info1 <- do.call(rbind, lapply(dual_df$list1, parse_dual_name))
info2 <- do.call(rbind, lapply(dual_df$list2, parse_dual_name))

dual_df$lineage1 <- info1$lineage
dual_df$dataset1 <- info1$dataset
dual_df$network1 <- info1$network

dual_df$lineage2 <- info2$lineage
dual_df$dataset2 <- info2$dataset
dual_df$network2 <- info2$network

dual_df$Group <- "Dual-pull"

## =========================
## 3) Combine and define comparison type
## =========================

summary_df <- rbind(
  dual_df[, c("list1", "list2", "jaccard", "Group", "dataset1", "dataset2", "network1", "network2")],
  lean_df[, c("list1", "list2", "jaccard", "Group", "dataset1", "dataset2", "network1", "network2")]
)

summary_df$Comparison <- ifelse(
  summary_df$dataset1 == summary_df$dataset2,
  "Within-dataset",
  "Cross-dataset"
)

summary_df$Group <- factor(
  summary_df$Group,
  levels = c("Dual-pull", "CM-lean", "CF-lean")
)

summary_df$Comparison <- factor(
  summary_df$Comparison,
  levels = c("Within-dataset", "Cross-dataset")
)

summary_df$Dataset_pair <- paste(summary_df$dataset1, "vs", summary_df$dataset2)
summary_df$Dataset_pair <- factor(
  summary_df$Dataset_pair,
  levels = unique(summary_df$Dataset_pair)
)

print(summary_df)

write.csv(summary_df, "summary_jaccard_dual_pull_lineage_lean.csv", row.names = FALSE)

## =========================
## 4) Mean/median summary
## =========================

mean_df <- aggregate(
  jaccard ~ Group + Comparison,
  data = summary_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)

median_df <- aggregate(
  jaccard ~ Group + Comparison,
  data = summary_df,
  FUN = function(x) median(x, na.rm = TRUE)
)

print(mean_df)
print(median_df)

write.csv(mean_df, "summary_jaccard_mean_dual_pull_lineage_lean.csv", row.names = FALSE)
write.csv(median_df, "summary_jaccard_median_dual_pull_lineage_lean.csv", row.names = FALSE)

## =========================
## 5) Main summary plot
## =========================
p_summary <- ggplot(summary_df, aes(x = Group, y = jaccard)) +
  
  ## CHANGED: color points by dataset pair
  geom_jitter(
    aes(color = Dataset_pair),
    width = 0.12,
    height = 0,
    size = 2.8,
    alpha = 0.85
  ) +
  
  ## CHANGED: mean shown as black diamond
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 5,
    color = "black"
  ) +
  
  facet_wrap(~ Comparison, nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  
  ## CHANGED: add color legend title
  labs(
    title = "Summary of Jaccard overlap for dual-pull and lineage-lean",
    x = NULL,
    y = "Jaccard score",
    color = "Dataset pair"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "grey95", colour = "grey50"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_summary)

ggsave(
  filename = "summary_jaccard_dual_pull_lineage_lean.pdf",  ## Fig 5E !!!
  plot = p_summary,
  width = 7,
  height = 4,
  dpi = 300
)
