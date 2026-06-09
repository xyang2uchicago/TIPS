rebuild_mat <- TRUE
source("/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/GSE87038_IID/code/24.0_acat_load_input_clean.R")

CTS_ID # [1] "8"
length(CHD) # 295
celltype_col # [1] "label"

updir <- "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/GSE87038_IID/results/PijuanSara2019_vsCF18/GSE181346_heart_scATAC"
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

library(igraph)
library(dplyr)
library(igraph)

file <- paste0(db_specifc_input_path, "GSE87038_IID_graph_perState_simplified_combinedweighted.rds")
graph_list <- readRDS(file)
names(graph_list)
#  [1] "HiG_1"     "HiG_2"     "HiG_3"     "HiG_4"     "HiG_5"     "HiG_6"    
#  [7] "HiG_9"     "HiG_10"    "HiG_12"    "HiG_14"    "HiG_17"    "HiG_18"   
# [13] "HiG_19"    "HiG_7"     "HiG_11"    "HiG_15"    "HiG_16"    "HiG_13"   
# [19] "HiG_8"     "HiGCTS_15" "CTS_7"     "CTS_11"    "CTS_15"    "CTS_16"   
# [25] "CTS_16.1"  "CTS_8"   


## extract the predicted TF regulators
(files <- list.files(pattern = "final.rds"))
# [1] "PPI_graph_HMGA2_GRN_prediction_CTS_8_CF_final.rds"
# [2] "PPI_graph_HMGA2_GRN_prediction_CTS_8_CM_final.rds"
# [3] "PPI_graph_KLF6_GRN_prediction_CTS_8_CM_final.rds" 
# [4] "PPI_graph_RARB_GRN_prediction_CTS_8_CM_final.rds" 
key_TFs <- lapply(files, function(x) strsplit(x, "_")[[1]][3] %>% unlist()) %>%
    unlist() %>%
    unique()
key_TFs # [1] "HMGA2" "KLF6"  "RARB"


####### generateing merged graph final table of the predicted subnetwork #################
(file <- list.files(pattern = "final_table.tsv"))
#  "PPI_graph_GRN_prediction_CTS_8_dualpull_final_table.tsv"


final_table <- read.table(file = file, header = TRUE, sep = "\t")
head(final_table)
#   linkeage   from     to         w1          w2      delta abs_delta direction
# 1       CF   DAB2  HMGA2  0.0000000 -0.12320331 -0.1232033 0.1232033  decrease
# 2       CM  HMGA2 MPPED2 -0.1528846  0.14451668  0.2974013 0.2974013  increase
# 3       CM   FGF8  HMGA2 -0.1243642  0.01338294  0.1377472 0.1377472  increase
# 4       CM   KLF6 MPPED2  0.0000000  0.12347294  0.1234729 0.1234729  increase
# 5       CM MPPED2   RARB  0.0000000  0.11809334  0.1180933 0.1180933  increase
#    status rank
# 1  gained    1
# 2 changed    1
# 3 changed    2
# 4  gained    1
# 5  gained    1
#                                                                                                                                                                                                   TF_highConf
# 1                                                                                                                                                                           HMGA1; HMGA2 (directAnnotation). 
# 2                                                                                                                                                                           HMGA1; HMGA2 (directAnnotation). 
# 3                                                                                                                                                                           HMGA1; HMGA2 (directAnnotation). 
# 4 DPF2; GTF2F1; IKZF3; MAZ; NFXL1; PRDM9; RBAK; SIN3A; TAF1; VEZF1; VEZF1; WRNIP1; ZBTB5; ZNF263; ZNF263; ZNF263; ZNF263; ZNF341; ZNF341; ZNF444; ZNF467; ZNF496; ZNF596; ZNF701; ZNF875 (directAnnotation). 
# 5 DPF2; GTF2F1; IKZF3; MAZ; NFXL1; PRDM9; RBAK; SIN3A; TAF1; VEZF1; VEZF1; WRNIP1; ZBTB5; ZNF263; ZNF263; ZNF263; ZNF263; ZNF341; ZNF341; ZNF444; ZNF467; ZNF496; ZNF596; ZNF701; ZNF875 (directAnnotation). 
#                  motif  NES
# 1 transfac_pro__M07320 3.47
# 2 transfac_pro__M07320 3.47
# 3 transfac_pro__M07320 3.47
# 4     metacluster_24.4 3.21
# 5     metacluster_24.4 3.21

source("/Users/felixyu/Documents/GitHub/TIPS/R/celltype_specific_weight_v10.R")

# g = graph_from_data_frame(final_table[,c('from','to')], directed = FALSE)


g_merged <- make_merged_TIPS_graph(subset(final_table, linkeage == "CM"),
    CHD = CHD,
    added_TF = setdiff(key_TFs, V(graph_list[[paste0("CTS_", CTS_ID)]])$name %>% toupper()), top_n_label = 5,
    g_string = graph_list[[paste0("CTS_", CTS_ID)]]
)

set.seed(2)
plot(
    g_merged,
    layout = layout_with_fr(g_merged, weights = NA),
    edge.curved = 0.15,
    vertex.size = 22,
    vertex.label.cex = 0.9,
    main = "Merged CMvsCP TIPS delta-edge reweighting"
)
mtext("CMvsCP edges labeled by delta (top abs_delta)", side = 1, line = -1, cex = 1.2)

dev.copy2pdf(file = "PPI_graph_merged_GRN_prediction_CTS_lateral_plate_mesoderm_CM_final.pdf")

g_merged <- make_merged_TIPS_graph(subset(final_table, linkeage == "CF"),
    CHD = CHD,
    added_TF = setdiff(key_TFs, V(graph_list[[paste0("CTS_", CTS_ID)]])$name %>% toupper()), top_n_label = 5,
    g_string = graph_list[[paste0("CTS_", CTS_ID)]]
)

set.seed(2)
plot(
    g_merged,
    layout = layout_with_fr(g_merged, weights = NA),
    edge.curved = 0.15,
    vertex.size = 22,
    vertex.label.cex = 0.9,
    main = "Merged CFvsCP TIPS delta-edge reweighting"
)
mtext("CFvsCP edges labeled by delta (top abs_delta)", side = 1, line = -1, cex = 1.2)

dev.copy2pdf(file = "PPI_graph_merged_GRN_prediction_CTS_lateral_plate_mesoderm_CF_final.pdf")


######## query the HMGA2 appeaing often in the results #########
## 1) Is HMGA2 actually highly expressed in the relevant cells? == yes, increased along cardiac trajectory from a->c
# p <- plotExpression(sce, features = c("HMGA2", "HMGA1"), x = celltype_col)
# p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
# dev.copy2pdf(file = "HMGA2_HMGA1_expression_across_clusters.pdf", width = 7, height = 5)

# cisTarget.res <- readRDS(file = "../cisTarget_targets_in_all_CTS.rds")
# cisTarget.res_HiG <- readRDS(file = "../cisTarget_targets_in_two_HiGs.rds")

# ## 2) How often does HMGA2 appear in TF_highConf among DEGs and CTS ?
# cisTarget.res[grepl("HMGA2", cisTarget.res$TF_highConf) | grepl("HMGA2", cisTarget.res$motif), ]$geneSet
# # [1] "endothelial.b"        "endothelial.b"        "cardiac.a"
# # [4] "presomiticMesoderm.a"
# cisTarget.res_HiG[grepl("HMGA2", cisTarget.res_HiG$TF_highConf) | grepl("HMGA2", cisTarget.res_HiG$motif), ]$geneSet
# # 0

# ## 3) 3) Is HMGA2 supported by many distinct motifs or only one recurring composite?
# hits_hmga <- subset(final_table, linkage == "CM" & grepl("HMGA2", TF_highConf))
# sort(tapply(hits_hmga$NES, hits_hmga$motif, max), decreasing = TRUE)
# # transfac_pro__M07320
# #                 3.81
