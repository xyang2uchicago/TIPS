rebuild_mat <- TRUE
source(here::here("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID", "code", "24.0_acat_load_input_clean.R"))

CTS_ID # 'cardiac.a'
length(CHD) # 295
celltype_col # [1] "subcelltype"

updir <- here::here("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID", "results", "GSE181346_heart_scATAC")
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

library(igraph)
library(dplyr)

db_specifc_input_path <- here::here("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID", "results/")
file <- paste0(db_specifc_input_path, "PPI_weight/IbarraSoria2018_IID_graph_perState_simplified_combinedweighted.rds")
graph_list <- readRDS(file)
names(graph_list)


## extract the predicted TF regulators
(files <- list.files(pattern = "final.rds"))
# [1] "PPI_graph_CBFB_GRN_prediction_CTS_cardiac.a_CM_final.rds"
# [2] "PPI_graph_ETV2_GRN_prediction_CTS_cardiac.a_CM_final.rds"
# [3] "PPI_graph_GATA4_GRN_prediction_CTS_cardiac.a_CM_final.rds"
# [4] "PPI_graph_HMGA2_GRN_prediction_CTS_cardiac.a_CF_final.rds"
# [5] "PPI_graph_HMGA2_GRN_prediction_CTS_cardiac.a_CM_final.rds"
# [6] "PPI_graph_PRDM6_GRN_prediction_CTS_cardiac.a_CM_final.rds"
# [7] "PPI_graph_RBFOX2_GRN_prediction_CTS_cardiac.a_CM_final.rds"
key_TFs <- lapply(files, function(x) strsplit(x, "_")[[1]][3] %>% unlist()) %>%
  unlist() %>%
  unique()
key_TFs # [1] "CBFB"   "ETV2"   "GATA4"  "HMGA2"  "PRDM6"  "RBFOX2"


####### generateing merged graph final table of the predicted subnetwork #################
(file <- list.files(pattern = "final_table.tsv"))
#  "PPI_graph_GRN_prediction_CTS_cardiac.a_dualpull_final_table.tsv"


final_table <- read.table(file = file, header = TRUE, sep = "\t")
head(final_table)
#   linkeage  from    to          w1            w2        delta   abs_delta
# 1       CM  CBFB GATA6  0.03651142 -0.0056903109 -0.042201727 0.042201727
# 2       CM  CBFB  TBX5 -0.01286985 -0.0002101566  0.012659694 0.012659694
# 3       CM GATA6  TBX5  0.01000000  0.0007727979 -0.009227202 0.009227202
# 4       CM  ETV2 SFRP5  0.00000000 -0.0404144174 -0.040414417 0.040414417
# 5       CM  ETV2 SFRP5  0.00000000 -0.0404144174 -0.040414417 0.040414417
# 6       CM GATA4 MYOCD  0.03826759  0.0001234696 -0.038144118 0.038144118
#   direction  status rank                          TF_highConf
# 1  decrease changed    1            CBFB (directAnnotation).
# 2  increase changed    2            CBFB (directAnnotation).
# 3  decrease changed    3
# 4  decrease  gained    1       ERF; ETV7 (directAnnotation).
# 5  decrease  gained    1      ETV2; ETV7 (directAnnotation).
# 6  decrease changed    1 EP300; GTF2IRD1 (directAnnotation).
#                                                    motif  NES
# 1                                             hdpi__CBFB 4.67
# 2                                             hdpi__CBFB 4.67
# 3                                                          NA
# 4 taipale_tf_pairs__ERF_ETV7_NSMGGACGGAYNTCCKSN_CAP_repr 3.18
# 5     taipale_tf_pairs__ETV2_ETV7_NSMGGACGGAYNTCCKSN_CAP 3.08
# 6                                      tfdimers__MD00175 4.36

source(here::here("R", "celltype_specific_weight_v10.R"))

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

dev.copy2pdf(file = "PPI_graph_merged_GRN_prediction_CTS_cardiac.a_CM_final.pdf")


######## query the HMGA2 appeaing often in the results #########
## 1) Is HMGA2 actually highly expressed in the relevant cells? == yes, increased along cardiac trajectory from a->c
p <- plotExpression(sce, features = c("HMGA2", "HMGA1"), x = celltype_col)
p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
dev.copy2pdf(file = "HMGA2_HMGA1_expression_across_clusters.pdf", width = 7, height = 5)

## 2) How often does HMGA2 appear in TF_highConf among DEGs and CTS ?
cisTarget.res[grepl("HMGA2", cisTarget.res$TF_highConf) | grepl("HMGA2", cisTarget.res$motif), ]$geneSet
# [1] "cardiac.a"     "endothelial.b"
cisTarget.res_HiG[grepl("HMGA2", cisTarget.res_HiG$TF_highConf) | grepl("HMGA2", cisTarget.res_HiG$motif), ]$geneSet
# 0

## 3) 3) Is HMGA2 supported by many distinct motifs or only one recurring composite?
hits_hmga <- subset(final_table, linkeage == "CM" & grepl("HMGA2", TF_highConf))
sort(tapply(hits_hmga$NES, hits_hmga$motif, max), decreasing = TRUE)
# transfac_pro__M07320
#                 3.81


###############################################################
## overview figure for publication
library(scater)

markers <- toupper(c("Mesp1", "Hand1", "Smarcd3", "Isl1", "Tbx18", "Krt19", "Col1a2", "Tnnt2"))
setdiff(markers, rownames(sce))

#' Tagln'
plotReducedDim(sce, "UMAP", colour = markers[[i]], by_exprs_values = "logcounts")
