rebuild_mat = FALSE
source("/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/code_3/24.0_acat_load_input_clean.R")

CTS_ID = '8'
celltype_col  = "label"

updir <- "/Users/felixyu/Documents/GitHub/TIPS/examples/GSE87038/results/GSE181346_heart_scATAC_C1"
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

library(igraph)
library(dplyr)
library(igraph)

# db_specifc_input_path <- "/Users/felixyu/Documents/GitHub/TIPS/examples/IbarraSoria2018/results"
# file <- paste0(db_specifc_input_path, "/PPI_weight/IbarraSoria2018_STRING_graph_perState_simplified_combinedweighted.rds")
# graph_list <- readRDS(file)
names(graph_list)


## extract the predicted TF regulators
(files <- list.files(pattern = "final.rds"))
# [1] "PPI_graph_HMGA2_GRN_prediction_CTS_8_CF_final.rds"
# [2] "PPI_graph_HMGA2_GRN_prediction_CTS_8_CM_final.rds"
# [3] "PPI_graph_RARB_GRN_prediction_CTS_8_CM_final.rds" 
key_TFs <- lapply(files, function(x) strsplit(x, "_")[[1]][3] %>% unlist()) %>%
    unlist() %>%
    unique()
key_TFs
# [1] "HMGA2" "RARB" 


####### generateing merged graph final table of the predicted subnetwork #################
(file <- list.files(pattern = "final_table.tsv"))
#  "PPI_graph_GRN_prediction_CTS_8_dualpull_final_table.tsv"


final_table <- read.table(file = file, header = TRUE, sep = "\t")
head(final_table)
#   linkage   from     to         w1          w2     delta abs_delta direction
# 1      CF   DAB2  HMGA2  0.0000000 0.201972551 0.2019726 0.2019726  increase
# 2      CM  HMGA2 MPPED2 -0.1528846 0.205327212 0.3582118 0.3582118  increase
# 3      CM   FGF8  HMGA2 -0.1243642 0.002014388 0.1263786 0.1263786  increase
# 4      CM MPPED2   RARB  0.0000000 0.118093337 0.1180933 0.1180933  increase
#    status rank
# 1  gained    1
# 2 changed    1
# 3 changed    2
# 4  gained    1
#                                                                                                                                                                                                   TF_highConf
# 1                                                                                                                                                                           HMGA1; HMGA2 (directAnnotation). 
# 2                                                                                                                                                                           HMGA1; HMGA2 (directAnnotation). 
# 3                                                                                                                                                                           HMGA1; HMGA2 (directAnnotation). 
# 4 DPF2; GTF2F1; IKZF3; MAZ; NFXL1; PRDM9; RBAK; SIN3A; TAF1; VEZF1; VEZF1; WRNIP1; ZBTB5; ZNF263; ZNF263; ZNF263; ZNF263; ZNF341; ZNF341; ZNF444; ZNF467; ZNF496; ZNF596; ZNF701; ZNF875 (directAnnotation). 
#                  motif  NES
# 1 transfac_pro__M07320 3.47
# 2 transfac_pro__M07320 3.47
# 3 transfac_pro__M07320 3.47
# 4     metacluster_24.4 3.21

source("../../../../../R/celltype_specific_weight_v10.R")

# g = graph_from_data_frame(final_table[,c('from','to')], directed = FALSE)

g_merged <- make_merged_TIPS_graph(subset(final_table, linkage == "CM"),
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

dev.copy2pdf(file="PPI_graph_merged_GRN_prediction_CTS_8_CM_final.pdf")



# ###############################################################
# ## overview figure for publication 
# library(scater)

# markers = toupper(c('Mesp1', 'Hand1', 'Smarcd3', 'Isl1', 'Tbx18','Krt19', 'Col1a2', 'Tnnt2' ) ) 
# setdiff(markers, rownames(sce))

# #'Tagln'
# plotReducedDim(sce, 'UMAP', colour = markers[[i]], by_exprs_values = "logcounts")
