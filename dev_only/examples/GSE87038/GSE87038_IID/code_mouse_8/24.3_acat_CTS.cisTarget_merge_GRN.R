## 24.3_acat_CTS.cisTarget_merge_GRN.R — code_mouse_8 (IID, full ATAC/ChIP pipeline)
##
## Rewritten from GSE87038_IID/code/24.3...R: sources code_mouse_8/24.0 (mm10, not the
## original's broken-path code/24.0 -- see that file's header), drops the redundant
## graph_list re-read (24.0 already loads it correctly), fixes "linkeage" -> "linkage"
## to match the corrected column name written by code_mouse_8/24.1, uses standard
## "CTS_8" filenames (original used a stale "CTS_lateral_plate_mesoderm" name), and uses
## pdf()+dev.off() instead of plot()+dev.copy2pdf() (needs a screen device, unavailable
## under headless Rscript). graph_list is already uppercase (code_mouse_8/24.0 keeps the
## uppercase-mouse-as-proxy convention throughout).

rebuild_mat <- FALSE
source(here::here("examples", "cardiac", "GSE87038", "GSE87038_IID", "code_mouse_8", "24.0_acat_load_input_clean.R"))

CTS_ID # [1] "8"
length(CHD) # 295
celltype_col # [1] "label"

updir <- here::here("examples", "cardiac", "GSE87038", "GSE87038_IID", "results_mouse_8", "GSE181346_heart_scATAC")
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

library(igraph)
library(dplyr)

names(graph_list)

## extract the predicted TF regulators
(files <- list.files(pattern = "final.rds"))
key_TFs <- lapply(files, function(x) strsplit(x, "_")[[1]][3] %>% unlist()) %>%
  unlist() %>%
  unique()
key_TFs

####### generateing merged graph final table of the predicted subnetwork #################
(file <- list.files(pattern = "final_table.tsv"))

final_table <- read.table(file = file, header = TRUE, sep = "\t")
head(final_table)

source(here::here("R", "celltype_specific_weight_v10.R"))

g_merged <- make_merged_TIPS_graph(subset(final_table, linkage == "CM"),
  CHD = CHD,
  added_TF = setdiff(key_TFs, V(graph_list[[paste0("CTS_", CTS_ID)]])$name %>% toupper()), top_n_label = 5,
  g_string = graph_list[[paste0("CTS_", CTS_ID)]]
)

set.seed(2)
pdf(file = "PPI_graph_merged_GRN_prediction_CTS_8_CM_final.pdf")
plot(
  g_merged,
  layout = layout_with_fr(g_merged, weights = NA),
  edge.curved = 0.15,
  vertex.size = 22,
  vertex.label.cex = 0.9,
  main = "Merged CMvsCP TIPS delta-edge reweighting"
)
mtext("CMvsCP edges labeled by delta (top abs_delta)", side = 1, line = -1, cex = 1.2)
dev.off()

final_table_cf <- subset(final_table, linkage == "CF")
if (nrow(final_table_cf) > 0) {
  g_merged_cf <- make_merged_TIPS_graph(final_table_cf,
    CHD = CHD,
    added_TF = setdiff(key_TFs, V(graph_list[[paste0("CTS_", CTS_ID)]])$name %>% toupper()), top_n_label = 5,
    g_string = graph_list[[paste0("CTS_", CTS_ID)]]
  )

  set.seed(2)
  pdf(file = "PPI_graph_merged_GRN_prediction_CTS_8_CF_final.pdf")
  plot(
    g_merged_cf,
    layout = layout_with_fr(g_merged_cf, weights = NA),
    edge.curved = 0.15,
    vertex.size = 22,
    vertex.label.cex = 0.9,
    main = "Merged CFvsCP TIPS delta-edge reweighting"
  )
  mtext("CFvsCP edges labeled by delta (top abs_delta)", side = 1, line = -1, cex = 1.2)
  dev.off()
} else {
  message("[24.3] No CF rows in final_table — skipping CF merge")
}
