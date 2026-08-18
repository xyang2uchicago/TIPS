## 24.3_acat_CTS.cisTarget_merge_GRN.R — code_core_mouse_8 (CORE: no ATAC / no ChIP)
##
## Adapted from GSE87038_STRING/code/24.3...R (original hg38/ATAC pipeline) to run against
## the mm10 core results: sources code_core_mouse_8/24.0 (not code/24.0), reads/writes
## results_core_mouse_8/cisTarget_predicted_8 (not the original's shared results dir),
## and uppercases graph_list to match final_table's gene names (24.1's GRN-extraction
## section uppercases everything from that point on -- see that file's header note).

rebuild_mat <- FALSE
source(here::here("examples", "cardiac", "GSE87038", "GSE87038_STRING", "code_core_mouse_8", "24.0_acat_load_input_clean.R"))

CTS_ID <- "8"
celltype_col <- "label"

updir <- here::here("examples", "cardiac", "GSE87038", "GSE87038_STRING", "results_core_mouse_8")
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

library(igraph)
library(dplyr)

source(here::here("R", "celltype_specific_weight_v10.R"))

## final_table (written by 24.1) has uppercase gene names -- match graph_list to it.
for (i in seq_along(graph_list)) V(graph_list[[i]])$name <- toupper(V(graph_list[[i]])$name)

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

g_merged <- make_merged_TIPS_graph(subset(final_table, linkage == "CM"),
    CHD = CHD,
    added_TF = setdiff(key_TFs, V(graph_list[[paste0("CTS_", CTS_ID)]])$name), top_n_label = 5,
    g_string = graph_list[[paste0("CTS_", CTS_ID)]]
)

## pdf()+dev.off() instead of plot()+dev.copy2pdf() -- the latter needs an active screen
## graphics device, which doesn't exist under headless Rscript.
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

## CF merge, mirroring the CM merge above (original only did CM)
final_table_cf <- subset(final_table, linkage == "CF")
if (nrow(final_table_cf) > 0) {
    g_merged_cf <- make_merged_TIPS_graph(final_table_cf,
        CHD = CHD,
        added_TF = setdiff(key_TFs, V(graph_list[[paste0("CTS_", CTS_ID)]])$name), top_n_label = 5,
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
