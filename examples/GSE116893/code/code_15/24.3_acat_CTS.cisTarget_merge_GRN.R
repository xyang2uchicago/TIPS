rebuild_mat <- FALSE
source("C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/code/code_15/24.0_acat_load_input_clean.R")

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

library(igraph)
library(dplyr)

key_TFs <- c("E2F3", "MYB", "ZNF730", "MTF2")

(updir <- getwd())
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

(file <- list.files(pattern = "final_table.tsv"))
# [1] "PPI_graph_GRN_prediction_CTS_15_dualpull_final_table.tsv"

final_table <- read.table(file = file, header = TRUE, sep = "\t")

# unique(final_table$linkage)  # [1] "CM"  (CTS_15 has no CF direction: all CF vcounts were 0 in 24.1)

for (lk in unique(final_table$linkage)) {
    g_merged <- make_merged_TIPS_graph(
        subset(final_table, linkage == lk),
        CHD = CHD,
        added_TF = setdiff(key_TFs, V(graph_list[[paste0("CTS_", CTS_ID)]])$name %>% toupper()),
        top_n_label = 5,
        g_string = graph_list[[paste0("CTS_", CTS_ID)]]
    )
    cat(lk, "merged graph: vcount =", vcount(g_merged), "ecount =", ecount(g_merged), "\n")
    # CM merged graph: vcount = 11  ecount = 16

    pdf(file = paste0("PPI_graph_merged_GRN_prediction_", CTS_name, "_", lk, "_final.pdf"))
    set.seed(2)
    plot(
        g_merged,
        layout = layout_with_fr(g_merged, weights = NA),
        edge.curved = 0.15,
        vertex.size = 22,
        vertex.label.cex = 0.9,
        main = paste0("Merged ", lk, "vsCP TIPS delta-edge reweighting")
    )
    mtext(paste0(lk, "vsCP edges labeled by delta (top abs_delta)"), side = 1, line = -1, cex = 1.2)
    dev.off()
    cat("Saved:", paste0("PPI_graph_merged_GRN_prediction_", CTS_name, "_", lk, "_final.pdf"), "\n")
}
