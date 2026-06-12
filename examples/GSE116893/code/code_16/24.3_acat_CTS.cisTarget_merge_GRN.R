rebuild_mat <- FALSE
source("C:/Users/felix/Documents/GitHub/TIPS/examples/GSE116893/code/code_16/24.0_acat_load_input_clean.R")

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

library(igraph)
library(dplyr)

key_TFs <- c("CHD1", "MYCN", "EBF1", "JUND")

(updir <- getwd())
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

(file <- list.files(pattern = "final_table.tsv"))
# [1] "PPI_graph_GRN_prediction_CTS_16_dualpull_final_table.tsv"

final_table <- read.table(file = file, header = TRUE, sep = "\t")

# unique(final_table$linkage)  # [1] "CM" "CF"  (both directions present; final_table: 31x13)

# deduplicate: same pair can appear under multiple TF sections with identical delta
final_table <- final_table %>%
    mutate(pair_key = paste(pmin(from, to), pmax(from, to), sep = "||")) %>%
    group_by(linkage, pair_key) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    select(-pair_key)

for (lk in unique(final_table$linkage)) {
    g_merged <- make_merged_TIPS_graph(
        subset(final_table, linkage == lk),
        CHD = CHD,
        added_TF = setdiff(key_TFs, V(graph_list[[paste0("CTS_", CTS_ID)]])$name %>% toupper()),
        top_n_label = 5,
        g_string = graph_list[[paste0("CTS_", CTS_ID)]]
    )
    cat(lk, "merged graph: vcount =", vcount(g_merged), "ecount =", ecount(g_merged), "\n")
    # CF merged graph: vcount = 10  ecount = 19
    # CM merged graph: vcount = 8   ecount = 9

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
