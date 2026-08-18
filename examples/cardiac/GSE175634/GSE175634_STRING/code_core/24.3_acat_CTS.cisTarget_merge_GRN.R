rebuild_mat <- FALSE
source(here::here("examples", "cardiac", "GSE175634", "GSE175634_STRING", "code_core", "24.0_acat_load_input_clean.R"))

source(paste0('https://raw.githubusercontent.com/xyang2uchicago/TIPS/refs/heads/main/R/celltype_specific_weight_v', celltype_specific_weight_version, '.R'))

library(igraph)
library(dplyr)

(updir <- getwd())
setwd(paste0(updir, "/cisTarget_predicted_", CTS_ID))

(file <- list.files(pattern = "final_table.tsv"))
# [1] "PPI_graph_GRN_prediction_CTS_CP_dualpull_final_table.tsv"

final_table <- read.table(file = file, header = TRUE, sep = "\t")

key_TFs <- lapply(
    list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final.rds"),
    function(x) strsplit(x, "_")[[1]][3] %>% unlist()
) %>% unlist() %>% unique()
key_TFs
# [1] "KLF6"  "MEIS1"

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
    # CF merged graph: vcount = 18 ecount = 20
    # CM merged graph: vcount = 6 ecount = 4

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
    # Saved: PPI_graph_merged_GRN_prediction_CTS_CP_CF_final.pdf
}
