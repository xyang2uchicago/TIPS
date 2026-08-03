## 24.3_acat_CTS.cisTarget_merge_GRN.R — Felix STRING (Larry MuTrans)

code_dir <- get0("code_dir", ifnotfound = "F:/projects/TIPS/source/GSE140802_lineage_tracking/7_data_MuTrans_TIPS_STRING/code_core_4_9vs11")
if (!exists("rebuild_mat")) rebuild_mat <- FALSE
source(file.path(code_dir, "24.0_acat_load_input_clean.R"))
source(file.path(tips_r_dir, paste0("celltype_specific_weight_v", celltype_specific_weight_version, ".R")))

suppressPackageStartupMessages({ library(igraph); library(dplyr) })

updir <- getwd()
setwd(file.path(updir, paste0("cisTarget_predicted_", CTS_ID)))

file <- list.files(pattern = "dualpull_final_table\\.tsv$")
if (!length(file)) stop("Run 24.1 first")
final_table <- read.table(file[1], header = TRUE, sep = "\t")

key_TFs <- lapply(
  list.files(pattern = "PPI_graph_.*_GRN_prediction_.*_final\\.rds$"),
  function(x) strsplit(x, "_")[[1]][3]
) %>% unlist() %>% unique()

final_table <- final_table %>%
  mutate(pair_key = paste(pmin(from, to), pmax(from, to), sep = "||")) %>%
  group_by(linkage, pair_key) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  select(-pair_key)

tips_real_name <- function(role, real_names = get0("Real_names")) {
  role <- as.character(role)
  if (is.null(real_names)) return(role)
  nms <- names(real_names)
  if (length(nms) && role %in% nms) return(as.character(real_names[[role]]))
  idx <- match(role, c("CP", "CM", "CF"))
  if (!is.na(idx) && idx <= length(real_names)) return(as.character(real_names[[idx]]))
  role
}

for (lk in unique(final_table$linkage)) {
  g_merged <- make_merged_TIPS_graph(
    subset(final_table, linkage == lk), CHD = CHD,
    added_TF = setdiff(key_TFs, toupper(V(graph_list[[CTS_name]])$name)),
    top_n_label = 5, g_string = graph_list[[CTS_name]]
  )
  message(lk, " merged graph: vcount = ", vcount(g_merged), " ecount = ", ecount(g_merged))
  desc_name <- tips_real_name(lk)
  cp_name <- tips_real_name("CP")
  pdf(paste0("PPI_graph_merged_GRN_prediction_", CTS_name, "_", lk, "_final.pdf"))
  set.seed(2)
  plot(
    g_merged, layout = layout_with_fr(g_merged, weights = NA), edge.curved = 0.15,
    vertex.size = 22, vertex.label.cex = 0.9,
    main = paste0("Merged ", desc_name, " vs ", cp_name, " TIPS delta-edge reweighting")
  )
  mtext(
    paste0(desc_name, " vs ", cp_name, " edges labeled by delta (top abs_delta)"),
    side = 1, line = -1, cex = 1.2
  )
  dev.off()
}
message("[24.3] done.")
