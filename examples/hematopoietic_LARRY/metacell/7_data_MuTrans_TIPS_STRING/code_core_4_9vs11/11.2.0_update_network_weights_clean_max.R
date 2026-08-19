## 11.2.0_update_network_weights_clean_max.R — Felix code_core_13 (Larry paths only)

library(SingleCellExperiment)
library(igraph)

########## BEGINNING OF USER INPUT ##########
code_dir <- get0("code_dir", ifnotfound = file.path(here::here("examples", "hematopoietic_LARRY", "metacell", "7_data_MuTrans_TIPS_STRING"), "code_core_4_9vs11"))
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
source(file.path(code_dir, "mutrans_sce_xy.R"))
setwd(ppi_path)

specificity_methods <- c("combined")
core_count <- 1L ## Felix: use 1 on Windows
step1 <- TRUE
step2 <- TRUE
########## END OF USER INPUT ##########

source(file.path(tips_r_dir, paste0("celltype_specific_weight_v", celltype_specific_weight_version, ".R")))

Sys.setenv(SEURAT_RDS = file.path(here::here("examples", "hematopoietic_LARRY", "data"), "seu_attractor_MuTrans_HVG.rds"))
Sys.unsetenv("SCE_RDS")

sce <- load_larry_sce()
sce <- sync_sce_cluster_labels(sce, group_col = group_col)
validate_tips_sce(sce, context = "11.2.0")
assayName <- assert_tips_sce_assay(sce)

graph_list <- readRDS(file = graph_notsimp_file)

if (file.exists(correct_edges_file)) {
  correct_n_edges <- readRDS(correct_edges_file)
  for (g_name in unique(correct_n_edges$graph_id)) {
    rows <- subset(correct_n_edges, graph_id == g_name)
    vertices_str <- rows$vertex_index_to_remove
    vertices_str <- vertices_str[!is.na(vertices_str)]
    vertices_to_remove <- unlist(lapply(vertices_str, function(s) {
      as.numeric(strsplit(s, ",")[[1]])
    }))
    vertices_to_remove <- sort(unique(vertices_to_remove), decreasing = TRUE)
    if (length(vertices_to_remove) > 0 && g_name %in% names(graph_list)) {
      graph_list[[g_name]] <- delete_vertices(graph_list[[g_name]], vertices_to_remove)
    }
  }
}

graph_list <- lapply(graph_list, igraph::simplify, edge.attr.comb = "max")

## Felix IID 11.2 / STRING 11.1: STRING aliases vs SCE rownames (e.g. GATA2 vs Gata2)
graph_list <- harmonize_graph_names_to_sce(graph_list, sce)
report_hig_sce_overlap(graph_list, sce)
if ("HiG_10" %in% names(graph_list)) {
  assert_sce_gene_overlap(
    sce, V(graph_list[["HiG_10"]])$name,
    min_overlap = 100L, label = "HiG_10 graph vertices"
  )
}

if (step1) {
  colData(sce)$cluster <- colData(sce)[[cluster_labels]]
  network_specificity_list <- calculate_network_specificity(
    sce, graph_list, assayName,
    celltype_col = "cluster", method = "pearson",
    cores = core_count, shrink = TRUE, verbose = FALSE
  )
  saveRDS(network_specificity_list, "network_specificity_list.rds")
}

if (step2) {
  network_specificity_list <- readRDS("network_specificity_list.rds")
  for (net in names(network_specificity_list)) {
    network_specificity_list[[net]]$corexp_sign <-
      network_specificity_list[[net]]$corexp_sign
  }
  for (s in specificity_methods) {
    weighted_graph_list <- update_network_weights(
      graph_list, network_specificity_list,
      specificity_method = s, verbose = FALSE, cores = 1L
    )
    saveRDS(
      weighted_graph_list,
      file = paste0(db, "_", ppi_tag, "_graph_perState_simplified_", s, "weighted.rds")
    )
  }
}

message('To see your top-3 TFs after rerunning 12.0:','\n')
message( 'res_pr[[paste0("HiGCTS_", CTS_ID)]]  ')