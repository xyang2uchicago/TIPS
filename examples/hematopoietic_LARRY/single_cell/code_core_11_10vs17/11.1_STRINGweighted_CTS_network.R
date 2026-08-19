## 11.1_STRINGweighted_CTS_network.R — Felix TIPS_core (LARRY MuTrans, attractor labels)
##
## setwd() not required — code_dir is set below in USER INPUT.

suppressPackageStartupMessages({
  library(dplyr)
  library(scran)
  library(SingleCellExperiment)
  library(igraph)
  library(STRINGdb)
})

########## BEGINNING OF USER INPUT ##########
code_dir <- get0("code_dir", ifnotfound = file.path(here::here("examples", "hematopoietic_LARRY", "single_cell"), "code_core_11_10vs17"))
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
setwd(results_dir)

## ---- DEG cutoffs: edit here (output filenames follow) ----
logFC.cut <- 0.5
fdr.cut   <- 0.05

deg_paths <- tips_deg_paths(logFC.cut, fdr.cut, min_prop, data_dir)
deg_tag <- deg_paths$deg_tag
deg_rdata <- deg_paths$deg_rdata
markers_rdata <- deg_paths$markers_rdata
message(
  "[11.1] DEG cutoffs: min.prop=", min_prop,
  " logFC>", logFC.cut, " FDR<", fdr.cut,
  " -> ", basename(deg_rdata)
)

## ---- BioTIP CTS significance (CTS_summary_data.csv) ----
## Filter order (all must pass; focal CTS_ID always kept):
##   1. MCI_P              < cts_mci_p_cut
##   2. IC_g_local_p       < cts_ic_g_p_cut
##   3. IC_s_local_p       < cts_ic_s_p_cut  (NA passes)
##   4. IC_s_delta_p       < cts_ic_s_delta_p_cut  (NA passes)
##   5. localmax           == "yes"  (lastly)
cts_mci_p_cut <- 0.05
cts_ic_g_p_cut <- 0.05
cts_ic_s_p_cut <- 0.05
cts_ic_s_delta_p_cut <- 0.05
cts_require_localmax <- TRUE
########## END OF USER INPUT ##########

Sys.setenv(SEURAT_RDS = file.path(here::here("examples", "hematopoietic_LARRY", "data"), "seu_attractor_MuTrans_HVG.rds"))

sce <- load_larry_sce()
sce <- sync_sce_cluster_labels(sce, group_col = group_col)
assert_tips_sce_assay(sce)
message("[11.1] sce assay=", tips_sce_assay, " (Seurat RNA data layer)")
message("[11.1] sce attractors (n=", length(unique(sce$label)), "): ",
        paste(sort(unique(sce$label)), collapse = ", "))

if (!file.exists(biotip_cts_rdata)) {
  stop("Missing BioTIP CTS: ", biotip_cts_rdata)
}
load(biotip_cts_rdata)
stopifnot(exists("CTS.Lib.Symbol"))
CTS_all <- CTS.Lib.Symbol
message("[11.1] CTS modules (all): ", paste(names(CTS_all), collapse = ", "))
# [11.1] CTS modules: 12, 6, 10, 11, 4, 1, 9, 13

CTS <- filter_cts_by_summary(
  CTS_all,
  summary_path = biotip_cts_summary_csv,
  mci_p_cut = cts_mci_p_cut,
  ic_g_local_p_cut = cts_ic_g_p_cut,
  ic_s_local_p_cut = cts_ic_s_p_cut,
  ic_s_delta_p_cut = cts_ic_s_delta_p_cut,
  require_localmax = cts_require_localmax,
  always_keep = CTS_ID
)
if (!CTS_ID %in% names(CTS)) {
  stop(
    "Focal CTS_", CTS_ID, " failed BioTIP significance filter. ",
    "Relax cts_*_p_cut / cts_require_localmax in 11.1 or check ", biotip_cts_summary_csv,
    call. = FALSE
  )
}
message("[11.1] CTS modules (significant): ", paste(names(CTS), collapse = ", "))
#[11.1] CTS modules (significant): 12, 6, 4
 

required_attractors <- normalize_cluster_labels(c(CP_cluster, CM_cluster, CF_cluster))
sce_labels <- normalize_cluster_labels(sce$label)
if (!all(required_attractors %in% unique(sce_labels))) {
  stop(
    "Missing attractors in sce$label: ",
    paste(setdiff(required_attractors, unique(sce_labels)), collapse = ", "),
    "\n  Found: ", paste(sort(unique(sce_labels)), collapse = ", "),
    "\n  Re-run from: sce <- load_larry_sce(); sce <- sync_sce_cluster_labels(sce)",
    call. = FALSE
  )
}
colData(sce)$label <- sce_labels

message("[11.1] findMarkers: min.prop=", min_prop, " logFC>", logFC.cut, " FDR<", fdr.cut)

markers.up <- findMarkers(
  sce, test = "t", groups = sce$label,
  min.prop = min_prop, direction = "up"
)

DEG <- list()
unique_CTS_ID <- names(CTS)
unique_CTS_ID <- unique_CTS_ID[!(grepl("\\.", unique_CTS_ID) &
  grepl("^[0-9]+$", sub("^[^.]*\\.", "", unique_CTS_ID)))]

for (i in c(setdiff(names(markers.up), names(CTS)), unique_CTS_ID)) {
  interesting.up <- markers.up[[i]]
  DEG[[i]] <- subset(interesting.up, summary.logFC > logFC.cut & FDR < fdr.cut) %>% rownames()
}

saveRDS(DEG, file = deg_rdata)
saveRDS(markers.up, file = markers_rdata)
deg_cutoffs_info <- list(
  logFC.cut = logFC.cut, fdr.cut = fdr.cut, min.prop = min_prop,
  deg_tag = deg_tag, deg_file = deg_rdata, markers_file = markers_rdata
)
saveRDS(deg_cutoffs_info, file.path(results_dir, paste0("deg_cutoffs_", deg_tag, ".rds")))
saveRDS(deg_cutoffs_info, file.path(results_dir, "deg_cutoffs_latest.rds"))
saveRDS(
  list(
    mci_p_cut = cts_mci_p_cut,
    ic_g_local_p_cut = cts_ic_g_p_cut,
    ic_s_local_p_cut = cts_ic_s_p_cut,
    ic_s_delta_p_cut = cts_ic_s_delta_p_cut,
    require_localmax = cts_require_localmax,
    summary_path = biotip_cts_summary_csv,
    sig_ids = names(CTS)
  ),
  file.path(results_dir, "cts_filter_latest.rds")
)
message("[11.1] DEG -> ", deg_rdata)
message("[11.1] markers.up -> ", markers_rdata)
message("[11.1] key states HiG: ",
        paste(sprintf("%s=%d", c(CP_cluster, CM_cluster, CF_cluster),
                      lengths(DEG[c(CP_cluster, CM_cluster, CF_cluster)])),
              collapse = "  "))
# key states HiG: 4=575  10=515  11=479

if (!dir.exists(string_ppin_dir)) {
  stop("STRING PPIN directory missing: ", string_ppin_dir)
}

string_db <- STRINGdb$new(
  version = "12.0", species = db_species,
  score_threshold = 200, network_type = "full",
  input_directory = string_ppin_dir
)

map_string_subgraph <- function(diff_exp, string_db) {
  mapped <- string_db$map(diff_exp, "symbol", removeUnmappedRows = TRUE)
  hits <- mapped$STRING_id
  graph <- string_db$get_subnetwork(hits)
  if (vcount(graph) < length(hits)) {
    missing_node <- setdiff(hits, V(graph)$name)
    if (length(missing_node) > 0) {
      mapped <- mapped[!mapped$STRING_id %in% missing_node, , drop = FALSE]
    }
  }
  stopifnot(all(mapped[match(V(graph)$name, mapped$STRING_id), ]$STRING_id == V(graph)$name))
  V(graph)$name <- mapped[match(V(graph)$name, mapped$STRING_id), ]$symbol
  ## Felix STRING 11.1: STRING upper-case aliases -> restore SCE/DEG symbol casing
  graph <- remap_graph_symbols_to_reference(graph, diff_exp$symbol)
  V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$summary.logFC
  V(graph)$FDR <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$FDR
  E(graph)$weight <- E(graph)$combined_score / 1000
  graph <- delete_edge_attr(graph, "combined_score")
  list(graph = graph, mapped = mapped)
}

graph_list <- list()
for (i in names(DEG)) {
  diff_exp <- markers.up[[i]]
  diff_exp$symbol <- rownames(diff_exp)
  diff_exp <- subset(diff_exp, summary.logFC > logFC.cut & FDR < fdr.cut)
  res <- map_string_subgraph(diff_exp, string_db)
  graph_list[[paste0("HiG_", i)]] <- res$graph
}

for (i in names(CTS)) {
  j <- if (grepl("\\.", i) && grepl("^[0-9]+$", sub("^[^.]*\\.", "", i))) {
    sub("\\..*$", "", i)
  } else {
    i
  }
  deg_table <- markers.up[[j]]
  if (is.null(deg_table)) {
    graph_list[[paste0("HiGCTS_", i)]] <- make_empty_graph()
    next
  }
  deg_table$symbol <- rownames(deg_table)
  deg_table <- deg_table[deg_table$symbol %in% CTS[[i]], ]
  diff_exp <- subset(deg_table, summary.logFC > logFC.cut & FDR < fdr.cut)
  if (nrow(diff_exp) == 0L) {
    graph_list[[paste0("HiGCTS_", i)]] <- make_empty_graph()
    next
  }
  res <- map_string_subgraph(diff_exp, string_db)
  graph_list[[paste0("HiGCTS_", i)]] <- res$graph
}

for (i in names(CTS)) {
  CTS_genes <- data.frame(symbol = CTS[[i]])
  mapped <- string_db$map(CTS_genes, "symbol", removeUnmappedRows = TRUE)
  graph <- string_db$get_subnetwork(mapped$STRING_id)
  stopifnot(all(mapped[match(V(graph)$name, mapped$STRING_id), ]$STRING_id == V(graph)$name))
  V(graph)$name <- mapped[match(V(graph)$name, mapped$STRING_id), ]$symbol
  graph <- remap_graph_symbols_to_reference(graph, CTS[[i]])
  E(graph)$weight <- E(graph)$combined_score / 1000
  graph <- delete_edge_attr(graph, "combined_score")
  graph_list[[paste0("CTS_", i)]] <- graph
}
# Warning:  we couldn't map to STRING 4% of your identifiersWarning:  we couldn't map to STRING 7% of your identifiersWarning:  we couldn't map to STRING 17% of your identifiers> 

df_graph_info <- data.frame(
  name = names(graph_list),
  vcount = vapply(graph_list, vcount, numeric(1)),
  ecount = vapply(graph_list, ecount, numeric(1)),
  stringsAsFactors = FALSE
)
print(df_graph_info, row.names = FALSE)
#       name vcount ecount
#      HiG_0    664  26030
#      HiG_1    706  30776
#     HiG_10    504  14968
#     HiG_11    465  13260
#     HiG_13    701  28232
#      HiG_2    704  35768
#      HiG_3    517  16422
#      HiG_5    572  31352
#      HiG_7    409  11376
#      HiG_8    650  24642
#      HiG_9    526  18308
#     HiG_12    162   3278
#      HiG_6    570  17858
#      HiG_4    555  14906
#  HiGCTS_12     28     96
#   HiGCTS_6     15     16
#   HiGCTS_4     32    216
#     CTS_12     48    186
#      CTS_6     38     54
#      CTS_4     73    476

saveRDS(graph_list, file = graph_notsimp_file)
graph_list <- readRDS(graph_notsimp_file)
graph_simp <- lapply(graph_list, igraph::simplify, edge.attr.comb = "max")
dups <- names(which(sapply(graph_simp, function(g) any(duplicated(V(g)$name)))))
if (length(dups)) {
  warning(
    "[11.1] duplicate vertex names after simplify: ", paste(dups, collapse = ", "),
    "\n  Run source('11.1.1_check_vertex_duplication.R') before 11.2.0"
  )
}

message("[11.1] saved -> ", graph_notsimp_file)
message("[11.1] next: source('11.1.1_check_vertex_duplication.R')")
message("[11.1] Venn QC: source('25_venn_CTS_MuTrans_pub_markers.R')")
