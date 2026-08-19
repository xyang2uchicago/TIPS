## 11.1_IIDweighted_CTS_cardiac_network.R — Felix TIPS_core (LARRY MuTrans, IID PPIN)
##
## IID global graph: iid2025_human_mouse_conserved_global.rds (see 00_configuration.R: iid_global_rds)

suppressPackageStartupMessages({
  library(dplyr)
  library(scran)
  library(SingleCellExperiment)
  library(igraph)
})

########## BEGINNING OF USER INPUT ##########
code_dir <- get0("code_dir", ifnotfound = file.path(here::here("examples", "hematopoietic_LARRY", "metacell", "7_data_MuTrans_TIPS_IID"), "code_core_4_9vs11"))
source(file.path(code_dir, "00_configuration.R"))
ensure_tips_configured(code_dir)
source(file.path(code_dir, "mutrans_sce_xy.R"))
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

if (!file.exists(iid_global_rds)) {
  stop("Missing IID global PPIN: ", iid_global_rds, call. = FALSE)
}
message("[11.1] IID global PPIN: ", iid_global_rds)
g_iid_global <- readRDS(iid_global_rds)
message("[11.1] g_iid_global: ", vcount(g_iid_global), " nodes, ", ecount(g_iid_global), " edges")

get_iid_subnetwork <- function(g_global, hits) {
  hits <- toupper(hits)
  hits <- unique(hits[!is.na(hits) & hits != ""])
  hits <- intersect(hits, V(g_global)$name)
  if (length(hits) < 2) {
    return(make_empty_graph())
  }
  induced_subgraph(g_global, vids = hits)
}

map_iid_subgraph <- function(diff_exp, g_global) {
  graph <- get_iid_subnetwork(g_global, diff_exp$symbol)
  if (vcount(graph) > 0) {
    graph <- remap_graph_symbols_to_reference(graph, diff_exp$symbol)
    V(graph)$weight <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$summary.logFC
    V(graph)$FDR <- diff_exp[match(V(graph)$name, diff_exp$symbol), ]$FDR
  }
  graph
}

graph_list <- list()
for (i in names(DEG)) {
  diff_exp <- markers.up[[i]]
  diff_exp$symbol <- rownames(diff_exp)
  diff_exp <- subset(diff_exp, summary.logFC > logFC.cut & FDR < fdr.cut)
  graph_list[[paste0("HiG_", i)]] <- map_iid_subgraph(diff_exp, g_iid_global)
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
  graph_list[[paste0("HiGCTS_", i)]] <- map_iid_subgraph(diff_exp, g_iid_global)
}

for (i in names(CTS)) {
  hits <- CTS[[i]]
  graph <- get_iid_subnetwork(g_iid_global, hits)
  if (vcount(graph) > 0) {
    graph <- remap_graph_symbols_to_reference(graph, hits)
  }
  graph_list[[paste0("CTS_", i)]] <- graph
}

df_graph_info <- data.frame(
  name = names(graph_list),
  vcount = vapply(graph_list, vcount, numeric(1)),
  ecount = vapply(graph_list, ecount, numeric(1)),
  stringsAsFactors = FALSE
)
print(df_graph_info, row.names = FALSE)

saveRDS(graph_list, file = graph_notsimp_file)
graph_list <- readRDS(graph_notsimp_file)
graph_simp <- lapply(graph_list, igraph::simplify, edge.attr.comb = "max")
dups <- names(which(sapply(graph_simp, function(g) any(duplicated(V(g)$name)))))
if (length(dups)) {
  warning(
    "[11.1] duplicate vertex names after simplify: ", paste(dups, collapse = ", "),
    "\n  IID graphs are usually unique — inspect HiG graphs before 11.2.0"
  )
}

message("[11.1] saved -> ", graph_notsimp_file)
message("[11.1] next: source('11.2.0_update_network_weights_clean_max.R')")
message("[11.1] Venn QC: source('25_venn_CTS_MuTrans_pub_markers.R')")
