## Larry GSE140802 MuTrans TIPS — STRING, TAG 4_9vs11 (CP4 Baso / CF11 Meg)
## Full docs: README.md in this folder
## USER PARAMETERS below → tips_configure() → 11.1 … 24.3 → 25_venn

TAG <- "4_9vs11"
db_species <- 10090L  # 10090 mouse, 9606 human
group_col <- Sys.getenv("TIPS_GROUP_COL", "attractor")  # Seurat meta.data column (e.g. attractor, leiden_r0_6)
CP_cluster <- Sys.getenv("TIPS_CP_CLUSTER", "4")   # source
CM_cluster <- Sys.getenv("TIPS_CM_CLUSTER", "9")  # Baso
CF_cluster <- Sys.getenv("TIPS_CF_CLUSTER", "11")  # Meg
Real_names <- c(CP = "progenitor", CM = "Baso", CF = "Meg")
motif_target_strategy <- "merge"  # c("solo", "fam", "merge"); or TIPS_MOTIF_TARGET_STRATEGY env

tips_wd <- here::here("examples", "hematopoietic_LARRY", "metacell", "7_data_MuTrans_TIPS_STRING/")

code_dir <- paste0(tips_wd, "code_core_", TAG)

source(file.path(code_dir, "00_configuration.R"))
tips_configure(
  TAG = TAG,
  CP_cluster = CP_cluster,
  CM_cluster = CM_cluster,
  CF_cluster = CF_cluster,
  group_col = group_col,
  db_species = db_species,
  Real_names = Real_names,
  motif_target_strategy = motif_target_strategy,
  wd = tips_wd
)

## Log console output + messages to results folder (after tips_configure)
log_file <- file.path(results_dir, "running_message.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)

.old_message <- getOption("message")
options(message = function(...) {
  m <- if (exists(".makeMessage", mode = "function", where = asNamespace("base"))) {
    get(".makeMessage", asNamespace("base"))(..., domain = NULL)
  } else {
    paste(..., collapse = "")
  }
  cat(m, file = log_con, append = TRUE)
  if (is.function(.old_message)) {
    .old_message(...)
  } else {
    cat(m, file = stderr())
  }
})

on.exit({
  options(message = .old_message)
  sink()
  close(log_con)
}, add = TRUE)
message("[log] writing to ", log_file)

rebuildgraph <- TRUE
if (rebuildgraph) {
  source(file.path(code_dir, "mutrans_sce_xy.R"))   # reload fixed helpers

  source(file.path(code_dir, "11.1_STRINGweighted_CTS_network.R"))
  source(file.path(code_dir, "11.1.1_check_vertex_duplication.R"))
  source(file.path(code_dir, "11.2.0_update_network_weights_clean_max.R"))
  source(file.path(code_dir, "11.3_CTS_network_ANND_pagerank.R"))
}

source(file.path(code_dir, "12.0_rank_by_PageRank_BC.R"))
# To see your top-3 TFs after rerunning 12.0:
res_pr[["HiGCTS_4"]]   # PageRank: top 3 TFs
# # A tibble: 1 × 19
  # signature gene  PageRank EigenCentrality PPI_cat strength normalized.strength
  # <chr>     <chr>    <dbl>           <dbl> <fct>      <dbl>               <dbl>
# 1 HiGCTS_4  Klf1    0.0594           0.776 HiGCTS      4.11               0.133
#
source(file.path(code_dir, "24.1_acat_CTS.cisTarget_dualpull_clean.R"))
source(file.path(code_dir, "24.3_acat_CTS.cisTarget_merge_GRN.R"))
# [filter_cts_by_summary] CTS_summary_data.csv MCI_P<0.05 IC_g_local_p<0.05 IC_s_local_p<0.05 (NA passes) IC_s_delta_p<0.05 (NA passes) localmax=yes (last) -> keep 3/8: 12, 6, 4
# [filter_cts_by_summary] dropped: 10, 11, 1, 9, 13
# [24.0] uppercase mouse symbols done
# CF merged graph: vcount = 7 ecount = 11
# [24.3] done.

message("[run_TIPS_core] done -> ", results_dir)

source(file.path(code_dir, "25_venn_CTS_MuTrans_pub_markers.R"))
#   Basophil:
#   Mega:
#   Mast cell:
# [25] CTS ∩ TD:
#   A9: Ces2g
#   A11: Atp1b2, Hbb-bs, Sptb, Tmem27
#   A10: Fgf3, Tmem27

savehistory(file = file.path(code_dir, "R_run1_out.txt"))
