## leiden_r0_8  CP=11 (Undiff, Meg-primed) / CM=10 (Baso) / CF=17 (Meg)
## Her code_core files byte-identical; only wd / inputs / group column differ.
TAG <- "11_10vs17"
db_species <- 10090L
group_col  <- Sys.getenv("TIPS_GROUP_COL", "leiden_r0_8")
CP_cluster <- Sys.getenv("TIPS_CP_CLUSTER", "11")
CM_cluster <- Sys.getenv("TIPS_CM_CLUSTER", "10")
CF_cluster <- Sys.getenv("TIPS_CF_CLUSTER", "17")
Real_names <- c(CP = "progenitor", CM = "Baso", CF = "Meg")
motif_target_strategy <- "merge"

tips_wd  <- Sys.getenv("TIPS_WD")
code_dir <- file.path(tips_wd, paste0("code_core_", TAG))

source(file.path(code_dir, "00_configuration.R"))
tips_configure(
  TAG = TAG, CP_cluster = CP_cluster, CM_cluster = CM_cluster,
  CF_cluster = CF_cluster, group_col = group_col, db_species = db_species,
  Real_names = Real_names, motif_target_strategy = motif_target_strategy,
  wd = tips_wd
)

rebuild_mat <- TRUE
rebuild_mat <- TRUE
source(file.path(code_dir, "24.1_acat_CTS.cisTarget_dualpull_clean.R"))
source(file.path(code_dir, "24.3_acat_CTS.cisTarget_merge_GRN.R"))
