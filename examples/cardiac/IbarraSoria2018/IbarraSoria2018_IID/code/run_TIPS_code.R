## cardiac IbarraSoria2018 (E8.25 mouse gastrulation) — IID arm, code/ (non-centralized)
##
## These are the PI's original scripts, untouched, each
## with its own embedded USER INPUT block (species/CP_cluster/etc. hardcoded
## per-file). This driver only sequences them in the documented order;
## nothing about the sourced files themselves is edited.
##
## 24.1 and 24.3 each source 24.0 internally (with their own rebuild_mat
## flag), so 24.0 is not sourced separately here. No 11.1.1/24.2 in this
## folder.

code_dir <- here::here("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID", "code")

source(file.path(code_dir, "11.1_IIDweighted_CTS_cardiac_network.R"))
source(file.path(code_dir, "11.2.0_update_network_weights_clean_max.R"))
source(file.path(code_dir, "11.3_CTS_cardiac_network_ANND_pagerank.R"))
source(file.path(code_dir, "12.0_rank_by_PageRank_BC.R"))
source(file.path(code_dir, "24.1_acat_CTS.cisTarget_dualpull_clean.R"))
source(file.path(code_dir, "24.3_acat_CTS.cisTarget_merge_GRN.R"))

message("[run_TIPS_code] done -> ", here::here("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_IID", "results"))
