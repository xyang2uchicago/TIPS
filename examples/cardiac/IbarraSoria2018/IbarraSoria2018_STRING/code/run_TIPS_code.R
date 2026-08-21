## cardiac IbarraSoria2018 (E8.25 mouse gastrulation) — STRING arm, code/ (non-centralized)
##
## Unlike code_core/run_TIPS_core.R, this driver does NOT go through
## tips_configure() -- these are the PI's original scripts, untouched, each
## with its own embedded USER INPUT block (species/CP_cluster/etc. hardcoded
## per-file). This driver only sequences them in the documented order;
## nothing about the sourced files themselves is edited.
##
## 24.1, 24.2, and 24.3 each source 24.0 internally (with their own
## rebuild_mat flag), so 24.0 is not sourced separately here. No 11.1.1 in
## this folder.
##
## 24.4.1_overlap_with_publications.R is intentionally excluded -- it's a
## publications-comparison script, not part of the core cisTarget sequence.

code_dir <- here::here("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_STRING", "code")

source(file.path(code_dir, "11.1_STRINGweighted_CTS_cardiac_network.R"))
source(file.path(code_dir, "11.2.0_update_network_weights_clean_max.R"))
source(file.path(code_dir, "11.3_CTS_cardiac_network_ANND_pagerank.R"))
source(file.path(code_dir, "12.0_rank_by_PageRank_BC.R"))
source(file.path(code_dir, "24.1_acat_CTS.cisTarget_dualpull_clean.R"))
source(file.path(code_dir, "24.2_acat_CTS.ChIPseq_dualpull_clean.R"))
source(file.path(code_dir, "24.3_acat_CTS.cisTarget_merge_GRN.R"))

message("[run_TIPS_code] done -> ", here::here("examples", "cardiac", "IbarraSoria2018", "IbarraSoria2018_STRING", "results"))
