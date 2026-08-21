## cardiac GSE175634 (human iPSC-CM differentiation) — IID arm, code/ (non-centralized)
##
## These are the PI's original scripts, untouched, each
## with its own embedded USER INPUT block.

code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_IID", "code")

source(file.path(code_dir, "11.1_IIDweighted_CTS_cardiac_network.R"))
source(file.path(code_dir, "11.2.0_update_network_weights_clean_max.R"))
source(file.path(code_dir, "11.3_CTS_cardiac_network_ANND_pagerank.R"))
source(file.path(code_dir, "24.0_acat_load_input_clean_v2.R"))
source(file.path(code_dir, "24.1.1_acat_CTS.cisTarget_dualpull_clean_CP_NEScut4.5.R"))
source(file.path(code_dir, "24.2_acat_CTS.ChIPseq_dualpull_clean.R"))

message("[run_TIPS_code] done -> ", here::here("examples", "cardiac", "GSE175634", "GSE175634_IID", "results"))
