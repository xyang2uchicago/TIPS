## cardiac GSE175634 (human iPSC-CM differentiation) — IID arm, code/ (non-centralized)
##
## Unlike code_core/run_TIPS_core.R, this driver does NOT go through
## tips_configure() -- these are the PI's original scripts, untouched, each
## with its own embedded USER INPUT block. This driver only sequences them;
## nothing about the sourced files themselves is edited.
##
## HPC-ONLY -- DO NOT RUN LOCALLY. 11.2.0 here hardcodes cores=4 against the
## same 230,786-cell sce.rds that crashed an 18GB local machine when the
## code_core/ equivalent was run with core_count=4 forking (see the
## GSE175634 HPC-only project memory). This file is left as the PI wrote it
## per the "no edits to sourced files" scope for this driver, so that
## hardcoded cores=4 is NOT fixed here the way it was in code_core/ -- run
## only on an HPC bigmem/gpu-for-memory node, same as code_core/run_TIPS_core.R.
##
## This folder has no 12.0 or 24.3-equivalent step -- see code_core/12.0 and
## code_core/24.3 instead for those. 24.0_acat_load_input_clean_v2.R and
## 24.1.1_..._CP_NEScut4.5.R stand in for 24.0/24.1 (later revisions of the
## same steps); 24.2_acat_CTS.ChIPseq_dualpull_clean.R is a genuinely
## different ChIP-seq validation analysis, not a cisTarget step -- included
## here because it's what this folder actually has under "24.2", not because
## it plays the same role as code_core's 24.x cisTarget sequence.

code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_IID", "code")

source(file.path(code_dir, "11.1_IIDweighted_CTS_cardiac_network.R"))
source(file.path(code_dir, "11.2.0_update_network_weights_clean_max.R"))
source(file.path(code_dir, "11.3_CTS_cardiac_network_ANND_pagerank.R"))
source(file.path(code_dir, "24.0_acat_load_input_clean_v2.R"))
source(file.path(code_dir, "24.1.1_acat_CTS.cisTarget_dualpull_clean_CP_NEScut4.5.R"))
source(file.path(code_dir, "24.2_acat_CTS.ChIPseq_dualpull_clean.R"))

message("[run_TIPS_code] done -> ", here::here("examples", "cardiac", "GSE175634", "GSE175634_IID", "results"))
