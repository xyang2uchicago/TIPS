## cardiac GSE175634 (human iPSC-CM differentiation) — STRING arm, code/ (non-centralized)
##
## Unlike code_core/run_TIPS_core.R, this driver does NOT go through
## tips_configure() -- these are the PI's original scripts, untouched, each
## with its own embedded USER INPUT block. This driver only sequences them;
## nothing about the sourced files themselves is edited.
##
## HPC-ONLY -- DO NOT RUN LOCALLY. 11.2.0 here hardcodes cores=4 against the
## same 230,786-cell sce.rds that crashed an 18GB local machine when the
## code_core/ equivalent was run with core_count=4 forking (see the
## GSE175634 HPC-only project memory). Left as the PI wrote it per the
## "no edits to sourced files" scope for this driver -- run only on an HPC
## bigmem/gpu-for-memory node, same as code_core/run_TIPS_core.R.
##
## This folder has no 12.0 step and no cisTarget 24.x step at all -- only
## 24.4.1_overlap_with_publications.R (a publications-comparison script) and
## the 24_acat_CTS.CP_pioneer_heatmap*.R variants (a different, separate
## analysis), neither of which plays the role of code_core's 24.x sequence.
## See code_core/12.0 and code_core/24.1 / 24.3 instead. Driver stops after
## 11.3, honestly reflecting what this folder actually has.

code_dir <- here::here("examples", "cardiac", "GSE175634", "GSE175634_STRING", "code")

source(file.path(code_dir, "11.1_STRINGweighted_CTS_cardiac_network.R"))
source(file.path(code_dir, "11.1.1_check_vertex_duplication.R"))
source(file.path(code_dir, "11.2.0_update_network_weights_clean_max.R"))
source(file.path(code_dir, "11.3_CTS_cardiac_network_ANND_pagerank.R"))

message("[run_TIPS_code] done (no 12.0/24.x cisTarget step in this folder; see code_core/ for those) -> ",
        here::here("examples", "cardiac", "GSE175634", "GSE175634_STRING", "results"))
