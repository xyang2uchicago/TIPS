#!/bin/bash
#SBATCH --job-name=GSE175634_STRING_bigmem
#SBATCH --partition=bigmem
#SBATCH --account=pi-xyang2
#SBATCH --mem=500G
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=/project/xyang2/felixy/TIPS/logs/GSE175634_STRING_bigmem_%j.out
#SBATCH --error=/project/xyang2/felixy/TIPS/logs/GSE175634_STRING_bigmem_%j.err

# bigmem alternative to midway3_submit.sbatch, matching the PI's own
# documented fallback config (mem=500G, cores=1). TIPS_CORE_COUNT=1 makes
# run_TIPS_core.R pass core_count=1 into calculate_network_specificity/
# update_network_weights instead of the default 4 used on the gpu partition.
#
# TIPS_RESULTS_TAG routes this run's output to results_core_bigmem/ instead
# of results_core/, so this is safe to submit and run fully concurrently
# with midway3_submit.sbatch (the gpu variant) for the same arm -- no
# queue-watching or scancel-the-loser required. Compare the two
# results_core*/ folders afterward; keep whichever, or both.

set -euo pipefail
cd /project/xyang2/felixy/TIPS

module load R/4.3.1

export TIPS_CORE_COUNT=1
export TIPS_RESULTS_TAG=_bigmem
Rscript examples/cardiac/GSE175634/GSE175634_STRING/code_core/run_TIPS_core.R
