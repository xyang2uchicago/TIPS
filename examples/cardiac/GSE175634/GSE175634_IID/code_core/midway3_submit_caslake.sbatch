#!/bin/bash
#SBATCH --job-name=GSE175634_IID_caslake
#SBATCH --partition=caslake
#SBATCH --account=pi-xyang2
#SBATCH --mem=180G
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=/project/xyang2/felixy/TIPS/logs/GSE175634_IID_caslake_%j.out
#SBATCH --error=/project/xyang2/felixy/TIPS/logs/GSE175634_IID_caslake_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL

# Third hedge alongside midway3_submit.sbatch (gpu) and midway3_submit_bigmem.sbatch
# (bigmem). caslake is RCC's large general-purpose partition -- far more
# node capacity than gpu/bigmem, so a job stuck on REASON=Priority there is
# likely to clear faster. Config matches the PI's own documented (if not
# explicitly "(used)") comment in 11.2.0: -p caslake --mem=180G -c 1.
#
# TIPS_RESULTS_TAG routes output to results_core_caslake/, separate from the
# gpu (results_core/) and bigmem (results_core_bigmem/) runs -- safe to run
# all three concurrently with zero collision risk.
#
# Set --mail-user below (or via `scontrol update JobId=<id> MailUser=...`
# after submission) to get notified without polling squeue.

set -euo pipefail
cd /project/xyang2/felixy/TIPS

module load R/4.3.1

export TIPS_CORE_COUNT=1
export TIPS_RESULTS_TAG=_caslake

# One-time resume flag: this arm already completed 11.2.0's ~1hr step1/step2
# (calculate_network_specificity/update_network_weights) on a prior caslake
# run that then failed later (missing brainGraph package, now installed).
# Skips redoing that step and reuses the already-saved
# network_specificity_list.rds / *_combinedweighted.rds from
# results_core_caslake/. Remove this line (or set to 0) for a genuinely
# fresh run where those cached files don't exist yet.
export TIPS_SKIP_SPECIFICITY=1

Rscript examples/cardiac/GSE175634/GSE175634_IID/code_core/run_TIPS_core.R
