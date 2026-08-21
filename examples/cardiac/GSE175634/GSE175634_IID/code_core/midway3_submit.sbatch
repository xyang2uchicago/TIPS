#!/bin/bash
#SBATCH --job-name=GSE175634_IID
#SBATCH --partition=gpu
#SBATCH --account=pi-xyang2
#SBATCH --gres=gpu:1
#SBATCH --mem=180G
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --output=/project/xyang2/felixy/TIPS/logs/GSE175634_IID_%j.out
#SBATCH --error=/project/xyang2/felixy/TIPS/logs/GSE175634_IID_%j.err

# Matches the PI's own validated invocation for this step (see the comment
# block at the top of 11.2.0_update_network_weights_clean_max.R, marked
# "(used)"). Adjust --account/--partition if felixy isn't on the pi-xyang2
# allocation, or switch to the bigmem alternative documented in the README
# this script ships next to (bigmem needs core_count=1 in run_TIPS_core.R
# instead of the default 4).

set -euo pipefail
cd /project/xyang2/felixy/TIPS

module load R/4.3.1

Rscript examples/cardiac/GSE175634/GSE175634_IID/code_core/run_TIPS_core.R
