#!/bin/bash
#SBATCH --job-name=network_robustness
#SBATCH --partition=bigmem
#SBATCH --account=pi-xyang2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=16:00:00
#SBATCH --mem=180G
#SBATCH --output=job_output_%j.txt

# Load R module
module load R/4.3.1

# Run your R script
Rscript /project/xyang2/felixy/2018_weighted/code/11.2.1_weighted_network_robustness_simulation.R
