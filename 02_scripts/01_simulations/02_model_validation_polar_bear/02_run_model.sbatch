#!/bin/bash
#
#SBATCH --job-name=CR_sim
#SBATCH --output=02_run_model.out
#SBATCH --error=02_run_model.err
#SBATCH --account account_name_here
#SBATCH --partition fast
#SBATCH --mem 50GB
#SBATCH --cpus-per-task 10
#SBATCH --time=0-5:00:00 

## load R
module load r/3.6.3
## run R script
Rscript 02_run_model.R