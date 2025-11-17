#!/bin/bash
#
#SBATCH --job-name=PB_no_biologging
#SBATCH --output=polar_bear_no_biologging.out
#SBATCH --error=polar_bear_no_biologging.err
#SBATCH --account polar_bear_sem_cr
#SBATCH --partition fast
#SBATCH --mem 25GB
#SBATCH --cpus-per-task 4
#SBATCH --time=0-3:00:00 

## load R
module load r/4.4.1
## run R script
Rscript polar_bear_no_biologging_run_model.R