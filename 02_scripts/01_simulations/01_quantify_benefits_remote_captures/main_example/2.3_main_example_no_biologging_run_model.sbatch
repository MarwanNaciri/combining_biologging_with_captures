#!/bin/bash
#
#SBATCH --job-name=0
#SBATCH --output=main_example_no_GPS.out
#SBATCH --error=main_example_no_GPS.err
#SBATCH --account polar_bear_sem_cr
#SBATCH --partition fast
#SBATCH --mem 5GB
#SBATCH --cpus-per-task 2
#SBATCH --time=0-0:15:00 

## load R
module load r/4.4.1
## run R script
Rscript main_example_no_biologging_run_model.R