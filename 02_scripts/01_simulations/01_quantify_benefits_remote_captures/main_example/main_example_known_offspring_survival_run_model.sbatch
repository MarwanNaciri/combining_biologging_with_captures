#!/bin/bash
#
#SBATCH --job-name=0
#SBATCH --output=main_example_offspring_survival.out
#SBATCH --error=main_example_offspring_survival.err
#SBATCH --account polar_bear_sem_cr
#SBATCH --partition fast
#SBATCH --mem 25GB
#SBATCH --cpus-per-task 4
#SBATCH --time=0-0:15:00 

## load R
module load r/4.4.1
## run R script
Rscript main_example_known_offspring_survival_run_model.R