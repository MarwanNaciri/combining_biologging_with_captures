#!/bin/bash
#
#SBATCH --job-name=0
#SBATCH --output=partridge_example_no_biologging.out
#SBATCH --error=partridge_example_no_biologging.err
#SBATCH --account polar_bear_sem_cr
#SBATCH --partition fast
#SBATCH --mem 25GB
#SBATCH --cpus-per-task 4
#SBATCH --time=0-0:10:00 

## load R
module load r/4.4.1
## run R script
Rscript partridge_example_no_biologging_run_model.R