#!/bin/bash
#
#SBATCH --job-name=CR
#SBATCH --output=CR_model_without_biologging.out
#SBATCH --error=CR_model_without_biologging.err
#SBATCH --account account_name_here
#SBATCH --partition long
#SBATCH --mem 50GB
#SBATCH --cpus-per-task 15
#SBATCH --time=1-12:00:00 

## load R
module load r/3.6.3
## run R script
Rscript CR_model_without_biologging.R