#!/bin/bash
#SBATCH --account=co_biostat
#SBATCH --job-name=hal_test
#SBATCH --partition=savio3
#SBATCH --nodes=8
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

#SBATCH --mail-user=haodong_li@berkeley.edu
#SBATCH --mail-type=ALL

module load r
Rscript run_para_boothal.R
