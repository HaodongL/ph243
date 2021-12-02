#!/bin/bash
#SBATCH --account=co_biostat
#SBATCH --job-name=hal_test
#SBATCH --partition=savio2
#SBATCH --nodes=4
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00

#SBATCH --mail-user=haodong_li@berkeley.edu
#SBATCH --mail-type=ALL

module load r
Rscript test.R
