#!/bin/bash
#SBATCH --account=co_biostat
#SBATCH --job-name=para_test
#SBATCH --partition=savio2
#SBATCH --nodes=2
#SBATCH --cpus-per-task=10
#SBATCH --time=10:00:00

#SBATCH --mail-user=haodong_li@berkeley.edu
#SBATCH --mail-type=ALL
module load r
Rscript test.R
