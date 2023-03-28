#!/bin/bash
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --job-name="Q3"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=18arg5@queensu.ca
#SBATCH --output=Q3.o
#SBATCH --error=Q3.e

module load StdEnv/2020 r/4.2.2
Rscript Q3.R

