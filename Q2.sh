#!/bin/bash
#SBATCH -c 1
#SBATCH --time=0:30:00
#SBATCH --mem=10G
#SBATCH --job-name="Q2"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=18arg5@queensu.ca
#SBATCH --output=Q2.o
#SBATCH --error=Q2.e

module load StdEnv/2020 r/4.2.2
Rscript Q2.R

