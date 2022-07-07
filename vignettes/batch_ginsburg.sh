#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH -A iicd
#SBATCH -J KNDinh
#SBATCH -t 1:00:00
#     SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=knd2127@columbia.edu

module load R

pwd

echo "Launching an R run"
date

R CMD BATCH --no-save --vanilla ginsburg.r routput
