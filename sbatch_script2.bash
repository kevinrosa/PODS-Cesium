#!/bin/bash

# SBATCH -p ivy-batch
# SBATCH -p sandy-batch
# SBATCH --account=epscor-condo
#SBATCH -p bigmem
#SBATCH -t 40:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=16G
#SBATCH --job-name=matlabbin
# SBATCH --mail-type=ALL # type of email notification BEGIN,END,FAIL,ALL
# SBATCH --mail-user=kevin_rosa@uri.edu


matlab < script2.m; exit
