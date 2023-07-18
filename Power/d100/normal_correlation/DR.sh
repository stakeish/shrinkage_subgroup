#!/bin/bash
#SBATCH --job-name=shrinkage
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=60:00:00
#SBATCH --account=stats_dept1
#SBATCH --partition=standard

module load R
LC_ALL=C.UTF-8 R CMD BATCH shrinkage_simulation_1000.R simulation.out