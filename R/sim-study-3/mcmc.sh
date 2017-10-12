#! /bin/bash
#SBATCH --time=4-0:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=gpu
#SBATCH --error=info.err
#SBATCH --output=info.out
#SBATCH --mail-user=emittman@iastate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load R
module load cuda

R --vanilla --no-save <  run1.R #run an R script using R
