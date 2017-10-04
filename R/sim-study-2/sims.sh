#! /bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --partition=gpu
#SBATCH --ntasks-per-node=1
#SBATCH --error=info.err
#SBATCH --output=info.out
#SBATCH --mail-user=emittman@iastate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load R
module load cuda

R --vanilla --no-save <  sims.R #run an R script using R
