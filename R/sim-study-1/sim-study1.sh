#! /bin/bash
#SBATCH --time=0-10:00:00
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

R --vanilla --no-save < small_test.R #run an R script using R
