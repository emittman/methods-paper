#! /bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=gpu
#SBATCH --error=msg/array%a.err
#SBATCH --output=msg/array%a.out
#SBATCH --mail-user=emittman@iastate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load R
module load cuda

cd /ptmp/emittman/methods-paper/R/sim-study-1/

R CMD BATCH --vanilla --no-save small_test.R #run an R script using R
