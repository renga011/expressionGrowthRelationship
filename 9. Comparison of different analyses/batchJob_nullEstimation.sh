#!/bin/bash -l        
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=60g
#SBATCH --tmp=4g
#SBATCH --array=1-46
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=renga011@umn.edu 

cd ~/ExpressionPhenotypeProject #change home directory to project directory
module load R/3.5.2_mkl #load R module
mkdir output_nullEstimation #make directory to collect all the output files
Rscript --vanilla computeNulls.R "$SLURM_ARRAY_TASK_ID" > output_nullEstimation/output_nullEstimation"$SLURM_ARRAY_TASK_ID".txt #run the computeNulls R code
