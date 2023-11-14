#!/bin/bash -l        
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=60g
#SBATCH --tmp=4g
#SBATCH --array=1-10
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=renga011@umn.edu 

cd ~/ExpressionPhenotypeProject #change home directory to project directory
module load R/3.5.2_mkl #load R module
mkdir output_randomH2Estimates #make directory to collect all the colocalization output files
Rscript --vanilla getRandomH2Estimates.R "$SLURM_ARRAY_TASK_ID" > output_randomH2Estimates/output_randomH2Estimates"$SLURM_ARRAY_TASK_ID".txt #run the randomH2Estimates R code
