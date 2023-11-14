#!/bin/bash -l        
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=60g
#SBATCH --tmp=4g
#SBATCH --array=1-1000
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=renga011@umn.edu 

cd ~/ExpressionPhenotypeProject #change home directory to project directory
module load R/3.5.2_mkl #load R module
mkdir output_colocalization_LD95 #make directory to collect all the colocalization output files
Rscript --vanilla cisEQTLs_LD95.R "$SLURM_ARRAY_TASK_ID" > output_colocTest/output_colocTest"$SLURM_ARRAY_TASK_ID".txt #run the colocalization R code
