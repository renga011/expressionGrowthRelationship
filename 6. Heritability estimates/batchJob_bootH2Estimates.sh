#!/bin/bash -l        
#SBATCH --time=48:00:00
#SBATCH --ntasks=128
#SBATCH --mem=100g
#SBATCH --tmp=4g
#SBATCH --array=1-50
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=renga011@umn.edu 

cd /home/albertf/renga011/ExpressionPhenotypeProject/ #change home directory to project directory
module load R/4.4.0-openblas-rocky8 #load R module
mkdir output_bootH2Estimates #make directory to collect all the colocalization output files
Rscript --vanilla 6.1_getH2ExplainedByHotspots_bootH2Estimates.R "$SLURM_ARRAY_TASK_ID" > output_bootH2Estimates/output_bootH2Estimates"$SLURM_ARRAY_TASK_ID".txt #run the bootH2Estimates R code
