#!/bin/bash -l        
#SBATCH --time=24:00:00
#SBATCH --ntasks=128
#SBATCH --mem=100g
#SBATCH --tmp=4g
#SBATCH --array=1-50
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=renga011@umn.edu 

cd /home/albertf/renga011/ExpressionPhenotypeProject/ #change home directory to project directory
module load R/4.4.0-openblas-rocky8 #load R module
mkdir output_GCHEC5050 #make directory to collect all the output files
Rscript --vanilla R2C10_HECVsGC5050.R "$SLURM_ARRAY_TASK_ID" > output_GCHEC5050/output_GCHEC5050"$SLURM_ARRAY_TASK_ID".txt #run the GCHEC5050 R code
