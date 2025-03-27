#!/bin/bash

cd /home/albertf/renga011/ExpressionPhenotypeProject/ #change home directory to project directory

module load R/4.4.0-openblas-rocky8 #load R module

mkdir output_GCQEC5050

# Loop from 1 to 100
for i in {1..100}
do
  # Call the R script with the current number as argument
  Rscript --vanilla R2C10_QECVsGC5050.R "$i" 
done