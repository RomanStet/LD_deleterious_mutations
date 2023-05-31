#!/bin/bash
#SBATCH --partition fast

module load r/4.1.0

srun Rscript ../7_SV_Analysis/5_SV_Analysis/script_match_duplic.r $1 $2 $3 $4 $5






