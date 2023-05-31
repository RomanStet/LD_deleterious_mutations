#!/bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 20GB

module load r/4.1.0
srun Rscript launch_Plot_Variant_Scores.r
