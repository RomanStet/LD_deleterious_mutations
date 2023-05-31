#!/bin/bash
#SBATCH --partition fast

module load r/4.1.0

srun Rscript .../7_SV_Analysis/5_SV_Analysis/launch_SV_analysis.r
