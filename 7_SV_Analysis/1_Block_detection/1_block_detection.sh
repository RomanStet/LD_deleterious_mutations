#!/bin/bash
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 2
#SBATCH --mem 20GB

module load r/4.1.0

Rscript launch_block_detection.r
