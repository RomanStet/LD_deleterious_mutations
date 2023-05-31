#! /bin/bash
#SBATCH --partition fast
#SBATCH --cpus-per-task 2
#SBATCH --mem 6GB

module load r/4.1.0
srun Rscript .../6_LD_Analysis/launch_LD_Analysis_4.r

