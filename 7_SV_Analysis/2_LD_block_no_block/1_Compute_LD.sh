#! /bin/bash
#SBATCH --partition fast
#SBATCH --cpus-per-task 2
#SBATCH --mem 6GB

module load r/4.1.0
srun Rscript .../7_SV_Analysis/2_LD_block_no_block/1_Compute_LD.r

