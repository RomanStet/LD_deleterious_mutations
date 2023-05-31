#! /bin/bash
#SBATCH --cpus-per-task 2
#SBATCH --mem 5GB

module load sra-tools/2.10.3

fasterq-dump $1