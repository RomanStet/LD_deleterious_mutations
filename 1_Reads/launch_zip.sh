#! /bin/bash
#SBATCH --cpus-per-task 2
#SBATCH --mem 5GB

gzip -k $1

