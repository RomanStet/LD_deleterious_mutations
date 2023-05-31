#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 5GB

git clone --recursive https://github.com/rvaser/sift4g.git sift4g
cd sift4g/
make