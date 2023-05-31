#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 2
#SBATCH --mem 2GB

cp ${2}${1}.g.vcf.gz ${3}${1}.g.vcf.gz
cp ${2}${1}.g.vcf.gz.tbi ${3}${1}.g.vcf.gz.tbi