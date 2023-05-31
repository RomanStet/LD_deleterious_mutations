#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 5GB

cd ${4}

gunzip cr145_${1}.vcf

cd ${2}

module load gatk4/4.2.6.1

gatk SelectVariants \
     -R ${3} \
     -V ${4}cr145_${1}.vcf \
     --select-type-to-include SNP \
     -O ${2}cr145_${1}_snps_raw.vcf