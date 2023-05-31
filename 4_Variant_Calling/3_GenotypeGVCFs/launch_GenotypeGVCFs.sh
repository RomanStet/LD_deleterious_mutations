#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 15GB

module load gatk4/4.2.6.1

gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R ${3} \
    -V gendb://${4}/database_${1} \
    -new-qual \
    --tmp-dir /scratch2/umi3614/bedim/rstetsenko/these/tmp \
    -O ${2}/cr145_${1}.vcf.gz
