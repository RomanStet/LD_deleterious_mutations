#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 15GB

module load gatk4/4.2.6.1

gatk --java-options "-Xmx27g" HaplotypeCaller  \
-R ${2} \
-I ${3} \
-O ${4}${1}.g.vcf.gz \
-ERC GVCF

