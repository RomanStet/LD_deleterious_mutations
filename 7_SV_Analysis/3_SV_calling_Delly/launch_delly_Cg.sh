#! /bin/bash
#SBATCH --cpus-per-task 5
#SBATCH --mem 2GB

module load delly/0.8.3
module load bcftools/1.15.1

delly call -x .../Reference_genome/not_chromosomes.bed -o ${1}.bcf -g $2 $3

#bcftools view ${1}.bcf > ${1}.vcf

bcftools query -f '%CHROM %ALT [%GT] %POS %CIPOS{0} %END %CIEND{1} %QUAL %FILTER \n' ${1}.bcf > ${1}.table

grep '<DUP>' ${1}.table | grep PASS > ${1}_DUP.table

