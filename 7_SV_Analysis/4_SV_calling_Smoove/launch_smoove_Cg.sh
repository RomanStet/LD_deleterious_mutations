#! /bin/bash
#SBATCH --cpus-per-task 5
#SBATCH --mem 5GB

module load smoove/0.2.6
module load bcftools/1.15.1

cd ${1}_results

smoove call -x --exclude .../Reference_genome/not_chromosomes.bed --name ${1} --fasta ${2} -p 5 --genotype ${3}

rm ${1}-smoove.genotyped.vcf

gunzip ${1}-smoove.genotyped.vcf.gz

bcftools query -f '%CHROM %ALT [%GT] %POS %CIPOS{0} %END %CIEND{1} %QUAL %FILTER \n' ${1}-smoove.genotyped.vcf > ${1}.table

grep '<DUP>' ${1}.table > ${1}_DUP.table





