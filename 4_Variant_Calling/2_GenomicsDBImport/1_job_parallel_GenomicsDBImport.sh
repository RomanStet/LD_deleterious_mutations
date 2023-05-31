#! /bin/bash

GENOME=".../Reference_genome/Cr145.fasta"
GENOME_FOLDER=".../Reference_genome/"
VCF_FOLDER=".../4_Variant_Calling/1_Haplotype_Caller/"
GDB_FOLDER=".../4_Variant_Calling/2_GenomicsDBImport/"

module load bedtools/2.30.0
module load samtools/1.15.1

cd ${GENOME_FOLDER}

samtools faidx Cr145.fasta

awk '{FS="\t"}; {print $1 FS "0" FS $2}' Cr145.fasta.fai | head -8 > chromosomes.bed

awk '{FS="\t"}; {print $1 FS "0" FS $2}' Cr145.fasta.fai | tail +9 > not_chromosomes.bed

cut -f1 chromosomes.bed > ${GDB_FOLDER}scaffold_file

cd ${VCF_FOLDER}

(cat bam_files; echo) | while read line
do
      name=$(echo $line | awk '{print $1}' | cut -d'.' -f2 | tr -d '/' )
      cd ${GDB_FOLDER}
      sbatch -p fast launch_cp_vcf.sh $name $VCF_FOLDER $GDB_FOLDER
done

cd ${GDB_FOLDER}

(cat scaffold_file; echo) | while read line
do
      scaffold=$(echo $line | awk '{print $1}')
      sbatch -p long launch_GenomicsDBImport.sh $scaffold $GDB_FOLDER
done
