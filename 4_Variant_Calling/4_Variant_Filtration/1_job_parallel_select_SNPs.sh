#! /bin/bash

FOLDER_VF=".../4_Variant_Calling/4_Variant_Filtration/"
FOLDER_GenotypeGVCFs=".../4_Variant_Calling/3_GenotypeGVCFs/"
GENOME=".../Reference_genome/Cr145.fasta"

cp ${FOLDER_GenotypeGVCFs}scaffold_file ${FOLDER_VF}scaffold_file

(cat scaffold_file; echo) | while read line
do
      scaffold=$(echo $line | awk '{print $1}')
      sbatch -p long launch_select_SNPs.sh $scaffold $FOLDER_VF $GENOME $FOLDER_GenotypeGVCFs
done
