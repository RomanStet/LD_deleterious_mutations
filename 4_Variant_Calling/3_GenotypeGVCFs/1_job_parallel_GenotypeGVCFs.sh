#! /bin/bash

GVCF_FOLDER=".../4_Variant_Calling/3_GenotypeGVCFs/"
DBImport_FOLDER=".../4_Variant_Calling/2_GenomicsDBImport/"
GENOME=".../Reference_genome/Cr145.fasta"
VCF_FOLDER=".../4_Variant_Calling/1_Haplotype_Caller/"

cp ${DBImport_FOLDER}scaffold_file ${GVCF_FOLDER}scaffold_file

cd ${VCF_FOLDER}

(cat bam_files; echo) | while read line
do
      name=$(echo $line | awk '{print $1}' | cut -d'.' -f2 | tr -d '/' )
      cd ${path}
      sbatch -p fast launch_cp_vcf.sh $name $VCF_FOLDER $DBImport_FOLDER 
done

cd ${GVCF_FOLDER}

(cat scaffold_file; echo) | while read line
do
     scaffold=$(echo $line | awk '{print $1}')
     sbatch -p long launch_GenotypeGVCFs.sh $scaffold $GVCF_FOLDER $GENOME $FOLDER_DBImport
done
