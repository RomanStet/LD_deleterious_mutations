#! /bin/bash

GENOME=".../Reference_genome/Cr145.fasta"
BAM_FOLDER=".../3_Mapping/"
HAP_CALL_FOLDER=".../4_Variant_Calling/1_Haplotype_Caller/"

cd ${BAM_FOLDER}

find . -name '*readgroup.bam.bai*' > ${HAP_CALL_FOLDER}/bam_files

cd ${HAP_CALL_FOLDER}

(cat bam_files; echo) | while read line
do
      name=$(echo $line | awk '{print $1}' | cut -d'.' -f2 | tr -d '/' )
      BAM_FILE="${BAM_FOLDER}${name}.fastq.gz.sorted_duplMarked_readgroup.bam"
      sbatch -p long launch_HaplotypeCaller.sh $name $GENOME $BAM_FILE $HAP_CALL_FOLDER
done




