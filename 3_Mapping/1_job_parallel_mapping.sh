#! /bin/bash

GENOME=".../Reference_Genome/Cr145.fasta"
TRIMMED_FOLDER=".../2_Trimming"
BAM_FOLDER=".../3_Mapping"
ID="@RG\tID:ind\tSM:ind\tPL:Illumina"
NCPU="5"

cd $TRIMMEDFOLDER

find . -name '*_trimmed_1.fastq.gz*' > .../1_Reads/fastq_files

cd $BAMFOLDER

(cat fastq_files; echo) | while read line
do
      name=$(echo $line | awk '{print $1}' | cut -d'_' -f1 | tr -d '/.' )
      sbatch -p long launch_mapping.sh $name $GENOME $TRIMMED_FOLDER $BAM_FOLDER $ID $NCPU
done




