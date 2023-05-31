#! /bin/bash

BAM_FOLDER=".../3_Mapping"

(cat fastq_files; echo) | while read line
do
      name=$(echo $line | awk '{print $1}' | cut -d'_' -f1 | tr -d '/.' )
      sbatch -p fast launch_mark_duplicates.sh $name $BAM_FOLDER
done

