#! /bin/bash

cd .../1_Reads/

find . -name '*_1.fastq.gz*' > .../2_Trimming/fastq_files

cd .../2_Trimming/

(cat fastq_files; echo) | while read line
do
      name=$(echo $line | awk '{print $1}' | rev | cut -d . -f 3- | tr -d '/.' | rev | cut -d'_' -f1)
      sbatch -p long launch_trimming.sh $name
done