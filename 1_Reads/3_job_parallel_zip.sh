#! /bin/bash

find . -name '*.fastq*' > .../1_Reads/fastq_files

(cat fastq_files; echo) | while read line
do
	FIRST=$(echo $line | awk '{print $1}')
	sbatch -p long launch_zip.sh $FIRST
done


