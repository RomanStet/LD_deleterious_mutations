#! /bin/bash

FASTQFOLDER=".../1_Reads"
LINES=$(cat double_runs.txt | wc -l)

for i in $(seq 1 $LINES)
do
	if (( $i % 2 == 1 ))
	then
	SRA1=$(sed "${i}q;d" double_runs)
	fi
	
	if (( $i % 2 == 0 ))
	then
	SRA2=$(sed "${i}q;d" double_runs)
	rm ${SRA2}_1.fastq.gz ${SRA2}_2.fastq.gz
	sbatch -p fast launch_merge_double_runs.sh $i $FASTQFOLDER $SRA1 $SRA2
	fi
done
