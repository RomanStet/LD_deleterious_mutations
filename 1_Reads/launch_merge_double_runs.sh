#! /bin/bash
#SBATCH --cpus-per-task 5
#SBATCH --mem 20GB

SRA_1_1="${2}${3}_1.fastq.gz"
SRA_1_2="${2}${3}_2.fastq.gz"
SRA_2_1="${2}${4}_1.fastq.gz"
SRA_2_2="${2}${4}_2.fastq.gz"
	
cat ${SRA_1_1} ${SRA_2_1} >${2}${3}_1_new.fastq.gz
cat ${SRA_1_2} ${SRA_2_2} >${2}${3}_2_new.fastq.gz

cat ${2}${3}_1_new.fastq.gz >${2}${3}_1.fastq.gz
cat ${2}${3}_2_new.fastq.gz >${2}${3}_2.fastq.gz

rm ${2}${3}_1_new.fastq.gz ${2}${3}_2_new.fastq.gz


