#! /bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem 2GB

module load trimmomatic/0.39

trimmomatic PE -threads 6 -phred33 .../1_Reads/${1}_1.fastq.gz .../1_Reads/${1}_2.fastq.gz .../2_Trimming/${1}_trimmed_1.fastq.gz .../2_Trimming/${1}_single_1.fastq.gz .../2_Trimming/${1}_trimmed_2.fastq.gz .../2_Trimming/${1}_single_2.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fasta:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

