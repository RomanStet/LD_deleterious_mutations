#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 5GB

module load samtools/1.9

module load bwa/0.7.17

bwa mem -t $6 -R $5 $2 ${3}/${1}_trimmed_1.fastq.gz ${3}/${1}_trimmed_2.fastq.gz | samtools view -Sb -q 10 - > ${4}/${1}.fastq.gz.bam

samtools sort --threads $6 ${4}/${1}.fastq.gz.bam > ${4}/${1}.fastq.gz.sorted.bam

samtools index ${4}/${1}.fastq.gz.sorted.bam

rm ${4}/${1}.fastq.gz.bam
