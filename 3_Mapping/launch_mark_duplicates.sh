#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 20GB

TMPDIR=".../tmp"
#mkdir -p "${TMPDIR}"
TMP="${TMPDIR}"
TEMP="${TMPDIR}"
export TMPDIR TMP TEMP

module load picard/2.23.5

PICARD_FILE=".../Picard/picard.jar"

java -Xmx8g -jar ${PICARD_FILE} MarkDuplicates -I ${2}/${1}.fastq.gz.sorted_2.bam -O ${2}/${1}.fastq.gz.sorted_duplMarked.bam -M ${2}/${1}_duplMarked.metrics -MAX_FILE_HANDLES 1000 --TMP_DIR ${TMPDIR}  

java -Xmx4g -jar ${PICARD_FILE} AddOrReplaceReadGroups -I ${2}/${1}.fastq.gz.sorted_duplMarked.bam -O ${2}/${1}.fastq.gz.sorted_duplMarked_readgroup.bam --LB library --PL illumina --PU barcode --SM sample

module load samtools/1.6

samtools index ${2}/${1}.fastq.gz.sorted_duplMarked_readgroup.bam

rm ${2}/${1}.fastq.gz.sorted_duplMarked.bam