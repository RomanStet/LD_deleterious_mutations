#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 8
#SBATCH --mem 25GB

TMPDIR=".../tmp"
mkdir -p "${TMPDIR}"
TMP="${TMPDIR}"
TEMP="${TMPDIR}"
export TMPDIR TMP TEMP

module load gatk4/4.2.6.1

gatk --java-options "-Xmx25g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path ${2}database_${1} \
       --batch-size 60 \
       -L ${1} \
       --sample-name-map ${2}ID_vcf.sample_map \
       --tmp-dir ${TMPDIR} \
       --reader-threads 8
