#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 50GB

TMPDIR=".../tmp"
mkdir -p "${TMPDIR}"
TMP="${TMPDIR}"
TEMP="${TMPDIR}"
export TMPDIR TMP TEMP

#### ######## create sift database
module load gffread/0.12.7
gffread Cr145.gff -T -o Cr145.gtf # convert gff to gtf
gzip Cr145.gtf 
gzip Cr145.fa # must end with fa, fas, fast
mkdir chr-src
mkdir gene-annotation-src
mv Cr145.fa.gz chr-src/
mv Cr145.gtf.gz gene-annotation-src/
module load perl/5.26.2 
module load sift4g/2.0.0
module load python/2.7
gunzip .../5_Building_Sift_DB/scripts_to_build_SIFT_db/protein_DB/uniref100.fasta.gz 
perl make-SIFT-db-all.pl -config capsellarubella.config # modify config file and run
gzip .../5_Building_Sift_DB/scripts_to_build_SIFT_db/protein_DB/uniref100.fasta
