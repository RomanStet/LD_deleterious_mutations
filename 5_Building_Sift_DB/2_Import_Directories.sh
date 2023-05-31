#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 50GB

TMPDIR=".../tmp"
mkdir -p "${TMPDIR}"
TMP="${TMPDIR}"
TEMP="${TMPDIR}"
export TMPDIR TMP TEMP

git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git scripts_to_build_SIFT_db

mv .../5_Building_Sift_DB/3_Run_SIFT.sh .../5_Building_Sift_DB/scripts_to_build_SIFT_db/

mv .../Reference_genome/Cr145.gff .../5_Building_Sift_DB/scripts_to_build_SIFT_db/

mv .../Reference_genome/Cr145.fa .../5_Building_Sift_DB/scripts_to_build_SIFT_db/

mv capsellarubella.config .../5_Building_Sift_DB/scripts_to_build_SIFT_db/

mv .../Protein_DB/uniref100.fasta .../5_Building_Sift_DB/scripts_to_build_SIFT_db/protein_DB/
