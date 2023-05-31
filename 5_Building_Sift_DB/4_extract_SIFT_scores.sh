#! /bin/bash
#SBATCH --partition fast
#SBATCH --cpus-per-task 5
#SBATCH --mem 15GB

OUTFOLDER=".../Sift_annotation"
SIFT_FOLDER=".../5_Building_Sift_DB/scripts_to_build_SIFT_db/Cr145/"

for i in $(seq 1 8)
do
	gunzip ${SIFT_FOLDER}SCF_${i}.gz

	head -1 ${SIFT_FOLDER}SCF_${i} > ${OUTFOLDER}Sift${i}
	grep CDS ${SIFT_FOLDER}SCF_${i} >> ${OUTFOLDER}Sift${i}

	cut -f1-3,11-13 ${OUTFOLDER}Sift${i} > ${OUTFOLDER}temp
	mv ${OUTFOLDER}temp ${OUTFOLDER}Sift${i}
done



