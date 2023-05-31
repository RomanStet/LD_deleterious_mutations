#! /bin/bash
#SBATCH --partition fast
#SBATCH --cpus-per-task 2
#SBATCH --mem 1GB

FOLDER_LD_ANALYSIS=".../6_LD_Analysis"

cd ${FOLDER_LD_ANALYSIS}

for i in $(seq 1 8)
do
	headers="POS CHROM_POS CHROM REF ALT"
	numcol=$(head -1 Cg_${i}_Sift.txt | wc -w)
	for((c=0; c<($numcol-5)*2; c++)) 
	do
		headers="$headers i"
	done
	headCg=($headers)

	tr '|' '/' < Cg_${i}_Sift.txt > Cg_${i}_Sift

	awk -F"/" '$1=$1' OFS="\t" Cg_${i}_Sift > temp
	tail -n +2 temp > temp_2
	( IFS=$'\t'; echo "${headCg[*]}"; cat temp_2 ) > Cg_${i}_Sift
done

for i in $(seq 1 8)
do
	headers="POS CHROM_POS CHROM REF ALT"
	numcol=$(head -1 Co_${i}_Sift.txt | wc -w)
	for((c=0; c<($numcol-5)*2; c++)) 
	do
		headers="$headers i"
	done
	headCo=($headers)

	tr '|' '/' < Co_${i}_Sift.txt > Co_${i}_Sift

	awk -F"/" '$1=$1' OFS="\t" Co_${i}_Sift > temp
	tail -n +2 temp > temp_2
	( IFS=$'\t'; echo "${headCo[*]}"; cat temp_2 ) > Co_${i}_Sift
done
rm temp temp_2



