#!/bin/bash

ID_IND=".../Accession_SRA_Number/ID_ind"
BLOCK_FOLDER="../7_SV_Analysis/1_Block_detection/"
SV_ANALYSIS_FOLDER="../7_SV_Analysis/5_SV_Analysis/"

species="Cg"

for i in {1..8..1}
do
	mkdir match_duplic_${species}${i}
	cd match_duplic_${species}${i}
	cat ${BLOCK_FOLDER}${species}${i}_bloc_stats_2.txt > temp_${species}${i}
	sed -i '1d' temp_${species}${i}
	
	(cat temp_${species}${i}; echo) | while read line
	do
		start=$(echo $line | awk '{print $3}')
		end=$(echo $line | awk '{print $4}')
		num_bloc=$(echo $line | awk '{print $2}')
		sbatch -p fast ${SV_ANALYSIS_FOLDER}launch_extract_duplic.sh $species $i $start $end $num_bloc
	done
	cd ..
done

