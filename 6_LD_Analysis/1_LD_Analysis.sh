#! /bin/bash
#SBATCH --partition fast
#SBATCH --cpus-per-task 2
#SBATCH --mem 1GB

FOLDER_LD_ANALYSIS=".../6_LD_Analysis"
FOLDER_SIFT=".../Sift_annotation/"

cd ${FOLDER_LD_ANALYSIS}

tr '|' '/' < cr145_Cg_SNPs.table > temp1
sed 's/\.GT//g' temp1 > temp2
sed 's/\.\/\././g' temp2 > temp
sed 's/\./\.\/\./g' temp > temp2
mv temp2 cr145_Cg_SNPs.table
rm temp1 temp
	
tr '|' '/' < cr145_Co_SNPs.table > temp1
sed 's/\.GT//g' temp1 > temp2
sed 's/\.\/\././g' temp2 > temp
sed 's/\./\.\/\./g' temp > temp2
mv temp2 cr145_Co_SNPs.table
rm temp1 temp

cut -f1 cr145_Cg_SNPs.table > col1
cut -f2 cr145_Cg_SNPs.table > col2
paste -d"_" col1 col2 > pos
paste pos cr145_Cg_SNPs.table > temp
mv temp SNP_Cg
rm col1 col2 pos

cut -f1 cr145_Co_SNPs.table > col1
cut -f2 cr145_Co_SNPs.table > col2
paste -d"_" col1 col2 > pos
paste pos cr145_Co_SNPs.table > temp
mv temp SNP_Co
rm col1 col2 pos

for i in $(seq 1 8)
do
	head -1 SNP_Cg > Cg${i}
	grep SCF_${i}_ SNP_Cg >> Cg${i}
done

for i in $(seq 1 8)
do
	head -1 SNP_Co > Co${i}
	grep SCF_${i}_ SNP_Co >> Co${i}
done

cd ${FOLDER_SIFT}

for i in $(seq 1 8)
do
	sed -i -e '1s/#Position/POS/' Sift${i}
done






