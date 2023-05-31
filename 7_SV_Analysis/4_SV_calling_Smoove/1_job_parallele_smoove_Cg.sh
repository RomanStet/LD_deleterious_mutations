#! /bin/bash

GENOME=".../Reference_genome/Cr145.fasta"
ID_IND=".../Accession_SRA_Number/ID_ind"

grep Cg $ID_IND | awk '{print $2}' > SRA_Cg

(cat SRA_Cg; echo) | while read line
do
      name=$(echo $line | awk '{print $1}' )
      ID=$(grep $name $ID_IND | awk '{print $3}')
      BAMFILE="${BAMFOLDER}${name}.fastq.gz.sorted_duplMarked_readgroup.bam"
      mkdir ${ID}_results
      sbatch -p fast launch_smoove_Cg.sh $ID $GENOME $BAMFILE
done

(cat SRA_Cg; echo) | while read line
do
	name=$(echo $line | awk '{print $1}' )
	ID=$(grep $name $ID_IND | awk '{print $3}')
	cat ${ID}_results/${ID}_DUP.table > temp
	cut -d' ' -f4,6,8 temp >> QUAL_DUP_smoove
done