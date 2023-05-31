#! /bin/bash

(cat SRA_List.txt; echo) | while read line
do
	FIRST=$(echo $line | awk '{print $1}')
	sbatch -p long launch_prefetch.sh $FIRST  
done
