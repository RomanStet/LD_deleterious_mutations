#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 20GB

VF_FOLDER=".../4_Variant_Calling/4_Variant_Filtration/"
GENOME=".../Reference_genome/Cr145.fasta"

module load gatk4/4.2.6.1

gatk VariantsToTable \
-R $GENOME \
-V ${VF_FOLDER}cr145_allScaffolds_snps_raw.vcf \
-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
-O ${VF_FOLDER}cr145_allScaffolds_snps_raw.table