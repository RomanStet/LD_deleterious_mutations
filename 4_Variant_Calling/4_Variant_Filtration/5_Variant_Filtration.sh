#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 20GB

VF_FOLDER=".../4_Variant_Calling/4_Variant_Filtration/"
GENOME=".../Reference_genome/Cr145.fasta"

module load gatk4/4.2.6.1

grep -vc "^#" ${VF_FOLDER}cr145_allScaffolds_snps_raw.vcf

gatk VariantFiltration \
   -R ${GENOME} \
   -V ${VF_FOLDER}cr145_allScaffolds_snps_raw.vcf \
   -O ${VF_FOLDER}cr145_allScaffolds_snps_filter_1.vcf \
    -filter "QD<2.0" --filter-name "QD" \
    -filter "QUAL<0" --filter-name "QUAL" \
    -filter "SOR>3.0" --filter-name "SOR" \
    -filter "FS>60.0" --filter-name "FS" \
    -filter "MQ<50.0" --filter-name "MQ" \
    -filter "MQRankSum<-5.0" --filter-name "MQRS" \
    -filter "ReadPosRankSum<-8.0" --filter-name "RPRS" \
    -filter "ReadPosRankSum>8.0" --filter-name "RPRS" \


grep -E '^#|PASS' ${VF_FOLDER}cr145_allScaffolds_snps_filter_1.vcf > ${VF_FOLDER}cr145_allScaffolds_snps_filter_1PASSED.vcf

grep -vc "^#" ${VF_FOLDER}cr145_allScaffolds_snps_filter_1PASSED.vcf
