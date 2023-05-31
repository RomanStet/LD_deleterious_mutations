#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 10GB

VF_FOLDER=".../4_Variant_Calling/4_Variant_Filtration/"

cd ${VF_FOLDER}

grep '^#' -h cr145_SCF_1_snps_raw.vcf > cr145_allScaffolds_snps_raw.vcf
grep -v '^#' -h cr145_SCF_1_snps_raw.vcf cr145_SCF_2_snps_raw.vcf cr145_SCF_3_snps_raw.vcf \
             cr145_SCF_4_snps_raw.vcf cr145_SCF_5_snps_raw.vcf cr145_SCF_6_snps_raw.vcf \
             cr145_SCF_7_snps_raw.vcf cr145_SCF_8_snps_raw.vcf >> cr145_allScaffolds_snps_raw.vcf

for i in {1..8..1}
do
	gzip cr145_SCF_${i}_snps_raw.vcf
done