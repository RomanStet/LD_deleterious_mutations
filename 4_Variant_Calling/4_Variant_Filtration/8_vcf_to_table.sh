#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 10GB

VF_FOLDER=".../4_Variant_Calling/4_Variant_Filtration/"
LD_ANALYSIS_FOLDER=".../5_LD_Analysis/"
GENOME=".../Reference_genome/Cr145.fasta"

module load gatk4/4.2.6.1

cd ${VF_FOLDER}

gatk VariantsToTable \
 -R ${GENOME} \
 -V cr145_Co_allScaffolds_snps_filter_1PASSED_DP_filter.vcf.recode.vcf \
 -F CHROM -F POS -F REF -F ALT -GF GT \
 -O ${LD_ANALYSIS_FOLDER}cr145_Co_SNPs.table

gatk VariantsToTable \
 -R ${GENOME} \
 -V cr145_Cg_allScaffolds_snps_filter_1PASSED_DP_filter.vcf.recode.vcf \
 -F CHROM -F POS -F REF -F ALT -GF GT \
 -O ${LD_ANALYSIS_FOLDER}cr145_Cg_SNPs.table
