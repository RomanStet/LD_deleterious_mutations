#! /bin/bash
#SBATCH --cpus-per-task 5
#SBATCH --mem 20GB

module load gatk4/4.2.6.1
module load vcftools/0.1.16

VF_FOLDER=".../4_Variant_Calling/4_Variant_Filtration/"

cd ${VF_FOLDER}

VCF_CG="cr145_Cg_allScaffolds_snps_filter_1PASSED.vcf"
VCF_CO="cr145_Co_allScaffolds_snps_filter_1PASSED.vcf"

vcftools --vcf $VCF_CG --site-mean-depth --out Cg_mean_depth_site
vcftools --vcf $VCF_CO --site-mean-depth --out Co_mean_depth_site

# R script to visualise the distribution of mean depth and choose a threshold for the maximum mean depth per site

MIN_DEPTH=5
MAX_DEPTH_MEAN_CG=71.65
MAX_DEPTH_MEAN_CO=45

vcftools --vcf $VCF_CG \
--max-meanDP $MAX_DEPTH_MEAN_CG \
--minDP $MIN_DEPTH --max-missing 0.5 --min-alleles 2 --max-alleles 2 --maf 0.002 --recode --out cr145_Cg_allScaffolds_snps_filter_1PASSED_DP_filter.vcf

vcftools --vcf $VCF_CO \
--max-meanDP $MAX_DEPTH_MEAN_CO \
--minDP $MIN_DEPTH --max-missing 0.5 --min-alleles 2 --max-alleles 2 --maf 0.01 --recode --out cr145_Co_allScaffolds_snps_filter_1PASSED_DP_filter.vcf