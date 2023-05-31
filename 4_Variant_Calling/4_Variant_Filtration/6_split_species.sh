#! /bin/bash
#SBATCH --partition long
#SBATCH --cpus-per-task 5
#SBATCH --mem 50GB

VF_FOLDER=".../4_Variant_Calling/4_Variant_Filtration/"

cd ${VF_FOLDER}

grep '^##' -h cr145_allScaffolds_snps_filter_1PASSED.vcf > cr145_allScaffolds_snps_filter_1PASSED_headers
grep -v '^##' -h cr145_allScaffolds_snps_filter_1PASSED.vcf > temp
cut -f1-9,10-78,81-103,106,110-114,118-201 temp > temp_Cg
cut -f1-9,79-80,104-105,107-109,115-117,202-224 temp > temp_Co
rm temp

cat cr145_allScaffolds_snps_filter_1PASSED_headers temp_Cg > cr145_Cg_allScaffolds_snps_filter_1PASSED.vcf
cat cr145_allScaffolds_snps_filter_1PASSED_headers temp_Co > cr145_Co_allScaffolds_snps_filter_1PASSED.vcf 

rm temp_Cg
rm temp_Co