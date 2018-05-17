#!/bin/bash


counter=0
for i in ${chr[@]}
do
  shapeit -B ${prefix[$counter]}_${data}_filtered \
        -M genetic_map_chr${chr[$counter]}_combined_b37.txt \
        -O ${prefix[$counter]}_${data}_filtered_phased \
        --thread $cpus
        
  shapeit -convert \
        --input-haps ${prefix[$counter]}_${data}_filtered_phased \
        --output-vcf ${prefix[$counter]}_${data}_filtered_phased.vcf
        
  counter=$counter+1
done 


. ./minimac3_imputation.sh

