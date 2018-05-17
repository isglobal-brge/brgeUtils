#!/bin/bash


if [[ $format == plink ]] || [[ $format == PLINK ]]
then
  counter=0
  for i in ${chr[@]}
  do
    plink --bfile ${data} \
        --make-bed \
        --mind 0.05 \
        --geno 0.05 \
        --chr ${chr[$counter]} \
        --out ${prefix[$counter]}_${data}_filtered
    counter=$counter+1
  done  
elif [[ $format == vcf ]] || [[ $format == VCF ]]
then
  counter=0
  for i in ${chr[@]}
  do
    plink --vcf ${data}.vcf.gz \
        --make-bed \
        --mind 0.05 \
        --geno 0.05 \
        --chr ${chr[$counter]} \
        --out ${prefix[$counter]}_${data}_filtered
    counter=$counter+1
  done 
fi 


. ./phasing_shapeit.sh








