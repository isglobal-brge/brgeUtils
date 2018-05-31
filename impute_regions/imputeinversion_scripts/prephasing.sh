#!/bin/bash


#if [[ $format == vcf ]] || [[ $format == VCF ]]
#then
#  plink --vcf ${data}.vcf --biallelic-only strict --double-id --out ${data}
#fi
#vcftools --vcf ${data}.vcf --plink --out ${data}
#plink --file ${data} --make-bed --out ${data}


# Create files for imputation (Imputation preparation)
mkdir prephasing_files

plink --bfile $data --freq --out $data
perl HRC-1000G-check-bim.pl -b $data.bim -f $data.frq -r /home/isglobal.lan/itolosana/homews/reference_panels/1000GP_Phase3_combined.legend -g

plink --bfile $data --exclude Exclude-$data-1000G.txt --make-bed --out TEMP1
plink --bfile TEMP1 --update-map Chromosome-$data-1000G.txt --update-chr --make-bed --out TEMP2
plink --bfile TEMP2 --update-map Position-$data-1000G.txt --make-bed --out TEMP3
plink --bfile TEMP3 --flip Strand-Flip-$data-1000G.txt --make-bed --out TEMP4
plink --bfile TEMP4 --reference-allele Force-Allele1-$data-1000G.txt --make-bed --out prephasing_files/${data}_updated
rm TEMP*


if [[ $format == plink ]] || [[ $format == PLINK ]]
then
  for i in ${unique_chr[@]}
  do
    mkdir prephasing_files/${i}
    plink --bfile prephasing_files/${data}_updated \
        --reference-allele Force-Allele1-$data-1000G.txt \
        --make-bed \
        --mind 0.05 \
        --geno 0.05 \
        --chr $i \
        --out prephasing_files/${i}/${i}_${data}_filtered
  done  
elif [[ $format == vcf ]] || [[ $format == VCF ]]
then
  counter=0
  for i in ${unique_chr[@]}
  do
    mkdir prephasing_files/${i}
    plink --vcf ${data}.vcf.gz \
        --make-bed \
        --mind 0.05 \
        --geno 0.05 \
        --chr $i \
        --out prephasing_files/${i}/${i}_${data}_filtered
    counter=$counter+1
  done 
fi 

mv *-1000G.txt prephasing_files/


. ./phasing_shapeit.sh







