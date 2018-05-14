#!/bin/bash

# Create files for imputation (Imputation preparation)
plink --bfile $data --freq --out $data
perl HRC-1000G-check-bim-v4.2.pl -b $data.bim -f $data.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

plink --bfile $data --exclude Exclude-$data-HRC.txt --make-bed --out TEMP1
plink --bfile TEMP1 --update-map Chromosome-$data-HRC.txt --update-chr --make-bed --out TEMP2
plink --bfile TEMP2 --update-map Position-$data-HRC.txt --make-bed --out TEMP3
plink --bfile TEMP3 --flip Strand-Flip-$data-HRC.txt --make-bed --out TEMP4
plink --bfile TEMP4 --reference-allele Force-Allele1-$data-HRC.txt --make-bed --out $data_updated
rm TEMP*


counter=0

for i in ${chr[@]}
do
  plink --bfile $data_updated \
      --reference-allele Force-Allele1-$data-HRC.txt \
      --make-bed \
      --mind 0.05 \
      --geno 0.05 \
      --chr ${chr[$counter]} \
      --out $data_chr${chr[$counter]}

 # plink --bfile $data_chr${chr[$counter]} \
 #     --mind 0.05 \
 #     --geno 0.05 \
 #     --make-bed \
 #     --out $data_chr${chr[$counter]}_filtered
  
  counter=$counter+1
done  









