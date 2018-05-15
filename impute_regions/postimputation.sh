#!/bin/bash

counter=0

for i in ${chr[@]}
do
  
  zcat ${chr[$counter]}_${data}_imputed.dose.vcf.gz | bgzip -c > ${chr[$counter]}_${data}_imputed_bgzip.vcf.gz && tabix ${chr[$counter]}_${data}_imputed_bgzip.vcf.gz
  
  bcftools view -h ${chr[$counter]}_${data}_imputed_bgzip.vcf.gz > hdr.txt
  sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' hdr.txt
  
  bcftools reheader -h hdr.txt --output TEMP0.vcf.gz ${chr[$counter]}_${data}_imputed_bgzip.vcf.gz
  tabix -p vcf TEMP0.vcf.gz
  
  bcftools view -r ${chr[$counter]}:${start[$counter]}-${end[$counter]} --output-type z --output-file TEMP1.vcf.gz TEMP0.vcf.gz # THIS STEP WON'T BE NECESSARY IF IMPUTING JUST THE INVERSION REGION
  tabix -p vcf TEMP1.vcf.gz
  
  bcftools annotate --remove ID --set-id +'%CHROM:%POS:%REF:%ALT' --threads $cpus --output-type z --output TEMP2.vcf.gz TEMP1.vcf.gz
  tabix -p vcf TEMP2.vcf.gz
  
  bcftools annotate --annotations All_new.vcf.gz --columns ID --threads $cpus --output-type z --output ${chr[$counter]}_${data}_imputed_final.vcf.gz TEMP2.vcf.gz
  tabix -p vcf ${chr[$counter]}_${data}_imputed_final.vcf.gz
  
  counter=$counter+1
done

rm TEMP*
rm hdr.txt

echo -e "\n\nPostimputation finished\n\n"

###################### CHANGE NAME OF THE VERY LAST FILE