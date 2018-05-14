#!/bin/bash

for i in ${chr[@]}
do
  
  zcat $data_chr${chr[$counter]}_imputed.dose.vcf.gz | bgzip -c > $data_chr${chr[$counter]}_imputed_bgzip.vcf.gz && tabix $data_chr${chr[$counter]}_imputed_bgzip.vcf.gz
  
  bcftools view -h $data_chr${chr[$counter]}_imputed_bgzip.vcf.gz > hdr.txt
  sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' hdr.txt
  
  bcftools reheader -h hdr.txt --output TEMP0.vcf.gz $data_chr${chr[$counter]}_imputed_bgzip.vcf.gz
  tabix -p vcf TEMP0.vcf.gz
  
  bcftools view -r ${chr[$counter]}:${start[$counter]}-${end[$counter]} --output-type z --output-file TEMP1.vcf.gz TEMP0.vcf.gz
  tabix -p vcf TEMP1.vcf.gz
  
  bcftools annotate --remove ID --set-id +'%CHROM:%POS:%REF:%ALT' --output-type z --output TEMP2.vcf.gz TEMP1.vcf.gz
  tabix -p vcf TEMP2.vcf.gz
  
  bcftools annotate --annotations All_new.vcf.gz --columns ID --output-type z --output $data_final_minimac_chr${chr[$counter]}.vcf.gz TEMP2.vcf.gz
  tabix -p vcf $data_final_minimac_chr${chr[$counter]}.vcf.gz
  
  rm TEMP*
  rm hdr.txt
  counter=$counter+1
done

###################### CHANGE NAME OF THE VERY LAST FILE