#!/bin/bash


counter=0
for i in ${chr[@]}
do
  
  zcat ${prefix[$counter]}_${data}_imputed.dose.vcf.gz | bgzip -c > ${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz && tabix ${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz
  
  bcftools view -h ${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz > hdr.txt
  sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' hdr.txt
  
  bcftools reheader -h hdr.txt --output TEMP1.vcf.gz ${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz 
  tabix -p vcf TEMP1.vcf.gz
  
  bcftools annotate --remove ID --set-id +'%CHROM:%POS:%REF:%ALT' --threads $cpus --output-type z --output TEMP2.vcf.gz TEMP1.vcf.gz
  tabix -p vcf TEMP2.vcf.gz
  
  bcftools annotate --annotations All_new.vcf.gz --columns ID --threads $cpus --output-type z --output ${prefix[$counter]}_${data}_imputed_final.vcf.gz TEMP2.vcf.gz
  tabix -p vcf ${prefix[$counter]}_${data}_imputed_final.vcf.gz
  
  counter=$counter+1
done

rm TEMP*
rm hdr.txt




