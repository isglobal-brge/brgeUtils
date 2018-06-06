#!/bin/bash

mkdir postimputation_files
mkdir ${data}_imputed_files

counter=0
for i in ${chr[@]}
do
  if [[ $i == X ]]
  then
    zcat pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed.dose.vcf.gz | bgzip -c > postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed_bgzip.vcf.gz && tabix postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed_bgzip.vcf.gz
    
    bcftools view -h postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed_bgzip.vcf.gz > postimputation_files/hdrmales.txt
    sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' postimputation_files/hdrmales.txt
  
    bcftools reheader -h postimputation_files/hdrmales.txt --output postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_hdr_imputed_bgzip.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed_bgzip.vcf.gz
    tabix -p vcf postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_hdr_imputed_bgzip.vcf.gz
    
    zcat pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed.dose.vcf.gz | bgzip -c > postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed_bgzip.vcf.gz && tabix postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed_bgzip.vcf.gz
    
    bcftools view -h postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed_bgzip.vcf.gz > postimputation_files/hdrfemales.txt
    sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' postimputation_files/hdrfemales.txt
  
    bcftools reheader -h postimputation_files/hdrfemales.txt --output postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_hdr_imputed_bgzip.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed_bgzip.vcf.gz
    tabix -p vcf postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_hdr_imputed_bgzip.vcf.gz
        
    bcftools merge --force-samples --output-type z --output postimputation_files/TEMP1.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_hdr_imputed_bgzip.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_hdr_imputed_bgzip.vcf.gz
    tabix -p vcf postimputation_files/TEMP1.vcf.gz
    
  else
    zcat pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed.dose.vcf.gz | bgzip -c > postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz && tabix postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz    
    
    bcftools view -h postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz > postimputation_files/hdr.txt
    sed -i '/^##FORMAT=<ID=GT/i ##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\n##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' postimputation_files/hdr.txt
  
    bcftools reheader -h postimputation_files/hdr.txt --output postimputation_files/TEMP1.vcf.gz postimputation_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_bgzip.vcf.gz
    tabix -p vcf postimputation_files/TEMP1.vcf.gz
  fi

  bcftools annotate --remove ID --set-id +'%CHROM:%POS:%REF:%ALT' --threads $cpus --output-type z --output postimputation_files/TEMP2.vcf.gz postimputation_files/TEMP1.vcf.gz
  tabix -p vcf postimputation_files/TEMP2.vcf.gz
  
  bcftools annotate --annotations /home/isglobal.lan/itolosana/homews/reference_panels/All_new.vcf.gz --columns ID --threads $cpus --output-type z --output ${data}_imputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_final.vcf.gz postimputation_files/TEMP2.vcf.gz
  tabix -p vcf ${data}_imputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed_final.vcf.gz
  
  counter=$counter+1
done

rm postimputation_files/TEMP*
rm postimputation_files/hdr*

if [[ -z "$keep_files" ]] || [[ $keep_files != Yes ]] && [[ $keep_files != YES ]] && [[ $keep_files != yes ]] # By default, the intermediate files will be deleted
then
  rm -rf pimputed_files
  rm -rf postimputation_files
fi

#bcftools view -r ${chr[$counter]}:${start[$counter]}-${end[$counter]} --output-type z --output-file TEMP1.vcf.gz TEMP0.vcf.gz # THIS STEP WON'T BE NECESSARY IF IMPUTING JUST THE INVERSION REGION
#tabix -p vcf TEMP1.vcf.gz
