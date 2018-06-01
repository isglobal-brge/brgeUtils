#!/bin/bash

mkdir pimputed_files

counter=0
for i in ${chr[@]}
do
  mkdir pimputed_files/${prefix[$counter]}
  if [[ $i == X ]]
  then
    sed 's/^23/X/' phased_files/${i}/${i}_${data}_filtered_phased.vcf  > phased_files/${i}/${i}_${data}_filtered_phased_X.vcf
    awk '$5 == 2 { print $2 > "females.txt"}' ${data}.fam
    awk '$5 == 1 { print $2 > "males.txt"}' ${data}.fam
    
    vcftools --vcf phased_files/${i}/${i}_${data}_filtered_phased_X.vcf \
            --keep males.txt\
            --recode \
            --out phased_files/${i}/${i}_${data}_males_filtered_phased_X

    vcftools --vcf phased_files/${i}/${i}_${data}_filtered_phased_X.vcf \
            --keep females.txt\
            --recode \
            --out phased_files/${i}/${i}_${data}_females_filtered_phased_X
    
    Minimac3-omp --refHaps /home/isglobal.lan/itolosana/homews/reference_panels/ALL.chrX.Non.Pseudo.Auto.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
            --haps phased_files/${i}/${i}_${data}_males_filtered_phased_X.recode.vcf \
            --rsid \
            --format GT,DS,GP \
            --chr X \
            --start ${start[$counter]} \
            --end ${end[$counter]} \
            --prefix pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_males_imputed \
            --cpus $cpus

    Minimac3-omp --refHaps /home/isglobal.lan/itolosana/homews/reference_panels/ALL.chrX.Non.Pseudo.Auto.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
            --haps phased_files/${i}/${i}_${data}_females_filtered_phased_X.recode.vcf \
            --rsid \
            --format GT,DS,GP \
            --chr X \
            --start ${start[$counter]} \
            --end ${end[$counter]} \
            --prefix pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_females_imputed \
            --cpus $cpus
          
  else
    Minimac3-omp --refHaps /home/isglobal.lan/itolosana/homews/reference_panels/ALL.chr${i}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
            --haps phased_files/${i}/${i}_${data}_filtered_phased.vcf \
            --rsid \
            --format GT,DS,GP \
            --chr $i \
            --start ${start[$counter]} \
            --end ${end[$counter]} \
            --prefix pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed \
            --cpus $cpus
  fi    
  counter=$counter+1
done  

if [[ -z "$keep_files" ]] || [[ $keep_files != Yes ]] && [[ $keep_files != YES ]] && [[ $keep_files != yes ]] # By default, the intermediate files will be deleted
then
  rm -rf phased_files
fi

. ./postimputation.sh



#### PONER OPCION DE IMPUTAR EL CROMOSOMA ENTERO