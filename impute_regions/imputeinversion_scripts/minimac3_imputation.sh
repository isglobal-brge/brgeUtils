#!/bin/bash

mkdir pimputed_files

counter=0
for i in ${chr[@]}
do
  mkdir pimputed_files/${prefix[$counter]}
  Minimac3-omp --refHaps /home/isglobal.lan/itolosana/homews/reference_panels/ALL.chr${i}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
          --haps phased_files/${i}/${i}_${data}_filtered_phased.vcf \
          --rsid \
          --format GT,DS,GP \
          --chr $i \
          --start ${start[$counter]} \
          --end ${end[$counter]} \
          --prefix pimputed_files/${prefix[$counter]}/${prefix[$counter]}_${data}_imputed \
          --cpus $cpus
        
  counter=$counter+1
done  


. ./postimputation.sh



#### PONER OPCION DE IMPUTAR EL CROMOSOMA ENTERO